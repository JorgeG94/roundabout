!! Shear-driven interior turbulence (kappa-shear) for the ocean core.
module rdb_ocean_kappa_shear
   !! Prognostic interior shear-mixing closure: shear-driven turbulence
   !! is modelled as a diffusivity field kappa(z) and a TKE field Q(z)
   !! at layer interfaces, coupled through two steady-state vertical
   !! diffusion-reaction equations solved per column, iteratively to
   !! convergence, with internal adaptive time-substepping as the
   !! column re-stratifies within one model step.  Unlike algebraic
   !! Richardson-number schemes (PP81, LMD94) the diffusivity diffuses
   !! in z with a stratification/rotation/boundary-limited decay length,
   !! so resolved shear layers entrain at the right rate even when the
   !! Richardson number is marginal; the closure is self-limiting and
   !! relaxes the column toward Ri >~ Ri_c.
   !!
   !! Reference: Jackson, Hallberg & Legg (2008), "A Parameterization
   !! of Shear-Driven Turbulence for Ocean Climate Models", J. Phys.
   !! Oceanogr. 38, 1033-1053.  Knob table: `docs/generated_nml_knobs.md`.
   !!
   !! Place in the stack: an INTERIOR closure.  It coexists with the
   !! surface boundary-layer schemes (KPP or EPBL) and with
   !! PP81/background; its kappa is ADDED to the other interior
   !! diffusivities (`kt += kd_int`) and its viscosity (`prandtl_turb *
   !! kd_int`) added to `kv`.  Not mutually exclusive with anything.
   !! Runs at thermo cadence (`kappa_shear_compute`); the merge into
   !! `vmix%kv` / `vmix%kt` runs every RK2 stage
   !! (`kappa_shear_merge_into_kv_kt`).
   !!
   !! Interface convention (same as `rdb_ocean_vmix` / EPBL):
   !! `kd_int(:,:,K)` lives at the bottom interface of layer K; global
   !! `kd_int(:,:,1)` is the bed and `kd_int(:,:,nz+1)` the free
   !! surface, both forced to exactly 0.  Layers are bottom-up: k=1 bed,
   !! k=nz surface.  The column solver itself runs surface-down (local
   !! k=1 = surface); the gather/scatter loops carry the index flip
   !! (global K = nz+2 - K_local).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, H_DIV_EPS
#endif
   use rdb_massless, only: massless_build_maps, massless_merge_fields, &
                           massless_interp_back
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_specvol_derivs
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ocean_kappa_shear_t
   public :: kappa_shear_compute
   public :: kappa_shear_merge_into_kv_kt
   public :: kappa_shear_vertex_scatter
   public :: ks_gather_corner

   ! Local-column array dimensions.  Layers run 1..NZL, interfaces
   ! 1..NZL+1 (surface-down); arrays are sized to the interface count.
   integer, parameter :: NZL = NZ_STACK_MAX
      !! Maximum number of layers a single column kernel can solve.
   integer, parameter :: NZLI = NZ_STACK_MAX + 1
      !! Interface array dimension (= NZL + 1).

   ! Corner-gather regularisers (vertex form).  Roles mirror MOM6's
   ! `H_tiny` / `H_subroundoff` / mask-sum `1e-36`:
   real(wp), parameter :: H_TINY_CORNER = 0.5_wp*H_DIV_EPS
      !! Sub-roundoff thickness added to the 2-point thickness-weight
      !! denominator of the corner u/v average (pure 1/0 armour).
   real(wp), parameter :: MASK_SUM_EPS = 1.0e-36_wp
      !! Armour added to a MASK sum (a nondimensional wet-cell count),
      !! so an all-land denominator gives 0/1e-36 = 0 rather than 0/0.

   ! Gather floor / merge threshold for vanished layers.  D4: this is
   ! the gather-time enforcement of the dynamic-vanish role, so it now
   ! uses the constant of record `H_VANISHED` (rdb_constants).  Bitwise
   ! no-op (1.5e-4 == the old private H_FLOOR).  When `massless_merge`
   ! is on, the blunt `max(h, H_VANISHED)` clamp is REPLACED by the
   ! merge (rdb_massless); when off, the clamp path is unchanged.

   type :: ocean_kappa_shear_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.

      ! ---- Scheme selection + master switch ----
      logical :: enable = .false.
         !! Master switch.  Default off — existing namelists and tests
         !! stay bit-identical.  Requires `vmix%use_closure` +
         !! thermodynamics (validated at configure).

      ! ---- JHL08 knobs (defaults = paper / OM4 production) ----
      real(wp) :: ri_crit = 0.25_wp
         !! Critical Richardson number (MOM6 RINO_CRIT).
      real(wp) :: shearmix_rate = 0.089_wp
         !! Source-rate coefficient (SHEARMIX_RATE).
      real(wp) :: fri_curvature = -0.97_wp
         !! Ri-function curvature (FRI_CURVATURE).
      real(wp) :: c_n = 0.24_wp
         !! TKE decay vs N (TKE_N_DECAY_CONST).
      real(wp) :: c_s = 0.14_wp
         !! TKE decay vs shear (TKE_SHEAR_DECAY_CONST).
      real(wp) :: lambda = 0.82_wp
         !! Buoyancy length-scale coefficient (KAPPA_BUOY_SCALE_COEF).
      real(wp) :: lz_rescale = 1.0_wp
         !! Boundary-distance length-scale rescale (LZ_RESCALE).
      real(wp) :: kappa_0 = 1.0e-7_wp
         !! Background diffusivity (m^2/s); also the pre-step kappa.
      real(wp) :: kappa_seed = 1.0_wp
         !! Iteration seed diffusivity (m^2/s).
      real(wp) :: kappa_trunc = 1.0e-9_wp
         !! Diffusivity below this -> 0 (m^2/s).
      real(wp) :: tke_bg = 0.0_wp
         !! Background TKE (m^2/s^2); Q is a denominator, floored.
      real(wp) :: tol_err = 0.1_wp
         !! Picard convergence tolerance (KAPPA_SHEAR_TOL_ERR).
      integer :: max_inner_it = 50
         !! Inner Picard iteration cap (MAX_RINO_IT).
      integer :: max_substep_it = 13
         !! Outer adaptive substep cap (MAX_KAPPA_SHEAR_IT).
      real(wp) :: src_max_chg = 10.0_wp
         !! Adaptive-dt source-change tolerance band.
      real(wp) :: prandtl_turb = 1.0_wp
         !! Kv = prandtl_turb * Kd into the momentum solve.
      real(wp) :: vel_underflow = 0.0_wp
         !! Velocity snap-to-zero magnitude (m/s) in the projection.
      logical :: massless_merge = .false.
         !! D4: when on, a column carrying vanished (< H_VANISHED) layers
         !! is merged onto its massive sub-grid (rdb_massless), solved on
         !! `nzc <= nz` layers, and the kappa/TKE interface fields are
         !! interpolated back — replacing the blunt `max(h, H_VANISHED)`
         !! gather floor.  Default off (bit-identical).  A per-column
         !! `any(h < H_VANISHED)` precheck makes healthy columns bypass
         !! the merge machinery entirely, so knob-on stays bit-identical
         !! on healthy envelopes (invariant I1).

      ! ---- Vertex (corner) form — MOM6 VERTEX_SHEAR ----
      logical :: at_vertex = .false.
         !! Solve the JHL08 columns at C-grid CORNERS (vorticity points)
         !! instead of tracer points, then average the corner Kd back to
         !! tracer points (MOM6 `VERTEX_SHEAR`; the OM5-class production
         !! setting).  The corner column sees the native face velocities
         !! without the u_h/v_h centre average, so resolved shear is not
         !! damped before the solve.  Default off — bit-identical.
         !! Kd: corner solve -> Pass-C average -> `kd_int` at tracer
         !! points.  Kv: routed corner->face (MOM6 `Kv_shear_Bu`
         !! consumed in vertvisc) — `kd_corner` feeds
         !! `vdiff_apply_momentum`'s `kv_corner_source` seam scaled by
         !! `prandtl_turb`, and the cell-centred kv merge is suppressed
         !! (no corner->centre->face smoothing, no double-count).
         !! Corner TKE is not carried (`tke_int` is zeroed in vertex
         !! mode).
      logical :: vertex_geometric_mean = .false.
         !! Corner->centre averaging: geometric mean of the 4 corner Kd
         !! (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN`) instead of the plain
         !! arithmetic mean.  A geometric mean is 0 if ANY corner is 0 —
         !! pair with `vertex_geomean_kdmin` (see below).
      real(wp) :: vertex_geomean_kdmin = 0.0_wp
         !! Floor (m^2/s) applied to each corner Kd BEFORE the geometric
         !! mean (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN_KDMIN`; inert unless
         !! `vertex_geometric_mean`).  With 0 the geometric mean hard-
         !! zeros Kd at every shear-zone edge; OM5 configs use 1e-9.

      ! ---- EOS hookup (shared handle from the eos slot) ----
      ! Value copy of the flat-POD `eos_t` set at configure
      ! from `ocean_state%eos` — the buoyancy derivatives use the
      ! SAME EOS the dyn-core runs.  Maps onto the device with the
      ! parent for free; one source of truth (no drifting copies).
      type(eos_t) :: eos
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).

      ! ---- Persistent fields ----
      real(wp), allocatable :: f_centre(:, :)
         !! |f| at cell centres (1/s); filled by `set_f_centre`.
      real(wp), allocatable :: kd_int(:, :, :)
         !! Kappa-shear diffusivity at interfaces (m^2/s), (nx,ny,nz+1),
         !! global bottom-up: zero at bed (K=1) and surface (K=nz+1).
      real(wp), allocatable :: tke_int(:, :, :)
         !! Time-mean TKE at interfaces (m^2/s^2), same shape/convention
         !! — diagnostic (currently filled to 0; reserved for the TKE
         !! budget diag).  Vertex mode zeroes it (corner TKE not carried).

      ! ---- Vertex-mode persistent fields (allocated by `init_vertex`,
      ! only when `at_vertex` — the corner carrier is nz+1 full planes) ----
      real(wp), allocatable :: f_corner(:, :)
         !! SIGNED Coriolis f at C-grid corners (1/s), `(nx+1,ny+1)`;
         !! corner (i,j) is the SW corner of cell (i,j).  The kernel
         !! squares it (MOM6 vertex form takes f^2 straight at the
         !! corner, no 4-point average).  Filled by `set_f_corner`
         !! (beta-plane) or `fill_coriolis_corner` at configure.
      real(wp), allocatable :: kd_corner(:, :, :)
         !! Corner diffusivity at interfaces (m^2/s), `(nx+1,ny+1,nz+1)`,
         !! global bottom-up.  The Pass-B -> Pass-C carrier (MOM6
         !! `kappa_vertex`): the corner solve writes it, the scatter
         !! kernel averages it to tracer points, and the momentum vdiff
         !! reads it as the corner Kv source (`kv_corner_source`,
         !! scaled by `prandtl_turb` — the corner->face viscosity seam,
         !! MOM6 `Kv_shear_Bu`).  A REQUIRED snapshot —
         !! fusing solve+scatter would be a read-neighbour/write-own
         !! `do concurrent` race.  Ring corners (ic=1, ic=nx+1, jc=1,
         !! jc=ny+1) are never solved and stay exactly 0.
   contains
      procedure, non_overridable :: init => ocean_kappa_shear_init
      procedure, non_overridable :: init_vertex => ocean_kappa_shear_init_vertex
      procedure, non_overridable :: destroy => ocean_kappa_shear_destroy
      procedure, non_overridable :: enter_data => ocean_kappa_shear_enter_data
      procedure, non_overridable :: exit_data => ocean_kappa_shear_exit_data
      procedure, non_overridable :: set_f_centre => ocean_kappa_shear_set_f_centre
      procedure, non_overridable :: set_f_corner => ocean_kappa_shear_set_f_corner
      procedure, non_overridable :: bytes => ocean_kappa_shear_bytes
   end type ocean_kappa_shear_t

contains

   ! =================================================================
   ! Lifecycle
   ! =================================================================

   subroutine ocean_kappa_shear_init(this, grid, nz_ml)
      !! Allocate the persistent fields.  Always allocates (configure
      !! runs after init, so `enable` is not known yet); the off-cost
      !! is the f_centre + two interface fields.
      class(ocean_kappa_shear_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      allocate (this%f_centre(nx, ny), source=0.0_wp)
      allocate (this%kd_int(nx, ny, nz + 1), source=0.0_wp)
      allocate (this%tke_int(nx, ny, nz + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_kappa_shear_init

   subroutine ocean_kappa_shear_init_vertex(this, grid, nz_ml)
      !! Allocate the vertex-mode corner fields and set `at_vertex`.
      !! Called at CONFIGURE time (after `init`, before `enter_data`) —
      !! deliberately NOT from `init`, so the (nx+1,ny+1,nz+1) corner
      !! carrier is only ever allocated when the vertex form is actually
      !! selected (~327 MB at 1000x800x50).
      class(ocean_kappa_shear_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_ml
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      this%at_vertex = .true.
      if (.not. allocated(this%f_corner)) then
         allocate (this%f_corner(nx + 1, ny + 1), source=0.0_wp)
      end if
      if (.not. allocated(this%kd_corner)) then
         allocate (this%kd_corner(nx + 1, ny + 1, nz_ml + 1), source=0.0_wp)
      end if
   end subroutine ocean_kappa_shear_init_vertex

   subroutine ocean_kappa_shear_destroy(this)
      class(ocean_kappa_shear_t), intent(inout) :: this
      this%is_init = .false.
      ! Reset the vertex arming with the fields it gates — a destroy ->
      ! re-init resurrection must not leave `at_vertex` pointing at
      ! deallocated corner fields.
      this%at_vertex = .false.
      if (allocated(this%f_centre)) deallocate (this%f_centre)
      if (allocated(this%kd_int)) deallocate (this%kd_int)
      if (allocated(this%tke_int)) deallocate (this%tke_int)
      if (allocated(this%f_corner)) deallocate (this%f_corner)
      if (allocated(this%kd_corner)) deallocate (this%kd_corner)
   end subroutine ocean_kappa_shear_destroy

   subroutine ocean_kappa_shear_enter_data(this)
      class(ocean_kappa_shear_t), intent(inout) :: this
      select type (this)
      type is (ocean_kappa_shear_t)
         call ocean_kappa_shear_enter_data_impl(this)
      end select
   end subroutine ocean_kappa_shear_enter_data

   subroutine ocean_kappa_shear_enter_data_impl(this)
      type(ocean_kappa_shear_t), intent(inout) :: this
      if (allocated(this%f_centre)) then
         !$acc enter data copyin(this%f_centre)
      end if
      if (allocated(this%kd_int)) then
         !$acc enter data copyin(this%kd_int)
      end if
      if (allocated(this%tke_int)) then
         !$acc enter data copyin(this%tke_int)
      end if
      if (allocated(this%f_corner)) then
         !$acc enter data copyin(this%f_corner)
      end if
      if (allocated(this%kd_corner)) then
         !$acc enter data copyin(this%kd_corner)
      end if
   end subroutine ocean_kappa_shear_enter_data_impl

   subroutine ocean_kappa_shear_exit_data(this)
      class(ocean_kappa_shear_t), intent(inout) :: this
      select type (this)
      type is (ocean_kappa_shear_t)
         call ocean_kappa_shear_exit_data_impl(this)
      end select
   end subroutine ocean_kappa_shear_exit_data

   subroutine ocean_kappa_shear_exit_data_impl(this)
      type(ocean_kappa_shear_t), intent(inout) :: this
      if (allocated(this%kd_corner)) then
         !$acc exit data delete(this%kd_corner)
      end if
      if (allocated(this%f_corner)) then
         !$acc exit data delete(this%f_corner)
      end if
      if (allocated(this%tke_int)) then
         !$acc exit data delete(this%tke_int)
      end if
      if (allocated(this%kd_int)) then
         !$acc exit data delete(this%kd_int)
      end if
      if (allocated(this%f_centre)) then
         !$acc exit data delete(this%f_centre)
      end if
   end subroutine ocean_kappa_shear_exit_data_impl

   subroutine ocean_kappa_shear_set_f_centre(this, grid, f_0, beta, y_ref)
      !! Fill `f_centre` with the beta-plane Coriolis magnitude at cell
      !! centres: |f_0 + beta*(y - y_ref)|.  Mirrors EPBL's
      !! `set_f_centre`.  Call after `init`, before `enter_data`.
      class(ocean_kappa_shear_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_0, beta, y_ref
      integer :: i, j, ng
      real(wp) :: y

      ng = grid%nghost
      do j = 1, size(this%f_centre, 2)
         y = (real(j + grid%j_offset_global - ng, wp) - 0.5_wp)*grid%dy
         do i = 1, size(this%f_centre, 1)
            this%f_centre(i, j) = abs(f_0 + beta*(y - y_ref))
         end do
      end do
   end subroutine ocean_kappa_shear_set_f_centre

   subroutine ocean_kappa_shear_set_f_corner(this, grid, f_0, beta, y_ref)
      !! Fill `f_corner` with the SIGNED beta-plane Coriolis at C-grid
      !! corners: f_0 + beta*(y - y_ref), corner row j at
      !! y = (j-1-nghost)*dy (half a cell below centre row j — corner
      !! (i,j) is the SW corner of cell (i,j)).  Bit-identical to
      !! `metrics_fill_coriolis`'s beta-plane corner fill.  Call after
      !! `init_vertex`, before `enter_data`.
      class(ocean_kappa_shear_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_0, beta, y_ref
      integer :: i, j, ng
      real(wp) :: y

      ng = grid%nghost
      do j = 1, size(this%f_corner, 2)
         y = real(j + grid%j_offset_global - 1 - ng, wp)*grid%dy
         do i = 1, size(this%f_corner, 1)
            this%f_corner(i, j) = f_0 + beta*(y - y_ref)
         end do
      end do
   end subroutine ocean_kappa_shear_set_f_corner

   ! =================================================================
   ! Merge into the vmix interface fields
   ! =================================================================

   pure subroutine kappa_shear_merge_into_kv_kt(this, nx, ny, nzp1, kv, kt)
      !! Fold the kappa-shear diffusivity into the vmix interface
      !! fields, ADDITIVELY (MOM6 interior-diffusivity semantics).
      !! Called EVERY stage (the interior closure rewrites kv/kt each
      !! stage; `kd_int` itself refreshes at thermo cadence).  Interior
      !! interfaces only — K=1 (bed) and K=nz+1 (surface) stay zero in
      !! both source and target.
      !!
      !! Vertex mode (`at_vertex`): the CELL-CENTRED kv add is
      !! suppressed (kv_fac = 0) — the shear viscosity instead reaches
      !! the momentum solve as the corner field `kd_corner`, routed
      !! corner→face inside `vdiff_apply_momentum`
      !! (`kv_corner_source`), never via a tracer point.  Adding it
      !! here too would double-count (the reference implementation
      !! zeroes the tracer-point `Kv_shear` when the vertex form is
      !! selected, for exactly this reason).  `kt` still receives the
      !! Pass-C tracer-point `kd_int` in both modes.
      type(ocean_kappa_shear_t), intent(in) :: this
      integer, intent(in) :: nx, ny, nzp1
         !! Interface-field extents (explicit shape: assumed-shape
         !! dummies in a `do concurrent` kernel make NVHPC walk the
         !! descriptor with per-launch memcpys — this runs every stage).
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
         !! Momentum viscosity at interfaces; gets prandtl_turb*kd
         !! (column mode) or nothing (vertex mode — see above).
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
         !! Tracer diffusivity at interfaces; gets kd.

      integer :: i, j, k
      real(wp) :: kv_fac

      ! Host-hoisted scalar (bit-identical in column mode: kv_fac IS
      ! prandtl_turb, same multiply).
      kv_fac = this%prandtl_turb
      if (this%at_vertex) kv_fac = 0.0_wp

      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx)
         kt(i, j, k) = kt(i, j, k) + this%kd_int(i, j, k)
         kv(i, j, k) = kv(i, j, k) + kv_fac*this%kd_int(i, j, k)
      end do
   end subroutine kappa_shear_merge_into_kv_kt

   ! =================================================================
   ! Main compute (outer shim)
   ! =================================================================

   pure subroutine kappa_shear_compute(grid, this, ms, dt, wet_t, wet_u, wet_v)
      !! Run kappa-shear over the domain: fill `this%kd_int` (interface
      !! diffusivity).  Call at thermo cadence with the thermo dt.
      !! Outer shim: dereferences the tracer-registry hTr arrays on the
      !! host (the array-of-DT indirection blocks NVHPC device codegen),
      !! then forwards to the column kernel.  Velocities reach the
      !! kernel as the C-grid face arrays, face-averaged to tracer
      !! points inside (same source as the PP81 interior shear).
      !!
      !! Vertex mode (`at_vertex`): solves the columns at C-grid corners
      !! from the NATIVE face velocities (no centre average), then
      !! averages the corner Kd back to tracer points — two kernels with
      !! `kd_corner` as the required carrier.  Needs the metrics wet
      !! masks (halo-valid `wet_T` + face masks); fail-loud if absent.
      type(hgrid_t), intent(in) :: grid
      type(ocean_kappa_shear_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
      ! assumed-shape-ok: host shim only — forwarded to the flat kernels
      ! as explicit-shape dummies; thermo cadence.
      real(wp), intent(in), optional :: wet_t(:, :)
         !! Halo-valid tracer-cell wet mask (`metrics%wet_T`), (nx,ny).
      real(wp), intent(in), optional :: wet_u(:, :)  ! assumed-shape-ok: host shim; thermo cadence
         !! u-face open mask (`metrics%wet_u`), (nx+1,ny).
      real(wp), intent(in), optional :: wet_v(:, :)  ! assumed-shape-ok: host shim; thermo cadence
         !! v-face open mask (`metrics%wet_v`), (nx,ny+1).

      if (ms%idx_temperature <= 0 .or. ms%idx_salinity <= 0) return

      if (this%at_vertex) then
         if (.not. (present(wet_t) .and. present(wet_u) .and. present(wet_v))) then
            error stop "kappa_shear_compute: at_vertex requires the wet_t/wet_u/wet_v masks"
         end if
         call kappa_shear_vertex_kernel(this, grid%nx_total, grid%ny_total, &
                                        ms%nz_ml, dt, &
                                        ms%h_layer, ms%u_face_x_layer, &
                                        ms%v_face_y_layer, &
                                        ms%tracers(ms%idx_temperature)%hTr, &
                                        ms%tracers(ms%idx_salinity)%hTr, &
                                        wet_t, wet_u, wet_v)
         call kappa_shear_vertex_scatter(grid%nx_total, grid%ny_total, &
                                         ms%nz_ml + 1, &
                                         this%vertex_geometric_mean, &
                                         this%vertex_geomean_kdmin, &
                                         wet_t, this%kd_corner, &
                                         this%kd_int, this%tke_int)
      else
         call kappa_shear_column_kernel(grid, this, ms, &
                                        ms%tracers(ms%idx_temperature)%hTr, &
                                        ms%tracers(ms%idx_salinity)%hTr, dt)
      end if
   end subroutine kappa_shear_compute

   pure subroutine kappa_shear_column_kernel(grid, this, ms, hT, hS, dt)
      !! Per-column JHL08 solve.  One `do concurrent (j, i)` over owned
      !! cells; every column is solved serially in surface-down order.
      !! Per-thread work is fixed-size `local()` arrays (L1 layout —
      !! all live in GPU registers, no shared memory); the column solve
      !! bodies are the same-module pure `!$acc routine seq` helpers below.
      !!
      !! Pipeline:
      !!   gather (FLIP global bottom-up -> local surface-down): h
      !!          (floored), u,v face-averaged to centre, T,S = hTr/h.
      !!   D4 (massless_merge on + column has a vanished layer): build the
      !!          merge maps, fold vanished layers onto the massive
      !!          sub-grid (nzc<=nz), solve there, interp kappa/TKE back.
      !!          Identity columns (none vanished) bypass the merge.
      !!   precompute: grids, h_Int, boundary length scale, background
      !!          kappa_0 pre-step tridiagonal, frozen EOS buoyancy
      !!          derivatives, initial N^2/S^2.
      !!   outer:  adaptive substepping with predictor-corrector and the
      !!          Picard inner (kappa,Q) solve.
      !!   scatter (FLIP back): kappa_avg -> kd_int, 0 at bed/surface.
      type(hgrid_t), intent(in) :: grid
      type(ocean_kappa_shear_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      ! assumed-shape-ok: tracer registry outer-shim — caller host-dereferences
      ! ms%tracers(idx)%hTr before passing; size varies per tracer slot
      ! (see CLAUDE.md "outer-shim + flat-impl" pattern); thermo cadence.
      real(wp), intent(in) :: hT(:, :, :)
         !! Temperature tracer hTr (degC*m), host-dereferenced.
      real(wp), intent(in) :: hS(:, :, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
         !! Salinity tracer hTr (PSU*m), host-dereferenced.
      real(wp), intent(in) :: dt

      integer :: i, j, k, kg, nx, ny, nz, nzc
      real(wp) :: f2_val, hk, inv_h
      logical :: merge_on, do_merge
      ! gathered surface-down column inputs
      real(wp) :: h_sd(NZL), u_sd(NZL), v_sd(NZL), t_sd(NZL), s_sd(NZL)
      ! precomputed interface/layer grids
      real(wp) :: idz_s(NZL), idz_int_s(NZLI), hint_s(NZLI), il2_s(NZLI)
      real(wp) :: kappa_avg_sd(NZLI), tke_avg_sd(NZLI)
      ! D4 massless-merge per-column scratch (+9 arrays: 8 real, 1 int).
      ! kc/kf are the merge maps; hc/uc/vc/tc/sc the merged column;
      ! kappa_c/tke_c the merged-grid interface outputs before interp.
      real(wp) :: hc(NZL), uc(NZL), vc(NZL), tc(NZL), sc(NZL)
      real(wp) :: kappa_c(NZLI), tke_c(NZLI), kf(NZL + 1)
      integer :: kc(NZL + 1)

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      merge_on = this%massless_merge

      do concurrent(j=1:ny, i=1:nx) &
         local(k, kg, nzc, f2_val, hk, inv_h, do_merge, &
               h_sd, u_sd, v_sd, t_sd, s_sd, &
               idz_s, idz_int_s, hint_s, il2_s, &
               kappa_avg_sd, tke_avg_sd, &
               hc, uc, vc, tc, sc, kappa_c, tke_c, kc, kf)

         ! ---- Default outputs (overwritten for wet columns) ----
         do k = 1, nz + 1
            kappa_avg_sd(k) = 0.0_wp
            tke_avg_sd(k) = 0.0_wp
         end do

         if (ms%wet_mask(i, j) > 0.0_wp) then
            ! ---- Gather with the index flip (local k=1 = surface) ----
            ! global layer (kg = nz+1-k) -> local layer k.  Gather RAW
            ! thickness; the I1 precheck below decides floor vs merge.
            do k = 1, nz
               kg = nz + 1 - k
               h_sd(k) = ms%h_layer(i, j, kg)
               u_sd(k) = 0.5_wp*(ms%u_face_x_layer(i, j, kg) + &
                                 ms%u_face_x_layer(i + 1, j, kg))
               v_sd(k) = 0.5_wp*(ms%v_face_y_layer(i, j, kg) + &
                                 ms%v_face_y_layer(i, j + 1, kg))
            end do

            ! ---- I1 precheck (cheap, per column): merge only when the
            ! knob is on AND the column actually carries a vanished layer.
            ! Identity columns bypass the merge machinery entirely, so
            ! knob-on is bit-identical to knob-off on healthy envelopes.
            do_merge = .false.
            if (merge_on) then
               do k = 1, nz
                  if (h_sd(k) < H_VANISHED) then
                     do_merge = .true.
                     exit
                  end if
               end do
            end if

            f2_val = this%f_centre(i, j)*this%f_centre(i, j)

            if (do_merge) then
               ! ---- D4 merge path: raw thickness, T/S back-out with the
               ! div-eps armour (hT = T*h, so hT/h recovers the layer mean
               ! even for a vanished layer).
               do k = 1, nz
                  kg = nz + 1 - k
                  inv_h = 1.0_wp/max(h_sd(k), H_DIV_EPS)
                  t_sd(k) = hT(i, j, kg)*inv_h
                  s_sd(k) = hS(i, j, kg)*inv_h
               end do

               ! Build maps -> merge fields -> solve on nzc -> interp back.
               call massless_build_maps(h_sd, nz, H_VANISHED, nzc, hc, kc, kf)
               call massless_merge_fields(h_sd, kc, nz, nzc, &
                                          u_sd, v_sd, t_sd, s_sd, &
                                          uc, vc, tc, sc)
               call ks_precompute(nzc, this%lz_rescale, hc, &
                                  idz_s, idz_int_s, hint_s, il2_s)
               call ks_solve_column(nzc, dt, f2_val, this%rho0, &
                                    this%ri_crit, this%shearmix_rate, &
                                    this%fri_curvature, this%c_n, this%c_s, &
                                    this%lambda, this%kappa_0, this%kappa_seed, &
                                    this%kappa_trunc, this%tke_bg, this%tol_err, &
                                    this%max_inner_it, this%max_substep_it, &
                                    this%src_max_chg, this%vel_underflow, &
                                    this%eos, &
                                    hc, uc, vc, tc, sc, &
                                    idz_s, idz_int_s, hint_s, il2_s, &
                                    kappa_c, tke_c)
               call massless_interp_back(kappa_c, kc, kf, nz, kappa_avg_sd)
               call massless_interp_back(tke_c, kc, kf, nz, tke_avg_sd)
            else
               ! ---- Existing path (knob off, or knob on + healthy column):
               ! blunt gather floor at H_VANISHED, solve on nz.  Bitwise
               ! identical to pre-D4 behaviour.
               do k = 1, nz
                  kg = nz + 1 - k
                  hk = max(h_sd(k), H_VANISHED)
                  inv_h = 1.0_wp/hk
                  h_sd(k) = hk
                  t_sd(k) = hT(i, j, kg)*inv_h
                  s_sd(k) = hS(i, j, kg)*inv_h
               end do

               ! ---- Precompute the thickness grids (h_Int, 1/h, L_bdry) ----
               call ks_precompute(nz, this%lz_rescale, h_sd, &
                                  idz_s, idz_int_s, hint_s, il2_s)

               ! ---- Adaptive outer solve (background pre-step + EOS
               !      buoyancy derivatives are done inside) ----
               call ks_solve_column(nz, dt, f2_val, this%rho0, &
                                    this%ri_crit, this%shearmix_rate, &
                                    this%fri_curvature, this%c_n, this%c_s, &
                                    this%lambda, this%kappa_0, this%kappa_seed, &
                                    this%kappa_trunc, this%tke_bg, this%tol_err, &
                                    this%max_inner_it, this%max_substep_it, &
                                    this%src_max_chg, this%vel_underflow, &
                                    this%eos, &
                                    h_sd, u_sd, v_sd, t_sd, s_sd, &
                                    idz_s, idz_int_s, hint_s, il2_s, &
                                    kappa_avg_sd, tke_avg_sd)
            end if
         end if

         ! ---- Scatter with the flip; force exact 0 at bed + surface ----
         ! local interface K -> global interface (nz+2-K).
         do k = 1, nz + 1
            kg = nz + 2 - k
            this%kd_int(i, j, kg) = kappa_avg_sd(k)
            this%tke_int(i, j, kg) = tke_avg_sd(k)
         end do
         this%kd_int(i, j, 1) = 0.0_wp
         this%kd_int(i, j, nz + 1) = 0.0_wp
         this%tke_int(i, j, 1) = 0.0_wp
         this%tke_int(i, j, nz + 1) = 0.0_wp
      end do
   end subroutine kappa_shear_column_kernel

   ! =================================================================
   ! Vertex (corner) form — MOM6 VERTEX_SHEAR.
   ! Corner (ic,jc) is the SW corner of cell (ic,jc): its four cells
   ! are SW=(ic-1,jc-1), SE=(ic,jc-1), NW=(ic-1,jc), NE=(ic,jc); its
   ! four faces are u-faces (ic,jc-1)/(ic,jc) and v-faces
   ! (ic-1,jc)/(ic,jc).  (MOM6's corner (I,J) is the NE corner of cell
   ! (i,j); the translation is (I,J) -> (ic,jc) = (I+1,J+1).)
   ! =================================================================

   pure subroutine ks_gather_corner(nx, ny, nz, ic, jc, h_layer, u_face, &
                                    v_face, hT, hS, wet_t, wet_u, wet_v, &
                                    h_sd, u_sd, v_sd, t_sd, s_sd)
      !! Assemble the surface-down column at corner (ic,jc) from the
      !! 2x2 cell patch + 4 adjacent faces (JHL08 vertex form; the
      !! interpolation recipes of the reference implementation):
      !!   u,v : 2-point THICKNESS-weighted average across the corner,
      !!         with the face thickness itself a mask-weighted 2-cell
      !!         average (recomputed inline — deterministic, so the
      !!         repeated evaluation is bitwise identical to MOM6's
      !!         precomputed h_at_u/h_at_v Pass A, non-OBC-bug form).
      !!   T,S : 4-cell mask-AND-thickness-weighted average.  The
      !!         registry stores hTr = h*T, which is exactly the
      !!         weighted quantity, so we sum wet*hTr directly.
      !!   h   : 4-cell mask-weighted average (no thickness weight —
      !!         it IS the thickness).
      !! Returns RAW h (no floor) — the caller decides floor vs
      !! massless-merge exactly as the column path does.  The deliberate
      !! (SW+NE)+(SE+NW) bracketing is reproducible-sum ordering — do
      !! not reassociate.
      !$acc routine seq
      integer, intent(in) :: nx, ny, nz, ic, jc
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: u_face(nx + 1, ny, nz)
      real(wp), intent(in) :: v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: hT(nx, ny, nz)
      real(wp), intent(in) :: hS(nx, ny, nz)
      real(wp), intent(in) :: wet_t(nx, ny)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(in) :: wet_v(nx, ny + 1)
      real(wp), intent(out) :: h_sd(NZL), u_sd(NZL), v_sd(NZL)
      real(wp), intent(out) :: t_sd(NZL), s_sd(NZL)

      integer :: k, kg
      real(wp) :: w_sw, w_se, w_nw, w_ne
      real(wp) :: h_sw, h_se, h_nw, h_ne
      real(wp) :: hu_s, hu_n, hv_w, hv_e
      real(wp) :: hwt, i_hwt

      w_sw = wet_t(ic - 1, jc - 1)
      w_se = wet_t(ic, jc - 1)
      w_nw = wet_t(ic - 1, jc)
      w_ne = wet_t(ic, jc)

      do k = 1, nz
         kg = nz + 1 - k  ! global bottom-up layer -> local surface-down
         h_sw = h_layer(ic - 1, jc - 1, kg)
         h_se = h_layer(ic, jc - 1, kg)
         h_nw = h_layer(ic - 1, jc, kg)
         h_ne = h_layer(ic, jc, kg)

         ! Face thicknesses (mask-weighted 2-cell averages, inline).
         hu_s = wet_u(ic, jc - 1)*(w_sw*h_sw + w_se*h_se)/ &
                (w_sw + w_se + MASK_SUM_EPS)
         hu_n = wet_u(ic, jc)*(w_nw*h_nw + w_ne*h_ne)/ &
                (w_nw + w_ne + MASK_SUM_EPS)
         hv_w = wet_v(ic - 1, jc)*(w_sw*h_sw + w_nw*h_nw)/ &
                (w_sw + w_nw + MASK_SUM_EPS)
         hv_e = wet_v(ic, jc)*(w_se*h_se + w_ne*h_ne)/ &
                (w_se + w_ne + MASK_SUM_EPS)

         ! Thickness-weighted 2-point transverse velocity averages.
         u_sd(k) = ((u_face(ic, jc - 1, kg)*hu_s) + &
                    (u_face(ic, jc, kg)*hu_n))/ &
                   ((hu_s + hu_n) + H_TINY_CORNER)
         v_sd(k) = ((v_face(ic - 1, jc, kg)*hv_w) + &
                    (v_face(ic, jc, kg)*hv_e))/ &
                   ((hv_w + hv_e) + H_TINY_CORNER)

         ! 4-cell mask*thickness weight (diagonal pairs first).
         hwt = ((w_sw*h_sw + w_ne*h_ne) + (w_se*h_se + w_nw*h_nw))
         i_hwt = 1.0_wp/(hwt + H_DIV_EPS)
         t_sd(k) = ((w_sw*hT(ic - 1, jc - 1, kg) + w_ne*hT(ic, jc, kg)) + &
                    (w_se*hT(ic, jc - 1, kg) + w_nw*hT(ic - 1, jc, kg)))*i_hwt
         s_sd(k) = ((w_sw*hS(ic - 1, jc - 1, kg) + w_ne*hS(ic, jc, kg)) + &
                    (w_se*hS(ic, jc - 1, kg) + w_nw*hS(ic - 1, jc, kg)))*i_hwt

         ! 4-cell mask-weighted thickness (mean of the WET cells only —
         ! a partially-wet corner is not diluted by land zeros).
         h_sd(k) = hwt/((w_sw + w_ne) + (w_se + w_nw) + MASK_SUM_EPS)
      end do
   end subroutine ks_gather_corner

   pure subroutine kappa_shear_vertex_kernel(this, nx, ny, nz, dt, h_layer, &
                                             u_face, v_face, hT, hS, &
                                             wet_t, wet_u, wet_v)
      !! Per-CORNER JHL08 solve (MOM6 vertex form, Pass B).  One
      !! `do concurrent (jc, ic)` over the interior corners
      !! [2,nx]x[2,ny] — the set whose full 2x2 cell patch exists
      !! in-array, which covers every corner any owned tracer cell
      !! needs (nghost >= 1).  Ring corners stay 0.
      !!
      !! Per corner: activity test = ANY adjacent velocity face open
      !! (NOT the 4-cell corner mask — a coastline corner with 1-3 wet
      !! cells IS solved, from the mask-weighted average of the wet
      !! cells); gather (`ks_gather_corner`); non-finite guard (an
      !! unfilled OBC ghost corner must yield kd=0, not a laundered
      !! clamp value); floor-or-merge exactly as the column path; the
      !! SAME shared column solve (`ks_precompute` + `ks_solve_column`);
      !! scatter into `kd_corner` with the surface-down -> bottom-up
      !! flip.  f^2 is taken straight at the corner (no averaging).
      type(ocean_kappa_shear_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: u_face(nx + 1, ny, nz)
      real(wp), intent(in) :: v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: hT(nx, ny, nz)
      real(wp), intent(in) :: hS(nx, ny, nz)
      real(wp), intent(in) :: wet_t(nx, ny)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(in) :: wet_v(nx, ny + 1)

      integer :: ic, jc, k, kg, nzc
      real(wp) :: f2_val
      logical :: merge_on, do_merge, col_ok
      ! gathered surface-down corner column
      real(wp) :: h_sd(NZL), u_sd(NZL), v_sd(NZL), t_sd(NZL), s_sd(NZL)
      ! precomputed interface/layer grids
      real(wp) :: idz_s(NZL), idz_int_s(NZLI), hint_s(NZLI), il2_s(NZLI)
      real(wp) :: kappa_avg_sd(NZLI), tke_avg_sd(NZLI)
      ! D4 massless-merge per-column scratch (same set as the column kernel)
      real(wp) :: hc(NZL), uc(NZL), vc(NZL), tc(NZL), sc(NZL)
      real(wp) :: kappa_c(NZLI), tke_c(NZLI), kf(NZL + 1)
      integer :: kc(NZL + 1)

      merge_on = this%massless_merge

      do concurrent(jc=2:ny, ic=2:nx) &
         local(k, kg, nzc, f2_val, do_merge, col_ok, &
               h_sd, u_sd, v_sd, t_sd, s_sd, &
               idz_s, idz_int_s, hint_s, il2_s, &
               kappa_avg_sd, tke_avg_sd, &
               hc, uc, vc, tc, sc, kappa_c, tke_c, kc, kf)

         ! ---- Default outputs (kept for inactive / non-finite corners) ----
         do k = 1, nz + 1
            kappa_avg_sd(k) = 0.0_wp
            tke_avg_sd(k) = 0.0_wp
         end do

         ! ---- Activity test: any adjacent velocity face open ----
         if ((wet_u(ic, jc - 1) + wet_u(ic, jc)) + &
             (wet_v(ic - 1, jc) + wet_v(ic, jc)) > 0.0_wp) then

            call ks_gather_corner(nx, ny, nz, ic, jc, h_layer, u_face, &
                                  v_face, hT, hS, wet_t, wet_u, wet_v, &
                                  h_sd, u_sd, v_sd, t_sd, s_sd)

            ! ---- Non-finite guard (OBC ghost-corner armour).  max()
            ! and comparisons LAUNDER/skip NaN under relaxed FP, so an
            ! explicit finite check is the only reliable gate; a corrupt
            ! gather yields kd=0 for this corner, never a plausible
            ! laundered value. ----
            col_ok = .true.
            do k = 1, nz
               if (.not. (ieee_is_finite(h_sd(k)) .and. &
                          ieee_is_finite(u_sd(k)) .and. &
                          ieee_is_finite(v_sd(k)) .and. &
                          ieee_is_finite(t_sd(k)) .and. &
                          ieee_is_finite(s_sd(k)))) col_ok = .false.
            end do

            if (col_ok) then
               ! ---- I1 precheck: merge only when the knob is on AND the
               ! corner column actually carries a vanished layer.
               do_merge = .false.
               if (merge_on) then
                  do k = 1, nz
                     if (h_sd(k) < H_VANISHED) then
                        do_merge = .true.
                        exit
                     end if
                  end do
               end if

               ! f^2 straight at the corner — f is natively a corner
               ! quantity on the C-grid (signed; square it).
               f2_val = this%f_corner(ic, jc)*this%f_corner(ic, jc)

               if (do_merge) then
                  call massless_build_maps(h_sd, nz, H_VANISHED, nzc, hc, kc, kf)
                  call massless_merge_fields(h_sd, kc, nz, nzc, &
                                             u_sd, v_sd, t_sd, s_sd, &
                                             uc, vc, tc, sc)
                  call ks_precompute(nzc, this%lz_rescale, hc, &
                                     idz_s, idz_int_s, hint_s, il2_s)
                  call ks_solve_column(nzc, dt, f2_val, this%rho0, &
                                       this%ri_crit, this%shearmix_rate, &
                                       this%fri_curvature, this%c_n, this%c_s, &
                                       this%lambda, this%kappa_0, this%kappa_seed, &
                                       this%kappa_trunc, this%tke_bg, this%tol_err, &
                                       this%max_inner_it, this%max_substep_it, &
                                       this%src_max_chg, this%vel_underflow, &
                                       this%eos, &
                                       hc, uc, vc, tc, sc, &
                                       idz_s, idz_int_s, hint_s, il2_s, &
                                       kappa_c, tke_c)
                  call massless_interp_back(kappa_c, kc, kf, nz, kappa_avg_sd)
                  call massless_interp_back(tke_c, kc, kf, nz, tke_avg_sd)
               else
                  ! Blunt gather floor at H_VANISHED, solve on nz —
                  ! mirrors the column path (T,S need no back-out: the
                  ! gather already returns layer values).
                  do k = 1, nz
                     h_sd(k) = max(h_sd(k), H_VANISHED)
                  end do
                  call ks_precompute(nz, this%lz_rescale, h_sd, &
                                     idz_s, idz_int_s, hint_s, il2_s)
                  call ks_solve_column(nz, dt, f2_val, this%rho0, &
                                       this%ri_crit, this%shearmix_rate, &
                                       this%fri_curvature, this%c_n, this%c_s, &
                                       this%lambda, this%kappa_0, this%kappa_seed, &
                                       this%kappa_trunc, this%tke_bg, this%tol_err, &
                                       this%max_inner_it, this%max_substep_it, &
                                       this%src_max_chg, this%vel_underflow, &
                                       this%eos, &
                                       h_sd, u_sd, v_sd, t_sd, s_sd, &
                                       idz_s, idz_int_s, hint_s, il2_s, &
                                       kappa_avg_sd, tke_avg_sd)
               end if
            end if
         end if

         ! ---- Scatter with the flip; force exact 0 at bed + surface.
         ! No mask multiply here (the non-bug store): land is handled by
         ! the activity test, and Pass C masks the tracer-point output. ----
         do k = 1, nz + 1
            kg = nz + 2 - k
            this%kd_corner(ic, jc, kg) = kappa_avg_sd(k)
         end do
         this%kd_corner(ic, jc, 1) = 0.0_wp
         this%kd_corner(ic, jc, nz + 1) = 0.0_wp
      end do
   end subroutine kappa_shear_vertex_kernel

   pure subroutine kappa_shear_vertex_scatter(nx, ny, nzp1, geometric, kdmin, &
                                              wet_t, kd_corner, kd_int, tke_int)
      !! Corner -> tracer-point averaging (MOM6 vertex form, Pass C).
      !! Cell (i,j) reads its four corners SW=(i,j), SE=(i+1,j),
      !! NW=(i,j+1), NE=(i+1,j+1) — a read-only corner stencil into an
      !! own-cell write, safe as its own `do concurrent` but NEVER
      !! fusable with the corner solve (`kd_corner` is the required
      !! snapshot).  Two modes:
      !!   arithmetic (default): 0.25 * ((SW+NE) + (NW+SE))
      !!   geometric:  4th root of the product of the four corner
      !!     values, each floored at `kdmin` first — a geometric mean is
      !!     0 if ANY factor is 0, which would otherwise blank Kd along
      !!     every shear-zone edge.  The floor applies to the CORNER
      !!     values only, never the output: a land cell still gets
      !!     exactly 0 via the wet_t multiply.
      !! Endpoints (bed K=1, surface K=nzp1) are forced to exactly 0.
      !! `tke_int` is zeroed — the vertex form does not carry a
      !! cell-centred TKE (corner TKE deliberately not materialised).
      !! Bracketing is reproducible-sum ordering; keep literal.
      integer, intent(in) :: nx, ny, nzp1
      logical, intent(in) :: geometric
      real(wp), intent(in) :: kdmin
      real(wp), intent(in) :: wet_t(nx, ny)
      real(wp), intent(in) :: kd_corner(nx + 1, ny + 1, nzp1)
      real(wp), intent(out) :: kd_int(nx, ny, nzp1)
      real(wp), intent(out) :: tke_int(nx, ny, nzp1)

      integer :: i, j, k
      real(wp) :: c_sw, c_se, c_nw, c_ne

      do concurrent(k=1:nzp1, j=1:ny, i=1:nx) local(c_sw, c_se, c_nw, c_ne)
         if (k == 1 .or. k == nzp1) then
            kd_int(i, j, k) = 0.0_wp
         else
            c_sw = kd_corner(i, j, k)
            c_se = kd_corner(i + 1, j, k)
            c_nw = kd_corner(i, j + 1, k)
            c_ne = kd_corner(i + 1, j + 1, k)
            if (geometric) then
               kd_int(i, j, k) = wet_t(i, j)*sqrt(sqrt( &
                                                  (max(c_sw, kdmin)*max(c_ne, kdmin))* &
                                                  (max(c_nw, kdmin)*max(c_se, kdmin))))
            else
               kd_int(i, j, k) = wet_t(i, j)*0.25_wp* &
                                 ((c_sw + c_ne) + (c_nw + c_se))
            end if
         end if
         tke_int(i, j, k) = 0.0_wp
      end do
   end subroutine kappa_shear_vertex_scatter

   ! =================================================================
   ! Column solve helpers (pure, device-callable, same-module)
   ! All work in LOCAL surface-down indices (k=1 surface, k=nz bed;
   ! interface K=1 surface, K=nz+1 bed).  Reference: JHL08.
   ! =================================================================

   pure function ks_src_func(ri_crit, shearmix_rate, fri_curvature, n2, s2) &
      result(ksrc)
      !! Shear-source function K_src at one interface (JHL08 eq. for the
      !! source term): nonzero only where N^2 < Ri_c * S^2.
      !$acc routine seq
      real(wp), intent(in) :: ri_crit, shearmix_rate, fri_curvature, n2, s2
      real(wp) :: ksrc
      real(wp) :: dnom

      ksrc = 0.0_wp
      if (n2 < ri_crit*s2) then
         dnom = ri_crit*s2 + fri_curvature*n2
         if (dnom /= 0.0_wp .and. s2 > 0.0_wp) then
            ksrc = 2.0_wp*shearmix_rate*sqrt(s2)*(ri_crit*s2 - n2)/dnom
         end if
      end if
   end function ks_src_func

   pure subroutine ks_precompute(nz, lz_rescale, h_sd, &
                                 idz_o, idz_int_o, hint_o, il2_o)
      !! Build the thickness-derived interface grids that the iteration
      !! reuses: 1/h, the interface 1/dz, the harmonic-mean interface FV
      !! cell thicknesses h_Int (Sum h_Int = Sum h), and the inverse
      !! boundary length scale squared (design doc section 5.1).
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: lz_rescale
      real(wp), intent(in) :: h_sd(NZL)
      real(wp), intent(out) :: idz_o(NZL), idz_int_o(NZLI)
      real(wp), intent(out) :: hint_o(NZLI), il2_o(NZLI)

      integer :: k
      real(wp) :: hk, hkm1, hkp1, norm_l, wt_a, wt_b, i_lz2, dtop
      real(wp) :: dbot(NZLI)

      i_lz2 = 1.0_wp/(lz_rescale*lz_rescale)

      do k = 1, nz
         idz_o(k) = 1.0_wp/h_sd(k)
      end do

      idz_int_o(1) = 2.0_wp/h_sd(1)
      do k = 2, nz
         idz_int_o(k) = 2.0_wp/(h_sd(k - 1) + h_sd(k))
      end do
      idz_int_o(nz + 1) = 2.0_wp/h_sd(nz)

      hint_o(1) = 0.0_wp
      if (nz >= 2) then
         hint_o(2) = h_sd(1)
         do k = 2, nz - 1
            hk = h_sd(k)
            hkm1 = h_sd(k - 1)
            hkp1 = h_sd(k + 1)
            norm_l = 1.0_wp/(hk*(hkm1 + hkp1) + 2.0_wp*hkm1*hkp1)
            wt_a = (hk + hkp1)*hkm1*norm_l
            wt_b = (hkm1 + hk)*hkp1*norm_l
            hint_o(k) = hint_o(k) + hk*wt_a
            hint_o(k + 1) = hk*wt_b
         end do
         hint_o(nz) = hint_o(nz) + h_sd(nz)
      end if
      hint_o(nz + 1) = 0.0_wp

      dbot(nz + 1) = 0.0_wp
      do k = nz, 1, -1
         dbot(k) = dbot(k + 1) + h_sd(k)
      end do
      il2_o(1) = 0.0_wp
      il2_o(nz + 1) = 0.0_wp
      dtop = 0.0_wp
      do k = 2, nz
         dtop = dtop + h_sd(k - 1)
         if (dtop > 0.0_wp .and. dbot(k) > 0.0_wp) then
            il2_o(k) = i_lz2*(dtop + dbot(k))**2/((dtop*dbot(k))**2)
         else
            il2_o(k) = 0.0_wp
         end if
      end do
   end subroutine ks_precompute

   pure subroutine ks_projected_state(nz, dt_now, ks, ke, vel_underflow, &
                                      dbuoy_t, dbuoy_s, h_sd, idz_int_s, &
                                      u0, v0, t0, s0, kappa_ps, &
                                      u_o, v_o, t_o, s_o, c1_o, n2_o, s2_o)
      !! Mix (u0,v0,T0,S0) implicitly with `kappa_ps` over `dt_now`,
      !! restricted to the layer band [ks,ke], and recompute N^2/S^2 at
      !! interfaces (band edges blend mixed inside / original outside).
      !! Backward-Euler tridiagonal; no-slip for u,v iff the band
      !! reaches the bed (ke==nz), insulating T,S.  N^2 floored at 0
      !! (design doc section 5.5).
      !$acc routine seq
      integer, intent(in) :: nz, ks, ke
      real(wp), intent(in) :: dt_now, vel_underflow
      real(wp), intent(in) :: dbuoy_t(NZLI), dbuoy_s(NZLI)
      real(wp), intent(in) :: h_sd(NZL), idz_int_s(NZLI)
      real(wp), intent(in) :: u0(NZL), v0(NZL), t0(NZL), s0(NZL)
      real(wp), intent(in) :: kappa_ps(NZLI)
      real(wp), intent(out) :: u_o(NZL), v_o(NZL), t_o(NZL), s_o(NZL)
      real(wp), intent(out) :: c1_o(NZLI), n2_o(NZLI), s2_o(NZLI)

      integer :: k, kk
      real(wp) :: a_b, a_a, b1, b1nz, d1, bd1
      real(wp) :: ua, ub, va, vb, ta, tb, sa, sb, n2v

      do k = 1, nz
         u_o(k) = u0(k)
         v_o(k) = v0(k)
         t_o(k) = t0(k)
         s_o(k) = s0(k)
      end do
      do k = 1, nz + 1
         c1_o(k) = 0.0_wp
      end do

      if (ks <= ke .and. dt_now > 0.0_wp) then
         ! Forward sweep (top layer ks).
         a_b = dt_now*kappa_ps(ks + 1)*idz_int_s(ks + 1)
         b1 = 1.0_wp/(h_sd(ks) + a_b)
         c1_o(ks + 1) = a_b*b1
         d1 = h_sd(ks)*b1
         u_o(ks) = b1*h_sd(ks)*u0(ks)
         v_o(ks) = b1*h_sd(ks)*v0(ks)
         t_o(ks) = b1*h_sd(ks)*t0(ks)
         s_o(ks) = b1*h_sd(ks)*s0(ks)

         do k = ks + 1, ke - 1
            a_a = a_b
            a_b = dt_now*kappa_ps(k + 1)*idz_int_s(k + 1)
            bd1 = h_sd(k) + d1*a_a
            b1 = 1.0_wp/(bd1 + a_b)
            c1_o(k + 1) = a_b*b1
            d1 = bd1*b1
            u_o(k) = b1*(h_sd(k)*u0(k) + a_a*u_o(k - 1))
            v_o(k) = b1*(h_sd(k)*v0(k) + a_a*v_o(k - 1))
            t_o(k) = b1*(h_sd(k)*t0(k) + a_a*t_o(k - 1))
            s_o(k) = b1*(h_sd(k)*s0(k) + a_a*s_o(k - 1))
         end do

         ! Bottom layer of the band (ke).
         a_a = a_b
         if (ke > ks) then
            b1 = 1.0_wp/(h_sd(ke) + d1*a_a)
            t_o(ke) = b1*(h_sd(ke)*t0(ke) + a_a*t_o(ke - 1))
            s_o(ke) = b1*(h_sd(ke)*s0(ke) + a_a*s_o(ke - 1))
            if (ke == nz) then
               b1nz = 1.0_wp/((h_sd(ke) + d1*a_a) + &
                              dt_now*kappa_ps(nz + 1)*idz_int_s(nz + 1))
            else
               b1nz = b1
            end if
            u_o(ke) = b1nz*(h_sd(ke)*u0(ke) + a_a*u_o(ke - 1))
            v_o(ke) = b1nz*(h_sd(ke)*v0(ke) + a_a*v_o(ke - 1))
         else
            b1 = 1.0_wp/(h_sd(ke) + a_a)
            t_o(ke) = b1*h_sd(ke)*t0(ke)
            s_o(ke) = b1*h_sd(ke)*s0(ke)
            if (ke == nz) then
               b1nz = 1.0_wp/(h_sd(ke) + a_a + &
                              dt_now*kappa_ps(nz + 1)*idz_int_s(nz + 1))
            else
               b1nz = b1
            end if
            u_o(ke) = b1nz*h_sd(ke)*u0(ke)
            v_o(ke) = b1nz*h_sd(ke)*v0(ke)
         end if
         if (abs(u_o(ke)) < vel_underflow) u_o(ke) = 0.0_wp
         if (abs(v_o(ke)) < vel_underflow) v_o(ke) = 0.0_wp

         ! Back-substitution upward.
         do k = ke - 1, ks, -1
            u_o(k) = u_o(k) + c1_o(k + 1)*u_o(k + 1)
            v_o(k) = v_o(k) + c1_o(k + 1)*v_o(k + 1)
            t_o(k) = t_o(k) + c1_o(k + 1)*t_o(k + 1)
            s_o(k) = s_o(k) + c1_o(k + 1)*s_o(k + 1)
            if (abs(u_o(k)) < vel_underflow) u_o(k) = 0.0_wp
            if (abs(v_o(k)) < vel_underflow) v_o(k) = 0.0_wp
         end do
      else
         do k = 1, nz
            if (abs(u_o(k)) < vel_underflow) u_o(k) = 0.0_wp
            if (abs(v_o(k)) < vel_underflow) v_o(k) = 0.0_wp
         end do
      end if

      ! N^2, S^2 — mixed inside the band, original values outside.
      n2_o(1) = 0.0_wp
      n2_o(nz + 1) = 0.0_wp
      s2_o(1) = 0.0_wp
      s2_o(nz + 1) = 0.0_wp
      do kk = 2, nz
         if (kk - 1 >= ks .and. kk - 1 <= ke) then
            ua = u_o(kk - 1)
            va = v_o(kk - 1)
            ta = t_o(kk - 1)
            sa = s_o(kk - 1)
         else
            ua = u0(kk - 1)
            va = v0(kk - 1)
            ta = t0(kk - 1)
            sa = s0(kk - 1)
         end if
         if (kk >= ks .and. kk <= ke) then
            ub = u_o(kk)
            vb = v_o(kk)
            tb = t_o(kk)
            sb = s_o(kk)
         else
            ub = u0(kk)
            vb = v0(kk)
            tb = t0(kk)
            sb = s0(kk)
         end if
         n2v = idz_int_s(kk)*(dbuoy_t(kk)*(ta - tb) + dbuoy_s(kk)*(sa - sb))
         if (n2v < 0.0_wp) n2v = 0.0_wp
         n2_o(kk) = n2v
         s2_o(kk) = ((ua - ub)**2 + (va - vb)**2)*idz_int_s(kk)**2
      end do
   end subroutine ks_projected_state

   pure subroutine ks_find_kappa_tke(nz, tke_min, f2_val, &
                                     ri_crit, shearmix_rate, fri_curvature, &
                                     c_n2, c_s2, ilambda2, kappa_0, &
                                     kappa_trunc, tke_bg, tol_err, max_inner_it, &
                                     n2_in, s2_in, kappa_seed, k_q_io, &
                                     idz_s, hint_s, il2_s, e1_s, &
                                     tke_o, kappa_o, &
                                     ksrc_sc, tkedec_sc, aq_sc, dq_sc, cq_sc, &
                                     dk_sc, ck_sc, ild2_sc)
      !! The inner Picard solve (design doc section 5.4): alternate a TKE
      !! tridiagonal sweep (Dirichlet surface, e1 tail below the deepest
      !! active interface) with a kappa tridiagonal sweep (smooth
      !! truncation ramp + active-range tracking) until the Picard
      !! increment converges.  Scratch arrays are supplied by the caller
      !! to avoid double-allocating per-thread stack.
      !$acc routine seq
      integer, intent(in) :: nz, max_inner_it
      real(wp), intent(in) :: tke_min, f2_val
      real(wp), intent(in) :: ri_crit, shearmix_rate, fri_curvature
      real(wp), intent(in) :: c_n2, c_s2, ilambda2, kappa_0
      real(wp), intent(in) :: kappa_trunc, tke_bg, tol_err
      real(wp), intent(in) :: n2_in(NZLI), s2_in(NZLI)
      real(wp), intent(in) :: kappa_seed(NZLI)
      real(wp), intent(inout) :: k_q_io(NZLI)
      real(wp), intent(in) :: idz_s(NZL), hint_s(NZLI), il2_s(NZLI), e1_s(NZLI)
      real(wp), intent(out) :: tke_o(NZLI), kappa_o(NZLI)
      real(wp), intent(inout) :: ksrc_sc(NZLI), tkedec_sc(NZLI)
      real(wp), intent(inout) :: aq_sc(NZL), dq_sc(NZLI), cq_sc(NZLI)
      real(wp), intent(inout) :: dk_sc(NZLI), ck_sc(NZLI), ild2_sc(NZLI)

      integer :: kk, k, k2, it
      integer :: ks_src, ke_src, ks_kap, ke_kap, ks_kp, ke_kp, ke_tke
      integer :: ks_new, ke_new, k_lo, k_hi
      real(wp) :: cqc_l, ckc_l, bqd1_l, bq_l, tsrc_l, raw_l, dnom_l
      real(wp) :: bkd1_l, bk_l, trv, tr2v, lhs_l, rhs_l
      logical :: conv

      ks_src = nz + 2
      ke_src = 0
      ksrc_sc(1) = 0.0_wp
      ksrc_sc(nz + 1) = 0.0_wp
      do kk = 2, nz
         ksrc_sc(kk) = ks_src_func(ri_crit, shearmix_rate, fri_curvature, &
                                   n2_in(kk), s2_in(kk))
         if (ksrc_sc(kk) > 0.0_wp) then
            if (ks_src > kk) ks_src = kk
            ke_src = kk
         end if
      end do

      if (ks_src > ke_src) then
         do kk = 1, nz + 1
            tke_o(kk) = tke_min
            kappa_o(kk) = 0.0_wp
            k_q_io(kk) = 0.0_wp
         end do
         return
      end if

      do kk = 2, nz
         tkedec_sc(kk) = sqrt(c_n2*n2_in(kk) + c_s2*s2_in(kk))
      end do

      tke_o(1) = tke_bg
      do kk = 2, nz
         if (kappa_seed(kk) > 0.0_wp .and. k_q_io(kk) > 0.0_wp) then
            tke_o(kk) = kappa_seed(kk)/k_q_io(kk)
         else
            tke_o(kk) = tke_min
         end if
      end do
      tke_o(nz + 1) = tke_min

      do kk = 1, nz + 1
         kappa_o(kk) = kappa_seed(kk)
      end do
      kappa_o(1) = 0.0_wp
      kappa_o(nz + 1) = 0.0_wp

      ks_kap = 2
      ke_kap = nz
      ks_kp = 2
      ke_kp = nz

      do it = 1, max_inner_it
         ! (a) TKE tridiagonal sweep.  ke_tke is an INTERFACE index
         ! (1..nz+1; bed = nz+1), so the clamp is nz+1 and the bed
         ! Dirichlet branch fires at ke_tke == nz+1 (design doc 5.4a,
         ! ambiguity #6 — "0-based nz" = the bed interface = local nz+1).
         ke_tke = min(max(ke_kap, ke_kp) + 1, nz + 1)
         do k = 1, min(ke_tke, nz)
            aq_sc(k) = (0.5_wp*(kappa_o(k) + kappa_o(k + 1)) + kappa_0)*idz_s(k)
         end do

         dq_sc(1) = -tke_o(1)
         tke_o(1) = tke_bg
         cq_sc(2) = 0.0_wp
         cqc_l = 1.0_wp
         do kk = 2, ke_tke - 1
            dq_sc(kk) = -tke_o(kk)
            tsrc_l = (kappa_o(kk) + kappa_0)*s2_in(kk) + tke_bg*tkedec_sc(kk)
            bqd1_l = hint_s(kk)*(tkedec_sc(kk) + n2_in(kk)*k_q_io(kk)) + &
                     cqc_l*aq_sc(kk - 1)
            bq_l = 1.0_wp/(bqd1_l + aq_sc(kk))
            tke_o(kk) = bq_l*(hint_s(kk)*tsrc_l + aq_sc(kk - 1)*tke_o(kk - 1))
            cq_sc(kk + 1) = aq_sc(kk)*bq_l
            cqc_l = bqd1_l*bq_l
         end do

         if (ke_tke == nz + 1) then
            tke_o(nz + 1) = tke_min
            dq_sc(nz + 1) = 0.0_wp
         else
            kk = ke_tke
            tsrc_l = kappa_0*s2_in(kk) + tke_bg*tkedec_sc(kk)
            bq_l = 1.0_wp/(hint_s(kk)*tkedec_sc(kk) + cqc_l*aq_sc(kk - 1) + &
                           aq_sc(kk))
            cq_sc(kk + 1) = aq_sc(kk)*bq_l
            dq_sc(kk) = -tke_o(kk)
            raw_l = bq_l*(hint_s(kk)*tsrc_l + aq_sc(kk - 1)*tke_o(kk - 1))
            dnom_l = 1.0_wp - cq_sc(kk + 1)*e1_s(kk + 1)
            if (abs(dnom_l) > 1.0e-30_wp) then
               tke_o(kk) = max((raw_l + cq_sc(kk + 1)* &
                                (tke_o(kk + 1) - e1_s(kk + 1)*tke_o(kk)))/dnom_l, &
                               tke_min)
            else
               tke_o(kk) = max(raw_l, tke_min)
            end if
            dq_sc(kk) = tke_o(kk) + dq_sc(kk)
            do k2 = ke_tke + 1, nz + 1
               dq_sc(k2) = e1_s(k2)*dq_sc(k2 - 1)
               tke_o(k2) = max(tke_o(k2) + dq_sc(k2), tke_min)
               if (abs(dq_sc(k2)) < 1.0e-16_wp*tke_o(k2)) exit
            end do
         end if

         do kk = ke_tke - 1, 1, -1
            tke_o(kk) = max(tke_o(kk) + cq_sc(kk + 1)*tke_o(kk + 1), tke_min)
            dq_sc(kk) = tke_o(kk) + dq_sc(kk)
         end do

         ! (b) kappa tridiagonal sweep with truncation ramp + range track.
         ks_kp = ks_kap
         ke_kp = ke_kap
         do kk = 2, nz
            if (tke_o(kk) > 0.0_wp) then
               ild2_sc(kk) = (n2_in(kk)*ilambda2 + f2_val)/tke_o(kk) + il2_s(kk)
            else
               ild2_sc(kk) = 1.0e30_wp
            end if
         end do

         dk_sc(1) = 0.0_wp
         ck_sc(2) = 0.0_wp
         ckc_l = 1.0_wp
         ke_new = 0
         ks_new = nz
         do kk = 2, nz
            dk_sc(kk) = -kappa_o(kk)
            bkd1_l = hint_s(kk)*ild2_sc(kk) + ckc_l*idz_s(kk - 1)
            bk_l = 1.0_wp/(bkd1_l + idz_s(kk))
            kappa_o(kk) = bk_l*(idz_s(kk - 1)*kappa_o(kk - 1) + &
                                hint_s(kk)*ksrc_sc(kk))
            ck_sc(kk + 1) = idz_s(kk)*bk_l
            ckc_l = bkd1_l*bk_l

            trv = ckc_l*kappa_trunc
            tr2v = 2.0_wp*trv
            if (kappa_o(kk) < trv) then
               kappa_o(kk) = 0.0_wp
               if (kk > ke_src) then
                  ke_kap = kk - 1
                  k_q_io(kk) = 0.0_wp
                  exit
               end if
            else if (kappa_o(kk) < tr2v) then
               kappa_o(kk) = 2.0_wp*(kappa_o(kk) - trv)
            end if
            ke_new = kk
         end do
         if (ke_new > 0) ke_kap = ke_new

         if (ke_kap >= 1 .and. tke_o(ke_kap) > 0.0_wp) then
            k_q_io(ke_kap) = kappa_o(ke_kap)/tke_o(ke_kap)
         end if
         dk_sc(ke_kap) = dk_sc(ke_kap) + kappa_o(ke_kap)

         do kk = ke_kap + 2, ke_kp + 1
            if (kk > nz) exit
            dk_sc(kk) = -kappa_o(kk)
            kappa_o(kk) = 0.0_wp
            k_q_io(kk) = 0.0_wp
         end do

         ks_new = 2
         do kk = ke_kap - 1, 2, -1
            kappa_o(kk) = kappa_o(kk) + ck_sc(kk + 1)*kappa_o(kk + 1)
            if (kappa_o(kk) <= kappa_trunc) then
               kappa_o(kk) = 0.0_wp
               if (kk < ks_src) then
                  ks_kap = kk + 1
                  k_q_io(kk) = 0.0_wp
                  exit
               end if
            else if (kappa_o(kk) < 2.0_wp*kappa_trunc) then
               kappa_o(kk) = 2.0_wp*(kappa_o(kk) - kappa_trunc)
            end if
            dk_sc(kk) = dk_sc(kk) + kappa_o(kk)
            if (tke_o(kk) > 0.0_wp) then
               k_q_io(kk) = kappa_o(kk)/tke_o(kk)
            else
               k_q_io(kk) = 0.0_wp
            end if
            ks_new = kk
         end do
         ks_kap = max(ks_new, 2)

         do kk = ks_kp, ks_kap - 2
            if (kk < 2) cycle
            kappa_o(kk) = 0.0_wp
            k_q_io(kk) = 0.0_wp
         end do

         ! (c) Picard convergence test.
         k_lo = min(ks_kap, ks_kp)
         k_hi = max(ke_kap, ke_kp)
         conv = .true.
         do kk = k_lo, k_hi
            lhs_l = abs(dk_sc(kk))
            rhs_l = tol_err*(kappa_0 + kappa_o(kk) - 0.5_wp*dk_sc(kk))
            if (lhs_l > rhs_l) then
               conv = .false.
               exit
            end if
         end do
         if (conv) exit
      end do

      kappa_o(1) = 0.0_wp
      kappa_o(nz + 1) = 0.0_wp
   end subroutine ks_find_kappa_tke

   pure function ks_adaptive_dt(nz, dt_rem, itt_outer, max_substep_it, &
                                ri_crit, shearmix_rate, fri_curvature, &
                                src_max_chg, tol_err, vel_underflow, &
                                dbuoy_t, dbuoy_s, &
                                h_s, u_cur, v_cur, t_cur, s_cur, &
                                kappa_out_s, kappa_src_s, local_src_s, &
                                local_src_avg_s, ks_kap, ke_kap, &
                                idz_int_s) result(dt_now_r)
      !! Largest dt_test <= dt_rem such that, after mixing for
      !! 0.5*dt_test with `kappa_out_s`, the regenerated source stays
      !! within the tolerance bands of the accepted-state source
      !! (design doc section 5.3): a halving pass followed by a 5-step
      !! refinement pass.
      !$acc routine seq
      integer, intent(in) :: nz, itt_outer, max_substep_it, ks_kap, ke_kap
      real(wp), intent(in) :: dt_rem, ri_crit, shearmix_rate, fri_curvature
      real(wp), intent(in) :: src_max_chg, tol_err, vel_underflow
      real(wp), intent(in) :: dbuoy_t(NZLI), dbuoy_s(NZLI)
      real(wp), intent(in) :: h_s(NZL), u_cur(NZL), v_cur(NZL)
      real(wp), intent(in) :: t_cur(NZL), s_cur(NZL)
      real(wp), intent(in) :: kappa_out_s(NZLI), kappa_src_s(NZLI)
      real(wp), intent(in) :: local_src_s(NZLI), local_src_avg_s(NZLI)
      real(wp), intent(in) :: idz_int_s(NZLI)
      real(wp) :: dt_now_r

      real(wp) :: tol_max(NZLI), tol_min_a(NZLI), tol_chg_a(NZLI)
      real(wp) :: u_pr(NZL), v_pr(NZL), t_pr(NZL), s_pr(NZL)
      real(wp) :: c1_pr(NZLI), n2_pr(NZLI), s2_pr(NZLI)

      integer :: kk, k_lo, k_hi, ih, ir, max_halvings, ks_lyr, ke_lyr
      real(wp) :: dt_tst, dt_inc, dt_try, idtt, tol_dksrc_low
      real(wp) :: ksrc_tst, upper_t, lower_t, tol2_l
      logical :: valid_dt, valid_try

      if (src_max_chg == 10.0_wp) then
         tol_dksrc_low = 0.95_wp
      else
         tol_dksrc_low = (src_max_chg - 0.5_wp)/src_max_chg
      end if
      tol2_l = 2.0_wp*tol_err

      do kk = 1, nz + 1
         tol_max(kk) = kappa_src_s(kk) + src_max_chg*local_src_s(kk)
         tol_min_a(kk) = kappa_src_s(kk) - tol_dksrc_low*local_src_s(kk)
         tol_chg_a(kk) = tol2_l*local_src_avg_s(kk)
      end do

      k_lo = max(ks_kap - 1, 2)
      k_hi = min(ke_kap + 1, nz)
      ks_lyr = max(ks_kap - 1, 1)
      ke_lyr = min(ke_kap, nz)

      dt_tst = dt_rem
      valid_dt = .false.
      max_halvings = (max_substep_it + 1 - itt_outer)/2

      ! Halving pass.
      do ih = 1, max(max_halvings, 1)
         call ks_projected_state(nz, 0.5_wp*dt_tst, ks_lyr, ke_lyr, &
                                 vel_underflow, dbuoy_t, dbuoy_s, h_s, &
                                 idz_int_s, u_cur, v_cur, t_cur, s_cur, &
                                 kappa_out_s, u_pr, v_pr, t_pr, s_pr, &
                                 c1_pr, n2_pr, s2_pr)
         idtt = 0.0_wp
         if (dt_tst > 0.0_wp) idtt = 1.0_wp/dt_tst
         valid_dt = .true.
         do kk = k_lo, k_hi
            if (n2_pr(kk) < ri_crit*s2_pr(kk)) then
               ksrc_tst = ks_src_func(ri_crit, shearmix_rate, fri_curvature, &
                                      n2_pr(kk), s2_pr(kk))
               upper_t = max(tol_max(kk), kappa_src_s(kk) + idtt*tol_chg_a(kk))
               lower_t = min(tol_min_a(kk), kappa_src_s(kk) - idtt*tol_chg_a(kk))
               if (ksrc_tst > upper_t .or. ksrc_tst < lower_t) then
                  valid_dt = .false.
                  exit
               end if
            else
               lower_t = min(tol_min_a(kk), kappa_src_s(kk) - idtt*tol_chg_a(kk))
               if (0.0_wp < lower_t) then
                  valid_dt = .false.
                  exit
               end if
            end if
         end do
         if (valid_dt) exit
         dt_tst = 0.5_wp*dt_tst
      end do

      ! Refinement pass.
      dt_inc = 0.0_wp
      if (dt_tst < dt_rem .and. valid_dt) then
         dt_inc = 0.5_wp*dt_tst
         do ir = 1, 5
            dt_try = dt_tst + dt_inc
            call ks_projected_state(nz, 0.5_wp*dt_try, ks_lyr, ke_lyr, &
                                    vel_underflow, dbuoy_t, dbuoy_s, h_s, &
                                    idz_int_s, u_cur, v_cur, t_cur, s_cur, &
                                    kappa_out_s, u_pr, v_pr, t_pr, s_pr, &
                                    c1_pr, n2_pr, s2_pr)
            idtt = 0.0_wp
            if (dt_try > 0.0_wp) idtt = 1.0_wp/dt_try
            valid_try = .true.
            do kk = k_lo, k_hi
               if (n2_pr(kk) < ri_crit*s2_pr(kk)) then
                  ksrc_tst = ks_src_func(ri_crit, shearmix_rate, &
                                         fri_curvature, n2_pr(kk), s2_pr(kk))
                  upper_t = max(tol_max(kk), kappa_src_s(kk) + idtt*tol_chg_a(kk))
                  lower_t = min(tol_min_a(kk), &
                                kappa_src_s(kk) - idtt*tol_chg_a(kk))
                  if (ksrc_tst > upper_t .or. ksrc_tst < lower_t) then
                     valid_try = .false.
                     exit
                  end if
               else
                  lower_t = min(tol_min_a(kk), &
                                kappa_src_s(kk) - idtt*tol_chg_a(kk))
                  if (0.0_wp < lower_t) then
                     valid_try = .false.
                     exit
                  end if
               end if
            end do
            if (valid_try) dt_tst = dt_try
            dt_inc = 0.5_wp*dt_inc
         end do
      end if

      dt_now_r = min(dt_tst*(1.0_wp + tol_err) + dt_inc, dt_rem)
   end function ks_adaptive_dt

   pure subroutine ks_solve_column(nz, dt, f2_val, rho0, &
                                   ri_crit, shearmix_rate, fri_curvature, &
                                   c_n, c_s, lambda, kappa_0, kappa_seed_in, &
                                   kappa_trunc, tke_bg, tol_err, max_inner_it, &
                                   max_substep_it, src_max_chg, vel_underflow, &
                                   eos, &
                                   h_sd, u_sd, v_sd, t_sd, s_sd, &
                                   idz_s, idz_int_s, hint_s, il2_s, &
                                   kappa_avg_sd, tke_avg_sd)
      !! Full JHL08 column solve in surface-down indices (design doc
      !! section 5.1-5.2): background kappa_0 pre-step (no-slip bed for
      !! u,v; insulating T,S), frozen interface buoyancy derivatives,
      !! e1 tail recursion, then the adaptive predictor-corrector outer
      !! loop driving the Picard inner solve.  Returns the time-mean
      !! diffusivity `kappa_avg_sd` and TKE `tke_avg_sd` over dt.
      !$acc routine seq
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle (by value) for the buoyancy derivatives.
      integer, intent(in) :: nz, max_inner_it, max_substep_it
      real(wp), intent(in) :: dt, f2_val, rho0
      real(wp), intent(in) :: ri_crit, shearmix_rate, fri_curvature
      real(wp), intent(in) :: c_n, c_s, lambda, kappa_0, kappa_seed_in
      real(wp), intent(in) :: kappa_trunc, tke_bg, tol_err, src_max_chg
      real(wp), intent(in) :: vel_underflow
      real(wp), intent(in) :: h_sd(NZL), u_sd(NZL), v_sd(NZL)
      real(wp), intent(in) :: t_sd(NZL), s_sd(NZL)
      real(wp), intent(in) :: idz_s(NZL), idz_int_s(NZLI)
      real(wp), intent(in) :: hint_s(NZLI), il2_s(NZLI)
      real(wp), intent(out) :: kappa_avg_sd(NZLI), tke_avg_sd(NZLI)

      ! Per-thread work arrays (fixed size, surface-down).
      real(wp) :: u_c(NZL), v_c(NZL), t_c(NZL), s_c(NZL)
      real(wp) :: dbuoy_t(NZLI), dbuoy_s(NZLI)
      real(wp) :: e1(NZLI)
      real(wp) :: kappa(NZLI), k_q(NZLI), kappa_avg(NZLI), tke_avg(NZLI)
      real(wp) :: n2(NZLI), s2(NZLI), tke(NZLI)
      real(wp) :: ksrc(NZLI), tkedec(NZLI)
      real(wp) :: aq(NZL), dq(NZLI), cq(NZLI)
      real(wp) :: dk(NZLI), ck(NZLI), ild2(NZLI)
      real(wp) :: cqsav(NZLI)
      real(wp) :: u_ps(NZL), v_ps(NZL), t_ps(NZL), s_ps(NZL), c1_ps(NZLI)
      real(wp) :: n2p(NZLI), s2p(NZLI), n2c(NZLI), s2c(NZLI)
      real(wp) :: kappa_out(NZLI), kq_tmp(NZLI)
      real(wp) :: tke_pred(NZLI), kappa_pred(NZLI), kappa_mid(NZLI)
      real(wp) :: tke_fin(NZLI), kappa_pred2(NZLI)
      real(wp) :: local_src_avg(NZLI), kappa_src(NZLI), local_src(NZLI)

      integer :: k, kk, io, ii
      real(wp) :: k0dt, tke_min, c_n2, c_s2, ilambda2
      real(wp) :: ome_l, eden1_l, eden2_l, i_eden_l
      real(wp) :: a1n, a1c, b1_l, d1_l, bd1_l, base_l, b1ns, b1in
      integer :: ks_src, ke_src, ks_kap, ke_kap
      real(wp) :: dsv_dt_k, dsv_ds_k, t_int, s_int, p_int, dpres
      real(wp) :: n2v, dt_rem, dt_now, dt_wt
      integer :: ks_ps, ke_ps, ks_mid, ke_mid
      logical :: no_mixing

      k0dt = dt*kappa_0
      tke_min = max(tke_bg, 1.0e-20_wp)
      c_n2 = c_n*c_n
      c_s2 = c_s*c_s
      ilambda2 = 1.0_wp/(lambda*lambda)

      ! ---- Background-diffusion pre-step (kappa_0) ----
      if (nz == 1) then
         b1_l = 1.0_wp/(h_sd(1) + k0dt*idz_int_s(2))
         u_c(1) = b1_l*h_sd(1)*u_sd(1)
         v_c(1) = b1_l*h_sd(1)*v_sd(1)
         t_c(1) = t_sd(1)
         s_c(1) = s_sd(1)
      else
         a1n = k0dt*idz_int_s(2)
         b1_l = 1.0_wp/(h_sd(1) + a1n)
         cqsav(2) = a1n*b1_l
         d1_l = h_sd(1)*b1_l
         u_c(1) = b1_l*h_sd(1)*u_sd(1)
         v_c(1) = b1_l*h_sd(1)*v_sd(1)
         t_c(1) = b1_l*h_sd(1)*t_sd(1)
         s_c(1) = b1_l*h_sd(1)*s_sd(1)
         do k = 2, nz - 1
            a1c = a1n
            a1n = k0dt*idz_int_s(k + 1)
            bd1_l = h_sd(k) + d1_l*a1c
            b1_l = 1.0_wp/(bd1_l + a1n)
            u_c(k) = b1_l*(h_sd(k)*u_sd(k) + a1c*u_c(k - 1))
            v_c(k) = b1_l*(h_sd(k)*v_sd(k) + a1c*v_c(k - 1))
            t_c(k) = b1_l*(h_sd(k)*t_sd(k) + a1c*t_c(k - 1))
            s_c(k) = b1_l*(h_sd(k)*s_sd(k) + a1c*s_c(k - 1))
            cqsav(k + 1) = a1n*b1_l
            d1_l = bd1_l*b1_l
         end do
         a1c = a1n
         base_l = h_sd(nz) + d1_l*a1c
         b1ns = 1.0_wp/(base_l + k0dt*idz_int_s(nz + 1))
         b1in = 1.0_wp/base_l
         u_c(nz) = b1ns*(h_sd(nz)*u_sd(nz) + a1c*u_c(nz - 1))
         v_c(nz) = b1ns*(h_sd(nz)*v_sd(nz) + a1c*v_c(nz - 1))
         t_c(nz) = b1in*(h_sd(nz)*t_sd(nz) + a1c*t_c(nz - 1))
         s_c(nz) = b1in*(h_sd(nz)*s_sd(nz) + a1c*s_c(nz - 1))
         cqsav(nz + 1) = 0.0_wp
         do k = nz - 1, 1, -1
            u_c(k) = u_c(k) + cqsav(k + 1)*u_c(k + 1)
            v_c(k) = v_c(k) + cqsav(k + 1)*v_c(k + 1)
            t_c(k) = t_c(k) + cqsav(k + 1)*t_c(k + 1)
            s_c(k) = s_c(k) + cqsav(k + 1)*s_c(k + 1)
         end do
      end if

      ! ---- Frozen interface buoyancy derivatives (design 5.1) ----
      ! Pressure accumulates downward (surface-relative, p(1)=0 -> D10);
      ! dbuoy_dX = -(g/rho0) drho_dX = g*rho0*dSV_dX.
      dbuoy_t(1) = 0.0_wp
      dbuoy_s(1) = 0.0_wp
      dbuoy_t(nz + 1) = 0.0_wp
      dbuoy_s(nz + 1) = 0.0_wp
      p_int = 0.0_wp
      do kk = 2, nz
         dpres = GRAVITY*rho0*h_sd(kk - 1)
         p_int = p_int + dpres
         t_int = 0.5_wp*(t_c(kk - 1) + t_c(kk))
         s_int = 0.5_wp*(s_c(kk - 1) + s_c(kk))
         call eos_specvol_derivs(eos, t_int, s_int, p_int, &
                                 dsv_dt_k, dsv_ds_k)
         dbuoy_t(kk) = GRAVITY*rho0*dsv_dt_k
         dbuoy_s(kk) = GRAVITY*rho0*dsv_ds_k
      end do

      ! ---- Initial N^2, S^2 ----
      n2(1) = 0.0_wp
      n2(nz + 1) = 0.0_wp
      s2(1) = 0.0_wp
      s2(nz + 1) = 0.0_wp
      do kk = 2, nz
         n2v = idz_int_s(kk)*(dbuoy_t(kk)*(t_c(kk - 1) - t_c(kk)) + &
                              dbuoy_s(kk)*(s_c(kk - 1) - s_c(kk)))
         if (n2v < 0.0_wp) n2v = 0.0_wp
         n2(kk) = n2v
         s2(kk) = ((u_c(kk - 1) - u_c(kk))**2 + (v_c(kk - 1) - v_c(kk))**2)* &
                  idz_int_s(kk)**2
      end do

      ! ---- e1 tail recursion ----
      e1(nz + 1) = 0.0_wp
      ome_l = 1.0_wp
      eden2_l = kappa_0*idz_s(nz)
      do kk = nz, 2, -1
         eden1_l = hint_s(kk)*sqrt(c_n2*n2(kk) + c_s2*s2(kk)) + ome_l*eden2_l
         eden2_l = kappa_0*idz_s(kk - 1)
         i_eden_l = 1.0_wp/(eden2_l + eden1_l)
         e1(kk) = eden2_l*i_eden_l
         ome_l = eden1_l*i_eden_l
      end do
      e1(1) = 0.0_wp

      ! ---- Outer-loop init ----
      do kk = 1, nz + 1
         k_q(kk) = 0.0_wp
         kappa_avg(kk) = 0.0_wp
         tke_avg(kk) = 0.0_wp
         kappa(kk) = kappa_seed_in
      end do
      kappa(1) = 0.0_wp
      kappa(nz + 1) = 0.0_wp
      dt_rem = dt
      local_src_avg(1) = 0.0_wp
      local_src_avg(nz + 1) = 0.0_wp
      do kk = 2, nz
         if (hint_s(kk) > 0.0_wp) then
            local_src_avg(kk) = 0.1_wp*k0dt*idz_int_s(kk)/hint_s(kk)
         else
            local_src_avg(kk) = 0.0_wp
         end if
      end do

      ! ---- Adaptive outer substepping ----
      do io = 1, max_substep_it
         ! Step 1: K_src + seed TKE from previous kappa/K_Q.
         ks_src = nz + 2
         ke_src = 0
         ksrc(1) = 0.0_wp
         ksrc(nz + 1) = 0.0_wp
         do kk = 2, nz
            ksrc(kk) = ks_src_func(ri_crit, shearmix_rate, fri_curvature, &
                                   n2(kk), s2(kk))
            if (ksrc(kk) > 0.0_wp) then
               if (ks_src > kk) ks_src = kk
               ke_src = kk
            end if
         end do
         do kk = 1, nz + 1
            kappa_src(kk) = ksrc(kk)
         end do

         do kk = 2, nz
            tkedec(kk) = sqrt(c_n2*n2(kk) + c_s2*s2(kk))
         end do
         tke(1) = tke_bg
         do kk = 2, nz
            if (kappa(kk) > 0.0_wp .and. k_q(kk) > 0.0_wp) then
               tke(kk) = kappa(kk)/k_q(kk)
            else
               tke(kk) = tke_min
            end if
         end do
         tke(nz + 1) = tke_min

         kappa_out = kappa

         if (ks_src > ke_src) then
            do kk = 1, nz + 1
               kappa_out(kk) = 0.0_wp
            end do
            do kk = 2, nz
               ild2(kk) = 0.0_wp
            end do
         else
            kq_tmp = k_q
            call ks_find_kappa_tke(nz, tke_min, f2_val, ri_crit, &
                                   shearmix_rate, fri_curvature, c_n2, c_s2, &
                                   ilambda2, kappa_0, kappa_trunc, tke_bg, &
                                   tol_err, max_inner_it, n2, s2, kappa, &
                                   kq_tmp, idz_s, hint_s, il2_s, e1, tke, &
                                   kappa_out, ksrc, tkedec, aq, dq, cq, dk, &
                                   ck, ild2)
            do kk = 1, nz + 1
               k_q(kk) = kq_tmp(kk)
            end do
         end if

         ! local_src for the adaptive-dt bands (design doc 5.4d): the
         ! K_src term, the kappa_0 background-change term, and the
         ! kappa diffusive-spreading term (added only when it is
         ! positive — a net inflow into the interface).
         local_src(1) = 0.0_wp
         local_src(nz + 1) = 0.0_wp
         do kk = 2, nz
            if (hint_s(kk) > 0.0_wp) then
               local_src(kk) = ksrc(kk) + kappa_0* &
                               ((idz_s(kk - 1) + idz_s(kk))/ &
                                max(hint_s(kk), 1.0e-30_wp) + ild2(kk))
               n2v = idz_s(kk - 1)*(kappa_out(kk - 1) - kappa_out(kk)) + &
                     idz_s(kk)*(kappa_out(kk + 1) - kappa_out(kk))
               if (n2v > 0.0_wp) then
                  local_src(kk) = local_src(kk) + n2v/max(hint_s(kk), 1.0e-30_wp)
               end if
            else
               local_src(kk) = ksrc(kk)
            end if
         end do

         ! Step 2: active range of kappa_out.
         ks_kap = nz + 2
         ke_kap = 0
         do kk = 2, nz
            if (kappa_out(kk) > 0.0_wp) then
               if (ks_kap > kk) ks_kap = kk
               ke_kap = kk
            end if
         end do
         if (ke_kap == nz) kappa_out(nz + 1) = 0.0_wp
         no_mixing = (ke_kap < ks_kap)

         ! Step 3: choose dt_now.
         if (no_mixing .or. io == max_substep_it) then
            dt_now = dt_rem
         else
            dt_now = ks_adaptive_dt(nz, dt_rem, io, max_substep_it, ri_crit, &
                                    shearmix_rate, fri_curvature, src_max_chg, &
                                    tol_err, vel_underflow, dbuoy_t, dbuoy_s, &
                                    h_sd, u_c, v_c, t_c, s_c, kappa_out, &
                                    kappa_src, local_src, local_src_avg, &
                                    ks_kap, ke_kap, idz_int_s)
         end if
         do kk = 2, nz
            local_src_avg(kk) = local_src_avg(kk) + dt_now*local_src(kk)
         end do
         dt_wt = dt_now/dt

         if (no_mixing) then
            ! No source reappears; remaining kappa stays 0.
            do kk = 1, nz + 1
               tke_avg(kk) = tke_avg(kk) + dt_wt*tke(kk)
            end do
            dt_rem = 0.0_wp
         else
            ! Predictor.
            ks_ps = max(ks_kap - 1, 1)
            ke_ps = min(ke_kap, nz)
            call ks_projected_state(nz, dt_now, ks_ps, ke_ps, vel_underflow, &
                                    dbuoy_t, dbuoy_s, h_sd, idz_int_s, u_c, &
                                    v_c, t_c, s_c, kappa_out, u_ps, v_ps, &
                                    t_ps, s_ps, c1_ps, n2p, s2p)
            kq_tmp = k_q
            call ks_find_kappa_tke(nz, tke_min, f2_val, ri_crit, &
                                   shearmix_rate, fri_curvature, c_n2, c_s2, &
                                   ilambda2, kappa_0, kappa_trunc, tke_bg, &
                                   tol_err, max_inner_it, n2p, s2p, kappa_out, &
                                   kq_tmp, idz_s, hint_s, il2_s, e1, tke_pred, &
                                   kappa_pred, ksrc, tkedec, aq, dq, cq, dk, &
                                   ck, ild2)
            do kk = 1, nz + 1
               kappa_mid(kk) = 0.5_wp*(kappa_out(kk) + kappa_pred(kk))
            end do
            ks_mid = nz + 2
            ke_mid = 0
            do kk = 1, nz + 1
               if (kappa_mid(kk) > 0.0_wp) then
                  if (ks_mid > kk) ks_mid = kk
                  ke_mid = kk
               end if
            end do
            ks_ps = max(ks_mid - 1, 1)
            ke_ps = min(ke_mid, nz)

            ! Corrector (real K_Q now).
            call ks_projected_state(nz, dt_now, ks_ps, ke_ps, vel_underflow, &
                                    dbuoy_t, dbuoy_s, h_sd, idz_int_s, u_c, &
                                    v_c, t_c, s_c, kappa_mid, u_ps, v_ps, &
                                    t_ps, s_ps, c1_ps, n2c, s2c)
            call ks_find_kappa_tke(nz, tke_min, f2_val, ri_crit, &
                                   shearmix_rate, fri_curvature, c_n2, c_s2, &
                                   ilambda2, kappa_0, kappa_trunc, tke_bg, &
                                   tol_err, max_inner_it, n2c, s2c, kappa_out, &
                                   k_q, idz_s, hint_s, il2_s, e1, tke_fin, &
                                   kappa_pred2, ksrc, tkedec, aq, dq, cq, dk, &
                                   ck, ild2)
            dt_rem = dt_rem - dt_now
            do kk = 1, nz + 1
               kappa_avg(kk) = kappa_avg(kk) + &
                               dt_wt*0.5_wp*(kappa_out(kk) + kappa_pred2(kk))
               tke_avg(kk) = tke_avg(kk) + &
                             dt_wt*0.5_wp*(tke_pred(kk) + tke_fin(kk))
               kappa(kk) = kappa_pred2(kk)
            end do
            kappa(1) = 0.0_wp
            kappa(nz + 1) = 0.0_wp

            ! Step 5: full-column real-state advance for the next substep.
            if (dt_rem > 0.0_wp) then
               do kk = 1, nz + 1
                  kappa_mid(kk) = 0.5_wp*(kappa_out(kk) + kappa_pred2(kk))
               end do
               call ks_projected_state(nz, dt_now, 1, nz, vel_underflow, &
                                       dbuoy_t, dbuoy_s, h_sd, idz_int_s, u_c, &
                                       v_c, t_c, s_c, kappa_mid, u_ps, v_ps, &
                                       t_ps, s_ps, c1_ps, n2, s2)
               do k = 1, nz
                  u_c(k) = u_ps(k)
                  v_c(k) = v_ps(k)
                  t_c(k) = t_ps(k)
                  s_c(k) = s_ps(k)
               end do
            end if
         end if

         if (dt_rem <= 0.0_wp) exit
      end do

      do kk = 1, nz + 1
         kappa_avg_sd(kk) = kappa_avg(kk)
         tke_avg_sd(kk) = tke_avg(kk)
      end do
      kappa_avg_sd(1) = 0.0_wp
      kappa_avg_sd(nz + 1) = 0.0_wp
   end subroutine ks_solve_column

   pure function ocean_kappa_shear_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the kappa-shear slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_kappa_shear_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_centre) &
               + arr_bytes(this%kd_int) &
               + arr_bytes(this%tke_int) &
               + arr_bytes(this%f_corner) &
               + arr_bytes(this%kd_corner)
   end function ocean_kappa_shear_bytes

end module rdb_ocean_kappa_shear
