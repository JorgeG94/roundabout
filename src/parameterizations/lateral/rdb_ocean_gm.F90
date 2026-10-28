!! Gent-McWilliams thickness-diffusion (eddy-induced bolus transport) slot.
module rdb_ocean_gm
   !! Gent & McWilliams (1990) / Griffies (1998) skew-flux thickness
   !! diffusion.  The eddy-induced ("bolus") overturning is realized as a
   !! layer THICKNESS flux (never an explicit velocity), applied as its OWN
   !! sequential operator on the CURRENT thickness, after the dynamics of
   !! every outer step (`rdb_continuity :: continuity_gm_apply`, driven by
   !! `rdb_ocean_dyn :: run_gm_step`) — mirroring how MOM6 applies its
   !! `thickness_diffuse` step after the split-RK2 dynamics step, updating
   !! `h` in place.
   !!
   !! ### Why sequential, not folded
   !!
   !! The per-face availability cap `h_avail = A·(h − H_VANISHED)/(4·dt)`
   !! bounds four faces' outflow by
   !! `h − H_VANISHED`, so `h_new >= H_VANISHED` — but ONLY for the `h` it
   !! was computed from.  Until 2026-10 this slot computed `uhD` from the
   !! stage-entry `h` and FOLDED it into the resolved continuity sweeps
   !! (`gm_fold_x/y`), so the cap bounded the stage-entry thickness and the
   !! resolved outflow came on top: on the 1-degree Southern Ocean (z*,
   !! open steps) an 8.6 cm partial bed cell with two open faces onto
   !! fillers was drained by GM at the cap on all four faces and the
   !! resolved flux took it to −8.2e-4 m.  Computed here from the thickness
   !! the dynamics LEFT, the cap bounds what is actually there (prototype:
   !! `python_prototypes/gm_sequential/gm_sequential.py`).
   !!
   !! At each u/v face the eddy streamfunction is
   !!     Sfn_unlim = -(KhTh * dy_Cu) * Slope
   !! using the STORED isopycnal slope from the slopes slot.  The per-layer
   !! bolus transport is the vertical difference of the limited
   !! streamfunction, formed by the column recurrence (surface->bed,
   !! uhtot=0 at surface):
   !!     uhD(k) = max( min(Sfn_in_H - uhtot, h_avail_i), -h_avail_{i+1} )
   !!     uhtot  = uhtot + uhD(k)
   !! so `Sum_k uhD = 0` (closes at the bed) — mass/tracer conservative.
   !!
   !! Bottom-blocking (MOM6's rule for the GM streamfunction) acts on the
   !! UNLIMITED streamfunction first: no transport from a donor layer that
   !! lies entirely below the receiving column's bed, and a share scaled by
   !! the fraction above it for a donor layer that straddles it
   !! (`gm_block_below_bed`).  On an open z-level step this is what keeps GM
   !! from pouring deep water through the step into the fillers below the
   !! shallow column's bed.
   !!
   !! Three limiters, in order: (1) safe-streamfunction blend toward a
   !! column-spread return flow where slope > slope_max; (2) mass-
   !! availability rsum bound (the conservation guard keeping each layer
   !! >= H_VANISHED); (3) per-layer donor cap.  The donor side is keyed on
   !! the SIGN of `uhtot` (column i when uhtot<=0, else i+1).
   !!
   !! KhTh is a 2D face field, constant-filled from `&ocean_gm_nml khth`
   !! (VarMix/MEKE seam populates it later via +=), CFL-clamped per face.
   !! `gm_src(nx,ny)` carries the GM PE release (>= 0 for a stable tilted
   !! column) for future MEKE coupling.
   !!
   !! Bottom-up convention (k=1 bed, k=nz surface; interface K=1 bed,
   !! K=nz+1 surface); the slopes slot zeroes the slope at both caps so the
   !! streamfunction vanishes there and the recurrence closes.
   !!
   !! ### Partial-step z-level faces (`&vcoord_nml zfixed_closed_faces`)
   !!
   !! Under `z_fixed` a face column is not the whole water column: a layer
   !! that is an inert FILLER on either side (inside the bed or the ice
   !! draft) is a WALL for that layer at that face (`metrics%open_u/open_v
   !! == 0`).  Continuity applies that mask to the resolved flux, and GM
   !! runs on its own (it does not pass through that mask), so GM must
   !! build its overturning on the OPEN part of each face column itself.  With the knob on, each face
   !! first marks its open layers,
   !!     ok(k) = open(k) .and. live(h_W(k)) .and. live(h_E(k))
   !! (`rdb_vl_is_live`, the one vanished-layer predicate), and the
   !! recurrence then:
   !!
   !!   * gives a not-`ok` layer ZERO transport and zero availability
   !!     (`h_avail = 0`, so it adds nothing to `rsum`): the streamfunction
   !!     is carried UNCHANGED across it, so it is 0 at the BOTTOM of the
   !!     open column (it starts at 0 at the bed, and the bed-side closed /
   !!     filler layers cannot change it);
   !!   * closes the column into the TOPMOST open layer `k_top_open`
   !!     (`uhD(k_top_open) = -uhtot`) instead of `k = nz`, so the
   !!     streamfunction is also 0 at the TOP of the open column — the free
   !!     surface, or the ice base with fillers above it.  The full-column
   !!     closure into `nz` would push the residual through a filler under
   !!     the ice, or through a face whose top layer is closed.
   !!
   !! Hence `uhD == 0` on every closed face-layer and every filler, and
   !! `Sum_k uhD = 0` at every face (the open-column integral IS the column
   !! integral).  With nothing closed `ok` is all-true, `k_top_open = nz`
   !! and the arithmetic is the full-column form operation for operation.
   !! The slopes slot zeroes slope / N^2 at every interface not strictly
   !! inside a face's open column on this path, so `gm_src` (the MEKE
   !! source) sees only the open column as well.  Knob OFF ⇒ `use_open =
   !! .false.`, `ok` all-true, the masks are never indexed ⇒ byte-identical.
   !!
   !! Non-finite guard: the streamfunction limiters are `min`/`max`
   !! clamps, which nvfortran's relaxed FP lowers NaN-blind (CLAUDE.md) —
   !! a NaN slope would come out as a plausible `±rsum` transport.  A
   !! non-finite slope is therefore read as ZERO slope (no GM at that
   !! interface) before any clamp sees it; finite inputs are untouched.
   !!
   !! Default off (`&ocean_gm_nml enable=.false.`) => slot allocated but
   !! `gm_compute_transports` no-ops => bit-identical.
   !!
   !! ### FGNV streamfunction (`use_fgnv`)
   !!
   !! With `use_fgnv` the local `Sfn_unlim` is replaced, per face column, by
   !! the solution of the Ferrari, Griffies, Nurser & Vallis (2010) vertical
   !! boundary-value problem
   !!     c^2 d2(psi)/dz2 - N^2 psi = -N^2 (1+s) psi_unlim,  psi = 0 at surface and bed,
   !! discretised as a symmetric positive-definite tridiagonal system over
   !! the interior interfaces and solved by Thomas (`gm_fgnv_solve_column`,
   !! called from the existing face-column kernels).  `psi_unlim` is the
   !! (bottom-blocked) local streamfunction, built from the same per-interface
   !! diffusivity as the local form (so the EBT structure of `use_ebt_struct`
   !! enters it); the limiter chain then runs unchanged on the SOLVED psi, so
   !! `Sum_k uhD = 0` still holds by construction.
   !! `c^2 = s * c_face^2` with `s = fgnv_filter_scale` and
   !! `c_face = (cg1_L + cg1_R)/2`: the filter scale multiplies the smoothing
   !! TERM, so the effective speed is `sqrt(s)*cg1` (OM3, 0.1: `0.316 cg1`).
   !! `psi_unlim` is amplified by `(1 + s)` before the solve, so the small-`c`
   !! limit is `(1+s)*psi_unlim` and `c^2 = 0` short-circuits to exactly that
   !! (no solve).  `N^2 >= (fgnv_strat_floor*Omega)^2`.  `fgnv_c_min` is
   !! accepted for namelist parity but not applied.  Away from the caps
   !! (`c/N` small against the vertical scale of the slope) psi follows
   !! psi_unlim; near the surface and bed it is smoothed over the deformation
   !! depth `c/N` and forced to zero, replacing ad-hoc tapering.
   !!
   !! `gm_src` defaults to the PE release `rho0*KH*S^2*N^2*h/4` of the
   !! LOCAL form (`gm_pe_release`), which does not see FGNV, the EBT
   !! structure or the limiters, so with FGNV on it overstates the release
   !! actually performed near the caps.  `src_work`
   !! (`&ocean_gm_nml gm_src_form = "work"`) instead has `gm_work_source`
   !! replace it by the rate of potential-energy release of the FINAL
   !! transports, face by face (see `gm_work_face`).
   !!
   !! Refs: Gent & McWilliams (1990); Gent et al. (1995); Griffies (1998);
   !! Ferrari, Griffies, Nurser & Vallis (2010), Ocean Modelling 32, 143-156.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS, GRAVITY, OMEGA_EARTH
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, H_VANISHED, H_DIV_EPS, GRAVITY, OMEGA_EARTH
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t
   use rdb_eos, only: eos_t, eos_density_specvol_derivs
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ocean_gm_t
   public :: gm_compute_transports
   public :: gm_fgnv_solve_column
   public :: gm_work_source
   public :: gm_work_source_impl

   type :: ocean_gm_t
      !! Gent-McWilliams thickness-diffusion state.  All fields default to
      !! the inert (`enable=.false.`) configuration so an ocean run that
      !! never sets `&ocean_gm_nml` is bit-identical.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; gate on this (never on
         !! `allocated`, which misses the GPU mapping).
      logical :: enable = .false.
         !! Master switch.  Off => `gm_compute_transports` and
         !! `continuity_gm_apply` are no-ops => bit-identity preserved.  Requires the
         !! slopes slot (`&ocean_slopes_nml enable`); the loud invariant is
         !! checked at configure (`configure_ocean_gm`).
      real(wp) :: khth = 0.0_wp
         !! Thickness-diffusion coefficient KhTh (m^2/s); constant-filled
         !! into the 2D face fields (production 1e2-1e3).
      real(wp) :: khth_max_cfl = 0.1_wp
         !! Fraction of the diffusive CFL the face KH may use:
         !! `KH = min(khth, 0.25*max_cfl/(dt*(idxCu^2+idyCv^2)))`.
      real(wp) :: khth_slope_max = 0.01_wp
         !! Slope magnitude (nondim) above which the safe-streamfunction
         !! blend takes over (`slope_max`).
      real(wp) :: rho0 = 1035.0_wp
         !! Reference density (kg/m^3) for the `gm_src` PE-release scaling.
      logical :: use_fgnv = .false.
         !! FGNV streamfunction BVP on (see module docstring).  Needs the
         !! `cg1` argument of `gm_compute_transports`.
      real(wp) :: fgnv_filter_scale = 1.0_wp
         !! Smoothing-term scale `s`: `c^2 = s*c_face^2`, `psi_unlim` gain `1+s`.
      real(wp) :: fgnv_strat_floor = 1.0e-15_wp
         !! `N^2 >= (floor*Omega)^2` (fraction of the rotation rate).
      logical :: use_ebt_struct = .false.
         !! Give the thickness diffusivity the vertical shape of the
         !! equivalent-barotropic mode (Hallberg 2013; Killworth & Blundell
         !! 2003): at each interior interface the 2-D face KhTh (after the
         !! CFL clamp) is multiplied by the mean of the two straddling
         !! columns' `ebt_struct` in the layer just above that interface.
         !! Needs `ebt_struct` passed to `gm_compute_transports`; off ⇒ the
         !! kernels are byte-identical.
      logical :: src_work = .false.
         !! `&ocean_gm_nml gm_src_form = "work"`: `gm_work_source` REPLACES
         !! `gm_src` by the PE-release rate of the final bolus transports
         !! (MOM6 `MEKE_GM_SRC_ALT = False`).  Off ⇒ the local form,
         !! bit-identical.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- 2D face KhTh fields (CFL-clamped) ----
      real(wp), allocatable :: khth_u(:, :)
         !! Thickness diffusivity at u-faces (m^2/s), `(nx+1,ny)`.
      real(wp), allocatable :: khth_v(:, :)
         !! Thickness diffusivity at v-faces (m^2/s), `(nx,ny+1)`.

      ! ---- Per-layer GM bolus transports (faces) ----
      ! k=1 bed, k=nz surface.  Sum_k over each face is ~0 (closed
      ! overturning).  Step-internal: filled from the current h and spent
      ! by `continuity_gm_apply` in the same outer step (not restart state).
      real(wp), allocatable :: uhD(:, :, :)
         !! GM x-face transport (m^3/s), `(nx+1,ny,nz)`.  Filled every outer
         !! step from the post-dynamics thickness and applied at once by
         !! `continuity_gm_apply`.
      real(wp), allocatable :: vhD(:, :, :)
         !! GM y-face transport (m^3/s), `(nx,ny+1,nz)`.

      ! ---- 2D GM PE-release source (for MEKE) ----
      real(wp), allocatable :: gm_src(:, :)
         !! Potential-energy release `-1/4 Sum_k rho0 KH Slope^2 N^2 h`
         !! (W/m^2-ish; >= 0 for a stable tilted column), `(nx,ny)`.
         !! CARRIED state: filled after the dynamics of step n, read by
         !! `meke_step` at the top of step n+1, so it is in the restart
         !! registry (`gm_src`).
   contains
      procedure, non_overridable :: init => ocean_gm_init
      procedure, non_overridable :: destroy => ocean_gm_destroy
      procedure, non_overridable :: enter_data => ocean_gm_enter_data
      procedure, non_overridable :: exit_data => ocean_gm_exit_data
      procedure, non_overridable :: bytes => ocean_gm_bytes
   end type ocean_gm_t

contains

   subroutine ocean_gm_init(this, grid, nz_ml)
      !! Allocate the 2D face KhTh fields, per-layer face transports, and
      !! the `gm_src` PE-release diagnostic.  Always allocates (configure
      !! runs after init); host allocation (no `do concurrent` before
      !! enter_data).
      class(ocean_gm_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml
      if (nz < 1) nz = 1
      ! Fail loud: the column recurrence uses NZ_STACK_MAX-sized locals.
      if (nz > NZ_STACK_MAX) then
         error stop "ocean_gm_init: nz_ml exceeds NZ_STACK_MAX "// &
            "(raise NZ_STACK_MAX in rdb_constants)"
      end if
      this%nx_total = nx
      this%ny_total = ny
      this%nz_ml = nz

      allocate (this%khth_u(nx + 1, ny), source=0.0_wp)
      allocate (this%khth_v(nx, ny + 1), source=0.0_wp)
      allocate (this%uhD(nx + 1, ny, nz), source=0.0_wp)
      allocate (this%vhD(nx, ny + 1, nz), source=0.0_wp)
      allocate (this%gm_src(nx, ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_gm_init

   subroutine ocean_gm_destroy(this)
      class(ocean_gm_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%khth_u)) deallocate (this%khth_u)
      if (allocated(this%khth_v)) deallocate (this%khth_v)
      if (allocated(this%uhD)) deallocate (this%uhD)
      if (allocated(this%vhD)) deallocate (this%vhD)
      if (allocated(this%gm_src)) deallocate (this%gm_src)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_gm_destroy

   subroutine ocean_gm_enter_data(this)
      ! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule:
      ! bare polymorphic `copyin(this)` maps the stack descriptor).
      class(ocean_gm_t), intent(inout) :: this
      select type (this)
      type is (ocean_gm_t)
         call ocean_gm_enter_data_impl(this)
      end select
   end subroutine ocean_gm_enter_data

   subroutine ocean_gm_enter_data_impl(this)
      type(ocean_gm_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%khth_u, this%khth_v)
      !$acc enter data copyin(this%uhD, this%vhD, this%gm_src)
   end subroutine ocean_gm_enter_data_impl

   subroutine ocean_gm_exit_data(this)
      class(ocean_gm_t), intent(inout) :: this
      select type (this)
      type is (ocean_gm_t)
         call ocean_gm_exit_data_impl(this)
      end select
   end subroutine ocean_gm_exit_data

   subroutine ocean_gm_exit_data_impl(this)
      type(ocean_gm_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%uhD, this%vhD, this%gm_src)
      !$acc exit data delete(this%khth_u, this%khth_v)
   end subroutine ocean_gm_exit_data_impl

   ! =================================================================
   ! Compute kernel
   ! =================================================================

   subroutine gm_compute_transports(grid, metrics, this, slopes, ms, dt, &
                                    khth_ext_u, khth_ext_v, ebt_struct, cg1)
      !! Fill `uhD`/`vhD` (m^3/s) with the GM bolus thickness transport and
      !! `gm_src` with the PE release, from the CURRENT `ms%h_layer` and the
      !! stored slopes.  Run every outer step AFTER the dynamics, and
      !! immediately followed by `continuity_gm_apply` with the SAME `dt`
      !! and the same, untouched `h_layer`: the availability cap
      !! `A·(h − H_VANISHED)/(4·dt)` bounds that `h` and no other (see the
      !! module docstring).  No-op when disabled, uninitialised, or the
      !! slopes slot is absent/disabled.
      !!
      !! `khth_ext_u/khth_ext_v` (optional): spatially-varying PRE-CFL base
      !! KhTh face field from VarMix.  When present the per-face CFL clamp
      !! uses it instead of the scalar `khth`; when absent, falls back to
      !! the constant `khth` (byte-identical to the pre-VarMix path).
      !!
      !! `ebt_struct` (with `this%use_ebt_struct`): the vertical shape of
      !! KhTh.  At interior interface `K` of a face the 2-D face KhTh is
      !! multiplied by `½(ebt_struct(L,ka) + ebt_struct(R,ka))`, `ka = K`
      !! the layer just ABOVE the interface (bottom-up: layer `k` lies
      !! between interfaces `k` and `k+1`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_gm_t), intent(inout) :: this
      type(ocean_slopes_t), intent(in) :: slopes
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: khth_ext_u(:, :)
      real(wp), intent(in), optional :: khth_ext_v(:, :)
      real(wp), intent(in), optional, contiguous :: ebt_struct(:, :, :)
         !! `ocean_wave_speed_t%ebt_struct` `(nx,ny,nz)`; read only when
         !! `this%use_ebt_struct` (required then).
      real(wp), intent(in), optional, contiguous :: cg1(:, :)
         !! First-baroclinic wave speed (m/s, cell centres) from the wave-speed
         !! slot; read only when `this%use_fgnv` (required then).

      integer :: nx, ny, nz
      logical :: use_ext, use_open, use_struct, use_fg

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. slopes%is_init) return
      if (.not. allocated(ms%h_layer)) return
      if (.not. allocated(slopes%slope_x)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      if (this%nz_ml /= nz) return
      if (slopes%nz_ml /= nz) return

      use_ext = present(khth_ext_u) .and. present(khth_ext_v)
      use_open = metrics%use_closed_faces
      use_fg = this%use_fgnv
      if (use_fg .and. .not. present(cg1)) then
         error stop "gm_compute_transports: use_fgnv requires the cg1 argument"
      end if

      use_struct = this%use_ebt_struct .and. present(ebt_struct)
      if (this%use_ebt_struct .and. .not. present(ebt_struct)) then
         error stop "gm_compute_transports: use_ebt_struct needs the ebt_struct field"
      end if
      if (use_struct) then
         if (size(ebt_struct, 1) /= nx .or. size(ebt_struct, 2) /= ny .or. &
             size(ebt_struct, 3) /= nz) then
            error stop "gm_compute_transports: ebt_struct extent mismatch (needs use_ebt)"
         end if
         if (use_fg) then
            call gm_dispatch(.true., ebt_struct, cg1)
         else
            call gm_dispatch(.true., ebt_struct, metrics%areaT)
         end if
      else
         ! Inert, device-present stand-in with the dummy's shape; never read.
         if (use_fg) then
            call gm_dispatch(.false., ms%h_layer, cg1)
         else
            call gm_dispatch(.false., ms%h_layer, metrics%areaT)
         end if
      end if

   contains

      subroutine gm_dispatch(use_struct_in, ebt, cgd)
         !! The four `gm_compute_impl` call shapes (open/ext), with the
         !! structure field handed over as an explicit-shape actual.
         logical, intent(in) :: use_struct_in
         real(wp), intent(in) :: ebt(nx, ny, nz)
         real(wp), intent(in) :: cgd(nx, ny)
            !! `cg1` when FGNV is on, else an inert `(nx,ny)` stand-in.
         ! Explicit-shape flat-impl: pass the top-level allocatables so NVHPC
         ! does not descriptor-walk per launch.  When the external VarMix base
         ! is absent, pass `this%khth_u/khth_v` as a harmless placeholder for
         ! the `khth_ext_*` dummy and gate it off with `use_ext`.
         !
         ! The z-level open masks follow the continuity convention (the
         ! `open_f` dummy of the barotropic renormaliser): with the knob OFF
         ! they are `(1,1,1)` placeholders, and NVHPC derives an explicit-shape
         ! dummy's present-check extent from its DECLARED bounds, so handing
         ! over the placeholder aborts `mem:separate` ("partially present")
         ! even though `use_open = .false.` never indexes it.  The knob-off
         ! calls therefore pass the slopes' own read-only `slope_x`/`slope_y`
         ! as the inert stand-in: `(nx+1,ny,nz+1)`/`(nx,ny+1,nz+1)`, i.e. at
         ! least the `(..,nz)` the dummy declares, already device-mapped, and
         ! only ever read here.
         if (use_open) then
            if (use_ext) then
               call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                                    this%khth_slope_max, this%rho0, use_ext, use_open, &
                                    metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                                    metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                                    metrics%areaT, metrics%wet_u, metrics%wet_v, &
                                    ms%h_layer, slopes%bathy, slopes%slope_x, slopes%slope_y, &
                                    slopes%n2_u, slopes%n2_v, khth_ext_u, khth_ext_v, &
                                    metrics%open_u, metrics%open_v, &
                                    this%khth_u, this%khth_v, &
                                    this%uhD, this%vhD, this%gm_src, &
                                    use_struct_in, ebt, &
                                    use_fg, this%fgnv_filter_scale, this%fgnv_strat_floor, cgd)
            else
               call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                                    this%khth_slope_max, this%rho0, use_ext, use_open, &
                                    metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                                    metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                                    metrics%areaT, metrics%wet_u, metrics%wet_v, &
                                    ms%h_layer, slopes%bathy, slopes%slope_x, slopes%slope_y, &
                                    slopes%n2_u, slopes%n2_v, this%khth_u, this%khth_v, &
                                    metrics%open_u, metrics%open_v, &
                                    this%khth_u, this%khth_v, &
                                    this%uhD, this%vhD, this%gm_src, &
                                    use_struct_in, ebt, &
                                    use_fg, this%fgnv_filter_scale, this%fgnv_strat_floor, cgd)
            end if
         else if (use_ext) then
            call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                                 this%khth_slope_max, this%rho0, use_ext, use_open, &
                                 metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                                 metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                                 metrics%areaT, metrics%wet_u, metrics%wet_v, &
                                 ms%h_layer, slopes%bathy, slopes%slope_x, slopes%slope_y, &
                                 slopes%n2_u, slopes%n2_v, khth_ext_u, khth_ext_v, &
                                 slopes%slope_x, slopes%slope_y, &
                                 this%khth_u, this%khth_v, &
                                 this%uhD, this%vhD, this%gm_src, &
                                 use_struct_in, ebt, &
                                 use_fg, this%fgnv_filter_scale, this%fgnv_strat_floor, cgd)
         else
            call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                                 this%khth_slope_max, this%rho0, use_ext, use_open, &
                                 metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                                 metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                                 metrics%areaT, metrics%wet_u, metrics%wet_v, &
                                 ms%h_layer, slopes%bathy, slopes%slope_x, slopes%slope_y, &
                                 slopes%n2_u, slopes%n2_v, this%khth_u, this%khth_v, &
                                 slopes%slope_x, slopes%slope_y, &
                                 this%khth_u, this%khth_v, &
                                 this%uhD, this%vhD, this%gm_src, &
                                 use_struct_in, ebt, &
                                 use_fg, this%fgnv_filter_scale, this%fgnv_strat_floor, cgd)
         end if
      end subroutine gm_dispatch
   end subroutine gm_compute_transports

   subroutine gm_compute_impl(nx, ny, nz, dt, khth, khth_max_cfl, slope_max, &
                              rho0, use_ext, use_open, dy_cu, dx_cv, idxCu, &
                              idyCv, idyCu, idxCv, areaT, wet_u, wet_v, h_layer, &
                              bathy, slope_x, slope_y, &
                              n2_u, n2_v, khth_ext_u, khth_ext_v, open_u, open_v, &
                              khth_u, khth_v, uhD, vhD, gm_src, &
                              use_struct, ebt_struct, &
                              use_fgnv, fgnv_filter_scale, fgnv_strat_floor, cg1)
      !! Flat-impl GM kernel.  Passes: CFL-clamp the 2D face KhTh, the
      !! u-face and v-face column recurrences into `uhD`/`vhD`, the
      !! `gm_src` PE release.  Each face's column sweep is serial in k (the
      !! `uhtot` recurrence) but faces parallelise over (i,j).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt, khth, khth_max_cfl, slope_max, rho0
      logical, intent(in) :: use_ext
      logical, intent(in) :: use_open
         !! z-level closed faces active (`metrics%use_closed_faces`).
         !! `.false.` ⇒ `open_u`/`open_v` are never indexed.
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny)
      real(wp), intent(in) :: idyCv(nx, ny + 1)
      real(wp), intent(in) :: idyCu(nx + 1, ny)
      real(wp), intent(in) :: idxCv(nx, ny + 1)
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(in) :: wet_v(nx, ny + 1)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: bathy(nx, ny)
         !! Bed depth `D` (m, positive down; `slopes%bathy`): the bed of each
         !! column is at `z = −D` for the bottom-blocking limiter.
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: khth_ext_u(nx + 1, ny)
      real(wp), intent(in) :: khth_ext_v(nx, ny + 1)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz)
         !! Per-layer 0/1 u-face open mask (`metrics%open_u`), read only
         !! when `use_open`; an inert device-present stand-in otherwise.
      real(wp), intent(in) :: open_v(nx, ny + 1, nz)
         !! v-face twin of `open_u`.
      real(wp), intent(out) :: khth_u(nx + 1, ny)
      real(wp), intent(out) :: khth_v(nx, ny + 1)
      real(wp), intent(out) :: uhD(nx + 1, ny, nz)
      real(wp), intent(out) :: vhD(nx, ny + 1, nz)
      real(wp), intent(out) :: gm_src(nx, ny)
      logical, intent(in) :: use_struct
      logical, intent(in) :: use_fgnv
         !! FGNV BVP on.  `.false.` ⇒ `cg1` and the two scalars below are
         !! never read (`cg1` is then an inert device-present stand-in).
      real(wp), intent(in) :: fgnv_filter_scale, fgnv_strat_floor
      real(wp), intent(in) :: cg1(nx, ny)
         !! EBT vertical structure of KhTh on (`&ocean_gm_nml use_ebt_struct`).
      real(wp), intent(in) :: ebt_struct(nx, ny, nz)
         !! T-cell EBT structure; an inert device-present stand-in when
         !! `use_struct = .false.` (never indexed).

      real(wp) :: i_smax2, i4dt, n2_floor
      integer :: i, j, k

      i_smax2 = 1.0_wp/(slope_max*slope_max)
      i4dt = 1.0_wp/(4.0_wp*dt)
      n2_floor = (fgnv_strat_floor*OMEGA_EARTH)**2

      ! ---- 1. CFL-clamp the constant KhTh onto the 2D face fields, using
      ! native face metrics (idxCu/idyCu at u-faces, idxCv/idyCv at
      ! v-faces) — exact on curvilinear grids.  Wall faces -> 0.
      call gm_clamp_khth(nx, ny, dt, khth, khth_max_cfl, use_ext, &
                         khth_ext_u, khth_ext_v, idxCu, idyCu, &
                         idxCv, idyCv, wet_u, wet_v, khth_u, khth_v)

      ! ---- 2. u-face bolus transport (column recurrence).
      do concurrent(k=1:nz, j=1:ny)
         uhD(1, j, k) = 0.0_wp
         uhD(nx + 1, j, k) = 0.0_wp
      end do
      call gm_column_x(nx, ny, nz, i_smax2, i4dt, use_open, &
                       dy_cu, areaT, h_layer, bathy, slope_x, khth_u, open_u, uhD, &
                       use_struct, ebt_struct, &
                       use_fgnv, fgnv_filter_scale, n2_floor, cg1, n2_u)

      ! ---- 3. v-face bolus transport (mirror).
      do concurrent(k=1:nz, i=1:nx)
         vhD(i, 1, k) = 0.0_wp
         vhD(i, ny + 1, k) = 0.0_wp
      end do
      call gm_column_y(nx, ny, nz, i_smax2, i4dt, use_open, &
                       dx_cv, areaT, h_layer, bathy, slope_y, khth_v, open_v, vhD, &
                       use_struct, ebt_struct, &
                       use_fgnv, fgnv_filter_scale, n2_floor, cg1, n2_v)

      ! ---- 4. gm_src PE release (cell centres).
      call gm_pe_release(nx, ny, nz, slope_max, rho0, h_layer, &
                         slope_x, slope_y, n2_u, n2_v, khth_u, khth_v, gm_src, &
                         use_struct, ebt_struct)
   end subroutine gm_compute_impl

   pure subroutine gm_clamp_khth(nx, ny, dt, khth, khth_max_cfl, use_ext, &
                                 khth_ext_u, khth_ext_v, idxCu, idyCu, &
                                 idxCv, idyCv, wet_u, wet_v, khth_u, khth_v)
      !! Fill the 2D face KhTh fields from the per-face base, CFL-clamped
      !! per face (native u-face idxCu/idyCu, v-face idxCv/idyCv) and zeroed
      !! on wall faces.  Base = VarMix `khth_ext_*` when `use_ext`, else the
      !! scalar `khth`; the same diffusive-CFL `min` is applied either way.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dt, khth, khth_max_cfl
      logical, intent(in) :: use_ext
      real(wp), intent(in) :: khth_ext_u(nx + 1, ny)
      real(wp), intent(in) :: khth_ext_v(nx, ny + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny)
      real(wp), intent(in) :: idyCu(nx + 1, ny)
      real(wp), intent(in) :: idxCv(nx, ny + 1)
      real(wp), intent(in) :: idyCv(nx, ny + 1)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(in) :: wet_v(nx, ny + 1)
      real(wp), intent(out) :: khth_u(nx + 1, ny)
      real(wp), intent(out) :: khth_v(nx, ny + 1)
      integer :: i, j
      real(wp) :: idy2, idx2, kh_cfl, denom, base

      ! u-faces: native idxCu/idyCu at the u-point.
      do concurrent(j=1:ny, i=1:nx + 1) local(idy2, idx2, kh_cfl, denom, base)
         khth_u(i, j) = 0.0_wp
         if (i >= 2 .and. i <= nx .and. wet_u(i, j) > 0.0_wp) then
            base = khth
            if (use_ext) base = khth_ext_u(i, j)
            idx2 = idxCu(i, j)*idxCu(i, j)
            idy2 = idyCu(i, j)*idyCu(i, j)
            denom = dt*(idx2 + idy2)
            kh_cfl = base
            if (denom > 0.0_wp) kh_cfl = min(base, 0.25_wp*khth_max_cfl/denom)
            khth_u(i, j) = kh_cfl
         end if
      end do
      ! v-faces: native idxCv/idyCv at the v-point.
      do concurrent(j=1:ny + 1, i=1:nx) local(idy2, idx2, kh_cfl, denom, base)
         khth_v(i, j) = 0.0_wp
         if (j >= 2 .and. j <= ny .and. wet_v(i, j) > 0.0_wp) then
            base = khth
            if (use_ext) base = khth_ext_v(i, j)
            idy2 = idyCv(i, j)*idyCv(i, j)
            idx2 = idxCv(i, j)*idxCv(i, j)
            denom = dt*(idx2 + idy2)
            kh_cfl = base
            if (denom > 0.0_wp) kh_cfl = min(base, 0.25_wp*khth_max_cfl/denom)
            khth_v(i, j) = kh_cfl
         end if
      end do
   end subroutine gm_clamp_khth

   pure function gm_h_frac(h_avail_k, rsum_k) result(hf)
      !$acc routine seq
      !! Donor mass fraction `h_avail(k)/rsum_above(k)` (0 when no mass is
      !! available above).  `rsum_above(k)` is the cumulative availability
      !! from the surface down to and including layer k, so `hf in [0,1]`.
      real(wp), intent(in) :: h_avail_k, rsum_k
      real(wp) :: hf
      hf = 0.0_wp
      if (h_avail_k > 0.0_wp) hf = h_avail_k/(rsum_k + H_DIV_EPS)
   end function gm_h_frac

   pure subroutine gm_column_x(nx, ny, nz, i_smax2, i4dt, use_open, &
                               dy_cu, areaT, h_layer, bathy, slope_x, khth_u, open_u, uhD, &
                               use_struct, ebt_struct, &
                               use_fgnv, fgnv_filter_scale, n2_floor, cg1, n2_u)
      !! u-face GM streamfunction + bolus-transport column recurrence.
      !! Interior u-face (i=2..nx) pairs columns iw=i-1 (west) and i (east).
      !! Bottom-up sweep: interior interfaces Kr=2 (bed-most) -> nz
      !! (surface-most), uhtot=0 at the bed; the surface BC (Sfn=0 at
      !! Kr=nz+1) is closed after the loop by `uhD(nz)=-uhtot`, giving
      !! `Sum_k uhD=0` exactly.  Interface Kr straddles ka=Kr (above) and
      !! kb=Kr-1 (below) and fills LAYER kb; the rsum bound keys on ka, the
      !! donor `h_frac` and per-layer cap on kb.
      !!
      !! OPEN COLUMN (`use_open`, `&vcoord_nml zfixed_closed_faces`): a
      !! layer is part of the face column only if `ok(k) = open_u(i,j,k)`
      !! and it is live on both sides.  A not-`ok` layer gets `h_avail = 0`
      !! and `uhD = 0`, leaving `uhtot` (the streamfunction at its top)
      !! unchanged, and the closure lands in the topmost `ok` layer
      !! `ktop` instead of `nz` — Sfn = 0 at the bottom AND the top of the
      !! open column (see the module docstring).  `use_open = .false.` ⇒
      !! `ok` all-true, `ktop = nz`: the full-column recurrence verbatim.
      !!
      !! DIVERGENCE (MOM6 nk_linear): MOM6's `thickness_diffuse_full`
      !! sets `nk_linear = max(GV%nkml, 1)`; in
      !! ALE mode `GV%nkml = 0`, so MOM6 always runs `nk_linear = 1` and it
      !! is NOT a namelist parameter — its top layer always takes a linear
      !! return-flow closure (comment: "Balance the deeper flow with a
      !! return flow uniformly distributed though the remaining
      !! near-surface layers"; MOM6 is top-down, k=1=surface, so
      !! `k <= nk_linear` selects the SURFACE region).  Roundabout is
      !! bottom-up (k=1=bed, CLAUDE.md); this column's `kb = k - 1` is the
      !! layer BELOW the interface (bed side), so `kb <= nk_linear` would
      !! select the BED-most layers — porting MOM6's `nk_linear=1` as
      !! written would put its SURFACE return-flow region at the SEABED.
      !! The dead `nk_linear` field + its `else` branch here were removed
      !! (PR-8) rather than wired: with the field permanently 0 (the only
      !! value ever set), `kb > nk_linear` was always true for every valid
      !! `kb >= 1`, so the limited path below ran unconditionally and the
      !! deleted branch never executed.  Adding a correct surface linear
      !! return-flow region is a live, unrecorded physics divergence from
      !! MOM6's default GM — it changes GM answers and needs its own PR +
      !! an analytical test (Psi -> 0 through the top layer), not a
      !! same-PR wire-up.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: i_smax2, i4dt
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: bathy(nx, ny)
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: khth_u(nx + 1, ny)
      logical, intent(in) :: use_open
      real(wp), intent(in) :: open_u(nx + 1, ny, nz)
      real(wp), intent(inout) :: uhD(nx + 1, ny, nz)
      logical, intent(in) :: use_struct
      logical, intent(in) :: use_fgnv
      real(wp), intent(in) :: fgnv_filter_scale, n2_floor
      real(wp), intent(in) :: cg1(nx, ny)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: ebt_struct(nx, ny, nz)
         !! EBT structure; read only when `use_struct`.  Interior interface
         !! `Kr` (loop index `k`) has layer `ka = Kr` above it, so its
         !! diffusivity is `kh*½(ebt_struct(iw,j,ka) + ebt_struct(i,j,ka))`
         !! (Hallberg 2013).  `use_struct = .false.` ⇒ the original
         !! expression, byte for byte.

      integer :: i, j, k, iw, ka, kb, ktop
      real(wp) :: havL(NZ_STACK_MAX), havR(NZ_STACK_MAX)
      real(wp) :: rsumL(NZ_STACK_MAX + 1), rsumR(NZ_STACK_MAX + 1)
      real(wp) :: eL(NZ_STACK_MAX + 1), eR(NZ_STACK_MAX + 1)
      logical :: ok(NZ_STACK_MAX)
      real(wp) :: uhtot, slope, s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h
      real(wp) :: h_frac_d, uhd_k, kh, khk
      real(wp) :: hh(NZ_STACK_MAX), n2k(NZ_STACK_MAX + 1)
      real(wp) :: psiu(NZ_STACK_MAX + 1), psis(NZ_STACK_MAX + 1)
      real(wp) :: cgf, c2
      logical :: fg

      do concurrent(j=1:ny, i=2:nx) &
         local(k, iw, ka, kb, ktop, havL, havR, rsumL, rsumR, eL, eR, ok, uhtot, slope, &
               s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h, h_frac_d, uhd_k, kh, khk, &
               hh, n2k, psiu, psis, cgf, c2, fg)
         iw = i - 1
         kh = khth_u(i, j)
         ! Interface heights (bed-up from z = -D) of the two columns, for
         ! the bottom-blocking limiter (`gm_block_below_bed`).
         eL(1) = -bathy(iw, j)
         eR(1) = -bathy(i, j)
         do k = 1, nz
            eL(k + 1) = eL(k) + h_layer(iw, j, k)
            eR(k + 1) = eR(k) + h_layer(i, j, k)
         end do

         ! The face's OPEN column (all-true without closed faces).
         do k = 1, nz
            ok(k) = .true.
         end do
         if (use_open) then
            do k = 1, nz
               ok(k) = open_u(i, j, k) > 0.5_wp .and. &
                       rdb_vl_is_live(h_layer(iw, j, k)) .and. &
                       rdb_vl_is_live(h_layer(i, j, k))
            end do
         end if
         ktop = 0
         do k = nz, 1, -1
            if (ok(k)) then
               ktop = k
               exit
            end if
         end do

         ! Per-layer availability + cumulative rsum from the SURFACE down:
         ! rsum*(k) = Sum_{k'=k}^{nz} h_avail(k') (mass above interface k).
         ! A layer outside the open column has nothing available here.
         rsumL(nz + 1) = 0.0_wp
         rsumR(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            if (ok(k)) then
               havL(k) = max(i4dt*areaT(iw, j)*(h_layer(iw, j, k) - H_VANISHED), 0.0_wp)
               havR(k) = max(i4dt*areaT(i, j)*(h_layer(i, j, k) - H_VANISHED), 0.0_wp)
            else
               havL(k) = 0.0_wp
               havR(k) = 0.0_wp
            end if
            rsumL(k) = rsumL(k + 1) + havL(k)
            rsumR(k) = rsumR(k + 1) + havR(k)
         end do

         ! FGNV: replace the local streamfunction by the BVP solution on the
         ! open column (interfaces 2..ktop; zero at the bed and above ktop).
         ! `psiu` is built from the SAME per-interface diffusivity as the
         ! local form (`khk` under the EBT structure, else `kh`).
         fg = .false.
         if (use_fgnv .and. ktop >= 2) then
            cgf = 0.5_wp*(cg1(iw, j) + cg1(i, j))
            if (.not. ieee_is_finite(cgf)) cgf = 0.0_wp
            c2 = fgnv_filter_scale*cgf*cgf
            fg = .true.
            do k = 1, ktop
               hh(k) = max(2.0_wp*h_layer(iw, j, k)*h_layer(i, j, k)/ &
                           (h_layer(iw, j, k) + h_layer(i, j, k) + H_DIV_EPS), H_VANISHED)
            end do
            psiu(1) = 0.0_wp
            psiu(ktop + 1) = 0.0_wp
            n2k(1) = n2_floor
            n2k(ktop + 1) = n2_floor
            do k = 2, ktop
               slope = slope_x(i, j, k)
               if (.not. ieee_is_finite(slope)) slope = 0.0_wp
               if (use_struct) then
                  khk = kh*(0.5_wp*(ebt_struct(iw, j, k) + ebt_struct(i, j, k)))
                  psiu(k) = gm_block_below_bed(-(khk*dy_cu(i, j))*slope, &
                                               eL(k), eL(k - 1), eR(1), eR(k), eR(k - 1), eL(1))
               else
                  psiu(k) = gm_block_below_bed(-(kh*dy_cu(i, j))*slope, &
                                               eL(k), eL(k - 1), eR(1), eR(k), eR(k - 1), eL(1))
               end if
               n2k(k) = gm_fgnv_n2(n2_u(i, j, k), n2_floor)
            end do
            call gm_fgnv_solve_column(ktop, c2, 1.0_wp + fgnv_filter_scale, hh, n2k, psiu, psis)
         end if

         ! Sweep interior interfaces bed-most (Kr=2) -> surface-most (Kr=nz).
         uhtot = 0.0_wp
         uhD(i, j, nz) = 0.0_wp
         do k = 2, nz       ! k is the interface index Kr
            ka = k           ! layer above interface (surface side)
            kb = k - 1       ! layer below interface (bed side) -> uhD(kb)
            if (.not. ok(kb) .or. kb >= ktop) then
               ! Outside the open column, or its top layer (closed below):
               ! no transport, Sfn carried unchanged across it.
               uhD(i, j, kb) = 0.0_wp
               cycle
            end if
            slope = slope_x(i, j, k)
            ! NaN-safe: a non-finite slope must not reach the min/max
            ! limiters below (they launder NaN into a bound, CLAUDE.md).
            if (.not. ieee_is_finite(slope)) slope = 0.0_wp
            s2r = slope*slope*i_smax2
            if (fg) then
               sfn_unlim = psis(k)
            else
               if (use_struct) then
                  khk = kh*(0.5_wp*(ebt_struct(iw, j, ka) + ebt_struct(i, j, ka)))
                  sfn_unlim = gm_block_below_bed(-(khk*dy_cu(i, j))*slope, &
                                                 eL(k), eL(k - 1), eR(1), eR(k), eR(k - 1), eL(1))
               else
                  sfn_unlim = gm_block_below_bed(-(kh*dy_cu(i, j))*slope, &
                                                 eL(k), eL(k - 1), eR(1), eR(k), eR(k - 1), eL(1))
               end if
            end if
            if (uhtot <= 0.0_wp) then
               h_frac_d = gm_h_frac(havL(kb), rsumL(kb))
            else
               h_frac_d = gm_h_frac(havR(kb), rsumR(kb))
            end if
            sfn_safe = uhtot*(1.0_wp - h_frac_d)
            sfn_est = (sfn_unlim + s2r*sfn_safe)/(1.0_wp + s2r)
            ! Mass above interface Kr bounds the streamfunction.
            sfn_in_h = min(max(sfn_est, -rsumL(ka)), rsumR(ka))
            uhd_k = max(min(sfn_in_h - uhtot, havL(kb)), -havR(kb))
            uhD(i, j, kb) = uhd_k
            uhtot = uhtot + uhd_k
         end do
         ! Top BC (Sfn=0 above the topmost open layer, `nz` without closed
         ! faces): close the column so Sum=0.  No open layer ⇒ all zero.
         if (ktop > 0) uhD(i, j, ktop) = -uhtot
      end do
   end subroutine gm_column_x

   pure subroutine gm_column_y(nx, ny, nz, i_smax2, i4dt, use_open, &
                               dx_cv, areaT, h_layer, bathy, slope_y, khth_v, open_v, vhD, &
                               use_struct, ebt_struct, &
                               use_fgnv, fgnv_filter_scale, n2_floor, cg1, n2_v)
      !! v-face GM column recurrence — mirror of `gm_column_x` with the
      !! v-stagger.  Interior v-face (j=2..ny) pairs columns js=j-1 (south)
      !! and j (north).  See `gm_column_x` for the open-column (`use_open`)
      !! rule and the MOM6 `nk_linear` divergence note.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: i_smax2, i4dt
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: bathy(nx, ny)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: khth_v(nx, ny + 1)
      logical, intent(in) :: use_open
      real(wp), intent(in) :: open_v(nx, ny + 1, nz)
      real(wp), intent(inout) :: vhD(nx, ny + 1, nz)
      logical, intent(in) :: use_struct
      logical, intent(in) :: use_fgnv
      real(wp), intent(in) :: fgnv_filter_scale, n2_floor
      real(wp), intent(in) :: cg1(nx, ny)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: ebt_struct(nx, ny, nz)
         !! v-face twin of the `gm_column_x` structure field.

      integer :: i, j, k, js, ka, kb, ktop
      real(wp) :: havS(NZ_STACK_MAX), havN(NZ_STACK_MAX)
      real(wp) :: rsumS(NZ_STACK_MAX + 1), rsumN(NZ_STACK_MAX + 1)
      real(wp) :: eS(NZ_STACK_MAX + 1), eN(NZ_STACK_MAX + 1)
      logical :: ok(NZ_STACK_MAX)
      real(wp) :: vhtot, slope, s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h
      real(wp) :: h_frac_d, vhd_k, kh, khk
      real(wp) :: hh(NZ_STACK_MAX), n2k(NZ_STACK_MAX + 1)
      real(wp) :: psiu(NZ_STACK_MAX + 1), psis(NZ_STACK_MAX + 1)
      real(wp) :: cgf, c2
      logical :: fg

      do concurrent(j=2:ny, i=1:nx) &
         local(k, js, ka, kb, ktop, havS, havN, rsumS, rsumN, eS, eN, ok, vhtot, slope, &
               s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h, h_frac_d, vhd_k, kh, khk, &
               hh, n2k, psiu, psis, cgf, c2, fg)
         js = j - 1
         kh = khth_v(i, j)
         eS(1) = -bathy(i, js)
         eN(1) = -bathy(i, j)
         do k = 1, nz
            eS(k + 1) = eS(k) + h_layer(i, js, k)
            eN(k + 1) = eN(k) + h_layer(i, j, k)
         end do

         do k = 1, nz
            ok(k) = .true.
         end do
         if (use_open) then
            do k = 1, nz
               ok(k) = open_v(i, j, k) > 0.5_wp .and. &
                       rdb_vl_is_live(h_layer(i, js, k)) .and. &
                       rdb_vl_is_live(h_layer(i, j, k))
            end do
         end if
         ktop = 0
         do k = nz, 1, -1
            if (ok(k)) then
               ktop = k
               exit
            end if
         end do

         rsumS(nz + 1) = 0.0_wp
         rsumN(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            if (ok(k)) then
               havS(k) = max(i4dt*areaT(i, js)*(h_layer(i, js, k) - H_VANISHED), 0.0_wp)
               havN(k) = max(i4dt*areaT(i, j)*(h_layer(i, j, k) - H_VANISHED), 0.0_wp)
            else
               havS(k) = 0.0_wp
               havN(k) = 0.0_wp
            end if
            rsumS(k) = rsumS(k + 1) + havS(k)
            rsumN(k) = rsumN(k + 1) + havN(k)
         end do

         ! FGNV: replace the local streamfunction by the BVP solution on the
         ! open column (interfaces 2..ktop; zero at the bed and above ktop).
         ! `psiu` is built from the SAME per-interface diffusivity as the
         ! local form (`khk` under the EBT structure, else `kh`).
         fg = .false.
         if (use_fgnv .and. ktop >= 2) then
            cgf = 0.5_wp*(cg1(i, js) + cg1(i, j))
            if (.not. ieee_is_finite(cgf)) cgf = 0.0_wp
            c2 = fgnv_filter_scale*cgf*cgf
            fg = .true.
            do k = 1, ktop
               hh(k) = max(2.0_wp*h_layer(i, js, k)*h_layer(i, j, k)/ &
                           (h_layer(i, js, k) + h_layer(i, j, k) + H_DIV_EPS), H_VANISHED)
            end do
            psiu(1) = 0.0_wp
            psiu(ktop + 1) = 0.0_wp
            n2k(1) = n2_floor
            n2k(ktop + 1) = n2_floor
            do k = 2, ktop
               slope = slope_y(i, j, k)
               if (.not. ieee_is_finite(slope)) slope = 0.0_wp
               if (use_struct) then
                  khk = kh*(0.5_wp*(ebt_struct(i, js, k) + ebt_struct(i, j, k)))
                  psiu(k) = gm_block_below_bed(-(khk*dx_cv(i, j))*slope, &
                                               eS(k), eS(k - 1), eN(1), eN(k), eN(k - 1), eS(1))
               else
                  psiu(k) = gm_block_below_bed(-(kh*dx_cv(i, j))*slope, &
                                               eS(k), eS(k - 1), eN(1), eN(k), eN(k - 1), eS(1))
               end if
               n2k(k) = gm_fgnv_n2(n2_v(i, j, k), n2_floor)
            end do
            call gm_fgnv_solve_column(ktop, c2, 1.0_wp + fgnv_filter_scale, hh, n2k, psiu, psis)
         end if

         vhtot = 0.0_wp
         vhD(i, j, nz) = 0.0_wp
         do k = 2, nz
            ka = k
            kb = k - 1
            if (.not. ok(kb) .or. kb >= ktop) then
               vhD(i, j, kb) = 0.0_wp
               cycle
            end if
            slope = slope_y(i, j, k)
            if (.not. ieee_is_finite(slope)) slope = 0.0_wp
            s2r = slope*slope*i_smax2
            if (fg) then
               sfn_unlim = psis(k)
            else
               if (use_struct) then
                  khk = kh*(0.5_wp*(ebt_struct(i, js, ka) + ebt_struct(i, j, ka)))
                  sfn_unlim = gm_block_below_bed(-(khk*dx_cv(i, j))*slope, &
                                                 eS(k), eS(k - 1), eN(1), eN(k), eN(k - 1), eS(1))
               else
                  sfn_unlim = gm_block_below_bed(-(kh*dx_cv(i, j))*slope, &
                                                 eS(k), eS(k - 1), eN(1), eN(k), eN(k - 1), eS(1))
               end if
            end if
            if (vhtot <= 0.0_wp) then
               h_frac_d = gm_h_frac(havS(kb), rsumS(kb))
            else
               h_frac_d = gm_h_frac(havN(kb), rsumN(kb))
            end if
            sfn_safe = vhtot*(1.0_wp - h_frac_d)
            sfn_est = (sfn_unlim + s2r*sfn_safe)/(1.0_wp + s2r)
            sfn_in_h = min(max(sfn_est, -rsumS(ka)), rsumN(ka))
            vhd_k = max(min(sfn_in_h - vhtot, havS(kb)), -havN(kb))
            vhD(i, j, kb) = vhd_k
            vhtot = vhtot + vhd_k
         end do
         if (ktop > 0) vhD(i, j, ktop) = -vhtot
      end do
   end subroutine gm_column_y

   pure function gm_fgnv_n2(n2, n2_floor) result(nf)
      !$acc routine seq
      !! FGNV stratification: `max(N^2, floor)` with a non-finite `N^2` read
      !! as the floor (a bare `max` is NaN-blind under relaxed FP).
      real(wp), intent(in) :: n2, n2_floor
      real(wp) :: nf
      nf = n2_floor
      if (ieee_is_finite(n2)) nf = max(n2, n2_floor)
   end function gm_fgnv_n2

   pure subroutine gm_fgnv_solve_column(nl, c2, gain, h, n2, psi_unlim, psi)
      !$acc routine seq
      !! Ferrari, Griffies, Nurser & Vallis (2010) streamfunction column
      !! solve.  Interfaces `K = 1..nl+1` are numbered bed -> surface
      !! (`K = 1` the bed, `K = nl+1` the surface; layer `k` lies between
      !! interfaces `k` and `k+1`).  Solves
      !!     c2 d2(psi)/dz2 - N2 psi = -N2 gain psi_unlim,  psi(1) = psi(nl+1) = 0
      !! at the interior interfaces `K = 2..nl` with the symmetric form
      !!     -c2/h(K-1) psi(K-1) + [c2/h(K-1) + c2/h(K) + N2(K) dz(K)] psi(K)
      !!     - c2/h(K) psi(K+1) = N2(K) dz(K) gain psi_unlim(K),
      !!     dz(K) = (h(K-1) + h(K))/2,
      !! (second-order accurate, SPD so the Thomas sweep needs no pivoting).
      !! `gain = 1 + fgnv_filter_scale` amplifies `psi_unlim` before the
      !! solve, so the small-`c` limit is `gain*psi_unlim` and the large-`c/N`
      !! (capped) limit is unchanged.  `h` must already be floored
      !! (`H_VANISHED`) and `n2` floored positive; `c2 >= 0` is the squared
      !! filter speed (`c2 = 0` returns `psi = gain*psi_unlim` exactly, no
      !! solve).  `nl <= NZ_STACK_MAX`.
      integer, intent(in) :: nl
      real(wp), intent(in) :: c2, gain
      real(wp), intent(in) :: h(nl)
      real(wp), intent(in) :: n2(nl + 1)
      real(wp), intent(in) :: psi_unlim(nl + 1)
      real(wp), intent(out) :: psi(nl + 1)

      real(wp) :: cp(NZ_STACK_MAX), dp(NZ_STACK_MAX)
      real(wp) :: lo, up, dz, dg, den
      integer :: k

      psi(1) = 0.0_wp
      psi(nl + 1) = 0.0_wp
      if (nl < 2) return
      if (c2 <= 0.0_wp) then
         do k = 2, nl
            psi(k) = gain*psi_unlim(k)
         end do
         return
      end if
      cp(1) = 0.0_wp
      dp(1) = 0.0_wp
      do k = 2, nl
         dz = 0.5_wp*(h(k - 1) + h(k))
         lo = c2/h(k - 1)
         up = c2/h(k)
         dg = lo + up + n2(k)*dz
         den = dg + lo*cp(k - 1)
         cp(k) = -up/den
         dp(k) = (n2(k)*dz*gain*psi_unlim(k) + lo*dp(k - 1))/den
      end do
      psi(nl) = dp(nl)
      do k = nl - 1, 2, -1
         psi(k) = dp(k) - cp(k)*psi(k + 1)
      end do
   end subroutine gm_fgnv_solve_column

   pure function gm_block_below_bed(sfn, e_top_l, e_bot_l, bed_r, e_top_r, e_bot_r, bed_l) &
      result(sfn_b)
      !$acc routine seq
      !! MOM6 bottom-blocking ("Avoid moving dense water upslope from below the
      !! level of the bottom on the receiving side").  `sfn` is the unlimited
      !! streamfunction at an interface: the transport of everything BELOW
      !! it, `> 0` from L to R.  Its donor layer is the one just below the
      !! interface on the donor side (`[e_bot, e_top]`).
      !!
      !!   * donor layer entirely below the RECEIVING column's bed
      !!     (`e_top < bed`): zero — GM may not push water into the rock (on
      !!     an open z-level step, into the fillers below the shallow column's
      !!     bed, which the remap then hands to its bottom live layer: a
      !!     transport THROUGH the step);
      !!   * donor layer straddling that bed: scaled by the fraction of the
      !!     donor layer above it.
      !!
      !! Bottom-up heights (m, +up): `e_top_*`/`e_bot_*` the top/bottom of the
      !! donor layer on each side, `bed_*` each column's bed `−D`.
      real(wp), intent(in) :: sfn, e_top_l, e_bot_l, bed_r, e_top_r, e_bot_r, bed_l
      real(wp) :: sfn_b
      sfn_b = sfn
      if (sfn > 0.0_wp) then
         if (e_top_l < bed_r) then
            sfn_b = 0.0_wp
         else if (bed_r > e_bot_l) then
            sfn_b = sfn*((e_top_l - bed_r)/((e_top_l - e_bot_l) + H_DIV_EPS))
         end if
      else
         if (e_top_r < bed_l) then
            sfn_b = 0.0_wp
         else if (bed_l > e_bot_r) then
            sfn_b = sfn*((e_top_r - bed_l)/((e_top_r - e_bot_r) + H_DIV_EPS))
         end if
      end if
   end function gm_block_below_bed

   pure function gm_clamp_slope(s, smax) result(sc)
      !$acc routine seq
      !! Clamp a slope to +/- smax (bounded slope for the PE release).
      !! A non-finite slope reads as 0 (no release) rather than being
      !! laundered into `±smax` by the NaN-blind clamp (CLAUDE.md).
      real(wp), intent(in) :: s, smax
      real(wp) :: sc
      sc = s
      if (.not. ieee_is_finite(sc)) sc = 0.0_wp
      if (sc > smax) sc = smax
      if (sc < -smax) sc = -smax
   end function gm_clamp_slope

   pure function gm_pos_n2(n2) result(np)
      !$acc routine seq
      !! `max(N^2, 0)` for the PE release, with a non-finite N^2 read as 0
      !! (a bare `max` is NaN-blind under relaxed FP, CLAUDE.md).
      real(wp), intent(in) :: n2
      real(wp) :: np
      np = 0.0_wp
      if (ieee_is_finite(n2)) np = max(n2, 0.0_wp)
   end function gm_pos_n2

   pure subroutine gm_pe_release(nx, ny, nz, slope_max, rho0, h_layer, &
                                 slope_x, slope_y, n2_u, n2_v, &
                                 khth_u, khth_v, gm_src, &
                                 use_struct, ebt_struct)
      !! GM potential-energy release at cell centres for the MEKE seam:
      !!   gm_src = 1/4 * Sum_k rho0 * (KH*Slope^2*N^2) * h
      !! over the four straddling faces, summed over interior interfaces;
      !! slope clamped to `slope_max`, N^2 floored at 0.  `gm_src >= 0` for
      !! a stable tilted column.  An interface value is attributed to the
      !! layer below it (kb=K-1); bed + surface carry zero slope/N^2.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: slope_max, rho0
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: khth_u(nx + 1, ny)
      real(wp), intent(in) :: khth_v(nx, ny + 1)
      real(wp), intent(out) :: gm_src(nx, ny)
      logical, intent(in) :: use_struct
      real(wp), intent(in) :: ebt_struct(nx, ny, nz)
         !! EBT structure; with `use_struct` each face KhTh at interface
         !! `Kr` is scaled by the mean structure of its two columns in the
         !! layer above (`ka = Kr`), exactly as in `gm_column_x/y`, so MEKE
         !! sees the structured diffusivity.  A domain-edge face has one
         !! real column: its structure is that column's own (its `khth` is
         !! the wall value 0 anyway).

      integer :: i, j, k, kb, iw0, ie0, js0, jn0
      real(wp) :: acc, fsum, sx_w, sx_e, sy_s, sy_n
      real(wp) :: n2w, n2e, n2s, n2n, kuw, kue, kvs, kvn

      do concurrent(j=1:ny, i=1:nx) &
         local(k, kb, acc, fsum, sx_w, sx_e, sy_s, sy_n, n2w, n2e, n2s, n2n, &
               kuw, kue, kvs, kvn, iw0, ie0, js0, jn0)
         acc = 0.0_wp
         iw0 = max(i - 1, 1)
         ie0 = min(i + 1, nx)
         js0 = max(j - 1, 1)
         jn0 = min(j + 1, ny)
         do k = 2, nz                 ! interior interface index Kr
            kb = k - 1                 ! layer below the interface
            sx_w = gm_clamp_slope(slope_x(i, j, k), slope_max)
            sx_e = gm_clamp_slope(slope_x(i + 1, j, k), slope_max)
            sy_s = gm_clamp_slope(slope_y(i, j, k), slope_max)
            sy_n = gm_clamp_slope(slope_y(i, j + 1, k), slope_max)
            n2w = gm_pos_n2(n2_u(i, j, k))
            n2e = gm_pos_n2(n2_u(i + 1, j, k))
            n2s = gm_pos_n2(n2_v(i, j, k))
            n2n = gm_pos_n2(n2_v(i, j + 1, k))
            if (use_struct) then
               kuw = khth_u(i, j)*(0.5_wp*(ebt_struct(iw0, j, k) + ebt_struct(i, j, k)))
               kue = khth_u(i + 1, j)*(0.5_wp*(ebt_struct(i, j, k) + ebt_struct(ie0, j, k)))
               kvs = khth_v(i, j)*(0.5_wp*(ebt_struct(i, js0, k) + ebt_struct(i, j, k)))
               kvn = khth_v(i, j + 1)*(0.5_wp*(ebt_struct(i, j, k) + ebt_struct(i, jn0, k)))
               fsum = (kuw*sx_w*sx_w*n2w + kue*sx_e*sx_e*n2e) + &
                      (kvs*sy_s*sy_s*n2s + kvn*sy_n*sy_n*n2n)
            else
               fsum = (khth_u(i, j)*sx_w*sx_w*n2w + khth_u(i + 1, j)*sx_e*sx_e*n2e) + &
                      (khth_v(i, j)*sy_s*sy_s*n2s + khth_v(i, j + 1)*sy_n*sy_n*n2n)
            end if
            acc = acc + fsum*h_layer(i, j, kb)
         end do
         gm_src(i, j) = 0.25_wp*rho0*acc
      end do
   end subroutine gm_pe_release

   subroutine gm_work_source(grid, metrics, eos, this, slopes, ms)
      !! Replace `gm_src` by the potential-energy release rate of the FINAL
      !! bolus transports (`&ocean_gm_nml gm_src_form = "work"`).  Run after
      !! `gm_compute_transports` and BEFORE `continuity_gm_apply`, so
      !! `uhD`/`vhD` (limited, FGNV, EBT-structured) and `h_layer` are the
      !! pair the operator will apply.  T/S come from the slopes slot's
      !! vert-filled copies (refreshed at thermo cadence); interface heights
      !! are rebuilt from the CURRENT `h_layer`.  No-op unless `src_work`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(eos_t), intent(in) :: eos
      type(ocean_gm_t), intent(inout) :: this
      type(ocean_slopes_t), intent(in) :: slopes
      type(multilayer_state_t), intent(in) :: ms

      integer :: nx, ny, nz

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. this%src_work) return
      if (.not. slopes%is_init) return
      if (.not. allocated(ms%h_layer)) return
      if (.not. allocated(slopes%t_fill)) return
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      if (this%nz_ml /= nz) return
      if (slopes%nz_ml /= nz) return

      ! The z-level open masks follow `gm_dispatch`: knob off ⇒ the slopes'
      ! own (nx+1,ny,nz+1) arrays stand in (device-present, never read).
      if (metrics%use_closed_faces) then
         call gm_work_source_impl(nx, ny, nz, eos, this%rho0, .true., metrics%areaT, &
                                  ms%h_layer, slopes%bathy, slopes%t_fill, slopes%s_fill, &
                                  this%uhD, this%vhD, metrics%open_u, metrics%open_v, this%gm_src)
      else
         call gm_work_source_impl(nx, ny, nz, eos, this%rho0, .false., metrics%areaT, &
                                  ms%h_layer, slopes%bathy, slopes%t_fill, slopes%s_fill, &
                                  this%uhD, this%vhD, slopes%slope_x, slopes%slope_y, this%gm_src)
      end if
   end subroutine gm_work_source

   subroutine gm_work_source_impl(nx, ny, nz, eos, rho0, use_open, areaT, h_layer, bathy, &
                                  t_fill, s_fill, uhD, vhD, open_u, open_v, gm_src)
      !! Flat-impl: `gm_src(i,j) = (1/2A) (W_u(i) + W_u(i+1) + W_v(j) + W_v(j+1))`
      !! (W m^-2), each face's work from `gm_work_face`.  A face is shared by
      !! its two cells, hence the 1/2.  Wall faces carry no transport.
      integer, intent(in) :: nx, ny, nz
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rho0
      logical, intent(in) :: use_open
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: bathy(nx, ny)
      real(wp), intent(in) :: t_fill(nx, ny, nz)
      real(wp), intent(in) :: s_fill(nx, ny, nz)
      real(wp), intent(in) :: uhD(nx + 1, ny, nz)
      real(wp), intent(in) :: vhD(nx, ny + 1, nz)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz)
         !! Per-layer u-face open mask, read only when `use_open`.
      real(wp), intent(in) :: open_v(nx, ny + 1, nz)
      real(wp), intent(out) :: gm_src(nx, ny)

      integer :: i, j, k
      real(wp) :: wsum
      real(wp) :: uc(NZ_STACK_MAX), oc(NZ_STACK_MAX)

      do concurrent(j=1:ny, i=1:nx) local(k, wsum, uc, oc)
         wsum = 0.0_wp
         if (i >= 2) then               ! west face, columns (i-1, i)
            do k = 1, nz
               uc(k) = uhD(i, j, k)
               oc(k) = 1.0_wp
               if (use_open) oc(k) = open_u(i, j, k)
            end do
            wsum = wsum + gm_work_face(nx, ny, nz, eos, rho0, h_layer, bathy, t_fill, s_fill, &
                                       i - 1, j, i, j, uc, oc, use_open)
         end if
         if (i + 1 <= nx) then          ! east face, columns (i, i+1)
            do k = 1, nz
               uc(k) = uhD(i + 1, j, k)
               oc(k) = 1.0_wp
               if (use_open) oc(k) = open_u(i + 1, j, k)
            end do
            wsum = wsum + gm_work_face(nx, ny, nz, eos, rho0, h_layer, bathy, t_fill, s_fill, &
                                       i, j, i + 1, j, uc, oc, use_open)
         end if
         if (j >= 2) then               ! south face, columns (j-1, j)
            do k = 1, nz
               uc(k) = vhD(i, j, k)
               oc(k) = 1.0_wp
               if (use_open) oc(k) = open_v(i, j, k)
            end do
            wsum = wsum + gm_work_face(nx, ny, nz, eos, rho0, h_layer, bathy, t_fill, s_fill, &
                                       i, j - 1, i, j, uc, oc, use_open)
         end if
         if (j + 1 <= ny) then          ! north face, columns (j, j+1)
            do k = 1, nz
               uc(k) = vhD(i, j + 1, k)
               oc(k) = 1.0_wp
               if (use_open) oc(k) = open_v(i, j + 1, k)
            end do
            wsum = wsum + gm_work_face(nx, ny, nz, eos, rho0, h_layer, bathy, t_fill, s_fill, &
                                       i, j, i, j + 1, uc, oc, use_open)
         end if
         gm_src(i, j) = 0.5_wp*wsum/areaT(i, j)
      end do
   end subroutine gm_work_source_impl

   pure function gm_work_face(nx, ny, nz, eos, rho0, h_layer, bathy, t_fill, s_fill, &
                              il, jl, ir, jr, uhf, opf, use_open) result(w)
      !$acc routine seq
      !! Potential-energy release rate (W) of the bolus transports `uhf(k)`
      !! (m^3/s, positive from column L to column R, bottom-up layers)
      !! through one face.  Boussinesq PE per area is
      !! `g Sum_k rho_k (e_{k+1}^2 - e_k^2)/2`, `e` the geopotential interface
      !! height (surface datum z = 0, bed at `-bathy`).  The layer mass flux
      !! `uhD_k` moves the interfaces (`de_K/dt = div U_K`, `U_K = Sum_{k>=K}
      !! uhD_k` the transport ABOVE interface K, 0 at the bed and the
      !! surface because `Sum uhD = 0`) and advects layer density at the
      !! face-mean value, `h_k d rho_k/dt = -q_k . grad rho_k`.  Then
      !!
      !!   W = -dPE/dt = g Sum_k uhD_k dRho_k ebar_k + g Sum_{K=2}^{nz} U_K (D_R - D_L)
      !!
      !! with `dRho_k = rho_R - rho_L` along layer k, `ebar_k` the mean of the
      !! four interface heights bounding layer k in both columns and
      !! `D = (rho_below - rho_above) e_K` per column (the second sum is the
      !! interface-motion term after one horizontal integration by parts).
      !! Both sums are positive for GM flattening a tilted front.  Density
      !! differences use `drho/dT`, `drho/dS` at the face-mean T/S and the
      !! hydrostatic pressure of the point (`rho0 g` times the water above),
      !! as the slope builder does.  Layer sweep runs surface -> bed.
      integer, intent(in) :: nx, ny, nz, il, jl, ir, jr
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rho0
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: bathy(nx, ny)
      real(wp), intent(in) :: t_fill(nx, ny, nz)
      real(wp), intent(in) :: s_fill(nx, ny, nz)
      real(wp), intent(in) :: uhf(nz)
      real(wp), intent(in) :: opf(nz)
      logical, intent(in) :: use_open
      real(wp) :: w

      integer :: k
      real(wp) :: el(NZ_STACK_MAX + 1), er(NZ_STACK_MAX + 1)
      real(wp) :: ubig, p, tm, sm, rho, dsv_dt, dsv_ds, drdt, drds
      real(wp) :: drl, drr, drho_layer, ebar

      el(1) = -bathy(il, jl)
      er(1) = -bathy(ir, jr)
      do k = 1, nz
         el(k + 1) = el(k) + h_layer(il, jl, k)
         er(k + 1) = er(k) + h_layer(ir, jr, k)
      end do

      w = 0.0_wp
      ubig = 0.0_wp
      do k = nz, 1, -1
         if (uhf(k) /= 0.0_wp) then
            ! Layer term at the layer-mid pressure and face-mean T/S.
            p = 0.5_wp*GRAVITY*rho0*((el(nz + 1) - 0.5_wp*(el(k) + el(k + 1))) + &
                                     (er(nz + 1) - 0.5_wp*(er(k) + er(k + 1))))
            tm = 0.5_wp*(t_fill(il, jl, k) + t_fill(ir, jr, k))
            sm = 0.5_wp*(s_fill(il, jl, k) + s_fill(ir, jr, k))
            call eos_density_specvol_derivs(eos, tm, sm, p, rho, dsv_dt, dsv_ds)
            drdt = -(rho*rho)*dsv_dt
            drds = -(rho*rho)*dsv_ds
            drho_layer = drdt*(t_fill(ir, jr, k) - t_fill(il, jl, k)) + &
                         drds*(s_fill(ir, jr, k) - s_fill(il, jl, k))
            ebar = 0.25_wp*((el(k) + el(k + 1)) + (er(k) + er(k + 1)))
            w = w + GRAVITY*uhf(k)*drho_layer*ebar
         end if
         ubig = ubig + uhf(k)          ! transport above interface K = k
         if (k >= 2 .and. ubig /= 0.0_wp) then
            if ((.not. use_open) .or. (opf(k) > 0.0_wp .and. opf(k - 1) > 0.0_wp)) then
               ! Interface term: kb = k-1 below, ka = k above.
               p = 0.5_wp*GRAVITY*rho0*((el(nz + 1) - el(k)) + (er(nz + 1) - er(k)))
               tm = 0.25_wp*((t_fill(il, jl, k) + t_fill(ir, jr, k)) + &
                             (t_fill(il, jl, k - 1) + t_fill(ir, jr, k - 1)))
               sm = 0.25_wp*((s_fill(il, jl, k) + s_fill(ir, jr, k)) + &
                             (s_fill(il, jl, k - 1) + s_fill(ir, jr, k - 1)))
               call eos_density_specvol_derivs(eos, tm, sm, p, rho, dsv_dt, dsv_ds)
               drdt = -(rho*rho)*dsv_dt
               drds = -(rho*rho)*dsv_ds
               drl = drdt*(t_fill(il, jl, k - 1) - t_fill(il, jl, k)) + &
                     drds*(s_fill(il, jl, k - 1) - s_fill(il, jl, k))
               drr = drdt*(t_fill(ir, jr, k - 1) - t_fill(ir, jr, k)) + &
                     drds*(s_fill(ir, jr, k - 1) - s_fill(ir, jr, k))
               w = w + GRAVITY*ubig*(drr*er(k) - drl*el(k))
            end if
         end if
      end do
   end function gm_work_face

   pure function ocean_gm_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the GM slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_gm_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%khth_u) &
               + arr_bytes(this%khth_v) &
               + arr_bytes(this%uhD) &
               + arr_bytes(this%vhD) &
               + arr_bytes(this%gm_src)
   end function ocean_gm_bytes

#include "rdb_vanished_layer.inc"

end module rdb_ocean_gm
