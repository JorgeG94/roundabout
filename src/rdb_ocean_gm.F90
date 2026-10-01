!! Gent-McWilliams thickness-diffusion (eddy-induced bolus transport) slot.
module rdb_ocean_gm
   !! Gent & McWilliams (1990) / Griffies (1998) skew-flux thickness
   !! diffusion.  The eddy-induced ("bolus") overturning is realized as a
   !! layer THICKNESS flux (never an explicit velocity), folded into the
   !! continuity / windowed-PPM tracer path like the Fox-Kemper MLE slot.
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
   !! Default off (`&ocean_gm_nml enable=.false.`) => slot allocated but
   !! `gm_compute_transports` no-ops => bit-identical.
   !! Refs: Gent & McWilliams (1990); Gent et al. (1995); Griffies (1998).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS, GRAVITY
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, H_VANISHED, H_DIV_EPS, GRAVITY
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t
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

   public :: ocean_gm_t
   public :: gm_compute_transports
   public :: gm_fold_x, gm_fold_y

   type :: ocean_gm_t
      !! Gent-McWilliams thickness-diffusion state.  All fields default to
      !! the inert (`enable=.false.`) configuration so an ocean run that
      !! never sets `&ocean_gm_nml` is bit-identical.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; gate on this (never on
         !! `allocated`, which misses the GPU mapping).
      logical :: enable = .false.
         !! Master switch.  Off => `gm_compute_transports` is a no-op and
         !! the folds add nothing => bit-identity preserved.  Requires the
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
      ! k=1 bed, k=nz surface.  Added to the continuity mass fluxes before
      ! the divergence; Sum_k over each face is ~0 (closed overturning).
      real(wp), allocatable :: uhD(:, :, :)
         !! GM x-face transport (m^3/s), `(nx+1,ny,nz)`.  Filled once per
         !! outer step at thermo cadence, folded into continuity every
         !! dynamics call via `gm_fold_x`.
      real(wp), allocatable :: vhD(:, :, :)
         !! GM y-face transport (m^3/s), `(nx,ny+1,nz)`.

      ! ---- 2D GM PE-release source (for MEKE) ----
      real(wp), allocatable :: gm_src(:, :)
         !! Potential-energy release `-1/4 Sum_k rho0 KH Slope^2 N^2 h`
         !! (W/m^2-ish; >= 0 for a stable tilted column), `(nx,ny)`.
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
                                    khth_ext_u, khth_ext_v)
      !! Fill `uhD`/`vhD` (m^3/s) with the GM bolus thickness transport and
      !! `gm_src` with the PE release.  Run once per outer step at THERMO
      !! cadence, AFTER `ocean_slopes_compute` and BEFORE the continuity
      !! divergence.  No-op when disabled, uninitialised, or the slopes
      !! slot is absent/disabled.
      !!
      !! `khth_ext_u/khth_ext_v` (optional): spatially-varying PRE-CFL base
      !! KhTh face field from VarMix.  When present the per-face CFL clamp
      !! uses it instead of the scalar `khth`; when absent, falls back to
      !! the constant `khth` (byte-identical to the pre-VarMix path).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_gm_t), intent(inout) :: this
      type(ocean_slopes_t), intent(in) :: slopes
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: khth_ext_u(:, :)
      real(wp), intent(in), optional :: khth_ext_v(:, :)

      integer :: nx, ny, nz
      logical :: use_ext

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

      ! Explicit-shape flat-impl: pass the top-level allocatables so NVHPC
      ! does not descriptor-walk per launch.  When the external VarMix base
      ! is absent, pass `this%khth_u/khth_v` as a harmless placeholder for
      ! the `khth_ext_*` dummy and gate it off with `use_ext`.
      if (use_ext) then
         call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                              this%khth_slope_max, this%rho0, &
                              use_ext, &
                              metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                              metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                              metrics%areaT, metrics%wet_u, &
                              metrics%wet_v, ms%h_layer, &
                              slopes%slope_x, slopes%slope_y, &
                              slopes%n2_u, slopes%n2_v, &
                              khth_ext_u, khth_ext_v, &
                              this%khth_u, this%khth_v, &
                              this%uhD, this%vhD, this%gm_src)
      else
         call gm_compute_impl(nx, ny, nz, dt, this%khth, this%khth_max_cfl, &
                              this%khth_slope_max, this%rho0, &
                              use_ext, &
                              metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                              metrics%idyCv, metrics%idyCu, metrics%idxCv, &
                              metrics%areaT, metrics%wet_u, &
                              metrics%wet_v, ms%h_layer, &
                              slopes%slope_x, slopes%slope_y, &
                              slopes%n2_u, slopes%n2_v, &
                              this%khth_u, this%khth_v, &
                              this%khth_u, this%khth_v, &
                              this%uhD, this%vhD, this%gm_src)
      end if
   end subroutine gm_compute_transports

   subroutine gm_compute_impl(nx, ny, nz, dt, khth, khth_max_cfl, slope_max, &
                              rho0, use_ext, dy_cu, dx_cv, idxCu, &
                              idyCv, idyCu, idxCv, areaT, wet_u, wet_v, h_layer, &
                              slope_x, slope_y, &
                              n2_u, n2_v, khth_ext_u, khth_ext_v, khth_u, khth_v, &
                              uhD, vhD, gm_src)
      !! Flat-impl GM kernel.  Passes: CFL-clamp the 2D face KhTh, the
      !! u-face and v-face column recurrences into `uhD`/`vhD`, the
      !! `gm_src` PE release.  Each face's column sweep is serial in k (the
      !! `uhtot` recurrence) but faces parallelise over (i,j).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt, khth, khth_max_cfl, slope_max, rho0
      logical, intent(in) :: use_ext
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
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: khth_ext_u(nx + 1, ny)
      real(wp), intent(in) :: khth_ext_v(nx, ny + 1)
      real(wp), intent(out) :: khth_u(nx + 1, ny)
      real(wp), intent(out) :: khth_v(nx, ny + 1)
      real(wp), intent(out) :: uhD(nx + 1, ny, nz)
      real(wp), intent(out) :: vhD(nx, ny + 1, nz)
      real(wp), intent(out) :: gm_src(nx, ny)

      real(wp) :: i_smax2, i4dt
      integer :: i, j, k

      i_smax2 = 1.0_wp/(slope_max*slope_max)
      i4dt = 1.0_wp/(4.0_wp*dt)

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
      call gm_column_x(nx, ny, nz, i_smax2, i4dt, &
                       dy_cu, areaT, h_layer, slope_x, khth_u, uhD)

      ! ---- 3. v-face bolus transport (mirror).
      do concurrent(k=1:nz, i=1:nx)
         vhD(i, 1, k) = 0.0_wp
         vhD(i, ny + 1, k) = 0.0_wp
      end do
      call gm_column_y(nx, ny, nz, i_smax2, i4dt, &
                       dx_cv, areaT, h_layer, slope_y, khth_v, vhD)

      ! ---- 4. gm_src PE release (cell centres).
      call gm_pe_release(nx, ny, nz, slope_max, rho0, h_layer, &
                         slope_x, slope_y, n2_u, n2_v, khth_u, khth_v, gm_src)
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

   pure subroutine gm_column_x(nx, ny, nz, i_smax2, i4dt, &
                               dy_cu, areaT, h_layer, slope_x, khth_u, uhD)
      !! u-face GM streamfunction + bolus-transport column recurrence.
      !! Interior u-face (i=2..nx) pairs columns iw=i-1 (west) and i (east).
      !! Bottom-up sweep: interior interfaces Kr=2 (bed-most) -> nz
      !! (surface-most), uhtot=0 at the bed; the surface BC (Sfn=0 at
      !! Kr=nz+1) is closed after the loop by `uhD(nz)=-uhtot`, giving
      !! `Sum_k uhD=0` exactly.  Interface Kr straddles ka=Kr (above) and
      !! kb=Kr-1 (below) and fills LAYER kb; the rsum bound keys on ka, the
      !! donor `h_frac` and per-layer cap on kb.
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
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: khth_u(nx + 1, ny)
      real(wp), intent(inout) :: uhD(nx + 1, ny, nz)

      integer :: i, j, k, iw, ka, kb
      real(wp) :: havL(NZ_STACK_MAX), havR(NZ_STACK_MAX)
      real(wp) :: rsumL(NZ_STACK_MAX + 1), rsumR(NZ_STACK_MAX + 1)
      real(wp) :: uhtot, slope, s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h
      real(wp) :: h_frac_d, uhd_k, kh

      do concurrent(j=1:ny, i=2:nx) &
         local(k, iw, ka, kb, havL, havR, rsumL, rsumR, uhtot, slope, s2r, &
               sfn_unlim, sfn_safe, sfn_est, sfn_in_h, h_frac_d, uhd_k, kh)
         iw = i - 1
         kh = khth_u(i, j)

         ! Per-layer availability + cumulative rsum from the SURFACE down:
         ! rsum*(k) = Sum_{k'=k}^{nz} h_avail(k') (mass above interface k).
         rsumL(nz + 1) = 0.0_wp
         rsumR(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            havL(k) = max(i4dt*areaT(iw, j)*(h_layer(iw, j, k) - H_VANISHED), 0.0_wp)
            havR(k) = max(i4dt*areaT(i, j)*(h_layer(i, j, k) - H_VANISHED), 0.0_wp)
            rsumL(k) = rsumL(k + 1) + havL(k)
            rsumR(k) = rsumR(k + 1) + havR(k)
         end do

         ! Sweep interior interfaces bed-most (Kr=2) -> surface-most (Kr=nz).
         uhtot = 0.0_wp
         do k = 2, nz       ! k is the interface index Kr
            ka = k           ! layer above interface (surface side)
            kb = k - 1       ! layer below interface (bed side) -> uhD(kb)
            slope = slope_x(i, j, k)
            s2r = slope*slope*i_smax2
            sfn_unlim = -(kh*dy_cu(i, j))*slope
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
         ! Surface BC (Sfn=0 above layer nz): close the column so Sum=0.
         uhD(i, j, nz) = -uhtot
      end do
   end subroutine gm_column_x

   pure subroutine gm_column_y(nx, ny, nz, i_smax2, i4dt, &
                               dx_cv, areaT, h_layer, slope_y, khth_v, vhD)
      !! v-face GM column recurrence — mirror of `gm_column_x` with the
      !! v-stagger.  Interior v-face (j=2..ny) pairs columns js=j-1 (south)
      !! and j (north).  See `gm_column_x` for the MOM6 `nk_linear`
      !! divergence note.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: i_smax2, i4dt
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: khth_v(nx, ny + 1)
      real(wp), intent(inout) :: vhD(nx, ny + 1, nz)

      integer :: i, j, k, js, ka, kb
      real(wp) :: havS(NZ_STACK_MAX), havN(NZ_STACK_MAX)
      real(wp) :: rsumS(NZ_STACK_MAX + 1), rsumN(NZ_STACK_MAX + 1)
      real(wp) :: vhtot, slope, s2r, sfn_unlim, sfn_safe, sfn_est, sfn_in_h
      real(wp) :: h_frac_d, vhd_k, kh

      do concurrent(j=2:ny, i=1:nx) &
         local(k, js, ka, kb, havS, havN, rsumS, rsumN, vhtot, slope, s2r, &
               sfn_unlim, sfn_safe, sfn_est, sfn_in_h, h_frac_d, vhd_k, kh)
         js = j - 1
         kh = khth_v(i, j)

         rsumS(nz + 1) = 0.0_wp
         rsumN(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            havS(k) = max(i4dt*areaT(i, js)*(h_layer(i, js, k) - H_VANISHED), 0.0_wp)
            havN(k) = max(i4dt*areaT(i, j)*(h_layer(i, j, k) - H_VANISHED), 0.0_wp)
            rsumS(k) = rsumS(k + 1) + havS(k)
            rsumN(k) = rsumN(k + 1) + havN(k)
         end do

         vhtot = 0.0_wp
         do k = 2, nz
            ka = k
            kb = k - 1
            slope = slope_y(i, j, k)
            s2r = slope*slope*i_smax2
            sfn_unlim = -(kh*dx_cv(i, j))*slope
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
         vhD(i, j, nz) = -vhtot
      end do
   end subroutine gm_column_y

   pure function gm_clamp_slope(s, smax) result(sc)
      !$acc routine seq
      !! Clamp a slope to +/- smax (bounded slope for the PE release).
      real(wp), intent(in) :: s, smax
      real(wp) :: sc
      sc = s
      if (sc > smax) sc = smax
      if (sc < -smax) sc = -smax
   end function gm_clamp_slope

   pure subroutine gm_pe_release(nx, ny, nz, slope_max, rho0, h_layer, &
                                 slope_x, slope_y, n2_u, n2_v, &
                                 khth_u, khth_v, gm_src)
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

      integer :: i, j, k, kb
      real(wp) :: acc, fsum, sx_w, sx_e, sy_s, sy_n
      real(wp) :: n2w, n2e, n2s, n2n

      do concurrent(j=1:ny, i=1:nx) &
         local(k, kb, acc, fsum, sx_w, sx_e, sy_s, sy_n, n2w, n2e, n2s, n2n)
         acc = 0.0_wp
         do k = 2, nz                 ! interior interface index Kr
            kb = k - 1                 ! layer below the interface
            sx_w = gm_clamp_slope(slope_x(i, j, k), slope_max)
            sx_e = gm_clamp_slope(slope_x(i + 1, j, k), slope_max)
            sy_s = gm_clamp_slope(slope_y(i, j, k), slope_max)
            sy_n = gm_clamp_slope(slope_y(i, j + 1, k), slope_max)
            n2w = max(n2_u(i, j, k), 0.0_wp)
            n2e = max(n2_u(i + 1, j, k), 0.0_wp)
            n2s = max(n2_v(i, j, k), 0.0_wp)
            n2n = max(n2_v(i, j + 1, k), 0.0_wp)
            fsum = (khth_u(i, j)*sx_w*sx_w*n2w + khth_u(i + 1, j)*sx_e*sx_e*n2e) + &
                   (khth_v(i, j)*sy_s*sy_s*n2s + khth_v(i, j + 1)*sy_n*sy_n*n2n)
            acc = acc + fsum*h_layer(i, j, kb)
         end do
         gm_src(i, j) = 0.25_wp*rho0*acc
      end do
   end subroutine gm_pe_release

   ! =================================================================
   ! Fold into continuity mass fluxes (called inside the split step)
   ! =================================================================

   pure subroutine gm_fold_x(this, mass_flux_x_layer, nx1, ny, nz)
      !! Add the GM x-face bolus transport into the per-layer zonal mass
      !! flux, AFTER `continuity_zonal_flux` fills it and BEFORE the zonal
      !! tracer advect / divergence — so the augmented flux transports both
      !! h and tracers (conservative; no velocity touched).  No-op when
      !! disabled.
      type(ocean_gm_t), intent(in) :: this
      integer, intent(in) :: nx1, ny, nz
      real(wp), intent(inout) :: mass_flux_x_layer(nx1, ny, nz)
      integer :: i, j, k
      if (.not. this%is_init) return
      if (.not. this%enable) return
      do concurrent(k=1:nz, j=1:ny, i=1:nx1)
         mass_flux_x_layer(i, j, k) = mass_flux_x_layer(i, j, k) + this%uhD(i, j, k)
      end do
   end subroutine gm_fold_x

   pure subroutine gm_fold_y(this, mass_flux_y_layer, nx, ny1, nz)
      !! Add the GM y-face bolus transport into the per-layer meridional
      !! mass flux.  Mirror of `gm_fold_x`.  No-op when disabled.
      type(ocean_gm_t), intent(in) :: this
      integer, intent(in) :: nx, ny1, nz
      real(wp), intent(inout) :: mass_flux_y_layer(nx, ny1, nz)
      integer :: i, j, k
      if (.not. this%is_init) return
      if (.not. this%enable) return
      do concurrent(k=1:nz, j=1:ny1, i=1:nx)
         mass_flux_y_layer(i, j, k) = mass_flux_y_layer(i, j, k) + this%vhD(i, j, k)
      end do
   end subroutine gm_fold_y

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

end module rdb_ocean_gm
