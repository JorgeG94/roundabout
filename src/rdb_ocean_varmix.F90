!! Spatially-varying GM/Redi lateral-diffusivity coefficients (VarMix).
module rdb_ocean_varmix
   !! VarMix capability [4]: produces the spatially-varying thickness- and
   !! tracer-diffusion coefficient face fields (`khth_u/khth_v`,
   !! `khtr_u/khtr_v`, m^2/s) that GM (`rdb_ocean_gm`) and the future Redi
   !! path consume — replacing the constant `khth` they fill from for v1.
   !!
   !! Two ingredients, both clean-room from the literature:
   !!
   !!   * A **resolution function** `Res_fn in [0,1]` (Hallberg 2013) that
   !!     turns the parameterization OFF where the deformation radius is
   !!     well resolved (`Ld >> dx`) and ON where it is not.  In the
   !!     divide-free form for power `p`:
   !!
   !!        f2_dx2   = (dx^2 + dy^2) * max(f^2, eps^2)
   !!        beta_dx2 = oneOrTwo * (dx^2 + dy^2) * |grad f|
   !!        dx_term  = f2_dx2 + cg1 * beta_dx2
   !!        Res_fn   = dx_term^(p/2) / (dx_term^(p/2) + (alpha*cg1)^p)
   !!
   !!     where `cg1` is the first baroclinic gravity-wave speed (the
   !!     `wavespeed` slot), `|grad f|` the FULL discrete 2D Coriolis-
   !!     gradient magnitude (so the equatorial `f -> 0` limit stays finite
   !!     via the beta term — no divide-by-zero), `oneOrTwo = 2` under the
   !!     Gill (1982) equatorial-Ld convention (default).  `f2_dx2` and
   !!     `beta_dx2` are STATIC (functions of metrics + f only) — precomputed
   !!     once at init with host loops; only `Res_fn` is rebuilt per thermo
   !!     step from the (slowly varying) `cg1`.  `eps = VERY_SMALL_FREQUENCY`.
   !!
   !!   * A **Visbeck (1997) / Eady (1949)** baroclinicity scaling
   !!     `KhTh += khth_slope_cff * L2 * SN` where `SN` is the Eady growth
   !!     rate (thickness-weighted vertical average of `sqrt(S^2 N^2)`) and
   !!     `L2 = visbeck_l_scale^2` (or `visbeck_l_scale^2 * areaCu` if the
   !!     knob is negative ⇒ a nondimensional scale times the cell area).
   !!
   !! `SN_u` (calc_Visbeck_coeffs thickness-weighted path):
   !!
   !!     SN_u = sum_k sqrt(S2 * N2) * H_geom / sum_k H_geom    (k = 2..nz)
   !!     H_geom = sqrt( sqrt(h(i,k)*h(i+1,k)) * sqrt(h(i,k-1)*h(i+1,k-1)) )
   !!     S2     = slope_x(i,k)^2 + (h-weighted avg of the 4 corner slope_y^2)
   !!     N2     = max(0, n2_u(i,k))
   !!     (optional S2 limit S2 = S2*S2max/(S2+S2max) only if visbeck_max_slope>0)
   !! The orthogonal slope is already folded into S2 (the h-weighted 4-corner
   !! slope_y^2 term), so SN_u is FINAL after the thickness-weighted sum — no
   !! separate SN_v combine (MOM6 calc_Visbeck_coeffs_old; the SN_v 4-corner
   !! combine belongs to the distinct calc_Eady_growth_rate_2D path).
   !!
   !! **Assembly order (Res_fn BEFORE the clamp)** per face:
   !!     Kh  = khth
   !!     Kh += khth_slope_cff * L2 * SN
   !!     Kh *= Res_fn                              (when resoln_scaled_khth)
   !!     Kh  = max(khth_min, min(Kh, khth_max))    (min only when khth_max>0)
   !! This is the **pre-CFL base** KhTh face field.  The diffusive-CFL cap is
   !! owned by GM (it holds dt + the native face metrics): GM does
   !! `KH = min(KH_CFL, base)`.  When VarMix is OFF, GM falls back to its
   !! constant `khth` ⇒ byte-identical.
   !!
   !! Bottom-up convention (k=1 bed, k=nz surface; interface K=1 bed,
   !! K=nz+1 surface — both carry zero slope from the slopes slot).
   !!
   !! Deferred (documented): EBT/SQG vertical-structure functions (KhTh stays
   !! 2D), the MEKE additive term, depth tapering, the
   !! `calc_Eady_growth_rate_2D` outcrop-cropping filter, and the OBC-aware
   !! face interpolation (interior straddle-average only for now).
   !!
   !! Default off (`&ocean_varmix_nml enable = .false.`) ⇒ `varmix_compute`
   !! is never called and GM keeps its constant `khth` ⇒ bit-identical.
   !!
   !! References: Visbeck, Marshall, Haine & Spall (1997) JPO 27, 381-402;
   !! Hallberg (2013) Ocean Modelling 72, 92-103; Eady (1949) Tellus 1, 33-52;
   !! Chelton, deSzoeke, Schlax, El Naggar & Siwertz (1998) JPO 28, 433-460;
   !! Gill (1982) "Atmosphere-Ocean Dynamics".  No source ported.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp
#else
   use rdb_constants, only: NZ_STACK_MAX, wp
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t
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

   public :: ocean_varmix_t
   public :: varmix_compute

   real(wp), parameter :: VERY_SMALL_FREQUENCY = 1.0e-17_wp
      !! Floor on f (1/s) inside `f2_dx2` so the resolution function stays
      !! finite at the equator (where f = 0; the beta term then dominates).
   real(wp), parameter :: H_SUBROUNDOFF4 = 1.0e-40_wp
      !! Tiny denominator armour for the 4-corner slope_y weighted average
      !! (mirrors the MOM6 `H_subroundoff^4` guard).

   type :: ocean_varmix_t
      !! Spatially-varying lateral-diffusivity-coefficient state.  All
      !! fields default to the inert (`enable=.false.`) configuration so an
      !! ocean run that never sets `&ocean_varmix_nml` is bit-identical.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; gate on this (never on
         !! `allocated`, which misses the GPU mapping).
      logical :: enable = .false.
         !! Master switch.  Off ⇒ `varmix_compute` is never called and GM
         !! uses the scalar `khth` ⇒ bit-identity.  Requires the slopes slot
         !! AND the wavespeed slot (loud configure invariant).
      logical :: use_visbeck = .false.
         !! Add the Visbeck/Eady `khth_slope_cff * L2 * SN` term.
      logical :: resoln_scaled_khth = .false.
         !! Multiply the assembled KhTh by `Res_fn`.
      logical :: resoln_scaled_khtr = .false.
         !! Multiply the assembled KhTr by `Res_fn`.
      logical :: gill_equatorial_ld = .true.
         !! Gill (1982) equatorial-Ld convention ⇒ `oneOrTwo = 2` in
         !! `beta_dx2` (else Pedlosky ⇒ 1).  Static; folded into the
         !! precomputed `beta_dx2_*` at init.
      logical :: interpolate_res_fn = .false.
         !! `.true.`: build `Res_fn` at centres then 2-pt-average to faces.
         !! `.false.` (default, MOM6 default): interpolate `cg1` to faces,
         !! then recompute `Res_fn` from the face `f2_dx2/beta_dx2/cg1`.
      integer :: kh_res_fn_power = 2
         !! Resolution-function power `p` (even; 2 is the production form).
      real(wp) :: kh_res_scale_coef = 1.0_wp
         !! Resolution-function `alpha` (the `(alpha*cg1)^p` denom coef).
      real(wp) :: khth = 0.0_wp
         !! Background thickness diffusivity KhTh (m^2/s) — the constant the
         !! Visbeck term and Res_fn scale.  Mirrors `&ocean_gm_nml khth`.
      real(wp) :: khtr = 0.0_wp
         !! Background tracer diffusivity KhTr (m^2/s) for the future Redi.
      real(wp) :: khth_slope_cff = 0.0_wp
         !! Visbeck coefficient `alpha_s` for the KhTh chain.
      real(wp) :: khtr_slope_cff = 0.0_wp
         !! Visbeck coefficient for the KhTr chain.
      real(wp) :: khth_min = 0.0_wp
         !! Lower clamp on the assembled KhTh (m^2/s).
      real(wp) :: khth_max = 0.0_wp
         !! Upper clamp on KhTh (m^2/s); <= 0 ⇒ no upper cap.
      real(wp) :: khtr_min = 0.0_wp
         !! Lower clamp on KhTr (m^2/s).
      real(wp) :: khtr_max = 0.0_wp
         !! Upper clamp on KhTr (m^2/s); <= 0 ⇒ no upper cap.
      real(wp) :: visbeck_l_scale = 0.0_wp
         !! Visbeck length scale L (m); if < 0, |L|^2 * areaCu/areaCv is used
         !! (a nondimensional scale times the local cell area).
      real(wp) :: visbeck_max_slope = 0.0_wp
         !! S^2 limiter scale; <= 0 ⇒ no S^2 limit.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- Static grid terms (precomputed once at init) ----
      real(wp), allocatable :: f2_dx2_u(:, :)
         !! `(dx^2+dy^2) max(f^2,eps^2)` at u-faces, `(nx+1,ny)`.  Static.
      real(wp), allocatable :: f2_dx2_v(:, :)
         !! Same at v-faces, `(nx,ny+1)`.
      real(wp), allocatable :: beta_dx2_u(:, :)
         !! `oneOrTwo*(dx^2+dy^2)*|grad f|` at u-faces, `(nx+1,ny)`.  Static.
      real(wp), allocatable :: beta_dx2_v(:, :)
         !! Same at v-faces, `(nx,ny+1)`.
      real(wp), allocatable :: l2_u(:, :)
         !! Visbeck `L^2` at u-faces (m^2), `(nx+1,ny)`.  Static.
      real(wp), allocatable :: l2_v(:, :)
         !! Visbeck `L^2` at v-faces (m^2), `(nx,ny+1)`.

      ! ---- Per-step diagnostics + outputs ----
      real(wp), allocatable :: res_fn_u(:, :)
         !! Resolution function at u-faces (nondim, [0,1]), `(nx+1,ny)`.
      real(wp), allocatable :: res_fn_v(:, :)
         !! Resolution function at v-faces, `(nx,ny+1)`.
      real(wp), allocatable :: sn_u(:, :)
         !! Eady growth rate `S*N` at u-faces (1/s), `(nx+1,ny)`.
      real(wp), allocatable :: sn_v(:, :)
         !! Eady growth rate at v-faces, `(nx,ny+1)`.
      real(wp), allocatable :: khth_u(:, :)
         !! Pre-CFL base thickness diffusivity at u-faces (m^2/s),
         !! `(nx+1,ny)` — threaded into GM as the optional external base.
      real(wp), allocatable :: khth_v(:, :)
         !! Pre-CFL base KhTh at v-faces, `(nx,ny+1)`.
      real(wp), allocatable :: khtr_u(:, :)
         !! Pre-CFL base tracer diffusivity at u-faces (m^2/s), `(nx+1,ny)`
         !! — consumed by the future Redi path.
      real(wp), allocatable :: khtr_v(:, :)
         !! Pre-CFL base KhTr at v-faces, `(nx,ny+1)`.
   contains
      procedure, non_overridable :: init => ocean_varmix_init
      procedure, non_overridable :: destroy => ocean_varmix_destroy
      procedure, non_overridable :: build_static => ocean_varmix_build_static
      procedure, non_overridable :: enter_data => ocean_varmix_enter_data
      procedure, non_overridable :: exit_data => ocean_varmix_exit_data
      procedure, non_overridable :: bytes => ocean_varmix_bytes
   end type ocean_varmix_t

contains

   subroutine ocean_varmix_init(this, grid, nz_ml)
      !! Allocate the static grid terms, the per-step diagnostics, and the
      !! KhTh/KhTr base face fields.  Always allocates (configure runs after
      !! init); the static `f2_dx2_*` / `beta_dx2_*` / `l2_*` are filled by
      !! `build_static` once the metrics + f_centre are known.  Setup uses
      !! plain host allocation (no `do concurrent` before enter_data).
      class(ocean_varmix_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml
      if (nz < 1) nz = 1
      if (nz > NZ_STACK_MAX) then
         error stop "ocean_varmix_init: nz_ml exceeds NZ_STACK_MAX "// &
            "(raise NZ_STACK_MAX in rdb_constants)"
      end if
      this%nx_total = nx
      this%ny_total = ny
      this%nz_ml = nz

      allocate (this%f2_dx2_u(nx + 1, ny), source=0.0_wp)
      allocate (this%f2_dx2_v(nx, ny + 1), source=0.0_wp)
      allocate (this%beta_dx2_u(nx + 1, ny), source=0.0_wp)
      allocate (this%beta_dx2_v(nx, ny + 1), source=0.0_wp)
      allocate (this%l2_u(nx + 1, ny), source=0.0_wp)
      allocate (this%l2_v(nx, ny + 1), source=0.0_wp)
      allocate (this%res_fn_u(nx + 1, ny), source=0.0_wp)
      allocate (this%res_fn_v(nx, ny + 1), source=0.0_wp)
      allocate (this%sn_u(nx + 1, ny), source=0.0_wp)
      allocate (this%sn_v(nx, ny + 1), source=0.0_wp)
      allocate (this%khth_u(nx + 1, ny), source=0.0_wp)
      allocate (this%khth_v(nx, ny + 1), source=0.0_wp)
      allocate (this%khtr_u(nx + 1, ny), source=0.0_wp)
      allocate (this%khtr_v(nx, ny + 1), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_varmix_init

   subroutine ocean_varmix_build_static(this, metrics, f_centre)
      !! Fill the static `f2_dx2_*`, `beta_dx2_*`, and `l2_*` face fields
      !! from the (curvilinear) metrics + the cell-centre Coriolis magnitude
      !! `f_centre`.  Called ONCE after init + configure + metrics fill (so
      !! `oneOrTwo`, `visbeck_l_scale`, and the device-resident metrics are
      !! known), BEFORE the device map.  HOST loops only — these are static
      !! (functions of geometry + planetary f) and never change with time.
      !!
      !! `|grad f|` is the FULL discrete 2D Coriolis-gradient magnitude built
      !! from `f_centre` (cell-centre) differences interpolated to each face;
      !! on an f-plane it is identically 0 ⇒ beta_dx2 = 0 (and the resolution
      !! function reduces to the midlatitude `f^2`-only form).
      class(ocean_varmix_t), intent(inout) :: this
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_centre(this%nx_total, this%ny_total)
      integer :: nx, ny, i, j
      real(wp) :: one_or_two, d2, fu, fv
      real(wp) :: dfdx, dfdy

      if (.not. this%is_init) return
      nx = this%nx_total
      ny = this%ny_total
      one_or_two = 1.0_wp
      if (this%gill_equatorial_ld) one_or_two = 2.0_wp

      ! ---- u-faces (i = 2..nx interior; edges left 0 ⇒ Res_fn = 0 there). ----
      do j = 1, ny
         do i = 2, nx
            d2 = metrics%dxCu(i, j)*metrics%dxCu(i, j) + &
                 metrics%dyCu(i, j)*metrics%dyCu(i, j)
            ! f at the u-face = straddle-average of the two centre f's.
            fu = 0.5_wp*(f_centre(i - 1, j) + f_centre(i, j))
            this%f2_dx2_u(i, j) = d2*max(fu*fu, &
                                         VERY_SMALL_FREQUENCY*VERY_SMALL_FREQUENCY)
            ! |grad f| at the u-face: df/dx across the face, df/dy averaged
            ! from the two straddling centre rows.
            dfdx = (f_centre(i, j) - f_centre(i - 1, j))*metrics%idxCu(i, j)
            dfdy = varmix_dfdy_uface(f_centre, nx, ny, i, j, metrics%idyCu(i, j))
            this%beta_dx2_u(i, j) = one_or_two*d2*sqrt(dfdx*dfdx + dfdy*dfdy)
            if (this%visbeck_l_scale < 0.0_wp) then
               this%l2_u(i, j) = this%visbeck_l_scale*this%visbeck_l_scale* &
                                 metrics%areaCu(i, j)
            else
               this%l2_u(i, j) = this%visbeck_l_scale*this%visbeck_l_scale
            end if
         end do
      end do

      ! ---- v-faces (j = 2..ny interior). ----
      do j = 2, ny
         do i = 1, nx
            d2 = metrics%dxCv(i, j)*metrics%dxCv(i, j) + &
                 metrics%dyCv(i, j)*metrics%dyCv(i, j)
            fv = 0.5_wp*(f_centre(i, j - 1) + f_centre(i, j))
            this%f2_dx2_v(i, j) = d2*max(fv*fv, &
                                         VERY_SMALL_FREQUENCY*VERY_SMALL_FREQUENCY)
            dfdy = (f_centre(i, j) - f_centre(i, j - 1))*metrics%idyCv(i, j)
            dfdx = varmix_dfdx_vface(f_centre, nx, ny, i, j, metrics%idxCv(i, j))
            this%beta_dx2_v(i, j) = one_or_two*d2*sqrt(dfdx*dfdx + dfdy*dfdy)
            if (this%visbeck_l_scale < 0.0_wp) then
               this%l2_v(i, j) = this%visbeck_l_scale*this%visbeck_l_scale* &
                                 metrics%areaCv(i, j)
            else
               this%l2_v(i, j) = this%visbeck_l_scale*this%visbeck_l_scale
            end if
         end do
      end do
   end subroutine ocean_varmix_build_static

   pure function varmix_dfdy_uface(f_centre, nx, ny, i, j, idy) result(dfdy)
      !! Cross-face df/dy at a u-face (i interior): average of the two
      !! centred y-derivatives in the west (i-1) and east (i) columns,
      !! clamped at the j edges (one-sided / zero there).
      integer, intent(in) :: nx, ny, i, j
      real(wp), intent(in) :: f_centre(nx, ny)
      real(wp), intent(in) :: idy
      real(wp) :: dfdy, dw, de
      integer :: jm, jp
      jm = max(j - 1, 1)
      jp = min(j + 1, ny)
      if (jp == jm) then
         dfdy = 0.0_wp
         return
      end if
      dw = (f_centre(i - 1, jp) - f_centre(i - 1, jm))/real(jp - jm, wp)*idy
      de = (f_centre(i, jp) - f_centre(i, jm))/real(jp - jm, wp)*idy
      dfdy = 0.5_wp*(dw + de)
   end function varmix_dfdy_uface

   pure function varmix_dfdx_vface(f_centre, nx, ny, i, j, idx) result(dfdx)
      !! Cross-face df/dx at a v-face (j interior): average of the centred
      !! x-derivatives in the south (j-1) and north (j) rows.
      integer, intent(in) :: nx, ny, i, j
      real(wp), intent(in) :: f_centre(nx, ny)
      real(wp), intent(in) :: idx
      real(wp) :: dfdx, ds, dn
      integer :: im, ip
      im = max(i - 1, 1)
      ip = min(i + 1, nx)
      if (ip == im) then
         dfdx = 0.0_wp
         return
      end if
      ds = (f_centre(ip, j - 1) - f_centre(im, j - 1))/real(ip - im, wp)*idx
      dn = (f_centre(ip, j) - f_centre(im, j))/real(ip - im, wp)*idx
      dfdx = 0.5_wp*(ds + dn)
   end function varmix_dfdx_vface

   subroutine ocean_varmix_destroy(this)
      class(ocean_varmix_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%f2_dx2_u)) deallocate (this%f2_dx2_u)
      if (allocated(this%f2_dx2_v)) deallocate (this%f2_dx2_v)
      if (allocated(this%beta_dx2_u)) deallocate (this%beta_dx2_u)
      if (allocated(this%beta_dx2_v)) deallocate (this%beta_dx2_v)
      if (allocated(this%l2_u)) deallocate (this%l2_u)
      if (allocated(this%l2_v)) deallocate (this%l2_v)
      if (allocated(this%res_fn_u)) deallocate (this%res_fn_u)
      if (allocated(this%res_fn_v)) deallocate (this%res_fn_v)
      if (allocated(this%sn_u)) deallocate (this%sn_u)
      if (allocated(this%sn_v)) deallocate (this%sn_v)
      if (allocated(this%khth_u)) deallocate (this%khth_u)
      if (allocated(this%khth_v)) deallocate (this%khth_v)
      if (allocated(this%khtr_u)) deallocate (this%khtr_u)
      if (allocated(this%khtr_v)) deallocate (this%khtr_v)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_varmix_destroy

   subroutine ocean_varmix_enter_data(this)
      ! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule:
      ! the slot header + scalars ride the orchestrator root copyin(state)).
      class(ocean_varmix_t), intent(inout) :: this
      select type (this)
      type is (ocean_varmix_t)
         call ocean_varmix_enter_data_impl(this)
      end select
   end subroutine ocean_varmix_enter_data

   subroutine ocean_varmix_enter_data_impl(this)
      type(ocean_varmix_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%f2_dx2_u, this%f2_dx2_v)
      !$acc enter data copyin(this%beta_dx2_u, this%beta_dx2_v)
      !$acc enter data copyin(this%l2_u, this%l2_v)
      !$acc enter data copyin(this%res_fn_u, this%res_fn_v)
      !$acc enter data copyin(this%sn_u, this%sn_v)
      !$acc enter data copyin(this%khth_u, this%khth_v)
      !$acc enter data copyin(this%khtr_u, this%khtr_v)
   end subroutine ocean_varmix_enter_data_impl

   subroutine ocean_varmix_exit_data(this)
      class(ocean_varmix_t), intent(inout) :: this
      select type (this)
      type is (ocean_varmix_t)
         call ocean_varmix_exit_data_impl(this)
      end select
   end subroutine ocean_varmix_exit_data

   subroutine ocean_varmix_exit_data_impl(this)
      type(ocean_varmix_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%khtr_u, this%khtr_v)
      !$acc exit data delete(this%khth_u, this%khth_v)
      !$acc exit data delete(this%sn_u, this%sn_v)
      !$acc exit data delete(this%res_fn_u, this%res_fn_v)
      !$acc exit data delete(this%l2_u, this%l2_v)
      !$acc exit data delete(this%beta_dx2_u, this%beta_dx2_v)
      !$acc exit data delete(this%f2_dx2_u, this%f2_dx2_v)
   end subroutine ocean_varmix_exit_data_impl

   ! =================================================================
   ! Compute entry — fill the per-step Res_fn, SN, and the assembled
   ! KhTh/KhTr base face fields.  Run once per outer step at THERMO
   ! cadence, BEFORE gm_compute_transports.
   ! =================================================================

   subroutine varmix_compute(grid, metrics, this, slopes, wavespeed, ms)
      !! Fill `res_fn_*`, `sn_*`, `khth_*`, `khtr_*` (the pre-CFL base face
      !! coefficients) from the static grid terms, the slopes/N^2 slot, and
      !! the first-mode wave speed `cg1`.  No-op when disabled / the deps are
      !! absent.  The CFL cap is applied downstream by GM.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_varmix_t), intent(inout) :: this
      type(ocean_slopes_t), intent(in) :: slopes
      type(ocean_wave_speed_t), intent(in) :: wavespeed
      type(multilayer_state_t), intent(in) :: ms
      integer :: nx, ny, nz
      logical :: do_visbeck

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. slopes%is_init) return
      if (.not. wavespeed%is_init) return
      if (.not. allocated(wavespeed%cg1)) return
      if (.not. allocated(ms%h_layer)) return
      if (.not. allocated(slopes%slope_x)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      if (this%nz_ml /= nz) return
      if (slopes%nz_ml /= nz) return

      ! Eady SN only when a Visbeck coefficient is on (item 8).
      do_visbeck = this%use_visbeck .and. &
                   (this%khth_slope_cff > 0.0_wp .or. this%khtr_slope_cff > 0.0_wp)

      call varmix_compute_impl(nx, ny, nz, this%kh_res_fn_power, &
                               this%kh_res_scale_coef, this%resoln_scaled_khth, &
                               this%resoln_scaled_khtr, this%interpolate_res_fn, &
                               do_visbeck, this%visbeck_max_slope, &
                               this%khth, this%khtr, this%khth_slope_cff, &
                               this%khtr_slope_cff, this%khth_min, this%khth_max, &
                               this%khtr_min, this%khtr_max, &
                               this%f2_dx2_u, this%f2_dx2_v, &
                               this%beta_dx2_u, this%beta_dx2_v, &
                               this%l2_u, this%l2_v, &
                               wavespeed%cg1, ms%h_layer, &
                               slopes%slope_x, slopes%slope_y, &
                               slopes%n2_u, slopes%n2_v, &
                               this%res_fn_u, this%res_fn_v, &
                               this%sn_u, this%sn_v, &
                               this%khth_u, this%khth_v, &
                               this%khtr_u, this%khtr_v)
   end subroutine varmix_compute

   subroutine varmix_compute_impl(nx, ny, nz, p, alpha, resoln_khth, &
                                  resoln_khtr, interp_res, do_visbeck, s2max, &
                                  khth, khtr, khth_cff, khtr_cff, khth_min, &
                                  khth_max, khtr_min, khtr_max, &
                                  f2_dx2_u, f2_dx2_v, beta_dx2_u, beta_dx2_v, &
                                  l2_u, l2_v, cg1, h_layer, slope_x, slope_y, &
                                  n2_u, n2_v, res_fn_u, res_fn_v, sn_u, sn_v, &
                                  khth_u, khth_v, khtr_u, khtr_v)
      !! Flat-impl VarMix kernel (explicit-shape; NVHPC descriptor-walk-free).
      !! Three phases: (1) Res_fn at faces (cg1 interpolated to faces or the
      !! centre-Res_fn averaged, per `interp_res`), (2) Eady SN at u/v faces
      !! (thickness-weighted column reductions with the orthogonal slope folded
      !! into S^2, scalar accumulators — SN is final, no SN_v combine), (3) the
      !! assembly (Visbeck addend, Res_fn scale, clamp) into the KhTh/KhTr base.
      integer, intent(in) :: nx, ny, nz, p
      real(wp), intent(in) :: alpha, s2max
      logical, intent(in) :: resoln_khth, resoln_khtr, interp_res, do_visbeck
      real(wp), intent(in) :: khth, khtr, khth_cff, khtr_cff
      real(wp), intent(in) :: khth_min, khth_max, khtr_min, khtr_max
      real(wp), intent(in) :: f2_dx2_u(nx + 1, ny)
      real(wp), intent(in) :: f2_dx2_v(nx, ny + 1)
      real(wp), intent(in) :: beta_dx2_u(nx + 1, ny)
      real(wp), intent(in) :: beta_dx2_v(nx, ny + 1)
      real(wp), intent(in) :: l2_u(nx + 1, ny)
      real(wp), intent(in) :: l2_v(nx, ny + 1)
      real(wp), intent(in) :: cg1(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(out) :: res_fn_u(nx + 1, ny)
      real(wp), intent(out) :: res_fn_v(nx, ny + 1)
      real(wp), intent(out) :: sn_u(nx + 1, ny)
      real(wp), intent(out) :: sn_v(nx, ny + 1)
      real(wp), intent(out) :: khth_u(nx + 1, ny)
      real(wp), intent(out) :: khth_v(nx, ny + 1)
      real(wp), intent(out) :: khtr_u(nx + 1, ny)
      real(wp), intent(out) :: khtr_v(nx, ny + 1)

      integer :: i, j
      real(wp) :: cg1u, cg1v

      ! ---- 1. Resolution function at faces. ----
      do concurrent(j=1:ny, i=1:nx + 1) local(cg1u)
         res_fn_u(i, j) = 0.0_wp
         if (i >= 2 .and. i <= nx) then
            if (interp_res) then
               ! Centre Res_fn (own + west) then 2-pt average.
               res_fn_u(i, j) = 0.5_wp*( &
                                varmix_res_fn(f2_dx2_u(i, j), beta_dx2_u(i, j), &
                                              cg1(i - 1, j), alpha, p) + &
                                varmix_res_fn(f2_dx2_u(i, j), beta_dx2_u(i, j), &
                                              cg1(i, j), alpha, p))
            else
               cg1u = 0.5_wp*(cg1(i - 1, j) + cg1(i, j))
               res_fn_u(i, j) = varmix_res_fn(f2_dx2_u(i, j), beta_dx2_u(i, j), &
                                              cg1u, alpha, p)
            end if
         end if
      end do
      do concurrent(j=1:ny + 1, i=1:nx) local(cg1v)
         res_fn_v(i, j) = 0.0_wp
         if (j >= 2 .and. j <= ny) then
            if (interp_res) then
               res_fn_v(i, j) = 0.5_wp*( &
                                varmix_res_fn(f2_dx2_v(i, j), beta_dx2_v(i, j), &
                                              cg1(i, j - 1), alpha, p) + &
                                varmix_res_fn(f2_dx2_v(i, j), beta_dx2_v(i, j), &
                                              cg1(i, j), alpha, p))
            else
               cg1v = 0.5_wp*(cg1(i, j - 1) + cg1(i, j))
               res_fn_v(i, j) = varmix_res_fn(f2_dx2_v(i, j), beta_dx2_v(i, j), &
                                              cg1v, alpha, p)
            end if
         end if
      end do

      ! ---- 2. Eady SN (thickness-weighted) at u/v faces.  These store the
      !         FINAL SN: S^2 = slope_x^2 + (h-weighted 4-corner slope_y^2) at
      !         u (mirror at v), so the orthogonal slope is already included —
      !         SN_u/SN_v are complete after the thickness-weighted sum (MOM6
      !         calc_Visbeck_coeffs_old; no separate SN_v combine). ----
      call varmix_sn_u(nx, ny, nz, do_visbeck, s2max, h_layer, slope_x, &
                       slope_y, n2_u, sn_u)
      call varmix_sn_v(nx, ny, nz, do_visbeck, s2max, h_layer, slope_x, &
                       slope_y, n2_v, sn_v)

      ! ---- 3. Assembly (Visbeck addend, Res_fn scale, clamp) into the pre-CFL
      !         base KhTh/KhTr face fields, consuming the final SN_u/SN_v.
      ! Assembly: SN_u/SN_v ALREADY carry the orthogonal slope (h4-weighted
      ! into S^2 in varmix_sn_*, per MOM6 calc_Visbeck_coeffs_old), so they
      ! are the FINAL Eady growth rate — NO extra 4-corner SN_v combine (that
      ! belongs to the separate calc_Eady_growth_rate_2D path and here would
      ! double-count the orthogonal slope).
      do concurrent(j=1:ny, i=1:nx + 1)
         khth_u(i, j) = varmix_assemble(khth, khth_cff, l2_u(i, j), sn_u(i, j), &
                                        res_fn_u(i, j), resoln_khth, khth_min, &
                                        khth_max, do_visbeck)
         khtr_u(i, j) = varmix_assemble(khtr, khtr_cff, l2_u(i, j), sn_u(i, j), &
                                        res_fn_u(i, j), resoln_khtr, khtr_min, &
                                        khtr_max, do_visbeck)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         khth_v(i, j) = varmix_assemble(khth, khth_cff, l2_v(i, j), sn_v(i, j), &
                                        res_fn_v(i, j), resoln_khth, khth_min, &
                                        khth_max, do_visbeck)
         khtr_v(i, j) = varmix_assemble(khtr, khtr_cff, l2_v(i, j), sn_v(i, j), &
                                        res_fn_v(i, j), resoln_khtr, khtr_min, &
                                        khtr_max, do_visbeck)
      end do
   end subroutine varmix_compute_impl

   pure function varmix_res_fn(f2_dx2, beta_dx2, cg1, alpha, p) result(r)
      !$acc routine seq
      !! Divide-free resolution function for power `p` (even).  p=2:
      !! `dx_term/(dx_term + (alpha*cg1)^2)`; general even p:
      !! `dx_term^(p/2)/(dx_term^(p/2) + (alpha*cg1)^p)`.  `dx_term =
      !! f2_dx2 + cg1*beta_dx2`.  -> 1 where unresolved, -> 0 where Ld>>dx.
      real(wp), intent(in) :: f2_dx2, beta_dx2, cg1, alpha
      integer, intent(in) :: p
      real(wp) :: r, dx_term, num, den_add
      integer :: ph
      dx_term = f2_dx2 + cg1*beta_dx2
      if (p == 2) then
         num = dx_term
         den_add = (alpha*cg1)*(alpha*cg1)
      else
         ph = p/2
         num = dx_term**ph
         den_add = (alpha*cg1)**p
      end if
      r = num/(num + den_add)
   end function varmix_res_fn

   pure function varmix_assemble(kh_bg, cff, l2, sn, res_fn, resoln, &
                                 kh_min, kh_max, do_visbeck) result(kh)
      !$acc routine seq
      !! Assembly chain in the load-bearing order: background + Visbeck
      !! addend, THEN Res_fn scale, THEN clamp.  `kh_max <= 0` ⇒ no upper cap.
      real(wp), intent(in) :: kh_bg, cff, l2, sn, res_fn, kh_min, kh_max
      logical, intent(in) :: resoln, do_visbeck
      real(wp) :: kh
      kh = kh_bg
      if (do_visbeck) kh = kh + cff*l2*sn
      if (resoln) kh = kh*res_fn
      if (kh_max > 0.0_wp) then
         kh = max(kh_min, min(kh, kh_max))
      else
         kh = max(kh_min, kh)
      end if
   end function varmix_assemble

   pure subroutine varmix_sn_u(nx, ny, nz, do_visbeck, s2max, h_layer, &
                               slope_x, slope_y, n2_u, sn_u)
      !! Thickness-weighted Eady growth rate at u-faces (own component).
      !! Interior u-face (i=2..nx) pairs centre columns iw=i-1 (west) and i
      !! (east).  Interior interfaces K=2..nz; `H_geom = sqrt(sqrt(h_iw,k *
      !! h_i,k) * sqrt(h_iw,k-1 * h_i,k-1))`.  S2 = slope_x^2 + the four
      !! corner slope_y^2 h-weighted to the u-face; S2 optionally limited.
      !! `SN_u = sum sqrt(S2*N2)*H_geom / sum H_geom`.
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: do_visbeck
      real(wp), intent(in) :: s2max
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_u(nx + 1, ny, nz + 1)
      real(wp), intent(out) :: sn_u(nx + 1, ny)

      integer :: i, j, k, iw, jm, jp
      real(wp) :: hgeom, hdn, hup, s2, n2, sn_acc, h_acc
      real(wp) :: wsw, wse, wnw, wne, sy2, denom

      do concurrent(j=1:ny, i=1:nx + 1) &
         local(k, iw, jm, jp, hgeom, hdn, hup, s2, n2, sn_acc, h_acc, &
               wsw, wse, wnw, wne, sy2, denom)
         sn_u(i, j) = 0.0_wp
         if (do_visbeck .and. i >= 2 .and. i <= nx) then
            iw = i - 1
            jm = max(1, j - 1)        ! south cell-row (array-edge clamp)
            jp = min(ny, j + 1)       ! north cell-row (array-edge clamp)
            sn_acc = 0.0_wp
            h_acc = 0.0_wp
            do k = 2, nz       ! interior interface index
               hdn = sqrt(max(h_layer(iw, j, k)*h_layer(i, j, k), 0.0_wp))
               hup = sqrt(max(h_layer(iw, j, k - 1)*h_layer(i, j, k - 1), 0.0_wp))
               hgeom = sqrt(hdn*hup)
               ! 4 corner slope_y values around the u-face, weighted by the
               ! MOM6 h4_v product co-located with each corner's slope_y
               ! v-point: the 4 thicknesses straddling that v-face (the two
               ! cell-rows it separates) at the two layers (k, k-1) the
               ! interface separates.  Under uniform h all four are equal ⇒
               ! bit-identical to a single-thickness weight.
               wnw = (h_layer(iw, j, k)*h_layer(iw, jp, k)) &
                     *(h_layer(iw, j, k - 1)*h_layer(iw, jp, k - 1))
               wne = (h_layer(i, j, k)*h_layer(i, jp, k)) &
                     *(h_layer(i, j, k - 1)*h_layer(i, jp, k - 1))
               wsw = (h_layer(iw, jm, k)*h_layer(iw, j, k)) &
                     *(h_layer(iw, jm, k - 1)*h_layer(iw, j, k - 1))
               wse = (h_layer(i, jm, k)*h_layer(i, j, k)) &
                     *(h_layer(i, jm, k - 1)*h_layer(i, j, k - 1))
               denom = ((wse + wnw) + (wne + wsw)) + H_SUBROUNDOFF4
               sy2 = (((wnw*slope_y(iw, j + 1, k)*slope_y(iw, j + 1, k)) + &
                       (wse*slope_y(i, j, k)*slope_y(i, j, k))) + &
                      ((wne*slope_y(i, j + 1, k)*slope_y(i, j + 1, k)) + &
                       (wsw*slope_y(iw, j, k)*slope_y(iw, j, k))))/denom
               s2 = slope_x(i, j, k)*slope_x(i, j, k) + sy2
               if (s2max > 0.0_wp) s2 = s2*s2max/(s2 + s2max)
               n2 = max(n2_u(i, j, k), 0.0_wp)
               sn_acc = sn_acc + sqrt(s2*n2)*hgeom
               h_acc = h_acc + hgeom
            end do
            if (h_acc > 0.0_wp) sn_u(i, j) = sn_acc/h_acc
         end if
      end do
   end subroutine varmix_sn_u

   pure subroutine varmix_sn_v(nx, ny, nz, do_visbeck, s2max, h_layer, &
                               slope_x, slope_y, n2_v, sn_v)
      !! Thickness-weighted Eady growth rate at v-faces (mirror of
      !! `varmix_sn_u`).  Interior v-face (j=2..ny) pairs js=j-1 (south) + j.
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: do_visbeck
      real(wp), intent(in) :: s2max
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(in) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(in) :: n2_v(nx, ny + 1, nz + 1)
      real(wp), intent(out) :: sn_v(nx, ny + 1)

      integer :: i, j, k, js, im, ip
      real(wp) :: hgeom, hdn, hup, s2, n2, sn_acc, h_acc
      real(wp) :: wsw, wse, wnw, wne, sx2, denom

      do concurrent(j=1:ny + 1, i=1:nx) &
         local(k, js, im, ip, hgeom, hdn, hup, s2, n2, sn_acc, h_acc, &
               wsw, wse, wnw, wne, sx2, denom)
         sn_v(i, j) = 0.0_wp
         if (do_visbeck .and. j >= 2 .and. j <= ny) then
            js = j - 1
            im = max(1, i - 1)        ! west cell-column (array-edge clamp)
            ip = min(nx, i + 1)       ! east cell-column (array-edge clamp)
            sn_acc = 0.0_wp
            h_acc = 0.0_wp
            do k = 2, nz
               hdn = sqrt(max(h_layer(i, js, k)*h_layer(i, j, k), 0.0_wp))
               hup = sqrt(max(h_layer(i, js, k - 1)*h_layer(i, j, k - 1), 0.0_wp))
               hgeom = sqrt(hdn*hup)
               ! 4 corner slope_x values around the v-face, weighted by the
               ! MOM6 h4_u product co-located with each corner's slope_x
               ! u-point: the 4 thicknesses straddling that u-face (the two
               ! cell-cols it separates) at the two layers (k, k-1) the
               ! interface separates.  Uniform h ⇒ bit-identical.
               wse = (h_layer(i, js, k)*h_layer(ip, js, k)) &
                     *(h_layer(i, js, k - 1)*h_layer(ip, js, k - 1))
               wnw = (h_layer(im, j, k)*h_layer(i, j, k)) &
                     *(h_layer(im, j, k - 1)*h_layer(i, j, k - 1))
               wne = (h_layer(i, j, k)*h_layer(ip, j, k)) &
                     *(h_layer(i, j, k - 1)*h_layer(ip, j, k - 1))
               wsw = (h_layer(im, js, k)*h_layer(i, js, k)) &
                     *(h_layer(im, js, k - 1)*h_layer(i, js, k - 1))
               denom = ((wse + wnw) + (wne + wsw)) + H_SUBROUNDOFF4
               sx2 = (((wse*slope_x(i + 1, js, k)*slope_x(i + 1, js, k)) + &
                       (wnw*slope_x(i, j, k)*slope_x(i, j, k))) + &
                      ((wne*slope_x(i + 1, j, k)*slope_x(i + 1, j, k)) + &
                       (wsw*slope_x(i, js, k)*slope_x(i, js, k))))/denom
               s2 = slope_y(i, j, k)*slope_y(i, j, k) + sx2
               if (s2max > 0.0_wp) s2 = s2*s2max/(s2 + s2max)
               n2 = max(n2_v(i, j, k), 0.0_wp)
               sn_acc = sn_acc + sqrt(s2*n2)*hgeom
               h_acc = h_acc + hgeom
            end do
            if (h_acc > 0.0_wp) sn_v(i, j) = sn_acc/h_acc
         end if
      end do
   end subroutine varmix_sn_v

   pure function ocean_varmix_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the VarMix slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_varmix_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f2_dx2_u) &
               + arr_bytes(this%f2_dx2_v) &
               + arr_bytes(this%beta_dx2_u) &
               + arr_bytes(this%beta_dx2_v) &
               + arr_bytes(this%l2_u) &
               + arr_bytes(this%l2_v) &
               + arr_bytes(this%res_fn_u) &
               + arr_bytes(this%res_fn_v) &
               + arr_bytes(this%sn_u) &
               + arr_bytes(this%sn_v) &
               + arr_bytes(this%khth_u) &
               + arr_bytes(this%khth_v) &
               + arr_bytes(this%khtr_u) &
               + arr_bytes(this%khtr_v)
   end function ocean_varmix_bytes

end module rdb_ocean_varmix
