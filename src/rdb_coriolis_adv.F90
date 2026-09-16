!! PV-conserving Coriolis-advection kernel state + barotropic step.
module rdb_coriolis_adv
   !! Holds variant flags + PV-at-corner workspace for the combined
   !! Coriolis + horizontal-momentum-advection kernel.  We follow the
   !! Sadourny (1975) energy/enstrophy-conserving form, optionally
   !! upgraded with the Hollingsworth-Källén-Arakawa correction that
   !! removes the spurious "Hollingsworth instability" on C-grids at
   !! eddy-resolving resolutions.
   !!
   !! The energy-conserving form adds the relative-vorticity flux +
   !! KE-gradient advection terms on top of the Coriolis force.  For
   !! uniform velocity fields zeta=0 and grad(KE)=0, so the kernel
   !! reduces to plain Coriolis.
   use rdb_constants, only: wp, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_ocean_porous, only: porous_narrow_3d
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: coriolis_adv_t
   public :: coriolis_adv_compute_tendencies_barotropic
   public :: coriolis_adv_apply_tendencies_barotropic
   public :: coriolis_adv_compute_tendencies
   public :: coriolis_adv_apply_tendencies
   public :: coriolis_adv_compute_tendencies_hk
   public :: coriolis_adv_compute_tendencies_sadourny_energy
   public :: parse_pv_variant
   public :: pv_variant_is_implemented
   public :: parse_pv_adv_scheme
   public :: pv_adv_scheme_is_implemented
   public :: pv_adv_required_nghost
   public :: weno3_recon, weno5_recon, weno7_recon

   ! PV-scheme variant tags.  Both `SADOURNY` and `SADOURNY_HK` ship;
   ! `AL81` is reserved.  The default below points at the safe
   ! enstrophy-conserving branch — flip to `SADOURNY_HK` for eddy-
   ! resolving runs via `&ocean_setup_nml ocean_coriolis_form`.
   integer, parameter, public :: PV_VARIANT_SADOURNY = 1
      !! Sadourny enstrophy-conserving — the default.
   integer, parameter, public :: PV_VARIANT_AL81 = 2
      !! Arakawa-Lamb 1981 PV-conserving (future; not yet wired).
   integer, parameter, public :: PV_VARIANT_SADOURNY_HK = 3
      !! Sadourny + Arakawa-Hsu (1990) "HK correction" — wider
      !! 3-corner stencil; suppresses the Hollingsworth-Källén
      !! instability that biases vanilla Sadourny at eddy-resolving
      !! resolutions.  Implemented in
      !! `coriolis_adv_compute_tendencies_hk`.
   integer, parameter, public :: PV_VARIANT_SADOURNY_ENERGY = 4
      !! Sadourny 1975 ENERGY-conserving (MOM6's `SADOURNY75_ENERGY`).
      !! Same stencil as the enstrophy form but the PV is grouped
      !! per-corner with its own v-fluxes rather than averaged:
      !!
      !!   ENSTRO (ours, default):
      !!     CAu = 0.5·avg(q_N, q_S) · avg(v_NW, v_NE, v_SW, v_SE)
      !!   ENERGY (MOM6):
      !!     CAu = 0.5·[q_N · 0.5·(v_NW + v_NE)
      !!                + q_S · 0.5·(v_SW + v_SE)]
      !!
      !! Conserves total kinetic energy; enstrophy form conserves
      !! squared vorticity instead.  ENSTRO dissipates v specifically
      !! at asymmetric WBC fronts; ENERGY preserves it (matches MOM6
      !! double_gyre at NK=1 within geostrophy on v).  Engaged via
      !! `ocean_coriolis_form = "sadourny_energy"`.
   integer, parameter, public :: PV_VARIANT_INVALID = -1
      !! Sentinel returned by `parse_pv_variant` for an unrecognised
      !! string (PR-6 fail-loud).  Distinct from `PV_VARIANT_AL81`
      !! (reserved-but-unwired): both are rejected by
      !! `pv_variant_is_implemented`, but a caller/validate_config can
      !! give a "typo" a different message than a "not-yet-implemented".

   ! PV face-interpolation scheme (orthogonal to PV_VARIANT_*): how the
   ! corner absolute vorticity is interpolated onto the velocity faces in
   ! the Sadourny enstrophy path.  `centered` (default) is the classic
   ! 2-point corner average -> bit-identical.  `weno{3,5,7}` replace it with
   ! an essentially-non-oscillatory upwind-biased reconstruction (MOM6
   ! CoriolisAdv WENO-VI: Large et al. WENO-Z weights) that sharpens PV
   ! fronts without the global dissipation the centred form leaks.  weno3
   ! (radius-2 stencil) fits nghost>=2; weno5/weno7 (radius 3/4) require
   ! nghost>=3/4 (fail-loud at configure via `pv_adv_required_nghost`).
   integer, parameter, public :: PV_ADV_CENTERED = 0
      !! 2-point corner average (default; bit-identical to pre-F1).
   integer, parameter, public :: PV_ADV_WENO3 = 1
      !! 3rd-order WENO-Z PV reconstruction (MOM6 WENOVI3RD, radius 2).
   integer, parameter, public :: PV_ADV_WENO5 = 2
      !! 5th-order WENO-Z PV reconstruction (MOM6 WENOVI5TH, radius 3 =>
      !! nghost>=3).
   integer, parameter, public :: PV_ADV_WENO7 = 3
      !! 7th-order WENO-Z PV reconstruction (MOM6 WENOVI7TH, radius 4 =>
      !! nghost>=4).
   integer, parameter, public :: PV_ADV_INVALID = -1
      !! Sentinel for an unrecognised string (fail-loud, mirrors
      !! `PV_VARIANT_INVALID`).

   real(wp), parameter :: PV_WENO_EPS_REL = 1.0e-20_wp
      !! MOM6 `fac_fn` guard threshold: a smoothness indicator `b` with
      !! `|b| <= PV_WENO_EPS_REL*tau` is treated as degenerate (that stencil
      !! dominates, factor -> PV_WENO_FAC_DEGEN).  Matches MOM6's implicit
      !! (no additive-epsilon) WENO-Z regulariser exactly.
   real(wp), parameter :: PV_WENO_FAC_DEGEN = 1.0e40_wp
      !! Degenerate-stencil nonlinear factor (MOM6's literal `1.0e40`).

   real(wp), parameter :: CORIOLIS_H_MIN_PV = 1.0e-12_wp
      !! Floor for the corner-h divide in the PV construction (shared by the
      !! energy/hk impls and the BOUND_CORIOLIS abs_vort recovery so the
      !! recomputed `h_corner` matches the pass-2 value bit-for-bit).
   real(wp), parameter :: PV_VOL_NEGLECT = 1.0e-20_wp
      !! Roundabout analogue of MOM6's `vol_neglect`
      !! (`H_subroundoff·(1e-4 m)²`): a corner VOLUME so small it is pure 1/0
      !! armor, NOT a physical floor.  Used only under `corner_h="mom6_area"`
      !! in `q = abs_vort·Area_q/(hArea_q + PV_VOL_NEGLECT)`.  Unlike the
      !! `cell_mean` path's `CORIOLIS_H_MIN_PV` thickness cap, this does NOT
      !! bound `q` at a vanishing corner — matching MOM6.
   integer, parameter, public :: CORNER_H_CELL_MEAN = 0
      !! `corner_h="cell_mean"` (default): wet-area 4-cell mean, `H_MIN_PV`-capped.
   integer, parameter, public :: CORNER_H_MOM6_AREA = 1
      !! `corner_h="mom6_area"`: MOM6 `q = abs_vort·Area_q/(hArea_q + vol_neglect)`.

   type :: coriolis_adv_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      integer :: pv_variant = PV_VARIANT_SADOURNY
         !! Active PV/Coriolis scheme variant.  Defaults to Sadourny;
         !! set to `PV_VARIANT_SADOURNY_HK` to engage the Arakawa-Hsu
         !! correction.  Driven from `&ocean_setup_nml
         !! ocean_coriolis_form` via `parse_pv_variant`.
      logical :: use_hk_correction = .false.
         !! Convenience flag mirroring `pv_variant == SADOURNY_HK`.
         !! Diagnostic only — the dispatch in `ocean_dyn` reads
         !! `pv_variant` directly.
      integer :: pv_adv_scheme = PV_ADV_CENTERED
         !! PV face-interpolation scheme (orthogonal to `pv_variant`):
         !! `PV_ADV_CENTERED` (default, 2-pt corner average -> bit-identical)
         !! or `PV_ADV_WENO{3,5,7}` (upwind-biased WENO-Z reconstruction of the
         !! corner absolute vorticity onto the faces in the Sadourny path).
         !! Driven from `&ocean_coriolis_nml pv_adv_scheme` via
         !! `parse_pv_adv_scheme`; weno5/weno7 require nghost>=3/4.
      logical :: weno_velocity_smooth = .false.
         !! `&ocean_coriolis_nml weno_velocity_smooth` (MOM6
         !! `WENO_VELOCITY_SMOOTH`, default off): when on, the WENO
         !! smoothness indicators are computed from the tangential velocity
         !! rather than the vorticity.  Off = vorticity-based (MOM6 default).
      logical :: state_fluxes = .false.
         !! `&ocean_coriolis_nml use_state_fluxes` — mass-consistent
         !! CorAdCalc.  Config guarantees `form="sadourny_energy"` +
         !! `split_scheme="pred_corr"`; the dyn driver forwards it to
         !! `coriolis_adv_compute_tendencies(use_state_fluxes=)` in the
         !! CORRECTOR stage only (the predictor keeps the recompute —
         !! its continuity has not run yet this step).
      logical :: no_slip = .false.
         !! Lateral boundary condition at coasts (spec §14 C1).
         !! `.false.` (default) = FREE-SLIP: the corner relative
         !! vorticity is multiplied by `wet_q` so a land corner
         !! contributes zero rel-vort (tangential velocity free at the
         !! wall, MOM6's free-slip closure).  `.true.` = NO-SLIP:
         !! the factor becomes `2 - wet_q` (image vorticity, MOM6's no-slip).
         !! Set from `&ocean_hvisc_nml no_slip` (shared with the lateral
         !! strain).  All-wet ⇒ `wet_q≡1` ⇒ factor `≡1` ⇒ bit-identical.
      integer :: corner_h_variant = CORNER_H_CELL_MEAN
         !! PV corner-thickness construction (energy scheme).
         !! `CORNER_H_CELL_MEAN` (default, bit-identical) — wet-area 4-cell
         !! mean, `CORIOLIS_H_MIN_PV`-capped.  `CORNER_H_MOM6_AREA` — MOM6's
         !! `q = abs_vort·Area_q/(hArea_q + PV_VOL_NEGLECT)`.  The two are
         !! algebraically identical above the floor; they differ only in the
         !! vanishing-thickness guard.  Set from `&ocean_coriolis_nml corner_h`
         !! at setup (config fail-loud on non-energy forms); read host-side
         !! into a local flag in the kernel.
      logical :: bound_coriolis = .false.
         !! MOM6 `BOUND_CORIOLIS`.  When
         !! `.true.` (energy scheme only — `&ocean_coriolis_nml bound_coriolis`,
         !! config fail-loud on any other form) the energy-scheme PV flux is
         !! clamped into the range of the four neighbouring `(f+ζ)·v`
         !! velocity-form estimates BEFORE the KE-gradient subtraction, capping
         !! the thin-layer `q·vh` blow-up.  Read host-side into a local flag in
         !! the kernel; default `.false.` ⇒ untaken branch ⇒ bit-identical.

      ! ---- Coriolis parameter field ----
      ! `f_corner` holds f at C-grid corners (SW corner of cell
      ! (i, j) sits at position (i-1/2, j-1/2)).  Both kernels read
      ! it here and average to the relevant face.  Shape (nx+1, ny+1).
      !
      ! Initialised to `f_0` everywhere by `init`; switch to a
      ! `f = f_0 + beta * (y - y_ref)` profile via
      ! `coriolis_adv_set_beta_plane`.  For an f-plane just leave
      ! `f_0` at the desired value and skip the beta call.
      real(wp), allocatable :: f_corner(:, :)
         !! Coriolis parameter (1/s) at C-grid corners, shape
         !! (nx+1, ny+1).
      real(wp) :: f_0 = 0.0_wp
         !! f-plane baseline.  `init` fills `f_corner = f_0`; later
         !! re-population (`set_beta_plane`) overrides on a per-
         !! corner basis.
      real(wp) :: beta = 0.0_wp
         !! Meridional gradient `df/dy` (1/(s·m)).  Diagnostic /
         !! convenience storage — the field-of-record is `f_corner`.

      ! ---- Per-step scratch ----
      type(scratch_3d_buffer_t) :: q_corner
         !! Sadourny path: relative vorticity ζ at corners.
         !! HK path: per-mass PV q = (f + ζ) / h_at_corner at corners.
         !! Both kernels write-then-read this in a single call so the
         !! repurposing is safe — they never coexist within one stage.
      type(scratch_3d_buffer_t) :: ke_centre
         !! Kinetic energy at cell centres for the gradient(KE) form.
      type(scratch_3d_buffer_t) :: pv_flux_x
         !! du/dt at east faces.
      type(scratch_3d_buffer_t) :: pv_flux_y
         !! dv/dt at north faces.
      type(scratch_3d_buffer_t) :: mass_flux_u
         !! HK path: u_face_x * h_at_u_face (mass flux per face) at
         !! east faces.  Same shape as `pv_flux_x` — (nx+1, ny, nz).
         !! Unused by the Sadourny path; allocated regardless so the
         !! lifecycle is uniform.
      type(scratch_3d_buffer_t) :: mass_flux_v
         !! HK path: v_face_y * h_at_v_face at north faces.  Shape
         !! (nx, ny+1, nz).  Same usage caveat as `mass_flux_u`.
   contains
      procedure, non_overridable :: init => coriolis_adv_init
      procedure, non_overridable :: destroy => coriolis_adv_destroy
      procedure, non_overridable :: enter_data => coriolis_adv_enter_data
      procedure, non_overridable :: exit_data => coriolis_adv_exit_data
      procedure, non_overridable :: set_beta_plane => coriolis_adv_set_beta_plane
      procedure, non_overridable :: bytes => coriolis_adv_bytes
   end type coriolis_adv_t

contains

   subroutine coriolis_adv_init(this, grid, nz_ml)
      !! Allocate the 4 scratch buffers sized at
      !! (nx_face / ny_face / corner, nz).  Default nz=1 covers the
      !! barotropic kernel; passing `nz_ml` sizes them for the
      !! multilayer kernel.  Same backward-compatible pattern as
      !! `continuity_init`.
      class(coriolis_adv_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      ! Corner-located: (nx+1, ny+1, nz)
      call this%q_corner%init(nx + 1, ny + 1, nz, "coriolis_adv_q_corner")
      ! Cell-centred: (nx, ny, nz)
      call this%ke_centre%init(nx, ny, nz, "coriolis_adv_ke_centre")
      ! East-face: (nx+1, ny, nz) — same shape as u_face_x_layer
      call this%pv_flux_x%init(nx + 1, ny, nz, "coriolis_adv_pv_flux_x")
      ! North-face: (nx, ny+1, nz) — same shape as v_face_y_layer
      call this%pv_flux_y%init(nx, ny + 1, nz, "coriolis_adv_pv_flux_y")
      ! HK mass fluxes — same shapes as pv_flux_x / pv_flux_y.
      call this%mass_flux_u%init(nx + 1, ny, nz, "coriolis_adv_mass_flux_u")
      call this%mass_flux_v%init(nx, ny + 1, nz, "coriolis_adv_mass_flux_v")

      ! Coriolis parameter at corners, initialised to f_0 everywhere.
      allocate (this%f_corner(nx + 1, ny + 1), source=this%f_0)

      this%is_init = .true.
   end subroutine coriolis_adv_init

   subroutine coriolis_adv_destroy(this)
      class(coriolis_adv_t), intent(inout) :: this
      this%is_init = .false.
      call this%q_corner%destroy()
      call this%ke_centre%destroy()
      call this%pv_flux_x%destroy()
      call this%pv_flux_y%destroy()
      call this%mass_flux_u%destroy()
      call this%mass_flux_v%destroy()
      if (allocated(this%f_corner)) deallocate (this%f_corner)
   end subroutine coriolis_adv_destroy

   subroutine coriolis_adv_enter_data(this)
      ! Bare `copyin(this)` removed (stack-descriptor map → AMD cross-slot
      ! overlap; see ocean_surfstress_enter_data).  f_corner + the scratch
      ! buffers attach below; cor-descriptor presence (so DCs reading
      ! cor%f_corner(i,j) don't per-launch memcpy) comes from the root
      ! copyin(state) in ocean_state_enter_data.  A V100 A/B with copyin(this)
      ! gone is bit-identical and faster overall, so the root copy covers it.
      class(coriolis_adv_t), intent(inout) :: this
      select type (this)
      type is (coriolis_adv_t)
         call coriolis_adv_enter_data_impl(this)
      end select
   end subroutine coriolis_adv_enter_data

   subroutine coriolis_adv_enter_data_impl(this)
      type(coriolis_adv_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%q_corner)
      call scratch_3d_buffer_enter_data_impl(this%ke_centre)
      call scratch_3d_buffer_enter_data_impl(this%pv_flux_x)
      call scratch_3d_buffer_enter_data_impl(this%pv_flux_y)
      call scratch_3d_buffer_enter_data_impl(this%mass_flux_u)
      call scratch_3d_buffer_enter_data_impl(this%mass_flux_v)
      !$acc enter data copyin(this%f_corner)
   end subroutine coriolis_adv_enter_data_impl

   subroutine coriolis_adv_exit_data(this)
      class(coriolis_adv_t), intent(inout) :: this
      select type (this)
      type is (coriolis_adv_t)
         call coriolis_adv_exit_data_impl(this)
      end select
   end subroutine coriolis_adv_exit_data

   subroutine coriolis_adv_exit_data_impl(this)
      type(coriolis_adv_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%q_corner)
      call scratch_3d_buffer_exit_data_impl(this%ke_centre)
      call scratch_3d_buffer_exit_data_impl(this%pv_flux_x)
      call scratch_3d_buffer_exit_data_impl(this%pv_flux_y)
      call scratch_3d_buffer_exit_data_impl(this%mass_flux_u)
      call scratch_3d_buffer_exit_data_impl(this%mass_flux_v)
      !$acc exit data delete(this%f_corner)
   end subroutine coriolis_adv_exit_data_impl

   subroutine coriolis_adv_set_beta_plane(this, grid, f_0, beta, y_ref)
      !! Populate `f_corner` with a beta-plane profile
      !!   f(y) = f_0 + beta * (y - y_ref)
      !!
      !! The C-grid SW corner of cell (i, j) sits at physical position
      !! `y = (j - 1 - nghost) * dy` under our convention (j=1+nghost
      !! is the first interior corner row).  Caller picks `y_ref` —
      !! typically the centre of the physical domain.
      !!
      !! Also records `f_0` and `beta` on the type for diagnostics.
      !! Must be called *after* `init` (so `f_corner` is allocated)
      !! and *before* `enter_data` if running on GPU (or follow with
      !! a host→device update if already mapped).
      class(coriolis_adv_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_0, beta, y_ref

      integer :: i, j, nx, ny
      real(wp) :: y

      nx = grid%nx_total
      ny = grid%ny_total
      this%f_0 = f_0
      this%beta = beta

      do j = 1, ny + 1
         y = real(j - 1 - grid%nghost, wp)*grid%dy
         do i = 1, nx + 1
            this%f_corner(i, j) = f_0 + beta*(y - y_ref)
         end do
      end do
   end subroutine coriolis_adv_set_beta_plane

   pure subroutine coriolis_adv_compute_tendencies_barotropic(grid, metrics, this, bs)
      !! Sadourny (1975) energy-conserving Coriolis + horizontal-
      !! momentum-advection form on the barotropic C-grid state:
      !!
      !!   du/dt = +(zeta + f) * v_at_u_face - d/dx(KE)
      !!   dv/dt = -(zeta + f) * u_at_v_face - d/dy(KE)
      !!
      !! where zeta = dv/dx - du/dy is the relative vorticity at
      !! cell corners (south-west corner of cell (i, j) sits at
      !! position (i-1/2, j-1/2)) and KE = (1/2)*(u^2 + v^2) is the
      !! kinetic energy per unit mass, averaged from the surrounding
      !! face velocities at each cell centre.
      !!
      !! For spatially uniform u, v this reduces to plain Coriolis:
      !! zeta vanishes and d/dx(KE) vanishes, leaving the f*v / -f*u
      !! pair (the Phase 3a kernel).
      !!
      !! Three passes:
      !!   1. zeta at corners (q_corner buffer; reused as "vorticity"
      !!      until Phase 5 adds the q = (zeta+f)/h division for
      !!      layered PV).
      !!   2. KE at cell centres (ke_centre buffer), using
      !!      KE = (1/4)*(u_W^2 + u_E^2 + v_S^2 + v_N^2) for the
      !!      energy-consistent form.
      !!   3. Tendencies at faces, combining the vorticity-flux and
      !!      KE-gradient terms.
      !!
      !! Walls: zeta uses zero fallback at outer corners (where the
      !! 4-point stencil falls off the grid); face tendencies at
      !! domain walls are zeroed (closed-wall, no momentum at the
      !! boundary).
      !!
      !! Loop order: j-then-i for NVHPC GPU coalescing.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(coriolis_adv_t), intent(inout) :: this
      type(barotropic_state_t), intent(in) :: bs

      integer :: i, j, nx, ny
      real(wp) :: v_at_u, u_at_v, zeta_at_u, zeta_at_v
      real(wp) :: f_at_u, f_at_v, ke_grad_x, ke_grad_y
      real(wp) :: h_vf_SW, h_vf_NW, h_vf_SE, h_vf_NE
      real(wp) :: h_uf_SW, h_uf_NW, h_uf_SE, h_uf_NE
      real(wp) :: vh_sum, uh_sum, h_eff_sum
      real(wp) :: ns

      nx = grid%nx_total
      ny = grid%ny_total
      ! Slip selector (C1): ns=0 ⇒ factor=wet_q (free-slip); ns=1 ⇒
      ! factor=2-wet_q (no-slip image vorticity).  Branchless per corner.
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)

      ! ---- Pass 1: relative vorticity at SW corners (circulation/area) ----
      ! zeta_corner(i, j) sits at position (i-1/2, j-1/2).  Curvilinear
      ! circulation form (design §2):
      !   zeta = ( v(i,j)·dyCv(i,j) - v(i-1,j)·dyCv(i-1,j)
      !          - (u(i,j)·dxCu(i,j) - u(i,j-1)·dxCu(i,j-1)) ) · iareaBu(i,j)
      ! On uniform square metrics this is `(Δv)/dx - (Δu)/dy` bitwise.
      ! Interior corners: i=2..nx, j=2..ny.  Outer corners
      ! (i=1, j=1, i=nx+1, j=ny+1) get zero — no neighbouring cell
      ! across the wall.  The rel-vort is multiplied by the slip factor
      ! `(1-2·ns)·wet_q + 2·ns` (C1, MOM6 free-slip/no-slip); planetary f added
      ! later stays UNMASKED.  Interior land corners zero exactly as the
      ! domain-wall corners do.
      do concurrent(j=2:ny, i=2:nx)
         this%q_corner%data(i, j, 1) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((bs%v_face_y(i, j)*metrics%dyCv(i, j) - bs%v_face_y(i - 1, j)*metrics%dyCv(i - 1, j)) - &
             (bs%u_face_x(i, j)*metrics%dxCu(i, j) - bs%u_face_x(i, j - 1)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      ! Outer-corner fallback: zero vorticity at the boundary.
      do concurrent(j=1:ny + 1)
         this%q_corner%data(1, j, 1) = 0.0_wp
         this%q_corner%data(nx + 1, j, 1) = 0.0_wp
      end do
      do concurrent(i=1:nx + 1)
         this%q_corner%data(i, 1, 1) = 0.0_wp
         this%q_corner%data(i, ny + 1, 1) = 0.0_wp
      end do

      ! ---- Pass 2: kinetic energy at cell centres (area-weighted) ----
      ! KE_centre(i,j) = 0.25·iareaT·( areaCu(i)·u(i)² + areaCu(i+1)·u(i+1)²
      !                              + areaCv(j)·v(j)² + areaCv(j+1)·v(j+1)² )
      ! On uniform metrics areaCu=areaCv=areaT, iareaT=1/areaT, so this
      ! reduces to the simple 0.25·(u²+u²+v²+v²) form bitwise.
      do concurrent(j=1:ny, i=1:nx)
         this%ke_centre%data(i, j, 1) = 0.25_wp*metrics%iareaT(i, j)*( &
                                        metrics%areaCu(i, j)*bs%u_face_x(i, j)**2 + &
                                        metrics%areaCu(i + 1, j)*bs%u_face_x(i + 1, j)**2 + &
                                        metrics%areaCv(i, j)*bs%v_face_y(i, j)**2 + &
                                        metrics%areaCv(i, j + 1)*bs%v_face_y(i, j + 1)**2)
      end do

      ! ---- Pass 3a: du/dt at interior east faces ----
      ! Thickness-weighted v at the u-face: average `bs%mass_flux_y`
      ! (= v · h_face) over the 4 abutting v-faces, divide by the
      ! sum of the v-face thicknesses.  See the multilayer kernel
      ! comment for the rationale — same form, dropped k axis.
      ! For uniform h this is bit-identical to the simple 4-point
      ! velocity average.
      do concurrent(j=1:ny, i=2:nx) &
         local(f_at_u, h_vf_SW, h_vf_NW, h_vf_SE, h_vf_NE, &
               vh_sum, h_eff_sum)
         h_vf_SW = 0.5_wp*(bs%h(i - 1, max(1, j - 1)) + bs%h(i - 1, j))
         h_vf_NW = 0.5_wp*(bs%h(i - 1, j) + bs%h(i - 1, min(ny, j + 1)))
         h_vf_SE = 0.5_wp*(bs%h(i, max(1, j - 1)) + bs%h(i, j))
         h_vf_NE = 0.5_wp*(bs%h(i, j) + bs%h(i, min(ny, j + 1)))
         vh_sum = (bs%v_face_y(i - 1, j)*h_vf_SW + bs%v_face_y(i - 1, j + 1)*h_vf_NW) + &
                  (bs%v_face_y(i, j)*h_vf_SE + bs%v_face_y(i, j + 1)*h_vf_NE)
         h_eff_sum = (h_vf_SW + h_vf_NW) + (h_vf_SE + h_vf_NE)
         if (h_eff_sum > 0.0_wp) then
            v_at_u = vh_sum/h_eff_sum
         else
            v_at_u = 0.0_wp
         end if
         zeta_at_u = 0.5_wp*(this%q_corner%data(i, j, 1) + &
                             this%q_corner%data(i, j + 1, 1))
         f_at_u = 0.5_wp*(this%f_corner(i, j) + this%f_corner(i, j + 1))
         ke_grad_x = (this%ke_centre%data(i, j, 1) - &
                      this%ke_centre%data(i - 1, j, 1))*metrics%idxCu(i, j)
         this%pv_flux_x%data(i, j, 1) = (zeta_at_u + f_at_u)*v_at_u - ke_grad_x
      end do
      do concurrent(j=1:ny)
         this%pv_flux_x%data(1, j, 1) = 0.0_wp
         this%pv_flux_x%data(nx + 1, j, 1) = 0.0_wp
      end do

      ! ---- Pass 3b: dv/dt at interior north faces ----
      do concurrent(j=2:ny, i=1:nx) &
         local(f_at_v, h_uf_SW, h_uf_NW, h_uf_SE, h_uf_NE, &
               uh_sum, h_eff_sum)
         h_uf_SW = 0.5_wp*(bs%h(max(1, i - 1), j - 1) + bs%h(i, j - 1))
         h_uf_SE = 0.5_wp*(bs%h(i, j - 1) + bs%h(min(nx, i + 1), j - 1))
         h_uf_NW = 0.5_wp*(bs%h(max(1, i - 1), j) + bs%h(i, j))
         h_uf_NE = 0.5_wp*(bs%h(i, j) + bs%h(min(nx, i + 1), j))
         uh_sum = (bs%u_face_x(i, j - 1)*h_uf_SW + bs%u_face_x(i + 1, j - 1)*h_uf_SE) + &
                  (bs%u_face_x(i, j)*h_uf_NW + bs%u_face_x(i + 1, j)*h_uf_NE)
         h_eff_sum = (h_uf_SW + h_uf_SE) + (h_uf_NW + h_uf_NE)
         if (h_eff_sum > 0.0_wp) then
            u_at_v = uh_sum/h_eff_sum
         else
            u_at_v = 0.0_wp
         end if
         zeta_at_v = 0.5_wp*(this%q_corner%data(i, j, 1) + &
                             this%q_corner%data(i + 1, j, 1))
         f_at_v = 0.5_wp*(this%f_corner(i, j) + this%f_corner(i + 1, j))
         ke_grad_y = (this%ke_centre%data(i, j, 1) - &
                      this%ke_centre%data(i, j - 1, 1))*metrics%idyCv(i, j)
         this%pv_flux_y%data(i, j, 1) = -(zeta_at_v + f_at_v)*u_at_v - ke_grad_y
      end do
      do concurrent(i=1:nx)
         this%pv_flux_y%data(i, 1, 1) = 0.0_wp
         this%pv_flux_y%data(i, ny + 1, 1) = 0.0_wp
      end do
   end subroutine coriolis_adv_compute_tendencies_barotropic

   subroutine coriolis_adv_compute_tendencies(grid, metrics, this, ms, u_src, v_src, h_src, &
                                              use_state_fluxes)
      !! Per-layer Coriolis + advection dispatcher.  Reads
      !! `this%pv_variant` and routes to the matching kernel body:
      !!
      !!   PV_VARIANT_SADOURNY_HK     → coriolis_adv_compute_tendencies_hk
      !!     (Arakawa-Hsu 1990 — wider 3-corner stencil, suppresses
      !!      the Hollingsworth-Källén instability)
      !!   PV_VARIANT_SADOURNY_ENERGY → coriolis_adv_compute_tendencies_sadourny_energy
      !!     (Sadourny 1975 energy-conserving transport form, q·vh —
      !!      MOM6 SADOURNY75_ENERGY)
      !!   otherwise                  → coriolis_adv_compute_tendencies_sadourny
      !!     (classical Sadourny 1975 enstrophy-conserving + al81)
      !!
      !! Callers that want to pin to a specific variant can call
      !! `_hk` / `_sadourny_energy` directly; tests do this.
      !!
      !! `use_state_fluxes` (optional, default .false.): consume the
      !! continuity-renormalised `ms%mass_flux_*_layer` transports in
      !! place of the kernel-internal `u·h_face` recompute — the MOM6
      !! mass-consistent CorAdCalc.  Only the `sadourny_energy`
      !! transport form honours it (config fail-loud enforces this).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(coriolis_adv_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      ! assumed-shape-ok: outer-shim passthrough, resolved to explicit-shape
      ! dummies in the variant bodies; never indexed here.
      real(wp), intent(in), optional :: u_src(:, :, :), v_src(:, :, :), h_src(:, :, :)
         !! Optional velocity/thickness source override (SPEC §4 S3,
         !! `split_scheme = "pred_corr"`): the driver passes the step
         !! time-means `ms%u_av_layer` / `v_av_layer` / `h_av_layer` so the
         !! Coriolis-advection tendency is evaluated on the MOM6 `u_av`
         !! family, never the prognostic (MOM6's `CorAdCalc(u_av, v_av,
         !! h_av, ...)`).  All three must be
         !! passed together.  Absent ⇒ prognostic arrays, bit-identical.
      logical, intent(in), optional :: use_state_fluxes

      logical :: usf

      usf = .false.
      if (present(use_state_fluxes)) usf = use_state_fluxes

      if (present(u_src)) then
         if (this%pv_variant == PV_VARIANT_SADOURNY_HK) then
            call coriolis_adv_compute_tendencies_hk(grid, metrics, this, ms, &
                                                    u_src, v_src, h_src)
         else if (this%pv_variant == PV_VARIANT_SADOURNY_ENERGY) then
            call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, this, ms, &
                                                                 u_src, v_src, h_src, &
                                                                 use_state_fluxes=usf)
         else
            call coriolis_adv_compute_tendencies_sadourny(grid, metrics, this, ms, &
                                                          u_src, v_src, h_src)
         end if
      else if (this%pv_variant == PV_VARIANT_SADOURNY_HK) then
         call coriolis_adv_compute_tendencies_hk(grid, metrics, this, ms, &
                                                 ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      else if (this%pv_variant == PV_VARIANT_SADOURNY_ENERGY) then
         call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, this, ms, &
                                                              ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, &
                                                              use_state_fluxes=usf)
      else
         call coriolis_adv_compute_tendencies_sadourny(grid, metrics, this, ms, &
                                                       ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      end if
   end subroutine coriolis_adv_compute_tendencies

   pure subroutine coriolis_adv_compute_tendencies_sadourny(grid, metrics, this, ms, u, v, h)
      !! Per-layer Sadourny Coriolis + advection tendency.  Same
      !! algorithm as the barotropic counterpart, lifted with a
      !! k-axis on every loop.  Each k-slice is independent (ζ
      !! stencil only reads same-k velocities; KE at centre only
      !! reads same-k face values), so the do-concurrent kernels
      !! parallelise over (k, j, i) for full GPU occupancy.
      !!
      !! The Coriolis parameter is read from `this%f_corner` (2D
      !! field at C-grid corners, shared with the barotropic
      !! kernel).  Uniform `f_0` is the f-plane default; call
      !! `this%set_beta_plane(grid, f_0, beta, y_ref)` to switch to
      !! a `f = f_0 + beta*(y - y_ref)` profile.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(coriolis_adv_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: u(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
         !! Face-velocity / thickness source arrays (outer-shim; the
         !! dispatcher forwards either the prognostic components or the
         !! `u_av` time-mean family under `split_scheme = "pred_corr"`).
      real(wp), intent(in) :: v(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total, ms%nz_ml)

      integer :: i, j, k, nx, ny, nz, pv_scheme
      real(wp) :: v_at_u, u_at_v, zeta_at_u, zeta_at_v
      real(wp) :: f_at_u, f_at_v, ke_grad_x, ke_grad_y
      real(wp) :: h_vf_SW, h_vf_NW, h_vf_SE, h_vf_NE
      real(wp) :: h_uf_SW, h_uf_NW, h_uf_SE, h_uf_NE
      real(wp) :: vh_sum, uh_sum, h_eff_sum
      real(wp) :: ns

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      ! Hoist the loop-invariant PV face-interp scheme to a plain scalar so
      ! the do-concurrent kernels never touch a derived-type component.
      pv_scheme = this%pv_adv_scheme
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)

      ! ---- Pass 1: relative vorticity at SW corners, per layer ----
      ! Circulation/area (design §2) — see the barotropic kernel for the
      ! reduction to `(Δv)/dx - (Δu)/dy` on uniform square metrics.
      ! Slip factor `(1-2·ns)·wet_q + 2·ns` masks the rel-vort at land
      ! corners (C1, free-slip default); planetary f stays unmasked.
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         this%q_corner%data(i, j, k) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((v(i, j, k)*metrics%dyCv(i, j) - &
              v(i - 1, j, k)*metrics%dyCv(i - 1, j)) - &
             (u(i, j, k)*metrics%dxCu(i, j) - &
              u(i, j - 1, k)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         this%q_corner%data(1, j, k) = 0.0_wp
         this%q_corner%data(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         this%q_corner%data(i, 1, k) = 0.0_wp
         this%q_corner%data(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 2: KE at cell centres, per layer (area-weighted) ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         this%ke_centre%data(i, j, k) = 0.25_wp*metrics%iareaT(i, j)*( &
                                        metrics%areaCu(i, j)*u(i, j, k)**2 + &
                                        metrics%areaCu(i + 1, j)*u(i + 1, j, k)**2 + &
                                        metrics%areaCv(i, j)*v(i, j, k)**2 + &
                                        metrics%areaCv(i, j + 1)*v(i, j + 1, k)**2)
      end do

      ! ---- Pass 3a: du/dt at interior east faces ----
      ! Thickness-weighted v at the u-face:
      !   v_at_u = Σ(v_face_y · h_at_v_face) / Σ(h_at_v_face)
      ! over the four abutting v-faces.  For uniform `h_layer` this
      ! reduces to the simple 4-point velocity average and is
      ! bit-identical to the previous code.  Under variable
      ! thickness it captures the mass-weighted advection the
      ! centred PV-advection form requires.  We compute `v · h_face`
      ! inline rather than reading `mass_flux_y_layer` so the kernel
      ! is self-contained — it doesn't matter whether continuity has
      ! run yet in this stage.
      ! ---- Pass 3a: u-face Coriolis-advection term (ζ+f)·v_at_u ----
      ! Writes (ζ+f)·v_at_u into pv_flux_x; the −∇KE term is subtracted
      ! in Pass 3c BELOW.  The split keeps the −∇KE subtraction in one
      ! place and future-proofs a BOUND_CORIOLIS clip on `CAu` BEFORE
      ! `−KEx` (mirroring MOM6's MOM_CoriolisAdv ordering — not yet ported).
      !
      ! This is the classical Sadourny (1975) ENSTROPHY-conserving form
      ! `CAu = (zeta_at_u + f_at_u)·v_at_u` (the default `form="sadourny"`,
      ! and `al81`).  The energy-conserving transport form
      ! (`form="sadourny_energy"`, MOM6 SADOURNY75_ENERGY) is a separate
      ! kernel, `coriolis_adv_compute_tendencies_sadourny_energy`.
      associate (q_corner => this%q_corner%data, f_corner => this%f_corner)
         do concurrent(k=1:nz, j=1:ny, i=2:nx) &
            local(f_at_u, h_vf_SW, h_vf_NW, h_vf_SE, h_vf_NE, &
                  vh_sum, h_eff_sum)
            h_vf_SW = 0.5_wp*(h(i - 1, max(1, j - 1), k) + h(i - 1, j, k))
            h_vf_NW = 0.5_wp*(h(i - 1, j, k) + h(i - 1, min(ny, j + 1), k))
            h_vf_SE = 0.5_wp*(h(i, max(1, j - 1), k) + h(i, j, k))
            h_vf_NE = 0.5_wp*(h(i, j, k) + h(i, min(ny, j + 1), k))
            vh_sum = (v(i - 1, j, k)*h_vf_SW + &
                      v(i - 1, j + 1, k)*h_vf_NW) + &
                     (v(i, j, k)*h_vf_SE + &
                      v(i, j + 1, k)*h_vf_NE)
            h_eff_sum = (h_vf_SW + h_vf_NW) + (h_vf_SE + h_vf_NE)
            if (h_eff_sum > 0.0_wp) then
               v_at_u = vh_sum/h_eff_sum
            else
               v_at_u = 0.0_wp
            end if
            ! Absolute vorticity (f+zeta) interpolated onto the u-face along j,
            ! upwind on v_at_u.  WENO reconstructs it directly (f baked into the
            ! stencil, MOM6 reconstructs f+zeta); centred = the 2-point average.
            ! Each order falls back to centred within its stencil radius of the
            ! j=1 / j=ny array edges (the nghost gate keeps every PHYSICAL face
            ! inside the band, so only ghost faces degrade).
            if (pv_scheme == PV_ADV_WENO7 .and. j >= 4 .and. j <= ny - 3) then
               zeta_at_u = weno7_recon( &
                           q_corner(i, j - 3, k) + f_corner(i, j - 3), &
                           q_corner(i, j - 2, k) + f_corner(i, j - 2), &
                           q_corner(i, j - 1, k) + f_corner(i, j - 1), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i, j + 1, k) + f_corner(i, j + 1), &
                           q_corner(i, j + 2, k) + f_corner(i, j + 2), &
                           q_corner(i, j + 3, k) + f_corner(i, j + 3), &
                           q_corner(i, j + 4, k) + f_corner(i, j + 4), v_at_u)
            else if (pv_scheme == PV_ADV_WENO5 .and. j >= 3 .and. j <= ny - 2) then
               zeta_at_u = weno5_recon( &
                           q_corner(i, j - 2, k) + f_corner(i, j - 2), &
                           q_corner(i, j - 1, k) + f_corner(i, j - 1), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i, j + 1, k) + f_corner(i, j + 1), &
                           q_corner(i, j + 2, k) + f_corner(i, j + 2), &
                           q_corner(i, j + 3, k) + f_corner(i, j + 3), v_at_u)
            else if (pv_scheme == PV_ADV_WENO3 .and. j >= 2 .and. j <= ny - 1) then
               zeta_at_u = weno3_recon( &
                           q_corner(i, j - 1, k) + f_corner(i, j - 1), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i, j + 1, k) + f_corner(i, j + 1), &
                           q_corner(i, j + 2, k) + f_corner(i, j + 2), v_at_u)
            else
               zeta_at_u = 0.5_wp*(q_corner(i, j, k) + q_corner(i, j + 1, k)) + &
                           0.5_wp*(f_corner(i, j) + f_corner(i, j + 1))
            end if
            this%pv_flux_x%data(i, j, k) = zeta_at_u*v_at_u
         end do
         ! ---- Pass 3c: subtract −∇KE from u-tendency ----
         do concurrent(k=1:nz, j=1:ny, i=2:nx) local(ke_grad_x)
            ke_grad_x = (this%ke_centre%data(i, j, k) - &
                         this%ke_centre%data(i - 1, j, k))*metrics%idxCu(i, j)
            this%pv_flux_x%data(i, j, k) = this%pv_flux_x%data(i, j, k) - ke_grad_x
         end do
         do concurrent(k=1:nz, j=1:ny)
            this%pv_flux_x%data(1, j, k) = 0.0_wp
            this%pv_flux_x%data(nx + 1, j, k) = 0.0_wp
         end do

         ! ---- Pass 4a: v-face Coriolis-advection term −(ζ+f)·u_at_v ----
         ! Mirror of Pass 3a.  Writes −(ζ+f)·u_at_v ONLY; ∇KE handled in 4c.
         do concurrent(k=1:nz, j=2:ny, i=1:nx) &
            local(f_at_v, h_uf_SW, h_uf_NW, h_uf_SE, h_uf_NE, &
                  uh_sum, h_eff_sum)
            h_uf_SW = 0.5_wp*(h(max(1, i - 1), j - 1, k) + h(i, j - 1, k))
            h_uf_SE = 0.5_wp*(h(i, j - 1, k) + h(min(nx, i + 1), j - 1, k))
            h_uf_NW = 0.5_wp*(h(max(1, i - 1), j, k) + h(i, j, k))
            h_uf_NE = 0.5_wp*(h(i, j, k) + h(min(nx, i + 1), j, k))
            uh_sum = (u(i, j - 1, k)*h_uf_SW + &
                      u(i + 1, j - 1, k)*h_uf_SE) + &
                     (u(i, j, k)*h_uf_NW + &
                      u(i + 1, j, k)*h_uf_NE)
            h_eff_sum = (h_uf_SW + h_uf_SE) + (h_uf_NW + h_uf_NE)
            if (h_eff_sum > 0.0_wp) then
               u_at_v = uh_sum/h_eff_sum
            else
               u_at_v = 0.0_wp
            end if
            ! Absolute vorticity onto the v-face along i, upwind on u_at_v.
            ! Sign mirrors the centred form (CAv = -(f+zeta)*u_at_v).  Same
            ! per-order boundary fallback as the u-face.
            if (pv_scheme == PV_ADV_WENO7 .and. i >= 4 .and. i <= nx - 3) then
               zeta_at_v = weno7_recon( &
                           q_corner(i - 3, j, k) + f_corner(i - 3, j), &
                           q_corner(i - 2, j, k) + f_corner(i - 2, j), &
                           q_corner(i - 1, j, k) + f_corner(i - 1, j), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i + 1, j, k) + f_corner(i + 1, j), &
                           q_corner(i + 2, j, k) + f_corner(i + 2, j), &
                           q_corner(i + 3, j, k) + f_corner(i + 3, j), &
                           q_corner(i + 4, j, k) + f_corner(i + 4, j), u_at_v)
            else if (pv_scheme == PV_ADV_WENO5 .and. i >= 3 .and. i <= nx - 2) then
               zeta_at_v = weno5_recon( &
                           q_corner(i - 2, j, k) + f_corner(i - 2, j), &
                           q_corner(i - 1, j, k) + f_corner(i - 1, j), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i + 1, j, k) + f_corner(i + 1, j), &
                           q_corner(i + 2, j, k) + f_corner(i + 2, j), &
                           q_corner(i + 3, j, k) + f_corner(i + 3, j), u_at_v)
            else if (pv_scheme == PV_ADV_WENO3 .and. i >= 2 .and. i <= nx - 1) then
               zeta_at_v = weno3_recon( &
                           q_corner(i - 1, j, k) + f_corner(i - 1, j), &
                           q_corner(i, j, k) + f_corner(i, j), &
                           q_corner(i + 1, j, k) + f_corner(i + 1, j), &
                           q_corner(i + 2, j, k) + f_corner(i + 2, j), u_at_v)
            else
               zeta_at_v = 0.5_wp*(q_corner(i, j, k) + q_corner(i + 1, j, k)) + &
                           0.5_wp*(f_corner(i, j) + f_corner(i + 1, j))
            end if
            this%pv_flux_y%data(i, j, k) = -zeta_at_v*u_at_v
         end do
      end associate

      ! ---- Pass 4c: subtract −∇KE from v-tendency ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(ke_grad_y)
         ke_grad_y = (this%ke_centre%data(i, j, k) - &
                      this%ke_centre%data(i, j - 1, k))*metrics%idyCv(i, j)
         this%pv_flux_y%data(i, j, k) = this%pv_flux_y%data(i, j, k) - ke_grad_y
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%pv_flux_y%data(i, 1, k) = 0.0_wp
         this%pv_flux_y%data(i, ny + 1, k) = 0.0_wp
      end do
   end subroutine coriolis_adv_compute_tendencies_sadourny

   pure subroutine coriolis_adv_compute_tendencies_hk(grid, metrics, this, ms, u, v, h)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Per-layer PV-conserving Coriolis + horizontal-advection tendency
      !! in the Arakawa-Hsu (1990) form ("HK correction").  The wider
      !! 3-corner PV stencil at each face suppresses the spurious
      !! Hollingsworth-Källén instability that biases the simpler
      !! Sadourny 2-corner form at eddy-resolving resolutions.
      !!
      !! Algorithm — 6 passes per call:
      !!   1. Pass 1: relative vorticity ζ at interior corners (same
      !!      stencil as the Sadourny multilayer kernel) into `q_corner`.
      !!      Wall corners get ζ = 0 (free-slip BC).
      !!   2. Pass 2: per-mass PV `q = (f + ζ) / h_at_corner` rewritten
      !!      into `q_corner`.  `h_at_corner` is the AREA-WEIGHTED 4-cell
      !!      mean (each T thickness weighted by its `areaT`, normalised
      !!      by the summed areas), with `min/max` clamps so wall corners
      !!      collapse onto the available cells.  Reduces to the plain
      !!      4-cell mean on uniform Cartesian (equal areas).
      !!   3. Pass 3: per-face mass fluxes — `mass_flux_u = u_face_x ·
      !!      h_at_u_face` (face-averaged thickness) and `mass_flux_v`
      !!      symmetrically.  Wall faces fall through to the single
      !!      available cell (velocities are zero there anyway).
      !!   4. Pass 4: KE at cell centres — identical to Sadourny.
      !!   5. Pass 5: u-face tendency
      !!        CAu(i,j) = a · mass_flux_v(i,   j+1)   (NE)
      !!                 + b · mass_flux_v(i-1, j+1)   (NW)
      !!                 + c · mass_flux_v(i-1, j)     (SW)
      !!                 + d · mass_flux_v(i,   j)     (SE)
      !!                 - grad_KE_x
      !!      where each coefficient combines the face's two end
      !!      corners + one diagonal corner:
      !!        a = (q(i,j+1) + q(i+1,j+1) + q(i,j))   / 12
      !!        b = (q(i,j+1) + q(i-1,j+1) + q(i,j))   / 12
      !!        c = (q(i,j+1) + q(i-1,j)   + q(i,j))   / 12
      !!        d = (q(i,j+1) + q(i+1,j)   + q(i,j))   / 12
      !!   6. Pass 6: v-face tendency — symmetric construction; overall
      !!      minus sign on the q-stencil sum (Coriolis on v is `-f·u`):
      !!        CAv(i,j) = -[a' · mass_flux_u(i+1, j)
      !!                   + b' · mass_flux_u(i,   j)
      !!                   + c' · mass_flux_u(i,   j-1)
      !!                   + d' · mass_flux_u(i+1, j-1)]
      !!                   - grad_KE_y
      !!
      !! Reduction property (uniform h, uniform v): each coefficient
      !! evaluates to `q/4`, so the sum of 4 mass-flux terms is `q · vh`
      !! and the kernel collapses to the Sadourny `(f+ζ)·v` form
      !! bit-identically.  This is the basis for the regression test.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(coriolis_adv_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: u(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
         !! Face-velocity / thickness source arrays (outer-shim; the
         !! dispatcher forwards either the prognostic components or the
         !! `u_av` time-mean family under `split_scheme = "pred_corr"`).
      real(wp), intent(in) :: v(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total, ms%nz_ml)

      integer :: i, j, k, nx, ny, nz, nu, nv
      real(wp) :: zeta_corner, h_corner, h_face
      real(wp) :: aSW, aSE, aNW, aNE, hm_num, hm_den
      integer :: iw, ie, js, jn
      real(wp) :: q_S, q_N, q_W, q_E, q_NE, q_NW, q_SE, q_SW
      real(wp) :: a_NE, b_NW, c_SW, d_SE
      real(wp) :: ke_grad_x, ke_grad_y
      real(wp) :: ns
      real(wp), parameter :: C1_12 = 1.0_wp/12.0_wp
      real(wp), parameter :: H_MIN_PV = 1.0e-12_wp
         !! Floor for the corner-h divide; vanishing-layer columns
         !! get q ≈ (f+ζ)/H_MIN_PV which is large but finite — paired
         !! with `mass_flux ≈ 0` at the same column so the product
         !! decays cleanly toward zero rather than blowing up.

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nu = size(u, 1)
      nv = size(v, 2)
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)

      ! ---- Pass 1: relative vorticity at corners (circulation/area) ----
      ! Slip factor masks the rel-vort at land corners (C1); Pass 2 reads
      ! this back as `zeta_corner` and adds the UNMASKED planetary f.
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         this%q_corner%data(i, j, k) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((v(i, j, k)*metrics%dyCv(i, j) - &
              v(i - 1, j, k)*metrics%dyCv(i - 1, j)) - &
             (u(i, j, k)*metrics%dxCu(i, j) - &
              u(i, j - 1, k)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         this%q_corner%data(1, j, k) = 0.0_wp
         this%q_corner%data(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         this%q_corner%data(i, 1, k) = 0.0_wp
         this%q_corner%data(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 2: PV q = (f + ζ) / h_at_corner ----
      ! h_at_corner is the AREA-WEIGHTED mean of the four surrounding
      ! T cells, with min/max wall clamps — wall corners collapse onto
      ! the available cells (single cell at the four grid corners;
      ! two-cell mean along an edge).  Area weighting (each T thickness
      ! weighted by its own areaT, normalised by the summed areas) is the
      ! conservative corner thickness on curvilinear grids where
      ! areaT varies cell-to-cell; on uniform Cartesian every areaT is
      ! equal so it reduces to the plain 4-cell mean (same value).
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx + 1) &
         local(zeta_corner, h_corner, iw, ie, js, jn, &
               aSW, aSE, aNW, aNE, hm_num, hm_den)
         zeta_corner = this%q_corner%data(i, j, k)
         iw = max(1, i - 1)
         ie = min(nx, i)
         js = max(1, j - 1)
         jn = min(ny, j)
         ! Mask-weighted corner area (C1, MOM6 `Area_h = mask2dT·areaT`):
         ! a blocked T-column contributes zero area + zero thickness, so
         ! `h_corner` is the wet-column mean only.  All-wet ⇒ `wet_T≡1` ⇒
         ! plain area-weighted mean (bit-identical).  `hm_den` floored at
         ! H_DIV_EPS for a fully-land corner (zeta=0 there, mass fluxes 0).
         !
         ! This is the MOM6 SADOURNY75_ENERGY (transport/energy-form) land
         ! treatment, NOT a divergence from it.  MOM6 masks corner area
         ! UNCONDITIONALLY at land (`Area_h = mask2dT·areaT`, zero on a land
         ! T-column) then forms `q = abs_vort·Area_q/(hArea_q + vol_neglect)`
         ! with `Area_q = Σ Area_h` over the four corners — algebraically
         ! `abs_vort·Σ(wet·areaT) / Σ(wet·areaT·h)`, IDENTICAL to our
         ! `(f+ζ)/h_corner` with `h_corner = Σ(wet·areaT·h)/Σ(wet·areaT)`.
         ! MOM6's *additional* `Area_h` area-mirroring across a boundary is
         ! OBC-SEGMENT-ONLY (inside `if (associated(OBC))`), not a land-coast
         ! rule — it does not apply here.  No HK-specific corner liberty.
         aSW = metrics%wet_T(iw, js)*metrics%areaT(iw, js)
         aSE = metrics%wet_T(ie, js)*metrics%areaT(ie, js)
         aNW = metrics%wet_T(iw, jn)*metrics%areaT(iw, jn)
         aNE = metrics%wet_T(ie, jn)*metrics%areaT(ie, jn)
         hm_num = aSW*h(iw, js, k) + aSE*h(ie, js, k) + &
                  aNW*h(iw, jn, k) + aNE*h(ie, jn, k)
         hm_den = aSW + aSE + aNW + aNE
         h_corner = hm_num/max(hm_den, H_DIV_EPS)
         h_corner = max(h_corner, H_MIN_PV)
         this%q_corner%data(i, j, k) = (this%f_corner(i, j) + zeta_corner)/h_corner
      end do

      ! ---- Pass 3a: u-face transport uh = u·h·dy_cu ----
      do concurrent(k=1:nz, j=1:ny, i=1:nu) local(h_face)
         if (i == 1) then
            h_face = h(1, j, k)
         else if (i == nu) then
            h_face = h(nx, j, k)
         else
            h_face = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
         end if
         this%mass_flux_u%data(i, j, k) = u(i, j, k)*h_face*metrics%dy_cu(i, j)
      end do

      ! ---- Pass 3b: v-face transport vh = v·h·dx_cv ----
      do concurrent(k=1:nz, j=1:nv, i=1:nx) local(h_face)
         if (j == 1) then
            h_face = h(i, 1, k)
         else if (j == nv) then
            h_face = h(i, ny, k)
         else
            h_face = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
         end if
         this%mass_flux_v%data(i, j, k) = v(i, j, k)*h_face*metrics%dx_cv(i, j)
      end do

      ! ---- Porous barriers (Adcroft 2013) ----
      ! The PV/advection transports must carry the SAME narrowed face
      ! width continuity uses, or the two mass-flux definitions disagree.
      ! Host-side gate => no kernel launch and no textual change to the
      ! passes above when the knob is off (byte-identical).
      if (metrics%use_porous) then
         call porous_narrow_3d(nu, ny, nz, metrics%por_face_area_u, &
                               this%mass_flux_u%data)
         call porous_narrow_3d(nx, nv, nz, metrics%por_face_area_v, &
                               this%mass_flux_v%data)
      end if

      ! ---- Pass 4: KE at cell centres (area-weighted) ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         this%ke_centre%data(i, j, k) = 0.25_wp*metrics%iareaT(i, j)*( &
                                        metrics%areaCu(i, j)*u(i, j, k)**2 + &
                                        metrics%areaCu(i + 1, j)*u(i + 1, j, k)**2 + &
                                        metrics%areaCv(i, j)*v(i, j, k)**2 + &
                                        metrics%areaCv(i, j + 1)*v(i, j + 1, k)**2)
      end do

      ! ---- Pass 5: u-face HK tendency ----
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(q_S, q_N, q_NE, q_NW, q_SE, q_SW, &
               a_NE, b_NW, c_SW, d_SE, ke_grad_x)
         q_S = this%q_corner%data(i, j, k)
         q_N = this%q_corner%data(i, j + 1, k)
         q_NE = this%q_corner%data(i + 1, j + 1, k)
         q_NW = this%q_corner%data(i - 1, j + 1, k)
         q_SE = this%q_corner%data(i + 1, j, k)
         q_SW = this%q_corner%data(i - 1, j, k)
         a_NE = (q_N + q_NE + q_S)*C1_12
         b_NW = (q_N + q_NW + q_S)*C1_12
         c_SW = (q_N + q_SW + q_S)*C1_12
         d_SE = (q_N + q_SE + q_S)*C1_12
         ke_grad_x = (this%ke_centre%data(i, j, k) - &
                      this%ke_centre%data(i - 1, j, k))*metrics%idxCu(i, j)
         ! q·vh sum is a transport-weighted PV flux (m³/s); the u-face
         ! IdxCu closes it to a per-length acceleration (= /dx on uniform).
         this%pv_flux_x%data(i, j, k) = &
            (a_NE*this%mass_flux_v%data(i, j + 1, k) + &
             b_NW*this%mass_flux_v%data(i - 1, j + 1, k) + &
             c_SW*this%mass_flux_v%data(i - 1, j, k) + &
             d_SE*this%mass_flux_v%data(i, j, k))*metrics%idxCu(i, j) - ke_grad_x
      end do
      do concurrent(k=1:nz, j=1:ny)
         this%pv_flux_x%data(1, j, k) = 0.0_wp
         this%pv_flux_x%data(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 6: v-face HK tendency ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(q_W, q_E, q_NE, q_NW, q_SE, q_SW, &
               a_NE, b_NW, c_SW, d_SE, ke_grad_y)
         q_W = this%q_corner%data(i, j, k)
         q_E = this%q_corner%data(i + 1, j, k)
         q_NE = this%q_corner%data(i + 1, j + 1, k)
         q_NW = this%q_corner%data(i, j + 1, k)
         q_SE = this%q_corner%data(i + 1, j - 1, k)
         q_SW = this%q_corner%data(i, j - 1, k)
         a_NE = (q_W + q_NE + q_E)*C1_12
         b_NW = (q_W + q_NW + q_E)*C1_12
         c_SW = (q_W + q_SW + q_E)*C1_12
         d_SE = (q_W + q_SE + q_E)*C1_12
         ke_grad_y = (this%ke_centre%data(i, j, k) - &
                      this%ke_centre%data(i, j - 1, k))*metrics%idyCv(i, j)
         this%pv_flux_y%data(i, j, k) = &
            -(a_NE*this%mass_flux_u%data(i + 1, j, k) + &
              b_NW*this%mass_flux_u%data(i, j, k) + &
              c_SW*this%mass_flux_u%data(i, j - 1, k) + &
              d_SE*this%mass_flux_u%data(i + 1, j - 1, k))*metrics%idyCv(i, j) - ke_grad_y
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%pv_flux_y%data(i, 1, k) = 0.0_wp
         this%pv_flux_y%data(i, ny + 1, k) = 0.0_wp
      end do
   end subroutine coriolis_adv_compute_tendencies_hk

   pure subroutine coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, this, ms, &
                                                                   u, v, h, use_state_fluxes)
      !! Faithful MOM6 SADOURNY75_ENERGY (Sadourny 1975 energy-conserving)
      !! per-layer Coriolis + horizontal-advection tendency.  This is the
      !! TRANSPORT form: the absolute-vorticity flux is the potential
      !! vorticity `q = (f + ζ)/h_at_corner` times the layer MASS TRANSPORT
      !! (`vh`/`uh`), so the discrete Coriolis term produces zero net domain
      !! kinetic energy (energy-conserving).  The default enstrophy form
      !! (`_sadourny`, `(f+ζ)·v`) only matches this under uniform thickness.
      !!
      !!   CAu(i,j) = 0.25·( q_N·(vh_NW+vh_NE) + q_S·(vh_SW+vh_SE) )·idxCu − ∂x KE
      !!   CAv(i,j) = −0.25·( q_E·(uh_NE+uh_SE) + q_W·(uh_NW+uh_SW) )·idyCv − ∂y KE
      !!
      !! with q_N/q_S the north/south corner PVs of the u-face (q_E/q_W the
      !! east/west corner PVs of the v-face) and vh/uh the per-face mass
      !! transports.  Reduction: uniform thickness ⇒ each q·(Σvh) collapses
      !! to `(f+ζ)·v`, bit-identical to the enstrophy form (regression test).
      !!
      !! Passes 1-4 (ζ at corners → PV `q = (f+ζ)/h_corner` with the
      !! wet-area-weighted corner thickness → mass transports → centre KE)
      !! mirror `coriolis_adv_compute_tendencies_hk` exactly; only the final
      !! corner→face stencil differs (Sadourny 2-corner vs HK 12-point).
      !! KEEP THE PREP PASSES IN SYNC with `_hk` (shared-prep extraction is
      !! tracked as a follow-up cleanup).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(coriolis_adv_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: u(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
         !! Face-velocity / thickness source arrays (outer-shim; the
         !! dispatcher forwards either the prognostic components or the
         !! `u_av` time-mean family under `split_scheme = "pred_corr"`).
      real(wp), intent(in) :: v(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total, ms%nz_ml)
      logical, intent(in), optional :: use_state_fluxes
         !! Mass-consistent CorAdCalc (MOM6 parity): fill the transport
         !! buffers from `ms%mass_flux_*_layer` — the continuity solve's
         !! renormalised uh/vh (same `u·h_face·dy_cu` m³/s convention,
         !! same shape, physical walls already zeroed) — instead of
         !! recomputing from the u/h source arrays.  In the pred_corr
         !! CORRECTOR those are the PREDICTOR chain's fluxes, i.e. exactly
         !! the transport field that produced the `u_av` evaluation state,
         !! so the q·vh product is energy-consistent on rim columns where
         !! the renorm/wall-zero and the naive `u·h_face` recompute
         !! disagree.  Absent / `.false.` ⇒ bit-identical recompute path.

      integer :: i, j, k, nx, ny, nz, nu, nv
      real(wp) :: zeta_corner, h_corner, h_face
      real(wp) :: aSW, aSE, aNW, aNE, hm_num, hm_den
      integer :: iw, ie, js, jn
      real(wp) :: q_S, q_N, q_W, q_E, ke_grad_x, ke_grad_y
      real(wp) :: ns
      logical :: usf
      real(wp) :: pv_part, av_a, av_b, fv1, fv2, fv3, fv4
      logical :: do_bound, use_mom6_ch
      real(wp), parameter :: H_MIN_PV = CORIOLIS_H_MIN_PV
         !! Floor for the corner-h divide; vanishing-layer columns get
         !! q ≈ (f+ζ)/H_MIN_PV (large but finite) paired with mass_flux ≈ 0
         !! at the same column so the product decays toward zero.  (MOM6
         !! instead floors the denominator, Area_q/(hArea_q+vol_neglect),
         !! so its q→0 as h→0; the forms differ only in the vanishing-layer
         !! limit, immaterial for non-vanishing envelopes.)

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nu = size(u, 1)
      nv = size(v, 2)
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)
      usf = .false.
      if (present(use_state_fluxes)) usf = use_state_fluxes
      ! BOUND_CORIOLIS host flag (read once; the clamp branch is untaken when
      ! off ⇒ bit-identical).  Energy scheme only (config-guaranteed).
      do_bound = this%bound_coriolis
      ! corner_h host flag: MOM6 area-weighted PV corner thickness.  Read once;
      ! the mom6 branch is untaken (and the div-then-cap path bit-identical to
      ! pre-knob) when off.  Energy scheme only (config-guaranteed).
      use_mom6_ch = this%corner_h_variant == CORNER_H_MOM6_AREA

      ! ---- Pass 1: relative vorticity at corners (circulation/area) ----
      ! Slip factor masks the rel-vort at land corners (C1); Pass 2 reads
      ! this back and adds the UNMASKED planetary f.
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         this%q_corner%data(i, j, k) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((v(i, j, k)*metrics%dyCv(i, j) - &
              v(i - 1, j, k)*metrics%dyCv(i - 1, j)) - &
             (u(i, j, k)*metrics%dxCu(i, j) - &
              u(i, j - 1, k)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         this%q_corner%data(1, j, k) = 0.0_wp
         this%q_corner%data(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         this%q_corner%data(i, 1, k) = 0.0_wp
         this%q_corner%data(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 2: PV q = (f + ζ) / h_at_corner ----
      ! h_at_corner is the wet-area-weighted 4-cell mean (MOM6 `Area_h =
      ! mask2dT·areaT`; land cells contribute zero area + thickness).  q =
      ! abs_vort·Area_q/hArea_q = abs_vort/h_corner.  All-wet ⇒ plain mean.
      ! `hm_num` ≡ MOM6 `hArea_q` (Σ area·h), `hm_den` ≡ MOM6 `Area_q` (Σ area)
      ! — the two constructions share these EXACTLY; they differ only in the
      ! vanishing-thickness guard (cell_mean caps h_corner at H_MIN_PV;
      ! mom6_area's PV_VOL_NEGLECT is pure 1/0 armor, matching MOM6).
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx + 1) &
         local(zeta_corner, h_corner, iw, ie, js, jn, &
               aSW, aSE, aNW, aNE, hm_num, hm_den)
         zeta_corner = this%q_corner%data(i, j, k)
         iw = max(1, i - 1)
         ie = min(nx, i)
         js = max(1, j - 1)
         jn = min(ny, j)
         aSW = metrics%wet_T(iw, js)*metrics%areaT(iw, js)
         aSE = metrics%wet_T(ie, js)*metrics%areaT(ie, js)
         aNW = metrics%wet_T(iw, jn)*metrics%areaT(iw, jn)
         aNE = metrics%wet_T(ie, jn)*metrics%areaT(ie, jn)
         hm_num = aSW*h(iw, js, k) + aSE*h(ie, js, k) + &
                  aNW*h(iw, jn, k) + aNE*h(ie, jn, k)
         hm_den = aSW + aSE + aNW + aNE
         if (use_mom6_ch) then
            ! MOM6 area form: q = abs_vort·Area_q/(hArea_q + vol_neglect).
            ! No thickness cap — vol_neglect is pure 1/0 armor.
            this%q_corner%data(i, j, k) = (this%f_corner(i, j) + zeta_corner)* &
                                          hm_den/(hm_num + PV_VOL_NEGLECT)
         else
            h_corner = hm_num/max(hm_den, H_DIV_EPS)
            h_corner = max(h_corner, H_MIN_PV)
            this%q_corner%data(i, j, k) = (this%f_corner(i, j) + zeta_corner)/h_corner
         end if
      end do

      ! ---- Pass 3a/3b: face transports uh / vh ----
      ! Two sources, same convention (u·h_face·dy_cu, m³/s):
      !   recompute (default) — self-contained `u·h_face` from the u/h
      !     source arrays; order-independent of continuity.
      !   state fluxes (`use_state_fluxes`) — copy the continuity
      !     solve's renormalised `ms%mass_flux_*_layer` (MOM6
      !     mass-consistent CorAdCalc; pred_corr-corrector stage only,
      !     where they still hold the predictor chain's fluxes — the
      !     transport field that produced the `u_av` evaluation state).
      if (usf) then
         do concurrent(k=1:nz, j=1:ny, i=1:nu)
            this%mass_flux_u%data(i, j, k) = ms%mass_flux_x_layer(i, j, k)
         end do
         do concurrent(k=1:nz, j=1:nv, i=1:nx)
            this%mass_flux_v%data(i, j, k) = ms%mass_flux_y_layer(i, j, k)
         end do
      else
         ! Pass 3a: u-face transport uh = u·h·dy_cu
         do concurrent(k=1:nz, j=1:ny, i=1:nu) local(h_face)
            if (i == 1) then
               h_face = h(1, j, k)
            else if (i == nu) then
               h_face = h(nx, j, k)
            else
               h_face = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
            end if
            this%mass_flux_u%data(i, j, k) = u(i, j, k)*h_face*metrics%dy_cu(i, j)
         end do

         ! Pass 3b: v-face transport vh = v·h·dx_cv
         do concurrent(k=1:nz, j=1:nv, i=1:nx) local(h_face)
            if (j == 1) then
               h_face = h(i, 1, k)
            else if (j == nv) then
               h_face = h(i, ny, k)
            else
               h_face = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
            end if
            this%mass_flux_v%data(i, j, k) = v(i, j, k)*h_face*metrics%dx_cv(i, j)
         end do

         ! ---- Porous barriers (Adcroft 2013) ----
         ! INSIDE the `else` only: the `usf` branch above copies
         ! continuity's mass fluxes, which are ALREADY narrowed, so
         ! applying the fraction again would square it.
         if (metrics%use_porous) then
            call porous_narrow_3d(nu, ny, nz, metrics%por_face_area_u, &
                                  this%mass_flux_u%data)
            call porous_narrow_3d(nx, nv, nz, metrics%por_face_area_v, &
                                  this%mass_flux_v%data)
         end if
      end if

      ! ---- Pass 4: KE at cell centres (area-weighted) ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         this%ke_centre%data(i, j, k) = 0.25_wp*metrics%iareaT(i, j)*( &
                                        metrics%areaCu(i, j)*u(i, j, k)**2 + &
                                        metrics%areaCu(i + 1, j)*u(i + 1, j, k)**2 + &
                                        metrics%areaCv(i, j)*v(i, j, k)**2 + &
                                        metrics%areaCv(i, j + 1)*v(i, j + 1, k)**2)
      end do

      ! ---- Pass 5: u-face energy tendency (2-corner Sadourny, q·vh) ----
      ! q·vh sum is a transport-weighted PV flux (m³/s); idxCu closes it to
      ! a per-length acceleration (= /dx on uniform).
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(q_S, q_N, ke_grad_x, pv_part, av_a, av_b, fv1, fv2, fv3, fv4)
         q_S = this%q_corner%data(i, j, k)
         q_N = this%q_corner%data(i, j + 1, k)
         ke_grad_x = (this%ke_centre%data(i, j, k) - &
                      this%ke_centre%data(i - 1, j, k))*metrics%idxCu(i, j)
         pv_part = 0.25_wp*( &
                   q_N*(this%mass_flux_v%data(i - 1, j + 1, k) + &
                        this%mass_flux_v%data(i, j + 1, k)) &
                   + q_S*(this%mass_flux_v%data(i - 1, j, k) + &
                          this%mass_flux_v%data(i, j, k)))* &
                   metrics%idxCu(i, j)
         if (do_bound) then
            ! BOUND_CORIOLIS (MOM6): clamp the PV flux into the range
            ! of the four neighbour (f+ζ)·v velocity-form estimates — north
            ! corner (i,j+1) × v(i-1/i,j+1); south corner (i,j) × v(i-1/i,j) —
            ! BEFORE subtracting the KE gradient.  abs_vort = q·h_corner.
            av_b = corner_abs_vort(i, j + 1, k, nx, ny, nz, q_N, h, &
                                   metrics%wet_T, metrics%areaT, use_mom6_ch)
            av_a = corner_abs_vort(i, j, k, nx, ny, nz, q_S, h, &
                                   metrics%wet_T, metrics%areaT, use_mom6_ch)
            fv1 = av_b*v(i - 1, j + 1, k)
            fv2 = av_b*v(i, j + 1, k)
            fv3 = av_a*v(i - 1, j, k)
            fv4 = av_a*v(i, j, k)
            pv_part = min(pv_part, max(max(fv1, fv2), max(fv3, fv4)))
            pv_part = max(pv_part, min(min(fv1, fv2), min(fv3, fv4)))
         end if
         this%pv_flux_x%data(i, j, k) = pv_part - ke_grad_x
      end do
      do concurrent(k=1:nz, j=1:ny)
         this%pv_flux_x%data(1, j, k) = 0.0_wp
         this%pv_flux_x%data(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 6: v-face energy tendency (2-corner Sadourny, −q·uh) ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(q_W, q_E, ke_grad_y, pv_part, av_a, av_b, fv1, fv2, fv3, fv4)
         q_W = this%q_corner%data(i, j, k)
         q_E = this%q_corner%data(i + 1, j, k)
         ke_grad_y = (this%ke_centre%data(i, j, k) - &
                      this%ke_centre%data(i, j - 1, k))*metrics%idyCv(i, j)
         pv_part = -0.25_wp*( &
                   q_E*(this%mass_flux_u%data(i + 1, j, k) + &
                        this%mass_flux_u%data(i + 1, j - 1, k)) &
                   + q_W*(this%mass_flux_u%data(i, j, k) + &
                          this%mass_flux_u%data(i, j - 1, k)))* &
                   metrics%idyCv(i, j)
         if (do_bound) then
            ! BOUND_CORIOLIS (MOM6): clamp into the four neighbour
            ! −(f+ζ)·u estimates — east corner (i+1,j) × u(i+1,j/j-1); west
            ! corner (i,j) × u(i,j/j-1) — BEFORE subtracting the KE gradient.
            av_b = corner_abs_vort(i + 1, j, k, nx, ny, nz, q_E, h, &
                                   metrics%wet_T, metrics%areaT, use_mom6_ch)
            av_a = corner_abs_vort(i, j, k, nx, ny, nz, q_W, h, &
                                   metrics%wet_T, metrics%areaT, use_mom6_ch)
            fv1 = -av_b*u(i + 1, j, k)
            fv2 = -av_b*u(i + 1, j - 1, k)
            fv3 = -av_a*u(i, j, k)
            fv4 = -av_a*u(i, j - 1, k)
            pv_part = min(pv_part, max(max(fv1, fv2), max(fv3, fv4)))
            pv_part = max(pv_part, min(min(fv1, fv2), min(fv3, fv4)))
         end if
         this%pv_flux_y%data(i, j, k) = pv_part - ke_grad_y
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%pv_flux_y%data(i, 1, k) = 0.0_wp
         this%pv_flux_y%data(i, ny + 1, k) = 0.0_wp
      end do
   end subroutine coriolis_adv_compute_tendencies_sadourny_energy

   pure function corner_abs_vort(ic, jc, k, nx, ny, nz, q_val, h, wet_T, areaT, &
                                 use_mom6_ch) result(av)
      !$acc routine seq
      !! BOUND_CORIOLIS abs_vort recovery: return `(f+ζ)` at corner `(ic,jc)`
      !! by multiplying the PV `q_val` back by the corner thickness Pass 2
      !! divided abs_vort by — recovering `(f+ζ)` to round-off.  Recomputes
      !! the SAME wet-area-weighted 4-cell `hm_num`/`hm_den`, then reconstructs
      !! the effective corner thickness with the SAME formula (and floors) the
      !! active `corner_h` variant used in Pass 2:
      !!   `cell_mean` (default): h_corner = max(hm_num/max(hm_den,H_DIV_EPS),
      !!     CORIOLIS_H_MIN_PV); av = q·h_corner.
      !!   `mom6_area` (`use_mom6_ch`): Pass 2 formed q = abs_vort·hm_den/
      !!     (hm_num + PV_VOL_NEGLECT), so the consistent inverse is
      !!     av = q·(hm_num + PV_VOL_NEGLECT)/hm_den (hm_den guarded by
      !!     H_DIV_EPS for the fully-land corner, where q≡0 ⇒ av≡0 anyway).
      !! Without matching the variant the recovered abs_vort would be slightly
      !! inconsistent with how q was made whenever BOTH bound_coriolis and
      !! corner_h="mom6_area" are on — visible only in the truly-vanishing-
      !! thickness limit.  Chosen over a persistent abs_vort buffer so the
      !! default-off knob costs ZERO memory.
      integer, intent(in) :: ic, jc, k, nx, ny, nz
      real(wp), intent(in) :: q_val
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: wet_T(nx, ny), areaT(nx, ny)
      logical, intent(in) :: use_mom6_ch
      real(wp) :: av
      integer :: iw, ie, js, jn
      real(wp) :: aSW, aSE, aNW, aNE, hm_num, hm_den, h_corner
      iw = max(1, ic - 1)
      ie = min(nx, ic)
      js = max(1, jc - 1)
      jn = min(ny, jc)
      aSW = wet_T(iw, js)*areaT(iw, js)
      aSE = wet_T(ie, js)*areaT(ie, js)
      aNW = wet_T(iw, jn)*areaT(iw, jn)
      aNE = wet_T(ie, jn)*areaT(ie, jn)
      hm_num = aSW*h(iw, js, k) + aSE*h(ie, js, k) + &
               aNW*h(iw, jn, k) + aNE*h(ie, jn, k)
      hm_den = aSW + aSE + aNW + aNE
      if (use_mom6_ch) then
         ! Algebraic inverse of Pass 2's mom6_area form (no thickness cap;
         ! PV_VOL_NEGLECT is the same 1/0 armor Pass 2 added to hm_num).
         h_corner = (hm_num + PV_VOL_NEGLECT)/max(hm_den, H_DIV_EPS)
      else
         h_corner = hm_num/max(hm_den, H_DIV_EPS)
         h_corner = max(h_corner, CORIOLIS_H_MIN_PV)
      end if
      av = q_val*h_corner
   end function corner_abs_vort

   subroutine coriolis_adv_apply_tendencies(this, ms, dt, no_wait)
      !! Per-layer forward-Euler velocity update.
      !! `no_wait` (optional, default .false.): when .true. the apply DC
      !! loops are issued on OpenACC queue 1 and the routine returns WITHOUT
      !! syncing, so a batched caller (`run_stage_split` velocity-apply chain)
      !! can pipeline the whole additive apply sequence and `!$acc wait(1)`
      !! ONCE.  Default ⇒ self-contained blocking apply (historical, safe for
      !! non-batched callers — e.g. the unsplit `run_stage`).  Not `pure`
      !! because of the async/wait directives; still functionally pure.
      type(coriolis_adv_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: no_wait
      integer :: i, j, k, nx_face, ny_uface, nx_vface, ny_face, nz
      logical :: lwait

      lwait = .true.
      if (present(no_wait)) lwait = .not. no_wait

      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      nz = ms%nz_ml

      !$acc kernels async(1)
      do concurrent(k=1:nz, j=1:ny_uface, i=1:nx_face)
         ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + &
                                      dt*this%pv_flux_x%data(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + &
                                      dt*this%pv_flux_y%data(i, j, k)
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine coriolis_adv_apply_tendencies

   pure subroutine coriolis_adv_apply_tendencies_barotropic(this, bs, dt)
      !! Forward-Euler velocity update from the tendencies the
      !! compute step wrote into pv_flux_x / pv_flux_y.
      !!   u_face_x(i, j) <- u_face_x(i, j) + dt * pv_flux_x(i, j)
      !!   v_face_y(i, j) <- v_face_y(i, j) + dt * pv_flux_y(i, j)
      !! Split-explicit RK2 (Phase 4) wraps a pair of these around an
      !! RK2 averaging pass.
      type(coriolis_adv_t), intent(in) :: this
      type(barotropic_state_t), intent(inout) :: bs
      real(wp), intent(in) :: dt
      integer :: i, j, nx_face, ny_uface, nx_vface, ny_face

      nx_face = size(bs%u_face_x, 1)
      ny_uface = size(bs%u_face_x, 2)
      do concurrent(j=1:ny_uface, i=1:nx_face)
         bs%u_face_x(i, j) = bs%u_face_x(i, j) + dt*this%pv_flux_x%data(i, j, 1)
      end do

      nx_vface = size(bs%v_face_y, 1)
      ny_face = size(bs%v_face_y, 2)
      do concurrent(j=1:ny_face, i=1:nx_vface)
         bs%v_face_y(i, j) = bs%v_face_y(i, j) + dt*this%pv_flux_y%data(i, j, 1)
      end do
   end subroutine coriolis_adv_apply_tendencies_barotropic

   pure function parse_pv_variant(name) result(code)
      !! Translate a namelist string into a `PV_VARIANT_*` code.
      !! An unrecognised string returns `PV_VARIANT_INVALID` (PR-6:
      !! fail-loud — a typo must NOT silently degrade to `SADOURNY`,
      !! which is a materially different conservation law).  `"al81"`
      !! still maps to `PV_VARIANT_AL81` (the reservation), but that
      !! code is rejected by `pv_variant_is_implemented` at configure.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("sadourny", "SADOURNY")
         code = PV_VARIANT_SADOURNY
      case ("sadourny_hk", "SADOURNY_HK", "hk", "HK")
         code = PV_VARIANT_SADOURNY_HK
      case ("sadourny_energy", "SADOURNY_ENERGY", "energy", "ENERGY")
         code = PV_VARIANT_SADOURNY_ENERGY
      case ("al81", "AL81")
         code = PV_VARIANT_AL81
      case default
         code = PV_VARIANT_INVALID
      end select
   end function parse_pv_variant

   pure function pv_variant_is_implemented(code) result(ok)
      !! `.true.` only for a Coriolis-advection variant that has a real
      !! kernel wired into `coriolis_adv_compute_tendencies`
      !! (SADOURNY / SADOURNY_HK / SADOURNY_ENERGY).  `PV_VARIANT_AL81`
      !! returns `.false.` — the constant is reserved but the
      !! Arakawa-Lamb kernel is not yet written, and the AL81 promise
      !! (simultaneous energy + enstrophy conservation) must not be
      !! silently substituted by the enstrophy-only Sadourny kernel.
      !! `PV_VARIANT_INVALID` also returns `.false.`.  The predicate is
      !! the single gate `validate_config` consumes (PR-6 fail-loud).
      integer, intent(in) :: code
      logical :: ok
      ok = (code == PV_VARIANT_SADOURNY) .or. &
           (code == PV_VARIANT_SADOURNY_HK) .or. &
           (code == PV_VARIANT_SADOURNY_ENERGY)
   end function pv_variant_is_implemented

   pure function parse_pv_adv_scheme(name) result(code)
      !! Translate a namelist string into a `PV_ADV_*` code.  An
      !! unrecognised string returns `PV_ADV_INVALID` (fail-loud — a typo
      !! must not silently degrade the PV interpolation).  `weno5`/`weno7`
      !! are implemented; they additionally require `nghost >= 3`/`4`
      !! (`pv_adv_required_nghost`), checked at configure.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("centered", "CENTERED", "centred", "center")
         code = PV_ADV_CENTERED
      case ("weno3", "WENO3")
         code = PV_ADV_WENO3
      case ("weno5", "WENO5")
         code = PV_ADV_WENO5
      case ("weno7", "WENO7")
         code = PV_ADV_WENO7
      case default
         code = PV_ADV_INVALID
      end select
   end function parse_pv_adv_scheme

   pure function pv_adv_scheme_is_implemented(code) result(ok)
      !! `.true.` for every wired PV face-interpolation scheme: `PV_ADV_CENTERED`
      !! + `PV_ADV_WENO3`/`WENO5`/`WENO7`.  weno5/weno7 additionally require a
      !! wider halo (`pv_adv_required_nghost`), checked separately at configure.
      !! `PV_ADV_INVALID` returns `.false.` (fail-loud on a typo).
      integer, intent(in) :: code
      logical :: ok
      ok = (code == PV_ADV_CENTERED) .or. (code == PV_ADV_WENO3) .or. &
           (code == PV_ADV_WENO5) .or. (code == PV_ADV_WENO7)
   end function pv_adv_scheme_is_implemented

   pure function pv_adv_required_nghost(code) result(ng)
      !! Minimum `nghost` for a PV face-interp scheme's stencil radius:
      !! weno5 (radius 3) -> 3, weno7 (radius 4) -> 4; centered/weno3 fit the
      !! nghost>=2 baseline.  Mirrors the tracer-WENO ladder's per-rung gate;
      !! `configure` fail-loud rejects an under-provisioned halo.
      integer, intent(in) :: code
      integer :: ng
      select case (code)
      case (PV_ADV_WENO5)
         ng = 3
      case (PV_ADV_WENO7)
         ng = 4
      case default
         ng = 2
      end select
   end function pv_adv_required_nghost

   pure function weno3_recon(qm1, q0, qp1, qp2, adv_vel) result(qf)
      !! 3rd-order WENO-Z reconstruction of a corner quantity onto the face
      !! between `q0` and `qp1`, upwind-biased on the sign of the advecting
      !! velocity `adv_vel` (MOM6 `weno_three_h_weight_reconstruction`).
      !! Blends a central candidate `c0` (ideal weight 2/3) with an
      !! upwind-side linear extrapolation `c1` (1/3); the WENO-Z nonlinear
      !! factor `(1+tau/b)^2` collapses the weight of whichever candidate
      !! straddles a PV front.  In smooth flow -> the fixed upwind-biased
      !! 3rd-order stencil; across a jump -> the ENO (non-oscillatory)
      !! branch.  `f`-baked absolute vorticity is passed in (rdb's
      !! vector-invariant form multiplies the result by the thickness-
      !! weighted face velocity, so there is no separate `h`-divide -- the
      !! mass-weighting lives in that velocity, not in a PV*vh product).
      !!
      !! Branchless apart from the upwind sign pick and MOM6's exact
      !! divide-guard (`|b| <= eps*tau` -> degenerate factor), both scalar
      !! per-lane predicates -- GPU-safe in `do concurrent`.
      !$acc routine seq
      real(wp), intent(in) :: qm1, q0, qp1, qp2
      real(wp), intent(in) :: adv_vel
      real(wp) :: qf
      real(wp) :: c0, c1, b0, b1, tau, w0, w1, f0, f1, sinv

      if (adv_vel > 0.0_wp) then
         c0 = 0.5_wp*(q0 + qp1)
         c1 = 0.5_wp*(-qm1 + 3.0_wp*q0)
         b0 = (q0 - qp1)**2
         b1 = (qm1 - q0)**2
      else
         c0 = 0.5_wp*(qp1 + q0)
         c1 = 0.5_wp*(-qp2 + 3.0_wp*qp1)
         b0 = (qp1 - q0)**2
         b1 = (qp2 - qp1)**2
      end if

      tau = abs(b0 - b1)
      f0 = fac_weno(tau, b0)
      f1 = fac_weno(tau, b1)
      w0 = (2.0_wp/3.0_wp)*f0
      w1 = (1.0_wp/3.0_wp)*f1
      sinv = 1.0_wp/(w0 + w1)
      qf = (w0*c0 + w1*c1)*sinv
   end function weno3_recon

   pure function fac_weno(tau, b) result(fac)
      !! MOM6 `fac_fn`: the WENO-Z nonlinear factor `(1+tau/b)^2`, EXACT-clamped
      !! (not an additive epsilon) to a dominant value when the smoothness
      !! indicator `b` is degenerate (`|b| <= eps*tau`), so a locally constant
      !! stencil takes over without a 0/0.  Shared by weno3/5/7.
      !$acc routine seq
      real(wp), intent(in) :: tau, b
      real(wp) :: fac
      if (abs(b) > PV_WENO_EPS_REL*tau) then
         fac = (1.0_wp + tau/b)**2
      else
         fac = PV_WENO_FAC_DEGEN
      end if
   end function fac_weno

   pure function beta5_0(a, b, c) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c
      real(wp) :: w
      w = a*(10.0_wp*a - 31.0_wp*b + 11.0_wp*c) + b*(25.0_wp*b - 19.0_wp*c) + 4.0_wp*c*c
   end function beta5_0

   pure function beta5_1(a, b, c) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c
      real(wp) :: w
      w = a*(4.0_wp*a - 13.0_wp*b + 5.0_wp*c) + b*(13.0_wp*b - 13.0_wp*c) + 4.0_wp*c*c
   end function beta5_1

   pure function beta5_2(a, b, c) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c
      real(wp) :: w
      w = a*(4.0_wp*a - 19.0_wp*b + 11.0_wp*c) + b*(25.0_wp*b - 31.0_wp*c) + 10.0_wp*c*c
   end function beta5_2

   pure function weno5_recon(q1, q2, q3, q4, q5, q6, adv_vel) result(qf)
      !! 5th-order WENO-Z reconstruction (MOM6 `weno_five_h_weight_reconstruction`)
      !! of a 6-point corner stencil onto the face between the two central points
      !! `q3,q4`, upwind-biased on `adv_vel`.  Three 3-point candidates blended by
      !! WENO-Z (ideal weights 3/10, 3/5, 1/10; tau = |b0-b2|).  Applied to the
      !! absolute vorticity directly (see `weno3_recon`).  Radius 3 (nghost>=3).
      !$acc routine seq
      real(wp), intent(in) :: q1, q2, q3, q4, q5, q6
      real(wp), intent(in) :: adv_vel
      real(wp) :: qf
      real(wp) :: c0, c1, c2, b0, b1, b2, tau, w0, w1, w2, sinv

      if (adv_vel > 0.0_wp) then
         c0 = (2.0_wp*q3 + 5.0_wp*q4 - q5)/6.0_wp
         b0 = beta5_0(q3, q4, q5)
         c1 = (-q2 + 5.0_wp*q3 + 2.0_wp*q4)/6.0_wp
         b1 = beta5_1(q2, q3, q4)
         c2 = (2.0_wp*q1 - 7.0_wp*q2 + 11.0_wp*q3)/6.0_wp
         b2 = beta5_2(q1, q2, q3)
      else
         c0 = (2.0_wp*q4 + 5.0_wp*q3 - q2)/6.0_wp
         b0 = beta5_0(q4, q3, q2)
         c1 = (-q5 + 5.0_wp*q4 + 2.0_wp*q3)/6.0_wp
         b1 = beta5_1(q5, q4, q3)
         c2 = (2.0_wp*q6 - 7.0_wp*q5 + 11.0_wp*q4)/6.0_wp
         b2 = beta5_2(q6, q5, q4)
      end if

      tau = abs(b0 - b2)
      w0 = 0.3_wp*fac_weno(tau, b0)
      w1 = 0.6_wp*fac_weno(tau, b1)
      w2 = 0.1_wp*fac_weno(tau, b2)
      sinv = 1.0_wp/((w0 + w1) + w2)
      qf = ((w0*c0) + (w1*c1) + (w2*c2))*sinv
   end function weno5_recon

   pure function beta7_0(a, b, c, d) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c, d
      real(wp) :: w
      w = a*(2.107_wp*a - 9.402_wp*b + 7.042_wp*c - 1.854_wp*d) &
          + b*(11.003_wp*b - 17.246_wp*c + 4.642_wp*d) &
          + c*(7.043_wp*c - 3.882_wp*d) + 0.547_wp*d*d
   end function beta7_0

   pure function beta7_1(a, b, c, d) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c, d
      real(wp) :: w
      w = a*(0.547_wp*a - 2.522_wp*b + 1.922_wp*c - 0.494_wp*d) &
          + b*(3.443_wp*b - 5.966_wp*c + 1.602_wp*d) &
          + c*(2.843_wp*c - 1.642_wp*d) + 0.267_wp*d*d
   end function beta7_1

   pure function beta7_2(a, b, c, d) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c, d
      real(wp) :: w
      w = a*(0.267_wp*a - 1.642_wp*b + 1.602_wp*c - 0.494_wp*d) &
          + b*(2.843_wp*b - 5.966_wp*c + 1.922_wp*d) &
          + c*(3.443_wp*c - 2.522_wp*d) + 0.547_wp*d*d
   end function beta7_2

   pure function beta7_3(a, b, c, d) result(w)
      !$acc routine seq
      real(wp), intent(in) :: a, b, c, d
      real(wp) :: w
      w = a*(0.547_wp*a - 3.882_wp*b + 4.642_wp*c - 1.854_wp*d) &
          + b*(7.043_wp*b - 17.246_wp*c + 7.042_wp*d) &
          + c*(11.003_wp*c - 9.402_wp*d) + 2.107_wp*d*d
   end function beta7_3

   pure function weno7_recon(q1, q2, q3, q4, q5, q6, q7, q8, adv_vel) result(qf)
      !! 7th-order WENO-Z reconstruction (MOM6 `weno_seven_h_weight_reconstruction`)
      !! of an 8-point corner stencil onto the face between the two central points
      !! `q4,q5`, upwind-biased on `adv_vel`.  Four 4-point candidates blended by
      !! WENO-Z (ideal weights 4/35, 18/35, 12/35, 1/35; Balsara-Shu smoothness;
      !! tau = |(b0-b3) + 3(b1-b2)|).  Applied to the absolute vorticity directly.
      !! Radius 4 (nghost>=4).
      !$acc routine seq
      real(wp), intent(in) :: q1, q2, q3, q4, q5, q6, q7, q8
      real(wp), intent(in) :: adv_vel
      real(wp) :: qf
      real(wp) :: c0, c1, c2, c3, b0, b1, b2, b3, tau, w0, w1, w2, w3, sinv

      if (adv_vel > 0.0_wp) then
         c0 = (6.0_wp*q4 + 26.0_wp*q5 - 10.0_wp*q6 + 2.0_wp*q7)/24.0_wp
         b0 = beta7_0(q4, q5, q6, q7)
         c1 = (-2.0_wp*q3 + 14.0_wp*q4 + 14.0_wp*q5 - 2.0_wp*q6)/24.0_wp
         b1 = beta7_1(q3, q4, q5, q6)
         c2 = (2.0_wp*q2 - 10.0_wp*q3 + 26.0_wp*q4 + 6.0_wp*q5)/24.0_wp
         b2 = beta7_2(q2, q3, q4, q5)
         c3 = (-6.0_wp*q1 + 26.0_wp*q2 - 46.0_wp*q3 + 50.0_wp*q4)/24.0_wp
         b3 = beta7_3(q1, q2, q3, q4)
      else
         c0 = (6.0_wp*q5 + 26.0_wp*q4 - 10.0_wp*q3 + 2.0_wp*q2)/24.0_wp
         b0 = beta7_0(q5, q4, q3, q2)
         c1 = (-2.0_wp*q6 + 14.0_wp*q5 + 14.0_wp*q4 - 2.0_wp*q3)/24.0_wp
         b1 = beta7_1(q6, q5, q4, q3)
         c2 = (2.0_wp*q7 - 10.0_wp*q6 + 26.0_wp*q5 + 6.0_wp*q4)/24.0_wp
         b2 = beta7_2(q7, q6, q5, q4)
         c3 = (-6.0_wp*q8 + 26.0_wp*q7 - 46.0_wp*q6 + 50.0_wp*q5)/24.0_wp
         b3 = beta7_3(q8, q7, q6, q5)
      end if

      tau = abs((b0 - b3) + 3.0_wp*(b1 - b2))
      w0 = (4.0_wp/35.0_wp)*fac_weno(tau, b0)
      w1 = (18.0_wp/35.0_wp)*fac_weno(tau, b1)
      w2 = (12.0_wp/35.0_wp)*fac_weno(tau, b2)
      w3 = (1.0_wp/35.0_wp)*fac_weno(tau, b3)
      sinv = 1.0_wp/((w0 + w1) + (w2 + w3))
      qf = ((w0*c0) + (w1*c1) + (w2*c2) + (w3*c3))*sinv
   end function weno7_recon

   pure function coriolis_adv_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the Coriolis-advection slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(coriolis_adv_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_corner) &
               + this%q_corner%bytes() &
               + this%ke_centre%bytes() &
               + this%pv_flux_x%bytes() &
               + this%pv_flux_y%bytes() &
               + this%mass_flux_u%bytes() &
               + this%mass_flux_v%bytes()
   end function coriolis_adv_bytes

end module rdb_coriolis_adv
