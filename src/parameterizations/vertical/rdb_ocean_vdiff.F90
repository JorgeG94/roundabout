!! Backward-Euler vertical diffusion solver for the ocean multilayer
!! C-grid.  Two public entry points:
!!
!!   * `vdiff_apply_momentum` — implicit vertical viscosity
!!     on `u_face_x_layer` and `v_face_y_layer`.
!!   * `vdiff_apply_tracers` — implicit vertical diffusivity on every
!!     registered tracer (operates on the concentration `c`;
!!     `h_layer` is untouched).  Takes TWO diffusivity
!!     sources, `kt_source` (temperature) and `ks_source` (salinity +
!!     every passive tracer, MOM6's `Kd_salt` convention) — see the
!!     routine's own docstring for the three-case dispatch.
!!
!! Both routines build a tridiagonal system per cell column and
!! solve in-place via the Thomas algorithm:
!!
!!   (1 + α_k + β_k) T_k^{n+1} - α_k T_{k+1}^{n+1} - β_k T_{k-1}^{n+1} = T_k^n
!!
!! with
!!   α_k = dt * K_v / (h_k * dz_face_{k+1})
!!   β_k = dt * K_v / (h_k * dz_face_k)
!!   dz_face_k = 0.5 * (h_{k-1} + h_k)
!!
!! Bed (`k = 1`) sets `β_1 = 0`; surface (`k = nz`) sets `α_nz = 0`
!! → closed top and bottom (no surface heat / wind-stress flux here;
!! those are separate kernels that add to the right-hand side).
!!
!! Backward-Euler is unconditionally stable, so the per-step
!! constraint is purely accuracy — `dt * K_v / h² ≲ 1` for the
!! diffusion-time scale to be resolved.  KPP-style large K_v in the
!! surface boundary layer (~1e-1 m²/s) over a 1 m layer at dt = 10 min
!! ⇒ `dt*K_v/h² ~ 60` — still stable, but the resolved decay rate
!! is heavily damped (which is the right thing under strong mixing).
!!
!! Per-column work is serial in k (Thomas recurrence); parallelism
!! is `do concurrent (j, i)`.
module rdb_ocean_vdiff
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: H_DIV_EPS, wp, H_VANISHED, GRAVITY
#else
   use rdb_constants, only: H_DIV_EPS, wp, H_VANISHED, NZ_STACK_MAX, GRAVITY
#endif
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t, eos_density_derivs
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_bbl_mixing, only: ocean_bbl_mixing_t
   use pic_logger, only: logger => global_logger
   implicit none
   private
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ocean_vdiff_t
   public :: vdiff_apply_momentum
   public :: vdiff_apply_tracers
   public :: vdiff_set_viscous_bbl
   public :: vdiff_bbl_configure
   public :: face_thick
   public :: BBL_FORM_LINEAR, BBL_FORM_QUADRATIC

   integer, parameter :: BBL_FORM_LINEAR = 1
      !! `&ocean_bdrag_nml form = "linear"` under the MOM6 BBL: MOM6
      !! `LINEAR_DRAG`, `u* = sqrt(cd)·DRAG_BG_VEL` (see `bbl_cd`).
   integer, parameter :: BBL_FORM_QUADRATIC = 2
      !! `form = "quadratic"`: `u* = sqrt(cd)·u_bbl`, `u_bbl` the HBBL-mean
      !! speed with the background velocity added in quadrature.

   type :: ocean_vdiff_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      real(wp) :: K_v_momentum = 0.0_wp
         !! Constant momentum vertical viscosity (m^2/s).  Zero is a
         !! no-op.
      real(wp) :: K_v_tracer = 0.0_wp
         !! Constant tracer vertical diffusivity (m^2/s).  Zero is a
         !! no-op.
      logical :: implicit_stress = .false.
         !! Fold the surface wind stress into the backward-Euler vertical-
         !! friction tridiagonal as a Neumann top-BC right-hand-side term
         !! (`&ocean_vdiff_nml implicit_stress`).  Roundabout is bottom-up, so
         !! the SURFACE is `k = nz`: the kinematic stress `τ/ρ₀` adds to the
         !! `rhs(:, :, nz)` row only (the diagonal is unchanged), filtered
         !! through the implicit operator together with the interior shear.
         !! When `.true.` the explicit surface-stress pre-solve apply in the
         !! dyn run-stage is suppressed (the driver gates it) so the forcing
         !! is not double-counted.  Default `.false.` ⇒ matrix + RHS built
         !! exactly as before ⇒ bit-identical to the explicit path.
      logical :: implicit_drag = .false.
         !! Fold the quadratic / linear bottom drag into the backward-Euler
         !! vertical-friction tridiagonal as a stress bottom-BC diagonal
         !! coupling (`&ocean_vdiff_nml implicit_drag`).  Roundabout bed is
         !! `k = 1`: the drag rate `λ_bot` adds to the `b_diag(:, :, 1)`
         !! diagonal only (a positive add ⇒ unconditionally stable on thin
         !! bottom layers where the explicit `u·(1 − dt·c_d|U|/h)` reverses
         !! sign once `dt·c_d|U|/h > 2`).  `λ_bot` is the SAME Rayleigh rate
         !! (`c_d·|U_bbl|/h_1` quadratic, `r` linear, `|U|` frozen at uⁿ)
         !! the bottom-drag slot already forms; the row is already
         !! normalized by `h_1` so the diagonal increment is `dt·λ_bot`
         !! directly.  Mutually exclusive with `&ocean_bdrag_nml implicit`
         !! (configure fails loud); the explicit drag apply is gated off
         !! when on.  Default `.false.` ⇒ bit-identical.
      logical :: implicit_top_drag = .false.
         !! Fold the ICE-SHELF TOP drag into the `k = nz` DIAGONAL
         !! (`&ocean_vdiff_nml implicit_top_drag`) instead of the explicit
         !! pre-solve add, and MASK the wind-stress RHS off on the faces
         !! the ice covers.  The mirror of `implicit_drag` at the other
         !! end of the column: a drag is a diagonal term and the wind is
         !! an RHS term, so the two share the surface row without either
         !! being approximated, and the cover mask is what makes the
         !! sharing physical (there is no atmosphere under a shelf).
         !! `λ_top` is the SAME Rayleigh rate (`C_d·|U|/h_nz` quadratic,
         !! `r` linear, `|U|` frozen at uⁿ) the top-drag slot forms, and
         !! it arrives already cover- and wet-masked.  Requires
         !! `&ocean_tdrag_nml enable`; mutually exclusive with
         !! `&ocean_tdrag_nml implicit` and with `htbl > 0`; all fail loud
         !! at configure.  The explicit top-drag apply is gated off when
         !! on.  Default `.false.` ⇒ bit-identical.
      logical :: hvel_mom6 = .false.
         !! MOM6 `HARMONIC_VISC = True` parity for the MOMENTUM face thickness
         !! (`hvel`).  Roundabout's `h_u` is the
         !! ARITHMETIC mean unconditionally while MOM6's is the HARMONIC mean
         !! blended back toward arithmetic near the bed when the flow runs
         !! thick->thin -- and MOM6's `h_shear` is then the ARITHMETIC mean of
         !! those hvels.  The two means are effectively SWAPPED relative to
         !! MOM6.  For a grounded sliver `harm(1e-10, 20) ~ 2e-10` against
         !! `arith ~ 10`: eleven orders, and that gap IS the ~1e19 suppression
         !! that renders the spurious sliver PGF harmless in MOM6.
         !!
         !! **Both halves must move together.** LAGRANGIAN_PGF_BUG.md 6.2
         !! records harmonic-`h_u` alone (and harmonic-`h_u` + arithmetic-`dz`)
         !! each making things WORSE -- because neither carried `botfn`.  Plain
         !! harmonic over-suppresses asymmetrically.
      real(wp) :: hbbl_visc = 10.0_wp
         !! Bottom-boundary-layer scale for the `botfn` blend (MOM6 `HBBL`,
         !! 10.0 m in the double_gyre reference) when the BBL glue is off.
         !! Under the per-face glue the blend is normalised by each face's
         !! `bbl_thick` instead (MOM6 `I_Hbbl = 1/bbl_thick`), and this is
         !! only the HBBL fallback for a drag configured with
         !! `&ocean_bdrag_nml hbbl = 0` (MOM6 has ONE `HBBL`).
      logical :: use_harmonic = .false.
         !! MOM6 HARMONIC_VISC analogue.  When `.true.` the implicit
         !! vdiff solver uses the harmonic mean of adjacent layer
         !! thicknesses in the face-thickness denominator instead of
         !! the arithmetic mean.  Better-conditioned at thin /
         !! vanishing layers — `harm = 2·h_a·h_b / (h_a + h_b)` stays
         !! small when one neighbour is small, whereas
         !! `arith = 0.5·(h_a + h_b)` is dominated by the thicker
         !! neighbour and produces stiff coefficients that the Thomas
         !! solve can't handle gracefully.  Default `.false.`
         !! preserves the arithmetic-mean behaviour bit-identically.
      logical :: bbl_glue = .false.
         !! MOM6 `bottomdraglaw` parity for the interface coupling
         !! (`&ocean_vdiff_nml bbl_glue`, PGF_BUG.md §9).  MOM6 holds a
         !! layered basin at rest NOT by computing a clean PGF — its PFu
         !! over grounded/sliver layers is identical to ours — but by
         !! absorbing it viscously every step (`du_dt_visc = −PFu` to
         !! machine precision, measured): `find_coupling_coef` raises the
         !! interface coupling to `a_cpl = kv_bbl/h_shear` within the
         !! bottom boundary layer (botfn weight, `z_i` accumulated from
         !! HARMONIC thicknesses so grounded stacks sit at z≈0), and the
         !! BED coupling is the piston `kv_bbl/(h₁/2)` — unbounded as the
         !! bottom layer thins (measured a_cpl 3e7–1.3e8 m/s on the
         !! rest-state reproducer).  This knob
         !! ports both halves into the momentum tridiagonal:
         !!   * interface kv → `kv + (kv_bbl − KV)·botfn` (`KV` =
         !!     `kv_bbl_bg`), with `h_shear` capped toward `bbl_thick` under
         !!     the same botfn weight;
         !!   * bed row → `dt·kv_bbl/(hf₁·(min(hvel₁/2, bbl_thick)))` — THE
         !!     bed sink whenever the glue is on: it replaces the
         !!     `dt·λ_bot` Rayleigh fold, and the driver skips the explicit
         !!     drag apply (the explicit `du_drag` still feeds the
         !!     barotropic `F_slow`, the split `implicit_drag` takes).
         !! `kv_bbl` and `bbl_thick` are PER FACE: `vdiff_set_viscous_bbl`
         !! (MOM6 `set_viscous_BBL`, quadratic or linear drag law, KW99
         !! rotation/stratification-limited thickness) fills them once per
         !! outer step when `bbl_per_face` (configure: a bottom drag is
         !! configured).  A hand-built slot without `bbl_per_face` gets the
         !! historical constants `bbl_piston·hbbl_visc` / `hbbl_visc`.
         !! Requires `hvel_mom6` (supplies the height-above-bed
         !! bookkeeping), enforced at configure.  With no bottom drag
         !! configured the configure latch turns it off (nothing to glue).
      real(wp) :: bbl_piston = 3.0e-4_wp
         !! BBL drag piston velocity u* (m/s) for `bbl_glue` —
         !! MOM6 `CDRAG·DRAG_BG_VEL` (0.003·0.1 with reference defaults).
         !! `kv_bbl = bbl_piston·hbbl_visc`.
      logical :: hvel_upwind = .true.
         !! Near-bed upwind (arithmetic-donor) blend inside the
         !! `hvel_mom6` face-thickness build.  `.true.` = MOM6 parity
         !! (the historical hvel_mom6 behaviour, bit-identical).
         !! `.false.` = pure harmonic hvel: the blend's u-sign test
         !! flip-flops on roundoff velocities at (near-)rest and
         !! collapses the BBL glue at flipped faces (PGF_BUG.md §9.8) —
         !! the rest-state envelope runs with it off.
      logical :: hvel_harmonic = .false.
         !! Which MOM6 face-thickness branch `hvel_mom6` builds.  `.false.`
         !! (default) is MOM6's DEFAULT, `HARMONIC_VISC = False`
         !! (`vertvisc_coef`): the face thickness is the ARITHMETIC mean,
         !! blended toward the harmonic mean near the bed when the flow runs
         !! from the thin side to the thick one, and the height above the
         !! bed is `max(zh, z_clear)` — `z_clear` the height of the higher of
         !! the two cells' interfaces above the SHALLOWER of the two beds.
         !! So every face layer below the shallower bed of a step sits at
         !! `z_i ≈ 0`, inside the bottom boundary layer, and the BBL glue
         !! couples it.  `.true.` is MOM6 `HARMONIC_VISC = True`: harmonic
         !! face thickness, blended toward arithmetic when the flow runs
         !! thick -> thin, height above bed from harmonic thicknesses alone
         !! (the historical roundabout `hvel_mom6`, gated by `hvel_upwind`).
      logical :: bbl_per_face = .false.
         !! The MOM6 bottom boundary layer (`set_viscous_BBL`) is computed
         !! per face once per outer step by `vdiff_set_viscous_bbl` into
         !! `kv_bbl_u/v` and `bbl_thick_u/v`, which the glue then reads.
         !! Set at configure when `bbl_glue` is on AND a bottom drag is
         !! configured with `&ocean_bdrag_nml hbbl > 0` (MOM6 requires
         !! HBBL).  `.false.` ⇒ the glue uses the historical constants
         !! `bbl_piston·hbbl_visc` / `hbbl_visc`, refilled each call (the
         !! hand-built test path).
      integer :: bbl_form = BBL_FORM_QUADRATIC
         !! Drag law of the BBL (`BBL_FORM_*`).
      real(wp) :: bbl_cd = 0.0_wp
         !! MOM6 `CDRAG`.  Quadratic: `&ocean_bdrag_nml cd`.  Linear: the
         !! equivalent `r·hbbl/bg_vel`, so that `CDRAG·DRAG_BG_VEL = r·hbbl`
         !! — the stress of the HBBL-distributed linear drag, `r·hbbl·u`.
      real(wp) :: bbl_hbbl = 0.0_wp
         !! MOM6 `HBBL` (`&ocean_bdrag_nml hbbl`): the thickness over which
         !! the near-bottom velocity is averaged for the drag law.
      real(wp) :: bbl_bg_vel = 0.0_wp
         !! MOM6 `DRAG_BG_VEL` (`&ocean_bdrag_nml bg_vel`).
      real(wp) :: bbl_thick_min = 0.0_wp
         !! MOM6 `BBL_THICK_MIN` (`&ocean_bdrag_nml bbl_thick_min`; MOM6
         !! default 0).
      logical :: bbl_rino_cap = .false.
         !! MOM6 `RiNo_mix` (= kappa-shear on): the BBL is capped at
         !! `0.5·HBBL`.
      real(wp) :: bbl_rho0 = 1035.0_wp
         !! Boussinesq reference density of the BBL stratification measure.
      real(wp) :: kv_bbl_bg = 1.0e-4_wp
         !! MOM6 `KV`: the background viscosity the BBL viscosity REPLACES
         !! within the boundary layer,
         !! `Kv_tot = Kv_tot + (kv_bbl − KV)·botfn` (`find_coupling_coef`) —
         !! the shear / boundary-layer contributions are kept.  Set from
         !! `&ocean_vmix_nml kv_bg`.
      real(wp), allocatable :: kv_bbl_u(:, :)
         !! BBL viscosity `kv_bbl` (m²/s) at u-faces, `(nx+1, ny)`.
      real(wp), allocatable :: kv_bbl_v(:, :)
         !! BBL viscosity at v-faces, `(nx, ny+1)`.
      real(wp), allocatable :: bbl_thick_u(:, :)
         !! BBL thickness (m) at u-faces, `(nx+1, ny)`.
      real(wp), allocatable :: bbl_thick_v(:, :)
         !! BBL thickness (m) at v-faces, `(nx, ny+1)`.
      real(wp), allocatable :: bbl_conc_t(:, :, :)
         !! Cell temperature CONCENTRATION `(nx, ny, nz)` for the BBL
         !! stratification measure, read through the I1′ column rule
         !! (`rdb_vl_column_conc`: a filler reads its donor).  Workspace of
         !! `vdiff_set_viscous_bbl`; `(1,1,1)` unless `bbl_per_face`.
      real(wp), allocatable :: bbl_conc_s(:, :, :)
         !! Cell salinity concentration, as `bbl_conc_t`.
      type(ocean_bbl_mixing_t) :: bbl_mix
         !! Law-of-the-wall BBL diapycnal mixing (`&ocean_bbl_mixing_nml`).
         !! Lives here because it consumes the per-face BBL (`kv_bbl_*`,
         !! `bbl_thick_*`, `bbl_cd`) this slot owns and rides every call
         !! path that already carries `vd`.  Allocated by its configure
         !! only when enabled (before `enter_data`).
      logical :: zlevel_faces = .false.
         !! `&vcoord_nml zfixed_closed_faces` — z-level partial steps.
         !! Latched by `configure_ocean_closed_faces`, forwarded as a
         !! plain scalar to `diffuse_velocity_columns_impl`, where it
         !! switches the face column to `min(h_L, h_R)` and CUTS the
         !! tridiagonal coupling across any interface touching a closed
         !! (inert-filler-on-one-side) layer.  See that routine's
         !! `zlevel_faces` docstring for the full argument.  Default
         !! `.false.` ⇒ bit-identical.

      ! ---- Tridiagonal scratch (cell-centred for tracers) ----
      ! Reused across all column solves in one tracer call.  The
      ! Thomas algorithm modifies `c_diag` and `rhs` in place during
      ! the forward sweep, so a separate `c_prime` / `d_prime` pair
      ! isn't needed.
      type(scratch_3d_buffer_t) :: a_diag_t
         !! Sub-diagonal for the tracer solve.  Shape (nx, ny, nz).
      type(scratch_3d_buffer_t) :: b_diag_t
         !! Main diagonal.
      type(scratch_3d_buffer_t) :: c_diag_t
         !! Super-diagonal (overwritten by `c_prime`).
      type(scratch_3d_buffer_t) :: rhs_t
         !! Right-hand side / solution.

      ! ---- Momentum tridiagonal scratch (face-located) ----
      ! u-face shape is (nx+1, ny, nz); v-face shape is (nx, ny+1, nz).
      ! Separate buffers keep the descriptor extents exact so the
      ! `do concurrent` indexing doesn't drift past the allocated
      ! footprint.
      type(scratch_3d_buffer_t) :: a_diag_u, b_diag_u, c_diag_u, rhs_u
         !! East-face (u) tridiagonal scratch.
      type(scratch_3d_buffer_t) :: a_diag_v, b_diag_v, c_diag_v, rhs_v
         !! North-face (v) tridiagonal scratch.

      ! ---- Cell-centred diffusivity workspace ----
      ! Shape (nx, ny, nz+1).  Filled with `K_v_*` (broadcast from
      ! the scalar) when the apply routines are called without a
      ! `kv_source` argument; left alone when the caller passes a
      ! closure-produced 3D field directly.  k=1 is the bed
      ! interface (forced to zero — closed bed BC), k=nz+1 is the
      ! free surface (also zero — closed top).
      type(scratch_3d_buffer_t) :: kv_scalar_buf
         !! Diffusivity workspace used when no 3D source is provided.
   contains
      procedure, non_overridable :: init => ocean_vdiff_init
      procedure, non_overridable :: destroy => ocean_vdiff_destroy
      procedure, non_overridable :: enter_data => ocean_vdiff_enter_data
      procedure, non_overridable :: exit_data => ocean_vdiff_exit_data
      procedure, non_overridable :: bytes => ocean_vdiff_bytes
   end type ocean_vdiff_t

contains

#include "rdb_vanished_layer.inc"
#include "rdb_tracer_conc.inc"

   subroutine ocean_vdiff_init(this, grid, nz_ml)
      class(ocean_vdiff_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      ! Tracer scratch: cell-centred (nx, ny, nz)
      call this%a_diag_t%init(nx, ny, nz, "ocean_vdiff_a_diag_t")
      call this%b_diag_t%init(nx, ny, nz, "ocean_vdiff_b_diag_t")
      call this%c_diag_t%init(nx, ny, nz, "ocean_vdiff_c_diag_t")
      call this%rhs_t%init(nx, ny, nz, "ocean_vdiff_rhs_t")

      ! u-face scratch: (nx+1, ny, nz)
      call this%a_diag_u%init(nx + 1, ny, nz, "ocean_vdiff_a_diag_u")
      call this%b_diag_u%init(nx + 1, ny, nz, "ocean_vdiff_b_diag_u")
      call this%c_diag_u%init(nx + 1, ny, nz, "ocean_vdiff_c_diag_u")
      call this%rhs_u%init(nx + 1, ny, nz, "ocean_vdiff_rhs_u")

      ! v-face scratch: (nx, ny+1, nz)
      call this%a_diag_v%init(nx, ny + 1, nz, "ocean_vdiff_a_diag_v")
      call this%b_diag_v%init(nx, ny + 1, nz, "ocean_vdiff_b_diag_v")
      call this%c_diag_v%init(nx, ny + 1, nz, "ocean_vdiff_c_diag_v")
      call this%rhs_v%init(nx, ny + 1, nz, "ocean_vdiff_rhs_v")

      ! Cell-centred diffusivity workspace: (nx, ny, nz+1).
      call this%kv_scalar_buf%init(nx, ny, nz + 1, "ocean_vdiff_kv_scalar_buf")

      ! Per-face bottom boundary layer (2-D, cheap): filled by
      ! `vdiff_set_viscous_bbl`, or with the glue constants.
      allocate (this%kv_bbl_u(nx + 1, ny), source=0.0_wp)
      allocate (this%kv_bbl_v(nx, ny + 1), source=0.0_wp)
      allocate (this%bbl_thick_u(nx + 1, ny), source=0.0_wp)
      allocate (this%bbl_thick_v(nx, ny + 1), source=0.0_wp)
      ! Concentration workspace: sized by `vdiff_bbl_configure` (before
      ! `enter_data`) only when the per-face BBL is on.
      allocate (this%bbl_conc_t(1, 1, 1), source=0.0_wp)
      allocate (this%bbl_conc_s(1, 1, 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_vdiff_init

   subroutine ocean_vdiff_destroy(this)
      class(ocean_vdiff_t), intent(inout) :: this
      this%is_init = .false.
      call this%a_diag_t%destroy()
      call this%b_diag_t%destroy()
      call this%c_diag_t%destroy()
      call this%rhs_t%destroy()
      call this%a_diag_u%destroy()
      call this%b_diag_u%destroy()
      call this%c_diag_u%destroy()
      call this%rhs_u%destroy()
      call this%a_diag_v%destroy()
      call this%b_diag_v%destroy()
      call this%c_diag_v%destroy()
      call this%rhs_v%destroy()
      call this%kv_scalar_buf%destroy()
      if (allocated(this%kv_bbl_u)) deallocate (this%kv_bbl_u)
      if (allocated(this%kv_bbl_v)) deallocate (this%kv_bbl_v)
      if (allocated(this%bbl_thick_u)) deallocate (this%bbl_thick_u)
      if (allocated(this%bbl_thick_v)) deallocate (this%bbl_thick_v)
      if (allocated(this%bbl_conc_t)) deallocate (this%bbl_conc_t)
      if (allocated(this%bbl_conc_s)) deallocate (this%bbl_conc_s)
      call this%bbl_mix%destroy()
   end subroutine ocean_vdiff_destroy

   subroutine ocean_vdiff_enter_data(this)
      class(ocean_vdiff_t), intent(inout) :: this
      select type (this)
      type is (ocean_vdiff_t)
         call ocean_vdiff_enter_data_impl(this)
      end select
   end subroutine ocean_vdiff_enter_data

   subroutine ocean_vdiff_enter_data_impl(this)
      type(ocean_vdiff_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%a_diag_t)
      call scratch_3d_buffer_enter_data_impl(this%b_diag_t)
      call scratch_3d_buffer_enter_data_impl(this%c_diag_t)
      call scratch_3d_buffer_enter_data_impl(this%rhs_t)
      call scratch_3d_buffer_enter_data_impl(this%a_diag_u)
      call scratch_3d_buffer_enter_data_impl(this%b_diag_u)
      call scratch_3d_buffer_enter_data_impl(this%c_diag_u)
      call scratch_3d_buffer_enter_data_impl(this%rhs_u)
      call scratch_3d_buffer_enter_data_impl(this%a_diag_v)
      call scratch_3d_buffer_enter_data_impl(this%b_diag_v)
      call scratch_3d_buffer_enter_data_impl(this%c_diag_v)
      call scratch_3d_buffer_enter_data_impl(this%rhs_v)
      call scratch_3d_buffer_enter_data_impl(this%kv_scalar_buf)
      !$acc enter data copyin(this%kv_bbl_u, this%kv_bbl_v, this%bbl_thick_u, this%bbl_thick_v)
      !$acc enter data copyin(this%bbl_conc_t, this%bbl_conc_s)
      if (this%bbl_mix%enable) call this%bbl_mix%enter_data()
   end subroutine ocean_vdiff_enter_data_impl

   subroutine ocean_vdiff_exit_data(this)
      class(ocean_vdiff_t), intent(inout) :: this
      select type (this)
      type is (ocean_vdiff_t)
         call ocean_vdiff_exit_data_impl(this)
      end select
   end subroutine ocean_vdiff_exit_data

   subroutine ocean_vdiff_exit_data_impl(this)
      type(ocean_vdiff_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%a_diag_t)
      call scratch_3d_buffer_exit_data_impl(this%b_diag_t)
      call scratch_3d_buffer_exit_data_impl(this%c_diag_t)
      call scratch_3d_buffer_exit_data_impl(this%rhs_t)
      call scratch_3d_buffer_exit_data_impl(this%a_diag_u)
      call scratch_3d_buffer_exit_data_impl(this%b_diag_u)
      call scratch_3d_buffer_exit_data_impl(this%c_diag_u)
      call scratch_3d_buffer_exit_data_impl(this%rhs_u)
      call scratch_3d_buffer_exit_data_impl(this%a_diag_v)
      call scratch_3d_buffer_exit_data_impl(this%b_diag_v)
      call scratch_3d_buffer_exit_data_impl(this%c_diag_v)
      call scratch_3d_buffer_exit_data_impl(this%rhs_v)
      call scratch_3d_buffer_exit_data_impl(this%kv_scalar_buf)
      !$acc exit data delete(this%kv_bbl_u, this%kv_bbl_v, this%bbl_thick_u, this%bbl_thick_v)
      if (this%bbl_mix%enable) call this%bbl_mix%exit_data()
      !$acc exit data delete(this%bbl_conc_t, this%bbl_conc_s)
   end subroutine ocean_vdiff_exit_data_impl

   subroutine vdiff_apply_momentum(grid, this, ms, dt, kv_source, &
                                   tau_u, tau_v, lambda_bot_u, lambda_bot_v, rho0, &
                                   visc_rem_u, visc_rem_v, remnant_only, &
                                   kv_corner_source, kv_corner_prandtl, &
                                   lambda_top_u, lambda_top_v, cover_u, cover_v)
      !! Backward-Euler vertical viscosity applied to
      !! `ms%u_face_x_layer` and `ms%v_face_y_layer`.  Per-face
      !! `h_face` averaged from the two abutting cell columns.
      !!
      !! Diffusivity source: either a 3D cell-centred field from a
      !! closure (PP81 / KPP) passed as `kv_source` (shape (nx, ny,
      !! nz+1) — interface-located, k=1 bed, k=nz+1 surface) OR the
      !! scalar `this%K_v_momentum` broadcast into the
      !! `kv_scalar_buf` workspace.  No-op when neither path
      !! supplies a non-zero diffusivity.
      !!
      !! Corner viscosity add-on (`kv_corner_source`, optional): a
      !! CORNER-staggered interface viscosity `(nx+1, ny+1, nz+1)`
      !! (corner (i,j) = SW corner of cell (i,j); k=1 bed interface,
      !! k=nz+1 surface), scaled by `kv_corner_prandtl` and ADDED to
      !! the face viscosity as the 2-point average of the face's two
      !! END corners — u-face (i,j) reads corners (i,j)/(i,j+1), v-face
      !! (i,j) reads corners (i,j)/(i+1,j).  This is the vertex
      !! kappa-shear `Kv` seam (JHL08 vertex form): the corner
      !! viscosity reaches the momentum solve WITHOUT passing through a
      !! tracer point, so it is NOT smoothed by a corner->centre->face
      !! round trip.  The add lands BEFORE the BBL-glue transform (the
      !! reference formulation folds all viscosity contributions into
      !! `Kv_tot` ahead of the coupling-coefficient bottom-boundary
      !! blend).  Absent ⇒ bit-identical.
      !!
      !! Implicit surface-stress / bottom-drag folding (optional, gated
      !! on `this%implicit_stress` / `this%implicit_drag`):
      !!   * `tau_u(nu, nv)` / `tau_v(nu, nv)` — wind stress (N/m²) on the
      !!     respective faces.  Divided by `rho0` and added to the surface
      !!     (`k = nz`) RHS row.  Only read when `implicit_stress`.
      !!   * `lambda_bot_u` / `lambda_bot_v` — bottom-drag Rayleigh rate
      !!     `λ` (1/s) on the respective faces, added to the bed row
      !!     (`k = kb_live`, the first LIVE row at/above the static
      !!     `k_bot_u/v` — #178) diagonal as `dt·λ`.  Only read when
      !!     `implicit_drag`.
      !!   * `lambda_top_u` / `lambda_top_v` + `cover_u` / `cover_v` —
      !!     the ICE-SHELF TOP-drag twin (`implicit_top_drag`): the rate
      !!     is added to the SURFACE (`k = nz`) diagonal as `dt·λ`, and
      !!     the wind-stress RHS on that same row is scaled by
      !!     `(1 − cover)` so an ice-covered face takes the drag and NOT
      !!     the wind.  All four are read only when `implicit_top_drag`
      !!     AND all four are supplied.
      !! All optional + device-resident; absent ⇒ the matrix is built
      !! exactly as before (bit-identical).
      !!
      !! `visc_rem_u(nu, nv, nz)` / `visc_rem_v(nu, nv, nz)` — optional
      !! viscous-remnant γ_k output (the `bt_work%visc_rem_u/v` seam,
      !! PR-19).  Both present ⇒ the SAME already-factorized tridiagonal
      !! solved a second time with RHS ≡ 1, written into these arrays
      !! (see `diffuse_velocity_columns_impl`).  Absent ⇒ no remnant
      !! work is done and the momentum answer is unperturbed
      !! (bit-identical) — this is a pure add-on, not a new physics
      !! path.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vdiff_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      ! EXPLICIT SHAPE — not `(:, :)` — on every optional array below, and
      ! the same on the two `rdb_ocean_dyn` frames that FORWARD them
      ! (`vmix_apply_in_stage`, `visc_rem_precompute`).  Load-bearing, not
      ! decoration.  Each of these is handed straight on to the
      ! EXPLICIT-SHAPE `(nu, nv)` / `(nx, ny, nz)` dummies of
      ! `diffuse_velocity_columns_impl`, i.e. sequence-associated.  An
      ! ASSUMED-SHAPE actual in that position forces the compiler to decide,
      ! at the call, whether the descriptor is contiguous or has to be
      ! packed into a temporary — and gfortran 15 computes that decision
      ! UNCONDITIONALLY, from descriptor extents/strides it only loaded
      ! under the `present()` guard.  An ABSENT optional therefore leaves
      ! the packing flag uninitialised, and the post-call cleanup
      ! (`if (.not. packed .and. ptr /= NULL) free(ptr)`) branches on it:
      ! 160 valgrind "conditional jump depends on uninitialised value(s)"
      ! from 10 contexts on `test_driver_ocean`, all rooted in this
      ! routine's stack frame.  Benign (the null check still gates the free)
      ! but real, and it is the compiler reading OUR uninitialised stack.
      !
      ! Declaring the shape removes the decision at its source: an
      ! explicit-shape dummy is contiguous by definition, so there is no
      ! flag, no temporary, no free, nothing to read — and, unlike the
      ! CONTIGUOUS attribute (tried and reverted during the cavity stack), it
      ! keeps working on nvfortran.  CONTIGUOUS on an assumed-shape OPTIONAL
      ! makes nvfortran
      ! 26.5 hand the callee an offset-adjusted, NON-NULL address for an
      ! ABSENT actual; the implicit present-or-copyin the compiler emits for
      ! the impl's explicit-shape dummy then loses its null short-circuit
      ! and issues a real `cuMemcpyHtoD` off that garbage pointer
      ! (`__pgi_uacc_dataonb ... 'tau_face(:,:)'` -> SIGSEGV, 44/220 tests
      ! on `-stdpar=gpu -gpu=cc70,mem:separate`).  An explicit-shape dummy
      ! carries no descriptor at all, so an absent actual stays a clean
      ! null and the runtime skips it, exactly as it does today.
      !
      ! The bounds are the ones the body already derives and the impl
      ! already assumes: `multilayer_state_t%init` allocates the u-face
      ! arrays `(nx_total + 1, ny_total, nz_ml)` and the v-face arrays
      ! `(nx_total, ny_total + 1, nz_ml)` off this same `grid`.  Every
      ! actual is a whole allocatable slot component or a whole local array,
      ! never a strided section (the GPU path already requires that: these
      ! arrays are device-mapped whole).  Bit-identical: the packing path
      ! was never taken, only speculated about.
      real(wp), intent(in), optional :: kv_source(grid%nx_total, grid%ny_total, ms%nz_ml + 1)
      real(wp), intent(in), optional :: tau_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in), optional :: tau_v(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in), optional :: lambda_bot_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in), optional :: lambda_bot_v(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in), optional :: lambda_top_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in), optional :: lambda_top_v(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in), optional :: cover_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in), optional :: cover_v(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in), optional :: rho0
         !! Boussinesq reference density for the implicit surface-stress
         !! fold, `dt*(tau/rho0)/h_nz`.  READ ONLY when that fold is active
         !! (`do_stress` below); every other path never reaches it.  There
         !! is deliberately NO defaulted value: the one rho0 of record is
         !! `&ocean_ic_nml rho_0` -> `eos%rho0`, which the production call
         !! site hands over as `surface_stress%rho0`
         !! (`configure_ocean_reference_density`).  Omitting it WITH the
         !! fold active is a wiring bug and fails loud — it used to fall
         !! back to a silent literal 1035, which is how a run configured at
         !! a different rho_0 could fold its wind stress in on the wrong
         !! one.
      real(wp), intent(inout), optional :: visc_rem_u(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
      real(wp), intent(inout), optional :: visc_rem_v(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
         !! Explicit-shape for the reason spelled out on the argument list
         !! above.
      logical, intent(in), optional :: remnant_only
         !! `.true.` = build the matrices and (re)fill `visc_rem_u/v`
         !! WITHOUT solving for / modifying the velocities — the
         !! pre-substep visc_rem refresh (PGF_BUG.md §9; MOM6 computes
         !! `vertvisc_coef` before `btstep`, so its BT weights never lag).
         !! Requires both `visc_rem_u/v` present.  Default `.false.`.
      real(wp), intent(in), optional :: kv_corner_source(grid%nx_total + 1, grid%ny_total + 1, ms%nz_ml + 1)
         !! Corner-staggered interface viscosity source, `(nx+1, ny+1,
         !! nz+1)` — see the corner add-on note above.  Typically the
         !! vertex kappa-shear `kd_corner` carrier (device-resident).
         !! Explicit-shape for the reason spelled out on the argument list
         !! above.
      real(wp), intent(in), optional :: kv_corner_prandtl
         !! Scale applied to `kv_corner_source` (Kv = Pr·Kd; the vertex
         !! kappa-shear `prandtl_turb`).  Default 1.

      integer :: nx, ny, nx_face, ny_uface, nx_vface, ny_face, nz
      logical :: do_stress, do_drag, do_remnant, solve_mom, do_corner, do_top
      real(wp) :: rho0_l, corner_prandtl_l

      nx = grid%nx_total
      ny = grid%ny_total
      nx_face = size(ms%u_face_x_layer, 1)
      ny_uface = size(ms%u_face_x_layer, 2)
      nx_vface = size(ms%v_face_y_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      nz = ms%nz_ml

      ! Only fold the boundary BCs when the slot knob is on AND the caller
      ! supplied the corresponding face field (driver wires both together).
      do_stress = this%implicit_stress .and. present(tau_u) .and. present(tau_v)
      do_drag = this%implicit_drag .and. present(lambda_bot_u) .and. present(lambda_bot_v)
      do_top = this%implicit_top_drag .and. present(lambda_top_u) .and. &
               present(lambda_top_v) .and. present(cover_u) .and. present(cover_v)
      do_remnant = present(visc_rem_u) .and. present(visc_rem_v)
      solve_mom = .true.
      if (present(remnant_only)) solve_mom = .not. remnant_only
      if (.not. solve_mom .and. .not. do_remnant) return   ! nothing to produce
      ! Reference density for the implicit stress fold.  The placeholder is
      ! inert — `rho0_l` only reaches an answer under `do_stress`, and that
      ! branch is guarded below — so it must never be a plausible-looking
      ! seawater density that a caller could silently inherit.
      rho0_l = 1.0_wp
      if (present(rho0)) then
         rho0_l = rho0
      else if (do_stress) then
         error stop "vdiff_apply_momentum: implicit_stress is on and tau_u/tau_v "// &
            "were supplied, but rho0 was not — pass the configured reference "// &
            "density (surface_stress%rho0, from &ocean_ic_nml rho_0)"
      end if
      do_corner = present(kv_corner_source)
      corner_prandtl_l = 1.0_wp
      if (present(kv_corner_prandtl)) corner_prandtl_l = kv_corner_prandtl
      ! Glue without the per-face MOM6 BBL (`bbl_per_face`, set at
      ! configure): the historical constant piston, `kv_bbl =
      ! bbl_piston·hbbl_visc` over `bbl_thick = hbbl_visc`, refilled every
      ! call so a hand-built slot that sets the knobs after `init` sees them.
      if (this%bbl_glue .and. .not. this%bbl_per_face) then
         call fill_bbl_constants(this%kv_bbl_u, this%bbl_thick_u, &
                                 this%bbl_piston*this%hbbl_visc, this%hbbl_visc, nx_face, ny_uface)
         call fill_bbl_constants(this%kv_bbl_v, this%bbl_thick_v, &
                                 this%bbl_piston*this%hbbl_visc, this%hbbl_visc, nx_vface, ny_face)
      end if

      if (present(kv_source)) then
         call diffuse_velocity_columns_impl( &
            nx_face, ny_uface, nz, dt, &
            ms%u_face_x_layer, ms%h_layer, kv_source, ms%wet_mask, &
            .true., nx, ny, this%use_harmonic, &
            this%a_diag_u%data, this%b_diag_u%data, &
            this%c_diag_u%data, this%rhs_u%data, &
            do_stress, do_drag, rho0_l, tau_u, lambda_bot_u, &
            do_top, lambda_top_u, cover_u, &
            solve_mom, do_remnant, visc_rem_u, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%hvel_upwind, this%hvel_harmonic, &
            this%kv_bbl_u, this%bbl_thick_u, this%kv_bbl_bg, &
            do_corner, corner_prandtl_l, kv_corner_source, &
            this%zlevel_faces, ms%k_top_u, ms%k_bot_u)
         call diffuse_velocity_columns_impl( &
            nx_vface, ny_face, nz, dt, &
            ms%v_face_y_layer, ms%h_layer, kv_source, ms%wet_mask, &
            .false., nx, ny, this%use_harmonic, &
            this%a_diag_v%data, this%b_diag_v%data, &
            this%c_diag_v%data, this%rhs_v%data, &
            do_stress, do_drag, rho0_l, tau_v, lambda_bot_v, &
            do_top, lambda_top_v, cover_v, &
            solve_mom, do_remnant, visc_rem_v, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%hvel_upwind, this%hvel_harmonic, &
            this%kv_bbl_v, this%bbl_thick_v, this%kv_bbl_bg, &
            do_corner, corner_prandtl_l, kv_corner_source, &
            this%zlevel_faces, ms%k_top_v, ms%k_bot_v)
      else
         ! Pure-vdiff no-op short-circuit ONLY when there is also no
         ! boundary forcing to fold; stress/drag/remnant must still be
         ! applied even at K_v = 0 (the surface Ekman + bed sink don't
         ! need interior viscosity to act, and an implicit_drag-only,
         ! K_v=0 config would otherwise silently leave visc_rem stale).
         ! A corner viscosity source likewise keeps the solve alive.
         if (this%K_v_momentum <= 0.0_wp .and. .not. do_stress .and. .not. do_drag &
             .and. .not. do_top .and. .not. do_remnant .and. .not. do_corner &
             .and. .not. this%bbl_glue) return
         call fill_kv_scalar_buf(this%kv_scalar_buf%data, &
                                 this%K_v_momentum, nx, ny, nz)
         call diffuse_velocity_columns_impl( &
            nx_face, ny_uface, nz, dt, &
            ms%u_face_x_layer, ms%h_layer, this%kv_scalar_buf%data, ms%wet_mask, &
            .true., nx, ny, this%use_harmonic, &
            this%a_diag_u%data, this%b_diag_u%data, &
            this%c_diag_u%data, this%rhs_u%data, &
            do_stress, do_drag, rho0_l, tau_u, lambda_bot_u, &
            do_top, lambda_top_u, cover_u, &
            solve_mom, do_remnant, visc_rem_u, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%hvel_upwind, this%hvel_harmonic, &
            this%kv_bbl_u, this%bbl_thick_u, this%kv_bbl_bg, &
            do_corner, corner_prandtl_l, kv_corner_source, &
            this%zlevel_faces, ms%k_top_u, ms%k_bot_u)
         call diffuse_velocity_columns_impl( &
            nx_vface, ny_face, nz, dt, &
            ms%v_face_y_layer, ms%h_layer, this%kv_scalar_buf%data, ms%wet_mask, &
            .false., nx, ny, this%use_harmonic, &
            this%a_diag_v%data, this%b_diag_v%data, &
            this%c_diag_v%data, this%rhs_v%data, &
            do_stress, do_drag, rho0_l, tau_v, lambda_bot_v, &
            do_top, lambda_top_v, cover_v, &
            solve_mom, do_remnant, visc_rem_v, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%hvel_upwind, this%hvel_harmonic, &
            this%kv_bbl_v, this%bbl_thick_v, this%kv_bbl_bg, &
            do_corner, corner_prandtl_l, kv_corner_source, &
            this%zlevel_faces, ms%k_top_v, ms%k_bot_v)
      end if
   end subroutine vdiff_apply_momentum

   subroutine vdiff_bbl_configure(this, nx, ny, nz, form, cd, r_linear, hbbl, bg_vel, &
                                  thick_min, rino_cap, rho0, kv_bg)
      !! Latch the MOM6 per-face bottom boundary layer (`bbl_per_face`) from
      !! the bottom-drag configuration and size its workspace.  Called once
      !! at configure, BEFORE `enter_data`, by `rdb_ocean_setup` when
      !! `bbl_glue` is on.  MOM6 `set_visc_init` is the parameter map:
      !!
      !!   * quadratic: `CDRAG = cd`, `LINEAR_DRAG = False`;
      !!   * linear: `LINEAR_DRAG = True` with `CDRAG·DRAG_BG_VEL = r·hbbl`,
      !!     the stress `r·hbbl·u` of roundabout's HBBL-distributed linear
      !!     drag, so `CDRAG = r·hbbl/bg_vel` (needs `bg_vel > 0`);
      !!   * `HBBL = hbbl`, `DRAG_BG_VEL = bg_vel`, `BBL_THICK_MIN =
      !!     thick_min`, `KV = kv_bg`, `RiNo_mix = rino_cap`;
      !!   * `BBL_USE_EOS = True` (mandatory in ALE mode),
      !!     `CORRECT_BBL_BOUNDS = False`, `BBL_USE_TIDAL_BG = False`,
      !!     `CHANNEL_DRAG = False`, `BODY_FORCE_DRAG = False`.
      !!
      !! A form with no drag (`cd <= 0` quadratic, `r <= 0` / `bg_vel <= 0`
      !! linear) or `hbbl <= 0` leaves the glue OFF (`bbl_glue = .false.`)
      !! and the caller's drag path untouched — there is no MOM6
      !! bottom-drag law to port without a drag, and MOM6 has no HBBL-less
      !! one.
      class(ocean_vdiff_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, nz
         !! Cell-centred extents (ghosts included) of the workspace.
      integer, intent(in) :: form
         !! `BBL_FORM_LINEAR` or `BBL_FORM_QUADRATIC`.
      real(wp), intent(in) :: cd
         !! Quadratic drag coefficient (dimensionless).
      real(wp), intent(in) :: r_linear
         !! Linear drag rate (1/s) of the HBBL-distributed linear drag.
      real(wp), intent(in) :: hbbl
         !! MOM6 `HBBL` (m).
      real(wp), intent(in) :: bg_vel
         !! MOM6 `DRAG_BG_VEL` (m/s).
      real(wp), intent(in) :: thick_min
         !! MOM6 `BBL_THICK_MIN` (m).
      logical, intent(in) :: rino_cap
         !! MOM6 `RiNo_mix` (kappa-shear on).
      real(wp), intent(in) :: rho0
         !! Boussinesq reference density (kg/m^3).
      real(wp), intent(in) :: kv_bg
         !! MOM6 `KV`, the background viscosity (m^2/s).

      this%bbl_per_face = .false.
      if (.not. this%bbl_glue) return
      this%bbl_form = form
      this%bbl_hbbl = hbbl
      this%bbl_bg_vel = bg_vel
      this%bbl_thick_min = thick_min
      this%bbl_rino_cap = rino_cap
      this%bbl_rho0 = rho0
      this%kv_bbl_bg = kv_bg
      if (form == BBL_FORM_LINEAR) then
         this%bbl_cd = 0.0_wp
         if (bg_vel > 0.0_wp) this%bbl_cd = r_linear*hbbl/bg_vel
      else
         this%bbl_cd = cd
      end if
      if (this%bbl_cd <= 0.0_wp .or. hbbl <= 0.0_wp) then
         this%bbl_glue = .false.
         return
      end if
      this%bbl_per_face = .true.
      if (allocated(this%bbl_conc_t)) deallocate (this%bbl_conc_t)
      if (allocated(this%bbl_conc_s)) deallocate (this%bbl_conc_s)
      allocate (this%bbl_conc_t(nx, ny, nz), source=0.0_wp)
      allocate (this%bbl_conc_s(nx, ny, nz), source=0.0_wp)
   end subroutine vdiff_bbl_configure

   subroutine vdiff_set_viscous_bbl(grid, this, ms, eos, f_corner)
      !! MOM6 `set_viscous_BBL` (the
      !! `BOTTOMDRAGLAW` branch): the per-face bottom-boundary-layer
      !! viscosity `kv_bbl_u/v` and thickness `bbl_thick_u/v` that the
      !! vertical-friction glue reads.  Called ONCE per outer step, before
      !! the stage loop, from the state at the start of the step — MOM6
      !! calls it once per step, from `step_MOM_dynamics`, before the
      !! predictor.  No-op unless `bbl_glue .and. bbl_per_face`.
      !!
      !! Per face (bed up; MOM6 counts from the surface, so its `k = nz`
      !! is our `k = 1`):
      !!
      !!   1. `h_at_vel`: harmonic mean of the two cells when the flow runs
      !!      thin -> thick (`u·Δh >= 0`), else arithmetic;
      !!      T and S at the face are plain averages of the two cells'
      !!      CONCENTRATIONS (I1′ column rule).
      !!   2. `u_bbl`: the `h_at_vel`-weighted mean, over the bottom `HBBL`,
      !!      of `sqrt(u² + v_at_u² + DRAG_BG_VEL²)` (`v_at_u` the masked
      !!      thickness-weighted average of the four transverse faces,
      !!      `set_v_at_u`); `u* = sqrt(CDRAG)·u_bbl`, or
      !!      `sqrt(CDRAG)·DRAG_BG_VEL` for the linear law or an empty
      !!      average.
      !!   3. `h_N`: the stratification-limited thickness, integrating the
      !!      density jump up from the bed until `h·Δρ` reaches
      !!      `400·ρ₀·u*²/g` (Killworth & Edwards 1999 eq. 2.22), with
      !!      `∂ρ/∂T`, `∂ρ/∂S` at the BBL-mean T/S and the
      !!      bottom pressure `ρ₀·g·Σh` (no surface-pressure term).
      !!   4. `bbl_thick = h_N/(1/2 + sqrt(1/4 + (2f·h_N/u*)²))` (rotation,
      !!      Killworth & Edwards 1999 eq. 2.20), floored at `BBL_THICK_MIN`,
      !!      capped at `HBBL/2` under kappa-shear.
      !!   5. `kv_bbl = sqrt(CDRAG)·u*·bbl_thick` — the viscosity whose
      !!      stress over `bbl_thick` is `CDRAG·u_bbl²`.
      !!
      !! A dry face gets `kv_bbl = 0`, `bbl_thick = HBBL`.  Open-boundary
      !! faces are NOT given MOM6's zero-gradient projection: they read
      !! their two cells like any interior face.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vdiff_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(eos_t), intent(in) :: eos
         !! The run's equation of state, for `∂ρ/∂T` and `∂ρ/∂S`.
      real(wp), intent(in) :: f_corner(grid%nx_total + 1, grid%ny_total + 1)
         !! Coriolis parameter at cell SW corners (1/s), `cor%f_corner`.

      integer :: nx, ny, nz
      logical :: use_eos

      if (.not. (this%bbl_glue .and. this%bbl_per_face)) return
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      use_eos = ms%idx_temperature > 0 .and. ms%idx_salinity > 0
      if (use_eos) then
         call bbl_column_conc_impl(nx, ny, nz, ms%h_layer, &
                                   ms%tracers(ms%idx_temperature)%c, this%bbl_conc_t)
         call bbl_column_conc_impl(nx, ny, nz, ms%h_layer, &
                                   ms%tracers(ms%idx_salinity)%c, this%bbl_conc_s)
      end if
      call bbl_faces_impl(nx + 1, ny, nx, ny, nz, .true., &
                          ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, ms%wet_mask, &
                          this%bbl_conc_t, this%bbl_conc_s, size(this%bbl_conc_t, 1), &
                          size(this%bbl_conc_t, 2), size(this%bbl_conc_t, 3), &
                          f_corner, use_eos, eos, this%bbl_form, this%bbl_cd, this%bbl_hbbl, &
                          this%bbl_bg_vel, this%bbl_thick_min, this%bbl_rino_cap, &
                          this%bbl_rho0, this%kv_bbl_u, this%bbl_thick_u)
      call bbl_faces_impl(nx, ny + 1, nx, ny, nz, .false., &
                          ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, ms%wet_mask, &
                          this%bbl_conc_t, this%bbl_conc_s, size(this%bbl_conc_t, 1), &
                          size(this%bbl_conc_t, 2), size(this%bbl_conc_t, 3), &
                          f_corner, use_eos, eos, this%bbl_form, this%bbl_cd, this%bbl_hbbl, &
                          this%bbl_bg_vel, this%bbl_thick_min, this%bbl_rino_cap, &
                          this%bbl_rho0, this%kv_bbl_v, this%bbl_thick_v)
   end subroutine vdiff_set_viscous_bbl

   pure subroutine bbl_column_conc_impl(nx, ny, nz, h, htr, conc)
      !! Cell concentrations of one tracer through THE vanished-layer column
      !! rule (`rdb_vl_column_conc`): `hTr/h` on a live layer, the donor's
      !! concentration on a filler, never a quotient by a near-zero or
      !! negative thickness.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: htr(nx, ny, nz)
      real(wp), intent(inout) :: conc(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: hc(NZ_STACK_MAX), qc(NZ_STACK_MAX), cc(NZ_STACK_MAX)

      do concurrent(j=1:ny, i=1:nx) local(hc, qc, cc, k)
         do k = 1, nz
            hc(k) = h(i, j, k)
            qc(k) = htr(i, j, k)
         end do
         call rdb_vl_column_conc(nz, hc, qc, cc)
         do k = 1, nz
            conc(i, j, k) = cc(k)
         end do
      end do
   end subroutine bbl_column_conc_impl

   pure subroutine bbl_faces_impl(nu, nv, nx, ny, nz, x_face, u_x, v_y, h, wet, &
                                  conc_t, conc_s, ncx, ncy, ncz, f_corner, use_eos, eos, &
                                  form, cd, hbbl, bg_vel, thick_min, rino_cap, rho0, &
                                  kv_bbl, bbl_thick)
      !! Per-face body of `vdiff_set_viscous_bbl` (see there).  `x_face`
      !! selects the u-faces `(nx+1, ny)` (cells `i-1`, `i`) or the v-faces
      !! `(nx, ny+1)` (cells `j-1`, `j`).
      integer, intent(in) :: nu, nv, nx, ny, nz, ncx, ncy, ncz
      logical, intent(in) :: x_face
      real(wp), intent(in) :: u_x(nx + 1, ny, nz)
      real(wp), intent(in) :: v_y(nx, ny + 1, nz)
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: conc_t(ncx, ncy, ncz)
      real(wp), intent(in) :: conc_s(ncx, ncy, ncz)
         !! Cell T/S concentrations; read only when `use_eos` (then
         !! `(nx, ny, nz)`; a `(1,1,1)` placeholder otherwise).
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1)
      logical, intent(in) :: use_eos
      type(eos_t), intent(in) :: eos
      integer, intent(in) :: form
      real(wp), intent(in) :: cd, hbbl, bg_vel, thick_min, rho0
      logical, intent(in) :: rino_cap
      real(wp), intent(inout) :: kv_bbl(nu, nv)
      real(wp), intent(inout) :: bbl_thick(nu, nv)

      real(wp), parameter :: ANGSTROM = 1.0e-10_wp
         !! MOM6 `ANGSTROM` default (m): the thinnest layer counted in the
         !! HBBL average (`hweight < 1.5·Angstrom` is skipped, line 666).
      real(wp), parameter :: EPS_NEG = 1.0e-30_wp
         !! MOM6 `h_neglect` / `dz_neglect` (Boussinesq, m).
      integer :: i, j, k, il, ir, jl, jr, im, ip, jm, jp
      real(wp) :: hav(NZ_STACK_MAX), tv(NZ_STACK_MAX), sv(NZ_STACK_MAX)
      real(wp) :: hl, hr, un, press, cd_sqrt, u2_bg, ustar, ustarsq
      real(wp) :: htot_vel, hwtot, hutot, thtot, shtot, hweight, vt
      real(wp) :: wa, wb, wc, wd, wsum, t_eos, s_eos
      real(wp) :: drdt, drds, htot, oldfn, dfn, dh, c2f, root, thick

      cd_sqrt = sqrt(cd)
      u2_bg = bg_vel*bg_vel
      do concurrent(j=1:nv, i=1:nu) &
         local(il, ir, jl, jr, im, ip, jm, jp, k, hav, tv, sv, hl, hr, un, press, &
               ustar, ustarsq, htot_vel, hwtot, hutot, thtot, shtot, hweight, vt, &
               wa, wb, wc, wd, wsum, t_eos, s_eos, drdt, drds, &
               htot, oldfn, dfn, dh, c2f, root, thick)
         if (x_face) then
            il = max(1, i - 1)
            ir = min(nx, i)
            jl = j
            jr = j
         else
            il = i
            ir = i
            jl = max(1, j - 1)
            jr = min(ny, j)
         end if
         kv_bbl(i, j) = 0.0_wp
         bbl_thick(i, j) = hbbl
         if (wet(il, jl)*wet(ir, jr) > 0.0_wp) then
            ! ---- 1. face thickness, T, S, bottom pressure (lines 491-516, 791-803)
            press = 0.0_wp
            do k = 1, nz
               hl = h(il, jl, k)
               hr = h(ir, jr, k)
               if (x_face) then
                  un = u_x(i, j, k)
               else
                  un = v_y(i, j, k)
               end if
               if (un*(hr - hl) >= 0.0_wp) then
                  hav(k) = 2.0_wp*hl*hr/(hl + hr + EPS_NEG)
               else
                  hav(k) = 0.5_wp*(hl + hr)
               end if
               press = press + rho0*GRAVITY*0.5_wp*(hl + hr)
               tv(k) = 0.0_wp
               sv(k) = 0.0_wp
               if (use_eos) then
                  tv(k) = 0.5_wp*(conc_t(il, jl, k) + conc_t(ir, jr, k))
                  sv(k) = 0.5_wp*(conc_s(il, jl, k) + conc_s(ir, jr, k))
               end if
            end do

            ! ---- 2. u_bbl and u* over the bottom HBBL ----
            htot_vel = 0.0_wp
            hwtot = 0.0_wp
            hutot = 0.0_wp
            thtot = 0.0_wp
            shtot = 0.0_wp
            do k = 1, nz
               if (htot_vel >= hbbl) exit
               hweight = min(hbbl - htot_vel, hav(k))
               if (hweight < 1.5_wp*ANGSTROM + EPS_NEG) cycle
               htot_vel = htot_vel + hav(k)
               hwtot = hwtot + hweight
               if (form == BBL_FORM_QUADRATIC) then
                  ! `set_v_at_u` / `set_u_at_v`: masked, thickness-weighted
                  ! mean of the four transverse faces at this layer.
                  if (x_face) then
                     un = u_x(i, j, k)
                     jm = max(1, j - 1)
                     jp = min(ny, j + 1)
                     wa = max(0.0_wp, h(il, jm, k) + h(il, j, k))*wet(il, jm)*wet(il, j)
                     wb = max(0.0_wp, h(ir, jm, k) + h(ir, j, k))*wet(ir, jm)*wet(ir, j)
                     wc = max(0.0_wp, h(il, j, k) + h(il, jp, k))*wet(il, j)*wet(il, jp)
                     wd = max(0.0_wp, h(ir, j, k) + h(ir, jp, k))*wet(ir, j)*wet(ir, jp)
                     wsum = (wa + wd) + (wb + wc)
                     vt = 0.0_wp
                     if (wsum > 0.0_wp) then
                        vt = ((wa*v_y(il, j, k) + wd*v_y(ir, j + 1, k)) + &
                              (wb*v_y(ir, j, k) + wc*v_y(il, j + 1, k)))/wsum
                     end if
                  else
                     un = v_y(i, j, k)
                     im = max(1, i - 1)
                     ip = min(nx, i + 1)
                     wa = max(0.0_wp, h(im, jl, k) + h(i, jl, k))*wet(im, jl)*wet(i, jl)
                     wb = max(0.0_wp, h(im, jr, k) + h(i, jr, k))*wet(im, jr)*wet(i, jr)
                     wc = max(0.0_wp, h(i, jl, k) + h(ip, jl, k))*wet(i, jl)*wet(ip, jl)
                     wd = max(0.0_wp, h(i, jr, k) + h(ip, jr, k))*wet(i, jr)*wet(ip, jr)
                     wsum = (wa + wd) + (wb + wc)
                     vt = 0.0_wp
                     if (wsum > 0.0_wp) then
                        vt = ((wa*u_x(i, jl, k) + wd*u_x(i + 1, jr, k)) + &
                              (wb*u_x(i, jr, k) + wc*u_x(i + 1, jl, k)))/wsum
                     end if
                  end if
                  hutot = hutot + hweight*sqrt(un*un + vt*vt + u2_bg)
               end if
               thtot = thtot + hweight*tv(k)
               shtot = shtot + hweight*sv(k)
            end do
            if (hwtot <= 0.0_wp .or. form == BBL_FORM_LINEAR) then
               ustar = cd_sqrt*bg_vel
            else
               ustar = cd_sqrt*hutot/hwtot
            end if

            ! ---- 3. stratification-limited thickness h_N ----
            drdt = 0.0_wp
            drds = 0.0_wp
            if (use_eos) then
               t_eos = 0.0_wp
               s_eos = 0.0_wp
               if (hwtot > 0.0_wp) then
                  t_eos = thtot/hwtot
                  s_eos = shtot/hwtot
               end if
               ! MOM6 `calculate_density_derivs` (closed-form ∂ρ/∂T, ∂ρ/∂S).
               call eos_density_derivs(eos, t_eos, s_eos, press, drdt, drds)
            end if
            ! The 400 is Ci² (KW99 eq. 2.22); Boussinesq `Rho0x400_G`.
            ustarsq = 400.0_wp*rho0/GRAVITY*ustar*ustar
            htot = 0.0_wp
            thtot = 0.0_wp
            shtot = 0.0_wp
            oldfn = 0.0_wp
            do k = 1, nz - 1
               if (hav(k) <= 0.0_wp) cycle
               ! Δρ·h of the BBL with everything below homogenised.
               oldfn = drdt*(thtot - tv(k)*htot) + drds*(shtot - sv(k)*htot)
               if (oldfn >= ustarsq) exit
               dfn = (drdt*(tv(k) - tv(k + 1)) + drds*(sv(k) - sv(k + 1)))*(hav(k) + htot)
               if (oldfn + dfn <= ustarsq) then
                  dh = hav(k)
               else
                  dh = hav(k)*sqrt((ustarsq - oldfn)/dfn)
               end if
               htot = htot + dh
               thtot = thtot + tv(k)*dh
               shtot = shtot + sv(k)*dh
            end do
            ! The top layer might be part of the BBL.
            if (oldfn < ustarsq .and. hav(nz) > 0.0_wp) then
               if (drdt*(thtot - tv(nz)*htot) + drds*(shtot - sv(nz)*htot) < ustarsq) then
                  htot = htot + hav(nz)
               end if
            end if

            ! ---- 4. rotation (Killworth & Edwards 1999 eq. 2.20) + caps ----
            if (x_face) then
               c2f = f_corner(i, j) + f_corner(i, j + 1)
            else
               c2f = f_corner(i, j) + f_corner(i + 1, j)
            end if
            if (cd*u2_bg <= 0.0_wp) then
               root = sqrt(0.25_wp*ustar*ustar + (htot*c2f)**2)
               if (htot*ustar <= (thick_min + EPS_NEG)*(0.5_wp*ustar + root)) then
                  thick = thick_min
               else
                  thick = (htot*ustar)/(0.5_wp*ustar + root)
               end if
            else
               thick = htot/(0.5_wp + sqrt(0.25_wp + htot*htot*c2f*c2f/(ustar*ustar)))
               if (thick < thick_min) thick = thick_min
            end if
            if (rino_cap .and. thick > 0.5_wp*hbbl) thick = 0.5_wp*hbbl

            ! ---- 5. the viscosity that carries CDRAG·u_bbl² ----
            kv_bbl(i, j) = cd_sqrt*ustar*thick
            bbl_thick(i, j) = thick
         end if
      end do
   end subroutine bbl_faces_impl

   pure subroutine fill_bbl_constants(kv_bbl, bbl_thick, kv_val, thick_val, nu, nv)
      !! The historical constant-piston glue (`bbl_per_face = .false.`).
      integer, intent(in) :: nu, nv
      real(wp), intent(inout) :: kv_bbl(nu, nv)
      real(wp), intent(inout) :: bbl_thick(nu, nv)
      real(wp), intent(in) :: kv_val, thick_val
      integer :: i, j

      do concurrent(j=1:nv, i=1:nu)
         kv_bbl(i, j) = kv_val
         bbl_thick(i, j) = thick_val
      end do
   end subroutine fill_bbl_constants

   subroutine vdiff_apply_tracers(grid, this, ms, dt, kt_source, ks_source)
      !! Backward-Euler vertical diffusivity on every registered
      !! tracer.  Each tracer is converted to `T = hTr/h`, the
      !! tridiagonal solve runs, then `hTr = T*h` is reconstituted.
      !! `h_layer` is untouched.  Per-tracer
      !! `do_vertical_diffusion` flag gates participation.
      !!
      !! Diffusivity dispatch — three cases:
      !!   1. `ks_source` absent -> single-source legacy path: one
      !!      factorize from `kt_source` (3D, interface-located) if
      !!      present, else the scalar `K_v_tracer` fallback; ALL
      !!      tracers use it.  Every existing caller lands here ⇒
      !!      bit-identical to the pre-PR-20 behaviour.
      !!   2. `ks_source` present, `kt_source` absent -> programming
      !!      error, fail loud (a caller that supplies salt but not
      !!      heat has a bug — MOM6 guards the same pairing).
      !!   3. Both present -> two passes, same buffers reused: pass 1
      !!      factorizes `kt_source` and applies it to temperature
      !!      only; pass 2 factorizes `ks_source` (overwriting the
      !!      same `a_diag_t`/`b_diag_t`/`c_diag_t` buffers — no new
      !!      scratch) and applies it to salinity AND every other
      !!      registered passive tracer.  This is MOM6's `Kd_salt`
      !!      convention — "the diapycnal diffusivity of salt AND
      !!      PASSIVE TRACERS" — so passives follow salt,
      !!      not heat.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vdiff_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: kt_source(:, :, :)
         !! Temperature diffusivity, (nx, ny, nz+1), interface-located.
      real(wp), intent(in), optional :: ks_source(:, :, :)
         !! Salinity + passive-tracer diffusivity, same shape.  MUST
         !! be accompanied by `kt_source` (case 2 fails loud).

      integer :: it, nx, ny, nz
      logical :: use_source, two_pass

      if (present(ks_source) .and. .not. present(kt_source)) then
         call logger%error("vdiff_apply_tracers: ks_source requires kt_source")
         error stop "vdiff_apply_tracers: ks_source requires kt_source"
      end if

      use_source = present(kt_source)
      two_pass = present(ks_source)
      if (.not. use_source .and. this%K_v_tracer <= 0.0_wp) return
      if (.not. allocated(ms%tracers)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (.not. use_source) then
         call fill_kv_scalar_buf(this%kv_scalar_buf%data, &
                                 this%K_v_tracer, nx, ny, nz)
      end if

      ! The tridiagonal is tracer-independent (depends only on kv/h/dt),
      ! so build + Thomas-factorize it ONCE per diffusivity field; every
      ! tracer sharing that field then reuses the factored coefficients.
      ! Saves rebuilding the matrix (the face_thick-heavy part) per
      ! tracer — the win grows with the tracer count.  Case 3 pays two
      ! factorizations (heat, then salt+passives) instead of one; case 1
      ! is the literal pre-PR-20 single-factorize path, unchanged.
      if (use_source) then
         call build_factorize_tracer_matrix(nx, ny, nz, this%use_harmonic, dt, &
                                            kt_source, ms%h_layer, &
                                            this%a_diag_t%data, this%b_diag_t%data, &
                                            this%c_diag_t%data)
      else
         call build_factorize_tracer_matrix(nx, ny, nz, this%use_harmonic, dt, &
                                            this%kv_scalar_buf%data, ms%h_layer, &
                                            this%a_diag_t%data, this%b_diag_t%data, &
                                            this%c_diag_t%data)
      end if

      do it = 1, size(ms%tracers)
         if (.not. ms%tracers(it)%do_vertical_diffusion) cycle
         if (two_pass .and. it /= ms%idx_temperature) cycle
         select case (ms%tracers(it)%budget_id)
         case (TRACER_BUDGET_HEAT)
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%c, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%heat_budget_vdiff)
         case (TRACER_BUDGET_SALT)
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%c, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%salt_budget_vdiff)
         case default
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%c, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data)
         end select
      end do

      if (.not. two_pass) return

      ! ---- Pass 2: ks_source -> salinity + every passive tracer ----
      call build_factorize_tracer_matrix(nx, ny, nz, this%use_harmonic, dt, &
                                         ks_source, ms%h_layer, &
                                         this%a_diag_t%data, this%b_diag_t%data, &
                                         this%c_diag_t%data)

      do it = 1, size(ms%tracers)
         if (.not. ms%tracers(it)%do_vertical_diffusion) cycle
         if (it == ms%idx_temperature) cycle
         if (it == ms%idx_salinity) then
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%c, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%salt_budget_vdiff)
         else
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%c, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data)
         end if
      end do
   end subroutine vdiff_apply_tracers

   pure subroutine fill_kv_scalar_buf(kv_buf, kappa, nx, ny, nz)
      !! Broadcast the scalar viscosity `kappa` into the
      !! interface-located workspace.  Boundary interfaces (bed at
      !! k=1, surface at k=nz+1) are forced to zero to match the
      !! closed BCs the column solve assumes; the original scalar
      !! kernel hard-coded those BCs via `α_1 = 0` and `β_nz = 0`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: kappa
      real(wp), intent(out) :: kv_buf(nx, ny, nz + 1)
      integer :: i, j, k

      do concurrent(k=1:nz + 1, j=1:ny, i=1:nx)
         if (k == 1 .or. k == nz + 1) then
            kv_buf(i, j, k) = 0.0_wp
         else
            kv_buf(i, j, k) = kappa
         end if
      end do
   end subroutine fill_kv_scalar_buf

   pure subroutine build_factorize_tracer_matrix(nx, ny, nz, use_harmonic, dt, kv_centre, &
                                                 h_layer, a_diag, b_diag, c_diag)
      !! Build the backward-Euler tridiagonal per cell column and run the
      !! Thomas forward factorization — the tracer-INDEPENDENT half of the
      !! vertical-diffusion solve (depends only on kv/h/dt).  Run once per
      !! stage; every registered tracer then reuses the factored
      !! coefficients via `apply_factored_tracer`.
      !!
      !! `kv_centre(:, :, k)` is the diffusivity at the BOTTOM interface
      !! of layer `k`; `kv_centre(:, :, nz+1)` is the surface interface.
      !! On exit:
      !!   a_diag = sub-diagonal (unchanged), for the per-tracer RHS sweep
      !!   b_diag = the Thomas pivots (b(1) = raw diagonal, b(k>1) = denom)
      !!   c_diag = super-diagonal already divided by its pivot
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: use_harmonic
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: kv_centre(nx, ny, nz + 1)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: a_diag(nx, ny, nz)
      real(wp), intent(inout) :: b_diag(nx, ny, nz)
      real(wp), intent(inout) :: c_diag(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: dz_face, alpha, beta, denom
      real(wp) :: hc, hm, hp

      ! Vanishing-layer handling (load-bearing for conservation under
      ! Fox-Kemper × windowed tracer advect): a thin z* surface layer can be
      ! driven to h ≤ 0 on an intermediate RK2 stage by the combined FK +
      ! resolved transport (Lagrangian, before the ALE remap).  The raw
      ! 1/h in α/β would then go negative/Inf and corrupt the WHOLE column's
      ! solve, dropping tracer mass.  We floor the per-layer thickness to
      ! H_VANISHED in the denominators (the SAME h̃ apply_factored_tracer
      ! uses) and ZERO the diffusive flux at any interface touching a
      ! vanishing layer (h ≤ H_VANISHED) — so a collapsed layer is decoupled
      ! (identity row) and its frozen tracer mass is preserved exactly,
      ! while its thick neighbours conserve among themselves.  For
      ! h ≫ H_VANISHED (sigma / double_gyre) h̃ = h and every interface is
      ! active ⇒ bit-identical.
      do concurrent(j=1:ny, i=1:nx) local(dz_face, alpha, beta, denom, k, hc, hm, hp)
         ! k = 1: bed BC (no flux below).  α uses the interface above
         ! layer 1 (= kv_centre(:, :, 2)).
         hc = max(h_layer(i, j, 1), H_VANISHED)
         ! SINGLE-LAYER COLUMN (nz = 1): the bed row IS the surface row and
         ! there is no interior interface, so alpha is identically zero.
         ! Without the gate this reads `h_layer(i, j, 2)` — past the end of
         ! a (nx, ny, 1) array — and the k = nz block below then reads
         ! `h_layer(i, j, 0)` and overwrites this row.  Same defect (and
         ! same shape of fix) as the velocity tridiagonal below.
         ! Loop-invariant gate INSIDE the single DC ⇒ nz >= 2 bit-identical.
         alpha = 0.0_wp
         if (nz > 1) then
            hp = max(h_layer(i, j, 2), H_VANISHED)
            dz_face = face_thick(hc, hp, use_harmonic)
            alpha = dt*kv_centre(i, j, 2)/(hc*dz_face)
            if (h_layer(i, j, 1) <= H_VANISHED .or. h_layer(i, j, 2) <= H_VANISHED) alpha = 0.0_wp
         end if
         a_diag(i, j, 1) = 0.0_wp
         c_diag(i, j, 1) = -alpha
         b_diag(i, j, 1) = 1.0_wp + alpha

         ! k = 2..nz-1: interior.  β = kv_centre(k), α = kv_centre(k+1).
         do k = 2, nz - 1
            hm = max(h_layer(i, j, k - 1), H_VANISHED)
            hc = max(h_layer(i, j, k), H_VANISHED)
            hp = max(h_layer(i, j, k + 1), H_VANISHED)
            beta = dt*kv_centre(i, j, k)/(hc*face_thick(hm, hc, use_harmonic))
            alpha = dt*kv_centre(i, j, k + 1)/(hc*face_thick(hc, hp, use_harmonic))
            if (h_layer(i, j, k) <= H_VANISHED .or. h_layer(i, j, k - 1) <= H_VANISHED) beta = 0.0_wp
            if (h_layer(i, j, k) <= H_VANISHED .or. h_layer(i, j, k + 1) <= H_VANISHED) alpha = 0.0_wp
            a_diag(i, j, k) = -beta
            c_diag(i, j, k) = -alpha
            b_diag(i, j, k) = 1.0_wp + alpha + beta
         end do

         ! k = nz: surface BC (no flux above).  β uses kv_centre(nz).
         ! SINGLE-LAYER COLUMN (nz = 1): already built as the bed row above;
         ! skip (see the k = 1 gate).  Otherwise this reads h_layer(:, :, 0).
         if (nz > 1) then
            hm = max(h_layer(i, j, nz - 1), H_VANISHED)
            hc = max(h_layer(i, j, nz), H_VANISHED)
            dz_face = face_thick(hm, hc, use_harmonic)
            beta = dt*kv_centre(i, j, nz)/(hc*dz_face)
            if (h_layer(i, j, nz) <= H_VANISHED .or. h_layer(i, j, nz - 1) <= H_VANISHED) beta = 0.0_wp
            a_diag(i, j, nz) = -beta
            c_diag(i, j, nz) = 0.0_wp
            b_diag(i, j, nz) = 1.0_wp + beta
         end if

         ! ---- Thomas forward factorization (matrix only) ----
         ! Store the pivot `denom` back into b_diag and the eliminated
         ! super-diagonal `c/denom` into c_diag.  b_diag(1) keeps the raw
         ! diagonal (its own pivot); a_diag stays the raw sub-diagonal.
         c_diag(i, j, 1) = c_diag(i, j, 1)/b_diag(i, j, 1)
         do k = 2, nz
            denom = b_diag(i, j, k) - a_diag(i, j, k)*c_diag(i, j, k - 1)
            c_diag(i, j, k) = c_diag(i, j, k)/denom
            b_diag(i, j, k) = denom
         end do
      end do
   end subroutine build_factorize_tracer_matrix

   pure subroutine apply_factored_tracer(nx, ny, nz, h_layer, c, &
                                         a_diag, b_diag, c_diag, rhs, budget)
      !! Apply the pre-factored tridiagonal (from
      !! `build_factorize_tracer_matrix`) to one tracer's CONCENTRATION: run
      !! the Thomas RHS sweep + back-substitution against the stored pivots
      !! and store `c_new`.  The matrix is weighted by the floored thickness
      !! `h̃ = max(h, H_VANISHED)`, so Σ h̃·c is conserved and a layer at or
      !! below `H_VANISHED` is an identity row that keeps its `c` exactly;
      !! `h_layer` and the factored `a/b/c` are read-only (shared across
      !! every tracer).
      !!
      !! When `budget` is present, accumulate the per-cell increment
      !! `h·c_new − h·c_old` into the contributor slot before overwriting
      !! (two standalone products so the budget telescopes with the console
      !! total).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: c(nx, ny, nz)
      real(wp), intent(in) :: a_diag(nx, ny, nz)
      real(wp), intent(in) :: b_diag(nx, ny, nz)
      real(wp), intent(in) :: c_diag(nx, ny, nz)
      real(wp), intent(inout) :: rhs(nx, ny, nz)
      real(wp), intent(inout), optional :: budget(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: p_new, p_old

      do concurrent(j=1:ny, i=1:nx) local(k, p_new, p_old)
         ! RHS = c.  A vanishing / collapsed layer (h ≤ H_VANISHED) is an
         ! identity row of the matrix (built with the same floored h̃), so its
         ! c is preserved EXACTLY.
         do k = 1, nz
            rhs(i, j, k) = c(i, j, k)
         end do

         ! ---- Thomas RHS forward sweep against the stored pivots ----
         rhs(i, j, 1) = rhs(i, j, 1)/b_diag(i, j, 1)
         do k = 2, nz
            rhs(i, j, k) = (rhs(i, j, k) - a_diag(i, j, k)*rhs(i, j, k - 1))/b_diag(i, j, k)
         end do

         ! ---- Back-substitution (rhs now holds T_new) ----
         do k = nz - 1, 1, -1
            rhs(i, j, k) = rhs(i, j, k) - c_diag(i, j, k)*rhs(i, j, k + 1)
         end do

         ! ---- Store c_new (budget = h·c_new − h·c_old) ----
         ! Budget write gated INSIDE the single do concurrent (splitting
         ! present() into two loops makes NVHPC emit a far slower kernel
         ! for one branch — see the remap fix).
         if (present(budget)) then
            do k = 1, nz
               p_new = h_layer(i, j, k)*rhs(i, j, k)
               p_old = h_layer(i, j, k)*c(i, j, k)
               budget(i, j, k) = budget(i, j, k) + (p_new - p_old)
               c(i, j, k) = rhs(i, j, k)
            end do
         else
            do k = 1, nz
               c(i, j, k) = rhs(i, j, k)
            end do
         end if
      end do
   end subroutine apply_factored_tracer

   pure subroutine diffuse_velocity_columns_impl(nu, nv, nz, dt, &
                                                 u_face, h_layer, kv_centre, &
                                                 wet_cell, &
                                                 x_face, nx_cells, ny_cells, &
                                                 use_harmonic, &
                                                 a_diag, b_diag, c_diag, rhs, &
                                                 do_stress, do_drag, rho0, &
                                                 tau_face, lambda_bot, &
                                                 do_top, lambda_top, cover_face, &
                                                 solve_momentum, &
                                                 do_remnant, visc_rem_out, &
                                                 hvel_mom6, hbbl_visc, &
                                                 bbl_glue, hvel_upwind, hvel_harmonic, &
                                                 kv_bbl_face, bbl_thick_face, kv_bbl_bg, &
                                                 do_corner, kv_prandtl, kv_corner, &
                                                 zlevel_faces, k_top_face, k_bot_face)
      !! Build + solve the tridiagonal system per face column.
      !! `u_face` is either `u_face_x_layer` (`x_face = .true.`,
      !! shape (nx+1, ny)) or `v_face_y_layer` (`x_face = .false.`,
      !! shape (nx, ny+1)).  `h_layer` is cell-centred (nx, ny, nz)
      !! and we average across the face direction to get the face
      !! thickness.  `kv_centre` is the same cell-centred diffusivity
      !! field used by the tracer kernel — we average it across the
      !! face to get a face-located value at each interface.
      !!
      !! For wall faces (i = 1 / i = nx_cells+1 for x_face;
      !! j = 1 / j = ny_cells+1 for y-face), there's no neighbour
      !! cell to average with — we fall back to the single available
      !! cell's thickness + diffusivity.  Mass through walls is
      !! forced to zero by the continuity kernel anyway, so the
      !! wall-face viscosity is only there to keep the system
      !! non-singular.
      integer, intent(in) :: nu, nv, nz
      integer, intent(in) :: nx_cells, ny_cells
      real(wp), intent(in) :: dt
      real(wp), intent(inout) :: u_face(nu, nv, nz)
      real(wp), intent(in) :: h_layer(nx_cells, ny_cells, nz)
      real(wp), intent(in) :: kv_centre(nx_cells, ny_cells, nz + 1)
      real(wp), intent(in) :: wet_cell(nx_cells, ny_cells)
         !! Cell-centred wet/dry land mask (1 wet, 0 land).  The wind-stress
         !! fold multiplies `tau_face` by the face mask `min` of the two
         !! bounding cells — matching the explicit surface-stress kernel —
         !! so no stress is injected at a no-normal-flow land face.  All-wet
         !! (mask ≡ 1) ⇒ bit-identical.  Only read when `do_stress`.
      logical, intent(in) :: x_face
      logical, intent(in) :: use_harmonic
      real(wp), intent(inout) :: a_diag(nu, nv, nz)
         !! Unused since the fused sweep (no row is stored before it is
         !! eliminated); kept so the call sites are unchanged.
      real(wp), intent(inout) :: b_diag(nu, nv, nz)
         !! Unused, as `a_diag`.
      real(wp), intent(inout) :: c_diag(nu, nv, nz)
         !! Scratch: ends holding the Thomas c'.
      real(wp), intent(inout) :: rhs(nu, nv, nz)
         !! Scratch: ends holding the Thomas d' (when `solve_momentum`).
      logical, intent(in) :: do_stress
         !! Fold the surface wind stress into the `k = nz` RHS row.
      logical, intent(in) :: do_drag
         !! Fold the bottom drag into the `k = k_bot_face` diagonal.
      real(wp), intent(in) :: rho0
         !! Boussinesq reference density for the `τ/ρ₀` stress conversion.
      real(wp), intent(in), optional :: tau_face(nu, nv)
         !! Wind stress (N/m²) on this face.  Present iff `do_stress`.
      real(wp), intent(in), optional :: lambda_bot(nu, nv)
         !! Bottom-drag Rayleigh rate λ (1/s) on this face.  Present iff
         !! `do_drag`.
      logical, intent(in) :: do_top
         !! Fold the ice-shelf top drag into the `k = nz` diagonal AND
         !! mask the wind-stress RHS by `(1 − cover_face)`.
      real(wp), intent(in), optional :: lambda_top(nu, nv)
         !! Top-drag Rayleigh rate λ (1/s) on this face.  Present iff
         !! `do_top`.
      real(wp), intent(in), optional :: cover_face(nu, nv)
         !! Face ice-cover mask (0 open, 1 under ice), the OR of the two
         !! abutting cells (`rdb_ocean_top_drag`).  Present iff `do_top`.
      logical, intent(in) :: hvel_mom6
         !! MOM6 HARMONIC_VISC parity for `h_u` + `h_shear` (see the slot-type
         !! docstring).  `.false.` => the historical arithmetic-`h_u` /
         !! `face_thick`-`dz` pair, bit-identical.
      real(wp), intent(in) :: hbbl_visc
         !! Bottom-layer scale for the `botfn` blend (MOM6 HBBL).
      logical, intent(in) :: bbl_glue
         !! MOM6 `bottomdraglaw` coupling parity (see the slot-type
         !! docstring).  Configure guarantees `hvel_mom6` when this is on;
         !! the bed piston then applies whether or not `do_drag`.
         !! `.false.` ⇒ bit-identical.
      logical, intent(in) :: hvel_upwind
         !! `.false.` = skip the near-bed upwind blend (pure harmonic
         !! hvel).  See the slot-type docstring.  `hvel_harmonic` branch only.
      logical, intent(in) :: hvel_harmonic
         !! MOM6 `HARMONIC_VISC`: `.false.` = MOM6's default `z_clear`
         !! face-thickness branch, `.true.` = the harmonic branch (slot-type
         !! docstring).  Only read when `hvel_mom6`.
      real(wp), intent(in) :: kv_bbl_face(nu, nv)
         !! BBL viscosity `kv_bbl` (m²/s) of this face (MOM6
         !! `visc%Kv_bbl_u`).  Only read when `bbl_glue`.
      real(wp), intent(in) :: bbl_thick_face(nu, nv)
         !! BBL thickness (m) of this face (MOM6 `visc%bbl_thick_u`); also
         !! the height-above-bed normaliser `1/I_Hbbl` under the glue, as in
         !! MOM6's `bottomdraglaw` branch.  Only read when `bbl_glue`.
      real(wp), intent(in) :: kv_bbl_bg
         !! MOM6 `KV`: the background viscosity the BBL viscosity replaces
         !! within botfn reach of the bed, `Kv_tot + (kv_bbl − KV)·botfn`.
      logical, intent(in) :: solve_momentum
         !! `.false.` = remnant-only mode: build the matrix and compute
         !! `visc_rem_out` WITHOUT touching `u_face` (the pre-substep
         !! visc_rem refresh, PGF_BUG.md §9 — MOM6 computes vertvisc_coef
         !! before btstep every stage, so its BT weights never lag; the
         !! stage-end producer alone leaves visc_rem ≡ 1 for the whole
         !! first stage, and the Δu corrector then deposits the spurious
         !! column-mean into wet layers).  `.true.` = the normal solve.
      logical, intent(in) :: do_remnant
         !! Fill `visc_rem_out` with the viscous remnant γ_k — the
         !! sensitivity of the post-friction layer velocity to a uniform
         !! barotropic acceleration, `γ_k ≡ (1/Δt)·∂u_k^{n+1}/∂Ā`.  γ
         !! solves the SAME tridiagonal system as the momentum solve
         !! (same `a_diag`/`b_diag`, and the already-factorized `c_diag`
         !! left behind by the momentum forward sweep) with RHS ≡ 1 —
         !! Roundabout's rows are pre-normalized by `h_k`, so the remnant RHS
         !! is the unit vector, NOT `h_u(k)` (MOM6's un-normalized
         !! convention).  Reusing the surviving factorization means γ
         !! cannot drift from the operator actually solved for momentum.
      real(wp), intent(inout), optional :: visc_rem_out(nu, nv, nz)
         !! Output γ_k, bottom-up (`k = 1` bed → γ smallest, `k = nz`
         !! surface → γ → 1).  Clamped to `min(γ, 1.0)` at production
         !! (cheap FP insurance on a quantity the maximum principle
         !! already bounds to (0, 1]).  Present iff `do_remnant`.
      integer, intent(in) :: k_top_face(nu, nv)
         !! `ms%k_top_u` / `k_top_v` — the index of the first layer that
         !! is LIVE on BOTH sides of this face, `nz` wherever nothing
         !! vanishes against the top.  The surface row of the column
         !! solve, the ice-shelf top-drag diagonal fold and the
         !! wind-stress RHS fold all sit on THIS row instead of on `nz`:
         !! under a quasi-geopotential coordinate beneath a shelf the
         !! rows above it are inert fillers, decoupled to the identity by
         !! the `zlevel_faces` gates, so a `dt*lambda_top` added at `nz`
         !! would damp a row that is not coupled to the ocean at all —
         !! and the drag rate `lambda_top` is itself captured at
         !! `k_top_face` by `rdb_ocean_top_drag`, so the two MUST agree.
         !! With `k_top_face ≡ nz` every expression below is the one it
         !! replaced, bit for bit.
      integer, intent(in) :: k_bot_face(nu, nv)
         !! `ms%k_bot_u` / `k_bot_v` — the first layer LIVE on BOTH sides
         !! of this face counting UP from the bed (`max` of its two
         !! columns), `1` wherever nothing vanishes against the bed.  The
         !! bed row of the column solve (no-flux-below BC, the implicit
         !! bottom-drag diagonal fold and the `bbl_glue` piston) sits on
         !! THIS row instead of on `k = 1`, and the rows below it — inert
         !! `z_fixed` bed fillers — are the identity.  The drag rate
         !! `lambda_bot` is itself captured at `k_bot_face` by
         !! `ocean_bottom_drag_compute_tendencies`, so the two MUST agree.
         !! The height-above-bed stack `zint` (BBL glue, upwind blend)
         !! also starts accumulating here.  With `k_bot_face ≡ 1` every
         !! expression below is the one it replaced, bit for bit.
      logical, intent(in) :: zlevel_faces
         !! `&vcoord_nml zfixed_closed_faces` — z-level partial steps.
         !!
         !! Changes TWO things about the face column, and only when on
         !! (`.false.` ⇒ every expression below is textually and
         !! bit-identically what it has always been):
         !!
         !! 1. the face thickness becomes `min(h_L, h_R)` instead of the
         !!    arithmetic (or MOM6-harmonic) mean — the standard
         !!    partial-cell choice, and the one that makes "this layer is
         !!    an inert FILLER on at least one side" visible locally as
         !!    `hvel(k) <= H_VANISHED`;
         !! 2. the tridiagonal coupling across an interface touching such
         !!    a layer is CUT (`alpha = beta = 0`), and the layer's own
         !!    row is forced to the identity.
         !!
         !! Without (2) an OPEN layer is frictionally coupled, every
         !! single step, to the CLOSED layer above or below it — whose
         !! velocity `mask_layer_velocities` has just zeroed — so the
         !! wall acts as a spurious side drag instead of as free-slip.
         !! The floor `max(hvel, H_VANISHED)` does not save it: it makes
         !! the coupling huge-but-finite (`dt·nu/1.5e-4`), i.e. a rigid
         !! glue to zero, which is the worst of the three options.
         !!
         !! This is the DYNAMIC twin of `metrics%open_u/open_v`: "filler
         !! on either side" is `min(h_L,h_R) <= H_VANISHED` here and
         !! `target_h <= H_VANISHED` on either side there, and under
         !! `z_fixed` the two agree (η is absorbed by the first LIVE
         !! layer).  It is derived locally rather than threaded as a 3-D
         !! mask because this kernel is reached through four flat-impl
         !! call sites and five dispatcher call sites, none of which sees
         !! `ocean_metrics_t`, and a scalar costs nothing.
         !! The equivalent tracer decoupling already exists, unconditionally,
         !! in `build_factorize_tracer_matrix`.
      logical, intent(in) :: do_corner
         !! Add the corner-staggered viscosity to every face interface.
         !! Gated INSIDE the single DC (no split loop — NVHPC penalty);
         !! `kv_corner` is guaranteed present when `do_corner`.
      real(wp), intent(in) :: kv_prandtl
         !! Scale on `kv_corner` (Kv = Pr·Kd).  Only read when
         !! `do_corner`.
      real(wp), intent(in), optional :: kv_corner(nx_cells + 1, ny_cells + 1, nz + 1)
         !! Corner-staggered interface viscosity (SW-corner convention:
         !! corner (i,j) is the SW corner of cell (i,j)).  A face reads
         !! the 2-point average of its two END corners: u-face (i,j) →
         !! corners (i,j)/(i,j+1); v-face (i,j) → corners (i,j)/(i+1,j)
         !! — the direct corner→face route, never via a tracer point.
         !! Present iff `do_corner`.

      integer :: i, j, k, ktop, kb, kb_live, i_left, i_right, j_below, j_above, i_c2, j_c2
      real(wp) :: zacc, zcol_l, zcol_r, d_min, i_hbbl, kv_bbl, bbl_thick
      ! Sliding window over the face-thickness / height-above-bed stack:
      ! `hv_0/hv_p` = hvel(k/k+1), `zi_0/zi_p` = zint(k/k+1).  It
      ! replaces the `hvel(NZ_STACK_MAX)` / `zint(NZ_STACK_MAX)` `local()`
      ! stacks, which at nz ~ 100 live in GPU local memory.  With the fused
      ! sweep below: 2082 -> 957 us per launch on the NK=100 double gyre
      ! (512x256, GH200).
      real(wp) :: hv_0, hv_p, zi_0, zi_p
      ! Interface couplings: `_b` = interface below row k, `_t` = above.  Each
      ! is formed ONCE, as `num = dt*nu_eff` and `dz_eff`, and read by both
      ! rows it couples (`alpha = num/(hf_k*dz)`), from the same operands
      ! either row used to form it from.
      real(wp) :: num_b, dz_b, num_t, dz_t, kv_cnr
      logical :: cut_b, cut_t
      real(wp) :: hf_k, alpha, beta, a_k, b_k, c_k, r_k, denom, wind_open
      real(wp) :: cp_prev, dp_prev, gp_prev, cp, dp, gp, u_next, g_next
      real(wp) :: inv_rho0
      real(wp), parameter :: EPS_HVEL = 1.0e-30_wp
         !! MOM6 `h_neglect` analogue in the harmonic mean / HBBL inverse.
      real(wp), parameter :: STRESS_H_MIN = 1.0e-3_wp
         !! Floor on the surface-layer face thickness in the `τ/(ρ₀·h)`
         !! stress conversion — matches `ocean_surface_stress_t%h_min` so
         !! the implicit fold reduces to the explicit kernel in the
         !! well-resolved limit.

      inv_rho0 = 1.0_wp/rho0

      ! ONE upward pass per face column: layer k+1's face thickness, row k's
      ! coefficients and row k's Thomas forward step (momentum c'/d' and the
      ! remnant g share the denominator) in the same iteration, then one
      ! back-substitution for both.  `a_diag` / `b_diag` are never stored —
      ! nothing reads them after the solve; `c_diag` / `rhs` end holding c'
      ! and d', as before.  The arithmetic is the separate-sweep version's,
      ! expression for expression (bit-identical).
      do concurrent(j=1:nv, i=1:nu) &
         local(i_left, i_right, j_below, j_above, i_c2, j_c2, k, ktop, kb, kb_live, &
               zacc, zcol_l, zcol_r, d_min, i_hbbl, kv_bbl, bbl_thick, &
               hv_0, hv_p, zi_0, zi_p, num_b, dz_b, num_t, dz_t, kv_cnr, &
               cut_b, cut_t, hf_k, alpha, beta, a_k, b_k, c_k, r_k, denom, wind_open, &
               cp_prev, dp_prev, gp_prev, cp, dp, gp, u_next, g_next)
         ! Neighbour cell indices for averaging h.  (i_c2, j_c2) is the
         ! SECOND end corner of this face for the corner-viscosity
         ! add-on (the first is (i, j) on both staggerings): a u-face
         ! runs south→north (corners (i,j)/(i,j+1)), a v-face west→east
         ! (corners (i,j)/(i+1,j)).
         if (x_face) then
            i_left = max(1, i - 1)
            i_right = min(nx_cells, i)
            j_below = j
            j_above = j
            i_c2 = i
            j_c2 = j + 1
         else
            i_left = i
            i_right = i
            j_below = max(1, j - 1)
            j_above = min(ny_cells, j)
            i_c2 = i + 1
            j_c2 = j
         end if

         ! First live layer of this face counting up from the bed (`1`
         ! everywhere but a `z_fixed` bed filler band), and the first live
         ! layer counting down from the top (`nz` off a rigid top).
         kb = k_bot_face(i, j)
         ktop = k_top_face(i, j)

         ! ---- This face's bottom boundary layer (MOM6 visc%Kv_bbl_u,
         ! visc%bbl_thick_u).  Under the glue the height above the bed is
         ! normalised by the BBL thickness, as MOM6's bottomdraglaw branch
         ! does (`I_Hbbl = 1/bbl_thick`, vertvisc_coef); otherwise by
         ! `hbbl_visc` (MOM6 HBBL).  The legacy glue fills `bbl_thick_face`
         ! with `hbbl_visc`, so its arithmetic is unchanged.
         kv_bbl = 0.0_wp
         bbl_thick = hbbl_visc
         if (bbl_glue) then
            kv_bbl = kv_bbl_face(i, j)
            bbl_thick = bbl_thick_face(i, j)
         end if
         i_hbbl = 1.0_wp/(bbl_thick + EPS_HVEL)

         ! MOM6 `HARMONIC_VISC = False` (`z_clear`) branch: the two columns'
         ! interface heights, from their beds up.  MOM6 starts them at the
         ! STATIC depth `-bathyT`; here the column sums stand in for it, which
         ! differs only by the free surface and cancels in `z_clear` wherever
         ! the two cells' eta agree (an O(d eta) difference across a face).
         zacc = 0.0_wp
         zcol_l = 0.0_wp
         zcol_r = 0.0_wp
         if (hvel_mom6 .and. .not. hvel_harmonic) then
            do k = 1, nz
               zcol_l = zcol_l - h_layer(i_left, j_below, k)
               zcol_r = zcol_r - h_layer(i_right, j_above, k)
            end do
         end if
         d_min = min(-zcol_l, -zcol_r)

         ! #178: the static bed index `kb` (`ms%k_bot_u/v`, filled once at
         ! configure from the η=0 target) can be a vanished z_fixed bed
         ! filler once a free-surface drawdown eats the nominal bed sliver
         ! on the live run (`h <= H_VANISHED` on one or both sides of the
         ! face) — `kb` itself never moves, so the bed BC / drag piston
         ! placed there would be glued to a massless row and decoupled
         ! from the live rows above by the `zlevel_faces` gates (visc_rem
         ! = 1 exactly, no bottom drag reaches the water).  `kb_live` is
         ! the first row at/above `kb` whose face thickness is LIVE by the
         ! shared `rdb_blf_is_live` criterion (`rdb_bed_live_face.inc` —
         ! `ocean_bottom_drag_compute_tendencies` walks the SAME search so
         ! the drag law samples the row the solver glues the BC to) —
         ! fall back to `kb` itself when every row at/above it is vanished
         ! (should not happen for a face with any live water, but keeps
         ! the row index in-bounds).  Only the BED-BC ROW PLACEMENT moves
         ! to `kb_live`; the height-above-bed accumulator (`zint`, gated
         ! on `k >= kb`) keeps referencing the true physical bed `kb`,
         ! unchanged.  `kb_live = kb` ⇒ bit-identical to the pre-#178 path.
         kb_live = kb
         do k = kb, nz
            if (rdb_blf_is_live(h_layer(i_left, j_below, k), h_layer(i_right, j_above, k), &
                                zlevel_faces)) then
               kb_live = k
               exit
            end if
         end do

         ! Prime the window with layer 1.
         call vdiff_face_layer(h_layer(i_left, j_below, 1), h_layer(i_right, j_above, 1), &
                               u_face(i, j, 1), 1, kb, hvel_mom6, hvel_harmonic, hvel_upwind, &
                               zlevel_faces, d_min, i_hbbl, zacc, zcol_l, zcol_r, hv_p, zi_p)
         hv_0 = 0.0_wp
         zi_0 = 0.0_wp
         num_t = 0.0_wp
         dz_t = 1.0_wp
         cut_t = .false.
         cp_prev = 0.0_wp
         dp_prev = 0.0_wp
         gp_prev = 0.0_wp

         do k = 1, nz
            hv_0 = hv_p
            zi_0 = zi_p
            num_b = num_t
            dz_b = dz_t
            cut_b = cut_t
            if (k < nz) then
               call vdiff_face_layer(h_layer(i_left, j_below, k + 1), h_layer(i_right, j_above, k + 1), &
                                     u_face(i, j, k + 1), k + 1, kb, hvel_mom6, hvel_harmonic, &
                                     hvel_upwind, zlevel_faces, d_min, i_hbbl, zacc, zcol_l, zcol_r, &
                                     hv_p, zi_p)
               ! Interface k+1/2 is read by row k (alpha: bed or interior
               ! row) and row k+1 (beta: interior or top row), so it is
               ! needed exactly for kb_live <= k <= ktop.
               if (k >= kb_live .and. k <= ktop) then
                  ! Corner-viscosity add-on (vertex kappa-shear Kv seam):
                  ! direct 2-point end-corner average onto this face, folded
                  ! in BEFORE the BBL glue.
                  kv_cnr = 0.0_wp
                  if (do_corner) kv_cnr = kv_prandtl*(0.5_wp*(kv_corner(i, j, k + 1) + &
                                                              kv_corner(i_c2, j_c2, k + 1)))
                  call vdiff_interface_coupling(hv_0, hv_p, zi_0, dt, hvel_mom6, use_harmonic, &
                                                bbl_glue, kv_bbl, kv_bbl_bg, bbl_thick, &
                                                zlevel_faces, do_corner, kv_cnr, &
                                                kv_centre(i_left, j_below, k + 1), &
                                                kv_centre(i_right, j_above, k + 1), &
                                                num_t, dz_t, cut_t)
               end if
            end if

            ! ---- Row k ----
            ! Rows below the first live bed row (`k < kb_live`: inert
            ! `z_fixed` bed fillers) and above the first live top row
            ! (`k > ktop`: inert fillers inside an ice draft) are the
            ! identity, `rhs = u^n`, so the sweep returns them unchanged.
            ! The row kinds are mutually exclusive `if` blocks over an
            ! init, not an if/else chain: gfortran 15.1 miscompiles a
            ! `local()` scalar reassigned across the arms of an if/else
            ! inside `do concurrent`.
            a_k = 0.0_wp
            b_k = 1.0_wp
            c_k = 0.0_wp
            ! Floor the face thickness at H_VANISHED before it enters the
            ! α/β denominators: an exactly-collapsed layer (hvel = 0) with
            ! kv > 0 would otherwise make α = Inf ⇒ NaN γ.  A no-op wherever
            ! hvel ≫ H_VANISHED.
            hf_k = max(hv_0, H_VANISHED)
            ! Bed row `k = kb_live` (when it is not above `ktop`): no flux
            ! below.  SINGLE-LAYER COLUMN (`kb_live = nz`): no interface
            ! above either, so `alpha = 0` and the surface is a pure
            ! stress-Neumann BC added below.
            if (k == kb_live .and. k <= ktop) then
               alpha = 0.0_wp
               if (kb_live < nz) then
                  alpha = num_t/(hf_k*dz_t)
                  if (cut_t) alpha = 0.0_wp
               end if
               c_k = -alpha
               b_k = 1.0_wp + alpha
               ! Bottom-drag stress BC (Roundabout bed = k=kb_live; MIRROR
               ! of MOM6 k=nz).  Under `bbl_glue` the drag is MOM6's viscous
               ! PISTON `kv_bbl/min(hvel/2, bbl_thick)` (`a_cpl(nz+1)`),
               ! which diverges as the bed layer thins and so anchors a
               ! grounded sliver to rest (PGF_BUG.md §9.4); it replaces the
               ! h-independent Rayleigh fold `dt·λ_bot`, used otherwise.  A
               ! drag is a SINK ⇒ POSITIVE diagonal add ⇒ unconditionally
               ! stable on thin shelf bottoms.
               if (bbl_glue) then
                  b_k = b_k + dt*kv_bbl/(hf_k*(min(0.5_wp*hv_0, bbl_thick) + EPS_HVEL))
               end if
               if (do_drag .and. .not. bbl_glue) then
                  b_k = b_k + dt*lambda_bot(i, j)
               end if
            end if
            ! Top row `k = ktop` above the bed row: no flux above.
            if (k == ktop .and. k > kb_live) then
               beta = num_b/(hf_k*dz_b)
               if (cut_b) beta = 0.0_wp
               a_k = -beta
               b_k = 1.0_wp + beta
            end if
            ! Interior rows.
            if (k > kb_live .and. k < ktop) then
               beta = num_b/(hf_k*dz_b)
               alpha = num_t/(hf_k*dz_t)
               if (cut_b) beta = 0.0_wp
               if (cut_t) alpha = 0.0_wp
               a_k = -beta
               c_k = -alpha
               b_k = 1.0_wp + alpha + beta
            end if

            r_k = u_face(i, j, k)
            if (k == ktop) then
               ! Ice-shelf top-drag stress BC (the MIRROR of the bed row's
               ! `dt·λ_bot` fold).  λ_top is already cover- and wet-masked,
               ! so an OPEN face adds an exact zero.  On a single-live-layer
               ! face this row is also the bed row and carries the bottom
               ! drag, the top drag and the wind together.
               if (do_top) b_k = b_k + dt*lambda_top(i, j)
               ! Surface wind-stress Neumann BC (MIRROR of MOM6 k=1): the
               ! kinematic stress τ/ρ₀ is a RHS source, pre-normalized by
               ! the face thickness, so with no interior coupling it reduces
               ! EXACTLY to the explicit `du = dt·τ/(ρ₀·h_top)`.  `tau_face`
               ! is masked by the face `min` of the two cells' wet mask (no
               ! stress at a land face) and by `(1 − cover_face)` under an
               ! ice shelf.  The cover factor is belt and braces — P2c
               ! already zeroes `tau` on every covered face at its source
               ! (`ocean_surface_stress_apply_cover`, face-for-face the same
               ! OR rule as `top_drag_fill_face_cover_impl`) — kept so the
               ! fold stays correct STANDALONE against a future forcing path
               ! that writes `tau` after the configure-time mask.
               if (do_stress) then
                  wind_open = 1.0_wp
                  if (do_top) wind_open = 1.0_wp - cover_face(i, j)
                  r_k = r_k + dt*tau_face(i, j)*wind_open &
                        *min(wet_cell(i_left, j_below), wet_cell(i_right, j_above)) &
                        *inv_rho0/max(hf_k, STRESS_H_MIN)
               end if
            end if

            ! ---- Thomas forward step ----
            ! remnant-only mode (`solve_momentum = .false.`) skips the d'
            ! leg; c' is always formed (the remnant back-substitution reads
            ! it).  The viscous remnant γ (`do_remnant`) solves the SAME
            ! operator with RHS ≡ 1 — Roundabout's rows are pre-normalized
            ! by h_k, so the unit vector, NOT MOM6's h_u(k) — and shares the
            ! denominator, so γ cannot drift from the operator actually
            ! solved for momentum.
            if (k == 1) then
               cp = c_k/b_k
               dp = r_k/b_k
               gp = 1.0_wp/b_k
            end if
            if (k > 1) then
               denom = b_k - a_k*cp_prev
               cp = c_k/denom
               dp = (r_k - a_k*dp_prev)/denom
               gp = (1.0_wp - a_k*gp_prev)/denom
            end if
            c_diag(i, j, k) = cp
            if (solve_momentum) rhs(i, j, k) = dp
            if (do_remnant) visc_rem_out(i, j, k) = gp
            cp_prev = cp
            dp_prev = dp
            gp_prev = gp
         end do

         ! ---- Back-substitution: momentum into u_face, γ into
         ! visc_rem_out, in one sweep (c' read once).  γ is clamped to
         ! `min(γ, 1)` as each element is finalized (MOM6 clamps at
         ! consumption; the maximum principle bounds γ to (0, 1], so this is
         ! cheap FP insurance). ----
         u_next = dp_prev
         g_next = min(gp_prev, 1.0_wp)
         if (solve_momentum) u_face(i, j, nz) = u_next
         if (do_remnant) visc_rem_out(i, j, nz) = g_next
         do k = nz - 1, 1, -1
            cp = c_diag(i, j, k)
            if (solve_momentum) then
               u_next = rhs(i, j, k) - cp*u_next
               u_face(i, j, k) = u_next
            end if
            if (do_remnant) then
               g_next = min(visc_rem_out(i, j, k) - cp*g_next, 1.0_wp)
               visc_rem_out(i, j, k) = g_next
            end if
         end do
      end do
   end subroutine diffuse_velocity_columns_impl

   pure function face_thick(h_a, h_b, use_harmonic) result(dz)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Face-thickness for the vdiff implicit operator.
      !!
      !! When `use_harmonic = .false.` (default): arithmetic mean
      !! `0.5 · (h_a + h_b)`.  Standard MOM6 / ROMS behaviour.
      !!
      !! When `use_harmonic = .true.`: harmonic mean
      !! `2 · h_a · h_b / max(h_a + h_b, eps)`.  Better-conditioned
      !! when one of `h_a`, `h_b` is small — arithmetic mean is
      !! dominated by the thicker neighbour, generating stiff
      !! tridiagonal coefficients at thin/vanishing layers.
      !! Equivalent to arithmetic when `h_a = h_b`.
      !!
      !! Marked `pure` + `!$acc routine seq` so NVHPC can inline the
      !! body into the `do concurrent` callers without descriptor
      !! marshalling.
      !$acc routine seq
      real(wp), intent(in) :: h_a, h_b
      logical, intent(in) :: use_harmonic
      real(wp) :: dz
      real(wp), parameter :: EPS = 1.0e-30_wp
      real(wp) :: sum_h
      if (use_harmonic) then
         sum_h = h_a + h_b
         if (sum_h > EPS) then
            dz = 2.0_wp*h_a*h_b/sum_h
         else
            dz = 0.5_wp*sum_h
         end if
      else
         dz = 0.5_wp*(h_a + h_b)
      end if
   end function face_thick

   pure subroutine vdiff_face_layer(hl_c, hr_c, u_k, k, kb, hvel_mom6, hvel_harmonic, &
                                    hvel_upwind, zlevel_faces, d_min, i_hbbl, &
                                    zacc, zcol_l, zcol_r, hvel, zint)
      !! Face thickness `hvel` (MOM6 `hvel`) and normalized height above the
      !! bed `zint` (MOM6 `z_i`, the TOP interface of layer k) of layer k of
      !! one face column, advancing the upward accumulators `zacc` /
      !! `zcol_l` / `zcol_r`.  Called for k = 1, 2, ..., nz in order by
      !! `diffuse_velocity_columns_impl`.  k = 1 is the BED (MOM6 counts
      !! from the surface), so the height accumulates UPWARD.
      !$acc routine seq
      real(wp), intent(in) :: hl_c, hr_c
         !! The two abutting cells' thickness of layer k.
      real(wp), intent(in) :: u_k
         !! Face velocity of layer k (u^n), for the upwind blend's sign test.
      integer, intent(in) :: k, kb
         !! Layer index and the face's static first-live bed layer.
      logical, intent(in) :: hvel_mom6, hvel_harmonic, hvel_upwind, zlevel_faces
      real(wp), intent(in) :: d_min
         !! Shallower of the two column depths (`z_clear` branch only).
      real(wp), intent(in) :: i_hbbl
         !! 1 / (height-above-bed normaliser).
      real(wp), intent(inout) :: zacc, zcol_l, zcol_r
         !! Upward accumulators: harmonic height above the live bed, and the
         !! two columns' interface heights (`z_clear` branch).
      real(wp), intent(out) :: hvel, zint
      real(wp), parameter :: EPS_HVEL = 1.0e-30_wp
      real(wp) :: h_harm, h_arith, h_delta, z_clear, z2, botfn

      ! Not read off `hvel_mom6` (configure guarantees `hvel_mom6`
      ! whenever `bbl_glue`, the only other reader).
      zint = 0.0_wp
      if (hvel_mom6 .and. .not. hvel_harmonic) then
         ! MOM6 vertvisc_coef, HARMONIC_VISC = False branch.  `zcol_*` is
         ! the height of the TOP of layer k in each cell; `z_clear` is the
         ! height of the higher of the two above the SHALLOWER bed
         ! (`-d_min`), so a face layer below the shallower bed of a step
         ! sits at z_clear <= 0 and only `zh` (harmonic, ~0 against a
         ! filler) lifts it.  `harm_BL_val = 0` (MOM6 HARMONIC_BL_SCALE
         ! default) ⇒ `z2_wt = 1`.
         h_harm = 2.0_wp*hl_c*hr_c/(hl_c + hr_c + EPS_HVEL)
         h_arith = 0.5_wp*(hl_c + hr_c)
         h_delta = hr_c - hl_c
         zcol_l = zcol_l + hl_c
         zcol_r = zcol_r + hr_c
         if (k >= kb) zacc = zacc + h_harm
         z_clear = max(zcol_l, zcol_r) + d_min
         zint = max(zacc, z_clear)*i_hbbl
         hvel = h_arith
         ! Flow from the THIN side to the THICK side: near the bed the
         ! face takes the harmonic (thin-side) thickness, so a
         ! near-massless layer cannot be drained through a thick face.
         if (u_k*h_delta > 0.0_wp) then
            z2 = zint
            botfn = 1.0_wp/(1.0_wp + 0.09_wp*z2*z2*z2*z2*z2*z2)
            hvel = (1.0_wp - botfn)*h_arith + botfn*h_harm
         end if
      else if (hvel_mom6) then
         h_harm = 2.0_wp*hl_c*hr_c/(hl_c + hr_c + EPS_HVEL)
         h_arith = 0.5_wp*(hl_c + hr_c)
         h_delta = hr_c - hl_c
         hvel = h_harm
         ! Upwind bias: only when the flow runs from the THICK side to the
         ! THIN side does the near-bed face take the arithmetic (donor)
         ! thickness.  Without this the harmonic mean over-suppresses
         ! asymmetrically -- the reason bare-harmonic attempts regressed
         ! (LAGRANGIAN_PGF_BUG.md 6.2).
         !
         ! `hvel_upwind = .false.` disables the blend (pure harmonic hvel).
         ! The sign test keys on u itself, so at (near-)rest it flip-flops
         ! faces between harmonic and arithmetic on ROUNDOFF-sign
         ! velocities, collapsing the BBL glue at whichever faces flip —
         ! measured ×200-800/stage residual amplification on the rest-state
         ! reproducer (PGF_BUG.md §9.8).  MOM6 carries the same test and the
         ! same hazard; it never bites there only because its state stays
         ! at 1e-14.
         if (hvel_upwind .and. u_k*h_delta < 0.0_wp) then
            z2 = zacc
            botfn = 1.0_wp/(1.0_wp + 0.09_wp*z2*z2*z2*z2*z2*z2)
            hvel = (1.0_wp - botfn)*h_harm + botfn*h_arith
         end if
         ! Height above the LIVE bed: the inert fillers below `kb` do not
         ! count (bit-identical at kb = 1).
         if (k >= kb) zacc = zacc + h_harm*i_hbbl
         zint = zacc
      else
         hvel = 0.5_wp*(hl_c + hr_c)
      end if
      ! z-level partial steps: the face column is the OVERLAP of the two
      ! cell columns, so its thickness is the min.  A filler on either side
      ! then reads at or below H_VANISHED and the coupling gates cut it out
      ! of the solve.
      if (zlevel_faces) hvel = min(hl_c, hr_c)
   end subroutine vdiff_face_layer

   pure subroutine vdiff_interface_coupling(hv_lo, hv_hi, zint_lo, dt, hvel_mom6, use_harmonic, &
                                            bbl_glue, kv_bbl, kv_bbl_bg, bbl_thick, &
                                            zlevel_faces, do_corner, kv_cnr, kv_l, kv_r, &
                                            num, dz, cut)
      !! Vertical-friction coupling across the interface between face layers
      !! `lo` (below) and `hi` (above): `num = dt·ν_eff` and the shear
      !! distance `dz`, so that each of the two rows it couples takes
      !! `num/(hf_k·dz)` as its alpha (row lo) / beta (row hi).  `cut` is the
      !! z-level partial-step decoupling.
      !$acc routine seq
      real(wp), intent(in) :: hv_lo, hv_hi
         !! Unfloored face thicknesses of the two layers.
      real(wp), intent(in) :: zint_lo
         !! Normalized height above the bed of the interface (layer lo's top).
      real(wp), intent(in) :: dt
      logical, intent(in) :: hvel_mom6, use_harmonic, bbl_glue, zlevel_faces, do_corner
      real(wp), intent(in) :: kv_bbl, kv_bbl_bg, bbl_thick
      real(wp), intent(in) :: kv_cnr
         !! Corner-viscosity add-on, `Pr·(corner average)`.  Only read when
         !! `do_corner`.
      real(wp), intent(in) :: kv_l, kv_r
         !! The two cells' interface diffusivity.
      real(wp), intent(out) :: num, dz
      logical, intent(out) :: cut
      real(wp) :: hf_lo, hf_hi, nu_f, botfn_int

      ! Floored face thicknesses (see the row build).
      hf_lo = max(hv_lo, H_VANISHED)
      hf_hi = max(hv_hi, H_VANISHED)
      if (hvel_mom6) then
         dz = 0.5_wp*(hf_lo + hf_hi)   ! MOM6 h_shear: arithmetic of hvels
      else
         dz = face_thick(hf_lo, hf_hi, use_harmonic)
      end if
      nu_f = 0.5_wp*(kv_l + kv_r)
      if (do_corner) nu_f = nu_f + kv_cnr
      ! BBL glue (MOM6 find_coupling_coef): within botfn reach of the bed
      ! the BBL viscosity REPLACES the background one, `Kv_tot + (kv_bbl −
      ! KV)·botfn` (shear / boundary-layer contributions kept; floored at 0
      ! for a `kv_bbl` below `KV`, where MOM6 relies on Kv_tot >= KV), and
      ! the shear distance is capped toward bbl_thick — grounded stacks
      ! (zint ≈ 0) become rigidly coupled, which is the mechanism that
      ! absorbs the spurious grounded-layer PGF every step (PGF_BUG.md §9).
      if (bbl_glue) then
         botfn_int = 1.0_wp/(1.0_wp + 0.09_wp*zint_lo*zint_lo*zint_lo*zint_lo*zint_lo*zint_lo)
         nu_f = max(0.0_wp, nu_f + (kv_bbl - kv_bbl_bg)*botfn_int)
         if (dz > bbl_thick) then
            dz = (1.0_wp - botfn_int)*dz + botfn_int*bbl_thick
         end if
      end if
      num = dt*nu_f
      ! z-level partial steps: no shear across an interface that touches a
      ! CLOSED face layer.  The tracer twin does exactly this,
      ! unconditionally, in `build_factorize_tracer_matrix`.
      ! vanished-ok: row decoupling across a z-level wall (the closed-face rule) —
      ! a momentum-grid question, not a tracer concentration.
      cut = zlevel_faces .and. (hv_lo <= H_VANISHED .or. hv_hi <= H_VANISHED)
   end subroutine vdiff_interface_coupling

   pure function ocean_vdiff_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the implicit vertical diffusion slot (0 when
      !! unallocated).
      class(ocean_vdiff_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%bbl_mix%bytes() &
               + this%a_diag_t%bytes() &
               + this%b_diag_t%bytes() &
               + this%c_diag_t%bytes() &
               + this%rhs_t%bytes() &
               + this%a_diag_u%bytes() &
               + this%b_diag_u%bytes() &
               + this%c_diag_u%bytes() &
               + this%rhs_u%bytes() &
               + this%a_diag_v%bytes() &
               + this%b_diag_v%bytes() &
               + this%c_diag_v%bytes() &
               + this%rhs_v%bytes() &
               + this%kv_scalar_buf%bytes() &
               + arr_bytes(this%kv_bbl_u) + arr_bytes(this%kv_bbl_v) &
               + arr_bytes(this%bbl_thick_u) + arr_bytes(this%bbl_thick_v) &
               + arr_bytes(this%bbl_conc_t) + arr_bytes(this%bbl_conc_s)
   end function ocean_vdiff_bytes

#include "rdb_bed_live_face.inc"

end module rdb_ocean_vdiff
