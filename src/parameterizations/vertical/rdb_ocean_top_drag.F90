!! Ice-shelf TOP drag: the momentum sink at the ocean's upper boundary
!! where an ice shelf sits on it.  The mirror image of
!! `rdb_ocean_bottom_drag` about the middle of the column — under the
!! ROMS-style bottom-up convention this codebase uses everywhere, the bed
!! is `k = 1` and the ice base is `k = nz`.
!!
!! ## Why a separate module and not a second mode of the bottom drag
!!
!! The two kernels share their ALGEBRA but nothing else.  The bottom-drag
!! slot already carries the channel (side-wall) drag, the `bed_factor`
!! bed-layer multiplier, and the `lambda_bot_u/v` fields the BBL-glue
!! bottom-BC consumes; bolting a second boundary onto it would make every
!! one of those fields ask "which end?".  Top drag also needs something
!! the bottom never does: a per-FACE ICE-COVER mask, because the upper
!! boundary is a free surface over open water and a no-slip-like wall
!! under a shelf, in the same domain, at the same time.  So: its own slot,
!! its own knob group (`&ocean_tdrag_nml`), and a deliberate, documented
!! one-to-one correspondence with the bottom-drag algebra that the
!! mirror-symmetry test (`test_ocean_top_drag`) pins numerically.
!!
!! ## The physics
!!
!! Under an ice shelf the ocean's top boundary is a solid, no-slip-like
!! wall.  The stress it exerts is quadratic in the boundary-layer flow,
!!
!!     tau_top = rho_0 * C_d * |u_top| * u_top          (opposing the flow)
!!
!! with `C_d = 2.5e-3` the ISOMIP+ protocol value (Asay-Davis et al.
!! (2016) Table 4, which prescribes the SAME quadratic law at the top and
!! the bottom), and the same coefficient Jenkins, Nicholls & Corr (2010)
!! use with the melt law.  A linear (Rayleigh) form is also provided for
!! analytic work — it is the form whose spin-down has a closed-form
!! exponential, which is how the distributed mode is verified.
!!
!! ## The FACE cover rule (stated once, here)
!!
!! `metrics%cover_frac` is a CELL-CENTRED binary 0/1 field.  Velocities
!! live on faces, so the mask has to be projected onto faces, and the
!! projection is a choice at the calving front where one neighbour is
!! covered and the other is open.  This module takes the **OR**:
!!
!!     cover_u(i,j) = max(cover_frac(i-1,j), cover_frac(i,j))
!!     cover_v(i,j) = max(cover_frac(i,j-1), cover_frac(i,j))
!!
!! so a face is "under ice" if EITHER of its two cells is.  The frontal
!! face therefore FEELS the drag.  The alternative (`min`, AND) leaves the
!! frontal face frictionless, which puts a slip line exactly where the
!! outflow jet leaves the cavity — the one place in the domain where the
!! top stress is largest and where a spurious free-slip band would be
!! systematically rectified into the overturning.  Erring toward too much
!! drag over one face-width is the conservative error; erring toward none
!! is not.  (A partial-cover area weighting is the v2 refinement;
!! `cover_frac` is binary in v1, so `max` IS the logical OR.)
!!
!! `cover_u`/`cover_v` are filled ONCE, at configure, because
!! `cover_frac` is static geometry (`&ocean_cavity_dyn_nml`, filled in the
!! IC seed and never touched again).  The per-step kernel reads the face
!! masks, never `metrics`.
!!
!! ## Distribution over a top boundary layer (`htbl`)
!!
!! Exactly the `hbbl` mode of the bottom drag, reflected: with
!! `htbl > 0` the stress is spread over the top `htbl` metres of the
!! column instead of being dumped into layer `nz` alone.  This matters
!! near a grounding line, where a sigma coordinate makes the top layer
!! arbitrarily thin and a single-layer explicit drag rate `C_d|U|/h_nz`
!! becomes arbitrarily stiff.  Killworth & Edwards (1999) is the reference
!! for a boundary layer of prescribed thickness in a layered model.
!!
!! ## Implicit forms — there are two, and they are different things
!!
!!   * `&ocean_tdrag_nml implicit` — backward-Euler in the DRAGGED
!!     velocity, formed inside this kernel: the tendency written is
!!     `-lambda*u/(1 + dt*lambda)` so the ordinary `u += dt*du_drag`
!!     apply reproduces `u^{n+1} = u/(1 + dt*lambda)` exactly.
!!     Unconditionally stable for any `h`.  Mirrors
!!     `&ocean_bdrag_nml implicit`.
!!   * `&ocean_vdiff_nml implicit_top_drag` — folds the rate into the
!!     vertical-friction tridiagonal's `k = nz` DIAGONAL, where the wind
!!     stress already owns the RHS.  Mirrors `&ocean_vdiff_nml
!!     implicit_drag`.  That path sets `implicit_fold` here, which fills
!!     `lambda_top_u/v` and makes the driver SKIP the explicit apply.
!!
!! The two are mutually exclusive at configure — both would damp the top
!! layer, which is a double count, not a stronger drag.
!!
!! ## What this module exports for the boundary-layer schemes
!!
!! `stress_top(nx,ny)` is the cell-centred magnitude of the top stress
!! (N/m^2), device-resident, refreshed by the same kernel that writes the
!! tendencies.  It is written and NOT read here: KPP/EPBL still take their
!! `u_*` from `stress_mag`, which this PR deliberately does not touch (a
!! later PR blends `ustar_shelf = sqrt(stress_top/rho_0)` into the
!! boundary-layer schemes).  Exposing it now means that PR is a consumer
!! change only.
!!
!! ## Citations (papers, never another model's source)
!!
!!   * Asay-Davis, X. S. et al. (2016): Geosci. Model Dev. 9, 2471-2497
!!     (ISOMIP+; quadratic top AND bottom drag, `C_d = 2.5e-3`).
!!   * Jenkins, A., Nicholls, K. W. and Corr, H. F. J. (2010): J. Phys.
!!     Oceanogr. 40, 2298-2312.
!!   * Killworth, P. D. and Edwards, N. R. (1999): J. Phys. Oceanogr. 29,
!!     1221-1238 (boundary layer of prescribed thickness).
module rdb_ocean_top_drag
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_top_drag_t
   public :: ocean_top_drag_compute_tendencies
   public :: ocean_top_drag_apply_tendencies
   public :: top_drag_fill_face_cover_impl
   public :: top_drag_tendencies_impl
   public :: top_drag_stress_mag_impl
   public :: parse_tdrag_variant
   public :: tdrag_variant_is_implemented

   integer, parameter, public :: TDRAG_LINEAR = 1
      !! Linear Rayleigh top drag, `du/dt = -r*u`.  The form with a
      !! closed-form spin-down; used by the analytic decay gate.
   integer, parameter, public :: TDRAG_QUADRATIC = 2
      !! Quadratic (log-layer) top drag, `du/dt = -C_d*|U|*u/h`.
      !! Production default when the group is enabled.
   integer, parameter, public :: TDRAG_INVALID = -1
      !! Sentinel for an unrecognised `form` string.  Linear and quadratic
      !! carry different coefficient DIMENSIONS and different decay laws,
      !! so a typo must abort rather than silently pick one.

   type :: ocean_top_drag_t
      !! Ice-shelf top-drag slot.  Every array is full size when
      !! `enable`, a `(1,1)` / `(1,1,1)` placeholder otherwise — the
      !! `ocean_cavity_flux_t` gating convention, latched in
      !! `ocean_state_init_from_config` BEFORE `init` so the allocation
      !! gate can read it.
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Always test this, never
         !! `allocated(...)`.
      logical :: enable = .false.
         !! `&ocean_tdrag_nml enable`, latched before `init`.  Off ⇒
         !! placeholders, no kernel launch, byte-identical.
      integer :: variant = TDRAG_QUADRATIC
         !! Active drag variant (`TDRAG_*`).
      real(wp) :: r_linear = 0.0_wp
         !! Linear Rayleigh coefficient (1/s).  Zero disables the linear
         !! branch even when the variant tag selects it.
      real(wp) :: c_drag = 0.0_wp
         !! Quadratic drag coefficient (dimensionless).  ISOMIP+ 2.5e-3.
         !! Zero disables the quadratic branch.
      real(wp) :: h_min = 1.0e-3_wp
         !! Floor on the top-layer face thickness inside the `u/h_top`
         !! division — keeps the kernel finite when the top layer pinches
         !! out.  Matches the bottom-drag `h_min`.
      real(wp) :: htbl = 0.0_wp
         !! Top-boundary-layer thickness (m) the stress is distributed
         !! over (the mirror of `hbbl`).  Zero (default) = layer-`nz`-only
         !! mode.
      real(wp) :: drag_bg_vel = 0.0_wp
         !! Background velocity floor (m/s) in the quadratic speed,
         !! `|U_eff| = max(drag_bg_vel, |U_tbl|)` (the mirror of MOM6
         !! DRAG_BG_VEL).  UNLIKE the bottom-drag slot, this is honoured
         !! in BOTH the layer-only and distributed modes — the two kernels
         !! are one code path here.  Zero (default) ⇒ no floor ⇒ the
         !! layer-only quadratic branch is the exact algebraic mirror of
         !! `ocean_bottom_drag_compute_tendencies`.
      real(wp) :: tbl_thick_min = 0.0_wp
         !! Minimum effective TBL thickness (m) in the `stress / h_tbl`
         !! denominator (mirror of BBL_THICK_MIN).  Zero (default) falls
         !! back to `h_min`.
      logical :: implicit = .false.
         !! Backward-Euler top drag in this kernel: the tendency is formed
         !! as `-lambda*u/(1 + dt*lambda)` so the standalone apply gives
         !! `u^{n+1} = u/(1 + dt*lambda)`.  Unconditionally stable for any
         !! layer thickness.  Default `.false.` = explicit forward Euler.
      logical :: implicit_fold = .false.
         !! `&ocean_vdiff_nml implicit_top_drag`: fill `lambda_top_u/v`
         !! so the vdiff solver can fold the drag into its `k = nz`
         !! diagonal, and let the driver SKIP the explicit apply.  Default
         !! `.false.` ⇒ the rate fields stay zero.
      real(wp) :: rho0 = 0.0_wp
         !! Boussinesq reference density (kg/m^3) — `eos%rho0` via
         !! `configure_ocean_reference_density`.  Used ONLY to turn the
         !! kinematic drag into the `stress_top` diagnostic; no dynamics
         !! reads it.

      type(scratch_3d_buffer_t) :: du_drag
         !! Top-drag tendency at east faces, shape (nx+1, ny, nz).  Only
         !! layers inside the top boundary layer carry a non-zero value.
      type(scratch_3d_buffer_t) :: dv_drag
         !! Top-drag tendency at north faces, shape (nx, ny+1, nz).

      real(wp), allocatable :: cover_u(:, :)
         !! Face ice-cover mask at east faces, shape (nx+1, ny): the OR of
         !! the two abutting cells' `cover_frac` (see the module
         !! docstring).  STATIC — filled once at configure.
      real(wp), allocatable :: cover_v(:, :)
         !! Face ice-cover mask at north faces, shape (nx, ny+1).
      real(wp), allocatable :: cover_t(:, :)
         !! CELL-CENTRED ice-cover mask, shape (nx, ny) — a configure-time
         !! copy of `metrics%cover_frac`.  Held on the slot (rather than
         !! reaching into `metrics` per step) so the kernel signature
         !! carries exactly the fields it reads.  Used only by the
         !! `stress_top` diagnostic, which is cell-centred.
      real(wp), allocatable :: lambda_top_u(:, :)
         !! Top-layer (k=nz) Rayleigh RATE lambda (1/s) at east faces,
         !! shape (nx+1, ny): `C_d*|U|/h_nz` (quadratic, `|U|` frozen at
         !! u^n) or `r` (linear), already cover- and wet-masked.  Consumed
         !! by `vdiff_apply_momentum` as the `+dt*lambda` add on the
         !! `k = nz` diagonal.  Zero unless `implicit_fold`.
      real(wp), allocatable :: lambda_top_v(:, :)
         !! Top-layer Rayleigh rate at north faces, shape (nx, ny+1).
      real(wp), allocatable :: stress_top(:, :)
         !! Cell-centred magnitude of the top stress (N/m^2),
         !! `rho_0 * |a_drag| * h_tbl` — see the module docstring.  Write
         !! only here; the KPP/EPBL `ustar_shelf` consumer is a later PR.
   contains
      procedure, non_overridable :: init => ocean_top_drag_init
      procedure, non_overridable :: destroy => ocean_top_drag_destroy
      procedure, non_overridable :: enter_data => ocean_top_drag_enter_data
      procedure, non_overridable :: exit_data => ocean_top_drag_exit_data
      procedure, non_overridable :: bytes => ocean_top_drag_bytes
   end type ocean_top_drag_t

contains

   ! ======================================================================
   ! Lifecycle
   ! ======================================================================

   subroutine ocean_top_drag_init(this, grid, nz_ml)
      !! Allocate the slot.  Gated on `enable` (latched before this runs),
      !! so a run without an ice shelf pays five `(1,1)` placeholders and
      !! two `(1,1,1)` scratch buffers.
      class(ocean_top_drag_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nz = 1
      if (this%enable) then
         nx = grid%nx_total
         ny = grid%ny_total
         if (present(nz_ml)) nz = nz_ml
      else
         nx = 1
         ny = 1
      end if

      call this%du_drag%init(nx + 1, ny, nz, "ocean_tdrag_du_drag")
      call this%dv_drag%init(nx, ny + 1, nz, "ocean_tdrag_dv_drag")
      allocate (this%cover_u(nx + 1, ny), source=0.0_wp)
      allocate (this%cover_v(nx, ny + 1), source=0.0_wp)
      allocate (this%cover_t(nx, ny), source=0.0_wp)
      allocate (this%lambda_top_u(nx + 1, ny), source=0.0_wp)
      allocate (this%lambda_top_v(nx, ny + 1), source=0.0_wp)
      allocate (this%stress_top(nx, ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_top_drag_init

   subroutine ocean_top_drag_destroy(this)
      !! Release the slot.  `is_init` is cleared FIRST.
      class(ocean_top_drag_t), intent(inout) :: this
      this%is_init = .false.
      call this%du_drag%destroy()
      call this%dv_drag%destroy()
      if (allocated(this%cover_u)) deallocate (this%cover_u)
      if (allocated(this%cover_v)) deallocate (this%cover_v)
      if (allocated(this%cover_t)) deallocate (this%cover_t)
      if (allocated(this%lambda_top_u)) deallocate (this%lambda_top_u)
      if (allocated(this%lambda_top_v)) deallocate (this%lambda_top_v)
      if (allocated(this%stress_top)) deallocate (this%stress_top)
   end subroutine ocean_top_drag_destroy

   subroutine ocean_top_drag_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! device-attach map base is the heap object, not a polymorphic box.
      class(ocean_top_drag_t), intent(inout) :: this
      select type (this)
      type is (ocean_top_drag_t)
         call ocean_top_drag_enter_data_impl(this)
      end select
   end subroutine ocean_top_drag_enter_data

   subroutine ocean_top_drag_enter_data_impl(this)
      !! `copyin` (not `create`) for the four host-filled 2-D fields —
      !! `cover_u`/`cover_v` are STATIC configure-time geometry and would
      !! otherwise reach the device as allocator leftovers
      !! (`mem:separate`).
      type(ocean_top_drag_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%du_drag)
      call scratch_3d_buffer_enter_data_impl(this%dv_drag)
      !$acc enter data copyin(this%cover_u, this%cover_v, this%cover_t, &
      !$acc                   this%lambda_top_u, this%lambda_top_v, &
      !$acc                   this%stress_top)
   end subroutine ocean_top_drag_enter_data_impl

   subroutine ocean_top_drag_exit_data(this)
      class(ocean_top_drag_t), intent(inout) :: this
      select type (this)
      type is (ocean_top_drag_t)
         call ocean_top_drag_exit_data_impl(this)
      end select
   end subroutine ocean_top_drag_exit_data

   subroutine ocean_top_drag_exit_data_impl(this)
      type(ocean_top_drag_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%du_drag)
      call scratch_3d_buffer_exit_data_impl(this%dv_drag)
      !$acc exit data delete(this%cover_u, this%cover_v, this%cover_t, &
      !$acc                  this%lambda_top_u, this%lambda_top_v, &
      !$acc                  this%stress_top)
   end subroutine ocean_top_drag_exit_data_impl

   ! ======================================================================
   ! Static face cover
   ! ======================================================================

   pure subroutine top_drag_fill_face_cover_impl(cover_u, cover_v, cover_frac, nx, ny)
      !! Project the cell-centred `cover_frac` onto velocity faces with
      !! the **OR** rule (see the module docstring for why OR and not AND
      !! at a calving front).  Host-side, once, at configure.
      !!
      !! Faces outside the drag stencil (`i = 1` and `i > nx` for u,
      !! `j = 1` and `j > ny` for v) are left at zero: they are wall /
      !! ghost faces that carry no prognostic velocity, so an exact zero
      !! there is not an approximation.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: cover_frac(nx, ny)
      real(wp), intent(out) :: cover_u(nx + 1, ny)
      real(wp), intent(out) :: cover_v(nx, ny + 1)

      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         cover_u(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         cover_v(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=2:nx)
         cover_u(i, j) = max(cover_frac(i - 1, j), cover_frac(i, j))
      end do
      do concurrent(j=2:ny, i=1:nx)
         cover_v(i, j) = max(cover_frac(i, j - 1), cover_frac(i, j))
      end do
   end subroutine top_drag_fill_face_cover_impl

   ! ======================================================================
   ! Kernel
   ! ======================================================================

   pure subroutine ocean_top_drag_compute_tendencies(this, ms, dt)
      !! Fill `du_drag` / `dv_drag` with the top-boundary drag
      !! acceleration, `stress_top` with the cell-centred stress
      !! magnitude, and (when `implicit_fold`) `lambda_top_u/v` with the
      !! `k = nz` Rayleigh rate the vdiff diagonal consumes.
      !!
      !! No-op (and no kernel launch) when the slot is disabled — the
      !! arrays are `(1,1[,1])` placeholders then and must not be indexed.
      !! `dt` is read only in the `implicit` branch.
      type(ocean_top_drag_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
         !! Outer-step length (s); only read when `this%implicit`.

      real(wp) :: tbl_min, dt_imp
      integer :: nx, ny, nz

      if (.not. this%enable) return

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml
      tbl_min = this%tbl_thick_min
      if (tbl_min <= 0.0_wp) tbl_min = this%h_min
      dt_imp = 0.0_wp
      if (this%implicit) dt_imp = dt

      call top_drag_tendencies_impl( &
         this%du_drag%data, this%dv_drag%data, &
         this%lambda_top_u, this%lambda_top_v, &
         ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, ms%wet_mask, &
         this%cover_u, this%cover_v, &
         this%variant, this%r_linear, this%c_drag, this%h_min, &
         this%htbl, this%drag_bg_vel, tbl_min, dt_imp, this%implicit_fold, &
         size(ms%u_face_x_layer, 1), size(ms%u_face_x_layer, 2), &
         size(ms%v_face_y_layer, 1), size(ms%v_face_y_layer, 2), &
         nx, ny, nz)

      call top_drag_stress_mag_impl( &
         this%stress_top, ms%u_face_x_layer, ms%v_face_y_layer, &
         ms%h_layer, ms%wet_mask, this%cover_t, &
         this%variant, this%r_linear, this%c_drag, this%h_min, &
         this%htbl, this%drag_bg_vel, tbl_min, this%rho0, &
         size(ms%u_face_x_layer, 1), size(ms%u_face_x_layer, 2), &
         size(ms%v_face_y_layer, 1), size(ms%v_face_y_layer, 2), &
         nx, ny, nz)
   end subroutine ocean_top_drag_compute_tendencies

   pure subroutine top_drag_tendencies_impl(du_drag, dv_drag, lambda_u, lambda_v, &
                                            u_face, v_face, h_layer, wet_mask, &
                                            cover_u, cover_v, &
                                            variant, r, c_d, h_floor, &
                                            htbl, bg_vel, tbl_min, dt_imp, fold, &
                                            nx_u, ny_u, nx_v, ny_v, &
                                            nx, ny, nz)
      !! Flat device kernel: explicit-shape dummies, no derived-type
      !! dereference inside the `do concurrent`.
      !!
      !! One code path covers both modes.  `htbl <= 0` is the LAYER-ONLY
      !! mode: the band is layer `nz` alone and `h_in/h_face == 1`, which
      !! reduces the formulae below to the exact algebraic mirror of
      !! `ocean_bottom_drag_compute_tendencies`' bed-only branch (at the
      !! default `bg_vel = 0`).  `htbl > 0` spreads the stress over the
      !! top `htbl` metres, the mirror of `compute_distributed_drag`.
      !!
      !! Per face, two sequential passes over `k = nz` downward:
      !!   1. band-mean velocity `U_tbl = sum_k u_k*h_in_k / max(sum_k
      !!      h_in_k, tbl_min)` and the band thickness;
      !!   2. per-layer rate and tendency.
      !!
      !! Linear:     `rate_k = r * (h_in_k/h_face_k)`
      !! Quadratic:  `rate_k = C_d * |U_eff| * (h_in_k/h_face_k) / h_tbl`
      !! with `|U_eff| = max(bg_vel, |U_tbl|)`, and the tendency
      !! `-rate_k*u_k/(1 + dt_imp*rate_k)` — `dt_imp = 0` gives the
      !! explicit form bit-identically.
      !!
      !! `h_in_k/h_face_k` is bounded by 1 by construction (`h_in_k =
      !! min(h_face_k, ...)`) and the loop exits on `h_face_k <= 0`, so
      !! the ratio needs no epsilon.
      integer, intent(in) :: nx_u, ny_u, nx_v, ny_v, nx, ny, nz, variant
      real(wp), intent(in) :: r, c_d, h_floor, htbl, bg_vel, tbl_min, dt_imp
      logical, intent(in) :: fold
      real(wp), intent(in) :: u_face(nx_u, ny_u, nz), v_face(nx_v, ny_v, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: cover_u(nx + 1, ny), cover_v(nx, ny + 1)
      real(wp), intent(out) :: du_drag(nx_u, ny_u, nz), dv_drag(nx_v, ny_v, nz)
      real(wp), intent(out) :: lambda_u(nx + 1, ny), lambda_v(nx, ny + 1)

      integer :: i, j, k
      logical :: layer_only, quad
      real(wp) :: cumul_h, h_face_k, h_in, mask_face, frac
      real(wp) :: h_in_total, u_int, v_int, u_tbl, v_tbl, u_at_v, v_at_u
      real(wp) :: abs_u_eff, h_eff, rate

      layer_only = (htbl <= 0.0_wp)
      quad = (variant == TDRAG_QUADRATIC)

      ! ---- Zero every level + both rate fields first ----
      do concurrent(k=1:nz, j=1:ny_u, i=1:nx_u)
         du_drag(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny_v, i=1:nx_v)
         dv_drag(i, j, k) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         lambda_u(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         lambda_v(i, j) = 0.0_wp
      end do

      if (variant == TDRAG_LINEAR .and. r <= 0.0_wp) return
      if (quad .and. c_d <= 0.0_wp) return
      if (.not. (quad .or. variant == TDRAG_LINEAR)) return

      ! ---- East (u) faces ----
      do concurrent(j=1:ny, i=2:nx) &
         local(k, cumul_h, h_face_k, h_in, mask_face, frac, &
               h_in_total, u_int, v_int, u_tbl, v_tbl, v_at_u, &
               abs_u_eff, h_eff, rate)
         mask_face = min(wet_mask(i - 1, j), wet_mask(i, j))*cover_u(i, j)
         ! Pass 1: band mean.
         cumul_h = 0.0_wp
         u_int = 0.0_wp
         v_int = 0.0_wp
         h_in_total = 0.0_wp
         do k = nz, 1, -1
            if (layer_only .and. k < nz) exit
            if ((.not. layer_only) .and. cumul_h >= htbl) exit
            h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            if (layer_only) then
               h_in = h_face_k
            else
               h_in = max(0.0_wp, min(h_face_k, htbl - cumul_h))
            end if
            u_int = u_int + u_face(i, j, k)*h_in
            v_at_u = 0.25_wp*( &
                     v_face(i - 1, j, k) + v_face(i, j, k) + &
                     v_face(i - 1, j + 1, k) + v_face(i, j + 1, k))
            v_int = v_int + v_at_u*h_in
            h_in_total = h_in_total + h_in
            cumul_h = cumul_h + h_face_k
         end do
         h_eff = max(h_in_total, tbl_min)
         h_eff = max(h_eff, h_floor)
         u_tbl = u_int/h_eff
         v_tbl = v_int/h_eff
         abs_u_eff = max(bg_vel, sqrt(u_tbl*u_tbl + v_tbl*v_tbl))
         ! Pass 2: per-layer rate + tendency.
         cumul_h = 0.0_wp
         do k = nz, 1, -1
            if (layer_only .and. k < nz) exit
            if ((.not. layer_only) .and. cumul_h >= htbl) exit
            h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            if (layer_only) then
               h_in = h_face_k
            else
               h_in = max(0.0_wp, min(h_face_k, htbl - cumul_h))
            end if
            frac = h_in/h_face_k
            if (quad) then
               rate = mask_face*c_d*abs_u_eff*frac/h_eff
            else
               rate = mask_face*r*frac
            end if
            du_drag(i, j, k) = -rate*u_face(i, j, k)/(1.0_wp + dt_imp*rate)
            if (fold .and. k == nz) lambda_u(i, j) = rate
            cumul_h = cumul_h + h_face_k
         end do
      end do

      ! ---- North (v) faces ----
      do concurrent(j=2:ny, i=1:nx) &
         local(k, cumul_h, h_face_k, h_in, mask_face, frac, &
               h_in_total, u_int, v_int, u_tbl, v_tbl, u_at_v, &
               abs_u_eff, h_eff, rate)
         mask_face = min(wet_mask(i, j - 1), wet_mask(i, j))*cover_v(i, j)
         cumul_h = 0.0_wp
         u_int = 0.0_wp
         v_int = 0.0_wp
         h_in_total = 0.0_wp
         do k = nz, 1, -1
            if (layer_only .and. k < nz) exit
            if ((.not. layer_only) .and. cumul_h >= htbl) exit
            h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            if (layer_only) then
               h_in = h_face_k
            else
               h_in = max(0.0_wp, min(h_face_k, htbl - cumul_h))
            end if
            v_int = v_int + v_face(i, j, k)*h_in
            u_at_v = 0.25_wp*( &
                     u_face(i, j - 1, k) + u_face(i + 1, j - 1, k) + &
                     u_face(i, j, k) + u_face(i + 1, j, k))
            u_int = u_int + u_at_v*h_in
            h_in_total = h_in_total + h_in
            cumul_h = cumul_h + h_face_k
         end do
         h_eff = max(h_in_total, tbl_min)
         h_eff = max(h_eff, h_floor)
         u_tbl = u_int/h_eff
         v_tbl = v_int/h_eff
         abs_u_eff = max(bg_vel, sqrt(u_tbl*u_tbl + v_tbl*v_tbl))
         cumul_h = 0.0_wp
         do k = nz, 1, -1
            if (layer_only .and. k < nz) exit
            if ((.not. layer_only) .and. cumul_h >= htbl) exit
            h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            if (layer_only) then
               h_in = h_face_k
            else
               h_in = max(0.0_wp, min(h_face_k, htbl - cumul_h))
            end if
            frac = h_in/h_face_k
            if (quad) then
               rate = mask_face*c_d*abs_u_eff*frac/h_eff
            else
               rate = mask_face*r*frac
            end if
            dv_drag(i, j, k) = -rate*v_face(i, j, k)/(1.0_wp + dt_imp*rate)
            if (fold .and. k == nz) lambda_v(i, j) = rate
            cumul_h = cumul_h + h_face_k
         end do
      end do
   end subroutine top_drag_tendencies_impl

   pure subroutine top_drag_stress_mag_impl(stress_top, u_face, v_face, &
                                            h_layer, wet_mask, cover_frac, &
                                            variant, r, c_d, h_floor, &
                                            htbl, bg_vel, tbl_min, rho0, &
                                            nx_u, ny_u, nx_v, ny_v, nx, ny, nz)
      !! Cell-centred magnitude of the top stress (N/m^2), for the
      !! later `ustar_shelf` consumer:
      !!
      !!   quadratic  `|tau| = rho_0 * C_d * |U_tbl_eff|^2`
      !!   linear     `|tau| = rho_0 * r * |U_tbl| * h_tbl`
      !!
      !! Both are `rho_0 * (drag acceleration) * (band thickness)`, so the
      !! two forms are one definition, and `sqrt(|tau|/rho_0)` is the
      !! friction velocity either way.  Cell-centred velocities come from
      !! the ordinary 2-point face averages; no cover interpolation is
      !! needed because `cover_frac` IS cell-centred here.
      integer, intent(in) :: nx_u, ny_u, nx_v, ny_v, nx, ny, nz, variant
      real(wp), intent(in) :: r, c_d, h_floor, htbl, bg_vel, tbl_min, rho0
      real(wp), intent(in) :: u_face(nx_u, ny_u, nz), v_face(nx_v, ny_v, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny), cover_frac(nx, ny)
      real(wp), intent(out) :: stress_top(nx, ny)

      integer :: i, j, k
      logical :: layer_only, quad
      real(wp) :: cumul_h, h_k, h_in, u_int, v_int, h_in_total
      real(wp) :: u_c, v_c, u_tbl, v_tbl, h_eff, abs_u_eff, spd

      layer_only = (htbl <= 0.0_wp)
      quad = (variant == TDRAG_QUADRATIC)

      do concurrent(j=1:ny, i=1:nx) &
         local(k, cumul_h, h_k, h_in, u_int, v_int, h_in_total, &
               u_c, v_c, u_tbl, v_tbl, h_eff, abs_u_eff, spd)
         cumul_h = 0.0_wp
         u_int = 0.0_wp
         v_int = 0.0_wp
         h_in_total = 0.0_wp
         do k = nz, 1, -1
            if (layer_only .and. k < nz) exit
            if ((.not. layer_only) .and. cumul_h >= htbl) exit
            h_k = h_layer(i, j, k)
            if (h_k <= 0.0_wp) exit
            if (layer_only) then
               h_in = h_k
            else
               h_in = max(0.0_wp, min(h_k, htbl - cumul_h))
            end if
            u_c = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
            v_c = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
            u_int = u_int + u_c*h_in
            v_int = v_int + v_c*h_in
            h_in_total = h_in_total + h_in
            cumul_h = cumul_h + h_k
         end do
         h_eff = max(max(h_in_total, tbl_min), h_floor)
         u_tbl = u_int/h_eff
         v_tbl = v_int/h_eff
         spd = sqrt(u_tbl*u_tbl + v_tbl*v_tbl)
         abs_u_eff = max(bg_vel, spd)
         if (quad) then
            stress_top(i, j) = rho0*c_d*abs_u_eff*abs_u_eff* &
                               wet_mask(i, j)*cover_frac(i, j)
         else
            stress_top(i, j) = rho0*r*spd*h_eff* &
                               wet_mask(i, j)*cover_frac(i, j)
         end if
      end do
   end subroutine top_drag_stress_mag_impl

   subroutine ocean_top_drag_apply_tendencies(this, ms, dt, no_wait)
      !! `u += dt*du_drag`, `v += dt*dv_drag` over the whole face array —
      !! layers outside the top boundary layer carry an exact zero.
      !!
      !! No-op when the slot is disabled (the buffers are `(1,1,1)`
      !! placeholders then).  `no_wait` semantics mirror
      !! `ocean_bottom_drag_apply_tendencies`: `.true.` runs the apply on
      !! OpenACC queue 1 and returns WITHOUT syncing, so the batched
      !! velocity-apply chain waits once.  Not `pure` (async/wait
      !! directives).
      type(ocean_top_drag_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: no_wait
      integer :: i, j, k, nx_face, ny_uface, nx_vface, ny_face, nz
      logical :: lwait

      if (.not. this%enable) return

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
                                      dt*this%du_drag%data(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + &
                                      dt*this%dv_drag%data(i, j, k)
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_top_drag_apply_tendencies

   ! ======================================================================
   ! Parsing + accounting
   ! ======================================================================

   pure function parse_tdrag_variant(name) result(code)
      !! Translate a namelist string into a `TDRAG_*` code.  An
      !! unrecognised string returns `TDRAG_INVALID` — a typo must abort
      !! rather than silently select one of two laws with different
      !! coefficient dimensions.  Accepts the same spellings as
      !! `parse_bdrag_variant`, deliberately: the two groups mirror.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("quadratic", "QUADRATIC", "cd", "CD")
         code = TDRAG_QUADRATIC
      case ("linear", "LINEAR", "rayleigh", "RAYLEIGH")
         code = TDRAG_LINEAR
      case default
         code = TDRAG_INVALID
      end select
   end function parse_tdrag_variant

   pure function tdrag_variant_is_implemented(code) result(ok)
      !! `.true.` only for a top-drag variant with a real kernel.  The
      !! single gate `validate_config` consumes.
      integer, intent(in) :: code
      logical :: ok
      ok = (code == TDRAG_LINEAR) .or. (code == TDRAG_QUADRATIC)
   end function tdrag_variant_is_implemented

   pure function ocean_top_drag_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the top-drag slot (0 when
      !! unallocated).  One `arr_bytes` term per array — add a term here
      !! when a new allocatable joins the type.
      class(ocean_top_drag_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%du_drag%bytes() &
               + this%dv_drag%bytes() &
               + arr_bytes(this%cover_u) &
               + arr_bytes(this%cover_v) &
               + arr_bytes(this%cover_t) &
               + arr_bytes(this%lambda_top_u) &
               + arr_bytes(this%lambda_top_v) &
               + arr_bytes(this%stress_top)
   end function ocean_top_drag_bytes

end module rdb_ocean_top_drag
