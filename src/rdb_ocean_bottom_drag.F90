!! Bottom drag for the ocean dynamical core.  Applies a momentum
!! sink at the bottom-most layer (`k = 1` under the ROMS-style
!! convention this codebase uses everywhere — bed at k=1, surface
!! at k=nz).  Two variants supported via a tag:
!!
!!   BDRAG_LINEAR    — Rayleigh-style: du/dt = -r * u, where `r`
!!                     has units of 1/s.
!!   BDRAG_QUADRATIC — log-layer drag: du/dt = -C_d * |U_bot| *
!!                     u / h_bot.  Closure for the unresolved
!!                     bottom-boundary-layer drag.  C_d ~ 2.5e-3
!!                     is the standard MOM6 / ROMS default.
!!
!! Only the bottom layer feels the drag — the kernel writes into
!! its own tendency buffer (`du_drag`, `dv_drag`) and the apply
!! step adds `dt * tendency` to `u_face_x_layer` / `v_face_y_layer`,
!! matching the additive-tendency pattern of Coriolis, PGF, and hvisc.
module rdb_ocean_bottom_drag
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_bottom_drag_t
   public :: ocean_bottom_drag_compute_tendencies
   public :: ocean_bottom_drag_apply_tendencies
   public :: ocean_channel_drag_compute_tendencies
   public :: ocean_channel_drag_apply_tendencies
   public :: parse_bdrag_variant
   public :: bdrag_variant_is_implemented

   real(wp), parameter :: SIDE_H_VANISH = 1.5e-4_wp
      !! Layer-thickness threshold below which a cross-stream neighbour
      !! layer counts as "blocked by sloping bathymetry" for the
      !! channel-drag perimeter fraction (matches the `H_VANISHED`
      !! dynamic-vanish role used elsewhere in the ocean core).

   integer, parameter, public :: BDRAG_LINEAR = 1
      !! Linear Rayleigh drag.  Use for analytic / regression tests
      !! where energy decay rate must be predictable.
   integer, parameter, public :: BDRAG_QUADRATIC = 2
      !! Log-layer / quadratic drag.  Production default.
   integer, parameter, public :: BDRAG_INVALID = -1
      !! Sentinel returned by `parse_bdrag_variant` for an unrecognised
      !! string (PR-6 fail-loud).  Linear (`τ = ρ·r·u`) and quadratic
      !! (`τ = ρ·C_d·|u|·u`) have different coefficient dimensions and
      !! energy-decay laws, so a typo must abort rather than silently
      !! pick one.

   type :: ocean_bottom_drag_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      integer :: variant = BDRAG_QUADRATIC
         !! Active drag variant.
      real(wp) :: r_linear = 0.0_wp
         !! Linear Rayleigh coefficient (1/s).  Zero disables the
         !! linear branch even if the variant tag selects it.
      real(wp) :: c_drag = 0.0_wp
         !! Quadratic-drag coefficient (dimensionless).  Zero
         !! disables the quadratic branch.  MOM6 / ROMS default
         !! is 2.5e-3.
      real(wp) :: h_min = 1.0e-3_wp
         !! Floor on the bottom-layer thickness inside the
         !! `u / h_bot` division — keeps the kernel finite when
         !! the bottom layer pinches out under ZSTAR_FULL.
      real(wp) :: hbbl = 0.0_wp
         !! Bottom-boundary-layer thickness (m) over which the drag
         !! is distributed (MOM6 HBBL).  When zero (default), the
         !! drag is applied to the bed layer only — historical
         !! behaviour, bit-identical to prior production.  Positive
         !! values activate the distributed-drag form: the stress
         !! is spread across the bottom `hbbl` metres so very thin
         !! bed layers don't see all the damping at once.
      real(wp) :: drag_bg_vel = 0.0_wp
         !! Background velocity floor (m/s) for the quadratic and
         !! distributed forms (MOM6 DRAG_BG_VEL).  The effective
         !! bottom speed used in the stress formula is
         !! `max(drag_bg_vel, |u_bbl|)` so laminar-bottom cells
         !! still see damping.  MOM6 production default 0.1 m/s.
         !! Zero (default) disables the floor.
      real(wp) :: bbl_thick_min = 0.0_wp
         !! Minimum effective BBL thickness (m) used in the
         !! `stress / h_bbl` denominator (MOM6 BBL_THICK_MIN).
         !! Guards against `hbbl - sum(h_bottom_layers)` going to
         !! zero when the bottom layers themselves are very thin.
         !! Zero (default) leaves no floor (falls back to `h_min`).
      real(wp) :: bed_factor = 1.0_wp
         !! Multiplier on the bed-layer (k=1) drag tendency only —
         !! layers `k >= 2` are unchanged.  Lets the bed get
         !! `bed_factor · r` while the surface keeps the nominal `r`
         !! (and via HBBL distribution typically gets none).  Default
         !! `1.0` keeps the historical drag bit-identical.  Driver
         !! writes from `cfg%ocean%bdrag%bed_factor` at init.  See
         !! `compute_distributed_drag` for the kernel-side application.
      logical :: implicit = .false.
         !! Backward-Euler (implicit) bottom drag when `.true.`:
         !! `u^{n+1} = u/(1 + dt·λ)`, unconditionally stable for any
         !! bottom-layer thickness (matches MOM6's implicit bottom-BC
         !! drag).  Default `.false.` = explicit forward-Euler tendency,
         !! bit-identical to the historical path but conditionally
         !! unstable on thin shelf bottom layers (`dt·λ > 1`).  Driver
         !! writes from `cfg%ocean%bdrag%implicit` at init.

      logical :: channel_drag = .false.
         !! Enable the per-layer lateral side-wall (channel) Rayleigh
         !! drag (MOM6 CHANNEL_DRAG analogue).  Default `.false.` ⇒ the
         !! channel-drag kernels are no-ops (bit-identical).  Driver
         !! writes from `cfg%ocean%bdrag%channel_drag` at init.
      real(wp) :: cdrag_side = 0.0_wp
         !! Side-wall drag coefficient (dimensionless).  Zero (default)
         !! ⇒ zero side-drag rate even when `channel_drag = .true.`.

      type(scratch_3d_buffer_t) :: du_drag
         !! Drag tendency at east faces, shape (nx+1, ny, nz).
         !! Only k=1 (bed) carries a non-zero value; k>=2 stays at
         !! 0 because the kernel only writes there.
      type(scratch_3d_buffer_t) :: dv_drag
         !! Drag tendency at north faces, shape (nx, ny+1, nz).

      logical :: implicit_fold = .false.
         !! When `.true.`, also fill the bed-layer (k=1) Rayleigh RATE
         !! fields `lambda_bot_u/v` (1/s) in `compute_tendencies` so the
         !! vdiff solver can fold the drag into its tridiagonal diagonal
         !! (`&ocean_vdiff_nml implicit_drag`).  The explicit `du_drag`
         !! tendency is then NOT applied by the driver (gated off) to avoid
         !! double-counting.  Default `.false.` ⇒ the rate fields stay zero
         !! and the explicit path is unchanged (bit-identical).  Driver sets
         !! this from `cfg%ocean%vdiff%implicit_drag` at configure.

      real(wp), allocatable :: lambda_bot_u(:, :)
         !! Bed-layer (k=1) bottom-drag Rayleigh RATE λ (1/s) at east
         !! faces, shape (nx+1, ny).  `λ = c_d·|U_bbl|/h_1` (quadratic) or
         !! `r` (linear), `|U|` frozen at uⁿ.  Consumed by `vdiff_apply_
         !! momentum` as the `+dt·λ` diagonal add when `implicit_drag`.
         !! Zero unless `implicit_fold = .true.`.
      real(wp), allocatable :: lambda_bot_v(:, :)
         !! Bed-layer bottom-drag Rayleigh rate λ (1/s) at north faces,
         !! shape (nx, ny+1).

      type(scratch_3d_buffer_t) :: lambda_side_u
         !! Per-layer side-drag Rayleigh RATE (1/s) at east faces,
         !! shape (nx+1, ny, nz).  Stored (not a tendency) so the apply
         !! can use the implicit form `u/(1+dt*lambda)`.  Zero where no
         !! lateral perimeter is blocked (all-wet / flat-bottom).
      type(scratch_3d_buffer_t) :: lambda_side_v
         !! Per-layer side-drag Rayleigh rate (1/s) at north faces,
         !! shape (nx, ny+1, nz).
   contains
      procedure, non_overridable :: init => ocean_bdrag_init
      procedure, non_overridable :: destroy => ocean_bdrag_destroy
      procedure, non_overridable :: enter_data => ocean_bdrag_enter_data
      procedure, non_overridable :: exit_data => ocean_bdrag_exit_data
      procedure, non_overridable :: bytes => ocean_bottom_drag_bytes
   end type ocean_bottom_drag_t

contains

   subroutine ocean_bdrag_init(this, grid, nz_ml)
      class(ocean_bottom_drag_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      call this%du_drag%init(nx + 1, ny, nz, "ocean_bdrag_du_drag")
      call this%dv_drag%init(nx, ny + 1, nz, "ocean_bdrag_dv_drag")
      allocate (this%lambda_bot_u(nx + 1, ny), source=0.0_wp)
      allocate (this%lambda_bot_v(nx, ny + 1), source=0.0_wp)
      call this%lambda_side_u%init(nx + 1, ny, nz, "ocean_bdrag_lambda_side_u")
      call this%lambda_side_v%init(nx, ny + 1, nz, "ocean_bdrag_lambda_side_v")
      this%is_init = .true.
   end subroutine ocean_bdrag_init

   subroutine ocean_bdrag_destroy(this)
      class(ocean_bottom_drag_t), intent(inout) :: this
      this%is_init = .false.
      call this%du_drag%destroy()
      call this%dv_drag%destroy()
      if (allocated(this%lambda_bot_u)) deallocate (this%lambda_bot_u)
      if (allocated(this%lambda_bot_v)) deallocate (this%lambda_bot_v)
      call this%lambda_side_u%destroy()
      call this%lambda_side_v%destroy()
   end subroutine ocean_bdrag_destroy

   subroutine ocean_bdrag_enter_data(this)
      class(ocean_bottom_drag_t), intent(inout) :: this
      select type (this)
      type is (ocean_bottom_drag_t)
         call ocean_bdrag_enter_data_impl(this)
      end select
   end subroutine ocean_bdrag_enter_data

   subroutine ocean_bdrag_enter_data_impl(this)
      type(ocean_bottom_drag_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%du_drag)
      call scratch_3d_buffer_enter_data_impl(this%dv_drag)
      !$acc enter data copyin(this%lambda_bot_u, this%lambda_bot_v)
      call scratch_3d_buffer_enter_data_impl(this%lambda_side_u)
      call scratch_3d_buffer_enter_data_impl(this%lambda_side_v)
   end subroutine ocean_bdrag_enter_data_impl

   subroutine ocean_bdrag_exit_data(this)
      class(ocean_bottom_drag_t), intent(inout) :: this
      select type (this)
      type is (ocean_bottom_drag_t)
         call ocean_bdrag_exit_data_impl(this)
      end select
   end subroutine ocean_bdrag_exit_data

   subroutine ocean_bdrag_exit_data_impl(this)
      type(ocean_bottom_drag_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%du_drag)
      call scratch_3d_buffer_exit_data_impl(this%dv_drag)
      !$acc exit data delete(this%lambda_bot_u, this%lambda_bot_v)
      call scratch_3d_buffer_exit_data_impl(this%lambda_side_u)
      call scratch_3d_buffer_exit_data_impl(this%lambda_side_v)
   end subroutine ocean_bdrag_exit_data_impl

   pure subroutine ocean_bottom_drag_compute_tendencies(grid, this, ms, dt)
      !! Fill `du_drag` / `dv_drag` with the bottom-layer drag
      !! acceleration.  Layers k >= 2 get zero in the bed-only mode
      !! (default).  When `hbbl > 0` the stress is distributed
      !! across the bottom-most `hbbl` metres — every layer with
      !! `cumulative_depth_from_bed_top ≤ hbbl` gets a proportional
      !! share of the drag tendency.
      !!
      !! `implicit=.true.` (knob, default .false.) makes the drag
      !! BACKWARD-EULER in the dragged velocity: a per-face/-layer rate
      !! `λ` (= `c_d·|U|/h` etc.) gives `u^{n+1} = u/(1+dt·λ)`, formed by
      !! the tendency `-λ·u/(1+dt·λ)` so the standalone `u += dt·du_drag`
      !! apply reproduces it exactly.  Unconditionally stable for ANY h
      !! (matches MOM6's implicit bottom-BC drag); the explicit form
      !! (default) is conditionally unstable on thin bottom layers
      !! (`λ·dt > 1`).  `dt` is unused in the explicit branch.
      !!
      !! Wall faces (i=1, i=nx+1 for u; j=1, j=ny+1 for v) get a
      !! zero tendency — they don't move under drag because they
      !! don't move at all under any kernel here.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bottom_drag_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
         !! Outer-step length (s); only read when `this%implicit`.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: r, c_d, h_floor, u_bot, v_bot, h_face, u_at_v, v_at_u
      real(wp) :: speed_at_u, speed_at_v, lam, dt_imp
      real(wp) :: hbbl, bg_vel, bbl_min
      logical :: implicit_drag, fold

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      r = this%r_linear
      c_d = this%c_drag
      h_floor = this%h_min
      hbbl = this%hbbl
      bg_vel = this%drag_bg_vel
      bbl_min = this%bbl_thick_min
      if (bbl_min <= 0.0_wp) bbl_min = h_floor
      implicit_drag = this%implicit
      dt_imp = merge(dt, 0.0_wp, implicit_drag)  ! 0 ⇒ explicit (bit-identical)
      fold = this%implicit_fold

      ! ---- Zero every level first; only k=1 will get filled ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         this%du_drag%data(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         this%dv_drag%data(i, j, k) = 0.0_wp
      end do
      ! Always re-zero the implicit-fold bed-rate fields so a face that
      ! left the wet stencil between steps doesn't carry a stale rate.
      do concurrent(j=1:ny, i=1:nx + 1)
         this%lambda_bot_u(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         this%lambda_bot_v(i, j) = 0.0_wp
      end do

      ! ---- Bed-only implicit-fold rate fields (k=1) ----
      ! λ is the Rayleigh RATE the vdiff diagonal consumes (`+dt·λ`): r
      ! (linear) or c_d·|U_bbl|/h_1 (quadratic, |U| frozen at uⁿ).  Filled
      ! here from the SAME stencil that forms du_drag so there is one drag
      ! algebra.  Configure forbids hbbl>0 + fold, so the bed-only 2D field
      ! is sufficient.  Zero where the face touches land (mask).
      if (fold) then
         if (this%variant == BDRAG_LINEAR .and. r > 0.0_wp) then
            do concurrent(j=1:ny, i=2:nx)
               this%lambda_bot_u(i, j) = min(ms%wet_mask(i - 1, j), ms%wet_mask(i, j))*r
            end do
            do concurrent(j=2:ny, i=1:nx)
               this%lambda_bot_v(i, j) = min(ms%wet_mask(i, j - 1), ms%wet_mask(i, j))*r
            end do
         else if (this%variant == BDRAG_QUADRATIC .and. c_d > 0.0_wp) then
            do concurrent(j=1:ny, i=2:nx) local(u_bot, v_at_u, h_face, speed_at_u)
               u_bot = ms%u_face_x_layer(i, j, 1)
               v_at_u = 0.25_wp*( &
                        ms%v_face_y_layer(i - 1, j, 1) + ms%v_face_y_layer(i, j, 1) + &
                        ms%v_face_y_layer(i - 1, j + 1, 1) + ms%v_face_y_layer(i, j + 1, 1))
               h_face = max(0.5_wp*(ms%h_layer(i - 1, j, 1) + ms%h_layer(i, j, 1)), h_floor)
               speed_at_u = sqrt(u_bot*u_bot + v_at_u*v_at_u)
               this%lambda_bot_u(i, j) = min(ms%wet_mask(i - 1, j), ms%wet_mask(i, j))* &
                                         c_d*speed_at_u/h_face
            end do
            do concurrent(j=2:ny, i=1:nx) local(v_bot, u_at_v, h_face, speed_at_v)
               v_bot = ms%v_face_y_layer(i, j, 1)
               u_at_v = 0.25_wp*( &
                        ms%u_face_x_layer(i, j - 1, 1) + ms%u_face_x_layer(i + 1, j - 1, 1) + &
                        ms%u_face_x_layer(i, j, 1) + ms%u_face_x_layer(i + 1, j, 1))
               h_face = max(0.5_wp*(ms%h_layer(i, j - 1, 1) + ms%h_layer(i, j, 1)), h_floor)
               speed_at_v = sqrt(v_bot*v_bot + u_at_v*u_at_v)
               this%lambda_bot_v(i, j) = min(ms%wet_mask(i, j - 1), ms%wet_mask(i, j))* &
                                         c_d*speed_at_v/h_face
            end do
         end if
      end if

      ! ---- HBBL-distributed branch (MOM6 LINEAR_DRAG / BBL_THICK_MIN) ----
      if (hbbl > 0.0_wp) then
         call compute_distributed_drag(this%du_drag%data, this%dv_drag%data, &
                                       ms%u_face_x_layer, ms%v_face_y_layer, &
                                       ms%h_layer, ms%wet_mask, &
                                       this%variant, r, c_d, hbbl, bg_vel, bbl_min, &
                                       this%bed_factor, dt_imp, &
                                       size(ms%u_face_x_layer, 1), size(ms%u_face_x_layer, 2), &
                                       size(ms%v_face_y_layer, 1), size(ms%v_face_y_layer, 2), &
                                       nx, ny, nz)
         return
      end if

      ! Face wet-mask: drag only fires at faces between two ocean cells.
      ! `min(wet_left, wet_right)` zeros the tendency at any face that
      ! touches land.  All-1.0 mask (analytical tests) is a no-op.

      ! ---- Linear branch ----
      ! Implicit (dt_imp>0): u^{n+1}=u/(1+dt·r) via tendency -r·u/(1+dt·r).
      if (this%variant == BDRAG_LINEAR .and. r > 0.0_wp) then
         do concurrent(j=1:ny, i=2:nx) local(u_bot)
            u_bot = ms%u_face_x_layer(i, j, 1)
            this%du_drag%data(i, j, 1) = &
               min(ms%wet_mask(i - 1, j), ms%wet_mask(i, j))*(-r*u_bot/(1.0_wp + dt_imp*r))
         end do
         do concurrent(j=2:ny, i=1:nx) local(v_bot)
            v_bot = ms%v_face_y_layer(i, j, 1)
            this%dv_drag%data(i, j, 1) = &
               min(ms%wet_mask(i, j - 1), ms%wet_mask(i, j))*(-r*v_bot/(1.0_wp + dt_imp*r))
         end do
         return
      end if

      ! ---- Quadratic branch ----
      if (this%variant == BDRAG_QUADRATIC .and. c_d > 0.0_wp) then
         ! du/dt = -C_d * |U| * u / h_bot, where |U| = sqrt(u^2 + v^2)
         ! evaluated at the same face.  For u-faces we average v from
         ! the four surrounding v-faces (standard C-grid stencil); for
         ! v-faces we average u from the four surrounding u-faces.
         do concurrent(j=1:ny, i=2:nx) &
            local(u_bot, v_at_u, h_face, speed_at_u)
            u_bot = ms%u_face_x_layer(i, j, 1)
            v_at_u = 0.25_wp*( &
                     ms%v_face_y_layer(i - 1, j, 1) + ms%v_face_y_layer(i, j, 1) + &
                     ms%v_face_y_layer(i - 1, j + 1, 1) + ms%v_face_y_layer(i, j + 1, 1))
            h_face = 0.5_wp*(ms%h_layer(i - 1, j, 1) + ms%h_layer(i, j, 1))
            h_face = max(h_face, h_floor)
            speed_at_u = sqrt(u_bot*u_bot + v_at_u*v_at_u)
            ! Implicit: denom h_face → h_face + dt·c_d·|U| ⇒ u/(1+dt·c_d·|U|/h).
            this%du_drag%data(i, j, 1) = &
               min(ms%wet_mask(i - 1, j), ms%wet_mask(i, j))* &
               (-c_d*speed_at_u*u_bot/(h_face + dt_imp*c_d*speed_at_u))
         end do
         do concurrent(j=2:ny, i=1:nx) &
            local(v_bot, u_at_v, h_face, speed_at_v)
            v_bot = ms%v_face_y_layer(i, j, 1)
            u_at_v = 0.25_wp*( &
                     ms%u_face_x_layer(i, j - 1, 1) + ms%u_face_x_layer(i + 1, j - 1, 1) + &
                     ms%u_face_x_layer(i, j, 1) + ms%u_face_x_layer(i + 1, j, 1))
            h_face = 0.5_wp*(ms%h_layer(i, j - 1, 1) + ms%h_layer(i, j, 1))
            h_face = max(h_face, h_floor)
            speed_at_v = sqrt(v_bot*v_bot + u_at_v*u_at_v)
            this%dv_drag%data(i, j, 1) = &
               min(ms%wet_mask(i, j - 1), ms%wet_mask(i, j))* &
               (-c_d*speed_at_v*v_bot/(h_face + dt_imp*c_d*speed_at_v))
         end do
      end if
   end subroutine ocean_bottom_drag_compute_tendencies

   pure subroutine compute_distributed_drag(du_drag, dv_drag, u_face, v_face, &
                                            h_layer, wet_mask, &
                                            variant, r, c_d, hbbl, bg_vel, bbl_min, &
                                            bed_factor, dt_imp, &
                                            nx_u, ny_u, nx_v, ny_v, nx, ny, nz)
      !! HBBL-distributed bottom drag.  Mirrors MOM6's LINEAR_DRAG
      !! and quadratic-with-HBBL formulations: the drag stress is
      !! spread across the bottom `hbbl` metres rather than dumped
      !! into the bed-most layer.
      !!
      !! Linear branch (`variant == BDRAG_LINEAR, r > 0`):
      !!   du_k/dt = -r · u_k · (h_in_bbl_k / h_face_k) · f_k
      !! Each layer that overlaps the BBL band gets a damping rate
      !! proportional to its fractional BBL coverage.  Integrating
      !! over k recovers the bulk MOM6 result `r · U_bbl`.
      !!
      !! Quadratic branch (`variant == BDRAG_QUADRATIC, c_d > 0`):
      !!   1. First pass per face: compute BBL-mean velocity
      !!      `U_bbl = Σ_k u_k · h_in_bbl_k / max(Σ_k h_in_bbl_k, bbl_min)`.
      !!   2. Effective speed `|U_eff| = max(bg_vel, |U_bbl|)`.
      !!   3. Per-layer apply:
      !!      du_k/dt = -c_d · |U_eff| · u_k · (h_in_bbl_k / h_face_k) · f_k /
      !!                 max(h_in_bbl_total, bbl_min)
      !!
      !! `f_k = bed_factor` for `k = 1` and `f_k = 1` otherwise — lets
      !! the bed layer carry stronger drag than the rest of the BBL
      !! while preserving HBBL-distribution shape for layers k>=2.
      !! Default `bed_factor = 1.0` ⇒ `f_k ≡ 1` ⇒ bit-identical to
      !! the pre-knob path.
      !!
      !! Flat-impl: explicit-shape dummies, no derived-type derefs in
      !! the device kernels.
      integer, intent(in) :: nx_u, ny_u, nx_v, ny_v, nx, ny, nz, variant
      real(wp), intent(in) :: r, c_d, hbbl, bg_vel, bbl_min, bed_factor
      real(wp), intent(in) :: dt_imp
         !! Implicit timestep: dt for backward-Euler drag, 0 for explicit.
      real(wp), intent(in)    :: u_face(nx_u, ny_u, nz), v_face(nx_v, ny_v, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(inout) :: du_drag(nx_u, ny_u, nz), dv_drag(nx_v, ny_v, nz)

      integer :: i, j, k
      real(wp) :: cumul_h, h_face_k, h_in_bbl, mask_face, f_k
      real(wp) :: h_in_bbl_total, u_bbl_int, v_bbl_int, u_bbl, v_bbl
      real(wp) :: u_at_v, v_at_u, abs_U_eff, h_eff_denom

      ! Linear: per-layer independent — single pass walks k=1 upward,
      ! accumulating BBL thickness and applying the per-layer drag.
      ! Bed layer (k=1) scaled by `bed_factor`; layers k>=2 unchanged.
      if (variant == BDRAG_LINEAR .and. r > 0.0_wp) then
         do concurrent(j=1:ny, i=2:nx) &
            local(k, cumul_h, h_face_k, h_in_bbl, mask_face, f_k)
            mask_face = min(wet_mask(i - 1, j), wet_mask(i, j))
            cumul_h = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
               if (h_face_k > 0.0_wp) then
                  h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
                  if (k == 1) then
                     f_k = bed_factor
                  else
                     f_k = 1.0_wp
                  end if
                  du_drag(i, j, k) = mask_face*(-r*u_face(i, j, k)*h_in_bbl/h_face_k)*f_k &
                                     /(1.0_wp + dt_imp*r*(h_in_bbl/h_face_k)*f_k)
                  cumul_h = cumul_h + h_face_k
               else
                  exit
               end if
            end do
         end do
         do concurrent(j=2:ny, i=1:nx) &
            local(k, cumul_h, h_face_k, h_in_bbl, mask_face, f_k)
            mask_face = min(wet_mask(i, j - 1), wet_mask(i, j))
            cumul_h = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
               if (h_face_k > 0.0_wp) then
                  h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
                  if (k == 1) then
                     f_k = bed_factor
                  else
                     f_k = 1.0_wp
                  end if
                  dv_drag(i, j, k) = mask_face*(-r*v_face(i, j, k)*h_in_bbl/h_face_k)*f_k &
                                     /(1.0_wp + dt_imp*r*(h_in_bbl/h_face_k)*f_k)
                  cumul_h = cumul_h + h_face_k
               else
                  exit
               end if
            end do
         end do
         return
      end if

      ! Quadratic with HBBL: need U_bbl, then apply.  Two-pass per
      ! face encoded as one DC with a sequential inner k-loop that
      ! does both passes (sum then write).
      if (variant == BDRAG_QUADRATIC .and. c_d > 0.0_wp) then
         do concurrent(j=2:ny, i=2:nx) &
            local(k, cumul_h, h_face_k, h_in_bbl, mask_face, f_k, &
                  h_in_bbl_total, u_bbl_int, v_bbl_int, u_bbl, v_bbl, &
                  v_at_u, abs_U_eff, h_eff_denom)
            mask_face = min(wet_mask(i - 1, j), wet_mask(i, j))
            ! Pass 1: integrate u_face * h_in_bbl over the BBL band.
            cumul_h = 0.0_wp
            u_bbl_int = 0.0_wp
            v_bbl_int = 0.0_wp
            h_in_bbl_total = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
               if (h_face_k <= 0.0_wp) exit
               h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
               u_bbl_int = u_bbl_int + u_face(i, j, k)*h_in_bbl
               v_at_u = 0.25_wp*( &
                        v_face(i - 1, j, k) + v_face(i, j, k) + &
                        v_face(i - 1, j + 1, k) + v_face(i, j + 1, k))
               v_bbl_int = v_bbl_int + v_at_u*h_in_bbl
               h_in_bbl_total = h_in_bbl_total + h_in_bbl
               cumul_h = cumul_h + h_face_k
            end do
            h_eff_denom = max(h_in_bbl_total, bbl_min)
            u_bbl = u_bbl_int/h_eff_denom
            v_bbl = v_bbl_int/h_eff_denom
            abs_U_eff = max(bg_vel, sqrt(u_bbl*u_bbl + v_bbl*v_bbl))
            ! Pass 2: apply stress per-layer.  Bed layer (k=1) scaled
            ! by `bed_factor`; layers k>=2 unchanged.
            cumul_h = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
               if (h_face_k <= 0.0_wp) exit
               h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
               if (k == 1) then
                  f_k = bed_factor
               else
                  f_k = 1.0_wp
               end if
               du_drag(i, j, k) = mask_face*( &
                                  -c_d*abs_U_eff*u_face(i, j, k)* &
                                  (h_in_bbl/h_face_k)/h_eff_denom)*f_k &
                                  /(1.0_wp + dt_imp*c_d*abs_U_eff*(h_in_bbl/h_face_k)/h_eff_denom*f_k)
               cumul_h = cumul_h + h_face_k
            end do
         end do
         do concurrent(j=2:ny, i=2:nx) &
            local(k, cumul_h, h_face_k, h_in_bbl, mask_face, f_k, &
                  h_in_bbl_total, u_bbl_int, v_bbl_int, u_bbl, v_bbl, &
                  u_at_v, abs_U_eff, h_eff_denom)
            mask_face = min(wet_mask(i, j - 1), wet_mask(i, j))
            cumul_h = 0.0_wp
            u_bbl_int = 0.0_wp
            v_bbl_int = 0.0_wp
            h_in_bbl_total = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
               if (h_face_k <= 0.0_wp) exit
               h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
               v_bbl_int = v_bbl_int + v_face(i, j, k)*h_in_bbl
               u_at_v = 0.25_wp*( &
                        u_face(i, j - 1, k) + u_face(i + 1, j - 1, k) + &
                        u_face(i, j, k) + u_face(i + 1, j, k))
               u_bbl_int = u_bbl_int + u_at_v*h_in_bbl
               h_in_bbl_total = h_in_bbl_total + h_in_bbl
               cumul_h = cumul_h + h_face_k
            end do
            h_eff_denom = max(h_in_bbl_total, bbl_min)
            u_bbl = u_bbl_int/h_eff_denom
            v_bbl = v_bbl_int/h_eff_denom
            abs_U_eff = max(bg_vel, sqrt(u_bbl*u_bbl + v_bbl*v_bbl))
            cumul_h = 0.0_wp
            do k = 1, nz
               if (cumul_h >= hbbl) exit
               h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
               if (h_face_k <= 0.0_wp) exit
               h_in_bbl = max(0.0_wp, min(h_face_k, hbbl - cumul_h))
               if (k == 1) then
                  f_k = bed_factor
               else
                  f_k = 1.0_wp
               end if
               dv_drag(i, j, k) = mask_face*( &
                                  -c_d*abs_U_eff*v_face(i, j, k)* &
                                  (h_in_bbl/h_face_k)/h_eff_denom)*f_k &
                                  /(1.0_wp + dt_imp*c_d*abs_U_eff*(h_in_bbl/h_face_k)/h_eff_denom*f_k)
               cumul_h = cumul_h + h_face_k
            end do
         end do
      end if
   end subroutine compute_distributed_drag

   subroutine ocean_bottom_drag_apply_tendencies(this, ms, dt, no_wait)
      !! `no_wait` (optional, default .false.): when .true. the apply DC
      !! loops run on OpenACC queue 1 and the routine returns WITHOUT
      !! syncing, so the batched velocity-apply chain in `run_stage_split`
      !! `!$acc wait(1)`s ONCE.  Default ⇒ blocking (safe for the unsplit
      !! `run_stage`).  Not `pure` because of the async/wait directives.
      type(ocean_bottom_drag_t), intent(in) :: this
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
   end subroutine ocean_bottom_drag_apply_tendencies

   ! =================================================================
   ! Channel (side-wall) drag — per-layer lateral Rayleigh drag
   ! =================================================================

   pure subroutine ocean_channel_drag_compute_tendencies(grid, metrics, this, ms)
      !! Per-layer lateral side-wall (channel) Rayleigh RATE
      !! (`lambda_side_u/v`, 1/s) for every velocity face whose
      !! cross-stream perimeter is partially blocked by land or a vanished
      !! (sloping-bathymetry) neighbour layer — fires at EVERY layer k that
      !! intersects the obstruction, not just the bed.
      !!   `f_blocked` = perimeter fraction blocked, summed over the two
      !!     flanking `wet_q` corners (each weighted 0.5), with a corner
      !!     also counted blocked when its two cross-stream cells' min
      !!     `h_layer < SIDE_H_VANISH`.
      !!   `lambda = cdrag_side·|U_face|·f_blocked/max(W, eps)`, W the
      !!     cross-stream face length (`dyCu` u-face, `dxCv` v-face).
      !! All-wet / flat-bottom ⇒ `f_blocked ≡ 0` ⇒ `lambda ≡ 0` (no-op);
      !! `channel_drag=.false.` or `cdrag_side=0` short-circuits to zero.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_bottom_drag_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call compute_channel_drag_rates(this%lambda_side_u%data, this%lambda_side_v%data, &
                                      ms%u_face_x_layer, ms%v_face_y_layer, &
                                      ms%h_layer, metrics%wet_q, &
                                      metrics%dyCu, metrics%dxCv, &
                                      this%channel_drag, this%cdrag_side, &
                                      size(ms%u_face_x_layer, 1), size(ms%u_face_x_layer, 2), &
                                      size(ms%v_face_y_layer, 1), size(ms%v_face_y_layer, 2), &
                                      nx, ny, nz)
   end subroutine ocean_channel_drag_compute_tendencies

   pure subroutine compute_channel_drag_rates(lambda_u, lambda_v, u_face, v_face, &
                                              h_layer, wet_q, dyCu, dxCv, &
                                              channel_drag, cdrag_side, &
                                              nx_u, ny_u, nx_v, ny_v, nx, ny, nz)
      !! Device kernel for the per-layer side-drag Rayleigh rate.  See
      !! `ocean_channel_drag_compute_tendencies` for the derivation.
      !! `wet_q` is `(nx+1, ny+1)`; `dyCu` is `(nx+1, ny)` (u-face length
      !! normal to the zonal flow); `dxCv` is `(nx, ny+1)`.
      integer, intent(in) :: nx_u, ny_u, nx_v, ny_v, nx, ny, nz
      logical, intent(in) :: channel_drag
      real(wp), intent(in) :: cdrag_side
      real(wp), intent(in) :: u_face(nx_u, ny_u, nz), v_face(nx_v, ny_v, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_q(nx + 1, ny + 1)
      real(wp), intent(in) :: dyCu(nx + 1, ny), dxCv(nx, ny + 1)
      real(wp), intent(out) :: lambda_u(nx_u, ny_u, nz), lambda_v(nx_v, ny_v, nz)

      integer :: i, j, k
      real(wp) :: f_s, f_n, f_blocked, speed, width, v_at_u, u_at_v
      real(wp), parameter :: WIDTH_EPS = 1.0e-20_wp

      ! ---- Zero everywhere first (also the all-disabled short-circuit) ----
      do concurrent(k=1:nz, j=1:ny_u, i=1:nx_u)
         lambda_u(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny_v, i=1:nx_v)
         lambda_v(i, j, k) = 0.0_wp
      end do
      if (.not. channel_drag) return
      if (cdrag_side <= 0.0_wp) return

      ! ---- u-faces: blocked by the two flanking corners (i,j),(i,j+1) ----
      ! Cross-stream cells at corner (i,j): T(i-1,j-1),T(i,j-1).
      ! Cross-stream cells at corner (i,j+1): T(i-1,j),T(i,j).
      ! Face range `j=2:ny` keeps the diagonal T reads in bounds: the
      ! south corner reads row `j-1` (>=1) and the north corner reads
      ! row `j` (<=ny); `wet_q(i,j+1)` reaches `j+1<=ny+1`, in bounds.
      ! Row `j=1` is a ghost row with no prognostic velocity, so leaving
      ! its rate at 0 is exact.
      do concurrent(k=1:nz, j=2:ny, i=2:nx) &
         local(f_s, f_n, f_blocked, speed, width, v_at_u)
         ! South corner (i,j): land OR layer vanished in T(i-1,j-1),T(i,j-1).
         f_s = 1.0_wp - wet_q(i, j)
         if (min(h_layer(i - 1, j - 1, k), h_layer(i, j - 1, k)) < SIDE_H_VANISH) then
            f_s = 1.0_wp
         end if
         ! North corner (i,j+1): land OR layer vanished in T(i-1,j),T(i,j).
         f_n = 1.0_wp - wet_q(i, j + 1)
         if (min(h_layer(i - 1, j, k), h_layer(i, j, k)) < SIDE_H_VANISH) then
            f_n = 1.0_wp
         end if
         f_blocked = 0.5_wp*f_s + 0.5_wp*f_n
         v_at_u = 0.25_wp*( &
                  v_face(i - 1, j, k) + v_face(i, j, k) + &
                  v_face(i - 1, j + 1, k) + v_face(i, j + 1, k))
         speed = sqrt(u_face(i, j, k)*u_face(i, j, k) + v_at_u*v_at_u)
         width = max(dyCu(i, j), WIDTH_EPS)
         lambda_u(i, j, k) = cdrag_side*speed*f_blocked/width
      end do

      ! ---- v-faces: blocked by the two flanking corners (i,j),(i+1,j) ----
      ! Cross-stream cells at corner (i,j): T(i-1,j-1),T(i-1,j).
      ! Cross-stream cells at corner (i+1,j): T(i,j-1),T(i,j).
      ! Face range `i=2:nx` keeps the diagonal T reads in bounds: the
      ! west corner reads column `i-1` (>=1) and the east corner reads
      ! column `i` (<=nx); `wet_q(i+1,j)` reaches `i+1<=nx+1`, in bounds.
      ! Column `i=1` is a ghost column with no prognostic velocity.
      do concurrent(k=1:nz, j=2:ny, i=2:nx) &
         local(f_s, f_n, f_blocked, speed, width, u_at_v)
         ! West corner (i,j): land OR layer vanished in T(i-1,j-1),T(i-1,j).
         f_s = 1.0_wp - wet_q(i, j)
         if (min(h_layer(i - 1, j - 1, k), h_layer(i - 1, j, k)) < SIDE_H_VANISH) then
            f_s = 1.0_wp
         end if
         ! East corner (i+1,j): land OR layer vanished in T(i,j-1),T(i,j).
         f_n = 1.0_wp - wet_q(i + 1, j)
         if (min(h_layer(i, j - 1, k), h_layer(i, j, k)) < SIDE_H_VANISH) then
            f_n = 1.0_wp
         end if
         f_blocked = 0.5_wp*f_s + 0.5_wp*f_n
         u_at_v = 0.25_wp*( &
                  u_face(i, j - 1, k) + u_face(i + 1, j - 1, k) + &
                  u_face(i, j, k) + u_face(i + 1, j, k))
         speed = sqrt(v_face(i, j, k)*v_face(i, j, k) + u_at_v*u_at_v)
         width = max(dxCv(i, j), WIDTH_EPS)
         lambda_v(i, j, k) = cdrag_side*speed*f_blocked/width
      end do
   end subroutine compute_channel_drag_rates

   subroutine ocean_channel_drag_apply_tendencies(this, ms, dt, no_wait)
      !! Apply the per-layer side drag IMPLICITLY:
      !!   `u <- u / (1 + dt·lambda_side_u)`,
      !!   `v <- v / (1 + dt·lambda_side_v)`.
      !! The implicit (backward-Euler) form is unconditionally stable on
      !! thin layers where an explicit `u - dt·lambda·u` would overshoot.
      !! `lambda ≡ 0` (default-off / all-wet / flat-bottom) ⇒ division by
      !! `1` ⇒ exact no-op.  `no_wait` semantics mirror
      !! `ocean_bottom_drag_apply_tendencies`.  Not `pure` (async/wait).
      type(ocean_bottom_drag_t), intent(in) :: this
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
         ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k)/ &
                                      (1.0_wp + dt*this%lambda_side_u%data(i, j, k))
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k)/ &
                                      (1.0_wp + dt*this%lambda_side_v%data(i, j, k))
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_channel_drag_apply_tendencies

   pure function parse_bdrag_variant(name) result(code)
      !! Translate a namelist string into a `BDRAG_*` code.  An
      !! unrecognised string returns `BDRAG_INVALID` (PR-6 fail-loud —
      !! a typo must not silently select the quadratic default over the
      !! linear form or vice versa; the two obey different physics).
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("quadratic", "QUADRATIC", "cd", "CD")
         code = BDRAG_QUADRATIC
      case ("linear", "LINEAR", "rayleigh", "RAYLEIGH")
         code = BDRAG_LINEAR
      case default
         code = BDRAG_INVALID
      end select
   end function parse_bdrag_variant

   pure function bdrag_variant_is_implemented(code) result(ok)
      !! `.true.` only for a bottom-drag variant with a real kernel
      !! (`BDRAG_LINEAR` / `BDRAG_QUADRATIC`).  `BDRAG_INVALID` returns
      !! `.false.`.  The single gate `validate_config` consumes (PR-6).
      integer, intent(in) :: code
      logical :: ok
      ok = (code == BDRAG_LINEAR) .or. (code == BDRAG_QUADRATIC)
   end function bdrag_variant_is_implemented

   pure function ocean_bottom_drag_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the bottom drag slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_bottom_drag_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%du_drag%bytes() &
               + this%dv_drag%bytes() &
               + arr_bytes(this%lambda_bot_u) &
               + arr_bytes(this%lambda_bot_v) &
               + this%lambda_side_u%bytes() &
               + this%lambda_side_v%bytes()
   end function ocean_bottom_drag_bytes

end module rdb_ocean_bottom_drag
