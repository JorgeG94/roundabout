!! Ocean lateral mixing parameterisation state.
module rdb_ocean_lateral_mix
   !! Flow-aware harmonic / biharmonic viscosity for the ocean dyn-core:
   !! per-face `ah_face_*` (m^2/s) and `nu4_face_*` (m^4/s), recomputed
   !! each step from the local flow and read by the horizontal-viscosity
   !! kernel in place of the scalar `nu_h`/`nu_4` (closure `LMIX_NONE`
   !! ⇒ scalar fallback ⇒ bit-identical).  The coastal path uses
   !! Smagorinsky in `rdb_ml_horizontal_viscosity`; the ocean path
   !! defaults to Leith, which scales with vorticity gradient and avoids
   !! over-damping coherent eddies.  Leith (1968); Smagorinsky (1963);
   !! Fox-Kemper & Menemenlis (2008); Griffies & Hallberg (2000).
   use rdb_constants, only: wp, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_lateral_mix_t
   public :: ocean_lateral_mix_compute_leith
   public :: ocean_lateral_mix_compute_smag
   public :: ocean_lateral_mix_compute_smag_ah
   public :: ocean_lateral_mix_compute_leith_biharm
   public :: ocean_lateral_mix_compute_vel_scale
   public :: ocean_lateral_mix_compute
   public :: parse_lateral_closure
   public :: lateral_closure_is_implemented
   public :: lateral_closure_conflicts_smag_ah
   public :: leith_biharm_is_inert
   public :: has_biharmonic_backstop
   public :: LMIX_NONE, LMIX_LEITH, LMIX_SMAGORINSKY, LMIX_BIHARMONIC, LMIX_LEITH_BIHARM
   public :: LMIX_INVALID

   ! Lateral closure tags.
   integer, parameter :: LMIX_INVALID = -1
      !! Sentinel for an unrecognised namelist string; aborts loudly at
      !! configure rather than silently falling back to background-only.
   integer, parameter :: LMIX_NONE = 0
      !! No flow-aware closure — falls back to the scalar `nu_h` field.
   integer, parameter :: LMIX_LEITH = 1
      !! Leith vorticity-gradient closure (default eddy-resolving).
   integer, parameter :: LMIX_SMAGORINSKY = 2
      !! Smagorinsky strain-rate closure.
   integer, parameter :: LMIX_BIHARMONIC = 3
      !! Constant-coefficient biharmonic.
   integer, parameter :: LMIX_LEITH_BIHARM = 4
      !! Leith-scaled biharmonic (Griffies & Hallberg 2000).

   type :: ocean_lateral_mix_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      integer  :: closure = LMIX_NONE
         !! Active closure tag.  Default `LMIX_NONE` ⇒ scalar-`nu_h`
         !! behaviour, bit-identical.
      logical  :: no_slip = .false.
         !! Lateral BC at coasts (shared with Coriolis).  `.false.`
         !! (default) = free-slip: corner shear strain `sh_xy` (and the
         !! Leith corner vorticity) is multiplied by `wet_q` so a land
         !! corner adds nothing.  `.true.` = no-slip: factor `2 - wet_q`.
         !! All-wet ⇒ `wet_q≡1` ⇒ bit-identical.
      real(wp) :: c_leith = 1.0_wp
         !! Leith dimensionless coefficient.  Typical 1.0–2.0.
      real(wp) :: c_smag = 0.15_wp
         !! Smagorinsky dimensionless coefficient (fallback).
      real(wp) :: ah_bg = 0.0_wp
         !! Background harmonic viscosity (m^2/s), floored beneath the
         !! closure to avoid zero damping in laminar patches.
      real(wp) :: ah_max = 1.0e4_wp
         !! Upper clip on harmonic viscosity (m^2/s).  Caps Leith spikes
         !! and enforces the viscous-CFL bound (nu·dt/dx² ≤ 0.5).

      real(wp) :: kh_vel_scale_live = 0.0_wp
         !! Live velocity-scale viscosity coefficient (m/s).  When
         !! positive, `A_vel = kh_vel_scale_live · L_grid · |u|`
         !! (`L_grid = sqrt(dxT·dyT)`) is `max`-combined into the
         !! per-face harmonic viscosity every step.  Default 0 ⇒ never
         !! computed ⇒ bit-identical.  State-dependent (evaluated per
         !! step); distinct from the `kh_vel_scale` background-floor knob
         !! set once at configure.  MOM6 `KH_VEL_SCALE` (Kh = U·Δ).

      ! ---- Biharmonic Smagorinsky ----
      logical  :: smag_ah_active = .false.
         !! When true, `compute_smag_ah` fills `nu4_face_x/y` each step
         !! from the local strain rate; the biharmonic kernel reads them
         !! instead of the scalar `nu_4`.  Independent of `closure` —
         !! Smag_KH (Laplacian) and Smag_AH (biharmonic) can both be on.
      logical  :: resoln_scaled_visc = .false.
         !! Hallberg (2013) resolution scaling.  When `.true.` AND the
         !! optional VarMix `res_fn_u/v` face fields are passed, the
         !! dynamic coefficients are multiplied by `Res_fn ∈ [0,1]`
         !! before the clamps (suppressed where the deformation radius
         !! is resolved).  Default `.false.` ⇒ unscaled ⇒ bit-identical.
      real(wp) :: smag_bi_const = 0.06_wp
         !! Nondimensional biharmonic Smagorinsky constant (typical
         !! 0.015–0.06).
      real(wp) :: c_leith_bi = 0.0_wp
         !! Nondimensional biharmonic Leith constant for
         !! `LMIX_LEITH_BIHARM` (Griffies & Hallberg 2000).  Default 0.0
         !! is a no-op; selecting the closure with it 0.0 warns at startup.
      real(wp) :: nu4_bg = 0.0_wp
         !! Background biharmonic viscosity floor (m⁴/s).
      real(wp) :: nu4_max = 1.0e12_wp
         !! Static upper clip on biharmonic viscosity (m⁴/s) applied when
         !! filling `nu4_face_*` — a cheap ceiling on the strain term.
         !! The true stability guard is the per-cell biharmonic-CFL clamp
         !! applied downstream in the biharmonic kernel.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- Face-located viscosity coefficients ----
      ! Recomputed every outer step (or per RK2 stage) from the local
      ! flow.  Shape matches u_face_x_layer / v_face_y_layer.
      real(wp), allocatable :: ah_face_x(:, :, :)
         !! Harmonic viscosity at east faces (m^2/s), shape
         !! `(nx+1, ny, nz_ml)`.
      real(wp), allocatable :: ah_face_y(:, :, :)
         !! Harmonic viscosity at north faces (m^2/s), shape
         !! `(nx, ny+1, nz_ml)`.
      real(wp), allocatable :: nu4_face_x(:, :, :)
         !! Biharmonic viscosity at east faces (m⁴/s), shape
         !! `(nx+1, ny, nz_ml)`.  Populated only when
         !! `smag_ah_active = .true.`.
      real(wp), allocatable :: nu4_face_y(:, :, :)
         !! Biharmonic viscosity at north faces (m⁴/s), shape
         !! `(nx, ny+1, nz_ml)`.  Populated only when
         !! `smag_ah_active = .true.`.

      ! ---- Vorticity-gradient scratch (Leith) ----
      ! ζ at corners, computed once per call from the face velocities.
      ! Shape `(nx+1, ny+1, nz_ml)`.
      type(scratch_3d_buffer_t) :: vort_corner
         !! Relative vorticity at C-grid corners.
   contains
      procedure, non_overridable :: init => ocean_lateral_mix_init
      procedure, non_overridable :: destroy => ocean_lateral_mix_destroy
      procedure, non_overridable :: enter_data => ocean_lateral_mix_enter_data
      procedure, non_overridable :: exit_data => ocean_lateral_mix_exit_data
      procedure, non_overridable :: bytes => ocean_lateral_mix_bytes
   end type ocean_lateral_mix_t

contains

   subroutine ocean_lateral_mix_init(this, grid, nz_ml)
      !! Allocate the face viscosity coefficients + corner-vorticity
      !! scratch.  Default `nz_ml = 1` preserves the barotropic-only
      !! constructor; pass `nz_ml = ms%nz_ml` for the multilayer driver.
      class(ocean_lateral_mix_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml
      if (nz < 1) nz = 1
      this%nx_total = nx
      this%ny_total = ny
      this%nz_ml = nz

      allocate (this%ah_face_x(nx + 1, ny, nz), source=0.0_wp)
      allocate (this%ah_face_y(nx, ny + 1, nz), source=0.0_wp)
      allocate (this%nu4_face_x(nx + 1, ny, nz), source=0.0_wp)
      allocate (this%nu4_face_y(nx, ny + 1, nz), source=0.0_wp)
      call this%vort_corner%init(nx + 1, ny + 1, nz, "lateral_mix_vort_corner")
      this%is_init = .true.
   end subroutine ocean_lateral_mix_init

   subroutine ocean_lateral_mix_destroy(this)
      class(ocean_lateral_mix_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%ah_face_x)) deallocate (this%ah_face_x)
      if (allocated(this%ah_face_y)) deallocate (this%ah_face_y)
      if (allocated(this%nu4_face_x)) deallocate (this%nu4_face_x)
      if (allocated(this%nu4_face_y)) deallocate (this%nu4_face_y)
      call this%vort_corner%destroy()
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_lateral_mix_destroy

   subroutine ocean_lateral_mix_enter_data(this)
      ! No bare `copyin(this)`: mapping the polymorphic dummy's stack
      ! descriptor makes a deep struct mapper on a recycled stack address
      ! that AMD libomptarget rejects (cross-slot overlap).  Slot header +
      ! scalar presence come from the orchestrator's root copyin(state)
      ! (every slot is inline in ocean_state_t).
      class(ocean_lateral_mix_t), intent(inout) :: this
      select type (this)
      type is (ocean_lateral_mix_t)
         call ocean_lateral_mix_enter_data_impl(this)
      end select
   end subroutine ocean_lateral_mix_enter_data

   subroutine ocean_lateral_mix_enter_data_impl(this)
      type(ocean_lateral_mix_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%ah_face_x, this%ah_face_y)
      !$acc enter data copyin(this%nu4_face_x, this%nu4_face_y)
      call scratch_3d_buffer_enter_data_impl(this%vort_corner)
   end subroutine ocean_lateral_mix_enter_data_impl

   subroutine ocean_lateral_mix_exit_data(this)
      class(ocean_lateral_mix_t), intent(inout) :: this
      select type (this)
      type is (ocean_lateral_mix_t)
         call ocean_lateral_mix_exit_data_impl(this)
      end select
   end subroutine ocean_lateral_mix_exit_data

   subroutine ocean_lateral_mix_exit_data_impl(this)
      type(ocean_lateral_mix_t), intent(inout) :: this
      if (.not. this%is_init) return
      call scratch_3d_buffer_exit_data_impl(this%vort_corner)
      !$acc exit data delete(this%nu4_face_y, this%nu4_face_x)
      !$acc exit data delete(this%ah_face_y, this%ah_face_x)
   end subroutine ocean_lateral_mix_exit_data_impl

   pure subroutine ocean_lateral_mix_compute_leith(grid, metrics, this, ms, &
                                                   res_fn_u, res_fn_v)
      !! Public only for the unit-test suite; ignore in production code.
      !! Populate `ah_face_x`/`ah_face_y` (m^2/s) with the Leith viscosity
      !!     A_h(face) = max(ah_bg, min(ah_max, (C_L · dx)^3 · |∇ζ|))
      !! where ζ is relative vorticity at C-grid corners (pass 1) and |∇ζ|
      !! the 2D gradient magnitude at each face (pass 2).  Wall faces get
      !! the background viscosity.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in), optional :: res_fn_u(grid%nx_total + 1, grid%ny_total)
         !! VarMix resolution function at u-faces (nondim, [0,1]).  When
         !! present AND `this%resoln_scaled_visc`, scales `A_h` before clamp.
      real(wp), intent(in), optional :: res_fn_v(grid%nx_total, grid%ny_total + 1)
         !! VarMix resolution function at v-faces.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: c_leith_local, ah_bg_local, ah_max_local
      real(wp) :: dzeta_dx, dzeta_dy, grad_mag, A_raw, leith_scale
      real(wp) :: ns
      logical :: do_resoln

      if (.not. this%is_init) return
      if (.not. allocated(ms%u_face_x_layer)) return
      if (.not. allocated(ms%v_face_y_layer)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Hoist scalar fields off `this` to locals — `do concurrent`
      ! bodies see plain real(wp) instead of a derived-type deref.
      ah_bg_local = this%ah_bg
      ah_max_local = this%ah_max
      c_leith_local = this%c_leith
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)   ! free-slip(0)/no-slip(1) selector (C1)
      ! Resolution-function scaling active only when the knob is on AND the
      ! VarMix face fields were supplied (Gap 1, Hallberg 2013); gated INSIDE
      ! the face loops so the optional is referenced only when present.
      do_resoln = this%resoln_scaled_visc .and. present(res_fn_u) .and. &
                  present(res_fn_v)

      ! Leith dimensionful prefactor: (C_L · L_grid)^3 with the grid
      ! scale `L_grid = sqrt(dxT·dyT)` evaluated per cell (design §2;
      ! = dx on uniform square metrics, so bit-reducing).  |∇ζ| has
      ! units 1/(m·s); A ~ L³·|∇ζ| is the right order for mesoscale
      ! closures.  The per-face scale below picks the adjacent T cell.

      ! ---- Pass 1: relative vorticity at SW corners, per layer ----
      ! Circulation/area form (consistent with the Coriolis kernel's
      ! converted zeta): ζ = (Δ(v·dyCv) − Δ(u·dxCu))·iareaBu.  Reduces
      ! to (Δv)/dx − (Δu)/dy on uniform square metrics.  Outer-most
      ! corners stay zero (closed-wall convention).
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         this%vort_corner%data(i, j, k) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((ms%v_face_y_layer(i, j, k)*metrics%dyCv(i, j) - &
              ms%v_face_y_layer(i - 1, j, k)*metrics%dyCv(i - 1, j)) - &
             (ms%u_face_x_layer(i, j, k)*metrics%dxCu(i, j) - &
              ms%u_face_x_layer(i, j - 1, k)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         this%vort_corner%data(1, j, k) = 0.0_wp
         this%vort_corner%data(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         this%vort_corner%data(i, 1, k) = 0.0_wp
         this%vort_corner%data(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 2a: A_h at u-faces (i-1/2, j) ----
      ! Adjacent corners: (i, j) at (i-1/2, j-1/2) and (i, j+1) at
      ! (i-1/2, j+1/2).  Across-face corners (one cell west/east):
      ! (i-1, j), (i-1, j+1), (i+1, j), (i+1, j+1).
      !
      ! ∂ζ/∂y at u-face: (ζ(i, j+1) - ζ(i, j))·idyCu  (along the face).
      ! ∂ζ/∂x at u-face: idxCu·[(ζ_E_S + ζ_E_N) - (ζ_W_S + ζ_W_N)]/4
      ! leith_scale = (C_L·sqrt(dxT·dyT))³ at the adjacent T cell
      ! (design §2; = (C_L·dx)³ on uniform square metrics).
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(dzeta_dx, dzeta_dy, grad_mag, A_raw, leith_scale)
         leith_scale = (c_leith_local*sqrt(metrics%dxT(i, j)*metrics%dyT(i, j)))**3
         dzeta_dy = (this%vort_corner%data(i, j + 1, k) - &
                     this%vort_corner%data(i, j, k))*metrics%idyCu(i, j)
         dzeta_dx = 0.25_wp*((this%vort_corner%data(i + 1, j, k) - &
                              this%vort_corner%data(i - 1, j, k)) + &
                             (this%vort_corner%data(i + 1, j + 1, k) - &
                              this%vort_corner%data(i - 1, j + 1, k)))*metrics%idxCu(i, j)
         grad_mag = sqrt(dzeta_dx*dzeta_dx + dzeta_dy*dzeta_dy)
         ! Resolution scaling applied BEFORE the clamp, as one assignment to
         ! the `local()` var `A_raw` per `do_resoln` branch (assigning a
         ! `local()` var once on each path; a conditional REASSIGN of a
         ! `do concurrent local()` var miscompiles on gfortran 15.1).
         if (do_resoln) then
            A_raw = leith_scale*grad_mag*res_fn_u(i, j)
         else
            A_raw = leith_scale*grad_mag
         end if
         this%ah_face_x(i, j, k) = min(ah_max_local, max(ah_bg_local, A_raw))
      end do
      ! Wall faces and i=1 edge: use background viscosity.
      do concurrent(k=1:nz, j=1:ny)
         this%ah_face_x(1, j, k) = ah_bg_local
         this%ah_face_x(nx + 1, j, k) = ah_bg_local
      end do

      ! ---- Pass 2b: A_h at v-faces (i, j-1/2) ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(dzeta_dx, dzeta_dy, grad_mag, A_raw, leith_scale)
         leith_scale = (c_leith_local*sqrt(metrics%dxT(i, j)*metrics%dyT(i, j)))**3
         dzeta_dx = (this%vort_corner%data(i + 1, j, k) - &
                     this%vort_corner%data(i, j, k))*metrics%idxCv(i, j)
         dzeta_dy = 0.25_wp*((this%vort_corner%data(i, j + 1, k) - &
                              this%vort_corner%data(i, j - 1, k)) + &
                             (this%vort_corner%data(i + 1, j + 1, k) - &
                              this%vort_corner%data(i + 1, j - 1, k)))*metrics%idyCv(i, j)
         grad_mag = sqrt(dzeta_dx*dzeta_dx + dzeta_dy*dzeta_dy)
         if (do_resoln) then
            A_raw = leith_scale*grad_mag*res_fn_v(i, j)
         else
            A_raw = leith_scale*grad_mag
         end if
         this%ah_face_y(i, j, k) = min(ah_max_local, max(ah_bg_local, A_raw))
      end do
      ! Wall faces: background viscosity.
      do concurrent(k=1:nz, i=1:nx)
         this%ah_face_y(i, 1, k) = ah_bg_local
         this%ah_face_y(i, ny + 1, k) = ah_bg_local
      end do
   end subroutine ocean_lateral_mix_compute_leith

   pure subroutine ocean_lateral_mix_compute_smag(grid, metrics, this, ms, &
                                                  res_fn_u, res_fn_v)
      !! Populate `ah_face_x`/`ah_face_y` (m^2/s) with the Smagorinsky
      !! Laplacian viscosity
      !!     A_h(face) = max(ah_bg, min(ah_max, (C_S · dx)^2 · |D|))
      !! where `|D| = sqrt(D_T^2 + D_S^2)` is the deformation-tensor
      !! magnitude — tension `D_T = ∂u/∂x − ∂v/∂y` (cell centred) and
      !! shear `D_S = ∂v/∂x + ∂u/∂y` (corner) — averaged onto the face.
      !! Wall faces get the background viscosity (wall-adjacent rows
      !! re-use the next interior row).  Smagorinsky (1963); C_S ≈ 0.15–0.2.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in), optional :: res_fn_u(grid%nx_total + 1, grid%ny_total)
         !! VarMix resolution function at u-faces (nondim, [0,1]).  When
         !! present AND `this%resoln_scaled_visc`, scales `A_h` before clamp.
      real(wp), intent(in), optional :: res_fn_v(grid%nx_total, grid%ny_total + 1)
         !! VarMix resolution function at v-faces.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: c_smag_local, smag_scale, ah_bg_local, ah_max_local
      real(wp) :: D_T_W, D_T_E, D_T_S, D_T_N, D_T_face
      real(wp) :: D_S_S, D_S_N, D_S_W, D_S_E, D_S_face
      real(wp) :: strain_mag, A_raw, ns
      logical :: do_resoln

      if (.not. this%is_init) return
      if (.not. allocated(ms%u_face_x_layer)) return
      if (.not. allocated(ms%v_face_y_layer)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Hoist scalar reads off `this` to local variables (see Leith
      ! sibling) — guards against device-side descriptor walks.
      ah_bg_local = this%ah_bg
      ah_max_local = this%ah_max
      c_smag_local = this%c_smag
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)   ! free-slip(0)/no-slip(1) selector (C1)
      ! Resolution-function scaling (Gap 1, Hallberg 2013) — see compute_leith.
      do_resoln = this%resoln_scaled_visc .and. present(res_fn_u) .and. &
                  present(res_fn_v)

      ! Strain components use per-stagger metric inverses (design §2):
      ! tension ∂u/∂x, ∂v/∂y on cell (idxT/idyT); shear ∂v/∂x, ∂u/∂y
      ! on the corner (idxBu/idyBu).  smag_scale = (C_S·sqrt(dxT·dyT))²
      ! per cell (= (C_S·dx)² on uniform square metrics, bit-reducing).

      ! ---- u-face viscosity (i-1/2, j) ----
      ! D_T at the face = mean of the two adjacent cell-centred values:
      !   D_T(i-1, j) and D_T(i, j).
      ! D_S at the face = mean of the two adjacent corner values along
      ! the face: D_S(i, j) at SW corner of (i, j) and D_S(i, j+1) at NW.
      !
      ! Curvilinear D_S at corner Bu(i,j) (= SW corner of T(i,j)):
      !   dvdx = dy_dxBu · (v(i,j)·idyCv(i,j)  - v(i-1,j)·idyCv(i-1,j))
      !   dudy = dx_dyBu · (u(i,j)·idxCu(i,j)  - u(i,j-1)·idxCu(i,j-1))
      !   D_S  = dvdx + dudy
      ! On uniform SQUARE metrics dy_dxBu=dx_dyBu=1 and idyCv=idxCu=1/dx,
      ! so this collapses to the old plain-difference form bit-for-bit.
      ! (design §2; mirrors MOM6 MOM_hor_visc shear-strain form)
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) &
         local(D_T_W, D_T_E, D_S_S, D_S_N, D_T_face, D_S_face, &
               strain_mag, A_raw, smag_scale)
         smag_scale = (c_smag_local*sqrt(metrics%dxT(i, j)*metrics%dyT(i, j)))**2
         ! Cell-centred D_T at (i-1, j) and (i, j)
         D_T_W = (ms%u_face_x_layer(i, j, k) - ms%u_face_x_layer(i - 1, j, k))*metrics%idxT(i - 1, j) &
                 - (ms%v_face_y_layer(i - 1, j + 1, k) - ms%v_face_y_layer(i - 1, j, k))*metrics%idyT(i - 1, j)
         D_T_E = (ms%u_face_x_layer(i + 1, j, k) - ms%u_face_x_layer(i, j, k))*metrics%idxT(i, j) &
                 - (ms%v_face_y_layer(i, j + 1, k) - ms%v_face_y_layer(i, j, k))*metrics%idyT(i, j)
         D_T_face = 0.5_wp*(D_T_W + D_T_E)
         ! Corner D_S: ratio-bundle form (design §2).  Bu(i,j) = SW corner
         ! of T(i,j); v at Cv(i,j) / Cv(i-1,j), u at Cu(i,j) / Cu(i,j-1).
         ! Corner shear strain sh_xy masked by slip factor (C1): free-slip
         ! ×wet_q / no-slip ×(2-wet_q).  Bit-identical for all-wet.
         D_S_S = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j)* &
                  (ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j) - &
                   ms%v_face_y_layer(i - 1, j, k)*metrics%idyCv(i - 1, j)) + &
                  metrics%dx_dyBu(i, j)* &
                  (ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j) - &
                   ms%u_face_x_layer(i, j - 1, k)*metrics%idxCu(i, j - 1)))
         D_S_N = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j + 1) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j + 1)* &
                  (ms%v_face_y_layer(i, j + 1, k)*metrics%idyCv(i, j + 1) - &
                   ms%v_face_y_layer(i - 1, j + 1, k)*metrics%idyCv(i - 1, j + 1)) + &
                  metrics%dx_dyBu(i, j + 1)* &
                  (ms%u_face_x_layer(i, j + 1, k)*metrics%idxCu(i, j + 1) - &
                   ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j)))
         D_S_face = 0.5_wp*(D_S_S + D_S_N)
         strain_mag = sqrt(D_T_face*D_T_face + D_S_face*D_S_face)
         if (do_resoln) then
            A_raw = smag_scale*strain_mag*res_fn_u(i, j)
         else
            A_raw = smag_scale*strain_mag
         end if
         this%ah_face_x(i, j, k) = min(ah_max_local, max(ah_bg_local, A_raw))
      end do
      ! j=1 and j=ny rows: re-use the j=2 / j=ny-1 values one row in.
      ! Avoids the j-1 / j+1 stencil walking into the wall.
      do concurrent(k=1:nz, i=2:nx)
         this%ah_face_x(i, 1, k) = this%ah_face_x(i, 2, k)
         this%ah_face_x(i, ny, k) = this%ah_face_x(i, ny - 1, k)
      end do
      ! Wall faces (i=1, i=nx+1): background.
      do concurrent(k=1:nz, j=1:ny)
         this%ah_face_x(1, j, k) = ah_bg_local
         this%ah_face_x(nx + 1, j, k) = ah_bg_local
      end do

      ! ---- v-face viscosity (i, j-1/2) ----
      ! Mirror of the u-face stencil.
      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) &
         local(D_T_S, D_T_N, D_S_W, D_S_E, D_T_face, D_S_face, &
               strain_mag, A_raw, smag_scale)
         smag_scale = (c_smag_local*sqrt(metrics%dxT(i, j)*metrics%dyT(i, j)))**2
         ! Cell-centred D_T at (i, j-1) and (i, j)
         D_T_S = (ms%u_face_x_layer(i + 1, j - 1, k) - ms%u_face_x_layer(i, j - 1, k))*metrics%idxT(i, j - 1) &
                 - (ms%v_face_y_layer(i, j, k) - ms%v_face_y_layer(i, j - 1, k))*metrics%idyT(i, j - 1)
         D_T_N = (ms%u_face_x_layer(i + 1, j, k) - ms%u_face_x_layer(i, j, k))*metrics%idxT(i, j) &
                 - (ms%v_face_y_layer(i, j + 1, k) - ms%v_face_y_layer(i, j, k))*metrics%idyT(i, j)
         D_T_face = 0.5_wp*(D_T_S + D_T_N)
         ! Corner D_S: ratio-bundle form (design §2).  Bu(i,j) = SW corner
         ! of T(i,j); v at Cv(i,j) / Cv(i-1,j), u at Cu(i,j) / Cu(i,j-1).
         ! Corner shear strain sh_xy masked by slip factor (C1).
         D_S_W = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j)* &
                  (ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j) - &
                   ms%v_face_y_layer(i - 1, j, k)*metrics%idyCv(i - 1, j)) + &
                  metrics%dx_dyBu(i, j)* &
                  (ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j) - &
                   ms%u_face_x_layer(i, j - 1, k)*metrics%idxCu(i, j - 1)))
         D_S_E = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i + 1, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i + 1, j)* &
                  (ms%v_face_y_layer(i + 1, j, k)*metrics%idyCv(i + 1, j) - &
                   ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j)) + &
                  metrics%dx_dyBu(i + 1, j)* &
                  (ms%u_face_x_layer(i + 1, j, k)*metrics%idxCu(i + 1, j) - &
                   ms%u_face_x_layer(i + 1, j - 1, k)*metrics%idxCu(i + 1, j - 1)))
         D_S_face = 0.5_wp*(D_S_W + D_S_E)
         strain_mag = sqrt(D_T_face*D_T_face + D_S_face*D_S_face)
         if (do_resoln) then
            A_raw = smag_scale*strain_mag*res_fn_v(i, j)
         else
            A_raw = smag_scale*strain_mag
         end if
         this%ah_face_y(i, j, k) = min(ah_max_local, max(ah_bg_local, A_raw))
      end do
      ! i=1, i=nx columns: re-use the i=2 / i=nx-1 values.
      do concurrent(k=1:nz, j=2:ny)
         this%ah_face_y(1, j, k) = this%ah_face_y(2, j, k)
         this%ah_face_y(nx, j, k) = this%ah_face_y(nx - 1, j, k)
      end do
      ! Wall faces (j=1, j=ny+1): background.
      do concurrent(k=1:nz, i=1:nx)
         this%ah_face_y(i, 1, k) = ah_bg_local
         this%ah_face_y(i, ny + 1, k) = ah_bg_local
      end do
   end subroutine ocean_lateral_mix_compute_smag

   pure subroutine ocean_lateral_mix_compute_smag_ah(grid, metrics, this, ms)
      !! Public only for the unit-test suite; ignore in production code.
      !! Populate `nu4_face_x`/`nu4_face_y` (m⁴/s) with the biharmonic
      !! Smagorinsky viscosity
      !!     A_4(face) = clamp(C_b · L⁴ · |D|, nu4_bg, nu4_max)
      !! where `L² = 2·dx²·dy²/(dx²+dy²)` (harmonic mean of dx²,dy²) and
      !! `|D|` is the strain-rate magnitude from `compute_smag`.  Wall
      !! faces get `nu4_bg`.  `SMAG_BI_CONST` ≈ 0.015–0.06.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: smag_bi_const_local, smag_bi_scale, nu4_bg_local, nu4_max_local
      real(wp) :: dx2, dy2, grid_sp_h2
      real(wp) :: D_T_W, D_T_E, D_T_S, D_T_N, D_T_face
      real(wp) :: D_S_S, D_S_N, D_S_W, D_S_E, D_S_face
      real(wp) :: strain_mag, A_raw, ns

      if (.not. this%is_init) return
      if (.not. allocated(ms%u_face_x_layer)) return
      if (.not. allocated(ms%v_face_y_layer)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)   ! free-slip(0)/no-slip(1) selector (C1)
      ! Biharmonic ν₄ is deliberately NOT resolution-scaled: the
      ! Hallberg (2013) resolution function suppresses only the
      ! scale-non-selective Laplacian, whereas the ∝k⁴ biharmonic already
      ! spares the resolved (large) scales and needs no suppression.

      ! Per-cell `C_b · (grid_sp_h2)^2` (MOM6 `Biharm_const_xx`), with
      ! grid_sp_h2 = 2·dx2h·dy2h/(dx2h+dy2h) the harmonic mean of the
      ! per-cell dxT²/dyT² (design §2; = the uniform value on square
      ! metrics, bit-reducing).  Strain inverses per stagger as in
      ! `compute_smag`.
      nu4_bg_local = this%nu4_bg
      nu4_max_local = this%nu4_max
      smag_bi_const_local = this%smag_bi_const

      ! ---- u-face viscosity (i-1/2, j) ----
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) &
         local(D_T_W, D_T_E, D_S_S, D_S_N, D_T_face, D_S_face, &
               strain_mag, A_raw, dx2, dy2, grid_sp_h2, smag_bi_scale)
         dx2 = metrics%dx2h(i, j)
         dy2 = metrics%dy2h(i, j)
         grid_sp_h2 = (2.0_wp*dx2*dy2)/(dx2 + dy2)
         smag_bi_scale = smag_bi_const_local*(grid_sp_h2*grid_sp_h2)
         D_T_W = (ms%u_face_x_layer(i, j, k) - ms%u_face_x_layer(i - 1, j, k))*metrics%idxT(i - 1, j) &
                 - (ms%v_face_y_layer(i - 1, j + 1, k) - ms%v_face_y_layer(i - 1, j, k))*metrics%idyT(i - 1, j)
         D_T_E = (ms%u_face_x_layer(i + 1, j, k) - ms%u_face_x_layer(i, j, k))*metrics%idxT(i, j) &
                 - (ms%v_face_y_layer(i, j + 1, k) - ms%v_face_y_layer(i, j, k))*metrics%idyT(i, j)
         D_T_face = 0.5_wp*(D_T_W + D_T_E)
         ! Corner D_S: ratio-bundle form (design §2).  Bu(i,j) = SW corner
         ! of T(i,j); v at Cv(i,j)/Cv(i-1,j), u at Cu(i,j)/Cu(i,j-1).
         ! Corner shear strain sh_xy masked by slip factor (C1): free-slip
         ! ×wet_q / no-slip ×(2-wet_q).  Bit-identical for all-wet.
         D_S_S = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j)* &
                  (ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j) - &
                   ms%v_face_y_layer(i - 1, j, k)*metrics%idyCv(i - 1, j)) + &
                  metrics%dx_dyBu(i, j)* &
                  (ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j) - &
                   ms%u_face_x_layer(i, j - 1, k)*metrics%idxCu(i, j - 1)))
         D_S_N = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j + 1) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j + 1)* &
                  (ms%v_face_y_layer(i, j + 1, k)*metrics%idyCv(i, j + 1) - &
                   ms%v_face_y_layer(i - 1, j + 1, k)*metrics%idyCv(i - 1, j + 1)) + &
                  metrics%dx_dyBu(i, j + 1)* &
                  (ms%u_face_x_layer(i, j + 1, k)*metrics%idxCu(i, j + 1) - &
                   ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j)))
         D_S_face = 0.5_wp*(D_S_S + D_S_N)
         strain_mag = sqrt(D_T_face*D_T_face + D_S_face*D_S_face)
         A_raw = smag_bi_scale*strain_mag
         this%nu4_face_x(i, j, k) = min(nu4_max_local, max(nu4_bg_local, A_raw))
      end do
      do concurrent(k=1:nz, i=2:nx)
         this%nu4_face_x(i, 1, k) = this%nu4_face_x(i, 2, k)
         this%nu4_face_x(i, ny, k) = this%nu4_face_x(i, ny - 1, k)
      end do
      do concurrent(k=1:nz, j=1:ny)
         this%nu4_face_x(1, j, k) = nu4_bg_local
         this%nu4_face_x(nx + 1, j, k) = nu4_bg_local
      end do

      ! ---- v-face viscosity (i, j-1/2) ----
      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) &
         local(D_T_S, D_T_N, D_S_W, D_S_E, D_T_face, D_S_face, &
               strain_mag, A_raw, dx2, dy2, grid_sp_h2, smag_bi_scale)
         dx2 = metrics%dx2h(i, j)
         dy2 = metrics%dy2h(i, j)
         grid_sp_h2 = (2.0_wp*dx2*dy2)/(dx2 + dy2)
         smag_bi_scale = smag_bi_const_local*(grid_sp_h2*grid_sp_h2)
         D_T_S = (ms%u_face_x_layer(i + 1, j - 1, k) - ms%u_face_x_layer(i, j - 1, k))*metrics%idxT(i, j - 1) &
                 - (ms%v_face_y_layer(i, j, k) - ms%v_face_y_layer(i, j - 1, k))*metrics%idyT(i, j - 1)
         D_T_N = (ms%u_face_x_layer(i + 1, j, k) - ms%u_face_x_layer(i, j, k))*metrics%idxT(i, j) &
                 - (ms%v_face_y_layer(i, j + 1, k) - ms%v_face_y_layer(i, j, k))*metrics%idyT(i, j)
         D_T_face = 0.5_wp*(D_T_S + D_T_N)
         ! Corner D_S: ratio-bundle form (design §2).  Bu(i,j) = SW corner
         ! of T(i,j); v at Cv(i,j)/Cv(i-1,j), u at Cu(i,j)/Cu(i,j-1).
         ! Corner shear strain sh_xy masked by slip factor (C1).
         D_S_W = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i, j)* &
                  (ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j) - &
                   ms%v_face_y_layer(i - 1, j, k)*metrics%idyCv(i - 1, j)) + &
                  metrics%dx_dyBu(i, j)* &
                  (ms%u_face_x_layer(i, j, k)*metrics%idxCu(i, j) - &
                   ms%u_face_x_layer(i, j - 1, k)*metrics%idxCu(i, j - 1)))
         D_S_E = ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i + 1, j) + 2.0_wp*ns)* &
                 (metrics%dy_dxBu(i + 1, j)* &
                  (ms%v_face_y_layer(i + 1, j, k)*metrics%idyCv(i + 1, j) - &
                   ms%v_face_y_layer(i, j, k)*metrics%idyCv(i, j)) + &
                  metrics%dx_dyBu(i + 1, j)* &
                  (ms%u_face_x_layer(i + 1, j, k)*metrics%idxCu(i + 1, j) - &
                   ms%u_face_x_layer(i + 1, j - 1, k)*metrics%idxCu(i + 1, j - 1)))
         D_S_face = 0.5_wp*(D_S_W + D_S_E)
         strain_mag = sqrt(D_T_face*D_T_face + D_S_face*D_S_face)
         A_raw = smag_bi_scale*strain_mag
         this%nu4_face_y(i, j, k) = min(nu4_max_local, max(nu4_bg_local, A_raw))
      end do
      do concurrent(k=1:nz, j=2:ny)
         this%nu4_face_y(1, j, k) = this%nu4_face_y(2, j, k)
         this%nu4_face_y(nx, j, k) = this%nu4_face_y(nx - 1, j, k)
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%nu4_face_y(i, 1, k) = nu4_bg_local
         this%nu4_face_y(i, ny + 1, k) = nu4_bg_local
      end do
   end subroutine ocean_lateral_mix_compute_smag_ah

   pure subroutine ocean_lateral_mix_compute_leith_biharm(grid, metrics, this, ms)
      !! Public only for the unit-test suite; ignore in production code.
      !! Populate `nu4_face_x`/`nu4_face_y` (m⁴/s) with the 2-D Leith
      !! biharmonic viscosity
      !!     A_4(face) = clamp(C_lb · grid_sp⁶ · inv_PI6 · |∇²ζ|,
      !!                       nu4_bg, nu4_max)
      !! where ζ is C-grid corner relative vorticity, `∇²ζ` its 5-point
      !! corner Laplacian, `grid_sp⁶ = grid_sp_h2³`, and `inv_PI6 = (1/π)⁶`.
      !! Per-face |∇²ζ| is the mean of the two adjacent corner Laplacians.
      !! Wall faces get `nu4_bg`.  Leith (1968); Griffies & Hallberg (2000).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: c_leith_bi_local, nu4_bg_local, nu4_max_local, ns, inv_pi6
      real(wp) :: dx2, dy2, grid_sp_h2, grid_sp6, leith_bi_scale
      real(wp) :: del2_a, del2_b, del2_face, A_raw

      if (.not. this%is_init) return
      if (.not. allocated(ms%u_face_x_layer)) return
      if (.not. allocated(ms%v_face_y_layer)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Hoist scalar fields off `this` to locals (see Leith/Smag
      ! siblings) — `do concurrent` bodies see plain real(wp).
      nu4_bg_local = this%nu4_bg
      nu4_max_local = this%nu4_max
      c_leith_bi_local = this%c_leith_bi
      ns = merge(1.0_wp, 0.0_wp, this%no_slip)   ! free-slip(0)/no-slip(1) selector (C1)
      inv_pi6 = (1.0_wp/PI)**6

      ! ---- Pass 1: relative vorticity at SW corners, per layer ----
      ! Identical circulation/area form to compute_leith (consistent
      ! with the Coriolis kernel's converted zeta).  Outer-most corners
      ! stay zero (closed-wall convention) so the corner Laplacian below
      ! sees a finite neighbourhood at the first interior corners.
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         this%vort_corner%data(i, j, k) = &
            ((1.0_wp - 2.0_wp*ns)*metrics%wet_q(i, j) + 2.0_wp*ns)* &
            ((ms%v_face_y_layer(i, j, k)*metrics%dyCv(i, j) - &
              ms%v_face_y_layer(i - 1, j, k)*metrics%dyCv(i - 1, j)) - &
             (ms%u_face_x_layer(i, j, k)*metrics%dxCu(i, j) - &
              ms%u_face_x_layer(i, j - 1, k)*metrics%dxCu(i, j - 1)))* &
            metrics%iareaBu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         this%vort_corner%data(1, j, k) = 0.0_wp
         this%vort_corner%data(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         this%vort_corner%data(i, 1, k) = 0.0_wp
         this%vort_corner%data(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 2a: A_4 at u-faces (i-1/2, j) ----
      ! ∇²ζ at the two adjacent corners Bu(i,j) (SW) and Bu(i,j+1) (NW)
      ! averaged onto the face; |∇²ζ| scaled by C_lb·grid_sp⁶·inv_PI6.
      ! Corner Laplacian needs j∈[2,ny-1] (j-1/j+1 in range for the NW
      ! corner at j+1); wall-adjacent rows are filled by row-copy below.
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) &
         local(dx2, dy2, grid_sp_h2, grid_sp6, leith_bi_scale, &
               del2_a, del2_b, del2_face, A_raw)
         dx2 = metrics%dx2h(i, j)
         dy2 = metrics%dy2h(i, j)
         grid_sp_h2 = (2.0_wp*dx2*dy2)/(dx2 + dy2)
         grid_sp6 = grid_sp_h2*grid_sp_h2*grid_sp_h2
         leith_bi_scale = c_leith_bi_local*grid_sp6*inv_pi6
         ! ∇²ζ at SW corner Bu(i,j): 5-point corner Laplacian.
         del2_a = (this%vort_corner%data(i + 1, j, k) - &
                   2.0_wp*this%vort_corner%data(i, j, k) + &
                   this%vort_corner%data(i - 1, j, k))/metrics%dx2q(i, j) + &
                  (this%vort_corner%data(i, j + 1, k) - &
                   2.0_wp*this%vort_corner%data(i, j, k) + &
                   this%vort_corner%data(i, j - 1, k))/metrics%dy2q(i, j)
         ! ∇²ζ at NW corner Bu(i,j+1).
         del2_b = (this%vort_corner%data(i + 1, j + 1, k) - &
                   2.0_wp*this%vort_corner%data(i, j + 1, k) + &
                   this%vort_corner%data(i - 1, j + 1, k))/metrics%dx2q(i, j + 1) + &
                  (this%vort_corner%data(i, j + 2, k) - &
                   2.0_wp*this%vort_corner%data(i, j + 1, k) + &
                   this%vort_corner%data(i, j, k))/metrics%dy2q(i, j + 1)
         del2_face = 0.5_wp*(del2_a + del2_b)
         A_raw = leith_bi_scale*abs(del2_face)
         this%nu4_face_x(i, j, k) = min(nu4_max_local, max(nu4_bg_local, A_raw))
      end do
      do concurrent(k=1:nz, i=2:nx)
         this%nu4_face_x(i, 1, k) = this%nu4_face_x(i, 2, k)
         this%nu4_face_x(i, ny, k) = this%nu4_face_x(i, ny - 1, k)
      end do
      do concurrent(k=1:nz, j=1:ny)
         this%nu4_face_x(1, j, k) = nu4_bg_local
         this%nu4_face_x(nx + 1, j, k) = nu4_bg_local
      end do

      ! ---- Pass 2b: A_4 at v-faces (i, j-1/2) ----
      ! Adjacent corners Bu(i,j) (SW) and Bu(i+1,j) (SE) averaged onto
      ! the face.  i∈[2,nx-1] keeps i-1/i+1 in range for the SE corner.
      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) &
         local(dx2, dy2, grid_sp_h2, grid_sp6, leith_bi_scale, &
               del2_a, del2_b, del2_face, A_raw)
         dx2 = metrics%dx2h(i, j)
         dy2 = metrics%dy2h(i, j)
         grid_sp_h2 = (2.0_wp*dx2*dy2)/(dx2 + dy2)
         grid_sp6 = grid_sp_h2*grid_sp_h2*grid_sp_h2
         leith_bi_scale = c_leith_bi_local*grid_sp6*inv_pi6
         ! ∇²ζ at SW corner Bu(i,j).
         del2_a = (this%vort_corner%data(i + 1, j, k) - &
                   2.0_wp*this%vort_corner%data(i, j, k) + &
                   this%vort_corner%data(i - 1, j, k))/metrics%dx2q(i, j) + &
                  (this%vort_corner%data(i, j + 1, k) - &
                   2.0_wp*this%vort_corner%data(i, j, k) + &
                   this%vort_corner%data(i, j - 1, k))/metrics%dy2q(i, j)
         ! ∇²ζ at SE corner Bu(i+1,j).
         del2_b = (this%vort_corner%data(i + 2, j, k) - &
                   2.0_wp*this%vort_corner%data(i + 1, j, k) + &
                   this%vort_corner%data(i, j, k))/metrics%dx2q(i + 1, j) + &
                  (this%vort_corner%data(i + 1, j + 1, k) - &
                   2.0_wp*this%vort_corner%data(i + 1, j, k) + &
                   this%vort_corner%data(i + 1, j - 1, k))/metrics%dy2q(i + 1, j)
         del2_face = 0.5_wp*(del2_a + del2_b)
         A_raw = leith_bi_scale*abs(del2_face)
         this%nu4_face_y(i, j, k) = min(nu4_max_local, max(nu4_bg_local, A_raw))
      end do
      do concurrent(k=1:nz, j=2:ny)
         this%nu4_face_y(1, j, k) = this%nu4_face_y(2, j, k)
         this%nu4_face_y(nx, j, k) = this%nu4_face_y(nx - 1, j, k)
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%nu4_face_y(i, 1, k) = nu4_bg_local
         this%nu4_face_y(i, ny + 1, k) = nu4_bg_local
      end do
   end subroutine ocean_lateral_mix_compute_leith_biharm

   pure subroutine ocean_lateral_mix_compute_vel_scale(grid, metrics, this, ms, seed_bg)
      !! Public only for the unit-test suite; ignore in production code.
      !! Live velocity-scale viscosity (MOM6 `KH_VEL_SCALE`, Kh = U·Δ):
      !! per face `A_vel = kh_vel_scale_live · L_grid · |u_face|`
      !! (`L_grid = sqrt(dxT·dyT)`), `max`-combined into `ah_face_*` so it
      !! floors — never reduces — the active closure.  `seed_bg = .true.`
      !! first fills every face with `ah_bg` (used when no closure ran);
      !! `.false.` only raises faces where `A_vel` exceeds the closure.
      !! No-op when `kh_vel_scale_live <= 0`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      logical, intent(in) :: seed_bg

      integer :: i, j, k, nx, ny, nz
      real(wp) :: vel_scale_local, ah_bg_local, ah_max_local
      real(wp) :: a_vel, l_grid

      if (.not. this%is_init) return
      if (this%kh_vel_scale_live <= 0.0_wp) return
      if (.not. allocated(ms%u_face_x_layer)) return
      if (.not. allocated(ms%v_face_y_layer)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      vel_scale_local = this%kh_vel_scale_live
      ah_bg_local = this%ah_bg
      ah_max_local = this%ah_max

      ! ---- u-faces (i-1/2, j): A_vel from |u_face_x| ----
      ! Interior faces i=2:nx (the closures' range); L_grid taken at
      ! the west-adjacent T cell (i-1).  Cap and floor as the closures
      ! do.  Wall faces (i=1, i=nx+1) get the background under seeding.
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(a_vel, l_grid)
         l_grid = sqrt(metrics%dxT(i - 1, j)*metrics%dyT(i - 1, j))
         a_vel = min(ah_max_local, vel_scale_local*l_grid*abs(ms%u_face_x_layer(i, j, k)))
         if (seed_bg) then
            this%ah_face_x(i, j, k) = max(ah_bg_local, a_vel)
         else
            this%ah_face_x(i, j, k) = max(this%ah_face_x(i, j, k), a_vel)
         end if
      end do
      do concurrent(k=1:nz, j=1:ny)
         if (seed_bg) then
            this%ah_face_x(1, j, k) = ah_bg_local
            this%ah_face_x(nx + 1, j, k) = ah_bg_local
         end if
      end do

      ! ---- v-faces (i, j-1/2): A_vel from |v_face_y| ----
      ! Interior faces j=2:ny; L_grid at the south-adjacent T cell.
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(a_vel, l_grid)
         l_grid = sqrt(metrics%dxT(i, j - 1)*metrics%dyT(i, j - 1))
         a_vel = min(ah_max_local, vel_scale_local*l_grid*abs(ms%v_face_y_layer(i, j, k)))
         if (seed_bg) then
            this%ah_face_y(i, j, k) = max(ah_bg_local, a_vel)
         else
            this%ah_face_y(i, j, k) = max(this%ah_face_y(i, j, k), a_vel)
         end if
      end do
      do concurrent(k=1:nz, i=1:nx)
         if (seed_bg) then
            this%ah_face_y(i, 1, k) = ah_bg_local
            this%ah_face_y(i, ny + 1, k) = ah_bg_local
         end if
      end do
   end subroutine ocean_lateral_mix_compute_vel_scale

   subroutine ocean_lateral_mix_compute(grid, metrics, this, ms, &
                                        res_fn_u, res_fn_v)
      !! Dispatcher — runs the compute kernel for the active closure tag.
      !! `LMIX_LEITH`/`LMIX_SMAGORINSKY` → harmonic `ah_face_*`;
      !! `LMIX_LEITH_BIHARM` → biharmonic `nu4_face_*`; `LMIX_BIHARMONIC`
      !! → scalar `nu_4` in the apply step (no per-face fill); `LMIX_NONE`
      !! → no-op.  The independent `smag_ah_active` switch separately
      !! fills `nu4_face_*` from the strain rate.  `this` is optional so
      !! the driver can call unconditionally.  `res_fn_u/v` (optional
      !! VarMix resolution-function fields): when supplied AND
      !! `resoln_scaled_visc`, scale the dynamic coefficients before the
      !! clamps; absent ⇒ unscaled ⇒ bit-identical.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_lateral_mix_t), intent(inout), optional :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in), optional :: res_fn_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in), optional :: res_fn_v(grid%nx_total, grid%ny_total + 1)
      logical :: ran_closure

      if (.not. present(this)) return
      if (.not. this%is_init) return
      ! `ran_closure` tracks whether a flow-aware closure populated the
      ! HARMONIC `ah_face_*` field this step — only LMIX_LEITH/SMAGORINSKY
      ! do (the biharmonic closures fill `nu4_face_*` instead).  It gates
      ! whether the live velocity-scale floor must first seed `ah_bg`.
      ran_closure = .false.
      ! Forward the resolution-function fields straight through; Fortran
      ! 2008+ propagates `present()` across optional dummies, so an absent
      ! actual stays absent in the callee (kernels self-gate on present).
      select case (this%closure)
      case (LMIX_LEITH)
         call ocean_lateral_mix_compute_leith(grid, metrics, this, ms, &
                                              res_fn_u=res_fn_u, res_fn_v=res_fn_v)
         ran_closure = .true.
      case (LMIX_SMAGORINSKY)
         call ocean_lateral_mix_compute_smag(grid, metrics, this, ms, &
                                             res_fn_u=res_fn_u, res_fn_v=res_fn_v)
         ran_closure = .true.
      case (LMIX_LEITH_BIHARM)
         call ocean_lateral_mix_compute_leith_biharm(grid, metrics, this, ms)
      case (LMIX_BIHARMONIC)
         ! Constant biharmonic — driven by the scalar `nu_4` in the
         ! apply step; no flow-aware per-face fill required.
      case default
         ! LMIX_NONE — scalar `nu_h` Laplacian, no flow-aware closure.
      end select
      ! Live velocity-scale floor — max-combined into `ah_face_*`; seed the
      ! background first only when no harmonic closure populated the field.
      if (this%kh_vel_scale_live > 0.0_wp) then
         call ocean_lateral_mix_compute_vel_scale(grid, metrics, this, ms, &
                                                  seed_bg=.not. ran_closure)
      end if
      if (this%smag_ah_active) then
         call ocean_lateral_mix_compute_smag_ah(grid, metrics, this, ms)
      end if
   end subroutine ocean_lateral_mix_compute

   pure function lateral_closure_is_implemented(code) result(ok)
      !! `.true.` iff the closure code has a working dispatcher path.
      !! Drives the configure-time fail-loud guard.  Keep in lock-step
      !! with the `select case` in `ocean_lateral_mix_compute`.
      integer, intent(in) :: code
      logical :: ok
      select case (code)
      case (LMIX_NONE, LMIX_LEITH, LMIX_SMAGORINSKY, &
            LMIX_BIHARMONIC, LMIX_LEITH_BIHARM)
         ok = .true.
      case default
         ok = .false.
      end select
   end function lateral_closure_is_implemented

   pure function lateral_closure_conflicts_smag_ah(code, smag_ah) result(conflict)
      !! `.true.` iff the closure is `LMIX_LEITH_BIHARM` and `smag_ah` is
      !! on — both would fill `nu4_face_*` and `smag_ah` runs last, so it
      !! would silently overwrite the Leith-biharmonic fill (we do not
      !! max-combine biharmonic closures).  Drives a fail-loud guard.
      integer, intent(in) :: code
      logical, intent(in) :: smag_ah
      logical :: conflict
      conflict = (code == LMIX_LEITH_BIHARM) .and. smag_ah
   end function lateral_closure_conflicts_smag_ah

   pure function leith_biharm_is_inert(code, c_leith_bi) result(inert)
      !! `.true.` iff the Leith-biharmonic closure is selected but its
      !! coefficient is <= 0 (PR-6 fail-loud).  ν₄ is linear in
      !! `c_leith_bi`, so `c_leith_bi <= 0` makes the whole closure a
      !! provable no-op — the user asked for biharmonic dissipation and
      !! got none.  Returns `.false.` for any other closure (they do not
      !! read `c_leith_bi`, so the guard must not fire on them).  Promotes
      !! the previous configure-time warning to an abort.
      integer, intent(in) :: code
      real(wp), intent(in) :: c_leith_bi
      logical :: inert
      inert = (code == LMIX_LEITH_BIHARM) .and. (c_leith_bi <= 0.0_wp)
   end function leith_biharm_is_inert

   pure function has_biharmonic_backstop(nu_4, smag_ah, smag_bi_const, code, &
                                         c_leith_bi, nu_4_bg) result(ok)
      !! `.true.` iff the configured biharmonic dispatch will produce a
      !! STRICTLY POSITIVE dissipation coefficient somewhere.  Drives the
      !! configure-time fail-loud guard requiring a biharmonic backstop when
      !! MEKE backscatter is on — the negative harmonic backscatter feeds a
      !! grid-scale mode that only a positive biharmonic can dissipate.
      !!
      !! Two failure modes closed here (both silently passed the guard
      !! before this signature): (1) a flow-aware closure SELECTED with a
      !! zero coefficient (`c_leith_bi` default 0, `smag_bi_const`
      !! user-settable to 0) still fills `nu4_face ≡ 0` — not a backstop
      !! unless the `nu_4_bg` floor clamp is positive; (2) the flow-aware
      !! face path, when engaged, `return`s before the scalar `nu_4` arm in
      !! `ocean_horizontal_viscosity_compute_tendencies` (the dispatch at
      !! `rdb_ocean_horizontal_viscosity.F90:470-497`), so `nu_4 > 0` is
      !! NOT a backstop whenever the face path is taken — only the scalar
      !! arm's own `nu_4` counts, and only when the face path is not.
      !!
      !! **Binding contract — keep in lock-step with the biharmonic
      !! dispatch** in `ocean_horizontal_viscosity_compute_tendencies`
      !! (`rdb_ocean_horizontal_viscosity.F90:470-497`): any PR that adds a
      !! biharmonic arm to that dispatch must extend this predicate in the
      !! same commit.
      real(wp), intent(in) :: nu_4
      logical, intent(in) :: smag_ah
      real(wp), intent(in) :: smag_bi_const
      integer, intent(in) :: code
      real(wp), intent(in) :: c_leith_bi
      real(wp), intent(in) :: nu_4_bg
      logical :: ok
      logical :: face_path
      real(wp) :: coeff

      ! Mirrors the dispatch guard at :470-473: the flow-aware face path
      ! is taken (and the scalar nu_4 arm below it never runs) whenever
      ! smag_ah is active or the closure is Leith-biharmonic. These two
      ! are already mutually exclusive at configure time
      ! (lateral_closure_conflicts_smag_ah).
      face_path = smag_ah .or. (code == LMIX_LEITH_BIHARM)
      if (face_path) then
         ! nu4_face = clamp(coeff * scale * |strain|, nu_4_bg, nu4_max)
         ! (rdb_ocean_lateral_mix.F90:599/643/743/781) — a positive
         ! coefficient OR a positive floor clamp both guarantee
         ! nu4_face > 0 somewhere.
         coeff = merge(smag_bi_const, c_leith_bi, smag_ah)
         ok = (coeff > 0.0_wp) .or. (nu_4_bg > 0.0_wp)
      else
         ok = (nu_4 > 0.0_wp)
      end if
   end function has_biharmonic_backstop

   pure function parse_lateral_closure(name) result(code)
      !! Translate a namelist string into an `LMIX_*` code.  Returns
      !! `LMIX_INVALID` on an unrecognised value (fail loud at configure);
      !! only "none"/"off"/"" map to `LMIX_NONE`.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("none", "NONE", "off", "OFF", "")
         code = LMIX_NONE
      case ("leith", "LEITH")
         code = LMIX_LEITH
      case ("smag", "SMAG", "smagorinsky", "SMAGORINSKY")
         code = LMIX_SMAGORINSKY
      case ("biharmonic", "BIHARMONIC")
         code = LMIX_BIHARMONIC
      case ("leith_biharm", "LEITH_BIHARM", "leith_biharmonic")
         code = LMIX_LEITH_BIHARM
      case default
         code = LMIX_INVALID
      end select
   end function parse_lateral_closure

   pure function ocean_lateral_mix_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the lateral viscosity slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_lateral_mix_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%ah_face_x) &
               + arr_bytes(this%ah_face_y) &
               + arr_bytes(this%nu4_face_x) &
               + arr_bytes(this%nu4_face_y) &
               + this%vort_corner%bytes()
   end function ocean_lateral_mix_bytes

end module rdb_ocean_lateral_mix
