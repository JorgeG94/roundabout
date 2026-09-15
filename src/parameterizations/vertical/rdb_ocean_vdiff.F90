!! Backward-Euler vertical diffusion solver for the ocean multilayer
!! C-grid.  Two public entry points:
!!
!!   * `vdiff_apply_momentum` — implicit vertical viscosity
!!     on `u_face_x_layer` and `v_face_y_layer`.
!!   * `vdiff_apply_tracers` — implicit vertical diffusivity on every
!!     registered tracer (operates on `T = hTr/h`, writes back
!!     `hTr = T*h`; `h_layer` is untouched).  Takes TWO diffusivity
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
   use rdb_constants, only: wp, H_VANISHED, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: ocean_vdiff_t
   public :: vdiff_apply_momentum
   public :: vdiff_apply_tracers
   public :: face_thick

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
         !! 10.0 m in the double_gyre reference).  Distinct from
         !! `&ocean_bdrag_nml hbbl`, which the implicit-drag fold forces to 0.
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
         !!   * interface kv → `kv + (kv_bbl − kv)·botfn`, with `h_shear`
         !!     capped toward `bbl_thick` under the same botfn weight;
         !!   * bed row → `dt·kv_bbl/(hf₁·(min(hvel₁/2, bbl_thick)))`
         !!     REPLACING the `dt·λ_bot` Rayleigh fold (same drag, stress
         !!     form — h-independent piston instead of h-proportional).
         !! `kv_bbl = bbl_piston·bbl_thick` with `bbl_thick = hbbl_visc`
         !! (v1: no rotational bbl_thick limit — MOM6's f-formula gave
         !! 21.4 m vs our 10 m on the reproducer; the glue magnitude is
         !! h_shear-dominated so the factor 2 is immaterial).  Requires
         !! `hvel_mom6` (supplies the harmonic z bookkeeping),
         !! `implicit_drag`, and `&ocean_bdrag_nml form="linear"` (the
         !! constant-piston parity; quadratic composition deferred) — all
         !! enforced at configure.  Default `.false.` ⇒ bit-identical.
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
   end subroutine ocean_vdiff_exit_data_impl

   subroutine vdiff_apply_momentum(grid, this, ms, dt, kv_source, &
                                   tau_u, tau_v, lambda_bot_u, lambda_bot_v, rho0, &
                                   visc_rem_u, visc_rem_v, remnant_only, &
                                   kv_corner_source, kv_corner_prandtl)
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
      !!     `λ` (1/s) on the respective faces, added to the bed (`k = 1`)
      !!     diagonal as `dt·λ`.  Only read when `implicit_drag`.
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
      real(wp), intent(in), optional :: kv_source(:, :, :)
      real(wp), intent(in), optional :: tau_u(:, :), tau_v(:, :)
      real(wp), intent(in), optional :: lambda_bot_u(:, :), lambda_bot_v(:, :)
      real(wp), intent(in), optional :: rho0
      real(wp), intent(inout), optional :: visc_rem_u(:, :, :), visc_rem_v(:, :, :)
      logical, intent(in), optional :: remnant_only
         !! `.true.` = build the matrices and (re)fill `visc_rem_u/v`
         !! WITHOUT solving for / modifying the velocities — the
         !! pre-substep visc_rem refresh (PGF_BUG.md §9; MOM6 computes
         !! `vertvisc_coef` before `btstep`, so its BT weights never lag).
         !! Requires both `visc_rem_u/v` present.  Default `.false.`.
      real(wp), intent(in), optional :: kv_corner_source(:, :, :)
         !! Corner-staggered interface viscosity source, `(nx+1, ny+1,
         !! nz+1)` — see the corner add-on note above.  Typically the
         !! vertex kappa-shear `kd_corner` carrier (device-resident).
      real(wp), intent(in), optional :: kv_corner_prandtl
         !! Scale applied to `kv_corner_source` (Kv = Pr·Kd; the vertex
         !! kappa-shear `prandtl_turb`).  Default 1.

      integer :: nx, ny, nx_face, ny_uface, nx_vface, ny_face, nz
      logical :: do_stress, do_drag, do_remnant, solve_mom, do_corner
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
      do_remnant = present(visc_rem_u) .and. present(visc_rem_v)
      solve_mom = .true.
      if (present(remnant_only)) solve_mom = .not. remnant_only
      if (.not. solve_mom .and. .not. do_remnant) return   ! nothing to produce
      rho0_l = 1035.0_wp
      if (present(rho0)) rho0_l = rho0
      do_corner = present(kv_corner_source)
      corner_prandtl_l = 1.0_wp
      if (present(kv_corner_prandtl)) corner_prandtl_l = kv_corner_prandtl

      if (present(kv_source)) then
         call diffuse_velocity_columns_impl( &
            nx_face, ny_uface, nz, dt, &
            ms%u_face_x_layer, ms%h_layer, kv_source, ms%wet_mask, &
            .true., nx, ny, this%use_harmonic, &
            this%a_diag_u%data, this%b_diag_u%data, &
            this%c_diag_u%data, this%rhs_u%data, &
            do_stress, do_drag, rho0_l, tau_u, lambda_bot_u, &
            solve_mom, do_remnant, visc_rem_u, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%bbl_piston, this%hvel_upwind, &
            do_corner, corner_prandtl_l, kv_corner_source)
         call diffuse_velocity_columns_impl( &
            nx_vface, ny_face, nz, dt, &
            ms%v_face_y_layer, ms%h_layer, kv_source, ms%wet_mask, &
            .false., nx, ny, this%use_harmonic, &
            this%a_diag_v%data, this%b_diag_v%data, &
            this%c_diag_v%data, this%rhs_v%data, &
            do_stress, do_drag, rho0_l, tau_v, lambda_bot_v, &
            solve_mom, do_remnant, visc_rem_v, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%bbl_piston, this%hvel_upwind, &
            do_corner, corner_prandtl_l, kv_corner_source)
      else
         ! Pure-vdiff no-op short-circuit ONLY when there is also no
         ! boundary forcing to fold; stress/drag/remnant must still be
         ! applied even at K_v = 0 (the surface Ekman + bed sink don't
         ! need interior viscosity to act, and an implicit_drag-only,
         ! K_v=0 config would otherwise silently leave visc_rem stale).
         ! A corner viscosity source likewise keeps the solve alive.
         if (this%K_v_momentum <= 0.0_wp .and. .not. do_stress .and. .not. do_drag &
             .and. .not. do_remnant .and. .not. do_corner) return
         call fill_kv_scalar_buf(this%kv_scalar_buf%data, &
                                 this%K_v_momentum, nx, ny, nz)
         call diffuse_velocity_columns_impl( &
            nx_face, ny_uface, nz, dt, &
            ms%u_face_x_layer, ms%h_layer, this%kv_scalar_buf%data, ms%wet_mask, &
            .true., nx, ny, this%use_harmonic, &
            this%a_diag_u%data, this%b_diag_u%data, &
            this%c_diag_u%data, this%rhs_u%data, &
            do_stress, do_drag, rho0_l, tau_u, lambda_bot_u, &
            solve_mom, do_remnant, visc_rem_u, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%bbl_piston, this%hvel_upwind, &
            do_corner, corner_prandtl_l, kv_corner_source)
         call diffuse_velocity_columns_impl( &
            nx_vface, ny_face, nz, dt, &
            ms%v_face_y_layer, ms%h_layer, this%kv_scalar_buf%data, ms%wet_mask, &
            .false., nx, ny, this%use_harmonic, &
            this%a_diag_v%data, this%b_diag_v%data, &
            this%c_diag_v%data, this%rhs_v%data, &
            do_stress, do_drag, rho0_l, tau_v, lambda_bot_v, &
            solve_mom, do_remnant, visc_rem_v, &
            this%hvel_mom6, this%hbbl_visc, &
            this%bbl_glue, this%bbl_piston, this%hvel_upwind, &
            do_corner, corner_prandtl_l, kv_corner_source)
      end if
   end subroutine vdiff_apply_momentum

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
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%hTr, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%heat_budget_vdiff)
         case (TRACER_BUDGET_SALT)
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%hTr, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%salt_budget_vdiff)
         case default
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%hTr, &
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
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%hTr, &
                                       this%a_diag_t%data, this%b_diag_t%data, &
                                       this%c_diag_t%data, this%rhs_t%data, &
                                       budget=ms%salt_budget_vdiff)
         else
            call apply_factored_tracer(nx, ny, nz, ms%h_layer, ms%tracers(it)%hTr, &
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

   pure subroutine apply_factored_tracer(nx, ny, nz, h_layer, hTr, &
                                         a_diag, b_diag, c_diag, rhs, budget)
      !! Apply the pre-factored tridiagonal (from
      !! `build_factorize_tracer_matrix`) to one tracer: form `T = hTr/h`,
      !! run the Thomas RHS sweep + back-substitution against the stored
      !! pivots, and reconstitute `hTr = T_new * h_layer`.  Operates on
      !! concentration so it conserves cell-centred T; `h_layer` and the
      !! factored `a/b/c` are read-only (shared across every tracer).
      !!
      !! When `budget` is present, accumulate the per-cell increment
      !! (T_new·h − hTr_old) into the contributor slot before overwriting.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: a_diag(nx, ny, nz)
      real(wp), intent(in) :: b_diag(nx, ny, nz)
      real(wp), intent(in) :: c_diag(nx, ny, nz)
      real(wp), intent(inout) :: rhs(nx, ny, nz)
      real(wp), intent(inout), optional :: budget(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: hTr_new

      do concurrent(j=1:ny, i=1:nx) local(k, hTr_new)
         ! RHS = T = hTr/h̃, where h̃ = max(h, H_VANISHED) is the SAME floored
         ! thickness used by build_factorize_tracer_matrix.  A vanishing /
         ! collapsed layer (h ≤ H_VANISHED, e.g. a thin z* surface layer
         ! driven ≤ 0 by combined Fox-Kemper + resolved transport on an
         ! intermediate RK2 stage) is decoupled (identity row in the matrix)
         ! and reconstituted with the same h̃ ⇒ hTr is preserved EXACTLY
         ! rather than zeroed (the previous `h ≤ 0 → T = 0 → hTr = 0` branch
         ! silently dropped its frozen tracer mass → ~5%/day leak under
         ! FK × windowed advect).  For h ≫ H_VANISHED (sigma / double_gyre)
         ! h̃ = h exactly ⇒ bit-identical.
         do k = 1, nz
            rhs(i, j, k) = hTr(i, j, k)/max(h_layer(i, j, k), H_VANISHED)
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

         ! ---- Reconstitute hTr = T_new * h ----
         ! Budget write gated INSIDE the single do concurrent (splitting
         ! present() into two loops makes NVHPC emit a far slower kernel
         ! for one branch — see the remap fix).
         if (present(budget)) then
            do k = 1, nz
               hTr_new = rhs(i, j, k)*max(h_layer(i, j, k), H_VANISHED)
               budget(i, j, k) = budget(i, j, k) + (hTr_new - hTr(i, j, k))
               hTr(i, j, k) = hTr_new
            end do
         else
            do k = 1, nz
               hTr(i, j, k) = rhs(i, j, k)*max(h_layer(i, j, k), H_VANISHED)
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
                                                 solve_momentum, &
                                                 do_remnant, visc_rem_out, &
                                                 hvel_mom6, hbbl_visc, &
                                                 bbl_glue, bbl_piston, hvel_upwind, &
                                                 do_corner, kv_prandtl, kv_corner)
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
      real(wp), intent(inout) :: b_diag(nu, nv, nz)
      real(wp), intent(inout) :: c_diag(nu, nv, nz)
      real(wp), intent(inout) :: rhs(nu, nv, nz)
      logical, intent(in) :: do_stress
         !! Fold the surface wind stress into the `k = nz` RHS row.
      logical, intent(in) :: do_drag
         !! Fold the bottom drag into the `k = 1` diagonal.
      real(wp), intent(in) :: rho0
         !! Boussinesq reference density for the `τ/ρ₀` stress conversion.
      real(wp), intent(in), optional :: tau_face(nu, nv)
         !! Wind stress (N/m²) on this face.  Present iff `do_stress`.
      real(wp), intent(in), optional :: lambda_bot(nu, nv)
         !! Bottom-drag Rayleigh rate λ (1/s) on this face.  Present iff
         !! `do_drag`.
      logical, intent(in) :: hvel_mom6
         !! MOM6 HARMONIC_VISC parity for `h_u` + `h_shear` (see the slot-type
         !! docstring).  `.false.` => the historical arithmetic-`h_u` /
         !! `face_thick`-`dz` pair, bit-identical.
      real(wp), intent(in) :: hbbl_visc
         !! Bottom-layer scale for the `botfn` blend (MOM6 HBBL).
      logical, intent(in) :: bbl_glue
         !! MOM6 `bottomdraglaw` coupling parity (see the slot-type
         !! docstring).  Configure guarantees `hvel_mom6` and `do_drag`
         !! are both on when this is.  `.false.` ⇒ bit-identical.
      real(wp), intent(in) :: bbl_piston
         !! BBL piston velocity u* (m/s); `kv_bbl = bbl_piston·hbbl_visc`.
      logical, intent(in) :: hvel_upwind
         !! `.false.` = skip the near-bed upwind blend (pure harmonic
         !! hvel).  See the slot-type docstring.
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

      integer :: i, j, k, i_left, i_right, j_below, j_above, i_c2, j_c2
      real(wp) :: hf_km1, hf_k, hf_kp1, dz_bot, dz_top
      real(wp) :: hvel(NZ_STACK_MAX)
      real(wp) :: zint(NZ_STACK_MAX)
         !! Normalized height (units of `hbbl_visc`) of the TOP interface
         !! of layer k above the bed, accumulated from HARMONIC face
         !! thicknesses — the mirror of MOM6's `z_i`.
         !! Grounded sliver stacks therefore sit at zint ≈ 0 no matter how
         !! thick their arithmetic-mean face layers are.  Only filled when
         !! `hvel_mom6` (configure guarantees that whenever `bbl_glue`).
      real(wp) :: zacc, z2, botfn, h_harm, h_arith, h_delta, hl_c, hr_c, i_hbbl
      real(wp) :: nu_face_k, nu_face_kp1, alpha, beta, denom
      real(wp) :: botfn_int, kv_bbl, bbl_thick
      real(wp) :: inv_rho0
      real(wp), parameter :: EPS_HVEL = 1.0e-30_wp
         !! MOM6 `h_neglect` analogue in the harmonic mean / HBBL inverse.
      real(wp), parameter :: STRESS_H_MIN = 1.0e-3_wp
         !! Floor on the surface-layer face thickness in the `τ/(ρ₀·h)`
         !! stress conversion — matches `ocean_surface_stress_t%h_min` so
         !! the implicit fold reduces to the explicit kernel in the
         !! well-resolved limit.

      inv_rho0 = 1.0_wp/rho0

      ! BBL-glue constants (MOM6 bottomdraglaw parity, PGF_BUG.md §9.5).
      ! v1 fixes bbl_thick at the botfn scale (no rotational limit) and
      ! forms kv_bbl = u*·bbl_thick as MOM6's linear-drag branch does
      ! (ustar = cdrag_sqrt·DRAG_BG_VEL).
      bbl_thick = hbbl_visc
      kv_bbl = bbl_piston*hbbl_visc

      do concurrent(j=1:nv, i=1:nu) &
         local(i_left, i_right, j_below, j_above, i_c2, j_c2, &
               hf_km1, hf_k, hf_kp1, dz_bot, dz_top, &
               nu_face_k, nu_face_kp1, alpha, beta, denom, k, &
               hvel, zint, zacc, z2, botfn, botfn_int, &
               h_harm, h_arith, h_delta, hl_c, hr_c)
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

         ! Seed RHS with u^n.
         do k = 1, nz
            rhs(i, j, k) = u_face(i, j, k)
         end do

         ! ---- Face thickness per layer (MOM6 `hvel`) ----
         ! k = 1 is the BED here (MOM6 counts from the surface), so the height
         ! above bed accumulates UPWARD and `z2` is the bottom interface of
         ! layer k -- the mirror of MOM6's `z_i(k+1)`.
         i_hbbl = 1.0_wp/(hbbl_visc + EPS_HVEL)
         zacc = 0.0_wp
         do k = 1, nz
            hl_c = h_layer(i_left, j_below, k)
            hr_c = h_layer(i_right, j_above, k)
            if (hvel_mom6) then
               h_harm = 2.0_wp*hl_c*hr_c/(hl_c + hr_c + EPS_HVEL)
               h_arith = 0.5_wp*(hl_c + hr_c)
               h_delta = hr_c - hl_c
               hvel(k) = h_harm
               ! Upwind bias: only when the flow runs from the THICK side to
               ! the THIN side does the near-bed face take the arithmetic
               ! (donor) thickness.  Without this the harmonic mean
               ! over-suppresses asymmetrically -- the reason bare-harmonic
               ! attempts regressed (LAGRANGIAN_PGF_BUG.md 6.2).
               !
               ! `hvel_upwind = .false.` disables the blend (pure harmonic
               ! hvel).  The sign test keys on u itself, so at (near-)rest it
               ! flip-flops faces between harmonic and arithmetic on
               ! ROUNDOFF-sign velocities, collapsing the BBL glue at
               ! whichever faces flip — measured ×200-800/stage residual
               ! amplification on the rest-state reproducer (PGF_BUG.md
               ! §9.8).  MOM6 carries the same test and the same hazard; it
               ! never bites there only because its state stays at 1e-14.
               if (hvel_upwind .and. u_face(i, j, k)*h_delta < 0.0_wp) then
                  z2 = zacc
                  botfn = 1.0_wp/(1.0_wp + 0.09_wp*z2*z2*z2*z2*z2*z2)
                  hvel(k) = (1.0_wp - botfn)*h_harm + botfn*h_arith
               end if
               zacc = zacc + h_harm*i_hbbl
               zint(k) = zacc
            else
               hvel(k) = 0.5_wp*(hl_c + hr_c)
            end if
         end do

         ! ---- k = 1: bed BC, no flux below ----
         ! Floor the face thicknesses at H_VANISHED before they enter the
         ! α/β denominators (dt·ν/(hf_k·dz)).  For an exactly-collapsed layer
         ! (both neighbour cells h=0 ⇒ hvel=0) with kv>0, a raw hf_k=0 (and the
         ! face_thick-derived dz=0) would make α=Inf ⇒ b_diag=Inf ⇒ NaN γ.
         ! Mirrors the tracer path's max(h, H_VANISHED) floor.  dz_top/dz_bot
         ! are built FROM these floored hf_* so they inherit the floor.  In
         ! every valid config hvel ≫ H_VANISHED ⇒ the max() is a no-op ⇒
         ! bit-identical.  Floor at the READ site (not by reassigning hvel(k))
         ! — gfortran mis-optimizes a local() array element reassigned across
         ! branches.
         hf_k = max(hvel(1), H_VANISHED)
         ! SINGLE-LAYER COLUMN (nz = 1).  The bed row IS the surface row:
         ! there is no interior interface above it, so the interior
         ! coupling is identically zero and the surface is a pure
         ! stress-Neumann BC (a RHS source, added below with the k = nz
         ! block's).  Without this gate the code reads `hvel(2)` — one
         ! past the layer loop that filled `hvel`, i.e. uninitialised
         ! `local()` stack — and the k = nz block below then OVERWRITES
         ! this row from `hvel(0)`/`zint(0)`, an out-of-bounds read that
         ! is a DETERMINISTIC `CUDA_ERROR_ILLEGAL_ADDRESS` on the NVHPC
         ! GPU build and silent stack garbage on the host.  Gated INSIDE
         ! the single `do concurrent` (a split loop costs an extra launch
         ! — see CLAUDE.md) and loop-invariant, so nz >= 2 is
         ! bit-identical.  Written as init-then-conditional-overwrite, not
         ! if/else: gfortran 15.1 miscompiles a `local()` scalar reassigned
         ! across the two arms of an if/else inside `do concurrent`.
         alpha = 0.0_wp
         if (nz > 1) then
            hf_kp1 = max(hvel(2), H_VANISHED)
            if (hvel_mom6) then
               dz_top = 0.5_wp*(hf_k + hf_kp1)   ! MOM6 h_shear: arithmetic of hvels
            else
               dz_top = face_thick(hf_k, hf_kp1, use_harmonic)
            end if
            nu_face_kp1 = 0.5_wp*(kv_centre(i_left, j_below, 2) + kv_centre(i_right, j_above, 2))
            ! Corner-viscosity add-on (vertex kappa-shear Kv seam): direct
            ! 2-point end-corner average onto this face, BEFORE the BBL
            ! glue (all viscosity contributions fold ahead of the
            ! coupling-coefficient blend).
            if (do_corner) then
               nu_face_kp1 = nu_face_kp1 + &
                             kv_prandtl*(0.5_wp*(kv_corner(i, j, 2) + kv_corner(i_c2, j_c2, 2)))
            end if
            ! BBL glue at the interface above the bed layer (MOM6
            ! find_coupling_coef, vert_friction:2214-2229): within botfn reach
            ! of the bed, the viscosity rises to kv_bbl and the shear distance
            ! is capped toward bbl_thick — grounded stacks (zint ≈ 0) become
            ! rigidly coupled, which is the mechanism that absorbs the
            ! spurious grounded-layer PGF every step (PGF_BUG.md §9).
            if (bbl_glue) then
               botfn_int = 1.0_wp/(1.0_wp + 0.09_wp*zint(1)*zint(1)*zint(1)* &
                                   zint(1)*zint(1)*zint(1))
               nu_face_kp1 = nu_face_kp1 + (kv_bbl - nu_face_kp1)*botfn_int
               if (dz_top > bbl_thick) then
                  dz_top = (1.0_wp - botfn_int)*dz_top + botfn_int*bbl_thick
               end if
            end if
            alpha = dt*nu_face_kp1/(hf_k*dz_top)
         end if
         a_diag(i, j, 1) = 0.0_wp
         c_diag(i, j, 1) = -alpha
         b_diag(i, j, 1) = 1.0_wp + alpha
         ! Bottom-drag stress BC (Roundabout bed = k=1; MIRROR of MOM6 k=nz).
         ! λ_bot is the Rayleigh RATE (c_d·|U|/h_1 quadratic, r linear) the
         ! bottom-drag slot already forms — the row is pre-normalized by
         ! h_1, so the MOM6 `h_1 + dt·a_bot` diagonal becomes `1 + dt·λ_bot`
         ! here.  A drag is a SINK ⇒ POSITIVE diagonal add ⇒ |amplification|
         ! ≤ 1 for ANY h/dt (unconditionally stable on thin shelf bottoms).
         ! Gated INSIDE the single DC (no split loop — NVHPC penalty); the
         ! optional `lambda_bot` is guaranteed present when `do_drag`.
         if (do_drag) then
            if (bbl_glue) then
               ! MOM6 bed coupling: the drag is a
               ! viscous PISTON `kv_bbl/(min(hvel₁/2, bbl_thick))` — it
               ! DIVERGES as the bottom layer thins (a sliver bed layer is
               ! anchored rigidly to rest), where the Rayleigh `dt·λ` fold
               ! below is h-independent and lets grounded stacks reach the
               ! ballistic balance u_eq = PGF·h/r (PGF_BUG.md §9.4).  Same
               ! drag physics for resolved columns: at h₁ = bbl_thick this
               ! equals dt·(bbl_piston/h₁), the distributed-linear-drag
               ! rate.  Replaces (not augments) the λ fold — one bed sink.
               b_diag(i, j, 1) = b_diag(i, j, 1) + dt*kv_bbl/ &
                                 (hf_k*(min(0.5_wp*hvel(1), bbl_thick) + EPS_HVEL))
            else
               b_diag(i, j, 1) = b_diag(i, j, 1) + dt*lambda_bot(i, j)
            end if
         end if

         ! ---- k = 2..nz-1 ----
         do k = 2, nz - 1
            ! Floor the face thicknesses (see the k=1 block) — keeps the α/β
            ! denominators non-zero for a collapsed interior layer; no-op for
            ! hvel ≫ H_VANISHED ⇒ bit-identical.
            hf_km1 = max(hvel(k - 1), H_VANISHED)
            hf_k = max(hvel(k), H_VANISHED)
            hf_kp1 = max(hvel(k + 1), H_VANISHED)
            if (hvel_mom6) then
               dz_bot = 0.5_wp*(hf_km1 + hf_k)
               dz_top = 0.5_wp*(hf_k + hf_kp1)
            else
               dz_bot = face_thick(hf_km1, hf_k, use_harmonic)
               dz_top = face_thick(hf_k, hf_kp1, use_harmonic)
            end if
            nu_face_k = 0.5_wp*(kv_centre(i_left, j_below, k) + kv_centre(i_right, j_above, k))
            nu_face_kp1 = 0.5_wp*(kv_centre(i_left, j_below, k + 1) + &
                                  kv_centre(i_right, j_above, k + 1))
            ! Corner-viscosity add-on (see the k=1 block).
            if (do_corner) then
               nu_face_k = nu_face_k + &
                           kv_prandtl*(0.5_wp*(kv_corner(i, j, k) + kv_corner(i_c2, j_c2, k)))
               nu_face_kp1 = nu_face_kp1 + &
                             kv_prandtl*(0.5_wp*(kv_corner(i, j, k + 1) + &
                                                 kv_corner(i_c2, j_c2, k + 1)))
            end if
            ! BBL glue (see the k=1 block).  The (nu, dz) transform uses the
            ! interface's own zint, so the SAME effective coupling lands in
            ! this row's alpha and the row-above's beta (symmetric matrix).
            if (bbl_glue) then
               botfn_int = 1.0_wp/(1.0_wp + 0.09_wp*zint(k - 1)*zint(k - 1)*zint(k - 1)* &
                                   zint(k - 1)*zint(k - 1)*zint(k - 1))
               nu_face_k = nu_face_k + (kv_bbl - nu_face_k)*botfn_int
               if (dz_bot > bbl_thick) then
                  dz_bot = (1.0_wp - botfn_int)*dz_bot + botfn_int*bbl_thick
               end if
               botfn_int = 1.0_wp/(1.0_wp + 0.09_wp*zint(k)*zint(k)*zint(k)* &
                                   zint(k)*zint(k)*zint(k))
               nu_face_kp1 = nu_face_kp1 + (kv_bbl - nu_face_kp1)*botfn_int
               if (dz_top > bbl_thick) then
                  dz_top = (1.0_wp - botfn_int)*dz_top + botfn_int*bbl_thick
               end if
            end if
            beta = dt*nu_face_k/(hf_k*dz_bot)
            alpha = dt*nu_face_kp1/(hf_k*dz_top)
            a_diag(i, j, k) = -beta
            c_diag(i, j, k) = -alpha
            b_diag(i, j, k) = 1.0_wp + alpha + beta
         end do

         ! ---- k = nz: surface BC, no flux above ----
         ! Floor the face thicknesses (see the k=1 block) — keeps the β
         ! denominator non-zero for a collapsed surface layer; no-op for
         ! hvel ≫ H_VANISHED ⇒ bit-identical.  (The stress-BC conversion
         ! below keeps its own STRESS_H_MIN floor on hf_k.)
         hf_k = max(hvel(nz), H_VANISHED)
         ! SINGLE-LAYER COLUMN (nz = 1): this row has already been built as
         ! the bed row above (a = c = 0, b = 1 + bottom drag).  Skip the
         ! interface-below build entirely — at nz = 1 it would read
         ! `hvel(0)` / `zint(0)` (out of bounds) and overwrite the bed
         ! row's drag.  Only the stress RHS below still applies, which is
         ! exactly right: one layer carries BOTH the wind stress and the
         ! bottom drag.  Loop-invariant gate ⇒ nz >= 2 bit-identical.
         if (nz > 1) then
            hf_km1 = max(hvel(nz - 1), H_VANISHED)
            if (hvel_mom6) then
               dz_bot = 0.5_wp*(hf_km1 + hf_k)
            else
               dz_bot = face_thick(hf_km1, hf_k, use_harmonic)
            end if
            nu_face_k = 0.5_wp*(kv_centre(i_left, j_below, nz) + kv_centre(i_right, j_above, nz))
            ! Corner-viscosity add-on (see the k=1 block).
            if (do_corner) then
               nu_face_k = nu_face_k + &
                           kv_prandtl*(0.5_wp*(kv_corner(i, j, nz) + kv_corner(i_c2, j_c2, nz)))
            end if
            ! BBL glue (see the k=1 block) — this row's beta is the twin of
            ! row nz-1's alpha and must see the same transform.
            if (bbl_glue) then
               botfn_int = 1.0_wp/(1.0_wp + 0.09_wp*zint(nz - 1)*zint(nz - 1)*zint(nz - 1)* &
                                   zint(nz - 1)*zint(nz - 1)*zint(nz - 1))
               nu_face_k = nu_face_k + (kv_bbl - nu_face_k)*botfn_int
               if (dz_bot > bbl_thick) then
                  dz_bot = (1.0_wp - botfn_int)*dz_bot + botfn_int*bbl_thick
               end if
            end if
            beta = dt*nu_face_k/(hf_k*dz_bot)
            a_diag(i, j, nz) = -beta
            c_diag(i, j, nz) = 0.0_wp
            b_diag(i, j, nz) = 1.0_wp + beta
         end if
         ! Surface wind-stress Neumann BC (Roundabout surface = k=nz; MIRROR of
         ! MOM6 k=1).  The free-surface flux is the prescribed kinematic
         ! stress τ/ρ₀, a RHS source — the DIAGONAL is unchanged.  The row
         ! is pre-normalized by h_nz (= hf_k here), so the MOM6 surface_stress
         ! term `dt·τ/ρ₀` becomes `dt·(τ/ρ₀)/h_nz_face`.  In the zero-
         ! interior-coupling limit this reduces EXACTLY to the explicit
         ! `du = dt·τ/(ρ₀·h_top)` increment.  `STRESS_H_MIN` matches the
         ! explicit surface-stress kernel's `h_min` so the two paths agree
         ! in the well-resolved limit.  `tau_face` is masked by the face
         ! `min` of the two bounding cells' wet/dry mask (as the explicit
         ! kernel does) so a land face receives no stress.  Gated INSIDE
         ! the single DC.
         if (do_stress) then
            rhs(i, j, nz) = rhs(i, j, nz) + &
                            dt*tau_face(i, j) &
                            *min(wet_cell(i_left, j_below), wet_cell(i_right, j_above)) &
                            *inv_rho0/max(hf_k, STRESS_H_MIN)
         end if

         ! ---- Thomas forward sweep ----
         ! remnant-only mode (`solve_momentum = .false.`) runs the c'
         ! recurrence WITHOUT the rhs leg and skips the back-substitution
         ! — u_face is untouched, but c_diag ends up holding exactly the
         ! c' values the remnant solve below needs.
         c_diag(i, j, 1) = c_diag(i, j, 1)/b_diag(i, j, 1)
         if (solve_momentum) then
            rhs(i, j, 1) = rhs(i, j, 1)/b_diag(i, j, 1)
            do k = 2, nz
               denom = b_diag(i, j, k) - a_diag(i, j, k)*c_diag(i, j, k - 1)
               c_diag(i, j, k) = c_diag(i, j, k)/denom
               rhs(i, j, k) = (rhs(i, j, k) - a_diag(i, j, k)*rhs(i, j, k - 1))/denom
            end do

            ! ---- Back-substitution: write directly into u_face ----
            u_face(i, j, nz) = rhs(i, j, nz)
            do k = nz - 1, 1, -1
               u_face(i, j, k) = rhs(i, j, k) - c_diag(i, j, k)*u_face(i, j, k + 1)
            end do
         else
            do k = 2, nz
               denom = b_diag(i, j, k) - a_diag(i, j, k)*c_diag(i, j, k - 1)
               c_diag(i, j, k) = c_diag(i, j, k)/denom
            end do
         end if

         ! ---- Viscous remnant γ_k: second solve, same operator, RHS ≡ 1 ----
         ! γ is DEFINED as the momentum solve's own sensitivity to a uniform
         ! barotropic acceleration, so it is built against the SAME
         ! factorization rather than a freshly-assembled matrix: `a_diag`/
         ! `b_diag` are untouched by the momentum forward sweep above,
         ! and `c_diag` already holds c' (the sweep overwrote it in place),
         ! so `denom = b_diag(k) - a_diag(k)*c_diag(k-1)` recomputes the
         ! IDENTICAL value the momentum sweep used at this k (bit-identical
         ! denominators — no new scratch, no re-derivation from h/kv).
         ! RHS is the unit vector, NOT h_u(k) — Roundabout's rows are already
         ! normalized by h_k (MOM6's un-normalized rows are why its RHS is
         ! h_u(k); dividing MOM6's h_1 remnant RHS by h_1 gives exactly 1).
         if (do_remnant) then
            visc_rem_out(i, j, 1) = 1.0_wp/b_diag(i, j, 1)
            do k = 2, nz
               denom = b_diag(i, j, k) - a_diag(i, j, k)*c_diag(i, j, k - 1)
               visc_rem_out(i, j, k) = (1.0_wp - a_diag(i, j, k)*visc_rem_out(i, j, k - 1))/denom
            end do
            ! Clamp at production (MOM6 clamps at
            ! consumption; the maximum principle bounds γ to (0, 1] so this
            ! is cheap FP insurance, applied as each element is finalized).
            visc_rem_out(i, j, nz) = min(visc_rem_out(i, j, nz), 1.0_wp)
            do k = nz - 1, 1, -1
               visc_rem_out(i, j, k) = visc_rem_out(i, j, k) &
                                       - c_diag(i, j, k)*visc_rem_out(i, j, k + 1)
               visc_rem_out(i, j, k) = min(visc_rem_out(i, j, k), 1.0_wp)
            end do
         end if
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

   pure function ocean_vdiff_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the implicit vertical diffusion slot (0 when
      !! unallocated).
      class(ocean_vdiff_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%a_diag_t%bytes() &
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
               + this%kv_scalar_buf%bytes()
   end function ocean_vdiff_bytes

end module rdb_ocean_vdiff
