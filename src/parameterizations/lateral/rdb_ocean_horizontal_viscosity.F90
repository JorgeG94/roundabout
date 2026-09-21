!! Constant-coefficient Laplacian horizontal viscosity on the ocean
!! C-grid.  Baseline closure for the dynamical core — gives a flat
!! `nu_h` momentum diffusion that damps grid-scale gravity-wave noise
!! and keeps the unsplit driver bounded.  Phase 5a / 5b will plug
!! Leith and Smagorinsky scaled `nu_h` fields into the same compute
!! kernel — the apply step is the same forward-Euler accumulation
!! regardless of closure.
module rdb_ocean_horizontal_viscosity
   !! Carries the closure parameters and per-step tendency workspace
   !! for the Laplacian horizontal-momentum-viscosity kernel.  Sits
   !! alongside `rdb_ocean_lateral_mix` (which holds the
   !! variable-coefficient state for Leith/Smag closures): this module
   !! is the *kernel*, that one is the *closure*.  Phase Tier-1 ships
   !! the constant-`nu_h` variant; later phases will read coefficients
   !! out of `ocean_lateral_mix_t%ah_face_*` and replace the scalar.
   !!
   !! Each compute pass writes a per-face tendency
   !! (`du_visc`, `dv_visc`); the apply step adds `dt * tendency` onto
   !! `u_face_x_layer` / `v_face_y_layer`.  The two-step pattern
   !! matches Coriolis and PGF — additive tendency buffers keep the
   !! SSP-RK2 driver simple (order between applies doesn't matter).
   !!
   !! Wall handling: at the four C-grid wall faces the tendency is
   !! forced to zero.  Interior faces use the standard 5-point
   !! Laplacian stencil on the face velocity itself (no thickness
   !! weighting — Phase 5+ adds the `h * A * grad u` flux-form once
   !! variable-thickness conservation matters).
   !!
   !! MOM6 stress-divergence path (`stress_tensor = .true.`, spec PR2):
   !! instead of the velocity Laplacian the kernel assembles a
   !! thickness-weighted stress and takes its divergence:
   !!   tension   `str_xx = A_T·(du/dx − dv/dy)·h_T`   (T-cell)
   !!   shear     `str_xy = A_q·(dv/dx + du/dy)·h_q·slip` (Bu corner)
   !!   `diffu = (1/(h_u + h_neglect))·iareaCu·∂(str)`,  `h_neglect`
   !!   an `H_VANISHED`-class floor.  Momentum-conserving and
   !!   down-weights vanishing layers.  A per-cell CFL viscosity
   !!   limiter (MOM6 BOUND_KH) clamps `A` from the actual discrete
   !!   stencil + `dt`, replacing the global `ah_max` cap; `wet_u`,
   !!   `wet_v`, `wet_q` mask the stress so momentum is not diffused
   !!   across coastlines.  On a uniform-grid + uniform-h + all-wet
   !!   column the cross terms cancel discretely and the operator
   !!   reduces to `A·∇²u` to round-off.
   !!
   !! Biharmonic composition: the constant-`nu_4` / flow-aware
   !! (`smag_ah`, `leith_biharm`) biharmonic add-on composes with
   !! **all three** harmonic operators — scalar Laplacian, face
   !! (Leith/Smagorinsky) Laplacian, and the stress-divergence path —
   !! matching MOM6 (`BIHARMONIC` "may be used with `LAPLACIAN`",
   !! default `.true.`).  `kh_aniso` therefore also composes with the
   !! biharmonic; it is no longer mutually exclusive.  The composition
   !! happens at the *tendency* level, not the *stress* level: under
   !! `stress_tensor`, the harmonic part is momentum-conserving and
   !! coast-masked while the velocity-form biharmonic part is neither
   !! (a documented fidelity divergence from MOM6, which sums both into
   !! one stress tensor before differencing — porting the biharmonic
   !! into the stress tensor is a follow-on PR, not this module today).
   use rdb_constants, only: wp, H_VANISHED, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, LMIX_NONE, &
                                    LMIX_LEITH, LMIX_SMAGORINSKY, LMIX_LEITH_BIHARM
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_horizontal_viscosity_t
   public :: ocean_horizontal_viscosity_compute_tendencies
   public :: ocean_horizontal_viscosity_apply_tendencies
   public :: ocean_horizontal_viscosity_compute_ke_diss
   public :: ocean_hvisc_set_aniso_direction
   public :: aniso_mode_is_implemented

   type :: ocean_horizontal_viscosity_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Tracks GPU device
         !! attachment too — prefer this to `allocated(...)`.
      real(wp) :: nu_h = 0.0_wp
         !! Constant Laplacian horizontal viscosity (m^2/s).  Phase
         !! Tier-1 default is zero (kernel becomes a no-op); set
         !! positive to enable damping.  Numerical-stability cap:
         !! `nu_h * dt * (1/dx^2 + 1/dy^2) <= 0.5`.
      logical :: stress_tensor = .false.
         !! When true, the compute step uses the MOM6-faithful
         !! thickness-weighted stress-divergence operator
         !! `diffu = (1/(h_u+h_neglect))·∇·(h·A·∇u)` with a per-cell
         !! CFL viscosity limiter and `wet_*` coast-masking instead of
         !! the velocity Laplacian `A·∇²u`.  Default false ⇒
         !! bit-identical to the historical kernel.  See the module
         !! header for the operator form.
      logical :: no_slip = .false.
         !! Coastal lateral BC selector (shared with the lateral-mix /
         !! Coriolis kernels).  Only consulted on the `stress_tensor`
         !! path: `.false.` (free-slip) masks the corner shear stress
         !! by `wet_q`; `.true.` (no-slip) by `2 - wet_q`.  All-wet ⇒
         !! factor ≡ 1 ⇒ bit-identical.
      real(wp) :: bound_coef = 0.8_wp
         !! CFL safety coefficient for the per-cell viscosity limiter
         !! (MOM6 `HORVISC_BOUND_COEF`).  Consulted on the
         !! `stress_tensor` path and, when `bound_kh` is set, on the
         !! velocity-Laplacian paths too.
      logical :: bound_kh = .false.
         !! MOM6 `BOUND_KH` analogue for the velocity-Laplacian paths
         !! (scalar `nu_h` and the flow-aware per-face closure).  When
         !! `.true.`, the per-face harmonic viscosity is clamped to
         !! `bound_coef·0.125/(dt·(idx²+idy²))` — ~0.25× the forward-
         !! Euler stability limit, matching MOM6's `Kh_Max_xx` margin.
         !! Load-bearing beyond simple FE stability: the barotropic
         !! mode receives the depth-mean viscous force FROZEN over the
         !! outer step (via `F_bt`), and for grid-scale gravity modes
         !! with `ω·dt ≳ 1` an unbounded `λ·dt = ν·k²·dt ≳ 0.6` frozen
         !! force is applied with reversed phase — ANTI-damping — which
         !! exponentially pumps rim-trapped barotropic modes through
         !! the Δu corrector (the 600² Lagrangian double-gyre h-guard
         !! blow-up; e-fold ~10 outer steps).  The clamp keeps
         !! `λ·dt ≤ 0.5·bound_coef` at every face so the corrector
         !! loop stays damped at any resolution.  Default `.false.`
         !! ⇒ bit-identical.

      ! ---- Anisotropic viscosity (Smith & McWilliams 2003) ----
      real(wp) :: kh_aniso = 0.0_wp
         !! Anisotropic Laplacian viscosity magnitude (m²/s).  When
         !! positive, a two-coefficient direction tensor splits the
         !! harmonic viscosity into along- and cross-direction parts.
         !! Only consulted on the `stress_tensor` path.  Default 0 ⇒
         !! isotropic, bit-identical.  MOM6 `KH_ANISO` analogue
         !! (Smith & McWilliams 2003, "Anisotropic horizontal
         !! viscosity for ocean models", Ocean Modelling 5(2), §2).
         !! Composes with the biharmonic add-on (`nu_4` / `smag_ah` /
         !! `leith_biharm`) — the two are independent linear operators
         !! that sum, matching MOM6's composition (`kh_aniso` inside
         !! the harmonic block, biharmonic added to the same tensor).
      real(wp) :: aniso_n1n2 = 0.0_wp
         !! Precomputed direction-tensor factor `2·n1·n2/(n1²+n2²)`
         !! (MOM6 `n1n2`).  Set by `ocean_hvisc_set_aniso_direction`
         !! from the `(n1,n2)` direction vector.  Default `(1,0)` =
         !! grid-i ⇒ `n1n2 = 0` (cross terms vanish; the operator just
         !! adds `kh_aniso` to the tension coefficient).
      real(wp) :: aniso_n1n1_m_n2n2 = 1.0_wp
         !! Precomputed direction-tensor factor `(n1²−n2²)/(n1²+n2²)`
         !! (MOM6 `n1n1_m_n2n2`).  Default `(1,0)` ⇒ 1.
      real(wp) :: nu_4 = 0.0_wp
         !! Constant biharmonic horizontal viscosity (m^4/s).  Scale-
         !! selective damping for stratified closed-basin runs: damps
         !! proportional to ν₄·k⁴, so it kills grid-scale baroclinic
         !! noise without bleeding into resolved scales the way a
         !! large Laplacian ν_h would.  Required to keep stratified
         !! Tasman-class runs bounded past day ~10 (Laplacian-only
         !! configurations grow exponentially via parametric
         !! amplification of roundoff seeds).  Numerical-stability
         !! cap: `ν₄ · dt · (1/dx² + 1/dy²)² <= 1/16`.  Typical 2 km
         !! resolution value ~1e10 m^4/s.  Composes with **all** of
         !! the harmonic dispatch arms, including `stress_tensor` — it
         !! is no longer silently disabled by it.

      type(scratch_3d_buffer_t) :: du_visc
         !! Per-step viscous tendency at east faces, shape
         !! (nx+1, ny, nz).  Filled by the compute step, consumed
         !! by the apply step.
      type(scratch_3d_buffer_t) :: dv_visc
         !! Per-step viscous tendency at north faces, shape
         !! (nx, ny+1, nz).
      type(scratch_3d_buffer_t) :: lap_u
         !! First-pass Laplacian buffer used by the biharmonic path —
         !! holds `∇²u_face` so the second Laplacian pass (`∇²(∇²u)`)
         !! can read it.  Same shape as `du_visc`.  Allocated
         !! unconditionally; sits inert when `nu_4 = 0`.
      type(scratch_3d_buffer_t) :: lap_v
         !! v-face counterpart of `lap_u`.

      ! ---- MOM6 stress-divergence scratch (stress_tensor path only) ----
      type(scratch_3d_buffer_t) :: str_xx
         !! Thickness-weighted tension stress `A_T·(du/dx−dv/dy)·h_T`
         !! at T-cell centres, shape `(nx, ny, nz)`.
      type(scratch_3d_buffer_t) :: str_xy
         !! Thickness-weighted shear stress
         !! `A_q·(dv/dx+du/dy)·h_q·slip` at Bu corners, shape
         !! `(nx+1, ny+1, nz)`.
      type(scratch_3d_buffer_t) :: ah_t
         !! Harmonic viscosity averaged onto T-cell centres (m²/s),
         !! shape `(nx, ny, nz)`.  CFL-clamped per cell.
      type(scratch_3d_buffer_t) :: ah_q
         !! Harmonic viscosity averaged onto Bu corners (m²/s),
         !! shape `(nx+1, ny+1, nz)`.  CFL-clamped per corner.

      ! ---- MEKE frictional-source coupling (capability [5] seam) ----
      logical :: compute_ke_diss = .false.
         !! When set (by configure when MEKE's frictional source is on),
         !! the apply step fills `ke_diss` with the lateral-viscosity KE
         !! dissipation rate.  Default off ⇒ no extra work, bit-identical.
      real(wp), allocatable :: ke_diss(:, :)
         !! Depth-integrated KE dissipation rate by the lateral viscosity,
         !! `Σ_k ρ_k h_k (u·du_visc + v·dv_visc)` (kg/s³, ≤0), at T-cell
         !! centres `(nx, ny)`.  MEKE consumes it as `-frcoeff·i_mass·ke_diss`
         !! (the mean→eddy frictional source).  Holds the most recent stage's
         !! rate (the quantity MEKE needs is a rate, so no accumulation).
   contains
      procedure, non_overridable :: init => ocean_hvisc_init
      procedure, non_overridable :: destroy => ocean_hvisc_destroy
      procedure, non_overridable :: enter_data => ocean_hvisc_enter_data
      procedure, non_overridable :: exit_data => ocean_hvisc_exit_data
      procedure, non_overridable :: bytes => ocean_horizontal_viscosity_bytes
   end type ocean_horizontal_viscosity_t

contains

   subroutine ocean_hvisc_init(this, grid, nz_ml)
      !! Allocate the two tendency scratch buffers.  Same
      !! optional-`nz_ml` pattern as the other ocean kernels — default
      !! 1 keeps the barotropic-only constructor valid; pass `nz_ml`
      !! to size for the multilayer driver.
      class(ocean_horizontal_viscosity_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      call this%du_visc%init(nx + 1, ny, nz, "ocean_hvisc_du_visc")
      call this%dv_visc%init(nx, ny + 1, nz, "ocean_hvisc_dv_visc")
      call this%lap_u%init(nx + 1, ny, nz, "ocean_hvisc_lap_u")
      call this%lap_v%init(nx, ny + 1, nz, "ocean_hvisc_lap_v")
      call this%str_xx%init(nx, ny, nz, "ocean_hvisc_str_xx")
      call this%str_xy%init(nx + 1, ny + 1, nz, "ocean_hvisc_str_xy")
      call this%ah_t%init(nx, ny, nz, "ocean_hvisc_ah_t")
      call this%ah_q%init(nx + 1, ny + 1, nz, "ocean_hvisc_ah_q")
      allocate (this%ke_diss(nx, ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_hvisc_init

   subroutine ocean_hvisc_destroy(this)
      class(ocean_horizontal_viscosity_t), intent(inout) :: this
      this%is_init = .false.
      call this%du_visc%destroy()
      call this%dv_visc%destroy()
      call this%lap_u%destroy()
      call this%lap_v%destroy()
      call this%str_xx%destroy()
      call this%str_xy%destroy()
      call this%ah_t%destroy()
      call this%ah_q%destroy()
      if (allocated(this%ke_diss)) deallocate (this%ke_diss)
   end subroutine ocean_hvisc_destroy

   pure subroutine ocean_hvisc_set_aniso_direction(this, n1, n2)
      !! Precompute the constant Smith & McWilliams (2003) direction-
      !! tensor factors from the anisotropy direction vector `(n1,n2)`
      !! (grid-relative i,j components).  Normalises by `n1²+n2²` so the
      !! caller need not pass a unit vector:
      !!
      !!     n1n2          = 2·n1·n2 / (n1²+n2²)
      !!     n1n1_m_n2n2   = (n1²−n2²) / (n1²+n2²)
      !!
      !! The default grid-i direction `(1,0)` gives `n1n2 = 0`,
      !! `n1n1_m_n2n2 = 1` — the cross terms vanish and only the tension
      !! coefficient picks up `kh_aniso`.  A degenerate `(0,0)` vector
      !! leaves the factors at their defaults (`n1n2=0`,
      !! `n1n1_m_n2n2=1`), i.e. grid-i.
      class(ocean_horizontal_viscosity_t), intent(inout) :: this
      real(wp), intent(in) :: n1, n2
      real(wp) :: recip_norm
      recip_norm = (n1*n1) + (n2*n2)
      if (recip_norm > 0.0_wp) then
         recip_norm = 1.0_wp/recip_norm
         this%aniso_n1n2 = 2.0_wp*(n1*n2)*recip_norm
         this%aniso_n1n1_m_n2n2 = ((n1*n1) - (n2*n2))*recip_norm
      else
         this%aniso_n1n2 = 0.0_wp
         this%aniso_n1n1_m_n2n2 = 1.0_wp
      end if
   end subroutine ocean_hvisc_set_aniso_direction

   pure function aniso_mode_is_implemented(mode) result(ok)
      !! `.true.` iff the anisotropy-direction mode has an implemented
      !! direction tensor.  Only mode 0 (grid-relative `(n1,n2) = aniso_dir`,
      !! a constant tensor) is built; the MOM6 flow-aligned modes are not
      !! ported.  Drives the configure-time fail-loud guard in
      !! `validate_config` so a requested-but-unimplemented mode aborts the
      !! run instead of silently falling back to the grid-i default.
      integer, intent(in) :: mode
      logical :: ok
      ok = (mode == 0)
   end function aniso_mode_is_implemented

   subroutine ocean_hvisc_enter_data(this)
      ! Bare `copyin(this)` removed (stack-descriptor map → AMD cross-slot
      ! overlap; see ocean_surfstress_enter_data).  Scratch buffers attach
      ! below; slot-descriptor presence comes from the root copyin(state)
      ! in ocean_state_enter_data.
      class(ocean_horizontal_viscosity_t), intent(inout) :: this
      select type (this)
      type is (ocean_horizontal_viscosity_t)
         call ocean_hvisc_enter_data_impl(this)
      end select
   end subroutine ocean_hvisc_enter_data

   subroutine ocean_hvisc_enter_data_impl(this)
      type(ocean_horizontal_viscosity_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%du_visc)
      call scratch_3d_buffer_enter_data_impl(this%dv_visc)
      call scratch_3d_buffer_enter_data_impl(this%lap_u)
      call scratch_3d_buffer_enter_data_impl(this%lap_v)
      call scratch_3d_buffer_enter_data_impl(this%str_xx)
      call scratch_3d_buffer_enter_data_impl(this%str_xy)
      call scratch_3d_buffer_enter_data_impl(this%ah_t)
      call scratch_3d_buffer_enter_data_impl(this%ah_q)
      !$acc enter data copyin(this%ke_diss)
   end subroutine ocean_hvisc_enter_data_impl

   subroutine ocean_hvisc_exit_data(this)
      class(ocean_horizontal_viscosity_t), intent(inout) :: this
      select type (this)
      type is (ocean_horizontal_viscosity_t)
         call ocean_hvisc_exit_data_impl(this)
      end select
   end subroutine ocean_hvisc_exit_data

   subroutine ocean_hvisc_exit_data_impl(this)
      type(ocean_horizontal_viscosity_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%du_visc)
      call scratch_3d_buffer_exit_data_impl(this%dv_visc)
      call scratch_3d_buffer_exit_data_impl(this%lap_u)
      call scratch_3d_buffer_exit_data_impl(this%lap_v)
      call scratch_3d_buffer_exit_data_impl(this%str_xx)
      call scratch_3d_buffer_exit_data_impl(this%str_xy)
      call scratch_3d_buffer_exit_data_impl(this%ah_t)
      call scratch_3d_buffer_exit_data_impl(this%ah_q)
      !$acc exit data delete(this%ke_diss)
   end subroutine ocean_hvisc_exit_data_impl

   subroutine ocean_horizontal_viscosity_compute_tendencies(grid, metrics, this, ms, lateral_mix, dt, u_src, v_src, h_src)
      !! Source-selecting shim over `ocean_horizontal_viscosity_compute_tendencies_on`:
      !! absent `u_src`/`v_src`/`h_src` (all callers today) forwards the
      !! prognostic components — bit-identical; the pred_corr driver passes
      !! the `u_av` time-mean family (SPEC §4 S3).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_horizontal_viscosity_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_lateral_mix_t), intent(in), optional :: lateral_mix
      real(wp), intent(in), optional :: dt
      ! assumed-shape-ok: pure passthrough to the _on shim.
      real(wp), intent(in), optional :: u_src(:, :, :), v_src(:, :, :), h_src(:, :, :)

      if (present(u_src)) then
         call ocean_horizontal_viscosity_compute_tendencies_on(grid, metrics, this, ms, &
                                                               u_src, v_src, h_src, lateral_mix=lateral_mix, dt=dt)
      else
         call ocean_horizontal_viscosity_compute_tendencies_on(grid, metrics, this, ms, &
                                                               ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, &
                                                               lateral_mix=lateral_mix, dt=dt)
      end if
   end subroutine ocean_horizontal_viscosity_compute_tendencies

   subroutine ocean_horizontal_viscosity_compute_tendencies_on(grid, metrics, this, ms, u, v, h, lateral_mix, dt)
      !! Fill `du_visc` and `dv_visc` with `nu_h * Laplacian` of the
      !! face velocities, per layer.  Closed-wall faces (i=1, i=nx+1
      !! for u; j=1, j=ny+1 for v) get zero tendency.  Interior y-
      !! boundary rows on u (j=1, j=ny) and interior x-boundary
      !! columns on v (i=1, i=nx) also get zero — equivalent to a
      !! free-slip wall condition on the tangential velocity.
      !!
      !! Loop order: (k, j, i) — `j` outermost-but-one for NVHPC GPU
      !! coalescing on the innermost-array-dimension `i`.
      !!
      !! Curvilinear (design §2): the face-velocity Laplacian is the
      !! finite-volume divergence of the velocity gradient over the
      !! face control volume — x-flux differenced across T-points
      !! (`dy_dxT` ratio), y-flux across Bu corners (`dx_dyBu`),
      !! normalised by `iareaCu`; the v-face mirror uses `dx_dyT` /
      !! `dy_dxBu` / `iareaCv`.  On uniform SQUARE metrics every ratio
      !! is 1 and `iareaCu = 1/(dx·dy)`, so the form collapses to the
      !! decoupled `Δ²u·(1/dx²) + Δ²u·(1/dy²)` to round-off (the only
      !! departure is FP reassociation of the 3-point grouping).
      !!
      !! Closure dispatch: when `lateral_mix` is present and its
      !! `closure` field is non-default (i.e. `/= LMIX_NONE`), the
      !! kernel reads per-face viscosity from
      !! `lateral_mix%ah_face_x` / `ah_face_y` instead of the scalar
      !! `nu_h`.  Caller must invoke
      !! `ocean_lateral_mix_compute_leith` (or equivalent) first to
      !! populate those fields.  Omitting the argument or leaving
      !! `closure = LMIX_NONE` preserves the scalar-`nu_h` behaviour
      !! bit-identically — existing call sites stay valid.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_horizontal_viscosity_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      ! assumed-shape-ok: forwarded straight to explicit-shape impl dummies,
      ! never indexed here (outer-shim source-selection layer, SPEC S3).
      real(wp), intent(in) :: u(:, :, :), v(:, :, :), h(:, :, :)
         !! Velocity/thickness source arrays — the prognostic components on
         !! the historical path, the `u_av` time-mean family under
         !! `split_scheme = "pred_corr"` (MOM6 evaluates horizontal_viscosity
         !! on `u_av`/`h_av`).
      type(ocean_lateral_mix_t), intent(in), optional :: lateral_mix
      real(wp), intent(in), optional :: dt
         !! Outer (or RK2-stage) time step.  Required when
         !! `this%stress_tensor` is true (drives the per-cell CFL
         !! viscosity limiter) and whenever a biharmonic add-on is
         !! active (drives its own per-face CFL clamp) — the velocity-
         !! Laplacian dispatch itself does not consume it, but the
         !! biharmonic block reached from every dispatch arm does.
         !! Omitting it defaults the biharmonic clamp to `dt_local = 1`
         !! (a huge, effectively-inactive bound); the caller is
         !! responsible for CFL in that case.

      integer :: nx, ny, nz
      real(wp) :: nu_h, nu_4, dt_local
      logical :: use_face_visc

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nu_h = this%nu_h
      nu_4 = this%nu_4

      ! dt for the per-cell CFL clamps (harmonic BOUND_KH on the
      ! stress-tensor path, biharmonic on the add-on below).  Defaults
      ! to 1 so the bound is huge when no dt was supplied — caller
      ! respects CFL.  Single hoisted site: both the stress branch and
      ! the biharmonic add-on below now read it.
      dt_local = 1.0_wp
      if (present(dt)) dt_local = dt

      ! Flow-aware HARMONIC face viscosity (`ah_face_*`) is read for the
      ! Laplacian closures (Leith / Smagorinsky) and for the live
      ! velocity-scale floor (`kh_vel_scale_live > 0`), which also
      ! populates `ah_face_*`.  The biharmonic closures
      ! (`LMIX_LEITH_BIHARM`, `smag_ah_active`) fill `nu4_face_*` instead
      ! and engage via the biharmonic add-on below; they must NOT pull the
      ! (unfilled) `ah_face_*` into the Laplacian — the scalar `nu_h`
      ! Laplacian stays in force underneath them.
      use_face_visc = .false.
      if (present(lateral_mix)) then
         if (lateral_mix%is_init .and. &
             (lateral_mix%closure == LMIX_LEITH .or. &
              lateral_mix%closure == LMIX_SMAGORINSKY .or. &
              lateral_mix%kh_vel_scale_live > 0.0_wp)) then
            use_face_visc = .true.
         end if
      end if

      ! ---- MOM6-faithful thickness-weighted stress-divergence path ----
      ! Replaces the velocity Laplacian with a single momentum-
      ! conserving stress assembly; per-cell CFL limiter + wet_*
      ! coast-masking baked in.  Composes with the biharmonic add-on
      ! below exactly as the Laplacian arms do (MOM6 `BIHARMONIC` "may
      ! be used with `LAPLACIAN`") — the harmonic part is momentum-
      ! conserving/coast-masked, the velocity-form biharmonic part is
      ! not (a documented fidelity divergence; a stress-level
      ! biharmonic is the follow-on).
      if (this%stress_tensor) then
         ! Step A: average the per-face harmonic A onto T-cells (ah_t)
         ! and Bu corners (ah_q), reading the flow-aware lateral-mix
         ! fields when active, else the scalar nu_h.  Then CFL-clamp
         ! each per its discrete stencil bound.
         if (use_face_visc) then
            call hvisc_avg_A_face( &
               lateral_mix%ah_face_x, lateral_mix%ah_face_y, &
               this%ah_t%data, this%ah_q%data, nx, ny, nz)
         else
            call hvisc_fill_A_scalar(this%ah_t%data, this%ah_q%data, nu_h, nx, ny, nz)
         end if

         ! Step A2 (anisotropic, Smith & McWilliams 2003): ADD the
         ! direction-tensor coefficients to the isotropic A before the
         ! CFL clamp (MOM6 adds anisotropy into Kh prior to BOUND_KH).
         ! Tension (T-cell): Kh += kh_aniso·(1−n1n2²);
         ! shear (Bu corner): Kh += kh_aniso·n1n2².  No-op when
         ! kh_aniso ≤ 0 ⇒ bit-identical isotropic path.
         if (this%kh_aniso > 0.0_wp) then
            call hvisc_add_aniso_coef( &
               this%ah_t%data, this%ah_q%data, &
               this%kh_aniso, this%aniso_n1n2, nx, ny, nz)
         end if

         call hvisc_clamp_A( &
            this%ah_t%data, this%ah_q%data, this%bound_coef, dt_local, &
            metrics%idxT, metrics%idyT, metrics%idxCu, metrics%idyCu, &
            metrics%idxCv, metrics%idyCv, metrics%iareaCu, metrics%iareaCv, &
            nx, ny, nz)

         ! Step B: stress assembly + thickness-weighted divergence.
         ! The anisotropic CROSS terms (kh_aniso·n1n2·(n1²−n2²)·strain)
         ! are folded into str_xx / str_xy inside the assembly when
         ! kh_aniso > 0; they vanish for the default grid-i direction
         ! (n1n2 = 0).
         call hvisc_compute_stress( &
            u, v, h, &
            this%str_xx%data, this%str_xy%data, &
            this%ah_t%data, this%ah_q%data, &
            this%du_visc%data, this%dv_visc%data, &
            merge(1.0_wp, 0.0_wp, this%no_slip), &
            this%kh_aniso, this%aniso_n1n2, this%aniso_n1n1_m_n2n2, &
            metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
            metrics%dx2h, metrics%dy2h, metrics%dx2q, metrics%dy2q, &
            metrics%dy_dxT, metrics%dx_dyT, metrics%dy_dxBu, metrics%dx_dyBu, &
            metrics%iareaCu, metrics%iareaCv, &
            metrics%wet_u, metrics%wet_v, metrics%wet_q, &
            nx, ny, nz)
         ! Outer shim: dereference the multilayer + lateral-mix + this
         ! components on host, then dispatch to a flat-impl with explicit-
         ! shape array dummies.  Without this layering NVHPC's `do
         ! concurrent` body sees `this%du_visc%data`, `lateral_mix%ah_face_x`,
         ! `u` as derived-type deep derefs and emits per-
         ! iteration descriptor-walk memcpys.
      else if (use_face_visc) then
         ! z-level closed faces: the mask actuals are ABSENT on the
         ! default path (the `(1,1,1)` placeholder must never reach an
         ! explicit-shape dummy), so the call is written twice rather
         ! than the argument once.  Cold dispatcher code.
         if (metrics%use_closed_faces) then
            call hvisc_compute_face_impl( &
               u, v, &
               lateral_mix%ah_face_x, lateral_mix%ah_face_y, &
               this%du_visc%data, this%dv_visc%data, &
               metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
               metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
               this%bound_kh, this%bound_coef, dt_local, &
               metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
               nx, ny, nz, metrics%open_u, metrics%open_v)
         else
            call hvisc_compute_face_impl( &
               u, v, &
               lateral_mix%ah_face_x, lateral_mix%ah_face_y, &
               this%du_visc%data, this%dv_visc%data, &
               metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
               metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
               this%bound_kh, this%bound_coef, dt_local, &
               metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
               nx, ny, nz)
         end if
      else
         if (metrics%use_closed_faces) then
            call hvisc_compute_scalar_impl( &
               u, v, &
               this%du_visc%data, this%dv_visc%data, &
               nu_h, &
               metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
               metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
               this%bound_kh, this%bound_coef, dt_local, &
               metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
               nx, ny, nz, metrics%open_u, metrics%open_v)
         else
            call hvisc_compute_scalar_impl( &
               u, v, &
               this%du_visc%data, this%dv_visc%data, &
               nu_h, &
               metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
               metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
               this%bound_kh, this%bound_coef, dt_local, &
               metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
               nx, ny, nz)
         end if
      end if

      ! Biharmonic add-on.  Independent of the harmonic dispatch above
      ! (stress-tensor, face, or scalar) — sums into the same tendency
      ! buffer whichever of the three wrote it.  When `nu_4 = 0` the
      ! kernel is a no-op (early return).  Scale-selective damping for
      ! stratified closed-basin runs: damps grid-scale modes
      ! O(ν_4·k⁴) without spilling into basin-scale modes the way a
      ! large Laplacian ν_h does.  Required for production Tasman-class
      ! runs where the dyn-core has a parametric instability that the
      ! standard Laplacian operator cannot reach.
      ! Flow-aware biharmonic ν₄ (`nu4_face_*`) is filled either by the
      ! strain-rate Smagorinsky_AH (`smag_ah_active`) or by the 2-D Leith
      ! biharmonic closure (`LMIX_LEITH_BIHARM`); both engage the same
      ! per-face biharmonic apply.  `dt_local` (hoisted above) drives the
      ! per-cell biharmonic CFL clamp.
      if (present(lateral_mix)) then
         if (lateral_mix%is_init .and. &
             (lateral_mix%smag_ah_active .or. &
              lateral_mix%closure == LMIX_LEITH_BIHARM)) then
            call hvisc_compute_biharmonic_face_impl( &
               u, v, &
               this%lap_u%data, this%lap_v%data, &
               this%du_visc%data, this%dv_visc%data, &
               lateral_mix%nu4_face_x, lateral_mix%nu4_face_y, &
               this%bound_coef, dt_local, &
               metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
               metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
               metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
               nx, ny, nz)
            return
         end if
      end if
      if (nu_4 > 0.0_wp) then
         call hvisc_compute_biharmonic_impl( &
            u, v, &
            this%lap_u%data, this%lap_v%data, &
            this%du_visc%data, this%dv_visc%data, &
            nu_4, this%bound_coef, dt_local, &
            metrics%idxCu, metrics%idyCu, metrics%idxCv, metrics%idyCv, &
            metrics%dy_dxT, metrics%dx_dyBu, metrics%iareaCu, &
            metrics%dx_dyT, metrics%dy_dxBu, metrics%iareaCv, &
            nx, ny, nz)
      end if
   end subroutine ocean_horizontal_viscosity_compute_tendencies_on

   pure subroutine hvisc_compute_face_impl(u_face, v_face, ah_face_x, ah_face_y, &
                                           du_visc, dv_visc, &
                                           dy_dxT, dx_dyBu, iareaCu, &
                                           dx_dyT, dy_dxBu, iareaCv, &
                                           bound_kh, bound_coef, dt, &
                                           idxCu, idyCu, idxCv, idyCv, &
                                           nx, ny, nz, open_u, open_v)
      !! Per-face metric Laplacian × spatially-varying viscosity.
      !! Explicit-shape dummies so NVHPC stdpar emits a device kernel
      !! without per-launch descriptor walks.  See
      !! `metric_lap_u`/`metric_lap_v` for the curvilinear FV form.
      !! `bound_kh` engages the per-face harmonic CFL ceiling
      !! (`hvisc_kh_cfl_bound`); `.false.` ⇒ bit-identical.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: ah_face_x(nx + 1, ny, nz), ah_face_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: dy_dxT(nx, ny), dx_dyBu(nx + 1, ny + 1), iareaCu(nx + 1, ny)
      real(wp), intent(in)    :: dx_dyT(nx, ny), dy_dxBu(nx + 1, ny + 1), iareaCv(nx, ny + 1)
      logical, intent(in)     :: bound_kh
      real(wp), intent(in)    :: bound_coef, dt
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in)    :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in), optional :: open_u(nx + 1, ny, nz)
         !! Per-layer 0/1 u-face open mask
         !! (`&vcoord_nml zfixed_closed_faces`).  ABSENT (the default
         !! path) => the interior loops below are textually the ones this
         !! routine has always run => bit-identical.
         !!
         !! PRESENT => the closed faces become FREE-SLIP walls, which is
         !! what a z-level partial step is (Adcroft, Hill & Marshall
         !! 1997).  Two things happen, and both are needed:
         !!
         !! 1. each neighbour difference is multiplied by the NEIGHBOUR
         !!    face's open flag, so a closed neighbour -- whose velocity
         !!    `mask_layer_velocities` has zeroed -- contributes nothing.
         !!    Without it the zero reads as a Dirichlet-0 boundary, i.e.
         !!    NO-SLIP: the wall would exert `nu*u/dx^2` of drag on the
         !!    live face beside it every step, which is the opposite of
         !!    the free-slip a partial step is supposed to be;
         !! 2. the whole tendency is multiplied by the face's OWN open
         !!    flag, so a closed face gets exactly zero viscous tendency
         !!    (it is a wall; `mask_layer_velocities` would zero it
         !!    anyway, but leaving a tendency there would make
         !!    `ke_diss` -- which MEKE reads -- account for work done on
         !!    water that is not there).
      real(wp), intent(in), optional :: open_v(nx, ny + 1, nz)
         !! v-face twin.  Present iff `open_u` is.
      integer :: i, j, k
      real(wp) :: lap_u, lap_v, nu_eff, idt

      idt = 0.0_wp
      if (dt > 0.0_wp) idt = 1.0_wp/dt

      ! u-face Laplacian interior + zero boundaries (curvilinear FV form)
      if (present(open_u)) then
         do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(lap_u, nu_eff)
            lap_u = iareaCu(i, j)*( &
                    (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k))*open_u(i + 1, j, k) - &
                     dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))*open_u(i - 1, j, k)) + &
                    (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k))*open_u(i, j + 1, k) - &
                     dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))*open_u(i, j - 1, k)))
            nu_eff = ah_face_x(i, j, k)
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
            end if
            du_visc(i, j, k) = nu_eff*lap_u*open_u(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(lap_u, nu_eff)
            lap_u = iareaCu(i, j)*( &
                    (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k)) - &
                     dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))) + &
                    (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k)) - &
                     dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))))
            nu_eff = ah_face_x(i, j, k)
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
            end if
            du_visc(i, j, k) = nu_eff*lap_u
         end do
      end if
      do concurrent(k=1:nz, j=1:ny)
         du_visc(1, j, k) = 0.0_wp
         du_visc(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         du_visc(i, 1, k) = 0.0_wp
         du_visc(i, ny, k) = 0.0_wp
      end do

      ! v-face Laplacian interior + zero boundaries (curvilinear FV form)
      if (present(open_v)) then
         do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(lap_v, nu_eff)
            lap_v = iareaCv(i, j)*( &
                    (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k))*open_v(i, j + 1, k) - &
                     dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))*open_v(i, j - 1, k)) + &
                    (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k))*open_v(i + 1, j, k) - &
                     dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))*open_v(i - 1, j, k)))
            nu_eff = ah_face_y(i, j, k)
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
            end if
            dv_visc(i, j, k) = nu_eff*lap_v*open_v(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(lap_v, nu_eff)
            lap_v = iareaCv(i, j)*( &
                    (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k)) - &
                     dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))) + &
                    (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k)) - &
                     dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))))
            nu_eff = ah_face_y(i, j, k)
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
            end if
            dv_visc(i, j, k) = nu_eff*lap_v
         end do
      end if
      do concurrent(k=1:nz, i=1:nx)
         dv_visc(i, 1, k) = 0.0_wp
         dv_visc(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         dv_visc(1, j, k) = 0.0_wp
         dv_visc(nx, j, k) = 0.0_wp
      end do
   end subroutine hvisc_compute_face_impl

   pure subroutine hvisc_compute_scalar_impl(u_face, v_face, du_visc, dv_visc, &
                                             nu_h, &
                                             dy_dxT, dx_dyBu, iareaCu, &
                                             dx_dyT, dy_dxBu, iareaCv, &
                                             bound_kh, bound_coef, dt, &
                                             idxCu, idyCu, idxCv, idyCv, &
                                             nx, ny, nz, open_u, open_v)
      !! Per-face metric Laplacian × scalar viscosity.  Used when no
      !! lateral-mix closure is active — falls back to constant
      !! `nu_h`.  When `nu_h = 0` the kernel still zeros all interior +
      !! boundary cells so the apply step sees a defined state.
      !! `bound_kh` engages the per-face harmonic CFL ceiling
      !! (`hvisc_kh_cfl_bound`); `.false.` ⇒ bit-identical.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(inout) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: nu_h
      real(wp), intent(in)    :: dy_dxT(nx, ny), dx_dyBu(nx + 1, ny + 1), iareaCu(nx + 1, ny)
      real(wp), intent(in)    :: dx_dyT(nx, ny), dy_dxBu(nx + 1, ny + 1), iareaCv(nx, ny + 1)
      logical, intent(in)     :: bound_kh
      real(wp), intent(in)    :: bound_coef, dt
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in)    :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in), optional :: open_u(nx + 1, ny, nz)
         !! Per-layer 0/1 u-face open mask
         !! (`&vcoord_nml zfixed_closed_faces`) — the FREE-SLIP closure
         !! of a z-level partial step.  See `hvisc_compute_face_impl`'s
         !! `open_u` docstring for the full argument; ABSENT (the default
         !! path) ⇒ the loops below are textually unchanged ⇒
         !! bit-identical.
      real(wp), intent(in), optional :: open_v(nx, ny + 1, nz)
         !! v-face twin.  Present iff `open_u` is.
      integer :: i, j, k
      real(wp) :: lap_u, lap_v, nu_eff, idt

      idt = 0.0_wp
      if (dt > 0.0_wp) idt = 1.0_wp/dt

      if (nu_h <= 0.0_wp) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
            du_visc(i, j, k) = 0.0_wp
         end do
         do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
            dv_visc(i, j, k) = 0.0_wp
         end do
         return
      end if

      if (present(open_u)) then
         do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(lap_u, nu_eff)
            lap_u = iareaCu(i, j)*( &
                    (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k))*open_u(i + 1, j, k) - &
                     dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))*open_u(i - 1, j, k)) + &
                    (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k))*open_u(i, j + 1, k) - &
                     dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))*open_u(i, j - 1, k)))
            nu_eff = nu_h
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
            end if
            du_visc(i, j, k) = nu_eff*lap_u*open_u(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(lap_u, nu_eff)
            lap_u = iareaCu(i, j)*( &
                    (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k)) - &
                     dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))) + &
                    (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k)) - &
                     dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))))
            nu_eff = nu_h
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
            end if
            du_visc(i, j, k) = nu_eff*lap_u
         end do
      end if
      do concurrent(k=1:nz, j=1:ny)
         du_visc(1, j, k) = 0.0_wp
         du_visc(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         du_visc(i, 1, k) = 0.0_wp
         du_visc(i, ny, k) = 0.0_wp
      end do

      if (present(open_v)) then
         do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(lap_v, nu_eff)
            lap_v = iareaCv(i, j)*( &
                    (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k))*open_v(i, j + 1, k) - &
                     dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))*open_v(i, j - 1, k)) + &
                    (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k))*open_v(i + 1, j, k) - &
                     dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))*open_v(i - 1, j, k)))
            nu_eff = nu_h
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
            end if
            dv_visc(i, j, k) = nu_eff*lap_v*open_v(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(lap_v, nu_eff)
            lap_v = iareaCv(i, j)*( &
                    (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k)) - &
                     dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))) + &
                    (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k)) - &
                     dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))))
            nu_eff = nu_h
            if (bound_kh) then
               nu_eff = min(nu_eff, hvisc_kh_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
            end if
            dv_visc(i, j, k) = nu_eff*lap_v
         end do
      end if
      do concurrent(k=1:nz, i=1:nx)
         dv_visc(i, 1, k) = 0.0_wp
         dv_visc(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         dv_visc(1, j, k) = 0.0_wp
         dv_visc(nx, j, k) = 0.0_wp
      end do
   end subroutine hvisc_compute_scalar_impl

   pure subroutine hvisc_compute_biharmonic_impl(u_face, v_face, lap_u, lap_v, &
                                                 du_visc, dv_visc, &
                                                 nu_4, bound_coef, dt, &
                                                 idxCu, idyCu, idxCv, idyCv, &
                                                 dy_dxT, dx_dyBu, iareaCu, &
                                                 dx_dyT, dy_dxBu, iareaCv, &
                                                 nx, ny, nz)
      !! Constant-coefficient biharmonic friction: applies
      !! `-ν₄ · ∇²(∇²u)` to the face velocities via two chained 5-point
      !! Laplacians.  Adds into the existing `du_visc` / `dv_visc`
      !! buffers (which Laplacian friction has already filled), so the
      !! caller can run with both `nu_h` and `nu_4` non-zero.
      !!
      !! Wall convention: the first Laplacian writes zero at every
      !! C-grid wall face (i=1, i=nx+1 on u; j=1, j=ny+1 on v) and at
      !! the interior tangential-boundary rows.  The second Laplacian
      !! reads `lap_u / lap_v` and sees those zeros — equivalent to a
      !! "no-stress on Δu" wall condition, which is the closure
      !! MOM6's `BIHARMONIC` block uses for closed boundaries.
      !!
      !! Sign convention: the operator is `du/dt = -ν₄·∇⁴u` so the
      !! sinusoid `u = e^{ikx}` has growth rate `-ν₄·k⁴` — damping for
      !! positive ν₄, scaling as `k⁴` rather than `k²` (Laplacian).
      !!
      !! Stability: forward-Euler biharmonic CFL is
      !! `ν₄ · dt · ((π/dx)² + (π/dy)²)² <= 2` (established project
      !! constant).  The constant scalar `nu_4` is clamped per face to
      !! this CFL ceiling (`hvisc_nu4_cfl_bound`, scaled by
      !! `bound_coef`) in Pass 2, so an over-large namelist `nu_4`
      !! cannot violate the local bound; below the ceiling the multiply
      !! is by `nu_4` bit-for-bit.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(inout) :: lap_u(nx + 1, ny, nz), lap_v(nx, ny + 1, nz)
      real(wp), intent(inout) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: nu_4
      real(wp), intent(in)    :: bound_coef, dt
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in)    :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in)    :: dy_dxT(nx, ny), dx_dyBu(nx + 1, ny + 1), iareaCu(nx + 1, ny)
      real(wp), intent(in)    :: dx_dyT(nx, ny), dy_dxBu(nx + 1, ny + 1), iareaCv(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: l_u, l_v, idt, nu4_u, nu4_v
      idt = 1.0_wp/dt

      ! ---- Pass 1: lap_u, lap_v at interior u-faces and v-faces ----
      ! Wall BC on the intermediate Laplacian: MIRROR (Neumann),
      ! not zero.  Hard-zero on lap_u/lap_v creates a step
      ! discontinuity that Pass 2 reads as a high-k mode at the
      ! wall-adjacent interior face → injects spurious work into the
      ! boundary band.  Mirror BC (copy the adjacent interior value)
      ! preserves the "smooth interior field ⇒ zero biharmonic"
      ! property — a smooth lap_u extended by mirror gives ∇² = 0 at
      ! the boundary too.  This is the equivalent of MOM6's "stress
      ! free at the boundary" closure on biharmonic friction.
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(l_u)
         l_u = iareaCu(i, j)*( &
               (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k)) - &
                dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))) + &
               (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k)) - &
                dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))))
         lap_u(i, j, k) = l_u
      end do
      do concurrent(k=1:nz, j=1:ny)
         lap_u(1, j, k) = lap_u(2, j, k)
         lap_u(nx + 1, j, k) = lap_u(nx, j, k)
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         lap_u(i, 1, k) = lap_u(i, 2, k)
         lap_u(i, ny, k) = lap_u(i, ny - 1, k)
      end do

      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(l_v)
         l_v = iareaCv(i, j)*( &
               (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k)) - &
                dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))) + &
               (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k)) - &
                dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))))
         lap_v(i, j, k) = l_v
      end do
      do concurrent(k=1:nz, i=1:nx)
         lap_v(i, 1, k) = lap_v(i, 2, k)
         lap_v(i, ny + 1, k) = lap_v(i, ny, k)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         lap_v(1, j, k) = lap_v(2, j, k)
         lap_v(nx, j, k) = lap_v(nx - 1, j, k)
      end do

      ! ---- Pass 2: -ν₄ · ∇²(lap_*) added into du_visc / dv_visc ----
      ! ν₄ is clamped per face to the local explicit-biharmonic CFL
      ! ceiling (bound_coef · 2/(dt·k⁴)).  Below the ceiling
      ! `min(nu_4, bound) == nu_4` ⇒ bit-identical to the unclamped add.
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(l_u, nu4_u)
         l_u = iareaCu(i, j)*( &
               (dy_dxT(i, j)*(lap_u(i + 1, j, k) - lap_u(i, j, k)) - &
                dy_dxT(i - 1, j)*(lap_u(i, j, k) - lap_u(i - 1, j, k))) + &
               (dx_dyBu(i, j + 1)*(lap_u(i, j + 1, k) - lap_u(i, j, k)) - &
                dx_dyBu(i, j)*(lap_u(i, j, k) - lap_u(i, j - 1, k))))
         nu4_u = min(nu_4, hvisc_nu4_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
         du_visc(i, j, k) = du_visc(i, j, k) - nu4_u*l_u
      end do

      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(l_v, nu4_v)
         l_v = iareaCv(i, j)*( &
               (dx_dyT(i, j)*(lap_v(i, j + 1, k) - lap_v(i, j, k)) - &
                dx_dyT(i, j - 1)*(lap_v(i, j, k) - lap_v(i, j - 1, k))) + &
               (dy_dxBu(i + 1, j)*(lap_v(i + 1, j, k) - lap_v(i, j, k)) - &
                dy_dxBu(i, j)*(lap_v(i, j, k) - lap_v(i - 1, j, k))))
         nu4_v = min(nu_4, hvisc_nu4_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
         dv_visc(i, j, k) = dv_visc(i, j, k) - nu4_v*l_v
      end do
   end subroutine hvisc_compute_biharmonic_impl

   pure subroutine hvisc_compute_biharmonic_face_impl(u_face, v_face, lap_u, lap_v, &
                                                      du_visc, dv_visc, &
                                                      nu4_face_x, nu4_face_y, &
                                                      bound_coef, dt, &
                                                      idxCu, idyCu, idxCv, idyCv, &
                                                      dy_dxT, dx_dyBu, iareaCu, &
                                                      dx_dyT, dy_dxBu, iareaCv, &
                                                      nx, ny, nz)
      !! Flow-aware biharmonic friction (MOM6 SMAGORINSKY_AH analogue).
      !! Identical to `hvisc_compute_biharmonic_impl` except Pass 2
      !! multiplies the second Laplacian by the per-face viscosity
      !! `nu4_face_x/y` instead of the scalar `nu_4`.  The face fields
      !! are filled upstream by `ocean_lateral_mix_compute_smag_ah`,
      !! which sets them to `C_b · L⁴ · |D|` clamped to
      !! `[nu4_bg, nu4_max]`.  Pass 2 additionally clamps each face
      !! coefficient to the per-face explicit-biharmonic CFL ceiling
      !! (`hvisc_nu4_cfl_bound`, scaled by `bound_coef`) on top of the
      !! static `nu4_max` floor/ceiling applied upstream — so a strain
      !! spike on a fine cell can never violate the local CFL bound.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(inout) :: lap_u(nx + 1, ny, nz), lap_v(nx, ny + 1, nz)
      real(wp), intent(inout) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: nu4_face_x(nx + 1, ny, nz), nu4_face_y(nx, ny + 1, nz)
      real(wp), intent(in)    :: bound_coef, dt
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in)    :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in)    :: dy_dxT(nx, ny), dx_dyBu(nx + 1, ny + 1), iareaCu(nx + 1, ny)
      real(wp), intent(in)    :: dx_dyT(nx, ny), dy_dxBu(nx + 1, ny + 1), iareaCv(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: l_u, l_v, idt, nu4_u, nu4_v
      idt = 1.0_wp/dt

      ! ---- Pass 1: lap_u, lap_v at interior u/v faces (metric FV) ----
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(l_u)
         l_u = iareaCu(i, j)*( &
               (dy_dxT(i, j)*(u_face(i + 1, j, k) - u_face(i, j, k)) - &
                dy_dxT(i - 1, j)*(u_face(i, j, k) - u_face(i - 1, j, k))) + &
               (dx_dyBu(i, j + 1)*(u_face(i, j + 1, k) - u_face(i, j, k)) - &
                dx_dyBu(i, j)*(u_face(i, j, k) - u_face(i, j - 1, k))))
         lap_u(i, j, k) = l_u
      end do
      do concurrent(k=1:nz, j=1:ny)
         lap_u(1, j, k) = lap_u(2, j, k)
         lap_u(nx + 1, j, k) = lap_u(nx, j, k)
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         lap_u(i, 1, k) = lap_u(i, 2, k)
         lap_u(i, ny, k) = lap_u(i, ny - 1, k)
      end do

      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(l_v)
         l_v = iareaCv(i, j)*( &
               (dx_dyT(i, j)*(v_face(i, j + 1, k) - v_face(i, j, k)) - &
                dx_dyT(i, j - 1)*(v_face(i, j, k) - v_face(i, j - 1, k))) + &
               (dy_dxBu(i + 1, j)*(v_face(i + 1, j, k) - v_face(i, j, k)) - &
                dy_dxBu(i, j)*(v_face(i, j, k) - v_face(i - 1, j, k))))
         lap_v(i, j, k) = l_v
      end do
      do concurrent(k=1:nz, i=1:nx)
         lap_v(i, 1, k) = lap_v(i, 2, k)
         lap_v(i, ny + 1, k) = lap_v(i, ny, k)
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         lap_v(1, j, k) = lap_v(2, j, k)
         lap_v(nx, j, k) = lap_v(nx - 1, j, k)
      end do

      ! ---- Pass 2: -nu4_face · ∇²(lap_*) added into du_visc / dv_visc ----
      ! Each per-face nu4 is clamped to the local CFL ceiling; below it
      ! `min(...) == nu4_face` ⇒ bit-identical to the unclamped add.
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx) local(l_u, nu4_u)
         l_u = iareaCu(i, j)*( &
               (dy_dxT(i, j)*(lap_u(i + 1, j, k) - lap_u(i, j, k)) - &
                dy_dxT(i - 1, j)*(lap_u(i, j, k) - lap_u(i - 1, j, k))) + &
               (dx_dyBu(i, j + 1)*(lap_u(i, j + 1, k) - lap_u(i, j, k)) - &
                dx_dyBu(i, j)*(lap_u(i, j, k) - lap_u(i, j - 1, k))))
         nu4_u = min(nu4_face_x(i, j, k), &
                     hvisc_nu4_cfl_bound(idxCu(i, j), idyCu(i, j), bound_coef, idt))
         du_visc(i, j, k) = du_visc(i, j, k) - nu4_u*l_u
      end do

      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1) local(l_v, nu4_v)
         l_v = iareaCv(i, j)*( &
               (dx_dyT(i, j)*(lap_v(i, j + 1, k) - lap_v(i, j, k)) - &
                dx_dyT(i, j - 1)*(lap_v(i, j, k) - lap_v(i, j - 1, k))) + &
               (dy_dxBu(i + 1, j)*(lap_v(i + 1, j, k) - lap_v(i, j, k)) - &
                dy_dxBu(i, j)*(lap_v(i, j, k) - lap_v(i - 1, j, k))))
         nu4_v = min(nu4_face_y(i, j, k), &
                     hvisc_nu4_cfl_bound(idxCv(i, j), idyCv(i, j), bound_coef, idt))
         dv_visc(i, j, k) = dv_visc(i, j, k) - nu4_v*l_v
      end do
   end subroutine hvisc_compute_biharmonic_face_impl

   pure subroutine hvisc_fill_A_scalar(ah_t, ah_q, nu_h, nx, ny, nz)
      !! Fill the T-cell and corner harmonic-viscosity fields with the
      !! scalar `nu_h` (no flow-aware closure active).  Constant ⇒ the
      !! T/corner averaging is exact, so the all-wet uniform-h reduction
      !! to the velocity Laplacian holds bit-for-bit.
      integer, intent(in)  :: nx, ny, nz
      real(wp), intent(out) :: ah_t(nx, ny, nz), ah_q(nx + 1, ny + 1, nz)
      real(wp), intent(in)  :: nu_h
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ah_t(i, j, k) = nu_h
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx + 1)
         ah_q(i, j, k) = nu_h
      end do
   end subroutine hvisc_fill_A_scalar

   pure subroutine hvisc_add_aniso_coef(ah_t, ah_q, kh_aniso, n1n2, nx, ny, nz)
      !! Add the Smith & McWilliams (2003) anisotropic direction-tensor
      !! coefficients onto the co-located isotropic viscosities.  The
      !! tension (T-cell) coefficient gains `kh_aniso·(1−n1n2²)` and the
      !! shear (Bu-corner) coefficient gains `kh_aniso·n1n2²`.  For the
      !! default grid-i direction `n1n2 = 0` ⇒ T gains `kh_aniso`, the
      !! corner gains nothing — stronger damping of along-i tension.
      !! The corner outer ring stays untouched (the divergence stencil
      !! never reads it; `hvisc_avg_A_face` already zeroed it).
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(inout) :: ah_t(nx, ny, nz), ah_q(nx + 1, ny + 1, nz)
      real(wp), intent(in)    :: kh_aniso, n1n2
      integer :: i, j, k
      real(wp) :: add_t, add_q
      add_t = kh_aniso*(1.0_wp - n1n2*n1n2)
      add_q = kh_aniso*(n1n2*n1n2)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ah_t(i, j, k) = ah_t(i, j, k) + add_t
      end do
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         ah_q(i, j, k) = ah_q(i, j, k) + add_q
      end do
   end subroutine hvisc_add_aniso_coef

   pure subroutine hvisc_avg_A_face(ah_face_x, ah_face_y, ah_t, ah_q, nx, ny, nz)
      !! Average the per-face harmonic viscosity (`ah_face_x` at
      !! u-faces, `ah_face_y` at v-faces) onto the T-cell centres
      !! (`ah_t`) and the Bu corners (`ah_q`).  The stress form needs A
      !! co-located with the tension (T-cell) and shear (corner)
      !! strains; the lateral-mix closure produces A at faces, so this
      !! is a 4-point face→cell / face→corner reduction.  On a uniform A
      !! field every average returns A, preserving the Laplacian
      !! reduction.
      integer, intent(in)  :: nx, ny, nz
      real(wp), intent(in)  :: ah_face_x(nx + 1, ny, nz), ah_face_y(nx, ny + 1, nz)
      real(wp), intent(out) :: ah_t(nx, ny, nz), ah_q(nx + 1, ny + 1, nz)
      integer :: i, j, k
      ! T-cell (i,j): mean of its west/east u-faces (i, i+1) and
      ! south/north v-faces (j, j+1).
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ah_t(i, j, k) = 0.25_wp*((ah_face_x(i, j, k) + ah_face_x(i + 1, j, k)) + &
                                  (ah_face_y(i, j, k) + ah_face_y(i, j + 1, k)))
      end do
      ! Bu corner (i,j) = SW corner of cell (i,j): mean of the two
      ! u-faces sharing it (i, j) & (i, j-1) and the two v-faces
      ! (i, j) & (i-1, j).  Interior corners only; the divergence
      ! stencil never reads the outermost corner ring.
      do concurrent(k=1:nz, j=2:ny, i=2:nx)
         ah_q(i, j, k) = 0.25_wp*((ah_face_x(i, j, k) + ah_face_x(i, j - 1, k)) + &
                                  (ah_face_y(i, j, k) + ah_face_y(i - 1, j, k)))
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         ah_q(1, j, k) = 0.0_wp
         ah_q(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         ah_q(i, 1, k) = 0.0_wp
         ah_q(i, ny + 1, k) = 0.0_wp
      end do
   end subroutine hvisc_avg_A_face

   pure subroutine hvisc_clamp_A(ah_t, ah_q, bound_coef, dt, &
                                 idxT, idyT, idxCu, idyCu, idxCv, idyCv, &
                                 iareaCu, iareaCv, nx, ny, nz)
      !! Per-cell CFL viscosity limiter (MOM6 `BOUND_KH`).  Clamps the
      !! T-cell viscosity to `Kh_Max_xx` and the corner viscosity to
      !! `Kh_Max_xy`, each derived from the actual discrete stress
      !! stencil metrics + `dt` so the explicit forward-Euler viscous
      !! update can never overshoot.  Replaces the global `ah_max` cap.
      !! On a uniform square grid `Kh_Max = bound_coef·0.25/(dt·(1/dx²+
      !! 1/dy²))`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: ah_t(nx, ny, nz), ah_q(nx + 1, ny + 1, nz)
      real(wp), intent(in) :: bound_coef, dt
      real(wp), intent(in) :: idxT(nx, ny), idyT(nx, ny)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in) :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in) :: iareaCu(nx + 1, ny), iareaCv(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: idt, kh_max, denom, tx, ty, dx2, dy2
      idt = 1.0_wp/dt
      ! T-cell bound: tension stress acts on the two u-faces (i, i+1)
      ! and two v-faces (j, j+1).  dx2/dy2 per cell from the metric
      ! inverses (dx2 = 1/idxT²).  Full T-cell range (incl. the edge
      ! rows whose stress feeds the interior divergence) — Cartesian /
      ! curvilinear metrics are valid in the ghosts, so the bound is
      ! defined everywhere; guard denom>0 for any zeroed land metric.
      do concurrent(k=1:nz, j=1:ny, i=1:nx) &
         local(kh_max, denom, tx, ty, dx2, dy2)
         dx2 = 1.0_wp/(idxT(i, j)*idxT(i, j))
         dy2 = 1.0_wp/(idyT(i, j)*idyT(i, j))
         tx = dy2*(idyT(i, j)/idxT(i, j))*(idyCu(i + 1, j) + idyCu(i, j))* &
              max(idyCu(i + 1, j)*iareaCu(i + 1, j), idyCu(i, j)*iareaCu(i, j))
         ty = dx2*(idxT(i, j)/idyT(i, j))*(idxCv(i, j + 1) + idxCv(i, j))* &
              max(idxCv(i, j + 1)*iareaCv(i, j + 1), idxCv(i, j)*iareaCv(i, j))
         denom = max(tx, ty)
         if (denom > 0.0_wp) then
            kh_max = bound_coef*0.25_wp*idt/denom
            if (ah_t(i, j, k) > kh_max) ah_t(i, j, k) = kh_max
         end if
      end do
      ! Corner bound: shear stress at Bu(i,j) acts on the u-faces
      ! (i, j-1)/(i, j) and v-faces (i-1, j)/(i, j).  Reuse the
      ! T-stencil metric magnitudes at the SW T-cell — the corner
      ! bound differs from the T-cell bound only by which faces it
      ! sums, and on a square grid both collapse to the same value.
      do concurrent(k=1:nz, j=2:ny, i=2:nx) &
         local(kh_max, denom, tx, ty, dx2, dy2)
         dx2 = 1.0_wp/(idxCu(i, j)*idxCu(i, j))
         dy2 = 1.0_wp/(idyCv(i, j)*idyCv(i, j))
         tx = dx2*(idxCu(i, j) + idxCu(i, j - 1))* &
              max(idxCu(i, j)*iareaCu(i, j), idxCu(i, j - 1)*iareaCu(i, j - 1))
         ty = dy2*(idyCv(i, j) + idyCv(i - 1, j))* &
              max(idyCv(i, j)*iareaCv(i, j), idyCv(i - 1, j)*iareaCv(i - 1, j))
         denom = max(tx, ty)
         if (denom > 0.0_wp) then
            kh_max = bound_coef*0.25_wp*idt/denom
            if (ah_q(i, j, k) > kh_max) ah_q(i, j, k) = kh_max
         end if
      end do
   end subroutine hvisc_clamp_A

   pure function hvisc_kh_cfl_bound(idx, idy, bound_coef, idt) result(kh_max_cfl)
      !! Per-face harmonic-viscosity ceiling for the velocity-Laplacian
      !! paths (MOM6 `Kh_Max_xx` analogue, uniform-grid reduction):
      !!     `ν_max = bound_coef · 0.125 / (dt · (idx² + idy²))`
      !! — one quarter of the forward-Euler stability limit
      !! `ν·dt·4·(idx²+idy²) ≤ 2`, the same margin MOM6's harmonic
      !! bound uses ("avoid overshoots when bound_coef < 1").  Keeps
      !! the FROZEN depth-mean viscous forcing on the barotropic mode
      !! (`F_bt`) out of the phase-reversed anti-damping regime for
      !! grid-scale gravity modes (see the `bound_kh` docstring).
      !! Returns huge (no clamp) for a fully-masked face.
      !! `!$acc routine seq` — called from the Laplacian `do concurrent`.
      real(wp), intent(in) :: idx, idy
         !! Metric inverses `1/dx`, `1/dy` at the face.
      real(wp), intent(in) :: bound_coef
         !! CFL safety coefficient (`this%bound_coef`).
      real(wp), intent(in) :: idt
         !! Reciprocal time step `1/dt`.
      real(wp) :: kh_max_cfl
      real(wp) :: k2
      !$acc routine seq
      k2 = idx*idx + idy*idy
      if (k2 > 0.0_wp) then
         kh_max_cfl = bound_coef*0.125_wp*idt/k2
      else
         kh_max_cfl = huge(1.0_wp)
      end if
   end function hvisc_kh_cfl_bound

   pure function hvisc_nu4_cfl_bound(idx, idy, bound_coef, idt) result(nu4_max_cfl)
      !! Per-face explicit-biharmonic CFL ceiling on the biharmonic
      !! viscosity `ν₄` (m⁴/s).  Forward-Euler stability for `−ν₄·∇⁴u`
      !! on a local cell of spacing `(dx, dy)` requires (established
      !! project constant)
      !!     `ν₄ · dt · ((π/dx)² + (π/dy)²)² ≤ 2`,
      !! so the per-face bound is
      !!     `ν₄_max = bound_coef · 2 / (dt · ((π·idx)² + (π·idy)²)²)`
      !! where `idx = 1/dx`, `idy = 1/dy` are the metric inverses at the
      !! face (`idxCu`/`idyCu` at u-faces, `idxCv`/`idyCv` at v-faces).
      !! `bound_coef` (MOM6 `HORVISC_BOUND_COEF`, default 0.8) is the
      !! CFL safety margin shared with the harmonic `hvisc_clamp_A`.
      !! Returns a huge value (no clamp) when the metric inverses are
      !! both zero (a fully-masked land face) so the caller's `min`
      !! leaves the coefficient untouched there.
      !! `!$acc routine seq` — called from the biharmonic `do concurrent`.
      real(wp), intent(in) :: idx, idy
         !! Metric inverses `1/dx`, `1/dy` at the face.
      real(wp), intent(in) :: bound_coef
         !! CFL safety coefficient (`this%bound_coef`).
      real(wp), intent(in) :: idt
         !! Reciprocal time step `1/dt`.
      real(wp) :: nu4_max_cfl
      real(wp) :: k2
      !$acc routine seq
      k2 = (PI*idx)*(PI*idx) + (PI*idy)*(PI*idy)
      if (k2 > 0.0_wp) then
         nu4_max_cfl = bound_coef*2.0_wp*idt/(k2*k2)
      else
         nu4_max_cfl = huge(1.0_wp)
      end if
   end function hvisc_nu4_cfl_bound

   pure function raw_sh_xx(u_face, v_face, idxCu, idyCu, idxCv, dy_dxT, dx_dyT, &
                           i, j, k, nx, ny, nz) result(sh_xx)
      !! Raw tension strain `sh_xx = du/dx − dv/dy` at T-cell (i,j),
      !! mirroring Phase-1's gradient form (unmasked — used only by the
      !! anisotropic cross term where the all-wet reduction is exact).
      !! `!$acc routine seq` so the stress-assembly `do concurrent` can
      !! call it on-device.
      !$acc routine seq
      integer, intent(in) :: i, j, k, nx, ny, nz
      real(wp), intent(in) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCu(nx + 1, ny), idxCv(nx, ny + 1)
      real(wp), intent(in) :: dy_dxT(nx, ny), dx_dyT(nx, ny)
      real(wp) :: sh_xx
      sh_xx = dy_dxT(i, j)*(idyCu(i + 1, j)*u_face(i + 1, j, k) - idyCu(i, j)*u_face(i, j, k)) &
              - dx_dyT(i, j)*(idxCv(i, j + 1)*v_face(i, j + 1, k) - idxCv(i, j)*v_face(i, j, k))
   end function raw_sh_xx

   pure function raw_sh_xy(u_face, v_face, idxCu, idyCv, dy_dxBu, dx_dyBu, &
                           i, j, k, nx, ny, nz) result(sh_xy)
      !! Raw shear strain `sh_xy = dv/dx + du/dy` at Bu corner (i,j),
      !! mirroring Phase-2's gradient form.  `!$acc routine seq` for the
      !! on-device cross-term loop.
      !$acc routine seq
      integer, intent(in) :: i, j, k, nx, ny, nz
      real(wp), intent(in) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: dy_dxBu(nx + 1, ny + 1), dx_dyBu(nx + 1, ny + 1)
      real(wp) :: sh_xy
      sh_xy = dy_dxBu(i, j)*(idyCv(i, j)*v_face(i, j, k) - idyCv(i - 1, j)*v_face(i - 1, j, k)) &
              + dx_dyBu(i, j)*(idxCu(i, j)*u_face(i, j, k) - idxCu(i, j - 1)*u_face(i, j - 1, k))
   end function raw_sh_xy

   pure subroutine hvisc_compute_stress(u_face, v_face, h_layer, &
                                        str_xx, str_xy, ah_t, ah_q, &
                                        du_visc, dv_visc, ns, &
                                        kh_aniso, n1n2, n1n1_m_n2n2, &
                                        idxCu, idyCu, idxCv, idyCv, &
                                        dx2h, dy2h, dx2q, dy2q, &
                                        dy_dxT, dx_dyT, dy_dxBu, dx_dyBu, &
                                        iareaCu, iareaCv, &
                                        wet_u, wet_v, wet_q, nx, ny, nz)
      !! MOM6 thickness-weighted stress-divergence operator.  Three
      !! phases: (1) tension `str_xx` at T-cells, (2) shear `str_xy` at
      !! Bu corners, (3) the divergence `(1/(h_u+h_neglect))·∂str`.
      !! `wet_u/wet_v/wet_q` mask the stress so no momentum is diffused
      !! across a coastline; `h_neglect = H_VANISHED` floors the
      !! velocity-point thickness.  See module header for the form +
      !! the all-wet uniform-h Laplacian reduction.
      !!
      !! Anisotropic cross terms (Smith & McWilliams 2003, §2): when
      !! `kh_aniso > 0` the tension stress at T-cells gains
      !! `+kh_aniso·n1n2·(n1²−n2²)·sh_xy_at_T·h_T` and the shear stress
      !! at corners gains `+kh_aniso·n1n2·(n1²−n2²)·sh_xx_at_q·h_q`,
      !! where the off-component strain is 4-point averaged onto the
      !! stress location.  (Sign is `+` here because Roundabout's stress
      !! convention `str = +A·strain·h` reduces the divergence to
      !! `+A·∇²u`, opposite to MOM6's `str = −Kh·sh`.)  For the default
      !! grid-i direction `n1n2 = 0` ⇒ the cross terms vanish and the
      !! assembly is bit-identical to the isotropic path.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)    :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: str_xx(nx, ny, nz), str_xy(nx + 1, ny + 1, nz)
      real(wp), intent(in)    :: ah_t(nx, ny, nz), ah_q(nx + 1, ny + 1, nz)
      real(wp), intent(inout) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: ns
      real(wp), intent(in)    :: kh_aniso, n1n2, n1n1_m_n2n2
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in)    :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(in)    :: dx2h(nx, ny), dy2h(nx, ny)
      real(wp), intent(in)    :: dx2q(nx + 1, ny + 1), dy2q(nx + 1, ny + 1)
      real(wp), intent(in)    :: dy_dxT(nx, ny), dx_dyT(nx, ny)
      real(wp), intent(in)    :: dy_dxBu(nx + 1, ny + 1), dx_dyBu(nx + 1, ny + 1)
      real(wp), intent(in)    :: iareaCu(nx + 1, ny), iareaCv(nx, ny + 1)
      real(wp), intent(in)    :: wet_u(nx + 1, ny), wet_v(nx, ny + 1)
      real(wp), intent(in)    :: wet_q(nx + 1, ny + 1)
      integer :: i, j, k
      real(wp) :: dudx, dvdy, dvdx, dudy, h_q, slip
      real(wp) :: aniso_cross, shxy_at_t, shxx_at_q
      real(wp), parameter :: H_NEGLECT = H_VANISHED

      aniso_cross = kh_aniso*n1n2*n1n1_m_n2n2

      ! ---- Phase 1: tension stress at T-cells (i,j) ----
      ! sh_xx = du/dx − dv/dy.  Masking the velocity gradients by the
      ! face wet masks keeps a land neighbour from contributing strain
      ! (MOM6 masks via reduction_xx + pre-masked metrics; the wet_u/
      ! wet_v product is the equivalent here).  All-wet ⇒ ×1.
      do concurrent(k=1:nz, j=1:ny, i=1:nx) local(dudx, dvdy)
         dudx = dy_dxT(i, j)*(idyCu(i + 1, j)*wet_u(i + 1, j)*u_face(i + 1, j, k) - &
                              idyCu(i, j)*wet_u(i, j)*u_face(i, j, k))
         dvdy = dx_dyT(i, j)*(idxCv(i, j + 1)*wet_v(i, j + 1)*v_face(i, j + 1, k) - &
                              idxCv(i, j)*wet_v(i, j)*v_face(i, j, k))
         str_xx(i, j, k) = ah_t(i, j, k)*(dudx - dvdy)*h_layer(i, j, k)
      end do

      ! ---- Phase 2: shear stress at Bu corners (i,j) = SW of cell (i,j) ----
      ! sh_xy = dv/dx + du/dy.  Corner thickness h_q = 4-pt mean of the
      ! surrounding cells (no halo issues — only interior corners are
      ! read by the divergence).  Slip factor: free-slip ×wet_q,
      ! no-slip ×(2−wet_q); all-wet ⇒ ×1.
      do concurrent(k=1:nz, j=2:ny, i=2:nx) local(dvdx, dudy, h_q, slip)
         dvdx = dy_dxBu(i, j)*(idyCv(i, j)*v_face(i, j, k) - idyCv(i - 1, j)*v_face(i - 1, j, k))
         dudy = dx_dyBu(i, j)*(idxCu(i, j)*u_face(i, j, k) - idxCu(i, j - 1)*u_face(i, j - 1, k))
         h_q = 0.25_wp*((h_layer(i - 1, j - 1, k) + h_layer(i, j, k)) + &
                        (h_layer(i - 1, j, k) + h_layer(i, j - 1, k)))
         slip = (1.0_wp - 2.0_wp*ns)*wet_q(i, j) + 2.0_wp*ns
         str_xy(i, j, k) = ah_q(i, j, k)*(dvdx + dudy)*h_q*slip
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         str_xy(1, j, k) = 0.0_wp
         str_xy(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         str_xy(i, 1, k) = 0.0_wp
         str_xy(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Phase 2.5: anisotropic cross terms (Smith & McWilliams 2003) ----
      ! Folded ONLY when kh_aniso·n1n2·(n1²−n2²) /= 0 — i.e. an
      ! off-grid-axis anisotropy direction.  The default grid-i
      ! direction (n1n2 = 0) skips both passes ⇒ bit-identical.  The
      ! cross terms couple the tension stress to the shear strain and
      ! vice versa, using the RAW velocity-gradient strains (not the
      ! already-formed stresses), recomputed inline:
      !   sh_xx(T)  = du/dx − dv/dy   (Phase-1 strain)
      !   sh_xy(Bu) = dv/dx + du/dy   (Phase-2 strain)
      ! str_xx(T) += aniso_cross · ⟨sh_xy⟩_corners→T · h_T
      ! str_xy(Bu)+= aniso_cross · ⟨sh_xx⟩_cells→Bu  · h_q · slip
      if (aniso_cross /= 0.0_wp) then
         ! str_xx gains the shear-strain contribution: average sh_xy
         ! from the 4 corners around T-cell (i,j): (i,j),(i+1,j),
         ! (i,j+1),(i+1,j+1).  Interior corners are valid for i=1:nx-1,
         ! j=1:ny-1; the boundary T-rows keep the isotropic value (the
         ! divergence weights them through the metric stencil).
         do concurrent(k=1:nz, j=2:ny - 1, i=2:nx - 1) local(shxy_at_t)
            shxy_at_t = 0.25_wp*( &
                        (raw_sh_xy(u_face, v_face, idxCu, idyCv, dy_dxBu, dx_dyBu, i, j, k, nx, ny, nz) + &
                         raw_sh_xy(u_face, v_face, idxCu, idyCv, dy_dxBu, dx_dyBu, i + 1, j + 1, k, nx, ny, nz)) + &
                        (raw_sh_xy(u_face, v_face, idxCu, idyCv, dy_dxBu, dx_dyBu, i + 1, j, k, nx, ny, nz) + &
                         raw_sh_xy(u_face, v_face, idxCu, idyCv, dy_dxBu, dx_dyBu, i, j + 1, k, nx, ny, nz)))
            str_xx(i, j, k) = str_xx(i, j, k) + aniso_cross*shxy_at_t*h_layer(i, j, k)
         end do
         ! str_xy gains the tension-strain contribution: average sh_xx
         ! from the 4 T-cells around corner (i,j): (i-1,j-1),(i,j-1),
         ! (i-1,j),(i,j).  Slip + h_q reuse Phase-2 forms.
         do concurrent(k=1:nz, j=2:ny, i=2:nx) local(shxx_at_q, h_q, slip)
            shxx_at_q = 0.25_wp*( &
                        (raw_sh_xx(u_face, v_face, idxCu, idyCu, idxCv, dy_dxT, dx_dyT, i - 1, j - 1, k, nx, ny, nz) + &
                         raw_sh_xx(u_face, v_face, idxCu, idyCu, idxCv, dy_dxT, dx_dyT, i, j, k, nx, ny, nz)) + &
                        (raw_sh_xx(u_face, v_face, idxCu, idyCu, idxCv, dy_dxT, dx_dyT, i, j - 1, k, nx, ny, nz) + &
                         raw_sh_xx(u_face, v_face, idxCu, idyCu, idxCv, dy_dxT, dx_dyT, i - 1, j, k, nx, ny, nz)))
            h_q = 0.25_wp*((h_layer(i - 1, j - 1, k) + h_layer(i, j, k)) + &
                           (h_layer(i - 1, j, k) + h_layer(i, j - 1, k)))
            slip = (1.0_wp - 2.0_wp*ns)*wet_q(i, j) + 2.0_wp*ns
            str_xy(i, j, k) = str_xy(i, j, k) + aniso_cross*shxx_at_q*h_q*slip
         end do
      end if

      ! ---- Phase 3a: u-face divergence ----
      ! diffu(i,j) = iareaCu·( idyCu·(dy2h·str_xx(i,j) − dy2h·str_xx(i-1,j))
      !                      + idxCu·(dx2q·str_xy(i,j+1) − dx2q·str_xy(i,j)) )
      !            / (h_u + h_neglect).   (signs verified to reduce to
      ! +A·∇²u on uniform-h all-wet; see module header.)
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx)
         du_visc(i, j, k) = wet_u(i, j)*iareaCu(i, j)* &
                            (idyCu(i, j)*(dy2h(i, j)*str_xx(i, j, k) - dy2h(i - 1, j)*str_xx(i - 1, j, k)) + &
                             idxCu(i, j)*(dx2q(i, j + 1)*str_xy(i, j + 1, k) - dx2q(i, j)*str_xy(i, j, k)))/ &
                            (0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k)) + H_NEGLECT)
      end do
      do concurrent(k=1:nz, j=1:ny)
         du_visc(1, j, k) = 0.0_wp
         du_visc(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx + 1)
         du_visc(i, 1, k) = 0.0_wp
         du_visc(i, ny, k) = 0.0_wp
      end do

      ! ---- Phase 3b: v-face divergence ----
      ! diffv(i,j) = iareaCv·( idxCv·(−dx2h·str_xx(i,j) + dx2h·str_xx(i,j-1))
      !                      + idyCv·(dy2q·str_xy(i+1,j) − dy2q·str_xy(i,j)) )
      !            / (h_v + h_neglect).  Tension enters with the
      ! opposite sign on the v-face (MOM6 diffv).
      do concurrent(k=1:nz, j=2:ny, i=2:nx - 1)
         dv_visc(i, j, k) = wet_v(i, j)*iareaCv(i, j)* &
                            (idyCv(i, j)*(dy2q(i + 1, j)*str_xy(i + 1, j, k) - dy2q(i, j)*str_xy(i, j, k)) - &
                             idxCv(i, j)*(dx2h(i, j)*str_xx(i, j, k) - dx2h(i, j - 1)*str_xx(i, j - 1, k)))/ &
                            (0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k)) + H_NEGLECT)
      end do
      do concurrent(k=1:nz, i=1:nx)
         dv_visc(i, 1, k) = 0.0_wp
         dv_visc(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1)
         dv_visc(1, j, k) = 0.0_wp
         dv_visc(nx, j, k) = 0.0_wp
      end do
   end subroutine hvisc_compute_stress

   subroutine ocean_horizontal_viscosity_apply_tendencies(this, ms, dt, no_wait)
      !! Forward-Euler accumulation of the viscous tendency onto the
      !! face velocities.  Shim — hoists `this%du_visc%data` etc. to
      !! the host before dispatching to the flat-impl.  Explicit shape
      !! dimensions are derived from `ms` here and passed as scalar
      !! args.
      !! `no_wait` (optional, default .false.): forwarded to the impl —
      !! when .true. the apply DC loops run on OpenACC queue 1 and the
      !! routine returns WITHOUT syncing, so the batched velocity-apply
      !! chain in `run_stage_split` `!$acc wait(1)`s ONCE.  Default ⇒
      !! blocking.  Not `pure` because of the async/wait directives.
      type(ocean_horizontal_viscosity_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: no_wait
      integer :: nx_cells, ny_cells, nz
      logical :: lwait
      lwait = .true.
      if (present(no_wait)) lwait = .not. no_wait
      nx_cells = size(ms%u_face_x_layer, 1) - 1
      ny_cells = size(ms%v_face_y_layer, 2) - 1
      nz = ms%nz_ml
      call hvisc_apply_impl(ms%u_face_x_layer, ms%v_face_y_layer, &
                            this%du_visc%data, this%dv_visc%data, dt, &
                            nx_cells, ny_cells, nz, lwait)
   end subroutine ocean_horizontal_viscosity_apply_tendencies

   subroutine hvisc_apply_impl(u_face, v_face, du_visc, dv_visc, dt, nx, ny, nz, lwait)
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(inout) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in)    :: dt
      logical, intent(in)    :: lwait
         !! .false. ⇒ leave the apply on queue 1 without syncing (batched).
      integer :: i, j, k
      !$acc kernels async(1)
      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         u_face(i, j, k) = u_face(i, j, k) + dt*du_visc(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         v_face(i, j, k) = v_face(i, j, k) + dt*dv_visc(i, j, k)
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine hvisc_apply_impl

   subroutine ocean_horizontal_viscosity_compute_ke_diss(this, ms)
      !! Fill `this%ke_diss` with the lateral-viscosity KE dissipation rate
      !! `Σ_k ρ_k h_k (u·du_visc + v·dv_visc)` at T-cell centres (kg/s³; ≤0
      !! where the viscosity removes KE).  MUST run AFTER
      !! `compute_tendencies` (du_visc fresh) and BEFORE the viscous apply,
      !! while `u_face`/`v_face` still hold the velocity the viscosity acted
      !! on.  No-op (and bit-identical) unless `compute_ke_diss` is set and
      !! the density field is live.  Feeds the MEKE frictional source.
      type(ocean_horizontal_viscosity_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      integer :: nx, ny, nz
      if (.not. this%compute_ke_diss) return
      if (.not. allocated(ms%rho_layer)) return
      nx = size(ms%u_face_x_layer, 1) - 1
      ny = size(ms%v_face_y_layer, 2) - 1
      nz = ms%nz_ml
      call hvisc_ke_diss_impl(nx, ny, nz, ms%u_face_x_layer, ms%v_face_y_layer, &
                              this%du_visc%data, this%dv_visc%data, &
                              ms%rho_layer, ms%h_layer, this%ke_diss)
   end subroutine ocean_horizontal_viscosity_compute_ke_diss

   pure subroutine hvisc_ke_diss_impl(nx, ny, nz, u_face, v_face, du_visc, dv_visc, &
                                      rho_layer, h_layer, ke_diss)
      !! C-grid KE budget: each face's `u·du_visc` rate is split half to each
      !! adjacent T-cell and depth-integrated with `ρ_k h_k`.  Race-free —
      !! every (i,j) writes only its own `ke_diss(i,j)`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: du_visc(nx + 1, ny, nz), dv_visc(nx, ny + 1, nz)
      real(wp), intent(in) :: rho_layer(nx, ny, nz), h_layer(nx, ny, nz)
      real(wp), intent(inout) :: ke_diss(nx, ny)
      integer :: i, j, k
      real(wp) :: acc, ke
      do concurrent(j=1:ny, i=1:nx) local(acc, ke, k)
         acc = 0.0_wp
         do k = 1, nz
            ke = 0.5_wp*(u_face(i, j, k)*du_visc(i, j, k) &
                         + u_face(i + 1, j, k)*du_visc(i + 1, j, k)) &
                 + 0.5_wp*(v_face(i, j, k)*dv_visc(i, j, k) &
                           + v_face(i, j + 1, k)*dv_visc(i, j + 1, k))
            acc = acc + rho_layer(i, j, k)*h_layer(i, j, k)*ke
         end do
         ke_diss(i, j) = acc
      end do
   end subroutine hvisc_ke_diss_impl

   pure function ocean_horizontal_viscosity_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the horizontal viscosity slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_horizontal_viscosity_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%du_visc%bytes() &
               + this%dv_visc%bytes() &
               + this%lap_u%bytes() &
               + this%lap_v%bytes() &
               + this%str_xx%bytes() &
               + this%str_xy%bytes() &
               + this%ah_t%bytes() &
               + this%ah_q%bytes() &
               + arr_bytes(this%ke_diss)
   end function ocean_horizontal_viscosity_bytes

end module rdb_ocean_horizontal_viscosity
