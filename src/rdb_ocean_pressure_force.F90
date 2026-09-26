!! Ocean pressure-force kernel state + the layered hydrostatic PGF variants.
module rdb_ocean_pressure_force
   !! Carries scheme-variant flags + reusable workspace for the
   !! layered hydrostatic pressure-force kernels on C-grid face arrays:
   !! Montgomery potential, finite-volume (lite / Wright / MOM6) and the
   !! NK=2 reduced-gravity form.
   !!
   !! The FV variants integrate pressure per column from the surface
   !! (p=0) down to the bed,
   !!   p_edge(k)   = p_edge(k+1) + g * rho_layer(k) * h_layer(k)
   !!   p_centre(k) = 0.5 * (p_edge(k) + p_edge(k+1))
   !! and difference cell-centred pressure across each face *at constant z*
   !! (pressure difference + `g*rho_face*dz_centre` correction) for the
   !! acceleration `du/dt = -(1/rho0) * dp/dx`.
   !!
   !! The Montgomery variant instead builds the Boussinesq Montgomery
   !! potential `M = p/rho0 + (g*rho/rho0)*z` — which is CONSTANT within a
   !! layer of uniform density — by a vertical recursion, and takes ONE
   !! horizontal difference of `M`.  See the OPGF_VARIANT_* constants below.
   !!
   !! MONT, FV_LITE and FV_WRIGHT are measured from the FREE SURFACE (z = 0
   !! at the surface, negative below), so they carry no barotropic
   !! `-g*grad(eta)` term: the split-explicit barotropic substep owns that.
   !! FV_MOM6 (`pa(nz+1) = rho_ref*g*eta_geo`) and GPRIME (`-g_FS*grad(eta)`
   !! in the top layer) DO carry it; the split sheds exactly that term from
   !! the barotropic forcing (`pgf_free_surface_gravity`,
   !! `set_fast_forcing_eta_pf`) and keeps the rest of the depth mean.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, H_DIV_EPS
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use rdb_eos, only: eos_wright_pgf_column_sweep_impl, eos_t, &
                      EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV
   use rdb_ocean_pgf_reconstruct, only: plm_edges_column, ppm_edges_column, &
                                        boole_dpa_intz_layer, boole_dpa_face, &
                                        boole_dpa_face_pcm, &
                                        PGF_RECON_PLM, PGF_RECON_PPM
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

   public :: ocean_pressure_force_t
   public :: ocean_pressure_force_compute
   public :: ocean_pressure_force_apply
   public :: parse_opgf_variant
   public :: gprime_nz_is_supported

   ! Pressure-force variant tags.
   integer, parameter, public :: OPGF_VARIANT_MONT = 1
      !! Boussinesq Montgomery-potential form (cheapest, layered, no EOS
      !! calls along the integration path, and no pressure stack at all).
      !!
      !! For a layer of horizontally uniform density the Boussinesq
      !! Montgomery potential
      !!   M = p/rho0 + rho_star*z,      rho_star = g*rho_layer/rho0
      !! is CONSTANT through the layer (moving up by `dz` costs
      !! `-g*rho*dz/rho0` of `p/rho0` and gains exactly `rho_star*dz`), and
      !!   -(1/rho0) * dp/dx|_z  ==  -dM/dx  along the layer.
      !! That is what makes ONE horizontal difference of `M` legitimate:
      !! `M` already carries the geopotential.
      !!
      !! Built by a vertical recursion, BOTTOM-UP here (k=1 bed, k=nz
      !! surface — MOM6 runs the same recursion top-down).  Seeded at the
      !! free surface, where the surface-relative interface height and the
      !! pressure are both zero:
      !!
      !!   e_edge(k)  = z_centre(k) + 0.5*h_layer(k)   ! TOP of layer k
      !!   M(nz)      = 0
      !!   M(k)       = M(k+1)
      !!              + (rho_star(k) - rho_star(k+1)) * e_edge(k+1)
      !!
      !! (`e_edge(k+1)` is the interface SHARED by layers k and k+1, where
      !! `p` and `z` agree, so the whole jump in `M` is the jump in
      !! `rho_star`.)  Then, per face,
      !!
      !!   PGF_x = -(M_R - M_L)*idxCu + (rho_star_R - rho_star_L)*z_eff*idxCu
      !!   z_eff = (e_L*h_R + e_R*h_L - h_L*h_R) / (h_L + h_R)
      !!
      !! There is NO `1/rho0` multiplying `dM`: `M` already has units of
      !! geopotential (m^2 s^-2), so its horizontal gradient IS an
      !! acceleration.  The `1/rho0` lives inside `rho_star`.
      !!
      !! The second term is load-bearing, not a refinement.  `-dM/dx` is
      !! the PGF only where `rho_layer` is horizontally uniform *within the
      !! layer*; where it is not, the exact relation picks up
      !! `+ z * d(rho_star)/dx`, and `z_eff` is the thickness-weighted
      !! height at which to evaluate it.  Drop it and the scheme gets the
      !! SIGN of a horizontal density contrast wrong.  On aligned columns
      !! (h_L == h_R) `z_eff` collapses to the arithmetic mean layer centre
      !! and the whole expression reduces ALGEBRAICALLY to FV_LITE.
      !!
      !! Exact at rest in an isopycnal (VCOORD_LAGRANGIAN) column over ANY
      !! bathymetry as long as the layer is present on both sides of the
      !! face: flat isopycnals make every `e_edge` horizontally uniform and
      !! every `rho_star` difference zero, so `M` is uniform and the PGF is
      !! zero to round-off.  Where an isopycnal layer has GROUNDED on one
      !! side the two columns share no common depth for that layer and the
      !! residual returns — `skip_nonoverlap` gates exactly that face, and
      !! honours `mont` for the same geometric reason it honours the FV
      !! variants.
   integer, parameter, public :: OPGF_VARIANT_FV_LITE = 2
      !! Finite-volume PGF with the z-position correction
      !!   PGF_x = -(1/rho_0) * [(p_R - p_L)/dx + g*rho_face*(z_R - z_L)/dx]
      !! that differences pressure *at constant z* rather than at constant
      !! layer index. Reduces to MONT for aligned columns; cancels the
      !! spurious bottom-current PGF when h_layer varies across columns.
      !! Reuses the layer-mean rho_layer from the EOS (no in-layer
      !! quadrature — use FV_WRIGHT for that).
   integer, parameter, public :: OPGF_VARIANT_GPRIME = 4
      !! Reduced-gravity / gprime PGF (reduced-gravity analogue).
      !! Per-layer acceleration
      !!   a_k = -g_FS · ∇η - Σ_{j>k} g'_j · ∇η_j_interface.
      !! NK = 2 only (asserted at init):
      !!   a_top = -g_FS · ∂(h_1 + h_2)/∂x
      !!   a_bot = -g_FS · ∂(h_1 + h_2)/∂x - g'_int · ∂h_1/∂x
      !! `pgf%gprime_gfs` is g_FS (m/s²), `pgf%gprime_gint` the single
      !! internal g'. No EOS, no Picard — densities fixed by GFS/GINT.
   integer, parameter, public :: OPGF_VARIANT_FV_WRIGHT = 3
      !! FV_LITE z-correction *plus* in-situ density along the
      !! integration path: each layer's centre density is re-evaluated by
      !! Wright at the pressure from a single Picard step (`ms%rho_layer`
      !! seeds the half-layer pressure). The pressure stack and `rho_face`
      !! both use the in-situ value, adding the compressibility piece
      !! FV_LITE misses; reduces to FV_LITE for incompressible water.
      !!
      !! Caller must have run `ocean_eos_compute` with
      !! `EOS_VARIANT_WRIGHT_97` first (the Picard seed is the nonlinear
      !! surface ρ, not a Boussinesq constant).
   integer, parameter, public :: OPGF_VARIANT_FV_MOM6 = 5
      !! Faithful layer-integrated FV-Bouss PGF (Boussinesq, per-layer
      !! `Rlay` density). Unlike FV_LITE (layer-centre pressure with a
      !! z-correction), uses layer-integrated pressure differences divided
      !! by face-averaged thickness:
      !!
      !!   PFu(I,j,k) = [ (pa·h + intz_dpa)_L − (pa·h + intz_dpa)_R
      !!                + (h_R − h_L) · intx_pa
      !!                − (e_bot_R − e_bot_L) · intx_dpa ]
      !!               · (2 · I_Rho0 · IdxCu) / (h_L + h_R + h_neglect)
      !!
      !! `pa` is the pressure anomaly stack relative to `rho_ref·g·z`,
      !! `intz_dpa` the in-layer vertical integral, `intx_pa/intx_dpa` the
      !! horizontal face integrals. The face-thickness divisor self-
      !! regulates bed-layer force at thin shelf-break cells (the mechanism
      !! FV_LITE+rho_init lacked, driving a 33× WBC bed overshoot).
      !!
      !! Density per layer from `ms%rho_layer`. `rho_ref` defaults to
      !! `rho0` so the surface layer's anomaly contribution vanishes.
      !! per-cell convention: our `k=1` is bed, MOM6's is surf.

   type :: ocean_pressure_force_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      integer :: variant = OPGF_VARIANT_MONT
         !! Active pressure-force scheme variant.  MONT is the default
         !! here AND the `&ocean_pgf_nml form` default, so a state built
         !! straight from this type and one built through the config
         !! agree.  MONT is general-purpose: valid over variable
         !! bathymetry and every vcoord.
      logical :: use_eos_along_path = .false.
         !! Convenience flag — re-evaluate density at each integration
         !! sample rather than at layer centres.  Phase 5d.
      logical :: scratch_gated = .false.
         !! ALLOCATION GATE for the variant-specific scratch buffers.
         !!
         !! `.false.` (default): `init` allocates all 16 buffers, whatever
         !! `variant` is — the historical behaviour every direct
         !! `pgf%init(grid, nz_ml=...)` call site (tests, benchmarks) relies
         !! on, since those set `variant` AFTER `init`.
         !!
         !! `.true.`: `init` allocates ONLY the buffers the active
         !! `variant` / `reconstruct_for_pressure` actually touch (see the
         !! per-buffer gate comments in `ocean_pressure_force_init`), saving
         !! up to 11 of the 16 (~3.5 GB at 1000x800x50).  Latched together
         !! with `variant` + `reconstruct_for_pressure` BEFORE `init(grid)`
         !! by `ocean_state_init_from_config` — the same conditional-
         !! allocation contract as the default-off closures and
         !! `continuity_t%windowed_advection`.
         !!
         !! Safety contract: every gated-off buffer is unreachable on its
         !! gated-off path because `ocean_pressure_force_compute` dispatches
         !! on `variant` at the TOP and `return`s out of each branch, and the
         !! one external `pgf%e_face` consumer (`compute_pbce`) `error stop`s
         !! unless `variant == OPGF_VARIANT_FV_MOM6`.  `enter_data` /
         !! `exit_data` / `bytes` all key off `allocated(...)`, so a gated-off
         !! buffer is neither mapped nor counted.
      real(wp) :: rho0 = 1035.0_wp
         !! BOUSSINESQ reference density (kg/m³) — the `ρ₀` that divides the
         !! pressure gradient into an acceleration, `du/dt = −(1/ρ₀)·∂p/∂x`.
         !! Read by EVERY variant (it is `inv_rho0` in the face passes and
         !! `g_over_rho0` in the Montgomery recursion), and by
         !! `compute_pbce` in the barotropic coupling.
         !!
         !! ASSIGNED FROM CONFIG by `configure_ocean_pgf`
         !! (`&ocean_ic_nml rho_0` → `eos%rho0`, the single ρ₀ of record —
         !! the EOS, the `eta_ib` surface-pressure seam, EPBL, kappa-shear,
         !! tidal mixing, MEKE/GM and the isopycnal slopes all take the same
         !! scalar).  The literal here is only the pre-configure default for
         !! direct `pgf%init(...)` call sites (tests, benchmarks) that never
         !! run the configure pass; it matches the `&ocean_ic_nml rho_0`
         !! default so a state built either way agrees.
         !!
         !! Host scalar: every read is host-side (into a local `inv_rho0`,
         !! or passed by value into a `*_impl`), so `configure_ocean_pgf`
         !! needs no `!$acc update device` — nothing reads it through the
         !! device-mapped `pgf` handle.

      ! ---- gprime / reduced-gravity knobs (OPGF_VARIANT_GPRIME) ----
      real(wp) :: gprime_gfs = 9.81_wp
         !! Free-surface reduced gravity (m/s²). Full physical g reduces to
         !! a standard free-surface model; reducing it slows the BT mode
         !! (wave speed `sqrt(gfs · H)`) for a slower BT CFL.
      real(wp) :: gprime_gint = 0.0098_wp
         !! Internal reduced gravity (m/s²) at the single layer-1/layer-2
         !! interface. Only used by `OPGF_VARIANT_GPRIME`.

      ! ---- FV_MOM6 reference state (OPGF_VARIANT_FV_MOM6) ----
      real(wp) :: rho_ref = 1035.0_wp
         !! ANOMALY reference density (kg/m³) subtracted from layer densities
         !! when building the `pa` pressure-anomaly stack (FV_MOM6), and the
         !! surface-layer `g·ρ_ref/ρ₀` in `compute_pbce`.
         !! `rho_ref = rho0 = ρ_surf` makes the surface layer's anomaly
         !! vanish; only denser bed layers contribute.
         !!
         !! A DISTINCT ROLE from `rho0` — `rho0` scales the gradient into an
         !! acceleration, `rho_ref` only shifts the baseline the anomaly is
         !! measured from — kept as a separate member so the two are never
         !! silently interchanged (MOM6 carries the same pair, `GV%Rho0` vs
         !! the `PressureForce_FV` `rho_ref`).  Both are nonetheless SOURCED
         !! FROM THE SAME CONFIGURED ρ₀ by `configure_ocean_pgf`
         !! (`&ocean_ic_nml rho_0` → `eos%rho0`): roundabout has no separate
         !! anomaly-reference knob, and `rho_ref ≠ rho0` would put a constant
         !! `g·(ρ_ref−ρ₀)/ρ₀` offset in `pa(top)` that the Boussinesq
         !! divisor no longer cancels.  Same host-scalar note as `rho0` —
         !! no device update needed.
      real(wp) :: h_neglect = 1.0e-10_wp
         !! Face-thickness floor in the FV_MOM6 denominator
         !! `(h_L + h_R + h_neglect)`. Prevents division by zero when both
         !! adjacent cells have vanishing bed layers.
      logical :: skip_nonoverlap = .false.
         !! Grounded-layer PGF gate (`&ocean_isopycnal_nml pgf_skip_nonoverlap`,
         !! set by the driver ONLY under VCOORD_LAGRANGIAN).  When `.true.` the
         !! face PGF is zeroed wherever the layer's z-extents in the two
         !! abutting columns do NOT overlap — a layer that has wedged out
         !! against the bed on one side, where the two-point Jacobian
         !! `Δp_centre + g·ρ_layer·Δz_centre` has no common depth to difference
         !! across and leaves `g·(ρ_layer − ρ̄_ambient)·∂z/∂x` of acceleration on
         !! a RESTING ocean.  `.false.` ⇒ the gate branch is never taken ⇒
         !! bit-identical.
         !!
         !! Honoured by `mont`, `fv_lite`, `fv_wright` and `fv_mom6`.  The gate
         !! reads `z_centre`; `ocean_pressure_force_init`'s allocation gate
         !! provides that buffer unconditionally for MONT and the two
         !! FV_LITE-family variants (all three also use it in the face passes)
         !! and, for FV_MOM6, exactly when this flag is set — so it MUST be
         !! latched before `init` (the ocean path does that in
         !! `ocean_state_init_from_config`, and `configure_ocean_pgf`
         !! re-checks that the latch did not drift).  `gprime` differences
         !! interface positions directly and has no such Jacobian, so it is
         !! the one variant left N/A (the driver warns).
         !!
         !! `mont` needs the gate for the SAME geometric reason the FV forms
         !! do, even though its face expression is not a two-point Jacobian:
         !! a grounded layer sits at the bed on the shallow side and at its
         !! flat-isopycnal height on the deep side, so the two `e_edge`
         !! values entering the `M` recursion are hundreds of metres apart
         !! and `M` stops being horizontally uniform at rest.
      logical :: mass_weight = .false.
         !! FV_MOM6 shelf-break `hWght` mass-weighting toggle. When `.true.`
         !! Pass-2's horizontal pressure integral biases the face density
         !! toward the thinner column at unequal-depth faces, cancelling the
         !! spurious bottom-layer shelf-break PGF. Equal-depth columns ⇒
         !! reduces exactly to the midpoint average ⇒ bit-identical.
      logical :: reconstruct_for_pressure = .false.
         !! FV_MOM6 in-layer T/S reconstruction toggle. `.false.` (default):
         !! layer-mean (PCM) density ⇒ bit-identical. `.true.`: per-layer
         !! `dpa` / `intz_dpa` from a 5-point Boole quadrature of a monotone
         !! PLM/PPM sub-layer T/S profile (Adcroft, Hallberg & Harrison 2008;
         !! White, Adcroft & Hallberg 2009), removing the spurious-PGF error
         !! on thick/sloped layers. Only consulted by FV_MOM6.
      integer :: recon_scheme = PGF_RECON_PLM
         !! In-layer reconstruction scheme: 1 = PLM, 2 = PPM. Only consulted
         !! when `reconstruct_for_pressure = .true.`.
      logical :: insitu_density = .true.
         !! FV_MOM6 constant-by-layer (PCM) density at its IN-SITU pressure
         !! (`&ocean_pgf_nml insitu_density`, MOM6 parity).  `.true.`
         !! (default): each layer's `dpa`/`intz_dpa` and the cross-face
         !! `intx_dpa`/`inty_dpa` are 5-point Boole quadratures of
         !! `EOS(T, S, p = -g·rho0·z)` with the layer-mean T/S — MOM6
         !! `int_density_dz_generic_pcm` (`compute_fv_mom6_insitu_pcm_impl`).
         !! `.false.`: the legacy PCM integral of `ms%rho_layer`, a
         !! POTENTIAL density at the single `&ocean_eos_nml p_ref`, which
         !! drops the pressure dependence of the horizontal density
         !! gradient below the reference level.  Consulted only by FV_MOM6
         !! with `reconstruct_for_pressure = .false.`, an EOS handle and
         !! T/S, and only for a PRESSURE-DEPENDENT EOS (Wright, Roquet):
         !! for the linear EOS in-situ and potential density coincide, so
         !! the legacy path runs and answers are bit-identical.
      logical :: p_top_in_bc = .false.
         !! FV_MOM6 top-of-column pressure in the surface boundary
         !! condition (`&ocean_pgf_nml p_top_in_bc`). `.false.` (default):
         !! `pa(nz+1) = rho_ref·g·eta_geo`, bit-identical. `.true.`: the
         !! load `multilayer_state_t%p_top` (Pa) is ADDED there, so the
         !! pressure stack measures down from the loaded surface —
         !! `pa(nz+1) = rho_ref·g·eta_geo + p_top`. Only consulted by
         !! FV_MOM6 (both the PCM and the `reconstruct_for_pressure`
         !! branch); fail-loud at configure for any other variant, which
         !! carries no injectable `pa` stack.
         !!
         !! A DEPTH-UNIFORM `p_top` perturbs every layer's `PFu` by the
         !! SAME `−(1/ρ₀)∇p_top` (Theorem 1 in the Pass-1 docstring
         !! below), and the split solver replaces the depth mean of the
         !! layer PGF with the barotropic solution, so the baroclinic
         !! operator does not see it and this does NOT double-count the
         !! `eta_forcing` seam. See the `p_top` seam contract in
         !! `src/core/ocean/README.md`.
      real(wp) :: gfs_scale = 1.0_wp
         !! Free-surface gravity scaling (= GFS / G_EARTH). Default 1.0 ⇒
         !! pure FV_MOM6, bit-identical. When < 1, Pass 5 applies the
         !! Montgomery `dM` correction subtracting
         !! `(1 - gfs_scale)·(g/ρ₀)·ρ_surf·∇η` from every layer's PGF
         !! (depth-independent, so the BT mass-flux invariant survives); the
         !! driver also drops `bt_work%g_bt` to `gfs_scale·GRAVITY`.
         !! Combined ⇒ wave speed `sqrt(gfs_scale · g · H)`.

      ! ---- Bathymetry (copy of state%barotropic%b) ----
      ! Used by the gprime PGF to recover ∇η from ∇(sum h_layer) - ∇b;
      ! without it variable-bathymetry runs see a spurious PGF dominated
      ! by ∇H_bathy (~10⁴× larger than ∇η at shelf-breaks). Default zero =
      ! flat bed. Set by the driver via `set_bathymetry()`.
      real(wp), allocatable :: b(:, :)

      ! ---- Workspace ----
      type(scratch_3d_buffer_t) :: p_edge
         !! Layer-edge pressure stack at cell centres.  Shape
         !! (nx, ny, nz_ml+1).  k=1 bed, k=nz_ml+1 surface.
         !! Reused across the two horizontal-gradient passes.
      type(scratch_3d_buffer_t) :: z_centre
         !! Layer-centre z at cell centres, from the free surface downward
         !! (z=0 surface, negative below). Shape (nx, ny, nz_ml). Filled for
         !! MONT and the FV_LITE-family branches. Surface-relative (NOT
         !! bed-relative) so the σ-coord Jacobian cancels the spurious
         !! cross-bathymetry pressure gradient — and so that neither form
         !! double-counts the barotropic `-g·∇η` the BT substep already
         !! carries. MONT reads it as the interface height
         !! `e_edge(k) = z_centre(k) + 0.5·h_layer(k)`.
      type(scratch_3d_buffer_t) :: mont_M
         !! Boussinesq Montgomery potential `M` at layer centres (m² s⁻²).
         !! Shape (nx, ny, nz_ml). Written by the MONT column recursion and
         !! read by its two face passes; no other variant touches it.
      type(scratch_3d_buffer_t) :: rho_insitu
         !! In-situ layer-centre density from the Wright
         !! column-sweep Picard step.  Shape (nx, ny, nz_ml).
         !! Populated only by the FV_WRIGHT branch; the other
         !! variants source ρ from `ms%rho_layer`.
      type(scratch_3d_buffer_t) :: dpdx_face
         !! East-face PGF acceleration: -(1/rho0) * dp/dx at
         !! u-face.  Shape (nx+1, ny, nz_ml).  Apply adds dt*dpdx
         !! to u_face_x_layer.
      type(scratch_3d_buffer_t) :: dpdy_face
         !! North-face PGF acceleration.  Shape (nx, ny+1, nz_ml).

      ! ---- FV_MOM6 scratch (OPGF_VARIANT_FV_MOM6) ----
      !! Interface heights `e` (positive-up).  `e(:,:,1)` is the bed
      !! (= -b), `e(:,:,nz+1)` is the free surface (= η).  Same
      !! sign convention as MOM6, but bottom-up indexing to match the
      !! Roundabout bed-up layer convention.
      type(scratch_3d_buffer_t) :: e_face
         !! Interface heights at cell centres, shape (nx, ny, nz+1).
      !! Pressure anomaly stack at interfaces, units Pa.  Relative to
      !! the `rho_ref · g · z` baseline so the surface-pressure
      !! contribution is just `pa(top) = rho_ref · g · η`.  Built by
      !! marching down from the surface; `pa(k) − pa(k+1) = dpa(k)`
      !! where `dpa(k) = (Rlay(k) − rho_ref) · g · h(k)`.
      type(scratch_3d_buffer_t) :: pa
         !! Pressure anomaly at interfaces, shape (nx, ny, nz+1).
      !! Per-layer vertical integral of `dpa` from layer top inward.
      !! For Boussinesq Rlay path: `intz_dpa(k) = 0.5 · (Rlay(k) −
      !! rho_ref) · g · h(k)²` (mid-point rule).
      type(scratch_3d_buffer_t) :: intz_dpa
         !! Per-layer ∫ dpa dz, shape (nx, ny, nz), units Pa·m.
      ! Horizontal integrals at faces — average of the two adjacent
      ! cells.  `intx_pa(K) = 0.5·(pa_L + pa_R)` at interface K;
      ! `intx_dpa(k) = 0.5·(Rlay(k) − rho_ref) · g · (h_L + h_R)`
      ! within layer k.  Surface BC fixes intx_pa at top; deeper
      ! values come from `intx_pa(K+1) = intx_pa(K) + intx_dpa(k)`.
      type(scratch_3d_buffer_t) :: intx_pa
         !! u-face ∫ pa across x, shape (nx+1, ny, nz+1), units Pa·m.
      type(scratch_3d_buffer_t) :: inty_pa
         !! v-face ∫ pa across y, shape (nx, ny+1, nz+1), units Pa·m.
      type(scratch_3d_buffer_t) :: intx_dpa
         !! u-face ∫ dpa, shape (nx+1, ny, nz), units Pa·m.
      type(scratch_3d_buffer_t) :: inty_dpa
         !! v-face ∫ dpa, shape (nx, ny+1, nz), units Pa·m.

      ! ---- In-layer reconstruction scratch (reconstruct_for_pressure) ----
      !! Per-column PLM/PPM top (shallower) and bottom (deeper) edge
      !! values of the layer-mean T and S, shape (nx, ny, nz).  Filled
      !! by `compute_fv_mom6_reconstruct_impl`'s edge-build pass; consumed
      !! by the 5-point Boole quadrature.  Allocated only when
      !! `reconstruct_for_pressure` (or `scratch_gated = .false.`).
      type(scratch_3d_buffer_t) :: recon_T_t
      type(scratch_3d_buffer_t) :: recon_T_b
      type(scratch_3d_buffer_t) :: recon_S_t
      type(scratch_3d_buffer_t) :: recon_S_b
   contains
      procedure, non_overridable :: init => ocean_pressure_force_init
      procedure, non_overridable :: destroy => ocean_pressure_force_destroy
      procedure, non_overridable :: enter_data => ocean_pressure_force_enter_data
      procedure, non_overridable :: exit_data => ocean_pressure_force_exit_data
      procedure, non_overridable :: set_bathymetry => ocean_pressure_force_set_bathymetry
      procedure, non_overridable :: bytes => ocean_pressure_force_bytes
   end type ocean_pressure_force_t

contains

   pure function gprime_nz_is_supported(variant, nz) result(ok)
      !! `.true.` unless `variant == OPGF_VARIANT_GPRIME` with `nz /= 2`
      !! (PR-6 fail-loud).  The reduced-gravity gprime PGF hard-writes
      !! ONLY k=1 (bottom) and k=2 (top) — with `nz > 2` layers k=3..nz
      !! carry zero pressure gradient AND the "top" branch lands on layer
      !! 2 of nz (the abyss under the bottom-up convention); with `nz < 2`
      !! the kernel early-returns leaving the whole PGF zero.  The guard
      !! is gprime-specific — every other variant supports general nz, so
      !! this returns `.true.` for them regardless of `nz`.  Wired into
      !! `validate_config` against `cfg%nz_layers` (config-time).
      integer, intent(in) :: variant, nz
      logical :: ok
      ok = (variant /= OPGF_VARIANT_GPRIME) .or. (nz == 2)
   end function gprime_nz_is_supported

   subroutine ocean_pressure_force_init(this, grid, nz_ml)
      !! Allocate the scratch buffers.  Default nz_ml=1 keeps the
      !! barotropic-only path constructible; passing nz_ml sizes them
      !! for the multilayer kernel.
      !!
      !! ALLOCATION GATE (`this%scratch_gated`, see the type docstring):
      !! when `.true.` only the buffers the active `variant` /
      !! `reconstruct_for_pressure` can actually reach are allocated.  The
      !! gate defaults `.false.` so a bare `pgf%init(...)` (every direct
      !! test/benchmark call site) keeps the historical allocate-everything
      !! behaviour.  Each gate's unreachability proof is stated inline.
      class(ocean_pressure_force_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz
      logical :: gate, need_p_edge, need_z_centre, need_rho_insitu
      logical :: need_mont_M, need_fv_mom6, need_recon

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      gate = this%scratch_gated
      ! `p_edge`: filled by Pass 1 and read by Pass 2/3.  The MONT, GPRIME and
      ! FV_MOM6 branches of `ocean_pressure_force_compute` `return` or branch
      ! away BEFORE Pass 1, so none of them ever touches it.  MONT builds the
      ! Montgomery potential straight from `rho_layer` + interface heights and
      ! needs no pressure stack at all.  No other module references it.
      need_p_edge = .not. gate .or. &
                    (this%variant == OPGF_VARIANT_FV_LITE .or. &
                     this%variant == OPGF_VARIANT_FV_WRIGHT)
      ! `z_centre`: written by Pass 1b and read by the MONT / FV_LITE /
      ! FV_WRIGHT face passes — every one of those sites sits inside a branch
      ! testing for exactly those three variants.  The Pass-4 grounded-layer
      ! gate reads it too, and that pass runs for FV_MOM6 as well — but ONLY
      ! when `skip_nonoverlap` is on, which is a pure config-time decision
      ! (vcoord + namelist) that the caller latches before `init`.  So FV_MOM6
      ! gets the buffer exactly when something will read it, and an ungated
      ! FV_MOM6 run still pays nothing.  No other module references it.
      need_z_centre = .not. gate .or. &
                      (this%variant == OPGF_VARIANT_MONT .or. &
                       this%variant == OPGF_VARIANT_FV_LITE .or. &
                       this%variant == OPGF_VARIANT_FV_WRIGHT .or. &
                       (this%variant == OPGF_VARIANT_FV_MOM6 .and. &
                        this%skip_nonoverlap))
      ! `mont_M`: written and read by the MONT branch alone.  No other
      ! variant, and no other module, references it.
      need_mont_M = .not. gate .or. (this%variant == OPGF_VARIANT_MONT)
      ! `rho_insitu`: written by the Pass-1 FV_WRIGHT branch (both the EOS
      ! column sweep and its no-tracer fallback) and read by the Pass 2/3
      ! FV_WRIGHT branches only.  No other module references it.
      need_rho_insitu = .not. gate .or. (this%variant == OPGF_VARIANT_FV_WRIGHT)
      ! FV_MOM6 stack: passed only to `compute_fv_mom6_impl` /
      ! `compute_fv_mom6_reconstruct_impl`, both inside the
      ! `if (variant == FV_MOM6)` branch that `return`s.  The one external
      ! reader, `compute_pbce` (rdb_barotropic_coupling), opens with an
      ! `error stop` unless `variant == OPGF_VARIANT_FV_MOM6`.
      need_fv_mom6 = .not. gate .or. (this%variant == OPGF_VARIANT_FV_MOM6)
      ! Reconstruction edge scratch: passed only to
      ! `compute_fv_mom6_reconstruct_impl`, whose call site additionally
      ! requires `reconstruct_for_pressure` (configure_ocean_pgf `error
      ! stop`s if that knob is set with any variant other than FV_MOM6).
      need_recon = .not. gate .or. this%reconstruct_for_pressure

      ! Cell-centred edge stack: (nx, ny, nz+1)
      if (need_p_edge) call this%p_edge%init(nx, ny, nz + 1, "ocean_pgf_p_edge")
      ! Layer-centre z-coordinate: (nx, ny, nz)
      if (need_z_centre) call this%z_centre%init(nx, ny, nz, "ocean_pgf_z_centre")
      ! Montgomery potential at layer centres: (nx, ny, nz)
      if (need_mont_M) call this%mont_M%init(nx, ny, nz, "ocean_pgf_mont_M")
      ! In-situ density at layer centres: (nx, ny, nz)
      if (need_rho_insitu) call this%rho_insitu%init(nx, ny, nz, "ocean_pgf_rho_insitu")
      ! East-face: (nx+1, ny, nz) — same shape as u_face_x_layer.  Every
      ! variant writes these (they ARE the PGF output), so never gated.
      call this%dpdx_face%init(nx + 1, ny, nz, "ocean_pgf_dpdx_face")
      ! North-face: (nx, ny+1, nz)
      call this%dpdy_face%init(nx, ny + 1, nz, "ocean_pgf_dpdy_face")

      ! FV_MOM6 scratch (interface heights, pa stack, per-layer
      ! integrals, per-face horizontal integrals).
      if (need_fv_mom6) then
         call this%e_face%init(nx, ny, nz + 1, "ocean_pgf_fv_mom6_e_face")
         call this%pa%init(nx, ny, nz + 1, "ocean_pgf_fv_mom6_pa")
         call this%intz_dpa%init(nx, ny, nz, "ocean_pgf_fv_mom6_intz_dpa")
         call this%intx_pa%init(nx + 1, ny, nz + 1, "ocean_pgf_fv_mom6_intx_pa")
         call this%inty_pa%init(nx, ny + 1, nz + 1, "ocean_pgf_fv_mom6_inty_pa")
         call this%intx_dpa%init(nx + 1, ny, nz, "ocean_pgf_fv_mom6_intx_dpa")
         call this%inty_dpa%init(nx, ny + 1, nz, "ocean_pgf_fv_mom6_inty_dpa")
      end if

      ! In-layer reconstruction edge scratch (nx, ny, nz).
      if (need_recon) then
         call this%recon_T_t%init(nx, ny, nz, "ocean_pgf_recon_T_t")
         call this%recon_T_b%init(nx, ny, nz, "ocean_pgf_recon_T_b")
         call this%recon_S_t%init(nx, ny, nz, "ocean_pgf_recon_S_t")
         call this%recon_S_b%init(nx, ny, nz, "ocean_pgf_recon_S_b")
      end if

      ! Bathymetry copy.  Default zero = flat bed.  Driver overwrites
      ! via `set_bathymetry` after `state%barotropic%b` is populated.
      allocate (this%b(nx, ny), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_pressure_force_init

   subroutine ocean_pressure_force_destroy(this)
      class(ocean_pressure_force_t), intent(inout) :: this
      this%is_init = .false.
      call this%p_edge%destroy()
      call this%z_centre%destroy()
      call this%mont_M%destroy()
      call this%rho_insitu%destroy()
      call this%dpdx_face%destroy()
      call this%dpdy_face%destroy()
      call this%e_face%destroy()
      call this%pa%destroy()
      call this%intz_dpa%destroy()
      call this%intx_pa%destroy()
      call this%inty_pa%destroy()
      call this%intx_dpa%destroy()
      call this%inty_dpa%destroy()
      call this%recon_T_t%destroy()
      call this%recon_T_b%destroy()
      call this%recon_S_t%destroy()
      call this%recon_S_b%destroy()
      if (allocated(this%b)) deallocate (this%b)
   end subroutine ocean_pressure_force_destroy

   subroutine ocean_pressure_force_set_bathymetry(this, b)
      !! Copy `b(:, :)` into `this%b` on the host. Must be called BEFORE
      !! `enter_data` (the device copy is taken from the host values).
      !! Issues no `update device`; to refresh post `enter_data` the
      !! caller must issue `!$acc update device(this%b)` itself.
      class(ocean_pressure_force_t), intent(inout) :: this
      real(wp), intent(in) :: b(:, :)
      integer :: nx, ny, i, j
      nx = size(this%b, 1)
      ny = size(this%b, 2)
      if (size(b, 1) /= nx .or. size(b, 2) /= ny) then
         error stop "ocean_pressure_force_set_bathymetry: shape mismatch"
      end if
      do j = 1, ny
         do i = 1, nx
            this%b(i, j) = b(i, j)
         end do
      end do
   end subroutine ocean_pressure_force_set_bathymetry

   subroutine ocean_pressure_force_enter_data(this)
      class(ocean_pressure_force_t), intent(inout) :: this
      select type (this)
      type is (ocean_pressure_force_t)
         call ocean_pressure_force_enter_data_impl(this)
      end select
   end subroutine ocean_pressure_force_enter_data

   subroutine ocean_pressure_force_enter_data_impl(this)
      type(ocean_pressure_force_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%p_edge)
      call scratch_3d_buffer_enter_data_impl(this%z_centre)
      call scratch_3d_buffer_enter_data_impl(this%mont_M)
      call scratch_3d_buffer_enter_data_impl(this%rho_insitu)
      call scratch_3d_buffer_enter_data_impl(this%dpdx_face)
      call scratch_3d_buffer_enter_data_impl(this%dpdy_face)
      call scratch_3d_buffer_enter_data_impl(this%e_face)
      call scratch_3d_buffer_enter_data_impl(this%pa)
      call scratch_3d_buffer_enter_data_impl(this%intz_dpa)
      call scratch_3d_buffer_enter_data_impl(this%intx_pa)
      call scratch_3d_buffer_enter_data_impl(this%inty_pa)
      call scratch_3d_buffer_enter_data_impl(this%intx_dpa)
      call scratch_3d_buffer_enter_data_impl(this%inty_dpa)
      call scratch_3d_buffer_enter_data_impl(this%recon_T_t)
      call scratch_3d_buffer_enter_data_impl(this%recon_T_b)
      call scratch_3d_buffer_enter_data_impl(this%recon_S_t)
      call scratch_3d_buffer_enter_data_impl(this%recon_S_b)
      !$acc enter data copyin(this%b)
   end subroutine ocean_pressure_force_enter_data_impl

   subroutine ocean_pressure_force_exit_data(this)
      class(ocean_pressure_force_t), intent(inout) :: this
      select type (this)
      type is (ocean_pressure_force_t)
         call ocean_pressure_force_exit_data_impl(this)
      end select
   end subroutine ocean_pressure_force_exit_data

   subroutine ocean_pressure_force_exit_data_impl(this)
      type(ocean_pressure_force_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%p_edge)
      call scratch_3d_buffer_exit_data_impl(this%z_centre)
      call scratch_3d_buffer_exit_data_impl(this%mont_M)
      call scratch_3d_buffer_exit_data_impl(this%rho_insitu)
      call scratch_3d_buffer_exit_data_impl(this%dpdx_face)
      call scratch_3d_buffer_exit_data_impl(this%dpdy_face)
      call scratch_3d_buffer_exit_data_impl(this%e_face)
      call scratch_3d_buffer_exit_data_impl(this%pa)
      call scratch_3d_buffer_exit_data_impl(this%intz_dpa)
      call scratch_3d_buffer_exit_data_impl(this%intx_pa)
      call scratch_3d_buffer_exit_data_impl(this%inty_pa)
      call scratch_3d_buffer_exit_data_impl(this%intx_dpa)
      call scratch_3d_buffer_exit_data_impl(this%inty_dpa)
      call scratch_3d_buffer_exit_data_impl(this%recon_T_t)
      call scratch_3d_buffer_exit_data_impl(this%recon_T_b)
      call scratch_3d_buffer_exit_data_impl(this%recon_S_t)
      call scratch_3d_buffer_exit_data_impl(this%recon_S_b)
      !$acc exit data delete(this%b)
   end subroutine ocean_pressure_force_exit_data_impl

   pure subroutine ocean_pressure_force_compute(grid, metrics, pgf, ms, eos)
      !! Compute the hydrostatic pressure-gradient acceleration at every
      !! C-grid face. Variants (see the OPGF_VARIANT_* / GPRIME / FV_MOM6
      !! constants for the per-variant formulas): MONT (layer-mean ρgh),
      !! FV_LITE (+ z-correction), FV_WRIGHT (+ in-situ Wright density),
      !! GPRIME, FV_MOM6.
      !!
      !! Pass layout: (1b) z_centre from h_layer (MONT + FV variants);
      !! then either (M1) the Montgomery column recursion and (M2/M3) its
      !! face passes, or (1) the column pressure sweep filling p_edge
      !! (+ rho_insitu for FV_WRIGHT) and (2/3) the FV east/north-face
      !! acceleration; finally (4) the grounded-layer gate.
      !!
      !! `ms%rho_layer` must be up to date — call the EOS kernel first.
      !! FV_WRIGHT needs `EOS_VARIANT_WRIGHT_97` for a good Picard seed.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics. The u-face gradient divides by
         !! `idxCu(i,j)`, the v-face by `idyCv(i,j)`. On uniform Cartesian
         !! `idxCu == 1/dx` bitwise (byte-identical to scalar inv_dx/inv_dy).
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(multilayer_state_t), intent(in) :: ms
      type(eos_t), intent(in), optional :: eos
         !! EOS handle — REQUIRED when `pgf%reconstruct_for_pressure` is
         !! on (the in-layer Boole quadrature evaluates the EOS at each
         !! sub-point).  Optional so the legacy PCM call sites (and the
         !! non-reconstruct variant tests) need not thread it through.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: inv_rho0, g_over_rho0
      real(wp) :: p_centre_left, p_centre_right
      real(wp) :: p_centre_below, p_centre_above
      real(wp) :: rho_face, z_correction, z_running
      real(wp) :: h_l, h_r, e_l, e_r, z_eff, drho_star

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      inv_rho0 = 1.0_wp/pgf%rho0

      ! ---- gprime / reduced-gravity branch (Tier-1: NK = 2 only) ----
      if (pgf%variant == OPGF_VARIANT_GPRIME) then
         call compute_gprime_impl(ms%h_layer, pgf%b, &
                                  pgf%dpdx_face%data, pgf%dpdy_face%data, &
                                  pgf%gprime_gfs, pgf%gprime_gint, &
                                  metrics%idxCu, metrics%idyCv, nx, ny, nz)
         return
      end if

      ! ---- Pass 1b: per-column z_centre from h_layer ----
      ! z_centre(:, :, k) = physical z of the layer-k mid-depth,
      ! measured from the free surface downward (z=0 at surface,
      ! z negative below).  Computed by summing h_layer top-down so
      ! the surface is the reference point regardless of column
      ! bathymetry.  This is the load-bearing change that makes the
      ! FV_LITE / FV_WRIGHT Jacobian cancel the cross-bathymetry
      ! pressure gradient: a bed-relative z_centre would inject a
      ! spurious `-g·ρ·dH/dx` residual at every shelf-break face.
      !
      ! Hoisted AHEAD of the variant branches because it depends on nothing
      ! but `h_layer`.  FV_LITE / FV_WRIGHT consume it in Passes 2/3; MONT
      ! recovers the layer-TOP interface height `z_centre + h/2` from it for
      ! both the `M` recursion and the face `z_eff`.  FV_MOM6's own face
      ! assembly never touches it, so on that path it is filled ONLY to feed
      ! the Pass-4 grounded-layer gate — which is also exactly when
      ! `ocean_pressure_force_init`'s allocation gate provides the buffer.
      if (pgf%variant == OPGF_VARIANT_MONT .or. &
          pgf%variant == OPGF_VARIANT_FV_LITE .or. &
          pgf%variant == OPGF_VARIANT_FV_WRIGHT .or. &
          (pgf%variant == OPGF_VARIANT_FV_MOM6 .and. pgf%skip_nonoverlap)) then
         do concurrent(j=1:ny, i=1:nx) local(z_running)
            z_running = 0.0_wp
            do k = nz, 1, -1
               pgf%z_centre%data(i, j, k) = z_running - 0.5_wp*ms%h_layer(i, j, k)
               z_running = z_running - ms%h_layer(i, j, k)
            end do
         end do
      end if

      ! ---- MONT branch — Boussinesq Montgomery potential ----------------
      ! Pass M1 (per column): the vertical recursion that builds `M`.  Pass
      ! M2/M3 (per face): ONE horizontal difference of `M`, plus the
      ! horizontal-density term.  No pressure stack is built on this path.
      ! See the OPGF_VARIANT_MONT docstring for the derivation.
      if (pgf%variant == OPGF_VARIANT_MONT) then
         g_over_rho0 = GRAVITY*inv_rho0

         ! ---- Pass M1: Montgomery potential per column -------------------
         ! Seeded at the free surface: `e_edge(nz+1) = 0` (z_centre is
         ! surface-relative) and the surface pressure is zero, so
         ! `M(nz) = p/rho0 + rho_star(nz)*e_edge(nz+1)` is identically zero
         ! in EVERY column.  That is not a loss: the barotropic `-g*grad(eta)`
         ! it would otherwise carry is the BT substep's job, and adding it
         ! here would double-count it.  Accumulated straight into the array
         ! (no `local` reassigned in the k-loop), mirroring the Pass-1
         ! pressure sweep.
         do concurrent(j=1:ny, i=1:nx)
            pgf%mont_M%data(i, j, nz) = 0.0_wp
            do k = nz - 1, 1, -1
               ! `e_edge(k+1)` — the interface SHARED by layers k and k+1,
               ! i.e. the TOP of layer k — is where `p` and `z` agree between
               ! the two layers, so the whole jump in M is the jump in
               ! rho_star.  Recovered from the layer-k centre.
               pgf%mont_M%data(i, j, k) = pgf%mont_M%data(i, j, k + 1) + &
                                          g_over_rho0*(ms%rho_layer(i, j, k) - &
                                                       ms%rho_layer(i, j, k + 1))* &
                                          (pgf%z_centre%data(i, j, k) + &
                                           0.5_wp*ms%h_layer(i, j, k))
            end do
         end do

         ! ---- Pass M2: east-face acceleration ----------------------------
         ! `-dM/dx` plus the horizontal-density term `+ z_eff * d(rho_star)/dx`.
         ! `z_eff` is the thickness-weighted height at which `M`'s two
         ! column anchors are reconciled; on aligned columns (h_L == h_R) it
         ! is exactly the mean layer centre and the pair reduces
         ! ALGEBRAICALLY to FV_LITE.  `H_DIV_EPS` is pure 1/0 armour for a
         ! face between two fully-vanished layers (numerator is then zero
         ! too, so the face value is a clean zero).
         do concurrent(k=1:nz, j=1:ny, i=2:nx) &
            local(h_l, h_r, e_l, e_r, z_eff, drho_star)
            h_l = ms%h_layer(i - 1, j, k)
            h_r = ms%h_layer(i, j, k)
            e_l = pgf%z_centre%data(i - 1, j, k) + 0.5_wp*h_l
            e_r = pgf%z_centre%data(i, j, k) + 0.5_wp*h_r
            z_eff = (e_l*h_r + e_r*h_l - h_l*h_r)/(h_l + h_r + H_DIV_EPS)
            drho_star = g_over_rho0*(ms%rho_layer(i, j, k) - ms%rho_layer(i - 1, j, k))
            pgf%dpdx_face%data(i, j, k) = &
               (-(pgf%mont_M%data(i, j, k) - pgf%mont_M%data(i - 1, j, k)) &
                + drho_star*z_eff)*metrics%idxCu(i, j)
         end do
         do concurrent(k=1:nz, j=1:ny)
            pgf%dpdx_face%data(1, j, k) = 0.0_wp
            pgf%dpdx_face%data(nx + 1, j, k) = 0.0_wp
         end do

         ! ---- Pass M3: north-face acceleration ---------------------------
         do concurrent(k=1:nz, j=2:ny, i=1:nx) &
            local(h_l, h_r, e_l, e_r, z_eff, drho_star)
            h_l = ms%h_layer(i, j - 1, k)
            h_r = ms%h_layer(i, j, k)
            e_l = pgf%z_centre%data(i, j - 1, k) + 0.5_wp*h_l
            e_r = pgf%z_centre%data(i, j, k) + 0.5_wp*h_r
            z_eff = (e_l*h_r + e_r*h_l - h_l*h_r)/(h_l + h_r + H_DIV_EPS)
            drho_star = g_over_rho0*(ms%rho_layer(i, j, k) - ms%rho_layer(i, j - 1, k))
            pgf%dpdy_face%data(i, j, k) = &
               (-(pgf%mont_M%data(i, j, k) - pgf%mont_M%data(i, j - 1, k)) &
                + drho_star*z_eff)*metrics%idyCv(i, j)
         end do
         do concurrent(k=1:nz, i=1:nx)
            pgf%dpdy_face%data(i, 1, k) = 0.0_wp
            pgf%dpdy_face%data(i, ny + 1, k) = 0.0_wp
         end do

         ! ---- FV_MOM6 branch — faithful port of MOM6 PressureForce_FV_Bouss ----
         ! Layer-integrated pressure differences with face-thickness
         ! divisor.  See OPGF_VARIANT_FV_MOM6 doc above.
      else if (pgf%variant == OPGF_VARIANT_FV_MOM6) then
         if (pgf%reconstruct_for_pressure .and. present(eos) .and. &
             ms%idx_salinity > 0 .and. ms%idx_temperature > 0) then
            ! In-layer PLM/PPM reconstruction: build per-layer Boole
            ! `dpa`/`intz_dpa` from a monotone sub-layer T/S profile, then
            ! reuse the unchanged FV_MOM6 face assembly.
            call compute_fv_mom6_reconstruct_impl(ms%h_layer, &
                                                  ms%tracers(ms%idx_salinity)%hTr, &
                                                  ms%tracers(ms%idx_temperature)%hTr, &
                                                  pgf%b, eos, &
                                                  pgf%recon_S_t%data, pgf%recon_S_b%data, &
                                                  pgf%recon_T_t%data, pgf%recon_T_b%data, &
                                                  pgf%e_face%data, pgf%pa%data, &
                                                  pgf%intz_dpa%data, &
                                                  pgf%intx_pa%data, pgf%inty_pa%data, &
                                                  pgf%intx_dpa%data, pgf%inty_dpa%data, &
                                                  pgf%dpdx_face%data, pgf%dpdy_face%data, &
                                                  pgf%rho0, pgf%rho_ref, pgf%h_neglect, &
                                                  pgf%gfs_scale, pgf%recon_scheme, &
                                                  ms%p_top, pgf%p_top_in_bc, &
                                                  metrics%idxCu, metrics%idyCv, nx, ny, nz)
         else if (use_insitu_pcm(pgf, ms, eos)) then
            ! Constant-by-layer T/S, density at the in-situ pressure
            ! (MOM6 `int_density_dz_generic_pcm`).  See `insitu_density`.
            call compute_fv_mom6_insitu_pcm_impl(ms%h_layer, &
                                                 ms%tracers(ms%idx_salinity)%hTr, &
                                                 ms%tracers(ms%idx_temperature)%hTr, &
                                                 pgf%b, eos, &
                                                 pgf%e_face%data, pgf%pa%data, &
                                                 pgf%intz_dpa%data, &
                                                 pgf%intx_pa%data, pgf%inty_pa%data, &
                                                 pgf%intx_dpa%data, pgf%inty_dpa%data, &
                                                 pgf%dpdx_face%data, pgf%dpdy_face%data, &
                                                 pgf%rho0, pgf%rho_ref, pgf%h_neglect, &
                                                 pgf%gfs_scale, pgf%mass_weight, &
                                                 ms%p_top, pgf%p_top_in_bc, &
                                                 metrics%idxCu, metrics%idyCv, nx, ny, nz)
         else
            call compute_fv_mom6_impl(ms%h_layer, ms%rho_layer, pgf%b, &
                                      pgf%e_face%data, pgf%pa%data, &
                                      pgf%intz_dpa%data, &
                                      pgf%intx_pa%data, pgf%inty_pa%data, &
                                      pgf%intx_dpa%data, pgf%inty_dpa%data, &
                                      pgf%dpdx_face%data, pgf%dpdy_face%data, &
                                      pgf%rho0, pgf%rho_ref, pgf%h_neglect, &
                                      pgf%gfs_scale, pgf%mass_weight, &
                                      ms%p_top, pgf%p_top_in_bc, &
                                      metrics%idxCu, metrics%idyCv, nx, ny, nz)
         end if
      else

         ! ---- Pass 1: hydrostatic integration per column ----
         ! Surface boundary condition: p_edge at the top of the
         ! water column is zero (atmospheric absorbed into Boussinesq).
         ! March down: each layer adds rho*g*h_layer to the pressure
         ! at the layer below.  FV_WRIGHT also computes rho_insitu
         ! per layer with a single Picard step.
         if (pgf%variant == OPGF_VARIANT_FV_WRIGHT) then
            if (ms%idx_salinity > 0 .and. ms%idx_temperature > 0) then
               call eos_wright_pgf_column_sweep_impl( &
                  ms%h_layer, &
                  ms%tracers(ms%idx_salinity)%hTr, &
                  ms%tracers(ms%idx_temperature)%hTr, &
                  ms%rho_layer, &
                  ms%p_top, &
                  pgf%p_edge%data, &
                  pgf%rho_insitu%data, &
                  GRAVITY, pgf%rho0, &
                  nx, ny, nz)
            else
               ! No S, T registered: fall through to rho_layer.  Keeps
               ! the test scaffolding (which sets rho_layer directly
               ! without registering tracers) workable.
               do concurrent(j=1:ny, i=1:nx)
                  pgf%p_edge%data(i, j, nz + 1) = 0.0_wp
                  do k = nz, 1, -1
                     pgf%p_edge%data(i, j, k) = pgf%p_edge%data(i, j, k + 1) + &
                                                GRAVITY*ms%rho_layer(i, j, k)*ms%h_layer(i, j, k)
                     pgf%rho_insitu%data(i, j, k) = ms%rho_layer(i, j, k)
                  end do
               end do
            end if
         else
            do concurrent(j=1:ny, i=1:nx)
               pgf%p_edge%data(i, j, nz + 1) = 0.0_wp
               do k = nz, 1, -1
                  pgf%p_edge%data(i, j, k) = pgf%p_edge%data(i, j, k + 1) + &
                                             GRAVITY*ms%rho_layer(i, j, k)*ms%h_layer(i, j, k)
               end do
            end do
         end if

         ! ---- Pass 2: east-face acceleration ----
         ! `idxCu(i,j)` replaces the scalar `inv_dx` (D4) — bit-identical on
         ! uniform Cartesian.
         if (pgf%variant == OPGF_VARIANT_FV_WRIGHT) then
            do concurrent(k=1:nz, j=1:ny, i=2:nx) &
               local(p_centre_left, p_centre_right, rho_face, z_correction)
               p_centre_right = 0.5_wp*(pgf%p_edge%data(i, j, k) + pgf%p_edge%data(i, j, k + 1))
               p_centre_left = 0.5_wp*(pgf%p_edge%data(i - 1, j, k) + pgf%p_edge%data(i - 1, j, k + 1))
               rho_face = 0.5_wp*(pgf%rho_insitu%data(i - 1, j, k) + pgf%rho_insitu%data(i, j, k))
               z_correction = GRAVITY*rho_face* &
                              (pgf%z_centre%data(i, j, k) - pgf%z_centre%data(i - 1, j, k))*metrics%idxCu(i, j)
               pgf%dpdx_face%data(i, j, k) = -inv_rho0*( &
                                             (p_centre_right - p_centre_left)*metrics%idxCu(i, j) + z_correction)
            end do
         else
            ! FV_LITE — the only variant that still reaches this pass
            ! (GPRIME returned, MONT and FV_MOM6 branched above).
            do concurrent(k=1:nz, j=1:ny, i=2:nx) &
               local(p_centre_left, p_centre_right, rho_face, z_correction)
               p_centre_right = 0.5_wp*(pgf%p_edge%data(i, j, k) + pgf%p_edge%data(i, j, k + 1))
               p_centre_left = 0.5_wp*(pgf%p_edge%data(i - 1, j, k) + pgf%p_edge%data(i - 1, j, k + 1))
               rho_face = 0.5_wp*(ms%rho_layer(i - 1, j, k) + ms%rho_layer(i, j, k))
               z_correction = GRAVITY*rho_face* &
                              (pgf%z_centre%data(i, j, k) - pgf%z_centre%data(i - 1, j, k))*metrics%idxCu(i, j)
               pgf%dpdx_face%data(i, j, k) = -inv_rho0*( &
                                             (p_centre_right - p_centre_left)*metrics%idxCu(i, j) + z_correction)
            end do
         end if
         do concurrent(k=1:nz, j=1:ny)
            pgf%dpdx_face%data(1, j, k) = 0.0_wp
            pgf%dpdx_face%data(nx + 1, j, k) = 0.0_wp
         end do

         ! ---- Pass 3: north-face acceleration ----
         ! `idyCv(i,j)` replaces the scalar `inv_dy` (D4).
         if (pgf%variant == OPGF_VARIANT_FV_WRIGHT) then
            do concurrent(k=1:nz, j=2:ny, i=1:nx) &
               local(p_centre_below, p_centre_above, rho_face, z_correction)
               p_centre_above = 0.5_wp*(pgf%p_edge%data(i, j, k) + pgf%p_edge%data(i, j, k + 1))
               p_centre_below = 0.5_wp*(pgf%p_edge%data(i, j - 1, k) + pgf%p_edge%data(i, j - 1, k + 1))
               rho_face = 0.5_wp*(pgf%rho_insitu%data(i, j - 1, k) + pgf%rho_insitu%data(i, j, k))
               z_correction = GRAVITY*rho_face* &
                              (pgf%z_centre%data(i, j, k) - pgf%z_centre%data(i, j - 1, k))*metrics%idyCv(i, j)
               pgf%dpdy_face%data(i, j, k) = -inv_rho0*( &
                                             (p_centre_above - p_centre_below)*metrics%idyCv(i, j) + z_correction)
            end do
         else
            ! FV_LITE — see the Pass-2 comment.
            do concurrent(k=1:nz, j=2:ny, i=1:nx) &
               local(p_centre_below, p_centre_above, rho_face, z_correction)
               p_centre_above = 0.5_wp*(pgf%p_edge%data(i, j, k) + pgf%p_edge%data(i, j, k + 1))
               p_centre_below = 0.5_wp*(pgf%p_edge%data(i, j - 1, k) + pgf%p_edge%data(i, j - 1, k + 1))
               rho_face = 0.5_wp*(ms%rho_layer(i, j - 1, k) + ms%rho_layer(i, j, k))
               z_correction = GRAVITY*rho_face* &
                              (pgf%z_centre%data(i, j, k) - pgf%z_centre%data(i, j - 1, k))*metrics%idyCv(i, j)
               pgf%dpdy_face%data(i, j, k) = -inv_rho0*( &
                                             (p_centre_above - p_centre_below)*metrics%idyCv(i, j) + z_correction)
            end do
         end if
         do concurrent(k=1:nz, i=1:nx)
            pgf%dpdy_face%data(i, 1, k) = 0.0_wp
            pgf%dpdy_face%data(i, ny + 1, k) = 0.0_wp
         end do
      end if

      ! ---- Pass 4: grounded-layer gate (VCOORD_LAGRANGIAN only) ----
      ! Passes 2/3 form a two-point Jacobian: the `Δp_centre` and
      ! `g·ρ_layer·Δz_centre` terms cancel AT REST only while the two abutting
      ! layer centres lie in a common z-interval whose ambient density IS
      ! `ρ_layer`.  Where an isopycnal layer has wedged out against the bed on
      ! one side, the centres are hundreds of metres apart, the interval
      ! between them holds OTHER density classes, and what survives is
      ! `g·(ρ_layer − ρ̄_ambient)·∂z/∂x` — a pressure gradient on a motionless
      ! ocean.  There is no common depth to difference the pressure across
      ! there, so the honest face value is ZERO; the layer is still free to be
      ! re-wetted by continuity's upwind flux and the barotropic correction.
      !
      ! FV_MOM6 is gated by the SAME test.  Its face assembly is not the
      ! two-point Jacobian but the layer-integrated FV-Bouss form, yet the
      ! defect is the geometry, not the quadrature: where the layer occupies
      ! disjoint z-intervals in the two columns the layer-integrated pressure
      ! difference is likewise being taken between depths that share no water
      ! of that density class, and the `1/(h_L + h_R + h_neglect)` divisor
      ! does NOT suppress it — the deep side keeps the thickness up while the
      ! grounded side contributes the whole `e_bot` offset.  Measured on
      ! `seamount_conservative_floor.nml` (form='fv_mom6'): En 2.17e-03 →
      ! 1.39e-26 at day 2.
      !
      ! MONT is gated by the same test for the same reason.  Its face
      ! expression is neither the two-point Jacobian nor the layer-integrated
      ! form, but a grounded layer still puts the two columns' `e_edge` values
      ! hundreds of metres apart, so the `M` recursion stops producing a
      ! horizontally uniform potential at rest and the residual reappears.
      !
      ! A separate guarded pass on purpose: when the gate is off (every vcoord
      ! but LAGRANGIAN) not one extra load is issued ⇒ bit-identical.
      ! Reads `z_centre`, so the driver only ever sets the flag for the four
      ! variants whose allocation gate provides it (MONT / FV_LITE /
      ! FV_WRIGHT, and FV_MOM6 — where `z_centre` is allocated + filled for
      ! this gate alone).
      if (pgf%skip_nonoverlap) then
         do concurrent(k=1:nz, j=1:ny, i=2:nx)
            if (min(pgf%z_centre%data(i - 1, j, k) + 0.5_wp*ms%h_layer(i - 1, j, k), &
                    pgf%z_centre%data(i, j, k) + 0.5_wp*ms%h_layer(i, j, k)) <= &
                max(pgf%z_centre%data(i - 1, j, k) - 0.5_wp*ms%h_layer(i - 1, j, k), &
                    pgf%z_centre%data(i, j, k) - 0.5_wp*ms%h_layer(i, j, k))) then
               pgf%dpdx_face%data(i, j, k) = 0.0_wp
            end if
         end do
         do concurrent(k=1:nz, j=2:ny, i=1:nx)
            if (min(pgf%z_centre%data(i, j - 1, k) + 0.5_wp*ms%h_layer(i, j - 1, k), &
                    pgf%z_centre%data(i, j, k) + 0.5_wp*ms%h_layer(i, j, k)) <= &
                max(pgf%z_centre%data(i, j - 1, k) - 0.5_wp*ms%h_layer(i, j - 1, k), &
                    pgf%z_centre%data(i, j, k) - 0.5_wp*ms%h_layer(i, j, k))) then
               pgf%dpdy_face%data(i, j, k) = 0.0_wp
            end if
         end do
      end if
   end subroutine ocean_pressure_force_compute

   subroutine ocean_pressure_force_apply(pgf, ms, dt, no_wait)
      !! Forward-Euler accumulation of the PGF acceleration onto the face
      !! velocities. Additive (not overwriting), so apply ordering vs the
      !! Coriolis apply doesn't matter before the next tendency-compute.
      !! `no_wait` (optional, default .false. ⇒ blocking): when .true. the
      !! apply loops run on OpenACC queue 1 and return WITHOUT syncing, so
      !! the batched velocity-apply chain can `!$acc wait(1)` ONCE. Not
      !! `pure` (async/wait directives); still functionally pure.
      type(ocean_pressure_force_t), intent(in) :: pgf
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
                                      dt*pgf%dpdx_face%data(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx_vface)
         ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + &
                                      dt*pgf%dpdy_face%data(i, j, k)
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_pressure_force_apply

   pure subroutine compute_gprime_impl(h_layer, b, dpdx_face, dpdy_face, &
                                       gfs, gint, idxCu, idyCv, nx, ny, nz)
      !! Reduced-gravity / gprime PGF for NK = 2.
      !!
      !! Convention: k = 1 bottom (heavier), k = nz = 2 surface (lighter).
      !! `b(i, j)` is bathymetric depth (positive-down), used to recover ∇η:
      !!   sum_h = h_1 + h_2 = b + η;  η = sum_h - b;  ∇η = ∇(h_1+h_2) - ∇b.
      !! The ∇b subtraction matters: without it ∇H_bathy dominates ∇η by
      !! 10³–10⁴ on shelf-break/spoon configs, over-driving the gyre.
      !!
      !!   a_top = -g_FS · ∇η
      !!   a_bot = -g_FS · ∇η - g'_int · ∇h_1
      !!
      !! u-face gradient `(f(i,j) - f(i-1,j)) · idxCu(i,j)`, mirror for
      !! v-face. Higher k stays zero ⇒ no-op on NK > 2.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: gfs, gint
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: b(nx, ny)
      real(wp), intent(inout) :: dpdx_face(nx + 1, ny, nz)
      real(wp), intent(inout) :: dpdy_face(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: sum_h_left, sum_h_right, sum_h_below, sum_h_above
      real(wp) :: db_dx, db_dy
      real(wp) :: dssh_dx, dssh_dy, dhbot_dx, dhbot_dy

      ! Zero everything first.  Then fill k = 1 and k = 2 explicitly.
      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         dpdx_face(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         dpdy_face(i, j, k) = 0.0_wp
      end do

      if (nz < 2) return

      ! u-face (east): two adjacent cells (i-1, j) and (i, j).
      do concurrent(j=1:ny, i=2:nx) &
         local(sum_h_left, sum_h_right, db_dx, dssh_dx, dhbot_dx)
         sum_h_left = h_layer(i - 1, j, 1) + h_layer(i - 1, j, 2)
         sum_h_right = h_layer(i, j, 1) + h_layer(i, j, 2)
         db_dx = (b(i, j) - b(i - 1, j))*idxCu(i, j)
         dssh_dx = (sum_h_right - sum_h_left)*idxCu(i, j) - db_dx
         dhbot_dx = (h_layer(i, j, 1) - h_layer(i - 1, j, 1))*idxCu(i, j)
         dpdx_face(i, j, 2) = -gfs*dssh_dx                 ! top layer (k=nz=2)
         dpdx_face(i, j, 1) = -gfs*dssh_dx - gint*dhbot_dx  ! bottom layer (k=1)
      end do

      ! v-face (north): two adjacent cells (i, j-1) and (i, j).
      do concurrent(j=2:ny, i=1:nx) &
         local(sum_h_below, sum_h_above, db_dy, dssh_dy, dhbot_dy)
         sum_h_below = h_layer(i, j - 1, 1) + h_layer(i, j - 1, 2)
         sum_h_above = h_layer(i, j, 1) + h_layer(i, j, 2)
         db_dy = (b(i, j) - b(i, j - 1))*idyCv(i, j)
         dssh_dy = (sum_h_above - sum_h_below)*idyCv(i, j) - db_dy
         dhbot_dy = (h_layer(i, j, 1) - h_layer(i, j - 1, 1))*idyCv(i, j)
         dpdy_face(i, j, 2) = -gfs*dssh_dy
         dpdy_face(i, j, 1) = -gfs*dssh_dy - gint*dhbot_dy
      end do
   end subroutine compute_gprime_impl

   pure subroutine compute_fv_mom6_impl(h_layer, rho_layer, b, &
                                        e_face, pa, intz_dpa, &
                                        intx_pa, inty_pa, &
                                        intx_dpa, inty_dpa, &
                                        dpdx_face, dpdy_face, &
                                        rho0, rho_ref, h_neglect, &
                                        gfs_scale, mass_weight, &
                                        p_top, p_top_in_bc, &
                                        idxCu, idyCv, nx, ny, nz)
      !! Faithful port of MOM6's `PressureForce_FV_Bouss` per-layer PGF
      !! for the Boussinesq + per-layer Rlay path.
      !!
      !! Pass layout (matches MOM6's `PressureForce_FV_Bouss`):
      !!
      !!   Pass 1 (per column):
      !!     e_face(i, j, k_face) — interface heights, positive up.
      !!       e_face(1) = -b (bed); e_face(k+1) = e_face(k) + h_layer(k);
      !!       e_face(nz+1) = -b + sum(h_layer) = η (free surface).
      !!     pa(i, j, nz+1) = rho_ref · g · η  (surface BC), plus
      !!       `p_top(i, j)` when `p_top_in_bc` (see the theorem below).
      !!     pa(i, j, k) = pa(i, j, k+1) + (rho_layer(k) − rho_ref) · g · h(k)
      !!       (marching down).
      !!     intz_dpa(i, j, k) = 0.5 · (rho_layer(k) − rho_ref) · g · h(k)²
      !!       (mid-point rule).
      !!
      !!   Pass 2 (per face, march down):
      !!     intx_pa(i, j, nz+1) = 0.5 · (pa(i-1, ., nz+1) + pa(i, ., nz+1))
      !!       (surface BC).
      !!     intx_dpa(i, j, k) = 0.5 · g · ((rho(i-1, k) − rho_ref) · h(i-1, k)
      !!                                    + (rho(i, k) − rho_ref) · h(i, k))
      !!       (generalised to per-cell rho).
      !!     intx_pa(i, j, k) = intx_pa(i, j, k+1) + intx_dpa(i, j, k)
      !!       (face pressure recurrence).
      !!     Symmetric on v-face.
      !!
      !!     When `mass_weight = .true.`, at hydrostatically-inconsistent
      !!     unequal-depth faces (`hWght > 0`,
      !!     `hWght = max(0, e_bed_R − e_top_L,k, e_bed_L − e_top_R,k)`)
      !!     the layer density entering `dpa_L`/`dpa_R` is replaced by
      !!     the MOM6 `hWt_LL/LR/RR/RL` blend biased toward the
      !!     thinner column, AND the integral uses the face-interpolated
      !!     thickness `dz = 0.5·(h_L + h_R)` for BOTH samples (MOM6
      !!     `dz_x · rho_anom`).  The interpolated
      !!     thickness is load-bearing: the plain per-cell form
      !!     `0.5·g·(ρ_L'·h_L + ρ_R'·h_R)` is invariant under the
      !!     (Σρh-conserving) hWt blend, so the shelf-break cancellation
      !!     only appears when one `dz` multiplies both blended
      !!     densities.  `hWght = 0` (aligned / equal-depth) ⇒ the exact
      !!     per-cell layer-midpoint average ⇒ bit-identical.
      !!
      !!   Pass 3 (PFu/PFv assembly):
      !!     numer = ((pa(L, k+1) · h(L, k) + intz_dpa(L, k))
      !!              − (pa(R, k+1) · h(R, k) + intz_dpa(R, k)))
      !!           + (h(R, k) − h(L, k)) · intx_pa(face, k+1)
      !!           − (e_face(R, k) − e_face(L, k)) · intx_dpa(face, k)
      !!     denom = h(L, k) + h(R, k) + h_neglect
      !!     PFu(face, k) = numer · (2 · I_Rho0 · IdxCu) / denom
      !!
      !! Convention map (MOM6 → Roundabout):
      !!   MOM6 K (top of layer k_mom6)        → ours k+1 (top of layer k)
      !!   MOM6 K+1 (bottom of layer k_mom6)   → ours k   (bottom of layer k)
      !!   MOM6 i (left of u-face I)           → ours i-1 (west of u-face i)
      !!   MOM6 i+1 (right of u-face I)        → ours i   (east of u-face i)
      !!   MOM6 k_mom6 = 1 (surface layer)     → ours k = nz
      !!   MOM6 k_mom6 = nz_mom6 (bed layer)   → ours k = 1
      !!
      !! Wall faces (face index 1 and N+1) get zero by convention — the
      !! BT-substep / slow continuity already enforces u=0 there.
      !!
      !! ## Theorem — a depth-uniform top load is baroclinically inert here
      !!
      !! Perturb the TOP boundary condition only: `pa(·,nz+1) → pa(·,nz+1)
      !! + δp` with `δp(i,j)` independent of `k`. Every `dpa`, `intz_dpa`
      !! and `intx_dpa` is unchanged, and the Pass-2 recurrence shifts
      !! `intx_pa(K) → intx_pa(K) + ½(δp_L + δp_R)` for EVERY `K`. The
      !! Pass-3 numerator therefore moves by
      !!
      !!   δnumer = δp_L·h_L − δp_R·h_R + (h_R − h_L)·½(δp_L + δp_R)
      !!          = ½(h_L + h_R)·(δp_L − δp_R)
      !!   δPFu(k) = −(1/ρ₀)·(δp_R − δp_L)·IdxCu · (h_L+h_R)/(h_L+h_R+h_neglect)
      !!
      !! — i.e. exactly `−(1/ρ₀)·∂δp/∂x`, **the same in every layer**, up
      !! to the `h_neglect` divisor (a relative `h_n/h_k ≈ 1e-10` for a
      !! metre-thick layer, `≈7e-7` for one at `H_VANISHED`).
      !!
      !! Consequence for the SPLIT solver (`&ocean_bt_nml bc_pgf_forcing`,
      !! default, MOM6 `BT_force`): the depth mean of the layer PGF FORCES
      !! the barotropic substep, so a depth-uniform `δPFu` reaches the
      !! barotropic mode.  The part of `p_top` that is the atmospheric /
      !! anomaly load `sf%p_surf` is ALSO on the `eta_forcing` seam, and
      !! `set_fast_forcing_eta_pf` sheds `g·∇η_ib` from the forcing so it
      !! is counted once; the static ice load `p_ice_ref` cancels inside
      !! `pa(nz+1)` against the datum-shifted `η_geo`, and what survives
      !! of it is the physical reference-density shortfall.  (The legacy
      !! split, `bc_pgf_forcing = .false.`, subtracted the whole depth
      !! mean, so there a depth-uniform `δPFu` cancelled identically.)
      !! The UNSPLIT driver has no seam, so there this term is the load's
      !! only path into the momentum.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: rho_layer(nx, ny, nz)
      real(wp), intent(in)    :: b(nx, ny)
      real(wp), intent(inout) :: e_face(nx, ny, nz + 1)
      real(wp), intent(inout) :: pa(nx, ny, nz + 1)
      real(wp), intent(inout) :: intz_dpa(nx, ny, nz)
      real(wp), intent(inout) :: intx_pa(nx + 1, ny, nz + 1)
      real(wp), intent(inout) :: inty_pa(nx, ny + 1, nz + 1)
      real(wp), intent(inout) :: intx_dpa(nx + 1, ny, nz)
      real(wp), intent(inout) :: inty_dpa(nx, ny + 1, nz)
      real(wp), intent(inout) :: dpdx_face(nx + 1, ny, nz)
      real(wp), intent(inout) :: dpdy_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: rho0, rho_ref, h_neglect, gfs_scale
      logical, intent(in)    :: mass_weight
      real(wp), intent(in)    :: p_top(nx, ny)
         !! Top-of-column pressure (Pa, `>= 0`), `multilayer_state_t%p_top`.
         !! Consulted only when `p_top_in_bc`; the zero array otherwise.
      logical, intent(in)    :: p_top_in_bc
         !! Add `p_top` to the Pass-1 surface BC. `.false.` ⇒ the
         !! assignment is character-for-character the pre-knob one ⇒
         !! bit-identical (same branch-on-a-scalar-knob shape as
         !! `mass_weight` in Pass 2).
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)

      integer  :: i, j, k
      real(wp) :: inv_rho0, eta, rho_anom, dpa_kk
      real(wp) :: dpa_L, dpa_R, h_L, h_R, e_bot_L, e_bot_R
      real(wp) :: pa_h_intz_L, pa_h_intz_R, numer, denom
      real(wp) :: dM_coeff, ddM_dx, ddM_dy
      real(wp) :: hwght, hwl, hwr, idenom_hw, hwt_ll, hwt_lr, hwt_rr, hwt_rl
      real(wp) :: rho_face_l, rho_face_r

      inv_rho0 = 1.0_wp/rho0

      ! ---- Pass 1: per-column build of e_face, pa, intz_dpa ----
      do concurrent(j=1:ny, i=1:nx) local(k, eta, rho_anom, dpa_kk)
         e_face(i, j, 1) = -b(i, j)
         do k = 1, nz
            e_face(i, j, k + 1) = e_face(i, j, k) + h_layer(i, j, k)
         end do
         eta = e_face(i, j, nz + 1)
         if (p_top_in_bc) then
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta + p_top(i, j)
         else
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta
         end if
         do k = nz, 1, -1
            rho_anom = rho_layer(i, j, k) - rho_ref
            dpa_kk = rho_anom*GRAVITY*h_layer(i, j, k)
            pa(i, j, k) = pa(i, j, k + 1) + dpa_kk
            intz_dpa(i, j, k) = 0.5_wp*dpa_kk*h_layer(i, j, k)
         end do
      end do

      ! ---- Pass 2a: u-face horizontal integrals ----
      do concurrent(j=1:ny, i=2:nx) local(k, dpa_L, dpa_R, hwght, hwl, hwr, &
                                          idenom_hw, hwt_ll, hwt_lr, hwt_rr, hwt_rl, &
                                          rho_face_l, rho_face_r)
         intx_pa(i, j, nz + 1) = 0.5_wp*(pa(i - 1, j, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            ! hWght: hydrostatic-inconsistency measure at this u-face for
            ! layer k.  Bed height = e_face(.,1); layer-k top = e_face(.,k+1).
            ! MOM6 form: max(0, bed_R - top_L, bed_L - top_R).
            hwght = 0.0_wp
            if (mass_weight) then
               hwght = max(0.0_wp, &
                           e_face(i, j, 1) - e_face(i - 1, j, k + 1), &
                           e_face(i - 1, j, 1) - e_face(i, j, k + 1))
            end if
            if (hwght > 0.0_wp) then
               ! Hydrostatically-inconsistent face: blend the layer
               ! density toward the thinner column (MOM6 hWt_*) and
               ! integrate with the face-interpolated thickness
               ! dz = 0.5·(h_L + h_R).  The interpolated thickness is
               ! what makes the blend non-trivial — the plain per-cell
               ! form `0.5·g·(ρ_L'·h_L + ρ_R'·h_R)` is invariant under
               ! the (Σρh-conserving) hWt blend, so the cancellation
               ! only appears when the SAME dz multiplies both samples
               ! (MOM6 `dz_x · rho_anom`).
               hwl = h_layer(i - 1, j, k) + h_neglect
               hwr = h_layer(i, j, k) + h_neglect
               hwght = hwght*((hwl - hwr)/(hwl + hwr))**2
               idenom_hw = 1.0_wp/(hwght*(hwr + hwl) + hwl*hwr)
               hwt_ll = (hwght*hwl + hwr*hwl)*idenom_hw
               hwt_lr = (hwght*hwr)*idenom_hw
               hwt_rr = (hwght*hwr + hwr*hwl)*idenom_hw
               hwt_rl = (hwght*hwl)*idenom_hw
               rho_face_l = hwt_ll*rho_layer(i - 1, j, k) + hwt_lr*rho_layer(i, j, k)
               rho_face_r = hwt_rl*rho_layer(i - 1, j, k) + hwt_rr*rho_layer(i, j, k)
               dpa_L = (rho_face_l - rho_ref)*GRAVITY* &
                       (0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k)))
               dpa_R = (rho_face_r - rho_ref)*GRAVITY* &
                       (0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k)))
            else
               ! Aligned / equal-depth face: exact layer-midpoint
               ! average (bit-identical to the pre-knob FV_MOM6 form).
               dpa_L = (rho_layer(i - 1, j, k) - rho_ref)*GRAVITY*h_layer(i - 1, j, k)
               dpa_R = (rho_layer(i, j, k) - rho_ref)*GRAVITY*h_layer(i, j, k)
            end if
            intx_dpa(i, j, k) = 0.5_wp*(dpa_L + dpa_R)
            intx_pa(i, j, k) = intx_pa(i, j, k + 1) + intx_dpa(i, j, k)
         end do
      end do
      do concurrent(k=1:nz, j=1:ny)
         intx_dpa(1, j, k) = 0.0_wp
         intx_dpa(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, j=1:ny)
         intx_pa(1, j, k) = 0.0_wp
         intx_pa(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 2b: v-face horizontal integrals ----
      do concurrent(j=2:ny, i=1:nx) local(k, dpa_L, dpa_R, hwght, hwl, hwr, &
                                          idenom_hw, hwt_ll, hwt_lr, hwt_rr, hwt_rl, &
                                          rho_face_l, rho_face_r)
         inty_pa(i, j, nz + 1) = 0.5_wp*(pa(i, j - 1, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            hwght = 0.0_wp
            if (mass_weight) then
               hwght = max(0.0_wp, &
                           e_face(i, j, 1) - e_face(i, j - 1, k + 1), &
                           e_face(i, j - 1, 1) - e_face(i, j, k + 1))
            end if
            if (hwght > 0.0_wp) then
               ! Hydrostatically-inconsistent face: hWt density blend +
               ! face-interpolated thickness (see Pass 2a comment).
               hwl = h_layer(i, j - 1, k) + h_neglect
               hwr = h_layer(i, j, k) + h_neglect
               hwght = hwght*((hwl - hwr)/(hwl + hwr))**2
               idenom_hw = 1.0_wp/(hwght*(hwr + hwl) + hwl*hwr)
               hwt_ll = (hwght*hwl + hwr*hwl)*idenom_hw
               hwt_lr = (hwght*hwr)*idenom_hw
               hwt_rr = (hwght*hwr + hwr*hwl)*idenom_hw
               hwt_rl = (hwght*hwl)*idenom_hw
               rho_face_l = hwt_ll*rho_layer(i, j - 1, k) + hwt_lr*rho_layer(i, j, k)
               rho_face_r = hwt_rl*rho_layer(i, j - 1, k) + hwt_rr*rho_layer(i, j, k)
               dpa_L = (rho_face_l - rho_ref)*GRAVITY* &
                       (0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k)))
               dpa_R = (rho_face_r - rho_ref)*GRAVITY* &
                       (0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k)))
            else
               ! Aligned / equal-depth face: exact layer-midpoint average.
               dpa_L = (rho_layer(i, j - 1, k) - rho_ref)*GRAVITY*h_layer(i, j - 1, k)
               dpa_R = (rho_layer(i, j, k) - rho_ref)*GRAVITY*h_layer(i, j, k)
            end if
            inty_dpa(i, j, k) = 0.5_wp*(dpa_L + dpa_R)
            inty_pa(i, j, k) = inty_pa(i, j, k + 1) + inty_dpa(i, j, k)
         end do
      end do
      do concurrent(k=1:nz, i=1:nx)
         inty_dpa(i, 1, k) = 0.0_wp
         inty_dpa(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, i=1:nx)
         inty_pa(i, 1, k) = 0.0_wp
         inty_pa(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 3: PFu assembly ----
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i - 1, j, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i - 1, j, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i - 1, j, k + 1)*h_L + intz_dpa(i - 1, j, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*intx_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*intx_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdx_face(i, j, k) = numer*(2.0_wp*inv_rho0*idxCu(i, j))/denom
      end do
      do concurrent(k=1:nz, j=1:ny)
         dpdx_face(1, j, k) = 0.0_wp
         dpdx_face(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 4: PFv assembly ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i, j - 1, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i, j - 1, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i, j - 1, k + 1)*h_L + intz_dpa(i, j - 1, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*inty_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*inty_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdy_face(i, j, k) = numer*(2.0_wp*inv_rho0*idyCv(i, j))/denom
      end do
      do concurrent(k=1:nz, i=1:nx)
         dpdy_face(i, 1, k) = 0.0_wp
         dpdy_face(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 5: Montgomery dM correction (MOM6 GFS_scale) ----
      ! Subtracts `(1 - gfs_scale)·(g/ρ₀)·ρ_surf·∇η` from every layer's
      ! PGF so the slow tendency only carries `gfs_scale · g · ∇η`.
      ! The BT substep evolves η with `g_bt = gfs_scale · GRAVITY`
      ! (set by the driver); combined they reproduce the physical
      ! surface-gravity coupling at the chosen reduced value.
      ! Depth-independent → BT mass-flux invariant intact.  Mirrors
      ! MOM6's Boussinesq non-EOS branch
      ! (`rho_surf = Rlay(top)`).  No-op when `gfs_scale = 1`
      ! (modulo the small tolerance below) so the default
      ! configuration is bit-identical to the pre-knob FV_MOM6 port.
      if (gfs_scale < 1.0_wp - 1.0e-12_wp) then
         dM_coeff = (gfs_scale - 1.0_wp)*GRAVITY*inv_rho0
         do concurrent(k=1:nz, j=1:ny, i=2:nx) local(ddM_dx)
            ddM_dx = dM_coeff*(rho_layer(i, j, nz)*e_face(i, j, nz + 1) &
                               - rho_layer(i - 1, j, nz)*e_face(i - 1, j, nz + 1))*idxCu(i, j)
            dpdx_face(i, j, k) = dpdx_face(i, j, k) - ddM_dx
         end do
         do concurrent(k=1:nz, j=2:ny, i=1:nx) local(ddM_dy)
            ddM_dy = dM_coeff*(rho_layer(i, j, nz)*e_face(i, j, nz + 1) &
                               - rho_layer(i, j - 1, nz)*e_face(i, j - 1, nz + 1))*idyCv(i, j)
            dpdy_face(i, j, k) = dpdy_face(i, j, k) - ddM_dy
         end do
      end if
   end subroutine compute_fv_mom6_impl

   pure subroutine compute_fv_mom6_reconstruct_impl(h_layer, hS, hT, b, eos, &
                                                    S_t, S_b, T_t, T_b, &
                                                    e_face, pa, intz_dpa, &
                                                    intx_pa, inty_pa, &
                                                    intx_dpa, inty_dpa, &
                                                    dpdx_face, dpdy_face, &
                                                    rho0, rho_ref, h_neglect, &
                                                    gfs_scale, recon_scheme, &
                                                    p_top, p_top_in_bc, &
                                                    idxCu, idyCv, nx, ny, nz)
      !! FV_MOM6 pressure-gradient with in-layer T/S reconstruction.
      !!
      !! Same Pass 3-5 face assembly as `compute_fv_mom6_impl`, but the
      !! two integrals that assembly consumes are BOTH taken from the
      !! reconstructed sub-layer T/S profile rather than a layer mean:
      !!
      !!   * Pass 1 (per column) replaces the PCM `dpa(k)` / `intz_dpa(k)`
      !!     with the 5-point VERTICAL Boole quadrature of the monotone
      !!     PLM/PPM profile — the side integrals of the control volume.
      !!   * Pass 2 (per face) replaces the two-column trapezoid
      !!     `0.5*(dpa_L + dpa_R)` with the 5-point HORIZONTAL Boole
      !!     quadrature `boole_dpa_face` — the top/bottom (tilted) edges.
      !!
      !! Both are required for the defining property: with a linear EOS
      !! and T/S linear in z, the PGF then vanishes to round-off for ANY
      !! layer geometry (Adcroft, Hallberg & Harrison 2008; Yung,
      !! Hallberg, Adcroft & Morrison 2026 §2.4).  Correcting the vertical
      !! integral alone leaves the horizontal trapezoid's curvature
      !! residual `g*(-drho/dz)*Delta_e^2/12` at every tilted interface,
      !! which is the sigma "second-kind" pressure-gradient error.
      !!
      !! `mass_weight` (hWght blend) is NOT applied here: it needs a
      !! per-cell density, whereas reconstruction works on column T/S
      !! edges.  Boundary layers take the linear-exact one-sided edge pair
      !! in the edge helper (`boundary_edges_linear`).
      !!
      !! `p_top_in_bc` injects the top load into the SAME Pass-1 surface
      !! BC as the PCM twin, and the Theorem in `compute_fv_mom6_impl`
      !! carries over verbatim: the reconstruction only changes `dpa` /
      !! `intz_dpa`, never the `pa(nz+1)` seed or the `intx_pa`
      !! recurrence, so a depth-uniform `p_top` still perturbs every
      !! layer's `PFu` by the same `−(1/ρ₀)∇p_top`. NOTE that this branch
      !! builds its OWN in-layer EOS pressure inside
      !! `boole_dpa_intz_layer` (`p = −g·ρ₀·z` from the surface-relative
      !! interface height) and that one is NOT offset by `p_top` — which
      !! is exactly why `validate_config` refuses
      !! `&ocean_psurf_nml in_eos` together with
      !! `reconstruct_for_pressure`. The BC injection here is a PRESSURE
      !! boundary condition, not an EOS argument; the two are independent
      !! seams.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: hS(nx, ny, nz)
         !! Salinity * thickness (PSU*m) — layer-mean S = hS / h.
      real(wp), intent(in)    :: hT(nx, ny, nz)
         !! Temperature * thickness (degC*m) — layer-mean T = hT / h.
      real(wp), intent(in)    :: b(nx, ny)
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: S_t(nx, ny, nz), S_b(nx, ny, nz)
      real(wp), intent(inout) :: T_t(nx, ny, nz), T_b(nx, ny, nz)
      real(wp), intent(inout) :: e_face(nx, ny, nz + 1)
      real(wp), intent(inout) :: pa(nx, ny, nz + 1)
      real(wp), intent(inout) :: intz_dpa(nx, ny, nz)
      real(wp), intent(inout) :: intx_pa(nx + 1, ny, nz + 1)
      real(wp), intent(inout) :: inty_pa(nx, ny + 1, nz + 1)
      real(wp), intent(inout) :: intx_dpa(nx + 1, ny, nz)
      real(wp), intent(inout) :: inty_dpa(nx, ny + 1, nz)
      real(wp), intent(inout) :: dpdx_face(nx + 1, ny, nz)
      real(wp), intent(inout) :: dpdy_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: rho0, rho_ref, h_neglect, gfs_scale
      integer, intent(in)    :: recon_scheme
      real(wp), intent(in)    :: p_top(nx, ny)
         !! Top-of-column pressure (Pa, `>= 0`), `multilayer_state_t%p_top`.
      logical, intent(in)    :: p_top_in_bc
         !! Add `p_top` to the Pass-1 surface BC (`.false.` ⇒ bit-identical).
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)

      integer  :: i, j, k
      real(wp) :: inv_rho0, eta, dpa_kk, intz_kk
      real(wp) :: h_L, h_R, e_bot_L, e_bot_R
      real(wp) :: t_m_L, t_m_R, s_m_L, s_m_R
      real(wp) :: pa_h_intz_L, pa_h_intz_R, numer, denom
      real(wp) :: dM_coeff, ddM_dx, ddM_dy
      logical  :: parabolic
      ! Per-column edge-build stacks (fixed-size for local()).
      real(wp) :: h_col(NZ_STACK_MAX), s_col(NZ_STACK_MAX), t_col(NZ_STACK_MAX)
      real(wp) :: st_col(NZ_STACK_MAX), sb_col(NZ_STACK_MAX)
      real(wp) :: tt_col(NZ_STACK_MAX), tb_col(NZ_STACK_MAX)
      ! Gate the layer-mean T/S recovery at H_VANISHED (the D4 vanished-layer
      ! role), not the old 1e-10: during an active drain the PPM limiter only
      ! guarantees h >= 0, so a layer in (0, H_VANISHED] would otherwise pass
      ! `h > 1e-10` and feed hS/h ≈ hS/1e-8 into the reconstruction.  Bit-
      ! identical for any layer with h > H_VANISHED (the `hS/h` branch is
      ! selected either way); only extreme-thin (0, H_VANISHED] layers switch
      ! to the floored fallback.  Reconstruct-for-pressure is default-off so
      ! no shipped anchor exercises this path.
      real(wp), parameter :: H_FLOOR = H_VANISHED

      inv_rho0 = 1.0_wp/rho0
      parabolic = (recon_scheme == PGF_RECON_PPM)

      ! ---- Pass 0: per-column PLM/PPM T/S edge values ----
      ! Build the layer-mean T,S column (= hTr/h, floored), call the
      ! edge helper, write the four edge slots.  Done as its own DC pass
      ! so the quadrature pass below reads clean edge stacks.
      do concurrent(j=1:ny, i=1:nx) local(k, h_col, s_col, t_col, &
                                          st_col, sb_col, tt_col, tb_col)
         do k = 1, nz
            h_col(k) = h_layer(i, j, k)
            if (h_col(k) > H_FLOOR) then
               s_col(k) = hS(i, j, k)/h_col(k)
               t_col(k) = hT(i, j, k)/h_col(k)
            else
               s_col(k) = hS(i, j, k)/H_FLOOR
               t_col(k) = hT(i, j, k)/H_FLOOR
            end if
         end do
         if (parabolic) then
            call ppm_edges_column(nz, h_col, s_col, st_col, sb_col)
            call ppm_edges_column(nz, h_col, t_col, tt_col, tb_col)
         else
            call plm_edges_column(nz, h_col, s_col, st_col, sb_col)
            call plm_edges_column(nz, h_col, t_col, tt_col, tb_col)
         end if
         do k = 1, nz
            S_t(i, j, k) = st_col(k)
            S_b(i, j, k) = sb_col(k)
            T_t(i, j, k) = tt_col(k)
            T_b(i, j, k) = tb_col(k)
         end do
      end do

      ! ---- Pass 1: per-column e_face, pa, intz_dpa via Boole quadrature ----
      ! e_top = e_face(k+1) is the shallower interface of layer k.  The
      ! reconstructed dpa(k) marches the pa stack; intz_dpa(k) is the
      ! first-moment piece.  Both replace the PCM forms.
      do concurrent(j=1:ny, i=1:nx) &
         local(k, eta, dpa_kk, intz_kk, s_col, t_col)
         e_face(i, j, 1) = -b(i, j)
         do k = 1, nz
            e_face(i, j, k + 1) = e_face(i, j, k) + h_layer(i, j, k)
            ! Layer mean (for the PPM parabolic curvature term).
            if (h_layer(i, j, k) > H_FLOOR) then
               s_col(k) = hS(i, j, k)/h_layer(i, j, k)
               t_col(k) = hT(i, j, k)/h_layer(i, j, k)
            else
               s_col(k) = hS(i, j, k)/H_FLOOR
               t_col(k) = hT(i, j, k)/H_FLOOR
            end if
         end do
         eta = e_face(i, j, nz + 1)
         if (p_top_in_bc) then
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta + p_top(i, j)
         else
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta
         end if
         do k = nz, 1, -1
            call boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                      e_face(i, j, k + 1), h_layer(i, j, k), &
                                      T_t(i, j, k), T_b(i, j, k), t_col(k), &
                                      S_t(i, j, k), S_b(i, j, k), s_col(k), &
                                      parabolic, dpa_kk, intz_kk)
            pa(i, j, k) = pa(i, j, k + 1) + dpa_kk
            intz_dpa(i, j, k) = intz_kk
         end do
      end do

      ! ---- Pass 2a: u-face horizontal integrals ----
      ! The along-face mean of the layer pressure increment, by the 5-point
      ! cross-face Boole quadrature of `boole_dpa_face` (sub-columns at the
      ! INTERPOLATED interface height with interpolated T/S).  The
      ! two-column trapezoid `0.5*(dpa_L + dpa_R)` this replaces is exact
      ! only for a pressure linear in x along the edge; under a tilted
      ! interface it leaves the sigma second-kind curvature residual at
      ! every interface — see the `boole_dpa_face` docstring.
      do concurrent(j=1:ny, i=2:nx) local(k, dpa_kk, t_m_L, t_m_R, s_m_L, s_m_R)
         intx_pa(i, j, nz + 1) = 0.5_wp*(pa(i - 1, j, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            t_m_L = recon_layer_mean(hT(i - 1, j, k), h_layer(i - 1, j, k))
            t_m_R = recon_layer_mean(hT(i, j, k), h_layer(i, j, k))
            s_m_L = recon_layer_mean(hS(i - 1, j, k), h_layer(i - 1, j, k))
            s_m_R = recon_layer_mean(hS(i, j, k), h_layer(i, j, k))
            call boole_dpa_face(eos, rho0, rho_ref, &
                                e_face(i - 1, j, k + 1), e_face(i, j, k + 1), &
                                h_layer(i - 1, j, k), h_layer(i, j, k), &
                                T_t(i - 1, j, k), T_b(i - 1, j, k), t_m_L, &
                                T_t(i, j, k), T_b(i, j, k), t_m_R, &
                                S_t(i - 1, j, k), S_b(i - 1, j, k), s_m_L, &
                                S_t(i, j, k), S_b(i, j, k), s_m_R, &
                                parabolic, dpa_kk)
            intx_dpa(i, j, k) = dpa_kk
            intx_pa(i, j, k) = intx_pa(i, j, k + 1) + dpa_kk
         end do
      end do
      do concurrent(k=1:nz, j=1:ny)
         intx_dpa(1, j, k) = 0.0_wp
         intx_dpa(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, j=1:ny)
         intx_pa(1, j, k) = 0.0_wp
         intx_pa(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 2b: v-face horizontal integrals ----
      do concurrent(j=2:ny, i=1:nx) local(k, dpa_kk, t_m_L, t_m_R, s_m_L, s_m_R)
         inty_pa(i, j, nz + 1) = 0.5_wp*(pa(i, j - 1, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            t_m_L = recon_layer_mean(hT(i, j - 1, k), h_layer(i, j - 1, k))
            t_m_R = recon_layer_mean(hT(i, j, k), h_layer(i, j, k))
            s_m_L = recon_layer_mean(hS(i, j - 1, k), h_layer(i, j - 1, k))
            s_m_R = recon_layer_mean(hS(i, j, k), h_layer(i, j, k))
            call boole_dpa_face(eos, rho0, rho_ref, &
                                e_face(i, j - 1, k + 1), e_face(i, j, k + 1), &
                                h_layer(i, j - 1, k), h_layer(i, j, k), &
                                T_t(i, j - 1, k), T_b(i, j - 1, k), t_m_L, &
                                T_t(i, j, k), T_b(i, j, k), t_m_R, &
                                S_t(i, j - 1, k), S_b(i, j - 1, k), s_m_L, &
                                S_t(i, j, k), S_b(i, j, k), s_m_R, &
                                parabolic, dpa_kk)
            inty_dpa(i, j, k) = dpa_kk
            inty_pa(i, j, k) = inty_pa(i, j, k + 1) + dpa_kk
         end do
      end do
      do concurrent(k=1:nz, i=1:nx)
         inty_dpa(i, 1, k) = 0.0_wp
         inty_dpa(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, i=1:nx)
         inty_pa(i, 1, k) = 0.0_wp
         inty_pa(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 3: PFu assembly (identical to compute_fv_mom6_impl) ----
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i - 1, j, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i - 1, j, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i - 1, j, k + 1)*h_L + intz_dpa(i - 1, j, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*intx_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*intx_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdx_face(i, j, k) = numer*(2.0_wp*inv_rho0*idxCu(i, j))/denom
      end do
      do concurrent(k=1:nz, j=1:ny)
         dpdx_face(1, j, k) = 0.0_wp
         dpdx_face(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 4: PFv assembly ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i, j - 1, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i, j - 1, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i, j - 1, k + 1)*h_L + intz_dpa(i, j - 1, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*inty_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*inty_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdy_face(i, j, k) = numer*(2.0_wp*inv_rho0*idyCv(i, j))/denom
      end do
      do concurrent(k=1:nz, i=1:nx)
         dpdy_face(i, 1, k) = 0.0_wp
         dpdy_face(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 5: Montgomery dM correction (MOM6 GFS_scale) ----
      ! rho_surf for the dM term is the reconstructed top-edge density of
      ! the surface layer; here we reuse the layer-mean surface density
      ! recovered from dpa(nz) = pa(nz) - pa(nz+1) divided by g*h, which
      ! equals (rho_surf - rho_ref).  Keep the same depth-independent form.
      if (gfs_scale < 1.0_wp - 1.0e-12_wp) then
         dM_coeff = (gfs_scale - 1.0_wp)*GRAVITY*inv_rho0
         do concurrent(k=1:nz, j=1:ny, i=2:nx) local(ddM_dx)
            ddM_dx = dM_coeff*(recon_rho_surf(pa(i, j, nz), pa(i, j, nz + 1), &
                                              h_layer(i, j, nz), rho_ref) &
                               *e_face(i, j, nz + 1) &
                               - recon_rho_surf(pa(i - 1, j, nz), pa(i - 1, j, nz + 1), &
                                                h_layer(i - 1, j, nz), rho_ref) &
                               *e_face(i - 1, j, nz + 1))*idxCu(i, j)
            dpdx_face(i, j, k) = dpdx_face(i, j, k) - ddM_dx
         end do
         do concurrent(k=1:nz, j=2:ny, i=1:nx) local(ddM_dy)
            ddM_dy = dM_coeff*(recon_rho_surf(pa(i, j, nz), pa(i, j, nz + 1), &
                                              h_layer(i, j, nz), rho_ref) &
                               *e_face(i, j, nz + 1) &
                               - recon_rho_surf(pa(i, j - 1, nz), pa(i, j - 1, nz + 1), &
                                                h_layer(i, j - 1, nz), rho_ref) &
                               *e_face(i, j - 1, nz + 1))*idyCv(i, j)
            dpdy_face(i, j, k) = dpdy_face(i, j, k) - ddM_dy
         end do
      end if
   end subroutine compute_fv_mom6_reconstruct_impl

   pure function use_insitu_pcm(pgf, ms, eos) result(yes)
      !! Does the FV_MOM6 constant-by-layer branch take the IN-SITU
      !! density path (`compute_fv_mom6_insitu_pcm_impl`)?  Only when it
      !! can change the answer: the knob is on, there is an EOS handle and
      !! T/S to evaluate it on, and the EOS depends on pressure.  For the
      !! linear EOS `ms%rho_layer` already IS the in-situ density, so the
      !! legacy path runs, bit-identical.
      type(ocean_pressure_force_t), intent(in) :: pgf
      type(multilayer_state_t), intent(in) :: ms
      type(eos_t), intent(in), optional :: eos
      logical :: yes
      yes = .false.
      if (.not. pgf%insitu_density) return
      if (.not. present(eos)) return
      if (ms%idx_salinity <= 0 .or. ms%idx_temperature <= 0) return
      yes = eos%variant == EOS_VARIANT_WRIGHT_97 .or. &
            eos%variant == EOS_VARIANT_ROQUET_SPV
   end function use_insitu_pcm

   pure subroutine compute_fv_mom6_insitu_pcm_impl(h_layer, hS, hT, b, eos, &
                                                   e_face, pa, intz_dpa, &
                                                   intx_pa, inty_pa, &
                                                   intx_dpa, inty_dpa, &
                                                   dpdx_face, dpdy_face, &
                                                   rho0, rho_ref, h_neglect, &
                                                   gfs_scale, mass_weight, &
                                                   p_top, p_top_in_bc, &
                                                   idxCu, idyCv, nx, ny, nz)
      !! FV_MOM6 pressure gradient, constant-by-layer (PCM) T/S, density at
      !! the IN-SITU pressure — MOM6 `PressureForce_FV_Bouss` with
      !! `RECONSTRUCT_FOR_PRESSURE = False` (`int_density_dz_generic_pcm`).
      !!
      !! The PCM twin `compute_fv_mom6_impl` integrates `ms%rho_layer`, a
      !! POTENTIAL density at the one horizontally uniform `p_ref`.  Its
      !! horizontal difference at depth is then the difference at the
      !! REFERENCE pressure, not at the local one: the thermal expansion
      !! coefficient roughly doubles between the surface and 4000 dbar
      !! (thermobaricity), so with `p_ref = 0` the deep baroclinic
      !! pressure gradient — the bottom-pressure gradient that forces the
      !! barotropic mode over topography — is systematically too weak.
      !! On the global 1-degree WOA13 spin-up it held Drake Passage at
      !! ~80 Sv where MOM6 on the same protocol adjusts to ~155 Sv, and
      !! the transport tracked `p_ref` (0 / 2000 / 4000 dbar: 83 / 143 /
      !! 203 Sv) — the tell of a reference-pressure artefact.
      !!
      !! Here every density is `EOS(T, S, p = −g·rho0·z)` at the point it
      !! is used, integrated by the same 5-point Boole rules as the
      !! reconstruction branch, with the sub-layer profile flat:
      !!
      !!   * Pass 1 (per column): `dpa(k)`, `intz_dpa(k)` from
      !!     `boole_dpa_intz_layer` with top = bottom = mean T/S.
      !!   * Pass 2 (per face): `intx_dpa` / `inty_dpa` from
      !!     `boole_dpa_face_pcm` — end points are the columns' own `dpa`,
      !!     the three interior sub-columns interpolate `z` linearly and
      !!     T/S with MOM6's near-bottom mass weighting (`hWght`, the same
      !!     measure and blend as `compute_fv_mom6_impl`) when
      !!     `mass_weight`.
      !!   * Passes 3–5: the face assembly, identical to the other two
      !!     FV_MOM6 branches.
      !!
      !! The trapezoid `0.5·(dpa_L + dpa_R)` of the potential-density twin
      !! is NOT kept: an in-situ density carries the compressibility
      !! gradient (`~4.4e-3 kg m⁻⁴`), and the trapezoid's curvature
      !! residual `g·(∂ρ/∂z)·Δe²/12` at a tilted interface (a partial-cell
      !! bed step) would be of the size of the signal.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: hS(nx, ny, nz)
         !! Salinity * thickness (PSU*m) — layer-mean S = hS / h.
      real(wp), intent(in)    :: hT(nx, ny, nz)
         !! Temperature * thickness (degC*m) — layer-mean T = hT / h.
      real(wp), intent(in)    :: b(nx, ny)
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: e_face(nx, ny, nz + 1)
      real(wp), intent(inout) :: pa(nx, ny, nz + 1)
      real(wp), intent(inout) :: intz_dpa(nx, ny, nz)
      real(wp), intent(inout) :: intx_pa(nx + 1, ny, nz + 1)
      real(wp), intent(inout) :: inty_pa(nx, ny + 1, nz + 1)
      real(wp), intent(inout) :: intx_dpa(nx + 1, ny, nz)
      real(wp), intent(inout) :: inty_dpa(nx, ny + 1, nz)
      real(wp), intent(inout) :: dpdx_face(nx + 1, ny, nz)
      real(wp), intent(inout) :: dpdy_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: rho0, rho_ref, h_neglect, gfs_scale
      logical, intent(in)    :: mass_weight
         !! MOM6 `MASS_WEIGHT_IN_PRESSURE_GRADIENT` (near-bottom `hWght`).
      real(wp), intent(in)    :: p_top(nx, ny)
         !! Top-of-column pressure (Pa, `>= 0`), `multilayer_state_t%p_top`.
      logical, intent(in)    :: p_top_in_bc
         !! Add `p_top` to the Pass-1 surface BC (`.false.` ⇒ the plain
         !! `rho_ref·g·eta` seed).
      real(wp), intent(in)    :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)

      integer  :: i, j, k
      real(wp) :: inv_rho0, eta, dpa_kk, intz_kk, t_m, s_m
      real(wp) :: t_m_L, t_m_R, s_m_L, s_m_R, dpa_L, dpa_R
      real(wp) :: hwght, hwl, hwr, idenom_hw, hwt_ll, hwt_lr, hwt_rr, hwt_rl
      real(wp) :: h_L, h_R, e_bot_L, e_bot_R
      real(wp) :: pa_h_intz_L, pa_h_intz_R, numer, denom
      real(wp) :: dM_coeff, ddM_dx, ddM_dy

      inv_rho0 = 1.0_wp/rho0

      ! ---- Pass 1: per-column e_face, pa, intz_dpa (in-situ Boole) ----
      do concurrent(j=1:ny, i=1:nx) local(k, eta, dpa_kk, intz_kk, t_m, s_m)
         e_face(i, j, 1) = -b(i, j)
         do k = 1, nz
            e_face(i, j, k + 1) = e_face(i, j, k) + h_layer(i, j, k)
         end do
         eta = e_face(i, j, nz + 1)
         if (p_top_in_bc) then
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta + p_top(i, j)
         else
            pa(i, j, nz + 1) = rho_ref*GRAVITY*eta
         end if
         do k = nz, 1, -1
            t_m = recon_layer_mean(hT(i, j, k), h_layer(i, j, k))
            s_m = recon_layer_mean(hS(i, j, k), h_layer(i, j, k))
            call boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                      e_face(i, j, k + 1), h_layer(i, j, k), &
                                      t_m, t_m, t_m, s_m, s_m, s_m, &
                                      .false., dpa_kk, intz_kk)
            pa(i, j, k) = pa(i, j, k + 1) + dpa_kk
            intz_dpa(i, j, k) = intz_kk
         end do
      end do

      ! ---- Pass 2a: u-face horizontal integrals ----
      do concurrent(j=1:ny, i=2:nx) local(k, dpa_kk, t_m_L, t_m_R, s_m_L, s_m_R, &
                                          dpa_L, dpa_R, hwght, hwl, hwr, idenom_hw, &
                                          hwt_ll, hwt_lr, hwt_rr, hwt_rl)
         intx_pa(i, j, nz + 1) = 0.5_wp*(pa(i - 1, j, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            hwght = 0.0_wp
            if (mass_weight) then
               hwght = max(0.0_wp, &
                           e_face(i, j, 1) - e_face(i - 1, j, k + 1), &
                           e_face(i - 1, j, 1) - e_face(i, j, k + 1))
            end if
            hwt_ll = 1.0_wp
            hwt_lr = 0.0_wp
            hwt_rr = 1.0_wp
            hwt_rl = 0.0_wp
            if (hwght > 0.0_wp) then
               hwl = h_layer(i - 1, j, k) + h_neglect
               hwr = h_layer(i, j, k) + h_neglect
               hwght = hwght*((hwl - hwr)/(hwl + hwr))**2
               idenom_hw = 1.0_wp/(hwght*(hwr + hwl) + hwl*hwr)
               hwt_ll = (hwght*hwl + hwr*hwl)*idenom_hw
               hwt_lr = (hwght*hwr)*idenom_hw
               hwt_rr = (hwght*hwr + hwr*hwl)*idenom_hw
               hwt_rl = (hwght*hwl)*idenom_hw
            end if
            t_m_L = recon_layer_mean(hT(i - 1, j, k), h_layer(i - 1, j, k))
            t_m_R = recon_layer_mean(hT(i, j, k), h_layer(i, j, k))
            s_m_L = recon_layer_mean(hS(i - 1, j, k), h_layer(i - 1, j, k))
            s_m_R = recon_layer_mean(hS(i, j, k), h_layer(i, j, k))
            dpa_L = pa(i - 1, j, k) - pa(i - 1, j, k + 1)
            dpa_R = pa(i, j, k) - pa(i, j, k + 1)
            call boole_dpa_face_pcm(eos, rho0, rho_ref, &
                                    e_face(i - 1, j, k + 1), e_face(i, j, k + 1), &
                                    h_layer(i - 1, j, k), h_layer(i, j, k), &
                                    t_m_L, t_m_R, s_m_L, s_m_R, dpa_L, dpa_R, &
                                    hwt_ll, hwt_lr, hwt_rr, hwt_rl, dpa_kk)
            intx_dpa(i, j, k) = dpa_kk
            intx_pa(i, j, k) = intx_pa(i, j, k + 1) + dpa_kk
         end do
      end do
      do concurrent(k=1:nz, j=1:ny)
         intx_dpa(1, j, k) = 0.0_wp
         intx_dpa(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, j=1:ny)
         intx_pa(1, j, k) = 0.0_wp
         intx_pa(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 2b: v-face horizontal integrals ----
      do concurrent(j=2:ny, i=1:nx) local(k, dpa_kk, t_m_L, t_m_R, s_m_L, s_m_R, &
                                          dpa_L, dpa_R, hwght, hwl, hwr, idenom_hw, &
                                          hwt_ll, hwt_lr, hwt_rr, hwt_rl)
         inty_pa(i, j, nz + 1) = 0.5_wp*(pa(i, j - 1, nz + 1) + pa(i, j, nz + 1))
         do k = nz, 1, -1
            hwght = 0.0_wp
            if (mass_weight) then
               hwght = max(0.0_wp, &
                           e_face(i, j, 1) - e_face(i, j - 1, k + 1), &
                           e_face(i, j - 1, 1) - e_face(i, j, k + 1))
            end if
            hwt_ll = 1.0_wp
            hwt_lr = 0.0_wp
            hwt_rr = 1.0_wp
            hwt_rl = 0.0_wp
            if (hwght > 0.0_wp) then
               hwl = h_layer(i, j - 1, k) + h_neglect
               hwr = h_layer(i, j, k) + h_neglect
               hwght = hwght*((hwl - hwr)/(hwl + hwr))**2
               idenom_hw = 1.0_wp/(hwght*(hwr + hwl) + hwl*hwr)
               hwt_ll = (hwght*hwl + hwr*hwl)*idenom_hw
               hwt_lr = (hwght*hwr)*idenom_hw
               hwt_rr = (hwght*hwr + hwr*hwl)*idenom_hw
               hwt_rl = (hwght*hwl)*idenom_hw
            end if
            t_m_L = recon_layer_mean(hT(i, j - 1, k), h_layer(i, j - 1, k))
            t_m_R = recon_layer_mean(hT(i, j, k), h_layer(i, j, k))
            s_m_L = recon_layer_mean(hS(i, j - 1, k), h_layer(i, j - 1, k))
            s_m_R = recon_layer_mean(hS(i, j, k), h_layer(i, j, k))
            dpa_L = pa(i, j - 1, k) - pa(i, j - 1, k + 1)
            dpa_R = pa(i, j, k) - pa(i, j, k + 1)
            call boole_dpa_face_pcm(eos, rho0, rho_ref, &
                                    e_face(i, j - 1, k + 1), e_face(i, j, k + 1), &
                                    h_layer(i, j - 1, k), h_layer(i, j, k), &
                                    t_m_L, t_m_R, s_m_L, s_m_R, dpa_L, dpa_R, &
                                    hwt_ll, hwt_lr, hwt_rr, hwt_rl, dpa_kk)
            inty_dpa(i, j, k) = dpa_kk
            inty_pa(i, j, k) = inty_pa(i, j, k + 1) + dpa_kk
         end do
      end do
      do concurrent(k=1:nz, i=1:nx)
         inty_dpa(i, 1, k) = 0.0_wp
         inty_dpa(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, i=1:nx)
         inty_pa(i, 1, k) = 0.0_wp
         inty_pa(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 3: PFu assembly (identical to compute_fv_mom6_impl) ----
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i - 1, j, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i - 1, j, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i - 1, j, k + 1)*h_L + intz_dpa(i - 1, j, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*intx_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*intx_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdx_face(i, j, k) = numer*(2.0_wp*inv_rho0*idxCu(i, j))/denom
      end do
      do concurrent(k=1:nz, j=1:ny)
         dpdx_face(1, j, k) = 0.0_wp
         dpdx_face(nx + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 4: PFv assembly ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(h_L, h_R, e_bot_L, e_bot_R, pa_h_intz_L, pa_h_intz_R, numer, denom)
         h_L = h_layer(i, j - 1, k)
         h_R = h_layer(i, j, k)
         e_bot_L = e_face(i, j - 1, k)
         e_bot_R = e_face(i, j, k)
         pa_h_intz_L = pa(i, j - 1, k + 1)*h_L + intz_dpa(i, j - 1, k)
         pa_h_intz_R = pa(i, j, k + 1)*h_R + intz_dpa(i, j, k)
         numer = (pa_h_intz_L - pa_h_intz_R) &
                 + (h_R - h_L)*inty_pa(i, j, k + 1) &
                 - (e_bot_R - e_bot_L)*inty_dpa(i, j, k)
         denom = h_L + h_R + h_neglect
         dpdy_face(i, j, k) = numer*(2.0_wp*inv_rho0*idyCv(i, j))/denom
      end do
      do concurrent(k=1:nz, i=1:nx)
         dpdy_face(i, 1, k) = 0.0_wp
         dpdy_face(i, ny + 1, k) = 0.0_wp
      end do

      ! ---- Pass 5: Montgomery dM correction (MOM6 GFS_scale) ----
      ! Same depth-independent form as the reconstruction branch, with the
      ! surface layer's mean in-situ density recovered from its `dpa`.
      if (gfs_scale < 1.0_wp - 1.0e-12_wp) then
         dM_coeff = (gfs_scale - 1.0_wp)*GRAVITY*inv_rho0
         do concurrent(k=1:nz, j=1:ny, i=2:nx) local(ddM_dx)
            ddM_dx = dM_coeff*(recon_rho_surf(pa(i, j, nz), pa(i, j, nz + 1), &
                                              h_layer(i, j, nz), rho_ref) &
                               *e_face(i, j, nz + 1) &
                               - recon_rho_surf(pa(i - 1, j, nz), pa(i - 1, j, nz + 1), &
                                                h_layer(i - 1, j, nz), rho_ref) &
                               *e_face(i - 1, j, nz + 1))*idxCu(i, j)
            dpdx_face(i, j, k) = dpdx_face(i, j, k) - ddM_dx
         end do
         do concurrent(k=1:nz, j=2:ny, i=1:nx) local(ddM_dy)
            ddM_dy = dM_coeff*(recon_rho_surf(pa(i, j, nz), pa(i, j, nz + 1), &
                                              h_layer(i, j, nz), rho_ref) &
                               *e_face(i, j, nz + 1) &
                               - recon_rho_surf(pa(i, j - 1, nz), pa(i, j - 1, nz + 1), &
                                                h_layer(i, j - 1, nz), rho_ref) &
                               *e_face(i, j - 1, nz + 1))*idyCv(i, j)
            dpdy_face(i, j, k) = dpdy_face(i, j, k) - ddM_dy
         end do
      end if
   end subroutine compute_fv_mom6_insitu_pcm_impl

   pure function recon_layer_mean(hq, h) result(q)
      !$acc routine seq
      !! Layer-mean tracer from the thickness-weighted prognostic,
      !! `q = hq/h`, with the D4 vanished-layer floor.  Same gate as the
      !! reconstruct kernel's Pass 0/1 (`H_VANISHED`, not `1e-10`): during
      !! an active drain the PPM positivity limiter only guarantees
      !! `h >= 0`, so a layer in `(0, H_VANISHED]` would otherwise divide
      !! by a near-zero thickness.
      real(wp), intent(in) :: hq
         !! Thickness-weighted tracer (e.g. `S*h`).
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      real(wp) :: q
      if (h > H_VANISHED) then
         q = hq/h
      else
         q = hq/H_VANISHED
      end if
   end function recon_layer_mean

   pure function recon_rho_surf(pa_k, pa_kp1, h_surf, rho_ref) result(rho_surf)
      !$acc routine seq
      !! Recover the layer-mean surface density from the reconstructed
      !! pressure-anomaly stack: dpa(nz) = pa(nz) - pa(nz+1) =
      !! (rho_surf - rho_ref)*g*h_surf, so rho_surf = rho_ref + dpa/(g*h).
      !! Used only by the gfs_scale Montgomery correction (Pass 5).
      real(wp), intent(in) :: pa_k, pa_kp1, h_surf, rho_ref
      real(wp) :: rho_surf
      real(wp), parameter :: H_FLOOR = 1.0e-10_wp
      real(wp) :: hh
      hh = max(h_surf, H_FLOOR)
      rho_surf = rho_ref + (pa_k - pa_kp1)/(GRAVITY*hh)
   end function recon_rho_surf

   pure function parse_opgf_variant(name) result(code)
      !! Translate a namelist string into an `OPGF_VARIANT_*` code.
      !! Unrecognised values fall back to `OPGF_VARIANT_FV_LITE` (the
      !! production default).
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("mont", "MONT", "montgomery")
         code = OPGF_VARIANT_MONT
      case ("fv_lite", "FV_LITE", "fv-lite", "")
         code = OPGF_VARIANT_FV_LITE
      case ("fv_wright", "FV_WRIGHT", "fv-wright", "wright")
         code = OPGF_VARIANT_FV_WRIGHT
      case ("gprime", "GPRIME", "reduced_gravity")
         code = OPGF_VARIANT_GPRIME
      case ("fv_mom6", "FV_MOM6", "fv-mom6", "mom6_fv", "MOM6_FV")
         code = OPGF_VARIANT_FV_MOM6
      case default
         code = OPGF_VARIANT_FV_LITE
      end select
   end function parse_opgf_variant

   pure function ocean_pressure_force_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the pressure force slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_pressure_force_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%b) &
               + this%p_edge%bytes() &
               + this%z_centre%bytes() &
               + this%mont_M%bytes() &
               + this%rho_insitu%bytes() &
               + this%dpdx_face%bytes() &
               + this%dpdy_face%bytes() &
               + this%e_face%bytes() &
               + this%pa%bytes() &
               + this%intz_dpa%bytes() &
               + this%intx_pa%bytes() &
               + this%inty_pa%bytes() &
               + this%intx_dpa%bytes() &
               + this%inty_dpa%bytes() &
               + this%recon_T_t%bytes() &
               + this%recon_T_b%bytes() &
               + this%recon_S_t%bytes() &
               + this%recon_S_b%bytes()
   end function ocean_pressure_force_bytes

end module rdb_ocean_pressure_force
