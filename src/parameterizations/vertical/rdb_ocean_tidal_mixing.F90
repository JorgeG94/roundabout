!! Internal-tide-driven interior diapycnal mixing (St-Laurent/Simmons) for the ocean core.
module rdb_ocean_tidal_mixing
   !! Bottom-intensified internal-tide diapycnal diffusivity.  A fixed
   !! fraction of the barotropic-to-baroclinic tidal energy conversion
   !! `E(x,y)` [W m-2] dissipates locally above rough topography; the
   !! resulting turbulent diffusivity decays exponentially upward from
   !! the bed with scale `zeta`, converted to a per-layer `Kd` through
   !! the stratification (`1/(dz*(N^2+Omega^2))`).  An INTERIOR closure:
   !! it is ADDED to the other interior diffusivities (KPP/EPBL/PP81/
   !! kappa-shear) and goes through the single `vmix_assemble` gate.
   !!
   !! References:
   !!   * Jayne & St Laurent (2001), GRL 28, "Parameterizing tidal
   !!     dissipation over rough topography" — the parameterized
   !!     conversion `E = 0.5*rho*kappa*<h^2>*U_tide^2`.
   !!   * St Laurent, Simmons & Garrett (2002), GRL 29, "Buoyancy forcing
   !!     by turbulence above rough topography" — the local-dissipation
   !!     fraction `q` and the bottom-anchored exponential vertical
   !!     structure F(z) with decay scale `zeta`.
   !!   * Simmons, Jayne, St Laurent & Weaver (2004), Ocean Modelling 6 —
   !!     the global implementation that combines the two.
   !!   * Melet, Hallberg, Legg & Polzin (2013), J. Phys. Oceanogr. 43 —
   !!     the near-bed N^2 averaged over the roughness height
   !!     (`n2_bot_rough_avg`) and the Polzin vertical structure.
   !!   * Polzin (2009), Ocean Modelling 30, "An abyssal recipe" — the
   !!     algebraic WKB-stretched profile `eps ~ (1 + z*/z_p*)^-2`
   !!     (`profile = "polzin"`).
   !! Knob table: `docs/generated_nml_knobs.md`.
   !!
   !! Static fields: `tideamp_source`/`h2_source = "file"` read time-less
   !! 2-D fields (`load_static_2d_field`); `tidal_mixing_fold_static_fields`
   !! folds the wet mask, `min_zbot` and the `(frac_rough*D)^2` roughness cap
   !! into `c_e(x,y)` once at setup.  `energy_limiter = "reinject"` and
   !! `bed_share = "lump"` make the power bookkeeping exact (OM3 parity);
   !! the defaults keep the discarding clamp and the dropped bed share.
   !!
   !! Polzin profile (`profile = "polzin"`): the remaining upward flux is
   !!   Phi(z*) = P_bot * I_p * z_p* / (z_p* + z*),   I_p = 1 + z_p*/H,
   !! in the WKB-stretched height `z*(z) = int_0^z N^2 dz' / Nbar^2`
   !! (`Nbar^2` the depth mean of `N^2` over the live column, so
   !! `z*(H) = H`), with the stretched scale
   !!   z_p* = f_d*nu*N_ref^2*U / (kappa^2*<h^2>*N_b*Nbar^2),
   !! capped at `f_max*H`; U = 0, <h^2> = 0, N_b = 0, Nbar^2 = 0 or a
   !! capped scale give `z_p* = f_max*H`.  The same deposit/limiter sweep
   !! as the exponential runs on `Phi`; `Phi(0) - Phi(H) = P_bot`.
   !!
   !! Discretization (conservative flux bookkeeping; St Laurent 2002):
   !! a downward TKE flux is tracked and the power dissipated IN each
   !! layer is converted to `Kd = TKE_lay / (rho*dz*(N^2+Omega^2))`.
   !! With the Adcroft-reciprocal normalization `Inv_int =
   !! 1/(1-exp(-H/zeta))` the column-integrated deposited power is
   !! exactly `q*mu*E` (energy conservation, to round-off).  The
   !! continuum equivalent is
   !!   Kd(z) = [q*mu*E/(rho*(N^2+Omega^2))] *
   !!           exp(-z/zeta) / (zeta*(1-exp(-H/zeta))),
   !! z measured UPWARD from the bed.
   !!
   !! Interface convention (same as `rdb_ocean_vmix` / EPBL /
   !! kappa-shear): `kd_int(:,:,K)` lives at the bottom interface of
   !! layer K; global `kd_int(:,:,1)` is the bed and `kd_int(:,:,nz+1)`
   !! the free surface, both forced to exactly 0.  Layers are bottom-up:
   !! k=1 bed, k=nz surface.  The exp-decay anchors at the BED (k=1) and
   !! the sweep runs k=1 -> nz.  The per-layer `Kd_add(k)` is split 50/50
   !! across the layer's two bounding interfaces (St Laurent 2002
   !! interface deposition; see divergence note D2).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, OMEGA_EARTH, &
                            H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY, OMEGA_EARTH, &
                            H_VANISHED, H_DIV_EPS
#endif
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_specvol_derivs
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

   public :: ocean_tidal_mixing_t
   public :: tidal_mixing_compute
   public :: tidal_mixing_merge_into_kt
   public :: tidal_mixing_is_inert
   public :: tidal_mixing_fold_static_fields
   public :: polzin_stretched_scale

   integer, parameter :: NZL = NZ_STACK_MAX
      !! Maximum number of layers a single column kernel can solve.
   integer, parameter :: NZLI = NZ_STACK_MAX + 1
      !! Interface array dimension (= NZL + 1).

   type :: ocean_tidal_mixing_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.

      ! ---- Master switch ----
      logical :: enable = .false.
         !! Master switch.  Default off — existing namelists and tests
         !! stay bit-identical.  Requires `vmix%use_closure` +
         !! thermodynamics (validated at configure).

      ! ---- St-Laurent / Simmons knobs (defaults = paper / MOM6) ----
      real(wp) :: gamma = 0.3333_wp
         !! Local-dissipation fraction `q` (GAMMA_ITIDES).
      real(wp) :: mu = 0.2_wp
         !! Mixing efficiency Gamma_mix (MU_ITIDES).
      real(wp) :: zeta = 500.0_wp
         !! Bottom decay scale (m) (INT_TIDE_DECAY_SCALE).
      real(wp) :: kd_max = 1.0e-2_wp
         !! Per-layer physical Kd cap (m^2/s); < 0 => no cap.  Distinct
         !! from the `vmix_assemble` ceiling (the final clip).
      real(wp) :: prandtl_tidal = 1.0_wp
         !! Kv = prandtl_tidal * Kd into the momentum solve.
      real(wp) :: min_zbot = 0.0_wp
         !! Mask off where column depth H < min_zbot (m).
      real(wp) :: omega2 = OMEGA_EARTH**2
         !! Rotation floor Omega^2 (s^-2) on N^2; physics, not merely a
         !! 1/0 guard.  The floor is plain Omega^2 (Melet et al. 2013
         !! efficiency rescaling N^2/(N^2+Omega^2)), not (2*Omega)^2.

      ! ---- Energy source (v1: prescribed scalar/uniform field) ----
      real(wp) :: e_uniform = 0.0_wp
         !! Uniform bottom energy input E (W m-2); default 0 => inert
         !! even when enabled (structural path live, physics zero).
      logical :: e_compute = .false.
         !! v1.1 state-dependent E: when on, `e_in = min(TKE_coef*N_bot,
         !! e_max)`.  Default off (uses the prescribed `e_uniform`).
      real(wp) :: kappa_itides = 6.2832e-4_wp
         !! Topographic wavenumber kappa (m^-1) for the v1.1 E recompute.
      real(wp) :: kappa_h2 = 1.0_wp
         !! KAPPA_H2_FACTOR for the v1.1 E recompute.
      real(wp) :: utide = 0.0_wp
         !! RMS barotropic tidal velocity amplitude (m/s) for v1.1 E.
      real(wp) :: h2_rough = 0.0_wp
         !! Sub-grid topographic roughness variance <h^2> (m^2) for v1.1 E.
      real(wp) :: frac_rough = 0.1_wp
         !! Roughness clamp: <h^2> <= (frac_rough*H)^2.
      real(wp) :: e_max = 1.0e3_wp
         !! TKE_itide_max cap on E (W m-2).

      ! ---- Static 2-D tide fields and OM3 refinements (all default off) ----
      logical :: use_fields = .false.
         !! True when `tideamp` and/or `h2` come from a file.  Then the
         !! `e_compute` coefficient is read from the setup-folded `c_e`
         !! (and the clamped `h2`) instead of the scalars; false keeps the
         !! scalar path bit-identical.
      logical :: n2_bot_rough_avg = .false.
         !! Near-bed N^2 = roughness-height average, imposed over that span
         !! (Melet et al. 2013).  Default off = raw two-deepest-layer
         !! gradient.  OM3 parity value: `.true.`.
      logical :: limiter_reinject = .false.
         !! `energy_limiter = "reinject"`: cap each interior layer's TKE at
         !! `rho0*h*(N^2+Omega^2)*kd_max` and offer the excess to the layer
         !! above (energy-conserving).  Default off = clip Kd and DISCARD
         !! the clipped power.  OM3 parity value: `.true.` (`kd_max = 0.1`).
      logical :: bed_lump = .false.
         !! `bed_share = "lump"`: the bed layer's share of the dissipated
         !! power stays in the upward flux and is deposited in the first
         !! interior layer.  Default off = the share is discarded (the bed
         !! layer is owned by the bottom boundary layer).  OM3 parity
         !! value: `.true.`.

      logical :: n2_mean_interfaces = .false.
         !! `n2_layer_form = "mean_interfaces"`: each layer's `N^2` is the
         !! mean of its two bounding interface density jumps over its own
         !! thickness (zero jump at the bed and surface interfaces, so the
         !! bed and surface layers get half their one jump).  Default off =
         !! the interface `N^2` is used as the layer value (bit-identical).
         !! Always used under `polzin`.

      ! ---- Vertical structure (default St Laurent exponential) ----
      logical :: polzin = .false.
         !! `profile = "polzin"`: Polzin (2009) algebraic profile in
         !! WKB-stretched height instead of the St Laurent exponential.
         !! Needs `e_compute` (the scale uses U and <h^2>).
      real(wp) :: nu_polzin = 0.0697_wp
         !! NU_POLZIN: dimensionless constant in the decay scale.
      real(wp) :: nbotref_polzin = 9.61e-4_wp
         !! NBOTREF_POLZIN: reference bottom buoyancy frequency (s^-1).
      real(wp) :: polzin_decay_scale_factor = 1.0_wp
         !! POLZIN_DECAY_SCALE_FACTOR `f_d` on the decay scale.
      real(wp) :: polzin_scale_max_factor = 1.0_wp
         !! POLZIN_SCALE_MAX_FACTOR `f_max`: scale capped at `f_max*H`.
      real(wp) :: polzin_min_decay_scale = 0.0_wp
         !! POLZIN_MIN_DECAY_SCALE: floor on the unstretched scale (m).

      ! ---- EOS hookup (shared flat-POD handle from the eos slot) ----
      type(eos_t) :: eos
         !! Value copy of `ocean_state%eos`, set at configure — the
         !! buoyancy derivatives use the SAME EOS the dyn-core runs.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).

      ! ---- Persistent fields ----
      real(wp), allocatable :: e_in(:, :)
         !! Bottom internal-tide energy input E (W m-2), (nx,ny).
      real(wp), allocatable :: tideamp(:, :)
         !! RMS barotropic tidal speed U (m/s), (nx,ny): the file field, or
         !! the scalar `utide` spread uniformly.  Filled at setup.
      real(wp), allocatable :: h2(:, :)
         !! Sub-grid roughness variance <h^2> (m^2), (nx,ny): the file field
         !! (after the `(frac_rough*D)^2` clamp, zeroed on dry cells) or the
         !! scalar `h2_rough` spread uniformly.  Filled at setup.
      real(wp), allocatable :: c_e(:, :)
         !! Folded Jayne & St Laurent coefficient `0.5*rho0*kappa_h2*
         !! kappa_itides*<h^2>*U^2` (kg s-2), (nx,ny), zero on dry cells and
         !! where the resting depth is below `min_zbot`.  Read by the kernel
         !! only when `use_fields`; the per-step energy is `c_e*N_b`.
      real(wp), allocatable :: kd_int(:, :, :)
         !! Tidal diffusivity at interfaces (m^2/s), (nx,ny,nz+1),
         !! global bottom-up: zero at bed (K=1) and surface (K=nz+1).
   contains
      procedure, non_overridable :: init => ocean_tidal_mixing_init
      procedure, non_overridable :: destroy => ocean_tidal_mixing_destroy
      procedure, non_overridable :: enter_data => ocean_tidal_mixing_enter_data
      procedure, non_overridable :: exit_data => ocean_tidal_mixing_exit_data
      procedure, non_overridable :: set_e_uniform => ocean_tidal_mixing_set_e_uniform
      procedure, non_overridable :: bytes => ocean_tidal_mixing_bytes
   end type ocean_tidal_mixing_t

contains

   pure function tidal_mixing_is_inert(enable, e_uniform, e_compute) result(inert)
      !! `.true.` iff tidal mixing is enabled but has no energy source
      !! (PR-6 fail-loud): `enable .and. e_uniform <= 0 .and. .not.
      !! e_compute`.  The bottom-intensified diffusivity `Kd ∝ E`, so a
      !! zero prescribed `e_uniform` with the `e_compute` estimator off
      !! makes `Kd ≡ 0` exactly — the closure runs its whole flux sweep
      !! for nothing.  `e_compute=.true.` supplies E internally, so that
      !! is NOT inert; a disabled closure is not "inert" either (it is
      !! simply off).  Drives a configure abort.
      logical, intent(in) :: enable, e_compute
      real(wp), intent(in) :: e_uniform
      logical :: inert
      inert = enable .and. (e_uniform <= 0.0_wp) .and. (.not. e_compute)
   end function tidal_mixing_is_inert

   pure function polzin_stretched_scale(u, h2, n_b, nbar2, h, nu, n_ref, kappa, &
                                        f_d, f_max, z_min) result(zp_st)
      !! Polzin (2009) stretched decay scale
      !!   z_p* = z_p * N_b^2/Nbar^2,
      !!   z_p  = max(f_d*nu*N_ref^2*U / (kappa^2*<h^2>*N_b^3), z_min),
      !! capped at `f_max*h` (`h` the live column depth).  A degenerate
      !! column (U, <h^2>, N_b or Nbar^2 not positive) returns the cap.
      !! `!$acc routine seq` so the column kernel may call it on device.
      !$acc routine seq
      real(wp), intent(in) :: u, h2, n_b, nbar2, h
         !! RMS tidal speed (m/s), roughness variance (m^2), near-bed N
         !! (s^-1), column-mean N^2 (s^-2), column depth (m).
      real(wp), intent(in) :: nu, n_ref, kappa, f_d, f_max, z_min
         !! NU_POLZIN, NBOTREF_POLZIN, topographic wavenumber, decay-scale
         !! factor, cap factor, floor on the unstretched scale.
      real(wp) :: zp_st
      real(wp) :: zp_max, zp_un

      zp_max = f_max*h
      zp_st = zp_max
      if (u > 0.0_wp .and. h2 > 0.0_wp .and. n_b > 0.0_wp .and. nbar2 > 0.0_wp) then
         zp_un = f_d*nu*n_ref**2*u/(kappa**2*h2*n_b**3)
         zp_un = max(zp_un, z_min)
         zp_st = min(zp_un*n_b**2/nbar2, zp_max)
      end if
      zp_st = max(zp_st, H_DIV_EPS)
   end function polzin_stretched_scale

   ! =================================================================
   ! Lifecycle
   ! =================================================================

   subroutine ocean_tidal_mixing_init(this, grid, nz_ml)
      !! Allocate the persistent fields.  Always allocates (configure
      !! runs after init, so `enable` is not known yet); the off-cost is
      !! one 2D field + one interface field.
      class(ocean_tidal_mixing_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      allocate (this%e_in(nx, ny), source=0.0_wp)
      allocate (this%tideamp(nx, ny), source=0.0_wp)
      allocate (this%h2(nx, ny), source=0.0_wp)
      allocate (this%c_e(nx, ny), source=0.0_wp)
      allocate (this%kd_int(nx, ny, nz + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_tidal_mixing_init

   subroutine ocean_tidal_mixing_destroy(this)
      class(ocean_tidal_mixing_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%e_in)) deallocate (this%e_in)
      if (allocated(this%tideamp)) deallocate (this%tideamp)
      if (allocated(this%h2)) deallocate (this%h2)
      if (allocated(this%c_e)) deallocate (this%c_e)
      if (allocated(this%kd_int)) deallocate (this%kd_int)
   end subroutine ocean_tidal_mixing_destroy

   subroutine ocean_tidal_mixing_enter_data(this)
      class(ocean_tidal_mixing_t), intent(inout) :: this
      select type (this)
      type is (ocean_tidal_mixing_t)
         call ocean_tidal_mixing_enter_data_impl(this)
      end select
   end subroutine ocean_tidal_mixing_enter_data

   subroutine ocean_tidal_mixing_enter_data_impl(this)
      type(ocean_tidal_mixing_t), intent(inout) :: this
      if (allocated(this%e_in)) then
         !$acc enter data copyin(this%e_in)
      end if
      if (allocated(this%tideamp)) then
         !$acc enter data copyin(this%tideamp)
      end if
      if (allocated(this%h2)) then
         !$acc enter data copyin(this%h2)
      end if
      if (allocated(this%c_e)) then
         !$acc enter data copyin(this%c_e)
      end if
      if (allocated(this%kd_int)) then
         !$acc enter data copyin(this%kd_int)
      end if
   end subroutine ocean_tidal_mixing_enter_data_impl

   subroutine ocean_tidal_mixing_exit_data(this)
      class(ocean_tidal_mixing_t), intent(inout) :: this
      select type (this)
      type is (ocean_tidal_mixing_t)
         call ocean_tidal_mixing_exit_data_impl(this)
      end select
   end subroutine ocean_tidal_mixing_exit_data

   subroutine ocean_tidal_mixing_exit_data_impl(this)
      type(ocean_tidal_mixing_t), intent(inout) :: this
      if (allocated(this%kd_int)) then
         !$acc exit data delete(this%kd_int)
      end if
      if (allocated(this%c_e)) then
         !$acc exit data delete(this%c_e)
      end if
      if (allocated(this%h2)) then
         !$acc exit data delete(this%h2)
      end if
      if (allocated(this%tideamp)) then
         !$acc exit data delete(this%tideamp)
      end if
      if (allocated(this%e_in)) then
         !$acc exit data delete(this%e_in)
      end if
   end subroutine ocean_tidal_mixing_exit_data_impl

   subroutine ocean_tidal_mixing_set_e_uniform(this, e_value)
      !! Fill the bottom energy field `e_in` with a uniform value.
      !! Host loop (setup phase); call after `init`, before `enter_data`.
      class(ocean_tidal_mixing_t), intent(inout) :: this
      real(wp), intent(in) :: e_value
      integer :: i, j

      do j = 1, size(this%e_in, 2)
         do i = 1, size(this%e_in, 1)
            this%e_in(i, j) = e_value
         end do
      end do
   end subroutine ocean_tidal_mixing_set_e_uniform

   pure subroutine tidal_mixing_fold_static_fields(this, nx, ny, depth, wet)
      !! Fold the setup-time masks and clamps into the static 2-D fields.
      !! Call once on the host after `tideamp` and `h2` hold the raw
      !! (file or uniform) values and BEFORE `enter_data`.  Per cell, with
      !! `D` the resting column depth:
      !!   * `U = tideamp` where wet and `D >= min_zbot`, else 0;
      !!   * `<h^2> <- min(<h^2>, (frac_rough*D)^2)`, 0 on dry cells
      !!     (Jayne & St Laurent 2001 roughness cap);
      !!   * `c_e = 0.5*rho0*kappa_h2*kappa_itides*<h^2>*U^2`  [kg s-2],
      !!     so the energy input is `E = min(c_e*N_b, e_max)`.
      type(ocean_tidal_mixing_t), intent(inout) :: this
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: depth(nx, ny)
         !! Resting column depth D (m), positive.
      real(wp), intent(in) :: wet(nx, ny)
         !! Wet mask (1 wet, 0 land).
      integer :: i, j
      real(wp) :: dd, h2c, u

      do j = 1, ny
         do i = 1, nx
            dd = max(depth(i, j), 0.0_wp)
            u = this%tideamp(i, j)
            h2c = this%h2(i, j)
            if (wet(i, j) <= 0.0_wp .or. dd < this%min_zbot) then
               u = 0.0_wp
               h2c = 0.0_wp
            end if
            h2c = min(h2c, (this%frac_rough*dd)**2)
            this%tideamp(i, j) = u
            this%h2(i, j) = h2c
            this%c_e(i, j) = 0.5_wp*this%rho0*this%kappa_h2*this%kappa_itides*h2c*u**2
         end do
      end do
   end subroutine tidal_mixing_fold_static_fields

   ! =================================================================
   ! Merge into the vmix interface fields
   ! =================================================================

   pure subroutine tidal_mixing_merge_into_kt(this, nx, ny, nzp1, kv, kt)
      !! Fold the tidal diffusivity into the vmix interface fields,
      !! ADDITIVELY (interior-diffusivity semantics; St Laurent tidal
      !! mixing is an interior source like kappa-shear).  Called EVERY
      !! stage (PP81 rewrites kv/kt each stage; `kd_int` itself refreshes
      !! at thermo cadence).  Interior interfaces only — K=1 (bed) and
      !! K=nzp1 (surface) stay zero in both source and target.
      type(ocean_tidal_mixing_t), intent(in) :: this
      integer, intent(in) :: nx, ny, nzp1
         !! Interface-field extents (explicit shape: assumed-shape
         !! dummies in a `do concurrent` kernel make NVHPC walk the
         !! descriptor with per-launch memcpys — this runs every stage).
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
         !! Momentum viscosity at interfaces; gets prandtl_tidal*kd.
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
         !! Tracer diffusivity at interfaces; gets kd.

      integer :: i, j, k

      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx)
         kt(i, j, k) = kt(i, j, k) + this%kd_int(i, j, k)
         kv(i, j, k) = kv(i, j, k) + this%prandtl_tidal*this%kd_int(i, j, k)
      end do
   end subroutine tidal_mixing_merge_into_kt

   ! =================================================================
   ! Main compute (outer shim)
   ! =================================================================

   pure subroutine tidal_mixing_compute(grid, this, ms, dt)
      !! Run tidal mixing over the domain: fill `this%kd_int` (interface
      !! diffusivity).  Call at thermo cadence with the thermo dt.
      !! Outer shim: dereferences the tracer-registry hTr arrays on the
      !! host (the array-of-DT indirection blocks NVHPC device codegen),
      !! then forwards to the column kernel.  `dt` is unused (the steady
      !! diagnostic Kd does not integrate in time) but carried for
      !! signature parity with the other interior closures.
      type(hgrid_t), intent(in) :: grid
      type(ocean_tidal_mixing_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt

      if (dt < 0.0_wp) return
      if (ms%idx_temperature <= 0 .or. ms%idx_salinity <= 0) return

      call tidal_mixing_column_kernel(grid, this, ms, &
                                      ms%tracers(ms%idx_temperature)%hTr, &
                                      ms%tracers(ms%idx_salinity)%hTr)
   end subroutine tidal_mixing_compute

   pure subroutine tidal_mixing_column_kernel(grid, this, ms, hT, hS)
      !! Per-column St-Laurent flux-bookkeeping sweep.  One
      !! `do concurrent (j, i)` over owned cells; each column is a serial
      !! upward sweep k=1(bed) -> nz(surface) over fixed-size `local()`
      !! arrays (register-resident).  Bottom-up global indexing
      !! throughout — no surface-down flip (the decay anchors at the bed).
      !!
      !! "The bed" is `kbed = ms%k_bot(i,j)`, the first LIVE layer counting
      !! up (`1` off `z_fixed`, ⇒ bit-identical there).  Under `z_fixed`
      !! the layers below it are inert fillers carrying the donor's T/S, so
      !! anchoring at `k = 1` read `N_bot = 0` across two fillers (killing
      !! the `e_compute` energy input on every column shallower than the
      !! nominal stack), exempted a FILLER from deposition instead of the
      !! bed layer, and handed the fillers a TKE share.  The sweep, the
      !! `N_bot` sample, the bed-layer exclusion and the bed end-cap all
      !! start at `kbed`; layers below it get no `Kd`.
      !!
      !! Pipeline per column:
      !!   gather: h (floored), per-layer T,S = hTr/h.
      !!   N^2(k): layer-centred buoyancy gradient (St-Laurent #4 —
      !!           N2_lay, not interface N^2), clamped >= 0.
      !!   E:      prescribed `e_in(i,j)` (v1) or state-dependent
      !!           TKE_coef*N_bot (v1.1, `e_compute`), masked where
      !!           H < min_zbot.
      !!   sweep:  Inv_int normalization, TKE flux bookkeeping, per-layer
      !!           Kd_add capped at kd_max, deposited 50/50 at the two
      !!           bounding interfaces.
      type(hgrid_t), intent(in) :: grid
      type(ocean_tidal_mixing_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      ! assumed-shape-ok: tracer registry outer-shim — caller host-dereferences
      ! ms%tracers(idx)%hTr before passing; size varies per tracer slot
      ! (see CLAUDE.md "outer-shim + flat-impl" pattern); thermo cadence.
      real(wp), intent(in) :: hT(:, :, :)
         !! Temperature tracer hTr (degC*m), host-dereferenced.
      real(wp), intent(in) :: hS(:, :, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
         !! Salinity tracer hTr (PSU*m), host-dereferenced.

      integer :: i, j, k, nx, ny, nz, kbed
      real(wp) :: e_col, h_tot, n_bot, tke_coef, h2c, e_cap
      real(wp) :: inv_int, hz, tke_bot, tke_rem, z_top, frac_top, tke_lay
      real(wp) :: dz_eff, denom, kd_lay
      real(wp) :: h_col(NZL), t_col(NZL), s_col(NZL)
      real(wp) :: n2_col(NZL), kd_lay_arr(NZL)
      real(wp) :: dbuoy_t, dbuoy_s, p_int
      real(wp) :: dsv_dt_k, dsv_ds_k, t_int, s_int
      real(wp) :: gamma_l, mu_l, zeta_l, kd_max_l, omega2_l, rho0_l, izeta
      real(wp) :: tke_used, tke_cap, h_amp, sum_db, sum_dz, zc_k, span_k, n2_b
      real(wp) :: u_col, nbar2, zp_st, i_p, zs_top, zs_max
      integer :: klast
      logical :: pz_col, mean_if
      real(wp) :: nb_int
      real(wp) :: n2_lay(NZL)

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      gamma_l = this%gamma
      mu_l = this%mu
      zeta_l = this%zeta
      kd_max_l = this%kd_max
      omega2_l = this%omega2
      rho0_l = this%rho0
      izeta = 1.0_wp/zeta_l

      do concurrent(j=1:ny, i=1:nx) &
         local(k, kbed, e_col, h_tot, n_bot, tke_coef, h2c, e_cap, &
               inv_int, hz, tke_bot, tke_rem, z_top, frac_top, tke_lay, &
               dz_eff, denom, kd_lay, h_col, t_col, s_col, n2_col, kd_lay_arr, &
               dbuoy_t, dbuoy_s, p_int, dsv_dt_k, dsv_ds_k, t_int, s_int, &
               tke_used, tke_cap, h_amp, sum_db, sum_dz, zc_k, span_k, n2_b, klast, &
               u_col, nbar2, zp_st, i_p, zs_top, zs_max, pz_col, mean_if, n2_lay, nb_int)

         ! ---- Default output: zero contribution (overwritten if wet) ----
         do k = 1, nz + 1
            this%kd_int(i, j, k) = 0.0_wp
         end do

         if (ms%wet_mask(i, j) > 0.0_wp) then
            kbed = ms%k_bot(i, j)
            ! ---- Gather (bottom-up; floor thickness, back out T,S) ----
            ! `h_tot` is the LIVE column (the fillers below `kbed` are not
            ! water the internal tide can mix).
            h_tot = 0.0_wp
            do k = 1, nz
               dz_eff = max(ms%h_layer(i, j, k), H_VANISHED)
               h_col(k) = dz_eff
               denom = 1.0_wp/max(ms%h_layer(i, j, k), H_DIV_EPS)
               t_col(k) = hT(i, j, k)*denom
               s_col(k) = hS(i, j, k)*denom
               if (k >= kbed) h_tot = h_tot + ms%h_layer(i, j, k)
            end do

            ! ---- Layer-centred N^2 (St-Laurent #4: N2_lay) ----
            ! Buoyancy gradient between the centres of layer k and the
            ! layer above (k+1, bottom-up), over their centre spacing
            ! 0.5*(h(k)+h(k+1)).  Pressure accumulates downward from the
            ! surface for the EOS in-situ derivatives.  Surface layer
            ! (k=nz) has no overlying layer -> N2=0 there.
            p_int = 0.0_wp
            n2_col(nz) = 0.0_wp
            do k = nz - 1, 1, -1
               ! Walk pressure down from the surface to the k/k+1 interface.
               p_int = p_int + GRAVITY*rho0_l*h_col(k + 1)
               t_int = 0.5_wp*(t_col(k) + t_col(k + 1))
               s_int = 0.5_wp*(s_col(k) + s_col(k + 1))
               call eos_specvol_derivs(this%eos, t_int, s_int, p_int, &
                                       dsv_dt_k, dsv_ds_k)
               dbuoy_t = GRAVITY*rho0_l*dsv_dt_k
               dbuoy_s = GRAVITY*rho0_l*dsv_ds_k
               dz_eff = 0.5_wp*(h_col(k) + h_col(k + 1))
               n2_col(k) = (dbuoy_t*(t_col(k + 1) - t_col(k)) + &
                            dbuoy_s*(s_col(k + 1) - s_col(k)))/ &
                           max(dz_eff, H_VANISHED)
               if (n2_col(k) < 0.0_wp) n2_col(k) = 0.0_wp
            end do

            ! ---- Energy input E(x,y) ----
            ! v1: prescribed field.  v1.1 (`e_compute`): state-dependent
            ! internal-tide generation, Jayne & St Laurent (2001):
            !   E = min( 0.5*rho0*kappa_h2*kappa_itides*<h^2>*U_tide^2 * N_bot,
            !            e_max )
            ! with the roughness clamp <h^2> <= (frac_rough*H)^2.  The rho0
            ! factor sets E in W/m^2 (the same units as the prescribed e_in),
            ! so it lands correctly in the rho0*dz*(N^2+Omega^2) Kd divisor.
            ! N_bot = sqrt(N^2) at the bed-most interior interface
            ! (n2_col(kbed), the gradient across the two deepest LIVE layers).
            if (this%e_compute .or. this%n2_bot_rough_avg .or. this%polzin) then
               if (this%use_fields) then
                  h2c = this%h2(i, j)   ! setup-clamped file/uniform field
               else
                  h2c = min(this%h2_rough, (this%frac_rough*h_tot)**2)
               end if
            end if

            ! ---- Near-bed N^2: roughness-height average (opt-in) ----
            ! Melet et al. (2013).  Walk up from the bed interface by
            ! interface, accumulating the buoyancy jump and the centre
            ! spacing; z is the centre height of the layer above the
            ! interface just taken.  Once z exceeds h_amp = sqrt(<h^2>),
            ! ONE MORE interface (if any) is accumulated and the walk stops,
            ! so even h_amp = 0 averages two interfaces (hence not
            ! bit-identical to the raw gradient; gated by the knob).
            ! N_b^2 = sum jumps / sum spans replaces the layer N^2 over
            ! every interface taken.
            mean_if = this%polzin .or. this%n2_mean_interfaces
            if (mean_if) then
               ! Layer N^2 from the RAW interface jumps (db = n2*span):
               ! mean of the two bounding jumps over the layer thickness;
               ! the bed and surface interfaces carry no jump.
               do k = 1, nz
                  n2_lay(k) = 0.0_wp
               end do
               do k = kbed, nz
                  sum_db = 0.0_wp
                  if (k > kbed) sum_db = sum_db + n2_col(k - 1)*0.5_wp*(h_col(k - 1) + h_col(k))
                  if (k < nz) sum_db = sum_db + n2_col(k)*0.5_wp*(h_col(k) + h_col(k + 1))
                  n2_lay(k) = 0.5_wp*sum_db/max(h_col(k), H_VANISHED)
               end do
            end if
            if (this%n2_bot_rough_avg) then
               h_amp = sqrt(max(h2c, 0.0_wp))
               sum_db = 0.0_wp
               sum_dz = 0.0_wp
               zc_k = 0.5_wp*h_col(kbed)
               klast = kbed - 1
               do k = kbed, nz - 1
                  span_k = 0.5_wp*(h_col(k) + h_col(k + 1))
                  sum_db = sum_db + n2_col(k)*span_k
                  sum_dz = sum_dz + span_k
                  klast = k
                  zc_k = zc_k + span_k
                  if (zc_k > h_amp) then
                     if (k + 1 <= nz - 1) then
                        span_k = 0.5_wp*(h_col(k + 1) + h_col(k + 2))
                        sum_db = sum_db + n2_col(k + 1)*span_k
                        sum_dz = sum_dz + span_k
                        klast = k + 1
                     end if
                     exit
                  end if
               end do
               if (klast >= kbed) then
                  n2_b = sum_db/max(sum_dz, H_VANISHED)
                  do k = kbed, klast
                     n2_col(k) = n2_b
                     if (mean_if) n2_lay(k) = n2_b
                  end do
               end if
            end if

            if (this%e_compute) then
               if (this%use_fields) then
                  tke_coef = this%c_e(i, j)
               else
                  tke_coef = 0.5_wp*rho0_l*this%kappa_h2*this%kappa_itides*h2c*this%utide**2
               end if
               n_bot = sqrt(max(n2_col(kbed), 0.0_wp))
               e_col = tke_coef*n_bot
               e_cap = this%e_max
               if (e_col > e_cap) e_col = e_cap
            else
               e_col = this%e_in(i, j)
            end if
            ! Mask off shallow columns (H < min_zbot).
            if (h_tot < this%min_zbot) e_col = 0.0_wp
            ! Layer-form N^2 feeds the profile and TKE->Kd conversion; the
            ! interface near-bed value (`nb_int`) already set the energy.
            nb_int = sqrt(max(n2_col(kbed), 0.0_wp))
            if (mean_if) then
               do k = kbed, nz
                  n2_col(k) = n2_lay(k)
               end do
            end if

            ! ---- Inv_int normalization (Adcroft reciprocal) ----
            ! L'Hospital degenerate-thin-column guard: H/zeta -> 0 => 1.
            hz = h_tot*izeta
            if (hz < 1.0e-14_wp) then
               inv_int = 1.0_wp
            else
               inv_int = 1.0_wp/(1.0_wp - exp(-hz))
            end if

            ! ---- Polzin (2009) scale: per-column scalars ----
            ! Replaces `inv_int` by `I_p` and the exponential fraction by
            ! the stretched algebraic one in the sweep below.  `nbar2` is
            ! the live-column mean of N^2 (after any near-bed replacement).
            pz_col = this%polzin
            if (pz_col) then
               if (this%use_fields) then
                  u_col = this%tideamp(i, j)
               else
                  u_col = this%utide
               end if
               nbar2 = 0.0_wp
               do k = kbed, nz
                  nbar2 = nbar2 + n2_col(k)*h_col(k)
               end do
               nbar2 = nbar2/max(h_tot, H_DIV_EPS)
               n_bot = nb_int
               zp_st = polzin_stretched_scale(u_col, h2c, n_bot, nbar2, h_tot, &
                                              this%nu_polzin, this%nbotref_polzin, &
                                              this%kappa_itides, &
                                              this%polzin_decay_scale_factor, &
                                              this%polzin_scale_max_factor, &
                                              this%polzin_min_decay_scale)
               i_p = 1.0_wp + zp_st/max(h_tot, H_DIV_EPS)
               zs_max = 1.0e14_wp*max(h_tot, H_DIV_EPS)*nbar2
               zs_top = 0.0_wp
               inv_int = i_p
            end if

            ! ---- Flux bookkeeping sweep (St-Laurent 2002), bed -> surf ----
            ! TKE_bot = q*mu*E is the total locally-dissipated power; the
            ! 1/rho0 Boussinesq factor lives in the TKE->Kd divisor
            ! (rho0*dz*(N^2+Omega^2)), so Sum(TKE_lay) == q*mu*E exactly
            ! (energy conservation is independent of the rho placement).
            tke_bot = gamma_l*mu_l*e_col
            tke_rem = inv_int*tke_bot
            z_top = 0.0_wp
            do k = 1, kbed - 1
               kd_lay_arr(k) = 0.0_wp   ! inert bed filler: no share of the TKE
            end do
            do k = kbed, nz
               z_top = z_top + h_col(k)
               if (pz_col) then
                  ! Stretched height: skip the increment when N^2 ~ 0 (the
                  ! column is then degenerate anyway).
                  if (nbar2 <= 0.0_wp) then
                     zs_top = z_top   ! unstratified: z* = z
                  else if (h_col(k)*n2_col(k) < zs_max) then
                     zs_top = zs_top + h_col(k)*n2_col(k)/nbar2
                  end if
                  frac_top = inv_int*zp_st/(zp_st + zs_top)
               else
                  frac_top = inv_int*exp(-z_top*izeta)
               end if
               tke_lay = tke_rem - tke_bot*frac_top

               ! TKE_to_Kd = 1/(rho0*dz*(N^2+Omega^2)).  dz floored at
               ! H_VANISHED (the vdiff vanishing-layer mass-drop
               ! precedent) and the divisor armoured with H_DIV_EPS.
               dz_eff = max(h_col(k), H_VANISHED)
               denom = rho0_l*dz_eff*(n2_col(k) + omega2_l)
               denom = max(denom, H_DIV_EPS)

               ! Power taken out of the upward flux by this layer.  Default:
               ! the whole natural share.  `bed_lump`: the bed layer takes
               ! none (its share rides up into the first interior layer).
               ! `limiter_reinject`: an interior layer takes at most
               ! `denom*kd_max` (the TKE that gives exactly Kd = kd_max); the
               ! excess stays in the flux for the layer above, and a
               ! vanished layer takes nothing.
               tke_used = tke_lay
               if (this%bed_lump .and. k == kbed) then
                  tke_used = 0.0_wp
               else if (this%limiter_reinject .and. k > kbed .and. k < nz) then
                  ! vanished-ok: a vanished layer is skipped (takes no power)
                  if (ms%h_layer(i, j, k) <= H_VANISHED) then
                     tke_used = 0.0_wp
                  else if (kd_max_l >= 0.0_wp) then
                     tke_cap = denom*kd_max_l
                     if (tke_used > tke_cap) tke_used = tke_cap
                  end if
               end if
               tke_rem = tke_rem - tke_used
               kd_lay = tke_used/denom
               ! Per-layer physical cap (St-Laurent).  Default
               ! (`energy_limiter = "clamp"`): a bare Kd_max clip that
               ! DISCARDS the clipped power; the energy-conservation
               ! invariant is asserted on the PRE-clip TKE_lay accordingly.
               ! Under `reinject` the layer TKE is already capped above, so
               ! this clip only guards round-off.
               if (kd_max_l >= 0.0_wp .and. kd_lay > kd_max_l) kd_lay = kd_max_l
               kd_lay_arr(k) = kd_lay
            end do

            ! ---- D1: exclude the BED and SURFACE LAYERS from deposition ----
            ! MOM6 (St Laurent/Simmons; MOM_set_diffusivity find_TKE_to_Kd,
            ! MOM-inspired) sets TKE_to_Kd(i,1)=TKE_to_Kd(i,nz)=0 — BOTH
            ! endpoint layers get zero Kd_add, and the deposit loop runs the
            ! interior only (top-down do k=nz-1,2,-1).  In MOM6 top-down
            ! indexing layer 1 is the surface and layer nz the bottom; in
            ! Roundabout bottom-up k=nz is the surface and k=1 the bed, so we
            ! zero kd_lay_arr at BOTH ends.  The bed layer (`kbed`: k=1, or
            ! the first LIVE layer `k_bot` under z_fixed) is owned by the BBL
            ! drag; the surface layer (k=nz) is the mixed layer,
            ! owned downstream by EPBL/KPP, and additionally has no overlying
            ! layer so its kernel N^2 is 0 (its Omega^2-only inflated Kd would
            ! over-mix interface K=nz on shallow energetic columns).  Zeroing
            ! the deposited Kd here does NOT touch the TKE energy bookkeeping
            ! above (tke_lay/tke_rem are unchanged) — the energy invariant is
            ! asserted on the pre-deposit TKE_lay, mirroring the original
            ! surface-only exclusion.  Zero both before the 50/50 split.
            kd_lay_arr(kbed) = 0.0_wp
            kd_lay_arr(nz) = 0.0_wp

            ! ---- Deposit per-layer Kd_add 50/50 at the two bounding
            ! interfaces (St-Laurent 2002).  DIVERGENCE D2: MOM6 splits
            ! the layer power 0.5/0.5 across its bounding interfaces; we
            ! mirror that interface deposition.  Global interface K is the
            ! BOTTOM interface of layer K, so layer k contributes to
            ! interfaces k (its bed) and k+1 (its top).
            do k = 1, nz
               this%kd_int(i, j, k) = this%kd_int(i, j, k) + 0.5_wp*kd_lay_arr(k)
               this%kd_int(i, j, k + 1) = this%kd_int(i, j, k + 1) + 0.5_wp*kd_lay_arr(k)
            end do
            ! D1 (cont.): zero the bed (K=kbed) and surface (K=nz+1) interface
            ! end-caps — `vmix_assemble` / EPBL / BBL own them.  With BOTH
            ! the bed LAYER (k=kbed) and surface LAYER (k=nz) now excluded
            ! above, interface K=kbed+1 receives nothing from k=kbed (only
            ! k=kbed+1's bed half)
            ! and interface K=nz nothing from k=nz (only k=nz-1's top half) —
            ! the full endpoint exclusion.  The near-bed N^2 is the
            ! raw per-layer value unless `n2_bot_rough_avg` (above).
            this%kd_int(i, j, kbed) = 0.0_wp
            this%kd_int(i, j, nz + 1) = 0.0_wp
         end if
      end do
   end subroutine tidal_mixing_column_kernel

   pure function ocean_tidal_mixing_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the tidal mixing slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_tidal_mixing_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%e_in) &
               + arr_bytes(this%tideamp) &
               + arr_bytes(this%h2) &
               + arr_bytes(this%c_e) &
               + arr_bytes(this%kd_int)
   end function ocean_tidal_mixing_bytes

end module rdb_ocean_tidal_mixing
