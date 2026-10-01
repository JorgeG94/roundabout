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
   !! Knob table: `docs/generated_nml_knobs.md`.
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

      ! ---- EOS hookup (shared flat-POD handle from the eos slot) ----
      type(eos_t) :: eos
         !! Value copy of `ocean_state%eos`, set at configure — the
         !! buoyancy derivatives use the SAME EOS the dyn-core runs.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).

      ! ---- Persistent fields ----
      real(wp), allocatable :: e_in(:, :)
         !! Bottom internal-tide energy input E (W m-2), (nx,ny).
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
      allocate (this%kd_int(nx, ny, nz + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_tidal_mixing_init

   subroutine ocean_tidal_mixing_destroy(this)
      class(ocean_tidal_mixing_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%e_in)) deallocate (this%e_in)
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
      !! throughout — no surface-down flip (the decay anchors at the bed,
      !! which IS k=1).
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

      integer :: i, j, k, nx, ny, nz
      real(wp) :: e_col, h_tot, n_bot, tke_coef, h2c, e_cap
      real(wp) :: inv_int, hz, tke_bot, tke_rem, z_top, frac_top, tke_lay
      real(wp) :: dz_eff, denom, kd_lay
      real(wp) :: h_col(NZL), t_col(NZL), s_col(NZL)
      real(wp) :: n2_col(NZL), kd_lay_arr(NZL)
      real(wp) :: dbuoy_t, dbuoy_s, p_int
      real(wp) :: dsv_dt_k, dsv_ds_k, t_int, s_int
      real(wp) :: gamma_l, mu_l, zeta_l, kd_max_l, omega2_l, rho0_l, izeta

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
         local(k, e_col, h_tot, n_bot, tke_coef, h2c, e_cap, &
               inv_int, hz, tke_bot, tke_rem, z_top, frac_top, tke_lay, &
               dz_eff, denom, kd_lay, h_col, t_col, s_col, n2_col, kd_lay_arr, &
               dbuoy_t, dbuoy_s, p_int, dsv_dt_k, dsv_ds_k, t_int, s_int)

         ! ---- Default output: zero contribution (overwritten if wet) ----
         do k = 1, nz + 1
            this%kd_int(i, j, k) = 0.0_wp
         end do

         if (ms%wet_mask(i, j) > 0.0_wp) then
            ! ---- Gather (bottom-up; floor thickness, back out T,S) ----
            h_tot = 0.0_wp
            do k = 1, nz
               dz_eff = max(ms%h_layer(i, j, k), H_VANISHED)
               h_col(k) = dz_eff
               denom = 1.0_wp/max(ms%h_layer(i, j, k), H_DIV_EPS)
               t_col(k) = hT(i, j, k)*denom
               s_col(k) = hS(i, j, k)*denom
               h_tot = h_tot + ms%h_layer(i, j, k)
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
            ! (n2_col(1), the gradient across the two deepest layers).
            if (this%e_compute) then
               h2c = min(this%h2_rough, (this%frac_rough*h_tot)**2)
               tke_coef = 0.5_wp*rho0_l*this%kappa_h2*this%kappa_itides*h2c*this%utide**2
               n_bot = sqrt(max(n2_col(1), 0.0_wp))
               e_col = tke_coef*n_bot
               e_cap = this%e_max
               if (e_col > e_cap) e_col = e_cap
            else
               e_col = this%e_in(i, j)
            end if
            ! Mask off shallow columns (H < min_zbot).
            if (h_tot < this%min_zbot) e_col = 0.0_wp

            ! ---- Inv_int normalization (Adcroft reciprocal) ----
            ! L'Hospital degenerate-thin-column guard: H/zeta -> 0 => 1.
            hz = h_tot*izeta
            if (hz < 1.0e-14_wp) then
               inv_int = 1.0_wp
            else
               inv_int = 1.0_wp/(1.0_wp - exp(-hz))
            end if

            ! ---- Flux bookkeeping sweep (St-Laurent 2002), bed -> surf ----
            ! TKE_bot = q*mu*E is the total locally-dissipated power; the
            ! 1/rho0 Boussinesq factor lives in the TKE->Kd divisor
            ! (rho0*dz*(N^2+Omega^2)), so Sum(TKE_lay) == q*mu*E exactly
            ! (energy conservation is independent of the rho placement).
            tke_bot = gamma_l*mu_l*e_col
            tke_rem = inv_int*tke_bot
            z_top = 0.0_wp
            do k = 1, nz
               z_top = z_top + h_col(k)
               frac_top = inv_int*exp(-z_top*izeta)
               tke_lay = tke_rem - tke_bot*frac_top
               tke_rem = tke_rem - tke_lay

               ! TKE_to_Kd = 1/(rho0*dz*(N^2+Omega^2)).  dz floored at
               ! H_VANISHED (the vdiff vanishing-layer mass-drop
               ! precedent) and the divisor armoured with H_DIV_EPS.
               dz_eff = max(h_col(k), H_VANISHED)
               denom = rho0_l*dz_eff*(n2_col(k) + omega2_l)
               denom = max(denom, H_DIV_EPS)
               kd_lay = tke_lay/denom
               ! Per-layer physical cap (St-Laurent).  DIVERGENCE D3:
               ! MOM6's max_TKE limiter re-injects the un-used power into
               ! the layer above (energy-conserving under clipping); this
               ! bare Kd_max clamp DISCARDS the clipped power.  v1 accepts
               ! the discard — it is a rarely-active safety cap, not the
               ! deposition mechanism.  The energy-conservation invariant
               ! is asserted on the PRE-clip TKE_lay accordingly.
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
            ! zero kd_lay_arr at BOTH ends.  The bed layer (k=1) is owned by
            ! the BBL drag; the surface layer (k=nz) is the mixed layer,
            ! owned downstream by EPBL/KPP, and additionally has no overlying
            ! layer so its kernel N^2 is 0 (its Omega^2-only inflated Kd would
            ! over-mix interface K=nz on shallow energetic columns).  Zeroing
            ! the deposited Kd here does NOT touch the TKE energy bookkeeping
            ! above (tke_lay/tke_rem are unchanged) — the energy invariant is
            ! asserted on the pre-deposit TKE_lay, mirroring the original
            ! surface-only exclusion.  Zero both before the 50/50 split.
            kd_lay_arr(1) = 0.0_wp
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
            ! D1 (cont.): zero the bed (K=1) and surface (K=nz+1) interface
            ! end-caps — `vmix_assemble` / EPBL / BBL own them.  With BOTH
            ! the bed LAYER (k=1) and surface LAYER (k=nz) now excluded above,
            ! interface K=2 receives nothing from k=1 (only k=2's bed half)
            ! and interface K=nz nothing from k=nz (only k=nz-1's top half) —
            ! the full MOM6 endpoint exclusion.  DIVERGENCE D4
            ! (BBL N^2 override): MOM6 replaces the near-bed N^2 with a
            ! roughness-height BBL average; v1 uses the raw per-layer N^2
            ! (a v1.1 refinement).
            this%kd_int(i, j, 1) = 0.0_wp
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
               + arr_bytes(this%kd_int)
   end function ocean_tidal_mixing_bytes

end module rdb_ocean_tidal_mixing
