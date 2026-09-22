!! Energetics-based planetary boundary layer (EPBL) for the ocean core.
module rdb_ocean_epbl
   !! Prognostic-energy surface boundary layer: each thermo step the
   !! wind supplies mechanical TKE (mstar·rho0·u*^3·dt) and surface
   !! buoyancy loss supplies convective PE (nstar-weighted); the
   !! scheme spends that energy interface by interface — the
   !! diffusivity at each interface is the largest value whose
   !! implicit-diffusion PE cost the remaining TKE can pay.  Mixing
   !! stops where the energy runs out, which IS the mixed-layer
   !! depth.  Energetically closed by construction; unconditionally
   !! stable (the sweep performs the forward elimination of the
   !! backward-Euler vertical-diffusion solve it feeds).
   !!
   !! References: Reichl & Hallberg (2018), Ocean Modelling 132
   !! ("ePBL"); the gravity-wave column-height correction follows the
   !! same paper's available-PE bookkeeping.  Knob table:
   !! `docs/generated_nml_knobs.md`.
   !!
   !! Scope (PR 1): surface boundary layer only — no bottom-boundary
   !! ePBL, no mean-KE→TKE conversion (so the closed-form direct
   !! energy solve applies; no inner Newton loop), no Langmuir
   !! enhancement (PR 2).  Boussinesq SI.
   !!
   !! Interface convention (same as `rdb_ocean_vmix`): `kd_int(:,:,k)`
   !! lives at the bottom interface of layer k; `kd_int(:,:,1)` is the
   !! bed and `kd_int(:,:,nz+1)` the free surface, both forced to 0.
   !! Layers are bottom-up: k=1 bed, k=nz surface.
   !!
   !! EPBL replaces the KPP overlay (mutually exclusive); the PP81
   !! interior closure + background diffusivity continue to run
   !! underneath, and `epbl_merge_into_kv_kt` folds `kd_int` into
   !! `vmix%kv` / `vmix%kt` additively (MOM6 EPBL_IS_ADDITIVE) or by
   !! max, every stage; `epbl_compute` itself runs at thermo cadence.
   use rdb_constants, only: wp, GRAVITY, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     sw_transmission, sw_pe_cost_shape
   use rdb_eos, only: eos_t, eos_specvol_derivs
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_epbl_t
   public :: epbl_compute
   public :: epbl_merge_into_kv_kt
   public :: epbl_find_mstar
   public :: epbl_mixlen_shape
   public :: epbl_lf17_wave_state
   public :: epbl_lf17_la
   public :: epbl_lt_enhance
   public :: parse_epbl_mstar_scheme
   public :: parse_epbl_vstar_scheme
   public :: parse_epbl_combine
   public :: parse_epbl_lt_scheme

   ! mstar (mechanical TKE / u*^3) scheme tags.
   integer, parameter, public :: EPBL_MSTAR_CONSTANT = 1
      !! Fixed mstar (the simplest credible configuration).
   integer, parameter, public :: EPBL_MSTAR_OM4 = 2
      !! Ekman/Obukhov-balance form used by MOM6's OM4 production
      !! config (Reichl & Hallberg 2018 Appendix).
   integer, parameter, public :: EPBL_MSTAR_RH18 = 3
      !! Reichl & Hallberg (2018) eq. fits (cN1..cS2 coefficients).

   ! Turbulent-velocity-scale scheme tags.
   integer, parameter, public :: EPBL_VSTAR_CUBE_ROOT = 1
      !! vstar = vstar_scale_fac * (TKE / (dt rho0))^(1/3).
   integer, parameter, public :: EPBL_VSTAR_RH18 = 2
      !! Separate mechanical (u*-proportional, surface-decaying) +
      !! convective cube-root contributions (RH18).

   ! Combination of kd_epbl with the interior closure's kv/kt.
   integer, parameter, public :: EPBL_COMBINE_ADD = 1
      !! kt += kd ; kv += prandtl*kd  (MOM6 EPBL_IS_ADDITIVE=.true.)
   integer, parameter, public :: EPBL_COMBINE_MAX = 2
      !! kt = max(kt, kd) ; kv = max(kv, prandtl*kd)

   ! Langmuir-enhancement scheme tags (Reichl & Li 2019 forms).
   integer, parameter, public :: EPBL_LT_NONE = 0
   integer, parameter, public :: EPBL_LT_RESCALE = 1
      !! mstar *= min(max_enh, 1 + coef * La^exp)  (multiplicative)
   integer, parameter, public :: EPBL_LT_ADDITIVE = 2
      !! mstar += coef * La^exp

   ! Thickness regulariser added to every layer thickness used inside
   ! the energy solve (MOM6 H_subroundoff role): keeps the pivot
   ! products bdt1 = hp_a*hp_b representable when a layer vanishes.
   ! Division-safety only (D4 role 2): aliases the constant of record
   ! `H_DIV_EPS` from rdb_constants — the local name is kept so the many
   ! use sites below read unchanged (bitwise no-op: 1e-20 == 1e-20).
   real(wp), parameter :: H_NEGLECT = H_DIV_EPS

   ! ---- LF17 / COARE 3.5 empirical constants ----
   ! Parts of the cited formulas, NOT namelist knobs.  Li & Fox-Kemper
   ! (2017) statistical wave state; Edson et al. (2013) COARE 3.5
   ! drag; Webb & Fox-Kemper (2015) / Breivik et al. (2016)
   ! surface-layer-averaged Stokes drift over a Phillips spectrum.
   real(wp), parameter :: LT_VONKAR_WAVES = 0.40_wp
      !! von Karman constant used by the COARE u*->U10 inversion.
   real(wp), parameter :: LT_NU_AIR = 1.0e-6_wp
      !! Kinematic viscosity of air (m^2/s) for the smooth-flow z0.
   real(wp), parameter :: LT_RHO_AIR = 1.225_wp
      !! Air density (kg/m^3) for the water->air u* conversion.
   real(wp), parameter :: LT_CHARNOCK_MIN = 0.028_wp
      !! Cap on the Charnock parameter.
   real(wp), parameter :: LT_CHARNOCK_SLOPE = 0.0017_wp
      !! d(alpha_Charnock)/dU10 (s/m).
   real(wp), parameter :: LT_CHARNOCK_ICPT = -0.005_wp
      !! Charnock intercept at U10 = 0.
   real(wp), parameter :: LT_US_TO_U10 = 0.0162_wp
      !! Surface Stokes drift / U10 (Webb 2011).
   real(wp), parameter :: LT_SWH_FROM_U10SQ = 0.0246_wp
      !! Significant wave height = c * U10^2 (s^2/m).
   real(wp), parameter :: LT_U19P5_TO_U10 = 1.075_wp
      !! Pierson-Moskowitz U19.5/U10 ratio.
   real(wp), parameter :: LT_FM_INTO_FP = 1.296_wp
      !! Mean / peak frequency ratio (Webb 2011).
   real(wp), parameter :: LT_R_LOSS = 0.667_wp
      !! Stokes-transport loss ratio.
   real(wp), parameter :: LT_PI = 4.0_wp*atan(1.0_wp)

   type :: ocean_epbl_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.

      ! ---- Scheme selection + master switch ----
      logical :: enable = .false.
         !! Master switch.  Default off — all existing namelists and
         !! tests stay bit-identical.  Requires `vmix%use_closure`;
         !! disables the KPP overlay (configure logs the override).
      integer :: mstar_scheme = EPBL_MSTAR_OM4
      integer :: vstar_scheme = EPBL_VSTAR_CUBE_ROOT
      integer :: combine_mode = EPBL_COMBINE_ADD

      ! ---- mstar knobs ----
      real(wp) :: mstar_const = 1.2_wp
         !! Constant-scheme mstar (MOM6 MSTAR).
      real(wp) :: mstar_cap = -1.0_wp
         !! Cap on mstar for the OM4/RH18 schemes; off when < 0.
      real(wp) :: mstar_coef1 = 0.3_wp
         !! OM4 stabilizing-balance coefficient (MOM6 MSTAR2_COEF1).
      real(wp) :: c_ek = 0.085_wp
         !! OM4 Ekman-limit coefficient (MOM6 MSTAR2_COEF2).
      real(wp) :: mstar_conv_adj = 0.0_wp
         !! Reduce mechanical mstar when convection dominates, in
         !! [0,1]; 0 = off (MOM6 MSTAR_CONV_ADJ).
      real(wp) :: rh18_cn1 = 0.275_wp
      real(wp) :: rh18_cn2 = 8.0_wp
      real(wp) :: rh18_cn3 = -5.0_wp
      real(wp) :: rh18_cs1 = 0.2_wp
      real(wp) :: rh18_cs2 = 0.4_wp
         !! RH18 mstar fit coefficients (paper values).

      ! ---- Energetics knobs ----
      real(wp) :: nstar = 0.2_wp
         !! Fraction of convectively released PE that becomes
         !! entrainment-driving TKE (MOM6 NSTAR).
      real(wp) :: tke_decay = 2.5_wp
         !! Ratio of the natural Ekman depth to the mechanical-TKE
         !! decay scale (MOM6 TKE_DECAY).
      real(wp) :: wstar_ustar_coef = 1.0_wp
         !! Weight of the convective reservoir in the velocity scale.
      real(wp) :: vstar_scale_fac = 1.0_wp
         !! Overall vstar multiplier (∝ diffusivity).
      real(wp) :: vstar_surf_fac = 1.2_wp
         !! RH18 vstar scheme: mechanical surface vstar ∝ u* factor.
      real(wp) :: von_karman = 0.41_wp
         !! kappa in Kd = vstar * kappa * mixing_length.
      real(wp) :: ekman_scale_coef = 1.0_wp
         !! Rotational inhibition of the mixing length.
      real(wp) :: min_mix_len = 0.0_wp
         !! Mixing-length floor (m).
      real(wp) :: mixlen_exponent = 2.0_wp
         !! Shape-function exponent (2 = KPP-like; fast path).
      real(wp) :: translay_scale = 0.1_wp
         !! Transition-layer floor of the shape function; must be in
         !! [0,1) when `mld_iteration` is on.

      ! ---- MLD iteration knobs ----
      logical :: mld_iteration = .true.
         !! Iterate the boundary-layer depth to self-consistency with
         !! the mixing-length shape function (MOM6 USE_MLD_ITERATION).
      real(wp) :: mld_tol = 1.0_wp
         !! Convergence tolerance on the MLD root-find (m).
      integer :: mld_max_its = 20
      logical :: mld_bisection = .false.
         !! Plain bisection instead of false-position.
      logical :: mld_use_prev_guess = .false.
         !! Seed the iteration from the previous step's MLD (faster,
         !! ~1-3 iterations; off = always start at half depth).

      ! ---- Environment ----
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).
      real(wp) :: omega = 7.2921e-5_wp
         !! Earth rotation rate (1/s), for the omega_frac blend.
      real(wp) :: omega_frac = 0.0_wp
         !! Blend |f| with 2*Omega: absf = sqrt((1-of) f^2 + of 4 Om^2).
      real(wp) :: ustar_min = 1.0e-8_wp
         !! Floor on u* (deliberate rdb floor — pure denominator
         !! guard in the TKE decay scale; MOM6 derives ~1e-10).
      real(wp) :: prandtl = 1.0_wp
         !! Kv = prandtl * Kd into the momentum solve.
      logical :: in_eos = .false.
         !! `&ocean_psurf_nml in_eos` (E3): start this scheme's column
         !! pressure stack at `multilayer_state_t%p_top` (the ice-shelf
         !! load + surface pressure, Pa) instead of at 0 Pa.
         !!
         !! The stack has TWO consumers here and the knob moves BOTH,
         !! deliberately: the IN-SITU pressure argument of
         !! `eos_specvol_derivs`, and the PE weight
         !! `dpe_* = dmass*p_mid*dsv_*`.  They are the same pressure --
         !! the weight is the hydrostatic load the layer's centre of mass
         !! has to lift, and under a floating shelf the ice is part of
         !! that load.  Splitting them would put two pressure conventions
         !! in one column, which is the failure the `p_top` seam contract
         !! exists to prevent.
         !!
         !! Host scalar, assigned by `configure_ocean_epbl` BEFORE
         !! `enter_data`, and passed BY VALUE into the column kernel --
         !! never read through the device-mapped handle.
         !!
         !! `.false.` (default) ⇒ the stack starts at 0 Pa, character for
         !! character the pre-E3 arithmetic.  It is NOT enough that
         !! `p_top` be the zero array: under a cavity `p_top` is the ice
         !! load whether or not `in_eos` is set, so this gate is what
         !! keeps an existing cavity + EPBL run bit-identical.
      logical :: tke_diags = .false.
         !! Compute + store the per-column TKE budget terms (W/m^2).
         !! The column ledger closes to round-off — see the design
         !! doc §7; `test_ocean_epbl` asserts it.

      ! ---- Langmuir turbulence (LF17 wind-only path) ----
      ! One scalar La per column per MLD iteration, computed from u*
      ! and the current MLD guess (Li & Fox-Kemper 2017 statistical
      ! waves; no wave model, no Stokes arrays).  Enhancement applies
      ! to mstar only (Reichl & Li 2019).  NOTE: the LF17 branch uses
      ! NEITHER MOM6's MIN_LANGMUIR nor LA_DEPTH_MIN floors — those
      ! belong to the profile-averaging wave paths only.
      logical :: use_lt = .false.
         !! Master Langmuir switch (default off — bit-identity).
      integer :: lt_scheme = EPBL_LT_RESCALE
         !! Enhancement form: rescale (multiplicative) or additive.
      real(wp) :: lt_enhance_coef = 0.447_wp
         !! Enhancement coefficient (MOM6 LT_ENHANCE_COEF).
      real(wp) :: lt_enhance_exp = -1.33_wp
         !! La exponent (MOM6 LT_ENHANCE_EXP).
      real(wp) :: lt_max_enhance = 5.0_wp
         !! Cap on the multiplicative enhancement factor.
      real(wp) :: la_frac_hbl = 0.04_wp
         !! Stokes surface-layer-average depth as a fraction of the
         !! BLD (MOM6 LA_DEPTH_RATIO).
      real(wp) :: lt_lac1 = -0.87_wp
      real(wp) :: lt_lac2 = 0.0_wp
      real(wp) :: lt_lac3 = 0.0_wp
      real(wp) :: lt_lac4 = 0.95_wp
      real(wp) :: lt_lac5 = 0.95_wp
         !! Stability-modified Langmuir number coefficients (MOM6
         !! LT_MOD_LAC1..5: MLD/Ekman, MLD/Obukhov stable, unstable,
         !! Ekman/Obukhov stable, unstable).  All-zero => La_mod = La.

      ! ---- EOS hookup (shared handle from the eos slot) ----
      ! The energy weights use the SAME EOS the dyn-core runs.  We
      ! carry a value copy of the flat-POD `eos_t` (no
      ! allocatable) set once at configure from `ocean_state%eos`;
      ! it maps onto the device with the parent `this` for free and
      ! the per-column `eos_specvol_derivs(this%eos, ...)` call reads
      ! it from registers.  One source of truth — no private scalar
      ! copies that can drift from the dyn-core EOS.
      type(eos_t) :: eos

      ! ---- Persistent fields ----
      real(wp), allocatable :: f_centre(:, :)
         !! |f| at cell centres (1/s); filled by `set_f_centre`.
      real(wp), allocatable :: mld(:, :)
         !! Converged active-mixing-layer depth (m) from the last
         !! call — the next step's first guess (when
         !! `mld_use_prev_guess`) and the primary diagnostic.
      real(wp), allocatable :: b0(:, :)
         !! Surface buoyancy flux (m^2/s^3, > 0 stabilizing) from the last
         !! solve — the same `B0 = g·rho0·(dSV/dT·q_T + dSV/dS·q_S)` the
         !! mstar scaling uses, persisted so the Bodner (2023) MLE
         !! restratification can form its convective velocity scale
         !! `w*^3 = max(0,-b0)·mld`.  Zero until the first EPBL step.
      real(wp), allocatable :: kd_int(:, :, :)
         !! EPBL diapycnal diffusivity at interfaces (m^2/s),
         !! (nx, ny, nz+1); zero at bed (k=1) and surface (k=nz+1).
      real(wp), allocatable :: la(:, :)
         !! Turbulent Langmuir number from the last call (diagnostic;
         !! 0 where `use_lt` is off or the column is dry).

      ! ---- Per-step TKE budget diagnostics (W/m^2) ----
      ! Overwritten each `epbl_compute` (final MLD iteration).  The
      ! ledger: wind + conv + forcing - mixing - mech_decay -
      ! conv_decay = 0 exactly (forcing <= 0, all others >= 0).
      real(wp), allocatable :: tke_wind(:, :)
      real(wp), allocatable :: tke_conv(:, :)
      real(wp), allocatable :: tke_forcing(:, :)
      real(wp), allocatable :: tke_mixing(:, :)
      real(wp), allocatable :: tke_mech_decay(:, :)
      real(wp), allocatable :: tke_conv_decay(:, :)

      ! ---- Iteration-invariant column workspaces ----
      ! Filled once per call by the column prep sweep; read by every
      ! MLD iteration.  All (nx, ny, nz).  Everything else the sweep
      ! needs is carried as scalars (see the design doc D6) — no
      ! per-column work arrays, no NZ_STACK_MAX locals.
      type(scratch_3d_buffer_t) :: t0
         !! Pre-mixing layer temperature (degC).
      type(scratch_3d_buffer_t) :: s0
         !! Pre-mixing layer salinity (PSU).
      type(scratch_3d_buffer_t) :: dpe_t
         !! d(column PE)/d(T_k): rho0*h_k * p_mid * dSV_dT (J/m^2/degC).
      type(scratch_3d_buffer_t) :: dpe_s
         !! d(column PE)/d(S_k) (J/m^2/PSU).
      type(scratch_3d_buffer_t) :: dcolht_t
         !! d(column height)/d(T_k): rho0*h_k * dSV_dT (m/degC) —
         !! steric sensitivity for the gravity-wave radiation term.
      type(scratch_3d_buffer_t) :: dcolht_s
         !! d(column height)/d(S_k) (m/PSU).
      type(scratch_3d_buffer_t) :: ctke_sw
         !! (PR-21) Per-layer TKE cost (J/m^2, <= 0 for solar heating) of
         !! homogenising the penetrating shortwave absorbed IN that layer.
         !! Filled once per call by the prep sweep from the shared two-band
         !! transmission + the in-layer PE-cost shape `Phi(h/zeta)`; the
         !! surface layer's share folds into `ctke_sfc`, the sub-surface
         !! layers drain the sweep's TKE reservoirs.  Zero (never filled)
         !! unless `epbl_sw_ctke .and. sf%has_sw`.  Same shape / lifetime /
         !! device contract as `dcolht_s`.
      logical :: epbl_sw_ctke = .true.
         !! (PR-21) Charge the EPBL TKE ledger for penetrating shortwave
         !! (`&ocean_thermo_nml epbl_sw_ctke`).  Default `.true.`; inert at
         !! `sw_pen_frac = 0` (⇒ `sf%has_sw = .false.`) ⇒ bit-identical.
   contains
      procedure, non_overridable :: init => ocean_epbl_init
      procedure, non_overridable :: destroy => ocean_epbl_destroy
      procedure, non_overridable :: enter_data => ocean_epbl_enter_data
      procedure, non_overridable :: exit_data => ocean_epbl_exit_data
      procedure, non_overridable :: set_f_centre => ocean_epbl_set_f_centre
      procedure, non_overridable :: bytes => ocean_epbl_bytes
   end type ocean_epbl_t

contains

   ! =================================================================
   ! Lifecycle
   ! =================================================================

   subroutine ocean_epbl_init(this, grid, nz_ml)
      !! Allocate the persistent fields + column workspaces.  Always
      !! allocates (configure runs after init, so `enable` isn't
      !! known yet); the memory cost when off is the same 7-field
      !! footprint the vmix slot already pays.
      class(ocean_epbl_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      allocate (this%f_centre(nx, ny), source=0.0_wp)
      allocate (this%mld(nx, ny), source=0.0_wp)
      allocate (this%b0(nx, ny), source=0.0_wp)
      allocate (this%kd_int(nx, ny, nz + 1), source=0.0_wp)
      allocate (this%la(nx, ny), source=0.0_wp)
      allocate (this%tke_wind(nx, ny), source=0.0_wp)
      allocate (this%tke_conv(nx, ny), source=0.0_wp)
      allocate (this%tke_forcing(nx, ny), source=0.0_wp)
      allocate (this%tke_mixing(nx, ny), source=0.0_wp)
      allocate (this%tke_mech_decay(nx, ny), source=0.0_wp)
      allocate (this%tke_conv_decay(nx, ny), source=0.0_wp)
      call this%t0%init(nx, ny, nz, "epbl_t0")
      call this%s0%init(nx, ny, nz, "epbl_s0")
      call this%dpe_t%init(nx, ny, nz, "epbl_dpe_t")
      call this%dpe_s%init(nx, ny, nz, "epbl_dpe_s")
      call this%dcolht_t%init(nx, ny, nz, "epbl_dcolht_t")
      call this%dcolht_s%init(nx, ny, nz, "epbl_dcolht_s")
      call this%ctke_sw%init(nx, ny, nz, "epbl_ctke_sw")

      this%is_init = .true.
   end subroutine ocean_epbl_init

   subroutine ocean_epbl_destroy(this)
      class(ocean_epbl_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%f_centre)) deallocate (this%f_centre)
      if (allocated(this%mld)) deallocate (this%mld)
      if (allocated(this%b0)) deallocate (this%b0)
      if (allocated(this%kd_int)) deallocate (this%kd_int)
      if (allocated(this%la)) deallocate (this%la)
      if (allocated(this%tke_wind)) deallocate (this%tke_wind)
      if (allocated(this%tke_conv)) deallocate (this%tke_conv)
      if (allocated(this%tke_forcing)) deallocate (this%tke_forcing)
      if (allocated(this%tke_mixing)) deallocate (this%tke_mixing)
      if (allocated(this%tke_mech_decay)) deallocate (this%tke_mech_decay)
      if (allocated(this%tke_conv_decay)) deallocate (this%tke_conv_decay)
      call this%t0%destroy()
      call this%s0%destroy()
      call this%dpe_t%destroy()
      call this%dpe_s%destroy()
      call this%dcolht_t%destroy()
      call this%dcolht_s%destroy()
      call this%ctke_sw%destroy()
   end subroutine ocean_epbl_destroy

   subroutine ocean_epbl_enter_data(this)
      class(ocean_epbl_t), intent(inout) :: this
      select type (this)
      type is (ocean_epbl_t)
         call ocean_epbl_enter_data_impl(this)
      end select
   end subroutine ocean_epbl_enter_data

   subroutine ocean_epbl_enter_data_impl(this)
      type(ocean_epbl_t), intent(inout) :: this
      if (allocated(this%f_centre)) then
         !$acc enter data copyin(this%f_centre)
      end if
      if (allocated(this%mld)) then
         !$acc enter data copyin(this%mld)
      end if
      if (allocated(this%b0)) then
         !$acc enter data copyin(this%b0)
      end if
      if (allocated(this%kd_int)) then
         !$acc enter data copyin(this%kd_int)
      end if
      if (allocated(this%la)) then
         !$acc enter data copyin(this%la)
      end if
      if (allocated(this%tke_wind)) then
         !$acc enter data copyin(this%tke_wind)
      end if
      if (allocated(this%tke_conv)) then
         !$acc enter data copyin(this%tke_conv)
      end if
      if (allocated(this%tke_forcing)) then
         !$acc enter data copyin(this%tke_forcing)
      end if
      if (allocated(this%tke_mixing)) then
         !$acc enter data copyin(this%tke_mixing)
      end if
      if (allocated(this%tke_mech_decay)) then
         !$acc enter data copyin(this%tke_mech_decay)
      end if
      if (allocated(this%tke_conv_decay)) then
         !$acc enter data copyin(this%tke_conv_decay)
      end if
      call scratch_3d_buffer_enter_data_impl(this%t0)
      call scratch_3d_buffer_enter_data_impl(this%s0)
      call scratch_3d_buffer_enter_data_impl(this%dpe_t)
      call scratch_3d_buffer_enter_data_impl(this%dpe_s)
      call scratch_3d_buffer_enter_data_impl(this%dcolht_t)
      call scratch_3d_buffer_enter_data_impl(this%dcolht_s)
      call scratch_3d_buffer_enter_data_impl(this%ctke_sw)
   end subroutine ocean_epbl_enter_data_impl

   subroutine ocean_epbl_exit_data(this)
      class(ocean_epbl_t), intent(inout) :: this
      select type (this)
      type is (ocean_epbl_t)
         call ocean_epbl_exit_data_impl(this)
      end select
   end subroutine ocean_epbl_exit_data

   subroutine ocean_epbl_exit_data_impl(this)
      type(ocean_epbl_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%ctke_sw)
      call scratch_3d_buffer_exit_data_impl(this%dcolht_s)
      call scratch_3d_buffer_exit_data_impl(this%dcolht_t)
      call scratch_3d_buffer_exit_data_impl(this%dpe_s)
      call scratch_3d_buffer_exit_data_impl(this%dpe_t)
      call scratch_3d_buffer_exit_data_impl(this%s0)
      call scratch_3d_buffer_exit_data_impl(this%t0)
      if (allocated(this%tke_conv_decay)) then
         !$acc exit data delete(this%tke_conv_decay)
      end if
      if (allocated(this%tke_mech_decay)) then
         !$acc exit data delete(this%tke_mech_decay)
      end if
      if (allocated(this%tke_mixing)) then
         !$acc exit data delete(this%tke_mixing)
      end if
      if (allocated(this%tke_forcing)) then
         !$acc exit data delete(this%tke_forcing)
      end if
      if (allocated(this%tke_conv)) then
         !$acc exit data delete(this%tke_conv)
      end if
      if (allocated(this%tke_wind)) then
         !$acc exit data delete(this%tke_wind)
      end if
      if (allocated(this%la)) then
         !$acc exit data delete(this%la)
      end if
      if (allocated(this%kd_int)) then
         !$acc exit data delete(this%kd_int)
      end if
      if (allocated(this%mld)) then
         !$acc exit data delete(this%mld)
      end if
      if (allocated(this%b0)) then
         !$acc exit data delete(this%b0)
      end if
      if (allocated(this%f_centre)) then
         !$acc exit data delete(this%f_centre)
      end if
   end subroutine ocean_epbl_exit_data_impl

   subroutine ocean_epbl_set_f_centre(this, grid, f_0, beta, y_ref)
      !! Fill `f_centre` with the beta-plane Coriolis magnitude at
      !! cell centres: |f_0 + beta*(y - y_ref)|.  Mirror of
      !! `coriolis_adv_set_beta_plane` (which fills corners).  Call
      !! after `init`, before `enter_data`.
      class(ocean_epbl_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_0, beta, y_ref
      integer :: i, j, ng
      real(wp) :: y

      ng = grid%nghost
      do j = 1, size(this%f_centre, 2)
         y = (real(j + grid%j_offset_global - ng, wp) - 0.5_wp)*grid%dy
         do i = 1, size(this%f_centre, 1)
            this%f_centre(i, j) = abs(f_0 + beta*(y - y_ref))
         end do
      end do
   end subroutine ocean_epbl_set_f_centre

   ! =================================================================
   ! Parsers (host-side, configure time)
   ! =================================================================

   pure function parse_epbl_mstar_scheme(name) result(tag)
      character(len=*), intent(in) :: name
      integer :: tag
      select case (trim(name))
      case ("constant")
         tag = EPBL_MSTAR_CONSTANT
      case ("om4", "OM4")
         tag = EPBL_MSTAR_OM4
      case ("rh18", "RH18", "reichl_h18")
         tag = EPBL_MSTAR_RH18
      case default
         tag = -1
      end select
   end function parse_epbl_mstar_scheme

   pure function parse_epbl_vstar_scheme(name) result(tag)
      character(len=*), intent(in) :: name
      integer :: tag
      select case (trim(name))
      case ("cube_root", "cube_root_tke")
         tag = EPBL_VSTAR_CUBE_ROOT
      case ("rh18", "RH18", "reichl_h18")
         tag = EPBL_VSTAR_RH18
      case default
         tag = -1
      end select
   end function parse_epbl_vstar_scheme

   pure function parse_epbl_combine(name) result(tag)
      character(len=*), intent(in) :: name
      integer :: tag
      select case (trim(name))
      case ("add", "additive")
         tag = EPBL_COMBINE_ADD
      case ("max")
         tag = EPBL_COMBINE_MAX
      case default
         tag = -1
      end select
   end function parse_epbl_combine

   pure function parse_epbl_lt_scheme(name) result(tag)
      character(len=*), intent(in) :: name
      integer :: tag
      select case (trim(name))
      case ("rescale")
         tag = EPBL_LT_RESCALE
      case ("additive", "add")
         tag = EPBL_LT_ADDITIVE
      case default
         tag = -1
      end select
   end function parse_epbl_lt_scheme

   ! =================================================================
   ! Pure helpers (device-callable; unit-tested directly)
   ! =================================================================

   pure subroutine epbl_find_mstar(scheme, mstar_const, mstar_cap, &
                                   mstar_coef1, c_ek, mstar_conv_adj, &
                                   cn1, cn2, cn3, cs1, cs2, &
                                   b0, ustar, bld, absf, mstar)
      !! mstar = (mechanical TKE available for entrainment) / u*^3.
      !! Schemes: constant; OM4 Ekman/Obukhov balance; RH18 fits.
      !! All followed by the optional convective reduction
      !! (mstar_conv_adj in [0,1]; the u*=0 corner multiplies by
      !! (1 - adj), matching the reference behaviour).
      !$acc routine seq
      integer, intent(in) :: scheme
      real(wp), intent(in) :: mstar_const, mstar_cap
      real(wp), intent(in) :: mstar_coef1, c_ek, mstar_conv_adj
      real(wp), intent(in) :: cn1, cn2, cn3, cs1, cs2
      real(wp), intent(in) :: b0
         !! Surface buoyancy flux (m^2/s^3); > 0 stabilizing.
      real(wp), intent(in) :: ustar
         !! Surface friction velocity (m/s), already floored.
      real(wp), intent(in) :: bld
         !! Boundary-layer depth guess (m).
      real(wp), intent(in) :: absf
         !! |Coriolis| (1/s), possibly omega-blended.
      real(wp), intent(out) :: mstar

      real(wp) :: mstar_n, mstar_s, msn_term, absf_floor
      real(wp) :: mscr_t1, mscr_t2

      absf_floor = max(absf, 1.0e-20_wp)
      select case (scheme)
      case (EPBL_MSTAR_OM4)
         mstar_s = mstar_coef1*sqrt(max(0.0_wp, b0)/(ustar*ustar*absf_floor))
         mstar_n = 0.0_wp
         if (ustar > absf_floor*bld) then
            mstar_n = c_ek*log(ustar/(absf_floor*max(bld, 1.0e-10_wp)))
         end if
         mstar = max(mstar_s, min(1.25_wp, mstar_n))
         if (mstar_cap > 0.0_wp) mstar = min(mstar_cap, mstar)
      case (EPBL_MSTAR_RH18)
         msn_term = cn2*exp(cn3*bld*absf/ustar)
         mstar_n = cn1*msn_term/(1.0_wp + msn_term)
         mstar_s = cs1*(max(0.0_wp, b0)**2*bld/(ustar**5*absf_floor))**cs2
         mstar = mstar_n + mstar_s
         if (mstar_cap > 0.0_wp) mstar = min(mstar_cap, mstar)
      case default  ! EPBL_MSTAR_CONSTANT
         mstar = mstar_const
      end select

      ! Convective reduction (no-op at the default adj = 0).
      if (mstar_conv_adj > 0.0_wp) then
         mscr_t1 = -bld*min(0.0_wp, b0)
         mscr_t2 = 2.0_wp*mstar*ustar**3
         if (mscr_t2 > 0.0_wp) then
            mstar = mstar*((1.0_wp - mstar_conv_adj)*mscr_t1 + mscr_t2)/ &
                    (mscr_t1 + mscr_t2)
         else
            mstar = mstar*(1.0_wp - mstar_conv_adj)
         end if
      end if
   end subroutine epbl_find_mstar

   pure function epbl_mixlen_shape(z_depth, mld_guess, translay_scale, &
                                   mixlen_exponent, shaped) result(shape_val)
      !! Mixing-length shape factor in [translay_scale, 1]: 1 at the
      !! surface, decaying to the transition-layer floor at the MLD.
      !! `shaped = .false.` (no MLD iteration) returns 1.
      !$acc routine seq
      real(wp), intent(in) :: z_depth
         !! Unconditional interface depth below the surface (m).
      real(wp), intent(in) :: mld_guess, translay_scale, mixlen_exponent
      logical, intent(in) :: shaped
      real(wp) :: shape_val
      real(wp) :: fr

      if ((.not. shaped) .or. translay_scale >= 1.0_wp .or. &
          translay_scale < 0.0_wp .or. mld_guess <= 0.0_wp) then
         shape_val = 1.0_wp
      else
         fr = max(0.0_wp, (mld_guess - z_depth)/mld_guess)
         if (mixlen_exponent == 2.0_wp) then
            shape_val = translay_scale + (1.0_wp - translay_scale)*fr*fr
         else
            shape_val = translay_scale + (1.0_wp - translay_scale)*fr**mixlen_exponent
         end if
      end if
   end function epbl_mixlen_shape

   pure subroutine epbl_lf17_wave_state(ustar_w, rho_ocn, u10, ustokes, kphil)
      !! LF17 statistical wave state from the water-side friction
      !! velocity alone: COARE 3.5 fixed-point inversion u* -> U10
      !! (Edson et al. 2013), then the Pierson-Moskowitz-based
      !! surface Stokes drift and Phillips peak wavenumber (Li &
      !! Fox-Kemper 2017; Breivik et al. 2016).  BLD-independent —
      !! call once per column, outside the MLD iteration.
      !$acc routine seq
      real(wp), intent(in) :: ustar_w
         !! Water-side u* (m/s), > 0 (caller floors it).
      real(wp), intent(in) :: rho_ocn
         !! Seawater reference density (kg/m^3).
      real(wp), intent(out) :: u10
         !! 10-m wind speed (m/s).
      real(wp), intent(out) :: ustokes
         !! Surface Stokes drift (m/s).
      real(wp), intent(out) :: kphil
         !! Phillips-spectrum peak wavenumber (1/m).

      real(wp) :: ustar_air, z0_smooth, z0w, alpha_ch, u10_new
      real(wp) :: hm0, f_mean, vstokes
      integer :: itt
      logical :: converged

      ustar_air = ustar_w*sqrt(rho_ocn/LT_RHO_AIR)
      z0_smooth = 0.11_wp*LT_NU_AIR/ustar_air
      u10 = ustar_air*sqrt(1000.0_wp)
      converged = .false.
      do itt = 1, 20
         alpha_ch = min(LT_CHARNOCK_MIN, LT_CHARNOCK_SLOPE*u10 + LT_CHARNOCK_ICPT)
         ! alpha can go negative below U10 ~ 3 m/s; clamp z0 positive
         ! before the log.
         z0w = max(z0_smooth + alpha_ch*ustar_air**2/GRAVITY, 1.0e-10_wp)
         u10_new = ustar_air*log(10.0_wp/z0w)/LT_VONKAR_WAVES
         if (abs(u10_new - u10) <= 1.0e-3_wp*u10) then
            u10 = u10_new
            converged = .true.
            exit
         end if
         u10 = u10_new
      end do
      if (.not. converged) u10 = 25.82_wp*ustar_air

      ustokes = LT_US_TO_U10*u10
      hm0 = LT_SWH_FROM_U10SQ*u10*u10
      f_mean = LT_FM_INTO_FP*0.877_wp*GRAVITY/ &
               (2.0_wp*LT_PI*LT_U19P5_TO_U10*u10)
      vstokes = 0.125_wp*LT_PI*LT_R_LOSS*f_mean*hm0*hm0
      kphil = 0.176_wp*ustokes/max(vstokes, 1.0e-30_wp)
   end subroutine epbl_lf17_wave_state

   pure function epbl_lf17_la(ustar_w, zsl, ustokes, kphil) result(la)
      !! Turbulent Langmuir number La = sqrt(u*/u_s_SL): the Stokes
      !! drift averaged over the surface layer of thickness `zsl`
      !! under the Phillips spectrum, in the singularity-safe form
      !! (Breivik et al. 2016 with Webb & Fox-Kemper 2015 directional
      !! spreading).  No MIN_LANGMUIR / LA_DEPTH_MIN floors — those
      !! belong to MOM6's profile-averaging wave paths, not LF17.
      !$acc routine seq
      real(wp), intent(in) :: ustar_w
         !! Water-side u* (m/s).
      real(wp), intent(in) :: zsl
         !! Surface-layer average depth (m) = la_frac_hbl * BLD.
      real(wp), intent(in) :: ustokes, kphil
         !! From `epbl_lf17_wave_state`.
      real(wp) :: la

      real(wp) :: xkz, root2kz, r1, r3, r5, ustokes_sl

      xkz = kphil*zsl
      root2kz = sqrt(2.0_wp*xkz)
      r1 = (0.302_wp - 1.68_wp*xkz)*one_m_exp_x(2.0_wp*xkz)
      r3 = (0.1264_wp + 0.64_wp*xkz)*one_m_exp_x(5.12_wp*xkz)
      if (root2kz > 1.0e-3_wp) then
         r5 = sqrt(LT_PI)*(root2kz*(-0.84_wp*erfc(root2kz) + &
                                    0.2_wp*erfc(1.6_wp*root2kz)) + &
                           0.1182_wp*(erfc(1.6_wp*root2kz) - erfc(root2kz))/root2kz)
      else
         r5 = -0.64_wp*sqrt(LT_PI)*root2kz + &
              (-0.14184_wp + 1.0839648_wp*root2kz*root2kz)
      end if
      ! Floor is a pure numerical guard (UStokes_sl -> 0+ for very
      ! deep surface layers); La then goes huge and the enhancement
      ! vanishes, which is the correct limit.
      ustokes_sl = max(ustokes*(0.715_wp + r1 + r3 + r5), 1.0e-10_wp)
      la = sqrt(ustar_w/ustokes_sl)
   end function epbl_lf17_la

   pure function one_m_exp_x(x) result(f)
      !! (1 - exp(-x)) / x, Taylor-safe at small x.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: f
      if (x < 1.0e-4_wp) then
         f = 1.0_wp - x*(0.5_wp - x/6.0_wp)
      else
         f = (1.0_wp - exp(-x))/x
      end if
   end function one_m_exp_x

   pure subroutine epbl_lt_enhance(scheme, coef, expo, max_enh, vonkar, &
                                   lac1, lac2, lac3, lac4, lac5, &
                                   la, b0, ustar, bld, absf, mstar)
      !! Apply the Langmuir enhancement to mstar (Reichl & Li 2019;
      !! Li et al. 2016 stability modification).  The modified
      !! Langmuir number folds the boundary-layer stability regime in
      !! via Ekman / Obukhov / MLD length-scale ratios, split by the
      !! sign of the surface buoyancy flux; all-zero lac coefficients
      !! give La_mod = La exactly.
      !$acc routine seq
      integer, intent(in) :: scheme
      real(wp), intent(in) :: coef, expo, max_enh, vonkar
      real(wp), intent(in) :: lac1, lac2, lac3, lac4, lac5
      real(wp), intent(in) :: la
         !! Raw turbulent Langmuir number (> 0).
      real(wp), intent(in) :: b0
         !! Surface buoyancy flux (m^2/s^3); > 0 stabilizing.
      real(wp), intent(in) :: ustar, bld, absf
      real(wp), intent(inout) :: mstar

      real(wp) :: mld_ek, ek_ob, mld_ob, la_mod, enh

      mld_ek = bld*absf/ustar
      ek_ob = abs(b0*vonkar)/(max(absf, 1.0e-20_wp)*ustar*ustar)
      mld_ob = abs(bld*b0*vonkar)/ustar**3
      if (b0 >= 0.0_wp) then
         la_mod = la*((1.0_wp + max(-0.5_wp, lac1*mld_ek)) + &
                      lac4*ek_ob + lac2*mld_ob)
      else
         la_mod = la*((1.0_wp + max(-0.5_wp, lac1*mld_ek)) + &
                      lac5*ek_ob + lac3*mld_ob)
      end if
      la_mod = max(la_mod, 1.0e-10_wp)

      if (scheme == EPBL_LT_ADDITIVE) then
         mstar = mstar + coef*la_mod**expo
      else
         enh = min(max_enh, 1.0_wp + coef*la_mod**expo)
         mstar = mstar*enh
      end if
   end subroutine epbl_lt_enhance

   ! =================================================================
   ! Main compute
   ! =================================================================

   pure subroutine epbl_compute(grid, this, ms, ss, dt, sf)
      !! Run EPBL over the domain: fill `this%kd_int` (interface
      !! diffusivity) and `this%mld`.  Call at thermo cadence with
      !! the thermo dt.  Outer shim: dereferences the tracer-registry
      !! hTr arrays + the 2D Q_heat/Q_salt forcing fields on the host
      !! (array-of-DT and allocatable indirection blocks NVHPC device
      !! codegen), then forwards to the column kernel which reads
      !! q_T_kin(i,j) / q_S_kin(i,j) per column inside the DC.
      !!
      !! `sf` (surface-flux slot) is REQUIRED — EPBL computes B_0 from
      !! the per-column kinematic fluxes.  The dispatch in
      !! `vmix_apply_in_stage` already errors before reaching here when
      !! `sf` is absent and EPBL is active.
      type(hgrid_t), intent(in) :: grid
      type(ocean_epbl_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_surface_stress_t), intent(in) :: ss
      real(wp), intent(in) :: dt
      type(ocean_surface_flux_t), intent(in) :: sf

      integer :: nx, ny
      logical :: sw_ctke_active

      if (ms%idx_temperature <= 0 .or. ms%idx_salinity <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total

      ! Kinematic scale factors (column-invariant multipliers):
      !   q_T_kin(i,j) = Q_heat(i,j) / (rho0 · cp)  [degC·m/s]
      !   q_S_kin(i,j) = Q_salt(i,j) / rho0           [PSU·m/s]
      ! The column kernel reads these at (i,j) inside the DC so that
      ! spatially-varying forcing (Area-A3/A4) is handled correctly.
      !
      ! (PR-21) Penetrating-SW TKE ledger: active only when the ledger is
      ! enabled AND penetrating SW is on.  The irradiance source is
      ! selected HOST-SIDE (`sf%sw_from_qsw`); the conditionally-allocated
      ! `q_sw` is passed only on the guarded branch (validate_config
      ! forces enable_components when sw_source="q_sw", making it total).
      ! When inactive, the source array is irrelevant (kernel never reads
      ! it), so the legacy `Q_heat` argument is passed on both branches.
      sw_ctke_active = this%epbl_sw_ctke .and. sf%has_sw
      if (sw_ctke_active .and. sf%sw_from_qsw) then
         call epbl_column_kernel(grid, this, ms, &
                                 ms%tracers(ms%idx_temperature)%hTr, &
                                 ms%tracers(ms%idx_salinity)%hTr, &
                                 ss, dt, &
                                 sf%Q_heat, 1.0_wp/(this%rho0*sf%cp), &
                                 sf%Q_salt, 1.0_wp/this%rho0, &
                                 sf%q_sw, sw_ctke_active, sf%sw_pen_frac, &
                                 sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                 ms%wet_mask, this%in_eos, nx, ny)
      else
         call epbl_column_kernel(grid, this, ms, &
                                 ms%tracers(ms%idx_temperature)%hTr, &
                                 ms%tracers(ms%idx_salinity)%hTr, &
                                 ss, dt, &
                                 sf%Q_heat, 1.0_wp/(this%rho0*sf%cp), &
                                 sf%Q_salt, 1.0_wp/this%rho0, &
                                 sf%Q_heat, sw_ctke_active, sf%sw_pen_frac, &
                                 sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                 ms%wet_mask, this%in_eos, nx, ny)
      end if
   end subroutine epbl_compute

   pure subroutine epbl_column_kernel(grid, this, ms, hT, hS, ss, dt, &
                                      Q_heat_field, inv_rho0_cp, &
                                      Q_salt_field, inv_rho0, &
                                      sw_src_field, sw_ctke_active, sw_pen_frac, &
                                      sw_R, sw_zeta1, sw_zeta2, wet_mask_field, &
                                      p_top_in_eos, nx_arg, ny_arg)
      !! Per-column EPBL solve.  One `do concurrent (j, i)` with the
      !! serial work in k inside (j -> i -> k ordering); ALL sweep
      !! state is carried in scalars (design doc D6) — the only
      !! column arrays are the six iteration-invariant workspaces
      !! filled by the prep sweep.
      !!
      !! Kinematic surface fluxes are read per-column inside the DC:
      !!   q_t_kin = Q_heat_field(i,j) / (rho0 · cp)  [degC·m/s]
      !!   q_s_kin = Q_salt_field(i,j) / rho0           [PSU·m/s]
      !! For the constant-fill default both fields are uniform, giving
      !! arithmetic identical to the former scalar-broadcast path.
      !!
      !! Algorithm:
      !!   prep:  T0, S0 and the pressure-weighted PE / steric
      !!          sensitivities per layer (downward pressure sum).
      !!   outer: MLD root-find (false position / bisection) —
      !!          mstar and the mixing-length shape depend on MLD.
      !!   sweep: interfaces Ki = nz..2 downward.  Decay mech TKE
      !!          across the layer above; rotation-reduce the
      !!          convective reservoir; closed-form energy solve for
      !!          the largest affordable Kd (the gravity-wave
      !!          column-height correction folds into PEc_core);
      !!          advance the embedded tridiagonal forward
      !!          elimination (hp_a / dX_to_dPE_a / Th_a recursions).
      type(hgrid_t), intent(in) :: grid
      type(ocean_epbl_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      ! assumed-shape-ok: tracer registry outer-shim — caller host-dereferences
      ! ms%tracers(idx)%hTr before passing; size varies per tracer slot
      ! (see CLAUDE.md "outer-shim + flat-impl" pattern); thermo cadence.
      real(wp), intent(in) :: hT(:, :, :)
         !! Temperature tracer hTr (degC*m), host-dereferenced.
      real(wp), intent(in) :: hS(:, :, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
         !! Salinity tracer hTr (PSU*m), host-dereferenced.
      type(ocean_surface_stress_t), intent(in) :: ss
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: Q_heat_field(:, :)   ! assumed-shape-ok: outer-shim pass; thermo cadence
         !! 2D heat-flux field (W/m²), host-dereferenced from sf%Q_heat.
      real(wp), intent(in) :: inv_rho0_cp
         !! Precomputed 1/(rho0·cp) multiplier.
      real(wp), intent(in) :: Q_salt_field(:, :)   ! assumed-shape-ok: outer-shim pass; thermo cadence
         !! 2D salt-flux field (kg/m²/s), host-dereferenced from sf%Q_salt.
      real(wp), intent(in) :: inv_rho0
         !! Precomputed 1/rho0 multiplier.
      real(wp), intent(in) :: sw_src_field(:, :)   ! assumed-shape-ok: outer-shim pass; thermo cadence
         !! (PR-21) 2D irradiance source (W/m², >= 0 on the q_sw path),
         !! host-selected from sf%Q_heat or sf%q_sw.  Read only when
         !! `sw_ctke_active`; otherwise the legacy sf%Q_heat is passed and
         !! ignored.
      logical, intent(in) :: sw_ctke_active
         !! (PR-21) Charge the TKE ledger for penetrating SW (host-side
         !! `epbl_sw_ctke .and. sf%has_sw`).  False ⇒ the surface
         !! energetics + sweep run the unmodified legacy lines.
      real(wp), intent(in) :: sw_pen_frac, sw_R, sw_zeta1, sw_zeta2
         !! (PR-21) Two-band SW parameters (from sf), by value.
      real(wp), intent(in) :: wet_mask_field(:, :)   ! assumed-shape-ok: outer-shim pass; thermo cadence
         !! (PR-21) Wet mask — mirrors the deposition kernel's I0 gate so
         !! the ledger and the tracer field account the same heat.
      logical, intent(in) :: p_top_in_eos
         !! (E3) Seed the column pressure stack at `ms%p_top(i,j)` rather
         !! than at 0 Pa — `&ocean_psurf_nml in_eos`, by value.  `.false.`
         !! ⇒ the pre-E3 arithmetic, character for character.
      integer, intent(in) :: nx_arg, ny_arg
         !! Grid extents — used to bounds-check Q_* indexing.

      integer :: i, j, nx, ny, nz
      ! prep locals
      integer :: k
      real(wp) :: q_t_kin, q_s_kin
      real(wp) :: hk, hk_eff, inv_h, t0k, s0k, dmass, dpres, p_mid, pres
      real(wp) :: dsv_dt_k, dsv_ds_k, dsv_dt_sfc, dsv_ds_sfc, h_sum
      ! forcing / environment locals
      real(wp) :: ustar, absf, idecay, mech_in
      real(wp) :: b0, ctke_sfc
      ! MLD iteration locals
      integer :: obl_it, n_its
      real(wp) :: min_mld, max_mld, mld_guess, mld_found
      real(wp) :: dmld_min, dmld_max
      logical :: have_min, have_max
      real(wp) :: mstar_val, mech_tke, conv_perel, forcing_clip
      ! sweep carried scalars
      integer :: ki, ka, kb
      real(wp) :: htot, z_int, pres_int, mld_output
      logical :: sfc_connected, sfc_disconnect
      real(wp) :: hp_a, dpe_t_a, dpe_s_a, dch_t_a, dch_s_a
      real(wp) :: th_a, sh_a, te_lag, se_lag, kddt_prev, kddt_cur
      real(wp) :: exp_kh, nstar_fc, tot_tke, dt_h, tke_here, vstar
      real(wp) :: h_ka, h_kb, hb_hs, shape_fn, hbs, mixlen, kd_g0
      real(wp) :: hp_b, th_b, sh_b, hps, bdt1, dt_c, ds_c
      real(wp) :: pec_core, colht_core, dkddt, pe_max
      real(wp) :: kd_val, tke_used, frac_bl, pe_g0, dpe_conv
      real(wp) :: b1, c1, r_reduc, te_new, se_new, surf_scale
      ! Langmuir (LF17) wave state + Langmuir number
      real(wp) :: lt_u10, lt_ustokes, lt_kphil, la_val
      ! TKE budget ledger (final iteration's values survive)
      real(wp) :: d_wind, d_conv, d_forcing, d_mixing, d_mdecay, d_cdecay
      ! (PR-21) penetrating-SW TKE ledger locals
      real(wp) :: i0_col, d_top_sw, d_bot_sw, heat_sw1, heat_sw2, q_nonpen_kin
      real(wp) :: ctke_sw_kb, r_sw

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      do concurrent(j=1:ny, i=1:nx) &
         local(k, hk, hk_eff, inv_h, t0k, s0k, dmass, dpres, p_mid, pres, &
               dsv_dt_k, dsv_ds_k, dsv_dt_sfc, dsv_ds_sfc, h_sum, &
               ustar, absf, idecay, mech_in, b0, ctke_sfc, &
               obl_it, n_its, min_mld, max_mld, mld_guess, mld_found, &
               dmld_min, dmld_max, have_min, have_max, &
               mstar_val, mech_tke, conv_perel, forcing_clip, &
               ki, ka, kb, htot, z_int, pres_int, mld_output, &
               sfc_connected, sfc_disconnect, &
               hp_a, dpe_t_a, dpe_s_a, dch_t_a, dch_s_a, &
               th_a, sh_a, te_lag, se_lag, kddt_prev, kddt_cur, &
               exp_kh, nstar_fc, tot_tke, dt_h, tke_here, vstar, &
               h_ka, h_kb, hb_hs, shape_fn, hbs, mixlen, kd_g0, &
               hp_b, th_b, sh_b, hps, bdt1, dt_c, ds_c, &
               pec_core, colht_core, dkddt, pe_max, &
               kd_val, tke_used, frac_bl, pe_g0, dpe_conv, &
               b1, c1, r_reduc, te_new, se_new, surf_scale, &
               lt_u10, lt_ustokes, lt_kphil, la_val, &
               d_wind, d_conv, d_forcing, d_mixing, d_mdecay, d_cdecay, &
               q_t_kin, q_s_kin, &
               i0_col, d_top_sw, d_bot_sw, heat_sw1, heat_sw2, q_nonpen_kin, &
               ctke_sw_kb, r_sw)

         ! ---- Column prep: T0/S0 + PE/steric weights (downward) ----
         ! (E3) The top of this column.  0 Pa at a free surface; the
         ! ice-shelf load + surface pressure `ms%p_top(i,j)` under a lid,
         ! when `&ocean_psurf_nml in_eos` is set.  Everything below
         ! accumulates `g*rho_0*h` downward from here, so `p_mid` is a
         ! true per-layer hydrostatic pressure either way — which is the
         ! test the `p_top` seam contract applies to a joining builder.
         ! It reaches BOTH consumers of the stack (the in-situ EOS
         ! argument and the PE weight `dmass*p_mid*dsv`), because they
         ! are the same pressure: the load a layer's centre of mass has
         ! to lift, ice included.
         pres = 0.0_wp
         if (p_top_in_eos) pres = ms%p_top(i, j)
         h_sum = 0.0_wp
         dsv_dt_sfc = 0.0_wp
         dsv_ds_sfc = 0.0_wp
         ! (PR-21) penetrating-SW column irradiance + running depth of the
         ! layer TOP below the free surface (0 at k=nz, grows downward to
         ! the opaque bed at k=1).  `i0_col` mirrors the deposition
         ! kernel's `I0 = sw_pen_frac·sw_src·wet_mask` exactly.
         if (sw_ctke_active) then
            i0_col = sw_pen_frac*sw_src_field(i, j)*wet_mask_field(i, j)
         else
            i0_col = 0.0_wp
         end if
         d_top_sw = 0.0_wp
         do k = nz, 1, -1
            hk = ms%h_layer(i, j, k)
            if (hk > 0.0_wp) then
               inv_h = 1.0_wp/(hk + H_NEGLECT)
               t0k = hT(i, j, k)*inv_h
               s0k = hS(i, j, k)*inv_h
            else
               t0k = 0.0_wp
               s0k = 0.0_wp
            end if
            dmass = this%rho0*hk
            dpres = GRAVITY*dmass
            p_mid = pres + 0.5_wp*dpres
            call eos_specvol_derivs(this%eos, t0k, s0k, p_mid, dsv_dt_k, dsv_ds_k)
            this%t0%data(i, j, k) = t0k
            this%s0%data(i, j, k) = s0k
            this%dpe_t%data(i, j, k) = dmass*p_mid*dsv_dt_k
            this%dpe_s%data(i, j, k) = dmass*p_mid*dsv_ds_k
            this%dcolht_t%data(i, j, k) = dmass*dsv_dt_k
            this%dcolht_s%data(i, j, k) = dmass*dsv_ds_k
            ! (PR-21) Per-layer penetrating-SW TKE cost.  Absorbed heat
            ! per band (degC·m) is the two-band difference form with an
            ! opaque bed at k=1 (identical to the deposition kernel, so
            ! Σ_k Σ_n heat = I0·dt/(rho0·cp) exactly).  The in-layer PE
            ! cost of homogenising that exponentially-distributed heating
            ! is `Phi(h/zeta) <= 1` of the skin cost; the skin identity
            ! `rho0²·h·dsv_dt ≡ rho0·dcolht_t` makes this division-free.
            if (sw_ctke_active) then
               d_bot_sw = d_top_sw + hk
               if (k > 1) then
                  heat_sw1 = i0_col*sw_R*inv_rho0_cp*dt* &
                             (exp(-d_top_sw/sw_zeta1) - exp(-d_bot_sw/sw_zeta1))
                  heat_sw2 = i0_col*(1.0_wp - sw_R)*inv_rho0_cp*dt* &
                             (exp(-d_top_sw/sw_zeta2) - exp(-d_bot_sw/sw_zeta2))
               else
                  ! opaque bed: absorb the whole remaining irradiance.
                  heat_sw1 = i0_col*sw_R*inv_rho0_cp*dt*exp(-d_top_sw/sw_zeta1)
                  heat_sw2 = i0_col*(1.0_wp - sw_R)*inv_rho0_cp*dt*exp(-d_top_sw/sw_zeta2)
               end if
               this%ctke_sw%data(i, j, k) = &
                  -0.5_wp*GRAVITY*this%rho0*this%dcolht_t%data(i, j, k)* &
                  (heat_sw1*sw_pe_cost_shape(hk/sw_zeta1) + &
                   heat_sw2*sw_pe_cost_shape(hk/sw_zeta2))
               d_top_sw = d_bot_sw
            end if
            if (k == nz) then
               dsv_dt_sfc = dsv_dt_k
               dsv_ds_sfc = dsv_ds_k
            end if
            pres = pres + dpres
            h_sum = h_sum + hk
         end do
         h_sum = h_sum + H_NEGLECT

         ! ---- Surface forcing energetics ----
         ! Kinematic fluxes at this column:
         !   q_T_kin = Q_heat(i,j) / (rho0·cp)  [degC·m/s]
         !   q_S_kin = Q_salt(i,j) / rho0         [PSU·m/s]
         q_t_kin = inv_rho0_cp*Q_heat_field(i, j)
         q_s_kin = inv_rho0*Q_salt_field(i, j)
         ! B0 = g rho0 (dSV_dT q_T + dSV_dS q_S); > 0 stabilizing.
         ! B0 drives mstar / Langmuir (surface-buoyancy scaling); it uses
         ! the FULL net heat flux, unchanged by this PR (MOM6 does the
         ! same — penetrating SW enters ONLY through the per-layer TKE
         ! ledger, not the mstar buoyancy scale).
         b0 = GRAVITY*this%rho0*(dsv_dt_sfc*q_t_kin + dsv_ds_sfc*q_s_kin)
         this%b0(i, j) = b0   ! persist for the Bodner MLE convective velocity scale
         ! cTKE_sfc: PE released (> 0) / required (< 0) to homogenize
         ! the freshly applied SKIN fluxes through the surface layer.
         ! (PR-21) When penetrating SW is charged to the ledger, the skin
         ! carries only the NON-penetrating heat `q_nonpen = q_heat - I0`;
         ! the SW absorbed within the surface layer is the separately
         ! computed `ctke_sw(nz)` (Phi-weighted, <= the all-skin charge).
         ! I0 = 0 ⇒ q_nonpen = q_heat and ctke_sw(nz) = 0 ⇒ the else
         ! branch is character-for-character the legacy line.
         if (sw_ctke_active) then
            q_nonpen_kin = inv_rho0_cp*(Q_heat_field(i, j) - i0_col)
            ctke_sfc = -0.5_wp*GRAVITY*this%rho0**2*ms%h_layer(i, j, nz)* &
                       (q_nonpen_kin*dt*dsv_dt_sfc + q_s_kin*dt*dsv_ds_sfc) &
                       + this%ctke_sw%data(i, j, nz)
         else
            ctke_sfc = -0.5_wp*GRAVITY*this%rho0**2*ms%h_layer(i, j, nz)* &
                       (q_t_kin*dt*dsv_dt_sfc + q_s_kin*dt*dsv_ds_sfc)
         end if

         ! PR-12 dedup: |tau| at cell centres is a shared field
         ! (ocean_surface_stress_set_derived) — the inner sqrt(tau_xc^2 +
         ! tau_yc^2) below IS stress_mag, computed with the identical FP
         ! op order, so this substitution is bit-identical (§7.5).
         ! Phase 4b: `stress_shelf` adds the ICE-SHELF base stress, which
         ! is not in `tau` (the cover mask zeroes the wind there).  The
         ! supports are disjoint, the field is always allocated, and it
         ! is the zero array without a cavity — `x + 0.0` is `x`.  See
         ! the contract in `rdb_ocean_surface_stress`.
         ustar = max(sqrt((ss%stress_mag(i, j) + ss%stress_shelf(i, j))/ &
                          this%rho0), this%ustar_min)
         absf = sqrt((1.0_wp - this%omega_frac)*this%f_centre(i, j)**2 + &
                     this%omega_frac*4.0_wp*this%omega**2)
         idecay = this%tke_decay*absf/ustar
         mech_in = dt*this%rho0*ustar**3

         ! LF17 wave state is BLD-independent: compute once per
         ! column; only the surface-layer average inside the MLD
         ! iteration depends on the guess.
         la_val = 0.0_wp
         lt_u10 = 0.0_wp
         lt_ustokes = 0.0_wp
         lt_kphil = 0.0_wp
         if (this%use_lt) then
            call epbl_lf17_wave_state(ustar, this%rho0, lt_u10, lt_ustokes, lt_kphil)
         end if

         d_wind = 0.0_wp
         d_conv = 0.0_wp
         d_forcing = 0.0_wp
         d_mixing = 0.0_wp
         d_mdecay = 0.0_wp
         d_cdecay = 0.0_wp
         mld_found = 0.0_wp

         if (ms%wet_mask(i, j) <= 0.0_wp .or. h_sum <= 2.0_wp*H_NEGLECT) then
            ! Dry / land column: no mixing.
            do k = 1, nz + 1
               this%kd_int(i, j, k) = 0.0_wp
            end do
         else

            ! ---- Outer MLD iteration ----
            min_mld = 0.0_wp
            max_mld = h_sum
            dmld_min = 0.0_wp
            dmld_max = 0.0_wp
            have_min = .false.
            have_max = .false.
            mld_guess = 0.5_wp*(min_mld + max_mld)
            if (this%mld_use_prev_guess .and. this%mld(i, j) > 0.0_wp) then
               mld_guess = min(this%mld(i, j), max_mld)
            end if
            n_its = 1
            if (this%mld_iteration) n_its = this%mld_max_its

            do obl_it = 1, n_its
               ! Budget ledger restarts each iteration; the last
               ! iteration's values are the ones reported.
               d_wind = 0.0_wp
               d_conv = 0.0_wp
               d_forcing = 0.0_wp
               d_mixing = 0.0_wp
               d_mdecay = 0.0_wp
               d_cdecay = 0.0_wp

               ! (A) mstar at the current MLD guess.
               call epbl_find_mstar(this%mstar_scheme, this%mstar_const, &
                                    this%mstar_cap, this%mstar_coef1, &
                                    this%c_ek, this%mstar_conv_adj, &
                                    this%rh18_cn1, this%rh18_cn2, this%rh18_cn3, &
                                    this%rh18_cs1, this%rh18_cs2, &
                                    b0, ustar, mld_guess, absf, mstar_val)
               if (this%use_lt) then
                  la_val = epbl_lf17_la(ustar, this%la_frac_hbl*mld_guess, &
                                        lt_ustokes, lt_kphil)
                  call epbl_lt_enhance(this%lt_scheme, this%lt_enhance_coef, &
                                       this%lt_enhance_exp, this%lt_max_enhance, &
                                       this%von_karman, &
                                       this%lt_lac1, this%lt_lac2, this%lt_lac3, &
                                       this%lt_lac4, this%lt_lac5, &
                                       la_val, b0, ustar, mld_guess, absf, mstar_val)
               end if
               mech_tke = mstar_val*mech_in
               d_wind = mech_tke

               ! (B) seed the reservoirs from the surface forcing.
               if (ctke_sfc <= 0.0_wp) then
                  forcing_clip = max(ctke_sfc, -mech_tke)
                  mech_tke = mech_tke + forcing_clip
                  conv_perel = 0.0_wp
                  d_forcing = forcing_clip
               else
                  conv_perel = ctke_sfc
                  d_conv = this%nstar*ctke_sfc
               end if

               ! (C) sweep initialization at the surface layer.
               h_ka = ms%h_layer(i, j, nz) + H_NEGLECT
               hp_a = h_ka
               dpe_t_a = this%dpe_t%data(i, j, nz)
               dpe_s_a = this%dpe_s%data(i, j, nz)
               dch_t_a = this%dcolht_t%data(i, j, nz)
               dch_s_a = this%dcolht_s%data(i, j, nz)
               th_a = h_ka*this%t0%data(i, j, nz)
               sh_a = h_ka*this%s0%data(i, j, nz)
               te_lag = 0.0_wp
               se_lag = 0.0_wp
               kddt_prev = 0.0_wp
               htot = ms%h_layer(i, j, nz)
               z_int = ms%h_layer(i, j, nz)
               pres_int = GRAVITY*this%rho0*ms%h_layer(i, j, nz)
               mld_output = ms%h_layer(i, j, nz)
               sfc_connected = .true.
               this%kd_int(i, j, nz + 1) = 0.0_wp

               ! (D) downward sweep over interfaces Ki = nz .. 2.
               ! Layer above the interface: ka = Ki; below: kb = Ki-1.
               do ki = nz, 2, -1
                  ka = ki
                  kb = ki - 1
                  h_ka = ms%h_layer(i, j, ka) + H_NEGLECT
                  h_kb = ms%h_layer(i, j, kb) + H_NEGLECT
                  sfc_disconnect = .false.

                  ! (1) mechanical TKE decays across the layer above
                  !     (Ekman-scale e-folding; no decay at f = 0).
                  exp_kh = exp(-h_ka*idecay)
                  d_mdecay = d_mdecay + (1.0_wp - exp_kh)*mech_tke
                  mech_tke = mech_tke*exp_kh

                  ! (2) per-layer convective forcing accrual (PR-21).
                  !     The surface layer's cTKE seeds the reservoirs in
                  !     (B); here the sub-surface layer `kb` accrues its
                  !     penetrating-SW cost.  A POSITIVE ctke_sw(kb) (SW
                  !     cooling: only reachable on the legacy net-heat
                  !     path at night) releases convective PE — into
                  !     conv_perel, posted to the convective ledger like
                  !     the surface seed.  A NEGATIVE ctke_sw(kb) (solar
                  !     heating: the physical case) is a TKE SINK and
                  !     drains the reservoirs in (3b), after tot_tke is
                  !     formed.
                  if (sw_ctke_active) then
                     ctke_sw_kb = this%ctke_sw%data(i, j, kb)
                     if (ctke_sw_kb > 0.0_wp) then
                        conv_perel = conv_perel + ctke_sw_kb
                        d_conv = d_conv + this%nstar*ctke_sw_kb
                     end if
                  else
                     ctke_sw_kb = 0.0_wp
                  end if

                  ! (3) rotation-reduced convective efficiency.
                  nstar_fc = this%nstar
                  if (conv_perel > 0.0_wp .and. absf > 0.0_wp) then
                     nstar_fc = this%nstar*conv_perel/ &
                                (conv_perel + 0.2_wp* &
                                 sqrt(0.5_wp*dt*this%rho0*(absf*htot)**3*conv_perel))
                  end if
                  tot_tke = mech_tke + nstar_fc*conv_perel

                  ! (3b) penetrating-SW TKE drain (PR-21).  A negative
                  !      ctke_sw(kb) (solar heating stratifies below the
                  !      interface) must be paid before any mixing here —
                  !      discard that TKE to homogenise the SW through the
                  !      next denser cell.  Mechanical + convective
                  !      reservoirs drain PROPORTIONATELY (Reichl &
                  !      Hallberg 2018).  The consumed energy is booked to
                  !      d_mixing (the PE-raising work the SW stratification
                  !      demands), and the nstar-vs-nstar_fc slop of the
                  !      drained convective reservoir to d_cdecay — mirroring
                  !      the proven (8b) closed-form drain below, so the
                  !      column TKE budget still closes to round-off.  (Note:
                  !      posting to d_mixing OR d_forcing balances the ledger;
                  !      posting to BOTH double-counts.)
                  if (sw_ctke_active) then
                     if (ctke_sw_kb < 0.0_wp) then
                        if (ctke_sw_kb + tot_tke < 0.0_wp) then
                           ! SW cost exhausts all the TKE at this interface.
                           d_mixing = d_mixing + tot_tke
                           d_cdecay = d_cdecay + (this%nstar - nstar_fc)*conv_perel
                           tot_tke = 0.0_wp
                           mech_tke = 0.0_wp
                           conv_perel = 0.0_wp
                        else
                           r_sw = (tot_tke + ctke_sw_kb)/tot_tke
                           d_mixing = d_mixing - ctke_sw_kb
                           d_cdecay = d_cdecay + &
                                      (1.0_wp - r_sw)*(this%nstar - nstar_fc)*conv_perel
                           tot_tke = r_sw*tot_tke
                           mech_tke = r_sw*mech_tke
                           conv_perel = r_sw*conv_perel
                        end if
                     end if
                  end if

                  ! (5) static-stability short-circuit: no energy and
                  !     a stable interface => no mixing here.
                  if (tot_tke <= 0.0_wp .and. &
                      0.0_wp <= (this%dcolht_t%data(i, j, kb) + this%dcolht_t%data(i, j, ka))* &
                      (this%t0%data(i, j, ka) - this%t0%data(i, j, kb)) + &
                      (this%dcolht_s%data(i, j, kb) + this%dcolht_s%data(i, j, ka))* &
                      (this%s0%data(i, j, ka) - this%s0%data(i, j, kb))) then
                     kd_val = 0.0_wp
                     sfc_disconnect = .true.
                  else
                     ! (6) velocity scale, mixing length, first-guess Kd.
                     dt_h = dt/max(0.5_wp*(h_ka + h_kb), 1.0e-15_wp*h_sum)
                     tke_here = mech_tke + this%wstar_ustar_coef*conv_perel
                     hb_hs = (h_sum - z_int)/h_sum
                     shape_fn = epbl_mixlen_shape(z_int, mld_guess, &
                                                  this%translay_scale, &
                                                  this%mixlen_exponent, &
                                                  this%mld_iteration)
                     hbs = min(hb_hs, shape_fn)
                     if (tke_here > 0.0_wp) then
                        if (this%vstar_scheme == EPBL_VSTAR_RH18) then
                           surf_scale = max(0.05_wp, 1.0_wp - htot/mld_guess)
                           vstar = this%vstar_scale_fac*surf_scale* &
                                   (this%vstar_surf_fac*ustar + &
                                    (this%wstar_ustar_coef*conv_perel/ &
                                     (dt*this%rho0))**(1.0_wp/3.0_wp))
                        else
                           vstar = this%vstar_scale_fac* &
                                   (tke_here/(dt*this%rho0))**(1.0_wp/3.0_wp)
                        end if
                        if (this%mld_iteration) then
                           mixlen = max(this%min_mix_len, &
                                        (htot*hbs*vstar)/ &
                                        (this%ekman_scale_coef*absf*htot*hbs + vstar))
                           kd_g0 = vstar*this%von_karman*mixlen
                        else
                           kd_g0 = vstar*this%von_karman*(htot*hbs*vstar)/ &
                                   (this%ekman_scale_coef*absf*htot*hbs + vstar)
                        end if
                     else
                        vstar = 0.0_wp
                        kd_g0 = 0.0_wp
                     end if

                     ! (7) pivot quantities for the layer below.
                     hp_b = h_kb
                     th_b = h_kb*this%t0%data(i, j, kb)
                     sh_b = h_kb*this%s0%data(i, j, kb)

                     ! (8) closed-form energy solve (direct path).
                     hps = hp_a + hp_b
                     bdt1 = hp_a*hp_b
                     dt_c = hp_a*th_b - hp_b*th_a
                     ds_c = hp_a*sh_b - hp_b*sh_a
                     pec_core = hp_b*(dpe_t_a*dt_c + dpe_s_a*ds_c) - &
                                hp_a*(this%dpe_t%data(i, j, kb)*dt_c + &
                                      this%dpe_s%data(i, j, kb)*ds_c)
                     colht_core = hp_b*(dch_t_a*dt_c + dch_s_a*ds_c) - &
                                  hp_a*(this%dcolht_t%data(i, j, kb)*dt_c + &
                                        this%dcolht_s%data(i, j, kb)*ds_c)
                     ! Gravity-wave radiation correction: a shrinking
                     ! column radiates energy that cannot drive mixing.
                     if (colht_core < 0.0_wp) then
                        pec_core = pec_core - pres_int*colht_core
                     end if
                     dkddt = kd_g0*dt_h
                     pe_max = pec_core/(bdt1*hps)

                     if (pe_max < 0.0_wp) then
                        ! (8a) convectively unstable: mixing RELEASES
                        ! PE.  Recompute vstar with the released
                        ! energy included; Kd from the mixing length
                        ! (not energy-limited); bank the release.
                        tke_here = mech_tke + &
                                   this%wstar_ustar_coef*(conv_perel - pe_max)
                        if (tke_here > 0.0_wp) then
                           if (this%vstar_scheme == EPBL_VSTAR_RH18) then
                              surf_scale = max(0.05_wp, 1.0_wp - htot/mld_guess)
                              vstar = this%vstar_scale_fac*surf_scale* &
                                      (this%vstar_surf_fac*ustar + &
                                       (this%wstar_ustar_coef*conv_perel/ &
                                        (dt*this%rho0))**(1.0_wp/3.0_wp))
                           else
                              vstar = this%vstar_scale_fac* &
                                      (tke_here/(dt*this%rho0))**(1.0_wp/3.0_wp)
                           end if
                           if (this%mld_iteration) then
                              mixlen = max(this%min_mix_len, &
                                           (htot*hbs*vstar)/ &
                                           (this%ekman_scale_coef*absf*htot*hbs + vstar))
                              kd_val = vstar*this%von_karman*mixlen
                           else
                              kd_val = vstar*this%von_karman*(htot*hbs*vstar)/ &
                                       (this%ekman_scale_coef*absf*htot*hbs + vstar)
                           end if
                        else
                           vstar = 0.0_wp
                           kd_val = 0.0_wp
                        end if
                        pe_g0 = pec_core*dkddt/(bdt1*(bdt1 + dkddt*hps))
                        dpe_conv = pec_core*(kd_val*dt_h)/ &
                                   (bdt1*(bdt1 + (kd_val*dt_h)*hps))
                        if (dpe_conv > 0.0_wp) then
                           kd_val = kd_g0
                           dpe_conv = pe_g0
                        end if
                        ! dpe_conv < 0 => the reservoir grows; the
                        ! reservoirs are NOT proportionally drained on
                        ! this branch.
                        conv_perel = conv_perel - dpe_conv
                        d_conv = d_conv - this%nstar*dpe_conv
                        if (sfc_connected) then
                           mld_output = mld_output + ms%h_layer(i, j, kb)
                        end if
                     else
                        ! (8b) stable: direct closed-form Kd from the
                        ! energy budget; drain the reservoirs.
                        if ((pec_core*dkddt <= &
                             tot_tke*(bdt1*(bdt1 + dkddt*hps))) .or. &
                            (pec_core <= 0.0_wp)) then
                           kd_val = kd_g0
                           tke_used = pec_core*dkddt/(bdt1*(bdt1 + dkddt*hps))
                           frac_bl = 1.0_wp
                        else
                           kd_val = (bdt1**2*tot_tke)/ &
                                    (dt_h*(pec_core - bdt1*hps*tot_tke))
                           tke_used = tot_tke
                           frac_bl = tot_tke*(bdt1*(bdt1 + dkddt*hps))/ &
                                     (pec_core*dkddt)
                        end if
                        if (sfc_connected) then
                           mld_output = mld_output + frac_bl*ms%h_layer(i, j, kb)
                        end if
                        if (frac_bl < 1.0_wp) sfc_disconnect = .true.
                        r_reduc = 0.0_wp
                        if (tot_tke > 0.0_wp .and. tot_tke > tke_used) then
                           r_reduc = (tot_tke - tke_used)/tot_tke
                        end if
                        d_mixing = d_mixing + tke_used
                        d_cdecay = d_cdecay + &
                                   (1.0_wp - r_reduc)*(this%nstar - nstar_fc)*conv_perel
                        mech_tke = r_reduc*mech_tke
                        conv_perel = r_reduc*conv_perel
                     end if
                  end if

                  this%kd_int(i, j, ki) = kd_val
                  kddt_cur = kd_val*dt_h

                  ! (9) advance the embedded forward elimination.
                  b1 = 1.0_wp/(hp_a + kddt_cur)
                  c1 = kddt_cur*b1
                  if (ki == nz) then
                     te_lag = b1*(h_ka*this%t0%data(i, j, ka))
                     se_lag = b1*(h_ka*this%s0%data(i, j, ka))
                  else
                     te_new = b1*(h_ka*this%t0%data(i, j, ka) + kddt_prev*te_lag)
                     se_new = b1*(h_ka*this%s0%data(i, j, ka) + kddt_prev*se_lag)
                     te_lag = te_new
                     se_lag = se_new
                  end if
                  hp_a = h_kb + (hp_a*b1)*kddt_cur
                  dpe_t_a = this%dpe_t%data(i, j, kb) + c1*dpe_t_a
                  dpe_s_a = this%dpe_s%data(i, j, kb) + c1*dpe_s_a
                  dch_t_a = this%dcolht_t%data(i, j, kb) + c1*dch_t_a
                  dch_s_a = this%dcolht_s%data(i, j, kb) + c1*dch_s_a
                  th_a = h_kb*this%t0%data(i, j, kb) + kddt_cur*te_lag
                  sh_a = h_kb*this%s0%data(i, j, kb) + kddt_cur*se_lag
                  kddt_prev = kddt_cur
                  if (sfc_disconnect) then
                     htot = ms%h_layer(i, j, kb)
                     sfc_connected = .false.
                  else
                     htot = htot + ms%h_layer(i, j, kb)
                  end if
                  z_int = z_int + ms%h_layer(i, j, kb)
                  pres_int = pres_int + GRAVITY*this%rho0*ms%h_layer(i, j, kb)
               end do
               this%kd_int(i, j, 1) = 0.0_wp

               ! Leftover stocks at the bed count as dissipated.
               d_mdecay = d_mdecay + mech_tke
               d_cdecay = d_cdecay + this%nstar*conv_perel

               mld_found = mld_output
               if (.not. this%mld_iteration) exit
               if (abs(mld_found - mld_guess) < this%mld_tol) exit
               if (obl_it == n_its) exit

               ! Bracket update + next guess (false position with a
               ! bisection fallback; bisection-only when configured).
               if (mld_found > mld_guess) then
                  min_mld = mld_guess
                  dmld_min = mld_found - mld_guess
                  have_min = .true.
               else
                  max_mld = mld_guess
                  dmld_max = mld_found - mld_guess
                  have_max = .true.
               end if
               if (this%mld_bisection) then
                  mld_guess = 0.5_wp*(min_mld + max_mld)
               else if (have_min .and. have_max .and. obl_it > 2 .and. &
                        mod(obl_it - 1, 4) > 0) then
                  mld_guess = min_mld + dmld_min*(max_mld - min_mld)/ &
                              (dmld_min - dmld_max)
               else if (mld_found > min_mld .and. mld_found < max_mld) then
                  mld_guess = mld_found
               else
                  mld_guess = 0.5_wp*(min_mld + max_mld)
               end if
            end do

         end if

         this%mld(i, j) = mld_found
         this%la(i, j) = la_val
         if (this%tke_diags) then
            this%tke_wind(i, j) = d_wind/dt
            this%tke_conv(i, j) = d_conv/dt
            this%tke_forcing(i, j) = d_forcing/dt
            this%tke_mixing(i, j) = d_mixing/dt
            this%tke_mech_decay(i, j) = d_mdecay/dt
            this%tke_conv_decay(i, j) = d_cdecay/dt
         end if
      end do
   end subroutine epbl_column_kernel

   pure subroutine epbl_merge_into_kv_kt(this, nx, ny, nzp1, kv, kt)
      !! Fold the EPBL diffusivity into the vmix interface fields.
      !! Called EVERY stage (the interior closure rewrites kv/kt each
      !! stage; `kd_int` itself refreshes at thermo cadence).
      !! Interior interfaces only — k=1 (bed) and k=nz+1 (surface)
      !! stay at the closed-boundary zero in both source and target.
      type(ocean_epbl_t), intent(in) :: this
      integer, intent(in) :: nx, ny, nzp1
         !! Interface-field extents (explicit shape: assumed-shape
         !! dummies in a `do concurrent` kernel make NVHPC walk the
         !! descriptor with per-launch memcpys — this runs every stage).
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
         !! Momentum viscosity at interfaces; gets prandtl*kd.
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
         !! Tracer diffusivity at interfaces; gets kd.

      integer :: i, j, k

      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx)
         if (this%combine_mode == EPBL_COMBINE_MAX) then
            kt(i, j, k) = max(kt(i, j, k), this%kd_int(i, j, k))
            kv(i, j, k) = max(kv(i, j, k), this%prandtl*this%kd_int(i, j, k))
         else
            kt(i, j, k) = kt(i, j, k) + this%kd_int(i, j, k)
            kv(i, j, k) = kv(i, j, k) + this%prandtl*this%kd_int(i, j, k)
         end if
      end do
   end subroutine epbl_merge_into_kv_kt

   pure function ocean_epbl_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the EPBL slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_epbl_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_centre) &
               + arr_bytes(this%mld) &
               + arr_bytes(this%b0) &
               + arr_bytes(this%kd_int) &
               + arr_bytes(this%la) &
               + arr_bytes(this%tke_wind) &
               + arr_bytes(this%tke_conv) &
               + arr_bytes(this%tke_forcing) &
               + arr_bytes(this%tke_mixing) &
               + arr_bytes(this%tke_mech_decay) &
               + arr_bytes(this%tke_conv_decay) &
               + this%t0%bytes() &
               + this%s0%bytes() &
               + this%dpe_t%bytes() &
               + this%dpe_s%bytes() &
               + this%dcolht_t%bytes() &
               + this%dcolht_s%bytes() &
               + this%ctke_sw%bytes()
   end function ocean_epbl_bytes

end module rdb_ocean_epbl
