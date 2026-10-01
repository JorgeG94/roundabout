!! Mesoscale eddy kinetic energy (MEKE) prognostic slot.
module rdb_ocean_meke
   !! Mesoscale eddy kinetic energy parameterization.  Carries a single 2D
   !! vertically-averaged eddy-energy field `meke(nx,ny)` [m^2/s^2], sourced
   !! by the GM potential-energy release, damped by an implicit
   !! (backward-Euler) bottom drag, transported laterally (harmonic-mass
   !! Laplacian + optional biharmonic + advection), and fed back as a
   !! thickness/tracer diffusivity `kh = khcoeff·sqrt(2·gamma_t²·E)·Lmix`
   !! added (geom mean) into VarMix's per-face KhTh before GM's CFL clamp —
   !! closing the GM↔eddy-energy loop.  Updated via a Strang split each
   !! thermo step.
   !!
   !! Clean-room (no source ported).  `meke` is PROGNOSTIC (restart-persistent).
   !! Default off (`&ocean_meke_nml enable = .false.`) ⇒ `meke_step` never
   !! called ⇒ bit-identical.  Bottom-up convention (k=1 bed, k=nz surface);
   !! only the column mass sum touches k and it is orientation-independent.
   !!
   !! References: Jansen, Adcroft, Hallberg & Held (2015); Eden & Greatbatch
   !! (2008); Marshall, Maddison & Berloff (2012).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, H_VANISHED, H_DIV_EPS, GRAVITY
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, H_VANISHED, H_DIV_EPS, GRAVITY
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_gm, only: ocean_gm_t
   use rdb_ocean_varmix, only: ocean_varmix_t
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t
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

   public :: ocean_meke_t
   public :: meke_step
   public :: meke_backscatter_apply

   real(wp), parameter :: MASS_NEGLECT = 1.0e-30_wp
      !! Floor in the harmonic-mass denominator (MOM6 `mass_neglect`).
   real(wp), parameter :: BACKSCATTER_CFL = 0.8_wp
      !! Forward-Euler viscous-CFL safety coefficient for the backscatter
      !! lower bound (MOM6 `BACKSCATTER_UNDERBOUND` analogue).  The NET
      !! (resolved − Ku) harmonic viscosity is floored at
      !! `−BACKSCATTER_CFL·0.5/(dt·(idx²+idy²))` so the negative mode's
      !! growth RATE is bounded — NOT so the operator is stable on its own
      !! (a negative Laplacian always amplifies; a positive biharmonic
      !! backstop, mandatory at configure, dissipates the fed grid-scale
      !! mode).  Matches the `bound_coef = 0.8` convention the resolved
      !! hvisc CFL limiter uses.

   type :: ocean_meke_t
      !! Mesoscale eddy kinetic energy state.  All fields default to the
      !! inert (`enable=.false.`) configuration so an ocean run that never
      !! sets `&ocean_meke_nml` is bit-identical.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; gate on this (never on
         !! `allocated`, which misses the GPU mapping).
      logical :: enable = .false.
         !! Master switch.  Off ⇒ `meke_step` is never called.  Requires
         !! `&ocean_gm_nml enable` (needs `gm%gm_src`); loud configure check.

      ! ---- Source / sink coefficients ----
      real(wp) :: gmcoeff = -1.0_wp
         !! Efficiency of PE->MEKE conversion (nondim).  < 0 ⇒ GM source off.
      real(wp) :: frcoeff = -1.0_wp
         !! Efficiency of mean->eddy frictional conversion (nondim); < 0
         !! (default) ⇒ off ⇒ bit-identical.  When >= 0, adds the frictional
         !! source `-frcoeff*i_mass*ke_diss` from the lateral-viscosity KE
         !! dissipation (`hvisc%ke_diss`).
      real(wp) :: bgsrc = 0.0_wp
         !! Background energy source (m^2/s^3).
      real(wp) :: damping = 0.0_wp
         !! Local depth-independent linear MEKE dissipation rate (1/s).

      ! ---- Lateral transport ----
      real(wp) :: kh = -1.0_wp
         !! Background lateral diffusion of MEKE (m^2/s).  < 0 ⇒ diffusion
         !! stage off.
      real(wp) :: k4 = -1.0_wp
         !! Background biharmonic diffusion of MEKE (m^4/s).  < 0 ⇒ off.
      real(wp) :: khmeke_fac = 0.0_wp
         !! Factor relating `meke%kh` (the derived diffusivity) to the
         !! diffusivity used for MEKE's own lateral spreading (nondim).
      real(wp) :: advection_factor = 0.0_wp
         !! Scaling on the barotropic-transport advection of MEKE (nondim);
         !! 0 (default) ⇒ the advection stage is skipped ⇒ bit-identical.

      ! ---- KhTh closure ----
      real(wp) :: khcoeff = 1.0_wp
         !! Scaling converting MEKE into Kh (nondim).  <= 0 ⇒ closure off.
      real(wp) :: cd_scale = 0.0_wp
         !! Ratio of bottom eddy velocity to column-mean eddy velocity
         !! (nondim); enters bottomFac2.
      real(wp) :: cb = 25.0_wp
         !! Coefficient in the gamma_bot (bottomFac2) expression (nondim).
      real(wp) :: ct = 50.0_wp
         !! Coefficient in the gamma_bt (barotrFac2) expression (nondim).
      real(wp) :: min_gamma2 = 1.0e-4_wp
         !! Floor on gamma_b^2 / gamma_t^2 (nondim).
      real(wp) :: uscale = 0.0_wp
         !! Background (e.g. tidal) eddy velocity scale for bottom drag (m/s).
      real(wp) :: dtscale = 1.0_wp
         !! Scale factor accelerating MEKE time-stepping (nondim).
      real(wp) :: cdrag = 2.5e-3_wp
         !! Bottom drag coefficient for MEKE (nondim).  Copied from
         !! `&ocean_bdrag_nml cdrag_side` at configure if that is > 0,
         !! else this default; enters drag_rate + Lfrict.
      logical :: use_bbl_drag = .false.
         !! Add the resolved bed-layer eddy velocity `|u_bed|²` to the MEKE
         !! bottom-drag rate `drag_rate = rho0·i_mass·sqrt(cdrag²·(2·bf2·E +
         !! |u_bed|² + uscale²))` (MOM6 `drag_rate_visc`).  Default off ⇒ the
         !! `u_bbl²` workspace stays 0 ⇒ bit-identical to the prior drag.

      ! ---- Length-scale alpha weights (default all 0 ⇒ Lmix path inert) ----
      real(wp) :: alpha_deform = 0.0_wp
         !! Weight on the deformation length scale Ldeform (nondim).
      real(wp) :: alpha_rhines = 0.0_wp
         !! Weight on the Rhines length scale Lrhines = sqrt(Ueddy/beta)
         !! (nondim).  Default 0 ⇒ beta term inert; > 0 activates the
         !! Rhines scale (needs `f_centre` filled, done at setup).
      real(wp) :: alpha_eady = 0.0_wp
         !! Weight on the Eady length scale Leady (needs VarMix SN) (nondim).
      real(wp) :: alpha_frict = 0.0_wp
         !! Weight on the frictional-arrest length scale Lfrict (nondim).
      real(wp) :: alpha_grid = 0.0_wp
         !! Weight on the grid length scale Lgrid (nondim).

      ! ---- Feedback seam factors ----
      real(wp) :: khth_fac = 0.0_wp
         !! Factor on the geometric-mean kh added into VarMix's KhTh face
         !! field (nondim).  0 (default) ⇒ feedback inert ⇒ bit-identical.
      real(wp) :: khtr_fac = 0.0_wp
         !! Factor on the geometric-mean kh added into VarMix's KhTr face
         !! field (nondim).  0 (default) ⇒ inert.

      ! ---- Backscatter (negative-viscosity momentum energy return, v1) ----
      logical :: backscatter = .false.
         !! Master switch for the MEKE → momentum negative-viscosity
         !! energy return (capability Gap 2).  Default `.false.` ⇒ the
         !! `ku` field stays 0 and `meke_backscatter_apply` adds nothing
         !! ⇒ bit-identical.  When `.true.`, `meke_step` fills `ku` and the
         !! driver subtracts it from the per-face harmonic viscosity,
         !! returning eddy energy to the resolved flow.  HARMONIC only in
         !! v1 (the biharmonic `Au` return + EBT/SQG vertical structure
         !! `BS_struct` are deferred).  Needs `enable=.true.` (so `meke_step`
         !! runs and refreshes `ku`) — inert otherwise.
      real(wp) :: visc_coeff_ku = 0.0_wp
         !! MOM6 `MEKE_VISCOSITY_COEFF_KU` — the nondimensional efficiency
         !! of the harmonic backscatter viscosity
         !! `Ku = visc_coeff_ku·sqrt(2·gamma_t²·E)·Lmix` (m²/s).  May be
         !! negative in MOM6 (negative viscosity); here it is the magnitude
         !! coefficient and the SIGN of the momentum effect is set by the
         !! SUBTRACTION in `meke_backscatter_apply` (`A_net = A − Ku`), so a
         !! positive `visc_coeff_ku` returns energy.  0 (default) ⇒ inert.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- Prognostic + derived 2D fields ----
      real(wp), allocatable :: meke(:, :)
         !! Eddy kinetic energy E (m^2/s^2), `(nx,ny)`.  PROGNOSTIC —
         !! restart-persistent.
      real(wp), allocatable :: kh_diff(:, :)
         !! Derived MEKE diffusivity kh (m^2/s), `(nx,ny)`; fed (geom-mean)
         !! into VarMix's face KhTh/KhTr.
      real(wp), allocatable :: le(:, :)
         !! Mixing length scale Lmix (m), `(nx,ny)` (diagnostic).
      real(wp), allocatable :: ku(:, :)
         !! Derived harmonic backscatter viscosity Ku (m²/s), `(nx,ny)`;
         !! `ku = visc_coeff_ku·sqrt(2·gamma_t²·E)·Lmix`.  Zero unless
         !! `backscatter`.  Subtracted (face-averaged, stability-floored)
         !! from the resolved harmonic viscosity by `meke_backscatter_apply`.

      ! ---- Strang-stage workspace ----
      real(wp), allocatable :: i_mass(:, :)
         !! 1 / column mass (m^2/kg), `(nx,ny)`.
      real(wp), allocatable :: depth_tot(:, :)
         !! Total column thickness (m), `(nx,ny)`; used in Lfrict.
      real(wp), allocatable :: bottom_fac2(:, :)
         !! gamma_b^2 (nondim), `(nx,ny)`.
      real(wp), allocatable :: barotr_fac2(:, :)
         !! gamma_t^2 (nondim), `(nx,ny)`.
      real(wp), allocatable :: src(:, :)
         !! Aggregate source (m^2/s^3), `(nx,ny)`.
      real(wp), allocatable :: uflux(:, :)
         !! u-face MEKE flux workspace, `(nx+1,ny)`.
      real(wp), allocatable :: vflux(:, :)
         !! v-face MEKE flux workspace, `(nx,ny+1)`.
      real(wp), allocatable :: del2(:, :)
         !! Laplacian of MEKE workspace (biharmonic), `(nx,ny)`.
      real(wp), allocatable :: mass_ws(:, :)
         !! Column mass (kg/m^2), `(nx,ny)`; the harmonic-mass input for the
         !! lateral flux (= 1/i_mass where i_mass>0).
      real(wp), allocatable :: u_bbl2(:, :)
         !! Resolved bed-layer speed² (m²/s²) at cell centres, `(nx,ny)`;
         !! filled from `ms` when `use_bbl_drag`, else 0.  Feeds `meke_drag`.
      real(wp), allocatable :: ke_diss_ws(:, :)
         !! Staged hvisc KE-dissipation rate (kg/s³, ≤0), `(nx,ny)`; copied
         !! from `hv%ke_diss` when the frictional source is wired, else 0.
         !! Feeds `meke_source` as `-frcoeff·i_mass·ke_diss`.
      real(wp), allocatable :: rd_ws(:, :)
         !! Staged copy of `wavespeed%rd_over_dx` (nondim), `(nx,ny)`; zero
         !! when wavespeed absent ⇒ Ldeform→0.
      real(wp), allocatable :: f_centre(:, :)
         !! |f| at cell centres (1/s), `(nx,ny)`; filled by `set_f_centre`
         !! from the same beta-plane the Coriolis slot uses (mirrors EPBL /
         !! kappa-shear).  Drives `beta = |grad f|` for the Rhines length.
         !! Zero until filled ⇒ Rhines weight inert (alpha_rhines default 0).
      real(wp), allocatable :: sn_u_ws(:, :)
         !! Staged copy of `varmix%sn_u` (1/s), `(nx+1,ny)`; zero when VarMix
         !! absent ⇒ Eady scale inert.
      real(wp), allocatable :: sn_v_ws(:, :)
         !! Staged copy of `varmix%sn_v` (1/s), `(nx,ny+1)`.
      real(wp), allocatable :: baro_hu(:, :)
         !! Depth-integrated u-face mass transport (kg/s), `(nx+1,ny)`;
         !! the barotropic transport that advects E.  Filled only when
         !! `advection_factor > 0`.
      real(wp), allocatable :: baro_hv(:, :)
         !! Depth-integrated v-face mass transport (kg/s), `(nx,ny+1)`.
   contains
      procedure, non_overridable :: init => ocean_meke_init
      procedure, non_overridable :: destroy => ocean_meke_destroy
      procedure, non_overridable :: enter_data => ocean_meke_enter_data
      procedure, non_overridable :: exit_data => ocean_meke_exit_data
      procedure, non_overridable :: set_f_centre => ocean_meke_set_f_centre
      procedure, non_overridable :: bytes => ocean_meke_bytes
   end type ocean_meke_t

contains

   subroutine ocean_meke_init(this, grid, nz_ml)
      !! Allocate the prognostic field, the derived diffusivity, and the
      !! Strang-stage workspaces.  Always allocates (configure runs after
      !! init); setup uses plain host allocation (no `do concurrent` before
      !! enter_data).
      class(ocean_meke_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml
      if (nz < 1) nz = 1
      if (nz > NZ_STACK_MAX) then
         error stop "ocean_meke_init: nz_ml exceeds NZ_STACK_MAX "// &
            "(raise NZ_STACK_MAX in rdb_constants)"
      end if
      this%nx_total = nx
      this%ny_total = ny
      this%nz_ml = nz

      allocate (this%meke(nx, ny), source=0.0_wp)
      allocate (this%kh_diff(nx, ny), source=0.0_wp)
      allocate (this%le(nx, ny), source=0.0_wp)
      allocate (this%ku(nx, ny), source=0.0_wp)
      allocate (this%i_mass(nx, ny), source=0.0_wp)
      allocate (this%depth_tot(nx, ny), source=0.0_wp)
      allocate (this%bottom_fac2(nx, ny), source=0.0_wp)
      allocate (this%barotr_fac2(nx, ny), source=0.0_wp)
      allocate (this%src(nx, ny), source=0.0_wp)
      allocate (this%uflux(nx + 1, ny), source=0.0_wp)
      allocate (this%vflux(nx, ny + 1), source=0.0_wp)
      allocate (this%del2(nx, ny), source=0.0_wp)
      allocate (this%mass_ws(nx, ny), source=0.0_wp)
      allocate (this%u_bbl2(nx, ny), source=0.0_wp)
      allocate (this%ke_diss_ws(nx, ny), source=0.0_wp)
      allocate (this%rd_ws(nx, ny), source=0.0_wp)
      allocate (this%f_centre(nx, ny), source=0.0_wp)
      allocate (this%sn_u_ws(nx + 1, ny), source=0.0_wp)
      allocate (this%sn_v_ws(nx, ny + 1), source=0.0_wp)
      allocate (this%baro_hu(nx + 1, ny), source=0.0_wp)
      allocate (this%baro_hv(nx, ny + 1), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_meke_init

   subroutine ocean_meke_destroy(this)
      class(ocean_meke_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%meke)) deallocate (this%meke)
      if (allocated(this%kh_diff)) deallocate (this%kh_diff)
      if (allocated(this%le)) deallocate (this%le)
      if (allocated(this%ku)) deallocate (this%ku)
      if (allocated(this%i_mass)) deallocate (this%i_mass)
      if (allocated(this%depth_tot)) deallocate (this%depth_tot)
      if (allocated(this%bottom_fac2)) deallocate (this%bottom_fac2)
      if (allocated(this%barotr_fac2)) deallocate (this%barotr_fac2)
      if (allocated(this%src)) deallocate (this%src)
      if (allocated(this%uflux)) deallocate (this%uflux)
      if (allocated(this%vflux)) deallocate (this%vflux)
      if (allocated(this%del2)) deallocate (this%del2)
      if (allocated(this%mass_ws)) deallocate (this%mass_ws)
      if (allocated(this%u_bbl2)) deallocate (this%u_bbl2)
      if (allocated(this%ke_diss_ws)) deallocate (this%ke_diss_ws)
      if (allocated(this%rd_ws)) deallocate (this%rd_ws)
      if (allocated(this%f_centre)) deallocate (this%f_centre)
      if (allocated(this%sn_u_ws)) deallocate (this%sn_u_ws)
      if (allocated(this%sn_v_ws)) deallocate (this%sn_v_ws)
      if (allocated(this%baro_hu)) deallocate (this%baro_hu)
      if (allocated(this%baro_hv)) deallocate (this%baro_hv)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_meke_destroy

   subroutine ocean_meke_enter_data(this)
      ! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule).
      class(ocean_meke_t), intent(inout) :: this
      select type (this)
      type is (ocean_meke_t)
         call ocean_meke_enter_data_impl(this)
      end select
   end subroutine ocean_meke_enter_data

   subroutine ocean_meke_enter_data_impl(this)
      type(ocean_meke_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%meke, this%kh_diff, this%le, this%ku)
      !$acc enter data copyin(this%i_mass, this%depth_tot)
      !$acc enter data copyin(this%bottom_fac2, this%barotr_fac2, this%src)
      !$acc enter data copyin(this%uflux, this%vflux, this%del2)
      !$acc enter data copyin(this%mass_ws, this%rd_ws, this%f_centre)
      !$acc enter data copyin(this%u_bbl2, this%ke_diss_ws)
      !$acc enter data copyin(this%sn_u_ws, this%sn_v_ws)
      !$acc enter data copyin(this%baro_hu, this%baro_hv)
   end subroutine ocean_meke_enter_data_impl

   subroutine ocean_meke_exit_data(this)
      class(ocean_meke_t), intent(inout) :: this
      select type (this)
      type is (ocean_meke_t)
         call ocean_meke_exit_data_impl(this)
      end select
   end subroutine ocean_meke_exit_data

   subroutine ocean_meke_exit_data_impl(this)
      type(ocean_meke_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%baro_hu, this%baro_hv)
      !$acc exit data delete(this%sn_u_ws, this%sn_v_ws)
      !$acc exit data delete(this%u_bbl2, this%ke_diss_ws)
      !$acc exit data delete(this%mass_ws, this%rd_ws, this%f_centre)
      !$acc exit data delete(this%uflux, this%vflux, this%del2)
      !$acc exit data delete(this%bottom_fac2, this%barotr_fac2, this%src)
      !$acc exit data delete(this%i_mass, this%depth_tot)
      !$acc exit data delete(this%meke, this%kh_diff, this%le, this%ku)
   end subroutine ocean_meke_exit_data_impl

   subroutine ocean_meke_set_f_centre(this, grid, f_centre)
      !! Copy a pre-filled cell-centre Coriolis magnitude |f| (1/s) onto the
      !! MEKE slot, so `beta = |grad f|` for the Rhines length is live.  The
      !! caller (setup) builds `f_centre` with the same `metrics_fill_coriolis`
      !! path the Coriolis / VarMix / EPBL slots use (handles beta-plane AND
      !! spherical).  Host loop — call after `init`, before `enter_data`.
      class(ocean_meke_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_centre(grid%nx_total, grid%ny_total)
      integer :: i, j

      if (.not. allocated(this%f_centre)) return
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            this%f_centre(i, j) = f_centre(i, j)
         end do
      end do
   end subroutine ocean_meke_set_f_centre

   ! =================================================================
   ! Compute driver
   ! =================================================================

   subroutine meke_step(grid, metrics, this, gm, varmix, wavespeed, ms, dt, ke_diss_ext)
      !! Advance the MEKE field one thermo step (Strang split), update the
      !! derived diffusivity `kh_diff`, and feed the geometric-mean kh into
      !! VarMix's per-face KhTh/KhTr (the GM↔MEKE feedback).  Run once per
      !! outer step at thermo cadence, after `varmix_compute` and before
      !! `gm_compute_transports` (MEKE reads `gm%gm_src` from the previous
      !! thermo step — a one-step lag).  No-op when `enable=.false.`,
      !! uninitialised, or the GM slot is absent.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_meke_t), intent(inout) :: this
      type(ocean_gm_t), intent(in) :: gm
      type(ocean_varmix_t), intent(inout), optional :: varmix
      type(ocean_wave_speed_t), intent(in), optional :: wavespeed
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: ke_diss_ext(:, :)
         !! hvisc KE-dissipation rate `(nx,ny)` for the frictional source;
         !! absent ⇒ the source is inert (staged to 0).

      integer :: nx, ny, nz
      real(wp) :: sdt, sdt_damp, damp_step
      logical :: have_rho, kh_flux_enabled, have_ws, have_vm

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. gm%is_init) return
      if (.not. allocated(ms%h_layer)) return
      if (.not. allocated(gm%gm_src)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      if (this%nz_ml /= nz) return

      sdt = dt*this%dtscale
      damp_step = 1.0_wp
      if (this%kh >= 0.0_wp .or. this%k4 >= 0.0_wp) damp_step = 0.5_wp
      sdt_damp = sdt*damp_step
      kh_flux_enabled = (this%kh >= 0.0_wp)
      have_rho = allocated(ms%rho_layer)
      have_ws = .false.
      if (present(wavespeed)) have_ws = wavespeed%is_init .and. allocated(wavespeed%rd_over_dx)
      have_vm = .false.
      if (present(varmix)) have_vm = varmix%is_init .and. allocated(varmix%sn_u)

      ! ---- 0. stage the optional inputs into device-resident workspaces.
      ! Copy ON-DEVICE from the slot arrays (which are device-resident) so the
      ! length-scale kernel reads explicit-shape slot args, never an
      ! optional-or-host-staged actual.  When the source slot is absent the
      ! workspace is zeroed (Ldeform->0 / Eady inert).  `f_centre` is the
      ! slot's own |f| (filled at setup via set_f_centre); when never filled
      ! it is zero ⇒ beta=0 ⇒ Rhines inert (alpha_rhines default 0).
      if (have_ws) then
         call meke_stage_rd(nx, ny, wavespeed%rd_over_dx, this%rd_ws)
      else
         call meke_zero_2d(nx, ny, this%rd_ws)
      end if
      if (have_vm) then
         call meke_stage_sn(nx, ny, varmix%sn_u, varmix%sn_v, &
                            this%sn_u_ws, this%sn_v_ws)
      else
         call meke_zero_2d(nx + 1, ny, this%sn_u_ws)
         call meke_zero_2d(nx, ny + 1, this%sn_v_ws)
      end if

      ! ---- 1. column mass + depth (I_mass, depth_tot, mass_ws). ----
      ! Pass rho_layer only when allocated; otherwise feed h_layer as a
      ! harmless placeholder (gated off by have_rho ⇒ mass uses gm%rho0).
      if (have_rho) then
         call meke_mass(nx, ny, nz, gm%rho0, .true., ms%h_layer, ms%rho_layer, &
                        this%i_mass, this%depth_tot, this%mass_ws)
      else
         call meke_mass(nx, ny, nz, gm%rho0, .false., ms%h_layer, ms%h_layer, &
                        this%i_mass, this%depth_tot, this%mass_ws)
      end if

      ! ---- 2. structure factors (bottomFac2, barotrFac2) + Lmix into le. ----
      call meke_length_scales(nx, ny, this%cd_scale, this%cb, this%ct, &
                              this%min_gamma2, this%cdrag, &
                              this%alpha_deform, this%alpha_rhines, &
                              this%alpha_eady, this%alpha_frict, this%alpha_grid, &
                              metrics%areaT, metrics%idxT, metrics%idyT, &
                              this%f_centre, &
                              this%depth_tot, this%meke, &
                              this%rd_ws, this%sn_u_ws, this%sn_v_ws, &
                              this%bottom_fac2, this%barotr_fac2, this%le)

      ! ---- 3. explicit source bump: E += sdt*src. ----
      ! Stage the hvisc KE-dissipation rate for the frictional source
      ! (-frcoeff·i_mass·ke_diss); 0 when the seam is not wired ⇒ inert.
      if (present(ke_diss_ext)) then
         call meke_stage_rd(nx, ny, ke_diss_ext, this%ke_diss_ws)
      else
         call meke_zero_2d(nx, ny, this%ke_diss_ws)
      end if
      call meke_source(nx, ny, this%bgsrc, this%gmcoeff, this%frcoeff, sdt, &
                       this%i_mass, gm%gm_src, this%ke_diss_ws, this%src, this%meke)

      ! ---- 4. implicit drag half: E <- E/(1+sdt_damp*damp_rate). ----
      ! Resolved bed-layer eddy velocity (MOM6 drag_rate_visc); 0 when off.
      if (this%use_bbl_drag) then
         call meke_bbl_speed2(nx, ny, nz, ms%u_face_x_layer, ms%v_face_y_layer, this%u_bbl2)
      end if
      call meke_drag(nx, ny, sdt_damp, this%damping, this%cdrag, this%uscale, gm%rho0, &
                     this%i_mass, this%bottom_fac2, this%u_bbl2, this%meke)

      ! ---- 5. lateral diffusion (+ biharmonic). ----
      if (kh_flux_enabled .or. this%k4 >= 0.0_wp) then
         call meke_lateral(nx, ny, sdt, this%kh, this%k4, this%khmeke_fac, &
                           kh_flux_enabled, &
                           metrics%dy_cu, metrics%dx_cv, metrics%idxCu, &
                           metrics%idyCv, metrics%iareaT, &
                           this%i_mass, this%mass_ws, &
                           this%kh_diff, this%uflux, this%vflux, this%del2, &
                           this%meke)
      end if

      ! ---- 5b. upwind advection of E by the barotropic transport. ----
      ! Self-contained transport stage: baroHu = Sum_k (mass-weighted face
      ! transport), upwind flux E*baroHu, single conservative divergence.
      ! `advection_factor = 0` (default) ⇒ no-op ⇒ bit-identical.
      if (this%advection_factor > 0.0_wp) then
         if (have_rho) then
            call meke_baro_transport(nx, ny, nz, gm%rho0, .true., &
                                     ms%mass_flux_x_layer, ms%mass_flux_y_layer, &
                                     ms%rho_layer, this%baro_hu, this%baro_hv)
         else
            call meke_baro_transport(nx, ny, nz, gm%rho0, .false., &
                                     ms%mass_flux_x_layer, ms%mass_flux_y_layer, &
                                     ms%rho_layer, this%baro_hu, this%baro_hv)
         end if
         call meke_advect(nx, ny, sdt, this%advection_factor, metrics%iareaT, &
                          this%i_mass, this%baro_hu, this%baro_hv, &
                          this%uflux, this%vflux, this%meke)
      end if

      ! ---- 6. implicit drag half (only when Strang-split, damp_step=0.5). ----
      ! `u_bbl2` already filled in step 4 (or held at 0 when off).
      if (this%kh >= 0.0_wp .or. this%k4 >= 0.0_wp) then
         call meke_drag(nx, ny, sdt_damp, this%damping, this%cdrag, this%uscale, gm%rho0, &
                        this%i_mass, this%bottom_fac2, this%u_bbl2, this%meke)
      end if

      ! ---- 7. MEKE -> KhTh closure: kh = khcoeff*sqrt(2*gamma_t2*E)*Lmix. ----
      call meke_kh_closure(nx, ny, this%khcoeff, &
                           this%barotr_fac2, this%meke, this%le, this%kh_diff)

      ! ---- 7b. MEKE -> Ku backscatter coefficient (v1, harmonic only). ----
      ! ku = visc_coeff_ku*sqrt(2*gamma_t2*E)*Lmix; same eddy-velocity ×
      ! mixing-length form as kh.  Filled only when `backscatter`; the
      ! driver subtracts a face-average of `ku` from the resolved harmonic
      ! viscosity (stability-floored) ⇒ a negative-viscosity energy return.
      ! Default off ⇒ ku stays 0 ⇒ bit-identical.
      call meke_ku_closure(nx, ny, this%backscatter, this%visc_coeff_ku, &
                           this%meke, this%le, this%ku)

      ! ---- 8. feedback seam: add geom-mean kh into VarMix face KhTh/KhTr. ----
      if (present(varmix)) then
         if (varmix%is_init .and. (this%khth_fac /= 0.0_wp .or. &
                                   this%khtr_fac /= 0.0_wp)) then
            call meke_feed_khth(nx, ny, this%khth_fac, this%khtr_fac, &
                                this%kh_diff, varmix%khth_u, varmix%khth_v, &
                                varmix%khtr_u, varmix%khtr_v)
         end if
      end if
   end subroutine meke_step

   pure subroutine meke_zero_2d(n1, n2, a)
      !! Zero a device-resident 2D workspace (absent-source fallback).
      integer, intent(in) :: n1, n2
      real(wp), intent(out) :: a(n1, n2)
      integer :: i, j
      do concurrent(j=1:n2, i=1:n1)
         a(i, j) = 0.0_wp
      end do
   end subroutine meke_zero_2d

   pure subroutine meke_stage_rd(nx, ny, rd_in, rd_ws)
      !! Copy rd_over_dx onto the device workspace.  Thermo-cadence;
      !! explicit-shape; on-device (source slot is device-resident).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: rd_in(nx, ny)
      real(wp), intent(out) :: rd_ws(nx, ny)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         rd_ws(i, j) = rd_in(i, j)
      end do
   end subroutine meke_stage_rd

   pure subroutine meke_stage_sn(nx, ny, sn_u_in, sn_v_in, sn_u_ws, sn_v_ws)
      !! Copy the SN faces onto the device workspaces.  Thermo-cadence;
      !! explicit-shape; on-device (VarMix slot is device-resident).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sn_u_in(nx + 1, ny)
      real(wp), intent(in) :: sn_v_in(nx, ny + 1)
      real(wp), intent(out) :: sn_u_ws(nx + 1, ny)
      real(wp), intent(out) :: sn_v_ws(nx, ny + 1)
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx + 1)
         sn_u_ws(i, j) = sn_u_in(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         sn_v_ws(i, j) = sn_v_in(i, j)
      end do
   end subroutine meke_stage_sn

   ! =================================================================
   ! Kernels
   ! =================================================================

   pure subroutine meke_mass(nx, ny, nz, rho0, have_rho, h_layer, rho_layer, &
                             i_mass, depth_tot, mass_ws)
      !! Column mass `mass = Sum_k rho*max(h,H_VANISHED)` (kg/m^2), its
      !! inverse `i_mass` (0 where mass<=0), `depth_tot = Sum_k h` (m), and
      !! `mass_ws = mass` (the harmonic-mass input for the lateral flux).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: rho0
      logical, intent(in) :: have_rho
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: rho_layer(nx, ny, nz)
      real(wp), intent(out) :: i_mass(nx, ny)
      real(wp), intent(out) :: depth_tot(nx, ny)
      real(wp), intent(out) :: mass_ws(nx, ny)
      integer :: i, j, k
      real(wp) :: mass, dsum, hk, rhok

      do concurrent(j=1:ny, i=1:nx) local(k, mass, dsum, hk, rhok)
         mass = 0.0_wp
         dsum = 0.0_wp
         do k = 1, nz
            hk = max(h_layer(i, j, k), H_VANISHED)
            rhok = rho0
            if (have_rho) rhok = rho_layer(i, j, k)
            mass = mass + rhok*hk
            dsum = dsum + h_layer(i, j, k)
         end do
         depth_tot(i, j) = dsum
         mass_ws(i, j) = mass
         if (mass > 0.0_wp) then
            i_mass(i, j) = 1.0_wp/mass
         else
            i_mass(i, j) = 0.0_wp
         end if
      end do
   end subroutine meke_mass

   pure subroutine meke_baro_transport(nx, ny, nz, rho0, have_rho, &
                                       mflux_x, mflux_y, rho_layer, &
                                       baro_hu, baro_hv)
      !! Depth-integrated, mass-weighted barotropic transport through each
      !! C-grid face: `baroHu(I,j) = Sum_k rho_face * mass_flux_x_layer`,
      !! where `mass_flux_*_layer` is the per-layer VOLUME transport (m^3/s,
      !! = u*h_face*dy) and `rho_face` the two-cell average density.  The
      !! result is a MASS transport (kg/s) so the advective divergence pairs
      !! exactly with `IareaT*I_mass` (1/(area*mass)) ⇒ Sum E*area*mass is
      !! conserved.  Array-edge faces carry zero transport (closed domain).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: rho0
      logical, intent(in) :: have_rho
      real(wp), intent(in) :: mflux_x(nx + 1, ny, nz)
      real(wp), intent(in) :: mflux_y(nx, ny + 1, nz)
      real(wp), intent(in) :: rho_layer(nx, ny, nz)
      real(wp), intent(out) :: baro_hu(nx + 1, ny)
      real(wp), intent(out) :: baro_hv(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: tsum, rho_face

      ! u-faces: interior I=2..nx between cells I-1 and I.
      do concurrent(j=1:ny, i=1:nx + 1) local(k, tsum, rho_face)
         tsum = 0.0_wp
         if (i >= 2 .and. i <= nx) then
            do k = 1, nz
               rho_face = rho0
               if (have_rho) rho_face = 0.5_wp*(rho_layer(i - 1, j, k) + rho_layer(i, j, k))
               tsum = tsum + rho_face*mflux_x(i, j, k)
            end do
         end if
         baro_hu(i, j) = tsum
      end do
      ! v-faces: interior J=2..ny between cells J-1 and J.
      do concurrent(j=1:ny + 1, i=1:nx) local(k, tsum, rho_face)
         tsum = 0.0_wp
         if (j >= 2 .and. j <= ny) then
            do k = 1, nz
               rho_face = rho0
               if (have_rho) rho_face = 0.5_wp*(rho_layer(i, j - 1, k) + rho_layer(i, j, k))
               tsum = tsum + rho_face*mflux_y(i, j, k)
            end do
         end if
         baro_hv(i, j) = tsum
      end do
   end subroutine meke_baro_transport

   pure subroutine meke_advect(nx, ny, sdt, adv_fac, iareaT, i_mass, &
                               baro_hu, baro_hv, uflux, vflux, meke)
      !! Upwind flux-form advection of E by the barotropic mass transport.
      !!   advFac = adv_fac/sdt
      !!   uflux(I) = baroHu(I)*advFac*E_upwind   (E_{I-1} if baroHu>0 else E_I)
      !!   E += sdt*IareaT*I_mass*((uflux_{i-1}-uflux_i)+(vflux_{j-1}-vflux_j))
      !! Conservative on a closed domain (interior faces only; the
      !! divergence telescopes ⇒ Sum E*area*mass conserved).  `adv_fac=0`
      !! never reaches here (gated by the caller) ⇒ default bit-identity.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sdt, adv_fac
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: i_mass(nx, ny)
      real(wp), intent(in) :: baro_hu(nx + 1, ny)
      real(wp), intent(in) :: baro_hv(nx, ny + 1)
      real(wp), intent(inout) :: uflux(nx + 1, ny)
      real(wp), intent(inout) :: vflux(nx, ny + 1)
      real(wp), intent(inout) :: meke(nx, ny)
      integer :: i, j
      real(wp) :: adv_per_t, bh, mke

      adv_per_t = 0.0_wp
      if (sdt > 0.0_wp) adv_per_t = adv_fac/sdt

      ! u-face upwind flux (array-edge faces carry zero ⇒ closed domain).
      do concurrent(j=1:ny, i=1:nx + 1)
         uflux(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=2:nx) local(bh)
         bh = baro_hu(i, j)
         if (bh > 0.0_wp) then
            uflux(i, j) = bh*adv_per_t*meke(i - 1, j)
         else if (bh < 0.0_wp) then
            uflux(i, j) = bh*adv_per_t*meke(i, j)
         end if
      end do
      ! v-face upwind flux.
      do concurrent(j=1:ny + 1, i=1:nx)
         vflux(i, j) = 0.0_wp
      end do
      do concurrent(j=2:ny, i=1:nx) local(bh)
         bh = baro_hv(i, j)
         if (bh > 0.0_wp) then
            vflux(i, j) = bh*adv_per_t*meke(i, j - 1)
         else if (bh < 0.0_wp) then
            vflux(i, j) = bh*adv_per_t*meke(i, j)
         end if
      end do
      ! conservative divergence (uflux at face i is the LEFT face of cell i;
      ! face i+1 the RIGHT face).  Inflow-left minus outflow-right.
      do concurrent(j=1:ny, i=1:nx) local(mke)
         mke = sdt*(iareaT(i, j)*i_mass(i, j))* &
               ((uflux(i, j) - uflux(i + 1, j)) + (vflux(i, j) - vflux(i, j + 1)))
         meke(i, j) = meke(i, j) + mke
      end do
   end subroutine meke_advect

   pure function meke_inv_lmix(ueddy, sn, beta, area, rd_over_dx, depth, cdrag, &
                               a_deform, a_rhines, a_eady, a_frict, a_grid) result(inv_l)
      !$acc routine seq
      !! Harmonic inverse mixing length `1/Lmix = Sum aX/LX` over the five
      !! length scales (deformation, frictional, Rhines, Eady, grid).  Each
      !! scale is gated `aX*LX > 0` so a zero weight or a degenerate scale
      !! contributes nothing.  Returns 1/Lmix (0 ⇒ Lmix degenerate).
      real(wp), intent(in) :: ueddy, sn, beta, area, rd_over_dx, depth, cdrag
      real(wp), intent(in) :: a_deform, a_rhines, a_eady, a_frict, a_grid
      real(wp) :: inv_l
      real(wp) :: lgrid, ldeform, lfrict, lrhines, leady

      lgrid = sqrt(max(area, 0.0_wp))
      ldeform = lgrid*rd_over_dx
      lfrict = 0.0_wp
      if (cdrag > 0.0_wp) lfrict = depth/cdrag
      lrhines = 0.0_wp
      if (beta > 0.0_wp) lrhines = sqrt(max(ueddy, 0.0_wp)/beta)
      leady = 0.0_wp
      if (sn > 1.0e-15_wp) leady = ueddy/sn

      inv_l = 0.0_wp
      if (a_deform*ldeform > 0.0_wp) inv_l = inv_l + 1.0_wp/(a_deform*ldeform)
      if (a_frict*lfrict > 0.0_wp) inv_l = inv_l + 1.0_wp/(a_frict*lfrict)
      if (a_rhines*lrhines > 0.0_wp) inv_l = inv_l + 1.0_wp/(a_rhines*lrhines)
      if (a_eady*leady > 0.0_wp) inv_l = inv_l + 1.0_wp/(a_eady*leady)
      if (a_grid*lgrid > 0.0_wp) inv_l = inv_l + 1.0_wp/(a_grid*lgrid)
   end function meke_inv_lmix

   pure subroutine meke_length_scales(nx, ny, cd_scale, cb, ct, min_gamma2, &
                                      cdrag, a_deform, a_rhines, &
                                      a_eady, a_frict, a_grid, &
                                      areaT, idxT, idyT, f_centre, &
                                      depth_tot, meke, &
                                      rd_over_dx, sn_u, sn_v, &
                                      bottom_fac2, barotr_fac2, le)
      !! Fill the structure factors gamma_b^2 (`bottom_fac2`) and gamma_t^2
      !! (`barotr_fac2`) plus the mixing length `le` (Lmix) at each cell
      !! centre.  `Ldeform/Lfrict` drives both gammas; `Lmix` is the harmonic
      !! sum of the alpha-weighted scales.  `beta = |grad f|` from centred
      !! `f_centre` differences scaled by `idxT`/`idyT` (zero when f_centre
      !! is unfilled ⇒ Rhines inert).  SN = 0.25*(sn_u(i)+sn_u(i-1)+
      !! sn_v(j)+sn_v(j-1)) only when aEady>0.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: cd_scale, cb, ct, min_gamma2, cdrag
      real(wp), intent(in) :: a_deform, a_rhines, a_eady, a_frict, a_grid
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: idxT(nx, ny)
      real(wp), intent(in) :: idyT(nx, ny)
      real(wp), intent(in) :: f_centre(nx, ny)
      real(wp), intent(in) :: depth_tot(nx, ny)
      real(wp), intent(in) :: meke(nx, ny)
      real(wp), intent(in) :: rd_over_dx(nx, ny)
      real(wp), intent(in) :: sn_u(nx + 1, ny)
      real(wp), intent(in) :: sn_v(nx, ny + 1)
      real(wp), intent(out) :: bottom_fac2(nx, ny)
      real(wp), intent(out) :: barotr_fac2(nx, ny)
      real(wp), intent(out) :: le(nx, ny)
      integer :: i, j
      real(wp) :: lgrid, ldeform, lfrict, ratio, bf2, tf2
      real(wp) :: ueddy, sn, beta, inv_l, rd

      do concurrent(j=1:ny, i=1:nx) &
         local(lgrid, ldeform, lfrict, ratio, bf2, tf2, ueddy, sn, beta, inv_l, rd)
         rd = rd_over_dx(i, j)
         lgrid = sqrt(max(areaT(i, j), 0.0_wp))
         ldeform = lgrid*rd
         lfrict = 0.0_wp
         if (cdrag > 0.0_wp) lfrict = depth_tot(i, j)/cdrag

         ! gamma_b^2 = cd_scale^2 + 1/(1+Cb*Ldeform/Lfrict)^0.8 (floor).
         bf2 = cd_scale*cd_scale
         if (lfrict*cb > 0.0_wp) then
            ratio = ldeform/lfrict
            bf2 = bf2 + 1.0_wp/(1.0_wp + cb*ratio)**0.8_wp
         end if
         bf2 = max(bf2, min_gamma2)
         bottom_fac2(i, j) = bf2

         ! gamma_t^2 = 1/(1+Ct*Ldeform/Lfrict)^0.25 (floor).
         tf2 = 1.0_wp
         if (lfrict*ct > 0.0_wp) then
            ratio = ldeform/lfrict
            tf2 = 1.0_wp/(1.0_wp + ct*ratio)**0.25_wp
         end if
         tf2 = max(tf2, min_gamma2)
         barotr_fac2(i, j) = tf2

         ! Mixing length: harmonic sum of alpha-weighted scales.
         ueddy = sqrt(2.0_wp*max(0.0_wp, tf2*meke(i, j)))
         ! beta = |grad f|; centred f_centre differences scaled to a true
         ! gradient by the cell-centre inverse metrics (df/dx ~ (f_{i+1} -
         ! f_{i-1})*idxT/2).  When f_centre is unfilled (default) beta=0 ⇒
         ! Rhines weight inert (alpha_rhines default 0).
         beta = 0.0_wp
         if (i > 1 .and. i < nx) then
            beta = beta + (0.5_wp*(f_centre(i + 1, j) - f_centre(i - 1, j))*idxT(i, j))**2
         end if
         if (j > 1 .and. j < ny) then
            beta = beta + (0.5_wp*(f_centre(i, j + 1) - f_centre(i, j - 1))*idyT(i, j))**2
         end if
         beta = sqrt(beta)
         sn = 0.0_wp
         if (a_eady > 0.0_wp) sn = 0.25_wp*((sn_u(i, j) + sn_u(i + 1, j)) + &
                                            (sn_v(i, j) + sn_v(i, j + 1)))

         inv_l = meke_inv_lmix(ueddy, sn, beta, areaT(i, j), rd, &
                               depth_tot(i, j), cdrag, &
                               a_deform, a_rhines, a_eady, a_frict, a_grid)
         if (inv_l > 0.0_wp) then
            le(i, j) = 1.0_wp/inv_l
         else
            le(i, j) = 0.0_wp
         end if
      end do
   end subroutine meke_length_scales

   pure subroutine meke_source(nx, ny, bgsrc, gmcoeff, frcoeff, sdt, i_mass, &
                               gm_src, ke_diss, src, meke)
      !! Aggregate source `src = bgsrc + gmcoeff*I_mass*gm_src
      !! - frcoeff*I_mass*ke_diss` and the explicit bump `E += sdt*src`.
      !! `gmcoeff<0` ⇒ GM source off; `frcoeff<0` ⇒ frictional source off.
      !! `ke_diss` is the lateral-viscosity KE dissipation rate (≤0), so
      !! `-frcoeff*I_mass*ke_diss ≥ 0` is a mean→eddy source (0 ⇒ inert).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: bgsrc, gmcoeff, frcoeff, sdt
      real(wp), intent(in) :: i_mass(nx, ny)
      real(wp), intent(in) :: gm_src(nx, ny)
      real(wp), intent(in) :: ke_diss(nx, ny)
      real(wp), intent(inout) :: src(nx, ny)
      real(wp), intent(inout) :: meke(nx, ny)
      integer :: i, j
      real(wp) :: s

      do concurrent(j=1:ny, i=1:nx) local(s)
         s = bgsrc
         if (gmcoeff >= 0.0_wp) s = s + gmcoeff*i_mass(i, j)*gm_src(i, j)
         if (frcoeff >= 0.0_wp) s = s - frcoeff*i_mass(i, j)*ke_diss(i, j)
         src(i, j) = s
         meke(i, j) = meke(i, j) + sdt*s
      end do
   end subroutine meke_source

   pure subroutine meke_drag(nx, ny, sdt_damp, damping, cdrag, uscale, rho0, &
                             i_mass, bottom_fac2, u_bbl2, meke)
      !! Implicit (backward-Euler) bottom-drag half-step.
      !!   drag_rate = rho0*i_mass*sqrt(cdrag^2*(max(0,2*bf2*E)+u_bbl2+uscale^2))  [1/s]
      !!   damp_rate = damping + drag_rate*bf2 ; =0 where E<0
      !!   E <- E/(1 + sdt_damp*damp_rate)
      !! `rho0*i_mass` = rho0/(Sum_k rho_k*h_k) ~= 1/depth_tot [1/m], so
      !! drag_rate ~= cdrag*|U_d|/H -- the MOM6 `GV%H_to_RZ * I_mass` factor.
      !! Without it drag_rate is m^3/(kg*s), not a
      !! rate. `i_mass=0` on dry columns still gives `drag_rate=0`.
      !! `u_bbl2` is the resolved bed-layer speed² (MOM6 `drag_rate_visc`);
      !! it is 0 unless `use_bbl_drag` is set, so the default is bit-identical.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sdt_damp, damping, cdrag, uscale, rho0
      real(wp), intent(in) :: i_mass(nx, ny)
      real(wp), intent(in) :: bottom_fac2(nx, ny)
      real(wp), intent(in) :: u_bbl2(nx, ny)
      real(wp), intent(inout) :: meke(nx, ny)
      integer :: i, j
      real(wp) :: drag_rate, damp_rate, cd2, e

      cd2 = cdrag*cdrag
      do concurrent(j=1:ny, i=1:nx) local(drag_rate, damp_rate, e)
         e = meke(i, j)
         drag_rate = (rho0*i_mass(i, j))*sqrt(cd2*(max(0.0_wp, 2.0_wp*bottom_fac2(i, j)*e) &
                                                   + u_bbl2(i, j) + uscale*uscale))
         damp_rate = damping + drag_rate*bottom_fac2(i, j)
         if (e < 0.0_wp) damp_rate = 0.0_wp
         meke(i, j) = e/(1.0_wp + sdt_damp*damp_rate)
      end do
   end subroutine meke_drag

   pure subroutine meke_bbl_speed2(nx, ny, nz, u_face, v_face, u_bbl2)
      !! Resolved bed-layer (k=1, bottom-up) speed² at cell centres:
      !! `u_bbl2 = u_c² + v_c²` with `u_c = ½(u_face(i)+u_face(i+1))`,
      !! `v_c = ½(v_face(j)+v_face(j+1))` from the BED layer.  The bottom
      !! eddy velocity the MEKE drag law needs (MOM6 `drag_rate_visc`).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: u_face(nx + 1, ny, nz)
      real(wp), intent(in) :: v_face(nx, ny + 1, nz)
      real(wp), intent(inout) :: u_bbl2(nx, ny)
      integer :: i, j
      real(wp) :: u_c, v_c
      do concurrent(j=1:ny, i=1:nx) local(u_c, v_c)
         u_c = 0.5_wp*(u_face(i, j, 1) + u_face(i + 1, j, 1))
         v_c = 0.5_wp*(v_face(i, j, 1) + v_face(i, j + 1, 1))
         u_bbl2(i, j) = u_c*u_c + v_c*v_c
      end do
   end subroutine meke_bbl_speed2

   pure subroutine meke_kh_closure(nx, ny, khcoeff, barotr_fac2, meke, le, kh_diff)
      !! Derived diffusivity `kh = khcoeff*sqrt(2*max(0,gamma_t2*E))*Lmix`.
      !! `khcoeff<=0` ⇒ kh left at 0.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: khcoeff
      real(wp), intent(in) :: barotr_fac2(nx, ny)
      real(wp), intent(in) :: meke(nx, ny)
      real(wp), intent(in) :: le(nx, ny)
      real(wp), intent(out) :: kh_diff(nx, ny)
      integer :: i, j
      real(wp) :: ueddy

      do concurrent(j=1:ny, i=1:nx) local(ueddy)
         if (khcoeff > 0.0_wp) then
            ueddy = sqrt(2.0_wp*max(0.0_wp, barotr_fac2(i, j)*meke(i, j)))
            kh_diff(i, j) = khcoeff*ueddy*le(i, j)
         else
            kh_diff(i, j) = 0.0_wp
         end if
      end do
   end subroutine meke_kh_closure

   pure subroutine meke_ku_closure(nx, ny, backscatter, visc_coeff_ku, &
                                   meke, le, ku)
      !! Derived harmonic backscatter viscosity
      !! `ku = visc_coeff_ku*sqrt(2*max(0,E))*Lmix` (m²/s), matching MOM6
      !! `MEKE%Ku = MEKE_VISCOSITY_COEFF_KU*sqrt(2*MEKE)*Lmix`.  Unlike the
      !! kh closure (which carries the barotropic-mode factor `gamma_t2`
      !! inside the eddy velocity), MOM6's Ku uses the PLAIN `sqrt(2*MEKE)` —
      !! no vertical-structure factor — so `gamma_t2` is deliberately absent
      !! here (vertical structure `BS_struct = 1`; EBT/SQG deferred).
      !! Off ⇒ ku left at 0 (bit-identical seam).  Always ≥ 0; the SIGN of
      !! the momentum effect is set by the subtraction downstream.
      integer, intent(in) :: nx, ny
      logical, intent(in) :: backscatter
      real(wp), intent(in) :: visc_coeff_ku
      real(wp), intent(in) :: meke(nx, ny)
      real(wp), intent(in) :: le(nx, ny)
      real(wp), intent(out) :: ku(nx, ny)
      integer :: i, j
      real(wp) :: ueddy

      do concurrent(j=1:ny, i=1:nx) local(ueddy)
         if (backscatter .and. visc_coeff_ku /= 0.0_wp) then
            ueddy = sqrt(2.0_wp*max(0.0_wp, meke(i, j)))
            ku(i, j) = visc_coeff_ku*ueddy*le(i, j)
         else
            ku(i, j) = 0.0_wp
         end if
      end do
   end subroutine meke_ku_closure

   pure subroutine meke_backscatter_apply(grid, metrics, this, dt, ah_face_x, ah_face_y)
      !! Inject the MEKE harmonic backscatter into the per-face resolved
      !! harmonic viscosity (capability Gap 2, v1).  Subtracts a face-average
      !! of the cell-centred `ku` field from `ah_face_x`/`ah_face_y` so the
      !! NET coefficient `A_net = A_resolved − Ku` can go NEGATIVE — that
      !! negative viscosity is the energy return into the momentum tendency
      !! (the hvisc Laplacian kernel reads these same face fields).
      !!
      !! STABILITY: a negative explicit Laplacian viscosity AMPLIFIES
      !! grid-scale modes — a pure negative Laplacian is unconditionally
      !! unstable (`g = 1 + |A|·dt·k² > 1`).  The floor does NOT by itself
      !! make the operator stable; it bounds the negative mode's GROWTH
      !! RATE so a co-present POSITIVE biharmonic (Smag_AH / `nu_4` /
      !! `leith_biharm`) can dissipate the grid-scale mode the backscatter
      !! feeds.  The net magnitude is floored at the forward-Euler bound
      !! `|A|·dt·(idx²+idy²) ≤ backscatter_cfl·0.5`, i.e.
      !! `A_floor = −backscatter_cfl·0.5/(dt·(idx²+idy²))` per face
      !! (`backscatter_cfl = 0.8` safety; MOM6 `BACKSCATTER_UNDERBOUND`
      !! analogue).  Result: `ah_face <- max(A_resolved − Ku, A_floor)`.
      !! A biharmonic backstop is therefore MANDATORY when `backscatter`
      !! is on — enforced fail-loud at configure (`validate_config`).
      !!
      !! No-op (bit-identical) when `backscatter` is off or the slot is
      !! uninitialised.  GPU-resident: explicit-shape flat-impl dispatch.
      !! Run AFTER `ocean_lateral_mix_compute` (which fills `ah_face_*`) and
      !! BEFORE `ocean_horizontal_viscosity_compute_tendencies` consumes them.
      !! Uses the PREVIOUS thermo step's `ku` (a one-step lag, like the
      !! kh→KhTh feedback) — `meke_step` refreshes `ku` once per outer step.
      !!
      !! REQUIRES a flow-aware lateral closure (`&ocean_hvisc_nml
      !! lateral_closure` = `leith`/`smag` ⇒ `LMIX_*`): the hvisc Laplacian
      !! only reads `ah_face_*` when `closure /= LMIX_NONE`.  With the
      !! default scalar-`nu_h` path the modified faces are never consumed,
      !! so backscatter has no effect (it is still bit-identical when off).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_meke_t), intent(in) :: this
      real(wp), intent(in) :: dt
      real(wp), intent(inout) :: ah_face_x(:, :, :)
      real(wp), intent(inout) :: ah_face_y(:, :, :)
      integer :: nx, ny, nz

      if (.not. this%is_init) return
      if (.not. this%backscatter) return
      if (.not. allocated(this%ku)) return
      nx = grid%nx_total
      ny = grid%ny_total
      nz = size(ah_face_x, 3)
      if (this%nx_total /= nx .or. this%ny_total /= ny) return

      call meke_backscatter_apply_impl(nx, ny, nz, dt, BACKSCATTER_CFL, this%ku, &
                                       metrics%idxCu, metrics%idyCu, &
                                       metrics%idxCv, metrics%idyCv, &
                                       ah_face_x, ah_face_y)
   end subroutine meke_backscatter_apply

   pure subroutine meke_backscatter_apply_impl(nx, ny, nz, dt, cfl_safety, ku, &
                                               idxCu, idyCu, idxCv, idyCv, &
                                               ah_face_x, ah_face_y)
      !! Flat-impl: subtract the face-averaged `ku` from each per-face
      !! harmonic viscosity and floor the net at the CFL-stable minimum.
      !! Explicit-shape dummies for NVHPC stdpar (no descriptor walk).
      !! The backscatter coefficient is z-independent in v1 (`BS_struct=1`),
      !! so the cell-centred `ku(i,j)` is broadcast to every layer.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt, cfl_safety
      real(wp), intent(in) :: ku(nx, ny)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCu(nx + 1, ny)
      real(wp), intent(in) :: idxCv(nx, ny + 1), idyCv(nx, ny + 1)
      real(wp), intent(inout) :: ah_face_x(nx + 1, ny, nz), ah_face_y(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: ku_face, denom, a_floor, a_net

      ! u-faces (i-1/2, j): Ku averaged from the two adjacent T-cells
      ! (i-1, i).  Interior faces only; wall faces (i=1, i=nx+1) keep the
      ! resolved background (no backscatter at the closed boundary, matching
      ! the lateral-mix wall convention).
      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(ku_face, denom, a_floor, a_net)
         ku_face = 0.5_wp*(ku(i - 1, j) + ku(i, j))
         denom = dt*(idxCu(i, j)*idxCu(i, j) + idyCu(i, j)*idyCu(i, j))
         a_net = ah_face_x(i, j, k) - ku_face
         if (denom > 0.0_wp) then
            a_floor = -cfl_safety*0.5_wp/denom
            if (a_net < a_floor) a_net = a_floor
         end if
         ah_face_x(i, j, k) = a_net
      end do

      ! v-faces (i, j-1/2): Ku averaged from the two adjacent T-cells
      ! (j-1, j).  Interior faces only.
      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(ku_face, denom, a_floor, a_net)
         ku_face = 0.5_wp*(ku(i, j - 1) + ku(i, j))
         denom = dt*(idxCv(i, j)*idxCv(i, j) + idyCv(i, j)*idyCv(i, j))
         a_net = ah_face_y(i, j, k) - ku_face
         if (denom > 0.0_wp) then
            a_floor = -cfl_safety*0.5_wp/denom
            if (a_net < a_floor) a_net = a_floor
         end if
         ah_face_y(i, j, k) = a_net
      end do
   end subroutine meke_backscatter_apply_impl

   pure subroutine meke_feed_khth(nx, ny, khth_fac, khtr_fac, kh_diff, &
                                  khth_u, khth_v, khtr_u, khtr_v)
      !! Add the geometric mean of neighbour `kh_diff` into VarMix's per-face
      !! KhTh (and KhTr) base BEFORE GM's CFL clamp:
      !!   khth_u(i,j) += khth_fac*sqrt(kh(i-1,j)*kh(i,j))
      !! `khth_fac=0` ⇒ nothing added ⇒ bit-identical seam.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: khth_fac, khtr_fac
      real(wp), intent(in) :: kh_diff(nx, ny)
      real(wp), intent(inout) :: khth_u(nx + 1, ny)
      real(wp), intent(inout) :: khth_v(nx, ny + 1)
      real(wp), intent(inout) :: khtr_u(nx + 1, ny)
      real(wp), intent(inout) :: khtr_v(nx, ny + 1)
      integer :: i, j
      real(wp) :: gm_u, gm_v

      ! u-faces: interior i=2..nx pairs (i-1,i).
      do concurrent(j=1:ny, i=2:nx) local(gm_u)
         gm_u = sqrt(max(0.0_wp, kh_diff(i - 1, j))*max(0.0_wp, kh_diff(i, j)))
         khth_u(i, j) = khth_u(i, j) + khth_fac*gm_u
         khtr_u(i, j) = khtr_u(i, j) + khtr_fac*gm_u
      end do
      ! v-faces: interior j=2..ny pairs (j-1,j).
      do concurrent(j=2:ny, i=1:nx) local(gm_v)
         gm_v = sqrt(max(0.0_wp, kh_diff(i, j - 1))*max(0.0_wp, kh_diff(i, j)))
         khth_v(i, j) = khth_v(i, j) + khth_fac*gm_v
         khtr_v(i, j) = khtr_v(i, j) + khtr_fac*gm_v
      end do
   end subroutine meke_feed_khth

   pure subroutine meke_lateral(nx, ny, sdt, kh_bg, k4, khmeke_fac, &
                                kh_flux_enabled, dy_cu, dx_cv, idxCu, idyCv, &
                                iareaT, i_mass, mass, kh_diff, &
                                uflux, vflux, del2, meke)
      !! Harmonic-mass Laplacian diffusion of MEKE (+ optional biharmonic).
      !! Flux-form, conservative on a closed domain (interior faces only;
      !! array-edge faces carry zero flux).
      !!   Kh_u = max(0,kh_bg) + khmeke_fac*0.5*(kh_i+kh_{i+1}), CFL-capped 0.25
      !!   uflux = Kh_u*(dy_cu*idxCu)*[2 m_i m_{i+1}/(m_i+m_{i+1}+eps)]*(E_i-E_{i+1})
      !!   E += sdt*iareaT*i_mass*((uflux_{i-1}-uflux_i)+(vflux_{j-1}-vflux_j))
      !! Biharmonic: del2 = iareaT*(d uflux' + d vflux') with the bare-gradient
      !! flux uflux' = (dy_cu*idxCu)*(E_{i+1}-E_i); then a harmonic-mass flux
      !! of del2 with CFL cap 0.3 and E += that divergence (additive).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sdt, kh_bg, k4, khmeke_fac
      logical, intent(in) :: kh_flux_enabled
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny)
      real(wp), intent(in) :: idyCv(nx, ny + 1)
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: i_mass(nx, ny)
      real(wp), intent(in) :: mass(nx, ny)
      real(wp), intent(in) :: kh_diff(nx, ny)
      real(wp), intent(inout) :: uflux(nx + 1, ny)
      real(wp), intent(inout) :: vflux(nx, ny + 1)
      real(wp), intent(inout) :: del2(nx, ny)
      real(wp), intent(inout) :: meke(nx, ny)
      integer :: i, j
      real(wp) :: kh_u, kh_v, hm, geo, inv_max, k4_u, k4_v, mke

      ! ---------- Biharmonic (computed first; tendency added after diffusion). ----------
      if (k4 >= 0.0_wp) then
         ! bare-gradient flux into uflux/vflux workspaces (units m^2/s^2).
         do concurrent(j=1:ny, i=1:nx + 1)
            uflux(i, j) = 0.0_wp
         end do
         do concurrent(j=1:ny, i=2:nx)
            uflux(i, j) = (dy_cu(i, j)*idxCu(i, j))*(meke(i, j) - meke(i - 1, j))
         end do
         do concurrent(j=1:ny + 1, i=1:nx)
            vflux(i, j) = 0.0_wp
         end do
         do concurrent(j=2:ny, i=1:nx)
            vflux(i, j) = (dx_cv(i, j)*idyCv(i, j))*(meke(i, j) - meke(i, j - 1))
         end do
         do concurrent(j=1:ny, i=1:nx)
            del2(i, j) = iareaT(i, j)*((uflux(i + 1, j) - uflux(i, j)) + &
                                       (vflux(i, j + 1) - vflux(i, j)))
         end do
         ! harmonic-mass flux of del2 with K4 (CFL cap 0.3).
         do concurrent(j=1:ny, i=1:nx + 1)
            uflux(i, j) = 0.0_wp
         end do
         do concurrent(j=1:ny, i=2:nx) local(k4_u, geo, hm, inv_max)
            geo = dy_cu(i, j)*idxCu(i, j)
            inv_max = 64.0_wp*sdt*(geo*max(iareaT(i - 1, j), iareaT(i, j)))**2
            k4_u = k4
            if (k4_u*inv_max > 0.3_wp) k4_u = 0.3_wp/inv_max
            hm = 2.0_wp*mass(i - 1, j)*mass(i, j)/((mass(i - 1, j) + mass(i, j)) + MASS_NEGLECT)
            uflux(i, j) = (k4_u*geo)*hm*(del2(i, j) - del2(i - 1, j))
         end do
         do concurrent(j=1:ny + 1, i=1:nx)
            vflux(i, j) = 0.0_wp
         end do
         do concurrent(j=2:ny, i=1:nx) local(k4_v, geo, hm, inv_max)
            geo = dx_cv(i, j)*idyCv(i, j)
            inv_max = 64.0_wp*sdt*(geo*max(iareaT(i, j - 1), iareaT(i, j)))**2
            k4_v = k4
            if (k4_v*inv_max > 0.3_wp) k4_v = 0.3_wp/inv_max
            hm = 2.0_wp*mass(i, j - 1)*mass(i, j)/((mass(i, j - 1) + mass(i, j)) + MASS_NEGLECT)
            vflux(i, j) = (k4_v*geo)*hm*(del2(i, j) - del2(i, j - 1))
         end do
         do concurrent(j=1:ny, i=1:nx) local(mke)
            mke = sdt*(iareaT(i, j)*i_mass(i, j))* &
                  ((uflux(i, j) - uflux(i + 1, j)) + (vflux(i, j) - vflux(i, j + 1)))
            del2(i, j) = mke   ! stash the biharmonic tendency in del2
         end do
      end if

      ! ---------- Laplacian (harmonic-mass) diffusion. ----------
      if (kh_flux_enabled) then
         do concurrent(j=1:ny, i=1:nx + 1)
            uflux(i, j) = 0.0_wp
         end do
         do concurrent(j=1:ny, i=2:nx) local(kh_u, geo, hm, inv_max)
            geo = dy_cu(i, j)*idxCu(i, j)
            kh_u = max(0.0_wp, kh_bg) + khmeke_fac*0.5_wp*(kh_diff(i - 1, j) + kh_diff(i, j))
            inv_max = 2.0_wp*sdt*(geo*max(iareaT(i - 1, j), iareaT(i, j)))
            if (kh_u*inv_max > 0.25_wp) kh_u = 0.25_wp/inv_max
            hm = 2.0_wp*mass(i - 1, j)*mass(i, j)/((mass(i - 1, j) + mass(i, j)) + MASS_NEGLECT)
            uflux(i, j) = (kh_u*geo)*hm*(meke(i - 1, j) - meke(i, j))
         end do
         do concurrent(j=1:ny + 1, i=1:nx)
            vflux(i, j) = 0.0_wp
         end do
         do concurrent(j=2:ny, i=1:nx) local(kh_v, geo, hm, inv_max)
            geo = dx_cv(i, j)*idyCv(i, j)
            kh_v = max(0.0_wp, kh_bg) + khmeke_fac*0.5_wp*(kh_diff(i, j - 1) + kh_diff(i, j))
            inv_max = 2.0_wp*sdt*(geo*max(iareaT(i, j - 1), iareaT(i, j)))
            if (kh_v*inv_max > 0.25_wp) kh_v = 0.25_wp/inv_max
            hm = 2.0_wp*mass(i, j - 1)*mass(i, j)/((mass(i, j - 1) + mass(i, j)) + MASS_NEGLECT)
            vflux(i, j) = (kh_v*geo)*hm*(meke(i, j - 1) - meke(i, j))
         end do
         do concurrent(j=1:ny, i=1:nx) local(mke)
            mke = sdt*(iareaT(i, j)*i_mass(i, j))* &
                  ((uflux(i, j) - uflux(i + 1, j)) + (vflux(i, j) - vflux(i, j + 1)))
            meke(i, j) = meke(i, j) + mke
         end do
      end if

      ! add the biharmonic tendency (computed above, stashed in del2).
      if (k4 >= 0.0_wp) then
         do concurrent(j=1:ny, i=1:nx)
            meke(i, j) = meke(i, j) + del2(i, j)
         end do
      end if
   end subroutine meke_lateral

   pure function ocean_meke_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the MEKE slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_meke_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%meke) &
               + arr_bytes(this%kh_diff) &
               + arr_bytes(this%le) &
               + arr_bytes(this%ku) &
               + arr_bytes(this%i_mass) &
               + arr_bytes(this%depth_tot) &
               + arr_bytes(this%bottom_fac2) &
               + arr_bytes(this%barotr_fac2) &
               + arr_bytes(this%src) &
               + arr_bytes(this%uflux) &
               + arr_bytes(this%vflux) &
               + arr_bytes(this%del2) &
               + arr_bytes(this%mass_ws) &
               + arr_bytes(this%u_bbl2) &
               + arr_bytes(this%ke_diss_ws) &
               + arr_bytes(this%rd_ws) &
               + arr_bytes(this%f_centre) &
               + arr_bytes(this%sn_u_ws) &
               + arr_bytes(this%sn_v_ws) &
               + arr_bytes(this%baro_hu) &
               + arr_bytes(this%baro_hv)
   end function ocean_meke_bytes

end module rdb_ocean_meke
