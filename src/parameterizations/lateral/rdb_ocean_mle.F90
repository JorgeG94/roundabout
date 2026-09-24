!! Ocean mixed-layer-eddy (Fox-Kemper) restratification slot.
module rdb_ocean_mle
   !! Fox-Kemper, Ferrari & Hallberg (2008) submesoscale mixed-layer-eddy
   !! (MLE) restratification.  Submesoscale eddies slump lateral buoyancy
   !! fronts in the surface mixed layer via an overturning streamfunction
   !!     Psi(z) = Ce * (H_ml^2 / |f|) * (grad b_bar x z_hat) * mu(z)
   !! injected as extra ML-confined per-layer mass transports `uhml`/`vhml`
   !! added to `ms%mass_flux_{x,y}_layer` BEFORE the continuity divergence +
   !! PPM tracer advection — never touching the velocity fields.  The
   !! per-layer weights `a(k)` satisfy `sum_k a(k) = mu(0) - mu(-1) = 0`
   !! exactly (a closed overturning cell ⇒ mass/tracer conservative).
   !!
   !! Vertical convention: k=1 bed, k=nz surface.  The ML band walks from
   !! the surface down; the FK sigma coordinate is 0 at the surface,
   !! -1 at the ML base.  `H_ml` is the EPBL `epbl%mld`.  Timescale forms
   !! (`use_mom_mixrate`): bare `Ce/max(|f|,f_floor)` (default), or FK11
   !! momentum-mixrate (production-recommended; suppresses restratification
   !! under vigorous mixing).  Default off ⇒ bit-identical.
   !!
   !! References: Fox-Kemper, Ferrari & Hallberg (2008); Fox-Kemper et al.
   !! (2011).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_epbl, only: ocean_epbl_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
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

   public :: ocean_mle_t
   public :: mle_compute_transports
   public :: mle_fold_x, mle_fold_y
   public :: mle_mu_shape, mle_layer_weights

   real(wp), parameter :: H_NEGLECT = 1.0e-30_wp
      !! Empty-column guard (spec 3.1 / 3.4).
   real(wp), parameter :: MLE_H_AVAIL_MIN = 1.0e-6_wp
      !! Floor for the per-layer transport availability cap (m): a donor
      !! layer may never be drained below this by the FK overturning.
      !! Matched to the continuity / windowed-drain `h_min` (1e-6 m) so the
      !! cap and the drain limiter protect the same minimum thickness.
   real(wp), parameter :: VONKAR = 0.41_wp
      !! von Karman constant for the FK11 momentum-mixrate form.
   real(wp), parameter :: PI = 3.14159265358979323846_wp

   type :: ocean_mle_t
      !! Fox-Kemper mixed-layer-eddy restratification state.  All fields
      !! default to the inert (`enable=.false.`) configuration so an
      !! ocean run that never sets `&ocean_foxkemper_nml` is bit-identical.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; tracks GPU attachment.
      logical :: enable = .false.
         !! Master switch.  Off => `mle_compute_transports` is a no-op and
         !! the folds add nothing => bit-identity preserved.
      real(wp) :: ce = 0.0625_wp
         !! FK08 coefficient Ce (typical 0.06-0.08).
      real(wp) :: f_floor = 1.0e-5_wp
         !! |f| regularisation floor (1/s); ~|f| at 4 degN.
      real(wp) :: mld_decay_time = 0.0_wp
         !! Running-mean MLD filter time-scale (s).  0 ⇒ filter off
         !! (instantaneous EPBL MLD, bit-identical).  Positive ⇒ the MLD
         !! resets instantly to a deeper value but decays over this scale
         !! when the diagnosed MLD retreats, bounding Psi ~ MLD² growth
         !! under a thinning/oscillating boundary layer.
      real(wp) :: tail_dh = 0.0_wp
         !! mu cubic-tail extension below the ML base.  Default 0 ⇒ exact
         !! FK08 mu (the tail itself is deferred).
      logical :: use_mom_mixrate = .false.
         !! Use the FK11 momentum-mixrate timescale instead of the bare
         !! Ce/|f| floor form.  Production-recommended; default off.
      logical :: resolution_taper = .false.
         !! Resolution-function taper hook.  Hard config error if on
         !! without the resolution function (deferred).
      logical :: use_bodner = .false.
         !! Bodner et al. (2023) frontogenesis-arrest variant.  When true the
         !! `Ce/|f|` timescale is replaced by `ts_bod = Cr·Δs·|f|·h/w'u'` (a
         !! timescale [s], so it drops into the same uDml formula); the
         !! frontal-arrest length enters inline as `|f|·h/w'u'`.  Default off.
      real(wp) :: cr = 0.0_wp
         !! Bodner efficiency coefficient `Cr` (nondim).
      real(wp) :: bodner_mstar = 0.5_wp
         !! Mechanical (u*) weight in `w'u' = (mstar·u*³+nstar·w*³)^{2/3}`.
      real(wp) :: bodner_nstar = 0.066_wp
         !! Convective (w*) weight; `w*³ = max(0,-b0)·mld` from `epbl%b0`.
      real(wp) :: min_wstar2 = 1.0e-24_wp
         !! Floor on `w'u'` (m²/s²), pure 1/0 armour.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- 2D ML diagnostics (cell centres) ----
      real(wp), allocatable :: b_ml(:, :)
         !! ML-averaged buoyancy b_bar (m/s^2), `(nx,ny)`.
      real(wp), allocatable :: htot_ml(:, :)
         !! Accumulated ML thickness on the grid (<= mld) (m), `(nx,ny)`.
      real(wp), allocatable :: mld_filtered(:, :)
         !! Running-mean filtered MLD (m), `(nx,ny)`; persistent across
         !! steps.  Only used when `mld_decay_time > 0`.  Lazily seeded
         !! from the first instantaneous MLD (the `< 0` sentinel marks the
         !! unseeded state so the running mean starts from the true value
         !! rather than a spurious zero — MOM6 seeds from a restart field).

      ! ---- Per-layer FK transports (faces) ----
      ! k=1 bed, k=nz surface.  These are added to the continuity mass
      ! fluxes before the divergence; sum_k over each face is ~0.
      real(wp), allocatable :: uhml(:, :, :)
         !! FK x-face transport (m^3/s), `(nx+1,ny,nz_ml)`.  Filled once per
         !! outer step at thermo cadence, then folded into continuity every
         !! dynamics call via `mle_fold_x`.  At `dt_therm_ratio > 1` the
         !! stale values are re-folded on non-thermo steps (over-applies FK;
         !! still conservative since sum_k a(k) = 0).
      real(wp), allocatable :: vhml(:, :, :)
         !! FK y-face transport (m^3/s), `(nx,ny+1,nz_ml)`; same limitation.
   contains
      procedure :: init => ocean_mle_init
      procedure :: destroy => ocean_mle_destroy
      procedure :: enter_data => ocean_mle_enter_data
      procedure :: exit_data => ocean_mle_exit_data
      procedure, non_overridable :: bytes => ocean_mle_bytes
   end type ocean_mle_t

contains

   subroutine ocean_mle_init(this, grid, nz_ml)
      !! Allocate the 2D ML diagnostics + per-layer face transports.
      !! Always allocates (configure runs after init); the off-state
      !! footprint is two 2D + two face-shaped 3D arrays.
      class(ocean_mle_t), intent(inout) :: this
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

      allocate (this%b_ml(nx, ny), source=0.0_wp)
      allocate (this%htot_ml(nx, ny), source=0.0_wp)
      ! -1 sentinel = unseeded; the first filtered step copies the
      ! instantaneous MLD into the running mean.
      allocate (this%mld_filtered(nx, ny), source=-1.0_wp)
      allocate (this%uhml(nx + 1, ny, nz), source=0.0_wp)
      allocate (this%vhml(nx, ny + 1, nz), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_mle_init

   subroutine ocean_mle_destroy(this)
      class(ocean_mle_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%b_ml)) deallocate (this%b_ml)
      if (allocated(this%htot_ml)) deallocate (this%htot_ml)
      if (allocated(this%mld_filtered)) deallocate (this%mld_filtered)
      if (allocated(this%uhml)) deallocate (this%uhml)
      if (allocated(this%vhml)) deallocate (this%vhml)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_mle_destroy

   subroutine ocean_mle_enter_data(this)
      ! select-type wrapper delegating to the non-poly `_impl` — the
      ! polymorphic dummy's stack descriptor must not be the map root
      ! (AMD libomptarget cross-slot overlap; see ocean_lateral_mix).
      ! Slot header + scalar knobs ride the orchestrator's root
      ! copyin(state) in ocean_state_enter_data.
      class(ocean_mle_t), intent(inout) :: this
      select type (this)
      type is (ocean_mle_t)
         call ocean_mle_enter_data_impl(this)
      end select
   end subroutine ocean_mle_enter_data

   subroutine ocean_mle_enter_data_impl(this)
      type(ocean_mle_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%b_ml, this%htot_ml, this%mld_filtered)
      !$acc enter data copyin(this%uhml, this%vhml)
   end subroutine ocean_mle_enter_data_impl

   subroutine ocean_mle_exit_data(this)
      class(ocean_mle_t), intent(inout) :: this
      select type (this)
      type is (ocean_mle_t)
         call ocean_mle_exit_data_impl(this)
      end select
   end subroutine ocean_mle_exit_data

   subroutine ocean_mle_exit_data_impl(this)
      type(ocean_mle_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%vhml, this%uhml)
      !$acc exit data delete(this%mld_filtered, this%htot_ml, this%b_ml)
   end subroutine ocean_mle_exit_data_impl

   ! =================================================================
   ! Pure device helpers (callable from do concurrent)
   ! =================================================================

   pure function mle_mu_shape(sigma) result(mu)
      !$acc routine seq
      !! FK08 second-order vertical structure function mu(sigma).
      !! sigma in [-1,0]: 0 = surface interface, -1 = ML base.  mu = 0 at
      !! sigma=0 (surface) and sigma <= -1 (ML base); peaks near sigma=-0.5.
      !!
      !!     mu(sigma) = max(0, (1-(2 sigma+1)^2) * (1 + (5/21)(2 sigma+1)^2))
      !!
      !! Reference: FK08 eq. 21 / FK11 eq. 5.
      real(wp), intent(in) :: sigma
      real(wp) :: mu, x
      x = 2.0_wp*sigma + 1.0_wp          ! maps sigma in [-1,0] -> x in [-1,1]
      mu = (1.0_wp - x*x)*(1.0_wp + (5.0_wp/21.0_wp)*x*x)
      mu = max(0.0_wp, mu)
   end function mle_mu_shape

   pure function mle_timescale(f_abs, ustar, h_vel, ce, f_floor, use_mom_mixrate) result(ts)
      !$acc routine seq
      !! FK restratification timescale [s].
      !!
      !! Bare (default, T1 analytic gate):
      !!     ts = Ce / max(|f|, f_floor)
      !! reproduces FK08 Psi_max = Ce*H^2*|grad b|/|f| exactly.
      !!
      !! FK11 momentum-mixrate (PRODUCTION-RECOMMENDED, `use_mom_mixrate`):
      !!     mr = (vonKar * pi^2) * u*^2 / (|f|*H_vel^2 + 4*(H_vel+h_neg)*u*)
      !!     ts = Ce * (|f| + 2*mr) / (|f|^2 + mr^2)
      !! The mr term regularises 1/|f| as f -> 0 and suppresses restrat
      !! under vigorous mixing.  The bare form is the mr -> 0 limit.
      real(wp), intent(in) :: f_abs, ustar, h_vel, ce, f_floor
      logical, intent(in) :: use_mom_mixrate
      real(wp) :: ts, mr, f_eff
      if (use_mom_mixrate) then
         mr = (VONKAR*PI*PI)*ustar*ustar/ &
              (f_abs*h_vel*h_vel + 4.0_wp*(h_vel + H_NEGLECT)*ustar + H_NEGLECT)
         ts = ce*(f_abs + 2.0_wp*mr)/(f_abs*f_abs + mr*mr + H_NEGLECT)
      else
         f_eff = max(f_abs, f_floor)
         ts = ce/f_eff
      end if
   end function mle_timescale

   pure function mle_bodner_timescale(f_abs, ustar, h_vel, b0_face, ds, &
                                      cr, mstar, nstar, min_wstar2) result(ts)
      !$acc routine seq
      !! Bodner et al. (2023) frontogenesis-arrest MLE timescale [s]:
      !!     ts = Cr · ds · |f| · h / w'u'
      !! where the frontal-arrest length enters inline as `|f|·h/w'u'` and
      !!     w'u' = max( (mstar·u*³ + nstar·w*³)^(2/3), min_wstar2 ),
      !!     w*³   = max(0, -b0)·h            (destabilizing buoyancy flux only)
      !! `ds = sqrt(0.5(dx²+dy²))` is the grid-scale front width.  Because the
      !! product is dimensionally a TIME, this drops straight into the FK
      !! transport form `uDml = ts·dyCu·idxCu·db·H²` (the swap that turns
      !! classic Fox-Kemper into Bodner).  `w'u'` is floored so a quiescent,
      !! unforced column (u*→0, b0→0) never divides by zero.
      real(wp), intent(in) :: f_abs, ustar, h_vel, b0_face, ds
      real(wp), intent(in) :: cr, mstar, nstar, min_wstar2
      real(wp) :: ts, wstar3, wpup
      wstar3 = max(0.0_wp, -b0_face)*h_vel
      wpup = max((mstar*ustar**3 + nstar*wstar3)**(2.0_wp/3.0_wp), min_wstar2)
      ts = cr*ds*f_abs*h_vel/wpup
   end function mle_bodner_timescale

   pure subroutine mle_layer_weights(h_face, nz, h_vel, a)
      !$acc routine seq
      !! Per-layer transport weights a(k) = mu(sigma_top) - mu(sigma_bot),
      !! walking surface (k=nz) -> bed (k=1).  `sum_k a(k) = mu(0)-mu(-1) = 0`
      !! (closed cell => conservation).  Layers below the ML base get
      !! sigma <= -1 => mu=0 on both interfaces => a(k)=0 (ML-confinement).
      !!
      !! `h_vel` is clamped to the column total inside the caller so the
      !! sigma walk reaches exactly -1 over the actual depth (guards the
      !! mld > column-depth case; mu(-1)=0 keeps sum a(k)=0 regardless).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_face(nz)   ! k=1 bed .. k=nz surface
      real(wp), intent(in) :: h_vel
      real(wp), intent(out) :: a(nz)       ! k=1 bed .. k=nz surface
      integer :: k
      real(wp) :: ih_tot, sigma_top, sigma_bot, mu_top, mu_bot
      ih_tot = 1.0_wp/(h_vel + H_NEGLECT)
      sigma_top = 0.0_wp                   ! sigma=0 at surface (above k=nz)
      ! Walk surface -> bed.
      do k = nz, 1, -1
         sigma_bot = sigma_top - h_face(k)*ih_tot   ! more negative going down
         mu_top = mle_mu_shape(sigma_top)
         mu_bot = mle_mu_shape(sigma_bot)
         a(k) = mu_top - mu_bot
         sigma_top = sigma_bot
      end do
   end subroutine mle_layer_weights

   ! =================================================================
   ! Compute kernel
   ! =================================================================

   subroutine mle_compute_transports(grid, metrics, this, ms, epbl, ss, dt_limit, bc)
      !! Fill `uhml`/`vhml` (m^3/s) with the FK MLE overturning transport.
      !! Run once per outer step at thermo cadence, before the continuity
      !! divergence.  Steps: (1) b_ml + htot_ml at cell centres (surface→bed
      !! band to mld, partial-weight the straddling layer); (2) uDml/vDml at
      !! faces from grad b_bar, timescale, H_vel²; (2b) optional per-layer
      !! availability cap (a scalar shrink keeping sum_k a(k)=0); (3) fold
      !! the mu profile a(k) → uhml/vhml.  No-op when `enable=.false.` or
      !! the slot / state arrays are absent.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_mle_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_epbl_t), intent(in) :: epbl
      type(ocean_surface_stress_t), intent(in), optional :: ss
      real(wp), intent(in), optional :: dt_limit
         !! Window (s) the FK transport is integrated over (= `dt_therm`).
         !! When present and > 0, the per-layer availability cap bounds the
         !! transport so the windowed drain's `hprev` reconstruction stays
         !! non-negative in thin surface layers.  Absent / ≤ 0 ⇒ no cap.
         !! PR-8: this cap is UNCONDITIONAL in production — both production
         !! call sites (`rdb_ocean_dyn.F90`) always pass
         !! `dt_limit=dyn%therm_dt(dt)` > 0, so `do_limit` is always
         !! `.true.` on the live path.  The former `&ocean_foxkemper_nml
         !! apply_cfl_limit` knob was deleted as dead/incoherent — it could
         !! only ever be set to a thing the code already always does; do
         !! not re-add a knob that toggles this cap.
      type(ocean_bc_state_t), intent(in), optional :: bc
         !! Per-edge OBC tags.  Masks the FK transport on closed
         !! (non-periodic) physical wall faces so no MLE overturning crosses
         !! a land boundary.  Absent ⇒ array-edge zeroing only.

      integer :: i, j, k, nx, ny, nz, nghost, nx_phys, ny_phys
      real(wp) :: ce_l, f_floor_l, rho0_l, g_over_rho0
      real(wp) :: cr_l, mstar_l, nstar_l, minw2_l
      logical :: use_mr_l, use_bodner_l, has_ustar, do_limit, do_filter, n_seam
      real(wp) :: h_remain, w, htot, rho_int
      real(wp) :: db, h_vel, f_abs, ustar, ts, uDml, vDml, i4dt, h_av
      real(wp) :: a_stack(NZ_STACK_MAX), hf_stack(NZ_STACK_MAX)
      real(wp) :: a_fac, b_fac, mld_inst, mld_use

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. allocated(ms%rho_layer)) return
      if (.not. allocated(ms%h_layer)) return
      if (.not. allocated(epbl%mld)) return
      if (.not. allocated(epbl%f_centre)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nghost = grid%nghost
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys

      ce_l = this%ce
      f_floor_l = this%f_floor
      use_mr_l = this%use_mom_mixrate
      use_bodner_l = this%use_bodner
      cr_l = this%cr
      mstar_l = this%bodner_mstar
      nstar_l = this%bodner_nstar
      minw2_l = this%min_wstar2
      ! Bodner needs the surface buoyancy flux (epbl%b0); if EPBL didn't
      ! persist it, there is nothing to restratify with -> no-op.
      if (use_bodner_l .and. .not. allocated(epbl%b0)) return
      rho0_l = epbl%rho0
      g_over_rho0 = GRAVITY/rho0_l
      ! Availability-cap setup (FK transport CFL limiter).  I4dt = 1/(4*dt)
      ! so a single donor face can evacuate at most 1/4 of the donor layer
      ! volume over the window — guarantees the windowed-drain hprev
      ! reconstruction stays positive in thin layers.
      do_limit = .false.
      i4dt = 0.0_wp
      if (present(dt_limit)) then
         if (dt_limit > 0.0_wp) then
            do_limit = .true.
            i4dt = 1.0_wp/(4.0_wp*dt_limit)
         end if
      end if
      ! u* only needed for the FK11 mixrate form; surface stress is a
      ! Pa wind stress tau, u* = sqrt(|tau|/rho0).  Bare form ignores it.
      has_ustar = .false.
      if ((use_mr_l .or. use_bodner_l) .and. present(ss)) has_ustar = .true.

      ! ---- Running-mean MLD filter ----------
      ! Psi ~ MLD^2, so a thinning/oscillating EPBL MLD injects spiky
      ! transport.  Damp with a running mean that resets instantly to a
      ! deeper MLD but decays over `mld_decay_time` when it retreats:
      !     aFac = T/(dt+T),  bFac = dt/(dt+T)
      !     MLD_filt = max( MLD, bFac*MLD + aFac*MLD_filt )
      ! `dt` = the FK call cadence (= `dt_limit`).  Filter off (default)
      ! reads the instantaneous MLD ⇒ bit-identical.  Fox-Kemper et al. (2011).
      do_filter = .false.
      if (this%mld_decay_time > 0.0_wp .and. present(dt_limit)) then
         if (dt_limit > 0.0_wp) do_filter = .true.
      end if
      if (do_filter) then
         a_fac = this%mld_decay_time/(dt_limit + this%mld_decay_time)
         b_fac = dt_limit/(dt_limit + this%mld_decay_time)
         do concurrent(j=1:ny, i=1:nx) local(mld_inst)
            mld_inst = epbl%mld(i, j)
            if (this%mld_filtered(i, j) < 0.0_wp) then
               ! Unseeded: start the running mean from the true MLD.
               this%mld_filtered(i, j) = mld_inst
            else
               this%mld_filtered(i, j) = max(mld_inst, &
                                             b_fac*mld_inst + a_fac*this%mld_filtered(i, j))
            end if
         end do
      end if

      ! ---- 1. ML-averaged buoyancy + clamped MLD at cell centres ----
      ! Walk k=nz (surface) -> bed; partial-weight the layer straddling
      ! the MLD base so htot reaches mld exactly.  b = -(g/rho0)*rho_bar.
      ! `mld_use` is the filtered MLD when the decay-time filter is on,
      ! else the instantaneous EPBL MLD (bit-identical legacy path).
      do concurrent(j=1:ny, i=1:nx) local(k, h_remain, w, htot, rho_int, mld_use)
         if (do_filter) then
            mld_use = this%mld_filtered(i, j)
         else
            mld_use = epbl%mld(i, j)
         end if
         htot = 0.0_wp
         rho_int = 0.0_wp
         do k = nz, 1, -1
            h_remain = mld_use - htot
            if (h_remain <= 0.0_wp) exit
            w = min(ms%h_layer(i, j, k), h_remain)   ! partial weight
            htot = htot + w
            rho_int = rho_int + ms%rho_layer(i, j, k)*w
         end do
         this%htot_ml(i, j) = htot
         this%b_ml(i, j) = -g_over_rho0*(rho_int/(htot + H_NEGLECT))
      end do

      ! ---- 2+3. u-face transports.  Interior faces i=2..nx; wall faces
      ! (i=1, i=nx+1) get zero transport so FK injects no flux through
      ! walls (mirrors the continuity wall convention).
      ! uDml = timescale * dy_cu * (b_E - b_W) * H_vel^2.
      ! Positive (b_E - b_W) (light to the east) => positive uDml; the
      ! surface a(k)<0 then drives SURFACE transport WESTWARD toward the
      ! dense column => the front slumps / restratifies (see module head).
      do concurrent(k=1:nz, j=1:ny)
         this%uhml(1, j, k) = 0.0_wp
         this%uhml(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=2:nx) &
         local(k, db, h_vel, f_abs, ustar, ts, uDml, h_av, a_stack, hf_stack)
         db = this%b_ml(i, j) - this%b_ml(i - 1, j)        ! b_E - b_W
         h_vel = 0.5_wp*(this%htot_ml(i - 1, j) + this%htot_ml(i, j))
         f_abs = abs(0.5_wp*(epbl%f_centre(i - 1, j) + epbl%f_centre(i, j)))
         ustar = 0.0_wp
         if (has_ustar) ustar = mle_face_ustar_x(ss, rho0_l, i, j)
         if (use_bodner_l) then
            ts = mle_bodner_timescale(f_abs, ustar, h_vel, &
                                      0.5_wp*(epbl%b0(i - 1, j) + epbl%b0(i, j)), &
                                      sqrt(0.5_wp*(metrics%dxCu(i, j)**2 + metrics%dyCu(i, j)**2)), &
                                      cr_l, mstar_l, nstar_l, minw2_l)
         else
            ts = mle_timescale(f_abs, ustar, h_vel, ce_l, f_floor_l, use_mr_l)
         end if
         ! db is the raw buoyancy DIFFERENCE b_E - b_W; idxCu = 1/dxCu turns
         ! it into the gradient db/dx, so dy_cu*idxCu = dyCu/dxCu is the
         ! transport aspect ratio (FK streamfunction * face width). Omitting
         ! idxCu made uDml ~dxCu too large (the FK over-amplification bug).
         uDml = ts*metrics%dy_cu(i, j)*metrics%idxCu(i, j)*db*h_vel*h_vel
         do k = 1, nz
            hf_stack(k) = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
         end do
         call mle_layer_weights(hf_stack, nz, h_vel, a_stack)
         ! Per-layer availability cap (MOM6).  Donor is the WEST cell (i-1)
         ! for positive transport a(k)*uDml > 0, the EAST cell (i) for
         ! negative.  Shrinking the scalar uDml preserves sum_k a(k)=0.
         if (do_limit) then
            do k = 1, nz
               if (a_stack(k)*uDml > 0.0_wp) then
                  h_av = max(i4dt*metrics%areaT(i - 1, j)* &
                             (ms%h_layer(i - 1, j, k) - MLE_H_AVAIL_MIN), 0.0_wp)
                  if (a_stack(k)*uDml > h_av) uDml = h_av/a_stack(k)
               else if (a_stack(k)*uDml < 0.0_wp) then
                  h_av = max(i4dt*metrics%areaT(i, j)* &
                             (ms%h_layer(i, j, k) - MLE_H_AVAIL_MIN), 0.0_wp)
                  if (-a_stack(k)*uDml > h_av) uDml = -h_av/a_stack(k)
               end if
            end do
         end if
         do k = 1, nz
            this%uhml(i, j, k) = a_stack(k)*uDml
         end do
      end do

      ! ---- 2+3. v-face transports.  Mirror; wall faces j=1, j=ny+1 = 0.
      do concurrent(k=1:nz, i=1:nx)
         this%vhml(i, 1, k) = 0.0_wp
         this%vhml(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(j=2:ny, i=1:nx) &
         local(k, db, h_vel, f_abs, ustar, ts, vDml, h_av, a_stack, hf_stack)
         db = this%b_ml(i, j) - this%b_ml(i, j - 1)        ! b_N - b_S
         h_vel = 0.5_wp*(this%htot_ml(i, j - 1) + this%htot_ml(i, j))
         f_abs = abs(0.5_wp*(epbl%f_centre(i, j - 1) + epbl%f_centre(i, j)))
         ustar = 0.0_wp
         if (has_ustar) ustar = mle_face_ustar_y(ss, rho0_l, i, j)
         if (use_bodner_l) then
            ts = mle_bodner_timescale(f_abs, ustar, h_vel, &
                                      0.5_wp*(epbl%b0(i, j - 1) + epbl%b0(i, j)), &
                                      sqrt(0.5_wp*(metrics%dxCv(i, j)**2 + metrics%dyCv(i, j)**2)), &
                                      cr_l, mstar_l, nstar_l, minw2_l)
         else
            ts = mle_timescale(f_abs, ustar, h_vel, ce_l, f_floor_l, use_mr_l)
         end if
         ! idyCv = 1/dyCv turns the raw db = b_N - b_S into db/dy; the
         ! dx_cv*idyCv = dxCv/dyCv aspect ratio mirrors the u-face above.
         vDml = ts*metrics%dx_cv(i, j)*metrics%idyCv(i, j)*db*h_vel*h_vel
         do k = 1, nz
            hf_stack(k) = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
         end do
         call mle_layer_weights(hf_stack, nz, h_vel, a_stack)
         ! Per-layer availability cap.  Donor is the SOUTH cell (i,j-1) for
         ! positive transport, the NORTH cell (i,j) for negative.
         if (do_limit) then
            do k = 1, nz
               if (a_stack(k)*vDml > 0.0_wp) then
                  h_av = max(i4dt*metrics%areaT(i, j - 1)* &
                             (ms%h_layer(i, j - 1, k) - MLE_H_AVAIL_MIN), 0.0_wp)
                  if (a_stack(k)*vDml > h_av) vDml = h_av/a_stack(k)
               else if (a_stack(k)*vDml < 0.0_wp) then
                  h_av = max(i4dt*metrics%areaT(i, j)* &
                             (ms%h_layer(i, j, k) - MLE_H_AVAIL_MIN), 0.0_wp)
                  if (-a_stack(k)*vDml > h_av) vDml = -h_av/a_stack(k)
               end if
            end do
         end if
         do k = 1, nz
            this%vhml(i, j, k) = a_stack(k)*vDml
         end do
      end do

      ! ---- Physical closed-wall face mask -------
      ! No MLE transport may cross a closed (land) boundary.  Array-edge
      ! zeroing above misses the PHYSICAL wall faces (at i=nghost+1 /
      ! i=nghost+nx_phys+1, interior to the array); a nonzero uhml/vhml
      ! there leaks tracer across the wall (worst under the windowed drain).
      ! Zero the FK transport on every non-periodic physical edge, mirroring
      ! the continuity wall convention.  Periodic edges are seams, not walls.
      if (present(bc)) then
         if (.not. bc%periodic_x) then
            do concurrent(k=1:nz, j=1:ny)
               this%uhml(nghost + 1, j, k) = 0.0_wp
               this%uhml(nghost + nx_phys + 1, j, k) = 0.0_wp
            end do
         end if
         if (.not. bc%periodic_y) then
            ! A tripolar north fold is a seam too: its fold-line face keeps
            ! the FK transport (projected antisymmetric with the resolved
            ! mass flux in `continuity_tracer_step_split`).
            n_seam = bc%north_fold   ! host scalar: never deref bc on device
            do concurrent(k=1:nz, i=1:nx)
               this%vhml(i, nghost + 1, k) = 0.0_wp
               if (.not. n_seam) this%vhml(i, nghost + ny_phys + 1, k) = 0.0_wp
            end do
         end if
      end if
   end subroutine mle_compute_transports

   pure function mle_face_ustar_x(ss, rho0, i, j) result(ustar)
      !$acc routine seq
      !! Friction velocity u* = sqrt(|tau|/rho0) at the u-face from the
      !! surface wind stress, averaged onto the face.  Only used by the
      !! FK11 mixrate form.  Returns 0 if stress fields are absent.
      type(ocean_surface_stress_t), intent(in) :: ss
      real(wp), intent(in) :: rho0
      integer, intent(in) :: i, j
      real(wp) :: ustar, tx, ty
      ustar = 0.0_wp
      if (.not. allocated(ss%tau_x) .or. .not. allocated(ss%tau_y)) return
      ! tau_x already lives on the u-face (i,j); tau_y on v-faces is
      ! averaged from the four neighbours straddling this u-face.
      tx = ss%tau_x(i, j)
      ty = 0.25_wp*(ss%tau_y(i - 1, j) + ss%tau_y(i, j) + &
                    ss%tau_y(i - 1, j + 1) + ss%tau_y(i, j + 1))
      ustar = sqrt(sqrt(tx*tx + ty*ty)/rho0)
   end function mle_face_ustar_x

   pure function mle_face_ustar_y(ss, rho0, i, j) result(ustar)
      !$acc routine seq
      !! u* at the v-face; mirror of `mle_face_ustar_x`.
      type(ocean_surface_stress_t), intent(in) :: ss
      real(wp), intent(in) :: rho0
      integer, intent(in) :: i, j
      real(wp) :: ustar, tx, ty
      ustar = 0.0_wp
      if (.not. allocated(ss%tau_x) .or. .not. allocated(ss%tau_y)) return
      ! tau_y already lives on the v-face (i,j); tau_x on u-faces is
      ! averaged from the four neighbours straddling this v-face.
      ty = ss%tau_y(i, j)
      tx = 0.25_wp*(ss%tau_x(i, j - 1) + ss%tau_x(i + 1, j - 1) + &
                    ss%tau_x(i, j) + ss%tau_x(i + 1, j))
      ustar = sqrt(sqrt(tx*tx + ty*ty)/rho0)
   end function mle_face_ustar_y

   ! =================================================================
   ! Fold into continuity mass fluxes (called inside the split step)
   ! =================================================================

   pure subroutine mle_fold_x(this, mass_flux_x_layer, nx1, ny, nz)
      !! Add the FK x-face transport into the per-layer zonal mass flux,
      !! AFTER `continuity_zonal_flux` fills it and BEFORE the zonal
      !! tracer advect / divergence — so the augmented flux transports
      !! both h and tracers (conservative; no velocity touched).
      !! No-op when disabled.
      type(ocean_mle_t), intent(in) :: this
      integer, intent(in) :: nx1, ny, nz
      real(wp), intent(inout) :: mass_flux_x_layer(nx1, ny, nz)
      integer :: i, j, k
      if (.not. this%is_init) return
      if (.not. this%enable) return
      do concurrent(k=1:nz, j=1:ny, i=1:nx1)
         mass_flux_x_layer(i, j, k) = mass_flux_x_layer(i, j, k) + this%uhml(i, j, k)
      end do
   end subroutine mle_fold_x

   pure subroutine mle_fold_y(this, mass_flux_y_layer, nx, ny1, nz)
      !! Add the FK y-face transport into the per-layer meridional mass
      !! flux.  Mirror of `mle_fold_x`.  No-op when disabled.
      type(ocean_mle_t), intent(in) :: this
      integer, intent(in) :: nx, ny1, nz
      real(wp), intent(inout) :: mass_flux_y_layer(nx, ny1, nz)
      integer :: i, j, k
      if (.not. this%is_init) return
      if (.not. this%enable) return
      do concurrent(k=1:nz, j=1:ny1, i=1:nx)
         mass_flux_y_layer(i, j, k) = mass_flux_y_layer(i, j, k) + this%vhml(i, j, k)
      end do
   end subroutine mle_fold_y

   pure function ocean_mle_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the MLE / Fox-Kemper slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_mle_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%b_ml) &
               + arr_bytes(this%htot_ml) &
               + arr_bytes(this%mld_filtered) &
               + arr_bytes(this%uhml) &
               + arr_bytes(this%vhml)
   end function ocean_mle_bytes

end module rdb_ocean_mle
