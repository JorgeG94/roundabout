!! First-baroclinic gravity-wave speed + Rossby deformation radius.
module rdb_ocean_wave_speed
   !! Per-column first-baroclinic internal gravity-wave speed `cg1`
   !! (m/s) and first-mode Rossby radius `Rd` (m), plus `Rd/dx` (the
   !! GM/Redi/MEKE resolution ratio).  Solves the rigid-lid
   !! Sturm-Liouville eigenproblem discretised from layer thicknesses
   !! and per-interface reduced gravities `gprime = (g/rho0)*max(0,drho)`
   !! (rho_layer authoritative); largest c^2 via a fixed-budget
   !! Sturm-count bisection (no early exit -> warp-divergence-free).
   !! Reference: Chelton et al. (1998).
   !!
   !! Vertical ordering (load-bearing): rdb is bottom-up (k=1 bed,
   !! k=nz surface); the eigensolve is surface-down so the column kernel
   !! FLIPS on gather (local k_loc=1 is the surface layer).
   !!
   !! Diagnostic, default off (`&ocean_wavespeed_nml enable=.false.`) ⇒
   !! kernel never called, bit-identical.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY
#else
   use rdb_constants, only: wp, GRAVITY, NZ_STACK_MAX
#endif
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_wave_speed_t
   public :: wavespeed_compute
   public :: wavespeed_cg1_column
   public :: wavespeed_rd

   ! ---- Fixed numerical constants (NOT namelist knobs) ----
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran workaround: kept module-local rather than imported from
      !! rdb_constants.  LFortran 0.64 lowers an *imported* parameter used to
      !! dimension an explicit-shape dummy at a call site inside a PURE
      !! procedure to an impure runtime getter (`__lcompilers_get_*`) and then
      !! rejects it — a module-local parameter of the same value folds cleanly.
      !! Must track rdb_constants%NZ_STACK_MAX (128).
#endif
   integer, parameter :: MAX_ITT = 40
      !! Fixed bisection budget — no early exit (GPU-divergence-free).
   integer, parameter :: MAX_DBL = 128
      !! Cap on the Sturm-count doubling prelude (data-dependent trip
      !! count, but O(1) det evals; the inner bisection stays fixed).
   real(wp), parameter :: C2_SCALE = 1.0_wp/(4096.0_wp*4096.0_wp)
      !! Per-row determinant rescale `s` (slows det growth between rows).
   real(wp), parameter :: RESCALE = 1024.0_wp**4
      !! Dynamic-rescale ceiling to keep `det` representable for large kc.
   real(wp), parameter :: I_RESCALE = 1.0_wp/(1024.0_wp**4)
   real(wp), parameter :: MIN_SPEED2 = 1.0e-8_wp
      !! `(1e-4 m/s)^2` floor: `speed2_tot <= MIN_SPEED2` => cg1 = 0
      !! (no resolvable first-baroclinic mode).
   real(wp), parameter :: LAM_SEED = 1.0e-20_wp
      !! Tiny lower seed for the doubling prelude (below any physical
      !! eigenvalue).
   real(wp), parameter :: TOL_MERGE = 0.001_wp
      !! Relative backtracking-merge threshold (MOM6 wave_speed_tol
      !! default).  The criterion is relative to the column's own
      !! stratification scale `drxh_sum`; no fixed absolute floor.
   real(wp), parameter :: F_DENOM_FLOOR = 1.0e-10_wp
      !! Floor on the Rd denominator (guards f = beta = 0).
   real(wp), parameter :: RD_GUARD = 1.0e-20_wp
      !! Inside-sqrt guard for the smooth equatorial Rd blend.

   type :: ocean_wave_speed_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.

      ! ---- Scheme selection + master switch ----
      logical :: enable = .false.
         !! Master switch.  Default off — bit-identity preserved.
      real(wp) :: mono_n2_depth = -1.0_wp
         !! DEFERRED: limit N^2 from increasing with depth below this
         !! depth.  `< 0` = off.
      logical :: use_ebt = .false.
         !! DEFERRED: pressure-Neumann / equivalent-barotropic variant.
      integer :: n_wavespeed = 1
         !! Cadence: recompute every `n_wavespeed` steps (slow
         !! diagnostic).  Default every step (cheap when default-off).

      ! ---- Environment (copied from the EOS slot at configure) ----
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3) for `gprime`.

      ! ---- Persistent fields ----
      real(wp), allocatable :: f_centre(:, :)
         !! |f| at cell centres (1/s); filled by `build_static` from
         !! `metrics_fill_coriolis` (planetary or beta-plane, per
         !! `&ocean_grid_nml coriolis_scheme`).
      real(wp), allocatable :: beta_centre(:, :)
         !! |grad f| at cell centres (1/(m*s)); static, filled at
         !! configure by `build_static` from the same Coriolis field
         !! (the `meke_length_scales` centred-difference idiom).  Feeds
         !! the equatorial branch of `wavespeed_rd` — replaces the old
         !! namelist-scalar `beta`.
      real(wp), allocatable :: cg1(:, :)
         !! First-baroclinic gravity-wave speed (m/s).
      real(wp), allocatable :: rd(:, :)
         !! First-mode Rossby deformation radius (m).
      real(wp), allocatable :: rd_over_dx(:, :)
         !! Rd / dx (nondim) — the B2 resolution ratio.  `dx` here is
         !! `metrics%dxT` (metres), NOT `grid%dx` (degrees on
         !! spherical/supergrid/tripolar).
   contains
      procedure :: init => ocean_wave_speed_init
      procedure :: destroy => ocean_wave_speed_destroy
      procedure :: enter_data => ocean_wave_speed_enter_data
      procedure :: exit_data => ocean_wave_speed_exit_data
      procedure :: build_static => ocean_wave_speed_build_static
      procedure, non_overridable :: bytes => ocean_wave_speed_bytes
   end type ocean_wave_speed_t

contains

   ! =================================================================
   ! Lifecycle
   ! =================================================================

   subroutine ocean_wave_speed_init(this, grid)
      !! Allocate the persistent (nx, ny) fields.  Always allocates
      !! (configure runs after init, so `enable` is not known yet);
      !! the off-state footprint is four 2D arrays.
      class(ocean_wave_speed_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      allocate (this%f_centre(nx, ny), source=0.0_wp)
      allocate (this%beta_centre(nx, ny), source=0.0_wp)
      allocate (this%cg1(nx, ny), source=0.0_wp)
      allocate (this%rd(nx, ny), source=0.0_wp)
      allocate (this%rd_over_dx(nx, ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_wave_speed_init

   subroutine ocean_wave_speed_destroy(this)
      class(ocean_wave_speed_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%f_centre)) deallocate (this%f_centre)
      if (allocated(this%beta_centre)) deallocate (this%beta_centre)
      if (allocated(this%cg1)) deallocate (this%cg1)
      if (allocated(this%rd)) deallocate (this%rd)
      if (allocated(this%rd_over_dx)) deallocate (this%rd_over_dx)
   end subroutine ocean_wave_speed_destroy

   subroutine ocean_wave_speed_enter_data(this)
      class(ocean_wave_speed_t), intent(inout) :: this
      select type (this)
      type is (ocean_wave_speed_t)
         call ocean_wave_speed_enter_data_impl(this)
      end select
   end subroutine ocean_wave_speed_enter_data

   subroutine ocean_wave_speed_enter_data_impl(this)
      type(ocean_wave_speed_t), intent(inout) :: this
      if (allocated(this%f_centre)) then
         !$acc enter data copyin(this%f_centre)
      end if
      if (allocated(this%beta_centre)) then
         !$acc enter data copyin(this%beta_centre)
      end if
      if (allocated(this%cg1)) then
         !$acc enter data copyin(this%cg1)
      end if
      if (allocated(this%rd)) then
         !$acc enter data copyin(this%rd)
      end if
      if (allocated(this%rd_over_dx)) then
         !$acc enter data copyin(this%rd_over_dx)
      end if
   end subroutine ocean_wave_speed_enter_data_impl

   subroutine ocean_wave_speed_exit_data(this)
      class(ocean_wave_speed_t), intent(inout) :: this
      select type (this)
      type is (ocean_wave_speed_t)
         call ocean_wave_speed_exit_data_impl(this)
      end select
   end subroutine ocean_wave_speed_exit_data

   subroutine ocean_wave_speed_exit_data_impl(this)
      type(ocean_wave_speed_t), intent(inout) :: this
      if (allocated(this%rd_over_dx)) then
         !$acc exit data delete(this%rd_over_dx)
      end if
      if (allocated(this%rd)) then
         !$acc exit data delete(this%rd)
      end if
      if (allocated(this%cg1)) then
         !$acc exit data delete(this%cg1)
      end if
      if (allocated(this%beta_centre)) then
         !$acc exit data delete(this%beta_centre)
      end if
      if (allocated(this%f_centre)) then
         !$acc exit data delete(this%f_centre)
      end if
   end subroutine ocean_wave_speed_exit_data_impl

   subroutine ocean_wave_speed_build_static(this, grid, metrics, f_centre)
      !! Copy a pre-filled cell-centre Coriolis magnitude |f| (1/s) onto
      !! the slot (mirror of `ocean_meke_set_f_centre` — the caller
      !! builds `f_centre` via `fill_coriolis_centre` /
      !! `metrics_fill_coriolis`, which handles beta-plane AND
      !! planetary/spherical), and fill the static `beta_centre =
      !! |grad f|` field with the SAME centred-difference stencil
      !! `meke_length_scales` uses (edge rows/columns left at 0 -> the
      !! extratropical `Rd = cg1/|f|` branch there).
      !!
      !! Divergence note (documented, not fixed here): this takes the
      !! gradient of `f_centre`, which is `|f|` (matching the VarMix /
      !! MEKE house idiom) rather than the signed Coriolis field, so
      !! `|grad|f||` has a kink at the equator where `|grad f|` would be
      !! smooth.  Shared pre-existing divergence with VarMix/MEKE — not
      !! fixed by this PR (see docs/CLOSURE_MATRIX.md).
      !!
      !! Host loops only — these are static (functions of geometry +
      !! planetary f) and never change with time.  Call after `init`,
      !! AFTER the metrics are filled + finalized, BEFORE `enter_data`.
      class(ocean_wave_speed_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_centre(grid%nx_total, grid%ny_total)
      integer :: i, j, nx, ny
      real(wp) :: dfdx, dfdy

      if (.not. allocated(this%f_centre)) return
      nx = grid%nx_total
      ny = grid%ny_total
      do j = 1, ny
         do i = 1, nx
            this%f_centre(i, j) = f_centre(i, j)
         end do
      end do

      do j = 1, ny
         do i = 1, nx
            dfdx = 0.0_wp
            dfdy = 0.0_wp
            if (i > 1 .and. i < nx) then
               dfdx = 0.5_wp*(f_centre(i + 1, j) - f_centre(i - 1, j))*metrics%idxT(i, j)
            end if
            if (j > 1 .and. j < ny) then
               dfdy = 0.5_wp*(f_centre(i, j + 1) - f_centre(i, j - 1))*metrics%idyT(i, j)
            end if
            this%beta_centre(i, j) = sqrt(dfdx*dfdx + dfdy*dfdy)
         end do
      end do
   end subroutine ocean_wave_speed_build_static

   ! =================================================================
   ! Pure device-callable column solver (unit-tested directly)
   ! =================================================================

   pure subroutine wavespeed_cg1_column(nz, h_rak, rho_rak, rho0, cg1)
      !! First-baroclinic wave speed for ONE column.  `h_rak`/`rho_rak`
      !! are in Roundabout ordering (k=1 bed, k=nz surface), fixed-size
      !! NZ_STACK_MAX arrays; only 1..nz are read.  Returns `cg1` (m/s),
      !! 0 for land / homogeneous / `kc<2` / sub-floor columns.
      !! Gathers+flips surface-down, backtracking convective merge,
      !! symmetric tridiag, fixed-budget Sturm-count bisection.
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_rak(NZ_STACK_MAX), rho_rak(NZ_STACK_MAX)
      real(wp), intent(in) :: rho0
      real(wp), intent(out) :: cg1

      ! column arrays (surface-down): 5 total (spec §7)
      real(wp) :: hc(NZ_STACK_MAX), rc(NZ_STACK_MAX)
      real(wp) :: gprime(NZ_STACK_MAX + 1)
      real(wp) :: igu(NZ_STACK_MAX + 1), igl(NZ_STACK_MAX + 1)
      ! scalar carries (names avoid intrinsic shadows: no count/sum/sign/mod/...)
      integer :: k, kg, kc, ki, it, n_chg, slo, smid
      real(wp) :: g_rho0, drxh, hnew, hnew2, rnew2, gp
      real(wp) :: s2tot, lam_lo, lam_hi, lam_mid, lam_probe
      logical :: do_weld, do_weld2

      cg1 = 0.0_wp
      if (nz < 2) return
      g_rho0 = GRAVITY/rho0

      ! ---- gather (FLIP): rdb k=1(bed)..nz(surf) -> surface-down ----
      do k = 1, nz
         kg = nz + 1 - k          ! k=1 -> surface layer
         hc(k) = h_rak(kg)
         rc(k) = rho_rak(kg)
      end do

      ! ---- drxh_sum on the RAW (unmerged) surface-down profile ----
      drxh = 0.0_wp
      do k = 2, nz
         drxh = drxh + 0.5_wp*(hc(k - 1) + hc(k))*max(0.0_wp, rc(k) - rc(k - 1))
      end do

      ! ---- backtracking convective merge (last-committed compare) ----
      ! kc counts committed layers; the committed layer kc lives at
      ! hc(kc)/rc(kc).  Walk the remaining raw layers (still in hc/rc at
      ! index k) and weld or commit.  We compact in place: the committed
      ! stack occupies hc(1..kc); the next raw layer is hc(k) (k>kc).
      kc = 1
      do k = 2, nz
         do_weld = ((rc(k) - rc(kc))*(hc(kc) + hc(k)) < 2.0_wp*TOL_MERGE*drxh)
         if (do_weld) then
            hnew = hc(kc) + hc(k)
            rc(kc) = (hc(kc)*rc(kc) + hc(k)*rc(k))/hnew
            hc(kc) = hnew
            ! backtrack: undo any inversion the weld created (looser tol)
            do
               if (kc < 2) exit
               do_weld2 = ((rc(kc) - rc(kc - 1))*(hc(kc) + hc(kc - 1)) &
                           < TOL_MERGE*drxh)
               if (.not. do_weld2) exit
               hnew2 = hc(kc) + hc(kc - 1)
               rnew2 = (hc(kc)*rc(kc) + hc(kc - 1)*rc(kc - 1))/hnew2
               rc(kc - 1) = rnew2
               hc(kc - 1) = hnew2
               kc = kc - 1
            end do
         else
            kc = kc + 1
            hc(kc) = hc(k)
            rc(kc) = rc(k)
         end if
      end do

      if (kc < 2) return

      ! ---- gprime, Igu, Igl at interfaces K=2..kc ----
      s2tot = 0.0_wp
      do ki = 2, kc
         gp = g_rho0*(rc(ki) - rc(ki - 1))
         gprime(ki) = gp
         if (gp > 0.0_wp) then
            igu(ki) = 1.0_wp/(gp*hc(ki - 1))
            igl(ki) = 1.0_wp/(gp*hc(ki))
            s2tot = s2tot + gp*(hc(ki - 1) + hc(ki))
         else
            igu(ki) = 0.0_wp
            igl(ki) = 0.0_wp
         end if
      end do

      if (s2tot <= MIN_SPEED2) return

      ! ---- bracket via Sturm-count doubling (robust for all nz) ----
      lam_probe = LAM_SEED
      do it = 1, MAX_DBL
         if (lam_probe >= 1.0_wp/MIN_SPEED2) then
            lam_probe = 1.0_wp/MIN_SPEED2
            exit
         end if
         n_chg = sturm_count(igu, igl, kc, lam_probe)
         if (n_chg >= 1) exit
         lam_probe = lam_probe*2.0_wp
      end do
      ! safety: no mode-1 found below the cap
      if (sturm_count(igu, igl, kc, lam_probe) < 1) return

      lam_lo = lam_probe*0.5_wp
      lam_hi = lam_probe

      ! ---- fixed-budget bisection on sign(det) (NO early exit) ----
      slo = det_sign(igu, igl, kc, lam_lo)
      do it = 1, MAX_ITT
         lam_mid = 0.5_wp*(lam_lo + lam_hi)
         smid = det_sign(igu, igl, kc, lam_mid)
         if (smid == slo) then
            lam_lo = lam_mid
         else
            lam_hi = lam_mid
         end if
      end do
      lam_mid = 0.5_wp*(lam_lo + lam_hi)
      if (lam_mid > 0.0_wp) cg1 = 1.0_wp/sqrt(lam_mid)
   end subroutine wavespeed_cg1_column

   pure integer function sturm_count(igu, igl, kc, lam) result(n_chg)
      !! Number of Sturm-sequence sign changes (eigenvalues < lam) via
      !! the three-term determinant recursion with dynamic rescaling.
      !$acc routine seq
      real(wp), intent(in) :: igu(NZ_STACK_MAX + 1), igl(NZ_STACK_MAX + 1)
      integer, intent(in) :: kc
      real(wp), intent(in) :: lam
      integer :: k, prev_s, cur_s
      real(wp) :: det, detm1, detm2, dval, offl, absd, sfac

      n_chg = 0
      detm1 = 1.0_wp
      det = (igu(2) + igl(2)) - lam
      prev_s = 1
      cur_s = merge(1, -1, det >= 0.0_wp)
      if (cur_s /= prev_s) n_chg = n_chg + 1
      prev_s = cur_s
      do k = 3, kc
         absd = abs(det)
         sfac = 1.0_wp
         if (absd > RESCALE) then
            sfac = I_RESCALE
         else if (absd < I_RESCALE .and. absd > 0.0_wp) then
            sfac = RESCALE
         end if
         detm2 = C2_SCALE*detm1*sfac
         detm1 = C2_SCALE*det*sfac
         dval = (igu(k) + igl(k)) - lam
         offl = igu(k)*igl(k - 1)
         det = dval*detm1 - offl*detm2
         cur_s = merge(1, -1, det >= 0.0_wp)
         if (cur_s /= prev_s) n_chg = n_chg + 1
         prev_s = cur_s
      end do
   end function sturm_count

   pure integer function det_sign(igu, igl, kc, lam) result(sgn)
      !! Sign of det(M(lam)) on rows 2..kc (Sturm/Hallberg recursion).
      !$acc routine seq
      real(wp), intent(in) :: igu(NZ_STACK_MAX + 1), igl(NZ_STACK_MAX + 1)
      integer, intent(in) :: kc
      real(wp), intent(in) :: lam
      integer :: k
      real(wp) :: det, detm1, detm2, dval, offl, absd, sfac

      detm1 = 1.0_wp
      det = (igu(2) + igl(2)) - lam
      do k = 3, kc
         absd = abs(det)
         sfac = 1.0_wp
         if (absd > RESCALE) then
            sfac = I_RESCALE
         else if (absd < I_RESCALE .and. absd > 0.0_wp) then
            sfac = RESCALE
         end if
         detm2 = C2_SCALE*detm1*sfac
         detm1 = C2_SCALE*det*sfac
         dval = (igu(k) + igl(k)) - lam
         offl = igu(k)*igl(k - 1)
         det = dval*detm1 - offl*detm2
      end do
      sgn = merge(1, -1, det >= 0.0_wp)
   end function det_sign

   pure function wavespeed_rd(cg1, fabs, beta) result(rd)
      !! Smooth equatorial Rd blend: Rd = cg1/sqrt(f^2 + 2*beta*cg1).
      !! Reduces to cg1/|f| away from the equator and sqrt(cg1/(2*beta))
      !! at f=0; a small inside-sqrt guard + denominator floor handle
      !! f = beta = 0.
      !$acc routine seq
      real(wp), intent(in) :: cg1, fabs, beta
      real(wp) :: rd
      real(wp) :: denom
      denom = sqrt(fabs*fabs + 2.0_wp*beta*cg1 + RD_GUARD)
      rd = cg1/max(denom, F_DENOM_FLOOR)
   end function wavespeed_rd

   ! =================================================================
   ! Main compute (outer shim + column kernel)
   ! =================================================================

   pure subroutine wavespeed_compute(grid, metrics, this, ms)
      !! Fill `this%cg1`, `this%rd`, `this%rd_over_dx` over the domain.
      !! `rho_layer` is a top-level allocatable that reaches the device
      !! directly, so no outer-shim tracer dereference is needed (unlike
      !! EPBL/kappa-shear).  Call at the `n_wavespeed` cadence.
      !! Host guards + dereference here; the explicit-shape `do
      !! concurrent` kernel lives in `wavespeed_compute_impl`
      !! (outer-shim + flat-impl pattern — mirrors `varmix_compute`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_wave_speed_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      integer :: nx, ny, nz

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. allocated(ms%rho_layer)) return
      if (.not. allocated(this%cg1)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call wavespeed_compute_impl(nx, ny, nz, this%rho0, &
                                  ms%h_layer, ms%rho_layer, ms%wet_mask, &
                                  this%f_centre, this%beta_centre, metrics%dxT, &
                                  this%cg1, this%rd, this%rd_over_dx)
   end subroutine wavespeed_compute

   pure subroutine wavespeed_compute_impl(nx, ny, nz, rho0, h_layer, rho_layer, &
                                          wet_mask, f_centre, beta_centre, dxT, &
                                          cg1, rd, rd_over_dx)
      !! Flat-impl wavespeed kernel (explicit-shape; NVHPC
      !! descriptor-walk-free).  Per-column Sturm-Liouville solve
      !! (`wavespeed_cg1_column`) + the deformation-radius blend
      !! (`wavespeed_rd`), then the metres-denominated resolution ratio
      !! `rd_over_dx = rd / dxT`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: rho0
      real(wp), intent(in) :: h_layer(nx, ny, nz), rho_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: f_centre(nx, ny), beta_centre(nx, ny), dxT(nx, ny)
      real(wp), intent(out) :: cg1(nx, ny), rd(nx, ny), rd_over_dx(nx, ny)

      integer :: i, j, k
      real(wp) :: h_col(NZ_STACK_MAX), rho_col(NZ_STACK_MAX)
      real(wp) :: cg1_v, rd_v

      do concurrent(j=1:ny, i=1:nx) &
         local(k, h_col, rho_col, cg1_v, rd_v)
         cg1_v = 0.0_wp
         rd_v = 0.0_wp
         if (wet_mask(i, j) > 0.0_wp .and. nz >= 2) then
            do k = 1, nz
               h_col(k) = h_layer(i, j, k)
               rho_col(k) = rho_layer(i, j, k)
            end do
            call wavespeed_cg1_column(nz, h_col, rho_col, rho0, cg1_v)
            rd_v = wavespeed_rd(cg1_v, f_centre(i, j), beta_centre(i, j))
         end if
         cg1(i, j) = cg1_v
         rd(i, j) = rd_v
         rd_over_dx(i, j) = rd_v/max(dxT(i, j), F_DENOM_FLOOR)
      end do
   end subroutine wavespeed_compute_impl

   pure function ocean_wave_speed_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the wave speed slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_wave_speed_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_centre) &
               + arr_bytes(this%beta_centre) &
               + arr_bytes(this%cg1) &
               + arr_bytes(this%rd) &
               + arr_bytes(this%rd_over_dx)
   end function ocean_wave_speed_bytes

end module rdb_ocean_wave_speed
