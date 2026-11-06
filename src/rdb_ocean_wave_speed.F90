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
   !! Optional equivalent-barotropic (EBT) mode (`use_ebt`): the gravest
   !! PRESSURE mode of the same layered Sturm-Liouville problem
   !! (Killworth, Chelton & de Szoeke 1997; de Szoeke & Chelton 1999;
   !! Killworth & Blundell 2003), used as the vertical shape of mesoscale
   !! eddy stirring (Hallberg 2013).  Unknowns live at the `kc` LAYER
   !! centres of the committed (merged) stack; surface Neumann, bed
   !! Dirichlet (p'=0 at the bed interface, half-cell ghost).  Gives the
   !! speed `cg1_ebt = 1/sqrt(lambda_min)` (-> 2NH/pi for constant N, twice
   !! the first-baroclinic NH/pi) and the one-signed structure
   !! `ebt_struct(i,j,k)`, normalised to 1 at the surface, remapped
   !! conservatively (PLM) from the merged stack onto the raw layers.
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
   use rdb_remap_column, only: remap_column_plm
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_wave_speed_t
   public :: wavespeed_compute
   public :: wavespeed_cg1_column
   public :: wavespeed_mode_column
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
   integer, parameter :: EBT_ITT = 64
      !! Fixed bisection budget for the EBT eigenvalue (bracket
      !! [0, Gershgorin] halves each pass, no early exit).
   integer, parameter :: EBT_NINV = 4
      !! Fixed inverse-iteration count for the EBT eigenvector.
   real(wp), parameter :: EBT_SHIFT = 1.0e-6_wp
      !! Inverse-iteration shift sits at `lambda*(1 - EBT_SHIFT)`.
   real(wp), parameter :: EBT_PIVMIN = 1.0e-100_wp
      !! Pivot floor in the EBT Sturm count.
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
         !! Limit N^2 from increasing with depth below this depth.
         !! `< 0` = off.  NOT IMPLEMENTED (inert); `validate_config`
         !! refuses it together with `use_ebt`.
      logical :: use_ebt = .false.
         !! Also compute the equivalent-barotropic speed `cg1_ebt` and
         !! vertical structure `ebt_struct` (computed in the same column
         !! launch as `cg1`).  `cg1` itself stays first-baroclinic.  GM does
         !! not consume the structure yet.  Default off: `cg1` bit-identical
         !! and `ebt_struct` stays 1.
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
      real(wp), allocatable :: cg1_ebt(:, :)
         !! Equivalent-barotropic gravity-wave speed (m/s); 0 on land and
         !! degenerate columns.  Only written when `use_ebt`.
      real(wp), allocatable :: ebt_struct(:, :, :)
         !! EBT vertical structure at T cells, bottom-up (k=nz surface),
         !! surface value 1.  Initialised to 1 ("no structure"); written
         !! only when `use_ebt` (0 on land / degenerate columns).
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

   subroutine ocean_wave_speed_init(this, grid, nz_ml)
      !! Allocate the persistent fields.  Always allocates (configure runs
      !! after init, so `enable` is not known yet); the off-state footprint
      !! is five 2D arrays plus the `(nx, ny, nz)` EBT structure (`nz_ml`
      !! defaults to 1 when absent: a column of ones).
      class(ocean_wave_speed_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nzs

      nx = grid%nx_total
      ny = grid%ny_total
      nzs = 1
      if (present(nz_ml)) nzs = max(1, nz_ml)
      allocate (this%f_centre(nx, ny), source=0.0_wp)
      allocate (this%beta_centre(nx, ny), source=0.0_wp)
      allocate (this%cg1(nx, ny), source=0.0_wp)
      allocate (this%cg1_ebt(nx, ny), source=0.0_wp)
      allocate (this%ebt_struct(nx, ny, nzs), source=1.0_wp)
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
      if (allocated(this%cg1_ebt)) deallocate (this%cg1_ebt)
      if (allocated(this%ebt_struct)) deallocate (this%ebt_struct)
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
      if (allocated(this%cg1_ebt)) then
         !$acc enter data copyin(this%cg1_ebt)
      end if
      if (allocated(this%ebt_struct)) then
         !$acc enter data copyin(this%ebt_struct)
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
      if (allocated(this%ebt_struct)) then
         !$acc exit data delete(this%ebt_struct)
      end if
      if (allocated(this%cg1_ebt)) then
         !$acc exit data delete(this%cg1_ebt)
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
      !! First-baroclinic wave speed for ONE column (no EBT output).
      !! Thin wrapper over `wavespeed_mode_column` with `use_ebt=.false.`.
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_rak(NZ_STACK_MAX), rho_rak(NZ_STACK_MAX)
      real(wp), intent(in) :: rho0
      real(wp), intent(out) :: cg1
      real(wp) :: cg1_ebt
      real(wp) :: phi(NZ_STACK_MAX)

      call wavespeed_mode_column(nz, h_rak, rho_rak, rho0, .false., cg1, cg1_ebt, phi)
   end subroutine wavespeed_cg1_column

   pure subroutine wavespeed_mode_column(nz, h_rak, rho_rak, rho0, use_ebt, &
                                         cg1, cg1_ebt, phi_rak)
      !! First-baroclinic wave speed `cg1` and (when `use_ebt`) the
      !! equivalent-barotropic speed `cg1_ebt` + structure `phi_rak` for ONE
      !! column.  `h_rak`/`rho_rak`/`phi_rak` are in Roundabout ordering
      !! (k=1 bed, k=nz surface), fixed-size NZ_STACK_MAX arrays; only 1..nz
      !! are read/written.  `cg1` is 0 for land / homogeneous / `kc<2` /
      !! sub-floor columns; those columns also get `cg1_ebt = 0` and
      !! `phi_rak = 0` (structure 0, not 1).  With `use_ebt = .false.`,
      !! `cg1_ebt`/`phi_rak` are set to 0 and the `cg1` arithmetic is
      !! untouched (bit-identical).
      !! Gathers+flips surface-down, backtracking convective merge,
      !! tridiag, fixed-budget Sturm-count bisection.
      !!
      !! EBT: with the merged stack surface-down (layer 1 = surface), layer
      !! `k` has upper interface `k` and lower interface `k+1`, `Igu(K) =
      !! 1/(g'_K h_{K-1})`, `Igl(K) = 1/(g'_K h_K)`.  Pressure rows,
      !!   lam p_1  = Igu(2) (p_1 - p_2)                      (surface Neumann)
      !!   lam p_k  = Igl(k) (p_k - p_{k-1}) + Igu(k+1) (p_k - p_{k+1})
      !!   lam p_kc = Igl(kc) (p_kc - p_{kc-1}) + 2 Igl(kc) p_kc  (bed Dirichlet)
      !! (the bed ghost is a half cell, hence the factor 2).  `lam =
      !! 1/c_ebt^2` is the smallest eigenvalue (Sturm-count bisection on the
      !! LDL pivots, fixed budget); the eigenvector is `EBT_NINV` inverse
      !! iterations at `lam (1 - EBT_SHIFT)`, normalised to 1 at the surface
      !! and remapped (conservative PLM, no boundary extrapolation) from the
      !! merged thicknesses to the raw ones.
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_rak(NZ_STACK_MAX), rho_rak(NZ_STACK_MAX)
      real(wp), intent(in) :: rho0
      logical, intent(in) :: use_ebt
      real(wp), intent(out) :: cg1, cg1_ebt
      real(wp), intent(out) :: phi_rak(NZ_STACK_MAX)

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
      cg1_ebt = 0.0_wp
      phi_rak = 0.0_wp
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

      if (use_ebt .and. cg1 > 0.0_wp) then
         call ebt_column(nz, kc, hc, h_rak, igu, igl, cg1_ebt, phi_rak)
      end if
   end subroutine wavespeed_mode_column

   pure subroutine ebt_column(nz, kc, hc, h_rak, igu, igl, cg1_ebt, phi_rak)
      !! EBT eigenvalue + structure on the committed stack (see
      !! `wavespeed_mode_column`).  `hc(1:kc)` are the merged thicknesses
      !! (surface-down); `h_rak(1:nz)` the raw bottom-up ones.
      !$acc routine seq
      integer, intent(in) :: nz, kc
      real(wp), intent(in) :: hc(NZ_STACK_MAX)
      real(wp), intent(in) :: h_rak(NZ_STACK_MAX)
      real(wp), intent(in) :: igu(NZ_STACK_MAX + 1), igl(NZ_STACK_MAX + 1)
      real(wp), intent(out) :: cg1_ebt
      real(wp), intent(out) :: phi_rak(NZ_STACK_MAX)

      real(wp) :: dg(NZ_STACK_MAX), bp(NZ_STACK_MAX), x(NZ_STACK_MAX)
      real(wp) :: dzo(NZ_STACK_MAX), dzn(NZ_STACK_MAX)
      real(wp) :: qo(NZ_STACK_MAX), qn(NZ_STACK_MAX)
      real(wp) :: lam_lo, lam_hi, lam_mid, lam, sig, q, xmax, rowsum
      integer :: k, it, nneg, kg

      cg1_ebt = 0.0_wp
      phi_rak = 0.0_wp

      ! ---- diagonal (surface-down layer rows) ----
      dg(1) = igu(2)
      do k = 2, kc - 1
         dg(k) = igl(k) + igu(k + 1)
      end do
      dg(kc) = 3.0_wp*igl(kc)

      ! ---- Gershgorin upper bound ----
      lam_hi = 0.0_wp
      do k = 1, kc
         rowsum = dg(k)
         if (k >= 2) rowsum = rowsum + igl(k)
         if (k < kc) rowsum = rowsum + igu(k + 1)
         lam_hi = max(lam_hi, rowsum)
      end do
      lam_lo = 0.0_wp

      ! ---- fixed-budget bisection on the Sturm count (pivots < 0) ----
      do it = 1, EBT_ITT
         lam_mid = 0.5_wp*(lam_lo + lam_hi)
         q = dg(1) - lam_mid
         nneg = 0
         if (q < 0.0_wp) nneg = 1
         do k = 2, kc
            if (abs(q) < EBT_PIVMIN) q = EBT_PIVMIN
            q = dg(k) - lam_mid - igu(k)*igl(k)/q
            if (q < 0.0_wp) nneg = nneg + 1
         end do
         if (nneg >= 1) then
            lam_hi = lam_mid
         else
            lam_lo = lam_mid
         end if
      end do
      lam = 0.5_wp*(lam_lo + lam_hi)
      if (lam <= 0.0_wp) return
      cg1_ebt = 1.0_wp/sqrt(lam)

      ! ---- inverse iteration (Thomas on the shifted rows) ----
      sig = lam*(1.0_wp - EBT_SHIFT)
      do k = 1, kc
         x(k) = 1.0_wp
      end do
      do it = 1, EBT_NINV
         bp(1) = dg(1) - sig
         do k = 2, kc
            bp(k) = dg(k) - sig - igu(k)*igl(k)/bp(k - 1)
            x(k) = x(k) + igl(k)/bp(k - 1)*x(k - 1)
         end do
         x(kc) = x(kc)/bp(kc)
         do k = kc - 1, 1, -1
            x(k) = (x(k) + igu(k + 1)*x(k + 1))/bp(k)
         end do
         xmax = 0.0_wp
         do k = 1, kc
            xmax = max(xmax, abs(x(k)))
         end do
         do k = 1, kc
            x(k) = x(k)/xmax
         end do
      end do

      ! ---- surface normalisation; zero surface value => structure 0 ----
      if (x(1) == 0.0_wp) return
      xmax = 1.0_wp/x(1)
      do k = 1, kc
         qo(k) = max(0.0_wp, x(k)*xmax)
         dzo(k) = hc(k)
      end do
      ! pad the committed stack to nz layers with zero-thickness copies of
      ! the bed value, then PLM-remap onto the raw (surface-down) thicknesses
      do k = kc + 1, nz
         qo(k) = qo(kc)
         dzo(k) = 0.0_wp
      end do
      do k = 1, nz
         dzn(k) = h_rak(nz + 1 - k)
      end do
      call remap_column_plm(nz, dzo, dzn, qo, qn)

      ! flip back to bottom-up (surface value pinned to exactly 1: the
      ! surface raw layer lies inside committed layer 1, whose PLM slope is
      ! flattened at the boundary); vanished raw layers inherit the committed
      ! value above them (remap returns 0 for dz_new <= 0)
      do k = 1, nz
         kg = nz + 1 - k
         if (dzn(k) > 0.0_wp) then
            phi_rak(kg) = max(0.0_wp, qn(k))
            if (k == 1) phi_rak(kg) = 1.0_wp
         else if (k > 1) then
            phi_rak(kg) = phi_rak(kg + 1)
         else
            phi_rak(kg) = 1.0_wp
         end if
      end do
   end subroutine ebt_column

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
      !! Fill `this%cg1`, `this%rd`, `this%rd_over_dx` over the domain
      !! (and `cg1_ebt` + `ebt_struct` when `use_ebt`).
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
      integer :: nx, ny, nz, nze

      if (.not. this%is_init) return
      if (.not. this%enable) return
      if (.not. allocated(ms%rho_layer)) return
      if (.not. allocated(this%cg1)) return
      if (.not. allocated(this%ebt_struct)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      ! The kernel's `ebt_struct` dummy is declared with the array's REAL
      ! third extent, not `nz`: an `init` without `nz_ml` allocates one
      ! layer, and declaring it `(nx, ny, nz)` would make the device map
      ! a partially-present section (an OpenACC abort) even with
      ! `use_ebt` off, when the kernel never touches it.  With `use_ebt`
      ! on, a mismatched extent is a setup bug, so fail loud rather than
      ! skip the whole solve (which would also leave `cg1` stale).
      nze = size(this%ebt_struct, 3)
      if (this%use_ebt .and. nze /= nz) then
         error stop "wavespeed_compute: use_ebt needs ebt_struct allocated with nz_ml layers "// &
            "(ocean_wave_speed_init nz_ml)"
      end if

      call wavespeed_compute_impl(nx, ny, nz, this%rho0, this%use_ebt, &
                                  ms%h_layer, ms%rho_layer, ms%wet_mask, &
                                  this%f_centre, this%beta_centre, metrics%dxT, &
                                  this%cg1, this%rd, this%rd_over_dx, &
                                  this%cg1_ebt, nze, this%ebt_struct)
   end subroutine wavespeed_compute

   pure subroutine wavespeed_compute_impl(nx, ny, nz, rho0, use_ebt, h_layer, &
                                          rho_layer, wet_mask, f_centre, beta_centre, &
                                          dxT, cg1, rd, rd_over_dx, cg1_ebt, nze, ebt_struct)
      !! Flat-impl wavespeed kernel (explicit-shape; NVHPC
      !! descriptor-walk-free).  Per-column Sturm-Liouville solve
      !! (`wavespeed_cg1_column`) + the deformation-radius blend
      !! (`wavespeed_rd`), then the metres-denominated resolution ratio
      !! `rd_over_dx = rd / dxT`.  With `use_ebt`, the same column launch
      !! also writes `cg1_ebt` and `ebt_struct` (0 on land / degenerate
      !! columns); with it off neither is touched.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: rho0
      logical, intent(in) :: use_ebt
      real(wp), intent(in) :: h_layer(nx, ny, nz), rho_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: f_centre(nx, ny), beta_centre(nx, ny), dxT(nx, ny)
      real(wp), intent(out) :: cg1(nx, ny), rd(nx, ny), rd_over_dx(nx, ny)
      integer, intent(in) :: nze
         !! Third extent of `ebt_struct` (== nz whenever `use_ebt`).
      real(wp), intent(inout) :: cg1_ebt(nx, ny), ebt_struct(nx, ny, nze)

      integer :: i, j, k
      real(wp) :: h_col(NZ_STACK_MAX), rho_col(NZ_STACK_MAX), phi_col(NZ_STACK_MAX)
      real(wp) :: cg1_v, rd_v, ebt_v

      do concurrent(j=1:ny, i=1:nx) &
         local(k, h_col, rho_col, phi_col, cg1_v, rd_v, ebt_v)
         cg1_v = 0.0_wp
         rd_v = 0.0_wp
         ebt_v = 0.0_wp
         phi_col = 0.0_wp
         if (wet_mask(i, j) > 0.0_wp .and. nz >= 2) then
            do k = 1, nz
               h_col(k) = h_layer(i, j, k)
               rho_col(k) = rho_layer(i, j, k)
            end do
            call wavespeed_mode_column(nz, h_col, rho_col, rho0, use_ebt, &
                                       cg1_v, ebt_v, phi_col)
            rd_v = wavespeed_rd(cg1_v, f_centre(i, j), beta_centre(i, j))
         end if
         cg1(i, j) = cg1_v
         rd(i, j) = rd_v
         rd_over_dx(i, j) = rd_v/max(dxT(i, j), F_DENOM_FLOOR)
         if (use_ebt) then
            cg1_ebt(i, j) = ebt_v
            do k = 1, nz
               ebt_struct(i, j, k) = phi_col(k)
            end do
         end if
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
               + arr_bytes(this%cg1_ebt) &
               + arr_bytes(this%ebt_struct) &
               + arr_bytes(this%rd) &
               + arr_bytes(this%rd_over_dx)
   end function ocean_wave_speed_bytes

end module rdb_ocean_wave_speed
