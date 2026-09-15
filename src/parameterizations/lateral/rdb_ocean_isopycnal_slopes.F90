!! Isopycnal (neutral) slope diagnostics for the ocean dynamical core.
module rdb_ocean_isopycnal_slopes
   !! Neutral-density slope `S = -∇ρ/∂_zρ` and interface stratification
   !! `N²` at C-grid layer interfaces (u-faces → `slope_x`/`n2_u`, v-faces
   !! → `slope_y`/`n2_v`).  Purely diagnostic — consumed by GM/Redi/VarMix/
   !! MLE; no flux consumer here.  Harmonic-thickness-weighted FD form
   !! (Griffies 1998).
   !!
   !! Bottom-up convention (k=1 bed, k=nz_ml surface).  Interface `K`
   !! 1..nz_ml+1; K=1 bed and K=nz_ml+1 surface forced to zero slope.
   !! Interior `K` straddles layer `k=K` (above, surface side) and `k=K-1`
   !! (below, bed side).  A vert-fill pre-pass diffuses T/S into massless
   !! layers so `T=S=0` ghosts don't corrupt gradients.
   !!
   !! Default off (`&ocean_slopes_nml enable=.false.`) ⇒ slot allocated but
   !! `ocean_slopes_compute` no-ops ⇒ bit-identical.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, H_DIV_EPS
#endif
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_specvol_derivs, eos_density_point
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

   public :: ocean_slopes_t
   public :: ocean_slopes_compute
   public :: ocean_slopes_vert_fill_ts
   public :: pressure_above_x   !! exposed for the interface-pressure unit test

   type :: ocean_slopes_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Guard on this, never on
         !! `allocated(...)` (host pointer only; misses GPU mapping).
      logical :: enable = .false.
         !! Master switch (`&ocean_slopes_nml enable`).  Default off ⇒
         !! `ocean_slopes_compute` no-ops ⇒ bit-identical.
      real(wp) :: kd_smooth = 1.0e-6_wp
         !! Vertical diffusivity (m²/s) used by `vert_fill_TS` to fill
         !! massless layers.  Multiplied by `dt` for the smoothing
         !! `kappa·dt`.
      real(wp) :: min_dz_for_n2 = 1.0_wp
         !! Minimum layer thickness (m) used to floor `h` in the N²
         !! vertical-difference denominator, so vanished layers don't
         !! produce a spurious N² spike.
      real(wp) :: rho0 = 1035.0_wp
         !! Reference density (kg/m³) for the N² scaling `g/ρ₀`.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
      integer :: nz_ml = 0

      ! ---- Outputs (interface-located) ----
      real(wp), allocatable :: slope_x(:, :, :)
         !! Isopycnal slope at u-faces, shape `(nx+1, ny, nz+1)`.
      real(wp), allocatable :: slope_y(:, :, :)
         !! Isopycnal slope at v-faces, shape `(nx, ny+1, nz+1)`.
      real(wp), allocatable :: n2_u(:, :, :)
         !! Brunt-Väisälä N² at u-faces (s⁻²), shape `(nx+1, ny, nz+1)`.
      real(wp), allocatable :: n2_v(:, :, :)
         !! Brunt-Väisälä N² at v-faces (s⁻²), shape `(nx, ny+1, nz+1)`.

      ! ---- Vert-filled T/S scratch + interface height ----
      real(wp), allocatable :: t_fill(:, :, :)
         !! Massless-layer-filled temperature scratch, `(nx, ny, nz)`.
      real(wp), allocatable :: s_fill(:, :, :)
         !! Massless-layer-filled salinity scratch, `(nx, ny, nz)`.
      real(wp), allocatable :: e_int(:, :, :)
         !! Interface height (m), bottom-up cumulative from bathy,
         !! `(nx, ny, nz+1)`; `e_int(:,:,1)` = bed, `(:,:,nz+1)` = surface.
   contains
      procedure, non_overridable :: init => ocean_slopes_init
      procedure, non_overridable :: destroy => ocean_slopes_destroy
      procedure, non_overridable :: enter_data => ocean_slopes_enter_data
      procedure, non_overridable :: exit_data => ocean_slopes_exit_data
      procedure, non_overridable :: bytes => ocean_slopes_bytes
   end type ocean_slopes_t

contains

   subroutine ocean_slopes_init(this, grid, nz_ml)
      !! Allocate the slope / N² outputs + the vert-fill T/S scratch +
      !! the interface-height buffer.  Default `nz_ml = 1` preserves the
      !! barotropic-only constructor; pass `nz_ml = ms%nz_ml` for the
      !! multilayer driver.  Setup code uses plain host allocation (no
      !! `do concurrent` before `enter_data`).
      class(ocean_slopes_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml
      if (nz < 1) nz = 1
      ! Fail loud at configure: vert_fill uses NZ_STACK_MAX-sized column
      ! locals; nz beyond that silently overruns the GPU stack.
      if (nz > NZ_STACK_MAX) then
         error stop "ocean_slopes_init: nz_ml exceeds NZ_STACK_MAX "// &
            "(raise NZ_STACK_MAX in rdb_constants)"
      end if
      this%nx_total = nx
      this%ny_total = ny
      this%nz_ml = nz

      allocate (this%slope_x(nx + 1, ny, nz + 1), source=0.0_wp)
      allocate (this%slope_y(nx, ny + 1, nz + 1), source=0.0_wp)
      allocate (this%n2_u(nx + 1, ny, nz + 1), source=0.0_wp)
      allocate (this%n2_v(nx, ny + 1, nz + 1), source=0.0_wp)
      allocate (this%t_fill(nx, ny, nz), source=0.0_wp)
      allocate (this%s_fill(nx, ny, nz), source=0.0_wp)
      allocate (this%e_int(nx, ny, nz + 1), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_slopes_init

   subroutine ocean_slopes_destroy(this)
      class(ocean_slopes_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%slope_x)) deallocate (this%slope_x)
      if (allocated(this%slope_y)) deallocate (this%slope_y)
      if (allocated(this%n2_u)) deallocate (this%n2_u)
      if (allocated(this%n2_v)) deallocate (this%n2_v)
      if (allocated(this%t_fill)) deallocate (this%t_fill)
      if (allocated(this%s_fill)) deallocate (this%s_fill)
      if (allocated(this%e_int)) deallocate (this%e_int)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_slopes_destroy

   subroutine ocean_slopes_enter_data(this)
      !! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule:
      !! bare polymorphic `copyin(this)` maps the stack descriptor → AMD
      !! libomptarget cross-slot overlap crash).
      class(ocean_slopes_t), intent(inout) :: this
      select type (this)
      type is (ocean_slopes_t)
         call ocean_slopes_enter_data_impl(this)
      end select
   end subroutine ocean_slopes_enter_data

   subroutine ocean_slopes_enter_data_impl(this)
      type(ocean_slopes_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%slope_x, this%slope_y)
      !$acc enter data copyin(this%n2_u, this%n2_v)
      !$acc enter data copyin(this%t_fill, this%s_fill, this%e_int)
   end subroutine ocean_slopes_enter_data_impl

   subroutine ocean_slopes_exit_data(this)
      class(ocean_slopes_t), intent(inout) :: this
      select type (this)
      type is (ocean_slopes_t)
         call ocean_slopes_exit_data_impl(this)
      end select
   end subroutine ocean_slopes_exit_data

   subroutine ocean_slopes_exit_data_impl(this)
      type(ocean_slopes_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%t_fill, this%s_fill, this%e_int)
      !$acc exit data delete(this%n2_u, this%n2_v)
      !$acc exit data delete(this%slope_x, this%slope_y)
   end subroutine ocean_slopes_exit_data_impl

   subroutine ocean_slopes_compute(grid, metrics, eos, slopes, ms, dt)
      !! Public entry point — fill `slope_x`/`slope_y` + `n2_u`/`n2_v` at
      !! all interfaces.  No-op if absent / uninitialised / disabled, so
      !! the driver can call it unconditionally.  Pipeline: vert-fill T/S
      !! → build interface heights `e_int` → u-face pass → v-face pass.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(eos_t), intent(in) :: eos
      type(ocean_slopes_t), intent(inout), optional :: slopes
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt

      if (.not. present(slopes)) return
      if (.not. slopes%is_init) return
      if (.not. slopes%enable) return
      if (ms%idx_temperature <= 0 .or. ms%idx_salinity <= 0) return
      if (.not. allocated(ms%h_layer)) return

      ! Outer shim: dereference the tracer-registry hTr arrays (array of
      ! derived types ⇒ device indirection) on the host, pass the flat
      ! top-level allocatables into the flat-impl kernel.
      call ocean_slopes_compute_impl(grid, metrics, eos, slopes, ms, &
                                     ms%h_layer, &
                                     ms%tracers(ms%idx_temperature)%hTr, &
                                     ms%tracers(ms%idx_salinity)%hTr, dt)
   end subroutine ocean_slopes_compute

   subroutine ocean_slopes_compute_impl(grid, metrics, eos, slopes, ms, &
                                        h_layer, t_htr, s_htr, dt)
      !! Flat-impl: explicit-shape dummies for the prognostic arrays so
      !! NVHPC doesn't descriptor-walk per launch.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(eos_t), intent(in) :: eos
      type(ocean_slopes_t), intent(inout) :: slopes
      type(multilayer_state_t), intent(in) :: ms
      integer :: nx, ny, nz
      real(wp), intent(in) :: h_layer(slopes%nx_total, slopes%ny_total, slopes%nz_ml)
      real(wp), intent(in) :: t_htr(slopes%nx_total, slopes%ny_total, slopes%nz_ml)
      real(wp), intent(in) :: s_htr(slopes%nx_total, slopes%ny_total, slopes%nz_ml)
      real(wp), intent(in) :: dt

      nx = slopes%nx_total
      ny = slopes%ny_total
      nz = slopes%nz_ml
      if (ms%nz_ml /= nz) return

      ! (1) Vert-fill T/S into the scratch (massless layers diffused).
      call ocean_slopes_vert_fill_ts(nx, ny, nz, h_layer, t_htr, s_htr, &
                                     slopes%kd_smooth, dt, &
                                     slopes%t_fill, slopes%s_fill)

      ! (2) Bottom-up interface heights from Σ h_layer (bed datum 0).
      call ocean_slopes_build_e(nx, ny, nz, h_layer, slopes%e_int)

      ! (3) u-face slopes + N².
      call ocean_slopes_pass_x(nx, ny, nz, eos, slopes%rho0, &
                               slopes%min_dz_for_n2, h_layer, &
                               slopes%t_fill, slopes%s_fill, slopes%e_int, &
                               metrics%idxCu, metrics%wet_u, &
                               slopes%slope_x, slopes%n2_u)

      ! (4) v-face slopes + N².
      call ocean_slopes_pass_y(nx, ny, nz, eos, slopes%rho0, &
                               slopes%min_dz_for_n2, h_layer, &
                               slopes%t_fill, slopes%s_fill, slopes%e_int, &
                               metrics%idyCv, metrics%wet_v, &
                               slopes%slope_y, slopes%n2_v)
   end subroutine ocean_slopes_compute_impl

   pure subroutine ocean_slopes_build_e(nx, ny, nz, h_layer, e_int)
      !! Build interface heights bottom-up: `e_int(:,:,1) = 0` (bed),
      !! `e_int(:,:,K+1) = e_int(:,:,K) + h_layer(:,:,K)`.  A per-column
      !! serial cumulative sum (parallel over i,j); only the across-face
      !! difference is consumed, so the absolute bed datum is irrelevant.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(out) :: e_int(nx, ny, nz + 1)
      integer :: i, j, k
      do concurrent(j=1:ny, i=1:nx)
         e_int(i, j, 1) = 0.0_wp
         do k = 1, nz
            e_int(i, j, k + 1) = e_int(i, j, k) + h_layer(i, j, k)
         end do
      end do
   end subroutine ocean_slopes_build_e

   subroutine ocean_slopes_vert_fill_ts(nx, ny, nz, h_layer, t_htr, s_htr, &
                                        kd_smooth, dt, t_fill, s_fill)
      !! Fill massless layers in T/S with sensible values via one pass of
      !! constant-`kappa·dt` vertical diffusion — a SINGLE forward-elim +
      !! back-sub Thomas sweep per column (no iteration).  Operates on the
      !! tracer-from-hTr conversion (`T = hTr/h`) and writes the scratch
      !! `t_fill`/`s_fill`; the prognostic tracers are untouched.
      !!
      !! `kap_dt_x2 = 2·kappa·dt`; the inter-layer entrainment is
      !! `ent(K) = kap_dt_x2 / ((h(k)+h(k+1)) + h_neglect)`.  Surface +
      !! bed boundary rows close the tridiagonal exactly.  Column locals
      !! are fixed-size (`NZ_STACK_MAX`) so the `local()` clause is legal
      !! on `-stdpar=gpu` (dummy-sized automatics crash).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: t_htr(nx, ny, nz)
      real(wp), intent(in) :: s_htr(nx, ny, nz)
      real(wp), intent(in) :: kd_smooth, dt
      real(wp), intent(out) :: t_fill(nx, ny, nz)
      real(wp), intent(out) :: s_fill(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: kap_dt_x2, h_neglect, h0c
      real(wp) :: ent(NZ_STACK_MAX + 1)
      real(wp) :: c1(NZ_STACK_MAX)
      real(wp) :: b1, d1, h_tr, h_eff
      real(wp) :: t_in(NZ_STACK_MAX), s_in(NZ_STACK_MAX)

      kap_dt_x2 = 2.0_wp*kd_smooth*dt
      h_neglect = H_DIV_EPS
      h0c = h_neglect

      if (kap_dt_x2 <= 0.0_wp .or. nz < 2) then
         ! No smoothing — pass the raw tracer-from-hTr through.
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            h_eff = max(h_layer(i, j, k), H_VANISHED)
            t_fill(i, j, k) = t_htr(i, j, k)/h_eff
            s_fill(i, j, k) = s_htr(i, j, k)/h_eff
         end do
         return
      end if

      ! Per-column Thomas sweep.  k=1 bed ... k=nz surface (bottom-up);
      ! the tridiagonal couples layer k to k±1 identically regardless of
      ! orientation, so the bottom-up index runs the published sweep with
      ! "k=1" as the first boundary row.
      do concurrent(j=1:ny, i=1:nx) &
         local(k, ent, c1, b1, d1, h_tr, h_eff, t_in, s_in)
         ! T,S from hTr/h (floor h consistently — vanished layers feed a
         ! near-zero raw T that the diffusion then overwrites).
         do k = 1, nz
            h_eff = max(h_layer(i, j, k), H_VANISHED)
            t_in(k) = t_htr(i, j, k)/h_eff
            s_in(k) = s_htr(i, j, k)/h_eff
         end do

         ! Forward elimination — first (bed) boundary row at k=1.
         ent(2) = kap_dt_x2/((h_layer(i, j, 1) + h_layer(i, j, 2)) + h0c)
         h_tr = h_layer(i, j, 1) + h_neglect
         b1 = 1.0_wp/(h_tr + ent(2))
         d1 = b1*h_tr
         t_fill(i, j, 1) = (b1*h_tr)*t_in(1)
         s_fill(i, j, 1) = (b1*h_tr)*s_in(1)
         do k = 2, nz - 1
            ent(k + 1) = kap_dt_x2/((h_layer(i, j, k) + h_layer(i, j, k + 1)) + h0c)
            h_tr = h_layer(i, j, k) + h_neglect
            c1(k) = ent(k)*b1
            b1 = 1.0_wp/((h_tr + d1*ent(k)) + ent(k + 1))
            d1 = b1*(h_tr + d1*ent(k))
            t_fill(i, j, k) = b1*(h_tr*t_in(k) + ent(k)*t_fill(i, j, k - 1))
            s_fill(i, j, k) = b1*(h_tr*s_in(k) + ent(k)*s_fill(i, j, k - 1))
         end do
         ! Last (surface) boundary row at k=nz.
         c1(nz) = ent(nz)*b1
         h_tr = h_layer(i, j, nz) + h_neglect
         b1 = 1.0_wp/(h_tr + d1*ent(nz))
         t_fill(i, j, nz) = b1*(h_tr*t_in(nz) + ent(nz)*t_fill(i, j, nz - 1))
         s_fill(i, j, nz) = b1*(h_tr*s_in(nz) + ent(nz)*s_fill(i, j, nz - 1))
         ! Back substitution.
         do k = nz - 1, 1, -1
            t_fill(i, j, k) = t_fill(i, j, k) + c1(k + 1)*t_fill(i, j, k + 1)
            s_fill(i, j, k) = s_fill(i, j, k) + c1(k + 1)*s_fill(i, j, k + 1)
         end do
      end do
   end subroutine ocean_slopes_vert_fill_ts

   pure subroutine ocean_slopes_pass_x(nx, ny, nz, eos, rho0, min_dz, &
                                       h_layer, t_fill, s_fill, e_int, &
                                       idxCu, wet_u, slope_x, n2_u)
      !! u-face slope + N² pass.  Interface `K` (interior 2..nz) straddles
      !! layer `k=K` (above, surface side) and `k=K-1` (below, bed side).
      !! Bed (K=1) + surface (K=nz+1) are forced to zero.  The u-face at
      !! (i,j) sits between cells (i-1,j) and (i,j); pairs columns
      !! `iw=i-1` (west) and `i` (east), so loop `i=2:nx`.
      integer, intent(in) :: nx, ny, nz
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rho0, min_dz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: t_fill(nx, ny, nz)
      real(wp), intent(in) :: s_fill(nx, ny, nz)
      real(wp), intent(in) :: e_int(nx, ny, nz + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny)
      real(wp), intent(in) :: wet_u(nx + 1, ny)
      real(wp), intent(out) :: slope_x(nx + 1, ny, nz + 1)
      real(wp), intent(out) :: n2_u(nx + 1, ny, nz + 1)

      integer :: i, j, k, iw, ka, kb
      real(wp) :: pres_u, t_u, s_u, rho_u, dsv_dt, dsv_ds, drdt, drds
      real(wp) :: drdiA, drdiB, drdkL, drdkR
      real(wp) :: hg2A, hg2B, hg2L, hg2R, haA, haB, haL, haR
      real(wp) :: dzaL, dzaR, wtA, wtB, wtL, wtR
      real(wp) :: drdx, drdz, mag2, slope, presL, presR
      real(wp) :: g_rho0, mask

      g_rho0 = GRAVITY/rho0

      ! Bed + surface interfaces: zero everywhere.
      do concurrent(j=1:ny, i=1:nx + 1)
         slope_x(i, j, 1) = 0.0_wp
         slope_x(i, j, nz + 1) = 0.0_wp
         n2_u(i, j, 1) = 0.0_wp
         n2_u(i, j, nz + 1) = 0.0_wp
      end do
      ! Wall faces (i=1, i=nx+1): zero at all interfaces (no interior pair).
      do concurrent(k=1:nz + 1, j=1:ny)
         slope_x(1, j, k) = 0.0_wp
         slope_x(nx + 1, j, k) = 0.0_wp
         n2_u(1, j, k) = 0.0_wp
         n2_u(nx + 1, j, k) = 0.0_wp
      end do

      ! Interior interfaces K = 2..nz, interior u-faces i = 2..nx.
      do concurrent(k=2:nz, j=1:ny, i=2:nx) &
         local(iw, ka, kb, pres_u, t_u, s_u, rho_u, dsv_dt, dsv_ds, &
               drdt, drds, drdiA, drdiB, drdkL, drdkR, &
               hg2A, hg2B, hg2L, hg2R, haA, haB, haL, haR, &
               dzaL, dzaR, wtA, wtB, wtL, wtR, drdx, drdz, &
               mag2, slope, presL, presR, mask)
         iw = i - 1
         ka = k          ! layer ABOVE the interface (surface side)
         kb = k - 1      ! layer BELOW the interface (bed side)

         ! Interface pressure: accumulate from the surface (k=nz) down to
         ! the layer above this interface.  Surface-relative hydrostatic
         ! pressure at the interface = g·ρ₀·Σ_{above} h.
         presL = pressure_above_x(nx, ny, nz, h_layer, iw, j, ka, rho0)
         presR = pressure_above_x(nx, ny, nz, h_layer, i, j, ka, rho0)
         pres_u = 0.5_wp*(presL + presR)

         ! 4-point interface T/S (two columns × two adjacent layers).
         t_u = 0.25_wp*((t_fill(iw, j, ka) + t_fill(i, j, ka)) + &
                        (t_fill(iw, j, kb) + t_fill(i, j, kb)))
         s_u = 0.25_wp*((s_fill(iw, j, ka) + s_fill(i, j, ka)) + &
                        (s_fill(iw, j, kb) + s_fill(i, j, kb)))

         ! Locally-referenced density derivatives: drho_dX = -ρ²·dSV/dX.
         rho_u = eos_density_point(eos, t_u, s_u, pres_u)
         call eos_specvol_derivs(eos, t_u, s_u, pres_u, dsv_dt, dsv_ds)
         drdt = -(rho_u*rho_u)*dsv_dt
         drds = -(rho_u*rho_u)*dsv_ds

         ! Along-layer horizontal ρ-gradients, above (A=ka) / below (B=kb).
         drdiA = drdt*(t_fill(i, j, ka) - t_fill(iw, j, ka)) + &
                 drds*(s_fill(i, j, ka) - s_fill(iw, j, ka))
         drdiB = drdt*(t_fill(i, j, kb) - t_fill(iw, j, kb)) + &
                 drds*(s_fill(i, j, kb) - s_fill(iw, j, kb))

         ! Vertical ρ-difference (below - above): drho_dX·(X[kb]-X[ka]).
         ! For stable stratification (lighter water above) this gives
         ! drdk>0 ⇒ drdz>0 ⇒ N²>0 (punch-list sign fix #2).
         drdkL = drdt*(t_fill(iw, j, kb) - t_fill(iw, j, ka)) + &
                 drds*(s_fill(iw, j, kb) - s_fill(iw, j, ka))
         drdkR = drdt*(t_fill(i, j, kb) - t_fill(i, j, ka)) + &
                 drds*(s_fill(i, j, kb) - s_fill(i, j, ka))

         ! Harmonic-mean thickness weights.
         hg2A = h_layer(iw, j, ka)*h_layer(i, j, ka) + H_DIV_EPS*H_DIV_EPS
         hg2B = h_layer(iw, j, kb)*h_layer(i, j, kb) + H_DIV_EPS*H_DIV_EPS
         hg2L = h_layer(iw, j, ka)*h_layer(iw, j, kb) + H_DIV_EPS*H_DIV_EPS
         hg2R = h_layer(i, j, ka)*h_layer(i, j, kb) + H_DIV_EPS*H_DIV_EPS
         haA = 0.5_wp*(h_layer(iw, j, ka) + h_layer(i, j, ka)) + H_DIV_EPS
         haB = 0.5_wp*(h_layer(iw, j, kb) + h_layer(i, j, kb)) + H_DIV_EPS
         haL = 0.5_wp*(h_layer(iw, j, ka) + h_layer(iw, j, kb)) + H_DIV_EPS
         haR = 0.5_wp*(h_layer(i, j, ka) + h_layer(i, j, kb)) + H_DIV_EPS
         ! Vertical centre spacing across the interface (floored).
         dzaL = max(haL, min_dz)
         dzaR = max(haR, min_dz)
         wtA = hg2A*haB
         wtB = hg2B*haA
         wtL = hg2L*(haR*dzaR)
         wtR = hg2R*(haL*dzaL)

         drdz = ((wtL*drdkL) + (wtR*drdkR))/((dzaL*wtL) + (dzaR*wtR))

         ! Interface-tilt rotation term + metric scaling.
         drdx = ((wtA*drdiA + wtB*drdiB)/(wtA + wtB) - &
                 drdz*(e_int(iw, j, k) - e_int(i, j, k)))*idxCu(i, j)

         mag2 = drdx*drdx + drdz*drdz
         if (mag2 > 0.0_wp) then
            slope = drdx/sqrt(mag2)
         else
            slope = 0.0_wp
         end if

         mask = wet_u(i, j)
         slope_x(i, j, k) = slope*mask
         n2_u(i, j, k) = g_rho0*drdz*mask
      end do
   end subroutine ocean_slopes_pass_x

   pure function pressure_above_x(nx, ny, nz, h_layer, ic, jc, ka, rho0) result(p)
      !! Surface-relative hydrostatic pressure at the interface K straddled
      !! by layer `ka` (above, surface-side) and `ka-1` (below): the
      !! interface sits at the BOTTOM of layer `ka`, so the water column
      !! above it is layers `ka..nz` (bottom-up, k=nz the surface).
      !! p = g·ρ₀·Σ_{k'=ka}^{nz} h(k') — the sum INCLUDES `ka` (the layer
      !! directly above the interface); omitting it shorts the pressure by
      !! one layer (~5e5 Pa) and biases pressure-dependent EOS derivatives.
      !$acc routine seq
      integer, intent(in) :: nx, ny, nz, ic, jc, ka
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: rho0
      real(wp) :: p
      integer :: kk
      p = 0.0_wp
      do kk = nz, ka, -1
         p = p + GRAVITY*rho0*h_layer(ic, jc, kk)
      end do
   end function pressure_above_x

   pure subroutine ocean_slopes_pass_y(nx, ny, nz, eos, rho0, min_dz, &
                                       h_layer, t_fill, s_fill, e_int, &
                                       idyCv, wet_v, slope_y, n2_v)
      !! v-face slope + N² pass — mirror of `pass_x` with v-staggering.
      !! The v-face at (i,j) sits between cells (i,j-1) and (i,j); pairs
      !! columns `js=j-1` (south) and `j` (north), loop `j=2:ny`.
      integer, intent(in) :: nx, ny, nz
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rho0, min_dz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: t_fill(nx, ny, nz)
      real(wp), intent(in) :: s_fill(nx, ny, nz)
      real(wp), intent(in) :: e_int(nx, ny, nz + 1)
      real(wp), intent(in) :: idyCv(nx, ny + 1)
      real(wp), intent(in) :: wet_v(nx, ny + 1)
      real(wp), intent(out) :: slope_y(nx, ny + 1, nz + 1)
      real(wp), intent(out) :: n2_v(nx, ny + 1, nz + 1)

      integer :: i, j, k, js, ka, kb
      real(wp) :: pres_v, t_v, s_v, rho_v, dsv_dt, dsv_ds, drdt, drds
      real(wp) :: drdjA, drdjB, drdkL, drdkR
      real(wp) :: hg2A, hg2B, hg2L, hg2R, haA, haB, haL, haR
      real(wp) :: dzaL, dzaR, wtA, wtB, wtL, wtR
      real(wp) :: drdy, drdz, mag2, slope, presS, presN
      real(wp) :: g_rho0, mask

      g_rho0 = GRAVITY/rho0

      do concurrent(i=1:nx, j=1:ny + 1)
         slope_y(i, j, 1) = 0.0_wp
         slope_y(i, j, nz + 1) = 0.0_wp
         n2_v(i, j, 1) = 0.0_wp
         n2_v(i, j, nz + 1) = 0.0_wp
      end do
      do concurrent(k=1:nz + 1, i=1:nx)
         slope_y(i, 1, k) = 0.0_wp
         slope_y(i, ny + 1, k) = 0.0_wp
         n2_v(i, 1, k) = 0.0_wp
         n2_v(i, ny + 1, k) = 0.0_wp
      end do

      do concurrent(k=2:nz, j=2:ny, i=1:nx) &
         local(js, ka, kb, pres_v, t_v, s_v, rho_v, dsv_dt, dsv_ds, &
               drdt, drds, drdjA, drdjB, drdkL, drdkR, &
               hg2A, hg2B, hg2L, hg2R, haA, haB, haL, haR, &
               dzaL, dzaR, wtA, wtB, wtL, wtR, drdy, drdz, &
               mag2, slope, presS, presN, mask)
         js = j - 1
         ka = k
         kb = k - 1

         presS = pressure_above_x(nx, ny, nz, h_layer, i, js, ka, rho0)
         presN = pressure_above_x(nx, ny, nz, h_layer, i, j, ka, rho0)
         pres_v = 0.5_wp*(presS + presN)

         t_v = 0.25_wp*((t_fill(i, js, ka) + t_fill(i, j, ka)) + &
                        (t_fill(i, js, kb) + t_fill(i, j, kb)))
         s_v = 0.25_wp*((s_fill(i, js, ka) + s_fill(i, j, ka)) + &
                        (s_fill(i, js, kb) + s_fill(i, j, kb)))

         rho_v = eos_density_point(eos, t_v, s_v, pres_v)
         call eos_specvol_derivs(eos, t_v, s_v, pres_v, dsv_dt, dsv_ds)
         drdt = -(rho_v*rho_v)*dsv_dt
         drds = -(rho_v*rho_v)*dsv_ds

         drdjA = drdt*(t_fill(i, j, ka) - t_fill(i, js, ka)) + &
                 drds*(s_fill(i, j, ka) - s_fill(i, js, ka))
         drdjB = drdt*(t_fill(i, j, kb) - t_fill(i, js, kb)) + &
                 drds*(s_fill(i, j, kb) - s_fill(i, js, kb))

         drdkL = drdt*(t_fill(i, js, kb) - t_fill(i, js, ka)) + &
                 drds*(s_fill(i, js, kb) - s_fill(i, js, ka))
         drdkR = drdt*(t_fill(i, j, kb) - t_fill(i, j, ka)) + &
                 drds*(s_fill(i, j, kb) - s_fill(i, j, ka))

         hg2A = h_layer(i, js, ka)*h_layer(i, j, ka) + H_DIV_EPS*H_DIV_EPS
         hg2B = h_layer(i, js, kb)*h_layer(i, j, kb) + H_DIV_EPS*H_DIV_EPS
         hg2L = h_layer(i, js, ka)*h_layer(i, js, kb) + H_DIV_EPS*H_DIV_EPS
         hg2R = h_layer(i, j, ka)*h_layer(i, j, kb) + H_DIV_EPS*H_DIV_EPS
         haA = 0.5_wp*(h_layer(i, js, ka) + h_layer(i, j, ka)) + H_DIV_EPS
         haB = 0.5_wp*(h_layer(i, js, kb) + h_layer(i, j, kb)) + H_DIV_EPS
         haL = 0.5_wp*(h_layer(i, js, ka) + h_layer(i, js, kb)) + H_DIV_EPS
         haR = 0.5_wp*(h_layer(i, j, ka) + h_layer(i, j, kb)) + H_DIV_EPS
         dzaL = max(haL, min_dz)
         dzaR = max(haR, min_dz)
         wtA = hg2A*haB
         wtB = hg2B*haA
         wtL = hg2L*(haR*dzaR)
         wtR = hg2R*(haL*dzaL)

         drdz = ((wtL*drdkL) + (wtR*drdkR))/((dzaL*wtL) + (dzaR*wtR))

         drdy = ((wtA*drdjA + wtB*drdjB)/(wtA + wtB) - &
                 drdz*(e_int(i, js, k) - e_int(i, j, k)))*idyCv(i, j)

         mag2 = drdy*drdy + drdz*drdz
         if (mag2 > 0.0_wp) then
            slope = drdy/sqrt(mag2)
         else
            slope = 0.0_wp
         end if

         mask = wet_v(i, j)
         slope_y(i, j, k) = slope*mask
         n2_v(i, j, k) = g_rho0*drdz*mask
      end do
   end subroutine ocean_slopes_pass_y

   pure function ocean_slopes_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the isopycnal slopes slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_slopes_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%slope_x) &
               + arr_bytes(this%slope_y) &
               + arr_bytes(this%n2_u) &
               + arr_bytes(this%n2_v) &
               + arr_bytes(this%t_fill) &
               + arr_bytes(this%s_fill) &
               + arr_bytes(this%e_int)
   end function ocean_slopes_bytes

end module rdb_ocean_isopycnal_slopes
