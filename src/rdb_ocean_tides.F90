!! Ocean tidal forcing state.
module rdb_ocean_tides
   !! Equilibrium (astronomical) body-force tidal forcing state for the
   !! ocean dyn-core (capability C1).  Fills a GPU-resident equilibrium
   !! tide elevation `eta_eq(x,y)` from a small set of harmonic
   !! constituents; the barotropic momentum solve then drives
   !! `-g grad(eta - eta_forcing)` (pure surface body force), where
   !! `eta_forcing = eta_eq + eta_sal` folds in the scalar self-attraction
   !! & loading (C2) surface elevation `eta_sal = beta_sal*eta` (Ray 1998;
   !! Accad & Pekeris 1978).  With `eta_sal = beta*eta` the surface term
   !! becomes the effective-gravity `-g(1-beta) grad(eta)`; `beta` is lagged
   !! one outer step (uses the stage-start barotropic `eta`).  SAL is
   !! opt-in (`use_sal`, default off) — off ⇒ `eta_forcing == eta_eq`,
   !! bit-identical to C1.
   !! MOM6 divergence (intentional): MOM6 scalar SAL scales the whole
   !! `(eta - eta_eq)` by `(1-beta)` (its `dgeo_de`), damping the body tide by
   !! `beta` too; we apply the Accad-Pekeris load `eta_sal = beta*eta` to the
   !! ocean surface only (body tide at full strength).  Both are valid scalar
   !! approximations; they differ by `g*beta*grad(eta_eq)` (~9% of the tidal
   !! forcing at `beta=0.09`).  Internal-tide drag (C4) and OBC-tide
   !! reconciliation (C3) remain out of scope; the dead `use_itd`/
   !! `itd_coeff`/`itd_global_scale` scaffolding for C4 was removed (PR-8)
   !! — PR-29 (barotropic linear wave drag) lands its own
   !! `lwd_drag_u/v` map on a different type instead.
   !!
   !! Per-outer-step update decomposes the sum-over-constituents into a
   !! host-side scalar update (`amp_cos`, `amp_sin`; `nconst` cos/sin
   !! calls) times a precomputed device-resident spatial structure
   !! (`cos_struct`, `sin_struct`, built once at init from lat/lon), so
   !! the device kernel does NO per-cell trigonometry and NO reduction —
   !! 2 mul + 1 add per constituent per cell.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_tide_astro, only: TIDE_SPECIES, TIDE_AMP, TIDE_LOVE, &
                                   TIDE_OMEGA, TIDES_CATALOG_SIZE, &
                                   equilibrium_arguments, nodal_fu, &
                                   TIDE_DEG2RAD
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_tides_t
   public :: tides_configure_astronomy, tides_build_struct
   public :: tides_update_eta_eq, tides_update_eta_sal

   integer, parameter, public :: TIDES_NCONST_DEFAULT = 8
      !! Standard constituents: M2 S2 N2 K2 K1 O1 P1 Q1.

   type :: ocean_tides_t
      logical :: is_init = .false.
         !! True between `init` and `destroy` (tracks GPU attachment too).

      ! ---- Master switches ----
      logical :: enable = .false.
         !! Master switch (default off => bit-identical).
      logical :: use_sal = .false.
         !! Apply scalar self-attraction & loading (C2).

      ! ---- Active constituent catalog (nconst <= TIDES_CATALOG_SIZE) ----
      integer :: nconst = 0
         !! Number of active harmonic constituents.
      integer, allocatable :: species_c(:)
         !! (nconst) structure-slice index 1/2/3 (diurnal/semidi/long-per).
      real(wp), allocatable :: omega_c(:)
         !! (nconst) angular frequencies (rad/s).
      real(wp), allocatable :: amp_c(:)
         !! (nconst) equilibrium amplitudes A (m).
      real(wp), allocatable :: love_c(:)
         !! (nconst) Love-number factors.
      real(wp), allocatable :: phase0(:)
         !! (nconst) equilibrium argument V_c at ref_date (rad).
      real(wp), allocatable :: f_nodal(:)
         !! (nconst) nodal amplitude factor (fixed at nodal_ref_date).
      real(wp), allocatable :: u_nodal(:)
         !! (nconst) nodal phase (rad).
      real(wp), allocatable :: amp_cos(:)
         !! (nconst) per-step scratch A*love*f*cos(omega*now+V+u).
      real(wp), allocatable :: amp_sin(:)
         !! (nconst) per-step scratch A*love*f*sin(...).

      ! ---- Spatial structure (built ONCE at init) ----
      real(wp), allocatable :: cos_struct(:, :, :)
         !! (nx,ny,3) cos-part spatial structure per species slice.
      real(wp), allocatable :: sin_struct(:, :, :)
         !! (nx,ny,3) sin-part spatial structure per species slice.

      ! ---- 2D forcing fields (filled once per OUTER step) ----
      real(wp), allocatable :: eta_eq(:, :)
         !! (nx,ny) equilibrium tide elevation (m) at cell centres.
      real(wp), allocatable :: eta_sal(:, :)
         !! (nx,ny) scalar-SAL elevation (m); `beta_sal*eta` (C2).
      real(wp), allocatable :: eta_forcing(:, :)
         !! (nx,ny) combined seam field `eta_eq + eta_sal` — the surface
         !! elevation the barotropic PGF drives `-g grad(eta - .)` against.
      real(wp) :: beta_sal = 0.0_wp
         !! Scalar SAL factor beta (~0.085-0.12; C2).

      ! ---- Time tracker ----
      real(wp) :: t_epoch = 0.0_wp
         !! Seconds from ref_date to the model's t=0 (v1: 0).
   contains
      procedure, non_overridable :: init => ocean_tides_init
      procedure, non_overridable :: destroy => ocean_tides_destroy
      procedure, non_overridable :: enter_data => ocean_tides_enter_data
      procedure, non_overridable :: exit_data => ocean_tides_exit_data
      procedure, non_overridable :: bytes => ocean_tides_bytes
   end type ocean_tides_t

contains

   subroutine ocean_tides_init(this, grid)
      !! Minimal init — the real allocation + astronomy fill happens in
      !! `tides_configure_astronomy` / `tides_build_struct` once the
      !! namelist + metrics are available (host, before enter_data).
      class(ocean_tides_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (.false.) this%beta_sal = real(grid%nx_total, wp)
      this%is_init = .true.
   end subroutine ocean_tides_init

   subroutine ocean_tides_destroy(this)
      class(ocean_tides_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%species_c)) deallocate (this%species_c)
      if (allocated(this%omega_c)) deallocate (this%omega_c)
      if (allocated(this%amp_c)) deallocate (this%amp_c)
      if (allocated(this%love_c)) deallocate (this%love_c)
      if (allocated(this%phase0)) deallocate (this%phase0)
      if (allocated(this%f_nodal)) deallocate (this%f_nodal)
      if (allocated(this%u_nodal)) deallocate (this%u_nodal)
      if (allocated(this%amp_cos)) deallocate (this%amp_cos)
      if (allocated(this%amp_sin)) deallocate (this%amp_sin)
      if (allocated(this%cos_struct)) deallocate (this%cos_struct)
      if (allocated(this%sin_struct)) deallocate (this%sin_struct)
      if (allocated(this%eta_eq)) deallocate (this%eta_eq)
      if (allocated(this%eta_sal)) deallocate (this%eta_sal)
      if (allocated(this%eta_forcing)) deallocate (this%eta_forcing)
   end subroutine ocean_tides_destroy

   subroutine tides_configure_astronomy(this, cat_idx, nconst, dref, dnodal, &
                                        add_nodal, nx, ny)
      !! Allocate the active-constituent arrays + 2D fields and fill the
      !! catalog copies + astronomy (phase0, nodal f/u) at the reference
      !! and nodal reference day numbers.  Host-side setup (before
      !! enter_data).  `cat_idx(1:nconst)` are catalog indices
      !! (1..TIDES_CATALOG_SIZE).
      class(ocean_tides_t), intent(inout) :: this
      integer, intent(in) :: nconst
      integer, intent(in) :: cat_idx(nconst)
      real(wp), intent(in) :: dref, dnodal
      logical, intent(in) :: add_nodal
      integer, intent(in) :: nx, ny
      real(wp) :: v_all(TIDES_CATALOG_SIZE)
      real(wp) :: f_all(TIDES_CATALOG_SIZE), u_all(TIDES_CATALOG_SIZE)
      integer :: c, ic

      this%nconst = nconst
      if (allocated(this%species_c)) deallocate (this%species_c)
      allocate (this%species_c(nconst))
      allocate (this%omega_c(nconst), this%amp_c(nconst), this%love_c(nconst))
      allocate (this%phase0(nconst), this%f_nodal(nconst), this%u_nodal(nconst))
      allocate (this%amp_cos(nconst), this%amp_sin(nconst))

      call equilibrium_arguments(dref, v_all)
      call nodal_fu(dnodal, add_nodal, f_all, u_all)

      do c = 1, nconst
         ic = cat_idx(c)
         this%species_c(c) = TIDE_SPECIES(ic)
         this%omega_c(c) = TIDE_OMEGA(ic)
         this%amp_c(c) = TIDE_AMP(ic)
         this%love_c(c) = TIDE_LOVE(ic)
         this%phase0(c) = v_all(ic)
         this%f_nodal(c) = f_all(ic)
         this%u_nodal(c) = u_all(ic)
         this%amp_cos(c) = 0.0_wp
         this%amp_sin(c) = 0.0_wp
      end do

      if (allocated(this%cos_struct)) deallocate (this%cos_struct)
      allocate (this%cos_struct(nx, ny, 3), source=0.0_wp)
      allocate (this%sin_struct(nx, ny, 3), source=0.0_wp)
      allocate (this%eta_eq(nx, ny), source=0.0_wp)
      allocate (this%eta_sal(nx, ny), source=0.0_wp)
      allocate (this%eta_forcing(nx, ny), source=0.0_wp)
   end subroutine tides_configure_astronomy

   subroutine tides_build_struct(this, geolat, geolon, nx, ny)
      !! Build the (nx,ny,3) cos/sin spatial-structure arrays from cell-
      !! centre latitude/longitude (degrees), via the angle-sum fold
      !! cos(theta + n*lambda) = cos(theta)cos(n*lambda) - sin(theta)sin(n*lambda).
      !! Plain host loop over all cells incl. ghosts (before enter_data).
      !!   slice 1 diurnal    (n=1): G1 = sin(2 phi)
      !!   slice 2 semidiurnal(n=2): G2 = cos^2 phi
      !!   slice 3 long-period(n=0): G0 = 1/2 - 3/2 sin^2 phi
      class(ocean_tides_t), intent(inout) :: this
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: geolat(nx, ny), geolon(nx, ny)
      integer :: i, j
      real(wp) :: phi, lam, s2phi, c2phi

      do j = 1, ny
         do i = 1, nx
            phi = geolat(i, j)*TIDE_DEG2RAD
            lam = geolon(i, j)*TIDE_DEG2RAD
            s2phi = sin(2.0_wp*phi)
            c2phi = cos(phi)**2
            ! slice 1 diurnal (n=1)
            this%cos_struct(i, j, 1) = s2phi*cos(lam)
            this%sin_struct(i, j, 1) = -s2phi*sin(lam)
            ! slice 2 semidiurnal (n=2)
            this%cos_struct(i, j, 2) = c2phi*cos(2.0_wp*lam)
            this%sin_struct(i, j, 2) = -c2phi*sin(2.0_wp*lam)
            ! slice 3 long-period (n=0)
            this%cos_struct(i, j, 3) = 0.5_wp - 1.5_wp*sin(phi)**2
            this%sin_struct(i, j, 3) = 0.0_wp
         end do
      end do
   end subroutine tides_build_struct

   subroutine tides_update_eta_eq(this, t)
      !! Refresh `eta_eq(x,y)` for the current outer-step time `t` (s).
      !! Host recomputes the `nconst` amplitude scalars, pushes them to
      !! the device, then a `do concurrent` fills `eta_eq` with no
      !! per-cell trig and no reduction.  Held static across the inner
      !! barotropic substep loop.
      type(ocean_tides_t), intent(inout) :: this
      real(wp), intent(in) :: t
      real(wp) :: now, ang, pre
      integer :: c, nx, ny

      now = t + this%t_epoch
      do c = 1, this%nconst
         ang = this%omega_c(c)*now + this%phase0(c) + this%u_nodal(c)
         pre = this%amp_c(c)*this%love_c(c)*this%f_nodal(c)
         this%amp_cos(c) = pre*cos(ang)
         this%amp_sin(c) = pre*sin(ang)
      end do
      !$acc update device(this%amp_cos, this%amp_sin)

      nx = size(this%eta_eq, 1)
      ny = size(this%eta_eq, 2)
      call tides_update_eta_eq_impl(nx, ny, this%nconst, this%species_c, &
                                    this%amp_cos, this%amp_sin, &
                                    this%cos_struct, this%sin_struct, this%eta_eq)
   end subroutine tides_update_eta_eq

   subroutine tides_update_eta_eq_impl(nx, ny, nconst, species_c, amp_cos, &
                                       amp_sin, cos_struct, sin_struct, eta_eq)
      !! Flat-impl device fill (explicit-shape dummies — no descriptor
      !! walk).  eta_eq(i,j) = sum_c amp_cos(c)*cos_struct(i,j,m)
      !!                            + amp_sin(c)*sin_struct(i,j,m),
      !! m = species_c(c).  Contiguous index (i) innermost.
      integer, intent(in) :: nx, ny, nconst
      integer, intent(in) :: species_c(nconst)
      real(wp), intent(in) :: amp_cos(nconst), amp_sin(nconst)
      real(wp), intent(in) :: cos_struct(nx, ny, 3), sin_struct(nx, ny, 3)
      real(wp), intent(out) :: eta_eq(nx, ny)
      integer :: i, j, c, m
      real(wp) :: acc

      do concurrent(j=1:ny, i=1:nx) local(acc, c, m)
         acc = 0.0_wp
         do c = 1, nconst
            m = species_c(c)
            acc = acc + amp_cos(c)*cos_struct(i, j, m) + amp_sin(c)*sin_struct(i, j, m)
         end do
         eta_eq(i, j) = acc
      end do
   end subroutine tides_update_eta_eq_impl

   subroutine tides_update_eta_sal(this, eta_current)
      !! Refresh the combined seam field `eta_forcing` for the current
      !! outer step.  Scalar self-attraction & loading (C2): when
      !! `use_sal`, `eta_sal = beta_sal*eta_current` and
      !! `eta_forcing = eta_eq + eta_sal`; otherwise `eta_forcing = eta_eq`
      !! (bit-identical to C1).  `eta_current` is the lagged (stage-start,
      !! previous outer step) barotropic surface elevation.  Host does no
      !! work; the fill is a single explicit-shape `do concurrent`.  Must
      !! be called AFTER `tides_update_eta_eq` (reads the fresh `eta_eq`).
      type(ocean_tides_t), intent(inout) :: this
      real(wp), intent(in) :: eta_current(:, :)
      integer :: nx, ny

      nx = size(this%eta_eq, 1)
      ny = size(this%eta_eq, 2)
      call tides_update_eta_sal_impl(nx, ny, this%use_sal, this%beta_sal, &
                                     eta_current, this%eta_eq, this%eta_sal, &
                                     this%eta_forcing)
   end subroutine tides_update_eta_sal

   subroutine tides_update_eta_sal_impl(nx, ny, use_sal, beta_sal, eta_current, &
                                        eta_eq, eta_sal, eta_forcing)
      !! Flat-impl device fill (explicit-shape dummies).  Loop-invariant
      !! `use_sal` branch kept INSIDE the single `do concurrent` (one
      !! launch, uniform branch is ~free).  Off ⇒ pure copy of `eta_eq`
      !! into `eta_forcing` ⇒ bit-identical.  Contiguous index (i) innermost.
      integer, intent(in) :: nx, ny
      logical, intent(in) :: use_sal
      real(wp), intent(in) :: beta_sal
      real(wp), intent(in) :: eta_current(nx, ny), eta_eq(nx, ny)
      real(wp), intent(inout) :: eta_sal(nx, ny)
      real(wp), intent(out) :: eta_forcing(nx, ny)
      integer :: i, j

      ! Write straight into the arrays (no `local` scalar): reading back
      ! `eta_sal(i,j)` in the same iteration is a within-iteration RAW
      ! (legal in do concurrent) and sidesteps the gfortran `local()`
      ! if/else codegen artefact that perturbs the product sub-ULP.
      do concurrent(j=1:ny, i=1:nx)
         if (use_sal) then
            eta_sal(i, j) = beta_sal*eta_current(i, j)
            eta_forcing(i, j) = eta_eq(i, j) + eta_sal(i, j)
         else
            eta_forcing(i, j) = eta_eq(i, j)
         end if
      end do
   end subroutine tides_update_eta_sal_impl

   subroutine ocean_tides_enter_data(this)
      !! Attach the device-resident tide arrays.  Only when enabled.
      !! select-type -> non-poly `_impl` (AMD libomptarget class-box rule).
      class(ocean_tides_t), intent(inout) :: this
      if (.not. this%enable) return
      select type (this)
      type is (ocean_tides_t)
         call ocean_tides_enter_data_impl(this)
      end select
   end subroutine ocean_tides_enter_data

   subroutine ocean_tides_enter_data_impl(this)
      type(ocean_tides_t), intent(inout) :: this
      if (allocated(this%species_c)) then
         !$acc enter data copyin(this%species_c, this%omega_c, this%amp_c, &
         !$acc                   this%love_c, this%phase0, this%f_nodal, &
         !$acc                   this%u_nodal, this%amp_cos, this%amp_sin, &
         !$acc                   this%cos_struct, this%sin_struct, this%eta_eq, &
         !$acc                   this%eta_sal, this%eta_forcing)
      end if
   end subroutine ocean_tides_enter_data_impl

   subroutine ocean_tides_exit_data(this)
      class(ocean_tides_t), intent(inout) :: this
      if (.not. this%enable) return
      select type (this)
      type is (ocean_tides_t)
         call ocean_tides_exit_data_impl(this)
      end select
   end subroutine ocean_tides_exit_data

   subroutine ocean_tides_exit_data_impl(this)
      type(ocean_tides_t), intent(inout) :: this
      if (allocated(this%species_c)) then
         !$acc exit data delete(this%eta_forcing, this%eta_sal, &
         !$acc                  this%eta_eq, this%sin_struct, this%cos_struct, &
         !$acc                  this%amp_sin, this%amp_cos, this%u_nodal, &
         !$acc                  this%f_nodal, this%phase0, this%love_c, &
         !$acc                  this%amp_c, this%omega_c, this%species_c)
      end if
   end subroutine ocean_tides_exit_data_impl

   pure function ocean_tides_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the tides slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_tides_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%species_c) &
               + arr_bytes(this%omega_c) &
               + arr_bytes(this%amp_c) &
               + arr_bytes(this%love_c) &
               + arr_bytes(this%phase0) &
               + arr_bytes(this%f_nodal) &
               + arr_bytes(this%u_nodal) &
               + arr_bytes(this%amp_cos) &
               + arr_bytes(this%amp_sin) &
               + arr_bytes(this%cos_struct) &
               + arr_bytes(this%sin_struct) &
               + arr_bytes(this%eta_eq) &
               + arr_bytes(this%eta_sal) &
               + arr_bytes(this%eta_forcing)
   end function ocean_tides_bytes

end module rdb_ocean_tides
