!! Ocean atmospheric surface-pressure loading state (inverse barometer).
module rdb_ocean_p_surf
   !! Surface-pressure loading for the ocean C-grid dyn-core (PR-17).  An
   !! atmospheric surface pressure `p_surf(x,y)` (Pa) adds a depth-uniform
   !! acceleration `-(1/rho0) grad(p_surf)` to the horizontal momentum
   !! (Wunsch & Stammer 1997).  Because the term is depth-independent it is
   !! a purely BAROTROPIC forcing: it drives the free surface and, at rest,
   !! integrates to the equilibrium inverse-barometer response
   !! `eta_ib = -p_surf/(rho0 g) + const` (~1 cm depression per hPa of
   !! atmospheric high; Ponte 2006 for the closed-domain mean removal).
   !!
   !! It lands on the SAME `eta_forcing` seam the equilibrium body tide and
   !! scalar SAL compose through: the barotropic substep drives
   !! `-g grad(eta - eta_forcing)`, so folding
   !! `eta_ib = -p_surf/(rho0 g_bt)` into the seam gives the momentum an
   !! extra `-g grad(eta) - (1/rho0) grad(p_surf)` exactly.  `eta_ib` and the
   !! combined seam field `eta_seam = eta_ib [+ eta_tide]` are cell-centred
   !! 2-D fields in metres, differenced on the identical two-point stencil
   !! the SSH uses, so the discrete inverse-barometer balance
   !! `eta = eta_ib + const` is an exact discrete steady state.
   !!
   !! Sign: `eta_ib = -p_surf/(rho0 g)` — a high depresses the surface, a low
   !! bulges it.  Gravity is `bt_work%g_bt` (what the barotropic substep that
   !! consumes the seam actually runs on), NOT `rdb_constants:GRAVITY`; the
   !! caller passes it in so `eta_ib` scales with the ocean's own dynamics.
   !!
   !! `p_surf` (the assembled total, Pa, `>= 0`) is consumed read-only from
   !! the surface-forcing type (`ocean_surface_flux_t`) — this module adds no
   !! forcing field of its own, and does not write `p_surf` (the dyn step
   !! holds `sf` `intent(in)`).  With no ice mass-loading the assembly
   !! `p_surf = p_surf_atm` is a configure-time seed (`set_p_surf_const`
   !! fills both) since `p_surf_atm` is static; when sea-ice loading lands
   !! (PR-18) its `ice_ocean_mass_load` overwrites `p_surf = p_surf_atm +
   !! g_load*mis` via its own `inout` access earlier in the same outer step,
   !! and this module reads whatever `p_surf` then holds.  Adcroft et al.
   !! (2019) levitating surface-load convention.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_p_surf_t
   public :: p_surf_configure, p_surf_update_seam

   type :: ocean_p_surf_t
      logical :: is_init = .false.
         !! True between `init` and `destroy` (tracks GPU attachment too).
      logical :: enable = .false.
         !! Master switch (default off => bit-identical).
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).  Assigned from
         !! `ocean_state%eos%rho0` at configure — the single ρ₀ of record.
      real(wp), allocatable :: eta_ib(:, :)
         !! (nx,ny) inverse-barometer elevation `-p_surf/(rho0 g_bt)`, metres.
      real(wp), allocatable :: eta_seam(:, :)
         !! (nx,ny) combined seam field `eta_ib [+ eta_tide]` — the surface
         !! elevation the barotropic PGF drives `-g grad(eta - .)` against.
   contains
      procedure, non_overridable :: init => ocean_p_surf_init
      procedure, non_overridable :: destroy => ocean_p_surf_destroy
      procedure, non_overridable :: enter_data => ocean_p_surf_enter_data
      procedure, non_overridable :: exit_data => ocean_p_surf_exit_data
      procedure, non_overridable :: bytes => ocean_p_surf_bytes
   end type ocean_p_surf_t

contains

   subroutine ocean_p_surf_init(this, grid)
      !! Minimal init — the real allocation happens in `p_surf_configure`
      !! once the namelist is available (host, before enter_data).
      class(ocean_p_surf_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (.false.) this%rho0 = real(grid%nx_total, wp)
      this%is_init = .true.
   end subroutine ocean_p_surf_init

   subroutine ocean_p_surf_destroy(this)
      class(ocean_p_surf_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%eta_ib)) deallocate (this%eta_ib)
      if (allocated(this%eta_seam)) deallocate (this%eta_seam)
   end subroutine ocean_p_surf_destroy

   subroutine p_surf_configure(this, nx, ny)
      !! Allocate `eta_ib`/`eta_seam` (host, before enter_data).  No-op when
      !! `.not. enable` (bit-identical; zero extra device memory).  Mirrors
      !! `tides_configure_astronomy`.
      class(ocean_p_surf_t), intent(inout) :: this
      integer, intent(in) :: nx, ny
      if (.not. this%enable) return
      if (allocated(this%eta_ib)) deallocate (this%eta_ib)
      if (allocated(this%eta_seam)) deallocate (this%eta_seam)
      allocate (this%eta_ib(nx, ny), source=0.0_wp)
      allocate (this%eta_seam(nx, ny), source=0.0_wp)
   end subroutine p_surf_configure

   subroutine p_surf_update_seam(this, p_surf, g_bt, eta_tide)
      !! Refresh the combined seam field `eta_seam` for the current outer
      !! step from the assembled total surface pressure `p_surf` (read-only;
      !! see the module header for where the assembly happens): fold
      !! `eta_ib = -p_surf/(rho0 g_bt)` into `eta_seam = eta_ib [+ eta_tide]`.
      !! `g_bt` is the barotropic-substep gravity (`bt_work%g_bt`), passed in
      !! at call time so the seam scales with the ocean's own dynamics.
      !! `eta_tide` (optional) is the equilibrium-tide + scalar-SAL seam,
      !! added when the tide is on.  Host does one reciprocal; the fill is a
      !! single explicit-shape `do concurrent`.  `p_surf` is the
      !! device-resident `ocean_surface_flux_t%p_surf` component array.
      class(ocean_p_surf_t), intent(inout) :: this
      real(wp), intent(in) :: p_surf(:, :)
      real(wp), intent(in) :: g_bt
      real(wp), intent(in), optional :: eta_tide(:, :)
      real(wp) :: i_rho0_g
      integer :: nx, ny

      nx = size(this%eta_ib, 1)
      ny = size(this%eta_ib, 2)
      i_rho0_g = 1.0_wp/(this%rho0*g_bt)
      if (present(eta_tide)) then
         call p_surf_update_seam_impl(nx, ny, i_rho0_g, p_surf, &
                                      this%eta_ib, this%eta_seam, &
                                      use_tide=.true., eta_tide=eta_tide)
      else
         call p_surf_update_seam_impl(nx, ny, i_rho0_g, p_surf, &
                                      this%eta_ib, this%eta_seam, &
                                      use_tide=.false.)
      end if
   end subroutine p_surf_update_seam

   subroutine p_surf_update_seam_impl(nx, ny, i_rho0_g, p_surf, &
                                      eta_ib, eta_seam, use_tide, eta_tide)
      !! Flat-impl device fill (explicit-shape dummies — no descriptor walk).
      !! `eta_ib(i,j) = -p_surf(i,j)*i_rho0_g` and
      !! `eta_seam(i,j) = eta_ib(i,j) [+ eta_tide(i,j)]`.  Loop-invariant
      !! `use_tide` branch kept INSIDE the single `do concurrent` (one
      !! launch, uniform branch is ~free).  Contiguous index (i) innermost.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: i_rho0_g
      real(wp), intent(in) :: p_surf(nx, ny)
      real(wp), intent(out) :: eta_ib(nx, ny)
      real(wp), intent(out) :: eta_seam(nx, ny)
      logical, intent(in) :: use_tide
      real(wp), intent(in), optional :: eta_tide(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         eta_ib(i, j) = -p_surf(i, j)*i_rho0_g
         if (use_tide) then
            eta_seam(i, j) = eta_ib(i, j) + eta_tide(i, j)
         else
            eta_seam(i, j) = eta_ib(i, j)
         end if
      end do
   end subroutine p_surf_update_seam_impl

   subroutine ocean_p_surf_enter_data(this)
      !! Attach the device-resident seam arrays.  Only when enabled.
      !! select-type -> non-poly `_impl` (AMD libomptarget class-box rule).
      class(ocean_p_surf_t), intent(inout) :: this
      if (.not. this%enable) return
      select type (this)
      type is (ocean_p_surf_t)
         call ocean_p_surf_enter_data_impl(this)
      end select
   end subroutine ocean_p_surf_enter_data

   subroutine ocean_p_surf_enter_data_impl(this)
      type(ocean_p_surf_t), intent(inout) :: this
      if (allocated(this%eta_ib)) then
         !$acc enter data copyin(this%eta_ib, this%eta_seam)
      end if
   end subroutine ocean_p_surf_enter_data_impl

   subroutine ocean_p_surf_exit_data(this)
      class(ocean_p_surf_t), intent(inout) :: this
      if (.not. this%enable) return
      select type (this)
      type is (ocean_p_surf_t)
         call ocean_p_surf_exit_data_impl(this)
      end select
   end subroutine ocean_p_surf_exit_data

   subroutine ocean_p_surf_exit_data_impl(this)
      type(ocean_p_surf_t), intent(inout) :: this
      if (allocated(this%eta_ib)) then
         !$acc exit data delete(this%eta_seam, this%eta_ib)
      end if
   end subroutine ocean_p_surf_exit_data_impl

   pure function ocean_p_surf_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the p_surf slot (0 when
      !! unallocated).  One arr_bytes term per array — add a term here when
      !! a new allocatable joins the type.
      class(ocean_p_surf_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%eta_ib) &
               + arr_bytes(this%eta_seam)
   end function ocean_p_surf_bytes

end module rdb_ocean_p_surf
