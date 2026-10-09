!! Barotropic (depth-integrated) solution state on an Arakawa C-grid.
module rdb_barotropic_state
   !! Depth-integrated state on an Arakawa C-grid for the ocean dynamical core.
   !! C-grid layout:
   !!   * Scalars (h, b, flux divergence) at cell centres, (nx_total, ny_total).
   !!   * x-velocity/momentum at east faces, (nx_total+1, ny_total).
   !!   * y-velocity/momentum at north faces, (nx_total, ny_total+1).
   !! Lives at `ocean_state%barotropic` for `sim_type='ocean'`.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: barotropic_state_t

   type :: barotropic_state_t
      !! Depth-integrated state on an Arakawa C-grid (face-located velocities as
      !! continuity-PPM produces/consumes; cell-centred h-update by flux div).

      logical :: is_init = .false.
         !! True between `init` and `destroy`. Prefer to `allocated(...)` — it
         !! also tracks GPU device attachment, which `allocated` cannot see.

      ! ---- Cell-centred scalars ----
      real(wp), allocatable :: h(:, :)
         !! Water depth at cell centre (m), shape (nx_total, ny_total)
      real(wp), allocatable :: b(:, :)
         !! Reference column DEPTH at cell centre (m, **positive down**):
         !! a 4000 m deep cell holds `b = +4000`, land holds `b = 0`, and a
         !! cell is land iff `b < LAND_DEPTH_THRESHOLD` (2 m,
         !! `rdb_constants`).  Fixed by the SSH identity
         !! `eta = sum_k h_layer - b` (`rdb_ocean_diag_fills.F90:157`) and by
         !! `barotropic%h = b + eta` (`rdb_ocean_state.F90:2482`).
         !!
         !! NOTE the opposite sign of the topographic HEIGHTS used by the
         !! porous-barrier statistics (`rdb_ocean_porous.F90:17-21`) and by
         !! `metrics%bed_height` (`rdb_ocean_metrics.F90:151`, `= -b`), and of
         !! GEBCO/ETOPO products, which give a negative height.  Handing a
         !! negative-height array in unflipped puts every cell below the land
         !! threshold and yields a silent, crash-free, all-land run.
         !! (This docstring previously read "positive up" — a coastal-path
         !! leftover; the coastal regime left this tree in the ocean-only
         !! split.)

      ! ---- Face-located velocities / momentum ----
      ! u on east face: u_face_x(i,j) at interface (i,j)|(i+1,j), shape
      ! (nx_total+1, ny_total) so face index runs 1 (west wall)..nx+1 (east wall).
      real(wp), allocatable :: u_face_x(:, :)
         !! x-velocity at east face of each cell (m/s)
      real(wp), allocatable :: hu_face_x(:, :)
         !! Face-thickness x momentum h*u (m^2/s)

      ! v on north face: v_face_y(i,j) at interface (i,j)|(i,j+1),
      ! shape (nx_total, ny_total+1).
      real(wp), allocatable :: v_face_y(:, :)
         !! y-velocity at north face of each cell (m/s)
      real(wp), allocatable :: hv_face_y(:, :)
         !! Face-thickness y momentum h*v (m^2/s)

      ! ---- Per-face mass fluxes (continuity-PPM primary output) ----
      ! PPM mass flux at east/north faces; consumed by tracer advection under
      ! CWC. Same shape as the face-velocity arrays.
      real(wp), allocatable :: mass_flux_x(:, :)
      real(wp), allocatable :: mass_flux_y(:, :)

      ! ---- Flux divergence (workspace) ----
      ! Cell-centred net divergence used to step h forward.
      real(wp), allocatable :: flux_h(:, :)

      ! ---- RK2 save buffers ----
      ! State at start of the outer RK step; the inner substeps run N times
      ! from this snapshot.
      real(wp), allocatable :: h0(:, :)
      real(wp), allocatable :: u_face_x0(:, :)
      real(wp), allocatable :: v_face_y0(:, :)

      ! ---- Physics scalars ----
      real(wp) :: manning_n = 0.0_wp
         !! Manning roughness for shelf-side bottom drag.
      real(wp) :: coriolis_f = 0.0_wp
         !! Beta-plane / f-plane Coriolis parameter (1/s).
   contains
      procedure, non_overridable :: init => barotropic_state_init
      procedure, non_overridable :: destroy => barotropic_state_destroy
      procedure, non_overridable :: enter_data => barotropic_state_enter_data
      procedure, non_overridable :: exit_data => barotropic_state_exit_data
      procedure, non_overridable :: bytes => barotropic_state_bytes
   end type barotropic_state_t

contains

   subroutine barotropic_state_init(this, grid)
      !! Allocate every C-grid barotropic array zero-filled. East-face arrays
      !! are (nx+1, ny), north-face (nx, ny+1). Scalars (manning_n, coriolis_f)
      !! are populated by `state_init_from_config` afterwards.
      class(barotropic_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      ! Cell-centred scalars
      allocate (this%h(nx, ny), source=0.0_wp)
      allocate (this%b(nx, ny), source=0.0_wp)
      allocate (this%flux_h(nx, ny), source=0.0_wp)

      ! East-face arrays: extra column, face index i runs 1..nx+1.
      allocate (this%u_face_x(nx + 1, ny), source=0.0_wp)
      allocate (this%hu_face_x(nx + 1, ny), source=0.0_wp)
      allocate (this%mass_flux_x(nx + 1, ny), source=0.0_wp)

      ! North-face arrays: extra row, face index j runs 1..ny+1.
      allocate (this%v_face_y(nx, ny + 1), source=0.0_wp)
      allocate (this%hv_face_y(nx, ny + 1), source=0.0_wp)
      allocate (this%mass_flux_y(nx, ny + 1), source=0.0_wp)

      ! RK2 save buffers
      allocate (this%h0(nx, ny), source=0.0_wp)
      allocate (this%u_face_x0(nx + 1, ny), source=0.0_wp)
      allocate (this%v_face_y0(nx, ny + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine barotropic_state_init

   subroutine barotropic_state_destroy(this)
      !! Tear down host allocations. Call `exit_data` first — destroying a
      !! still-device-mapped state leaks device memory silently.
      class(barotropic_state_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%h)) deallocate (this%h)
      if (allocated(this%b)) deallocate (this%b)
      if (allocated(this%u_face_x)) deallocate (this%u_face_x)
      if (allocated(this%hu_face_x)) deallocate (this%hu_face_x)
      if (allocated(this%v_face_y)) deallocate (this%v_face_y)
      if (allocated(this%hv_face_y)) deallocate (this%hv_face_y)
      if (allocated(this%mass_flux_x)) deallocate (this%mass_flux_x)
      if (allocated(this%mass_flux_y)) deallocate (this%mass_flux_y)
      if (allocated(this%flux_h)) deallocate (this%flux_h)
      if (allocated(this%h0)) deallocate (this%h0)
      if (allocated(this%u_face_x0)) deallocate (this%u_face_x0)
      if (allocated(this%v_face_y0)) deallocate (this%v_face_y0)
   end subroutine barotropic_state_destroy

   pure function barotropic_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the ocean C-grid barotropic slot.
      class(barotropic_state_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%h) + arr_bytes(this%b) &
               + arr_bytes(this%u_face_x) + arr_bytes(this%hu_face_x) &
               + arr_bytes(this%v_face_y) + arr_bytes(this%hv_face_y) &
               + arr_bytes(this%mass_flux_x) + arr_bytes(this%mass_flux_y) &
               + arr_bytes(this%flux_h) + arr_bytes(this%h0) &
               + arr_bytes(this%u_face_x0) + arr_bytes(this%v_face_y0)
   end function barotropic_state_bytes

   subroutine barotropic_state_enter_data(this)
      !! Attach the C-grid barotropic allocatables to the device. The parent
      !! struct is mapped by the orchestrating routine. Read+write fields use
      !! `copyin`; pure-workspace fields use `create`.
      class(barotropic_state_t), intent(inout) :: this
      select type (this)
      type is (barotropic_state_t)
         call barotropic_state_enter_data_impl(this)
      end select
   end subroutine barotropic_state_enter_data

   subroutine barotropic_state_enter_data_impl(this)
      type(barotropic_state_t), intent(inout) :: this

      !$acc enter data copyin(this%h, this%b, &
      !$acc&                  this%u_face_x, this%hu_face_x, &
      !$acc&                  this%v_face_y, this%hv_face_y, &
      !$acc&                  this%h0, this%u_face_x0, this%v_face_y0)
      !$acc enter data create(this%mass_flux_x, this%mass_flux_y, &
      !$acc&                  this%flux_h)
   end subroutine barotropic_state_enter_data_impl

   subroutine barotropic_state_exit_data(this)
      !! Reverse of `enter_data` — copy out prognostic fields, drop scratch.
      !! Caller detaches the parent struct after this returns.
      class(barotropic_state_t), intent(inout) :: this
      select type (this)
      type is (barotropic_state_t)
         call barotropic_state_exit_data_impl(this)
      end select
   end subroutine barotropic_state_exit_data

   subroutine barotropic_state_exit_data_impl(this)
      type(barotropic_state_t), intent(inout) :: this

      !$acc exit data copyout(this%h, this%u_face_x, this%v_face_y, &
      !$acc&                  this%hu_face_x, this%hv_face_y)
      !$acc exit data delete(this%b, this%h0, this%u_face_x0, this%v_face_y0, &
      !$acc&                 this%mass_flux_x, this%mass_flux_y, this%flux_h)
   end subroutine barotropic_state_exit_data_impl

end module rdb_barotropic_state
