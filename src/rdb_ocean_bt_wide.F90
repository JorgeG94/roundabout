!! Wide-halo BT march-in state (Phase 3c, D3 v1.1).
module rdb_ocean_bt_wide
   !! Shadow state for the wide-halo barotropic march-in.
   !! When `bt_halo > 0`, the BT fast loop runs on WIDE arrays with ghost
   !! width `ng_wide = nghost + bt_halo`.  Ghost cells outside the valid
   !! band evolve stale data that creeps INWARD at 2 cells/substep; one
   !! grouped exchange every `bt_halo/2` substeps keeps the physical interior
   !! clean.  The mid-substep u exchange is absorbed (within the 2-cell/substep
   !! stencil budget).
   !!
   !! Lifecycle:
   !!   1. `init`         — allocate wide arrays + wide grid/metrics/f_corner.
   !!   2. `enter_data`   — attach all wide arrays to the GPU present table.
   !!   3. Per outer step: `copy_in` -> (entry_exchange + substep) -> `copy_out`
   !!                      -> normal-width exit exchange (in the caller).
   !!   4. `exit_data`    — release GPU present table entries.
   !!   5. `destroy`      — free host memory.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, &
                                metrics_fill_cartesian, metrics_fill_spherical, &
                                metrics_finalize, &
                                GRID_CONFIG_CARTESIAN, GRID_CONFIG_SPHERICAL, &
                                CORIOLIS_SCHEME_BETA_PLANE, &
                                metrics_fill_coriolis
   use rdb_ocean_halo, only: ocean_halo_bt_group_2d_wide, &
                             ocean_halo_centre_2d_wide, &
                             ocean_halo_face_x_2d_wide, &
                             ocean_halo_face_y_2d_wide
   use rdb_profiler, only: profiler_start, profiler_stop
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_mem_report, only: arr_bytes
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: bt_wide_t
   public :: bt_wide_substep

   type :: bt_wide_t
      !! Wide-halo shadow state for the barotropic fast loop.
      logical :: is_init = .false.
         !! True after `init`, before `destroy`.
      integer :: bt_halo = 0
         !! Requested wide-halo width (cells; even, > 0).
      integer :: ng_wide = 0
         !! Effective ghost width: `nghost + bt_halo`.
      integer :: num_cycles = 0
         !! Substeps between grouped wide exchanges: `bt_halo / 2`.

      ! Wide grid and metrics.  Same `nx_phys/ny_phys` as the normal grid;
      ! only `nghost` differs (`ng_wide`).
      type(hgrid_t) :: grid_w
         !! Wide hgrid_t: nghost = ng_wide.
      type(ocean_metrics_t) :: metrics_w
         !! Metrics built on `grid_w` via the same formula generator.

      ! Wide Coriolis corners — filled by `metrics_fill_coriolis`.
      real(wp), allocatable :: f_corner_w(:, :)
         !! Coriolis at wide-grid C-grid corners (nx_w+1, ny_w+1).

      ! Wide shadow copies of the 22 fast-loop 2D arrays.
      ! Naming convention: w_<field> mirrors <field> on the wide grid.
      ! Centres (nx_w, ny_w):
      real(wp), allocatable :: w_eta(:, :)
      real(wp), allocatable :: w_H_ref(:, :)
      real(wp), allocatable :: w_eta_new(:, :)
      real(wp), allocatable :: w_ke(:, :)
      real(wp), allocatable :: w_eta_sum(:, :)
      real(wp), allocatable :: w_eta_end(:, :)
      ! East-face u (nx_w+1, ny_w):
      real(wp), allocatable :: w_ubt(:, :)
      real(wp), allocatable :: w_ubt_prev(:, :)
      real(wp), allocatable :: w_rem_u(:, :)
      real(wp), allocatable :: w_ubt_sum(:, :)
      real(wp), allocatable :: w_uhbt_sum(:, :)
      real(wp), allocatable :: w_uhbt(:, :)
      real(wp), allocatable :: w_ubt_end(:, :)
      ! North-face v (nx_w, ny_w+1):
      real(wp), allocatable :: w_vbt(:, :)
      real(wp), allocatable :: w_vbt_prev(:, :)
      real(wp), allocatable :: w_rem_v(:, :)
      real(wp), allocatable :: w_vbt_sum(:, :)
      real(wp), allocatable :: w_vhbt_sum(:, :)
      real(wp), allocatable :: w_vhbt(:, :)
      real(wp), allocatable :: w_vbt_end(:, :)
      ! Corner zeta (nx_w+1, ny_w+1):
      real(wp), allocatable :: w_zeta(:, :)
      ! Wide force arrays (copied in from the normal-width slow tendencies):
      real(wp), allocatable :: w_force_u(:, :)
      real(wp), allocatable :: w_force_v(:, :)

   contains
      procedure, non_overridable :: init => bt_wide_init
      procedure, non_overridable :: enter_data => bt_wide_enter_data
      procedure, non_overridable :: exit_data => bt_wide_exit_data
      procedure, non_overridable :: copy_in => bt_wide_copy_in
      procedure, non_overridable :: entry_exchange => bt_wide_entry_exchange
      procedure, non_overridable :: copy_out => bt_wide_copy_out
      procedure, non_overridable :: destroy => bt_wide_destroy
      procedure, non_overridable :: bytes => bt_wide_bytes
   end type bt_wide_t

contains

   pure function bt_wide_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the wide-halo BT shadow state
      !! (0 when unallocated, i.e. whenever `&ocean_bt_nml bt_halo = 0`).
      !!
      !! Every term here is device-mapped by `bt_wide_enter_data_impl` and
      !! was previously counted by nothing — `ocean_dyn_bytes` omitted the
      !! slot entirely because this function did not exist.  ~436 MB at
      !! 1000x800 with `bt_halo = 10, nghost = 4`.  One `arr_bytes` term
      !! per array: add one here when a new allocatable joins the type.
      class(bt_wide_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%metrics_w%bytes() &
               + arr_bytes(this%f_corner_w) &
               + arr_bytes(this%w_eta) &
               + arr_bytes(this%w_H_ref) &
               + arr_bytes(this%w_eta_new) &
               + arr_bytes(this%w_ke) &
               + arr_bytes(this%w_eta_sum) &
               + arr_bytes(this%w_eta_end) &
               + arr_bytes(this%w_ubt) &
               + arr_bytes(this%w_ubt_prev) &
               + arr_bytes(this%w_rem_u) &
               + arr_bytes(this%w_ubt_sum) &
               + arr_bytes(this%w_uhbt_sum) &
               + arr_bytes(this%w_uhbt) &
               + arr_bytes(this%w_ubt_end) &
               + arr_bytes(this%w_vbt) &
               + arr_bytes(this%w_vbt_prev) &
               + arr_bytes(this%w_rem_v) &
               + arr_bytes(this%w_vbt_sum) &
               + arr_bytes(this%w_vhbt_sum) &
               + arr_bytes(this%w_vhbt) &
               + arr_bytes(this%w_vbt_end) &
               + arr_bytes(this%w_zeta) &
               + arr_bytes(this%w_force_u) &
               + arr_bytes(this%w_force_v)
   end function bt_wide_bytes

   subroutine bt_wide_init(this, grid, dx, dy, lon_west, lat_south, &
                           rad_earth, grid_config, &
                           f_0, beta, y_ref, coriolis_scheme)
      !! Allocate the wide shadow state.  Builds `grid_w` (same nx_phys/ny_phys
      !! as `grid`, nghost = grid%nghost + bt_halo), fills wide metrics via the
      !! same formula generator, fills wide f_corner.
      !! `grid_config` must be GRID_CONFIG_CARTESIAN or GRID_CONFIG_SPHERICAL;
      !! supergrid/tripolar are excluded at configure time.
      class(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
         !! Normal-width grid descriptor for this subdomain.
      real(wp), intent(in) :: dx, dy
         !! Cell spacing (m for Cartesian; deg for spherical).
      real(wp), intent(in) :: lon_west, lat_south
         !! South-west corner (used only for spherical; ignored for Cartesian).
      real(wp), intent(in) :: rad_earth
         !! Earth radius (m; used only for spherical; ignored for Cartesian).
      integer, intent(in) :: grid_config
         !! GRID_CONFIG_CARTESIAN or GRID_CONFIG_SPHERICAL.
      real(wp), intent(in) :: f_0, beta, y_ref
         !! Beta-plane Coriolis parameters.
      integer, intent(in) :: coriolis_scheme
         !! CORIOLIS_SCHEME_BETA_PLANE or CORIOLIS_SCHEME_PLANETARY.

      integer :: nx_w, ny_w, bt_h

      bt_h = this%bt_halo
      this%ng_wide = grid%nghost + bt_h
      this%num_cycles = bt_h/2

      ! Build the wide grid descriptor.  Physical size identical to the normal
      ! grid; nghost increases by bt_halo.  Global offsets must match so the
      ! formula-metric and f_corner fills land on the correct physical rows.
      call this%grid_w%init(grid%nx_phys, grid%ny_phys, this%ng_wide, dx, dy)
      this%grid_w%i_offset_global = grid%i_offset_global
      this%grid_w%j_offset_global = grid%j_offset_global
      this%grid_w%nx_global = grid%nx_global
      this%grid_w%ny_global = grid%ny_global

      nx_w = this%grid_w%nx_total
      ny_w = this%grid_w%ny_total

      ! Build wide metrics via the same formula generator.
      call this%metrics_w%init(this%grid_w)
      select case (grid_config)
      case (GRID_CONFIG_SPHERICAL)
         call metrics_fill_spherical(this%metrics_w, this%grid_w, &
                                     lon_west, lat_south, dx, dy, rad_earth)
      case default   ! GRID_CONFIG_CARTESIAN
         call metrics_fill_cartesian(this%metrics_w, this%grid_w, dx, dy)
      end select
      call metrics_finalize(this%metrics_w)

      ! Fill wide f_corner via metrics_fill_coriolis.
      allocate (this%f_corner_w(nx_w + 1, ny_w + 1), source=0.0_wp)
      block
         real(wp), allocatable :: f_centre_scratch(:, :)
         allocate (f_centre_scratch(nx_w, ny_w), source=0.0_wp)
         call metrics_fill_coriolis(this%metrics_w, coriolis_scheme, &
                                    f_0, beta, y_ref, 0.0_wp, &
                                    this%grid_w, this%f_corner_w, f_centre_scratch)
         deallocate (f_centre_scratch)
      end block

      ! Allocate wide shadow arrays, initialised to zero.
      allocate (this%w_eta(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_H_ref(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_eta_new(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_ke(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_eta_sum(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_eta_end(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_ubt(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_ubt_prev(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_rem_u(nx_w + 1, ny_w), source=1.0_wp)
      allocate (this%w_ubt_sum(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_uhbt_sum(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_uhbt(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_ubt_end(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_vbt(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_vbt_prev(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_rem_v(nx_w, ny_w + 1), source=1.0_wp)
      allocate (this%w_vbt_sum(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_vhbt_sum(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_vhbt(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_vbt_end(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_zeta(nx_w + 1, ny_w + 1), source=0.0_wp)
      allocate (this%w_force_u(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_force_v(nx_w, ny_w + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine bt_wide_init

   subroutine bt_wide_enter_data(this)
      !! Attach all wide arrays (and wide metrics leaf arrays) to the GPU
      !! present table.  The containing `ocean_dyn_t` is already mapped by
      !! the caller; this routine attaches the components.
      class(bt_wide_t), intent(inout) :: this
      select type (this)
      type is (bt_wide_t)
         call bt_wide_enter_data_impl(this)
      end select
   end subroutine bt_wide_enter_data

   subroutine bt_wide_enter_data_impl(this)
      !! Non-polymorphic enter_data body (avoids class-box GPU descriptor issue).
      type(bt_wide_t), intent(inout) :: this
      ! Wide metrics: delegate to the ocean_metrics_t enter_data TBP —
      ! same attach set as the normal metrics, and keeps the acc
      ! directives on one-level names (associate-leaf/ifx rule; the
      ! openmp-portability hook rejects deep struct chains in data
      ! clauses).
      call this%metrics_w%enter_data()
      ! Wide Coriolis.
      !$acc enter data copyin(this%f_corner_w)
      ! Wide shadow arrays.
      !$acc enter data copyin(this%w_eta, this%w_H_ref, this%w_eta_new, this%w_ke)
      !$acc enter data copyin(this%w_eta_sum, this%w_eta_end)
      !$acc enter data copyin(this%w_ubt, this%w_ubt_prev, this%w_rem_u)
      !$acc enter data copyin(this%w_ubt_sum, this%w_uhbt_sum, this%w_uhbt, this%w_ubt_end)
      !$acc enter data copyin(this%w_vbt, this%w_vbt_prev, this%w_rem_v)
      !$acc enter data copyin(this%w_vbt_sum, this%w_vhbt_sum, this%w_vhbt, this%w_vbt_end)
      !$acc enter data copyin(this%w_zeta)
      !$acc enter data copyin(this%w_force_u, this%w_force_v)
   end subroutine bt_wide_enter_data_impl

   subroutine bt_wide_exit_data(this)
      !! Detach all wide arrays from the GPU present table.
      class(bt_wide_t), intent(inout) :: this
      select type (this)
      type is (bt_wide_t)
         call bt_wide_exit_data_impl(this)
      end select
   end subroutine bt_wide_exit_data

   subroutine bt_wide_exit_data_impl(this)
      !! Non-polymorphic exit_data body.
      type(bt_wide_t), intent(inout) :: this
      !$acc exit data delete(this%w_force_u, this%w_force_v)
      !$acc exit data delete(this%w_zeta)
      !$acc exit data delete(this%w_vbt_sum, this%w_vhbt_sum, this%w_vhbt, this%w_vbt_end)
      !$acc exit data delete(this%w_vbt, this%w_vbt_prev, this%w_rem_v)
      !$acc exit data delete(this%w_ubt_sum, this%w_uhbt_sum, this%w_uhbt, this%w_ubt_end)
      !$acc exit data delete(this%w_ubt, this%w_ubt_prev, this%w_rem_u)
      !$acc exit data delete(this%w_eta_sum, this%w_eta_end)
      !$acc exit data delete(this%w_eta, this%w_H_ref, this%w_eta_new, this%w_ke)
      !$acc exit data delete(this%f_corner_w)
      ! Wide metrics: delegate (see enter_data note — one-level names only).
      call this%metrics_w%exit_data()
   end subroutine bt_wide_exit_data_impl

   subroutine bt_wide_copy_in(this, grid, &
                              bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                              bt_ubt_prev, bt_vbt_prev, &
                              bt_rem_u, bt_rem_v, &
                              force_u, force_v)
      !! Offset-copy normal-width input arrays into the wide shadow arrays.
      !! Dispatches to the non-polymorphic `_impl` body to avoid the
      !! class-box GPU descriptor issue (same pattern as enter/exit_data).
      class(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: bt_eta(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: bt_H_ref(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: bt_ubt(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: bt_vbt(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: bt_ubt_prev(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: bt_vbt_prev(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: bt_rem_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: bt_rem_v(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: force_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: force_v(grid%nx_total, grid%ny_total + 1)
      select type (this)
      type is (bt_wide_t)
         call bt_wide_copy_in_impl(this, grid, &
                                   bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                                   bt_ubt_prev, bt_vbt_prev, &
                                   bt_rem_u, bt_rem_v, force_u, force_v)
      end select
   end subroutine bt_wide_copy_in

   subroutine bt_wide_copy_in_impl(this, grid, &
                                   bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                                   bt_ubt_prev, bt_vbt_prev, &
                                   bt_rem_u, bt_rem_v, &
                                   force_u, force_v)
      !! Non-polymorphic copy_in body.  Offset = bt_halo:
      !! w_X(iw, jw) = X(clamp(iw-off), clamp(jw-off)) over the FULL wide
      !! extent — the inner band is a direct offset copy; the outer bt_halo
      !! ring is a clamped-index (constant-extrapolation) fill.  The ring
      !! fill matters: without it the ring carries stale end-of-fast-loop
      !! values from the previous stage (H_ref = 0, eta from t-1), which at
      !! a PHYSICAL (non-seam) edge is never refreshed by any exchange and
      !! free-runs an inconsistent zero-depth integration that blows up in
      !! O(25) outer steps.  At an MPI seam the ring is immediately
      !! overwritten with true neighbour data by `entry_exchange`, so the
      !! clamped fill only governs physical edges — the same sane ghost-band
      !! construction the v1 normal-width path gets from its own ghosts.
      !! Scratch / accumulator arrays (w_eta_new, w_ke, w_eta_sum, …) do not
      !! need copy-in — the substep initialises them.
      type(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
         !! Normal-width grid (provides nx_total, ny_total for loop bounds).
      real(wp), intent(in) :: bt_eta(grid%nx_total, grid%ny_total)
         !! Barotropic SSH (cell centres, normal-width).
      real(wp), intent(in) :: bt_H_ref(grid%nx_total, grid%ny_total)
         !! Reference column depth (cell centres, normal-width).
      real(wp), intent(in) :: bt_ubt(grid%nx_total + 1, grid%ny_total)
         !! BT u (east faces, normal-width).
      real(wp), intent(in) :: bt_vbt(grid%nx_total, grid%ny_total + 1)
         !! BT v (north faces, normal-width).
      real(wp), intent(in) :: bt_ubt_prev(grid%nx_total + 1, grid%ny_total)
         !! BEBT u^{n-1} snapshot (east faces, normal-width).
      real(wp), intent(in) :: bt_vbt_prev(grid%nx_total, grid%ny_total + 1)
         !! BEBT v^{n-1} snapshot (north faces, normal-width).
      real(wp), intent(in) :: bt_rem_u(grid%nx_total + 1, grid%ny_total)
         !! Multiplicative drag factor for u (east faces, normal-width).
      real(wp), intent(in) :: bt_rem_v(grid%nx_total, grid%ny_total + 1)
         !! Multiplicative drag factor for v (north faces, normal-width).
      real(wp), intent(in) :: force_u(grid%nx_total + 1, grid%ny_total)
         !! BT slow forcing for u (east faces, normal-width).
      real(wp), intent(in) :: force_v(grid%nx_total, grid%ny_total + 1)
         !! BT slow forcing for v (north faces, normal-width).

      integer :: i, j, nx, ny, off
      integer :: i_src, j_src
      off = this%bt_halo
      nx = grid%nx_total
      ny = grid%ny_total

      ! Centre arrays: full wide extent (nx + 2*off, ny + 2*off), clamped
      ! source indices (constant extrapolation into the outer ring).
      do concurrent(j=1:ny + 2*off, i=1:nx + 2*off) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx)
         j_src = min(max(j - off, 1), ny)
         this%w_eta(i, j) = bt_eta(i_src, j_src)
         this%w_H_ref(i, j) = bt_H_ref(i_src, j_src)
      end do
      ! East-face arrays: wide extent (nx + 2*off + 1, ny + 2*off).
      do concurrent(j=1:ny + 2*off, i=1:nx + 2*off + 1) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx + 1)
         j_src = min(max(j - off, 1), ny)
         this%w_ubt(i, j) = bt_ubt(i_src, j_src)
         this%w_ubt_prev(i, j) = bt_ubt_prev(i_src, j_src)
         this%w_rem_u(i, j) = bt_rem_u(i_src, j_src)
         this%w_force_u(i, j) = force_u(i_src, j_src)
      end do
      ! North-face arrays: wide extent (nx + 2*off, ny + 2*off + 1).
      do concurrent(j=1:ny + 2*off + 1, i=1:nx + 2*off) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx)
         j_src = min(max(j - off, 1), ny + 1)
         this%w_vbt(i, j) = bt_vbt(i_src, j_src)
         this%w_vbt_prev(i, j) = bt_vbt_prev(i_src, j_src)
         this%w_rem_v(i, j) = bt_rem_v(i_src, j_src)
         this%w_force_v(i, j) = force_v(i_src, j_src)
      end do
   end subroutine bt_wide_copy_in_impl

   subroutine bt_wide_entry_exchange(this)
      !! One wide grouped exchange (eta+ubt+vbt) + wide singles for the other
      !! 7 input arrays (H_ref, ubt_prev, rem_u, force_u, vbt_prev, rem_v,
      !! force_v).  Fills the entire wide ghost band before the fast loop.
      !! Counter effect: +1 bt_group, +1 centre_2d, +3 face_x_2d, +3 face_y_2d.
      !! Dispatches to the non-polymorphic `_impl` body.
      class(bt_wide_t), intent(inout) :: this
      select type (this)
      type is (bt_wide_t)
         call bt_wide_entry_exchange_impl(this)
      end select
   end subroutine bt_wide_entry_exchange

   subroutine bt_wide_entry_exchange_impl(this)
      !! Non-polymorphic entry_exchange body.
      type(bt_wide_t), intent(inout) :: this
      integer :: ng_w
      ng_w = this%ng_wide
      call profiler_start("ocean_comms_bt")
      call ocean_halo_bt_group_2d_wide(this%w_eta, this%w_ubt, this%w_vbt, ng_w)
      call ocean_halo_centre_2d_wide(this%w_H_ref, ng_w)
      call ocean_halo_face_x_2d_wide(this%w_ubt_prev, ng_w)
      call ocean_halo_face_x_2d_wide(this%w_rem_u, ng_w)
      call ocean_halo_face_x_2d_wide(this%w_force_u, ng_w)
      call ocean_halo_face_y_2d_wide(this%w_vbt_prev, ng_w)
      call ocean_halo_face_y_2d_wide(this%w_rem_v, ng_w)
      call ocean_halo_face_y_2d_wide(this%w_force_v, ng_w)
      call profiler_stop("ocean_comms_bt")
   end subroutine bt_wide_entry_exchange_impl

   subroutine bt_wide_copy_out(this, grid, &
                               bt_eta, bt_ubt, bt_vbt, &
                               bt_uhbt, bt_vhbt, &
                               bt_eta_end, bt_ubt_end, bt_vbt_end)
      !! Offset-copy wide output arrays back to the normal-width arrays.
      !! Dispatches to the non-polymorphic `_impl` body.
      class(bt_wide_t), intent(in) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: bt_eta(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: bt_ubt(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(out) :: bt_vbt(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(out) :: bt_uhbt(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(out) :: bt_vhbt(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(out) :: bt_eta_end(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: bt_ubt_end(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(out) :: bt_vbt_end(grid%nx_total, grid%ny_total + 1)
      select type (this)
      type is (bt_wide_t)
         call bt_wide_copy_out_impl(this, grid, &
                                    bt_eta, bt_ubt, bt_vbt, &
                                    bt_uhbt, bt_vhbt, &
                                    bt_eta_end, bt_ubt_end, bt_vbt_end)
      end select
   end subroutine bt_wide_copy_out

   subroutine bt_wide_copy_out_impl(this, grid, &
                                    bt_eta, bt_ubt, bt_vbt, &
                                    bt_uhbt, bt_vhbt, &
                                    bt_eta_end, bt_ubt_end, bt_vbt_end)
      !! Non-polymorphic copy_out body.  X(i,j) = w_X(i+off, j+off)
      !! over the full normal index range.  Outputs: time-mean eta/ubt/vbt,
      !! uhbt/vhbt, *_end snapshots.
      type(bt_wide_t), intent(in) :: this
      type(hgrid_t), intent(in) :: grid
         !! Normal-width grid.
      real(wp), intent(out) :: bt_eta(grid%nx_total, grid%ny_total)
         !! Time-mean barotropic SSH (output).
      real(wp), intent(out) :: bt_ubt(grid%nx_total + 1, grid%ny_total)
         !! Time-mean BT u (output).
      real(wp), intent(out) :: bt_vbt(grid%nx_total, grid%ny_total + 1)
         !! Time-mean BT v (output).
      real(wp), intent(out) :: bt_uhbt(grid%nx_total + 1, grid%ny_total)
         !! Time-mean depth-integrated u transport (output).
      real(wp), intent(out) :: bt_vhbt(grid%nx_total, grid%ny_total + 1)
         !! Time-mean depth-integrated v transport (output).
      real(wp), intent(out) :: bt_eta_end(grid%nx_total, grid%ny_total)
         !! End-of-loop eta snapshot (output).
      real(wp), intent(out) :: bt_ubt_end(grid%nx_total + 1, grid%ny_total)
         !! End-of-loop u snapshot (output).
      real(wp), intent(out) :: bt_vbt_end(grid%nx_total, grid%ny_total + 1)
         !! End-of-loop v snapshot (output).

      integer :: i, j, nx, ny, off
      off = this%bt_halo
      nx = grid%nx_total
      ny = grid%ny_total

      do concurrent(j=1:ny, i=1:nx)
         bt_eta(i, j) = this%w_eta(i + off, j + off)
         bt_eta_end(i, j) = this%w_eta_end(i + off, j + off)
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_ubt(i, j) = this%w_ubt(i + off, j + off)
         bt_uhbt(i, j) = this%w_uhbt(i + off, j + off)
         bt_ubt_end(i, j) = this%w_ubt_end(i + off, j + off)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_vbt(i, j) = this%w_vbt(i + off, j + off)
         bt_vhbt(i, j) = this%w_vhbt(i + off, j + off)
         bt_vbt_end(i, j) = this%w_vbt_end(i + off, j + off)
      end do
   end subroutine bt_wide_copy_out_impl

   subroutine bt_wide_destroy(this)
      !! Deallocate all wide state.
      class(bt_wide_t), intent(inout) :: this
      this%is_init = .false.
      call this%metrics_w%destroy()
      if (allocated(this%f_corner_w)) deallocate (this%f_corner_w)
      if (allocated(this%w_eta)) deallocate (this%w_eta)
      if (allocated(this%w_H_ref)) deallocate (this%w_H_ref)
      if (allocated(this%w_eta_new)) deallocate (this%w_eta_new)
      if (allocated(this%w_ke)) deallocate (this%w_ke)
      if (allocated(this%w_eta_sum)) deallocate (this%w_eta_sum)
      if (allocated(this%w_eta_end)) deallocate (this%w_eta_end)
      if (allocated(this%w_ubt)) deallocate (this%w_ubt)
      if (allocated(this%w_ubt_prev)) deallocate (this%w_ubt_prev)
      if (allocated(this%w_rem_u)) deallocate (this%w_rem_u)
      if (allocated(this%w_ubt_sum)) deallocate (this%w_ubt_sum)
      if (allocated(this%w_uhbt_sum)) deallocate (this%w_uhbt_sum)
      if (allocated(this%w_uhbt)) deallocate (this%w_uhbt)
      if (allocated(this%w_ubt_end)) deallocate (this%w_ubt_end)
      if (allocated(this%w_vbt)) deallocate (this%w_vbt)
      if (allocated(this%w_vbt_prev)) deallocate (this%w_vbt_prev)
      if (allocated(this%w_rem_v)) deallocate (this%w_rem_v)
      if (allocated(this%w_vbt_sum)) deallocate (this%w_vbt_sum)
      if (allocated(this%w_vhbt_sum)) deallocate (this%w_vhbt_sum)
      if (allocated(this%w_vhbt)) deallocate (this%w_vhbt)
      if (allocated(this%w_vbt_end)) deallocate (this%w_vbt_end)
      if (allocated(this%w_zeta)) deallocate (this%w_zeta)
      if (allocated(this%w_force_u)) deallocate (this%w_force_u)
      if (allocated(this%w_force_v)) deallocate (this%w_force_v)
   end subroutine bt_wide_destroy

   subroutine bt_wide_substep(bt_wide, bt_work, n_steps, dt_inner, bc)
      !! Wide-halo (march-in) entry point for the nonlinear barotropic fast
      !! loop.  Unpacks the wide shadow arrays (`w_*`, wide grid/metrics, wide
      !! `f_corner`) and forwards them to `barotropic_substep_nonlinear` with
      !! `bt_halo = bt_wide%bt_halo`, so the ~20-array plumbing lives here once
      !! rather than at the call site.  `bc` propagates by absence.  Tides
      !! (`eta_forcing`) are a configure-time exclusion on the wide path, so
      !! none is forwarded.  Interior twin: `barotropic_substep_nonlinear_interior`.
      type(bt_wide_t), intent(inout) :: bt_wide
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      type(ocean_bc_state_t), intent(inout), optional :: bc

      call barotropic_substep_nonlinear(bt_wide%grid_w, &
                                        bt_work, &
                                        bt_wide%w_force_u, bt_wide%w_force_v, &
                                        n_steps, dt_inner, &
                                        bt_eta=bt_wide%w_eta, bt_H_ref=bt_wide%w_H_ref, &
                                        bt_eta_new=bt_wide%w_eta_new, &
                                        bt_ke_centre=bt_wide%w_ke, &
                                        eta_sum=bt_wide%w_eta_sum, &
                                        bt_eta_end=bt_wide%w_eta_end, &
                                        bt_ubt=bt_wide%w_ubt, &
                                        bt_ubt_prev=bt_wide%w_ubt_prev, &
                                        bt_rem_u=bt_wide%w_rem_u, &
                                        ubt_sum=bt_wide%w_ubt_sum, &
                                        uhbt_sum=bt_wide%w_uhbt_sum, &
                                        bt_uhbt=bt_wide%w_uhbt, &
                                        bt_ubt_end=bt_wide%w_ubt_end, &
                                        bt_vbt=bt_wide%w_vbt, &
                                        bt_vbt_prev=bt_wide%w_vbt_prev, &
                                        bt_rem_v=bt_wide%w_rem_v, &
                                        vbt_sum=bt_wide%w_vbt_sum, &
                                        vhbt_sum=bt_wide%w_vhbt_sum, &
                                        bt_vhbt=bt_wide%w_vhbt, &
                                        bt_vbt_end=bt_wide%w_vbt_end, &
                                        bt_zeta_corner=bt_wide%w_zeta, &
                                        f_corner=bt_wide%f_corner_w, &
                                        area_cu=bt_wide%metrics_w%areaCu, area_cv=bt_wide%metrics_w%areaCv, &
                                        dx_cu=bt_wide%metrics_w%dxCu, dx_cv=bt_wide%metrics_w%dx_cv, &
                                        dy_cu=bt_wide%metrics_w%dy_cu, dy_cv=bt_wide%metrics_w%dyCv, &
                                        iarea_bu=bt_wide%metrics_w%iareaBu, iarea_t=bt_wide%metrics_w%iareaT, &
                                        idx_cu=bt_wide%metrics_w%idxCu, idy_cv=bt_wide%metrics_w%idyCv, &
                                        bc=bc, bt_halo=bt_wide%bt_halo)
   end subroutine bt_wide_substep

end module rdb_ocean_bt_wide
