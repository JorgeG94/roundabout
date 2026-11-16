!! Wide-halo BT march-in state (Phase 3c, D3 v1.1).
module rdb_ocean_bt_wide
   !! Shadow state for the wide-halo barotropic march-in.
   !! When `bt_halo > 0`, the BT fast loop runs on WIDE arrays with ghost
   !! width `ng_wide = nghost + bt_halo`.  Ghost cells outside the valid
   !! band evolve stale data that creeps INWARD at one cell per substep; one
   !! PACKED exchange (a single message per neighbour per direction pass for
   !! eta, ubt, vbt and, with `bebt > 0`, `ubt_prev`/`vbt_prev`) every
   !! `bt_wide_cycle_len(ng_wide) = ng_wide - 2` substeps keeps the owned cells
   !! exact (proven limit `ng_wide - 1`; one substep of margin).  The per-solve
   !! entry set (eta, ubt, vbt, rem, force) is one packed exchange too, the
   !! static `H_ref` is exchanged once, and the normal-width time means,
   !! end-of-loop snapshots and transports leave in one packed exit exchange.
   !! The mid-substep u exchange is absorbed.  Proof: owned cells are bitwise
   !! equal to the one-rank run with every invalidated layer poisoned between
   !! exchanges (`RDB_BT_INHALO_POISON`, `test_ocean_dyn_mpi`), and the cadence
   !! one substep past the limit breaks them.
   !!
   !! Lifecycle:
   !!   1. `init`         — allocate wide arrays; fill the STATIC wide metrics and
   !!                      `f_corner` from the model's own arrays (host, once).
   !!   2. `enter_data`   — attach all wide arrays to the GPU present table.
   !!   3. Per outer step: `copy_in` -> (entry_exchange + substep) -> `copy_out`
   !!                      -> normal-width exit exchange (in the caller).
   !!   4. `exit_data`    — release GPU present table entries.
   !!   5. `destroy`      — free host memory.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_halo, only: ocean_halo_packed_wide, &
                             ocean_halo_centre_2d_wide, &
                             ocean_halo_face_x_2d_wide, &
                             ocean_halo_face_y_2d_wide
   use rdb_profiler, only: profiler_start, profiler_stop
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear, bt_wide_cycle_len
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_mem_report, only: arr_bytes
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: bt_wide_t
   public :: bt_wide_substep

   type :: bt_wide_t
      !! Wide-halo shadow state for the barotropic fast loop.
      logical :: is_init = .false.
         !! True after `init`, before `destroy`.
      logical :: h_ref_ready = .false.
         !! `w_H_ref` has been copied in and wide-exchanged (once, first solve).
      integer :: bt_halo = 0
         !! Requested wide-halo width (cells; even, > 0).
      integer :: ng_wide = 0
         !! Effective ghost width: `nghost + bt_halo`.
      integer :: num_cycles = 0
         !! Substeps between in-loop wide exchanges (`bt_wide_cycle_len`).

      ! Wide grid.  Same `nx_phys/ny_phys` as the normal grid; only `nghost`
      ! differs (`ng_wide`).
      type(hgrid_t) :: grid_w
         !! Wide hgrid_t: nghost = ng_wide.

      ! STATIC wide fields the fast loop reads.  They come from their OWNERS
      ! (the model's own, masked, possibly file-based, possibly closed-face
      ! narrowed arrays) by an offset copy plus one host-side wide exchange at
      ! setup (`bt_wide_init`) -- never from a grid formula, so land, the
      ! wall-ghost mask and file metrics are the numbers every rank and the
      ! one-rank `bt_halo = 0` run see.  Mapped once, never refilled
      ! (`w_dy_cu`/`w_dx_cv` excepted under porous barriers, per solve).
      real(wp), allocatable :: w_area_cu(:, :)
         !! `areaCu` (nx_w+1, ny_w).
      real(wp), allocatable :: w_area_cv(:, :)
         !! `areaCv` (nx_w, ny_w+1).
      real(wp), allocatable :: w_dx_cu(:, :)
         !! `dxCu` (nx_w+1, ny_w).
      real(wp), allocatable :: w_dx_cv(:, :)
         !! BT width `dx_cv_bt` (nx_w, ny_w+1).
      real(wp), allocatable :: w_dy_cu(:, :)
         !! BT width `dy_cu_bt` (nx_w+1, ny_w).
      real(wp), allocatable :: w_dy_cv(:, :)
         !! `dyCv` (nx_w, ny_w+1).
      real(wp), allocatable :: w_iarea_bu(:, :)
         !! `iareaBu` (nx_w+1, ny_w+1).
      real(wp), allocatable :: w_iarea_t(:, :)
         !! `iareaT` (nx_w, ny_w).
      real(wp), allocatable :: w_idx_cu(:, :)
         !! `idxCu` (nx_w+1, ny_w).
      real(wp), allocatable :: w_idy_cv(:, :)
         !! `idyCv` (nx_w, ny_w+1).
      real(wp), allocatable :: f_corner_w(:, :)
         !! Coriolis at the wide corners (nx_w+1, ny_w+1).
      logical :: refresh_widths = .false.
         !! Porous barriers or z-level closed faces on: `w_dy_cu`/`w_dx_cv`
         !! are re-copied from the model's `dy_cu_bt`/`dx_cv_bt` and
         !! re-exchanged once per solve (the model recomputes them once per
         !! outer step from the live `h`).
      integer :: ring_poison = 0
         !! Test hook (`RDB_BT_RING_POISON`, set at setup): units digit 1 =
         !! finite garbage, 2 = NaN; tens digit = EXTRA layers poisoned beyond
         !! the ring, into the normal-width ghost band the kernel does read
         !! (a liveness control: the owned result must then change).  The
         !! outer `bt_halo` ring of EVERY wide array (static and
         !! per-solve) is overwritten with it before the wide exchange, so
         !! that at a physical edge -- where no exchange refills it -- the ring
         !! stays garbage for the whole run.  Owned results must not change:
         !! the substep's insets never read it.  Seam / periodic rings are
         !! refilled by the exchange, so the leg is vacuous there by design.
      real(wp), allocatable :: w_eta_forcing(:, :)
         !! Equilibrium-tide / pressure-loading elevation (nx_w, ny_w); copied
         !! in and exchanged once per solve, only when `eta_forcing` is passed.

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
      procedure, non_overridable :: set_eta_forcing => bt_wide_set_eta_forcing
      procedure, non_overridable :: copy_out => bt_wide_copy_out
      procedure, non_overridable :: destroy => bt_wide_destroy
      procedure, non_overridable :: bytes => bt_wide_bytes
   end type bt_wide_t

contains

   pure function bt_wide_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the wide-halo BT shadow state
      !! (0 when unallocated, i.e. whenever `&ocean_bt_nml bt_halo = 0`).
      !!
      !! Every term here is device-mapped by `bt_wide_enter_data_impl`.  The
      !! static metrics are the 10 arrays the fast loop reads (plus
      !! `f_corner`), not a whole `ocean_metrics_t` clone.  One `arr_bytes`
      !! term per array: add one here when a new allocatable joins the type.
      class(bt_wide_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_corner_w) &
               + arr_bytes(this%w_area_cu) &
               + arr_bytes(this%w_area_cv) &
               + arr_bytes(this%w_dx_cu) &
               + arr_bytes(this%w_dx_cv) &
               + arr_bytes(this%w_dy_cu) &
               + arr_bytes(this%w_dy_cv) &
               + arr_bytes(this%w_iarea_bu) &
               + arr_bytes(this%w_iarea_t) &
               + arr_bytes(this%w_idx_cu) &
               + arr_bytes(this%w_idy_cv) &
               + arr_bytes(this%w_eta_forcing) &
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

   subroutine bt_wide_init(this, grid, metrics, f_corner, refresh_widths, ring_poison, &
                           eta_forcing_on)
      !! Allocate the wide shadow state and fill its STATIC fields from their
      !! owners.  Builds `grid_w` (same nx_phys/ny_phys as `grid`,
      !! nghost = grid%nghost + bt_halo).  Every metric array the fast loop
      !! reads, and `f_corner`, is the model's own array offset-copied into
      !! the wide extent (the outer `bt_halo` ring is a clamped copy, finite
      !! and never read at a physical edge: the substep's insets re-anchor
      !! every array-edge closure at the NORMAL array edge), then ONE
      !! host-side wide exchange gives every seam / periodic ghost its owner's
      !! value.  Land masks, wall-ghost zeros, file-based (supergrid) metrics
      !! and closed-face widths therefore arrive exactly as the one-rank
      !! `bt_halo = 0` kernel sees them.
      !!
      !! Must run after `ocean_halo_init` + `ocean_halo_reserve(ng_wide)` and
      !! BEFORE `enter_data` (host fill: the device copy is the copyin of the
      !! finished host arrays), with the model's metrics and `f_corner` final
      !! on the host.
      class(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
         !! Normal-width grid descriptor for this subdomain.
      type(ocean_metrics_t), intent(in) :: metrics
         !! The model's metrics (host copy, final).
      real(wp), intent(in) :: f_corner(:, :)
         !! The model's Coriolis at corners (host copy, final).
      logical, intent(in), optional :: refresh_widths
         !! Porous barriers or closed faces on (per-solve width refresh).
      integer, intent(in), optional :: ring_poison
         !! Test hook: see `bt_wide_t%ring_poison`.
      logical, intent(in), optional :: eta_forcing_on
         !! A body tide or surface-pressure load will be passed
         !! (`set_eta_forcing`): allocate the wide `eta_forcing`.  Absent /
         !! false => a 1 x 1 placeholder (still mapped, never read).

      integer :: nx_w, ny_w, bt_h

      bt_h = this%bt_halo
      this%ng_wide = grid%nghost + bt_h
      this%num_cycles = bt_wide_cycle_len(this%ng_wide)
      this%refresh_widths = .false.
      if (present(refresh_widths)) this%refresh_widths = refresh_widths
      this%ring_poison = 0
      if (present(ring_poison)) this%ring_poison = ring_poison

      ! Build the wide grid descriptor.  Physical size identical to the normal
      ! grid; nghost increases by bt_halo.
      call this%grid_w%init(grid%nx_phys, grid%ny_phys, this%ng_wide, grid%dx, grid%dy)
      this%grid_w%i_offset_global = grid%i_offset_global
      this%grid_w%j_offset_global = grid%j_offset_global
      this%grid_w%nx_global = grid%nx_global
      this%grid_w%ny_global = grid%ny_global

      nx_w = this%grid_w%nx_total
      ny_w = this%grid_w%ny_total

      allocate (this%w_area_cu(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_dx_cu(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_dy_cu(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_idx_cu(nx_w + 1, ny_w), source=0.0_wp)
      allocate (this%w_area_cv(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_dx_cv(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_dy_cv(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_idy_cv(nx_w, ny_w + 1), source=0.0_wp)
      allocate (this%w_iarea_t(nx_w, ny_w), source=0.0_wp)
      allocate (this%w_iarea_bu(nx_w + 1, ny_w + 1), source=0.0_wp)
      allocate (this%f_corner_w(nx_w + 1, ny_w + 1), source=0.0_wp)
      if (present(eta_forcing_on)) then
         if (eta_forcing_on) allocate (this%w_eta_forcing(nx_w, ny_w), source=0.0_wp)
      end if
      if (.not. allocated(this%w_eta_forcing)) allocate (this%w_eta_forcing(1, 1), source=0.0_wp)

      call offset_fill(this%w_area_cu, metrics%areaCu, bt_h)
      call offset_fill(this%w_dx_cu, metrics%dxCu, bt_h)
      call offset_fill(this%w_dy_cu, metrics%dy_cu_bt, bt_h)
      call offset_fill(this%w_idx_cu, metrics%idxCu, bt_h)
      call offset_fill(this%w_area_cv, metrics%areaCv, bt_h)
      call offset_fill(this%w_dx_cv, metrics%dx_cv_bt, bt_h)
      call offset_fill(this%w_dy_cv, metrics%dyCv, bt_h)
      call offset_fill(this%w_idy_cv, metrics%idyCv, bt_h)
      call offset_fill(this%w_iarea_t, metrics%iareaT, bt_h)
      call offset_fill(this%w_iarea_bu, metrics%iareaBu, bt_h)
      call offset_fill(this%f_corner_w, f_corner, bt_h)

      if (this%ring_poison > 0) then
         call poison_ring_host(this%w_area_cu, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_dx_cu, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_dy_cu, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_idx_cu, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_area_cv, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_dx_cv, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_dy_cv, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_idy_cv, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_iarea_t, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%w_iarea_bu, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
         call poison_ring_host(this%f_corner_w, bt_h + this%ring_poison/10, mod(this%ring_poison, 10))
      end if

      ! ONE host-side wide exchange of the static set: seam / periodic ghosts
      ! get their owner's value.  (Single-rank non-periodic: a no-op, the
      ! clamped fill stands.)  Corner fields have no primitive: see
      ! `exchange_corner_host`.
      call ocean_halo_face_x_2d_wide(this%w_area_cu, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_x_2d_wide(this%w_dx_cu, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_x_2d_wide(this%w_dy_cu, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_x_2d_wide(this%w_idx_cu, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_y_2d_wide(this%w_area_cv, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_y_2d_wide(this%w_dx_cv, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_y_2d_wide(this%w_dy_cv, this%ng_wide, device_resident=.false.)
      call ocean_halo_face_y_2d_wide(this%w_idy_cv, this%ng_wide, device_resident=.false.)
      call ocean_halo_centre_2d_wide(this%w_iarea_t, this%ng_wide, device_resident=.false.)
      call exchange_corner_host(this%w_iarea_bu, this%ng_wide)
      call exchange_corner_host(this%f_corner_w, this%ng_wide)

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

   pure function ring_value(mode, i, j) result(v)
      !! Garbage for the ring-poison hook: finite (mode 1) or NaN (mode 2).
      !$acc routine seq
      integer, intent(in) :: mode, i, j
      real(wp) :: v
      if (mode == 2) then
         v = ieee_value(1.0_wp, ieee_quiet_nan)
      else
         v = 1.0e3_wp*real(1 + mod(abs(7*i + 13*j), 11), wp)
      end if
   end function ring_value

   subroutine poison_ring_host(a, off, mode)
      !! Overwrite the outer `off` ring of `a` (any stagger) with garbage.  Host.
      real(wp), intent(inout) :: a(:, :)
      integer, intent(in) :: off, mode
      integer :: i, j, n1, n2
      n1 = size(a, 1)
      n2 = size(a, 2)
      do j = 1, n2
         do i = 1, n1
            if (i <= off .or. i > n1 - off .or. j <= off .or. j > n2 - off) then
               a(i, j) = ring_value(mode, i, j)
            end if
         end do
      end do
   end subroutine poison_ring_host

   subroutine poison_ring_dev(a, n1, n2, off, mode)
      !! Device twin of `poison_ring_host` (explicit shape, `do concurrent`).
      integer, intent(in) :: n1, n2, off, mode
      real(wp), intent(inout) :: a(n1, n2)
      integer :: i, j
      do concurrent(j=1:n2, i=1:n1)
         if (i <= off .or. i > n1 - off .or. j <= off .or. j > n2 - off) then
            a(i, j) = ring_value(mode, i, j)
         end if
      end do
   end subroutine poison_ring_dev

   subroutine offset_fill(w, a, off)
      !! `w(i,j) = a(clamp(i-off), clamp(j-off))` over the whole of `w`: the
      !! inner band is a direct offset copy of the owner's array, the outer
      !! `off` ring a clamped copy (finite; overwritten by the wide exchange
      !! at a seam and never read at a physical edge).  Host, setup only.
      real(wp), intent(inout) :: w(:, :)
      real(wp), intent(in) :: a(:, :)
      integer, intent(in) :: off
      integer :: i, j
      if (size(w, 1) /= size(a, 1) + 2*off .or. size(w, 2) /= size(a, 2) + 2*off) then
         error stop "bt_wide offset_fill: wide/normal extents disagree"
      end if
      do j = 1, size(w, 2)
         do i = 1, size(w, 1)
            w(i, j) = a(min(max(i - off, 1), size(a, 1)), min(max(j - off, 1), size(a, 2)))
         end do
      end do
   end subroutine offset_fill

   subroutine exchange_corner_host(c, ng)
      !! Wide host exchange of a CORNER-staggered field `c(nxt+1, nyt+1)`
      !! (no corner primitive exists in the halo module).  A corner array is
      !! x-face-like in x and, in y, differs from a cell-centred row index
      !! only by one extra row: its rows map to global rows exactly as
      !! cell rows do.  So the x-face exchange applied to rows 1:nyt, then to
      !! rows 2:nyt+1 (as two `(nxt+1, nyt)` temporaries), gives every row of
      !! `c` its owner's value (the second pass supplies the last row), with
      !! the seam column duplicated as for a face and the two-pass corner
      !! rule, periodic wrap included.  Host, setup only.
      real(wp), intent(inout) :: c(:, :)
      integer, intent(in) :: ng
      real(wp), allocatable :: t(:, :)
      integer :: ny1
      ny1 = size(c, 2)
      allocate (t(size(c, 1), ny1 - 1))
      t = c(:, 1:ny1 - 1)
      call ocean_halo_face_x_2d_wide(t, ng, device_resident=.false.)
      c(:, 1:ny1 - 1) = t
      t = c(:, 2:ny1)
      call ocean_halo_face_x_2d_wide(t, ng, device_resident=.false.)
      c(:, 2:ny1) = t
      deallocate (t)
   end subroutine exchange_corner_host

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
      ! Static wide metrics + Coriolis (host-filled from the owners by
      ! `bt_wide_init`: this copyin is what puts the owners' values on the
      ! device).
      !$acc enter data copyin(this%w_area_cu, this%w_area_cv, this%w_dx_cu, this%w_dx_cv)
      !$acc enter data copyin(this%w_dy_cu, this%w_dy_cv, this%w_iarea_bu, this%w_iarea_t)
      !$acc enter data copyin(this%w_idx_cu, this%w_idy_cv, this%f_corner_w)
      !$acc enter data copyin(this%w_eta_forcing)
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
      !$acc exit data delete(this%w_eta_forcing)
      !$acc exit data delete(this%w_idx_cu, this%w_idy_cv, this%f_corner_w)
      !$acc exit data delete(this%w_dy_cu, this%w_dy_cv, this%w_iarea_bu, this%w_iarea_t)
      !$acc exit data delete(this%w_area_cu, this%w_area_cv, this%w_dx_cu, this%w_dx_cv)
   end subroutine bt_wide_exit_data_impl

   subroutine bt_wide_copy_in(this, grid, &
                              bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                              bt_ubt_prev, bt_vbt_prev, &
                              bt_rem_u, bt_rem_v, &
                              force_u, force_v, dy_cu_bt, dx_cv_bt)
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
      real(wp), intent(in) :: dy_cu_bt(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: dx_cv_bt(grid%nx_total, grid%ny_total + 1)
      select type (this)
      type is (bt_wide_t)
         call bt_wide_copy_in_impl(this, grid, &
                                   bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                                   bt_ubt_prev, bt_vbt_prev, &
                                   bt_rem_u, bt_rem_v, force_u, force_v, &
                                   dy_cu_bt, dx_cv_bt)
      end select
   end subroutine bt_wide_copy_in

   subroutine bt_wide_copy_in_impl(this, grid, &
                                   bt_eta, bt_H_ref, bt_ubt, bt_vbt, &
                                   bt_ubt_prev, bt_vbt_prev, &
                                   bt_rem_u, bt_rem_v, &
                                   force_u, force_v, dy_cu_bt, dx_cv_bt)
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
      real(wp), intent(in) :: dy_cu_bt(grid%nx_total + 1, grid%ny_total)
         !! BT u-face width (normal-width); re-copied only under porous barriers.
      real(wp), intent(in) :: dx_cv_bt(grid%nx_total, grid%ny_total + 1)
         !! BT v-face width (normal-width); re-copied only under porous barriers.

      integer :: i, j, nx, ny, off
      logical :: refresh_widths
      integer :: i_src, j_src
      off = this%bt_halo
      nx = grid%nx_total
      ny = grid%ny_total
      refresh_widths = this%refresh_widths

      ! Centre arrays: full wide extent (nx + 2*off, ny + 2*off), clamped
      ! source indices (constant extrapolation into the outer ring).
      do concurrent(j=1:ny + 2*off, i=1:nx + 2*off) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx)
         j_src = min(max(j - off, 1), ny)
         this%w_eta(i, j) = bt_eta(i_src, j_src)
      end do
      ! `H_ref` is static (the datum `b`, set at configure; the cavity, which
      ! would move it, is refused with bt_halo > 0): copied and exchanged ONCE,
      ! by the first `entry_exchange`, and left alone afterwards.
      if (.not. this%h_ref_ready) then
         do concurrent(j=1:ny + 2*off, i=1:nx + 2*off) local(i_src, j_src)
            i_src = min(max(i - off, 1), nx)
            j_src = min(max(j - off, 1), ny)
            this%w_H_ref(i, j) = bt_H_ref(i_src, j_src)
         end do
      end if
      ! East-face arrays: wide extent (nx + 2*off + 1, ny + 2*off).
      do concurrent(j=1:ny + 2*off, i=1:nx + 2*off + 1) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx + 1)
         j_src = min(max(j - off, 1), ny)
         this%w_ubt(i, j) = bt_ubt(i_src, j_src)
         this%w_ubt_prev(i, j) = bt_ubt_prev(i_src, j_src)
         this%w_rem_u(i, j) = bt_rem_u(i_src, j_src)
         this%w_force_u(i, j) = force_u(i_src, j_src)
      end do
      ! Porous barriers: the BT widths are recomputed once per outer step on
      ! the device; this wide copy is re-exchanged in `entry_exchange`.
      if (refresh_widths) then
         do concurrent(j=1:ny + 2*off, i=1:nx + 2*off + 1) local(i_src, j_src)
            i_src = min(max(i - off, 1), nx + 1)
            j_src = min(max(j - off, 1), ny)
            this%w_dy_cu(i, j) = dy_cu_bt(i_src, j_src)
         end do
      end if
      ! North-face arrays: wide extent (nx + 2*off, ny + 2*off + 1).
      do concurrent(j=1:ny + 2*off + 1, i=1:nx + 2*off) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx)
         j_src = min(max(j - off, 1), ny + 1)
         this%w_vbt(i, j) = bt_vbt(i_src, j_src)
         this%w_vbt_prev(i, j) = bt_vbt_prev(i_src, j_src)
         this%w_rem_v(i, j) = bt_rem_v(i_src, j_src)
         this%w_force_v(i, j) = force_v(i_src, j_src)
      end do
      if (refresh_widths) then
         do concurrent(j=1:ny + 2*off + 1, i=1:nx + 2*off) local(i_src, j_src)
            i_src = min(max(i - off, 1), nx)
            j_src = min(max(j - off, 1), ny + 1)
            this%w_dx_cv(i, j) = dx_cv_bt(i_src, j_src)
         end do
      end if
      if (this%ring_poison > 0) call bt_wide_poison_dynamic_ring(this)
   end subroutine bt_wide_copy_in_impl

   subroutine bt_wide_poison_dynamic_ring(this)
      !! Ring-poison hook: garbage into the outer ring of every per-solve wide
      !! input (after the copy-in, before the entry exchange).
      type(bt_wide_t), intent(inout) :: this
      integer :: n1, n2, off, m
      off = this%bt_halo + this%ring_poison/10
      m = mod(this%ring_poison, 10)
      ! `ubt_prev` / `vbt_prev` are NOT poisoned: they are no longer part of
      ! the packed entry exchange (nothing refills a seam ring of them), they
      ! are re-seeded from `ubt` / `vbt` over the whole wide array when
      ! `bebt > 0`, and at `bebt = 0` the kernel still evaluates
      ! `(1+bebt)*ubt - bebt*ubt_prev`, so a NaN there would be `0*NaN` in a
      ! cell whose value is not a physical dependency.
      n1 = size(this%w_eta, 1)
      n2 = size(this%w_eta, 2)
      call poison_ring_dev(this%w_eta, n1, n2, off, m)
      ! `H_ref` is copied and exchanged once: poison it with the first solve
      ! only, or a later solve would garble a seam ring nothing refills.
      if (.not. this%h_ref_ready) call poison_ring_dev(this%w_H_ref, n1, n2, off, m)
      n1 = size(this%w_ubt, 1)
      n2 = size(this%w_ubt, 2)
      call poison_ring_dev(this%w_ubt, n1, n2, off, m)
      call poison_ring_dev(this%w_rem_u, n1, n2, off, m)
      call poison_ring_dev(this%w_force_u, n1, n2, off, m)
      if (this%refresh_widths) call poison_ring_dev(this%w_dy_cu, n1, n2, off, m)
      n1 = size(this%w_vbt, 1)
      n2 = size(this%w_vbt, 2)
      call poison_ring_dev(this%w_vbt, n1, n2, off, m)
      call poison_ring_dev(this%w_rem_v, n1, n2, off, m)
      call poison_ring_dev(this%w_force_v, n1, n2, off, m)
      if (this%refresh_widths) call poison_ring_dev(this%w_dx_cv, n1, n2, off, m)
   end subroutine bt_wide_poison_dynamic_ring

   subroutine bt_wide_entry_exchange(this)
      !! Fills the wide ghost band before the fast loop with ONE packed
      !! exchange (eta, ubt, vbt, rem_u, force_u, rem_v, force_v).  `H_ref` is
      !! exchanged once, at the first call.  The previous-substep velocities
      !! are not exchanged: the fast loop re-seeds them from `ubt`/`vbt` over
      !! the whole wide array when `bebt > 0` and never reads them otherwise
      !! (their coefficient is 0).  Under porous barriers / closed faces the
      !! two per-solve BT widths (`w_dy_cu`, `w_dx_cv`) follow in one wide
      !! face exchange each (the packed set is full at three x / y arrays).
      !! Counter effect: +1 bt_packed per call (+1 centre_2d on the first;
      !! +1 face_x_2d, +1 face_y_2d more when the widths are refreshed).
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
      if (.not. this%h_ref_ready) then
         call ocean_halo_centre_2d_wide(this%w_H_ref, ng_w)
         this%h_ref_ready = .true.
      end if
      call ocean_halo_packed_wide(ng_w, this%w_eta, this%w_ubt, this%w_vbt, &
                                  xa=this%w_rem_u, xb=this%w_force_u, &
                                  ya=this%w_rem_v, yb=this%w_force_v)
      if (this%refresh_widths) then
         call ocean_halo_face_x_2d_wide(this%w_dy_cu, ng_w)
         call ocean_halo_face_y_2d_wide(this%w_dx_cv, ng_w)
      end if
      call profiler_stop("ocean_comms_bt")
   end subroutine bt_wide_entry_exchange_impl

   subroutine bt_wide_set_eta_forcing(this, grid, eta_forcing)
      !! Per-solve wide copy of `eta_forcing` (equilibrium tide / pressure
      !! loading; the PGF drives `grad(eta - eta_forcing)`) plus ONE wide
      !! centre exchange, so a ghost cell carries its owner's value.  Called
      !! only when the solve has a forcing (absent otherwise: the wide loop
      !! then gets no forcing, as the one-rank `bt_halo = 0` path).
      class(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: eta_forcing(grid%nx_total, grid%ny_total)
      select type (this)
      type is (bt_wide_t)
         call bt_wide_set_eta_forcing_impl(this, grid, eta_forcing)
      end select
   end subroutine bt_wide_set_eta_forcing

   subroutine bt_wide_set_eta_forcing_impl(this, grid, eta_forcing)
      !! Non-polymorphic body of `bt_wide_set_eta_forcing`.
      type(bt_wide_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: eta_forcing(grid%nx_total, grid%ny_total)
      integer :: i, j, nx, ny, off, i_src, j_src
      off = this%bt_halo
      nx = grid%nx_total
      ny = grid%ny_total
      if (size(this%w_eta_forcing, 1) /= nx + 2*off .or. size(this%w_eta_forcing, 2) /= ny + 2*off) then
         error stop "bt_wide_set_eta_forcing: wide eta_forcing not allocated (init without eta_forcing_on)"
      end if
      do concurrent(j=1:ny + 2*off, i=1:nx + 2*off) local(i_src, j_src)
         i_src = min(max(i - off, 1), nx)
         j_src = min(max(j - off, 1), ny)
         this%w_eta_forcing(i, j) = eta_forcing(i_src, j_src)
      end do
      if (this%ring_poison > 0) then
         call poison_ring_dev(this%w_eta_forcing, size(this%w_eta_forcing, 1), &
                              size(this%w_eta_forcing, 2), off + this%ring_poison/10, &
                              mod(this%ring_poison, 10))
      end if
      call profiler_start("ocean_comms_bt")
      call ocean_halo_centre_2d_wide(this%w_eta_forcing, this%ng_wide)
      call profiler_stop("ocean_comms_bt")
   end subroutine bt_wide_set_eta_forcing_impl

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
      if (allocated(this%f_corner_w)) deallocate (this%f_corner_w)
      if (allocated(this%w_area_cu)) deallocate (this%w_area_cu)
      if (allocated(this%w_area_cv)) deallocate (this%w_area_cv)
      if (allocated(this%w_dx_cu)) deallocate (this%w_dx_cu)
      if (allocated(this%w_dx_cv)) deallocate (this%w_dx_cv)
      if (allocated(this%w_dy_cu)) deallocate (this%w_dy_cu)
      if (allocated(this%w_dy_cv)) deallocate (this%w_dy_cv)
      if (allocated(this%w_iarea_bu)) deallocate (this%w_iarea_bu)
      if (allocated(this%w_iarea_t)) deallocate (this%w_iarea_t)
      if (allocated(this%w_idx_cu)) deallocate (this%w_idx_cu)
      if (allocated(this%w_idy_cv)) deallocate (this%w_idy_cv)
      if (allocated(this%w_eta_forcing)) deallocate (this%w_eta_forcing)
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

   subroutine bt_wide_substep(bt_wide, bt_work, n_steps, dt_inner, bc, t, tide)
      !! Wide-halo (march-in) entry point for the nonlinear barotropic fast
      !! loop.  Unpacks the wide shadow arrays (`w_*`, wide grid, the static
      !! wide metrics and `f_corner_w`) and forwards them to
      !! `barotropic_substep_nonlinear` with `bt_halo = bt_wide%bt_halo`, so
      !! the ~20-array plumbing lives here once rather than at the call site.
      !! `bc` and `t` propagate by absence; `t` is the outer-step model time
      !! the tidal open-edge target is evaluated at (as on the one-rank
      !! `bt_halo = 0` path).  `tide` selects the wide `eta_forcing`
      !! (`bt_wide_set_eta_forcing` ran this solve).  Interior twin:
      !! `barotropic_substep_nonlinear_interior`.
      type(bt_wide_t), intent(inout) :: bt_wide
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      type(ocean_bc_state_t), intent(inout), optional :: bc
      real(wp), intent(in), optional :: t
      logical, intent(in), optional :: tide

      logical :: tide_on
      tide_on = .false.
      if (present(tide)) tide_on = tide

      if (tide_on) then
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
                                           area_cu=bt_wide%w_area_cu, area_cv=bt_wide%w_area_cv, &
                                           dx_cu=bt_wide%w_dx_cu, dx_cv=bt_wide%w_dx_cv, &
                                           dy_cu=bt_wide%w_dy_cu, dy_cv=bt_wide%w_dy_cv, &
                                           iarea_bu=bt_wide%w_iarea_bu, iarea_t=bt_wide%w_iarea_t, &
                                           idx_cu=bt_wide%w_idx_cu, idy_cv=bt_wide%w_idy_cv, &
                                           bc=bc, t=t, eta_forcing=bt_wide%w_eta_forcing, &
                                           bt_halo=bt_wide%bt_halo)
      else
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
                                           area_cu=bt_wide%w_area_cu, area_cv=bt_wide%w_area_cv, &
                                           dx_cu=bt_wide%w_dx_cu, dx_cv=bt_wide%w_dx_cv, &
                                           dy_cu=bt_wide%w_dy_cu, dy_cv=bt_wide%w_dy_cv, &
                                           iarea_bu=bt_wide%w_iarea_bu, iarea_t=bt_wide%w_iarea_t, &
                                           idx_cu=bt_wide%w_idx_cu, idy_cv=bt_wide%w_idy_cv, &
                                           bc=bc, t=t, bt_halo=bt_wide%bt_halo)
      end if
   end subroutine bt_wide_substep

end module rdb_ocean_bt_wide
