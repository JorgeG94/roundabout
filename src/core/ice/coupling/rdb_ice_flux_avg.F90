!! Window-mean flux averager between the ice step and the ocean thermo step (roundabout-original).
module rdb_ice_flux_avg
   !! The CICE-port ice steps EVERY outer step (ice dt = outer dt, `ndtd = 1`),
   !! while the ocean applies surface fluxes at the thermo step
   !! (`is_thermo_step()`, every `dt_therm_ratio` outer steps; MOM6
   !! `THERMO_SPANS_COUPLING`).  `ice_flux_avg_t` holds the TIME-WEIGHTED MEAN of
   !! the ice's cell fluxes `fresh, fsalt, fhocn, fswthru, fswthru_vis`
   !! (MASTER_PLAN C.4: per unit cell area, per second, Icepack signs) over the
   !! ice steps since the last hand-off.  The ocean applies that mean for the
   !! whole window, so `sum_cells A*applied*window = sum_steps sum_cells A*flux*dt`
   !! (closure B).
   !!
   !! **Running mean, not `acc += flux*dt`.**  The plan writes the sum form
   !! `acc += flux*dt; Q = acc/(n*dt)`.  This module keeps the mathematically
   !! identical running mean `acc <- acc + (flux - acc)*w`, `w = dt/(t_acc + dt)`,
   !! because the sum form is NOT bitwise the direct hand-off at ratio 1
   !! (`(f*dt)/dt /= f` for some `f`).  With the running mean the first step of a
   !! window has `t_acc = 0`, so `w = dt/dt = 1` exactly and
   !! `acc = 0 + (f - 0)*1 = f` bit for bit.  Later steps are weighted by their
   !! own `dt` (the weight is the ice step, so a changed `dt` mid-window stays
   !! exact); with equal steps `acc` is the arithmetic mean to round-off.  A
   !! non-finite flux enters `acc` as NaN (no clamp, no mask) and stays until
   !! `ice_flux_avg_reset`.
   !!
   !! **Window = the steps between two resets.**  The averager does not know the
   !! cadence: the caller accumulates after every ice step and, at
   !! `is_thermo_step()`, hands `acc_*` to the adapters
   !! (`rdb_ice_cice_adapters`, which also call `ice_flux_avg_reset`).
   !! `n_acc` counts the steps in the open window, `t_acc` the seconds.
   !!
   !! **Valid on OWNED cells only** (the merged ice fluxes are); ghosts are zero.
   !!
   !! **Restart.**  A checkpoint in the middle of a window must resume bitwise:
   !! `register_restart` adds the five mean arrays under the tags
   !! `ice_avg_{fresh,fsalt,fhocn,fswthru,fswthru_vis}` (owned slice, like every
   !! registry entry) and the host scalars `ice_avg_n`, `ice_avg_t`, `optional` by
   !! default and REQUIRED (`register_restart(reg, required=.true.)`) when the file was
   !! written by the new scheme, so a partial window cannot be dropped silently.  The slot is initialised, mapped and registered ONLY under the
   !! new scheme; under `scheme = "sis2"` nothing exists, so no tag is added
   !! (7 entries of `MAX_RESTART_ENTRIES = 256`).
   !!
   !! **Lifecycle** (`mem:separate`): `init` allocates zeros; `enter_data` maps the
   !! COMPONENT arrays with `copyin` (host zeros or restart-seeded values);
   !! `exit_data` deletes them; `bytes()` counts every array.
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_flux_state, only: ice_flux_t
   implicit none
   private

   public :: ice_flux_avg_t, ice_flux_avg_accumulate, ice_flux_avg_reset, ice_flux_avg_stash_pending
   public :: ICE_AVG_OK, ICE_AVG_ERR_SHAPE, ICE_AVG_ERR_DT

   integer, parameter :: ICE_AVG_OK = 0
      !! Call ran.
   integer, parameter :: ICE_AVG_ERR_SHAPE = 1
      !! Slots uninitialised or disagree on shape; nothing ran.
   integer, parameter :: ICE_AVG_ERR_DT = 2
      !! `dt` not finite and positive; nothing ran.

   type :: ice_flux_avg_t
      !! The window-mean accumulator slot.
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data`.
      integer :: nx = 0
         !! Total extent in x including ghosts.
      integer :: ny = 0
         !! Total extent in y including ghosts.
      integer :: nghost = 0
         !! Ghost width.
      integer :: nx_phys = 0
         !! Owned cells in x.
      integer :: ny_phys = 0
         !! Owned cells in y.
      real(wp) :: n_acc = 0.0_wp
         !! Ice steps in the open window (real so the restart registry can carry it).
      real(wp) :: t_acc = 0.0_wp
         !! Seconds in the open window.
      real(wp), allocatable :: acc_fresh(:, :)
         !! [kg m-2 s-1] window mean of `fresh`, + into the ocean.
      real(wp), allocatable :: acc_fsalt(:, :)
         !! [kg m-2 s-1] window mean of `fsalt`, + salinifies.
      real(wp), allocatable :: acc_fhocn(:, :)
         !! [W m-2] window mean of `fhocn`, + into the ocean (never frazil).
      real(wp), allocatable :: acc_fswthru(:, :)
         !! [W m-2] window mean of `fswthru`, + down.
      real(wp), allocatable :: acc_fswthru_vis(:, :)
         !! [W m-2] window mean of `fswthru_vis`, + down (visible part).
      real(wp) :: t_pend = 0.0_wp
         !! Seconds of the PENDING hand-over: the last window's mean `fhocn` that the ocean has not
         !! debited yet (it applies the hand-over of window k at the END of window k+1).  The melt cap
         !! subtracts it too.  0: nothing pending.
      real(wp), allocatable :: pend_fhocn(:, :)
         !! [W m-2] mean `fhocn` of the pending window (see `t_pend`).
   contains
      procedure, non_overridable :: init => avg_init
      procedure, non_overridable :: destroy => avg_destroy
      procedure, non_overridable :: enter_data => avg_enter_data
      procedure, non_overridable :: exit_data => avg_exit_data
      procedure, non_overridable :: bytes => avg_bytes
      procedure, non_overridable :: register_restart => avg_register_restart
   end type ice_flux_avg_t

contains

   subroutine avg_init(this, grid)
      !! Allocate the five mean arrays at the ghosted extent, zeroed; empty window.
      class(ice_flux_avg_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      if (this%is_init) call this%destroy()
      this%nx = grid%nx_total
      this%ny = grid%ny_total
      this%nghost = grid%nghost
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%n_acc = 0.0_wp
      this%t_acc = 0.0_wp
      allocate (this%acc_fresh(this%nx, this%ny), source=0.0_wp)
      allocate (this%acc_fsalt(this%nx, this%ny), source=0.0_wp)
      allocate (this%acc_fhocn(this%nx, this%ny), source=0.0_wp)
      allocate (this%acc_fswthru(this%nx, this%ny), source=0.0_wp)
      allocate (this%acc_fswthru_vis(this%nx, this%ny), source=0.0_wp)
      allocate (this%pend_fhocn(this%nx, this%ny), source=0.0_wp)
      this%t_pend = 0.0_wp
      this%is_init = .true.
      this%mapped = .false.
   end subroutine avg_init

   subroutine avg_destroy(this)
      !! Reverse of `init` (unmaps first).
      class(ice_flux_avg_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%acc_fresh)) deallocate (this%acc_fresh)
      if (allocated(this%acc_fsalt)) deallocate (this%acc_fsalt)
      if (allocated(this%acc_fhocn)) deallocate (this%acc_fhocn)
      if (allocated(this%acc_fswthru)) deallocate (this%acc_fswthru)
      if (allocated(this%acc_fswthru_vis)) deallocate (this%acc_fswthru_vis)
      if (allocated(this%pend_fhocn)) deallocate (this%pend_fhocn)
      this%nx = 0
      this%ny = 0
      this%nghost = 0
      this%nx_phys = 0
      this%ny_phys = 0
      this%n_acc = 0.0_wp
      this%t_acc = 0.0_wp
   end subroutine avg_destroy

   subroutine avg_enter_data(this)
      !! Map the component arrays (poly TBP delegating to a `type(...)` impl).
      class(ice_flux_avg_t), intent(inout) :: this
      select type (this)
      type is (ice_flux_avg_t)
         call avg_enter_data_impl(this)
      end select
   end subroutine avg_enter_data

   subroutine avg_enter_data_impl(this)
      !! `copyin` every component array (host zeros or restart-seeded values).
      type(ice_flux_avg_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%acc_fresh, this%acc_fsalt, this%acc_fhocn, this%acc_fswthru, this%acc_fswthru_vis)
      !$acc enter data copyin(this%pend_fhocn)
      this%mapped = .true.
   end subroutine avg_enter_data_impl

   subroutine avg_exit_data(this)
      !! Unmap the component arrays (poly TBP delegating to a `type(...)` impl).
      class(ice_flux_avg_t), intent(inout) :: this
      select type (this)
      type is (ice_flux_avg_t)
         call avg_exit_data_impl(this)
      end select
   end subroutine avg_exit_data

   subroutine avg_exit_data_impl(this)
      !! `copyout` is the restart path's `update self`; here the arrays are deleted.
      type(ice_flux_avg_t), intent(inout) :: this
      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%acc_fresh, this%acc_fsalt, this%acc_fhocn, this%acc_fswthru, this%acc_fswthru_vis)
      !$acc exit data delete(this%pend_fhocn)
      this%mapped = .false.
   end subroutine avg_exit_data_impl

   pure function avg_bytes(this) result(nbytes)
      !! Counted allocatable footprint (one `arr_bytes` term per array).
      class(ice_flux_avg_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = 0_int64 &
               + arr_bytes(this%acc_fresh) &
               + arr_bytes(this%acc_fsalt) &
               + arr_bytes(this%acc_fhocn) &
               + arr_bytes(this%acc_fswthru) &
               + arr_bytes(this%acc_fswthru_vis) &
               + arr_bytes(this%pend_fhocn) &
               + 0_int64
   end function avg_bytes

   subroutine avg_register_restart(this, reg, required)
      !! Register the mid-window state (see the module docstring), host scalars for the
      !! window counters.  Call ONLY when the new scheme is active.  No-op on an
      !! uninitialised slot.  `required = .true.` (the caller sets it when the file was
      !! written by the new scheme, i.e. carries the ice state version scalar) registers
      !! every entry NON-optional, so a checkpoint that lost part of a window fails to read
      !! instead of silently resuming with an empty one; default `.false.` (a file from the
      !! SIS2 scheme, converted, legitimately has no window).
      class(ice_flux_avg_t), intent(in), target :: this
      type(restart_registry_t), intent(inout) :: reg
      logical, intent(in), optional :: required
      logical :: opt

      if (.not. this%is_init) return
      opt = .true.
      if (present(required)) opt = .not. required
      call reg%register_2d("ice_avg_fresh", this%acc_fresh, this%nghost, this%nx_phys, this%ny_phys, optional=opt)
      call reg%register_2d("ice_avg_fsalt", this%acc_fsalt, this%nghost, this%nx_phys, this%ny_phys, optional=opt)
      call reg%register_2d("ice_avg_fhocn", this%acc_fhocn, this%nghost, this%nx_phys, this%ny_phys, optional=opt)
      call reg%register_2d("ice_avg_fswthru", this%acc_fswthru, this%nghost, this%nx_phys, this%ny_phys, &
                           optional=opt)
      call reg%register_2d("ice_avg_fswthru_vis", this%acc_fswthru_vis, this%nghost, this%nx_phys, &
                           this%ny_phys, optional=opt)
      call reg%register_2d("ice_avg_pend_fhocn", this%pend_fhocn, this%nghost, this%nx_phys, this%ny_phys, &
                           optional=.true.)
      call reg%register_scalar("ice_avg_pend_t", this%t_pend, optional=.true.)
      call reg%register_scalar("ice_avg_n", this%n_acc, optional=opt)
      call reg%register_scalar("ice_avg_t", this%t_acc, optional=opt)
   end subroutine avg_register_restart

   subroutine ice_flux_avg_accumulate(avg, flux, dt, ierr)
      !! Fold one ice step's cell fluxes into the window mean over the OWNED cells,
      !! ONE launch: `acc <- acc + (flux - acc)*w`, `w = dt/(t_acc + dt)`; then
      !! `t_acc += dt`, `n_acc += 1`.  Call after the ice step, with the merged
      !! fluxes of THAT step (before `ice_flux_zero` of the next).  `ierr =
      !! ICE_AVG_ERR_SHAPE` / `ICE_AVG_ERR_DT` (nothing ran) on a shape or `dt`
      !! error.  Both slots must be mapped on an offload build.
      type(ice_flux_avg_t), intent(inout) :: avg
      type(ice_flux_t), intent(in) :: flux
      real(wp), intent(in) :: dt
         !! The ice step [s].
      integer, intent(out) :: ierr

      real(wp) :: w

      ierr = ICE_AVG_OK
      if (.not. (avg%is_init .and. flux%is_init) .or. avg%nx /= flux%nx .or. avg%ny /= flux%ny .or. &
          avg%nghost /= flux%nghost) then
         ierr = ICE_AVG_ERR_SHAPE
         return
      end if
      if (.not. ieee_is_finite(dt)) then
         ierr = ICE_AVG_ERR_DT
         return
      end if
      if (dt <= 0.0_wp) then
         ierr = ICE_AVG_ERR_DT
         return
      end if
      w = dt/(avg%t_acc + dt)
      call ice_flux_avg_accumulate_impl(avg%nx, avg%ny, avg%nghost + 1, avg%nghost + avg%nx_phys, &
                                        avg%nghost + 1, avg%nghost + avg%ny_phys, w, &
                                        flux%fresh, flux%fsalt, flux%fhocn, flux%fswthru, flux%fswthru_vis, &
                                        avg%acc_fresh, avg%acc_fsalt, avg%acc_fhocn, avg%acc_fswthru, &
                                        avg%acc_fswthru_vis)
      avg%t_acc = avg%t_acc + dt
      avg%n_acc = avg%n_acc + 1.0_wp
   end subroutine ice_flux_avg_accumulate

   pure subroutine ice_flux_avg_accumulate_impl(nx, ny, is, ie, js, je, w, &
                                                fresh, fsalt, fhocn, fswthru, fswthru_vis, &
                                                a_fresh, a_fsalt, a_fhocn, a_fswthru, a_fswthru_vis)
      !! Device kernel: one `do concurrent (j, i)` over the owned cells updating
      !! the five running means (no loop index reaches a call).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: w
      real(wp), intent(in) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), fswthru(nx, ny), fswthru_vis(nx, ny)
      real(wp), intent(inout) :: a_fresh(nx, ny), a_fsalt(nx, ny), a_fhocn(nx, ny), a_fswthru(nx, ny)
      real(wp), intent(inout) :: a_fswthru_vis(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         a_fresh(i, j) = a_fresh(i, j) + (fresh(i, j) - a_fresh(i, j))*w
         a_fsalt(i, j) = a_fsalt(i, j) + (fsalt(i, j) - a_fsalt(i, j))*w
         a_fhocn(i, j) = a_fhocn(i, j) + (fhocn(i, j) - a_fhocn(i, j))*w
         a_fswthru(i, j) = a_fswthru(i, j) + (fswthru(i, j) - a_fswthru(i, j))*w
         a_fswthru_vis(i, j) = a_fswthru_vis(i, j) + (fswthru_vis(i, j) - a_fswthru_vis(i, j))*w
      end do
   end subroutine ice_flux_avg_accumulate_impl

   subroutine ice_flux_avg_stash_pending(avg)
      !! At the ocean thermo step, BEFORE the window is applied and reset: keep the window's mean
      !! `fhocn` and its length as the PENDING hand-over (the ocean debits it only at the end of the
      !! NEXT window), which the melt cap of the next window's imports subtracts.  Owned cells; one
      !! launch.
      type(ice_flux_avg_t), intent(inout) :: avg

      if (.not. avg%is_init) return
      call ice_flux_avg_copy_impl(avg%nx, avg%ny, avg%acc_fhocn, avg%pend_fhocn)
      avg%t_pend = avg%t_acc
   end subroutine ice_flux_avg_stash_pending

   subroutine ice_flux_avg_copy_impl(nx, ny, src, dst)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: src(nx, ny)
      real(wp), intent(inout) :: dst(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         dst(i, j) = src(i, j)
      end do
   end subroutine ice_flux_avg_copy_impl

   subroutine ice_flux_avg_reset(avg)
      !! Open a new window: zero the five means (ghosts included, ONE launch) and
      !! the window counters.  Called right after the hand-off.  The slot must be
      !! mapped on an offload build.
      type(ice_flux_avg_t), intent(inout) :: avg

      if (.not. avg%is_init) return
      call ice_flux_avg_reset_impl(avg%nx, avg%ny, avg%acc_fresh, avg%acc_fsalt, avg%acc_fhocn, &
                                   avg%acc_fswthru, avg%acc_fswthru_vis)
      avg%t_acc = 0.0_wp
      avg%n_acc = 0.0_wp
   end subroutine ice_flux_avg_reset

   pure subroutine ice_flux_avg_reset_impl(nx, ny, a_fresh, a_fsalt, a_fhocn, a_fswthru, a_fswthru_vis)
      !! Device kernel zeroing the five means.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: a_fresh(nx, ny), a_fsalt(nx, ny), a_fhocn(nx, ny), a_fswthru(nx, ny)
      real(wp), intent(inout) :: a_fswthru_vis(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         a_fresh(i, j) = 0.0_wp
         a_fsalt(i, j) = 0.0_wp
         a_fhocn(i, j) = 0.0_wp
         a_fswthru(i, j) = 0.0_wp
         a_fswthru_vis(i, j) = 0.0_wp
      end do
   end subroutine ice_flux_avg_reset_impl

end module rdb_ice_flux_avg
