module rdb_ice_dyn_c_state
   !! State slot of the C-grid EVP dynamics running under the SIS2-derived ice
   !! state (`&ocean_ice_dyn_nml dyn_scheme = "cice_c"`, CICE6 port D9.2).
   !!
   !! Under the default `dyn_scheme = "sis2_c"` this type is never initialised:
   !! nothing is allocated, mapped or registered for restart, so the default is
   !! bit-identical.  Under `"cice_c"` it holds, allocated once and mapped once
   !! (component arrays only), everything the C-grid solver of `rdb_ice_evp_c`
   !! carries or needs between two outer steps:
   !!
   !!   * the solver work/state type `w` (`ice_evp_c_work_t`: the T stresses,
   !!     the corner stress, the ice masks and the per-step scratch),
   !!   * (the C-grid metric set `ice_cgrid_metrics_t` is NOT here: its module
   !!     reaches the ocean halo and so the configuration, which this slot, a
   !!     component of the ice slot read by the configuration, cannot depend
   !!     on.  It lives in `ocean_state_t%ice_cgrid`.)
   !!   * the coupling inputs at E and N (`fcorE`, `strairxE`, `ss_tltxE`,
   !!     `uocnE`, `vocnE`, ... : rebuilt every step by `rdb_ice_dyn_c`),
   !!   * `strocnT_x/y`, the ice-ocean stress per ICE area at T (restart state,
   !!     read one step late by the basal heat flux) and the per-ice-area face
   !!     stress it is built from,
   !!   * `ice_err`, the per-cell error words of the step (`rdb_ice_errors`),
   !!   * `emask_r/nmask_r`, the real-valued mirror of the two face ice masks
   !!     (restart state: the logical masks of `w` decide the start velocity of a
   !!     newly covered face, so a bit-exact restart carries them).
   !!
   !! The ice velocity itself is NOT held here: it is `ocean_sea_ice_t%u_ice` /
   !! `v_ice` (the face arrays the transport already reads), which is why the
   !! category transport follows the new velocities with no further wiring.
   !!
   !! **Restart.**  Registered `optional` by `ocean_state_build_restart_registry`
   !! under `dyn_scheme = "cice_c"` only: `stresspT`, `stressmT`, `stress12U`
   !! (whole arrays, ghost rings included, so no exchange is needed after the
   !! read), the two mask mirrors and `strocnT_x/y`.  `zetax2T`, `etax2T` and
   !! every per-step scratch are recomputed.  The scalar `ice_dyn_scheme_id` (1) is
   !! registered with them: a checkpoint written under the other core (or before
   !! this one existed) is REFUSED on read instead of resuming half-mixed.
   !! `strocnT_x/y` is restart state that has NO CONSUMER until the basal-heat
   !! wiring (`ustar` of the next step) lands; it is kept so that wiring needs no
   !! restart change.
   !!
   !! **Notes.**  Open ocean boundaries are land for the ice (the dynamics
   !! envelope allows only wall or periodic edges, as the SIS2 dynamics does).
   !! The tilt force uses the free-surface height `sum_k h_layer - b`; when the
   !! ice mass loads the surface (not the case today) it must use the
   !! pressure-loaded surface instead.
   !! **Interim choice (logged).**  The ice strength is the Hibler (1979)
   !! strength of the aggregates (`kstrength = 0`, `strength = Pstar vice
   !! exp(-Cstar (1 - aice))`), because the SIS2-derived state has no CICE
   !! category arrays for the Rothrock strength; the ice and snow masses use the
   !! SIS2 densities so the mass seen by the dynamics is the mass of the state.
   !!
   !! Roundabout-original; no upstream code.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_mem_report, only: arr_bytes
   use rdb_ice_params, only: ice_params_t
   use rdb_ice_evp_c_prep, only: ice_evp_c_work_t
   implicit none
   private

   public :: ice_dyn_c_t

   type :: ice_dyn_c_t
      !! The slot.  `enable` is latched from the configuration BEFORE `init`.
      logical :: enable = .false.
         !! `dyn_scheme = "cice_c"` (and `ice%dynamics`).
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      integer :: nx = 0
      integer :: ny = 0
      real(wp) :: scheme_id = 1.0_wp
         !! Dynamics-scheme id carried in the checkpoint (`ice_dyn_scheme_id`):
         !! 1 = `"cice_c"`.  A file without the tag is the existing scheme (0).
      type(ice_params_t) :: p
         !! Parameter bundle (`ice_params_from_config`), mapped once.
      type(ice_evp_c_work_t) :: w
      real(wp), allocatable :: ssh(:, :)
         !! Sea-surface height `sum_k h_layer - b` at T [m], rebuilt every step.
      real(wp), allocatable :: fcorE(:, :), fcorN(:, :)
         !! Coriolis parameter at E and N [1/s].
      real(wp), allocatable :: strairxE(:, :), strairyN(:, :)
         !! Wind stress on the ice at E (x) and N (y): the ice-area-weighted
         !! snapshot of the wind stress felt by the ocean [N/m2].
      real(wp), allocatable :: ss_tltxE(:, :), ss_tltyN(:, :)
         !! Sea-surface slope along E (x) and N (y) [1].
      real(wp), allocatable :: uocnE(:, :), vocnE(:, :), uocnN(:, :), vocnN(:, :)
         !! Ocean surface velocity components at E and at N [m/s].
      real(wp), allocatable :: wE(:, :), wN(:, :)
         !! Face stress per ICE area, `strocnxE/aiE` and `strocnyN/aiN` (0 where
         !! the face concentration is below `puny`) [N/m2]; scratch of `strocnT`.
      real(wp), allocatable :: strocnT_x(:, :), strocnT_y(:, :)
         !! Ice-ocean stress on the ice per ICE area at T [N/m2] (restart state).
      real(wp), allocatable :: emask_r(:, :), nmask_r(:, :)
         !! 0/1 mirror of `w%iceEmask`, `w%iceNmask` (restart state).
      integer(int32), allocatable :: ice_err(:, :)
         !! Per-cell error words of the step.
   contains
      procedure :: init => dyn_c_init
      procedure :: destroy => dyn_c_destroy
      procedure :: enter_data => dyn_c_enter_data
      procedure :: exit_data => dyn_c_exit_data
      procedure :: bytes => dyn_c_bytes
   end type ice_dyn_c_t

contains

   subroutine dyn_c_init(this, nx, ny)
      !! Allocate and zero every array (host).  `this%p` must be set by the
      !! caller (`ice_params_from_config`); `kstrength = 0` is enforced
      !! at configure (Hibler strength, see the module note).
      class(ice_dyn_c_t), intent(inout) :: this
      integer, intent(in) :: nx, ny

      if (this%is_init) error stop "ice_dyn_c_t%init: already initialised"
      this%nx = nx
      this%ny = ny
      call this%w%init(nx, ny)
      allocate (this%ssh(nx, ny), source=0.0_wp)
      allocate (this%fcorE(nx, ny), this%fcorN(nx, ny), source=0.0_wp)
      allocate (this%strairxE(nx, ny), this%strairyN(nx, ny), source=0.0_wp)
      allocate (this%ss_tltxE(nx, ny), this%ss_tltyN(nx, ny), source=0.0_wp)
      allocate (this%uocnE(nx, ny), this%vocnE(nx, ny), source=0.0_wp)
      allocate (this%uocnN(nx, ny), this%vocnN(nx, ny), source=0.0_wp)
      allocate (this%wE(nx, ny), this%wN(nx, ny), source=0.0_wp)
      allocate (this%strocnT_x(nx, ny), this%strocnT_y(nx, ny), source=0.0_wp)
      allocate (this%emask_r(nx, ny), this%nmask_r(nx, ny), source=0.0_wp)
      allocate (this%ice_err(nx, ny), source=0_int32)
      this%is_init = .true.
   end subroutine dyn_c_init

   subroutine dyn_c_destroy(this)
      class(ice_dyn_c_t), intent(inout) :: this

      if (.not. this%is_init) return
      call this%w%destroy()
      deallocate (this%ssh, this%fcorE, this%fcorN, this%strairxE, this%strairyN)
      deallocate (this%ss_tltxE, this%ss_tltyN, this%uocnE, this%vocnE, this%uocnN, this%vocnN)
      deallocate (this%wE, this%wN, this%strocnT_x, this%strocnT_y, this%emask_r, this%nmask_r)
      deallocate (this%ice_err)
      this%is_init = .false.
   end subroutine dyn_c_destroy

   pure function dyn_c_bytes(this) result(nb)
      !! Counted allocatable footprint (0 when the slot is not initialised).
      class(ice_dyn_c_t), intent(in) :: this
      integer(int64) :: nb

      nb = 0_int64
      if (.not. this%is_init) return
      nb = this%w%bytes() &
           + arr_bytes(this%ssh) &
           + arr_bytes(this%fcorE) + arr_bytes(this%fcorN) &
           + arr_bytes(this%strairxE) + arr_bytes(this%strairyN) &
           + arr_bytes(this%ss_tltxE) + arr_bytes(this%ss_tltyN) &
           + arr_bytes(this%uocnE) + arr_bytes(this%vocnE) &
           + arr_bytes(this%uocnN) + arr_bytes(this%vocnN) &
           + arr_bytes(this%wE) + arr_bytes(this%wN) &
           + arr_bytes(this%strocnT_x) + arr_bytes(this%strocnT_y) &
           + arr_bytes(this%emask_r) + arr_bytes(this%nmask_r) &
           + int(size(this%ice_err), int64)*int(storage_size(this%ice_err)/8, int64)
   end function dyn_c_bytes

   subroutine dyn_c_enter_data(this)
      !! Map every array once.  The logical face masks of `w` are rebuilt from
      !! their restart mirrors on the host first (cold start: all false).
      class(ice_dyn_c_t), intent(inout) :: this

      if (.not. this%is_init) return
      this%w%iceEmask = this%emask_r > 0.5_wp
      this%w%iceNmask = this%nmask_r > 0.5_wp
      call this%w%enter_data()
      !$acc enter data copyin(this%p)
      !$acc enter data copyin(this%ssh, this%fcorE, this%fcorN, this%strairxE, this%strairyN)
      !$acc enter data copyin(this%ss_tltxE, this%ss_tltyN, this%uocnE, this%vocnE, this%uocnN, this%vocnN)
      !$acc enter data copyin(this%wE, this%wN, this%strocnT_x, this%strocnT_y, this%emask_r, this%nmask_r)
      !$acc enter data copyin(this%ice_err)
   end subroutine dyn_c_enter_data

   subroutine dyn_c_exit_data(this)
      !! Unmap.  The restart state is copied out; scratch is deleted.
      class(ice_dyn_c_t), intent(inout) :: this

      if (.not. this%is_init) return
      !$acc exit data delete(this%ssh, this%fcorE, this%fcorN, this%strairxE, this%strairyN)
      !$acc exit data delete(this%ss_tltxE, this%ss_tltyN, this%uocnE, this%vocnE, this%uocnN, this%vocnN)
      !$acc exit data delete(this%wE, this%wN, this%ice_err)
      !$acc exit data copyout(this%strocnT_x, this%strocnT_y, this%emask_r, this%nmask_r)
      !$acc exit data delete(this%p)
      !$acc update self(this%w%stresspT, this%w%stressmT, this%w%stress12U)
      call this%w%exit_data()
   end subroutine dyn_c_exit_data

end module rdb_ice_dyn_c_state
