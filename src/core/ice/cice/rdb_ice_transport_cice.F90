!! Sea-ice horizontal transport step of the CICE-port scheme: checks, pack,
!! upwind advection, conservation test, unpack, post-transport state update.
module rdb_ice_transport_cice
   !! Roundabout-original driver (no upstream ancestry; the kernels it calls
   !! carry their own attribution).  Donor-cell advection: Lipscomb and Hunke
   !! (2004, Mon. Wea. Rev. 132 1341) for the tracer hierarchy it carries.
   !!
   !! **What one call does** (`ice_step_transport`, in this order; every step
   !! is one device launch over the cells, no host loop over cells, no
   !! transfer of fields):
   !!
   !! 0. host: the configuration (`advection`, counts, range, tracer table) is
   !!    checked before any launch;
   !! 1. ALWAYS: face copy + CFL check on the velocities (a face displacement
   !!    above the narrower of its two cells, either sign, is an error, never
   !!    clamped);
   !! 2. pack the state into the workspace (`transport_pack`, with
   !!    `thick_floor = UPWIND_THICK_FLOOR`: the volume of a category of
   !!    negligible area is carried, as the reference's extensive fields do);
   !! 3. non-finite and negative-area entry checks on the packed fields (per
   !!    cell; always when the pack counted a non-finite value, and under
   !!    `transport_check`), and under `transport_check` the EFP conservation
   !!    sums "before" and the monotonicity bounds (only if the workspace was
   !!    built `want_bounds`);
   !! 4. if any entry check fired: STOP here (the state was never written);
   !! 5. upwind advection into the second workspace (`rdb_ice_upwind`);
   !! 6. the advection counts non-finite new products and new areas below
   !!    `-puny`; when it counted any, and under `transport_check`, the per-cell
   !!    checks run on the result; under `transport_check` also the EFP sums
   !!    "after", their comparison and the monotonicity check;
   !! 7. if any fatal bit fired: STOP (the state is untouched: the advection
   !!    wrote only the second workspace);
   !! 8. unpack over the owned cells, then `upwind_empty_tsfc`
   !!    (`Tsfc = Tf` on an emptied category, `0` on land);
   !! 9. the post step selected by `post`:
   !!    `ICE_TRANSPORT_POST_NONE` (state only; aggregates stale),
   !!    `ICE_TRANSPORT_POST_AGGREGATE` (`ice_aggregate_state_impl`: `aice`,
   !!    `vice`, `vsno`, and `aice0 = max(1 - aice, 0)`) or
   !!    `ICE_TRANSPORT_POST_CLEANUP` (the clean-up of the ridging step,
   !!    `ice_cleanup_itd_step` with rebin, `limit_aice` and the aggregate).
   !!
   !! **The reference sequence, and which post step the step driver uses.**
   !! The reference runs transport, then the ridging (which closes the ADVECTED
   !! open-water fraction `aice0` and whose own clean-up removes slivers), then
   !! the aggregate.  The step-sequencing PR must therefore call this driver
   !! with `ICE_TRANSPORT_POST_NONE`, run the ridging on the advected state
   !! (the advected open-water category stays in `aice0`), and aggregate after
   !! it.  `POST_AGGREGATE` before the ridging would replace the open water by
   !! `max(1 - aice, 0)` and so change the closing rate; `POST_CLEANUP` raises
   !! the FATAL `ICE_ERR_AICE_RANGE` bit on every convergent cell (aice above
   !! `1 + puny`, which only the ridging removes).  Both are for tests and
   !! diagnostics, not for the model step.
   !!
   !! **Fatal means "step not taken", and the cell error word says so.**  The
   !! error bits of `rdb_ice_errors` (`ICE_ERR_NONFINITE`,
   !! `_TRANSPORT_CFL`, `_TRANSPORT_NEG_AREA`, `_TRANSPORT_CONSERV`) are ORed
   !! into `ice_err(i,j)` AND into `res%ierr`, with `transport_check` on or
   !! off, so a caller that enforces through the cell error word (`ice_err`
   !! reduction and `ice_err_enforce`) never sees a silently skipped
   !! transport: a refused step has `res%stepped = .false.` and at least one
   !! cell bit.  With the checks off the attribution is by cell of the result
   !! (a NaN velocity marks the cells it reached), not of the first offence.
   !! A CFL violation is counted (`res%n_cfl` cells), never clamped.  A NaN
   !! face velocity refuses the step even on a masked (land) face: the
   !! dynamics must hand over finite velocities on every face.
   !!
   !! **What must be valid on entry** (owned range `is:ie, js:je`, one ghost
   !! ring around it, `is >= 2`, `ie <= nx - 1`, likewise in `j`):
   !! the state arrays (`aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon`,
   !! `aice0`) over `is-1:ie+1, js-1:je+1` (the ring is a copy of the
   !! neighbour's owned cells: after a halo exchange / periodic wrap /
   !! fold); `u(is:ie+1, js-1:je)` and `v(is-1:ie, js:je+1)` (the owned faces
   !! and the ring faces the CFL check reads; the seam face is stored on both
   !! sides, equal); `HTE, HTN, epm, npm` at `i = is-1:ie`, `j = js-1:je`;
   !! `tarea` over the owned cells; `dxT, dyT` over `is-1:ie+1, js-1:je+1`;
   !! `Tf` and `hm` on the owned cells.  Slivers (`0 < aicen <= puny`) are
   !! allowed: their volume is carried (the reference removes them in the
   !! ridging step's clean-up).  **On exit** the state is advanced on the OWNED
   !! cells only; its ghost ring is the old one and MUST be refreshed (halo
   !! exchange, periodic wrap, fold) before anything reads it.  Multi-rank
   !! exchange points, for the multi-rank PR: (a) state ring before the call,
   !! (b) face velocities (already needed by the dynamics), (c) the EFP sums
   !! `before`/`after` need a bin-wise `MPI_SUM` before the comparison
   !! (`efp_to_transport`), (d) the counts in `res` need a sum/max over ranks,
   !! (e) the state ring after the call.
   !!
   !! **Conservation** is to round-off (one division by the cell area, the
   !! packed-frame re-multiply, possible FMA asymmetry of a face flux on the
   !! device), checked under `transport_check` by
   !! `transport_conservation_compare_floor`: per quantity the change must stay
   !! within `puny` times the quantity, or, for a nearly empty category, within
   !! `puny * (total area) * (the category's mean of that field)` (the area
   !! tracer content of cells whose new area is `<= puny` is dropped by design,
   !! so a relative-only test would fail on categories that hold almost
   !! nothing).  A genuine loss of a non-negligible amount still trips it.
   !!
   !! **Memory.**  Two packed workspaces (`2 x (1 + ntrace*ncat/(ncat+1)) *
   !! (ncat + 1)` fields of the grid), the check scratch (`rows`, `tot`, and
   !! `tmin`/`tmax` only with `want_bounds`) and two face copies; mapped once
   !! (`enter_data`), never allocated per call.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_to_real
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, &
                             ICE_ERR_TRANSPORT_NEG_AREA, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_params, only: ice_params_t, ice_transport_params_valid, ICE_ADVECT_UPWIND
   use rdb_ice_aggregate, only: ice_aggregate_state_impl
   use rdb_ice_itd_cleanup, only: ice_cleanup_itd_step
   use rdb_ice_transport_pack, only: transport_workspace_t, transport_tracer_table_build, TT_OK, &
                                     transport_pack, transport_unpack
   use rdb_ice_transport_checks, only: transport_check_workspace_t, transport_nquantities, &
                                       transport_check_finite, transport_check_neg_area, transport_check_cfl, &
                                       transport_conserved_sums, transport_bounds, transport_check_monotone
   use rdb_ice_upwind, only: upwind_advect, upwind_empty_tsfc, upwind_table_ok, UPWIND_CFG_OK, UPWIND_CFG_DEPEND, &
                             UPWIND_THICK_FLOOR
   implicit none
   private

   public :: ICE_TRANSPORT_POST_NONE, ICE_TRANSPORT_POST_AGGREGATE, ICE_TRANSPORT_POST_CLEANUP
   public :: TRANSPORT_CFG_OK, TRANSPORT_CFG_PARAMS, TRANSPORT_CFG_WORKSPACE, TRANSPORT_CFG_TABLE, TRANSPORT_CFG_RANGE
   public :: transport_cice_ws_t, transport_result_t
   public :: ice_step_transport, transport_conservation_compare_floor

   integer, parameter :: ICE_TRANSPORT_POST_NONE = 0
      !! Unpack only.
   integer, parameter :: ICE_TRANSPORT_POST_AGGREGATE = 1
      !! Unpack, then `aice`, `vice`, `vsno`, `aice0 = max(1 - aice, 0)`.
   integer, parameter :: ICE_TRANSPORT_POST_CLEANUP = 2
      !! Unpack, then the small-area / rebin / limiter clean-up with the aggregates.

   integer, parameter :: TRANSPORT_CFG_OK = 0
      !! Call accepted.
   integer, parameter :: TRANSPORT_CFG_PARAMS = 1
      !! `advection` is not `ICE_ADVECT_UPWIND`, or `ncat`/`nilyr`/`nslyr` differ from the workspace's.
   integer, parameter :: TRANSPORT_CFG_WORKSPACE = 2
      !! The workspace is not initialised / not the shape of the arrays.
   integer, parameter :: TRANSPORT_CFG_TABLE = 3
      !! The tracer table carries a dependent of a parent tracer (not supported by the upwind kernel).
   integer, parameter :: TRANSPORT_CFG_RANGE = 4
      !! The range has no room for the one-cell ghost ring.

   type :: transport_cice_ws_t
      !! Persistent buffers of the driver (see the module docstring, Memory).
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0, ny = 0, ncat = 0, nilyr = 0, nslyr = 0, ntrace = 0, nq = 0
      type(transport_workspace_t) :: a
         !! Packed input.
      type(transport_workspace_t) :: b
         !! Advected output.
      type(transport_check_workspace_t) :: chk
      type(efp_t), allocatable :: before(:), after(:)
      real(wp), allocatable :: uE(:, :), vN(:, :)
         !! Face velocities shifted to the T-cell index (`uE(i,j) = u(i+1,j)`) for the CFL check.
   contains
      procedure :: init => tcw_init
      procedure :: enter_data => tcw_enter_data
      procedure :: exit_data => tcw_exit_data
      procedure :: destroy => tcw_destroy
      procedure :: bytes => tcw_bytes
   end type transport_cice_ws_t

   type :: transport_result_t
      !! What one call did (host scalars).
      integer(int32) :: ierr = ICE_ERR_NONE
         !! OR of the global bits raised by the transport checks.
      integer :: cfg = TRANSPORT_CFG_OK
         !! `TRANSPORT_CFG_*`: nonzero means nothing was run.
      logical :: stepped = .false.
         !! The state was advanced.
      integer :: n_cfl = 0
         !! Faces violating the CFL check (check on).
      integer :: n_slivers = 0
         !! Sliver categories at entry (check on).
      integer :: n_nonfinite = 0
         !! Non-finite packed values at entry + non-finite new products (always) plus entry cells (check on).
      integer :: n_neg = 0
         !! Categories with a new area below `-puny` (always) / entry cells (check on).
      integer :: n_mono = 0
         !! Cells outside the monotonicity bound (bounds on).
      real(wp) :: worst_rel = 0.0_wp
         !! Largest relative change of a conserved quantity (check on).
      integer :: worst_q = 0
         !! Index of that quantity.
      logical :: conserv_done = .false.
         !! The conservation comparison ran (check on and the range small enough for the EFP sums).
      integer :: n_launch = 0
         !! Device launches issued by the driver (kernels, not the one-off host work).
   end type transport_result_t

contains

   ! ------------------------------------------------------------------
   ! workspace
   ! ------------------------------------------------------------------

   subroutine tcw_init(this, nx, ny, ncat, nilyr, nslyr, want_bounds, ierr)
      !! Build the table for `nilyr`/`nslyr`, allocate both packed workspaces
      !! and the check scratch.  `ierr /= 0` (a `TT_ERR_*` or
      !! `TRANSPORT_CFG_*`): nothing is left allocated.
      class(transport_cice_ws_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      logical, intent(in) :: want_bounds
      integer, intent(out) :: ierr
      type(transport_workspace_t) :: tmp
      integer :: ierr2

      ierr = 0
      call this%destroy()
      if (nslyr /= 1) then
         ierr = TRANSPORT_CFG_PARAMS
         return
      end if
      call transport_tracer_table_build(nilyr, nslyr, tmp%tab)
      if (tmp%tab%ierr /= TT_OK) then
         ierr = tmp%tab%ierr
         return
      end if
      call this%a%init(nx, ny, ncat, tmp%tab, ierr2)
      if (ierr2 /= 0) then
         ierr = ierr2
         return
      end if
      call this%b%init(nx, ny, ncat, tmp%tab, ierr2)
      if (ierr2 /= 0) then
         call this%destroy()
         ierr = ierr2
         return
      end if
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%nilyr = nilyr
      this%nslyr = nslyr
      this%ntrace = tmp%tab%ntrace
      this%nq = transport_nquantities(ncat, this%ntrace)
      call this%chk%init(nx, ny, ncat, this%ntrace, want_bounds)
      allocate (this%before(this%nq), this%after(this%nq))
      allocate (this%uE(nx, ny), this%vN(nx, ny))
      this%uE = 0.0_wp
      this%vN = 0.0_wp
      this%is_init = .true.
   end subroutine tcw_init

   subroutine tcw_enter_data(this)
      !! Map the component arrays once (never the aggregate).
      class(transport_cice_ws_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      call this%a%enter_data()
      call this%b%enter_data()
      call this%chk%enter_data()
      !$acc enter data copyin(this%uE, this%vN)
      this%mapped = .true.
   end subroutine tcw_enter_data

   subroutine tcw_exit_data(this)
      class(transport_cice_ws_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%uE, this%vN)
      call this%chk%exit_data()
      call this%b%exit_data()
      call this%a%exit_data()
      this%mapped = .false.
   end subroutine tcw_exit_data

   subroutine tcw_destroy(this)
      class(transport_cice_ws_t), intent(inout) :: this
      call this%exit_data()
      call this%a%destroy()
      call this%b%destroy()
      call this%chk%destroy()
      if (allocated(this%before)) deallocate (this%before)
      if (allocated(this%after)) deallocate (this%after)
      if (allocated(this%uE)) deallocate (this%uE)
      if (allocated(this%vN)) deallocate (this%vN)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%nilyr = 0
      this%nslyr = 0
      this%ntrace = 0
      this%nq = 0
   end subroutine tcw_destroy

   pure function tcw_bytes(this) result(b)
      !! Bytes of the grid-sized buffers (host copy; the device mirror is the same).
      class(transport_cice_ws_t), intent(in) :: this
      integer(int64) :: b
      b = this%a%bytes() + this%b%bytes() + this%chk%bytes()
      if (allocated(this%uE)) b = b + int(size(this%uE), int64)*8_int64
      if (allocated(this%vN)) b = b + int(size(this%vN), int64)*8_int64
   end function tcw_bytes

   ! ------------------------------------------------------------------
   ! driver
   ! ------------------------------------------------------------------

   subroutine ice_step_transport(w, p, dt, post, is, ie, js, je, u, v, HTE, HTN, tarea, dxT, dyT, epm, npm, hm, Tf, &
                                 hin_max, aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                 aice, vice, vsno, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ice_err, res)
      !! One transport step of the owned cells `is:ie, js:je` over the
      !! SoA state (see the module docstring for the order, the entry/exit
      !! contract and what the post step is).  Every array must be
      !! device-present on an offload build; the parameter bundle `p` is read
      !! on the host and its scalars passed to the kernels.  `dt` is the
      !! dynamics step [s].  `hin_max(0:ncat)` and the flux/aggregate arrays are
      !! read or written only by `ICE_TRANSPORT_POST_CLEANUP`
      !! (`aice`, `vice`, `vsno` also by `_AGGREGATE`).  `ice_err` is ORed into,
      !! never reset.
      type(transport_cice_ws_t), intent(inout) :: w
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      integer, intent(in) :: post, is, ie, js, je
      real(wp), intent(in) :: u(w%nx + 1, w%ny), v(w%nx, w%ny + 1)
      real(wp), intent(in) :: HTE(w%nx, w%ny), HTN(w%nx, w%ny), tarea(w%nx, w%ny), dxT(w%nx, w%ny), dyT(w%nx, w%ny)
      real(wp), intent(in) :: epm(w%nx, w%ny), npm(w%nx, w%ny), hm(w%nx, w%ny), Tf(w%nx, w%ny)
      real(wp), intent(in) :: hin_max(0:w%ncat)
      real(wp), intent(inout) :: aice0(w%nx, w%ny), aicen(w%nx, w%ny, w%ncat), vicen(w%nx, w%ny, w%ncat)
      real(wp), intent(inout) :: vsnon(w%nx, w%ny, w%ncat), Tsfcn(w%nx, w%ny, w%ncat)
      real(wp), intent(inout) :: qicen(w%nx, w%ny, w%ncat, w%nilyr), sicen(w%nx, w%ny, w%ncat, w%nilyr)
      real(wp), intent(inout) :: qsnon(w%nx, w%ny, w%ncat, w%nslyr)
      real(wp), intent(inout) :: aice(w%nx, w%ny), vice(w%nx, w%ny), vsno(w%nx, w%ny)
      real(wp), intent(inout) :: Tsfc_c(w%nx, w%ny), qice_c(w%nx, w%ny, w%nilyr), sice_c(w%nx, w%ny, w%nilyr)
      real(wp), intent(inout) :: qsno_c(w%nx, w%ny, w%nslyr)
      real(wp), intent(inout) :: fresh(w%nx, w%ny), fsalt(w%nx, w%ny), fhocn(w%nx, w%ny)
      integer(int32), intent(inout) :: ice_err(w%nx, w%ny)
      type(transport_result_t), intent(out) :: res

      integer :: nxg, nyg, ncatg, nilg, nslg, ntrg
      integer :: ni, nviol, nbad, nneg, cfg, nbadv
      integer(int32) :: ierr_pack, ierr_cmp
      logical :: check, ok, okb
      real(wp) :: puny

      res%ierr = ICE_ERR_NONE
      res%cfg = TRANSPORT_CFG_OK
      if (.not. w%is_init) then
         res%cfg = TRANSPORT_CFG_WORKSPACE
         return
      end if
      nxg = w%nx
      nyg = w%ny
      ncatg = w%ncat
      nilg = w%nilyr
      nslg = w%nslyr
      ntrg = w%ntrace
      if (p%transport%advection /= ICE_ADVECT_UPWIND .or. .not. ice_transport_params_valid(p%transport) &
          .or. p%ncat /= ncatg .or. p%nilyr /= nilg .or. p%nslyr /= nslg) then
         res%cfg = TRANSPORT_CFG_PARAMS
         return
      end if
      if (.not. upwind_table_ok(w%a)) then
         res%cfg = TRANSPORT_CFG_TABLE
         return
      end if
      if (is < 2 .or. ie > nxg - 1 .or. js < 2 .or. je > nyg - 1 .or. ie < is .or. je < js) then
         res%cfg = TRANSPORT_CFG_RANGE
         return
      end if
      puny = p%c%puny
      check = p%transport%check
      okb = .false.

      ! 1. velocities: finite and within the CFL displacement (always)
      call shift_faces(nxg, nyg, is - 1, ie, js - 1, je, u, v, w%uE, w%vN, nbadv)
      call transport_check_cfl(nxg, nyg, is - 1, ie, js - 1, je, dt, w%uE, w%vN, dxT, dyT, ice_err, nviol)
      res%n_cfl = nviol
      res%n_launch = res%n_launch + 2
      if (nbadv > 0) res%ierr = ior(res%ierr, ICE_ERR_NONFINITE)
      if (nviol > 0 .and. nbadv == 0) res%ierr = ior(res%ierr, ICE_ERR_TRANSPORT_CFL)
      res%n_nonfinite = res%n_nonfinite + nbadv
      if (res%ierr /= ICE_ERR_NONE) return

      ! 2. pack over the owned cells plus one ring; volumes of negligible areas are carried
      call transport_pack(w%a, p%c, nilg, nslg, is - 1, ie + 1, js - 1, je + 1, aice0, aicen, vicen, vsnon, &
                          Tsfcn, qicen, sicen, qsnon, ierr_pack, thick_floor=UPWIND_THICK_FLOOR)
      res%n_launch = res%n_launch + 1
      res%ierr = ior(res%ierr, ierr_pack)

      ! 3. entry checks on the packed fields (per cell: always when the pack saw a non-finite value)
      if (check .or. ierr_pack /= ICE_ERR_NONE) then
         call transport_check_finite(nxg, nyg, ncatg, ntrg, is - 1, ie + 1, js - 1, je + 1, w%a%aim, w%a%trm, ice_err, nviol)
         if (nviol > 0) res%ierr = ior(res%ierr, ICE_ERR_NONFINITE)
         res%n_nonfinite = res%n_nonfinite + nviol
         res%n_launch = res%n_launch + 1
      end if
      if (check) then
         call transport_check_neg_area(nxg, nyg, ncatg, is - 1, ie + 1, js - 1, je + 1, puny, w%a%aim, ice_err, nviol)
         if (nviol > 0) res%ierr = ior(res%ierr, ICE_ERR_TRANSPORT_NEG_AREA)
         res%n_neg = res%n_neg + nviol
         res%n_launch = res%n_launch + 1
      end if
      if (res%ierr /= ICE_ERR_NONE) return

      if (check) then
         call transport_conserved_sums(w%chk, nxg, nyg, ncatg, ntrg, is, ie, js, je, w%a%tab%depend(1:ntrg), tarea, &
                                       w%a%aim, w%a%trm, w%before, ok)
         res%n_launch = res%n_launch + 2
         res%conserv_done = ok
         if (w%chk%has_bounds) then
            call transport_bounds(w%chk, nxg, nyg, ncatg, ntrg, is + 1, ie - 1, js + 1, je - 1, puny, &
                                  w%a%tab%tracer_type(1:ntrg), w%a%tab%depend(1:ntrg), w%a%aim, w%a%trm, okb)
            res%n_launch = res%n_launch + 1
         end if
      end if

      ! 5. advection (one launch)
      call upwind_advect(w%a, w%b, is, ie, js, je, dt, puny, u, v, HTE, HTN, tarea, epm, npm, nbad, nneg, cfg)
      if (cfg /= UPWIND_CFG_OK) then
         res%cfg = merge(TRANSPORT_CFG_TABLE, TRANSPORT_CFG_WORKSPACE, cfg == UPWIND_CFG_DEPEND)
         return
      end if
      res%n_launch = res%n_launch + 1
      res%n_nonfinite = res%n_nonfinite + nbad
      res%n_neg = res%n_neg + nneg
      if (nbad > 0) res%ierr = ior(res%ierr, ICE_ERR_NONFINITE)
      if (nneg > 0) res%ierr = ior(res%ierr, ICE_ERR_TRANSPORT_NEG_AREA)

      ! 6. per-cell bits on the result (always when the advection counted something) and the sums
      if (check .or. nbad > 0 .or. nneg > 0) then
         call transport_check_finite(nxg, nyg, ncatg, ntrg, is, ie, js, je, w%b%aim, w%b%trm, ice_err, nviol)
         if (nviol > 0) res%ierr = ior(res%ierr, ICE_ERR_NONFINITE)
         call transport_check_neg_area(nxg, nyg, ncatg, is, ie, js, je, puny, w%b%aim, ice_err, nviol)
         if (nviol > 0) res%ierr = ior(res%ierr, ICE_ERR_TRANSPORT_NEG_AREA)
         res%n_launch = res%n_launch + 2
      end if
      if (check .and. res%conserv_done .and. res%ierr == ICE_ERR_NONE) then
         call transport_conserved_sums(w%chk, nxg, nyg, ncatg, ntrg, is, ie, js, je, w%b%tab%depend(1:ntrg), tarea, &
                                       w%b%aim, w%b%trm, w%after, ok)
         res%n_launch = res%n_launch + 2
         call transport_conservation_compare_floor(ncatg, ntrg, w%before, w%after, puny, ierr_cmp, res%worst_rel, &
                                                   res%worst_q)
         res%ierr = ior(res%ierr, ierr_cmp)
         if (ierr_cmp /= ICE_ERR_NONE) call mark_all_conserv(ice_err, nxg, nyg, is, js)
      end if
      if (check .and. w%chk%has_bounds .and. okb .and. res%ierr == ICE_ERR_NONE) then
         call transport_check_monotone(w%chk, nxg, nyg, ncatg, ntrg, is + 1, ie - 1, js + 1, je - 1, puny, &
                                       w%b%tab%tracer_type(1:ntrg), w%b%tab%depend(1:ntrg), w%b%aim, w%b%trm, &
                                       ice_err, ni, ok)
         res%n_mono = ni
         res%n_launch = res%n_launch + 1
         if (ni > 0) res%ierr = ior(res%ierr, ICE_ERR_TRANSPORT_CONSERV)
      end if
      if (res%ierr /= ICE_ERR_NONE) return

      ! 8. unpack (owned cells) and the surface temperature of emptied categories
      call transport_unpack(w%b, p%c, nilg, nslg, is, ie, js, je, aice0, aicen, vicen, vsnon, &
                            Tsfcn, qicen, sicen, qsnon, ierr_pack)
      call upwind_empty_tsfc(nxg, nyg, ncatg, is, ie, js, je, puny, Tf, hm, aicen, Tsfcn)
      res%n_launch = res%n_launch + 2
      res%ierr = ior(res%ierr, ierr_pack)
      res%stepped = .true.
      if (ierr_pack /= ICE_ERR_NONE) return

      ! 9. post step
      select case (post)
      case (ICE_TRANSPORT_POST_AGGREGATE)
         call ice_aggregate_state_impl(nxg, nyg, ncatg, is, ie, js, je, aicen, vicen, vsnon, aice, vice, vsno, aice0, ice_err)
         res%n_launch = res%n_launch + 1
      case (ICE_TRANSPORT_POST_CLEANUP)
         call ice_cleanup_itd_step(nxg, nyg, ncatg, nilg, nslg, is, ie, js, je, p, hin_max, dt, Tf, &
                                   .true., .true., .true., aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                   aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ice_err)
         res%n_launch = res%n_launch + 1
      case default
      end select
   end subroutine ice_step_transport

   subroutine mark_all_conserv(ice_err, nx, ny, is, js)
      !! A global conservation failure has no cell: mark the first owned cell so the cell error word carries it.
      integer, intent(in) :: nx, ny, is, js
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j

      i = is
      j = js
      !$acc serial present(ice_err)
      ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_TRANSPORT_CONSERV)
      !$acc end serial
   end subroutine mark_all_conserv

   subroutine transport_conservation_compare_floor(ncat, ntrace, before, after, puny, ierr, worst_rel, worst_q)
      !! Conservation comparison with an absolute floor (host, a few scalars).  Quantity order as
      !! `transport_nquantities`: the first `ncat + 1` are the category areas (`n = 0..ncat`), then
      !! `ncat*ntrace` tracer contents.  With `Aref = sum of the category areas before` (the owned
      !! area, to the extent the fractions sum to 1) the change of quantity `q` must satisfy
      !! `|after - before| <= puny*max(|before|, floor_q)`, `floor_q = Aref` for an area and
      !! `Aref*|before|/|area of its category|` for a tracer content (`Aref` times the category's
      !! mean of that field): a category that holds almost nothing is allowed `puny*Aref*mean`,
      !! which is what dropping the area-tracer content of its sub-`puny` cells costs, and a real
      !! loss, far above that, still trips.  A non-finite total is a violation.  `worst_rel` is
      !! the largest `|change|/max(|before|, floor_q)` (NaN if poisoned), `worst_q` its quantity.
      !! The reference comparison (`transport_conservation_compare`) is unchanged for other users.
      integer, intent(in) :: ncat, ntrace
      type(efp_t), intent(in) :: before(:), after(:)
      real(wp), intent(in) :: puny
      integer(int32), intent(out) :: ierr
      real(wp), intent(out) :: worst_rel
      integer, intent(out) :: worst_q
      integer :: q, nq, n, qq
      real(wp) :: b, a, aref, narea, fl, ref, rel

      nq = transport_nquantities(ncat, ntrace)
      ierr = ICE_ERR_NONE
      worst_rel = 0.0_wp
      worst_q = 0
      if (size(before) < nq .or. size(after) < nq) then
         ierr = ICE_ERR_TRANSPORT_CONSERV
         worst_rel = ieee_value(worst_rel, ieee_quiet_nan)
         return
      end if
      aref = 0.0_wp
      do q = 1, ncat + 1
         aref = aref + efp_to_real(before(q))
      end do
      do q = 1, nq
         b = efp_to_real(before(q))
         a = efp_to_real(after(q))
         if (.not. (ieee_is_finite(b) .and. ieee_is_finite(a))) then
            ierr = ICE_ERR_TRANSPORT_CONSERV
            worst_rel = ieee_value(worst_rel, ieee_quiet_nan)
            worst_q = q
            return
         end if
         if (q <= ncat + 1) then
            fl = max(aref, 0.0_wp)
         else
            qq = q - ncat - 1
            n = (qq - 1)/ntrace + 1
            narea = efp_to_real(before(n + 1))
            fl = 0.0_wp
            if (narea > 0.0_wp) fl = max(aref, 0.0_wp)*abs(b)/narea
         end if
         ref = max(abs(b), fl)
         if (ref > 0.0_wp) then
            rel = abs(a - b)/ref
            if (rel > worst_rel) then
               worst_rel = rel
               worst_q = q
            end if
            if (rel > puny) ierr = ICE_ERR_TRANSPORT_CONSERV
         else if (a /= b) then
            ierr = ICE_ERR_TRANSPORT_CONSERV
         end if
      end do
   end subroutine transport_conservation_compare_floor

   subroutine shift_faces(nx, ny, is, ie, js, je, u, v, uE, vN, nbad)
      !! `uE(i,j) = u(i+1,j)`, `vN(i,j) = v(i,j+1)` over `is:ie, js:je` (one launch):
      !! the T-cell indexing the CFL check wants.  `nbad` counts non-finite faces.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(inout) :: uE(nx, ny), vN(nx, ny)
      integer, intent(out) :: nbad
      integer :: i, j

      nbad = 0
      !$acc data present(u, v, uE, vN)
      do concurrent(j=js:je, i=is:ie) reduce(+:nbad)
         uE(i, j) = u(i + 1, j)
         vN(i, j) = v(i, j + 1)
         if (.not. ieee_is_finite(uE(i, j))) nbad = nbad + 1
         if (.not. ieee_is_finite(vN(i, j))) nbad = nbad + 1
      end do
      !$acc end data
   end subroutine shift_faces

end module rdb_ice_transport_cice
