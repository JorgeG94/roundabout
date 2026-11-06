!! Vertical tracer advection on the ocean multilayer C-grid.
!!
!! Two-piece kernel:
!!
!!   1. `compute_w_from_continuity` diagnoses
!!      `ms%w_interface` from the horizontal-continuity residual
!!      `ms%flux_h_layer`.  Integrates from the bed upward under the
!!      no-flow boundary condition `w_interface(:, :, 1) = 0`.  The
!!      surface value `w_interface(:, :, nz+1)` then represents the
!!      free-surface displacement rate (Eulerian z).
!!
!!   2. `tracer_advect_vertical` consumes `w_interface`
!!      and applies a first-order upwind-in-z vertical advection
!!      step to every registered tracer.  Same outer-shim + flat-
!!      impl pattern as horizontal tracer advection, so the registry
!!      indirection only appears in the wrapper.
!!
!! Convention (ROMS-style, bottom-up):
!!   * `w_interface(:, :, 1)`     — bed BC, forced to zero
!!   * `w_interface(:, :, k)`     — interface between layer k-1
!!                                  (below) and layer k (above)
!!   * `w_interface(:, :, nz+1)`  — surface, equals the diagnosed
!!                                  free-surface displacement rate
!!
!! Tracer flux at interface k = `w_interface(k) * Tr_face_upwind(k)`.
!! Layer k's hTr update:
!!   ∂(hTr)/∂t = w_interface(k)   * Tr_face(k)
!!             - w_interface(k+1) * Tr_face(k+1)
!! (flux in at bottom minus flux out at top, positive w = upward).
!!
!! Per-tracer gating: respects `tracer_t%do_vertical_exchange`
!! (default true).  Forced / diagnostic tracers can opt out at the
!! registry level.
!!
!! Phase Tier-1 uses first-order upwind for robustness; Phase 5+
!! may switch to PPM-in-z once vertical resolution warrants the
!! higher-order treatment.
module rdb_ocean_vertical_advection
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_vertical_advection_t
   public :: compute_w_from_continuity
   public :: tracer_advect_vertical

   type :: ocean_vertical_advection_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: enforce_bed_bc = .true.
         !! When true the bed interface w is forced to zero before
         !! integration regardless of what `compute_w_from_continuity`
         !! is handed.  Disable only when the caller wants to inject
         !! a sub-bed flow source (test instrumentation).
      type(scratch_3d_buffer_t) :: F_face
         !! Per-step vertical tracer-flux at interfaces.  Shape
         !! (nx, ny, nz_ml+1).  Two-pass design: tracer impl first
         !! fills `F_face` from `hTr` + `w_interface` (read-only),
         !! then applies the divergence onto `hTr`.  This avoids
         !! the do-concurrent race that a single-pass version would
         !! have (each layer reads its neighbour's hTr while
         !! another iteration writes it).  Reused across all
         !! tracers in one call.
   contains
      procedure, non_overridable :: init => ocean_vert_adv_init
      procedure, non_overridable :: destroy => ocean_vert_adv_destroy
      procedure, non_overridable :: enter_data => ocean_vert_adv_enter_data
      procedure, non_overridable :: exit_data => ocean_vert_adv_exit_data
      procedure, non_overridable :: bytes => ocean_vertical_advection_bytes
   end type ocean_vertical_advection_t

contains

   subroutine ocean_vert_adv_init(this, grid, nz_ml)
      class(ocean_vertical_advection_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      call this%F_face%init(nx, ny, nz + 1, "ocean_vert_adv_F_face")
      this%is_init = .true.
   end subroutine ocean_vert_adv_init

   subroutine ocean_vert_adv_destroy(this)
      class(ocean_vertical_advection_t), intent(inout) :: this
      this%is_init = .false.
      call this%F_face%destroy()
   end subroutine ocean_vert_adv_destroy

   subroutine ocean_vert_adv_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! device-attach map base is the heap object, not a polymorphic stack
      !! box (AMD libomptarget cross-slot-overlap fix).
      class(ocean_vertical_advection_t), intent(inout) :: this
      select type (this)
      type is (ocean_vertical_advection_t)
         call ocean_vert_adv_enter_data_impl(this)
      end select
   end subroutine ocean_vert_adv_enter_data

   subroutine ocean_vert_adv_enter_data_impl(this)
      type(ocean_vertical_advection_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%F_face)
   end subroutine ocean_vert_adv_enter_data_impl

   subroutine ocean_vert_adv_exit_data(this)
      class(ocean_vertical_advection_t), intent(inout) :: this
      select type (this)
      type is (ocean_vertical_advection_t)
         call ocean_vert_adv_exit_data_impl(this)
      end select
   end subroutine ocean_vert_adv_exit_data

   subroutine ocean_vert_adv_exit_data_impl(this)
      type(ocean_vertical_advection_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%F_face)
   end subroutine ocean_vert_adv_exit_data_impl

   pure subroutine compute_w_from_continuity(grid, this, ms)
      !! Fill `ms%w_interface` by integrating the horizontal-
      !! continuity residual upward from the bed.  Eulerian z:
      !!
      !!   w_interface(:, :, 1)    = 0  (bed BC)
      !!   w_interface(:, :, k+1)  = w_interface(:, :, k) - flux_h_layer(k)
      !!
      !! With this w, ∂h/∂t = -horizontal_div + (w(k) - w(k+1)) = 0,
      !! so layer thicknesses stay at their initial Eulerian z
      !! positions.  `continuity_compute_fluxes` must run
      !! first — `flux_h_layer` is the input.
      !!
      !! Recurrence is data-dependent in k, so the outer loop is
      !! serial in k inside each column.  Columns are parallelised
      !! over `(j, i)`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vertical_advection_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms

      integer :: i, j, k, nx, ny, nz
      logical :: enforce_bed

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      enforce_bed = this%enforce_bed_bc

      do concurrent(j=1:ny, i=1:nx)
         if (enforce_bed) then
            ms%w_interface(i, j, 1) = 0.0_wp
         end if
         do k = 1, nz
            ms%w_interface(i, j, k + 1) = ms%w_interface(i, j, k) - &
                                          ms%flux_h_layer(i, j, k)
         end do
      end do
   end subroutine compute_w_from_continuity

   subroutine tracer_advect_vertical(grid, this, ms, dt, active)
      !! Apply first-order upwind-in-z vertical advection to every
      !! registered tracer, then update `h_layer` by the same
      !! vertical mass-flux divergence.  Updating h is what makes
      !! the kernel CWC-consistent — uniform `T = hTr/h` stays
      !! uniform regardless of how divergent the w field is.
      !!
      !! Order matters: tracers read the pre-update `h_layer` for
      !! their upwind concentration; `h_layer` is updated only
      !! after all tracers have run.
      !!
      !! Eulerian-z mode: the h update here exactly cancels the
      !! horizontal h update produced by
      !! `continuity_tracer_step_split` (both consume
      !! the same `flux_h_layer` — the total horizontal divergence
      !! summed over both substeps), so `h_layer` stays at its
      !! initial z-coordinate values.
      !!
      !! Per-tracer `do_vertical_exchange` gate lets diagnostic /
      !! forced tracers opt out without affecting horizontal
      !! advection.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vertical_advection_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional gate (thermo cadence).  Absent ⇒ kernel runs;
         !! present-and-false ⇒ early return.

      integer :: it

      if (present(active)) then
         if (.not. active) return
      end if
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. ms%tracers(it)%do_vertical_exchange) cycle
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call tracer_advect_vertical_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, dt, &
                  ms%h_layer, ms%w_interface, ms%tracers(it)%hTr, &
                  this%F_face%data, budget=ms%heat_budget_vert_adv)
            case (TRACER_BUDGET_SALT)
               call tracer_advect_vertical_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, dt, &
                  ms%h_layer, ms%w_interface, ms%tracers(it)%hTr, &
                  this%F_face%data, budget=ms%salt_budget_vert_adv)
            case default
               call tracer_advect_vertical_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, dt, &
                  ms%h_layer, ms%w_interface, ms%tracers(it)%hTr, &
                  this%F_face%data)
            end select
         end do
      end if

      call apply_w_to_h_layer(grid%nx_total, grid%ny_total, ms%nz_ml, dt, &
                              ms%w_interface, ms%h_layer)
   end subroutine tracer_advect_vertical

   pure subroutine apply_w_to_h_layer(nx, ny, nz, dt, w_interface, h_layer)
      !! `h_layer(k) += dt * (w(k) - w(k+1))` per cell.  Bed and
      !! surface interfaces feed through whatever w the caller set
      !! (zero for the bed by default after
      !! `compute_w_from_continuity` with `enforce_bed_bc = true`).
      !! No h-floor — plain accumulation; the caller is responsible
      !! for guarding against negative thickness if the prescribed
      !! w + dt is large.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: w_interface(nx, ny, nz + 1)
      real(wp), intent(inout) :: h_layer(nx, ny, nz)

      integer :: i, j, k

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         h_layer(i, j, k) = h_layer(i, j, k) + &
                            dt*(w_interface(i, j, k) - w_interface(i, j, k + 1))
      end do
   end subroutine apply_w_to_h_layer

   pure subroutine tracer_advect_vertical_one_impl(nx, ny, nz, dt, &
                                                   h, w_interface, hTr, F_face, budget)
      !! Flat-impl first-order upwind-in-z vertical advection for
      !! one tracer.  Two passes:
      !!
      !!   Pass 1: fill `F_face(:, :, k)` for k = 1..nz+1 with the
      !!     upwind tracer flux at interface k.  Bed (k=1) and
      !!     surface (k=nz+1) faces get zero.  Reads `hTr` only.
      !!
      !!   Pass 2: `hTr(k) += dt * (F_face(k) - F_face(k+1))`.
      !!     Reads `F_face` only.
      !!
      !! Two-pass avoids the do-concurrent race a single-pass
      !! version would have (each layer reads its neighbour's hTr
      !! while another iteration writes it).  The shared `F_face`
      !! scratch lives on `ocean_vertical_advection_t` and is
      !! reused across all tracers in one call.
      !!
      !! Vanishing-layer donor side: if h_donor <= 0, `T_donor`
      !! falls back to zero (no flux).
      !!
      !! `budget` (optional): Phase D v2 contributor slot.  When
      !! present, the per-cell `dt · (F(k) - F(k+1))` increment is
      !! also accumulated into the slot (alongside the hTr update)
      !! for budget closure.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: w_interface(nx, ny, nz + 1)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: F_face(nx, ny, nz + 1)
      real(wp), intent(inout), optional :: budget(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: w_at_k, T_donor, h_donor

      ! ---- Pass 1: per-interface upwind tracer flux ----
      do concurrent(k=2:nz, j=1:ny, i=1:nx) &
         local(w_at_k, T_donor, h_donor)
         w_at_k = w_interface(i, j, k)
         if (w_at_k >= 0.0_wp) then
            h_donor = h(i, j, k - 1)
            if (h_donor > 0.0_wp) then
               T_donor = hTr(i, j, k - 1)/h_donor
            else
               T_donor = 0.0_wp
            end if
         else
            h_donor = h(i, j, k)
            if (h_donor > 0.0_wp) then
               T_donor = hTr(i, j, k)/h_donor
            else
               T_donor = 0.0_wp
            end if
         end if
         F_face(i, j, k) = w_at_k*T_donor
      end do
      ! Bed and surface fluxes: zero by closed BC.
      do concurrent(j=1:ny, i=1:nx)
         F_face(i, j, 1) = 0.0_wp
         F_face(i, j, nz + 1) = 0.0_wp
      end do

      ! ---- Pass 2: apply flux divergence ----
      if (present(budget)) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            hTr(i, j, k) = hTr(i, j, k) + &
                           dt*(F_face(i, j, k) - F_face(i, j, k + 1))
            budget(i, j, k) = budget(i, j, k) + &
                              dt*(F_face(i, j, k) - F_face(i, j, k + 1))
         end do
      else
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            hTr(i, j, k) = hTr(i, j, k) + &
                           dt*(F_face(i, j, k) - F_face(i, j, k + 1))
         end do
      end if
   end subroutine tracer_advect_vertical_one_impl

   pure function ocean_vertical_advection_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the vertical advection slot (0 when
      !! unallocated).
      class(ocean_vertical_advection_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%F_face%bytes()
   end function ocean_vertical_advection_bytes

end module rdb_ocean_vertical_advection
