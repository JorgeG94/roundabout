!! Global-integral conservation-budget PRIMITIVES — a verification
!! instrument, not a production reporter.
module rdb_ocean_budgets
   !! Globally-integrated conservation scalars (total mass, KE, salt,
   !! heat) plus a per-kernel contributor registry, used to check a
   !! kernel's per-cell budget accumulator against the actual state
   !! change (LHS = current_total − values_init vs RHS =
   !! Σ contributors%total_integrated).  **No production caller by
   !! design (PR-8):** this module has no `ocean_state_t` slot — a
   !! consumer constructs a local `type(ocean_budgets_t)` (seven of
   !! Roundabout's ocean conservation tests already do this; see
   !! `tests/test_ocean_surface_flux.F90` for the idiom). It is
   !! single-rank with no MPI reduction; a production run's global
   !! conservation report is `rdb_ocean_console_stats.F90:
   !! ocean_console_stats_report`, which does allreduce correctly.
   !! `init_snapshot(state)` captures `values_init`;
   !! `evaluate(state)` recomputes `values`.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_diag_mask, only: diag_mask_t, diag_mask_global, diag_mask_destroy
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_budgets_t, budget_contributor_t
   public :: budget_total_mass, budget_total_tracer, budget_total_ke

   ! Integral tags — what `evaluate` recomputes.
   integer, parameter, public :: BUDGET_MASS = 1
   integer, parameter, public :: BUDGET_KE = 2  !! kinetic energy
   integer, parameter, public :: BUDGET_SALT_TOTAL = 3
   integer, parameter, public :: BUDGET_HEAT_TOTAL = 4
   integer, parameter, public :: BUDGET_N_DEFAULT = 4

   character(len=12), parameter :: BUDGET_NAMES(BUDGET_N_DEFAULT) = [ &
                                   "mass        ", "KE          ", &
                                   "salt_total  ", "heat_total  "]

   integer, parameter :: INITIAL_CONTRIBUTOR_CAPACITY = 16

   type :: budget_contributor_t
      !! One source/sink of a conserved quantity, populated by a physics
      !! kernel each step (adds `dt·delta` into `per_cell`) and drained at
      !! eval cadence (`drain_contributors` integrates `per_cell·mask·dA`
      !! over k into `total_integrated`, then zeroes `per_cell`).
      !! `total_integrated` is the cumulative RHS checked against
      !! LHS = current_total − values_init.
      character(len=64) :: name = ""
         !! Identifier — "continuity", "vert_diff_heat", "surface_heat_flux", ...
      integer :: quantity = 0
         !! Tag from BUDGET_* (MASS / SALT_TOTAL / HEAT_TOTAL / ...).
      real(wp), pointer :: per_cell(:, :, :) => null()
         !! Kernel-owned per-step accumulator; manager holds a non-owning
         !! pointer (kernel guarantees the array outlives registration).
      real(wp) :: total_integrated = 0.0_wp
         !! Cumulative spatial integral of contributions drained so
         !! far this run.  Sign convention: positive = source.
      logical :: is_active = .false.
      logical :: device_resident = .false.
         !! True when `per_cell` lives on the GPU (mapped via OpenACC by the
         !! owning kernel); `drain_contributors` then `acc update self`
         !! before integrating and `acc update device` after zeroing.
         !! Default false = host-only kernels / synthetic tests.
   end type budget_contributor_t

   type :: ocean_budgets_t
      logical :: is_init = .false.

      ! ---- Active integrals ----
      integer :: nbudgets = BUDGET_N_DEFAULT
      real(wp) :: values(BUDGET_N_DEFAULT) = 0.0_wp
      real(wp) :: values_init(BUDGET_N_DEFAULT) = 0.0_wp
      logical :: has_snapshot = .false.
         !! True once `init_snapshot` has populated `values_init`.

      ! ---- Geometry cache ----
      type(hgrid_t) :: grid
      real(wp), allocatable :: areaT(:, :)
         !! Per-cell T-area (m²) copied from the metrics slot via `set_area`.
         !! Physical integrals weight by `mask%weight·areaT`; falls back to
         !! `grid%dx·grid%dy` until set.

      ! ---- Region restriction ----
      type(diag_mask_t), allocatable :: mask
         !! Optional region mask.  Default = whole interior (built lazily);
         !! reset via `set_mask`.

      ! ---- Per-kernel contributors ----
      ! Populated via `register_contributor`, drained each eval cadence.
      ! Budget identity: LHS = values(q) − values_init(q),
      ! RHS = Σ contributors of quantity q %total_integrated, residual = LHS − RHS.
      type(budget_contributor_t), allocatable :: contributors(:)
      integer :: n_contributors = 0
      integer :: n_contributors_max = 0
   contains
      procedure, non_overridable :: init => ocean_budgets_init
      procedure, non_overridable :: destroy => ocean_budgets_destroy
      procedure, non_overridable :: set_area => ocean_budgets_set_area
      procedure, non_overridable :: set_mask => ocean_budgets_set_mask
      procedure, non_overridable :: init_snapshot => ocean_budgets_init_snapshot
      procedure, non_overridable :: evaluate => ocean_budgets_evaluate
      procedure, non_overridable :: register_contributor => ocean_budgets_register_contributor
      procedure, non_overridable :: drain_contributors => ocean_budgets_drain_contributors
      procedure, non_overridable :: bytes => ocean_budgets_bytes
   end type ocean_budgets_t

contains

   subroutine ocean_budgets_init(this, grid)
      class(ocean_budgets_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      this%grid = grid
      this%has_snapshot = .false.
      this%values = 0.0_wp
      this%values_init = 0.0_wp
      this%n_contributors = 0
      this%n_contributors_max = INITIAL_CONTRIBUTOR_CAPACITY
      allocate (this%contributors(this%n_contributors_max))
      this%is_init = .true.
   end subroutine ocean_budgets_init

   subroutine ocean_budgets_destroy(this)
      class(ocean_budgets_t), intent(inout) :: this
      integer :: i
      this%is_init = .false.
      this%nbudgets = 0
      this%values = 0.0_wp
      this%values_init = 0.0_wp
      this%has_snapshot = .false.
      if (allocated(this%areaT)) deallocate (this%areaT)
      if (allocated(this%mask)) then
         call diag_mask_destroy(this%mask)
         deallocate (this%mask)
      end if
      if (allocated(this%contributors)) then
         do i = 1, this%n_contributors
            nullify (this%contributors(i)%per_cell)
         end do
         deallocate (this%contributors)
      end if
      this%n_contributors = 0
      this%n_contributors_max = 0
   end subroutine ocean_budgets_destroy

   subroutine ocean_budgets_set_area(this, areaT)
      !! Cache the per-cell T-area (m²) from the metrics slot.  After this,
      !! physical integrals weight by `mask%weight·areaT`; before it they
      !! fall back to `grid%dx·grid%dy` (= areaT on uniform Cartesian).
      class(ocean_budgets_t), intent(inout) :: this
      real(wp), intent(in) :: areaT(:, :)
      if (allocated(this%areaT)) deallocate (this%areaT)
      allocate (this%areaT, source=areaT)
   end subroutine ocean_budgets_set_area

   subroutine ocean_budgets_set_mask(this, mask)
      !! Override the default global mask.  Useful for regional
      !! conservation diagnostics — "mass north of 30°S stays at FP."
      class(ocean_budgets_t), intent(inout) :: this
      type(diag_mask_t), intent(in) :: mask
      if (allocated(this%mask)) then
         call diag_mask_destroy(this%mask)
         deallocate (this%mask)
      end if
      allocate (this%mask, source=mask)
   end subroutine ocean_budgets_set_mask

   subroutine ocean_budgets_init_snapshot(this, ms)
      !! Snapshot `values_init` from the current state.  Call once
      !! after the IC is set + the EOS has run, before stepping.
      class(ocean_budgets_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      if (.not. allocated(this%mask)) then
         allocate (this%mask, source=diag_mask_global(this%grid))
      end if
      call this%evaluate(ms)
      this%values_init = this%values
      this%has_snapshot = .true.
   end subroutine ocean_budgets_init_snapshot

   subroutine ocean_budgets_evaluate(this, ms)
      !! Recompute `this%values(:)` from the current multilayer state.
      !! Honours `this%mask` — sums are weighted by `mask%weight · areaT`.
      class(ocean_budgets_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp) :: total_mass, total_KE, total_S, total_T
      real(wp), allocatable :: areaT(:, :)
      integer :: it_S, it_T

      if (.not. allocated(this%mask)) then
         allocate (this%mask, source=diag_mask_global(this%grid))
      end if

      call budget_area_field(this, areaT)

      total_mass = budget_total_mass(ms, this%mask, areaT)
      total_KE = budget_total_ke(ms, this%mask, areaT)

      it_S = ms%idx_salinity
      it_T = ms%idx_temperature
      total_S = 0.0_wp
      total_T = 0.0_wp
      if (it_S > 0) total_S = budget_total_tracer(ms, it_S, this%mask, areaT)
      if (it_T > 0) total_T = budget_total_tracer(ms, it_T, this%mask, areaT)

      this%values(BUDGET_MASS) = total_mass
      this%values(BUDGET_KE) = total_KE
      this%values(BUDGET_SALT_TOTAL) = total_S
      this%values(BUDGET_HEAT_TOTAL) = total_T
   end subroutine ocean_budgets_evaluate

   subroutine ocean_budgets_register_contributor(this, name, quantity, per_cell, device_resident)
      !! Register a per-kernel contributor: store a non-owning pointer to
      !! the kernel-owned `per_cell` + metadata.  Set `device_resident`
      !! when `per_cell` is GPU-mapped (drain syncs it).  Idempotent on
      !! `name` (re-registering updates the pointer, resets the integral).
      class(ocean_budgets_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: quantity
      real(wp), intent(in), target :: per_cell(:, :, :)
      logical, intent(in), optional :: device_resident
      type(budget_contributor_t), allocatable :: tmp(:)
      integer :: i, slot

      if (.not. this%is_init) return

      slot = 0
      do i = 1, this%n_contributors
         if (trim(this%contributors(i)%name) == trim(name)) then
            slot = i
            exit
         end if
      end do

      if (slot == 0) then
         if (this%n_contributors == this%n_contributors_max) then
            allocate (tmp(2*this%n_contributors_max))
            do i = 1, this%n_contributors
               tmp(i) = this%contributors(i)
            end do
            call move_alloc(tmp, this%contributors)
            this%n_contributors_max = size(this%contributors)
         end if
         this%n_contributors = this%n_contributors + 1
         slot = this%n_contributors
      end if

      associate (c => this%contributors(slot))
         c%name = name
         c%quantity = quantity
         c%per_cell => per_cell
         c%total_integrated = 0.0_wp
         c%is_active = .true.
         c%device_resident = .false.
         if (present(device_resident)) c%device_resident = device_resident
      end associate
   end subroutine ocean_budgets_register_contributor

   subroutine ocean_budgets_drain_contributors(this)
      !! For each contributor: integrate `per_cell·mask·areaT` (over k) into
      !! `total_integrated`, then zero `per_cell` for the next window.
      class(ocean_budgets_t), intent(inout) :: this
      integer :: ic, i, j, k, nx, ny, nz
      real(wp) :: w_dA, col_sum, increment
      real(wp), allocatable :: areaT(:, :)

      if (.not. this%is_init) return
      if (.not. allocated(this%mask)) then
         allocate (this%mask, source=diag_mask_global(this%grid))
      end if
      call budget_area_field(this, areaT)

      do ic = 1, this%n_contributors
         associate (c => this%contributors(ic))
            if (.not. c%is_active) cycle
            if (.not. associated(c%per_cell)) cycle
            if (c%device_resident) then
               !$acc update self(c%per_cell) if_present
            end if
            nx = min(this%mask%nx, size(c%per_cell, 1))
            ny = min(this%mask%ny, size(c%per_cell, 2))
            nz = size(c%per_cell, 3)
            increment = 0.0_wp
            do j = 1, ny
               do i = 1, nx
                  w_dA = this%mask%weight(i, j)*areaT(i, j)
                  if (w_dA <= 0.0_wp) cycle
                  col_sum = 0.0_wp
                  do k = 1, nz
                     col_sum = col_sum + c%per_cell(i, j, k)
                  end do
                  increment = increment + w_dA*col_sum
               end do
            end do
            c%total_integrated = c%total_integrated + increment
            c%per_cell = 0.0_wp
            if (c%device_resident) then
               !$acc update device(c%per_cell) if_present
            end if
         end associate
      end do
   end subroutine ocean_budgets_drain_contributors

   ! ---------------------------------------------------------------------
   ! Module-level integral helpers — public so tests can call them directly.
   ! ---------------------------------------------------------------------

   subroutine budget_area_field(this, areaT)
      !! Per-cell T-area (m²) for physical integrals: cached metrics `areaT`
      !! once `set_area` has run, else uniform `grid%dx·grid%dy`.
      type(ocean_budgets_t), intent(in) :: this
      real(wp), allocatable, intent(out) :: areaT(:, :)
      if (allocated(this%areaT)) then
         areaT = this%areaT
      else
         allocate (areaT(this%grid%nx_total, this%grid%ny_total), &
                   source=this%grid%dx*this%grid%dy)
      end if
   end subroutine budget_area_field

   pure function budget_total_mass(ms, mask, areaT) result(total)
      !! Total mass: Σ over masked interior of `h_layer·areaT·weight` over k.
      !! Units: m³.  Public only for the unit-test suite.
      type(multilayer_state_t), intent(in) :: ms
      type(diag_mask_t), intent(in) :: mask
      real(wp), intent(in) :: areaT(:, :)
      real(wp) :: total
      integer :: i, j, k, nx, ny, nz
      real(wp) :: w_dA, col_sum
      total = 0.0_wp
      nx = min(mask%nx, size(ms%h_layer, 1))
      ny = min(mask%ny, size(ms%h_layer, 2))
      nz = ms%nz_ml
      do j = 1, ny
         do i = 1, nx
            w_dA = mask%weight(i, j)*areaT(i, j)
            if (w_dA <= 0.0_wp) cycle
            col_sum = 0.0_wp
            do k = 1, nz
               col_sum = col_sum + ms%h_layer(i, j, k)
            end do
            total = total + w_dA*col_sum
         end do
      end do
   end function budget_total_mass

   pure function budget_total_tracer(ms, it, mask, areaT) result(total)
      !! Total tracer content: Σ over masked interior of `hTr·areaT·weight`
      !! over k.  `hTr` is thickness-weighted, so this is ∫(tracer·volume),
      !! the conserved quantity.  Public only for the unit-test suite.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: it
      type(diag_mask_t), intent(in) :: mask
      real(wp), intent(in) :: areaT(:, :)
      real(wp) :: total
      integer :: i, j, k, nx, ny, nz
      real(wp) :: w_dA, col_sum
      total = 0.0_wp
      if (it <= 0) return
      if (.not. allocated(ms%tracers)) return
      if (it > size(ms%tracers)) return
      if (.not. allocated(ms%tracers(it)%hTr)) return
      nx = min(mask%nx, size(ms%tracers(it)%hTr, 1))
      ny = min(mask%ny, size(ms%tracers(it)%hTr, 2))
      nz = ms%nz_ml
      do j = 1, ny
         do i = 1, nx
            w_dA = mask%weight(i, j)*areaT(i, j)
            if (w_dA <= 0.0_wp) cycle
            col_sum = 0.0_wp
            do k = 1, nz
               col_sum = col_sum + ms%tracers(it)%hTr(i, j, k)
            end do
            total = total + w_dA*col_sum
         end do
      end do
   end function budget_total_tracer

   pure function budget_total_ke(ms, mask, areaT) result(total)
      !! Total KE: Σ over masked interior of `0.5·h·(u²+v²)·areaT·weight`,
      !! u/v averaged from C-grid faces to centres.  Public only for tests.
      type(multilayer_state_t), intent(in) :: ms
      type(diag_mask_t), intent(in) :: mask
      real(wp), intent(in) :: areaT(:, :)
      real(wp) :: total
      integer :: i, j, k, nx, ny, nz
      real(wp) :: w_dA, uc, vc, ke_col
      total = 0.0_wp
      nx = min(mask%nx, size(ms%h_layer, 1), &
               size(ms%u_face_x_layer, 1) - 1, &
               size(ms%v_face_y_layer, 1))
      ny = min(mask%ny, size(ms%h_layer, 2), &
               size(ms%u_face_x_layer, 2), &
               size(ms%v_face_y_layer, 2) - 1)
      nz = ms%nz_ml
      do j = 1, ny
         do i = 1, nx
            w_dA = mask%weight(i, j)*areaT(i, j)
            if (w_dA <= 0.0_wp) cycle
            ke_col = 0.0_wp
            do k = 1, nz
               uc = 0.5_wp*(ms%u_face_x_layer(i, j, k) + ms%u_face_x_layer(i + 1, j, k))
               vc = 0.5_wp*(ms%v_face_y_layer(i, j, k) + ms%v_face_y_layer(i, j + 1, k))
               ke_col = ke_col + ms%h_layer(i, j, k)*0.5_wp*(uc*uc + vc*vc)
            end do
            total = total + w_dA*ke_col
         end do
      end do
   end function budget_total_ke

   pure function ocean_budgets_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the conservation budgets (own accumulator; registry contributors are counted by their owning slot) slot (0 when
      !! unallocated).
      class(ocean_budgets_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%areaT)
      ! Region mask — HOST-only (lazily allocated by `init_snapshot` /
      ! `evaluate`, never `!$acc enter data`'d).  Counted anyway now that
      ! `diag_mask_t` has a `bytes()`: the term is 0 on every shipped
      ! namelist (`ocean_budgets_t` has no production consumer) and keeping
      ! it makes the drift hook's "every countable component" rule hold.
      if (allocated(this%mask)) nbytes = nbytes + this%mask%bytes()
   end function ocean_budgets_bytes

end module rdb_ocean_budgets
