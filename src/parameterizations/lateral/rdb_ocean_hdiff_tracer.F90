!! Constant-coefficient horizontal tracer diffusion for the ocean
!! multilayer C-grid.  Conservative curvilinear flux-form Laplacian
!! on `T = hTr/h` (design §2, mirrors continuity):
!!
!!   F_x(i,j,k) = kappa * h_face_x * (T(i,j,k) - T(i-1,j,k)) * idxCu * dy_cu
!!   hTr(i,j,k) += dt * [(F_x(i+1)-F_x(i)) + (F_y(j+1)-F_y(j))] * iareaT
!!
!! The face flux carries the width-weighted transport [tracer·m³/s]
!! and the divergence is area-normalised by `iareaT`; on uniform
!! Cartesian this collapses to the `/dx`/`/dy` form to round-off.
!! `h_face_x = 0.5 * (h(i-1, j, k) + h(i, j, k))` keeps the kernel
!! conservative under variable thickness — the divergence of an
!! h-weighted flux integrates to zero over the closed-wall domain.
!!
!! Wall faces force `F = 0` (closed-wall, no tracer flux across the
!! boundary).  Vanishing-layer guards on the `T = hTr/h` division
!! mirror the horizontal advection kernel.
!!
!! Phase Tier-1 ships the constant-`kappa_h` variant; Phase 5+ may
!! drop in a Smagorinsky-style isotropic closure or a Redi-style
!! isopycnal rotation.
module rdb_ocean_hdiff_tracer
   !! Kernel state for the per-tracer horizontal-diffusion sweep.
   !! Sits alongside `rdb_ocean_horizontal_viscosity` (which does
   !! the momentum half of horizontal mixing); this module is the
   !! tracer-side counterpart.  Same compute-then-apply pattern,
   !! same scratch-buffer two-pass design for race-free
   !! `do concurrent` parallelism.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, OBC_WALL
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_hdiff_tracer_t
   public :: tracer_hdiff

   type :: ocean_hdiff_tracer_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      real(wp) :: kappa_h = 0.0_wp
         !! Constant horizontal tracer diffusivity (m^2/s).  Zero
         !! is a no-op — the kernel short-circuits without
         !! touching hTr.  Stability bound (explicit forward-Euler):
         !!   kappa_h * dt * (1/dx^2 + 1/dy^2) <= 0.5

      type(scratch_3d_buffer_t) :: T_centre
         !! `T = hTr/h` at cell centres.  Shape (nx, ny, nz_ml).
         !! Filled per tracer at the start of each impl call;
         !! reused across the F_x / F_y face passes so the
         !! divergence step reads a consistent snapshot.
      type(scratch_3d_buffer_t) :: F_x_face
         !! East-face tracer flux `kappa * h_face * dT/dx`.
         !! Shape (nx+1, ny, nz_ml).
      type(scratch_3d_buffer_t) :: F_y_face
         !! North-face counterpart.  Shape (nx, ny+1, nz_ml).
   contains
      procedure, non_overridable :: init => ocean_hdiff_tracer_init
      procedure, non_overridable :: destroy => ocean_hdiff_tracer_destroy
      procedure, non_overridable :: enter_data => ocean_hdiff_tracer_enter_data
      procedure, non_overridable :: exit_data => ocean_hdiff_tracer_exit_data
      procedure, non_overridable :: bytes => ocean_hdiff_tracer_bytes
   end type ocean_hdiff_tracer_t

contains

   subroutine ocean_hdiff_tracer_init(this, grid, nz_ml)
      class(ocean_hdiff_tracer_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      call this%T_centre%init(nx, ny, nz, "ocean_hdiff_T_centre")
      call this%F_x_face%init(nx + 1, ny, nz, "ocean_hdiff_F_x_face")
      call this%F_y_face%init(nx, ny + 1, nz, "ocean_hdiff_F_y_face")
      this%is_init = .true.
   end subroutine ocean_hdiff_tracer_init

   subroutine ocean_hdiff_tracer_destroy(this)
      class(ocean_hdiff_tracer_t), intent(inout) :: this
      this%is_init = .false.
      call this%T_centre%destroy()
      call this%F_x_face%destroy()
      call this%F_y_face%destroy()
   end subroutine ocean_hdiff_tracer_destroy

   subroutine ocean_hdiff_tracer_enter_data(this)
      class(ocean_hdiff_tracer_t), intent(inout) :: this
      select type (this)
      type is (ocean_hdiff_tracer_t)
         call ocean_hdiff_tracer_enter_data_impl(this)
      end select
   end subroutine ocean_hdiff_tracer_enter_data

   subroutine ocean_hdiff_tracer_enter_data_impl(this)
      type(ocean_hdiff_tracer_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%T_centre)
      call scratch_3d_buffer_enter_data_impl(this%F_x_face)
      call scratch_3d_buffer_enter_data_impl(this%F_y_face)
   end subroutine ocean_hdiff_tracer_enter_data_impl

   subroutine ocean_hdiff_tracer_exit_data(this)
      class(ocean_hdiff_tracer_t), intent(inout) :: this
      select type (this)
      type is (ocean_hdiff_tracer_t)
         call ocean_hdiff_tracer_exit_data_impl(this)
      end select
   end subroutine ocean_hdiff_tracer_exit_data

   subroutine ocean_hdiff_tracer_exit_data_impl(this)
      type(ocean_hdiff_tracer_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%T_centre)
      call scratch_3d_buffer_exit_data_impl(this%F_x_face)
      call scratch_3d_buffer_exit_data_impl(this%F_y_face)
   end subroutine ocean_hdiff_tracer_exit_data_impl

   subroutine tracer_hdiff(grid, metrics, this, ms, dt, active, bc)
      !! Iterate the tracer registry and apply constant-`kappa_h`
      !! horizontal Laplacian diffusion to every tracer whose
      !! `do_horizontal_diffusion` flag is set.  Outer-shim:
      !! forwards each tracer's hTr to the flat-impl below.
      !! Reuses the same scratch buffers across all tracers in one
      !! call.
      !!
      !! Optional `active` — when present and false the kernel is a
      !! no-op (used by the dyn step to gate on
      !! thermodynamics + thermo-substep cadence).
      !!
      !! Optional `bc` — per-edge OBC tags, resolved to plain host
      !! logicals before the kernel call (mirrors
      !! `continuity_zonal_flux`; a derived-type dummy must never
      !! reach a `do concurrent` kernel — mem:separate).  Absent =>
      !! WALL on every edge (single-rank default).  A non-WALL edge,
      !! or an MPI seam (`has_* = .false.`), leaves the physical-edge
      !! flux computed rather than hard-zeroed — the halo/periodic-wrap
      !! preamble the split-solver caller runs immediately before this
      !! call has already filled the ghost band with the correct
      !! neighbour/periodic value there.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_hdiff_tracer_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
      type(ocean_bc_state_t), intent(in), optional :: bc

      integer :: it
      logical :: wall_w, wall_e, wall_s, wall_n

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      if (this%kappa_h <= 0.0_wp) return

      ! Resolve the physical-edge wall flags on the host — see
      ! `continuity_zonal_flux` for the identical pattern.  Default
      ! .true. => single-rank, all-wall closed boundary (bit-identical
      ! to a run with no `bc`).
      wall_w = .true.
      wall_e = .true.
      wall_s = .true.
      wall_n = .true.
      if (present(bc)) then
         wall_w = (bc%west%bc_type == OBC_WALL) .and. bc%has_west
         wall_e = (bc%east%bc_type == OBC_WALL) .and. bc%has_east
         wall_s = (bc%south%bc_type == OBC_WALL) .and. bc%has_south
         wall_n = (bc%north%bc_type == OBC_WALL) .and. bc%has_north
      end if

      ! The mask actuals are chosen ONCE, outside the tracer loop: they are
      ! passed as ABSENT optionals on the default path (see the `open_u`
      ! docstring for why an inert stand-in is not available here), and
      ! an absent optional cannot be selected inside an expression, so
      ! the select-case is written twice rather than the arguments once.
      ! Thermo cadence, cold code.
      if (metrics%use_closed_faces) then
         do it = 1, size(ms%tracers)
            if (.not. ms%tracers(it)%do_horizontal_diffusion) cycle
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n, &
                  open_u=metrics%open_u, open_v=metrics%open_v, &
                  budget=ms%heat_budget_hdiff)
            case (TRACER_BUDGET_SALT)
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n, &
                  open_u=metrics%open_u, open_v=metrics%open_v, &
                  budget=ms%salt_budget_hdiff)
            case default
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n, &
                  open_u=metrics%open_u, open_v=metrics%open_v)
            end select
         end do
      else
         do it = 1, size(ms%tracers)
            if (.not. ms%tracers(it)%do_horizontal_diffusion) cycle
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n, &
                  budget=ms%heat_budget_hdiff)
            case (TRACER_BUDGET_SALT)
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n, &
                  budget=ms%salt_budget_hdiff)
            case default
               call tracer_hdiff_one_impl( &
                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                  dt, this%kappa_h, &
                  metrics%dy_cu, metrics%dx_cv, metrics%idxCu, metrics%idyCv, &
                  metrics%iareaT, &
                  ms%h_layer, ms%tracers(it)%hTr, &
                  this%T_centre%data, &
                  this%F_x_face%data, this%F_y_face%data, &
                  wall_w, wall_e, wall_s, wall_n)
            end select
         end do
      end if
   end subroutine tracer_hdiff

   pure subroutine tracer_hdiff_one_impl(nx, ny, nz, nghost, nx_phys, ny_phys, dt, kappa, &
                                         dy_cu, dx_cv, idxCu, idyCv, iareaT, &
                                         h, hTr, T_centre, F_x_face, F_y_face, &
                                         wall_w, wall_e, wall_s, wall_n, &
                                         open_u, open_v, budget)
      !! Four-pass flat-impl horizontal-Laplacian tracer diffusion in
      !! conservative curvilinear form (design §2, mirrors continuity):
      !!
      !!   1. Compute T = hTr / h at every cell centre (vanishing
      !!      layer → T = 0).
      !!   2. East-face TRANSPORT flux [tracer·m³/s]:
      !!        F_x(i) = kappa · 0.5·(h(i-1)+h(i)) · (T(i)-T(i-1))
      !!                 · idxCu(i) · dy_cu(i)
      !!      PHYSICAL wall faces (i = nghost+1, i = nghost+nx_phys+1):
      !!      zeroed when `wall_w`/`wall_e` (mirrors
      !!      `continuity_zonal_flux`'s host-resolved has_*/OBC_WALL
      !!      gate).  ARRAY-bound faces (i = 1, i = nx+1) are always
      !!      zeroed too — the pass-4 divergence reads F_x_face(1,...)
      !!      and F_x_face(nx+1,...) at the domain edge cells even
      !!      when the physical wall sits inside the ghost band, so
      !!      those two faces must stay initialised.
      !!   3. North-face flux, mirror of pass 2 (idyCv · dx_cv).
      !!   4. Divergence into hTr:
      !!        hTr(i,j,k) += dt · [(F_x(i+1)-F_x(i))
      !!                            + (F_y(j+1)-F_y(j))] · iareaT(i,j)
      !!
      !! On uniform Cartesian `dy_cu=dy`, `idxCu=1/dx`,
      !! `iareaT=1/(dx·dy)`, so the form collapses to the old
      !! `Δ(kappa·h·ΔT/dx)/dx` to round-off.  Each pass writes a
      !! different buffer than it reads, so `do concurrent` is
      !! race-free.  Constancy: uniform T → zero face fluxes → hTr
      !! unchanged.  Conservation over the PHYSICAL domain: closed
      !! physical-wall + flux-form divergence over the area-weighted
      !! cells integrates to zero exactly (see `tracer_hdiff`'s
      !! host-side wall-flag resolution for the OBC/MPI-seam gate).
      !!
      !! Stability (per cell, explicit forward-Euler):
      !!   kappa · dt · (idxT² + idyT²) <= 0.5.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      real(wp), intent(in) :: dt, kappa
      real(wp), intent(in) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: T_centre(nx, ny, nz)
      real(wp), intent(inout) :: F_x_face(nx + 1, ny, nz)
      real(wp), intent(inout) :: F_y_face(nx, ny + 1, nz)
      logical, intent(in) :: wall_w, wall_e, wall_s, wall_n
      real(wp), intent(in), optional :: open_u(nx + 1, ny, nz)
         !! Per-layer 0/1 u-face open mask
         !! (`&vcoord_nml zfixed_closed_faces`).  A CLOSED face is a
         !! z-level WALL for that layer, so it carries no diffusive
         !! tracer flux either — the same statement continuity makes
         !! about mass.
         !!
         !! ABSENT (the default path) ⇒ no masking pass is generated and
         !! every expression below is byte-identical to the un-masked
         !! form.  It is OPTIONAL rather than a `use_open` + inert
         !! stand-in pair precisely because the `(1,1,1)` placeholder
         !! must never reach an explicit-shape dummy, and this routine
         !! has no full-size read-only array of its own to lend (the
         !! `F_*_face` buffers it would otherwise borrow are its own
         !! `intent(inout)` scratch, so lending them would alias).
      real(wp), intent(in), optional :: open_v(nx, ny + 1, nz)
         !! v-face twin.  Present iff `open_u` is.
      real(wp), intent(inout), optional :: budget(nx, ny, nz)

      integer :: i, j, k
      real(wp) :: h_face, delta

      ! ---- Pass 1: T = hTr/h at centres ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         if (h(i, j, k) > 0.0_wp) then
            T_centre(i, j, k) = hTr(i, j, k)/h(i, j, k)
         else
            T_centre(i, j, k) = 0.0_wp
         end if
      end do

      ! ---- Pass 2: east-face transport flux ----
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(h_face)
         h_face = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
         F_x_face(i, j, k) = kappa*h_face* &
                             (T_centre(i, j, k) - T_centre(i - 1, j, k))* &
                             idxCu(i, j)*dy_cu(i, j)
      end do
      ! z-level closed faces: a separate host-gated pass so the loop above
      ! is textually unchanged with the knob off.
      if (present(open_u)) then
         do concurrent(k=1:nz, j=1:ny, i=2:nx)
            F_x_face(i, j, k) = F_x_face(i, j, k)*open_u(i, j, k)
         end do
      end if
      ! Array-bound faces: always zeroed (pass-4 divergence at the
      ! domain-edge cells reads them even when the physical wall is
      ! elsewhere inside the ghost band).
      do concurrent(k=1:nz, j=1:ny)
         F_x_face(1, j, k) = 0.0_wp
         F_x_face(nx + 1, j, k) = 0.0_wp
      end do
      ! Physical wall faces: zero only when the edge is actually a
      ! closed wall here (single-rank all-wall default; OBC/MPI-seam
      ! edges keep the computed flux read from a correctly-filled
      ! ghost column — see `tracer_hdiff`).
      do concurrent(k=1:nz, j=1:ny)
         if (wall_w) F_x_face(nghost + 1, j, k) = 0.0_wp
         if (wall_e) F_x_face(nghost + nx_phys + 1, j, k) = 0.0_wp
      end do

      ! ---- Pass 3: north-face transport flux ----
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(h_face)
         h_face = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
         F_y_face(i, j, k) = kappa*h_face* &
                             (T_centre(i, j, k) - T_centre(i, j - 1, k))* &
                             idyCv(i, j)*dx_cv(i, j)
      end do
      ! z-level closed faces: see the zonal twin.
      if (present(open_v)) then
         do concurrent(k=1:nz, j=2:ny, i=1:nx)
            F_y_face(i, j, k) = F_y_face(i, j, k)*open_v(i, j, k)
         end do
      end if
      do concurrent(k=1:nz, i=1:nx)
         F_y_face(i, 1, k) = 0.0_wp
         F_y_face(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, i=1:nx)
         if (wall_s) F_y_face(i, nghost + 1, k) = 0.0_wp
         if (wall_n) F_y_face(i, nghost + ny_phys + 1, k) = 0.0_wp
      end do

      ! ---- Pass 4: apply area-weighted divergence ----
      if (present(budget)) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx) local(delta)
            delta = dt*( &
                    (F_x_face(i + 1, j, k) - F_x_face(i, j, k)) + &
                    (F_y_face(i, j + 1, k) - F_y_face(i, j, k)))*iareaT(i, j)
            hTr(i, j, k) = hTr(i, j, k) + delta
            budget(i, j, k) = budget(i, j, k) + delta
         end do
      else
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            hTr(i, j, k) = hTr(i, j, k) + dt*( &
                           (F_x_face(i + 1, j, k) - F_x_face(i, j, k)) + &
                           (F_y_face(i, j + 1, k) - F_y_face(i, j, k)))*iareaT(i, j)
         end do
      end if
   end subroutine tracer_hdiff_one_impl

   pure function ocean_hdiff_tracer_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the horizontal tracer diffusion slot (0 when
      !! unallocated).
      class(ocean_hdiff_tracer_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%T_centre%bytes() &
               + this%F_x_face%bytes() &
               + this%F_y_face%bytes()
   end function ocean_hdiff_tracer_bytes

end module rdb_ocean_hdiff_tracer
