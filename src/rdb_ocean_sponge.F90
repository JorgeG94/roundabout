!! Sponge relaxation kernel.
module rdb_ocean_sponge
   !! Two sponge implementations share this module.
   !!
   !! LEGACY BAND PATH (`ocean_sponge_apply` / `ocean_sponge_apply_tracers`,
   !! `&ocean_bc_nml`): dissipates outgoing waves in a band of cells next to
   !! a sponge-tagged edge. Reads `bc%<edge>%bc_type == OBC_SPONGE` +
   !! per-edge `sponge_width` (cells) and `sponge_strength` (1/s) and
   !! applies `du/dt += -τ·u`, with `τ` a cosine ramp from 0 at the interior
   !! to `sponge_strength` at the outer edge, toward **zero** momentum and a
   !! **scalar** per-edge tracer target. State evolution (`hTr`, momentum)
   !! kept byte-identical — this is what `&ocean_sponge_nml enable=.false.`
   !! (the default) still runs. PR-23 additionally mirrors the tracer
   !! relaxation's salinity/temperature increment into `ms%salt_budget_sponge`
   !! / `heat_budget_sponge` (a pure side-channel add, new arrays that were
   !! always zero before — no change to `hTr` itself) so the console
   !! salt/heat budget closes with the legacy `sponge_relax_tracers = .true.`
   !! knob too, not just the map-driven path.
   !!
   !! MAP-DRIVEN PATH (`ocean_sponge_t` + `ocean_sponge_apply_maps`,
   !! `&ocean_sponge_nml enable=.true.`): a real sponge (MOM6's ALE sponge
   !! is the model in spirit, not in code — see below).
   !! A per-cell inverse-damping-time map (`idamp_h`/`idamp_u`/`idamp_v`,
   !! 1/s; `Idamp = 0` IS the sponge mask, no separate width/extent
   !! bookkeeping) relaxes momentum toward a 3-D `u_ref`/`v_ref` (not zero)
   !! and every registered tracer toward a 3-D `ref_tracer(i,j,k,it)`
   !! concentration field (not a scalar). `damp_source="band"` (the only
   !! implemented source in v1) fills the maps from the same cosine ramp as
   !! the legacy kernel, at the exact cell/u-face/v-face offsets the legacy
   !! kernel touches, so `enable=.true., damp_source="band"` is *physically*
   !! the same band as today (see `sponge_source_is_implemented` /
   !! `rdb_ocean_setup::configure_ocean_sponge`). `target_source="ic"` (the
   !! only implemented target in v1) snapshots the reference from the
   !! seeded initial condition via `ocean_sponge_snapshot_reference` —
   !! reachable as a "nudge toward a parent climatology" via
   !! `&ocean_zinit_nml` with no new file reader.
   !!
   !! Divergences from MOM6's ALE-sponge / sponge approach (deliberate,
   !! documented per CLAUDE.md "cite the paper not other codebases" — this
   !! list records WHERE we differ, not MOM6's implementation):
   !!   1. Exponential relaxation `decay = exp(-Idamp*dt); phi <-
   !!      decay*phi + (1-decay)*phi_ref` (exact for constant phi_ref over
   !!      the step, unconditionally stable, monotone) instead of MOM6's
   !!      backward-Euler `I1pdamp = 1/(1+Idamp*dt)` (physics-equivalent to
   !!      O(dt), less accurate).
   !!   2. Dense `(nx,ny)`/`(nx,ny,nz)` maps, not MOM6's compressed-column
   !!      `Iresttime_col`/`Ref_val%p(k,c)` sparse structure — a
   !!      gather/scatter indirection is worse than `idamp=0` (one multiply,
   !!      race-free `do concurrent`) on GPU.
   !!   3. The reference is snapshotted ONCE on model layers at seed time,
   !!      not held on a source z-grid and remapped to the live column every
   !!      apply (MOM6's `Ref_dz` + `remapping_core_h`). Valid because
   !!      Roundabout's ALE remap pins layer depths to the coordinate every
   !!      thermo step; the residual motion is far below the target's own
   !!      vertical resolution.
   !!   4. Applied at `dt` (every RK2 stage), not MOM6's diabatic/thermo
   !!      cadence.
   !!
   !! Restart contract (mirrors `rdb_ocean_surface_flux`'s docstring): the
   !! `idamp_*` maps and the `ref_*` fields are NOT registered in the
   !! restart registry. The maps are a pure function of config + geometry
   !! (restart-invariant, rebuilt every configure). The reference is
   !! re-derived from the re-seeded initial condition on every resume via
   !! `ocean_sponge_snapshot_reference`, called BEFORE the restart read (see
   !! `rdb_driver.F90`) — this is deliberate: `target_source="ic"` must mean
   !! the IC, not whatever state a resumed run happens to be in (see
   !! `docs/plans/PLAN_PR23_real_sponge.md` §13.1 item 4 for the reasoning;
   !! `test_ocean_sponge::sponge_reference_is_the_ic_not_the_restart` pins
   !! it).
   use rdb_constants, only: wp, H_DIV_EPS, H_VANISHED
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, OBC_SPONGE
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_sponge_apply
   public :: ocean_sponge_apply_tracers
   public :: ocean_sponge_t
   public :: ocean_sponge_apply_maps
   public :: ocean_sponge_snapshot_reference

   type :: ocean_sponge_t
      !! Map-driven sponge state: per-cell `Idamp` [1/s] + a 3-D reference
      !! state. See the module docstring for the physics + the MOM6
      !! divergences. Value-semantics slot (no `CS` pointer, no
      !! `associated(CS)` guards) per `src/core/ocean/README.md`.
      logical :: is_init = .false.
         !! True once `init` has run (only called when `enable`, mirroring
         !! the gated-closure convention — epbl/kshear/... — so a disabled
         !! sponge never allocates the (potentially large) `ref_tracer`).
      logical :: enable = .false.
         !! Master switch. Default `.false.` ⇒ the legacy band kernels run
         !! unchanged ⇒ existing nmls + tests are bit-identical.
      logical :: relax_uv = .true.
         !! Relax `u_face_x_layer`/`v_face_y_layer` toward `u_ref`/`v_ref`.
         !! Default `.true.` matches the legacy path (momentum is the one
         !! thing today's sponge always damps); MOM6's `SPONGE_UV` defaults
         !! `.false.` — a deliberate divergence to preserve legacy parity.
      logical :: relax_tracers = .true.
         !! Relax every registered tracer's concentration toward
         !! `ref_tracer`.
      logical :: relax_h = .false.
         !! Interior-interface thickness damping. NOT IMPLEMENTED in PR-23
         !! v1 (deferred to PR-23b alongside the file-backed targets —
         !! `docs/plans/PLAN_PR23_real_sponge.md` §14 Q1); `enable=.true.,
         !! relax_h=.true.` aborts fail-loud at `validate_config`.
      character(len=16) :: damp_source = "band"
         !! How `idamp_h`/`idamp_u`/`idamp_v` are filled. `"band"` (only
         !! value implemented in v1): cosine ramp from every
         !! `bc%<edge>%bc_type == OBC_SPONGE` edge tag, summed at overlaps
         !! (§3.2 of the plan). `"file"` recognised but aborts at
         !! `validate_config` (PR-23b, needs the PR-14 reader).
      character(len=16) :: target_source = "ic"
         !! Reference-state source. `"ic"` (only value implemented in v1):
         !! `ocean_sponge_snapshot_reference` copies the seeded initial
         !! condition. `"file"` recognised but aborts at `validate_config`
         !! (PR-23b).
      integer :: n_tracers = 0
         !! Registered tracer count — sizes `ref_tracer`'s 4th dimension.
      real(wp), allocatable :: idamp_h(:, :)
         !! Inverse damping time for tracers + (deferred) thickness, 1/s,
         !! shape `(nx_total, ny_total)`. Zero outside the sponge and in
         !! every ghost cell — `Idamp = 0` IS the sponge mask (no separate
         !! width/extent bookkeeping).
      real(wp), allocatable :: idamp_u(:, :)
         !! Inverse damping time for `u_face_x_layer`, 1/s, shape
         !! `(nx_total+1, ny_total)` (matches the u-face stagger).
      real(wp), allocatable :: idamp_v(:, :)
         !! Inverse damping time for `v_face_y_layer`, 1/s, shape
         !! `(nx_total, ny_total+1)` (matches the v-face stagger).
      real(wp), allocatable :: ref_tracer(:, :, :, :)
         !! Reference tracer CONCENTRATION (not `hTr`) — PSU, degC, ... per
         !! the registered tracer's own units — shape `(nx_total, ny_total,
         !! nz, n_tracers)`, indexed by the multilayer tracer-registry
         !! index (matches `ms%idx_salinity` / `ms%idx_temperature`).
      real(wp), allocatable :: u_ref(:, :, :)
         !! Reference x-velocity, shape `(nx_total+1, ny_total, nz)` —
         !! matches `ms%u_face_x_layer`.
      real(wp), allocatable :: v_ref(:, :, :)
         !! Reference y-velocity, shape `(nx_total, ny_total+1, nz)` —
         !! matches `ms%v_face_y_layer`.
   contains
      procedure, non_overridable :: init => ocean_sponge_init
      procedure, non_overridable :: destroy => ocean_sponge_destroy
      procedure, non_overridable :: enter_data => ocean_sponge_enter_data
      procedure, non_overridable :: exit_data => ocean_sponge_exit_data
      procedure, non_overridable :: bytes => ocean_sponge_bytes
   end type ocean_sponge_t

contains

   subroutine ocean_sponge_init(this, grid, nz_ml, n_tracers)
      !! Allocate the maps + reference-state arrays. Call ONLY when
      !! `this%enable` (the gated-closure convention, mirroring
      !! epbl/kshear/...) — the caller sets `enable` before calling, and
      !! `ocean_state_init` only reaches this when `this%enable` is
      !! `.true.` (see `rdb_ocean_state.F90`). `ref_tracer` is sized at
      !! `max(n_tracers, 1)` so an as-yet-empty tracer registry never
      !! trips a zero-extent allocate.
      class(ocean_sponge_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_ml
      integer, intent(in) :: n_tracers

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      this%n_tracers = n_tracers
      allocate (this%idamp_h(nx, ny), source=0.0_wp)
      allocate (this%idamp_u(nx + 1, ny), source=0.0_wp)
      allocate (this%idamp_v(nx, ny + 1), source=0.0_wp)
      allocate (this%ref_tracer(nx, ny, nz_ml, max(n_tracers, 1)), source=0.0_wp)
      allocate (this%u_ref(nx + 1, ny, nz_ml), source=0.0_wp)
      allocate (this%v_ref(nx, ny + 1, nz_ml), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_sponge_init

   subroutine ocean_sponge_destroy(this)
      class(ocean_sponge_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%idamp_h)) deallocate (this%idamp_h)
      if (allocated(this%idamp_u)) deallocate (this%idamp_u)
      if (allocated(this%idamp_v)) deallocate (this%idamp_v)
      if (allocated(this%ref_tracer)) deallocate (this%ref_tracer)
      if (allocated(this%u_ref)) deallocate (this%u_ref)
      if (allocated(this%v_ref)) deallocate (this%v_ref)
   end subroutine ocean_sponge_destroy

   subroutine ocean_sponge_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! device-attach map base is the heap object, not a polymorphic stack
      !! box (AMD libomptarget cross-slot-overlap fix; copied verbatim from
      !! `rdb_ocean_surface_flux.F90`).
      class(ocean_sponge_t), intent(inout) :: this
      select type (this)
      type is (ocean_sponge_t)
         call ocean_sponge_enter_data_impl(this)
      end select
   end subroutine ocean_sponge_enter_data

   subroutine ocean_sponge_enter_data_impl(this)
      type(ocean_sponge_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%idamp_h, this%idamp_u, this%idamp_v, &
      !$acc&                  this%ref_tracer, this%u_ref, this%v_ref)
      !$acc update device(this%idamp_h, this%idamp_u, this%idamp_v, &
      !$acc&               this%ref_tracer, this%u_ref, this%v_ref)
   end subroutine ocean_sponge_enter_data_impl

   subroutine ocean_sponge_exit_data(this)
      class(ocean_sponge_t), intent(inout) :: this
      select type (this)
      type is (ocean_sponge_t)
         call ocean_sponge_exit_data_impl(this)
      end select
   end subroutine ocean_sponge_exit_data

   subroutine ocean_sponge_exit_data_impl(this)
      type(ocean_sponge_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%idamp_h, this%idamp_u, this%idamp_v, &
      !$acc&                 this%ref_tracer, this%u_ref, this%v_ref)
   end subroutine ocean_sponge_exit_data_impl

   pure function ocean_sponge_bytes(this) result(nbytes)
      !! Counted allocatable footprint (0 when `enable=.false.` — every
      !! array is unallocated then). `tools/check_bytes_accounting.py`
      !! reconciles this against the measured device mapping.
      class(ocean_sponge_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%idamp_h) + arr_bytes(this%idamp_u) &
               + arr_bytes(this%idamp_v) + arr_bytes(this%ref_tracer) &
               + arr_bytes(this%u_ref) + arr_bytes(this%v_ref)
   end function ocean_sponge_bytes

   subroutine ocean_sponge_snapshot_reference(sp, grid, ms)
      !! Snapshot the seeded initial condition into `sp%ref_tracer` /
      !! `sp%u_ref` / `sp%v_ref` (`target_source = "ic"`, the only
      !! implemented source in v1). HOST-side, plain `do` loops (mirrors
      !! `seed_ts_from_zfile` — CLAUDE.md gotcha: this runs before
      !! `ocean_state_enter_data`).
      !!
      !! CALL-SITE CONTRACT (load-bearing, see the module docstring +
      !! `docs/plans/PLAN_PR23_real_sponge.md` §13.1 item 4): must run
      !! AFTER `ocean_state_seed_from_cfg` and BEFORE
      !! `ocean_state_restart_read` in `rdb_driver.F90`. Snapshotting any
      !! later would capture a warm-restarted run's mid-run state instead
      !! of the IC — a silent, resume-dependent physics change with a
      !! bit-identical first step and no error.
      !!
      !! Tracer concentration is `hTr / h_layer`, guarded by `H_DIV_EPS`
      !! (pure 1/0 armour, D4 taxonomy). Under `VCOORD_ZSTAR_FULL` bed-side
      !! layers can vanish (`h_layer <= H_VANISHED`); those columns fall
      !! back to the nearest massive layer's concentration (bed-up then
      !! surface-down fill), mirroring the `k_dep` scan in
      !! `apply_geothermal_src_impl`. No-op when `.not. sp%enable` (keeps
      !! the default-off path free of extra work) or `target_source /=
      !! "ic"` (the `"file"` source is PR-23b; `validate_config` has
      !! already aborted before this runs if requested in v1).
      type(ocean_sponge_t), intent(inout) :: sp
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, it, nx, ny, nz

      if (.not. sp%enable) return
      if (trim(sp%target_source) /= "ic") return
      if (.not. sp%is_init) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (allocated(ms%tracers)) then
         do it = 1, min(size(ms%tracers), sp%n_tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            do j = 1, ny
               do i = 1, nx
                  call snapshot_column_concentration( &
                     sp%ref_tracer(i, j, :, it), &
                     ms%tracers(it)%hTr(i, j, :), ms%h_layer(i, j, :), nz)
               end do
            end do
         end do
      end if

      if (allocated(ms%u_face_x_layer) .and. size(ms%u_face_x_layer, 1) == nx + 1 &
          .and. size(ms%u_face_x_layer, 2) == ny) then
         sp%u_ref = ms%u_face_x_layer
      end if
      if (allocated(ms%v_face_y_layer) .and. size(ms%v_face_y_layer, 1) == nx &
          .and. size(ms%v_face_y_layer, 2) == ny + 1) then
         sp%v_ref = ms%v_face_y_layer
      end if
   end subroutine ocean_sponge_snapshot_reference

   pure subroutine snapshot_column_concentration(conc, hTr_col, h_col, nz)
      !! Fill one water column's reference concentration from `hTr/h`,
      !! falling back to the nearest massive layer (`h > H_VANISHED`) for
      !! any vanished layer — bed-up pass first, then a surface-down pass
      !! to backfill any vanished layers below the first massive one.
      !! Plain host loop (called from a host `do i,j` loop, never `do
      !! concurrent` — this is configure-time setup, not a per-step kernel).
      integer, intent(in) :: nz
      real(wp), intent(out) :: conc(nz)
      real(wp), intent(in)  :: hTr_col(nz), h_col(nz)
      integer :: k
      real(wp) :: last_valid
      logical :: have_valid

      have_valid = .false.
      last_valid = 0.0_wp
      do k = 1, nz
         if (h_col(k) > H_VANISHED) then
            conc(k) = hTr_col(k)/max(h_col(k), H_DIV_EPS)
            last_valid = conc(k)
            have_valid = .true.
         else if (have_valid) then
            conc(k) = last_valid
         else
            conc(k) = 0.0_wp   ! filled by the backward pass below, if ever massive
         end if
      end do
      if (.not. have_valid) return
      last_valid = conc(nz)
      do k = nz, 1, -1
         if (h_col(k) > H_VANISHED) then
            last_valid = conc(k)
         else
            conc(k) = last_valid
         end if
      end do
   end subroutine snapshot_column_concentration

   subroutine ocean_sponge_apply_maps(grid, sp, ms, dt)
      !! Map-driven sponge dispatch (`&ocean_sponge_nml enable=.true.`).
      !! Relaxes momentum toward `u_ref`/`v_ref` (when `relax_uv`) and every
      !! registered tracer toward `ref_tracer` (when `relax_tracers`),
      !! mirroring the S/T tracer relaxation into `ms%salt_budget_sponge` /
      !! `ms%heat_budget_sponge` so the console salt/heat budget can close
      !! with the sponge on (`rdb_ocean_console_stats::ocean_heat_src_sum` /
      !! `ocean_salt_src_sum`). No-op when `.not. sp%is_init .or. .not.
      !! sp%enable`. Run from the same slot as the legacy path — see
      !! `rdb_ocean_dyn::run_stage_split`'s dispatch (exactly one of the
      !! legacy band / map-driven path runs).
      type(hgrid_t), intent(in) :: grid
      type(ocean_sponge_t), intent(in) :: sp
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      integer :: nx, ny, nz, it

      if (.not. sp%is_init .or. .not. sp%enable) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (sp%relax_uv) then
         call relax_map_u_impl(ms%u_face_x_layer, sp%u_ref, sp%idamp_u, nx + 1, ny, nz, dt)
         call relax_map_v_impl(ms%v_face_y_layer, sp%v_ref, sp%idamp_v, nx, ny + 1, nz, dt)
      end if

      if (sp%relax_tracers .and. allocated(ms%tracers)) then
         ! Outer-shim: the per-tracer `it` loop stays outside the `do
         ! concurrent` kernels (§6.4 of the plan) — `ms%tracers(it)%hTr` is
         ! dereferenced here on the host, never inside the device loop.
         do it = 1, min(size(ms%tracers), sp%n_tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (it == ms%idx_salinity) then
               call relax_map_tracer_budget_impl(ms%tracers(it)%hTr, ms%h_layer, &
                                                 sp%ref_tracer, it, sp%n_tracers, &
                                                 sp%idamp_h, ms%salt_budget_sponge, &
                                                 nx, ny, nz, dt)
            else if (it == ms%idx_temperature) then
               call relax_map_tracer_budget_impl(ms%tracers(it)%hTr, ms%h_layer, &
                                                 sp%ref_tracer, it, sp%n_tracers, &
                                                 sp%idamp_h, ms%heat_budget_sponge, &
                                                 nx, ny, nz, dt)
            else
               call relax_map_tracer_impl(ms%tracers(it)%hTr, ms%h_layer, &
                                          sp%ref_tracer, it, sp%n_tracers, &
                                          sp%idamp_h, nx, ny, nz, dt)
            end if
         end do
      end if
   end subroutine ocean_sponge_apply_maps

   pure subroutine relax_map_u_impl(u, u_ref, idamp_u, nxu, ny, nz, dt)
      !! Relax `u_face_x_layer` toward `u_ref` at rate `idamp_u`. Explicit-
      !! shape dummies, dims declared first (decl-order hook). `idamp_u(i,j)
      !! <= 0` is a bit-exact no-op (Idamp = 0 IS the sponge mask; also
      !! guarantees `sponge_idamp_zero_is_exact_identity`, since `decay =
      !! exp(-0*dt) = 1.0` exactly would already be a no-op algebraically —
      !! the early `if` skips the arithmetic entirely instead of relying on
      !! that).
      integer, intent(in)    :: nxu, ny, nz
      real(wp), intent(inout) :: u(nxu, ny, nz)
      real(wp), intent(in)    :: u_ref(nxu, ny, nz)
      real(wp), intent(in)    :: idamp_u(nxu, ny)
      real(wp), intent(in)    :: dt
      integer :: i, j, k
      real(wp) :: decay
      do concurrent(k=1:nz, j=1:ny, i=1:nxu) local(decay)
         if (idamp_u(i, j) > 0.0_wp) then
            decay = exp(-idamp_u(i, j)*dt)
            u(i, j, k) = relax_toward(u(i, j, k), u_ref(i, j, k), decay)
         end if
      end do
   end subroutine relax_map_u_impl

   pure subroutine relax_map_v_impl(v, v_ref, idamp_v, nx, nyv, nz, dt)
      !! Relax `v_face_y_layer` toward `v_ref` at rate `idamp_v`. Mirror of
      !! `relax_map_u_impl` for the y-direction.
      integer, intent(in)    :: nx, nyv, nz
      real(wp), intent(inout) :: v(nx, nyv, nz)
      real(wp), intent(in)    :: v_ref(nx, nyv, nz)
      real(wp), intent(in)    :: idamp_v(nx, nyv)
      real(wp), intent(in)    :: dt
      integer :: i, j, k
      real(wp) :: decay
      do concurrent(k=1:nz, j=1:nyv, i=1:nx) local(decay)
         if (idamp_v(i, j) > 0.0_wp) then
            decay = exp(-idamp_v(i, j)*dt)
            v(i, j, k) = relax_toward(v(i, j, k), v_ref(i, j, k), decay)
         end if
      end do
   end subroutine relax_map_v_impl

   pure subroutine relax_map_tracer_impl(hTr, h_layer, ref_tracer, it, n_tr, &
                                         idamp_h, nx, ny, nz, dt)
      !! Relax one tracer's `hTr` toward `ref_tracer(:,:,:,it)*h_layer` at
      !! rate `idamp_h`, no budget mirror (every tracer except S/T — see
      !! `relax_map_tracer_budget_impl`). `ref_tracer` is passed WHOLE +
      !! indexed by the scalar `it` inside the kernel — never sliced by the
      !! caller (§6.4: a device array section of a mapped array is the
      !! descriptor-walk trap the `dc-assumed-shape` hook exists to catch).
      integer, intent(in)    :: it, n_tr, nx, ny, nz
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: ref_tracer(nx, ny, nz, n_tr)
      real(wp), intent(in)    :: idamp_h(nx, ny)
      real(wp), intent(in)    :: dt
      integer :: i, j, k
      real(wp) :: decay, tgt
      do concurrent(k=1:nz, j=1:ny, i=1:nx) local(decay, tgt)
         if (idamp_h(i, j) > 0.0_wp) then
            decay = exp(-idamp_h(i, j)*dt)
            tgt = ref_tracer(i, j, k, it)*h_layer(i, j, k)
            hTr(i, j, k) = relax_toward(hTr(i, j, k), tgt, decay)
         end if
      end do
   end subroutine relax_map_tracer_impl

   pure subroutine relax_map_tracer_budget_impl(hTr, h_layer, ref_tracer, it, n_tr, &
                                                idamp_h, budget, nx, ny, nz, dt)
      !! `relax_map_tracer_impl` + mirror the per-cell increment into
      !! `budget` (salt or heat), the S/T budget-instrumented variant —
      !! mirrors `apply_geothermal_src_impl`'s host-shim + flat-impl +
      !! budget-mirror shape. `delta` is the ALGEBRAIC increment
      !! `hTr_new - hTr_old = (1-decay)*(tgt-hTr_old)`, written once and
      !! used for BOTH the state update and the budget mirror so the two
      !! stay exactly consistent (no independent recomputation to drift).
      integer, intent(in)    :: it, n_tr, nx, ny, nz
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: ref_tracer(nx, ny, nz, n_tr)
      real(wp), intent(in)    :: idamp_h(nx, ny)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: dt
      integer :: i, j, k
      real(wp) :: decay, tgt, delta
      do concurrent(k=1:nz, j=1:ny, i=1:nx) local(decay, tgt, delta)
         if (idamp_h(i, j) > 0.0_wp) then
            decay = exp(-idamp_h(i, j)*dt)
            tgt = ref_tracer(i, j, k, it)*h_layer(i, j, k)
            delta = (1.0_wp - decay)*(tgt - hTr(i, j, k))
            hTr(i, j, k) = hTr(i, j, k) + delta
            budget(i, j, k) = budget(i, j, k) + delta
         end if
      end do
   end subroutine relax_map_tracer_budget_impl

   pure real(wp) function relax_toward(cur, tgt, decay) result(res)
      !! Shared core algebra for BOTH sponge paths: `cur*decay +
      !! tgt*(1-decay)` — the exact solution of `dphi/dt = -Idamp*(phi -
      !! phi_ref)` over one step at constant `phi_ref`, for whatever
      !! `decay = exp(-rate*dt)` the caller derived (`relax_one` derives it
      !! from the legacy edge/band cosine ramp; the map kernels above derive
      !! it directly from a per-cell `Idamp`). "One home" for the algebra
      !! per the plan (§3.1 / step 10).
      !$acc routine seq
      real(wp), intent(in) :: cur, tgt, decay
      res = cur*decay + tgt*(1.0_wp - decay)
   end function relax_toward

   subroutine ocean_sponge_apply(grid, bc, ms, dt)
      !! Apply momentum relaxation in any sponge-tagged edge band.
      !! No-op when no edge is OBC_SPONGE.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      integer :: i, j, k, nz, d
      integer :: i0, i1, j0, j1
      integer :: wall_face, band
      integer :: nx_u, ny_u, nx_v, ny_v
      real(wp) :: strength, tau, decay, alpha
      real(wp), parameter :: PI = acos(-1.0_wp)

      nz = ms%nz_ml
      i0 = grid%nghost + 1
      i1 = grid%nghost + grid%nx_phys
      j0 = grid%nghost + 1
      j1 = grid%nghost + grid%ny_phys
      ! Hoisted array extents for the in-band bounds guards (don't call
      ! size() on a mapped array inside a do concurrent body).
      nx_u = size(ms%u_face_x_layer, 1)
      ny_u = size(ms%u_face_x_layer, 2)
      nx_v = size(ms%v_face_y_layer, 1)
      ny_v = size(ms%v_face_y_layer, 2)

      ! GPU note: `do concurrent` on the device-resident layer velocities.
      ! Within one edge block every (d, ·, k) tuple writes a distinct
      ! element, so the DC is race-free; overlapping corner faces between
      ! edge blocks get both decays applied sequentially.

      ! ---- West edge ----
      ! has_* gate: a subdomain seam never carries a sponge band (O0).
      if (bc%west%bc_type == OBC_SPONGE .and. bc%west%sponge_width > 0 .and. bc%has_west) then
         band = bc%west%sponge_width
         strength = bc%west%sponge_strength
         wall_face = grid%nghost + 1
         do concurrent(k=1:nz, j=j0:j1, d=0:band - 1) local(alpha, tau, decay, i)
            ! Cosine ramp: τ peaks at d=0 (outer), tapers to 0 at the
            ! interior edge of the band.
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            ! u-face just east of the west wall sits at (wall_face + d + 1)
            i = wall_face + d + 1
            if (i >= 1 .and. i <= nx_u) then
               ms%u_face_x_layer(i, j, k) = decay*ms%u_face_x_layer(i, j, k)
            end if
         end do
         do concurrent(k=1:nz, j=j0:j1 + 1, d=0:band - 1) local(alpha, tau, decay, i)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            i = wall_face + d
            if (i >= 1 .and. i <= nx_v) then
               ms%v_face_y_layer(i, j, k) = decay*ms%v_face_y_layer(i, j, k)
            end if
         end do
      end if

      ! ---- East edge ----
      if (bc%east%bc_type == OBC_SPONGE .and. bc%east%sponge_width > 0 .and. bc%has_east) then
         band = bc%east%sponge_width
         strength = bc%east%sponge_strength
         wall_face = grid%nghost + grid%nx_phys + 1
         do concurrent(k=1:nz, j=j0:j1, d=0:band - 1) local(alpha, tau, decay, i)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            i = wall_face - d - 1
            if (i >= 1 .and. i <= nx_u) then
               ms%u_face_x_layer(i, j, k) = decay*ms%u_face_x_layer(i, j, k)
            end if
         end do
         do concurrent(k=1:nz, j=j0:j1 + 1, d=0:band - 1) local(alpha, tau, decay, i)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            i = wall_face - d - 1
            if (i >= 1 .and. i <= nx_v) then
               ms%v_face_y_layer(i, j, k) = decay*ms%v_face_y_layer(i, j, k)
            end if
         end do
      end if

      ! ---- South edge ----
      if (bc%south%bc_type == OBC_SPONGE .and. bc%south%sponge_width > 0 .and. bc%has_south) then
         band = bc%south%sponge_width
         strength = bc%south%sponge_strength
         wall_face = grid%nghost + 1
         do concurrent(k=1:nz, i=i0:i1 + 1, d=0:band - 1) local(alpha, tau, decay, j)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            j = wall_face + d
            if (j >= 1 .and. j <= ny_u) then
               ms%u_face_x_layer(i, j, k) = decay*ms%u_face_x_layer(i, j, k)
            end if
         end do
         do concurrent(k=1:nz, i=i0:i1, d=0:band - 1) local(alpha, tau, decay, j)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            j = wall_face + d + 1
            if (j >= 1 .and. j <= ny_v) then
               ms%v_face_y_layer(i, j, k) = decay*ms%v_face_y_layer(i, j, k)
            end if
         end do
      end if

      ! ---- North edge ----
      if (bc%north%bc_type == OBC_SPONGE .and. bc%north%sponge_width > 0 .and. bc%has_north) then
         band = bc%north%sponge_width
         strength = bc%north%sponge_strength
         wall_face = grid%nghost + grid%ny_phys + 1
         do concurrent(k=1:nz, i=i0:i1 + 1, d=0:band - 1) local(alpha, tau, decay, j)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            j = wall_face - d - 1
            if (j >= 1 .and. j <= ny_u) then
               ms%u_face_x_layer(i, j, k) = decay*ms%u_face_x_layer(i, j, k)
            end if
         end do
         do concurrent(k=1:nz, i=i0:i1, d=0:band - 1) local(alpha, tau, decay, j)
            alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
            tau = strength*alpha
            decay = exp(-tau*dt)
            j = wall_face - d - 1
            if (j >= 1 .and. j <= ny_v) then
               ms%v_face_y_layer(i, j, k) = decay*ms%v_face_y_layer(i, j, k)
            end if
         end do
      end if
   end subroutine ocean_sponge_apply

   subroutine ocean_sponge_apply_tracers(grid, bc, ms, dt)
      !! Relax tracer concentrations (hTr/h) toward per-edge targets in any
      !! sponge band with `sponge_relax_tracers = .true.`. Same cosine ramp
      !! and `sponge_strength` as the momentum kernel. No-op otherwise.
      !! Per-tracer loop outside the inner loops (outer-shim for the
      !! array-of-derived-types registry).
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      integer :: it
      integer :: i0, i1, j0, j1
      integer :: wall_face, band, nz, nxt, nyt
      real(wp) :: strength, C_bc

      nz = ms%nz_ml
      nxt = grid%nx_total
      nyt = grid%ny_total
      i0 = grid%nghost + 1
      i1 = grid%nghost + grid%nx_phys
      j0 = grid%nghost + 1
      j1 = grid%nghost + grid%ny_phys

      if (.not. allocated(ms%tracers)) return

      ! Per-tracer `it` loop stays OUTSIDE the `do concurrent` kernels
      ! (outer-shim): the hTr slice + h_layer reach the device as top-level
      ! explicit-shape allocatables in the `_impl` routines, never
      ! dereferenced as `ms%tracers(it)%hTr(...)` inside a device loop.

      ! ---- West edge ----
      ! has_* gate: a subdomain seam never carries a sponge band (O0).
      if (bc%west%bc_type == OBC_SPONGE .and. &
          bc%west%sponge_relax_tracers .and. &
          bc%west%sponge_width > 0 .and. &
          allocated(bc%west%clamped_tracer) .and. &
          bc%has_west) then
         band = bc%west%sponge_width
         strength = bc%west%sponge_strength
         wall_face = grid%nghost + 1
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (it > size(bc%west%clamped_tracer)) cycle
            C_bc = bc%west%clamped_tracer(it)
            call sponge_relax_band_x_tracer(ms, it, nxt, nyt, nz, &
                                            j0, j1, wall_face, band, +1, strength, C_bc, dt)
         end do
      end if

      ! ---- East edge ----
      if (bc%east%bc_type == OBC_SPONGE .and. &
          bc%east%sponge_relax_tracers .and. &
          bc%east%sponge_width > 0 .and. &
          allocated(bc%east%clamped_tracer) .and. &
          bc%has_east) then
         band = bc%east%sponge_width
         strength = bc%east%sponge_strength
         wall_face = grid%nghost + grid%nx_phys + 1
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (it > size(bc%east%clamped_tracer)) cycle
            C_bc = bc%east%clamped_tracer(it)
            call sponge_relax_band_x_tracer(ms, it, nxt, nyt, nz, &
                                            j0, j1, wall_face, band, -1, strength, C_bc, dt)
         end do
      end if

      ! ---- South edge ----
      if (bc%south%bc_type == OBC_SPONGE .and. &
          bc%south%sponge_relax_tracers .and. &
          bc%south%sponge_width > 0 .and. &
          allocated(bc%south%clamped_tracer) .and. &
          bc%has_south) then
         band = bc%south%sponge_width
         strength = bc%south%sponge_strength
         wall_face = grid%nghost + 1
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (it > size(bc%south%clamped_tracer)) cycle
            C_bc = bc%south%clamped_tracer(it)
            call sponge_relax_band_y_tracer(ms, it, nxt, nyt, nz, &
                                            i0, i1, wall_face, band, +1, strength, C_bc, dt)
         end do
      end if

      ! ---- North edge ----
      if (bc%north%bc_type == OBC_SPONGE .and. &
          bc%north%sponge_relax_tracers .and. &
          bc%north%sponge_width > 0 .and. &
          allocated(bc%north%clamped_tracer) .and. &
          bc%has_north) then
         band = bc%north%sponge_width
         strength = bc%north%sponge_strength
         wall_face = grid%nghost + grid%ny_phys + 1
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (it > size(bc%north%clamped_tracer)) cycle
            C_bc = bc%north%clamped_tracer(it)
            call sponge_relax_band_y_tracer(ms, it, nxt, nyt, nz, &
                                            i0, i1, wall_face, band, -1, strength, C_bc, dt)
         end do
      end if
   end subroutine ocean_sponge_apply_tracers

   subroutine sponge_relax_band_x_tracer(ms, it, nx_total, ny_total, nz, &
                                         j0, j1, wall_face, band, side, strength, C_bc, dt)
      !! Host dispatch for one tracer's `relax_band_x_impl` call: routes the
      !! salinity/temperature budget mirror (PR-23) so the legacy band
      !! sponge's tracer sink is instrumented exactly like the map-driven
      !! path, without threading an `if`-branch into every one of the four
      !! edge blocks in `ocean_sponge_apply_tracers`.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: it, nx_total, ny_total, nz, j0, j1, wall_face, band, side
      real(wp), intent(in) :: strength, C_bc, dt

      if (it == ms%idx_salinity) then
         call relax_band_x_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                j0, j1, wall_face, band, side, strength, C_bc, dt, &
                                budget=ms%salt_budget_sponge)
      else if (it == ms%idx_temperature) then
         call relax_band_x_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                j0, j1, wall_face, band, side, strength, C_bc, dt, &
                                budget=ms%heat_budget_sponge)
      else
         call relax_band_x_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                j0, j1, wall_face, band, side, strength, C_bc, dt)
      end if
   end subroutine sponge_relax_band_x_tracer

   subroutine sponge_relax_band_y_tracer(ms, it, nx_total, ny_total, nz, &
                                         i0, i1, wall_face, band, side, strength, C_bc, dt)
      !! Mirror of `sponge_relax_band_x_tracer` for the y-direction.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: it, nx_total, ny_total, nz, i0, i1, wall_face, band, side
      real(wp), intent(in) :: strength, C_bc, dt

      if (it == ms%idx_salinity) then
         call relax_band_y_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                i0, i1, wall_face, band, side, strength, C_bc, dt, &
                                budget=ms%salt_budget_sponge)
      else if (it == ms%idx_temperature) then
         call relax_band_y_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                i0, i1, wall_face, band, side, strength, C_bc, dt, &
                                budget=ms%heat_budget_sponge)
      else
         call relax_band_y_impl(ms%tracers(it)%hTr, ms%h_layer, nx_total, ny_total, nz, &
                                i0, i1, wall_face, band, side, strength, C_bc, dt)
      end if
   end subroutine sponge_relax_band_y_tracer

   pure subroutine relax_band_x_impl(hTr, h_layer, nx_total, ny_total, nz, &
                                     j0, j1, wall_face, band, side, strength, C_bc, dt, budget)
      !! Cosine-ramp tracer relaxation in a west/east sponge band.
      !! Explicit-shape dummies (flat-impl + outer-shim — the per-tracer
      !! slice is passed by the caller, never dereferenced inside the DC).
      !! `side = +1` for the west edge (cells wall_face .. wall_face+band-1),
      !! `side = -1` for the east edge (cells wall_face-1 .. wall_face-band).
      !!
      !! `budget` (PR-23, optional): when present, mirrors the per-cell
      !! `hTr` increment into it (salt or heat) so the legacy band sponge's
      !! tracer relaxation closes the console budget exactly like the
      !! map-driven path (`relax_map_tracer_budget_impl`) — the caller
      !! passes it only for the salinity/temperature tracer indices.
      !! `has_budget` is hoisted to a plain host logical and read (never
      !! `present()`) INSIDE the single `do concurrent` — splitting on
      !! `present()` into two separate DC loops is the
      !! `nvhpc_split_optional_dc_codegen` trap (~25x slower); this is the
      !! one-loop fix.
      integer, intent(in) :: nx_total, ny_total, nz, j0, j1, wall_face, band, side
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: strength, C_bc, dt
      real(wp), intent(inout), optional :: budget(nx_total, ny_total, nz)

      integer :: k, j, d, ii
      real(wp) :: new_val
      logical :: has_budget

      has_budget = present(budget)
      do concurrent(k=1:nz, j=j0:j1, d=0:band - 1) local(ii, new_val)
         if (side > 0) then
            ii = wall_face + d
         else
            ii = wall_face - d - 1
         end if
         new_val = relax_one(hTr(ii, j, k), C_bc*h_layer(ii, j, k), &
                             strength, d, band, dt)
         if (has_budget) budget(ii, j, k) = budget(ii, j, k) + (new_val - hTr(ii, j, k))
         hTr(ii, j, k) = new_val
      end do
   end subroutine relax_band_x_impl

   pure subroutine relax_band_y_impl(hTr, h_layer, nx_total, ny_total, nz, &
                                     i0, i1, wall_face, band, side, strength, C_bc, dt, budget)
      !! Cosine-ramp tracer relaxation in a south/north sponge band.
      !! Mirror of `relax_band_x_impl` for the y-direction; see its
      !! docstring for the optional `budget` mirror (PR-23).
      integer, intent(in) :: nx_total, ny_total, nz, i0, i1, wall_face, band, side
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: strength, C_bc, dt
      real(wp), intent(inout), optional :: budget(nx_total, ny_total, nz)

      integer :: k, i, d, jj
      real(wp) :: new_val
      logical :: has_budget

      has_budget = present(budget)
      do concurrent(k=1:nz, i=i0:i1, d=0:band - 1) local(jj, new_val)
         if (side > 0) then
            jj = wall_face + d
         else
            jj = wall_face - d - 1
         end if
         new_val = relax_one(hTr(i, jj, k), C_bc*h_layer(i, jj, k), &
                             strength, d, band, dt)
         if (has_budget) budget(i, jj, k) = budget(i, jj, k) + (new_val - hTr(i, jj, k))
         hTr(i, jj, k) = new_val
      end do
   end subroutine relax_band_y_impl

   pure real(wp) function relax_one(hTr_cur, tgt, strength, d, band, dt) result(res)
      !! Cosine-ramp implicit relaxation of one cell toward `tgt`.
      !! `decay = exp(-strength*alpha*dt)`, alpha the cosine ramp at band
      !! position `d`; result = hTr_cur*decay + tgt*(1-decay). Funnels
      !! through `relax_toward` for the shared core algebra (§3.1 / step 10
      !! "one home" — the map-driven kernels derive `decay` from a per-cell
      !! `Idamp` instead of this band/strength ramp, but both land on the
      !! same `cur*decay + tgt*(1-decay)` update).
      !$acc routine seq
      real(wp), intent(in) :: hTr_cur, tgt, strength, dt
      integer, intent(in) :: d, band
      real(wp) :: alpha, decay
      real(wp), parameter :: PI = acos(-1.0_wp)
      alpha = 0.5_wp*(1.0_wp + cos(PI*real(d, wp)/real(band, wp)))
      decay = exp(-strength*alpha*dt)
      res = relax_toward(hTr_cur, tgt, decay)
   end function relax_one

end module rdb_ocean_sponge
