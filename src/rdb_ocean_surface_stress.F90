!! Surface momentum stress for the ocean dynamical core.  Carries
!! a 2D wind-stress field `tau_x, tau_y` (N/m^2) and applies it as
!! a momentum source at the surface-most layer (`k = nz` under the
!! ROMS-style k=1-bed convention).
!!
!! Momentum equation: `du/dt|_stress = tau_x / (rho_0 * h_top)`.
!! `h_top` is the top layer thickness, averaged onto the face from
!! the two abutting cells.  Bottom drag's mirror image — same
!! tendency-buffer pattern, opposite k.
!!
!! `tau_x` lives on east faces (shape `(nx+1, ny)`); `tau_y` lives
!! on north faces (shape `(nx, ny+1)`).  Use
!! `set_wind_stress_const(tau_x, tau_y)` to fill them uniformly
!! (spatially-constant wind — the Tier-1 default).  Spatially-
!! varying wind sets the arrays directly before `enter_data`.
!!
!! ## Ice-shelf cover (`&ocean_cavity_dyn_nml`)
!!
!! Under an ice shelf there is no atmosphere, so there is no wind.  The
!! mask is applied ONCE to the `tau` PAIR
!! (`ocean_surface_stress_apply_cover`) rather than to each derived
!! view, because `tau_x`/`tau_y` have THREE independent consumers and
!! only one of them lives in this module:
!! `ocean_surface_stress_compute_tendencies` (the explicit momentum
!! source), the implicit surface-stress fold in `rdb_ocean_vdiff`
!! (`tau_u=ss%tau_x`, nine call sites into a hot tridiagonal kernel)
!! and the MLE front-stress sampler in `rdb_ocean_mle`.  Masking the
!! derived views would leave the other two blowing an atmosphere
!! through several hundred metres of solid ice; masking the source
!! makes all three correct by construction and adds no branch to any
!! hot kernel.
!!
!! FACE RULE — it is `&ocean_tdrag_nml`'s, not a second one.  The
!! projection of the cell-centred `cover_frac` onto faces is stated once,
!! in `rdb_ocean_top_drag`'s module docstring: a face is under ice if
!! EITHER abutting cell is (`max(cover_L, cover_R)`), rim faces are left
!! alone, and the CALVING-FRONT face is therefore closed to the wind
!! exactly where the top drag closes on it.  The two implementations are
!! held identical face-for-face by
!! `cavity_cover_face_rule_matches_top_drag`.  Its cost is one face-wide
!! transition at the front, where the first open cell receives half the
!! zonal impulse and therefore a `stress_mag` of `|tau|/2` — the `u*`
!! consistent with the momentum it actually got.
!!
!! CONTRACT: every writer of the `tau` pair must re-apply the cover
!! before refreshing `stress_mag` — `ocean_surface_stress_apply_cover`
!! does both, and `ocean_surface_stress_set_derived`'s optional
!! `cover_frac` routes through it.  Absent ⇒ the original path,
!! byte-identical.  The mask is idempotent (it multiplies by 0 or 1),
!! so re-applying it after every bracket read costs nothing.
!!
!! ## The two upper-boundary stresses (`stress_mag` + `stress_shelf`)
!!
!! Masking the wind is only half the physics: under a shelf the
!! turbulent boundary layer is not unforced, it is forced by the
!! ICE-OCEAN stress instead, and that stress is not in `tau` — the
!! ice-shelf top drag is a separate momentum tendency with its own slot.
!! So this type carries TWO cell-centred stress magnitudes and they are
!! different things:
!!
!!   * `stress_mag`   — `|tau|`, and ONLY `|tau|`.  Purely DERIVED from
!!     the `tau` pair, rebuilt from scratch by every writer of it
!!     (`ocean_surfstress_refresh_stress_mag`).  Carries the wind, the
!!     sea-ice blend, and the ice-shelf cover mask, because all three
!!     are written into `tau`.
!!   * `stress_shelf` — `|tau_top|` at an ICE-SHELF BASE (N/m^2), the
!!     cell-centred magnitude published by `rdb_ocean_top_drag` as
!!     `stress_top`, or (top drag off, basal melt on) `rho_0*u_*^2` from
!!     the melt slot's own `u_*` (via
!!     `ocean_surface_stress_set_shelf_from_ustar`).  NOT derived from
!!     `tau`; refreshed by the driver, not by the `tau` refresh; exactly
!!     zero without a cavity.
!!
!! `u_*^2 = (stress_mag + stress_shelf) / rho_0` is the ONE definition
!! both boundary-layer schemes use (`rdb_ocean_vmix` KPP,
!! `rdb_ocean_epbl`).  The SUM is the area-weighted total upper-boundary
!! momentum flux, not a double count: `tau` is zeroed on every face
!! touching a covered cell, so `stress_mag` is exactly zero wherever
!! `cover_frac = 1`, and `stress_top` is multiplied by `cover_frac`, so
!! it is exactly zero wherever `cover_frac = 0`.  The two supports are
!! disjoint under the binary v1 cover, and each term already carries its
!! own area weight, so the sum generalises unchanged to a fractional
!! cover.
!!
!! `stress_shelf` is kept SEPARATE from `stress_mag` rather than blended
!! into it for three reasons, each of which a folded-in design gets
!! wrong: (1) `stress_mag` is rebuilt from `tau` by every `tau` writer,
!! which would silently wipe a folded contribution at the next
!! data-forcing bracket or sea-ice blend; (2) the top-drag stress is
!! recomputed EVERY RK2 stage while the `tau` refresh is per outer step,
!! so a fold would either accumulate across stages or go stale; (3)
!! keeping `stress_mag` a pure function of `tau` is what makes the
!! existing cover gate (`test_ocean_cavity_flux`'s "exactly zero under
!! cover") and the sea-ice gates (`test_ocean_ice_stress_mag`) still
!! mean what they say.
!!
!! `stress_shelf` is ALWAYS allocated and zero-filled (`init`), mapped
!! (`enter_data`), counted (`bytes`) — never conditionally — so both
!! boundary-layer kernels have ONE code path: no optional dummy, no
!! placeholder array reaching an explicit-shape device kernel.  With no
!! cavity it is the zero array and `stress_mag + 0.0` is `stress_mag`
!! bit-for-bit under IEEE-754 (both are finite and non-negative).
module rdb_ocean_surface_stress
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_surface_stress_t
   public :: ocean_surface_stress_compute_tendencies
   public :: ocean_surface_stress_apply_tendencies
   public :: ocean_surface_stress_set_derived
   public :: ocean_surface_stress_refresh_mag
   public :: ocean_surface_stress_apply_cover
   public :: ocean_surface_stress_set_shelf_from_ustar

   type :: ocean_surface_stress_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3) used in the
         !! `tau / (rho_0 * h)` acceleration (top-layer and DIRECT_STRESS
         !! distributed forms alike).
         !!
         !! ASSIGNED FROM CONFIG by `configure_ocean_reference_density`,
         !! which copies the single rho0 of record (`&ocean_ic_nml rho_0`
         !! -> `eos%rho0`).  The literal here is only the pre-configure
         !! type default.  Host scalar: passed by value into the
         !! `surfstress_*_impl` kernels, so no `!$acc update device`.
      real(wp) :: h_min = 1.0e-3_wp
         !! Floor on the surface-layer thickness in the `1/h_top`
         !! division — keeps the kernel finite when the top layer
         !! pinches out (rare; e.g., wave breaking under ZSTAR_FULL).
      logical :: direct_stress = .false.
         !! MOM6 DIRECT_STRESS analogue.  When `.true.` the wind
         !! stress is distributed across the top `hmix_stress` metres
         !! rather than concentrated in the surface-most layer.
         !! Mirrors HYCOM's approach: stress acts on a surface slab
         !! of fixed thickness, with each layer in that slab getting
         !! a proportional share.  Has no effect when the top layer
         !! is already thinner than `hmix_stress` (bit-identical to
         !! the bed-only branch in that limit).
      real(wp) :: hmix_stress = 0.0_wp
         !! Thickness (m) of the surface slab over which the stress
         !! is spread when `direct_stress = .true.`.  MOM6 production
         !! default 20 m.  Zero (default) keeps the surface-layer-only
         !! behaviour even if `direct_stress` is flipped on.

      real(wp), allocatable :: tau_x(:, :)
         !! Zonal wind stress (N/m^2) on east faces, shape `(nx+1, ny)`.
         !! Fill via `set_wind_stress_const` for spatially-uniform
         !! wind, or assign the array directly for spatially-varying
         !! input.  Default zero.
      real(wp), allocatable :: tau_y(:, :)
         !! Meridional wind stress (N/m^2) on north faces, shape
         !! `(nx, ny+1)`.
      logical :: has_stress_mag = .false.
         !! True once `ocean_surface_stress_set_derived` has filled
         !! `stress_mag` at least once.  Informational — no kernel gates
         !! on it (`stress_mag` is always allocated and zero-safe).
      real(wp), allocatable :: stress_mag(:, :)
         !! Cell-centred wind-stress magnitude `|tau|` (Pa,
         !! `sqrt(tau_x_cell^2 + tau_y_cell^2)`), shape `(nx_total,
         !! ny_total)`.  **Always allocated** (PR-12 §7.4: a pure τ
         !! property, unlike a shared `ustar` which would need a
         !! coherent `rho0` — deferred).  Refreshed in step by every
         !! writer of the `tau` pair: the `set_wind_stress_*` setters and
         !! `ocean_surface_stress_set_derived` at configure, the
         !! data-forcing reader's seam refresh per bracket, and the
         !! sea-ice stress coupler's on-device blend every outer step
         !! (all via `ocean_surface_stress_refresh_mag`).  KPP
         !! (`rdb_ocean_vmix`)
         !! and EPBL (`rdb_ocean_epbl`) both read it instead of
         !! re-deriving `tau_mag` inline (bit-identical dedup — same
         !! three lines, same FP op order, just computed once) — which is
         !! why a `tau` write that skips the refresh silently freezes
         !! BOTH schemes' `u_*` at the last refreshed stress.
         !!
         !! Under an ice-shelf cavity the refresh is not enough on its
         !! own: a `tau` writer must ALSO re-apply the cover mask, which
         !! is why `ocean_surface_stress_apply_cover` does both in one
         !! call (see the module docstring's cover contract).  What the
         !! mask leaves behind under the ice — the ICE-OCEAN stress — is
         !! `stress_shelf` below, deliberately not folded in here.
      real(wp), allocatable :: stress_shelf(:, :)
         !! Cell-centred magnitude of the stress an ICE-SHELF BASE
         !! exerts on the ocean (N/m^2, `>= 0`), shape `(nx_total,
         !! ny_total)`, valid including ghosts.  **Always allocated**,
         !! mapped and counted; exactly zero without a cavity, and
         !! exactly zero on every cell with `cover_frac = 0`.
         !!
         !! NOT derived from `tau` — the ice-ocean stress is a separate
         !! momentum tendency (`rdb_ocean_top_drag`), so this field is
         !! NOT touched by `ocean_surfstress_refresh_stress_mag` and NOT
         !! touched by `ocean_surface_stress_apply_cover`.  It is
         !! refreshed by the split/unsplit RK2 drivers, inline, in the
         !! same stage that recomputes the top drag and strictly before
         !! `vmix_apply_in_stage` reads it (no lag), from
         !! `ocean_top_drag_t%stress_top`.  When the top drag is OFF but
         !! `&ocean_cavity_melt_nml enable` is on, `engine_step_finalize`
         !! fills it instead from `rho_0*u_*^2` with the melt slot's own
         !! `u_*` — the SAME `C_d` under the one-drag-coefficient rule,
         !! but at the thermo cadence, so THAT path is lagged one outer
         !! step.  Both are documented in the module docstring.
         !!
         !! Consumers: KPP (`rdb_ocean_vmix`) and EPBL
         !! (`rdb_ocean_epbl`), both as
         !! `u_* = sqrt((stress_mag + stress_shelf)/rho_0)`.

      type(scratch_3d_buffer_t) :: du_stress
         !! Surface stress tendency at east faces, shape
         !! (nx+1, ny, nz).  Only k=nz carries a non-zero value;
         !! k<nz stays zero.
      type(scratch_3d_buffer_t) :: dv_stress
         !! Surface stress tendency at north faces, shape
         !! (nx, ny+1, nz).
   contains
      procedure, non_overridable :: init => ocean_surfstress_init
      procedure, non_overridable :: destroy => ocean_surfstress_destroy
      procedure, non_overridable :: enter_data => ocean_surfstress_enter_data
      procedure, non_overridable :: exit_data => ocean_surfstress_exit_data
      procedure, non_overridable :: set_wind_stress_const => ocean_surfstress_set_const
      procedure, non_overridable :: set_wind_stress_2gyre => ocean_surfstress_set_2gyre
      procedure, non_overridable :: set_wind_stress_neverworld2 => ocean_surfstress_set_neverworld2
      procedure, non_overridable :: bytes => ocean_surface_stress_bytes
   end type ocean_surface_stress_t

contains

   subroutine ocean_surfstress_init(this, grid, nz_ml)
      class(ocean_surface_stress_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      allocate (this%tau_x(nx + 1, ny), source=0.0_wp)
      allocate (this%tau_y(nx, ny + 1), source=0.0_wp)
      allocate (this%stress_mag(nx, ny), source=0.0_wp)
      allocate (this%stress_shelf(nx, ny), source=0.0_wp)
      call this%du_stress%init(nx + 1, ny, nz, "ocean_surfstress_du_stress")
      call this%dv_stress%init(nx, ny + 1, nz, "ocean_surfstress_dv_stress")
      this%is_init = .true.
   end subroutine ocean_surfstress_init

   subroutine ocean_surfstress_destroy(this)
      class(ocean_surface_stress_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%tau_x)) deallocate (this%tau_x)
      if (allocated(this%tau_y)) deallocate (this%tau_y)
      if (allocated(this%stress_mag)) deallocate (this%stress_mag)
      if (allocated(this%stress_shelf)) deallocate (this%stress_shelf)
      call this%du_stress%destroy()
      call this%dv_stress%destroy()
   end subroutine ocean_surfstress_destroy

   subroutine ocean_surfstress_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! device-attach map base is the heap object, not a polymorphic stack
      !! box (AMD libomptarget cross-slot-overlap fix).  Slot-header
      !! presence comes from the orchestrator's root copyin(state); only
      !! the leaf arrays + scratch are attached here.
      class(ocean_surface_stress_t), intent(inout) :: this
      select type (this)
      type is (ocean_surface_stress_t)
         call ocean_surfstress_enter_data_impl(this)
      end select
   end subroutine ocean_surfstress_enter_data

   subroutine ocean_surfstress_enter_data_impl(this)
      type(ocean_surface_stress_t), intent(inout) :: this
      !$acc enter data copyin(this%tau_x, this%tau_y, this%stress_mag, &
      !$acc                   this%stress_shelf)
      ! Force the host wind values onto the device.  On OpenMP the root
      ! map(to:state) can leave tau_x/tau_y already "present" (descriptors
      ! come over with the parent), making the copyin above a no-op copy —
      ! `update device` (-> omp target update to) pushes the values
      ! regardless of presence.  Harmless on OpenACC.  `stress_mag` is
      ! filled host-side at configure by `ocean_surface_stress_set_derived`
      ! BEFORE this call — the `update device` here is what pushes it (the
      ! mem:separate trap: a missed push here gives `stress_mag == 0` on
      ! device, silently killing KPP/EPBL wind mixing, CLAUDE.md:312).
      !$acc update device(this%tau_x, this%tau_y, this%stress_mag, &
      !$acc                this%stress_shelf)
      call scratch_3d_buffer_enter_data_impl(this%du_stress)
      call scratch_3d_buffer_enter_data_impl(this%dv_stress)
   end subroutine ocean_surfstress_enter_data_impl

   subroutine ocean_surfstress_exit_data(this)
      class(ocean_surface_stress_t), intent(inout) :: this
      select type (this)
      type is (ocean_surface_stress_t)
         call ocean_surfstress_exit_data_impl(this)
      end select
   end subroutine ocean_surfstress_exit_data

   subroutine ocean_surfstress_exit_data_impl(this)
      type(ocean_surface_stress_t), intent(inout) :: this
      call scratch_3d_buffer_exit_data_impl(this%du_stress)
      call scratch_3d_buffer_exit_data_impl(this%dv_stress)
      !$acc exit data delete(this%tau_x, this%tau_y, this%stress_mag, &
      !$acc                  this%stress_shelf)
   end subroutine ocean_surfstress_exit_data_impl

   subroutine ocean_surfstress_set_const(this, tau_x_val, tau_y_val)
      !! Fill `tau_x` / `tau_y` uniformly with the given scalar
      !! values.  Convenience helper for the spatially-constant
      !! wind case (the Tier-1 default and most unit tests).  Host
      !! only — call `enter_data` afterwards (or `!$acc update
      !! device` if already mapped) to sync to GPU.  Refreshes
      !! `stress_mag` in step (PR-12) so every existing caller —
      !! production and unit test alike — gets a consistent
      !! `stress_mag` with no separate call required.
      class(ocean_surface_stress_t), intent(inout) :: this
      real(wp), intent(in) :: tau_x_val, tau_y_val
      this%tau_x = tau_x_val
      this%tau_y = tau_y_val
      call ocean_surfstress_refresh_stress_mag(this)
   end subroutine ocean_surfstress_set_const

   subroutine ocean_surfstress_set_2gyre(this, grid, taux_mag, j_offset, ny_global)
      !! Fill `tau_x` with the MOM6 2gyre profile,
      !! `tau_x(i,j) = taux_mag · (1 − cos(2π · (y − y_south) / y_len))`,
      !! and zero `tau_y`.  In Cartesian terms `(y − y_south) / y_len`
      !! is the normalised position from the south wall of the physical
      !! domain (0 at south, 1 at north), so the formula reduces to
      !! `taux_mag · (1 − cos(2π · ((j_phys − 0.5) / ny_phys)))` with
      !! `j_phys = j − nghost`.  Physical-interior rows only; ghost rows
      !! stay at zero so wall faces see no spurious stress.  Host only —
      !! call `enter_data` afterwards (or `!$acc update device` if
      !! already mapped).
      !!
      !! `j_offset` is the global physical-index offset of this rank's
      !! first physical row (= `decomp%j_start - 1`).  `ny_global` is
      !! the global meridional physical extent.  Both default to the
      !! single-rank values (`0` and `grid%ny_phys`), preserving
      !! byte-identical single-rank behaviour.
      class(ocean_surface_stress_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: taux_mag
      integer, intent(in), optional :: j_offset, ny_global
      real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)
      real(wp) :: y_rel
      integer :: i, j, j_phys, ng
      integer :: joff, nyg

      joff = 0
      if (present(j_offset)) joff = j_offset
      nyg = grid%ny_phys
      if (present(ny_global)) nyg = ny_global

      this%tau_x = 0.0_wp
      this%tau_y = 0.0_wp
      ng = grid%nghost
      do j = ng + 1, ng + grid%ny_phys
         j_phys = j - ng
         ! Global meridional normalisation: (j_phys + joff) gives the
         ! global physical row index; ny_global is the global extent.
         y_rel = (real(j_phys + joff, wp) - 0.5_wp)/real(nyg, wp)
         do i = 1, size(this%tau_x, 1)
            this%tau_x(i, j) = taux_mag*(1.0_wp - cos(TWO_PI*y_rel))
         end do
      end do
      call ocean_surfstress_refresh_stress_mag(this)
   end subroutine ocean_surfstress_set_2gyre

   subroutine ocean_surfstress_set_neverworld2(this, grid, taux_mag, j_offset, ny_global)
      !! Fill `tau_x` with the **Neverworld2** zonal wind-stress profile
      !! (Marques et al. 2022, GMD; MOM6-inspired) and zero `tau_y`.  τ_x is a
      !! 3-band piecewise function of the normalized meridional position
      !! `y = (j_phys − 0.5)/ny_phys ∈ [0,1]` (which equals MOM6's
      !! `(lat − south)/len_lat` on a uniform grid), scaled by the peak stress
      !! `taux_mag` (Pa), with `off = 0.02`:
      !!
      !!   band 1  y ≤ 0.29:           τ = taux·[ (1/0.29)·y − (1/2π)·sin(2π·y/0.29) ]
      !!   band 2  0.29 < y ≤ 0.8−off: τ = taux·[ 0.35 + 0.65·cos(π·(y−0.29)/(0.51−off)) ]
      !!   band 3  0.8−off < y ≤ 1−off: τ = taux·[ 1.5·((y−1+off) − (0.1/π)·sin(10π·(y−0.8+off))) ]
      !!   else (polar):               τ = 0
      !!
      !! This reproduces the canonical southern-westerlies / trades /
      !! polar-easterlies pattern.  τ_y ≡ 0.  Physical-interior rows only;
      !! ghost rows stay at zero so wall faces see no spurious stress (land
      !! masking is applied downstream via `wet_mask` in the stress-acceleration
      !! kernel).  Host only — call `enter_data` afterwards.
      !!
      !! `j_offset` is the global physical-index offset of this rank's
      !! first physical row (= `decomp%j_start - 1`).  `ny_global` is
      !! the global meridional physical extent.  Both default to the
      !! single-rank values, preserving byte-identical behaviour.
      class(ocean_surface_stress_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: taux_mag
      integer, intent(in), optional :: j_offset, ny_global
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp), parameter :: TWO_PI = 8.0_wp*atan(1.0_wp)
      real(wp), parameter :: OFF = 0.02_wp
      real(wp) :: y
      integer :: i, j, j_phys, ng
      integer :: joff, nyg

      joff = 0
      if (present(j_offset)) joff = j_offset
      nyg = grid%ny_phys
      if (present(ny_global)) nyg = ny_global

      this%tau_x = 0.0_wp
      this%tau_y = 0.0_wp
      ng = grid%nghost
      do j = ng + 1, ng + grid%ny_phys
         j_phys = j - ng
         ! Global meridional normalisation: (j_phys + joff) gives the
         ! global physical row index; nyg is the global extent.
         y = (real(j_phys + joff, wp) - 0.5_wp)/real(nyg, wp)
         do i = 1, size(this%tau_x, 1)
            if (y <= 0.29_wp) then
               this%tau_x(i, j) = taux_mag*((1.0_wp/0.29_wp)*y &
                                            - (1.0_wp/TWO_PI)*sin(TWO_PI*y/0.29_wp))
            else if (y <= 0.8_wp - OFF) then
               this%tau_x(i, j) = taux_mag*(0.35_wp + 0.65_wp*cos(PI*(y - 0.29_wp)/(0.51_wp - OFF)))
            else if (y <= 1.0_wp - OFF) then
               this%tau_x(i, j) = taux_mag*(1.5_wp*((y - 1.0_wp + OFF) &
                                                    - (0.1_wp/PI)*sin(10.0_wp*PI*(y - 0.8_wp + OFF))))
            end if
         end do
      end do
      call ocean_surfstress_refresh_stress_mag(this)
   end subroutine ocean_surfstress_set_neverworld2

   subroutine ocean_surface_stress_compute_tendencies(grid, this, ms)
      !! Fill `du_stress` / `dv_stress` with the surface stress
      !! acceleration `tau / (rho_0 * h_top)` at k = nz.  Outer-shim:
      !! hoist the derived-type derefs (`this%tau_x`, `this%du_stress%data`,
      !! `ms%h_layer`, `ms%wet_mask`) to the host, dispatch to flat-impl.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      if (this%direct_stress .and. this%hmix_stress > 0.0_wp) then
         call surfstress_distributed_impl( &
            ms%h_layer, ms%wet_mask, &
            this%tau_x, this%tau_y, &
            this%du_stress%data, this%dv_stress%data, &
            this%rho0, this%h_min, this%hmix_stress, &
            grid%nx_total, grid%ny_total, ms%nz_ml)
      else
         call surfstress_compute_impl( &
            ms%h_layer, ms%wet_mask, &
            this%tau_x, this%tau_y, &
            this%du_stress%data, this%dv_stress%data, &
            this%rho0, this%h_min, &
            grid%nx_total, grid%ny_total, ms%nz_ml)
      end if
   end subroutine ocean_surface_stress_compute_tendencies

   pure subroutine surfstress_distributed_impl(h_layer, wet_mask, tau_x, tau_y, &
                                               du_stress, dv_stress, &
                                               rho0, h_min, hmix_stress, &
                                               nx, ny, nz)
      !! DIRECT_STRESS branch — distribute the wind stress across
      !! the top `hmix_stress` metres of the column.  For each face,
      !! walk layers from k = nz (surface) down to k = 1, accumulating
      !! thickness; the surface boundary layer (SBL) is the set of
      !! layers whose top sits within `hmix_stress` of the free
      !! surface.  The acceleration per layer is:
      !!
      !!     du_k/dt = (tau / (rho0 · hmix_stress)) · (h_in_sbl_k / h_face_k)
      !!
      !! Sum over k: total impulse per unit area = tau / rho0,
      !! matching the bed-only formulation.  When `hmix_stress` is
      !! smaller than the top layer thickness the SBL is contained in
      !! the surface-most layer and the formula collapses to
      !! `tau / (rho0 · h_face_nz)` — bit-identical to the existing
      !! kernel.  When `hmix_stress` straddles multiple layers, each
      !! gets its fractional acceleration.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz), wet_mask(nx, ny)
      real(wp), intent(in)    :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      real(wp), intent(inout) :: du_stress(nx + 1, ny, nz), dv_stress(nx, ny + 1, nz)
      real(wp), intent(in)    :: rho0, h_min, hmix_stress
      integer :: i, j, k
      real(wp) :: inv_rho_hmix, cumul_h, h_face_k, h_in_sbl, mask_face

      inv_rho_hmix = 1.0_wp/(rho0*hmix_stress)

      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         du_stress(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         dv_stress(i, j, k) = 0.0_wp
      end do

      do concurrent(j=1:ny, i=2:nx) &
         local(k, cumul_h, h_face_k, h_in_sbl, mask_face)
         mask_face = min(wet_mask(i - 1, j), wet_mask(i, j))
         cumul_h = 0.0_wp
         do k = nz, 1, -1
            if (cumul_h >= hmix_stress) exit
            h_face_k = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            h_in_sbl = max(0.0_wp, min(h_face_k, hmix_stress - cumul_h))
            du_stress(i, j, k) = mask_face* &
                                 tau_x(i, j)*inv_rho_hmix* &
                                 (h_in_sbl/max(h_face_k, h_min))
            cumul_h = cumul_h + h_face_k
         end do
      end do
      do concurrent(j=2:ny, i=1:nx) &
         local(k, cumul_h, h_face_k, h_in_sbl, mask_face)
         mask_face = min(wet_mask(i, j - 1), wet_mask(i, j))
         cumul_h = 0.0_wp
         do k = nz, 1, -1
            if (cumul_h >= hmix_stress) exit
            h_face_k = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            if (h_face_k <= 0.0_wp) exit
            h_in_sbl = max(0.0_wp, min(h_face_k, hmix_stress - cumul_h))
            dv_stress(i, j, k) = mask_face* &
                                 tau_y(i, j)*inv_rho_hmix* &
                                 (h_in_sbl/max(h_face_k, h_min))
            cumul_h = cumul_h + h_face_k
         end do
      end do
   end subroutine surfstress_distributed_impl

   pure subroutine surfstress_compute_impl(h_layer, wet_mask, tau_x, tau_y, &
                                           du_stress, dv_stress, &
                                           rho0, h_min, nx, ny, nz)
      !! Flat-array surface-stress kernel.  Explicit-shape dummies so
      !! NVHPC stdpar can compile the device kernel against static
      !! bounds.
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz), wet_mask(nx, ny)
      real(wp), intent(in)    :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      real(wp), intent(inout) :: du_stress(nx + 1, ny, nz), dv_stress(nx, ny + 1, nz)
      real(wp), intent(in)    :: rho0, h_min
      integer :: i, j, k
      real(wp) :: inv_rho0, h_top_face

      inv_rho0 = 1.0_wp/rho0

      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         du_stress(i, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         dv_stress(i, j, k) = 0.0_wp
      end do

      ! Face wet-mask: only apply stress at faces between two ocean
      ! cells.  `min(wet_left, wet_right)` zeros stress at any face
      ! that touches land.  Default mask (all 1.0) preserves the
      ! flat-bottom / analytical behaviour bit-identically.
      do concurrent(j=1:ny, i=2:nx) local(h_top_face)
         h_top_face = 0.5_wp*(h_layer(i - 1, j, nz) + h_layer(i, j, nz))
         h_top_face = max(h_top_face, h_min)
         du_stress(i, j, nz) = &
            min(wet_mask(i - 1, j), wet_mask(i, j))* &
            tau_x(i, j)*inv_rho0/h_top_face
      end do
      do concurrent(j=2:ny, i=1:nx) local(h_top_face)
         h_top_face = 0.5_wp*(h_layer(i, j - 1, nz) + h_layer(i, j, nz))
         h_top_face = max(h_top_face, h_min)
         dv_stress(i, j, nz) = &
            min(wet_mask(i, j - 1), wet_mask(i, j))* &
            tau_y(i, j)*inv_rho0/h_top_face
      end do
   end subroutine surfstress_compute_impl

   subroutine ocean_surface_stress_apply_tendencies(this, ms, dt, no_wait)
      !! Outer-shim — flattens the derived-type derefs before the
      !! `do concurrent` body sees them.  Explicit-shape dimensions
      !! derived from `ms` and passed to the impl as scalar args.
      !! `no_wait` (optional, default .false.): forwarded to the impl —
      !! when .true. the apply DC loops run on OpenACC queue 1 without a
      !! trailing sync, so the batched velocity-apply chain in
      !! `run_stage_split` `!$acc wait(1)`s ONCE.  Default ⇒ blocking.
      !! Not `pure` because of the async/wait directives.
      type(ocean_surface_stress_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: no_wait
      integer :: nx_cells, ny_cells, nz
      logical :: lwait
      lwait = .true.
      if (present(no_wait)) lwait = .not. no_wait
      nx_cells = size(ms%u_face_x_layer, 1) - 1
      ny_cells = size(ms%v_face_y_layer, 2) - 1
      nz = ms%nz_ml
      call surfstress_apply_impl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                 this%du_stress%data, this%dv_stress%data, &
                                 dt, nx_cells, ny_cells, nz, lwait)
   end subroutine ocean_surface_stress_apply_tendencies

   subroutine surfstress_apply_impl(u_face, v_face, du_stress, dv_stress, &
                                    dt, nx, ny, nz, lwait)
      integer, intent(in)    :: nx, ny, nz
      real(wp), intent(inout) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in)    :: du_stress(nx + 1, ny, nz), dv_stress(nx, ny + 1, nz)
      real(wp), intent(in)    :: dt
      logical, intent(in)    :: lwait
         !! .false. ⇒ leave the apply on queue 1 without syncing (batched).
      integer :: i, j, k
      !$acc kernels async(1)
      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
         u_face(i, j, k) = u_face(i, j, k) + dt*du_stress(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
         v_face(i, j, k) = v_face(i, j, k) + dt*dv_stress(i, j, k)
      end do
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine surfstress_apply_impl

   pure subroutine ocean_surface_stress_set_shelf_from_ustar(stress_shelf, ustar, &
                                                             rho0, nx, ny)
      !! Publish the ice-base stress from a FRICTION VELOCITY:
      !! `stress_shelf = rho_0 * u_*^2`.
      !!
      !! The melt-only route into the `stress_shelf` seam.  When
      !! `&ocean_tdrag_nml` is on, the RK2 stage drivers fill
      !! `stress_shelf` inline from `ocean_top_drag_t%stress_top` and
      !! this is never called; when the top drag is OFF but
      !! `&ocean_cavity_melt_nml` is on, `engine_step_finalize` calls it
      !! with the melt slot's own `u_*`, which was solved with the SAME
      !! `C_d` (the one-drag-coefficient rule refuses a disagreeing
      !! `cdrag_top` at configure).  So this is a change of variable, not
      !! a second drag law.
      !!
      !! **Cadence, stated because it is a real cost:** the melt `u_*` is
      !! refreshed at the THERMO cadence at the END of an outer step, so
      !! this path reaches KPP/EPBL ONE OUTER STEP LATE.  The top-drag
      !! path has no such lag.
      !!
      !! A CALL rather than an inline loop at the call site, deliberately:
      !! `engine_step_finalize` would have to walk `engine%state%...`
      !! inside a `do concurrent`, and `ocean_engine_t` is not a mapped
      !! object, so nvfortran emits a data clause for the whole engine and
      !! aborts with "partially present on the device".  Host-dereference
      !! at the call site, flat explicit-shape dummies here.  The
      !! escaping-actual pessimisation CLAUDE.md warns about does not
      !! apply: `engine_step_finalize` owns no `do concurrent` of its own.
      !!
      !! `u_*` is exactly zero on every uncovered column
      !! (`cavity_melt_columns_2d`), so the published field keeps
      !! `stress_shelf`'s "exactly zero off the cover" invariant and the
      !! disjoint-support argument in the module docstring still holds.
      integer, intent(in) :: nx, ny
         !! Extents of BOTH arrays.  The caller gates on the melt slot's
         !! `enable`, which is exactly when `ustar` is full size, so an
         !! `(1,1)` placeholder can never reach these dummies.
      real(wp), intent(in) :: rho0
         !! Reference density (kg/m^3) -- the single rho0 of record, by
         !! value from `ocean_surface_stress_t%rho0`.
      real(wp), intent(inout) :: stress_shelf(nx, ny)
         !! `ocean_surface_stress_t%stress_shelf` (N/m^2), overwritten.
      real(wp), intent(in) :: ustar(nx, ny)
         !! `ocean_cavity_flux_t%ustar` (m/s).

      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         stress_shelf(i, j) = rho0*ustar(i, j)*ustar(i, j)
      end do
   end subroutine ocean_surface_stress_set_shelf_from_ustar

   pure function ocean_surface_stress_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the surface stress slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_surface_stress_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%tau_x) &
               + arr_bytes(this%tau_y) &
               + arr_bytes(this%stress_mag) &
               + arr_bytes(this%stress_shelf) &
               + this%du_stress%bytes() &
               + this%dv_stress%bytes()
   end function ocean_surface_stress_bytes

   pure subroutine ocean_surface_stress_apply_cover(ss, cover_frac)
      !! Zero the wind-stress pair on every C-grid face that touches an
      !! ice-covered cell, then refresh `stress_mag` from the masked
      !! pair.  The two halves are ONE call on purpose: a masked `tau`
      !! with a stale `stress_mag` would leave KPP/EPBL mixing on a wind
      !! that no longer reaches the water.
      !!
      !! FACE RULE — see the module docstring.  `open_face =
      !! min(1 - cover_left, 1 - cover_right)`, i.e. a face is closed if
      !! EITHER neighbour is covered, so no wind acts on a covered cell.
      !! A rim face (`i = 1` / `i = nx+1`, `j = 1` / `j = ny+1`) has only
      !! one neighbour in range and takes that cell's cover.
      !!
      !! Idempotent (multiplies by 0 or 1), so the configure-time call
      !! and the per-bracket data-forcing seam can both run it.
      !!
      !! **Runs where the data lives** — a plain `do concurrent`, so it
      !! is a host loop before `enter_data` and a device kernel after,
      !! exactly like `ocean_surface_stress_refresh_mag`.
      type(ocean_surface_stress_t), intent(inout) :: ss
      real(wp), intent(in) :: cover_frac(:, :)
         !! Ice-cover fraction at cell centres (`metrics%cover_frac`,
         !! v1 binary 0/1), shape `(nx_total, ny_total)`.
      ! assumed-shape-ok: configure-time / per-forcing-bracket cadence,
      ! forwarded to an explicit-shape `_impl` before any device loop.
      integer :: nx, ny

      if (.not. allocated(ss%tau_x) .or. .not. allocated(ss%tau_y)) return
      nx = size(ss%tau_x, 1) - 1
      ny = size(ss%tau_x, 2)
      ! A placeholder-sized cover (the `use_cavity = .false.` (1,1)
      ! allocation) is not a mask — refuse it rather than mask a corner.
      if (size(cover_frac, 1) /= nx .or. size(cover_frac, 2) /= ny) return
      call ocean_surfstress_cover_impl(ss%tau_x, ss%tau_y, cover_frac, nx, ny)
      call ocean_surface_stress_refresh_mag(ss)
   end subroutine ocean_surface_stress_apply_cover

   pure subroutine ocean_surfstress_cover_impl(tau_x, tau_y, cover_frac, nx, ny)
      !! Flat `do concurrent` kernel behind `ocean_surface_stress_apply_cover`
      !! — explicit-shape dummies, integer dims first (decl-order, ifx
      !! #8586).
      !!
      !! THIS IS `top_drag_fill_face_cover_impl`'S RULE, APPLIED IN PLACE.
      !! The face projection is stated once, in `rdb_ocean_top_drag`'s
      !! module docstring (§ The FACE cover rule), and there are two
      !! implementations of it only because one fills face ARRAYS the drag
      !! kernel reads per step and the other multiplies `tau` in place at
      !! configure — materialising two more face arrays here, in a slot
      !! that exists whether or not there is a cavity, to then throw them
      !! away, is the worse trade.  They are kept identical by a TEST
      !! (`cavity_cover_face_rule_matches_top_drag`), which is the same
      !! way `mirror_of_bottom_drag` keeps the top and bottom drag one
      !! closure rather than two.  Both halves of the rule matter:
      !!
      !!   * interior faces take the OR, `max(cover(i-1,j), cover(i,j))`,
      !!     so the CALVING-FRONT face is closed to the wind exactly where
      !!     it is opened to the drag;
      !!   * RIM faces (`i = 1`, `i = nx+1`, `j = 1`, `j = ny+1`) are left
      !!     UNTOUCHED, because they have only one neighbour in range.
      !!     Inventing a one-sided cover there would (a) disagree with the
      !!     drag's `cover_u = 0` and (b) clobber a periodic `tau` ghost
      !!     that the halo seam had already wrapped — this routine runs
      !!     AFTER `ocean_seam_refresh_surface_stress`.  Those faces carry
      !!     no prognostic velocity, so leaving them alone is not an
      !!     approximation.
      integer, intent(in)    :: nx, ny
      real(wp), intent(inout) :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      real(wp), intent(in)    :: cover_frac(nx, ny)
      integer :: i, j
      real(wp) :: cov
      do concurrent(j=1:ny, i=2:nx) local(cov)
         cov = max(cover_frac(i - 1, j), cover_frac(i, j))
         tau_x(i, j) = tau_x(i, j)*(1.0_wp - cov)
      end do
      do concurrent(j=2:ny, i=1:nx) local(cov)
         cov = max(cover_frac(i, j - 1), cover_frac(i, j))
         tau_y(i, j) = tau_y(i, j)*(1.0_wp - cov)
      end do
   end subroutine ocean_surfstress_cover_impl

   subroutine ocean_surface_stress_set_derived(grid, ss, cover_frac)
      !! Fill `stress_mag` from the current `tau_x`/`tau_y` — the MOM6
      !! `set_derived_forcing_fields` analogue (PR-12), and the public
      !! entry point `configure_ocean_forcing` calls after the
      !! `wind_config` dispatch.  Host-side (the wind field is
      !! configure-static in v1, so one fill at configure suffices — a
      !! future time-varying wind reader re-calls this after each read).
      !! Not `pure`: writes into `ss`.  Call BEFORE `enter_data`, or
      !! follow with `!$acc update device(ss%stress_mag)` if already
      !! mapped.  In practice this is a defensive re-fill only: every
      !! `set_wind_stress_*` setter already refreshes `stress_mag` in
      !! step (`ocean_surfstress_refresh_stress_mag`) so `stress_mag` is
      !! never stale relative to `tau_x`/`tau_y` regardless of call site
      !! (production driver OR a unit test that never reaches
      !! `configure_ocean_forcing`).
      type(hgrid_t), intent(in) :: grid
         !! Kept in the public signature (PR-12 plan §5.2) though shape
         !! is now taken from `tau_x`/`tau_y` directly — see
         !! `ocean_surfstress_refresh_stress_mag`.
      type(ocean_surface_stress_t), intent(inout) :: ss
      real(wp), intent(in), optional :: cover_frac(:, :)
         !! Ice-shelf cover fraction at cell centres
         !! (`metrics%cover_frac`).  Present ⇒ the `tau` pair is masked
         !! on every face touching a covered cell BEFORE `stress_mag` is
         !! rebuilt (`ocean_surface_stress_apply_cover`).  Absent ⇒ the
         !! original refresh-only path, byte-identical.
      ! assumed-shape-ok: configure-time / per-forcing-bracket cadence.
      if (present(cover_frac)) then
         call ocean_surface_stress_apply_cover(ss, cover_frac)
      else
         call ocean_surface_stress_refresh_mag(ss)
      end if
   end subroutine ocean_surface_stress_set_derived

   subroutine ocean_surfstress_refresh_stress_mag(this)
      !! Type-bound-facing shim behind every `set_wind_stress_*` setter:
      !! strips the polymorphic box (same reason as `enter_data`) and
      !! delegates to `ocean_surface_stress_refresh_mag`.  Keeping this
      !! call inside each setter — rather than requiring a separate
      !! explicit call — is what keeps `stress_mag` correct for every
      !! existing caller, including unit tests that build
      !! `ocean_surface_stress_t` directly and never reach
      !! `configure_ocean_forcing`.
      class(ocean_surface_stress_t), intent(inout) :: this
      select type (this)
      type is (ocean_surface_stress_t)
         call ocean_surface_stress_refresh_mag(this)
      end select
   end subroutine ocean_surfstress_refresh_stress_mag

   pure subroutine ocean_surface_stress_refresh_mag(ss)
      !! Recompute `stress_mag` from the CURRENT `tau_x`/`tau_y`, shape
      !! taken from the already-allocated arrays (no `grid` needed — this
      !! is the grid-free twin of `ocean_surface_stress_set_derived`, for
      !! callers that hold the slot but not the grid).
      !!
      !! **Runs where the data lives.**  `ocean_surfstress_derived_impl` is
      !! a plain `do concurrent`, so at configure time (before
      !! `enter_data`) this is a host loop over host arrays, and once the
      !! slot is device-mapped the SAME call is a device kernel over the
      !! mapped `tau_x`/`tau_y`/`stress_mag`.  A per-step caller therefore
      !! pays no host round trip and no allocation — which is what lets the
      !! sea-ice stress coupler (`ice_ocean_stress_flux`,
      !! `rdb_ice_ocean_coupler`) refresh `stress_mag` in step with the
      !! ice-mediated `tau` it writes on the device every outer step.
      !! Without that refresh KPP (`rdb_ocean_vmix`) and EPBL
      !! (`rdb_ocean_epbl`) — whose only source of `u_*` is this field —
      !! keep mixing on the configure-time WIND under sea ice.
      type(ocean_surface_stress_t), intent(inout) :: ss
      integer :: nx, ny
      nx = size(ss%tau_x, 1) - 1
      ny = size(ss%tau_x, 2)
      call ocean_surfstress_derived_impl(ss%tau_x, ss%tau_y, ss%stress_mag, nx, ny)
      ss%has_stress_mag = .true.
   end subroutine ocean_surface_stress_refresh_mag

   pure subroutine ocean_surfstress_derived_impl(tau_x, tau_y, stress_mag, nx, ny)
      !! `stress_mag(i,j) = |tau|` at cell centres — literal copy of the
      !! three lines this dedups from `rdb_ocean_vmix.F90` (KPP,
      !! `:575-577`/`:641-643`, pre-PR-12) and `rdb_ocean_epbl.F90`
      !! (`:974-976`, pre-PR-12): SAME face-average op order, so the
      !! substitution at each call site is bit-identical (PR-12 §7.5).
      !! Explicit-shape dummies, integer dims first (decl-order).
      integer, intent(in)    :: nx, ny
      real(wp), intent(in)    :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      real(wp), intent(inout) :: stress_mag(nx, ny)
      integer :: i, j
      real(wp) :: tau_x_cell, tau_y_cell
      do concurrent(j=1:ny, i=1:nx) local(tau_x_cell, tau_y_cell)
         tau_x_cell = 0.5_wp*(tau_x(i, j) + tau_x(i + 1, j))
         tau_y_cell = 0.5_wp*(tau_y(i, j) + tau_y(i, j + 1))
         stress_mag(i, j) = sqrt(tau_x_cell*tau_x_cell + tau_y_cell*tau_y_cell)
      end do
   end subroutine ocean_surfstress_derived_impl

end module rdb_ocean_surface_stress
