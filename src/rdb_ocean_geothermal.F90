!! Geothermal bottom heat flux for the ocean core — bed-side analogue
!! of `rdb_ocean_surface_flux`.  Stamps a constant bottom heat flux
!! `Q_geo` (W/m^2, positive into the ocean from below) into the bottom
!! layer (`k=1`, ROMS-style k=1 bed / k=nz surface).  Heat only.
!!
!! In `hTr` space (`hTr = T*h`, units K*m) the increment is
!! thickness-independent: `d(hT_{k=1}) = Q_geo * dt / (rho0 * cp)`.
!! Under `VCOORD_ZSTAR_FULL` the bed layer can pinch to near-zero, so
!! the increment lands in the lowest massive layer (first `k` with
!! `h_layer > h_min`), `k=1` in the common case.
module rdb_ocean_geothermal
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   implicit none
   private

   public :: ocean_geothermal_t
   public :: ocean_geothermal_apply_tracers

   type :: ocean_geothermal_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: enable = .false.
         !! Master switch.  Default `.false.` — the kernel no-ops, so
         !! existing nmls + tests stay bit-identical.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).
      real(wp) :: cp = SEAWATER_CP
         !! Specific heat capacity (J/kg/K).
      real(wp) :: h_min = 1.0e-3_wp
         !! Thickness floor for the lowest-massive-layer scan.
      real(wp) :: q_geo_const = 0.0_wp
         !! Scalar constant bottom heat flux (W/m^2, positive into the
         !! ocean from below).  Typical geothermal ~0.05-0.1 W/m^2.
   contains
      procedure :: init => ocean_geothermal_init
      procedure :: destroy => ocean_geothermal_destroy
   end type ocean_geothermal_t

contains

   subroutine ocean_geothermal_init(this, grid)
      !! No-op shell mirroring the surface-flux init; sets `is_init`.
      class(ocean_geothermal_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (.false.) this%rho0 = real(grid%nx_total, wp)
      this%is_init = .true.
   end subroutine ocean_geothermal_init

   subroutine ocean_geothermal_destroy(this)
      !! No-op shell mirroring the surface-flux destroy.
      class(ocean_geothermal_t), intent(inout) :: this
      this%is_init = .false.
   end subroutine ocean_geothermal_destroy

   subroutine ocean_geothermal_apply_tracers(grid, geo, ms, dt, active)
      !! Add the geothermal bottom heat flux to the lowest massive
      !! tracer layer.  Operates in `hTr` space (concentration*
      !! thickness):
      !!   d(hT_{k=1})/dt = Q_geo / (rho0 * cp)
      !! (units (W/m^2)/(kg/m^3 * J/kg/K) = K*m/s, matching `hTr`).
      !!
      !! No-op when `geo` is absent, `.not. enable`, `q_geo_const == 0`,
      !! no temperature tracer is registered, or `ms%tracers` is
      !! unallocated — preserving the default-off bit-identity contract.
      type(hgrid_t), intent(in) :: grid
      type(ocean_geothermal_t), intent(in), optional :: geo
         !! Optional — when absent the kernel is a no-op (no geothermal
         !! forcing configured).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional gate (thermo cadence).  Absent => kernel runs;
         !! present-and-false => early return.

      integer :: nx, ny, nz, idx_T
      real(wp) :: src_T

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. present(geo)) return
      if (.not. geo%enable) return
      if (geo%q_geo_const == 0.0_wp) return
      if (.not. allocated(ms%tracers)) return

      idx_T = ms%idx_temperature
      if (idx_T <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      src_T = dt*geo%q_geo_const/(geo%rho0*geo%cp)

      ! Shim+_impl split keeps `tracers(idx)%hTr` deref on the host
      ! (array-of-DT registry indirection blocks NVHPC device codegen).
      call apply_geothermal_src_impl(ms%tracers(idx_T)%hTr, &
                                     ms%heat_budget_geothermal, &
                                     ms%wet_mask, ms%h_layer, src_T, nz, &
                                     geo%h_min)
   end subroutine ocean_geothermal_apply_tracers

   pure subroutine apply_geothermal_src_impl(hTr, budget, wet_mask, h_layer, src, nz, h_min)
      !! Stamp `src * wet_mask(i,j)` onto the lowest *massive* layer of
      !! a tracer's hTr array (first `k` with `h_layer > h_min`,
      !! scanning `k = 1..nz` from the bed up), mirror into the matching
      !! budget contributor.  Flat-impl over plain allocatables — the
      !! outer subroutine reaches `ms%tracers(idx)%hTr` on the host
      !! before calling this.
      ! assumed-shape-ok: tracer registry outer-shim — caller host-dereferences
      ! ms%tracers(idx)%hTr before passing; size varies per tracer slot;
      ! called once per tracer per thermo step (per CLAUDE.md outer-shim pattern).
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(inout) :: budget(:, :, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
      real(wp), intent(in)    :: wet_mask(:, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
      real(wp), intent(in)    :: h_layer(:, :, :)  ! assumed-shape-ok: tracer registry outer-shim; thermo cadence
      real(wp), intent(in)    :: src
      integer, intent(in)    :: nz
      real(wp), intent(in)    :: h_min
      integer :: i, j, nx, ny, k, k_dep
      real(wp) :: cell
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      do concurrent(j=1:ny, i=1:nx) local(cell, k, k_dep)
         ! Lowest massive layer: scan from the bed (k=1) up.  In the
         ! common case h_layer(i,j,1) > h_min and k_dep = 1.  Falls back
         ! to nz if every layer is below the floor (deposits at the
         ! surface rather than dropping the energy).
         k_dep = nz
         do k = 1, nz
            if (h_layer(i, j, k) > h_min) then
               k_dep = k
               exit
            end if
         end do
         cell = src*wet_mask(i, j)
         hTr(i, j, k_dep) = hTr(i, j, k_dep) + cell
         budget(i, j, k_dep) = budget(i, j, k_dep) + cell
      end do
   end subroutine apply_geothermal_src_impl

end module rdb_ocean_geothermal
