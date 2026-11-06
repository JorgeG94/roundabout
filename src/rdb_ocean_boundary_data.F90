!! Abstract ocean boundary data source.
module rdb_ocean_boundary_data
   !! Polymorphic interface for populating `ocean_bc_state_t`'s per-edge
   !! `data_*` buffers each outer step via a single `update(t, bc)` entry
   !! point. Concrete backends (file, callback, tidal table, constant)
   !! handle the data plumbing; the driver invokes `update` once per step.
   use rdb_constants, only: wp
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   implicit none
   private

   public :: ocean_boundary_data_source_t
   public :: ocean_boundary_data_constant_t

   type, abstract :: ocean_boundary_data_source_t
      !! Polymorphic base. Concrete backends extend this and provide
      !! `update` + `destroy`.
   contains
      procedure(boundary_data_update_iface), deferred :: update
      procedure(boundary_data_destroy_iface), deferred :: destroy
   end type ocean_boundary_data_source_t

   abstract interface
      subroutine boundary_data_update_iface(this, t, bc)
         !! Refresh `bc%data_*` (and any other time-varying fields)
         !! for the current outer-step wall time `t` (s).
         import :: wp, ocean_boundary_data_source_t, ocean_bc_state_t
         implicit none
         class(ocean_boundary_data_source_t), intent(inout) :: this
         real(wp), intent(in) :: t
         type(ocean_bc_state_t), intent(inout) :: bc
      end subroutine boundary_data_update_iface

      subroutine boundary_data_destroy_iface(this)
         import :: ocean_boundary_data_source_t
         implicit none
         class(ocean_boundary_data_source_t), intent(inout) :: this
      end subroutine boundary_data_destroy_iface
   end interface

   type, extends(ocean_boundary_data_source_t) :: ocean_boundary_data_constant_t
      !! Trivial backend: writes user-supplied constant scalars straight
      !! to `bc%west/east/south/north`'s `clamped_*` fields. For
      !! fixed-inflow tests and as a polymorphic-dispatch sanity check.
      real(wp) :: u_west = 0.0_wp, u_east = 0.0_wp
      real(wp) :: v_south = 0.0_wp, v_north = 0.0_wp
      real(wp) :: eta_west = 0.0_wp, eta_east = 0.0_wp
      real(wp) :: eta_south = 0.0_wp, eta_north = 0.0_wp
   contains
      procedure :: update => constant_update
      procedure :: destroy => constant_destroy
   end type ocean_boundary_data_constant_t

contains

   subroutine constant_update(this, t, bc)
      class(ocean_boundary_data_constant_t), intent(inout) :: this
      real(wp), intent(in) :: t
      type(ocean_bc_state_t), intent(inout) :: bc
      bc%west%clamped_u = this%u_west
      bc%east%clamped_u = this%u_east
      bc%south%clamped_v = this%v_south
      bc%north%clamped_v = this%v_north
      bc%west%clamped_eta = this%eta_west
      bc%east%clamped_eta = this%eta_east
      bc%south%clamped_eta = this%eta_south
      bc%north%clamped_eta = this%eta_north
      ! `t` unused for a constant source; kept to match the interface.
      if (.false.) then
         bc%west%clamped_u = t
         bc%west%clamped_u = this%u_west
      end if
   end subroutine constant_update

   subroutine constant_destroy(this)
      class(ocean_boundary_data_constant_t), intent(inout) :: this
      this%u_west = 0.0_wp
      this%u_east = 0.0_wp
      this%v_south = 0.0_wp
      this%v_north = 0.0_wp
      this%eta_west = 0.0_wp
      this%eta_east = 0.0_wp
      this%eta_south = 0.0_wp
      this%eta_north = 0.0_wp
   end subroutine constant_destroy

end module rdb_ocean_boundary_data
