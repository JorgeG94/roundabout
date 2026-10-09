!! Ocean open-boundary nesting state.
module rdb_ocean_obc
   !! Parent-state ingestion + interpolation state for nesting a regional
   !! ocean run inside a global/basin parent. Distinct from coastal
   !! `BC_NESTED`: parent state on a coarser mesh (ingest + remap onto our
   !! vertical grid), tides composed in via `rdb_ocean_tides`, mandatory
   !! T/S clamp, and a Flow-Relaxation-Scheme (FRS) zone of `nrelax` cells
   !! rather than a hard clamp. Scaffold (kernels not yet implemented).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_obc_t

   ! Edge-band tags.
   integer, parameter, public :: OBC_EDGE_WEST = 1
   integer, parameter, public :: OBC_EDGE_EAST = 2
   integer, parameter, public :: OBC_EDGE_SOUTH = 3
   integer, parameter, public :: OBC_EDGE_NORTH = 4

   type :: ocean_obc_t
      logical :: is_init = .false.
         !! True between `init` and `destroy` (tracks GPU device attachment).

      ! ---- Master switches ----
      logical :: enabled = .false.
         !! Master OBC switch.  Falls back to coastal BCs if false.
      logical :: clamp_tracers = .true.
         !! Apply the FRS clamp to T, S (not just SSH/UV).

      ! ---- FRS zone ----
      integer  :: nrelax = 8
         !! Width (in cells) of the Flow-Relaxation zone.
      real(wp) :: relax_inner = 0.0_wp
         !! Relaxation timescale (1/s) at the interior edge of FRS.
      real(wp) :: relax_outer = 1.0_wp/600.0_wp
         !! Relaxation timescale at the outer (boundary) edge.

      ! ---- Per-edge parent-state ring buffers ----
      ! Two-time-level buffers. Shape per edge:
      !   (nx, nz_ml [+1 for face arrays], 2 time levels)
      real(wp), allocatable :: parent_eta_west(:, :)
      real(wp), allocatable :: parent_eta_east(:, :)
      real(wp), allocatable :: parent_eta_south(:, :)
      real(wp), allocatable :: parent_eta_north(:, :)

      real(wp), allocatable :: parent_u_west(:, :, :)
      real(wp), allocatable :: parent_u_east(:, :, :)
      real(wp), allocatable :: parent_v_south(:, :, :)
      real(wp), allocatable :: parent_v_north(:, :, :)

      real(wp), allocatable :: parent_T_west(:, :, :)
      real(wp), allocatable :: parent_T_east(:, :, :)
      real(wp), allocatable :: parent_T_south(:, :, :)
      real(wp), allocatable :: parent_T_north(:, :, :)

      real(wp), allocatable :: parent_S_west(:, :, :)
      real(wp), allocatable :: parent_S_east(:, :, :)
      real(wp), allocatable :: parent_S_south(:, :, :)
      real(wp), allocatable :: parent_S_north(:, :, :)

      ! ---- Ingest timing ----
      real(wp) :: dt_ingest = 3600.0_wp
         !! Cadence at which the next parent slab is read (s).
      real(wp) :: t_next_ingest = 0.0_wp
   contains
      procedure, non_overridable :: init => ocean_obc_init
      procedure, non_overridable :: destroy => ocean_obc_destroy
      procedure, non_overridable :: bytes => ocean_obc_bytes
   end type ocean_obc_t

contains

   subroutine ocean_obc_init(this, grid)
      class(ocean_obc_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (.false.) this%nrelax = grid%nx_total
      this%is_init = .true.
   end subroutine ocean_obc_init

   subroutine ocean_obc_destroy(this)
      class(ocean_obc_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%parent_eta_west)) deallocate (this%parent_eta_west)
      if (allocated(this%parent_eta_east)) deallocate (this%parent_eta_east)
      if (allocated(this%parent_eta_south)) deallocate (this%parent_eta_south)
      if (allocated(this%parent_eta_north)) deallocate (this%parent_eta_north)
      if (allocated(this%parent_u_west)) deallocate (this%parent_u_west)
      if (allocated(this%parent_u_east)) deallocate (this%parent_u_east)
      if (allocated(this%parent_v_south)) deallocate (this%parent_v_south)
      if (allocated(this%parent_v_north)) deallocate (this%parent_v_north)
      if (allocated(this%parent_T_west)) deallocate (this%parent_T_west)
      if (allocated(this%parent_T_east)) deallocate (this%parent_T_east)
      if (allocated(this%parent_T_south)) deallocate (this%parent_T_south)
      if (allocated(this%parent_T_north)) deallocate (this%parent_T_north)
      if (allocated(this%parent_S_west)) deallocate (this%parent_S_west)
      if (allocated(this%parent_S_east)) deallocate (this%parent_S_east)
      if (allocated(this%parent_S_south)) deallocate (this%parent_S_south)
      if (allocated(this%parent_S_north)) deallocate (this%parent_S_north)
   end subroutine ocean_obc_destroy

   pure function ocean_obc_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the open boundary condition slot (0 when
      !! unallocated).
      class(ocean_obc_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%parent_eta_west) &
               + arr_bytes(this%parent_eta_east) &
               + arr_bytes(this%parent_eta_south) &
               + arr_bytes(this%parent_eta_north) &
               + arr_bytes(this%parent_u_west) &
               + arr_bytes(this%parent_u_east) &
               + arr_bytes(this%parent_v_south) &
               + arr_bytes(this%parent_v_north) &
               + arr_bytes(this%parent_T_west) &
               + arr_bytes(this%parent_T_east) &
               + arr_bytes(this%parent_T_south) &
               + arr_bytes(this%parent_T_north) &
               + arr_bytes(this%parent_S_west) &
               + arr_bytes(this%parent_S_east) &
               + arr_bytes(this%parent_S_south) &
               + arr_bytes(this%parent_S_north)
   end function ocean_obc_bytes

end module rdb_ocean_obc
