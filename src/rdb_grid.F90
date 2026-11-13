!! Structured grid metadata
module rdb_grid
   !! Pure metadata describing a structured Cartesian grid: physical and
   !! ghosted dimensions, ghost width, and cell sizes. No allocatables,
   !! no GPU mapping needed.
   use rdb_constants, only: wp
   implicit none
   private

   public :: hgrid_t

   type :: hgrid_t
      integer  :: nx_total = 0
         !! Total cells in x including ghosts (nx_phys + 2*nghost)
      integer  :: ny_total = 0
         !! Total cells in y including ghosts
      integer  :: nx_phys = 0
         !! Physical (interior) cells in x
      integer  :: ny_phys = 0
         !! Physical (interior) cells in y
      integer  :: nghost = 2
         !! Ghost cell width on each side
      real(wp) :: dx = 0.0_wp
         !! Cell size in x (m)
      real(wp) :: dy = 0.0_wp
         !! Cell size in y (m)
      integer  :: i_offset_global = 0
         !! Global index offset for this subdomain: i_global = i_phys_local + i_offset_global
         !! (0 on a single-rank / undecomposed grid — formula fills must add this
         !! when converting a local loop index to a physical coordinate)
      integer  :: j_offset_global = 0
         !! Global index offset in y (see i_offset_global)
      integer  :: nx_global = 0
         !! GLOBAL physical (interior) cells in x, i.e. the undecomposed
         !! domain extent.  Equals `nx_phys` on a single rank; under MPI it
         !! is the sum over the rank row.
         !!
         !! Why this lives on the grid: every formula fill (bathymetry
         !! setter, IC seeder) that normalises a position by the domain
         !! size, or centres a feature on the domain, needs the GLOBAL
         !! extent.  These used to take `nx_global`/`ny_global` (and the
         !! offsets) as OPTIONAL dummy arguments defaulting to the LOCAL
         !! `nx_phys`/`ny_phys` — so a fill that simply forgot to pass them
         !! compiled fine and silently rebuilt the whole pattern inside
         !! each rank's tile.  That defect shipped three times.  Carrying
         !! the global extents here, alongside `i_offset_global` /
         !! `j_offset_global`, makes the correct value the only value a
         !! fill can read: there is no argument left to forget.
      integer  :: ny_global = 0
         !! GLOBAL physical (interior) cells in y (see nx_global)
   contains
      procedure, non_overridable :: init => grid_init
   end type hgrid_t

contains

   subroutine grid_init(this, nx_phys, ny_phys, nghost, dx, dy)
      class(hgrid_t), intent(inout) :: this
      integer, intent(in) :: nx_phys, ny_phys, nghost
      real(wp), intent(in) :: dx, dy

      this%nx_phys = nx_phys
      this%ny_phys = ny_phys
      this%nghost = nghost
      this%dx = dx
      this%dy = dy
      this%nx_total = nx_phys + 2*nghost
      this%ny_total = ny_phys + 2*nghost
      ! Default the GLOBAL extents to the local ones — correct on a single
      ! rank / undecomposed grid.  The MPI driver overwrites both from
      ! `decomp%nx_global` / `decomp%ny_global` right after this call.
      this%nx_global = nx_phys
      this%ny_global = ny_phys
   end subroutine grid_init

end module rdb_grid
