!! Domain decomposition for MPI: splits a global structured grid across ranks.
module rdb_decomp
   !! Pure Fortran, no MPI dependency. Each rank calls decomp_init with its
   !! process-grid coordinates to learn its local subdomain size and offsets.
   use rdb_config, only: config_t
   implicit none
   private

   public :: decomp_t
   public :: decomp_init
   public :: decomp_init_from_config
   public :: decomp_global_to_local
   public :: decomp_local_to_global
   public :: decomp_rank_from_coords
   public :: decomp_auto_factor
   public :: decomp_log_summary

   type :: decomp_t
      !! Domain decomposition descriptor

      integer :: px = 1
         !! Number of processes in x-direction
      integer :: py = 1
         !! Number of processes in y-direction
      integer :: rx = 0
         !! This rank's x-coordinate in the process grid (0-based)
      integer :: ry = 0
         !! This rank's y-coordinate in the process grid (0-based)
      integer :: nx_global = 0
         !! Global physical cells in x
      integer :: ny_global = 0
         !! Global physical cells in y
      integer :: nx_local = 0
         !! Local physical cells in x for this rank
      integer :: ny_local = 0
         !! Local physical cells in y for this rank
      integer :: i_start = 1
         !! Global i-index of this rank's first physical cell (1-based)
      integer :: j_start = 1
         !! Global j-index of this rank's first physical cell (1-based)
      logical :: has_west = .true.
         !! True if this rank is on the west domain boundary
      logical :: has_east = .true.
         !! True if this rank is on the east domain boundary
      logical :: has_south = .true.
         !! True if this rank is on the south domain boundary
      logical :: has_north = .true.
         !! True if this rank is on the north domain boundary
   end type decomp_t

contains

   subroutine decomp_init(d, nx_global, ny_global, px, py, rank)
      !! Initialise decomposition for a given rank
      !!
      !! Distributes remainder cells to the first ranks in each direction.
      !! rank is mapped to (rx, ry) using row-major order: rank = ry * px + rx
      type(decomp_t), intent(out) :: d
         !! Decomposition descriptor to populate
      integer, intent(in) :: nx_global
         !! Total physical cells in x
      integer, intent(in) :: ny_global
         !! Total physical cells in y
      integer, intent(in) :: px
         !! Process grid size in x
      integer, intent(in) :: py
         !! Process grid size in y
      integer, intent(in) :: rank
         !! MPI rank (0-based)

      integer :: base_nx, rem_nx, base_ny, rem_ny

      d%px = px
      d%py = py
      d%nx_global = nx_global
      d%ny_global = ny_global

      ! Row-major mapping: rank = ry * px + rx
      d%rx = mod(rank, px)
      d%ry = rank/px

      ! Distribute cells with remainder going to first ranks
      base_nx = nx_global/px
      rem_nx = mod(nx_global, px)
      if (d%rx < rem_nx) then
         d%nx_local = base_nx + 1
         d%i_start = d%rx*(base_nx + 1) + 1
      else
         d%nx_local = base_nx
         d%i_start = rem_nx*(base_nx + 1) + (d%rx - rem_nx)*base_nx + 1
      end if

      base_ny = ny_global/py
      rem_ny = mod(ny_global, py)
      if (d%ry < rem_ny) then
         d%ny_local = base_ny + 1
         d%j_start = d%ry*(base_ny + 1) + 1
      else
         d%ny_local = base_ny
         d%j_start = rem_ny*(base_ny + 1) + (d%ry - rem_ny)*base_ny + 1
      end if

      ! Boundary flags
      d%has_west = (d%rx == 0)
      d%has_east = (d%rx == px - 1)
      d%has_south = (d%ry == 0)
      d%has_north = (d%ry == py - 1)

   end subroutine decomp_init

   subroutine decomp_init_from_config(d, cfg, nprocs, rank)
      !! Initialise decomposition from config, with optional auto-factoring
      !!
      !! If cfg%px and cfg%py are both 1 but nprocs > 1, automatically
      !! computes an optimal process grid. Updates cfg%px and cfg%py
      !! with the chosen values.
      type(decomp_t), intent(out) :: d
         !! Decomposition descriptor to populate
      type(config_t), intent(inout) :: cfg
         !! Configuration (px, py may be updated)
      integer, intent(in) :: nprocs
         !! Total number of MPI ranks
      integer, intent(in) :: rank
         !! This MPI rank (0-based)

      ! Auto-compute process grid if user didn't specify one
      if (nprocs > 1 .and. cfg%px == 1 .and. cfg%py == 1) then
         call decomp_auto_factor(nprocs, cfg%nx, cfg%ny, cfg%px, cfg%py)
      end if

      call decomp_init(d, cfg%nx, cfg%ny, cfg%px, cfg%py, rank)

   end subroutine decomp_init_from_config

   pure subroutine decomp_global_to_local(d, ig, jg, il, jl)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Convert global physical indices to local physical indices
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: ig
         !! Global i-index (1-based, physical cells)
      integer, intent(in) :: jg
         !! Global j-index (1-based, physical cells)
      integer, intent(out) :: il
         !! Local i-index (1-based, physical cells)
      integer, intent(out) :: jl
         !! Local j-index (1-based, physical cells)

      il = ig - d%i_start + 1
      jl = jg - d%j_start + 1

   end subroutine decomp_global_to_local

   pure subroutine decomp_local_to_global(d, il, jl, ig, jg)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Convert local physical indices to global physical indices
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: il
         !! Local i-index (1-based, physical cells)
      integer, intent(in) :: jl
         !! Local j-index (1-based, physical cells)
      integer, intent(out) :: ig
         !! Global i-index (1-based, physical cells)
      integer, intent(out) :: jg
         !! Global j-index (1-based, physical cells)

      ig = il + d%i_start - 1
      jg = jl + d%j_start - 1

   end subroutine decomp_local_to_global

   pure function decomp_rank_from_coords(px, rx, ry) result(rank)
      !! Compute rank from process grid coordinates (row-major)
      integer, intent(in) :: px
         !! Process grid size in x
      integer, intent(in) :: rx
         !! x-coordinate (0-based)
      integer, intent(in) :: ry
         !! y-coordinate (0-based)
      integer :: rank

      rank = ry*px + rx

   end function decomp_rank_from_coords

   pure subroutine decomp_auto_factor(nprocs, nx, ny, px, py)
      !! Choose px, py to minimise halo communication cost
      !!
      !! Tries all factorisations nprocs = px * py and picks the one
      !! that minimises the total halo perimeter: px*ny + py*nx.
      !! This produces roughly square subdomains matched to the grid
      !! aspect ratio.
      integer, intent(in) :: nprocs
         !! Total number of MPI ranks
      integer, intent(in) :: nx
         !! Global physical cells in x
      integer, intent(in) :: ny
         !! Global physical cells in y
      integer, intent(out) :: px
         !! Chosen process grid size in x
      integer, intent(out) :: py
         !! Chosen process grid size in y

      integer :: p, best_px
      integer :: cost, best_cost

      best_px = nprocs
      best_cost = huge(0)

      do p = 1, nprocs
         if (mod(nprocs, p) /= 0) cycle
         cost = p*ny + (nprocs/p)*nx
         if (cost < best_cost) then
            best_cost = cost
            best_px = p
         end if
      end do

      px = best_px
      py = nprocs/best_px

   end subroutine decomp_auto_factor

   subroutine decomp_log_summary(d, nprocs)
      !! Log the process grid + this rank's subdomain shape.  Call from
      !! rank 0 after decomp_init_from_config so a run's decomposition is
      !! visible at startup (the halo perimeter is what auto_factor
      !! minimises, printed here as the per-rank interior:ghost ratio cue).
      use pic_logger, only: logger => global_logger
      use pic_strings, only: to_string
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: nprocs
         !! Total compute ranks (px*py)

      character(len=16) :: shape_note

      if (d%px == 1) then
         shape_note = " (x-strips)"
      else if (d%py == 1) then
         shape_note = " (y-strips)"
      else
         shape_note = " (2-D tiles)"
      end if

      call logger%info("Domain decomposition: "//to_string(nprocs)// &
                       " ranks as "//to_string(d%px)//" x "//to_string(d%py)// &
                       " (px x py)"//trim(shape_note)//" over "// &
                       to_string(d%nx_global)//" x "//to_string(d%ny_global)//" cells")
      call logger%info("  rank-0 subdomain: "//to_string(d%nx_local)//" x "// &
                       to_string(d%ny_local)//" interior cells")
   end subroutine decomp_log_summary

end module rdb_decomp
