!! Configure-time global reductions over a tiled grid.
module rdb_grid_reduce
   !! Rank-uniform verdicts at configure time need GLOBAL quantities: a test
   !! on a tile's own cells lets one rank refuse (or return before a
   !! collective) while its neighbours carry on, and the job then hangs or
   !! reduces garbage.  These wrappers reduce over the compute communicator
   !! ONLY when the grid is a tile of a larger domain (`nx_phys < nx_global`
   !! or `ny_phys < ny_global`) and are the identity otherwise.
   !!
   !! The identity case is load-bearing twice: unit tests build a state with
   !! no communicator initialised (a reduction would fault), and the
   !! bit-identity gate runs a serial reference (`compute_size = 1`, a global
   !! grid) on EVERY physical rank, where a reduction over the real
   !! communicator would sum `nprocs` copies of one problem.
   !!
   !! The tiled test depends only on the decomposition, which every rank
   !! agrees on, so all ranks take the same branch.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_halo, only: halo_allreduce_max, halo_allreduce_min, halo_allreduce_sum_i8
   implicit none
   private

   public :: grid_is_tiled, grid_allreduce_max, grid_allreduce_min, grid_allreduce_sum_i8

contains

   pure function grid_is_tiled(grid) result(tiled)
      !! True when this rank holds only part of the global domain.
      type(hgrid_t), intent(in) :: grid
      logical :: tiled
      tiled = grid%nx_phys < grid%nx_global .or. grid%ny_phys < grid%ny_global
   end function grid_is_tiled

   subroutine grid_allreduce_max(grid, local_val, global_val)
      !! Global max over the tiles of `grid` (identity on an untiled grid).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: local_val
      real(wp), intent(out) :: global_val
      if (grid_is_tiled(grid)) then
         call halo_allreduce_max(local_val, global_val)
      else
         global_val = local_val
      end if
   end subroutine grid_allreduce_max

   subroutine grid_allreduce_min(grid, local_val, global_val)
      !! Global min over the tiles of `grid` (identity on an untiled grid).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: local_val
      real(wp), intent(out) :: global_val
      if (grid_is_tiled(grid)) then
         call halo_allreduce_min(local_val, global_val)
      else
         global_val = local_val
      end if
   end subroutine grid_allreduce_min

   subroutine grid_allreduce_sum_i8(grid, local_val, global_val)
      !! Exact global integer sum over the tiles of `grid` (identity on an
      !! untiled grid).
      type(hgrid_t), intent(in) :: grid
      integer(int64), intent(in) :: local_val
      integer(int64), intent(out) :: global_val
      if (grid_is_tiled(grid)) then
         call halo_allreduce_sum_i8(local_val, global_val)
      else
         global_val = local_val
      end if
   end subroutine grid_allreduce_sum_i8

end module rdb_grid_reduce
