!! Unit tests for domain decomposition (no MPI required)
module test_decomp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_decomp, only: decomp_t, decomp_init, decomp_global_to_local, &
                         decomp_local_to_global, decomp_rank_from_coords
   implicit none
   private

   public :: collect_decomp_tests

contains

   subroutine collect_decomp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("local_sizes_sum_to_global", test_sizes_sum), &
                  new_unittest("remainder_distribution", test_remainder), &
                  new_unittest("boundary_flags_corners", test_boundary_corners), &
                  new_unittest("boundary_flags_interior", test_boundary_interior), &
                  new_unittest("coord_mapping_roundtrip", test_coord_roundtrip), &
                  new_unittest("single_rank_is_identity", test_single_rank), &
                  new_unittest("rank_from_coords", test_rank_coords) &
                  ]
   end subroutine collect_decomp_tests

   subroutine test_sizes_sum(error)
      !! Local nx sizes across all x-ranks must sum to nx_global
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      integer :: px, py, rx, sum_nx, sum_ny, ry

      px = 3
      py = 2
      ry = 0

      ! Sum local x-sizes across x-ranks (ry=0)
      sum_nx = 0
      do rx = 0, px - 1
         call decomp_init(d, 100, 50, px, py, ry*px + rx)
         sum_nx = sum_nx + d%nx_local
      end do
      call check(error, sum_nx == 100, "sum of nx_local must equal nx_global")
      if (allocated(error)) return

      ! Sum local y-sizes across y-ranks (rx=0)
      sum_ny = 0
      do ry = 0, py - 1
         call decomp_init(d, 100, 50, px, py, ry*px)
         sum_ny = sum_ny + d%ny_local
      end do
      call check(error, sum_ny == 50, "sum of ny_local must equal ny_global")

   end subroutine test_sizes_sum

   subroutine test_remainder(error)
      !! 10 cells over 3 ranks: first rank gets 4, others get 3
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d

      call decomp_init(d, 10, 1, 3, 1, 0)
      call check(error, d%nx_local == 4, "rank 0 should get 4 cells (10/3 + remainder)")
      if (allocated(error)) return

      call decomp_init(d, 10, 1, 3, 1, 1)
      call check(error, d%nx_local == 3, "rank 1 should get 3 cells")
      if (allocated(error)) return

      call decomp_init(d, 10, 1, 3, 1, 2)
      call check(error, d%nx_local == 3, "rank 2 should get 3 cells")
      if (allocated(error)) return

      ! i_start offsets: rank0=1, rank1=5, rank2=8
      call decomp_init(d, 10, 1, 3, 1, 0)
      call check(error, d%i_start == 1, "rank 0 i_start should be 1")
      if (allocated(error)) return

      call decomp_init(d, 10, 1, 3, 1, 1)
      call check(error, d%i_start == 5, "rank 1 i_start should be 5")
      if (allocated(error)) return

      call decomp_init(d, 10, 1, 3, 1, 2)
      call check(error, d%i_start == 8, "rank 2 i_start should be 8")

   end subroutine test_remainder

   subroutine test_boundary_corners(error)
      !! 2x2 grid: corner ranks should have exactly 2 boundary flags
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      integer :: n_bnd

      ! rank 0 = (rx=0, ry=0) -> has_west, has_south
      call decomp_init(d, 100, 100, 2, 2, 0)
      n_bnd = 0
      if (d%has_west) n_bnd = n_bnd + 1
      if (d%has_east) n_bnd = n_bnd + 1
      if (d%has_south) n_bnd = n_bnd + 1
      if (d%has_north) n_bnd = n_bnd + 1
      call check(error, n_bnd == 2, "corner rank should have 2 boundary flags")
      if (allocated(error)) return
      call check(error, d%has_west .and. d%has_south, &
                 "rank 0 should be west+south")
      if (allocated(error)) return

      ! rank 3 = (rx=1, ry=1) -> has_east, has_north
      call decomp_init(d, 100, 100, 2, 2, 3)
      call check(error, d%has_east .and. d%has_north, &
                 "rank 3 should be east+north")

   end subroutine test_boundary_corners

   subroutine test_boundary_interior(error)
      !! 3x3 grid, rank 4 (centre) should have no boundary flags
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d

      call decomp_init(d, 90, 90, 3, 3, 4)
      call check(error,.not. d%has_west, "interior rank should not have west")
      if (allocated(error)) return
      call check(error,.not. d%has_east, "interior rank should not have east")
      if (allocated(error)) return
      call check(error,.not. d%has_south, "interior rank should not have south")
      if (allocated(error)) return
      call check(error,.not. d%has_north, "interior rank should not have north")

   end subroutine test_boundary_interior

   subroutine test_coord_roundtrip(error)
      !! global_to_local followed by local_to_global should be identity
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      integer :: il, jl, ig_out, jg_out

      call decomp_init(d, 100, 50, 2, 2, 3)

      ! Pick a global point inside this rank's domain
      call decomp_global_to_local(d, d%i_start + 5, d%j_start + 3, il, jl)
      call decomp_local_to_global(d, il, jl, ig_out, jg_out)

      call check(error, ig_out == d%i_start + 5, "i roundtrip should match")
      if (allocated(error)) return
      call check(error, jg_out == d%j_start + 3, "j roundtrip should match")

   end subroutine test_coord_roundtrip

   subroutine test_single_rank(error)
      !! px=1, py=1 should give full global domain
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d

      call decomp_init(d, 200, 100, 1, 1, 0)
      call check(error, d%nx_local == 200, "single rank nx_local == nx_global")
      if (allocated(error)) return
      call check(error, d%ny_local == 100, "single rank ny_local == ny_global")
      if (allocated(error)) return
      call check(error, d%i_start == 1, "single rank i_start == 1")
      if (allocated(error)) return
      call check(error, d%j_start == 1, "single rank j_start == 1")
      if (allocated(error)) return
      call check(error, d%has_west .and. d%has_east .and. &
                 d%has_south .and. d%has_north, &
                 "single rank should have all boundaries")

   end subroutine test_single_rank

   subroutine test_rank_coords(error)
      !! rank_from_coords should give correct row-major rank
      type(error_type), allocatable, intent(out) :: error

      call check(error, decomp_rank_from_coords(3, 0, 0) == 0, "rank(0,0)=0")
      if (allocated(error)) return
      call check(error, decomp_rank_from_coords(3, 2, 0) == 2, "rank(2,0)=2")
      if (allocated(error)) return
      call check(error, decomp_rank_from_coords(3, 0, 1) == 3, "rank(0,1)=3")
      if (allocated(error)) return
      call check(error, decomp_rank_from_coords(3, 1, 1) == 4, "rank(1,1)=4")

   end subroutine test_rank_coords

end module test_decomp
