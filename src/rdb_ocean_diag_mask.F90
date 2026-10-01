!! Region masks for ocean diagnostics.
!!
!! A mask is a per-cell weight in `[0, 1]` over the structured grid.
!! Diagnostic vars carry an optional `diag_mask_t`; when present, the
!! manager multiplies each sample by the mask weight at fold time so
!! cells outside the region contribute zero to the accumulator.
!! Builders: `diag_mask_global` (whole grid), `diag_mask_bbox`
!! (index rectangle), `diag_mask_h_section` / `diag_mask_v_section`
!! (single row/column strips for transport across a constant line).
!!
!! `total_area` caches the NOMINAL `Σ weight · dx · dy`; the physical
!! area-weighted integrals (budgets, console stats) do NOT use it — they
!! weight each cell by metric `areaT(i,j)` directly (`dx·dy` is dead on
!! curvilinear grids). On uniform Cartesian the two agree.
module rdb_ocean_diag_mask
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private

   public :: diag_mask_t
   public :: diag_mask_global, diag_mask_bbox
   public :: diag_mask_h_section, diag_mask_v_section
   public :: diag_mask_destroy

   type :: diag_mask_t
      character(len=64) :: name = ""
         !! Human-readable identifier ("global", "tasman_box", ...).
      integer :: nx = 0, ny = 0
         !! Mask shape — matches the grid's total extent.
      real(wp), allocatable :: weight(:, :)
         !! Per-cell weight in `[0, 1]`. 0 = excluded, 1 = full,
         !! fractional = partial coverage.
      real(wp) :: total_area = 0.0_wp
         !! Cached nominal `Σ weight · dx · dy`.
   contains
      procedure, non_overridable :: bytes => diag_mask_bytes
   end type diag_mask_t

contains

   pure function diag_mask_bytes(this) result(nbytes)
      !! Counted allocatable footprint of one region mask (0 when
      !! unallocated).  `weight` is device-mapped by
      !! `ocean_diag_enter_data_impl`, so it must carry a term — see
      !! `diag_var_bytes`, which folds this in per masked diagnostic.
      class(diag_mask_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%weight)
   end function diag_mask_bytes

   subroutine recompute_total_area(this, dx, dy)
      type(diag_mask_t), intent(inout) :: this
      real(wp), intent(in) :: dx, dy
      this%total_area = sum(this%weight)*dx*dy
   end subroutine recompute_total_area

   function diag_mask_global(grid) result(m)
      !! Whole-grid mask — all weights = 1. Equivalent to "no mask";
      !! exists for explicit registration and global-integral tests.
      type(hgrid_t), intent(in) :: grid
      type(diag_mask_t) :: m
      m%name = "global"
      m%nx = grid%nx_total
      m%ny = grid%ny_total
      allocate (m%weight(m%nx, m%ny), source=1.0_wp)
      call recompute_total_area(m, grid%dx, grid%dy)
   end function diag_mask_global

   function diag_mask_bbox(grid, i0, i1, j0, j1, name) result(m)
      !! Index-based bounding box (test-only). Weight = 1 inside the
      !! closed interval `[i0,i1] × [j0,j1]`, 0 outside. Indices 1-based
      !! and inclusive; ghost cells count if inside the box.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: i0, i1, j0, j1
      character(len=*), intent(in), optional :: name
      type(diag_mask_t) :: m
      integer :: i, j
      m%name = "bbox"
      if (present(name)) m%name = name
      m%nx = grid%nx_total
      m%ny = grid%ny_total
      allocate (m%weight(m%nx, m%ny), source=0.0_wp)
      do j = max(1, j0), min(m%ny, j1)
         do i = max(1, i0), min(m%nx, i1)
            m%weight(i, j) = 1.0_wp
         end do
      end do
      call recompute_total_area(m, grid%dx, grid%dy)
   end function diag_mask_bbox

   function diag_mask_h_section(grid, j_row, i0, i1, name) result(m)
      !! Horizontal section (test-only): one cell row at `j = j_row`,
      !! spanning `i ∈ [i0,i1]`. Integrates transport across a constant-y
      !! line (`v_face_y_layer` at `j_row` = meridional throughflow).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: j_row, i0, i1
      character(len=*), intent(in), optional :: name
      type(diag_mask_t) :: m
      integer :: i
      m%name = "h_section"
      if (present(name)) m%name = name
      m%nx = grid%nx_total
      m%ny = grid%ny_total
      allocate (m%weight(m%nx, m%ny), source=0.0_wp)
      if (j_row >= 1 .and. j_row <= m%ny) then
         do i = max(1, i0), min(m%nx, i1)
            m%weight(i, j_row) = 1.0_wp
         end do
      end if
      call recompute_total_area(m, grid%dx, grid%dy)
   end function diag_mask_h_section

   function diag_mask_v_section(grid, i_col, j0, j1, name) result(m)
      !! Vertical section (test-only): one cell column at `i = i_col`,
      !! spanning `j ∈ [j0,j1]`. Mirror of `h_section` for zonal throughflow.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: i_col, j0, j1
      character(len=*), intent(in), optional :: name
      type(diag_mask_t) :: m
      integer :: j
      m%name = "v_section"
      if (present(name)) m%name = name
      m%nx = grid%nx_total
      m%ny = grid%ny_total
      allocate (m%weight(m%nx, m%ny), source=0.0_wp)
      if (i_col >= 1 .and. i_col <= m%nx) then
         do j = max(1, j0), min(m%ny, j1)
            m%weight(i_col, j) = 1.0_wp
         end do
      end if
      call recompute_total_area(m, grid%dx, grid%dy)
   end function diag_mask_v_section

   subroutine diag_mask_destroy(this)
      type(diag_mask_t), intent(inout) :: this
      if (allocated(this%weight)) deallocate (this%weight)
      this%nx = 0
      this%ny = 0
      this%total_area = 0.0_wp
      this%name = ""
   end subroutine diag_mask_destroy

end module rdb_ocean_diag_mask
