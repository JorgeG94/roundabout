!! Round-trip tests for the generic static 2-D NetCDF reader
!! (`load_static_2d_field`): a time-less `(x, y)` field over the whole
!! grid, read into a ghosted tile array.
!!
!! The oracle is an analytic field `f(i, j) = 1000*j + i` written by the test
!! itself (so a transposition, an off-by-one in the ghost offset or a wrong
!! variable would all change the value).  Covered: both on-disk orientations
!! (Fortran `(x, y)` and C-order `(y, x)`), dimension names `x/y` and
!! `nx/ny`, the constant-extrapolation ghosts, a missing variable, a wrong
!! grid size, and the unchanged 3-name bathymetry lookup.  Files go to
!! `tmp_local_artifacts/` (ctest's working directory is the build's tests/).
module test_static_2d_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_bathymetry, only: load_static_2d_field, load_bathymetry_into_array
   use rdb_io_netcdf, only: nc_create_file, nc_def_dim, nc_def_var_2d, nc_enddef, &
                            nc_put_var_2d, nc_close
   implicit none
   private

   public :: collect_static_2d_reader_tests

   integer, parameter :: NX = 7, NY = 5, NG = 2

contains

   subroutine collect_static_2d_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("static2d_xy_order", test_xy_order), &
                  new_unittest("static2d_yx_order_nxny_names", test_yx_order), &
                  new_unittest("static2d_missing_and_wrong_dims", test_failures), &
                  new_unittest("static2d_bathymetry_wrapper", test_bathy_wrapper) &
                  ]
   end subroutine collect_static_2d_reader_tests

   pure function analytic(i, j) result(f)
      integer, intent(in) :: i, j
      real(wp) :: f
      f = 1000.0_wp*real(j, wp) + real(i, wp)
   end function analytic

   subroutine write_field(filename, varname, nm_a, nm_b, transposed, nxf, nyf)
      !! Write `varname` = analytic field.  Not transposed: dims
      !! `(nm_a=nxf, nm_b=nyf)`; transposed: array `(nyf, nxf)` with dims
      !! `(nm_b, nm_a)` (what C/Python writers produce).
      character(len=*), intent(in) :: filename, varname, nm_a, nm_b
      logical, intent(in) :: transposed
      integer, intent(in) :: nxf, nyf
      integer :: ncid, d_a, d_b, vid, i, j
      real(wp), allocatable :: f(:, :), ft(:, :)

      allocate (f(nxf, nyf))
      do j = 1, nyf
         do i = 1, nxf
            f(i, j) = analytic(i, j)
         end do
      end do
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, nm_a, nxf, d_a)
      call nc_def_dim(ncid, nm_b, nyf, d_b)
      if (transposed) then
         call nc_def_var_2d(ncid, varname, [d_b, d_a], vid)
      else
         call nc_def_var_2d(ncid, varname, [d_a, d_b], vid)
      end if
      call nc_enddef(ncid)
      if (transposed) then
         allocate (ft(nyf, nxf))
         ft = transpose(f)
         call nc_put_var_2d(ncid, vid, ft)
      else
         call nc_put_var_2d(ncid, vid, f)
      end if
      call nc_close(ncid)
   end subroutine write_field

   subroutine check_loaded(error, b, grid, label)
      !! Interior equals the analytic field; ghosts equal the nearest
      !! interior cell (constant extrapolation).
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      character(len=*), intent(in) :: label
      integer :: i, j, ic, jc
      real(wp) :: want, worst

      worst = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            ic = min(max(i - NG, 1), NX)
            jc = min(max(j - NG, 1), NY)
            want = analytic(ic, jc)
            worst = max(worst, abs(b(i, j) - want))
         end do
      end do
      call check(error, worst == 0.0_wp, label//": loaded field differs from analytic")
   end subroutine check_loaded

   subroutine test_xy_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ierr
      call grid%init(NX, NY, NG, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      call write_field("tmp_local_artifacts/s2d_xy.nc", "tideamp", "x", "y", .false., NX, NY)
      call load_static_2d_field("tmp_local_artifacts/s2d_xy.nc", "tideamp", b, grid, ierr=ierr)
      call check(error, ierr == 0, "x/y order: reader returned an error")
      if (allocated(error)) return
      call check_loaded(error, b, grid, "x/y order")
   end subroutine test_xy_order

   subroutine test_yx_order(error)
      !! C-order file with the OM3 dimension names `ny`, `nx`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ierr
      call grid%init(NX, NY, NG, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      call write_field("tmp_local_artifacts/s2d_yx.nc", "h2", "nx", "ny", .true., NX, NY)
      call load_static_2d_field("tmp_local_artifacts/s2d_yx.nc", "h2", b, grid, ierr=ierr)
      call check(error, ierr == 0, "(ny,nx) order: reader returned an error")
      if (allocated(error)) return
      call check_loaded(error, b, grid, "(ny,nx) order")
   end subroutine test_yx_order

   subroutine test_failures(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid, gwrong
      real(wp), allocatable :: b(:, :), bw(:, :)
      integer :: ierr
      call grid%init(NX, NY, NG, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      call write_field("tmp_local_artifacts/s2d_fail.nc", "tideamp", "nx", "ny", .true., NX, NY)

      ierr = 0
      call load_static_2d_field("tmp_local_artifacts/s2d_fail.nc", "no_such_var", b, grid, ierr=ierr)
      call check(error, ierr /= 0, "missing variable must fail loud")
      if (allocated(error)) return

      ierr = 0
      call load_static_2d_field("tmp_local_artifacts/s2d_does_not_exist.nc", "tideamp", b, grid, ierr=ierr)
      call check(error, ierr /= 0, "missing file must fail loud")
      if (allocated(error)) return

      ierr = 0
      call load_static_2d_field("tmp_local_artifacts/s2d_fail.nc", "", b, grid, ierr=ierr)
      call check(error, ierr /= 0, "blank variable name must fail loud")
      if (allocated(error)) return

      ! Grid of a different size than the file.
      call gwrong%init(NX + 1, NY, NG, 1.0_wp, 1.0_wp)
      allocate (bw(gwrong%nx_total, gwrong%ny_total), source=-1.0_wp)
      ierr = 0
      call load_static_2d_field("tmp_local_artifacts/s2d_fail.nc", "tideamp", bw, gwrong, ierr=ierr)
      call check(error, ierr /= 0, "wrong grid size must fail loud")
   end subroutine test_failures

   subroutine test_bathy_wrapper(error)
      !! The bathymetry loader still finds `depth` (3-name lookup) in an
      !! `x/y` file and returns the same field as the generic reader.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), c(:, :)
      integer :: ierr
      call grid%init(NX, NY, NG, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      allocate (c(grid%nx_total, grid%ny_total), source=-2.0_wp)
      call write_field("tmp_local_artifacts/s2d_bathy.nc", "depth", "x", "y", .false., NX, NY)
      call load_bathymetry_into_array("tmp_local_artifacts/s2d_bathy.nc", b, grid, ierr=ierr)
      call check(error, ierr == 0, "bathymetry wrapper: reader returned an error")
      if (allocated(error)) return
      call check_loaded(error, b, grid, "bathymetry wrapper")
      if (allocated(error)) return
      call load_static_2d_field("tmp_local_artifacts/s2d_bathy.nc", "depth", c, grid, ierr=ierr)
      call check(error, ierr == 0 .and. all(b == c), "wrapper and generic reader disagree")
   end subroutine test_bathy_wrapper

end module test_static_2d_reader
