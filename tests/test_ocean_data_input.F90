!! Tests for the shared time-varying NetCDF input reader
!! (`rdb_ocean_data_input`, PR-14).
!!
!! Tests (NetCDF I/O — must run with OMP_NUM_THREADS=1):
!!   data_input_linear_in_time     — T1: f(x,y,t) = a+b*x+c*y+d*t, queried
!!                                   OFF-GRID, recovered to 1e-13 relative.
!!   data_input_linear_in_time_3d  — T1's 3-D twin, f(x,y,z,t).
!!   data_input_on_grid_exact      — T2: query at exactly a node returns
!!                                   w=0 on THAT record, not w=1 on the
!!                                   previous bracket.
!!   data_input_cyclic_wraparound  — T3: the (nt,1) seam bracket + the
!!                                   second-cycle bit-identity.
!!   data_input_static_mode        — T4: DATA_TIME_STATIC reads once,
!!                                   ever; two widely-separated queries
!!                                   return the same value.
!!   data_input_scale_offset       — T5: scale/add_offset applied once at
!!                                   read, recovered at two steps in the
!!                                   same bracket.
!!   data_input_oor_behaviour      — T6: DATA_OOR_CLAMP does not abort and
!!                                   holds the end record; the ERROR
!!                                   default's abort condition is checked
!!                                   via the pure `out_of_range` flag (no
!!                                   in-process death test — house
!!                                   convention, see `zinit`).
!!   data_input_bad_mode_aborts    — T7: `data_input_time_mode_is_implemented`
!!                                   false for an unknown tag (the aborting
!!                                   `..._from_string` itself is not
!!                                   invoked — same convention).
!!   data_input_dim_mismatch       — T8: `data_input_dims_ok` predicate,
!!                                   mirrors `zinit_dims_ok`.
!!   data_input_locate_pure        — T9: `data_input_locate` swept across
!!                                   LINEAR/CYCLIC x interior/node/seam/
!!                                   out-of-range, incl. nt==1.
!!   data_input_gpu_resident       — T10: register + enter_data + TWO
!!                                   queries spanning a bracket advance,
!!                                   all within ONE mapped session — the
!!                                   one test that can only pass if
!!                                   `!$acc update device` fires on the
!!                                   advance (CLAUDE.md mem:separate
!!                                   gotcha).
!!   data_input_fill_static_host   — fill_static_host (host-only,
!!                                   pre-map) reproduces update_2d's
!!                                   (device) result bit-for-bit.
!!
!! Bit-identity note: the reader is only exercised when a consumer
!! registers a field; PR-14 ships no consumer, so the full `ctest -R
!! rdb` suite staying green with zero registrations is the
!! bit-identity guarantee (see `ocean_data_input_update_all`'s
!! `nfields == 0` early return).
module test_ocean_data_input
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_data_input, only: ocean_data_input_t, &
                                   ocean_data_input_register_2d, &
                                   ocean_data_input_register_3d, &
                                   ocean_data_input_update_2d, &
                                   ocean_data_input_update_3d, &
                                   ocean_data_input_update_all, &
                                   ocean_data_input_fill_static_host, &
                                   data_input_locate, &
                                   data_input_time_mode_is_implemented, &
                                   data_input_dims_ok, &
                                   DATA_TIME_LINEAR, DATA_TIME_CYCLIC, DATA_TIME_STATIC, &
                                   DATA_OOR_ERROR, DATA_OOR_CLAMP
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            nc_def_var_4d, rdb_def_var_1d, nc_enddef, rdb_put_var_1d, &
                            nc_put_att
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_data_input_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_data_input_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("data_input_linear_in_time", test_linear_in_time), &
                  new_unittest("data_input_linear_in_time_3d", test_linear_in_time_3d), &
                  new_unittest("data_input_on_grid_exact", test_on_grid_exact), &
                  new_unittest("data_input_cyclic_wraparound", test_cyclic_wraparound), &
                  new_unittest("data_input_static_mode", test_static_mode), &
                  new_unittest("data_input_scale_offset", test_scale_offset), &
                  new_unittest("data_input_oor_behaviour", test_oor_behaviour), &
                  new_unittest("data_input_bad_mode_aborts", test_bad_mode), &
                  new_unittest("data_input_dim_mismatch", test_dim_mismatch), &
                  new_unittest("data_input_locate_pure", test_locate_pure), &
                  new_unittest("data_input_gpu_resident", test_gpu_resident), &
                  new_unittest("data_input_fill_static_host", test_fill_static_host) &
                  ]
   end subroutine collect_ocean_data_input_tests

   ! =================================================================
   ! Test-support: synthetic Fortran-ordered (x, y[, z], time) files.
   ! =================================================================

   subroutine write_2d_time_file(filename, nx, ny, nt, t_axis, f, tunits)
      !! Write a `field(x, y, time)` NetCDF, Fortran dim order (no
      !! transpose needed on read back — mirrors `write_zfile`).
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nt
      real(wp), intent(in) :: t_axis(nt), f(nx, ny, nt)
      character(len=*), intent(in), optional :: tunits

      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_3d(ncid, "field", [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      if (present(tunits)) call nc_put_att(ncid, vid_t, "units", tunits)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_2d_time_file

   subroutine write_3d_time_file(filename, nx, ny, nz, nt, t_axis, f)
      !! `field(x, y, z, time)` twin of `write_2d_time_file`.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nz, nt
      real(wp), intent(in) :: t_axis(nt), f(nx, ny, nz, nt)

      integer :: ncid, dim_x, dim_y, dim_z, dim_t, vid_f, vid_t, ierr

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "z", nz, dim_z)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_4d(ncid, "field", [dim_x, dim_y, dim_z, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_3d_time_file

   function make_grid(nx, ny) result(g)
      integer, intent(in) :: nx, ny
      type(hgrid_t) :: g
      call g%init(nx, ny, NGHOST, 1.0_wp, 1.0_wp)
   end function make_grid

   ! =================================================================
   ! T1 — linear-in-time, off-grid query (2-D).
   ! =================================================================

   subroutine test_linear_in_time(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 2, NT = 6
      real(wp), parameter :: A = 1.0_wp, B = 2.0_wp, C = 3.0_wp, D = 0.5_wp
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest(NX, NY)
      integer :: i, j, k, id
      real(wp) :: tq, expect

      do k = 1, NT
         t_axis(k) = real(k - 1, wp)*100.0_wp
      end do
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = A + B*real(i, wp) + C*real(j, wp) + D*t_axis(k)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_linear.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)

      grid = make_grid(NX, NY)
      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id)

      tq = 137.0_wp
      dest = -999.0_wp
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_2d(reader, id, tq, NX, NY, dest)
      !$acc update self(dest)
      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = 1, NY
         do i = 1, NX
            expect = A + B*real(i, wp) + C*real(j, wp) + D*tq
            call check(error, abs(dest(i, j) - expect) < 1.0e-10_wp, &
                       "T1: off-grid linear-in-time recovery wrong")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
   end subroutine test_linear_in_time

   ! =================================================================
   ! T1's 3-D twin.
   ! =================================================================

   subroutine test_linear_in_time_3d(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 2, NZ = 3, NT = 4
      real(wp), parameter :: A = 1.0_wp, B = 2.0_wp, C = 3.0_wp, E = 4.0_wp, D = 0.25_wp
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NZ, NT), dest(NX, NY, NZ)
      integer :: i, j, k, m, id
      real(wp) :: tq, expect

      do m = 1, NT
         t_axis(m) = real(m - 1, wp)*150.0_wp
      end do
      do m = 1, NT
         do k = 1, NZ
            do j = 1, NY
               do i = 1, NX
                  f(i, j, k, m) = A + B*real(i, wp) + C*real(j, wp) + &
                                  E*real(k, wp) + D*t_axis(m)
               end do
            end do
         end do
      end do

      fname = "/tmp/test_data_input_linear_3d.nc"
      call write_3d_time_file(trim(fname), NX, NY, NZ, NT, t_axis, f)

      grid = make_grid(NX, NY)
      call reader%init()
      call ocean_data_input_register_3d(reader, trim(fname), "field", grid, NZ, &
                                        NX, NY, NZ, 1, 1, id)

      tq = 217.0_wp
      dest = -999.0_wp
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_3d(reader, id, tq, NX, NY, NZ, dest)
      !$acc update self(dest)
      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               expect = A + B*real(i, wp) + C*real(j, wp) + E*real(k, wp) + D*tq
               call check(error, abs(dest(i, j, k) - expect) < 1.0e-10_wp, &
                          "T1-3D: off-grid linear-in-time recovery wrong")
               if (allocated(error)) return
            end do
         end do
      end do

      call reader%destroy()
   end subroutine test_linear_in_time_3d

   ! =================================================================
   ! T2 — on-grid query returns w=0 at the node, not w=1 from the left.
   ! =================================================================

   subroutine test_on_grid_exact(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t_axis(4)
      integer :: n0, n1
      real(wp) :: w
      logical :: oor

      ! Half-open bracket convention: t == t_axis(3) exactly must land
      ! on the (3,4) bracket at w=0 (record 3 dominates exactly), NOT
      ! the (2,3) bracket at w=1 — the classic off-by-one at the node.
      ! Either labeling recovers record 3's VALUE exactly (w picks the
      ! side), but the bracket indices must be the (3,4)/w=0 pairing —
      ! that is the actual invariant `data_input_locate` promises.
      t_axis = [0.0_wp, 100.0_wp, 200.0_wp, 300.0_wp]
      call data_input_locate(t_axis, 4, DATA_TIME_LINEAR, 0.0_wp, 200.0_wp, n0, n1, w, oor)

      call check(error, n0 == 3, "T2: on-grid query bracket n0 wrong (off-by-one)")
      if (allocated(error)) return
      call check(error, n1 == 4, "T2: on-grid query bracket n1 wrong (off-by-one)")
      if (allocated(error)) return
      call check(error, abs(w) < 1.0e-14_wp, "T2: on-grid query weight must be exactly 0")
      if (allocated(error)) return
      call check(error,.not. oor, "T2: an in-range node must not be flagged out-of-range")
   end subroutine test_on_grid_exact

   ! =================================================================
   ! T3 — cyclic climatology: the (nt, 1) seam + second-cycle identity.
   ! =================================================================

   subroutine test_cyclic_wraparound(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 2, NT = 4
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp), parameter :: PERIOD = 360.0_wp
      real(wp), parameter :: BX = 0.01_wp, BY = 0.02_wp
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest1(NX, NY), dest2(NX, NY)
      integer :: i, j, k, id
      real(wp) :: expect

      t_axis = [0.0_wp, 90.0_wp, 180.0_wp, 270.0_wp]
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = sin(2.0_wp*PI*t_axis(k)/PERIOD) + BX*real(i, wp) + BY*real(j, wp)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_cyclic.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)
      grid = make_grid(NX, NY)

      ! First cycle: query the seam bracket at t = 315.
      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, time_mode="cyclic", &
                                        cycle_period=PERIOD)
      dest1 = -999.0_wp
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest1)
      call ocean_data_input_update_all(reader, 315.0_wp)
      call ocean_data_input_update_2d(reader, id, 315.0_wp, NX, NY, dest1)
      !$acc update self(dest1)
      !$acc exit data delete(dest1)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()

      ! Analytic seam value: 0.5*(f(270) + f(0)).
      do j = 1, NY
         do i = 1, NX
            expect = 0.5_wp*(sin(2.0_wp*PI*270.0_wp/PERIOD) + sin(0.0_wp)) + &
                     BX*real(i, wp) + BY*real(j, wp)
            call check(error, abs(dest1(i, j) - expect) < 1.0e-10_wp, &
                       "T3: seam-bracket (nt,1) value wrong")
            if (allocated(error)) return
         end do
      end do

      ! Second cycle: t = 315 + 360 = 675 must reproduce the SAME value.
      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, time_mode="cyclic", &
                                        cycle_period=PERIOD)
      dest2 = -999.0_wp
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest2)
      call ocean_data_input_update_all(reader, 675.0_wp)
      call ocean_data_input_update_2d(reader, id, 675.0_wp, NX, NY, dest2)
      !$acc update self(dest2)
      !$acc exit data delete(dest2)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()

      do j = 1, NY
         do i = 1, NX
            call check(error, abs(dest2(i, j) - dest1(i, j)) < 1.0e-12_wp, &
                       "T3: second-cycle value not bit-identical to the first")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_cyclic_wraparound

   ! =================================================================
   ! T4 — DATA_TIME_STATIC reads once, ever.
   ! =================================================================

   subroutine test_static_mode(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 2, NT = 3
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest(NX, NY)
      integer :: i, j, k, id
      real(wp) :: expect

      t_axis = [0.0_wp, 10.0_wp, 20.0_wp]
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               ! Records 2/3 are deliberately DIFFERENT from record 1 —
               ! if the reader ever re-reads, this test will fail.
               f(i, j, k) = real(i, wp) + 10.0_wp*real(j, wp) + 1000.0_wp*real(k, wp)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_static.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)
      grid = make_grid(NX, NY)

      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, time_mode="static")

      call check(error, reader%fields(id)%nreads == 1, &
                 "T4: registration of a STATIC field must read exactly once")
      if (allocated(error)) return

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)

      dest = -999.0_wp
      call ocean_data_input_update_all(reader, 0.0_wp)
      call ocean_data_input_update_2d(reader, id, 0.0_wp, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect = real(i, wp) + 10.0_wp*real(j, wp) + 1000.0_wp
            call check(error, abs(dest(i, j) - expect) < 1.0e-12_wp, &
                       "T4: static value at t=0 must be record 1")
            if (allocated(error)) return
         end do
      end do

      dest = -999.0_wp
      call ocean_data_input_update_all(reader, 1.0e9_wp)
      call ocean_data_input_update_2d(reader, id, 1.0e9_wp, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect = real(i, wp) + 10.0_wp*real(j, wp) + 1000.0_wp
            call check(error, abs(dest(i, j) - expect) < 1.0e-12_wp, &
                       "T4: static value at t=1e9 must STILL be record 1")
            if (allocated(error)) return
         end do
      end do

      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)

      call check(error, reader%fields(id)%nreads == 1, &
                 "T4: two per-step queries must not trigger a second slab read")
      if (allocated(error)) return

      call reader%destroy()
   end subroutine test_static_mode

   ! =================================================================
   ! T5 — scale/add_offset applied once at read.
   ! =================================================================

   subroutine test_scale_offset(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 1, NT = 4
      real(wp), parameter :: SCALE = 0.01_wp, ADD_OFFSET = 273.15_wp
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f_raw(NX, NY, NT), dest(NX, NY)
      integer :: i, j, k, id
      real(wp) :: tq, expect_raw

      do k = 1, NT
         t_axis(k) = real(k - 1, wp)*10.0_wp
      end do
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f_raw(i, j, k) = 100.0_wp*real(i, wp) + 5.0_wp*t_axis(k)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_scale.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f_raw)
      grid = make_grid(NX, NY)

      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, &
                                        scale=SCALE, add_offset=ADD_OFFSET)

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)

      ! Two steps in the SAME bracket [t_axis(1), t_axis(2)] = [0, 10].
      tq = 3.0_wp
      dest = -999.0_wp
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_2d(reader, id, tq, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect_raw = 100.0_wp*real(i, wp) + 5.0_wp*tq
            call check(error, abs(dest(i, j) - (SCALE*expect_raw + ADD_OFFSET)) < 1.0e-9_wp, &
                       "T5: scale/add_offset wrong at step 1")
            if (allocated(error)) return
         end do
      end do

      tq = 7.0_wp
      dest = -999.0_wp
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_2d(reader, id, tq, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect_raw = 100.0_wp*real(i, wp) + 5.0_wp*tq
            call check(error, abs(dest(i, j) - (SCALE*expect_raw + ADD_OFFSET)) < 1.0e-9_wp, &
                       "T5: scale/add_offset wrong at step 2 (same bracket)")
            if (allocated(error)) return
         end do
      end do

      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()
   end subroutine test_scale_offset

   ! =================================================================
   ! T6 — out-of-range: CLAMP does not abort; ERROR's condition checked
   ! via the pure `out_of_range` flag (no in-process death test).
   ! =================================================================

   subroutine test_oor_behaviour(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 1, NT = 3
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest(NX, NY)
      integer :: i, j, k, id
      real(wp) :: expect
      real(wp) :: w
      integer :: n0, n1
      logical :: oor

      ! -- Pure-helper half: an out-of-domain query IS flagged. --
      t_axis = [0.0_wp, 10.0_wp, 20.0_wp]
      call data_input_locate(t_axis, 3, DATA_TIME_LINEAR, 0.0_wp, 25.0_wp, n0, n1, w, oor)
      call check(error, oor, "T6: a query past t_axis(nt) must be flagged out_of_range")
      if (allocated(error)) return
      call check(error, n0 == 3 .and. n1 == 3, &
                 "T6: an out-of-range-past-the-end query clamps to the last record")
      if (allocated(error)) return

      ! -- Live path: DATA_OOR_CLAMP does NOT abort, holds the end record. --
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = real(i, wp) + 1000.0_wp*real(k, wp)
            end do
         end do
      end do
      fname = "/tmp/test_data_input_oor.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)
      grid = make_grid(NX, NY)

      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, oor=DATA_OOR_CLAMP)

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)
      dest = -999.0_wp
      call ocean_data_input_update_all(reader, 999.0_wp)
      call ocean_data_input_update_2d(reader, id, 999.0_wp, NX, NY, dest)
      !$acc update self(dest)
      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = 1, NY
         do i = 1, NX
            expect = real(i, wp) + 3000.0_wp
            call check(error, abs(dest(i, j) - expect) < 1.0e-10_wp, &
                       "T6: DATA_OOR_CLAMP must hold the last record's value")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
   end subroutine test_oor_behaviour

   ! =================================================================
   ! T7 — unimplemented time-mode tag: the fail-loud predicate.
   ! =================================================================

   subroutine test_bad_mode(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error,.not. data_input_time_mode_is_implemented("bilinear"), &
                 "T7: 'bilinear' must not be flagged implemented")
      if (allocated(error)) return
      call check(error, data_input_time_mode_is_implemented("linear"), &
                 "T7: 'linear' must be implemented")
      if (allocated(error)) return
      call check(error, data_input_time_mode_is_implemented("cyclic"), &
                 "T7: 'cyclic' must be implemented")
      if (allocated(error)) return
      call check(error, data_input_time_mode_is_implemented("static"), &
                 "T7: 'static' must be implemented")
   end subroutine test_bad_mode

   ! =================================================================
   ! T8 — dimension-mismatch predicate (mirrors zinit_dims_ok).
   ! =================================================================

   subroutine test_dim_mismatch(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 4, NY = 3, NT = 2
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT)

      t_axis = [0.0_wp, 10.0_wp]
      f = 1.0_wp

      fname = "/tmp/test_data_input_mismatch.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)

      call check(error, data_input_dims_ok(trim(fname), "field", NX, NY, .false.), &
                 "T8: correct dims flagged as mismatch")
      if (allocated(error)) return
      call check(error,.not. data_input_dims_ok(trim(fname), "field", NX + 1, NY, .false.), &
                 "T8: x mismatch not detected")
      if (allocated(error)) return
      call check(error,.not. data_input_dims_ok(trim(fname), "field", NX, NY + 2, .false.), &
                 "T8: y mismatch not detected")
      if (allocated(error)) return
      call check(error,.not. data_input_dims_ok(trim(fname), "nope", NX, NY, .false.), &
                 "T8: missing variable not detected")
      if (allocated(error)) return
      call check(error,.not. data_input_dims_ok(trim(fname), "field", NX, NY, .true.), &
                 "T8: rank mismatch (2-D file registered as 3-D) not detected")
   end subroutine test_dim_mismatch

   ! =================================================================
   ! T9 — data_input_locate, exhaustive sweep.
   ! =================================================================

   subroutine test_locate_pure(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t_axis(4), t_one(1)
      integer :: n0, n1
      real(wp) :: w
      logical :: oor

      t_axis = [0.0_wp, 100.0_wp, 200.0_wp, 300.0_wp]

      ! Interior, LINEAR.
      call data_input_locate(t_axis, 4, DATA_TIME_LINEAR, 0.0_wp, 150.0_wp, n0, n1, w, oor)
      call check(error, n0 == 2 .and. n1 == 3 .and. abs(w - 0.5_wp) < 1.0e-14_wp .and. .not. oor, &
                 "T9: LINEAR interior bracket/weight wrong")
      if (allocated(error)) return

      ! Below-range, LINEAR (clamped low + flagged).
      call data_input_locate(t_axis, 4, DATA_TIME_LINEAR, 0.0_wp, -50.0_wp, n0, n1, w, oor)
      call check(error, n0 == 1 .and. n1 == 1 .and. oor, &
                 "T9: LINEAR below-range must clamp to record 1 and flag out_of_range")
      if (allocated(error)) return

      ! Above-range, LINEAR (clamped high + flagged).
      call data_input_locate(t_axis, 4, DATA_TIME_LINEAR, 0.0_wp, 999.0_wp, n0, n1, w, oor)
      call check(error, n0 == 4 .and. n1 == 4 .and. oor, &
                 "T9: LINEAR above-range must clamp to record nt and flag out_of_range")
      if (allocated(error)) return

      ! CYCLIC interior (first period).
      call data_input_locate(t_axis, 4, DATA_TIME_CYCLIC, 400.0_wp, 250.0_wp, n0, n1, w, oor)
      call check(error, n0 == 3 .and. n1 == 4 .and. abs(w - 0.5_wp) < 1.0e-14_wp .and. .not. oor, &
                 "T9: CYCLIC interior bracket/weight wrong")
      if (allocated(error)) return

      ! CYCLIC seam bracket (nt, 1).
      call data_input_locate(t_axis, 4, DATA_TIME_CYCLIC, 400.0_wp, 350.0_wp, n0, n1, w, oor)
      call check(error, n0 == 4 .and. n1 == 1 .and. .not. oor, &
                 "T9: CYCLIC seam bracket wrong")
      if (allocated(error)) return

      ! CYCLIC never out-of-range, even far beyond the axis.
      call data_input_locate(t_axis, 4, DATA_TIME_CYCLIC, 400.0_wp, 1.0e6_wp, n0, n1, w, oor)
      call check(error,.not. oor, "T9: CYCLIC must never flag out_of_range")
      if (allocated(error)) return

      ! Degenerate nt == 1: never divides, never out of range.
      t_one = [42.0_wp]
      call data_input_locate(t_one, 1, DATA_TIME_LINEAR, 0.0_wp, -1.0e6_wp, n0, n1, w, oor)
      call check(error, n0 == 1 .and. n1 == 1 .and. abs(w) < 1.0e-14_wp .and. .not. oor, &
                 "T9: nt==1 must be a degenerate always-record-1 case")
   end subroutine test_locate_pure

   ! =================================================================
   ! T10 — GPU residency: a bracket ADVANCE must push new data to the
   ! device.  Both queries share ONE mapped session (no re-enter_data in
   ! between) — the only way this can fail silently on multicore and
   ! pass is if `!$acc update device` didn't fire on the advance.
   ! =================================================================

   subroutine test_gpu_resident(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 2, NT = 6
      real(wp), parameter :: A = 1.0_wp, B = 2.0_wp, C = 3.0_wp, D = 0.5_wp
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest(NX, NY)
      integer :: i, j, k, id
      real(wp) :: tq, expect

      do k = 1, NT
         t_axis(k) = real(k - 1, wp)*100.0_wp
      end do
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = A + B*real(i, wp) + C*real(j, wp) + D*t_axis(k)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_gpu_resident.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)
      grid = make_grid(NX, NY)

      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id)

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest)

      ! Query 1: bracket (1,2).
      tq = 137.0_wp
      dest = -999.0_wp
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_2d(reader, id, tq, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect = A + B*real(i, wp) + C*real(j, wp) + D*tq
            call check(error, abs(dest(i, j) - expect) < 1.0e-10_wp, &
                       "T10: first (pre-advance) query wrong")
            if (allocated(error)) return
         end do
      end do

      ! Query 2: advance PAST the current bracket to (4,5) — this is the
      ! query that can only be right if the device copy of f0/f1 was
      ! actually refreshed.  NO enter_data call between query 1 and 2.
      tq = 337.0_wp
      dest = -999.0_wp
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_input_update_2d(reader, id, tq, NX, NY, dest)
      !$acc update self(dest)
      do j = 1, NY
         do i = 1, NX
            expect = A + B*real(i, wp) + C*real(j, wp) + D*tq
            call check(error, abs(dest(i, j) - expect) < 1.0e-10_wp, &
                       "T10: post-advance query wrong — device slab not refreshed?")
            if (allocated(error)) return
         end do
      end do

      !$acc exit data delete(dest)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()
   end subroutine test_gpu_resident

   ! =================================================================
   ! fill_static_host: host-only route agrees bit-for-bit with the
   ! device blend route.
   ! =================================================================

   subroutine test_fill_static_host(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 3, NT = 2
      type(ocean_data_input_t) :: reader
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), f(NX, NY, NT), dest_host(NX, NY), dest_dev(NX, NY)
      integer :: i, j, k, id

      t_axis = [0.0_wp, 10.0_wp]
      do k = 1, NT
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = real(i, wp) + 100.0_wp*real(j, wp) + 1.0e6_wp*real(k, wp)
            end do
         end do
      end do

      fname = "/tmp/test_data_input_fill_static.nc"
      call write_2d_time_file(trim(fname), NX, NY, NT, t_axis, f)
      grid = make_grid(NX, NY)

      call reader%init()
      call ocean_data_input_register_2d(reader, trim(fname), "field", grid, &
                                        NX, NY, 1, 1, id, time_mode="static")

      ! Host-only route — BEFORE any device mapping exists.
      dest_host = -999.0_wp
      call ocean_data_input_fill_static_host(reader, id, NX, NY, dest_host)

      ! Device-blend route, for comparison.
      dest_dev = -999.0_wp
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(dest_dev)
      call ocean_data_input_update_all(reader, 0.0_wp)
      call ocean_data_input_update_2d(reader, id, 0.0_wp, NX, NY, dest_dev)
      !$acc update self(dest_dev)
      !$acc exit data delete(dest_dev)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = 1, NY
         do i = 1, NX
            call check(error, abs(dest_host(i, j) - dest_dev(i, j)) < 1.0e-12_wp, &
                       "fill_static_host disagrees with the device blend route")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
   end subroutine test_fill_static_host

end module test_ocean_data_input
