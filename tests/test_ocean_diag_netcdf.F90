!! Phase 6d serial-NetCDF emit tests.
!!
!! Opens a per-process NetCDF stream, fires the diag manager a few
!! times, closes, reopens via the bare NetCDF wrapper, and verifies
!! that the time dim grew and the data slices round-trip.  Skipped
!! when the build does not include NetCDF.
module test_ocean_diag_netcdf
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: ocean_diag_t, DIAG_OP_INSTANT
   use rdb_ocean_diag_fills, only: fill_ssh, fill_temperature, register_default_diags
   use rdb_ocean_diag_derived, only: apply_diag_selection
   use rdb_ocean_diag_netcdf, only: open_stream, close_stream
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, nc_get_dim_len, &
                            nc_get_varid
   use netcdf, only: nf90_get_var
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_diag_netcdf_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1.0_wp

contains

   subroutine collect_ocean_diag_netcdf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("netcdf_writes_2d_var_per_fire", test_nc_write_2d), &
                  new_unittest("netcdf_writes_3d_var_per_fire", test_nc_write_3d), &
                  new_unittest("netcdf_close_unbinds_emit_hook", test_nc_close_unbinds), &
                  new_unittest("diag_registered_set_unchanged", test_registered_set_unchanged) &
                  ]
   end subroutine collect_ocean_diag_netcdf_tests

   subroutine dummy_fill(state_handle, buf)
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      buf = 42.0_wp
      if (.false.) then
         select type (state_handle)
         class default
         end select
      end if
   end subroutine dummy_fill

   subroutine setup_state(grid, state)
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
   end subroutine setup_state

   subroutine cleanup_file(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=filename, exist=exists)
      if (exists) then
         open (newunit=u, file=filename, status='old')
         close (u, status='delete')
      end if
   end subroutine cleanup_file

   subroutine test_nc_write_2d(error)
      !! Open NetCDF stream, register a 2D var (SSH), fire twice via
      !! the diag manager, close, reopen via NetCDF and verify the time
      !! dim length and the stored slice values.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=*), parameter :: FN = "test_diag_nc_2d.nc"
      integer :: step, ncid, nfires, varid, i, j
      real(wp), allocatable :: read_buf(:, :, :)

      call cleanup_file(FN)
      call setup_state(grid, state)

      do j = 1, NY
         do i = 1, NX
            state%barotropic%h(i, j) = real(i + j, wp)
            state%barotropic%b(i, j) = 0.0_wp
            ! `fill_ssh` reads from `multilayer%h_layer` on the
            ! multilayer path; stamp the layered column-sum to match
            ! `barotropic%h` so the test asserts the same value.
            state%multilayer%h_layer(i, j, :) = real(i + j, wp)/real(NZ, wp)
         end do
      end do

      call state%diag%register("SSH", units="m", fill=fill_ssh, &
                               n1=NX, n2=NY, n3=1, &
                               long_name="sea_surface_height", &
                               time_op=DIAG_OP_INSTANT, dt_out=4.0_wp)
      call open_stream(state%diag, FN)
      checks: block
         call check(error, state%diag%nc_stream%is_open, "stream should be open")
         if (allocated(error)) exit checks
         call check(error, associated(state%diag%emit_post_fire), &
                    "emit_post_fire should be bound after open_stream")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         ! 8 steps × dt=1 → 2 fires (at steps 4 and 8).
         do step = 1, 8
            call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
         end do
         call ocean_state_exit_data(state)
         call close_stream(state%diag)

         ! Reopen and verify.  Nested blocks track the acquire/release
         ! pairs for `ncid` (open at nc_open_read, close at end) and
         ! `read_buf` (allocated then deallocated).
         call nc_open_read(FN, ncid)
         ncid_scope: block
            call nc_get_dim_len(ncid, "time_SSH", nfires)
            call check(error, nfires == 2, "time_SSH dim should have 2 entries")
            if (allocated(error)) exit ncid_scope

            call nc_get_varid(ncid, "SSH", varid)
            allocate (read_buf(NX, NY, 2))
            buf_scope: block
               call nc_check(nf90_get_var(ncid, varid, read_buf), "reading SSH back")
               call check(error, abs(read_buf(2, 2, 1) - state%barotropic%h(2, 2)) < 1.0e-10_wp, &
                          "SSH slice 1 should match h(2,2)")
               if (allocated(error)) exit buf_scope
               call check(error, abs(read_buf(2, 2, 2) - state%barotropic%h(2, 2)) < 1.0e-10_wp, &
                          "SSH slice 2 should match h(2,2)")
            end block buf_scope
            deallocate (read_buf)
         end block ncid_scope
         call nc_close(ncid)
      end block checks
      call state%destroy()
      call cleanup_file(FN)
   end subroutine test_nc_write_2d

   subroutine test_nc_write_3d(error)
      !! 3D var (temperature, NZ layers) exercises the 4D NetCDF path
      !! `(x, y, z_temperature, time)`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=*), parameter :: FN = "test_diag_nc_3d.nc"
      integer :: it, k, step, ncid, nz_dim, varid
      real(wp), allocatable :: read_buf(:, :, :, :)

      call cleanup_file(FN)
      call setup_state(grid, state)

      state%multilayer%h_layer = 5.0_wp
      it = state%multilayer%idx_temperature
      do k = 1, NZ
         state%multilayer%tracers(it)%hTr(:, :, k) = 5.0_wp*real(k*2, wp)  ! T = 2k
      end do

      call state%diag%register("temperature", units="degC", fill=fill_temperature, &
                               n1=NX, n2=NY, n3=NZ, &
                               long_name="sea_water_potential_temperature", &
                               time_op=DIAG_OP_INSTANT, dt_out=2.0_wp)
      call open_stream(state%diag, FN)
      call ocean_state_enter_data(state)
      do step = 1, 4
         call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
      end do
      call ocean_state_exit_data(state)
      call close_stream(state%diag)

      call nc_open_read(FN, ncid)
      ncid_scope: block
         call nc_get_dim_len(ncid, "z_temperature", nz_dim)
         call check(error, nz_dim == NZ, "z_temperature dim should equal NZ")
         if (allocated(error)) exit ncid_scope

         call nc_get_varid(ncid, "temperature", varid)
         allocate (read_buf(NX, NY, NZ, 2))
         buf_scope: block
            call nc_check(nf90_get_var(ncid, varid, read_buf), "reading temperature back")
            call check(error, abs(read_buf(2, 2, NZ, 1) - 6.0_wp) < 1.0e-10_wp, &
                       "T at (2,2,k=NZ,t=1) should be 6.0 (k=3 -> 2*3)")
            if (allocated(error)) exit buf_scope
            call check(error, abs(read_buf(2, 2, 1, 1) - 2.0_wp) < 1.0e-10_wp, &
                       "T at (2,2,k=1,t=1) should be 2.0 (k=1 -> 2*1)")
         end block buf_scope
         deallocate (read_buf)
      end block ncid_scope
      call nc_close(ncid)
      call state%destroy()
      call cleanup_file(FN)
   end subroutine test_nc_write_3d

   subroutine test_nc_close_unbinds(error)
      !! After close_stream, emit_post_fire must be null and is_open
      !! false.  Guards against second-open / use-after-close paths.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=*), parameter :: FN = "test_diag_nc_close.nc"
      checks: block

         call cleanup_file(FN)
         call setup_state(grid, state)
         call state%diag%register("foo", units="m", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=1)
         call open_stream(state%diag, FN)
         call close_stream(state%diag)

         call check(error,.not. state%diag%nc_stream%is_open, "stream should be closed")
         if (allocated(error)) exit checks
         call check(error,.not. associated(state%diag%emit_post_fire), &
                    "emit_post_fire should be nullified after close")

      end block checks
      call state%destroy()
      call cleanup_file(FN)
   end subroutine test_nc_close_unbinds

   subroutine test_registered_set_unchanged(error)
      !! PR-64 §9 case 6 — the bit-identity argument, made observable. For
      !! a spec with no gated-off entries (a blank `diags`), the registered
      !! name-set + ORDER + `nvars` must be exactly the pre-PR canonical
      !! six (SSH, temperature, salinity, u, v, KE) on a default state (ice
      !! / EPBL / kappa-shear / ideal-age all off). Order matters: it is
      !! the NetCDF variable write order. This is the unit-level proxy for
      !! "every shipped validation_examples/ocean/*.nml registers the same
      !! set before and after this PR" — the warn only fires on a path
      !! that today registers nothing at all, so it cannot perturb this.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=16), parameter :: expect(6) = [character(len=16) :: &
                                                   "SSH", "temperature", "salinity", &
                                                   "u", "v", "KE"]
      integer :: i
      checks: block
         call setup_state(grid, state)
         call apply_diag_selection(state, "", dt_out=3600.0_wp)

         call check(error, state%diag%nvars == 6, "registered set should have exactly 6 vars")
         if (allocated(error)) exit checks
         do i = 1, 6
            call check(error, trim(state%diag%vars(i)%name) == trim(expect(i)), &
                       "var "//trim(expect(i))//" should be at registry slot "// &
                       achar(48 + i))
            if (allocated(error)) exit checks
         end do
      end block checks
      call state%destroy()
   end subroutine test_registered_set_unchanged

end module test_ocean_diag_netcdf
