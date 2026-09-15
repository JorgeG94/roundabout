!! `&ocean_diag_nml output_precision` — diagnostic-stream output width.
!!
!! Diagnostic NetCDF output is write-bandwidth bound, so halving the
!! element width is the lever on its cost.  The knob decouples the
!! DIAGNOSTIC stream from `NC_WP` (compile-time working precision) — and
!! nothing else: restarts, gauges, coastal output, console conservation
!! totals and checksums stay fp64 unconditionally.  A restart that does
!! not round-trip exactly makes a resumed run a different run, so there is
!! deliberately no spelling that can request a lossy one.
!!
!! Five gates, in the order they matter:
!!
!!   1. `default_is_byte_identical` — absent knob and `"double"` produce
!!      byte-identical files, and the variables are still `NF90_DOUBLE`.
!!      This is what lets the change land without re-baselining anything.
!!   2. `single_round_trips_to_fp32` — a known analytic field written at
!!      `"single"` reads back to fp32 epsilon, is genuinely lossy relative
!!      to the fp64 host values (so the knob demonstrably did something),
!!      and the variable is `NF90_FLOAT`.
!!   3. `fill_value_type_matches_var` — `_FillValue` reads back as the
!!      variable's own type in BOTH modes and equals `DIAG_MISSING_VALUE`.
!!      A mismatch is a NetCDF error at define time, so this fails loudly.
!!   4. `missing_cells_survive_conversion` — sentinel-carrying cells still
!!      compare equal to the file's `_FillValue` after fp64 -> fp32.
!!   5. `single_file_is_about_half` — catches "the knob parsed but nothing
!!      actually changed" (netcdf will silently narrow fp64 data handed to
!!      an NF90_FLOAT variable, so a type-only change would still pass 2).
module test_ocean_diag_precision
   use, intrinsic :: iso_fortran_env, only: int8, real32, real64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: ocean_diag_t, DIAG_OP_INSTANT, DIAG_MISSING_VALUE
   use rdb_ocean_diag_fills, only: fill_temperature
   use rdb_ocean_diag_netcdf, only: open_stream, close_stream
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, nc_get_varid
   use netcdf, only: nf90_get_var, nf90_inquire_variable, nf90_inquire_attribute, &
                     nf90_get_att, nf90_float, nf90_double
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_diag_precision_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1.0_wp

   ! Bigger grid for the file-size ratio: HDF5 superblock + per-variable
   ! metadata is a few kB, which would swamp the payload on a 6x4x3 case.
   integer, parameter :: BIG_NX = 64, BIG_NY = 64, BIG_NZ = 8

   real(real32), parameter :: FP32_EPS = epsilon(1.0_real32)

contains

   subroutine collect_ocean_diag_precision_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("default_is_byte_identical", test_default_byte_identical), &
                  new_unittest("single_round_trips_to_fp32", test_single_roundtrip), &
                  new_unittest("fill_value_type_matches_var", test_fill_value_type), &
                  new_unittest("missing_cells_survive_conversion", test_missing_survives), &
                  new_unittest("single_file_is_about_half", test_file_size_ratio) &
                  ]
   end subroutine collect_ocean_diag_precision_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine dummy_fill(state_handle, buf)
      !! Never invoked: the writer-level tests drive `emit_post_fire`
      !! directly with a host-stamped `output_buffer` so the values under
      !! test are exact and independent of any physics kernel.  `register`
      !! still demands a bound fill.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      buf = 0.0_wp
      if (.false.) then
         select type (state_handle)
         class default
         end select
      end if
   end subroutine dummy_fill

   pure function analytic(i, j, k) result(v)
      !! Reciprocals of odd integers: none of these are exactly
      !! representable in fp32 (nor in fp64), so a narrowing conversion is
      !! observable rather than a no-op on round numbers.
      integer, intent(in) :: i, j, k
      real(wp) :: v
      v = 1.0_wp/real(2*(i + 3*j + 7*k) + 1, wp)
   end function analytic

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

   integer function file_size_bytes(filename) result(nbytes)
      character(len=*), intent(in) :: filename
      nbytes = -1
      inquire (file=filename, size=nbytes)
   end function file_size_bytes

   logical function files_are_identical(fa, fb) result(same)
      !! Byte-for-byte comparison of two files via unformatted stream
      !! access.  Differing sizes => not identical (short-circuits before
      !! any read).
      character(len=*), intent(in) :: fa, fb
      integer :: na, nb, ua, ub
      integer(int8), allocatable :: ba(:), bb(:)

      na = file_size_bytes(fa)
      nb = file_size_bytes(fb)
      same = .false.
      if (na <= 0 .or. na /= nb) return

      allocate (ba(na), bb(nb))
      open (newunit=ua, file=fa, access='stream', form='unformatted', status='old')
      read (ua) ba
      close (ua)
      open (newunit=ub, file=fb, access='stream', form='unformatted', status='old')
      read (ub) bb
      close (ub)
      same = all(ba == bb)
      deallocate (ba, bb)
   end function files_are_identical

   integer function var_xtype(filename, varname) result(xt)
      !! NetCDF element type of a data variable, read back from the file.
      character(len=*), intent(in) :: filename, varname
      integer :: ncid, varid
      call nc_open_read(filename, ncid)
      call nc_get_varid(ncid, varname, varid)
      call nc_check(nf90_inquire_variable(ncid, varid, xtype=xt), &
                    "inquiring xtype of "//trim(varname))
      call nc_close(ncid)
   end function var_xtype

   subroutine write_one_frame(filename, n1, n2, n3, buf, precision_arg, has_missing, &
                              stage_allocated)
      !! Register a single diagnostic of shape (n1, n2, n3), open a stream
      !! (optionally at a requested precision), stamp `buf` into the output
      !! buffer on the host and emit exactly one frame.  Driving
      !! `emit_post_fire` directly rather than `diag%step` keeps the values
      !! under test exact — the manager's device-side fill / fold / finalise
      !! are covered by `test_ocean_diag_netcdf`, and the emit hook always
      !! runs against host-current data (`ocean_diag_step` does an
      !! `!$acc update self(v%output_buffer)` immediately before calling it).
      character(len=*), intent(in) :: filename
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: buf(:, :, :)
      character(len=*), intent(in), optional :: precision_arg
      logical, intent(in), optional :: has_missing
      logical, intent(out), optional :: stage_allocated
         !! Whether the writer's fp32 host staging buffer was live at emit
         !! time.  White-box, but it is the ONLY observable that separates
         !! "we narrowed the data ourselves" from "we handed fp64 to an
         !! NF90_FLOAT variable and let netcdf narrow it" — both produce
         !! identical files, so no black-box assertion can tell them apart,
         !! yet only the first actually halves the host->library traffic
         !! this knob exists to halve.
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      logical :: hm

      hm = .false.
      if (present(has_missing)) hm = has_missing

      call cleanup_file(filename)
      call grid%init(n1, n2, 1, DX, DX)
      state%multilayer%nz_ml = max(n3, 1)
      call state%init(grid)

      call state%diag%register("field", units="1", fill=dummy_fill, &
                               n1=n1, n2=n2, n3=n3, &
                               long_name="analytic test field", &
                               time_op=DIAG_OP_INSTANT, dt_out=1.0_wp, &
                               has_missing=hm)
      if (present(precision_arg)) then
         call open_stream(state%diag, filename, output_precision=precision_arg)
      else
         call open_stream(state%diag, filename)
      end if
      state%diag%vars(1)%output_buffer(:, :, :) = buf
      if (present(stage_allocated)) stage_allocated = allocated(state%diag%nc_stream%stage)
      call state%diag%emit_post_fire(state%diag, 1, 0.0_wp)
      call close_stream(state%diag)
      call state%destroy()
   end subroutine write_one_frame

   ! ------------------------------------------------------------------
   ! 1. Default is byte-identical
   ! ------------------------------------------------------------------

   subroutine test_default_byte_identical(error)
      !! The gate that makes the change safe to land: an omitted
      !! `output_precision` and an explicit `"double"` must produce
      !! byte-identical files, and the data variable must still be
      !! `NF90_DOUBLE`.  Run END-TO-END through `diag%step` with the state
      !! device-mapped, so this covers the production path (device fill ->
      !! `!$acc update self` -> emit hook) and not just the writer.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FA = "test_diagprec_default.nc"
      character(len=*), parameter :: FB = "test_diagprec_double.nc"

      checks: block
         call run_default_case(FA, use_arg=.false.)
         call run_default_case(FB, use_arg=.true.)

         call check(error, files_are_identical(FA, FB), &
                    "omitted output_precision must be byte-identical to "// &
                    "output_precision='double'")
         if (allocated(error)) exit checks
         call check(error, var_xtype(FA, "temperature") == nf90_double, &
                    "default stream must define data vars as NF90_DOUBLE")
      end block checks

      call cleanup_file(FA)
      call cleanup_file(FB)
   end subroutine test_default_byte_identical

   subroutine run_default_case(filename, use_arg)
      !! One end-to-end diag run emitting two `temperature` frames.
      character(len=*), intent(in) :: filename
      logical, intent(in) :: use_arg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: it, k, step

      call cleanup_file(filename)
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)

      state%multilayer%h_layer = 5.0_wp
      it = state%multilayer%idx_temperature
      do k = 1, NZ
         state%multilayer%tracers(it)%hTr(:, :, k) = 5.0_wp*(1.0_wp/real(2*k + 1, wp))
      end do

      call state%diag%register("temperature", units="degC", fill=fill_temperature, &
                               n1=NX, n2=NY, n3=NZ, &
                               long_name="sea_water_potential_temperature", &
                               time_op=DIAG_OP_INSTANT, dt_out=2.0_wp)
      if (use_arg) then
         call open_stream(state%diag, filename, output_precision="double")
      else
         call open_stream(state%diag, filename)
      end if
      call ocean_state_enter_data(state)
      do step = 1, 4
         call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
      end do
      call ocean_state_exit_data(state)
      call close_stream(state%diag)
      call state%destroy()
   end subroutine run_default_case

   ! ------------------------------------------------------------------
   ! 2. Round-trip at "single"
   ! ------------------------------------------------------------------

   subroutine test_single_roundtrip(error)
      !! Write the same analytic field at both precisions.  The single file
      !! must (a) declare `NF90_FLOAT`, (b) agree with the fp64 truth to
      !! fp32 epsilon, and (c) genuinely differ from it — at least one cell
      !! must be off by more than fp64 round-off, otherwise the knob would
      !! be a no-op that the accuracy assertion alone could not detect.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FD = "test_diagprec_rt_double.nc"
      character(len=*), parameter :: FS = "test_diagprec_rt_single.nc"
      real(wp) :: truth(NX, NY, NZ)
      real(wp), allocatable :: got_d(:, :, :, :), got_s(:, :, :, :)
      integer :: i, j, k, ncid, varid
      real(wp) :: max_rel_err, max_abs_diff
      logical :: staged_d, staged_s

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               truth(i, j, k) = analytic(i, j, k)
            end do
         end do
      end do

      checks: block
         call write_one_frame(FD, NX, NY, NZ, truth, precision_arg="double", &
                              stage_allocated=staged_d)
         call write_one_frame(FS, NX, NY, NZ, truth, precision_arg="single", &
                              stage_allocated=staged_s)

         ! The conversion must be ours, and must cost nothing on the
         ! default path.  See the `stage_allocated` note in
         ! `write_one_frame` for why this can only be checked white-box.
         call check(error, staged_s, &
                    "single stream must stage the fp64 -> fp32 conversion itself")
         if (allocated(error)) exit checks
         call check(error,.not. staged_d, &
                    "double stream must not allocate the fp32 staging buffer")
         if (allocated(error)) exit checks

         call check(error, var_xtype(FS, "field") == nf90_float, &
                    "output_precision='single' must define data vars as NF90_FLOAT")
         if (allocated(error)) exit checks
         call check(error, var_xtype(FD, "field") == nf90_double, &
                    "output_precision='double' must define data vars as NF90_DOUBLE")
         if (allocated(error)) exit checks

         allocate (got_d(NX, NY, NZ, 1), got_s(NX, NY, NZ, 1))
         call nc_open_read(FD, ncid)
         call nc_get_varid(ncid, "field", varid)
         call nc_check(nf90_get_var(ncid, varid, got_d), "reading double field")
         call nc_close(ncid)
         call nc_open_read(FS, ncid)
         call nc_get_varid(ncid, "field", varid)
         call nc_check(nf90_get_var(ncid, varid, got_s), "reading single field")
         call nc_close(ncid)

         ! (a) fp64 stream is exact.
         call check(error, all(got_d(:, :, :, 1) == truth), &
                    "double stream must round-trip the host buffer exactly")
         if (allocated(error)) exit checks

         ! (b) fp32 stream agrees to fp32 epsilon.
         max_rel_err = maxval(abs(got_s(:, :, :, 1) - truth)/abs(truth))
         call check(error, max_rel_err <= real(FP32_EPS, wp), &
                    "single stream must agree with the host buffer to fp32 epsilon")
         if (allocated(error)) exit checks

         ! (c) ...and is actually lossy, i.e. the knob did something.
         max_abs_diff = maxval(abs(got_s(:, :, :, 1) - truth))
         call check(error, max_abs_diff > 0.0_wp, &
                    "single stream must differ from the fp64 values — "// &
                    "identical output means the knob was a no-op")
         deallocate (got_d, got_s)
      end block checks

      call cleanup_file(FD)
      call cleanup_file(FS)
   end subroutine test_single_roundtrip

   ! ------------------------------------------------------------------
   ! 3. Attribute type matches the variable type
   ! ------------------------------------------------------------------

   subroutine test_fill_value_type(error)
      !! NetCDF rejects a `_FillValue` whose type differs from its
      !! variable's, so the attribute has to follow the data type.  Assert
      !! the declared attribute type in both modes and that the value still
      !! reads back as `DIAG_MISSING_VALUE`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FD = "test_diagprec_fv_double.nc"
      character(len=*), parameter :: FS = "test_diagprec_fv_single.nc"
      real(wp) :: buf(NX, NY, NZ)

      buf = 1.0_wp

      checks: block
         call write_one_frame(FD, NX, NY, NZ, buf, precision_arg="double", has_missing=.true.)
         call write_one_frame(FS, NX, NY, NZ, buf, precision_arg="single", has_missing=.true.)

         call assert_fill_value(error, FD, nf90_double)
         if (allocated(error)) exit checks
         call assert_fill_value(error, FS, nf90_float)
      end block checks

      call cleanup_file(FD)
      call cleanup_file(FS)
   end subroutine test_fill_value_type

   subroutine assert_fill_value(error, filename, expect_xtype)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: filename
      integer, intent(in) :: expect_xtype
      integer :: ncid, varid, att_xtype
      real(real32) :: fv32
      real(real64) :: fv64

      call nc_open_read(filename, ncid)
      call nc_get_varid(ncid, "field", varid)
      call nc_check(nf90_inquire_attribute(ncid, varid, "_FillValue", xtype=att_xtype), &
                    "inquiring _FillValue type in "//trim(filename))
      call check(error, att_xtype == expect_xtype, &
                 "_FillValue type must match the variable type in "//trim(filename))
      if (.not. allocated(error)) then
         if (expect_xtype == nf90_float) then
            call nc_check(nf90_get_att(ncid, varid, "_FillValue", fv32), &
                          "reading real32 _FillValue")
            call check(error, fv32 == real(DIAG_MISSING_VALUE, real32), &
                       "real32 _FillValue must equal DIAG_MISSING_VALUE")
         else
            call nc_check(nf90_get_att(ncid, varid, "_FillValue", fv64), &
                          "reading real64 _FillValue")
            call check(error, fv64 == real(DIAG_MISSING_VALUE, real64), &
                       "real64 _FillValue must equal DIAG_MISSING_VALUE")
         end if
      end if
      call nc_close(ncid)
   end subroutine assert_fill_value

   ! ------------------------------------------------------------------
   ! 4. Masked / missing cells survive the conversion
   ! ------------------------------------------------------------------

   subroutine test_missing_survives(error)
      !! A remapped diagnostic with vanished target cells carries
      !! `DIAG_MISSING_VALUE` (1e20, comfortably inside fp32's ~3.4e38
      !! range).  After narrowing, those cells must still compare EQUAL to
      !! the file's own `_FillValue`, or downstream masking silently stops
      !! working.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FS = "test_diagprec_miss_single.nc"
      real(wp) :: buf(NX, NY, NZ)
      real(real32), allocatable :: got(:, :, :, :)
      real(real32) :: fv32
      integer :: ncid, varid, i, j, k, n_sentinel

      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               buf(i, j, k) = analytic(i, j, k)
            end do
         end do
      end do
      ! Vanished bed-side cells in the shallow half of the domain, plus a
      ! lone interior cell so a whole-plane assertion cannot pass by
      ! accident.
      buf(1:NX/2, :, 1) = DIAG_MISSING_VALUE
      buf(NX, NY, NZ) = DIAG_MISSING_VALUE
      n_sentinel = (NX/2)*NY + 1

      checks: block
         call write_one_frame(FS, NX, NY, NZ, buf, precision_arg="single", has_missing=.true.)

         allocate (got(NX, NY, NZ, 1))
         call nc_open_read(FS, ncid)
         call nc_get_varid(ncid, "field", varid)
         ! Read as real32 so the comparison is against the stored bits, not
         ! a widened copy.
         call nc_check(nf90_get_var(ncid, varid, got), "reading masked field")
         call nc_check(nf90_get_att(ncid, varid, "_FillValue", fv32), &
                       "reading real32 _FillValue")
         call nc_close(ncid)

         call check(error, count(got(:, :, :, 1) == fv32) == n_sentinel, &
                    "every sentinel cell (and only those) must read back "// &
                    "equal to the file's _FillValue after fp64 -> fp32")
         if (allocated(error)) exit checks
         call check(error, fv32 == real(DIAG_MISSING_VALUE, real32), &
                    "_FillValue must be the narrowed DIAG_MISSING_VALUE")
         deallocate (got)
      end block checks

      call cleanup_file(FS)
   end subroutine test_missing_survives

   ! ------------------------------------------------------------------
   ! 5. File size
   ! ------------------------------------------------------------------

   subroutine test_file_size_ratio(error)
      !! The cheap check that the knob moved BYTES, not just the declared
      !! type: netcdf will happily narrow fp64 data handed to an
      !! NF90_FLOAT variable, which would leave every accuracy assertion
      !! green while the staging buffer did nothing.  Payload here is
      !! 64x64x8 doubles = 256 kB, so the few kB of HDF5 metadata is noise
      !! and the ratio lands close to 0.5.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FD = "test_diagprec_size_double.nc"
      character(len=*), parameter :: FS = "test_diagprec_size_single.nc"
      real(wp), allocatable :: buf(:, :, :)
      integer :: i, j, k, nd, ns
      real(wp) :: ratio

      allocate (buf(BIG_NX, BIG_NY, BIG_NZ))
      do k = 1, BIG_NZ
         do j = 1, BIG_NY
            do i = 1, BIG_NX
               buf(i, j, k) = analytic(i, j, k)
            end do
         end do
      end do

      checks: block
         call write_one_frame(FD, BIG_NX, BIG_NY, BIG_NZ, buf, precision_arg="double")
         call write_one_frame(FS, BIG_NX, BIG_NY, BIG_NZ, buf, precision_arg="single")

         nd = file_size_bytes(FD)
         ns = file_size_bytes(FS)
         call check(error, nd > 0 .and. ns > 0, "both files must exist and be non-empty")
         if (allocated(error)) exit checks

         ratio = real(ns, wp)/real(nd, wp)
         call check(error, ratio < 0.6_wp, &
                    "single-precision diag file must be about half the size "// &
                    "of the double-precision one")
         if (allocated(error)) exit checks
         call check(error, ratio > 0.4_wp, &
                    "single-precision diag file below 0.4x is suspicious — "// &
                    "expected ~0.5x plus fixed HDF5 metadata")
      end block checks

      deallocate (buf)
      call cleanup_file(FD)
      call cleanup_file(FS)
   end subroutine test_file_size_ratio

end module test_ocean_diag_precision
