!! Tests for the A2 z-level T/S initial-condition reader
!! (`rdb_ocean_z_init:seed_ts_from_zfile`).
!!
!! Tests (NetCDF I/O — must run with OMP_NUM_THREADS=1):
!!   zinit_exact_recovery   — T1: source z-levels placed EXACTLY at the
!!                            model layer-centre depths ⇒ interpolated
!!                            T/S recovers the source values to roundoff.
!!   zinit_linear           — T2: known linear profile T(z)=a+b*z; the
!!                            mid-layer values equal the numpy-oracle
!!                            hand-computed a+b*z_ctr(k).
!!   zinit_extrapolation    — T3: layer centres shallower than z_src(1)
!!                            and deeper than z_src(nz_src) recover the
!!                            end source values exactly (no linear tails).
!!   zinit_dry_column       — T4: a wet_mask=0 column gets land_fill_t/_s,
!!                            not the interp result.
!!   zinit_dim_mismatch     — T5: a file with the wrong x/y dims trips the
!!                            zinit_dims_ok validator (no error stop).
!!
!! Oracle (numpy `local_archive/specs/a2_oracle.py`, not committed):
!!   Column: nz_ml=4 layers, H=400 m, uniform h_layer=100 m each.
!!   z_ctr(k=1..4, bed→surface) = [350, 250, 150, 50].
!!   T2: a=10, b=0.02, z_src=[25,75,...,375] ⇒ T=[17,15,13,11],
!!       S=35-0.001*z ⇒ S=[34.65,34.75,34.85,34.95].
!!   T3: z_src=[100,200,300], T_src=[12,14,16] ⇒
!!       T_ex=[16,15,13,12]  (k=1 z=350>300 clamp deep ⇒ 16;
!!       k=4 z=50<100 clamp shallow ⇒ 12; k=2 z=250 ⇒ 15; k=3 z=150 ⇒ 13).
!!
!! Bit-identity note: with `enable=.false.` (the config default), the
!! z-init branch in `ocean_state_seed_from_cfg` is never taken — the
!! analytical IC is produced unchanged.  The full suite staying green
!! with the default-off branch is the bit-identity guarantee; this test
!! exercises only the explicitly-enabled path.
module test_ocean_zinit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_config, only: ocean_zinit_config_t
   use rdb_ocean_z_init, only: seed_ts_from_zfile, interp_column_linear_z, &
                               zinit_dims_ok
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, &
                            nc_def_var_3d, rdb_def_var_1d, nc_enddef, &
                            rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_zinit_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_zinit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zinit_exact_recovery", test_exact_recovery), &
                  new_unittest("zinit_linear", test_linear), &
                  new_unittest("zinit_extrapolation", test_extrapolation), &
                  new_unittest("zinit_dry_column", test_dry_column), &
                  new_unittest("zinit_dim_mismatch", test_dim_mismatch) &
                  ]
   end subroutine collect_ocean_zinit_tests

   ! =================================================================
   ! Test-support: write a tiny (x, y, z) Fortran-ordered z-level file.
   ! =================================================================

   subroutine write_zfile(filename, nx, ny, nz_src, z_src, t_src, s_src)
      !! Write a model-grid T/S NetCDF with Fortran-order dims (x, y, z)
      !! so the reader's first dim name is "x" (no transpose needed).
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nz_src
      real(wp), intent(in) :: z_src(nz_src)
      real(wp), intent(in) :: t_src(nx, ny, nz_src), s_src(nx, ny, nz_src)

      integer :: ncid, dim_x, dim_y, dim_z
      integer :: vid_t, vid_s, vid_z

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "z", nz_src, dim_z)
      call nc_def_var_3d(ncid, "temp", [dim_x, dim_y, dim_z], vid_t)
      call nc_def_var_3d(ncid, "salt", [dim_x, dim_y, dim_z], vid_s)
      call rdb_def_var_1d(ncid, "z_src", dim_z, vid_z)
      call nc_enddef(ncid)
      call nc_put_3d(ncid, vid_t, t_src)
      call nc_put_3d(ncid, vid_s, s_src)
      call rdb_put_var_1d(ncid, vid_z, z_src)
      call nc_close(ncid)
   end subroutine write_zfile

   subroutine nc_put_3d(ncid, varid, data)
      !! Write a full 3D array (no per-time slicing helper exists).
      integer, intent(in) :: ncid, varid
      real(wp), intent(in) :: data(:, :, :)
      integer :: ierr
      ierr = nf90_put_var(ncid, varid, data)
   end subroutine nc_put_3d

   subroutine make_state(ms, nx_phys, ny_phys, nz_ml, depth)
      !! Build a tiny multilayer C-grid state with uniform bathymetry
      !! `depth`, uniform layers (h_layer = depth/nz_ml), and an all-wet
      !! mask.  Tracers (S, T) are registered by `init`.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx_phys, ny_phys, nz_ml
      real(wp), intent(in) :: depth
      type(hgrid_t) :: grid

      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_ml
      call ms%init(grid)
      ms%h_layer = depth/real(nz_ml, wp)
      ms%wet_mask = 1.0_wp
      ms%tracers(ms%idx_salinity)%hTr = 0.0_wp
      ms%tracers(ms%idx_temperature)%hTr = 0.0_wp
   end subroutine make_state

   function make_grid(nx_phys, ny_phys) result(g)
      integer, intent(in) :: nx_phys, ny_phys
      type(hgrid_t) :: g
      call g%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end function make_grid

   ! =================================================================
   ! T1 — exact recovery: source levels AT the layer centres.
   ! =================================================================

   subroutine test_exact_recovery(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 2, NZ_ML = 4, NZ_SRC = 4
      real(wp), parameter :: DEPTH = 400.0_wp
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      type(ocean_zinit_config_t) :: cfg
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)
      ! z_ctr for uniform 4x100m layers = [350,250,150,50] (bed→surf).
      ! Place source levels at the SAME depths (ascending order).
      real(wp), parameter :: ZC(NZ_ML) = [350.0_wp, 250.0_wp, 150.0_wp, 50.0_wp]
      integer :: i, j, k, ng
      real(wp) :: val, expect

      ! source levels ascending = [50,150,250,350]; value = depth itself
      z_src = [50.0_wp, 150.0_wp, 250.0_wp, 350.0_wp]
      do k = 1, NZ_SRC
         t_src(:, :, k) = z_src(k)            ! T = depth
         s_src(:, :, k) = 30.0_wp + z_src(k)  ! S = 30 + depth
      end do

      fname = "/tmp/test_zinit_exact.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      grid = make_grid(NX, NY)
      call make_state(ms, NX, NY, NZ_ML, DEPTH)
      cfg%file = trim(fname)
      call seed_ts_from_zfile(ms, grid, cfg)

      ng = NGHOST
      do j = 1, NY
         do i = 1, NX
            do k = 1, NZ_ML
               ! hTr / h_layer must equal the source value at z_ctr(k).
               val = ms%tracers(ms%idx_temperature)%hTr(ng + i, ng + j, k)/ &
                     ms%h_layer(ng + i, ng + j, k)
               expect = ZC(k)   ! T = depth
               call check(error, abs(val - expect) < 1.0e-10_wp, &
                          "T1 temperature recovery wrong at layer")
               if (allocated(error)) return
               val = ms%tracers(ms%idx_salinity)%hTr(ng + i, ng + j, k)/ &
                     ms%h_layer(ng + i, ng + j, k)
               expect = 30.0_wp + ZC(k)
               call check(error, abs(val - expect) < 1.0e-10_wp, &
                          "T1 salinity recovery wrong at layer")
               if (allocated(error)) return
            end do
         end do
      end do

      call ms%destroy()
   end subroutine test_exact_recovery

   ! =================================================================
   ! T2 — linear-in-z correctness vs the numpy oracle.
   ! =================================================================

   subroutine test_linear(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 2, NZ_ML = 4, NZ_SRC = 8
      real(wp), parameter :: DEPTH = 400.0_wp
      real(wp), parameter :: A = 10.0_wp, B = 0.02_wp
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      type(ocean_zinit_config_t) :: cfg
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)
      ! Oracle expected T/S at z_ctr=[350,250,150,50] (bed→surface):
      real(wp), parameter :: T_EXPECT(NZ_ML) = [17.0_wp, 15.0_wp, 13.0_wp, 11.0_wp]
      real(wp), parameter :: S_EXPECT(NZ_ML) = [34.65_wp, 34.75_wp, 34.85_wp, 34.95_wp]
      integer :: i, j, k, ng
      real(wp) :: val

      z_src = [25.0_wp, 75.0_wp, 125.0_wp, 175.0_wp, &
               225.0_wp, 275.0_wp, 325.0_wp, 375.0_wp]
      do k = 1, NZ_SRC
         t_src(:, :, k) = A + B*z_src(k)
         s_src(:, :, k) = 35.0_wp - 0.001_wp*z_src(k)
      end do

      fname = "/tmp/test_zinit_linear.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      grid = make_grid(NX, NY)
      call make_state(ms, NX, NY, NZ_ML, DEPTH)
      cfg%file = trim(fname)
      call seed_ts_from_zfile(ms, grid, cfg)

      ng = NGHOST
      do j = 1, NY
         do i = 1, NX
            do k = 1, NZ_ML
               val = ms%tracers(ms%idx_temperature)%hTr(ng + i, ng + j, k)/ &
                     ms%h_layer(ng + i, ng + j, k)
               call check(error, abs(val - T_EXPECT(k)) < 1.0e-10_wp, &
                          "T2 linear temperature mismatch vs oracle")
               if (allocated(error)) return
               val = ms%tracers(ms%idx_salinity)%hTr(ng + i, ng + j, k)/ &
                     ms%h_layer(ng + i, ng + j, k)
               call check(error, abs(val - S_EXPECT(k)) < 1.0e-10_wp, &
                          "T2 linear salinity mismatch vs oracle")
               if (allocated(error)) return
            end do
         end do
      end do

      call ms%destroy()
   end subroutine test_linear

   ! =================================================================
   ! T3 — constant extrapolation beyond the source range.
   ! =================================================================

   subroutine test_extrapolation(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 2, NY = 1, NZ_ML = 4, NZ_SRC = 3
      real(wp), parameter :: DEPTH = 400.0_wp
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      type(ocean_zinit_config_t) :: cfg
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)
      ! Oracle: z_src=[100,200,300], T_src=[12,14,16];
      ! z_ctr=[350,250,150,50] ⇒ T_ex=[16,15,13,12].
      real(wp), parameter :: T_EXPECT(NZ_ML) = [16.0_wp, 15.0_wp, 13.0_wp, 12.0_wp]
      integer :: i, j, k, ng
      real(wp) :: val

      z_src = [100.0_wp, 200.0_wp, 300.0_wp]
      do k = 1, NZ_SRC
         t_src(:, :, k) = 10.0_wp + 0.02_wp*z_src(k)
         s_src(:, :, k) = 35.0_wp
      end do

      fname = "/tmp/test_zinit_extrap.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      grid = make_grid(NX, NY)
      call make_state(ms, NX, NY, NZ_ML, DEPTH)
      cfg%file = trim(fname)
      call seed_ts_from_zfile(ms, grid, cfg)

      ng = NGHOST
      do j = 1, NY
         do i = 1, NX
            do k = 1, NZ_ML
               val = ms%tracers(ms%idx_temperature)%hTr(ng + i, ng + j, k)/ &
                     ms%h_layer(ng + i, ng + j, k)
               call check(error, abs(val - T_EXPECT(k)) < 1.0e-10_wp, &
                          "T3 constant-extrapolation mismatch vs oracle")
               if (allocated(error)) return
            end do
         end do
      end do

      call ms%destroy()
   end subroutine test_extrapolation

   ! =================================================================
   ! T4 — dry column skipped (gets land_fill_t/_s).
   ! =================================================================

   subroutine test_dry_column(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 2, NZ_ML = 4, NZ_SRC = 4
      real(wp), parameter :: DEPTH = 400.0_wp
      real(wp), parameter :: FILL_T = 7.5_wp, FILL_S = 33.0_wp
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      type(ocean_zinit_config_t) :: cfg
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)
      integer :: k, ng, di, dj
      real(wp) :: val

      z_src = [50.0_wp, 150.0_wp, 250.0_wp, 350.0_wp]
      do k = 1, NZ_SRC
         t_src(:, :, k) = 20.0_wp     ! distinct from FILL_T
         s_src(:, :, k) = 36.0_wp     ! distinct from FILL_S
      end do

      fname = "/tmp/test_zinit_dry.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      grid = make_grid(NX, NY)
      call make_state(ms, NX, NY, NZ_ML, DEPTH)
      ng = NGHOST
      ! Mark interior column (i=1, j=1) dry.
      di = ng + 1
      dj = ng + 1
      ms%wet_mask(di, dj) = 0.0_wp

      cfg%file = trim(fname)
      cfg%land_fill_t = FILL_T
      cfg%land_fill_s = FILL_S
      call seed_ts_from_zfile(ms, grid, cfg)

      ! Dry column must hold the land-fill constants.
      do k = 1, NZ_ML
         val = ms%tracers(ms%idx_temperature)%hTr(di, dj, k)/ms%h_layer(di, dj, k)
         call check(error, abs(val - FILL_T) < 1.0e-10_wp, &
                    "T4 dry column did not get land_fill_t")
         if (allocated(error)) return
         val = ms%tracers(ms%idx_salinity)%hTr(di, dj, k)/ms%h_layer(di, dj, k)
         call check(error, abs(val - FILL_S) < 1.0e-10_wp, &
                    "T4 dry column did not get land_fill_s")
         if (allocated(error)) return
      end do

      ! A neighbouring WET column must still get the interp result (20.0).
      val = ms%tracers(ms%idx_temperature)%hTr(ng + 2, dj, 1)/ms%h_layer(ng + 2, dj, 1)
      call check(error, abs(val - 20.0_wp) < 1.0e-10_wp, &
                 "T4 wet neighbour did not get the interp value")
      if (allocated(error)) return

      call ms%destroy()
   end subroutine test_dry_column

   ! =================================================================
   ! T5 — dim-mismatch caught by the in-process validator (no error stop).
   ! =================================================================

   subroutine test_dim_mismatch(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 4, NY = 3, NZ_SRC = 2
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)

      z_src = [50.0_wp, 150.0_wp]
      t_src = 10.0_wp
      s_src = 35.0_wp

      fname = "/tmp/test_zinit_mismatch.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      ! Correct query → ok.
      call check(error, zinit_dims_ok(trim(fname), "temp", NX, NY), &
                 "correct dims flagged as mismatch")
      if (allocated(error)) return
      ! x mismatch.
      call check(error,.not. zinit_dims_ok(trim(fname), "temp", NX + 1, NY), &
                 "x mismatch not detected")
      if (allocated(error)) return
      ! y mismatch.
      call check(error,.not. zinit_dims_ok(trim(fname), "temp", NX, NY + 2), &
                 "y mismatch not detected")
      if (allocated(error)) return
      ! missing variable name.
      call check(error,.not. zinit_dims_ok(trim(fname), "nope", NX, NY), &
                 "missing variable not detected")
   end subroutine test_dim_mismatch

end module test_ocean_zinit
