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
!!   zinit_draft_offset     — P5.3a: under a SLOPING draft, every layer
!!                            centre receives T at its TRUE geopotential
!!                            depth, so the isopycnals are flat in z and
!!                            not in the coordinate.
!!   zinit_draft_bitident   — P5.3b: `z_draft == 0` reproduces the
!!                            no-draft answer BIT-for-bit, on both the
!!                            file and the analytic source.
!!   zinit_linear_source    — P5.3c: the analytic `source='linear'` path
!!                            is exact on the profile it claims, fills
!!                            ghosts and land, and is exactly recovered
!!                            by the file reader fed the same profile.
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
                               seed_ts_linear_z, build_z_ctr, zinit_dims_ok
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, &
                            nc_def_var_3d, rdb_def_var_1d, nc_enddef, &
                            rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_zinit_tests

   integer, parameter :: NGHOST = 2

   ! Profile of record for these tests: T(z) = T_REF + DTDZ*z, z positive
   ! UP, so T at depth d is T_REF - DTDZ*d.  DTDZ > 0 is stable.
   real(wp), parameter :: T_REF = 1.0_wp, DTDZ = 2.0e-3_wp
   real(wp), parameter :: S_REF = 34.2_wp, DSDZ = -1.0e-3_wp
   ! Tolerances, derived rather than tuned.  Every assertion below
   ! compares two evaluations of the SAME affine formula (or a linear
   ! interpolation of it, which is exact for an affine function), so the
   ! only error is floating-point rounding: a handful of eps times the
   ! largest intermediate.  The largest intermediate is |DTDZ|*z_ctr <=
   ! 2e-3 * 1000 = 2, so a few eps is ~1e-15.  TOL_EXACT sits 3 decades
   ! above that, which is loose enough to survive any reassociation a
   ! compiler may do and 11 decades below the SIGNAL the test is looking
   ! for (see MISMATCH_FLOOR).
   real(wp), parameter :: TOL_EXACT = 1.0e-12_wp
   ! What the BUG this slice fixes would look like: without the draft
   ! offset, neighbouring columns disagree by DTDZ*(z_draft difference)
   ! at equal geopotential depth.  For the geometry below that is
   ! 2e-3 * 50 = 1e-1 degC per column step.  The "would have caught it"
   ! assertion demands at least a tenth of that.
   real(wp), parameter :: MISMATCH_FLOOR = 1.0e-2_wp

contains

   subroutine collect_ocean_zinit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zinit_exact_recovery", test_exact_recovery), &
                  new_unittest("zinit_linear", test_linear), &
                  new_unittest("zinit_extrapolation", test_extrapolation), &
                  new_unittest("zinit_dry_column", test_dry_column), &
                  new_unittest("zinit_dim_mismatch", test_dim_mismatch), &
                  new_unittest("zinit_draft_offset", test_draft_offset), &
                  new_unittest("zinit_draft_bitident", test_draft_bitident), &
                  new_unittest("zinit_linear_source", test_linear_source) &
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

   ! =================================================================
   ! P5.3 — the draft offset and the analytic linear source.
   !
   ! The geometry these three share: a FLAT bed at `BED`, an ice base
   ! that deepens linearly with i (`z_draft(i) = D0 + DSLOPE*(i-1)`),
   ! and sigma layers that split the LOCAL water column
   ! `(BED - z_draft(i))/nz` — i.e. exactly what
   ! `ocean_state_seed_from_cfg` lays down under `VCOORD_SIGMA` with a
   ! cavity.  The layer interfaces therefore TILT with the ice base, and
   ! a profile laid out per layer index tilts with them.  Only a profile
   ! sampled at the TRUE geopotential depth is flat in z, and only a
   ! stratification that is flat in z is a state of rest.
   ! =================================================================

   subroutine make_sloping_cavity(ms, z_draft, nx, ny, nz_ml, bed, d0, dslope)
      !! Tiny multilayer state under a linearly-deepening ice base, plus
      !! the matching FULL ghosted `z_draft`.  Sigma split of the LOCAL
      !! water column, all wet.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), allocatable, intent(out) :: z_draft(:, :)
      integer, intent(in) :: nx, ny, nz_ml
      real(wp), intent(in) :: bed, d0, dslope
      type(hgrid_t) :: grid
      integer :: i, j, nxt, nyt
      real(wp) :: water

      call grid%init(nx, ny, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_ml
      call ms%init(grid)
      nxt = size(ms%h_layer, 1)
      nyt = size(ms%h_layer, 2)
      allocate (z_draft(nxt, nyt))
      ! Analytic on the FULL array, ghosts included -- the formula
      ! setters' contract (CLAUDE.md ghost-fill gotcha).
      do j = 1, nyt
         do i = 1, nxt
            z_draft(i, j) = d0 + dslope*real(i - NGHOST - 1, wp)
            water = bed - z_draft(i, j)
            ms%h_layer(i, j, :) = water/real(nz_ml, wp)
         end do
      end do
      ms%wet_mask = 1.0_wp
      ms%tracers(ms%idx_salinity)%hTr = 0.0_wp
      ms%tracers(ms%idx_temperature)%hTr = 0.0_wp
   end subroutine make_sloping_cavity

   pure function tracer_value(ms, idx, i, j, k) result(v)
      !! `hTr / h_layer` -- the concentration the seeders wrote.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: idx, i, j, k
      real(wp) :: v
      v = ms%tracers(idx)%hTr(i, j, k)/ms%h_layer(i, j, k)
   end function tracer_value

   pure function interp_at_depth(z_ctr, t_col, nz_ml, z_target) result(v)
      !! Linear interpolation of a column's (depth, value) pairs at
      !! `z_target`, with `z_ctr` DESCENDING in k (k=1 is the deepest).
      !! Exact for an affine profile, which is the whole point: any
      !! discrepancy between two columns is then a geometry error, not
      !! an interpolation error.
      integer, intent(in) :: nz_ml
      real(wp), intent(in) :: z_ctr(nz_ml), t_col(nz_ml), z_target
      real(wp) :: v
      integer :: k
      real(wp) :: w

      v = t_col(nz_ml)
      do k = nz_ml, 2, -1
         ! Layer centres k (deeper) and k-1 ... walking DOWN in k means
         ! walking DEEPER, so the bracket is [z_ctr(k), z_ctr(k-1)].
         if (z_target >= z_ctr(k) .and. z_target <= z_ctr(k - 1)) then
            w = (z_target - z_ctr(k))/(z_ctr(k - 1) - z_ctr(k))
            v = (1.0_wp - w)*t_col(k) + w*t_col(k - 1)
            return
         end if
      end do
      if (z_target < z_ctr(nz_ml)) v = t_col(nz_ml)
      if (z_target > z_ctr(1)) v = t_col(1)
   end function interp_at_depth

   subroutine test_draft_offset(error)
      !! P5.3a.  Under a sloping ice base, every layer centre must
      !! receive `T` at its TRUE geopotential depth, so that isopycnals
      !! are flat in `z` and not in the terrain-following coordinate.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 6, NY = 2, NZ_ML = 5
      real(wp), parameter :: BED = 1000.0_wp, D0 = 200.0_wp, DSLOPE = 50.0_wp
      type(multilayer_state_t) :: ms
      type(ocean_zinit_config_t) :: cfg
      real(wp), allocatable :: z_draft(:, :)
      real(wp) :: zc(NZ_ML), t_col(NZ_ML), zc_r(NZ_ML), t_col_r(NZ_ML)
      real(wp) :: expect, got, z_probe, v_l, v_r, worst, worst_noff
      integer :: i, j, k, ng

      ng = NGHOST
      cfg%source = "linear"
      cfg%lin_t_ref = T_REF
      cfg%lin_dt_dz = DTDZ
      cfg%lin_s_ref = S_REF
      cfg%lin_ds_dz = DSDZ

      ! ---- (1) every layer centre gets T at its geopotential depth ----
      call make_sloping_cavity(ms, z_draft, NX, NY, NZ_ML, BED, D0, DSLOPE)
      call seed_ts_linear_z(ms, cfg, z_draft=z_draft)
      do j = 1, NY
         do i = 1, NX
            call build_z_ctr(ms%h_layer(ng + i, ng + j, :), NZ_ML, &
                             z_draft(ng + i, ng + j), zc)
            do k = 1, NZ_ML
               expect = T_REF - DTDZ*zc(k)
               got = tracer_value(ms, ms%idx_temperature, ng + i, ng + j, k)
               call check(error, abs(got - expect) < TOL_EXACT, &
                          "P5.3a: T must equal the profile at the layer centre's "// &
                          "GEOPOTENTIAL depth")
               if (allocated(error)) return
               expect = S_REF - DSDZ*zc(k)
               got = tracer_value(ms, ms%idx_salinity, ng + i, ng + j, k)
               call check(error, abs(got - expect) < TOL_EXACT, &
                          "P5.3a: S must equal the profile at the layer centre's "// &
                          "GEOPOTENTIAL depth")
               if (allocated(error)) return
            end do
         end do
      end do

      ! ---- (2) the headline: isopycnals FLAT IN z ---------------------
      ! Compare neighbouring columns at the SAME geopotential depth.  The
      ! layer centres do not coincide (the sigma layers tilt with the
      ! lid), so the comparison goes through a linear interpolation in
      ! each column -- exact for this affine profile, so anything left is
      ! a geometry error.  The probe depth is chosen inside BOTH columns'
      ! centre range: the deeper column's shallowest centre is the
      ! binding constraint.
      worst = 0.0_wp
      j = 1
      do i = 1, NX - 1
         call build_z_ctr(ms%h_layer(ng + i, ng + j, :), NZ_ML, &
                          z_draft(ng + i, ng + j), zc)
         call build_z_ctr(ms%h_layer(ng + i + 1, ng + j, :), NZ_ML, &
                          z_draft(ng + i + 1, ng + j), zc_r)
         do k = 1, NZ_ML
            t_col(k) = tracer_value(ms, ms%idx_temperature, ng + i, ng + j, k)
            t_col_r(k) = tracer_value(ms, ms%idx_temperature, ng + i + 1, ng + j, k)
         end do
         do k = 1, NZ_ML
            z_probe = max(zc(NZ_ML), zc_r(NZ_ML)) + &
                      (min(zc(1), zc_r(1)) - max(zc(NZ_ML), zc_r(NZ_ML)))* &
                      real(k - 1, wp)/real(NZ_ML - 1, wp)
            v_l = interp_at_depth(zc, t_col, NZ_ML, z_probe)
            v_r = interp_at_depth(zc_r, t_col_r, NZ_ML, z_probe)
            worst = max(worst, abs(v_l - v_r))
         end do
      end do
      call check(error, worst < TOL_EXACT, &
                 "P5.3a: neighbouring columns must agree at equal GEOPOTENTIAL "// &
                 "depth -- the isopycnals are flat in z, not in the coordinate")
      if (allocated(error)) return

      ! ---- (3) and the same measurement WITHOUT the offset fails ------
      ! This is what the code did before P5.3: depth measured from the
      ! column TOP.  The test must be able to SEE that, or it gates
      ! nothing.
      call seed_ts_linear_z(ms, cfg)
      worst_noff = 0.0_wp
      do i = 1, NX - 1
         call build_z_ctr(ms%h_layer(ng + i, ng + j, :), NZ_ML, &
                          z_draft(ng + i, ng + j), zc)
         call build_z_ctr(ms%h_layer(ng + i + 1, ng + j, :), NZ_ML, &
                          z_draft(ng + i + 1, ng + j), zc_r)
         do k = 1, NZ_ML
            t_col(k) = tracer_value(ms, ms%idx_temperature, ng + i, ng + j, k)
            t_col_r(k) = tracer_value(ms, ms%idx_temperature, ng + i + 1, ng + j, k)
         end do
         do k = 1, NZ_ML
            z_probe = max(zc(NZ_ML), zc_r(NZ_ML)) + &
                      (min(zc(1), zc_r(1)) - max(zc(NZ_ML), zc_r(NZ_ML)))* &
                      real(k - 1, wp)/real(NZ_ML - 1, wp)
            v_l = interp_at_depth(zc, t_col, NZ_ML, z_probe)
            v_r = interp_at_depth(zc_r, t_col_r, NZ_ML, z_probe)
            worst_noff = max(worst_noff, abs(v_l - v_r))
         end do
      end do
      call check(error, worst_noff > MISMATCH_FLOOR, &
                 "P5.3a: dropping the draft offset MUST tilt the isopycnals -- "// &
                 "if it does not, this test is not gating the fix")

      call ms%destroy()
   end subroutine test_draft_offset

   subroutine test_draft_bitident(error)
      !! P5.3b.  A `z_draft` of zeros is the open ocean, and it must
      !! reproduce the no-draft answer BIT-for-bit on both sources.
      !!
      !! Exact equality is the right assertion here and not a violation
      !! of the never-compare-across-paths rule: `build_z_ctr` starts its
      !! accumulator at `z_top`, so `z_top = 0` makes `above` hold the
      !! same IEEE value the original `above = 0.0_wp` held, and every
      !! subsequent operation is the identical expression on identical
      !! operands.  It is the SAME arithmetic, not an equivalent one --
      !! which is exactly the claim this slice has to make about every
      !! shipped non-cavity zinit case.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 2, NZ_ML = 4, NZ_SRC = 6
      real(wp), parameter :: DEPTH = 400.0_wp
      type(multilayer_state_t) :: ms_a, ms_b
      type(hgrid_t) :: grid
      type(ocean_zinit_config_t) :: cfg
      character(len=256) :: fname
      real(wp) :: z_src(NZ_SRC), t_src(NX, NY, NZ_SRC), s_src(NX, NY, NZ_SRC)
      real(wp), allocatable :: zero_draft(:, :)
      integer :: i, j, k, ng
      logical :: same

      ng = NGHOST
      grid = make_grid(NX, NY)

      ! ---- analytic source -------------------------------------------
      cfg%source = "linear"
      cfg%lin_t_ref = T_REF
      cfg%lin_dt_dz = DTDZ
      cfg%lin_s_ref = S_REF
      cfg%lin_ds_dz = DSDZ
      call make_state(ms_a, NX, NY, NZ_ML, DEPTH)
      call make_state(ms_b, NX, NY, NZ_ML, DEPTH)
      allocate (zero_draft(size(ms_b%h_layer, 1), size(ms_b%h_layer, 2)), &
                source=0.0_wp)
      call seed_ts_linear_z(ms_a, cfg)
      call seed_ts_linear_z(ms_b, cfg, z_draft=zero_draft)
      same = all(ms_a%tracers(ms_a%idx_temperature)%hTr == &
                 ms_b%tracers(ms_b%idx_temperature)%hTr) .and. &
         all(ms_a%tracers(ms_a%idx_salinity)%hTr == &
             ms_b%tracers(ms_b%idx_salinity)%hTr)
      call check(error, same, &
                 "P5.3b: a zero draft must be BIT-identical to no draft (analytic)")
      if (allocated(error)) then
         call ms_a%destroy(); call ms_b%destroy(); return
      end if
      call ms_a%destroy()
      call ms_b%destroy()
      deallocate (zero_draft)

      ! ---- file source -----------------------------------------------
      do k = 1, NZ_SRC
         z_src(k) = 25.0_wp + 70.0_wp*real(k - 1, wp)
         t_src(:, :, k) = T_REF - DTDZ*z_src(k)
         s_src(:, :, k) = S_REF - DSDZ*z_src(k)
      end do
      fname = "test_zinit_draft_bitident.nc"
      call write_zfile(trim(fname), NX, NY, NZ_SRC, z_src, t_src, s_src)

      cfg%source = "file"
      cfg%file = trim(fname)
      call make_state(ms_a, NX, NY, NZ_ML, DEPTH)
      call make_state(ms_b, NX, NY, NZ_ML, DEPTH)
      allocate (zero_draft(size(ms_b%h_layer, 1), size(ms_b%h_layer, 2)), &
                source=0.0_wp)
      call seed_ts_from_zfile(ms_a, grid, cfg)
      call seed_ts_from_zfile(ms_b, grid, cfg, z_draft=zero_draft)
      same = .true.
      do k = 1, NZ_ML
         do j = 1, NY
            do i = 1, NX
               if (ms_a%tracers(ms_a%idx_temperature)%hTr(ng + i, ng + j, k) /= &
                   ms_b%tracers(ms_b%idx_temperature)%hTr(ng + i, ng + j, k)) same = .false.
               if (ms_a%tracers(ms_a%idx_salinity)%hTr(ng + i, ng + j, k) /= &
                   ms_b%tracers(ms_b%idx_salinity)%hTr(ng + i, ng + j, k)) same = .false.
            end do
         end do
      end do
      call check(error, same, &
                 "P5.3b: a zero draft must be BIT-identical to no draft (file)")
      if (allocated(error)) then
         call ms_a%destroy(); call ms_b%destroy(); return
      end if

      ! ---- and build_z_ctr shifts by EXACTLY z_top --------------------
      block
         real(wp) :: zc0(NZ_ML), zcd(NZ_ML)
         real(wp), parameter :: Z_TOP = 137.5_wp
         call build_z_ctr(ms_a%h_layer(ng + 1, ng + 1, :), NZ_ML, 0.0_wp, zc0)
         call build_z_ctr(ms_a%h_layer(ng + 1, ng + 1, :), NZ_ML, Z_TOP, zcd)
         call check(error, maxval(abs((zcd - zc0) - Z_TOP)) < TOL_EXACT, &
                    "P5.3b: build_z_ctr must shift every layer centre by exactly z_top")
      end block

      call ms_a%destroy()
      call ms_b%destroy()
   end subroutine test_draft_bitident

   subroutine test_linear_source(error)
      !! P5.3c.  The analytic source is exact on the profile it claims,
      !! and -- unlike the file path, which has no data for them -- it
      !! fills the GHOST rows and the DRY columns too, because the
      !! formula is defined there and a ghost column left on some other
      !! profile is the wall-adjacent spurious-density-jump gotcha.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 3, NY = 3, NZ_ML = 4
      real(wp), parameter :: DEPTH = 400.0_wp
      type(multilayer_state_t) :: ms
      type(ocean_zinit_config_t) :: cfg
      real(wp) :: zc(NZ_ML), expect, got
      integer :: i, j, k, nxt, nyt, ng

      ng = NGHOST
      cfg%source = "linear"
      cfg%lin_t_ref = T_REF
      cfg%lin_dt_dz = DTDZ
      cfg%lin_s_ref = S_REF
      cfg%lin_ds_dz = DSDZ
      cfg%land_fill_t = -999.0_wp
      cfg%land_fill_s = -999.0_wp

      call make_state(ms, NX, NY, NZ_ML, DEPTH)
      ! One DRY column: the analytic path must still write the profile
      ! there (no land_fill), so a wall-adjacent EOS evaluation sees a
      ! continuous density field.
      ms%wet_mask(ng + 2, ng + 2) = 0.0_wp
      call seed_ts_linear_z(ms, cfg)

      nxt = size(ms%h_layer, 1)
      nyt = size(ms%h_layer, 2)
      do j = 1, nyt
         do i = 1, nxt
            call build_z_ctr(ms%h_layer(i, j, :), NZ_ML, 0.0_wp, zc)
            do k = 1, NZ_ML
               expect = T_REF - DTDZ*zc(k)
               got = tracer_value(ms, ms%idx_temperature, i, j, k)
               call check(error, abs(got - expect) < TOL_EXACT, &
                          "P5.3c: the analytic profile must be written on the FULL "// &
                          "array -- ghosts and dry columns included")
               if (allocated(error)) return
               expect = S_REF - DSDZ*zc(k)
               got = tracer_value(ms, ms%idx_salinity, i, j, k)
               call check(error, abs(got - expect) < TOL_EXACT, &
                          "P5.3c: the analytic S profile must be written on the "// &
                          "FULL array")
               if (allocated(error)) return
            end do
         end do
      end do

      ! A zero gradient is a uniform column, exactly.
      cfg%lin_dt_dz = 0.0_wp
      call seed_ts_linear_z(ms, cfg)
      call check(error, all(abs(ms%tracers(ms%idx_temperature)%hTr - &
                                T_REF*ms%h_layer) < TOL_EXACT), &
                 "P5.3c: lin_dt_dz = 0 must give a uniform column at lin_t_ref")

      call ms%destroy()
   end subroutine test_linear_source

end module test_ocean_zinit
