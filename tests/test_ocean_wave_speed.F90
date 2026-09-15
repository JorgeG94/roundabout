!! Unit + analytical tests for the first-baroclinic wave speed (B1,
!! `rdb_ocean_wave_speed`).  Golden values come from the verified
!! prototype `local_archive/prototypes/b1_wavespeed_prototype.py`
!! (run `--report` for the nz=4 worked example); spec is
!! `local_archive/specs/b1_wavespeed_spec.md`.
!!
!! Cases:
!!   T1 — two-layer reduced gravity: cg1 = sqrt(g' H1 H2 / H) < 1%
!!        across an H1/H2, H, drho sweep.
!!   T2 — flat-N^2 WKB limit: cg1 -> N H / pi as nz grows.
!!   T3 — monotone cg1/Rd in stratification + depth; equatorial Rd
!!        finite and == sqrt(cg1/(2 beta)).
!!   T4 — degenerate columns (homogeneous + statically unstable +
!!        single layer): cg1 = 0 exactly.
!!   T5 — golden nz=4 column == prototype's 3.266407 m/s to ~1e-6.
!!   T6 — partial mid-column inversion: cg1 matches the MOM6 merged
!!        value (2.574131), NOT 0 and NOT 1e10 — the backtracking-merge
!!        regression gate for the closed CRITICAL defect.
!!   T7 — wavespeed_compute (the wire): fills cg1 over a stratified
!!        Cartesian column to the T5 golden, cg1=0 at a land cell.  Runs
!!        the DEVICE kernel via the shim; the mem:separate canary.
!!   T8 — rd_over_dx uses metrics%dxT (metres), NOT grid%dx (which is
!!        degrees on spherical) — metric bug 1.
!!   T9 — f_centre / beta_centre from the planetary Coriolis path on an
!!        equator-crossing spherical sector; equatorial Rd matches the
!!        Gill form — metric bug 2 + its tail.
!!   T10 — beta-plane build_static is bit-identical to the legacy
!!        set_f_centre; beta_centre == the namelist beta in the interior.
module test_ocean_wave_speed
   use testdrive, only: new_unittest, unittest_type, error_type, check
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY
#endif
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_finalize, metrics_fill_coriolis, &
                                CORIOLIS_SCHEME_PLANETARY, CORIOLIS_SCHEME_BETA_PLANE
   use ocean_test_metrics, only: make_cartesian_metrics, make_spherical_metrics, &
                                 destroy_cartesian_metrics
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t, wavespeed_compute, &
                                   wavespeed_cg1_column, wavespeed_rd
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 128
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=128).
#endif

   public :: collect_ocean_wave_speed_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_wave_speed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wavespeed_two_layer_reduced_gravity", test_two_layer), &
                  new_unittest("wavespeed_flat_n2_wkb_limit", test_wkb), &
                  new_unittest("wavespeed_monotone_and_equatorial_rd", test_monotone_rd), &
                  new_unittest("wavespeed_degenerate_columns", test_degenerate), &
                  new_unittest("wavespeed_golden_nz4_column", test_golden_nz4), &
                  new_unittest("wavespeed_partial_inversion_regression", test_inversion), &
                  new_unittest("wavespeed_compute_fills_cg1_on_stratified_column", &
                               test_compute_fills_cg1), &
                  new_unittest("wavespeed_rd_over_dx_uses_metres_not_grid_units", &
                               test_rd_over_dx_metres), &
                  new_unittest("wavespeed_f_centre_and_beta_from_planetary_coriolis", &
                               test_f_centre_beta_planetary), &
                  new_unittest("wavespeed_beta_plane_build_static_bit_identical", &
                               test_beta_plane_bit_identical) &
                  ]
   end subroutine collect_ocean_wave_speed_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   !! Pack a Roundabout-ordered (k=1 bed) column into the NZ_STACK_MAX
   !! fixed-size arrays + call the column solver on the host.
   subroutine cg1_of(nz, h_rak, rho_rak, cg1)
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_rak(nz), rho_rak(nz)
      real(wp), intent(out) :: cg1
      real(wp) :: h_col(NZ_STACK_MAX), rho_col(NZ_STACK_MAX)
      integer :: k
      h_col = 0.0_wp
      rho_col = 0.0_wp
      do k = 1, nz
         h_col(k) = h_rak(k)
         rho_col(k) = rho_rak(k)
      end do
      call wavespeed_cg1_column(nz, h_col, rho_col, RHO0, cg1)
   end subroutine cg1_of

   ! -----------------------------------------------------------------
   ! T1 — two-layer reduced gravity
   ! -----------------------------------------------------------------
   subroutine test_two_layer(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h1h2_set(3), h_set(3), drho_set(3)
      real(wp) :: h1h2, h, drho, h1, h2, gp, c2_ex, cg1_ex, cg1_got, relerr
      real(wp) :: h_rak(2), rho_rak(2)
      integer :: a, b, c

      h1h2_set = [0.1_wp, 1.0_wp, 10.0_wp]
      h_set = [200.0_wp, 1000.0_wp, 4000.0_wp]
      drho_set = [0.5_wp, 2.0_wp, 8.0_wp]

      do a = 1, 3
         do b = 1, 3
            do c = 1, 3
               h1h2 = h1h2_set(a)
               h = h_set(b)
               drho = drho_set(c)
               h2 = h/(1.0_wp + h1h2)
               h1 = h - h2
               gp = GRAVITY*drho/RHO0
               c2_ex = gp*h1*h2/h
               cg1_ex = sqrt(c2_ex)
               ! Roundabout ordering: k=1=bed (H2, denser), k=2=surf (H1, lighter)
               h_rak = [h2, h1]
               rho_rak = [RHO0 + drho, RHO0]
               call cg1_of(2, h_rak, rho_rak, cg1_got)
               relerr = abs(cg1_got - cg1_ex)/max(abs(cg1_ex), 1.0e-300_wp)
               call check(error, relerr < 0.01_wp, &
                          "two-layer cg1 rel-err >= 1%")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_two_layer

   ! -----------------------------------------------------------------
   ! T2 — flat-N^2 WKB limit: cg1 -> N H / pi, error shrinks with nz
   ! -----------------------------------------------------------------
   subroutine test_wkb(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HCOL = 1000.0_wp, NBV = 0.01_wp
      real(wp) :: cg1_wkb, dz, cg1_got, err10, err40, pi
      real(wp) :: h_rak(NZ_STACK_MAX), rho_rak(NZ_STACK_MAX)
      integer :: nz, k

      pi = 4.0_wp*atan(1.0_wp)
      cg1_wkb = NBV*HCOL/pi

      ! nz = 10
      nz = 10
      dz = HCOL/real(nz, wp)
      do k = 1, nz
         h_rak(k) = dz
         ! surface-down rho linear; flip to rdb (k=1 bed = densest)
         rho_rak(nz + 1 - k) = RHO0 + (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
      end do
      call cg1_of(nz, h_rak(1:nz), rho_rak(1:nz), cg1_got)
      err10 = abs(cg1_got - cg1_wkb)/cg1_wkb

      ! nz = 40
      nz = 40
      dz = HCOL/real(nz, wp)
      do k = 1, nz
         h_rak(k) = dz
         rho_rak(nz + 1 - k) = RHO0 + (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
      end do
      call cg1_of(nz, h_rak(1:nz), rho_rak(1:nz), cg1_got)
      err40 = abs(cg1_got - cg1_wkb)/cg1_wkb

      call check(error, err40 < 0.05_wp, "WKB nz=40 error >= 5%")
      if (allocated(error)) return
      call check(error, err40 < err10, "WKB not converging with nz")
   end subroutine test_wkb

   ! -----------------------------------------------------------------
   ! T3 — monotone cg1/Rd + equatorial Rd branch
   ! -----------------------------------------------------------------
   subroutine test_monotone_rd(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp), parameter :: HCOL = 500.0_wp, F_MID = 1.0e-4_wp
      real(wp), parameter :: BETA = 2.3e-11_wp, NBV = 0.01_wp
      real(wp) :: drho_set(6), h_test_set(4)
      real(wp) :: dz, cg1, cg1_prev, rd
      real(wp) :: h_rak(NZ), rho_rak(NZ)
      real(wp) :: cg1_eq, rd_eq, rd_eq_expected
      integer :: a, k

      drho_set = [0.1_wp, 0.5_wp, 1.0_wp, 2.0_wp, 5.0_wp, 10.0_wp]
      h_test_set = [200.0_wp, 500.0_wp, 1000.0_wp, 2000.0_wp]
      dz = HCOL/real(NZ, wp)

      ! stratification sweep — cg1 strictly increasing
      cg1_prev = -1.0_wp
      do a = 1, 6
         do k = 1, NZ
            h_rak(k) = dz
            ! surface-down linspace(rho0, rho0+drho); flip to rdb
            rho_rak(NZ + 1 - k) = RHO0 + drho_set(a)*real(k - 1, wp)/real(NZ - 1, wp)
         end do
         call cg1_of(NZ, h_rak, rho_rak, cg1)
         call check(error, cg1 > cg1_prev, "cg1 not increasing with stratification")
         if (allocated(error)) return
         cg1_prev = cg1
      end do

      ! depth sweep — cg1 strictly increasing
      cg1_prev = -1.0_wp
      do a = 1, 4
         dz = h_test_set(a)/real(NZ, wp)
         do k = 1, NZ
            h_rak(k) = dz
            rho_rak(NZ + 1 - k) = RHO0 + (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
         end do
         call cg1_of(NZ, h_rak, rho_rak, cg1)
         call check(error, cg1 > cg1_prev, "cg1 not increasing with depth")
         if (allocated(error)) return
         cg1_prev = cg1
      end do

      ! equatorial Rd branch: f -> 0, Rd = sqrt(cg1/(2 beta)), finite
      dz = 100.0_wp
      do k = 1, NZ
         h_rak(k) = dz
         rho_rak(NZ + 1 - k) = RHO0 + (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
      end do
      call cg1_of(NZ, h_rak, rho_rak, cg1_eq)
      ! Rd = cg1 / sqrt(f^2 + 2 beta cg1 + guard), f = 0 — via the SHIPPED
      ! `wavespeed_rd` (now public), not a copy-pasted formula: this test
      ! pins the actual production function rather than another
      ! independent restatement of it.
      rd_eq = wavespeed_rd(cg1_eq, 0.0_wp, BETA)
      rd_eq_expected = sqrt(cg1_eq/(2.0_wp*BETA))
      call check(error, rd_eq > 0.0_wp .and. rd_eq < 1.0e9_wp, &
                 "equatorial Rd not finite")
      if (allocated(error)) return
      call check(error, abs(rd_eq - rd_eq_expected)/rd_eq_expected < 1.0e-3_wp, &
                 "equatorial Rd != sqrt(cg1/2beta)")
   end subroutine test_monotone_rd

   ! -----------------------------------------------------------------
   ! T4 — degenerate columns: cg1 = 0 exactly
   ! -----------------------------------------------------------------
   subroutine test_degenerate(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp) :: dz, cg1
      real(wp) :: h_rak(NZ), rho_rak(NZ)
      real(wp) :: h1(1), r1(1)
      integer :: k

      dz = 500.0_wp/real(NZ, wp)

      ! homogeneous
      do k = 1, NZ
         h_rak(k) = dz
         rho_rak(k) = RHO0
      end do
      call cg1_of(NZ, h_rak, rho_rak, cg1)
      call check(error, cg1 == 0.0_wp, "homogeneous cg1 /= 0")
      if (allocated(error)) return

      ! statically unstable: denser at surface (rho increases bed->surface)
      do k = 1, NZ
         h_rak(k) = dz
         rho_rak(k) = RHO0 + 2.0_wp*real(k - 1, wp)/real(NZ - 1, wp)
      end do
      call cg1_of(NZ, h_rak, rho_rak, cg1)
      call check(error, cg1 == 0.0_wp, "statically-unstable cg1 /= 0")
      if (allocated(error)) return

      ! single layer
      h1 = [500.0_wp]
      r1 = [RHO0]
      call cg1_of(1, h1, r1, cg1)
      call check(error, cg1 == 0.0_wp, "single-layer cg1 /= 0")
   end subroutine test_degenerate

   ! -----------------------------------------------------------------
   ! T5 — golden nz=4 column (prototype --report: cg1 = 3.266407 m/s)
   ! -----------------------------------------------------------------
   subroutine test_golden_nz4(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp), parameter :: HCOL = 1000.0_wp, NBV = 0.01_wp
      real(wp), parameter :: CG1_GOLD = 3.266407_wp
      real(wp) :: dz, cg1
      real(wp) :: h_rak(NZ), rho_rak(NZ)
      integer :: k

      dz = HCOL/real(NZ, wp)
      ! surface-down rho(k_loc) = rho0 + (N^2 rho0 / g)(k-0.5)dz; flip to rdb.
      ! Same GRAVITY used in the rho construction AND the gprime solve, so cg1
      ! matches the prototype golden built with G=9.80665 (== GRAVITY here).
      do k = 1, NZ
         h_rak(k) = dz
         rho_rak(NZ + 1 - k) = RHO0 + (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
      end do
      call cg1_of(NZ, h_rak, rho_rak, cg1)
      call check(error, abs(cg1 - CG1_GOLD) < 1.0e-5_wp, &
                 "golden nz=4 cg1 != 3.266407 m/s")
   end subroutine test_golden_nz4

   ! -----------------------------------------------------------------
   ! T6 — partial mid-column inversion (CRITICAL regression gate).
   ! Stable nz=11 ramp with ONE inverted interior interface; the
   ! backtracking merge must weld it so cg1 matches the MOM6 merged
   ! value 2.574131 — NOT 0 and NOT ~1e10.
   ! -----------------------------------------------------------------
   subroutine test_inversion(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 11
      real(wp), parameter :: CG1_MOM = 2.574131_wp
      real(wp) :: rho_sd(NZ), cg1, relerr
      real(wp) :: h_rak(NZ), rho_rak(NZ)
      integer :: k

      ! surface-down stable ramp 0..5, then invert layer index 5 (0-based)
      ! by making it denser than the layer below: rho_sd(6) = rho_sd(7)+1.5.
      do k = 1, NZ
         rho_sd(k) = RHO0 + 5.0_wp*real(k - 1, wp)/real(NZ - 1, wp)
      end do
      rho_sd(6) = rho_sd(7) + 1.5_wp        ! inv@5 (0-based) dip=1.5

      ! flip surface-down -> rdb (k=1 bed), uniform dz = 100 m
      do k = 1, NZ
         h_rak(k) = 100.0_wp
         rho_rak(NZ + 1 - k) = rho_sd(k)
      end do
      call cg1_of(NZ, h_rak, rho_rak, cg1)

      call check(error, cg1 > 0.0_wp .and. cg1 < 1.0e6_wp, &
                 "inversion cg1 blew up (zero-gprime-row defect) or went to 0")
      if (allocated(error)) return
      relerr = abs(cg1 - CG1_MOM)/CG1_MOM
      call check(error, relerr < 1.0e-3_wp, "inversion cg1 != MOM6 merged value")
   end subroutine test_inversion

   ! -----------------------------------------------------------------
   ! Device-map helpers for the T7-T10 integration tests (the
   ! `mem:separate` contract — every array `wavespeed_compute_impl`
   ! touches must be explicitly present on device).  Mirrors
   ! `test_ocean_varmix.F90:map_in`/`map_out`; `metrics` is mapped
   ! separately by `make_cartesian_metrics`/`make_spherical_metrics`.
   ! -----------------------------------------------------------------
   subroutine map_in_ws(ms, ws)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_wave_speed_t), intent(inout) :: ws
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ws)
      call ws%enter_data()
   end subroutine map_in_ws

   subroutine map_out_ws(ms, ws)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_wave_speed_t), intent(inout) :: ws
      !$acc update self(ws%cg1, ws%rd, ws%rd_over_dx)
      call ws%exit_data()
      !$acc exit data delete(ws)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_ws

   ! -----------------------------------------------------------------
   ! T7 — wavespeed_compute fills cg1 on a stratified column (the wire).
   ! Runs the DEVICE kernel; the mem:separate canary for a kernel that
   ! has never executed on a GPU before this PR.
   ! -----------------------------------------------------------------
   subroutine test_compute_fills_cg1(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_wave_speed_t) :: ws
      integer, parameter :: NX = 5, NY = 4, NZ = 4
      real(wp), parameter :: DX = 10000.0_wp, HCOL = 1000.0_wp, NBV = 0.01_wp
      real(wp), parameter :: CG1_GOLD = 3.266407_wp
      real(wp), allocatable :: f_centre(:, :)
      real(wp) :: dz
      integer :: i, j, k, land_i, land_j

      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call make_cartesian_metrics(metrics, grid)

         ms%nz_ml = NZ
         call ms%init(grid)
         dz = HCOL/real(NZ, wp)
         ms%h_layer = dz
         ms%wet_mask = 1.0_wp
         ! One land cell (interior, away from the ghost band): cg1 must
         ! stay exactly 0 there.
         land_i = grid%nghost + 2
         land_j = grid%nghost + 2
         ms%wet_mask(land_i, land_j) = 0.0_wp
         ! Same surface-down stratified column as T5's golden, flipped to
         ! rdb (k=1 bed), broadcast to every column.
         do k = 1, NZ
            ms%rho_layer(:, :, NZ + 1 - k) = RHO0 + &
                                             (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
         end do

         allocate (f_centre(grid%nx_total, grid%ny_total), source=1.0e-4_wp)
         call ws%init(grid)
         ws%enable = .true.
         ws%rho0 = RHO0
         call ws%build_static(grid, metrics, f_centre)

         call map_in_ws(ms, ws)
         call wavespeed_compute(grid, metrics, ws, ms)
         call map_out_ws(ms, ws)

         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               if (i == land_i .and. j == land_j) then
                  call check(error, ws%cg1(i, j) == 0.0_wp, &
                             "land cell cg1 must be exactly 0")
               else
                  call check(error, abs(ws%cg1(i, j) - CG1_GOLD) < 1.0e-5_wp, &
                             "wavespeed_compute cg1 != T5 golden 3.266407 m/s")
               end if
               if (allocated(error)) exit checks
            end do
         end do

         call ws%destroy()
         call ms%destroy()
         call destroy_cartesian_metrics(metrics)
         deallocate (f_centre)
      end block checks
   end subroutine test_compute_fills_cg1

   ! -----------------------------------------------------------------
   ! T8 — rd_over_dx must use metrics%dxT (metres), NOT grid%dx (metric
   ! bug 1: grid%dx is degrees on spherical/supergrid/tripolar).
   ! -----------------------------------------------------------------
   subroutine test_rd_over_dx_metres(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_wave_speed_t) :: ws
      integer, parameter :: NX = 4, NY = 4, NZ = 4
      real(wp), parameter :: DXT_M = 20000.0_wp, HCOL = 1000.0_wp, NBV = 0.01_wp
      real(wp), allocatable :: f_centre(:, :)
      real(wp) :: dz, rd_over_dx_a, rd_over_dx_b, rd_a, rd_b, dx_consistent
      integer :: i0, j0, k

      checks: block
         ! --- (a) consistent metres grid: grid%dx == dxT == 20000 m. ---
         call grid%init(NX, NY, NGHOST, DXT_M, DXT_M)
         call metrics%init(grid)
         call metrics_fill_cartesian(metrics, grid, DXT_M, DXT_M)
         call metrics_finalize(metrics)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()

         ms%nz_ml = NZ
         call ms%init(grid)
         dz = HCOL/real(NZ, wp)
         ms%h_layer = dz
         ms%wet_mask = 1.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, NZ + 1 - k) = RHO0 + &
                                             (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
         end do

         allocate (f_centre(grid%nx_total, grid%ny_total), source=1.0e-4_wp)
         call ws%init(grid)
         ws%enable = .true.
         ws%rho0 = RHO0
         call ws%build_static(grid, metrics, f_centre)

         call map_in_ws(ms, ws)
         call wavespeed_compute(grid, metrics, ws, ms)
         call map_out_ws(ms, ws)

         i0 = grid%nghost + 2
         j0 = grid%nghost + 2
         rd_over_dx_a = ws%rd_over_dx(i0, j0)
         rd_a = ws%rd(i0, j0)
         dx_consistent = grid%dx

         call ws%destroy()
         call ms%destroy()
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call metrics%destroy()
         deallocate (f_centre)

         ! --- (b) same physics, but grid%dx left at a DEGREES value
         ! (0.18) while metrics%dxT is still the real 20000 m — the
         ! units mismatch the spherical path actually produces (dx in
         ! degrees, e.g. acc_channel_kitchensink.nml dx=0.2). ---
         call grid%init(NX, NY, NGHOST, 0.18_wp, 0.18_wp)
         call metrics%init(grid)
         call metrics_fill_cartesian(metrics, grid, DXT_M, DXT_M)
         call metrics_finalize(metrics)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()

         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = dz
         ms%wet_mask = 1.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, NZ + 1 - k) = RHO0 + &
                                             (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
         end do

         allocate (f_centre(grid%nx_total, grid%ny_total), source=1.0e-4_wp)
         call ws%init(grid)
         ws%enable = .true.
         ws%rho0 = RHO0
         call ws%build_static(grid, metrics, f_centre)

         call map_in_ws(ms, ws)
         call wavespeed_compute(grid, metrics, ws, ms)
         call map_out_ws(ms, ws)

         rd_over_dx_b = ws%rd_over_dx(i0, j0)
         rd_b = ws%rd(i0, j0)

         call ws%destroy()
         call ms%destroy()
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call metrics%destroy()
         deallocate (f_centre)

         ! The ratio must depend ONLY on the metric, not on grid%dx: (a)
         ! and (b) share the same dxT=20000 m and must give the SAME
         ! rd_over_dx despite grid%dx differing by ~1.1e5x.
         call check(error, rd_over_dx_a == rd_over_dx_b, &
                    "rd_over_dx depends on grid%dx, not just metrics%dxT")
         if (allocated(error)) exit checks
         call check(error, rd_a == rd_b, "rd itself must be grid%dx-independent")
         if (allocated(error)) exit checks
         ! Bit-identity leg: on THIS Cartesian grid dxT == grid%dx (case
         ! a), so rd_over_dx must equal rd/grid%dx exactly.
         call check(error, rd_over_dx_a == rd_a/dx_consistent, &
                    "rd_over_dx != rd/grid%dx on a pure Cartesian grid")
      end block checks
   end subroutine test_rd_over_dx_metres

   ! -----------------------------------------------------------------
   ! T9 — f_centre / beta_centre from the planetary Coriolis path on an
   ! equator-crossing spherical sector (metric bug 2 + its tail).
   ! -----------------------------------------------------------------
   subroutine test_f_centre_beta_planetary(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_wave_speed_t) :: ws
      integer, parameter :: NX = 6, NY = 8, NZ = 6
      real(wp), parameter :: DLON = 0.5_wp, DLAT = 0.05_wp, LON_WEST = 0.0_wp
      real(wp), parameter :: LAT_SOUTH = -0.19_wp
         !! Fine dlat + an OFF-CENTRE south edge: no T-row sits exactly on
         !! lat=0 (deliberately — a row exactly astride the equator makes
         !! the |f|-gradient stencil see a symmetric cancellation, |f(+h)|
         !! == |f(-h)|, so beta_centre would spuriously read exactly 0 at
         !! that one row: an artifact of the documented |grad|f|| kink, NOT
         !! a bug).  The nearest row (physical row 4, lat=-0.015 deg) sits
         !! close enough to the equator that f_centre^2 is negligible next
         !! to 2*beta_centre*cg1, so the Gill-form self-consistency check
         !! below holds to high precision regardless of the kink.
      real(wp), parameter :: REARTH = 6371000.0_wp
      real(wp), parameter :: OMEGA = 7.292115e-5_wp
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: D2R = PI/180.0_wp
      real(wp), parameter :: HCOL = 1000.0_wp, NBV = 0.01_wp
      real(wp), allocatable :: f_centre(:, :), f_corner(:, :)
      real(wp) :: dz, expected_f, rel, expected_beta, gill, lat_deg
      real(wp) :: best_abs_lat, this_abs_lat
      integer :: i, j, ic, jc, jeq
      logical :: found_eq

      checks: block
         call grid%init(NX, NY, NGHOST, DLON, DLAT)
         call make_spherical_metrics(metrics, grid, LON_WEST, LAT_SOUTH, &
                                     DLON, DLAT, REARTH)

         allocate (f_centre(grid%nx_total, grid%ny_total))
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1))
         call metrics_fill_coriolis(metrics, CORIOLIS_SCHEME_PLANETARY, &
                                    0.0_wp, 0.0_wp, 0.0_wp, OMEGA, grid, &
                                    f_corner, f_centre)

         ms%nz_ml = NZ
         call ms%init(grid)
         dz = HCOL/real(NZ, wp)
         ms%h_layer = dz
         ms%wet_mask = 1.0_wp
         block
            integer :: k
            do k = 1, NZ
               ms%rho_layer(:, :, NZ + 1 - k) = RHO0 + &
                                                (NBV*NBV*RHO0/GRAVITY)*(real(k, wp) - 0.5_wp)*dz
            end do
         end block

         call ws%init(grid)
         ws%enable = .true.
         ws%rho0 = RHO0
         call ws%build_static(grid, metrics, f_centre)

         ! `metrics` is already device-mapped by `make_spherical_metrics`;
         ! only `ms`/`ws` need mapping here.
         call map_in_ws(ms, ws)
         call wavespeed_compute(grid, metrics, ws, ms)
         !$acc update self(ws%f_centre, ws%beta_centre)
         call map_out_ws(ms, ws)

         ! (i) planetary dispatch: f_centre == |2*Omega*sin(geolat)|.
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               expected_f = abs(2.0_wp*OMEGA*sin(metrics%geolatT(i, j)*D2R))
               call check(error, abs(ws%f_centre(i, j) - expected_f) < 1.0e-12_wp, &
                          "f_centre != |2*Omega*sin(geolat)| (planetary dispatch)")
               if (allocated(error)) exit checks
            end do
         end do

         ! (ii) beta_centre ~= 2*Omega*cos(geolat)/Rearth AWAY from the
         ! equator (physical row 2, lat=-1.0 deg).  NOT at the equator:
         ! beta_centre is the gradient of |f|, which has a kink there (the
         ! documented |grad|f|| divergence) — |f(+dlat)| - |f(-dlat)| = 0
         ! by symmetry at f=0, so an equatorial sample point would spuriously
         ! read beta_centre ~ 0 regardless of correctness.
         ic = grid%nghost + NX/2
         jc = grid%nghost + 2
         lat_deg = metrics%geolatT(ic, jc)
         expected_beta = 2.0_wp*OMEGA*cos(lat_deg*D2R)/REARTH
         rel = abs(ws%beta_centre(ic, jc) - expected_beta)/expected_beta
         call check(error, rel < 1.0e-3_wp, &
                    "beta_centre != 2*Omega*cos(geolat)/Rearth (not a real field)")
         if (allocated(error)) exit checks

         ! (iii) equatorial (nearest-to-equator) row: rd finite and matches
         ! the Gill form sqrt(cg1/(2*beta_centre)) — the leg that fails if
         ! f_centre is fixed but beta stays a namelist scalar (Rd would be
         ! ~3e10 m).  Deliberately the row NEAREST lat=0, not one forced
         ! exactly onto it (see the LAT_SOUTH/DLAT comment above).
         found_eq = .false.
         jeq = 0
         best_abs_lat = huge(1.0_wp)
         do j = 1, grid%ny_total
            this_abs_lat = abs(metrics%geolatT(grid%nghost + 1, j))
            if (this_abs_lat < best_abs_lat) then
               best_abs_lat = this_abs_lat
               jeq = j
               found_eq = .true.
            end if
         end do
         call check(error, found_eq, "no equatorial T-row found in test grid")
         if (allocated(error)) exit checks

         ic = grid%nghost + 1
         call check(error, ws%rd(ic, jeq) > 0.0_wp .and. ws%rd(ic, jeq) < 1.0e9_wp, &
                    "equatorial Rd not finite (namelist-scalar-beta regression)")
         if (allocated(error)) exit checks
         gill = sqrt(ws%cg1(ic, jeq)/(2.0_wp*ws%beta_centre(ic, jeq)))
         rel = abs(ws%rd(ic, jeq) - gill)/gill
         call check(error, rel < 1.0e-3_wp, &
                    "equatorial Rd != Gill form sqrt(cg1/2beta)")

         call ws%destroy()
         call ms%destroy()
         call destroy_cartesian_metrics(metrics)
         deallocate (f_centre, f_corner)
      end block checks
   end subroutine test_f_centre_beta_planetary

   ! -----------------------------------------------------------------
   ! T10 — beta-plane build_static bit-identity with the legacy
   ! set_f_centre; beta_centre == the namelist beta in the interior.
   ! -----------------------------------------------------------------
   subroutine test_beta_plane_bit_identical(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_wave_speed_t) :: ws
      integer, parameter :: NX = 6, NY = 6
      real(wp), parameter :: DX = 10000.0_wp
      real(wp), parameter :: F0 = 5.0e-5_wp, BETA = 2.0e-11_wp, YREF = 12345.0_wp
      real(wp), allocatable :: f_centre(:, :), f_corner(:, :)
      real(wp) :: y, expected_f, rel
      integer :: i, j, ng, ic, jc

      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call metrics%init(grid)
         call metrics_fill_cartesian(metrics, grid, DX, DX)
         call metrics_finalize(metrics)

         allocate (f_centre(grid%nx_total, grid%ny_total))
         allocate (f_corner(grid%nx_total + 1, grid%ny_total + 1))
         call metrics_fill_coriolis(metrics, CORIOLIS_SCHEME_BETA_PLANE, &
                                    F0, BETA, YREF, 0.0_wp, grid, f_corner, f_centre)

         call ws%init(grid)
         ws%enable = .true.
         call ws%build_static(grid, metrics, f_centre)

         ng = grid%nghost
         do j = 1, grid%ny_total
            y = (real(j - ng, wp) - 0.5_wp)*grid%dy
            do i = 1, grid%nx_total
               expected_f = abs(F0 + BETA*(y - YREF))
               call check(error, ws%f_centre(i, j) == expected_f, &
                          "beta-plane f_centre not bit-identical to legacy set_f_centre")
               if (allocated(error)) exit checks
            end do
         end do

         ic = grid%nghost + NX/2
         jc = grid%nghost + NY/2
         rel = abs(ws%beta_centre(ic, jc) - BETA)/BETA
         call check(error, rel < 1.0e-12_wp, &
                    "beta_centre != namelist beta on a linear beta-plane field")

         call ws%destroy()
         call metrics%destroy()
         deallocate (f_centre, f_corner)
      end block checks
   end subroutine test_beta_plane_bit_identical

end module test_ocean_wave_speed
