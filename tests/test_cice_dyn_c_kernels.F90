!! Coupling kernels of the C-grid EVP dynamics under the SIS2-derived ice state
!! (CICE6 port D9.2): `rdb_ice_dyn_c` (sea-surface height, the forcing at E and
!! N, the face-stress post-processing, `strocnT`) and the face blend of
!! `rdb_ice_ocean_coupler` (`ice_ocean_stress_face_impl`).
!!
!! All cases run on a NON-UNIFORM lon-lat sector (dx varies with latitude, the
!! areas differ face to face), so a kernel that secretly uses a uniform spacing
!! or `dx*dy` fails here.  Every array a kernel touches is device-mapped
!! (`mem:separate`) and results are pulled back with `update self` on the
!! component arrays.
!!
!! Cases:
!!   * `tilt_closed_form`        the sea-surface slope of a field that is linear
!!                               in the along-face distance is its coefficient
!!                               (x from E, y from N), and a field that depends
!!                               on the other axis only has zero slope; closed
!!                               faces carry zero slope.
!!   * `ocean_velocity_layers`   uniform per-layer ocean velocity: the face
!!                               takes the first layer live in BOTH abutting
!!                               columns (a vanished top layer is skipped), the
!!                               along component and the tangential mean.
!!   * `coriolis_wind_weight`    Coriolis mean of the two face end points on a
!!                               field linear in latitude; the wind stress on
!!                               the ice is the snapshot times the face mean
!!                               concentration, equal to the preparation's `aiE`.
!!   * `inputs_reference`        every output on every owned cell against an
!!                               independent host transcription, on a random
!!                               state with vanished layers.
!!   * `post_and_strocnT`        `strocn/ai` where `ai > puny`, the mask
!!                               mirrors, and `strocnT`: closed form for a
!!                               uniform face stress, flux conservation
!!                               `sum(A_T strocnT) = sum(A_E w_E)` for a random
!!                               one with the wall faces empty.
!!   * `stress_handoff_exact`    with no wind the ocean stress is exactly
!!                               `-strocn` and `sum(area tau) = -sum(area
!!                               strocn)` bitwise; with wind the ice part
!!                               matches to rounding; faces outside the ice get
!!                               the pure wind.
module test_cice_dyn_c_kernels
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init
   use rdb_ice_cgrid_metrics, only: ice_cgrid_metrics_t
   use rdb_ice_dyn_c, only: dync_inputs_impl, dync_post_impl, dync_strocnt_impl
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_face_impl
   implicit none
   private

   public :: collect_cice_dyn_c_kernels_tests

   integer, parameter :: NG = 3
   integer, parameter :: NXP = 10, NYP = 8, NZ = 3
   real(wp), parameter :: DLON = 2.0_wp
   real(wp), parameter :: DLAT = 1.5_wp
   real(wp), parameter :: REARTH = 6.378e6_wp
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: TOL = 1.0e-12_wp   ! FMA on host code too
#else
   real(wp), parameter :: TOL = 1.0e-13_wp
#endif

   type :: fx_t
      integer :: nx = 0, ny = 0, is = 0, ie = 0, js = 0, je = 0
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m
      type(ocean_bc_state_t) :: bc
      type(ice_cgrid_metrics_t) :: cm
      real(wp), allocatable :: h_layer(:, :, :), u_x(:, :, :), v_y(:, :, :)
      integer, allocatable :: k_top_u(:, :), k_top_v(:, :)
      real(wp), allocatable :: f_corner(:, :), ssh(:, :), aice(:, :), tau_a_x(:, :), tau_a_y(:, :)
      real(wp), allocatable :: fcorE(:, :), fcorN(:, :), strairxE(:, :), strairyN(:, :)
      real(wp), allocatable :: ss_tltxE(:, :), ss_tltyN(:, :), uocnE(:, :), vocnE(:, :), uocnN(:, :), vocnN(:, :)
   end type fx_t

contains

   subroutine collect_cice_dyn_c_kernels_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("tilt_closed_form", test_tilt), &
                  new_unittest("ocean_velocity_layers", test_layers), &
                  new_unittest("coriolis_wind_weight", test_cor_wind), &
                  new_unittest("inputs_reference", test_reference), &
                  new_unittest("post_and_strocnT", test_post), &
                  new_unittest("stress_handoff_exact", test_handoff) &
                  ]
   end subroutine collect_cice_dyn_c_kernels_tests

   ! ------------------------------------------------------------------
   ! fixture
   ! ------------------------------------------------------------------

   subroutine fx_build(fx)
      type(fx_t), intent(inout) :: fx
      integer :: nx, ny

      call fx%grid%init(NXP, NYP, NG, DLON, DLAT)
      call make_spherical_metrics(fx%m, fx%grid, 0.0_wp, 40.0_wp, DLON, DLAT, REARTH)
      call ocean_bc_state_init(fx%bc, fx%grid, nz_ml=0, n_tracers=0)
      call fx%cm%init(fx%grid, fx%m, fx%bc)
      nx = fx%grid%nx_total
      ny = fx%grid%ny_total
      fx%nx = nx
      fx%ny = ny
      fx%is = NG + 1
      fx%ie = NG + NXP
      fx%js = NG + 1
      fx%je = NG + NYP
      allocate (fx%h_layer(nx, ny, NZ), source=10.0_wp)
      allocate (fx%u_x(nx + 1, ny, NZ), source=0.0_wp)
      allocate (fx%v_y(nx, ny + 1, NZ), source=0.0_wp)
      allocate (fx%k_top_u(nx + 1, ny), source=NZ)
      allocate (fx%k_top_v(nx, ny + 1), source=NZ)
      allocate (fx%f_corner(nx + 1, ny + 1), source=1.0e-4_wp)
      allocate (fx%ssh(nx, ny), source=0.0_wp)
      allocate (fx%aice(nx, ny), source=0.0_wp)
      allocate (fx%tau_a_x(nx + 1, ny), source=0.0_wp)
      allocate (fx%tau_a_y(nx, ny + 1), source=0.0_wp)
      allocate (fx%fcorE(nx, ny), fx%fcorN(nx, ny), fx%strairxE(nx, ny), fx%strairyN(nx, ny), source=0.0_wp)
      allocate (fx%ss_tltxE(nx, ny), fx%ss_tltyN(nx, ny), source=0.0_wp)
      allocate (fx%uocnE(nx, ny), fx%vocnE(nx, ny), fx%uocnN(nx, ny), fx%vocnN(nx, ny), source=0.0_wp)
   end subroutine fx_build

   subroutine fx_map(fx)
      type(fx_t), intent(inout) :: fx
      call fx%cm%enter_data()
      !$acc enter data copyin(fx%h_layer, fx%u_x, fx%v_y, fx%k_top_u, fx%k_top_v, fx%f_corner, fx%ssh)
      !$acc enter data copyin(fx%aice, fx%tau_a_x, fx%tau_a_y)
      !$acc enter data copyin(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc enter data copyin(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN)
   end subroutine fx_map

   subroutine fx_push(fx)
      type(fx_t), intent(inout) :: fx
      !$acc update device(fx%h_layer, fx%u_x, fx%v_y, fx%k_top_u, fx%k_top_v, fx%f_corner, fx%ssh)
      !$acc update device(fx%aice, fx%tau_a_x, fx%tau_a_y)
   end subroutine fx_push

   subroutine fx_run(fx)
      !! Push the inputs, run the kernel, pull the outputs.
      type(fx_t), intent(inout) :: fx
      call fx_push(fx)
      call dync_inputs_impl(fx%nx, fx%ny, NZ, fx%is, fx%ie, fx%js, fx%je, fx%h_layer, fx%u_x, fx%v_y, &
                            fx%k_top_u, fx%k_top_v, fx%f_corner, fx%ssh, fx%m%idxCu, fx%m%idyCv, &
                            fx%cm%hm, fx%cm%tarea, fx%cm%earea, fx%cm%narea, fx%cm%epm, fx%cm%npm, &
                            fx%aice, fx%tau_a_x, fx%tau_a_y, fx%fcorE, fx%fcorN, fx%strairxE, &
                            fx%strairyN, fx%ss_tltxE, fx%ss_tltyN, fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN)
      !$acc update self(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc update self(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN)
   end subroutine fx_run

   subroutine fx_destroy(fx)
      type(fx_t), intent(inout) :: fx
      !$acc exit data delete(fx%h_layer, fx%u_x, fx%v_y, fx%k_top_u, fx%k_top_v, fx%f_corner, fx%ssh)
      !$acc exit data delete(fx%aice, fx%tau_a_x, fx%tau_a_y)
      !$acc exit data delete(fx%fcorE, fx%fcorN, fx%strairxE, fx%strairyN, fx%ss_tltxE, fx%ss_tltyN)
      !$acc exit data delete(fx%uocnE, fx%vocnE, fx%uocnN, fx%vocnN)
      call fx%cm%exit_data()
      call fx%cm%destroy()
      call destroy_cartesian_metrics(fx%m)
   end subroutine fx_destroy

   function frac(s) result(r)
      integer(int64), intent(inout) :: s
      real(wp) :: r
      s = lcg_next(s)
      r = lcg_unit(s)
   end function frac

   ! ------------------------------------------------------------------
   ! tests
   ! ------------------------------------------------------------------

   subroutine test_tilt(error)
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp), allocatable :: xc(:, :), yc(:, :)
      real(wp), parameter :: S0 = 3.0e-6_wp, S1 = -2.0e-6_wp
      integer :: i, j
      real(wp) :: worst, ref, scale

      call fx_build(fx)
      call fx_map(fx)
      allocate (xc(fx%nx, fx%ny), yc(fx%nx, fx%ny))
      ! Distance along the row / column between neighbouring T centres is the
      ! face length dxCu(i+1,j) / dyCv(i,j+1).
      do j = 1, fx%ny
         xc(1, j) = 0.0_wp
         do i = 1, fx%nx - 1
            xc(i + 1, j) = xc(i, j) + fx%m%dxCu(i + 1, j)
         end do
      end do
      do i = 1, fx%nx
         yc(i, 1) = 0.0_wp
         do j = 1, fx%ny - 1
            yc(i, j + 1) = yc(i, j) + fx%m%dyCv(i, j + 1)
         end do
      end do
      ! x-linear field: x slope is S0 on every open E face
      fx%ssh = S0*xc
      call fx_run(fx)
      worst = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            ref = S0*fx%cm%epm(i, j)
            if (.not. ieee_is_finite(fx%ss_tltxE(i, j))) then
               call check(error, .false., "non-finite x slope")
               return
            end if
            worst = max(worst, abs(fx%ss_tltxE(i, j) - ref))
         end do
      end do
      call check(error, worst <= 1.0e-12_wp*abs(S0), "x slope of an x-linear field is its coefficient")
      if (allocated(error)) return
      ! y-linear field: y slope is S1; the x slope of a field of j only is zero
      fx%ssh = S1*yc
      call fx_run(fx)
      worst = 0.0_wp
      scale = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            ref = S1*fx%cm%npm(i, j)
            worst = max(worst, abs(fx%ss_tltyN(i, j) - ref))
            scale = max(scale, abs(fx%ss_tltxE(i, j)))
         end do
      end do
      call check(error, worst <= 1.0e-12_wp*abs(S1), "y slope of a y-linear field is its coefficient")
      if (allocated(error)) return
      call check(error, scale <= 1.0e-12_wp*abs(S1), "x slope of a field of j only is zero")
      if (allocated(error)) return
      ! a closed face (the east wall face E(ie)) carries no slope whatever the field
      fx%ssh = 7.0_wp*xc
      call fx_run(fx)
      call check(error, all(fx%ss_tltxE(fx%ie, fx%js:fx%je) == 0.0_wp), "the wall face carries zero slope")
      call fx_destroy(fx)
   end subroutine test_tilt

   subroutine test_layers(error)
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer :: i, j, k, ik
      real(wp) :: uk(NZ), vk(NZ), expect
      logical :: ok

      call fx_build(fx)
      call fx_map(fx)
      uk = [1.0_wp, 2.0_wp, 3.0_wp]       ! bed .. surface
      vk = [-1.0_wp, -2.0_wp, -3.0_wp]
      do k = 1, NZ
         fx%u_x(:, :, k) = uk(k)
         fx%v_y(:, :, k) = vk(k)
      end do
      ! the surface layer vanishes in the columns i = is+2 (all rows)
      fx%h_layer(fx%is + 2, :, NZ) = 0.0_wp
      call fx_run(fx)
      ok = .true.
      do j = fx%js, fx%je
         do i = fx%is, fx%ie - 1
            ! E(i,j) is x-face i+1, between columns i and i+1
            ik = NZ
            if (i == fx%is + 2 .or. i + 1 == fx%is + 2) ik = NZ - 1
            expect = uk(ik)*fx%cm%epm(i, j)
            if (abs(fx%uocnE(i, j) - expect) > TOL) ok = .false.
            ! N(i,j) is y-face j+1 between rows j and j+1: the layers of the columns
            ! are the same in both rows, a vanished column is the column itself
            ik = NZ
            if (i == fx%is + 2) ik = NZ - 1
            expect = vk(ik)*fx%cm%npm(i, j)
            if (abs(fx%vocnN(i, j) - expect) > TOL) ok = .false.
         end do
      end do
      call check(error, ok, "the along-face component takes the first layer live on both sides")
      if (allocated(error)) return
      ! tangential: a face whose four neighbours are all in the same layer reads that layer
      ok = .true.
      do j = fx%js, fx%je - 1
         do i = fx%is, fx%is
            expect = vk(NZ)*fx%cm%epm(i, j)
            if (abs(fx%vocnE(i, j) - expect) > TOL) ok = .false.
         end do
      end do
      call check(error, ok, "the tangential mean of one layer is that layer")
      call fx_destroy(fx)
   end subroutine test_layers

   subroutine test_cor_wind(error)
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp), parameter :: F0 = 1.0e-4_wp, FY = 3.0e-6_wp, A0 = 0.62_wp, TX = 0.07_wp, TY = -0.03_wp
      integer :: i, j
      logical :: ok

      call fx_build(fx)
      call fx_map(fx)
      do j = 1, fx%ny + 1
         fx%f_corner(:, j) = F0 + FY*real(j, wp)
      end do
      fx%aice = A0
      fx%tau_a_x = TX
      fx%tau_a_y = TY
      call fx_run(fx)
      ok = .true.
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            ! E(i,j): corners (i+1,j) and (i+1,j+1): mean row j + 1/2.  N(i,j): corners row j+1.
            if (abs(fx%fcorE(i, j) - (F0 + FY*(real(j, wp) + 0.5_wp))) > TOL*F0) ok = .false.
            if (abs(fx%fcorN(i, j) - (F0 + FY*real(j + 1, wp))) > TOL*F0) ok = .false.
         end do
      end do
      call check(error, ok, "Coriolis is the mean of the two face end points")
      if (allocated(error)) return
      ok = .true.
      do j = fx%js, fx%je
         do i = fx%is, fx%ie - 1
            ! both cells wet and ice-covered: the face concentration is A0
            if (abs(fx%strairxE(i, j) - A0*TX) > TOL) ok = .false.
         end do
      end do
      do j = fx%js, fx%je - 1
         do i = fx%is, fx%ie
            if (abs(fx%strairyN(i, j) - A0*TY) > TOL) ok = .false.
         end do
      end do
      call check(error, ok, "the wind on the ice is the snapshot times the face concentration")
      call fx_destroy(fx)
   end subroutine test_cor_wind

   subroutine test_reference(error)
      !! Every output of `dync_inputs_impl` against an independent host
      !! transcription on a random state (vanished layers included).
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      integer(int64) :: s
      integer :: i, j, k, ka, kb
      real(wp) :: w1, w2, a, e1, e2, e3, e4, ex, ey, worst, wt, r1, r2

      call fx_build(fx)
      call fx_map(fx)
      s = 20260514_int64
      do k = 1, NZ
         do j = 1, fx%ny
            do i = 1, fx%nx
               r1 = frac(s)
               r2 = frac(s)
               fx%h_layer(i, j, k) = merge(0.0_wp, 5.0_wp + 10.0_wp*r2, r1 < 0.15_wp)
            end do
         end do
      end do
      ! the bed layer is always live so a column is never empty
      fx%h_layer(:, :, 1) = 20.0_wp
      do k = 1, NZ
         do j = 1, fx%ny
            do i = 1, fx%nx + 1
               fx%u_x(i, j, k) = frac(s) - 0.5_wp
            end do
         end do
         do j = 1, fx%ny + 1
            do i = 1, fx%nx
               fx%v_y(i, j, k) = frac(s) - 0.5_wp
            end do
         end do
      end do
      do j = 1, fx%ny
         do i = 1, fx%nx
            fx%ssh(i, j) = 0.1_wp*frac(s)
            fx%aice(i, j) = frac(s)
            fx%tau_a_x(i, j) = 0.1_wp*frac(s)
            fx%tau_a_y(i, j) = 0.1_wp*frac(s)
         end do
      end do
      do j = 1, fx%ny + 1
         do i = 1, fx%nx + 1
            fx%f_corner(i, j) = 1.0e-4_wp*(1.0_wp + frac(s))
         end do
      end do
      call fx_run(fx)
      worst = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            ex = 0.5_wp*(fx%f_corner(i + 1, j) + fx%f_corner(i + 1, j + 1))
            worst = max(worst, abs(fx%fcorE(i, j) - ex)/1.0e-4_wp)
            ey = 0.5_wp*(fx%f_corner(i, j + 1) + fx%f_corner(i + 1, j + 1))
            worst = max(worst, abs(fx%fcorN(i, j) - ey)/1.0e-4_wp)
            ex = (fx%ssh(i + 1, j) - fx%ssh(i, j))*fx%m%idxCu(i + 1, j)*fx%cm%epm(i, j)
            worst = max(worst, abs(fx%ss_tltxE(i, j) - ex)/max(abs(ex), 1.0e-12_wp))
            ey = (fx%ssh(i, j + 1) - fx%ssh(i, j))*fx%m%idyCv(i, j + 1)*fx%cm%npm(i, j)
            worst = max(worst, abs(fx%ss_tltyN(i, j) - ey)/max(abs(ey), 1.0e-12_wp))
            ! wind
            w1 = fx%cm%hm(i, j)*fx%cm%tarea(i, j)
            w2 = fx%cm%hm(i + 1, j)*fx%cm%tarea(i + 1, j)
            a = 0.0_wp
            if (w1 + w2 /= 0.0_wp) a = (w1*fx%aice(i, j) + w2*fx%aice(i + 1, j))/(w1 + w2)
            worst = max(worst, abs(fx%strairxE(i, j) - a*fx%tau_a_x(i + 1, j)))
            w2 = fx%cm%hm(i, j + 1)*fx%cm%tarea(i, j + 1)
            a = 0.0_wp
            if (w1 + w2 /= 0.0_wp) a = (w1*fx%aice(i, j) + w2*fx%aice(i, j + 1))/(w1 + w2)
            worst = max(worst, abs(fx%strairyN(i, j) - a*fx%tau_a_y(i, j + 1)))
            ! ocean velocity at E
            ka = ref_ku(fx, i + 1, j)
            worst = max(worst, abs(fx%uocnE(i, j) - fx%u_x(i + 1, j, ka)*fx%cm%epm(i, j)))
            e1 = fx%v_y(i, j, ref_kv(fx, i, j))
            e2 = fx%v_y(i + 1, j, ref_kv(fx, i + 1, j))
            e3 = fx%v_y(i, j + 1, ref_kv(fx, i, j + 1))
            e4 = fx%v_y(i + 1, j + 1, ref_kv(fx, i + 1, j + 1))
            wt = fx%cm%narea(i, j - 1) + fx%cm%narea(i + 1, j - 1) + fx%cm%narea(i, j) + fx%cm%narea(i + 1, j)
            ex = 0.0_wp
            if (wt /= 0.0_wp) ex = (e1*fx%cm%narea(i, j - 1) + e2*fx%cm%narea(i + 1, j - 1) &
                                    + e3*fx%cm%narea(i, j) + e4*fx%cm%narea(i + 1, j))/wt
            worst = max(worst, abs(fx%vocnE(i, j) - ex*fx%cm%epm(i, j)))
            ! ocean velocity at N
            kb = ref_kv(fx, i, j + 1)
            worst = max(worst, abs(fx%vocnN(i, j) - fx%v_y(i, j + 1, kb)*fx%cm%npm(i, j)))
            e1 = fx%u_x(i, j, ref_ku(fx, i, j))
            e2 = fx%u_x(i + 1, j, ref_ku(fx, i + 1, j))
            e3 = fx%u_x(i, j + 1, ref_ku(fx, i, j + 1))
            e4 = fx%u_x(i + 1, j + 1, ref_ku(fx, i + 1, j + 1))
            wt = fx%cm%earea(i - 1, j) + fx%cm%earea(i, j) + fx%cm%earea(i - 1, j + 1) + fx%cm%earea(i, j + 1)
            ey = 0.0_wp
            if (wt /= 0.0_wp) ey = (e1*fx%cm%earea(i - 1, j) + e2*fx%cm%earea(i, j) &
                                    + e3*fx%cm%earea(i - 1, j + 1) + e4*fx%cm%earea(i, j + 1))/wt
            worst = max(worst, abs(fx%uocnN(i, j) - ey*fx%cm%npm(i, j)))
         end do
      end do
      if (.not. ieee_is_finite(worst)) then
         call check(error, .false., "non-finite difference")
         return
      end if
      call check(error, worst <= 50.0_wp*TOL, "dync_inputs_impl differs from the host transcription")
      call fx_destroy(fx)
   contains
      integer function ref_ku(f, ii, jj)
         !! First layer live in both columns of x-face (ii,jj), from the surface down.
         type(fx_t), intent(in) :: f
         integer, intent(in) :: ii, jj
         integer :: kk
         ref_ku = NZ
         do kk = NZ, 1, -1
            if (f%h_layer(ii - 1, jj, kk) > H_VANISHED .and. f%h_layer(ii, jj, kk) > H_VANISHED) then
               ref_ku = kk
               return
            end if
         end do
      end function ref_ku
      integer function ref_kv(f, ii, jj)
         type(fx_t), intent(in) :: f
         integer, intent(in) :: ii, jj
         integer :: kk
         ref_kv = NZ
         do kk = NZ, 1, -1
            if (f%h_layer(ii, jj - 1, kk) > H_VANISHED .and. f%h_layer(ii, jj, kk) > H_VANISHED) then
               ref_kv = kk
               return
            end if
         end do
      end function ref_kv
   end subroutine test_reference

   subroutine test_post(error)
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp), allocatable :: sx(:, :), sy(:, :), aiE(:, :), aiN(:, :), wE(:, :), wN(:, :), er(:, :), nr(:, :)
      real(wp), allocatable :: tx(:, :), ty(:, :)
      logical, allocatable :: le(:, :), ln(:, :)
      integer(int64) :: s
      integer :: i, j, nx, ny
      real(wp) :: puny, lhs, rhs, worst, uw, r1, r2, rsum
      logical :: ok

      call fx_build(fx)
      call fx_map(fx)
      nx = fx%nx
      ny = fx%ny
      puny = 1.0e-11_wp
      allocate (sx(nx, ny), sy(nx, ny), aiE(nx, ny), aiN(nx, ny), wE(nx, ny), wN(nx, ny), er(nx, ny), nr(nx, ny))
      allocate (tx(nx, ny), ty(nx, ny), le(nx, ny), ln(nx, ny))
      s = 77_int64
      sx = 0.0_wp
      sy = 0.0_wp
      aiE = 0.0_wp
      aiN = 0.0_wp
      le = .false.
      ln = .false.
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            aiE(i, j) = frac(s)
            aiN(i, j) = frac(s)
            r1 = frac(s)
            r2 = frac(s)
            if (r1 < 0.2_wp) aiE(i, j) = 0.0_wp
            if (r2 < 0.2_wp) aiN(i, j) = 0.0_wp
            le(i, j) = aiE(i, j) > puny
            ln(i, j) = aiN(i, j) > puny
            sx(i, j) = aiE(i, j)*(frac(s) - 0.5_wp)*0.2_wp
            sy(i, j) = aiN(i, j)*(frac(s) - 0.5_wp)*0.2_wp
         end do
      end do
      ! the wall faces carry no stress (no ice on a closed face)
      aiE(fx%ie, :) = 0.0_wp
      sx(fx%ie, :) = 0.0_wp
      le(fx%ie, :) = .false.
      aiN(:, fx%je) = 0.0_wp
      sy(:, fx%je) = 0.0_wp
      ln(:, fx%je) = .false.
      !$acc enter data copyin(sx, sy, aiE, aiN, le, ln, wE, wN, er, nr, tx, ty)
      call dync_post_impl(nx, ny, fx%is, fx%ie, fx%js, fx%je, puny, sx, sy, aiE, aiN, le, ln, wE, wN, er, nr)
      !$acc update self(wE, wN, er, nr)
      ok = .true.
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            if (aiE(i, j) > puny) then
               if (abs(wE(i, j) - sx(i, j)/aiE(i, j)) > TOL) ok = .false.
            else
               if (wE(i, j) /= 0.0_wp) ok = .false.
            end if
            if (aiN(i, j) > puny) then
               if (abs(wN(i, j) - sy(i, j)/aiN(i, j)) > TOL) ok = .false.
            else
               if (wN(i, j) /= 0.0_wp) ok = .false.
            end if
            if ((er(i, j) > 0.5_wp) .neqv. le(i, j)) ok = .false.
            if ((nr(i, j) > 0.5_wp) .neqv. ln(i, j)) ok = .false.
         end do
      end do
      call check(error, ok, "strocn/ai, its zero below puny and the 0/1 mask mirrors")
      if (allocated(error)) return

      ! strocnT of a uniform face stress per ice area: every T cell between two
      ! covered faces takes that stress when the two face areas equal the cell's
      uw = 0.37_wp
      wE = 0.0_wp
      wN = 0.0_wp
      wE(fx%is:fx%ie - 1, fx%js:fx%je) = uw
      wN(fx%is:fx%ie, fx%js:fx%je - 1) = -uw
      !$acc update device(wE, wN)
      call dync_strocnt_impl(nx, ny, fx%is, fx%ie, fx%js, fx%je, wE, wN, fx%cm%earea, fx%cm%narea, &
                             fx%cm%tarear, tx, ty)
      !$acc update self(tx, ty)
      worst = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is + 1, fx%ie - 1
            ! interior cells: both E faces are covered; the mean of A_E over A_T is not 1 on a
            ! lon-lat grid, so compare with the closed form of the map itself
            rhs = 0.5_wp*uw*(fx%cm%earea(i - 1, j) + fx%cm%earea(i, j))*fx%cm%tarear(i, j)
            worst = max(worst, abs(tx(i, j) - rhs)/uw)
         end do
      end do
      call check(error, worst <= TOL, "strocnT_x of a uniform face stress is the area-weighted two-face mean")
      if (allocated(error)) return

      ! flux conservation of the E -> T map with the wall faces empty
      wE = 0.0_wp
      wN = 0.0_wp
      s = 5_int64
      do j = fx%js, fx%je
         do i = fx%is, fx%ie - 1
            wE(i, j) = frac(s) - 0.5_wp
         end do
      end do
      do j = fx%js, fx%je - 1
         do i = fx%is, fx%ie
            wN(i, j) = frac(s) - 0.5_wp
         end do
      end do
      !$acc update device(wE, wN)
      call dync_strocnt_impl(nx, ny, fx%is, fx%ie, fx%js, fx%je, wE, wN, fx%cm%earea, fx%cm%narea, &
                             fx%cm%tarear, tx, ty)
      !$acc update self(tx, ty)
      lhs = 0.0_wp
      rhs = 0.0_wp
      rsum = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            lhs = lhs + fx%cm%tarea(i, j)*tx(i, j)
            rhs = rhs + wE(i, j)*fx%cm%earea(i, j)
            rsum = rsum + abs(wE(i, j))*fx%cm%earea(i, j)
         end do
      end do
      call check(error, abs(lhs - rhs) <= 1.0e-13_wp*rsum, &
                 "sum(A_T strocnT_x) equals sum(A_E w_E) (flux-conserving map)")
      !$acc exit data delete(sx, sy, aiE, aiN, le, ln, wE, wN, er, nr, tx, ty)
      call fx_destroy(fx)
   end subroutine test_post

   subroutine test_handoff(error)
      type(error_type), allocatable, intent(out) :: error
      type(fx_t) :: fx
      real(wp), allocatable :: tau_x(:, :), tau_y(:, :), ta_x(:, :), ta_y(:, :), sx(:, :), sy(:, :)
      real(wp), allocatable :: aiE(:, :), aiN(:, :), zx(:, :), zy(:, :), wu(:, :), wv(:, :)
      integer(int64) :: s
      integer :: i, j, nx, ny
      real(wp) :: lhs, rhs, worst
      logical :: ok

      call fx_build(fx)
      nx = fx%nx
      ny = fx%ny
      allocate (tau_x(nx + 1, ny), tau_y(nx, ny + 1), ta_x(nx + 1, ny), ta_y(nx, ny + 1))
      allocate (sx(nx, ny), sy(nx, ny), aiE(nx, ny), aiN(nx, ny), zx(nx + 1, ny), zy(nx, ny + 1))
      allocate (wu(nx + 1, ny), wv(nx, ny + 1))
      wu = 1.0_wp
      wv = 1.0_wp
      s = 99_int64
      aiE = 0.0_wp
      aiN = 0.0_wp
      sx = 0.0_wp
      sy = 0.0_wp
      do j = fx%js, fx%je
         do i = fx%is, fx%ie
            aiE(i, j) = frac(s)
            aiN(i, j) = frac(s)
            sx(i, j) = aiE(i, j)*(frac(s) - 0.5_wp)
            sy(i, j) = aiN(i, j)*(frac(s) - 0.5_wp)
         end do
      end do
      zx = 0.0_wp
      zy = 0.0_wp
      do j = 1, ny
         do i = 1, nx + 1
            ta_x(i, j) = 0.05_wp*frac(s)
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            ta_y(i, j) = 0.05_wp*frac(s)
         end do
      end do
      !$acc enter data copyin(tau_x, tau_y, ta_x, ta_y, zx, zy, sx, sy, aiE, aiN, wu, wv)
      ! ---- no wind: the ocean receives exactly -strocn, bitwise ----
      tau_x = 99.0_wp
      tau_y = 99.0_wp
      !$acc update device(tau_x, tau_y)
      call ice_ocean_stress_face_impl(nx, ny, tau_x, tau_y, zx, zy, sx, sy, aiE, aiN, wu, wv)
      !$acc update self(tau_x, tau_y)
      lhs = 0.0_wp
      rhs = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            lhs = lhs + fx%cm%earea(i, j)*tau_x(i + 1, j)
            rhs = rhs - fx%cm%earea(i, j)*sx(i, j)
         end do
      end do
      call check(error, lhs == rhs, "sum(area tau_x) = -sum(area strocnxE) holds bitwise without wind")
      if (allocated(error)) return
      lhs = 0.0_wp
      rhs = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            lhs = lhs + fx%cm%narea(i, j)*tau_y(i, j + 1)
            rhs = rhs - fx%cm%narea(i, j)*sy(i, j)
         end do
      end do
      call check(error, lhs == rhs, "sum(area tau_y) = -sum(area strocnyN) holds bitwise without wind")
      if (allocated(error)) return
      ! ---- with wind: the ice part is -strocn, the rest the unintercepted wind ----
      !$acc update device(tau_x, tau_y)
      call ice_ocean_stress_face_impl(nx, ny, tau_x, tau_y, ta_x, ta_y, sx, sy, aiE, aiN, wu, wv)
      !$acc update self(tau_x, tau_y)
      worst = 0.0_wp
      ok = .true.
      do j = 1, ny
         do i = 1, nx
            worst = max(worst, abs((tau_x(i + 1, j) - (1.0_wp - aiE(i, j))*ta_x(i + 1, j)) + sx(i, j)))
            worst = max(worst, abs((tau_y(i, j + 1) - (1.0_wp - aiN(i, j))*ta_y(i, j + 1)) + sy(i, j)))
            if (aiE(i, j) == 0.0_wp .and. tau_x(i + 1, j) /= ta_x(i + 1, j)) ok = .false.
            if (aiN(i, j) == 0.0_wp .and. tau_y(i, j + 1) /= ta_y(i, j + 1)) ok = .false.
         end do
      end do
      do j = 1, ny
         if (tau_x(1, j) /= ta_x(1, j)) ok = .false.
      end do
      do i = 1, nx
         if (tau_y(i, 1) /= ta_y(i, 1)) ok = .false.
      end do
      call check(error, worst <= 1.0e-15_wp, "ice part of the ocean stress is -strocn to rounding")
      if (allocated(error)) return
      call check(error, ok, "faces without ice and the first faces take the pure wind")
      ! ---- closed faces carry no stress (full ice cover, full wind snapshot) ----
      wu(fx%is + 2, :) = 0.0_wp
      wv(:, fx%js + 2) = 0.0_wp
      aiE = 1.0_wp
      aiN = 1.0_wp
      !$acc update device(wu, wv, aiE, aiN)
      call ice_ocean_stress_face_impl(nx, ny, tau_x, tau_y, ta_x, ta_y, sx, sy, aiE, aiN, wu, wv)
      !$acc update self(tau_x, tau_y)
      ok = all(tau_x(fx%is + 2, :) == 0.0_wp) .and. all(tau_y(:, fx%js + 2) == 0.0_wp)
      call check(error, ok, "closed faces carry no stress")
      if (allocated(error)) return
      call check(error, tau_x(fx%is + 3, fx%js) == -sx(fx%is + 2, fx%js) .and. &
                 tau_y(fx%is, fx%js + 3) == -sy(fx%is, fx%js + 2), "open faces are unmasked")
      !$acc exit data delete(tau_x, tau_y, ta_x, ta_y, zx, zy, sx, sy, aiE, aiN, wu, wv)
      call fx_destroy_partial(fx)
   end subroutine test_handoff

   subroutine fx_destroy_partial(fx)
      !! `test_handoff` never maps the `fx` arrays; only the metric set.
      type(fx_t), intent(inout) :: fx
      call fx%cm%destroy()
      call destroy_cartesian_metrics(fx%m)
   end subroutine fx_destroy_partial

end module test_cice_dyn_c_kernels
