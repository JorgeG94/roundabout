!! The isopycnal-slope bed datum: `ocean_slopes_build_e` builds GEOPOTENTIAL
!! interface heights (bed at `z = −D`), so the interface-tilt rotation term
!! `−∂zρ·(e_W − e_E)` sees the real tilt of an interface and never the
!! bathymetry step.  Every case runs the production device kernels
!! (`ocean_slopes_compute` [+ `gm_compute_transports`]) over a bed that
!! steps down eastward AND northward (30 m per column, 20 m per row, so
!! both the u- and the v-face pass see a step), on three coordinates:
!!
!!   * `sigma`  — `h_k = (D + η)/nz` at `η = 0` (the production
!!                `compute_target_h` SIGMA/ZSTAR branch, uniform `dsig`);
!!   * `zstar`  — the same branch at a uniformly RAISED free surface
!!                `η = ETA_Z`: the column top sits at `+η`, not at the
!!                datum, and must still read flat;
!!   * `z_fixed` — the production `ocean_vcoord_z_fixed_target_uniform`
!!                staircase (partial bottom cells + bed fillers) with the
!!                production closed-face masks (`zfixed_closed_faces`).
!!
!! Bottom-up (k = 1 bed, k = nz top).  Linear EOS, T-only stratification,
!! `kd_smooth = 0` (the vert-fill regulariser perturbs the boundary rows
!! of a linear profile at O(κΔt/Δz²); it has its own test).
!!
!!  1. `flat_over_step_<coord>` — a horizontally uniform stratification at
!!     rest: slope ≡ 0 to round-off at every interface and GM moves
!!     nothing.  On sigma / zstar, T is linear in the GEOPOTENTIAL height
!!     of each layer centre (the layers themselves tilt with the bed); on
!!     z_fixed it is layer-uniform (a z-level stratification).  Fails on
!!     the zero-bed-datum slopes for sigma and zstar (measured max|S| =
!!     1.5e-3 — the 30 m / 20 km step read as a slope — and max|uhD| =
!!     3.0e4 m³/s; z_fixed the same unless the tilt term is dropped).
!!  2. `tilted_over_step_<coord>` — `T = T0 + GZ·z + GX·x + GY·y` over the
!!     same step: every interior u-/v-face reads the analytic bounded
!!     neutral slope `(−G_h/GZ)/√(1 + (G_h/GZ)²)`.  On sigma / zstar that
!!     is exact (uniform layers in each column make the along-layer
!!     weights equal, and the tilt term cancels the layers' own tilt
!!     exactly); on z_fixed it is checked at interfaces whose four cells
!!     are full nominal cells (a partial cell's centre is not at its
!!     nominal depth, which is truncation, not the datum).  The zero-datum
!!     answer misses by the step itself, 1.5e-3 (x) / 1.0e-3 (y).
module test_ocean_slopes_datum
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_closed_face_masks
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_slopes_datum_tests

   integer, parameter :: COORD_SIGMA = 1, COORD_ZSTAR = 2, COORD_ZFIXED = 3
   integer, parameter :: NG = 2
   integer, parameter :: NXP = 8, NYP = 4, NZ = 8
   real(wp), parameter :: DX = 20000.0_wp
   real(wp), parameter :: H_NOM = 50.0_wp
      !! Nominal z_fixed spacing (m); NZ*H_NOM = 400 m = the deepest bed.
   real(wp), parameter :: H_MIN = 1.0e-4_wp
      !! `zstar_h_min` at its default: the inert z_fixed filler thickness.
   real(wp), parameter :: ETA_Z = 0.7_wp
      !! Uniform free-surface height of the zstar case (m).
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: KHTH = 1000.0_wp
   real(wp), parameter :: GZ = 0.02_wp
      !! Vertical T gradient (K/m), warm above: stable.
   real(wp), parameter :: GX = 2.0e-6_wp, GY = -1.0e-6_wp
      !! Horizontal T gradients of case 2 (K/m): slopes ~1e-4, 5e-5.

contains

   subroutine collect_ocean_slopes_datum_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("flat_over_step_sigma", test_flat_sigma), &
                  new_unittest("flat_over_step_zstar", test_flat_zstar), &
                  new_unittest("flat_over_step_zfixed", test_flat_zfixed), &
                  new_unittest("tilted_over_step_sigma", test_tilted_sigma), &
                  new_unittest("tilted_over_step_zstar", test_tilted_zstar), &
                  new_unittest("tilted_over_step_zfixed", test_tilted_zfixed) &
                  ]
   end subroutine collect_ocean_slopes_datum_tests

   ! ------------------------------------------------------------------
   ! Setup
   ! ------------------------------------------------------------------

   pure real(wp) function bed_depth(ip, jp) result(d)
      !! The step: 400 m at the south-west corner, 30 m shallower per
      !! column eastward and 20 m per row northward (off the 50 m z_fixed
      !! spacing, so z_fixed carries partial cells and 0..3 bed fillers).
      integer, intent(in) :: ip, jp
      d = real(NZ, wp)*H_NOM - 30.0_wp*real(ip - 1, wp) - 20.0_wp*real(jp - 1, wp)
   end function bed_depth

   subroutine build_case(coord, gx, gy, grid, metrics, ms, sl, gm, eos, full)
      !! Grid + metrics + a resting state over the step on `coord`, with
      !! `T = T0 + GZ·z + gx·x + gy·y`.  Ghost columns replicate the
      !! nearest physical column (bed AND state), as a wall halo does.
      !! `full(i,j,k)` flags the cells that are full nominal layers
      !! (z_fixed) — every live cell on sigma / zstar.  The slopes slot's
      !! bed datum is set from the SAME `D` the thicknesses were built on.
      integer, intent(in) :: coord
      real(wp), intent(in) :: gx, gy
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(eos_t), intent(out) :: eos
      logical, allocatable, intent(out) :: full(:, :, :)
      real(wp), allocatable :: d(:, :), eta0(:, :), z_top(:, :), tgt(:, :, :)
      real(wp) :: eta, e_bot, zc, t_val, xc, yc
      integer :: i, j, k, ni, nj, ip, jp

      call grid%init(NXP, NYP, NG, DX, DX)
      ni = grid%nx_total
      nj = grid%ny_total

      allocate (d(ni, nj), eta0(ni, nj), z_top(ni, nj))
      allocate (tgt(ni, nj, NZ), source=0.0_wp)
      allocate (full(ni, nj, NZ), source=.false.)
      do j = 1, nj
         jp = min(max(j - NG, 1), NYP)
         do i = 1, ni
            ip = min(max(i - NG, 1), NXP)
            d(i, j) = bed_depth(ip, jp)
         end do
      end do
      eta = 0.0_wp
      if (coord == COORD_ZSTAR) eta = ETA_Z
      eta0 = 0.0_wp
      z_top = 0.0_wp

      select case (coord)
      case (COORD_SIGMA, COORD_ZSTAR)
         call make_cartesian_metrics(metrics, grid)
         do k = 1, NZ
            tgt(:, :, k) = (d + eta)/real(NZ, wp)
         end do
         full = .true.
      case (COORD_ZFIXED)
         ! `nz_closed` grows open_u/open_v BEFORE the device map.
         call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
         call ocean_vcoord_z_fixed_target_uniform(tgt, d, eta0, z_top, &
                                                  ni, nj, NZ, H_NOM, H_MIN)
         call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                             tgt, ni, nj, NZ, H_VANISHED)
         metrics%use_closed_faces = .true.
         !$acc update self(metrics%open_u, metrics%open_v)
         full = abs(tgt - H_NOM) <= 1.0e-9_wp*H_NOM
      end select

      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = tgt
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do j = 1, nj
         jp = min(max(j - NG, 1), NYP)
         yc = real(jp, wp)*DX
         do i = 1, ni
            ip = min(max(i - NG, 1), NXP)
            xc = real(ip, wp)*DX
            e_bot = -d(i, j)
            do k = 1, NZ
               if (coord == COORD_ZFIXED) then
                  ! z-level: the NOMINAL centre height, so a partial or
                  ! filler cell carries its level's value (layer-uniform).
                  zc = -(real(NZ - k, wp) + 0.5_wp)*H_NOM
               else
                  ! Terrain-following: the GEOPOTENTIAL centre height.
                  zc = e_bot + 0.5_wp*tgt(i, j, k)
                  e_bot = e_bot + tgt(i, j, k)
               end if
               t_val = T0 + GZ*zc + gx*xc + gy*yc
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*tgt(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*tgt(i, j, k)
            end do
         end do
      end do
      if (coord == COORD_ZFIXED) then
         ! Each bed filler takes its donor's T before the I1' sweep, so the
         ! sweep pools equal concentrations and leaves the donor's T
         ! unshifted (as in test_ocean_gm_zfixed).
         do j = 1, nj
            do i = 1, ni
               t_val = -huge(1.0_wp)
               do k = NZ, 1, -1
                  if (tgt(i, j, k) > H_VANISHED) then
                     t_val = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/tgt(i, j, k)
                  else if (t_val > -huge(1.0_wp)) then
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*tgt(i, j, k)
                  end if
               end do
            end do
         end do
         call ms%enforce_vanished_content_host(ni, nj)
      end if

      call sl%init(grid, nz_ml=NZ)
      call sl%set_bathymetry(d)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%kd_smooth = 0.0_wp
      sl%min_dz_for_n2 = 1.0_wp

      call gm%init(grid, nz_ml=NZ)
      gm%enable = .true.
      gm%khth = KHTH
      gm%khth_slope_max = 0.01_wp
      gm%rho0 = RHO0

      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine build_case

   subroutine run_kernels(grid, metrics, ms, sl, gm, eos)
      !! Map, run slopes + GM on the device, pull the outputs back, unmap.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      type(eos_t), intent(in) :: eos
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()

      call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT)

      !$acc update self(gm%uhD, gm%vhD, gm%gm_src)
      !$acc update self(sl%slope_x, sl%slope_y, sl%n2_u, sl%n2_v)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_kernels

   subroutine teardown(metrics, ms, sl, gm)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_gm_t), intent(inout) :: gm
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown

   ! ------------------------------------------------------------------
   ! Case 1: flat stratification at rest over the step
   ! ------------------------------------------------------------------

   subroutine check_flat(error, coord, label)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: coord
      character(len=*), intent(in) :: label
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      logical, allocatable :: full(:, :, :)
      real(wp) :: smax, fmax, n2min
      character(len=128) :: msg
      integer :: k

      call build_case(coord, 0.0_wp, 0.0_wp, grid, metrics, ms, sl, gm, eos, full)
      checks: block
         call run_kernels(grid, metrics, ms, sl, gm, eos)

         smax = max(maxval(abs(sl%slope_x)), maxval(abs(sl%slope_y)))
         fmax = max(maxval(abs(gm%uhD)), maxval(abs(gm%vhD)))
         write (msg, '(a,es10.3,a,es10.3)') "max|S| = ", smax, &
            ", max|uhD,vhD| = ", fmax
         ! Non-vacuous: the stratification must be live at the interior
         ! u-faces (N² > 0 somewhere strictly inside), or "S = 0" is the
         ! mag2 = 0 branch, not the datum.
         n2min = huge(1.0_wp)
         do k = 2, NZ
            n2min = min(n2min, maxval(sl%n2_u(NG + 2:NG + NXP, NG + 1:NG + NYP, k)))
         end do
         call check(error, n2min > 0.0_wp, label//": the interior N^2 must be live")
         if (allocated(error)) exit checks
         ! Round-off scale: the tilt term is drdz·(e_W − e_E) with e ~ 400 m,
         ! so a cancellation error ~ 1e-13 m over DX = 2e4 m reads as
         ! S ~ 1e-17.  The zero-datum answer was ~1e-3 (sigma / zstar).
         call check(error, smax <= 1.0e-14_wp, &
                    label//": a flat stratification over a bathymetry step has "// &
                    "ZERO isopycnal slope: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, fmax <= 1.0e-12_wp*KHTH*DX, &
                    label//": GM must move nothing at rest: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, all(ieee_is_finite(sl%n2_u)) .and. &
                    all(ieee_is_finite(sl%n2_v)), label//": N^2 finite everywhere")
      end block checks
      call teardown(metrics, ms, sl, gm)
   end subroutine check_flat

   subroutine test_flat_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      call check_flat(error, COORD_SIGMA, "sigma")
   end subroutine test_flat_sigma

   subroutine test_flat_zstar(error)
      type(error_type), allocatable, intent(out) :: error
      call check_flat(error, COORD_ZSTAR, "zstar")
   end subroutine test_flat_zstar

   subroutine test_flat_zfixed(error)
      type(error_type), allocatable, intent(out) :: error
      call check_flat(error, COORD_ZFIXED, "z_fixed")
   end subroutine test_flat_zfixed

   ! ------------------------------------------------------------------
   ! Case 2: uniformly tilted isopycnals over the step
   ! ------------------------------------------------------------------

   subroutine check_tilted(error, coord, label)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: coord
      character(len=*), intent(in) :: label
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      logical, allocatable :: full(:, :, :)
      real(wp) :: sx_ref, sy_ref, ex, ey
      integer :: i, j, k, nx_chk, ny_chk
      character(len=128) :: msg

      call build_case(coord, GX, GY, grid, metrics, ms, sl, gm, eos, full)
      checks: block
         call run_kernels(grid, metrics, ms, sl, gm, eos)

         sx_ref = (-GX/GZ)/sqrt(1.0_wp + (GX/GZ)**2)
         sy_ref = (-GY/GZ)/sqrt(1.0_wp + (GY/GZ)**2)
         ex = 0.0_wp
         ey = 0.0_wp
         nx_chk = 0
         ny_chk = 0
         ! Interior u-faces between two PHYSICAL columns, interior
         ! interfaces whose four cells are full layers.
         do k = 2, NZ
            do j = NG + 1, NG + NYP
               do i = NG + 2, NG + NXP
                  if (.not. (full(i - 1, j, k) .and. full(i, j, k) .and. &
                             full(i - 1, j, k - 1) .and. full(i, j, k - 1))) cycle
                  ex = max(ex, abs(sl%slope_x(i, j, k) - sx_ref))
                  nx_chk = nx_chk + 1
               end do
            end do
            do j = NG + 2, NG + NYP
               do i = NG + 1, NG + NXP
                  if (.not. (full(i, j - 1, k) .and. full(i, j, k) .and. &
                             full(i, j - 1, k - 1) .and. full(i, j, k - 1))) cycle
                  ey = max(ey, abs(sl%slope_y(i, j, k) - sy_ref))
                  ny_chk = ny_chk + 1
               end do
            end do
         end do
         write (msg, '(a,es10.3,a,es10.3,a,i0,a,i0)') "err_x = ", ex, ", err_y = ", ey, &
            ", n_x = ", nx_chk, ", n_y = ", ny_chk
         call check(error, nx_chk >= NXP .and. ny_chk >= NYP, &
                    label//": too few full-cell faces to test: "//trim(msg))
         if (allocated(error)) exit checks
         ! Relative to the slope itself (~1e-4): exact up to round-off.
         ! The zero-datum answer misses by ~ΔD/Δx ~ 1e-3 on sigma / zstar.
         call check(error, ex <= 1.0e-9_wp*abs(sx_ref) .and. &
                    ey <= 1.0e-9_wp*abs(sy_ref), &
                    label//": a uniformly tilted isopycnal over a step must read "// &
                    "the analytic slope: "//trim(msg))
      end block checks
      call teardown(metrics, ms, sl, gm)
   end subroutine check_tilted

   subroutine test_tilted_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      call check_tilted(error, COORD_SIGMA, "sigma")
   end subroutine test_tilted_sigma

   subroutine test_tilted_zstar(error)
      type(error_type), allocatable, intent(out) :: error
      call check_tilted(error, COORD_ZSTAR, "zstar")
   end subroutine test_tilted_zstar

   subroutine test_tilted_zfixed(error)
      type(error_type), allocatable, intent(out) :: error
      call check_tilted(error, COORD_ZFIXED, "z_fixed")
   end subroutine test_tilted_zfixed

end module test_ocean_slopes_datum
