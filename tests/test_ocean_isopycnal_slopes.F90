!! Analytical + device tests for the isopycnal-slope diagnostics
!! (`rdb_ocean_isopycnal_slopes`).  All cases RUN THE DEVICE KERNEL via
!! `ocean_slopes_compute` (the production entry) with a linear EOS so the
!! neutral slope has a closed form.
!!
!! Bottom-up convention: k=1 bed, k=nz surface; interface K=1 bed,
!! K=nz+1 surface (both forced 0).  Interior interface K straddles layer
!! ka=K (above) and kb=K-1 (below).
module test_ocean_isopycnal_slopes
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute, &
                                         pressure_above_x
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_isopycnal_slopes_tests

   integer, parameter :: NGHOST = 1
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   !! Linear-EOS coefficients in kg/m^3 per unit (Fortran convention):
   !!   rho = rho0 + beta_S*(S-S_ref) - alpha_T*(T-T_ref)
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3/degC
   real(wp), parameter :: BETA_S = 0.78_wp     ! kg/m^3/PSU
   real(wp), parameter :: DT = 1800.0_wp

contains

   subroutine collect_ocean_isopycnal_slopes_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("slopes_linear_tilt", test_linear_tilt), &
                  new_unittest("slopes_tilted_interface", test_tilted_interface), &
                  new_unittest("slopes_flat_zero", test_flat_zero), &
                  new_unittest("slopes_bound", test_bound), &
                  new_unittest("slopes_n2_positive", test_n2_positive), &
                  new_unittest("slopes_vert_fill", test_vert_fill), &
                  new_unittest("slopes_pressure_column", test_pressure_column) &
                  ]
   end subroutine collect_ocean_isopycnal_slopes_tests

   ! ------------------------------------------------------------------
   ! Test 7: interface pressure must INCLUDE the layer directly above.
   ! Guards the off-by-one fixed in pressure_above_x (the loop must run
   ! kk=nz..ka, not ka+1).  Invisible to the linear-EOS slope tests
   ! (pressure-independent derivs), so checked here directly — biases
   ! slopes/N2 under the production Wright/Roquet (pressure-dependent) EOS.
   ! ------------------------------------------------------------------
   subroutine test_pressure_column(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3, NZ = 8
      real(wp), parameter :: DZ = 50.0_wp
      real(wp) :: h(NX, NY, NZ), p, expect
      h = DZ
      ! Interface straddled by layer ka (above): water above = layers ka..nz.
      ! ka = NZ (top interior interface): only the surface layer is above.
      p = pressure_above_x(NX, NY, NZ, h, 2, 2, NZ, RHO0)
      expect = GRAVITY*RHO0*DZ            ! 1 layer, NOT 0 (old off-by-one)
      call check(error, abs(p - expect) < 1.0e-6_wp*expect, &
                 "pressure at interface ka=nz must include the surface layer")
      if (allocated(error)) return
      ! ka = 2: layers 2..NZ are above => (NZ-1) layers.
      p = pressure_above_x(NX, NY, NZ, h, 2, 2, 2, RHO0)
      expect = GRAVITY*RHO0*real(NZ - 1, wp)*DZ
      call check(error, abs(p - expect) < 1.0e-6_wp*expect, &
                 "pressure at interface ka=2 must sum layers 2..nz (incl. ka)")
   end subroutine test_pressure_column

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine setup_slopes(sl, grid, nz, kd_smooth)
      type(ocean_slopes_t), intent(inout) :: sl
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in), optional :: kd_smooth
      call sl%init(grid, nz_ml=nz)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%kd_smooth = 1.0e-6_wp
      if (present(kd_smooth)) sl%kd_smooth = kd_smooth
      sl%min_dz_for_n2 = 1.0_wp
   end subroutine setup_slopes

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine make_eos

   subroutine map_in(ms, metrics, grid, sl)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(ocean_slopes_t), intent(inout) :: sl
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(sl)
      call sl%enter_data()
   end subroutine map_in

   subroutine map_out(ms, metrics, sl)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_slopes_t), intent(inout) :: sl
      ! Pull device-side outputs back to host before unmapping.
      !$acc update self(sl%slope_x, sl%slope_y, sl%n2_u, sl%n2_v)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! Test 1: linear tilt — interior slope = -Gx/Gz
   ! ------------------------------------------------------------------

   subroutine test_linear_tilt(error)
      !! T(x,z) = T0 + Gz*z + Gx*x, S uniform, uniform thickness.
      !! Analytic neutral slope = -Gx/Gz (independent of alpha, rho0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 4, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: GZ = 0.02_wp, GX = 1.0e-4_wp
      real(wp) :: s_analytic, z_c, x_c, t_val, err, mx
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         ! kd_smooth=0 isolates the slope FORMULA from vert_fill_TS: with
         ! smoothing on, the regularizer perturbs the boundary-adjacent
         ! layers of a clean linear profile by O(kappa*dt/dz^2) (~7e-7 here);
         ! vert_fill has its own test (slopes_vert_fill).  Off => the slope
         ! formula is exact at every interior interface.
         call setup_slopes(sl, grid, NZ, kd_smooth=0.0_wp)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ      ! bottom-up centre height
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
               end do
            end do
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         ! Kernel stores the BOUNDED slope drdx/sqrt(drdx^2+drdz^2) (MOM6
         ! convention; |S|<=1).  For a linear field drdx/drdz = -Gx/Gz
         ! exactly, so the bounded value is (-Gx/Gz)/sqrt(1+(Gx/Gz)^2) — it
         ! differs from the raw ratio by O(S^2) (~6e-8 here), which is the
         ! physically-correct target, matched to machine precision.
         s_analytic = (-GX/GZ)/sqrt(1.0_wp + (GX/GZ)**2)
         mx = 0.0_wp
         ! Interior u-faces (i=2..nx) at interior interfaces (K=2..nz).
         do k = 2, NZ
            do j = 1, nj
               do i = 2, ni
                  err = abs(sl%slope_x(i, j, k) - s_analytic)
                  if (err > mx) mx = err
               end do
            end do
         end do
         call check(error, mx < 1.0e-9_wp, &
                    "linear tilt: interior slope_x must match bounded -Gx/Gz")
         if (allocated(error)) exit checks

         ! Bed + surface interfaces must be exactly 0.
         call check(error, all(sl%slope_x(:, :, 1) == 0.0_wp), &
                    "bed interface slope must be 0")
         if (allocated(error)) exit checks
         call check(error, all(sl%slope_x(:, :, NZ + 1) == 0.0_wp), &
                    "surface interface slope must be 0")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_linear_tilt

   ! ------------------------------------------------------------------
   ! Test 2: tilted interface (rotation term) — drdiA=drdiB=0
   ! ------------------------------------------------------------------

   subroutine test_tilted_interface(error)
      !! x-uniform stratification (no along-layer gradient ⇒ drdi=0) but
      !! sloped interfaces (thickness varies with i).  Then the only
      !! contribution is the rotation term:
      !!   slope = -(e_iw - e_i)/dx   (small-slope limit).
      !! Exercises the interface-tilt term the prototype could not.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 8, NY = 4, NZ = 6
      real(wp), parameter :: DX = 1000.0_wp, DZ0 = 50.0_wp
      real(wp), parameter :: GZ = 0.02_wp           ! stable strat (per m)
      real(wp), parameter :: DHDX = 0.5_wp          ! thickness grows with i (m per cell)
      real(wp) :: z_c, t_val, e_iw, e_i, s_expected, err, mx
      integer :: i, j, k, ni, nj, kk
      real(wp) :: hcol
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_slopes(sl, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Thickness varies with i (sloped interfaces), uniform across layers.
         ! T depends ONLY on the layer centre height within the column so
         ! that the along-layer (same-k) horizontal gradient is zero
         ! (drdiA=drdiB=0) but vertical stratification is present.
         do j = 1, nj
            do i = 1, ni
               hcol = DZ0 + DHDX*real(i, wp)
               do k = 1, NZ
                  ms%h_layer(i, j, k) = hcol
               end do
            end do
         end do
         do j = 1, nj
            do i = 1, ni
               do k = 1, NZ
                  ! layer-index-based T (NOT height) ⇒ same value across i
                  ! at fixed k ⇒ zero along-layer horizontal ρ-gradient.
                  z_c = (real(k, wp) - 0.5_wp)
                  t_val = T0 + GZ*DZ0*z_c
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     t_val*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     S0*ms%h_layer(i, j, k)
               end do
            end do
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         ! Expected slope at u-face (i) interface K: -(e_iw - e_i)/dx where
         ! e is the bottom-up cumulative interface height.  Compare at an
         ! interior probe (mid column / mid depth).
         mx = 0.0_wp
         do k = 2, NZ
            do i = 2, ni
               e_iw = 0.0_wp
               e_i = 0.0_wp
               do kk = 1, k - 1
                  e_iw = e_iw + ms%h_layer(i - 1, 2, kk)
                  e_i = e_i + ms%h_layer(i, 2, kk)
               end do
               s_expected = -(e_iw - e_i)/DX
               err = abs(sl%slope_x(i, 2, k) - s_expected)
               if (err > mx) mx = err
            end do
         end do
         ! Small-slope tolerance: the normalization 1/sqrt(1+S^2) differs
         ! from the linearized -(Δe)/dx at O(S^3); slopes here are ~1e-3.
         call check(error, mx < 1.0e-5_wp, &
                    "tilted interface: rotation term must give -(Δe)/dx")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_tilted_interface

   ! ------------------------------------------------------------------
   ! Test 3: flat — horizontally uniform ⇒ slope ≡ 0
   ! ------------------------------------------------------------------

   subroutine test_flat_zero(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 6, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp, GZ = 0.02_wp
      real(wp) :: z_c, t_val, mx
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_slopes(sl, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            t_val = T0 + GZ*z_c
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_val*DZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S0*DZ
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         mx = max(maxval(abs(sl%slope_x)), maxval(abs(sl%slope_y)))
         call check(error, mx < 1.0e-12_wp, &
                    "horizontally uniform column: slope must be identically 0")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_flat_zero

   ! ------------------------------------------------------------------
   ! Test 4: steep front ⇒ |slope| <= 1
   ! ------------------------------------------------------------------

   subroutine test_bound(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 8, NY = 4, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: GZ = 1.0e-4_wp, GX = 0.5_wp  ! near-vertical isopycnals
      real(wp) :: z_c, x_c, t_val, mx
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_slopes(sl, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
               end do
            end do
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         mx = max(maxval(abs(sl%slope_x)), maxval(abs(sl%slope_y)))
         call check(error, mx <= 1.0_wp, &
                    "steep front: |slope| must be bounded by 1")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_bound

   ! ------------------------------------------------------------------
   ! Test 5: stable column ⇒ N2_u > 0 (sign-fix guard)
   ! ------------------------------------------------------------------

   subroutine test_n2_positive(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 4, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp, GZ = 0.02_wp
      real(wp) :: z_c, t_val, n2_min
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_slopes(sl, grid, NZ)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            t_val = T0 + GZ*z_c     ! warmer up ⇒ stable
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = t_val*DZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S0*DZ
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         ! Interior u-faces (i=2..nx) interior interfaces (K=2..nz): N2 > 0.
         n2_min = huge(1.0_wp)
         do k = 2, NZ
            do j = 1, nj
               do i = 2, ni
                  n2_min = min(n2_min, sl%n2_u(i, j, k))
               end do
            end do
         end do
         call check(error, n2_min > 0.0_wp, &
                    "stable column: N2_u must be strictly positive")
         if (allocated(error)) exit checks
         ! Bed + surface N2 forced to 0.
         call check(error, all(sl%n2_u(:, :, 1) == 0.0_wp) .and. &
                    all(sl%n2_u(:, :, NZ + 1) == 0.0_wp), &
                    "bed + surface N2_u must be 0")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_n2_positive

   ! ------------------------------------------------------------------
   ! Test 6: vanished interior layer ⇒ filled value, smooth slope
   ! ------------------------------------------------------------------

   subroutine test_vert_fill(error)
      !! A column with a vanished interior layer (h<H_VANISHED, T/S=0)
      !! must NOT produce a spurious slope spike: vert_fill_TS diffuses
      !! the neighbouring T/S into the gap.  We check that the slope at
      !! the interfaces adjacent to the vanished layer stays close to the
      !! all-wet linear-tilt slope (no spike).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 4, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: GZ = 0.02_wp, GX = 1.0e-4_wp
      integer, parameter :: KBAD = 4
      real(wp) :: z_c, x_c, t_val, spike
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_slopes(sl, grid, NZ)
         ! Larger smoothing so the fill clearly bridges the gap.
         sl%kd_smooth = 1.0e-2_wp
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total

         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c
                  if (k == KBAD) then
                     ! Vanished layer: tiny thickness, garbage (0) tracer.
                     ms%h_layer(i, j, k) = 0.5_wp*H_VANISHED
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 0.0_wp
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 0.0_wp
                  else
                     ms%h_layer(i, j, k) = DZ
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
                  end if
               end do
            end do
         end do

         call map_in(ms, metrics, grid, sl)
         call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
         call map_out(ms, metrics, sl)

         ! No NaN / Inf and no huge spike anywhere (fill prevented the
         ! T=0 ghost from blowing up the gradient).  All slopes |.|<=1.
         spike = maxval(abs(sl%slope_x))
         call check(error, spike <= 1.0_wp, &
                    "vanished layer: filled slope must stay bounded (no spike)")
         if (allocated(error)) exit checks
         ! Sanity: interior slopes away from the gap still resolve the
         ! tilt to a loose tolerance (fill perturbs the gap rows only).
         spike = 0.0_wp
         do k = 2, NZ
            if (k == KBAD .or. k == KBAD + 1) cycle
            do j = 1, nj
               do i = 2, ni
                  spike = max(spike, abs(sl%slope_x(i, j, k) - (-GX/GZ)))
               end do
            end do
         end do
         call check(error, spike < 1.0e-3_wp, &
                    "vanished layer: slopes away from the gap match the tilt")
      end block checks
      call sl%destroy()
      call ms%destroy()
   end subroutine test_vert_fill

end module test_ocean_isopycnal_slopes
