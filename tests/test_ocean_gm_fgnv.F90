!! Analytical + device tests for the FGNV (Ferrari, Griffies, Nurser & Vallis
!! 2010, Ocean Modelling 32, 143-156) GM streamfunction, `&ocean_gm_nml
!! use_fgnv` (`rdb_ocean_gm :: gm_fgnv_solve_column` + the face-column
!! kernels).
!!
!! Bottom-up convention: interface K=1 is the bed (z = -H), K=nl+1 the
!! surface (z = 0).
!!
!!  1. fgnv_cosh_profile     — constant N, c, psi_unlim = P0, gain 1+s: the column solve
!!                             matches (1+s)*P0*[1 - cosh(N(z+H/2)/c)/cosh(NH/2c)]
!!                             (z from the surface, z in [-H,0]) to < 1e-3 of
!!                             P0 at nz=100, with O(dz^2) convergence over
!!                             nz = 25, 50, 100 (error ratio ~ 4).
!!  2. fgnv_cap_monotone     — psi beside the surface and the bed falls
!!                             monotonically as c/(N H) grows.
!!  3. fgnv_c0_bit_exact     — cg1 = 0, use_fgnv gives (1+s)*local uhD/vhD (limiters
!!                             inactive); knob off ignores cg1 bit for bit.
!!  4. fgnv_kernel_sums      — tilted-front column pairs, nz = 6: Sum_k uhD = 0
!!                             and Sum_k vhD = 0 per face to round-off,
!!                             gm_src >= 0, and the FGNV transports really
!!                             differ from the local ones (the test has power).
!!  5. fgnv_pe_decreases     — tilted 2-layer interface under FGNV GM: the
!!                             potential energy decreases monotonically.
!!  6. fgnv_ebt_struct_psiu  — use_ebt_struct (1 above k0, 0 at/below) + FGNV: the
!!                             unlimited psi fed to the solve is exactly 0 below
!!                             the step (cg1=0 makes psi = (1+s) psi_unlim, so
!!                             uhD = 0 under it), and with cg1 > 0 the column
!!                             still closes.
!!  7. fgnv_config_refusals  — validate_config refuses use_fgnv without GM,
!!                             without the wave-speed slot, and non-positive
!!                             fgnv_filter_scale / fgnv_strat_floor, negative
!!                             fgnv_c_min; a valid set passes.
module test_ocean_gm_fgnv
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports, gm_fgnv_solve_column
   use rdb_continuity, only: continuity_t, continuity_gm_apply, TR_MODE_ADVECT
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_gm_fgnv_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: KHTH = 1000.0_wp
   real(wp), parameter :: FSCALE = 0.1_wp    ! fgnv_filter_scale (OM3)
   real(wp), parameter :: GAIN = 1.0_wp + FSCALE

contains

   subroutine collect_ocean_gm_fgnv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fgnv_cosh_profile", test_cosh_profile), &
                  new_unittest("fgnv_cap_monotone", test_cap_monotone), &
                  new_unittest("fgnv_c0_bit_exact", test_c0_bit_exact), &
                  new_unittest("fgnv_kernel_sums", test_kernel_sums), &
                  new_unittest("fgnv_pe_decreases", test_pe_decreases), &
                  new_unittest("fgnv_ebt_struct_psiu", test_ebt_struct_psiu), &
                  new_unittest("fgnv_config_refusals", test_config_refusals) &
                  ]
   end subroutine collect_ocean_gm_fgnv_tests

   ! ------------------------------------------------------------------
   ! Analytic column helpers
   ! ------------------------------------------------------------------

   subroutine solve_uniform(nl, hh, n2c, cc, p0, psi)
      !! Uniform column of `nl` layers, total depth `hh`, constant N^2 and c,
      !! psi_unlim = p0 at the interior interfaces.
      integer, intent(in) :: nl
      real(wp), intent(in) :: hh, n2c, cc, p0
      real(wp), intent(out) :: psi(nl + 1)
      real(wp) :: h(nl), n2(nl + 1), pu(nl + 1)
      h = hh/real(nl, wp)
      n2 = n2c
      pu = p0
      pu(1) = 0.0_wp
      pu(nl + 1) = 0.0_wp
      call gm_fgnv_solve_column(nl, cc*cc, GAIN, h, n2, pu, psi)
   end subroutine solve_uniform

   function cosh_error(nl, hh, n2c, cc, p0) result(emax)
      !! Max |psi - analytic| / p0 over the interior interfaces.
      integer, intent(in) :: nl
      real(wp), intent(in) :: hh, n2c, cc, p0
      real(wp) :: emax
      real(wp) :: psi(nl + 1), z, xi, ana, nn
      integer :: k
      nn = sqrt(n2c)
      call solve_uniform(nl, hh, n2c, cc, p0, psi)
      emax = 0.0_wp
      do k = 2, nl
         z = -hh + real(k - 1, wp)*hh/real(nl, wp)   ! from the surface (0), bed at -H
         xi = z + 0.5_wp*hh
         ana = GAIN*p0*(1.0_wp - cosh(nn*xi/cc)/cosh(nn*hh/(2.0_wp*cc)))
         emax = max(emax, abs(psi(k) - ana)/p0)
      end do
   end function cosh_error

   ! ------------------------------------------------------------------
   ! Test 1: cosh profile + second-order convergence
   ! ------------------------------------------------------------------
   subroutine test_cosh_profile(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HH = 1000.0_wp, N2C = 1.0e-5_wp, P0 = 7.0_wp
      real(wp), parameter :: CC = 0.5_wp   ! c/N = 158 m: resolved at 10 m
      real(wp) :: e25, e50, e100, r1, r2
      checks: block
         e25 = cosh_error(25, HH, N2C, CC, P0)
         e50 = cosh_error(50, HH, N2C, CC, P0)
         e100 = cosh_error(100, HH, N2C, CC, P0)
         call check(error, e100 < 1.0e-3_wp, "cosh profile: max rel error at nz=100 must be < 1e-3")
         if (allocated(error)) exit checks
         call check(error, e100 > 0.0_wp, "cosh profile: discretisation error must be nonzero")
         if (allocated(error)) exit checks
         r1 = e25/e50
         r2 = e50/e100
         call check(error, abs(r1 - 4.0_wp) < 0.5_wp, "cosh profile: e(25)/e(50) must be ~4 (O(dz^2))")
         if (allocated(error)) exit checks
         call check(error, abs(r2 - 4.0_wp) < 0.5_wp, "cosh profile: e(50)/e(100) must be ~4 (O(dz^2))")
      end block checks
   end subroutine test_cosh_profile

   ! ------------------------------------------------------------------
   ! Test 2: caps pull psi to zero as c/(N H) grows
   ! ------------------------------------------------------------------
   subroutine test_cap_monotone(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NL = 50
      real(wp), parameter :: HH = 1000.0_wp, N2C = 1.0e-5_wp, P0 = 3.0_wp
      real(wp) :: ratio(3), cc, psi(NL + 1), top(3), bed(3), nh
      integer :: m
      checks: block
         nh = sqrt(N2C)*HH
         ratio = [0.02_wp, 0.1_wp, 0.5_wp]
         do m = 1, 3
            cc = ratio(m)*nh
            call solve_uniform(NL, HH, N2C, cc, P0, psi)
            top(m) = psi(NL)
            bed(m) = psi(2)
         end do
         call check(error, top(1) > 0.0_wp .and. bed(1) > 0.0_wp .and. top(1) < P0 .and. bed(1) < P0, &
                    "cap: psi beside the caps must lie in (0, psi_unlim)")
         if (allocated(error)) exit checks
         call check(error, top(1) > top(2) .and. top(2) > top(3), &
                    "cap: surface-adjacent psi must fall monotonically with c/(NH)")
         if (allocated(error)) exit checks
         call check(error, bed(1) > bed(2) .and. bed(2) > bed(3), &
                    "cap: bed-adjacent psi must fall monotonically with c/(NH)")
         if (allocated(error)) exit checks
         call check(error, top(3) < 0.1_wp*P0, "cap: strong filtering must pull psi well below psi_unlim")
      end block checks
   end subroutine test_cap_monotone

   ! ------------------------------------------------------------------
   ! Kernel harness helpers
   ! ------------------------------------------------------------------

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine make_eos

   subroutine setup_slopes(sl, grid, nz)
      type(ocean_slopes_t), intent(inout) :: sl
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call sl%init(grid, nz_ml=nz)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%kd_smooth = 0.0_wp
      sl%min_dz_for_n2 = 1.0_wp
   end subroutine setup_slopes

   subroutine setup_gm(gm, grid, nz)
      type(ocean_gm_t), intent(inout) :: gm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%khth = KHTH
      gm%khth_max_cfl = 1.0e6_wp
      gm%khth_slope_max = 0.01_wp
      gm%rho0 = RHO0
   end subroutine setup_gm

   subroutine fill_tilted_TS(ms, grid, nz, dz, gx, gy)
      !! T tilts in x (and y when gy /= 0), stable in z: tilted isopycnals.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz, gx, gy
      real(wp), parameter :: GZ = 0.02_wp
      integer :: i, j, k
      real(wp) :: t_val
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               t_val = T0 + GZ*(real(k, wp) - 0.5_wp)*dz + gx*real(i, wp)*grid%dx &
                       + gy*real(j, wp)*grid%dy
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*dz
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*dz
            end do
         end do
      end do
   end subroutine fill_tilted_TS

   subroutine run_pair(nz, cg1_val, fgnv_scale, slope_max, uh_loc, vh_loc, uh_fg, vh_fg, src_fg, uh_nocg)
      !! Tilted-front column pairs (x and y tilt); compute GM transports with
      !! use_fgnv off and on (same state) and return both.
      integer, intent(in) :: nz
      real(wp), intent(in) :: cg1_val, fgnv_scale, slope_max
      real(wp), allocatable, intent(out) :: uh_loc(:, :, :), vh_loc(:, :, :)
      real(wp), allocatable, intent(out) :: uh_fg(:, :, :), vh_fg(:, :, :), src_fg(:, :), uh_nocg(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      real(wp), allocatable :: cg1(:, :)
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      call grid%init(8, 6, NGHOST, DX, DX)
      ms%nz_ml = nz
      call ms%init(grid)
      call make_cartesian_metrics(metrics, grid)
      call setup_slopes(sl, grid, nz)
      call setup_gm(gm, grid, nz)
      call make_eos(eos)
      call fill_tilted_TS(ms, grid, nz, 20.0_wp, 2.0e-5_wp, 1.0e-5_wp)
      allocate (cg1(grid%nx_total, grid%ny_total), source=cg1_val)
      gm%fgnv_filter_scale = fgnv_scale
      gm%khth_slope_max = slope_max
      gm%khth = 100.0_wp

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call sl%set_bathymetry(sum(ms%h_layer, dim=3))
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
      !$acc enter data copyin(cg1)

      call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
      gm%use_fgnv = .false.
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT, cg1=cg1)
      !$acc update self(gm%uhD, gm%vhD)
      uh_loc = gm%uhD
      vh_loc = gm%vhD
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
      !$acc update self(gm%uhD)
      uh_nocg = gm%uhD
      gm%use_fgnv = .true.
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT, cg1=cg1)
      !$acc update self(gm%uhD, gm%vhD, gm%gm_src)
      uh_fg = gm%uhD
      vh_fg = gm%vhD
      src_fg = gm%gm_src

      !$acc exit data delete(cg1)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pair

   ! ------------------------------------------------------------------
   ! Test 3: c -> 0 limit is (1+s)*local GM; knob off is untouched
   ! ------------------------------------------------------------------
   subroutine test_c0_bit_exact(error)
      !! cg1 = 0: psi = (1+s) psi_unlim exactly.  With the safe-streamfunction
      !! blend and the availability bound inactive (huge slope_max, small
      !! KhTh) the limiter chain is the identity, so uhD/vhD are (1+s) times
      !! the local ones to round-off.  Separately, use_fgnv=.false. must not
      !! depend on whether a cg1 field is handed in at all (bit for bit).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ul(:, :, :), vl(:, :, :), uf(:, :, :), vf(:, :, :), sf(:, :), un(:, :, :)
      checks: block
         call run_pair(6, 0.0_wp, FSCALE, 1.0e6_wp, ul, vl, uf, vf, sf, un)
         call check(error, maxval(abs(ul)) > 0.0_wp, "c0: the local GM must move water (test has power)")
         if (allocated(error)) exit checks
         call check(error, all(un == ul), "off: uhD must not depend on the cg1 argument (bit-identical)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(uf - GAIN*ul)) <= 1.0e-12_wp*maxval(abs(ul)), &
                    "c0: uhD must be (1+s) times the local transport")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vf - GAIN*vl)) <= 1.0e-12_wp*maxval(abs(vl)), &
                    "c0: vhD must be (1+s) times the local transport")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(uf - ul)) > 0.05_wp*maxval(abs(ul)), &
                    "c0: the amplification must be visible (test has power)")
      end block checks
   end subroutine test_c0_bit_exact

   ! ------------------------------------------------------------------
   ! Test 4: column sums, gm_src, and the FGNV answer really differs
   ! ------------------------------------------------------------------
   subroutine test_kernel_sums(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ul(:, :, :), vl(:, :, :), uf(:, :, :), vf(:, :, :), sf(:, :)
      real(wp), allocatable :: un(:, :, :)
      real(wp) :: sumu, sumv, scale_u, dif
      integer :: i, j, nx, ny
      checks: block
         call run_pair(6, 0.4_wp, 1.0_wp, 0.01_wp, ul, vl, uf, vf, sf, un)
         nx = size(sf, 1)
         ny = size(sf, 2)
         scale_u = maxval(abs(uf))
         call check(error, scale_u > 0.0_wp, "sums: FGNV must move water")
         if (allocated(error)) exit checks
         sumu = 0.0_wp
         do j = 1, ny
            do i = 1, nx + 1
               sumu = max(sumu, abs(sum(uf(i, j, :))))
            end do
         end do
         sumv = 0.0_wp
         do j = 1, ny + 1
            do i = 1, nx
               sumv = max(sumv, abs(sum(vf(i, j, :))))
            end do
         end do
         call check(error, sumu <= 1.0e-12_wp*scale_u, "sums: Sum_k uhD must vanish to round-off")
         if (allocated(error)) exit checks
         call check(error, sumv <= 1.0e-12_wp*max(scale_u, maxval(abs(vf))), &
                    "sums: Sum_k vhD must vanish to round-off")
         if (allocated(error)) exit checks
         call check(error, minval(sf) >= 0.0_wp, "sums: gm_src must be >= 0")
         if (allocated(error)) exit checks
         dif = maxval(abs(uf - ul))
         call check(error, dif > 1.0e-3_wp*scale_u, &
                    "sums: FGNV transports must differ from the local GM (c = 0.4 m/s)")
      end block checks
   end subroutine test_kernel_sums

   ! ------------------------------------------------------------------
   ! Test 5: PE decreases under FGNV GM (tilted 2-layer interface)
   ! ------------------------------------------------------------------
   function column_pe(ms, grid, nz) result(pe)
      !! PE = Sum g * rho' * h * z_mid, z_mid measured up from the bed.
      type(multilayer_state_t), intent(in) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp) :: pe
      integer :: i, j, k
      real(wp) :: zb, tt, rho
      pe = 0.0_wp
      do j = NGHOST + 1, NGHOST + grid%ny_phys
         do i = NGHOST + 1, NGHOST + grid%nx_phys
            zb = 0.0_wp
            do k = 1, nz
               tt = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ms%h_layer(i, j, k)
               rho = RHO0 - ALPHA_T*(tt - T0)
               pe = pe + GRAVITY*rho*ms%h_layer(i, j, k)*(zb + 0.5_wp*ms%h_layer(i, j, k))
               zb = zb + ms%h_layer(i, j, k)
            end do
         end do
      end do
   end function column_pe

   subroutine test_pe_decreases(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(continuity_t) :: ct
      type(eos_t) :: eos
      real(wp), allocatable :: cg1(:, :)
      integer, parameter :: NX = 8, NY = 3, NZ = 2, NSTEP = 40
      real(wp), parameter :: DX = 2000.0_wp, HTOT = 100.0_wp, DTL = 600.0_wp
      real(wp) :: pe0, pe, prev, h1w
      integer :: i, j, it
      logical :: monotone
      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_slopes(sl, grid, NZ)
         call setup_gm(gm, grid, NZ)
         gm%use_fgnv = .true.
         gm%fgnv_filter_scale = 0.1_wp
         call make_eos(eos)
         allocate (cg1(grid%nx_total, grid%ny_total), source=0.5_wp)
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               h1w = 0.5_wp*HTOT + 20.0_wp*(real(grid%nx_phys + 1, wp)*0.5_wp - real(i - NGHOST, wp)) &
                     /real(grid%nx_phys, wp)
               h1w = min(max(h1w, 10.0_wp), HTOT - 10.0_wp)
               ms%h_layer(i, j, 1) = h1w
               ms%h_layer(i, j, 2) = HTOT - h1w
               ms%tracers(ms%idx_temperature)%hTr(i, j, 1) = 8.0_wp*h1w
               ms%tracers(ms%idx_temperature)%hTr(i, j, 2) = 16.0_wp*(HTOT - h1w)
               ms%tracers(ms%idx_salinity)%hTr(i, j, 1) = S0*h1w
               ms%tracers(ms%idx_salinity)%hTr(i, j, 2) = S0*(HTOT - h1w)
            end do
         end do

         !$acc enter data copyin(ms)
         call ms%enter_data()
         call sl%set_bathymetry(sum(ms%h_layer, dim=3))
         !$acc enter data copyin(sl)
         call sl%enter_data()
         !$acc enter data copyin(gm)
         call gm%enter_data()
         !$acc enter data copyin(ct)
         call ct%enter_data()
         !$acc enter data copyin(cg1)

         pe0 = column_pe(ms, grid, NZ)
         prev = pe0
         monotone = .true.
         pe = pe0
         do it = 1, NSTEP
            call ocean_slopes_compute(grid, metrics, eos, sl, ms, DTL)
            call gm_compute_transports(grid, metrics, gm, sl, ms, DTL, cg1=cg1)
            call continuity_gm_apply(grid, metrics, ct, ms, gm, DTL, 1.0_wp, TR_MODE_ADVECT)
            !$acc update self(ms%h_layer)
            !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
            pe = column_pe(ms, grid, NZ)
            if (pe > prev + 1.0e-12_wp*abs(pe0)) monotone = .false.
            prev = pe
         end do

         !$acc exit data delete(cg1)
         call ct%exit_data()
         !$acc exit data delete(ct)
         call gm%exit_data()
         !$acc exit data delete(gm)
         call sl%exit_data()
         !$acc exit data delete(sl)
         call ms%exit_data()
         !$acc exit data delete(ms)

         call check(error, pe < pe0 - 1.0e-9_wp*abs(pe0), "PE must decrease under FGNV GM")
         if (allocated(error)) exit checks
         call check(error, monotone, "PE must be non-increasing every step")
      end block checks
      call ct%destroy()
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_pe_decreases

   ! ------------------------------------------------------------------
   ! Test 6: EBT structure feeds the FGNV solve
   ! ------------------------------------------------------------------
   subroutine run_ebt(nz, k0, cg1_val, uh_ebt, uh_both, vh_both)
      !! use_ebt_struct with a step structure (1 for ka > k0, 0 for ka <= k0):
      !! returns the EBT-only and EBT+FGNV transports (same state).
      integer, intent(in) :: nz, k0
      real(wp), intent(in) :: cg1_val
      real(wp), allocatable, intent(out) :: uh_ebt(:, :, :), uh_both(:, :, :), vh_both(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      real(wp), allocatable :: cg1(:, :), ebt(:, :, :)
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      integer :: k
      call grid%init(8, 6, NGHOST, DX, DX)
      ms%nz_ml = nz
      call ms%init(grid)
      call make_cartesian_metrics(metrics, grid)
      call setup_slopes(sl, grid, nz)
      call setup_gm(gm, grid, nz)
      call make_eos(eos)
      call fill_tilted_TS(ms, grid, nz, 20.0_wp, 2.0e-5_wp, 1.0e-5_wp)
      gm%khth_slope_max = 1.0e6_wp
      gm%khth = 100.0_wp
      gm%fgnv_filter_scale = FSCALE
      gm%use_ebt_struct = .true.
      allocate (cg1(grid%nx_total, grid%ny_total), source=cg1_val)
      allocate (ebt(grid%nx_total, grid%ny_total, nz))
      do k = 1, nz
         ebt(:, :, k) = merge(1.0_wp, 0.0_wp, k > k0)
      end do
      !$acc enter data copyin(ms)
      call ms%enter_data()
      call sl%set_bathymetry(sum(ms%h_layer, dim=3))
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
      !$acc enter data copyin(cg1, ebt)
      call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
      gm%use_fgnv = .false.
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT, ebt_struct=ebt, cg1=cg1)
      !$acc update self(gm%uhD)
      uh_ebt = gm%uhD
      gm%use_fgnv = .true.
      call gm_compute_transports(grid, metrics, gm, sl, ms, DT, ebt_struct=ebt, cg1=cg1)
      !$acc update self(gm%uhD, gm%vhD)
      uh_both = gm%uhD
      vh_both = gm%vhD
      !$acc exit data delete(cg1, ebt)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call gm%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_ebt

   subroutine test_ebt_struct_psiu(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ue(:, :, :), ub(:, :, :), vb(:, :, :)
      integer, parameter :: NZ = 6, K0 = 3
      real(wp) :: scale_e, sm
      integer :: i, j
      checks: block
         ! cg1 = 0: psi = (1+s)*psi_unlim, no smoothing, so a zero input
         ! below the step shows up as exactly zero transport under it.
         call run_ebt(NZ, K0, 0.0_wp, ue, ub, vb)
         scale_e = maxval(abs(ue))
         call check(error, scale_e > 0.0_wp, "ebt+fgnv: the structured GM must move water")
         if (allocated(error)) exit checks
         call check(error, all(ub(:, :, 1:K0 - 1) == 0.0_wp), &
                    "ebt+fgnv: psi_unlim fed to the solve must be exactly 0 below the step")
         if (allocated(error)) exit checks
         call check(error, all(ue(:, :, 1:K0 - 1) == 0.0_wp), &
                    "ebt only: transport below the step must be exactly 0 (control)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ub - GAIN*ue)) <= 1.0e-12_wp*scale_e, &
                    "ebt+fgnv at c=0 must be (1+s) times the EBT-only transport")
         if (allocated(error)) exit checks
         ! cg1 > 0: the BVP smooths, so psi below the step need not vanish, but
         ! the column must still close.
         call run_ebt(NZ, K0, 0.4_wp, ue, ub, vb)
         sm = 0.0_wp
         do j = 1, size(ub, 2)
            do i = 1, size(ub, 1)
               sm = max(sm, abs(sum(ub(i, j, :))))
            end do
         end do
         call check(error, maxval(abs(ub)) > 0.0_wp, "ebt+fgnv (c>0): must move water")
         if (allocated(error)) exit checks
         call check(error, sm <= 1.0e-12_wp*maxval(abs(ub)), "ebt+fgnv (c>0): Sum_k uhD must vanish")
      end block checks
   end subroutine test_ebt_struct_psiu

   ! ------------------------------------------------------------------
   ! Test 7: configure refusals
   ! ------------------------------------------------------------------
   subroutine parse_case(cfg, extra)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: extra
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
   end subroutine parse_case

   function validates(extra) result(ok)
      character(len=*), intent(in) :: extra
      logical :: ok
      type(config_t) :: cfg
      integer :: ierr
      call parse_case(cfg, extra)
      call validate_config(cfg, ierr=ierr)
      ok = (ierr == OCEAN_STATUS_OK)
   end function validates

   subroutine test_config_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: SLO = "&ocean_slopes_nml enable = .true. /"
      character(len=*), parameter :: WS = "&ocean_wavespeed_nml enable = .true. /"
      character(len=*), parameter :: GM = "&ocean_gm_nml enable = .true., khth = 500.0, "
      checks: block
         call check(error, validates(SLO//WS//GM//"use_fgnv = .true. /"), &
                    "refusals: a valid FGNV set must pass (positive control)")
         if (allocated(error)) exit checks
         call check(error, validates(SLO//WS//GM//"use_fgnv = .true., fgnv_filter_scale = 0.1, "// &
                                     "fgnv_c_min = 0.01 /"), "refusals: OM3-style values must pass")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//WS//"&ocean_gm_nml enable = .false., use_fgnv = .true. /"), &
                    "refusals: use_fgnv without GM must be refused")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//GM//"use_fgnv = .true. /"), &
                    "refusals: use_fgnv without the wave-speed slot must be refused")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//WS//GM//"use_fgnv = .true., fgnv_filter_scale = 0.0 /"), &
                    "refusals: fgnv_filter_scale = 0 must be refused")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//WS//GM//"use_fgnv = .true., fgnv_filter_scale = -1.0 /"), &
                    "refusals: fgnv_filter_scale < 0 must be refused")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//WS//GM//"use_fgnv = .true., fgnv_strat_floor = 0.0 /"), &
                    "refusals: fgnv_strat_floor = 0 must be refused")
         if (allocated(error)) exit checks
         call check(error,.not. validates(SLO//WS//GM//"use_fgnv = .true., fgnv_c_min = -0.1 /"), &
                    "refusals: fgnv_c_min < 0 must be refused")
      end block checks
   end subroutine test_config_refusals

end module test_ocean_gm_fgnv
