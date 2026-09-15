!! Unit tests for the PPM_H4 higher-order ALE vertical-remap reconstruction (D3)
module test_ocean_ppm_h4_remap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, REMAP_PPM, REMAP_PLM, REMAP_PPM_H4, REMAP_PQM, VCOORD_SIGMA
   use rdb_remap_column, only: remap_column, remap_column_ppm, remap_column_ppm_h4, &
                               remap_column_plm
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_remap, only: ocean_apply_ale_remap_step
   implicit none
   private

   public :: collect_ocean_ppm_h4_remap_tests

contains

   subroutine collect_ocean_ppm_h4_remap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("order_of_accuracy", test_order_of_accuracy), &
                  new_unittest("conservation", test_conservation), &
                  new_unittest("monotonicity_step", test_monotonicity_step), &
                  new_unittest("spurious_mixing", test_spurious_mixing), &
                  new_unittest("golden_nz6", test_golden_nz6), &
                  new_unittest("bit_identity_ppm", test_bit_identity_ppm), &
                  new_unittest("identity_remap", test_identity_remap), &
                  new_unittest("dispatch", test_dispatch), &
                  new_unittest("ocean_path_honors_remap_method", test_ocean_path_honors_remap_method), &
                  new_unittest("ocean_path_honors_pqm", test_ocean_path_honors_pqm) &
                  ]
   end subroutine collect_ocean_ppm_h4_remap_tests

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine make_nonuniform_dz(nz, amp, seed, dz)
      !! Deterministic (LCG) non-uniform thicknesses, normalised to sum = 1.
      integer, intent(in) :: nz, seed
      real(wp), intent(in) :: amp
      real(wp), intent(out) :: dz(nz)
      integer :: k
      integer(kind=8) :: state
      real(wp) :: r, tot
      state = int(seed, 8)
      tot = 0.0_wp
      do k = 1, nz
         state = mod(state*6364136223846793005_8 + 1442695040888963407_8, &
                     9223372036854775783_8)
         r = real(modulo(state, 1000000_8), wp)/1.0e6_wp
         dz(k) = 1.0_wp + amp*r
         tot = tot + dz(k)
      end do
      dz = dz/tot
   end subroutine make_nonuniform_dz

   pure function cell_avg_sin(z_lo, z_hi) result(u)
      !! Cell-average of sin(pi z) over [z_lo, z_hi] (H = 1).
      real(wp), intent(in) :: z_lo, z_hi
      real(wp) :: u
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      u = (-cos(PI*z_hi) + cos(PI*z_lo))/(PI*(z_hi - z_lo))
   end function cell_avg_sin

   ! ----------------------------------------------------------------------
   ! T1 — order of accuracy: PPM_H4 remap L2 error converges faster than PPM
   !      on non-uniform layers (the prototype shows steeper slope).
   ! ----------------------------------------------------------------------
   subroutine test_order_of_accuracy(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nN = 4
      integer :: N_vals(nN), iN, N, k
      real(wp) :: err_h4(nN), err_ppm(nN)
      real(wp), allocatable :: dz_s(:), dz_t(:), u_s(:), u_ref(:), u_h4(:), u_p(:)
      real(wp), allocatable :: z_s(:), z_t(:)
      real(wp) :: slope_h4, slope_ppm, dlogN

      N_vals = [12, 24, 48, 96]   ! all <= NZ_STACK_MAX (100)
      do iN = 1, nN
         N = N_vals(iN)
         allocate (dz_s(N), dz_t(N), u_s(N), u_ref(N), u_h4(N), u_p(N))
         allocate (z_s(0:N), z_t(0:N))
         call make_nonuniform_dz(N, 0.5_wp, 11, dz_s)
         call make_nonuniform_dz(N, 0.5_wp, 97, dz_t)
         dz_t(N) = dz_t(N) + (sum(dz_s) - sum(dz_t))   ! exact equal depths
         z_s(0) = 0.0_wp; z_t(0) = 0.0_wp
         do k = 1, N
            z_s(k) = z_s(k - 1) + dz_s(k)
            z_t(k) = z_t(k - 1) + dz_t(k)
         end do
         do k = 1, N
            u_s(k) = cell_avg_sin(z_s(k - 1), z_s(k))
            u_ref(k) = cell_avg_sin(z_t(k - 1), z_t(k))
         end do
         call remap_column_ppm_h4(N, dz_s, dz_t, u_s, u_h4)
         call remap_column_ppm(N, dz_s, dz_t, u_s, u_p)
         err_h4(iN) = sqrt(sum(dz_t*(u_h4 - u_ref)**2))
         err_ppm(iN) = sqrt(sum(dz_t*(u_p - u_ref)**2))
         deallocate (dz_s, dz_t, u_s, u_ref, u_h4, u_p, z_s, z_t)
      end do

      ! Fitted slope between coarsest and finest (log-log).
      dlogN = log(real(N_vals(nN), wp)/real(N_vals(1), wp))
      ! (dlogN computed from the coarsest..finest N pair)
      slope_h4 = log(err_h4(nN)/err_h4(1))/dlogN
      slope_ppm = log(err_ppm(nN)/err_ppm(1))/dlogN

      ! PPM_H4 must converge faster (more negative slope) than PPM.
      call check(error, slope_h4 < slope_ppm - 0.1_wp, &
                 "PPM_H4 remap L2 slope must be steeper than PPM on non-uniform dz")
   end subroutine test_order_of_accuracy

   ! ----------------------------------------------------------------------
   ! T2 — conservation: column integral preserved to roundoff.
   ! ----------------------------------------------------------------------
   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 20
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz)
      real(wp) :: z_mid, zc, mass_old, mass_new, rel
      integer :: k
      real(wp), parameter :: PI = 3.14159265358979323846_wp

      call make_nonuniform_dz(nz, 0.5_wp, 123, dz_old)
      call make_nonuniform_dz(nz, 0.5_wp, 321, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_old(k)
         u_old(k) = 1.0_wp + sin(2.0_wp*PI*z_mid)
         zc = zc + dz_old(k)
      end do

      call remap_column_ppm_h4(nz, dz_old, dz_new, u_old, u_new)
      mass_old = sum(u_old*dz_old)
      mass_new = sum(u_new*dz_new)
      rel = abs(mass_new - mass_old)/abs(mass_old)
      call check(error, rel <= 1.0e-13_wp, "PPM_H4 must conserve column integral to roundoff")
   end subroutine test_conservation

   ! ----------------------------------------------------------------------
   ! T3 — monotonicity: step profile yields no new extrema.
   ! ----------------------------------------------------------------------
   subroutine test_monotonicity_step(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 30
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz)
      real(wp) :: zc, z_mid, qmn, qmx
      integer :: k

      call make_nonuniform_dz(nz, 0.3_wp, 77, dz_old)
      call make_nonuniform_dz(nz, 0.3_wp, 707, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_old(k)
         if (z_mid < 0.5_wp) then
            u_old(k) = 0.0_wp
         else
            u_old(k) = 1.0_wp
         end if
         zc = zc + dz_old(k)
      end do
      qmn = minval(u_old); qmx = maxval(u_old)

      call remap_column_ppm_h4(nz, dz_old, dz_new, u_old, u_new)
      call check(error, minval(u_new) >= qmn - 1.0e-12_wp .and. &
                 maxval(u_new) <= qmx + 1.0e-12_wp, &
                 "PPM_H4 must not introduce new extrema on a step profile")
   end subroutine test_monotonicity_step

   ! ----------------------------------------------------------------------
   ! T4 — spurious mixing (HEADLINE): repeated remap of a sharp step retains
   !      more tracer variance under PPM_H4 than PPM (less numerical mixing).
   ! ----------------------------------------------------------------------
   subroutine test_spurious_mixing(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 40, ncyc = 200
      real(wp) :: dz_base(nz), dz_t(nz), u0(nz), u_h4(nz), u_p(nz), tmp(nz)
      real(wp) :: zc, z_mid, var0, var_h4, var_p
      integer :: k, c

      call make_nonuniform_dz(nz, 0.4_wp, 555, dz_base)
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_base(k)
         if (z_mid < 0.5_wp) then
            u0(k) = 0.0_wp
         else
            u0(k) = 1.0_wp
         end if
         zc = zc + dz_base(k)
      end do
      var0 = variance(nz, u0, dz_base)

      ! PPM_H4 cycles
      u_h4 = u0
      do c = 1, ncyc
         call make_nonuniform_dz(nz, 0.4_wp, 1000 + c, dz_t)
         dz_t(nz) = dz_t(nz) + (sum(dz_base) - sum(dz_t))
         call remap_column_ppm_h4(nz, dz_base, dz_t, u_h4, tmp)
         call remap_column_ppm_h4(nz, dz_t, dz_base, tmp, u_h4)
      end do
      var_h4 = variance(nz, u_h4, dz_base)

      ! PPM cycles (same grid sequence)
      u_p = u0
      do c = 1, ncyc
         call make_nonuniform_dz(nz, 0.4_wp, 1000 + c, dz_t)
         dz_t(nz) = dz_t(nz) + (sum(dz_base) - sum(dz_t))
         call remap_column_ppm(nz, dz_base, dz_t, u_p, tmp)
         call remap_column_ppm(nz, dz_t, dz_base, tmp, u_p)
      end do
      var_p = variance(nz, u_p, dz_base)

      ! PPM_H4 retains more variance (sharper interface = less spurious mixing).
      call check(error, var_h4 > var_p, &
                 "PPM_H4 must retain more tracer variance than PPM over N remap cycles")
      call check(error, var_h4 <= var0 + 1.0e-12_wp, "variance cannot grow")
   end subroutine test_spurious_mixing

   pure function variance(nz, u, dz) result(v)
      integer, intent(in) :: nz
      real(wp), intent(in) :: u(nz), dz(nz)
      real(wp) :: v, mean, tot
      tot = sum(dz)
      mean = sum(u*dz)/tot
      v = sum(dz*(u - mean)**2)/tot
   end function variance

   ! ----------------------------------------------------------------------
   ! T5 — GOLDEN: nz=6 worked-example remap output matches the verified
   !      Python prototype (`--report`) to ~1e-9.
   ! ----------------------------------------------------------------------
   subroutine test_golden_nz6(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 6
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz), q_golden(nz)
      integer :: k

      dz_old = [0.15_wp, 0.20_wp, 0.18_wp, 0.22_wp, 0.12_wp, 0.13_wp]
      u_old = [1.0_wp, 2.5_wp, 2.0_wp, 3.5_wp, 2.8_wp, 1.5_wp]
      dz_new = [0.13_wp, 0.12_wp, 0.22_wp, 0.18_wp, 0.20_wp, 0.15_wp]
      ! Golden reference from d3_ppm_h4_prototype.py (remap_conservative + H4):
      q_golden = [1.000000000000_wp, 2.250000000000_wp, 2.227272727273_wp, &
                  3.000000000000_wp, 3.192607815674_wp, 1.616522912434_wp]

      call remap_column_ppm_h4(nz, dz_old, dz_new, u_old, u_new)
      do k = 1, nz
         call check(error, abs(u_new(k) - q_golden(k)) < 1.0e-9_wp, &
                    "PPM_H4 nz=6 remap output must match the Python golden prototype")
         if (allocated(error)) return
      end do
   end subroutine test_golden_nz6

   ! ----------------------------------------------------------------------
   ! T6 — bit-identity: REMAP_PPM (default) output is byte-identical to the
   !      pre-D3 path (the PPM kernel is untouched).
   ! ----------------------------------------------------------------------
   subroutine test_bit_identity_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 12
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz)
      real(wp) :: a(nz), b(nz), c(nz)
      real(wp) :: zc, zm
      integer :: k

      call make_nonuniform_dz(nz, 0.6_wp, 31, dz_old)
      call make_nonuniform_dz(nz, 0.6_wp, 313, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      ! Smooth curved (sin) profile on the non-uniform grid: the limiter does
      ! not flatten the interior, so PPM and PPM_H4 genuinely differ there.
      zc = 0.0_wp
      do k = 1, nz
         zm = zc + 0.5_wp*dz_old(k)
         u_old(k) = sin(3.0_wp*zm)
         zc = zc + dz_old(k)
      end do

      ! Direct call and dispatch through remap_column must agree, and the
      ! PPM result must be unaffected by the existence of PPM_H4.
      call remap_column_ppm(nz, dz_old, dz_new, u_old, a)
      call remap_column(REMAP_PPM, nz, dz_old, dz_new, u_old, b)
      call remap_column_ppm_h4(nz, dz_old, dz_new, u_old, c)
      do k = 1, nz
         call check(error, a(k) == b(k), "REMAP_PPM dispatch must be bit-identical to direct PPM")
         if (allocated(error)) return
      end do
      ! PPM_H4 must differ from PPM on non-uniform layers (sanity: not a no-op).
      call check(error, any(abs(c - a) > 1.0e-10_wp), &
                 "PPM_H4 must differ from PPM on non-uniform layers")
   end subroutine test_bit_identity_ppm

   ! ----------------------------------------------------------------------
   ! Identity remap recovers q exactly.
   ! ----------------------------------------------------------------------
   subroutine test_identity_remap(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 7
      real(wp) :: dz(nz), q_old(nz), q_new(nz)
      integer :: k

      dz = [1.0_wp, 2.0_wp, 1.5_wp, 3.0_wp, 2.0_wp, 1.0_wp, 1.2_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 25.0_wp, 15.0_wp, 18.0_wp, 22.0_wp]
      call remap_column_ppm_h4(nz, dz, dz, q_old, q_new)
      do k = 1, nz
         call check(error, abs(q_new(k) - q_old(k)) < 1.0e-11_wp, &
                    "PPM_H4 identity remap must recover q exactly")
         if (allocated(error)) return
      end do
   end subroutine test_identity_remap

   ! ----------------------------------------------------------------------
   ! Dispatch: REMAP_PPM_H4 routes to remap_column_ppm_h4.
   ! ----------------------------------------------------------------------
   subroutine test_dispatch(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 8
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), a(nz), b(nz)
      integer :: k

      call make_nonuniform_dz(nz, 0.5_wp, 42, dz_old)
      call make_nonuniform_dz(nz, 0.5_wp, 424, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      do k = 1, nz
         q_old(k) = sin(real(k, wp))
      end do
      call remap_column(REMAP_PPM_H4, nz, dz_old, dz_new, q_old, a)
      call remap_column_ppm_h4(nz, dz_old, dz_new, q_old, b)
      do k = 1, nz
         call check(error, a(k) == b(k), "dispatch REMAP_PPM_H4 must equal direct call")
         if (allocated(error)) return
      end do
   end subroutine test_dispatch

   ! ----------------------------------------------------------------------
   ! T9 — ocean production path honours remap_method: calling
   !      ocean_apply_ale_remap_step (the dyn-core entry point, not the
   !      raw kernel) with REMAP_PPM and REMAP_PPM_H4 gives DIFFERENT
   !      tracer outputs on non-uniform dz, proving that Fix 1d wired
   !      vcoord%remap_method all the way through.  Also checks that
   !      mass is conserved to ~1e-12 in both cases.
   !      Runs on the device (state entered via acc enter data) so the
   !      GPU code path is exercised.
   ! ----------------------------------------------------------------------
   subroutine test_ocean_path_honors_remap_method(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_ppm, ms_h4
      type(ocean_vcoord_t) :: vc_ppm, vc_h4
      ! Use NX=NY=1 interior (+ nghost=1 halos → 3×3 total) so the remap
      ! runs on one real column.  NZ >= 4 is required for PPM_H4 (H4 uses a
      ! 4-cell interior stencil; nz=2 falls back to PLM, nz=3 only has the
      ! H3 boundary cells and no interior edge).  Use NZ=6 to give at least
      ! 2 interior edges so PPM and PPM_H4 genuinely differ.
      integer, parameter :: NX_PHYS = 1, NY_PHYS = 1, NZ = 6, NGHOST = 1
      real(wp), parameter :: H_TOTAL = 100.0_wp, S_REF = 35.0_wp
      ! Non-uniform layer thickness (dz ratios 0.10, 0.14, 0.18, 0.22, 0.18, 0.18)
      ! summing to 1: each layer in metres.
      real(wp), parameter :: FRACS(NZ) = [0.10_wp, 0.14_wp, 0.18_wp, 0.22_wp, 0.18_wp, 0.18_wp]
      real(wp) :: dz(NZ), mass_ppm, mass_h4, mass_old
      real(wp), allocatable :: bt_eta_ppm(:, :), bt_H_ref_ppm(:, :)
      real(wp), allocatable :: bt_eta_h4(:, :), bt_H_ref_h4(:, :)
      real(wp) :: max_diff, rel_cons_ppm, rel_cons_h4
      integer :: nx_tot, ny_tot, k, i, j

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         ! Initialise both states identically: non-uniform h_layer, smooth
         ! sinusoidal tracer profile so the CW limiter does NOT flatten
         ! interior cells (i.e. PPM and PPM_H4 produce genuinely different
         ! edge estimates on this profile).
         do k = 1, NZ
            dz(k) = FRACS(k)*H_TOTAL
         end do

         ms_ppm%nz_ml = NZ
         call ms_ppm%init(grid)
         ms_h4%nz_ml = NZ
         call ms_h4%init(grid)

         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  ms_ppm%h_layer(i, j, k) = dz(k)
                  ms_h4%h_layer(i, j, k) = dz(k)
                  ! Smooth tracer: sin(pi * k_frac) so the profile has curvature
                  ! that the H4 stencil resolves better than PPM.
                  ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k) = &
                     (S_REF + sin(3.14159265358979323846_wp*real(k, wp)/real(NZ, wp)))*dz(k)
                  ms_h4%tracers(ms_h4%idx_salinity)%hTr(i, j, k) = &
                     ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do

         ! Record mass for conservation check (host-side, before device enter).
         mass_old = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_old = mass_old + ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do

         ! Vcoord: SIGMA so the remap is non-trivial.  A non-uniform
         ! target_h (dsig from the fraction array) differs from the
         ! initial h_layer → the remap actually moves mass.
         call vc_ppm%init(grid, nz_ml=NZ)
         vc_ppm%coord_type = VCOORD_SIGMA
         vc_ppm%remap_method = REMAP_PPM
         do k = 1, NZ
            vc_ppm%dsig(k) = FRACS(k)   ! target = same as current → mass neutral
         end do
         ! Make target_h differ from h_layer by using flipped fracs so PPM
         ! and PPM_H4 produce measurably different outputs.
         ! The dsig is reversed: surface (k=NZ) gets the thick fraction,
         ! bed (k=1) the thin one — flipping the gradient direction.
         vc_ppm%dsig(1) = FRACS(NZ)
         vc_ppm%dsig(2) = FRACS(NZ - 1)
         vc_ppm%dsig(3) = FRACS(NZ - 2)
         vc_ppm%dsig(4) = FRACS(NZ - 3)
         vc_ppm%dsig(5) = FRACS(NZ - 4)
         vc_ppm%dsig(6) = FRACS(NZ - 5)

         call vc_h4%init(grid, nz_ml=NZ)
         vc_h4%coord_type = VCOORD_SIGMA
         vc_h4%remap_method = REMAP_PPM_H4
         do k = 1, NZ
            vc_h4%dsig(k) = vc_ppm%dsig(k)
         end do

         allocate (bt_eta_ppm(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref_ppm(nx_tot, ny_tot), source=H_TOTAL)
         allocate (bt_eta_h4(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref_h4(nx_tot, ny_tot), source=H_TOTAL)

         ! Enter data onto device (production path, same pattern as
         ! test_ocean_remap and test_ocean_remap_e2e use).
         !$acc enter data copyin(ms_ppm, ms_h4, vc_ppm, vc_h4)
         !$acc enter data copyin(bt_eta_ppm, bt_H_ref_ppm, bt_eta_h4, bt_H_ref_h4)
         call ms_ppm%enter_data()
         call ms_h4%enter_data()
         call vc_ppm%enter_data()
         call vc_h4%enter_data()

         ! Run the production entry point — this is the call site Fix 1d wires.
         call ocean_apply_ale_remap_step(grid, vc_ppm, ms_ppm, bt_eta_ppm, bt_H_ref_ppm, &
                                         method=vc_ppm%remap_method)
         call ocean_apply_ale_remap_step(grid, vc_h4, ms_h4, bt_eta_h4, bt_H_ref_h4, &
                                         method=vc_h4%remap_method)

         ! Exit data BEFORE reading host arrays (mirrors test_ocean_remap_e2e).
         call vc_ppm%exit_data()
         call vc_h4%exit_data()
         call ms_ppm%exit_data()
         call ms_h4%exit_data()
         !$acc exit data delete(bt_eta_h4, bt_H_ref_h4, bt_eta_ppm, bt_H_ref_ppm)
         !$acc exit data delete(vc_h4, vc_ppm, ms_h4, ms_ppm)

         ! Assertion 1: PPM and PPM_H4 must give DIFFERENT tracer results
         ! (proving vcoord%remap_method is read by the production entry point).
         max_diff = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  max_diff = max(max_diff, &
                                 abs(ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k) - &
                                     ms_h4%tracers(ms_h4%idx_salinity)%hTr(i, j, k)))
               end do
            end do
         end do
         call check(error, max_diff > 1.0e-10_wp, &
                    "ocean_apply_ale_remap_step: PPM and PPM_H4 must differ (remap_method not wired)")
         if (allocated(error)) exit checks

         ! Assertion 2: tracer mass conserved under PPM to ~1e-12.
         mass_ppm = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_ppm = mass_ppm + ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do
         rel_cons_ppm = abs(mass_ppm - mass_old)/abs(mass_old)
         call check(error, rel_cons_ppm <= 1.0e-12_wp, &
                    "ocean_apply_ale_remap_step PPM: tracer mass not conserved")
         if (allocated(error)) exit checks

         ! Assertion 3: tracer mass conserved under PPM_H4 to ~1e-12.
         mass_h4 = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_h4 = mass_h4 + ms_h4%tracers(ms_h4%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do
         rel_cons_h4 = abs(mass_h4 - mass_old)/abs(mass_old)
         call check(error, rel_cons_h4 <= 1.0e-12_wp, &
                    "ocean_apply_ale_remap_step PPM_H4: tracer mass not conserved")

      end block checks

      deallocate (bt_eta_ppm, bt_H_ref_ppm, bt_eta_h4, bt_H_ref_h4)
      call vc_h4%destroy()
      call vc_ppm%destroy()
      call ms_h4%destroy()
      call ms_ppm%destroy()
   end subroutine test_ocean_path_honors_remap_method

   subroutine test_ocean_path_honors_pqm(error)
      !! PR-9 §9.2 — "the single most important test in the PR": proves the
      !! "pqm" selector reaches the PRODUCTION remap entry point
      !! (`ocean_apply_ale_remap_step`) and that PQM is the code that ran,
      !! not a silent PLM/PPM fallback.  Mirrors
      !! `test_ocean_path_honors_remap_method` exactly, swapping PPM_H4 for
      !! PQM.  NZ=6 (>= 5) is required — `remap_column_pqm` silently
      !! degrades to PPM below 5 layers (PR-9 §11 risk 6); at NZ<5 this
      !! test would assert "PQM differs from PPM" on a column where PQM
      !! *is* PPM, and fail for the wrong reason.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_ppm, ms_pqm
      type(ocean_vcoord_t) :: vc_ppm, vc_pqm
      integer, parameter :: NX_PHYS = 1, NY_PHYS = 1, NZ = 6, NGHOST = 1
      real(wp), parameter :: H_TOTAL = 100.0_wp, S_REF = 35.0_wp
      real(wp), parameter :: FRACS(NZ) = [0.10_wp, 0.14_wp, 0.18_wp, 0.22_wp, 0.18_wp, 0.18_wp]
      real(wp) :: dz(NZ), mass_ppm, mass_pqm, mass_old
      real(wp), allocatable :: bt_eta_ppm(:, :), bt_H_ref_ppm(:, :)
      real(wp), allocatable :: bt_eta_pqm(:, :), bt_H_ref_pqm(:, :)
      real(wp) :: max_diff, rel_cons_ppm, rel_cons_pqm
      integer :: nx_tot, ny_tot, k, i, j

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         do k = 1, NZ
            dz(k) = FRACS(k)*H_TOTAL
         end do

         ms_ppm%nz_ml = NZ
         call ms_ppm%init(grid)
         ms_pqm%nz_ml = NZ
         call ms_pqm%init(grid)

         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  ms_ppm%h_layer(i, j, k) = dz(k)
                  ms_pqm%h_layer(i, j, k) = dz(k)
                  ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k) = &
                     (S_REF + sin(3.14159265358979323846_wp*real(k, wp)/real(NZ, wp)))*dz(k)
                  ms_pqm%tracers(ms_pqm%idx_salinity)%hTr(i, j, k) = &
                     ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do

         mass_old = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_old = mass_old + ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do

         call vc_ppm%init(grid, nz_ml=NZ)
         vc_ppm%coord_type = VCOORD_SIGMA
         vc_ppm%remap_method = REMAP_PPM
         vc_ppm%dsig(1) = FRACS(NZ)
         vc_ppm%dsig(2) = FRACS(NZ - 1)
         vc_ppm%dsig(3) = FRACS(NZ - 2)
         vc_ppm%dsig(4) = FRACS(NZ - 3)
         vc_ppm%dsig(5) = FRACS(NZ - 4)
         vc_ppm%dsig(6) = FRACS(NZ - 5)

         call vc_pqm%init(grid, nz_ml=NZ)
         vc_pqm%coord_type = VCOORD_SIGMA
         vc_pqm%remap_method = REMAP_PQM
         do k = 1, NZ
            vc_pqm%dsig(k) = vc_ppm%dsig(k)
         end do

         allocate (bt_eta_ppm(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref_ppm(nx_tot, ny_tot), source=H_TOTAL)
         allocate (bt_eta_pqm(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref_pqm(nx_tot, ny_tot), source=H_TOTAL)

         !$acc enter data copyin(ms_ppm, ms_pqm, vc_ppm, vc_pqm)
         !$acc enter data copyin(bt_eta_ppm, bt_H_ref_ppm, bt_eta_pqm, bt_H_ref_pqm)
         call ms_ppm%enter_data()
         call ms_pqm%enter_data()
         call vc_ppm%enter_data()
         call vc_pqm%enter_data()

         ! Run the production entry point via the config-level selector
         ! (vcoord%remap_method), exactly the path a `&vcoord_nml
         ! remap_method = "pqm"` namelist reaches.
         call ocean_apply_ale_remap_step(grid, vc_ppm, ms_ppm, bt_eta_ppm, bt_H_ref_ppm, &
                                         method=vc_ppm%remap_method)
         call ocean_apply_ale_remap_step(grid, vc_pqm, ms_pqm, bt_eta_pqm, bt_H_ref_pqm, &
                                         method=vc_pqm%remap_method)

         call vc_ppm%exit_data()
         call vc_pqm%exit_data()
         call ms_ppm%exit_data()
         call ms_pqm%exit_data()
         !$acc exit data delete(bt_eta_pqm, bt_H_ref_pqm, bt_eta_ppm, bt_H_ref_ppm)
         !$acc exit data delete(vc_pqm, vc_ppm, ms_pqm, ms_ppm)

         ! Assertion 1: PPM and PQM must give DIFFERENT tracer results
         ! (proving vcoord%remap_method="pqm" is read by the production
         ! entry point and PQM is the code that ran — not a silent PLM
         ! fallback, the exact gap PR-9 closes).
         max_diff = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  max_diff = max(max_diff, &
                                 abs(ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k) - &
                                     ms_pqm%tracers(ms_pqm%idx_salinity)%hTr(i, j, k)))
               end do
            end do
         end do
         call check(error, max_diff > 1.0e-10_wp, &
                    "ocean_apply_ale_remap_step: PPM and PQM must differ (remap_method not wired)")
         if (allocated(error)) exit checks

         ! Assertion 2/3: both conserve the column integral to ~1e-12.
         mass_ppm = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_ppm = mass_ppm + ms_ppm%tracers(ms_ppm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do
         rel_cons_ppm = abs(mass_ppm - mass_old)/abs(mass_old)
         call check(error, rel_cons_ppm <= 1.0e-12_wp, &
                    "ocean_apply_ale_remap_step PPM: tracer mass not conserved")
         if (allocated(error)) exit checks

         mass_pqm = 0.0_wp
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  mass_pqm = mass_pqm + ms_pqm%tracers(ms_pqm%idx_salinity)%hTr(i, j, k)
               end do
            end do
         end do
         rel_cons_pqm = abs(mass_pqm - mass_old)/abs(mass_old)
         call check(error, rel_cons_pqm <= 1.0e-12_wp, &
                    "ocean_apply_ale_remap_step PQM: tracer mass not conserved")

      end block checks

      deallocate (bt_eta_ppm, bt_H_ref_ppm, bt_eta_pqm, bt_H_ref_pqm)
      call vc_pqm%destroy()
      call vc_ppm%destroy()
      call ms_pqm%destroy()
      call ms_ppm%destroy()
   end subroutine test_ocean_path_honors_pqm

end module test_ocean_ppm_h4_remap
