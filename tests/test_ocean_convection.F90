!! Unit tests for Brunt-Vaisala-triggered convective adjustment
!! (`vmix_apply_convection` in `rdb_ocean_vmix`).
!!
!! Convective adjustment is a CONTRIBUTOR into `kv`/`kt`: where the
!! interior N^2 < n2_thresh (dense-over-light), it raises
!! `kt -> max(kt, kd_conv)` and `kv -> max(kv, prandtl_conv*kd_conv)`,
!! strictly below the active surface boundary-layer depth.  `kd_conv`
!! defaults to 1 m^2/s -- ~100x PP81's own `Ri<0` clip ceiling
!! (`pp81_nu_bg + pp81_nu0` ~= 1.01e-2 m^2/s) -- which is why PP81's own
!! "convective" test case is NOT full mixing (see the renamed
!! `pp81_unstable_ri_clip` in `test_ocean_pp81`).
!!
!! Every case drives `vmix_compute_pp81 -> vmix_apply_convection ->
!! vmix_assemble` through the `run_conv` harness (mirrors
!! `test_ocean_pp81`'s `run_pp81`).  `run_no_conv` runs the same chain
!! WITHOUT the convection call, for the default-off / stable-column
!! bit-identity discriminators.  A single generic `bld(:, :)` argument
!! stands in for both KPP's `vmix%bl_depth` and EPBL's `epbl%mld` --
!! the kernel treats either identically (one generic dummy, no
!! special-casing), so a `vmix%bl_depth`-only harness exercises the same
!! code EPBL would drive.
!!
!! Cases:
!!   * conv_disabled_bit_identical      -- default-off contract.
!!   * conv_stable_column_bit_identical -- the highest-consequence sign
!!     check: a stable column must not be touched.
!!   * conv_unstable_sets_kd_conv       -- exact magnitude + Prandtl.
!!   * conv_below_bl_only               -- BL ownership (kOBL+1 analogue).
!!   * conv_n2_thresh_live              -- n2_thresh is a live knob.
!!   * conv_homogenises_inversion       -- analytical: vdiff integration
!!     removes an interior inversion at kd_conv, with a PP81-only
!!     control that (over the same integration) does not.
!!   * conv_vanished_interface_skipped  -- H_VANISHED gate, no NaN/1/0.
module test_ocean_convection
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81, vmix_apply_convection, &
                             vmix_assemble
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_tracers
   implicit none
   private

   public :: collect_ocean_convection_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_convection_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("conv_disabled_bit_identical", test_disabled_bit_identical), &
                  new_unittest("conv_stable_column_bit_identical", test_stable_bit_identical), &
                  new_unittest("conv_unstable_sets_kd_conv", test_unstable_sets_kd_conv), &
                  new_unittest("conv_below_bl_only", test_below_bl_only), &
                  new_unittest("conv_n2_thresh_live", test_n2_thresh_live), &
                  new_unittest("conv_homogenises_inversion", test_homogenises_inversion), &
                  new_unittest("conv_vanished_interface_skipped", test_vanished_interface_skipped) &
                  ]
   end subroutine collect_ocean_convection_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine run_conv(grid, ms, vmix, bld)
      !! Map ms + vmix to the device, run
      !! pp81 -> apply_convection -> assemble, pull kv/kt back to the
      !! host.  `bld` must be set on the host BEFORE this call so
      !! `copyin` carries it (matches CLAUDE.md's mem:separate
      !! contract -- no implicit host<->device copies).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      real(wp), intent(in) :: bld(:, :)
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_apply_convection(grid, vmix, ms, bld)
      call vmix_assemble(grid, vmix, ms)
      !$acc update self(vmix%kv, vmix%kt)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_conv

   subroutine run_no_conv(grid, ms, vmix)
      !! Same chain WITHOUT the convection call -- the "PP81 only"
      !! control used by the bit-identity + discriminator cases.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      call vmix_assemble(grid, vmix, ms)
      !$acc update self(vmix%kv, vmix%kt)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_no_conv

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_disabled_bit_identical(error)
      !! Default-off contract: with `conv_enable = .false.`, calling
      !! `vmix_apply_convection` (the shim early-returns) must be
      !! byte-identical to never calling it at all -- on an UNSTABLE
      !! column, so this is not vacuously true because nothing would
      !! have triggered anyway.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp) :: max_diff_kv, max_diff_kt
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms_a%nz_ml = NZ; call ms_a%init(grid)
         ms_b%nz_ml = NZ; call ms_b%init(grid)
         call vmix_a%init(grid, nz_ml=NZ)
         call vmix_b%init(grid, nz_ml=NZ)
         vmix_a%conv_enable = .false.
         vmix_a%conv_kd = 1.0_wp

         ms_a%h_layer = H_LAYER; ms_b%h_layer = H_LAYER
         ms_a%u_face_x_layer = 0.0_wp; ms_b%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms_a%rho_layer(:, :, k) = 1030.0_wp + 1.0_wp*real(k - 1, wp)
            ms_b%rho_layer(:, :, k) = 1030.0_wp + 1.0_wp*real(k - 1, wp)
         end do

         call run_conv(grid, ms_a, vmix_a, vmix_a%bl_depth)
         call run_no_conv(grid, ms_b, vmix_b)

         max_diff_kv = maxval(abs(vmix_a%kv - vmix_b%kv))
         max_diff_kt = maxval(abs(vmix_a%kt - vmix_b%kt))

         call check(error, max_diff_kv < 1.0e-15_wp, &
                    "conv disabled: kv not bit-identical to no-conv run")
         if (allocated(error)) exit checks
         call check(error, max_diff_kt < 1.0e-15_wp, &
                    "conv disabled: kt not bit-identical to no-conv run")

      end block checks
      call vmix_a%destroy(); call vmix_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_disabled_bit_identical

   subroutine test_stable_bit_identical(error)
      !! THE highest-consequence sign check: a monotonically stable
      !! column (N^2 > 0 at every interior interface) run with
      !! `conv_enable = .true.` must be byte-identical to the no-conv
      !! run.  A sign error here would apply 1 m^2/s to the whole
      !! stratified ocean.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp) :: max_diff_kv, max_diff_kt
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms_a%nz_ml = NZ; call ms_a%init(grid)
         ms_b%nz_ml = NZ; call ms_b%init(grid)
         call vmix_a%init(grid, nz_ml=NZ)
         call vmix_b%init(grid, nz_ml=NZ)
         vmix_a%conv_enable = .true.
         vmix_a%conv_kd = 1.0_wp

         ms_a%h_layer = H_LAYER; ms_b%h_layer = H_LAYER
         ms_a%u_face_x_layer = 0.0_wp; ms_b%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp
         ! Stable: heavier at bed, decreasing upward (test_quiescent's
         ! profile in test_ocean_pp81).
         do k = 1, NZ
            ms_a%rho_layer(:, :, k) = 1030.0_wp - 1.0_wp*real(k - 1, wp)
            ms_b%rho_layer(:, :, k) = 1030.0_wp - 1.0_wp*real(k - 1, wp)
         end do

         call run_conv(grid, ms_a, vmix_a, vmix_a%bl_depth)
         call run_no_conv(grid, ms_b, vmix_b)

         max_diff_kv = maxval(abs(vmix_a%kv - vmix_b%kv))
         max_diff_kt = maxval(abs(vmix_a%kt - vmix_b%kt))

         call check(error, max_diff_kv < 1.0e-15_wp, &
                    "conv on, stable column: kv not bit-identical to no-conv run")
         if (allocated(error)) exit checks
         call check(error, max_diff_kt < 1.0e-15_wp, &
                    "conv on, stable column: kt not bit-identical to no-conv run")

      end block checks
      call vmix_a%destroy(); call vmix_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_stable_bit_identical

   subroutine test_unstable_sets_kd_conv(error)
      !! Inverted column, no BL (`bld = 0`): every interior interface
      !! must hit EXACTLY `kt = kd_conv`, `kv = prandtl_conv*kd_conv`.
      !! `prandtl_conv = 2.0` (non-default) so a `kv = kd_conv`
      !! copy-paste bug fails.  Boundary interfaces (k=1, k=NZ+1) stay
      !! exactly zero (closed-BC contract, untouched by convection).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: KD_CONV = 1.0_wp
      real(wp), parameter :: PRANDTL_CONV = 2.0_wp
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%conv_enable = .true.
         vmix%conv_kd = KD_CONV
         vmix%conv_prandtl = PRANDTL_CONV

         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp + 1.0_wp*real(k - 1, wp)
         end do
         vmix%bl_depth = 0.0_wp

         call run_conv(grid, ms, vmix, vmix%bl_depth)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2

         do k = 2, NZ
            call check(error, abs(vmix%kt(i_probe, j_probe, k) - KD_CONV) < TOL, &
                       "unstable: kt did not hit kd_conv exactly")
            if (allocated(error)) exit checks
            call check(error, &
                       abs(vmix%kv(i_probe, j_probe, k) - PRANDTL_CONV*KD_CONV) < TOL, &
                       "unstable: kv did not hit prandtl_conv*kd_conv exactly")
            if (allocated(error)) exit checks
         end do

         call check(error, &
                    maxval(abs(vmix%kt(:, :, 1))) < 1.0e-15_wp .and. &
                    maxval(abs(vmix%kv(:, :, 1))) < 1.0e-15_wp, &
                    "unstable: bed interface not exactly zero")
         if (allocated(error)) exit checks
         call check(error, &
                    maxval(abs(vmix%kt(:, :, NZ + 1))) < 1.0e-15_wp .and. &
                    maxval(abs(vmix%kv(:, :, NZ + 1))) < 1.0e-15_wp, &
                    "unstable: surface interface not exactly zero")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_unstable_sets_kd_conv

   subroutine test_below_bl_only(error)
      !! nz=10, h=20m (H=200m), bld=100.0, unstable everywhere.
      !! Interfaces with z_int(k) < 100 (k=10..7; z_int = 20,40,60,80)
      !! must stay at the PP81 unstable value (kappa_bg+nu0); interfaces
      !! with z_int(k) >= 100 (k=6..2; z_int = 100,120,...,180) must hit
      !! kd_conv.  Exercises the ">=" boundary exactly at k=6
      !! (z_int(6) == bld == 100).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ_BL = 10
      real(wp), parameter :: H_LAYER = 20.0_wp
      real(wp), parameter :: BLD = 100.0_wp
      real(wp), parameter :: KD_CONV = 1.0_wp
      real(wp), parameter :: PP81_UNSTABLE_KT = 1.0e-5_wp + 1.0e-2_wp
      real(wp), parameter :: TOL = 1.0e-10_wp
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ_BL; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ_BL)
         vmix%conv_enable = .true.
         vmix%conv_kd = KD_CONV

         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ_BL
            ms%rho_layer(:, :, k) = 1028.0_wp + 1.0_wp*real(k - 1, wp)
         end do
         vmix%bl_depth = BLD

         call run_conv(grid, ms, vmix, vmix%bl_depth)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2

         ! k = 10..7 : z_int = 20,40,60,80 < bld -> untouched (PP81 value).
         do k = 7, NZ_BL
            call check(error, &
                       abs(vmix%kt(i_probe, j_probe, k) - PP81_UNSTABLE_KT) < TOL, &
                       "below_bl: interface inside BL was touched by convection")
            if (allocated(error)) exit checks
         end do

         ! k = 6..2 : z_int = 100,120,...,180 >= bld -> kd_conv.
         do k = 2, 6
            call check(error, abs(vmix%kt(i_probe, j_probe, k) - KD_CONV) < TOL, &
                       "below_bl: interface below BL was NOT set to kd_conv")
            if (allocated(error)) exit checks
         end do

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_below_bl_only

   subroutine test_n2_thresh_live(error)
      !! Weakly stable column tuned to N^2 = 5e-6 exactly at every
      !! interior interface (uniform gradient).  n2_thresh=0 and
      !! n2_thresh=1e-6 must NOT trigger (5e-6 not < either);ensures
      !! the strict "<" plus a live magnitude comparison.
      !! n2_thresh=1e-5 MUST trigger (5e-6 < 1e-5).  Proves n2_thresh
      !! actually reaches the kernel (not a dead knob) and that the
      !! sign/magnitude of the N^2 comparison is correct.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: N2_TARGET = 5.0e-6_wp
      real(wp), parameter :: KD_CONV = 1.0_wp
      real(wp), parameter :: PP81_BG_KT = 1.0e-5_wp
      real(wp) :: drho, rho0_val
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%conv_enable = .true.
         vmix%conv_kd = KD_CONV
         rho0_val = vmix%rho0
         ! N2 = -GRAVITY*(rho_k - rho_km1)/(rho0*dz_face); dz_face=H_LAYER
         ! (constant layers) => drho = -N2*rho0*H_LAYER/GRAVITY, applied
         ! at every interior interface via a constant per-layer slope.
         drho = N2_TARGET*rho0_val*H_LAYER/9.80665_wp

         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp - drho*real(k - 1, wp)
         end do
         vmix%bl_depth = 0.0_wp
         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2

         ! n2_thresh = 0: 5e-6 not < 0 -> untouched.
         vmix%conv_n2_thresh = 0.0_wp
         call run_conv(grid, ms, vmix, vmix%bl_depth)
         call check(error, abs(vmix%kt(i_probe, j_probe, 2) - PP81_BG_KT) < 1.0e-9_wp, &
                    "n2_thresh=0: interface wrongly triggered")
         if (allocated(error)) exit checks

         ! n2_thresh = 1e-6: 5e-6 not < 1e-6 -> untouched.
         vmix%conv_n2_thresh = 1.0e-6_wp
         call run_conv(grid, ms, vmix, vmix%bl_depth)
         call check(error, abs(vmix%kt(i_probe, j_probe, 2) - PP81_BG_KT) < 1.0e-9_wp, &
                    "n2_thresh=1e-6: interface wrongly triggered")
         if (allocated(error)) exit checks

         ! n2_thresh = 1e-5: 5e-6 < 1e-5 -> triggers.
         vmix%conv_n2_thresh = 1.0e-5_wp
         call run_conv(grid, ms, vmix, vmix%bl_depth)
         call check(error, abs(vmix%kt(i_probe, j_probe, 2) - KD_CONV) < 1.0e-9_wp, &
                    "n2_thresh=1e-5: interface did not trigger")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_n2_thresh_live

   subroutine test_homogenises_inversion(error)
      !! Analytical: a uniformly-unstable linear column (rho ramps
      !! 1028 -> 1030 over 20 layers of 10 m, H=200 m) is integrated
      !! through `vdiff_apply_tracers` for 7 steps of dt=1200s using
      !! the vmix snapshot from ONE `run_conv` call (kt held fixed
      !! across the integration, matching the per-stage contributor
      !! contract -- PP81 rewrites kv/kt fresh each stage in production;
      !! this test isolates the diffusive effect of a single stage's
      !! kt).  Asserts:
      !!   (a) the inversion is substantially removed with convection on
      !!       (kd_conv=1) -- column spread collapses from 2.0 to well
      !!       under half;
      !!   (b) column-integrated hT is conserved (backward-Euler vdiff
      !!       is conservative -- convection must not create/destroy
      !!       heat);
      !!   (c) the discriminator: the SAME column/integration with
      !!       convection OFF (PP81's ~1.001e-2 only) retains MOST of
      !!       the inversion -- proving this is not something PP81's
      !!       clip alone can do.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_on, ms_off
      type(ocean_vmix_t) :: vmix_on, vmix_off
      type(ocean_vdiff_t) :: vd_on, vd_off
      integer, parameter :: NZ_H = 20
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 1200.0_wp
      integer, parameter :: NSTEPS = 7
      real(wp) :: T0(NZ_H), spread_on, spread_off, total0, total_on, total_off
      integer :: i_probe, j_probe, k, step, nx, ny
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         nx = grid%nx_total; ny = grid%ny_total

         ms_on%nz_ml = NZ_H; call ms_on%init(grid)
         ms_off%nz_ml = NZ_H; call ms_off%init(grid)
         call vmix_on%init(grid, nz_ml=NZ_H)
         call vmix_off%init(grid, nz_ml=NZ_H)
         call vd_on%init(grid, nz_ml=NZ_H)
         call vd_off%init(grid, nz_ml=NZ_H)

         vmix_on%conv_enable = .true.
         vmix_on%conv_kd = 1.0_wp
         vmix_off%conv_enable = .false.

         ms_on%h_layer = H_LAYER; ms_off%h_layer = H_LAYER
         ms_on%u_face_x_layer = 0.0_wp; ms_off%u_face_x_layer = 0.0_wp
         ms_on%v_face_y_layer = 0.0_wp; ms_off%v_face_y_layer = 0.0_wp
         do k = 1, NZ_H
            T0(k) = 1028.0_wp + 2.0_wp*real(k - 1, wp)/real(NZ_H - 1, wp)
            ms_on%rho_layer(:, :, k) = T0(k)
            ms_off%rho_layer(:, :, k) = T0(k)
            ms_on%tracers(ms_on%idx_temperature)%hTr(:, :, k) = T0(k)*H_LAYER
            ms_off%tracers(ms_off%idx_temperature)%hTr(:, :, k) = T0(k)*H_LAYER
         end do
         vmix_on%bl_depth = 0.0_wp
         vmix_off%bl_depth = 0.0_wp

         ! One vmix snapshot each (kt held fixed across the integration).
         call run_conv(grid, ms_on, vmix_on, vmix_on%bl_depth)
         call run_no_conv(grid, ms_off, vmix_off)

         ! Integrate: map, step NSTEPS times with the fixed kt, pull back.
         !$acc enter data copyin(ms_on, vmix_on, vd_on)
         call ms_on%enter_data(); call vmix_on%enter_data(); call vd_on%enter_data()
         !$acc update device(vmix_on%kt)
         do step = 1, NSTEPS
            call vdiff_apply_tracers(grid, vd_on, ms_on, DT, kt_source=vmix_on%kt)
         end do
         !$acc update self(ms_on%tracers(ms_on%idx_temperature)%hTr)
         call vd_on%exit_data(); call vmix_on%exit_data(); call ms_on%exit_data()
         !$acc exit data delete(ms_on, vmix_on, vd_on)

         !$acc enter data copyin(ms_off, vmix_off, vd_off)
         call ms_off%enter_data(); call vmix_off%enter_data(); call vd_off%enter_data()
         !$acc update device(vmix_off%kt)
         do step = 1, NSTEPS
            call vdiff_apply_tracers(grid, vd_off, ms_off, DT, kt_source=vmix_off%kt)
         end do
         !$acc update self(ms_off%tracers(ms_off%idx_temperature)%hTr)
         call vd_off%exit_data(); call vmix_off%exit_data(); call ms_off%exit_data()
         !$acc exit data delete(ms_off, vmix_off, vd_off)

         i_probe = nx/2; j_probe = ny/2

         spread_on = maxval(ms_on%tracers(ms_on%idx_temperature)%hTr(i_probe, j_probe, :)/H_LAYER) - &
                     minval(ms_on%tracers(ms_on%idx_temperature)%hTr(i_probe, j_probe, :)/H_LAYER)
         spread_off = maxval(ms_off%tracers(ms_off%idx_temperature)%hTr(i_probe, j_probe, :)/H_LAYER) - &
                      minval(ms_off%tracers(ms_off%idx_temperature)%hTr(i_probe, j_probe, :)/H_LAYER)

         ! (a) convection substantially removes the inversion.
         call check(error, spread_on < 0.6_wp, &
                    "homogenises: convection-on spread did not collapse")
         if (allocated(error)) exit checks

         ! (c) the PP81-only control retains most of it (~100x slower).
         call check(error, spread_off > 1.5_wp, &
                    "homogenises: PP81-only control unexpectedly mixed away")
         if (allocated(error)) exit checks

         ! (b) column-integrated hT conserved (backward-Euler vdiff is
         ! conservative; convection must not create/destroy heat).
         total0 = sum(T0)*H_LAYER
         total_on = sum(ms_on%tracers(ms_on%idx_temperature)%hTr(i_probe, j_probe, :))
         total_off = sum(ms_off%tracers(ms_off%idx_temperature)%hTr(i_probe, j_probe, :))

         call check(error, abs(total_on - total0) < 1.0e-9_wp*abs(total0), &
                    "homogenises: conv-on hT not conserved")
         if (allocated(error)) exit checks
         call check(error, abs(total_off - total0) < 1.0e-9_wp*abs(total0), &
                    "homogenises: conv-off hT not conserved")

      end block checks
      call vd_on%destroy(); call vd_off%destroy()
      call vmix_on%destroy(); call vmix_off%destroy()
      call ms_on%destroy(); call ms_off%destroy()
   end subroutine test_homogenises_inversion

   subroutine test_vanished_interface_skipped(error)
      !! Two adjacent layers (k=2,3) pinched to 1e-6 m (< H_VANISHED =
      !! 1.5e-4 m) so their shared interface (k=3) has
      !! dz_face = 0.5*(1e-6+1e-6) = 1e-6 <= H_VANISHED, with a strong
      !! density inversion across it (rho jumps from 1028 to 1040).
      !! Pinching a SINGLE layer does not reproduce the gate (dz_face is
      !! an average, so one thin layer next to a normal 10 m neighbour
      !! still averages to ~5 m) -- two adjacent thin layers are the
      !! minimal reproduction of a genuinely vanished interface under
      !! this averaging convention.
      !!
      !! Asserts kt at that interface is FINITE and UNCHANGED from
      !! PP81's own value there (PP81 itself does not blow up at
      !! dz_face=1e-6 -- only exactly-zero trips its own `<=0.0` guard --
      !! so this isolates the convection-specific H_VANISHED gate: with
      !! it removed, convection would (wrongly) classify this
      !! numerically-degenerate interface as convecting and set
      !! kt = kd_conv instead of leaving it at the PP81 value).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: H_PINCH = 1.0e-6_wp
      real(wp), parameter :: KD_CONV = 1.0_wp
      real(wp), parameter :: PP81_UNSTABLE_KT = 1.0e-5_wp + 1.0e-2_wp
      integer :: i_probe, j_probe
      real(wp) :: kt_val
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%conv_enable = .true.
         vmix%conv_kd = KD_CONV

         ms%h_layer = H_LAYER
         ms%h_layer(:, :, 2) = H_PINCH
         ms%h_layer(:, :, 3) = H_PINCH
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%rho_layer(:, :, 1) = 1028.0_wp
         ms%rho_layer(:, :, 2) = 1028.0_wp
         ms%rho_layer(:, :, 3) = 1040.0_wp
         ms%rho_layer(:, :, 4) = 1040.0_wp
         vmix%bl_depth = 0.0_wp

         call run_conv(grid, ms, vmix, vmix%bl_depth)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         kt_val = vmix%kt(i_probe, j_probe, 3)

         call check(error, kt_val == kt_val, &
                    "vanished interface: kt is NaN")
         if (allocated(error)) exit checks
         call check(error, abs(kt_val) < huge(1.0_wp), &
                    "vanished interface: kt is Inf")
         if (allocated(error)) exit checks
         call check(error, abs(kt_val - PP81_UNSTABLE_KT) < 1.0e-9_wp, &
                    "vanished interface: convection touched a pinched interface")
         if (allocated(error)) exit checks
         ! Discriminator: an ungated kernel would have set kd_conv here.
         call check(error, abs(kt_val - KD_CONV) > 0.5_wp, &
                    "vanished interface: kt looks like it was set to kd_conv")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_vanished_interface_skipped

end module test_ocean_convection
