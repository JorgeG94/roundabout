!! Unit tests for `vmix_add_kv_ml_invz2`: MOM6's KV_ML_INVZ2 + HMIX_FIXED
!! surface-band viscosity augmentation.
!!
!! Cases:
!!   * No-op when `kv_ml_invz2 = 0` — `kv` field untouched.
!!   * Single interior interface inside HMIX — `kv` += predictable
!!     amount based on the 1/z² profile.
!!   * Interface below HMIX — `kv` unchanged.
!!   * Multi-layer column — profile decreases monotonically with depth
!!     and reaches the floor at z = HMIX.
module test_ocean_kv_ml_invz2
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_add_kv_ml_invz2
   implicit none
   private

   public :: collect_ocean_kv_ml_invz2_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 4
   integer, parameter :: NY_PHYS = 4

contains

   subroutine collect_ocean_kv_ml_invz2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zero_amplitude_is_noop", test_zero_amplitude), &
                  new_unittest("interface_below_hmix_unchanged", &
                               test_interface_below_hmix), &
                  new_unittest("profile_decreases_with_depth", &
                               test_profile_monotone), &
                  new_unittest("value_at_z_equals_hmix_is_baseline", &
                               test_value_at_hmix) &
                  ]
   end subroutine collect_ocean_kv_ml_invz2_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1000.0_wp, 1000.0_wp)
   end subroutine make_grid

   subroutine setup(grid, ms, vmix, nz, h_each)
      !! Build a uniform-h column.  `kv` initialised to zero so the
      !! augmentation reads cleanly.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_each

      call make_grid(grid)
      ms%nz_ml = nz
      call ms%init(grid)
      call vmix%init(grid, nz_ml=nz)
      ms%h_layer = h_each
      vmix%kv = 0.0_wp
   end subroutine setup

   subroutine test_zero_amplitude(error)
      !! `kv_ml_invz2 = 0` must leave `kv` untouched, regardless of
      !! `hmix_fixed` or the layer thicknesses.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp) :: max_kv

      call setup(grid, ms, vmix, nz=4, h_each=5.0_wp)
      vmix%kv_ml_invz2 = 0.0_wp
      vmix%hmix_fixed = 20.0_wp
      call vmix_add_kv_ml_invz2(grid, vmix, ms)
      max_kv = maxval(abs(vmix%kv))
      call check(error, max_kv < 1.0e-15_wp, &
                 "kv_ml_invz2 = 0 should be a no-op")
      call vmix%destroy(); call ms%destroy()
   end subroutine test_zero_amplitude

   subroutine test_interface_below_hmix(error)
      !! Single 100 m layer, hmix_fixed = 20 m.  The topmost interior
      !! interface (k = nz) sits at depth 100 m below the surface —
      !! outside HMIX.  The augmentation must skip it; `kv` stays at 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp) :: max_kv

      call setup(grid, ms, vmix, nz=2, h_each=100.0_wp)
      vmix%kv_ml_invz2 = 0.01_wp
      vmix%hmix_fixed = 20.0_wp
      call vmix_add_kv_ml_invz2(grid, vmix, ms)
      max_kv = maxval(abs(vmix%kv))
      call check(error, max_kv < 1.0e-15_wp, &
                 "Interfaces below HMIX should not be augmented")
      call vmix%destroy(); call ms%destroy()
   end subroutine test_interface_below_hmix

   subroutine test_profile_monotone(error)
      !! Thin layers (h_each = 1 m, nz = 50).  Compare `kv` at two
      !! interior interfaces inside HMIX: shallow z = 1 m and
      !! deeper z = 10 m.  Profile `(hmix/z)²` must give shallow > deep.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: KV0 = 0.01_wp, HMIX = 20.0_wp, HK = 1.0_wp
      integer, parameter :: NZ = 50
      integer :: nx, ny, nz_eff, k_shallow, k_deep, ig, jg
      real(wp) :: kv_shallow, kv_deep

      call setup(grid, ms, vmix, nz=NZ, h_each=HK)
      vmix%kv_ml_invz2 = KV0
      vmix%hmix_fixed = HMIX
      call vmix_add_kv_ml_invz2(grid, vmix, ms)
      nx = grid%nx_total
      ny = grid%ny_total
      nz_eff = NZ
      ig = nx/2
      jg = ny/2
      ! Surface-most interior interface = k = nz (depth = HK)
      k_shallow = nz_eff
      ! Deeper interface = k = nz - 9 (depth = 10 m)
      k_deep = nz_eff - 9
      kv_shallow = vmix%kv(ig, jg, k_shallow)
      kv_deep = vmix%kv(ig, jg, k_deep)
      call check(error, kv_shallow > kv_deep, &
                 "Profile not monotone-decreasing with depth")
      if (.not. allocated(error)) call check(error, kv_shallow > 0.0_wp .and. kv_deep > 0.0_wp, &
                                             "Both interfaces inside HMIX should be non-zero")
      call vmix%destroy(); call ms%destroy()
   end subroutine test_profile_monotone

   subroutine test_value_at_hmix(error)
      !! At z = HMIX the profile term `(hmix/z)² = 1`, so the
      !! augmentation should land at exactly `kv_ml_invz2`.  Set up
      !! layers so an interface sits at depth = HMIX.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: KV0 = 0.01_wp, HMIX = 20.0_wp, HK = 5.0_wp
      integer, parameter :: NZ = 10
      integer :: nx, ny, ig, jg, k_at_hmix
      real(wp) :: kv_obs

      ! 10 layers of 5 m each.  Interfaces at depths 5, 10, 15, 20, ...
      ! Interface at z = 20 m is k = nz - 3 = 7.  The do-loop in
      ! vmix_add_kv_ml_invz2 exits when `z >= hmix`, so the k whose
      ! cumulative-z hits HMIX exactly is NOT augmented (we stop
      ! before processing it).  Inspect the previous interface
      ! (k = nz - 2 = 8, depth = 15 m) instead; expected value is
      ! KV0 * (HMIX/15)² = 0.01 * (20/15)² = 0.01778.
      call setup(grid, ms, vmix, nz=NZ, h_each=HK)
      vmix%kv_ml_invz2 = KV0
      vmix%hmix_fixed = HMIX
      call vmix_add_kv_ml_invz2(grid, vmix, ms)
      nx = grid%nx_total
      ny = grid%ny_total
      ig = nx/2
      jg = ny/2
      k_at_hmix = NZ - 2  ! interface at depth = 3*HK = 15 m
      kv_obs = vmix%kv(ig, jg, k_at_hmix)
      call check(error, &
                 abs(kv_obs - KV0*(HMIX/15.0_wp)**2) < 1.0e-14_wp, &
                 "kv augmentation off the (hmix/z)² profile")
      call vmix%destroy(); call ms%destroy()
   end subroutine test_value_at_hmix

end module test_ocean_kv_ml_invz2
