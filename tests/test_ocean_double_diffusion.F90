!! Unit tests for double diffusion (`vmix_split_ddiff_impl`, dispatched by
!! `vmix_split_kd_heat_salt` when `vmix%ddiff_enable`).
!!
!! Double diffusion is NOT a kv/kt contributor -- it is folded INTO the
!! heat/salt split, producing an ASYMMETRIC divergence between the salt
!! diffusivity `ks` and the temperature diffusivity `kt`:
!!     ks = kt_pre + kd_extra_s ;  kt = kt_pre + kd_extra_t
!! both from the SAME pre-split kt.  Interior interfaces only; branched on
!! the SIGNED alpha*dT / beta*dS (never a pre-divided R_rho).  CVMix
!! `cvmix_coeffs_ddiff` algebra: Large, McWilliams & Doney (1994) fingering
!! + Marmorino-Caldwell (1976) / Kelley (1990) diffusive convection.
!!
!! Every case isolates the double-diffusion contribution by setting the
!! pre-split `kt` to zero on the device (so `ks == kd_extra_s`,
!! `kt == kd_extra_t`) and driving `vmix_split_kd_heat_salt` directly, via
!! the `run_ddiff` harness (mem:separate mapping mirrors `test_ocean_pp81`
!! / `test_ocean_convection`).  Column T/S are set as thickness-integrals
!! `hTr = T*h_layer` on the two default registry tracers.
!!
!! Golden K_S / K_T values are the closed CVMix forms, cross-checked in
!! Python (docs/ocean_double_diffusion_plan.md).  Bottom-up convention:
!! interface k is at the bottom of layer k, upper = layer k, lower = k-1.
!!
!! Cases:
!!   * dd_disabled_bit_identical -- default-off: ks := kt exactly.
!!   * dd_fingering_exact        -- warm/salty over cool/fresh (Rrho=1.3):
!!                                  K_S > K_T > 0, exact; K_T = 0.7*K_S.
!!   * dd_convection_exact       -- cool/fresh over warm/salty (Rrho=0.7):
!!                                  K_T > K_S > 0, exact MC76 + flux ratio.
!!   * dd_above_cutoff_zero      -- Rrho >= strat_param_max: both zero.
!!   * dd_doubly_stable_untouched-- neither regime fires: ks == kt.
!!   * dd_boundary_interfaces_zero-- k=1 (bed) / k=nz+1 (surface) stay zero.
module test_ocean_double_diffusion
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_split_kd_heat_salt
   implicit none
   private

   public :: collect_ocean_double_diffusion_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: H_LAYER = 10.0_wp

   ! EOS linear coefficients (eos_t defaults) -- the constant
   ! alpha_T / beta_S the v1 kernel uses.
   real(wp), parameter :: ALPHA_T = 1.7e-4_wp
   real(wp), parameter :: BETA_S = 7.6e-4_wp

   ! Golden diffusivities (Python-verified, docs plan).
   real(wp), parameter :: K_S_FINGER = 5.2448726125e-5_wp
   real(wp), parameter :: K_T_FINGER = 3.6714108288e-5_wp   ! = 0.7*K_S
   real(wp), parameter :: K_T_CONV = 5.2441171668e-5_wp
   real(wp), parameter :: K_S_CONV = 2.3336321392e-5_wp     ! = 0.445*K_T

   real(wp), parameter :: RTOL = 1.0e-7_wp   ! tolerates fast-math exp/pow

contains

   subroutine collect_ocean_double_diffusion_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dd_disabled_bit_identical", test_disabled_bit_identical), &
                  new_unittest("dd_fingering_exact", test_fingering_exact), &
                  new_unittest("dd_convection_exact", test_convection_exact), &
                  new_unittest("dd_above_cutoff_zero", test_above_cutoff_zero), &
                  new_unittest("dd_doubly_stable_untouched", test_doubly_stable_untouched), &
                  new_unittest("dd_boundary_interfaces_zero", test_boundary_interfaces_zero) &
                  ]
   end subroutine collect_ocean_double_diffusion_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine set_column(ms, tprof, sprof)
      !! Fill the temperature/salinity registry tracers with a per-layer
      !! profile (bottom-up: index 1 = bed, NZ = surface), as the
      !! thickness-integral hTr = T*h_layer.  Uniform h = H_LAYER.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: tprof(NZ), sprof(NZ)
      integer :: k
      ms%h_layer = H_LAYER
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = tprof(k)*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = sprof(k)*H_LAYER
      end do
   end subroutine set_column

   subroutine run_ddiff(grid, ms, vmix)
      !! Map ms + vmix, push the host-set pre-split kt/ks to the device,
      !! run the double-diffusion split, pull kt/ks back.  Caller sets
      !! `vmix%kt` / `vmix%ks` on the host BEFORE this call (mem:separate:
      !! the explicit `update device` carries them regardless of whether
      !! enter_data maps them copyin or create).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      !$acc update device(vmix%kt, vmix%ks)
      call vmix_split_kd_heat_salt(grid, vmix, ms)
      !$acc update self(vmix%kt, vmix%ks)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_ddiff

   logical function close_to(a, b)
      real(wp), intent(in) :: a, b
      close_to = abs(a - b) <= RTOL*abs(b) + 1.0e-18_wp
   end function close_to

   ! -----------------------------------------------------------------

   subroutine test_disabled_bit_identical(error)
      !! `ddiff_enable = .false.` => the split reduces to `ks := kt`
      !! exactly, for an arbitrary (nonzero, non-uniform) pre-split kt --
      !! and on a column that WOULD trigger fingering, so this is not
      !! vacuously true.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer :: k
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .false.
         call set_column(ms, [5.0_wp, 7.0_wp, 9.0_wp], &
                         [34.0_wp, 34.3441295547_wp, 34.6882591094_wp])
         vmix%ks = 0.0_wp
         do k = 1, NZ + 1
            vmix%kt(:, :, k) = 1.0e-3_wp*real(k, wp)   ! arbitrary nonzero pattern
         end do

         call run_ddiff(grid, ms, vmix)

         call check(error, maxval(abs(vmix%ks - vmix%kt)) < 1.0e-15_wp, &
                    "disabled: ks must equal kt bit-for-bit")
         if (allocated(error)) exit checks
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_disabled_bit_identical

   subroutine test_fingering_exact(error)
      !! Warm/salty over cool/fresh, uniform interface Rrho = 1.3.  Both
      !! interior interfaces (k=2,3) must hit the exact clamped-cubic
      !! K_S and K_T = 0.7*K_S; K_S > K_T (fingering enhances salt).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer :: ip, jp, k
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .true.
         ! T warmer / S saltier toward the surface (index NZ).
         call set_column(ms, [5.0_wp, 7.0_wp, 9.0_wp], &
                         [34.0_wp, 34.3441295547_wp, 34.6882591094_wp])
         vmix%kt = 0.0_wp; vmix%ks = 0.0_wp

         call run_ddiff(grid, ms, vmix)

         ip = grid%nx_total/2; jp = grid%ny_total/2
         do k = 2, NZ
            call check(error, close_to(vmix%ks(ip, jp, k), K_S_FINGER), &
                       "fingering: K_S wrong")
            if (allocated(error)) exit checks
            call check(error, close_to(vmix%kt(ip, jp, k), K_T_FINGER), &
                       "fingering: K_T wrong (must be 0.7*K_S)")
            if (allocated(error)) exit checks
            call check(error, vmix%ks(ip, jp, k) > vmix%kt(ip, jp, k), &
                       "fingering: K_S must exceed K_T")
            if (allocated(error)) exit checks
         end do
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_fingering_exact

   subroutine test_convection_exact(error)
      !! Cool/fresh over warm/salty, uniform interface Rrho = 0.7.  Both
      !! interior interfaces must hit the exact MC76 K_T and the flux-ratio
      !! K_S = (1.85*Rrho - 0.85)*K_T; K_T > K_S (convection enhances heat).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer :: ip, jp, k
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .true.
         ! T colder / S fresher toward the surface (index NZ).
         call set_column(ms, [9.0_wp, 7.0_wp, 5.0_wp], &
                         [35.0_wp, 34.3609022556_wp, 33.7218045112_wp])
         vmix%kt = 0.0_wp; vmix%ks = 0.0_wp

         call run_ddiff(grid, ms, vmix)

         ip = grid%nx_total/2; jp = grid%ny_total/2
         do k = 2, NZ
            call check(error, close_to(vmix%kt(ip, jp, k), K_T_CONV), &
                       "convection: K_T wrong")
            if (allocated(error)) exit checks
            call check(error, close_to(vmix%ks(ip, jp, k), K_S_CONV), &
                       "convection: K_S wrong (flux ratio)")
            if (allocated(error)) exit checks
            call check(error, vmix%kt(ip, jp, k) > vmix%ks(ip, jp, k), &
                       "convection: K_T must exceed K_S")
            if (allocated(error)) exit checks
         end do
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_convection_exact

   subroutine test_above_cutoff_zero(error)
      !! Fingering geometry but Rrho = 3.0 > strat_param_max (2.55): the
      !! clamped-cubic is zero, so K_S = K_T = 0 and ks == kt == 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp) :: dS
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .true.
         ! Rrho = (alpha*dT)/(beta*dS) = 3.0 with dT = 2.0.
         dS = ALPHA_T*2.0_wp/(BETA_S*3.0_wp)
         call set_column(ms, [5.0_wp, 7.0_wp, 9.0_wp], &
                         [34.0_wp, 34.0_wp + dS, 34.0_wp + 2.0_wp*dS])
         vmix%kt = 0.0_wp; vmix%ks = 0.0_wp

         call run_ddiff(grid, ms, vmix)

         call check(error, maxval(abs(vmix%ks)) < 1.0e-18_wp .and. &
                    maxval(abs(vmix%kt)) < 1.0e-18_wp, &
                    "above-cutoff: diffusivities must be exactly zero")
         if (allocated(error)) exit checks
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_above_cutoff_zero

   subroutine test_doubly_stable_untouched(error)
      !! Warm AND fresh toward the surface: alpha*dT > 0 but beta*dS < 0,
      !! so neither the fingering (needs beta*dS > 0) nor the convection
      !! (needs alpha*dT < 0) trigger fires.  ks == kt everywhere -- the
      !! ~most-of-the-ocean doubly-stable case.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .true.
         ! Warm (T up) + fresh (S down) toward the surface.
         call set_column(ms, [5.0_wp, 7.0_wp, 9.0_wp], &
                         [35.0_wp, 34.5_wp, 34.0_wp])
         vmix%kt = 0.0_wp; vmix%ks = 0.0_wp

         call run_ddiff(grid, ms, vmix)

         call check(error, maxval(abs(vmix%ks - vmix%kt)) < 1.0e-15_wp .and. &
                    maxval(abs(vmix%kt)) < 1.0e-18_wp, &
                    "doubly-stable: no regime fires, ks == kt == 0")
         if (allocated(error)) exit checks
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_doubly_stable_untouched

   subroutine test_boundary_interfaces_zero(error)
      !! Under an active fingering column, the bed (k=1) and surface
      !! (k=nz+1) interfaces carry no double-diffusion contribution and
      !! stay at the pre-split kt (here zero) -- the closed-BC invariant.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      checks: block
         call make_grid(grid)
         ms%nz_ml = NZ; call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         vmix%ddiff_enable = .true.
         call set_column(ms, [5.0_wp, 7.0_wp, 9.0_wp], &
                         [34.0_wp, 34.3441295547_wp, 34.6882591094_wp])
         vmix%kt = 0.0_wp; vmix%ks = 0.0_wp

         call run_ddiff(grid, ms, vmix)

         call check(error, maxval(abs(vmix%ks(:, :, 1))) < 1.0e-18_wp .and. &
                    maxval(abs(vmix%kt(:, :, 1))) < 1.0e-18_wp, &
                    "bed interface (k=1) must be exactly zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vmix%ks(:, :, NZ + 1))) < 1.0e-18_wp .and. &
                    maxval(abs(vmix%kt(:, :, NZ + 1))) < 1.0e-18_wp, &
                    "surface interface (k=nz+1) must be exactly zero")
         if (allocated(error)) exit checks
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_boundary_interfaces_zero

end module test_ocean_double_diffusion
