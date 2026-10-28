module test_ocean_epbl_max_shear
   !! `&ocean_epbl_nml combine = "max_shear"`: kappa-shear and ePBL combine by
   !! max with the background added outside,
   !!   kt = bg + max(kd_shear, kd_epbl),
   !!   kv = kv_bg + max(prandtl_turb*kd_shear, prandtl_epbl*kd_epbl),
   !! unlike "max" (max(bg, ePBL) + kappa-shear) and "add" (bg + shear + ePBL).
   !!
   !! One column, one 0.2 m/s velocity jump (kappa-shear peaks at the jump
   !! interface B, decays away from it).  ePBL kd is imposed: far above the
   !! local shear at interface A, far below it at B.  Production order on this
   !! branch: seed background -> kappa-shear merge -> max_shear ePBL merge.
   !!   * A: kt == bg + kd_epbl, kv == kv_bg + Pr_epbl*kd_epbl
   !!   * B: kt == bg + kd_shear, kv == kv_bg + prandtl_turb*kd_shear
   !!   (1e-15 relative; Pr_epbl is deliberately != prandtl_turb)
   !! and "add" differs at both interfaces and "max" at A (they coincide at B), so the
   !! assertions cannot pass by accident.  With ePBL's kd identical to a
   !! kappa-shear-free column, "max_shear" is bit-identical to "add".
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_seed_interior_none, VMIX_INTERIOR_NONE
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t, kappa_shear_compute, &
                                    kappa_shear_merge_into_kv_kt
   use rdb_ocean_epbl, only: ocean_epbl_t, epbl_merge_into_kv_kt, epbl_merge_max_shear, &
                             parse_epbl_combine, EPBL_COMBINE_ADD, EPBL_COMBINE_MAX, &
                             EPBL_COMBINE_MAX_SHEAR
   use rdb_eos, only: EOS_VARIANT_LINEAR
   implicit none
   private

   public :: collect_ocean_epbl_max_shear_tests

   integer, parameter :: NGHOST = 2, NZK = 12, KJ = 6
   integer, parameter :: KA = 3, KB = KJ + 1
   real(wp), parameter :: DZ = 5.0_wp, N2 = 3.0e-5_wp, DU = 0.2_wp
   real(wp), parameter :: KT_BG = 3.0e-5_wp, KV_BG = 2.0e-4_wp
   real(wp), parameter :: PR_EPBL = 1.5_wp
   real(wp), parameter :: KD_A = 0.05_wp, KD_B = 1.0e-4_wp

contains

   subroutine collect_ocean_epbl_max_shear_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("epbl_max_shear_two_interfaces", test_two_interfaces), &
                  new_unittest("epbl_max_shear_no_shear_is_add", test_no_shear_is_add), &
                  new_unittest("epbl_combine_parse", test_parse) &
                  ]
   end subroutine collect_ocean_epbl_max_shear_tests

   subroutine run_column(mode, with_shear, kv, kt, kd_shear, prandtl_turb)
      !! Seed -> [kappa-shear merge] -> ePBL merge of the requested mode, on
      !! the device.  Returns kv/kt and the kappa-shear kd_int at the probe.
      integer, intent(in) :: mode
      logical, intent(in) :: with_shear
      real(wp), intent(out) :: kv(NZK + 1), kt(NZK + 1), kd_shear(NZK + 1), prandtl_turb
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_kappa_shear_t) :: ks
      type(ocean_epbl_t) :: ep
      real(wp), parameter :: RHO0 = 1025.0_wp, ALPHA = RHO0*2.0e-4_wp, DT = 1800.0_wp
      real(wp) :: dtdz
      integer :: k, i, j

      call grid%init(4, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZK
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZK)
      vmix%use_closure = .true.
      vmix%interior_closure = VMIX_INTERIOR_NONE
      vmix%pp81_nu_bg = KV_BG
      vmix%pp81_kappa_bg = KT_BG
      call vmix%seed_backgrounds()
      call ks%init(grid, nz_ml=NZK)
      ks%enable = .true.
      ks%eos%variant = EOS_VARIANT_LINEAR
      ks%eos%alpha_T = ALPHA
      ks%eos%beta_S = RHO0*7.6e-4_wp
      ks%eos%rho0 = RHO0
      ks%rho0 = RHO0
      call ks%set_f_centre(grid, 1.0e-4_wp, 0.0_wp, 0.0_wp)
      call ep%init(grid, nz_ml=NZK)
      ep%prandtl = PR_EPBL
      ep%combine_mode = mode
      ep%kd_int = 0.0_wp
      ep%kd_int(:, :, KA) = KD_A
      ep%kd_int(:, :, KB) = KD_B

      dtdz = N2*RHO0/(GRAVITY*ALPHA)
      ms%h_layer = DZ
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZK
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (20.0_wp + dtdz*DZ*real(k - NZK, wp))*DZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*DZ
         if (k > KJ) ms%u_face_x_layer(:, :, k) = DU
      end do

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vmix)
      call vmix%enter_data()
      !$acc enter data copyin(ks)
      call ks%enter_data()
      !$acc enter data copyin(ep)
      call ep%enter_data()
      call vmix_seed_interior_none(grid, vmix, ms)
      if (with_shear) then
         call kappa_shear_compute(grid, ks, ms, DT)
         call kappa_shear_merge_into_kv_kt(ks, grid%nx_total, grid%ny_total, NZK + 1, vmix%kv, vmix%kt)
      end if
      if (mode == EPBL_COMBINE_MAX_SHEAR .and. with_shear) then
         call epbl_merge_max_shear(ep, grid%nx_total, grid%ny_total, NZK + 1, vmix%kv, vmix%kt, &
                                   ks%kd_int, ks%prandtl_turb)
      else
         call epbl_merge_into_kv_kt(ep, grid%nx_total, grid%ny_total, NZK + 1, vmix%kv, vmix%kt)
      end if
      !$acc update self(vmix%kv, vmix%kt, ks%kd_int)
      call ep%exit_data()
      !$acc exit data delete(ep)
      call ks%exit_data()
      !$acc exit data delete(ks)
      call vmix%exit_data()
      !$acc exit data delete(vmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      i = grid%nx_total/2; j = grid%ny_total/2
      kv = vmix%kv(i, j, :)
      kt = vmix%kt(i, j, :)
      kd_shear = ks%kd_int(i, j, :)
      prandtl_turb = ks%prandtl_turb
      call ep%destroy(); call ks%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine run_column

   pure function close_rel(a, b) result(ok)
      real(wp), intent(in) :: a, b
      logical :: ok
      ok = abs(a - b) <= 1.0e-15_wp*abs(b)
   end function close_rel

   subroutine test_two_interfaces(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kv(NZK + 1), kt(NZK + 1), kd_s(NZK + 1), pt
      real(wp) :: kv_a(NZK + 1), kt_a(NZK + 1), kv_m(NZK + 1), kt_m(NZK + 1)
      checks: block
         call run_column(EPBL_COMBINE_MAX_SHEAR, .true., kv, kt, kd_s, pt)
         call check(error, kd_s(KA) < 0.2_wp*KD_A .and. kd_s(KB) > 100.0_wp*KD_B, &
                    "test setup: need shear < ePBL at A and shear >> ePBL at B")
         if (allocated(error)) exit checks
         call check(error, close_rel(kt(KA), KT_BG + KD_A), "A: kt must be bg + kd_epbl")
         if (allocated(error)) exit checks
         call check(error, close_rel(kv(KA), KV_BG + PR_EPBL*KD_A), "A: kv must be kv_bg + Pr*kd_epbl")
         if (allocated(error)) exit checks
         call check(error, close_rel(kt(KB), KT_BG + kd_s(KB)), "B: kt must be bg + kd_shear")
         if (allocated(error)) exit checks
         call check(error, close_rel(kv(KB), KV_BG + pt*kd_s(KB)), "B: kv must be kv_bg + kv_shear")
         if (allocated(error)) exit checks
         ! Discrimination: "add" and "max" disagree with max_shear at both interfaces.
         call run_column(EPBL_COMBINE_ADD, .true., kv_a, kt_a, kd_s, pt)
         call run_column(EPBL_COMBINE_MAX, .true., kv_m, kt_m, kd_s, pt)
         call check(error,.not. close_rel(kt_a(KA), kt(KA)) .and. .not. close_rel(kt_a(KB), kt(KB)), &
                    "add must differ from max_shear at both interfaces")
         if (allocated(error)) exit checks
         call check(error,.not. close_rel(kt_m(KA), kt(KA)), &
                    "max must differ from max_shear at A (it drops the background there)")
      end block checks
   end subroutine test_two_interfaces

   subroutine test_no_shear_is_add(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kv1(NZK + 1), kt1(NZK + 1), kv2(NZK + 1), kt2(NZK + 1), kd_s(NZK + 1), pt
      checks: block
         call run_column(EPBL_COMBINE_MAX_SHEAR, .false., kv1, kt1, kd_s, pt)
         call run_column(EPBL_COMBINE_ADD, .false., kv2, kt2, kd_s, pt)
         call check(error, all(kt1 == kt2) .and. all(kv1 == kv2), &
                    "max_shear without kappa-shear must equal add exactly")
         if (allocated(error)) exit checks
         call check(error, close_rel(kt1(KA), KT_BG + KD_A) .and. close_rel(kv1(KB), KV_BG + PR_EPBL*KD_B), &
                    "no-shear column must be bg + ePBL (test is not vacuous)")
      end block checks
   end subroutine test_no_shear_is_add

   subroutine test_parse(error)
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_epbl_combine("max_shear") == EPBL_COMBINE_MAX_SHEAR .and. &
                    parse_epbl_combine("add") == EPBL_COMBINE_ADD .and. &
                    parse_epbl_combine("max") == EPBL_COMBINE_MAX .and. &
                    parse_epbl_combine("max_shea") < 0, "combine parser")
      end block checks
   end subroutine test_parse

end module test_ocean_epbl_max_shear
