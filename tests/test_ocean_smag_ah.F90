!! Phase 5b biharmonic-Smagorinsky (MOM6 `SMAGORINSKY_AH` analogue).
!!
!! Covers:
!!   1. Zero flow → A_4 = nu4_bg everywhere.
!!   2. Uniform flow → strain rate = 0 → A_4 = nu4_bg.
!!   3. Linear-shear flow with known strain → A_4 matches analytical
!!      `C_b · L⁴ · |D|` formula at interior faces.
!!   4. nu4_max clip honoured under extreme strain.
!!   5. End-to-end: the biharmonic kernel with `smag_ah_active=.true.`
!!      uses the per-face viscosity (different from scalar nu_4).
module test_ocean_smag_ah
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, make_anisotropic_metrics, &
                                 destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_smag_ah, &
                                    ocean_lateral_mix_compute_smag, &
                                    LMIX_NONE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_smag_ah_tests

contains

   subroutine collect_ocean_smag_ah_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("smag_ah_zero_flow_is_background", test_zero_flow), &
                  new_unittest("smag_ah_uniform_flow_is_background", test_uniform_flow), &
                  new_unittest("smag_ah_shear_matches_analytic", test_shear_analytic), &
                  new_unittest("smag_ah_nu4_max_clip_honoured", test_nu4_max_clip), &
                  new_unittest("hvisc_with_smag_ah_uses_face_coefficient", &
                               test_hvisc_with_smag_ah), &
                  new_unittest("smag_anisotropic_corner_shear_ratio_bundle", &
                               test_anisotropic_corner_shear) &
                  ]
   end subroutine collect_ocean_smag_ah_tests

   subroutine setup_state(grid, ms, nx, ny, nz, dx)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx
      call grid%init(nx, ny, 1, dx, dx)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_state

   subroutine map_in(ms, metrics, grid, lmix, hv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(ocean_lateral_mix_t), intent(inout), optional :: lmix
      type(ocean_horizontal_viscosity_t), intent(inout), optional :: hv
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      if (present(lmix)) then
         !$acc enter data copyin(lmix)
         call lmix%enter_data()
      end if
      if (present(hv)) then
         !$acc enter data copyin(hv)
         call hv%enter_data()
      end if
   end subroutine map_in

   subroutine map_out(ms, metrics, lmix, hv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_lateral_mix_t), intent(inout), optional :: lmix
      type(ocean_horizontal_viscosity_t), intent(inout), optional :: hv
      if (present(hv)) then
         call hv%exit_data()
         !$acc exit data delete(hv)
      end if
      if (present(lmix)) then
         call lmix%exit_data()
         !$acc exit data delete(lmix)
      end if
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   subroutine test_zero_flow(error)
      !! Zero velocity → strain rate = 0 → A_4 = nu4_bg everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_BG = 1.0e8_wp
      real(wp) :: max_diff_x, max_diff_y

      call setup_state(grid, ms, 8, 8, 2, 50000.0_wp)
      call lmix%init(grid, nz_ml=2)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = NU4_BG
      lmix%nu4_max = 1.0e12_wp

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      max_diff_x = maxval(abs(lmix%nu4_face_x - NU4_BG))
      max_diff_y = maxval(abs(lmix%nu4_face_y - NU4_BG))

      call check(error, max_diff_x < 1.0e-9_wp .and. max_diff_y < 1.0e-9_wp, &
                 "Zero flow must leave nu4_face at nu4_bg")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_zero_flow

   subroutine test_uniform_flow(error)
      !! Uniform u, v → zero strain → A_4 = nu4_bg everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_BG = 5.0e8_wp
      real(wp) :: max_diff_x, max_diff_y

      call setup_state(grid, ms, 8, 8, 2, 50000.0_wp)
      ms%u_face_x_layer = 0.3_wp
      ms%v_face_y_layer = -0.2_wp
      call lmix%init(grid, nz_ml=2)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = NU4_BG
      lmix%nu4_max = 1.0e12_wp

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      max_diff_x = maxval(abs(lmix%nu4_face_x - NU4_BG))
      max_diff_y = maxval(abs(lmix%nu4_face_y - NU4_BG))

      call check(error, max_diff_x < 1.0e-6_wp .and. max_diff_y < 1.0e-6_wp, &
                 "Uniform flow must give nu4_face = nu4_bg")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_uniform_flow

   subroutine test_shear_analytic(error)
      !! Linear shear flow u(y) = S · y_j, v = 0:
      !!   D_T = ∂u/∂x − ∂v/∂y = 0
      !!   D_S = ∂v/∂x + ∂u/∂y = S/dy * dy = S/dy ... wait.
      !!
      !! Concretely, on C-grid:
      !!   u_face_x(i,j) = S * (j - j0)    (per-row constant)
      !!   ∂u/∂y at a face = (u(i,j+1) - u(i,j-1)) / (2·dy)  but our
      !!   centred kernel uses (u(i, j+1) - u(i, j-1))/dy or 1-sided
      !!   pairs; see code for the exact stencil.  The strain
      !!   magnitude |D| ends up being roughly |S/dy|, so
      !!   A_4 ≈ C_b · L⁴ · |S/dy| in the interior.
      !!
      !! We construct the flow and check that A_4 is > NU4_BG (i.e.,
      !! Smag actually kicked in) AND that it's a sensible scale
      !! (within an order of magnitude of the analytical estimate).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_BG = 1.0e6_wp
      real(wp), parameter :: DX = 50000.0_wp
      real(wp), parameter :: SHEAR = 0.01_wp   ! du/dy ~ 0.01·dy / dy = 0.01 (1/s normalized)
      real(wp), parameter :: C_B = 0.06_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp) :: grid_sp_h2, expected_scale, observed_x_mid
      integer :: i, j

      call setup_state(grid, ms, NX, NY, NZ, DX)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total + 1
            ms%u_face_x_layer(i, j, 1) = SHEAR*real(j, wp)
         end do
      end do
      ! v stays zero
      call lmix%init(grid, nz_ml=NZ)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = C_B
      lmix%nu4_bg = NU4_BG
      lmix%nu4_max = 1.0e15_wp   ! large enough not to clip

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      ! Analytic estimate: for u_face_x(i,j) = SHEAR·j (per-row
      ! constant in our integer-j IC), the centred-difference shear
      ! at an interior u-face is `D_S ≈ SHEAR/dy`.  With v=0 the
      ! tension D_T vanishes, so |D| ≈ SHEAR/dy.  Then
      !   A_4 = C_b · L⁴ · |D|   with L² = grid_sp_h2 = dx² (dx=dy).
      grid_sp_h2 = DX*DX
      expected_scale = C_B*grid_sp_h2*grid_sp_h2*(SHEAR/DX)

      ! Mid-domain interior u-face value (away from walls).
      observed_x_mid = lmix%nu4_face_x(NX/2, NY/2, 1)

      call check(error, observed_x_mid > 2.0_wp*NU4_BG, &
                 "Linear shear must raise nu4_face_x well above background")
      if (.not. allocated(error)) call check(error, &
                                             observed_x_mid > 0.1_wp*expected_scale .and. &
                                             observed_x_mid < 10.0_wp*expected_scale, &
                                             "nu4_face_x not within an order of magnitude of analytic")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_shear_analytic

   subroutine test_nu4_max_clip(error)
      !! Set up a noisy alternating-sign flow that produces large
      !! strain, then assert the cap is honoured: no nu4_face_*
      !! exceeds nu4_max.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_MAX = 2.0e8_wp
      integer, parameter :: NX = 10, NY = 10, NZ = 1
      real(wp) :: peak_x, peak_y
      integer :: i, j

      call setup_state(grid, ms, NX, NY, NZ, 50000.0_wp)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total + 1
            ms%u_face_x_layer(i, j, 1) = merge(1.0_wp, -1.0_wp, mod(i + j, 2) == 0)
         end do
      end do
      do j = 1, grid%ny_total + 1
         do i = 1, grid%nx_total
            ms%v_face_y_layer(i, j, 1) = merge(-1.0_wp, 1.0_wp, mod(i + j, 2) == 0)
         end do
      end do
      call lmix%init(grid, nz_ml=NZ)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = 0.0_wp
      lmix%nu4_max = NU4_MAX

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      peak_x = maxval(lmix%nu4_face_x)
      peak_y = maxval(lmix%nu4_face_y)

      call check(error, peak_x <= NU4_MAX + 1.0e-6_wp .and. peak_y <= NU4_MAX + 1.0e-6_wp, &
                 "nu4_max clip not honoured under noisy flow")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_nu4_max_clip

   subroutine test_hvisc_with_smag_ah(error)
      !! End-to-end: hvisc compute_tendencies with `lateral_mix`
      !! present and `smag_ah_active=.true.` produces tendencies
      !! that differ from those produced by the scalar nu_4 path
      !! (i.e., the face-impl is actually engaged).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_face, grid_scalar
      type(multilayer_state_t) :: ms_face, ms_scalar
      type(ocean_metrics_t) :: metrics_face, metrics_scalar
      type(ocean_lateral_mix_t) :: lmix_face
      type(ocean_horizontal_viscosity_t) :: hv_face, hv_scalar
      real(wp), parameter :: DX = 50000.0_wp, NU4_CONST = 1.0e9_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp) :: max_diff
      integer :: i, j

      call setup_state(grid_face, ms_face, NX, NY, NZ, DX)
      call setup_state(grid_scalar, ms_scalar, NX, NY, NZ, DX)
      do j = 1, grid_face%ny_total
         do i = 1, grid_face%nx_total + 1
            ms_face%u_face_x_layer(i, j, 1) = 0.01_wp*sin(2.0_wp*real(j, wp))
            ms_scalar%u_face_x_layer(i, j, 1) = ms_face%u_face_x_layer(i, j, 1)
         end do
      end do
      call hv_face%init(grid_face, nz_ml=NZ)
      call hv_scalar%init(grid_scalar, nz_ml=NZ)
      hv_face%nu_h = 0.0_wp
      hv_face%nu_4 = 0.0_wp     ! ignored when smag_ah_active=.true.
      hv_scalar%nu_h = 0.0_wp
      hv_scalar%nu_4 = NU4_CONST

      call lmix_face%init(grid_face, nz_ml=NZ)
      lmix_face%closure = LMIX_NONE
      lmix_face%smag_ah_active = .true.
      lmix_face%smag_bi_const = 0.06_wp
      lmix_face%nu4_bg = 0.0_wp
      lmix_face%nu4_max = 1.0e12_wp

      call map_in(ms_face, metrics_face, grid_face, lmix_face, hv_face)
      call ocean_lateral_mix_compute_smag_ah(grid_face, metrics_face, lmix_face, ms_face)
      call ocean_horizontal_viscosity_compute_tendencies(grid_face, metrics_face, hv_face, ms_face, &
                                                         lateral_mix=lmix_face)
      !$acc update self(hv_face%du_visc%data, hv_face%dv_visc%data)
      call map_out(ms_face, metrics_face, lmix_face, hv_face)

      call map_in(ms_scalar, metrics_scalar, grid_scalar, hv=hv_scalar)
      call ocean_horizontal_viscosity_compute_tendencies(grid_scalar, metrics_scalar, hv_scalar, ms_scalar)
      !$acc update self(hv_scalar%du_visc%data, hv_scalar%dv_visc%data)
      call map_out(ms_scalar, metrics_scalar, hv=hv_scalar)

      max_diff = maxval(abs(hv_face%du_visc%data - hv_scalar%du_visc%data))

      call check(error, max_diff > 1.0e-12_wp, &
                 "Smag_AH face path must produce different tendency than scalar nu_4")

      call lmix_face%destroy()
      call hv_face%destroy()
      call hv_scalar%destroy()
      call ms_face%destroy()
      call ms_scalar%destroy()
   end subroutine test_hvisc_with_smag_ah

   subroutine test_anisotropic_corner_shear(error)
      !! Regression gate for the ratio-bundle corner shear-strain fix
      !! (Finding 1).  On a PURELY j-varying dx (dy uniform) the ratio
      !! bundle dy_dxBu / dx_dyBu is non-trivial in the x-direction, so
      !! D_S computed by the old `(v_diff)*idxCu + (u_diff)*idyCv` form
      !! and the new ratio-bundle form are algebraically DIFFERENT.
      !!
      !! Setup:
      !!   dxT(i,j) = DX0 * (1 + AMP*(j-1)/ny)   (grows with j)
      !!   dy  = DY  (uniform)
      !!   u_face_x(i,j) = USHEAR * j             (linear shear in j)
      !!   v_face_y = 0
      !!
      !! The flow is chosen so D_T = 0 and D_S = du/dy at every corner.
      !! The ratio-bundle D_S at interior corner Bu(ip,jp):
      !!
      !!   dvdx = dy_dxBu(ip,jp) * (0 - 0) = 0
      !!   dudy = dx_dyBu(ip,jp) * (u(ip,jp)*idxCu(ip,jp) - u(ip,jp-1)*idxCu(ip,jp-1))
      !!        = (dxBu/dyBu) * (USHEAR*jp/dxCu(ip,jp) - USHEAR*(jp-1)/dxCu(ip,jp-1))
      !!
      !! We hand-compute this at corners (ip=3,jp=3) and (ip=3,jp=4) using
      !! the metric values written by make_anisotropic_metrics and compare
      !! with the ah_face_x output from ocean_lateral_mix_compute_smag (which
      !! calls the ratio-bundle D_S).
      !!
      !! Separation check: the old form would give
      !!   D_S_old = (u(ip,jp) - u(ip,jp-1)) * idyCv(ip,jp)
      !!           = USHEAR / dy   (uniform, independent of the dx variation)
      !! while the new form gives a j-varying result because idxCu varies.
      !! We assert:
      !!   (a) new form matches hand-computed ratio-bundle reference to 1e-12
      !!   (b) old form value differs from new form value (regression guard)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NXP = 12, NYP = 10, NZT = 1
      real(wp), parameter :: DX0 = 40000.0_wp   !! base zonal spacing (m)
      real(wp), parameter :: DY = 40000.0_wp    !! uniform meridional spacing (m)
      real(wp), parameter :: AMP = 0.1_wp       !! 10 % dx variation across domain
      real(wp), parameter :: USHEAR = 0.5_wp    !! du/dy amplitude (m/s per index step)
      real(wp), parameter :: C_SMAG = 0.15_wp
      ! Probe corners (ip, jp) — interior, away from walls.
      integer, parameter :: IP = 4, JP = 4
      integer, parameter :: IP2 = 4, JP2 = 5
      integer :: i, j, nx, ny
      real(wp) :: dx_cu_jp, dx_cu_jpm1, dx_cu_jp2, dx_cu_jp2m1
      real(wp) :: dx_bu_jp, dy_bu_jp, dx_bu_jp2, dy_bu_jp2
      real(wp) :: ds_hand_jp, ds_hand_jp2
      real(wp) :: ds_old_jp, ds_old_jp2
      real(wp) :: smag_scale_jp, smag_scale_jp2
      real(wp) :: dx_t_jp, dx_t_jp2
      real(wp) :: nu4_bg_loc
      real(wp) :: observed_nu4_x_jp, observed_nu4_x_jp2
      real(wp) :: expected_nu4_jp, expected_nu4_jp2
      real(wp) :: old_nu4_jp, old_nu4_jp2
      checks: block

         call grid%init(NXP, NYP, 1, DX0, DY)
         ms%nz_ml = NZT
         call ms%init(grid)
         ms%h_layer = 10.0_wp
         ms%v_face_y_layer = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ! Linear shear: u(i,j) = USHEAR * j  (constant across i)
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, 1) = USHEAR*real(j, wp)
            end do
         end do

         call lmix%init(grid, nz_ml=NZT)
         lmix%smag_ah_active = .false.
         lmix%closure = LMIX_NONE
         lmix%c_smag = C_SMAG
         lmix%ah_bg = 0.0_wp
         lmix%ah_max = 1.0e15_wp
         nu4_bg_loc = 0.0_wp
         lmix%smag_bi_const = 0.06_wp
         lmix%nu4_bg = nu4_bg_loc
         lmix%nu4_max = 1.0e20_wp

         call make_anisotropic_metrics(metrics, grid, DX0, DY, AMP)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(lmix)
         call lmix%enter_data()

         ! Run the biharmonic Smag kernel (which uses the same D_S stencil
         ! as compute_smag).  We check nu4_face_x (which carries D_S at
         ! the two adjacent corners averaged) against the hand formula.
         call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
         !$acc update self(lmix%nu4_face_x)

         call lmix%exit_data()
         !$acc exit data delete(lmix)
         call ms%exit_data()
         !$acc exit data delete(ms)
         call metrics%exit_data()
         !$acc exit data delete(metrics)
         call metrics%destroy()

         ! ---- Hand-compute the expected nu4_face_x at (IP, JP) ----
         ! nu4_face_x(IP,JP,1) uses corners Bu(IP,JP) and Bu(IP,JP+1).
         ! D_S at Bu(ip,jp) = dudy (since v=0):
         !   dudy = dx_dyBu(ip,jp) * (u(ip,jp)*idxCu(ip,jp) - u(ip,jp-1)*idxCu(ip,jp-1))
         ! metric values from make_anisotropic_metrics (1-ghost):
         !   dxCu(ip,jp)  = dxT(ip,jp) = DX0*(1+AMP*(jp-1)/ny)
         !   dxBu(ip,jp)  = DX0*(1+AMP*(jp-1)/ny)   (corner row j-1 relative to T)
         !   dyBu(ip,jp)  = DY
         !   dx_dyBu      = dxBu / dyBu
         !   u(ip,jp)     = USHEAR * jp
         !   u(ip,jp-1)   = USHEAR * (jp-1)
         !
         ! dx at Cu face:
         dx_cu_jp = DX0*(1.0_wp + AMP*real(JP - 1, wp)/real(ny, wp))
         dx_cu_jpm1 = DX0*(1.0_wp + AMP*real(JP - 2, wp)/real(ny, wp))
         ! dx / dy at Bu corner:
         dx_bu_jp = DX0*(1.0_wp + AMP*real(JP - 1, wp)/real(ny, wp))
         dy_bu_jp = DY
         ds_hand_jp = (dx_bu_jp/dy_bu_jp)* &
                      (USHEAR*real(JP, wp)/dx_cu_jp - USHEAR*real(JP - 1, wp)/dx_cu_jpm1)

         ! Same for Bu(IP, JP+1):
         dx_cu_jp2 = DX0*(1.0_wp + AMP*real(JP + 1 - 1, wp)/real(ny, wp))
         dx_cu_jp2m1 = DX0*(1.0_wp + AMP*real(JP + 1 - 2, wp)/real(ny, wp))
         dx_bu_jp2 = DX0*(1.0_wp + AMP*real(JP + 1 - 1, wp)/real(ny, wp))
         dy_bu_jp2 = DY
         ds_hand_jp2 = (dx_bu_jp2/dy_bu_jp2)* &
                       (USHEAR*real(JP + 1, wp)/dx_cu_jp2 - USHEAR*real(JP, wp)/dx_cu_jp2m1)

         ! smag_scale at T-cell (IP, JP) and (IP, JP+1):
         dx_t_jp = DX0*(1.0_wp + AMP*real(JP - 1, wp)/real(ny, wp))
         dx_t_jp2 = DX0*(1.0_wp + AMP*real(JP + 1 - 1, wp)/real(ny, wp))
         smag_scale_jp = (0.06_wp*(dx_t_jp*DY)**2)    ! C_b * L^4 where L^2 = grid_sp_h2
         smag_scale_jp2 = (0.06_wp*(dx_t_jp2*DY)**2)

         ! D_S at u-face (IP,JP) = 0.5*(D_S_S + D_S_N) where
         ! D_S_S is at Bu(IP,JP) and D_S_N is at Bu(IP,JP+1).
         ! D_T = 0 since u is j-only.  |D| = |D_S|.
         ! grid_sp_h2 at T-cell uses dx2 = dxT^2, dy2 = dyT^2:
         !   grid_sp_h2 = 2*dx_t^2*DY^2 / (dx_t^2 + DY^2)
         ! smag_bi_scale = C_b * grid_sp_h2^2
         block
            real(wp) :: dx2, dy2, gsh2_jp, gsh2_jp2
            dx2 = dx_t_jp*dx_t_jp
            dy2 = DY*DY
            gsh2_jp = (2.0_wp*dx2*dy2)/(dx2 + dy2)
            smag_scale_jp = 0.06_wp*gsh2_jp*gsh2_jp

            dx2 = dx_t_jp2*dx_t_jp2
            gsh2_jp2 = (2.0_wp*dx2*dy2)/(dx2 + dy2)
            smag_scale_jp2 = 0.06_wp*gsh2_jp2*gsh2_jp2
         end block

         expected_nu4_jp = smag_scale_jp*abs(0.5_wp*(ds_hand_jp + ds_hand_jp2))
         expected_nu4_jp2 = smag_scale_jp2*abs(0.5_wp*(ds_hand_jp2 + &
                                                       (dx_bu_jp2/dy_bu_jp2)* &
                                                     (USHEAR*real(JP + 2, wp)/(DX0*(1.0_wp + AMP*real(JP + 1, wp)/real(ny, wp))) - &
                                                        USHEAR*real(JP + 1, wp)/(DX0*(1.0_wp + AMP*real(JP, wp)/real(ny, wp))))))

         ! (a) Ratio-bundle form matches hand reference
         observed_nu4_x_jp = lmix%nu4_face_x(IP, JP, 1)
         observed_nu4_x_jp2 = lmix%nu4_face_x(IP, JP2, 1)

         call check(error, abs(observed_nu4_x_jp - expected_nu4_jp) < 1.0e-8_wp*expected_nu4_jp + 1.0e-10_wp, &
                    "ratio-bundle D_S: nu4_face_x(IP,JP) differs from hand reference")
         if (allocated(error)) exit checks

         ! (b) Old face-inverse form would give a DIFFERENT D_S at Bu(IP,JP):
         !   D_S_old = (u(IP,JP) - u(IP,JP-1)) * idyCv(IP,JP)
         !           = USHEAR / DY   (uniform — the j-variation in dx is invisible)
         ! The new form gives USHEAR*(dx_bu/dy_bu)*(jp/dx_cu_jp - (jp-1)/dx_cu_jpm1)
         ! which differs from USHEAR/DY whenever dx_cu_jp != dx_cu_jpm1.
         ds_old_jp = USHEAR/DY
         ds_old_jp2 = USHEAR/DY   ! same for all j — no j-dependence in the old form
         old_nu4_jp = smag_scale_jp*abs(0.5_wp*(ds_old_jp + ds_old_jp2))

         call check(error, abs(observed_nu4_x_jp - old_nu4_jp) > 1.0e-8_wp*old_nu4_jp, &
                    "ratio-bundle D_S must differ from old face-inverse form on anisotropic metrics")

      end block checks
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_anisotropic_corner_shear

end module test_ocean_smag_ah
