!! Analytical + device tests for VarMix (capability [4], `rdb_ocean_varmix`):
!! spatially-varying GM/Redi lateral-diffusivity coefficient face fields.
!! All cases RUN THE DEVICE KERNELS via `varmix_compute` (or, for the seam
!! test, `gm_compute_transports` with the VarMix external base).
!!
!! Bottom-up convention: k=1 bed, k=nz surface; interface K=1 bed,
!! K=nz+1 surface (both zero slope).
!!
!!  1. varmix_resfn_limits     — Ld>>dx ⇒ Res_fn→0; Ld<<dx ⇒ Res_fn→1.
!!  2. varmix_resfn_equatorial — f_centre=0 region ⇒ Res_fn finite (the
!!                               beta_dx2 term keeps it non-singular).
!!  3. varmix_visbeck_scaling  — uniform tilted column ⇒ SN = sqrt(S²N²) and
!!                               KhTh base = KHTH + cff·L²·SN (to ~1e-9).
!!  4. varmix_assembly_order   — Res_fn applied BEFORE the clamp (hand value).
!!  5. varmix_feeds_gm         — GM's khth_u = min(CFL, varmix base) when on;
!!                               with VarMix off GM = min(CFL, const) byte-id.
module test_ocean_varmix
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t
   use rdb_ocean_varmix, only: ocean_varmix_t, varmix_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_varmix_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: KHTH = 1000.0_wp

contains

   subroutine collect_ocean_varmix_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("varmix_resfn_limits", test_resfn_limits), &
                  new_unittest("varmix_resfn_equatorial", test_resfn_equatorial), &
                  new_unittest("varmix_visbeck_scaling", test_visbeck_scaling), &
                  new_unittest("varmix_assembly_order", test_assembly_order), &
                  new_unittest("varmix_feeds_gm", test_feeds_gm), &
                  new_unittest("varmix_sn_orthogonal", test_sn_orthogonal), &
                  new_unittest("varmix_sn_h4_weighted", test_sn_h4_weighted) &
                  ]
   end subroutine collect_ocean_varmix_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine setup_ms(ms, grid, nz, dz)
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = dz
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_ms

   subroutine setup_slopes(sl, grid, nz)
      type(ocean_slopes_t), intent(inout) :: sl
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call sl%init(grid, nz_ml=nz)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%slope_x = 0.0_wp
      sl%slope_y = 0.0_wp
      sl%n2_u = 0.0_wp
      sl%n2_v = 0.0_wp
   end subroutine setup_slopes

   subroutine setup_ws(ws, grid, cg1_val, f_val)
      !! Init the wavespeed slot with a uniform cg1 + a uniform f_centre.
      type(ocean_wave_speed_t), intent(inout) :: ws
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: cg1_val, f_val
      call ws%init(grid)
      ws%enable = .true.
      ws%cg1 = cg1_val
      ws%f_centre = f_val
   end subroutine setup_ws

   subroutine map_in(ms, sl, ws, vm)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_wave_speed_t), intent(inout) :: ws
      type(ocean_varmix_t), intent(inout) :: vm
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(ws)
      call ws%enter_data()
      !$acc enter data copyin(vm)
      call vm%enter_data()
   end subroutine map_in

   subroutine map_out(ms, sl, ws, vm)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_wave_speed_t), intent(inout) :: ws
      type(ocean_varmix_t), intent(inout) :: vm
      !$acc update self(vm%res_fn_u, vm%res_fn_v, vm%sn_u, vm%sn_v)
      !$acc update self(vm%khth_u, vm%khth_v, vm%khtr_u, vm%khtr_v)
      call vm%exit_data()
      !$acc exit data delete(vm)
      call ws%exit_data()
      !$acc exit data delete(ws)
      call sl%exit_data()
      !$acc exit data delete(sl)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! Test 1: resolution-function limits.
   ! ------------------------------------------------------------------
   subroutine test_resfn_limits(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 6, NY = 5, NZ = 4
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp, FMID = 1.0e-4_wp
      real(wp), allocatable :: f_centre(:, :)
      real(wp) :: r_coarse, r_fine

      checks: block
         ! --- Coarse grid: Ld = cg1/f ~ 2 km << 10 km ⇒ Res_fn → 1. ---
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         call setup_ws(ws, grid, cg1_val=0.2_wp, f_val=FMID)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%resoln_scaled_khth = .true.
         vm%khth = KHTH
         allocate (f_centre(grid%nx_total, grid%ny_total), source=FMID)
         call vm%build_static(metrics, f_centre)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)
         r_coarse = vm%res_fn_u(3, 3)
         call destroy_all(grid, ms, metrics, sl, ws, vm, f_centre)

         ! --- Fine grid: Ld = cg1/f ~ 500 km >> 10 km ⇒ Res_fn → 0. ---
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         call setup_ws(ws, grid, cg1_val=50.0_wp, f_val=FMID)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%resoln_scaled_khth = .true.
         vm%khth = KHTH
         allocate (f_centre(grid%nx_total, grid%ny_total), source=FMID)
         call vm%build_static(metrics, f_centre)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)
         r_fine = vm%res_fn_u(3, 3)
         call destroy_all(grid, ms, metrics, sl, ws, vm, f_centre)

         call check(error, r_coarse > 0.9_wp, "coarse grid Res_fn must be ~1")
         if (allocated(error)) exit checks
         call check(error, r_fine < 0.1_wp, "fine grid Res_fn must be ~0")
         if (allocated(error)) exit checks
         call check(error, r_coarse >= 0.0_wp .and. r_coarse <= 1.0_wp, "Res_fn in [0,1]")
         if (allocated(error)) exit checks
         call check(error, r_fine >= 0.0_wp .and. r_fine <= 1.0_wp, "Res_fn in [0,1]")
      end block checks
   end subroutine test_resfn_limits

   ! ------------------------------------------------------------------
   ! Test 2: equatorial (f=0) finite Res_fn (no div-by-zero).
   ! ------------------------------------------------------------------
   subroutine test_resfn_equatorial(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 8, NY = 7, NZ = 4
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: BETA = 2.0e-11_wp
      real(wp), allocatable :: f_centre(:, :)
      integer :: i, j, ng
      real(wp) :: y, r_eq

      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         call setup_ws(ws, grid, cg1_val=2.0_wp, f_val=0.0_wp)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%resoln_scaled_khth = .true.
         vm%khth = KHTH
         ! True equatorial beta-plane f-field: f passes through 0 mid-domain.
         ng = grid%nghost
         allocate (f_centre(grid%nx_total, grid%ny_total))
         do j = 1, grid%ny_total
            y = (real(j - ng, wp) - 0.5_wp)*grid%dy
            do i = 1, grid%nx_total
               f_centre(i, j) = BETA*(y - 0.5_wp*real(NY, wp)*grid%dy)
            end do
         end do
         ws%f_centre = f_centre
         call vm%build_static(metrics, f_centre)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)
         ! Pick the interior row nearest the equator (f≈0).
         r_eq = vm%res_fn_u(4, ng + NY/2)
         call destroy_all(grid, ms, metrics, sl, ws, vm, f_centre)

         call check(error, r_eq == r_eq, "Res_fn must be finite (not NaN) at the equator")
         if (allocated(error)) exit checks
         call check(error, r_eq > 0.0_wp .and. r_eq <= 1.0_wp, &
                    "equatorial Res_fn must be in (0,1] (beta_dx2 keeps it non-singular)")
      end block checks
   end subroutine test_resfn_equatorial

   ! ------------------------------------------------------------------
   ! Test 3: Visbeck SN + KhTh base scaling for a uniform tilted column.
   ! ------------------------------------------------------------------
   subroutine test_visbeck_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 6, NY = 5, NZ = 6
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: SLOPE = 0.005_wp, N2 = 1.0e-5_wp
      real(wp), parameter :: LSCALE = 50000.0_wp, CFF = 0.1_wp
      real(wp) :: sn_expect, kh_expect, sn_got, kh_got

      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         ! Uniform slope_x, zero slope_y, uniform N² on interior interfaces.
         sl%slope_x = SLOPE
         sl%slope_y = 0.0_wp
         sl%n2_u = N2
         sl%n2_v = N2
         call setup_ws(ws, grid, cg1_val=2.0_wp, f_val=1.0e-4_wp)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%use_visbeck = .true.
         vm%resoln_scaled_khth = .false.   ! isolate the Visbeck term
         vm%khth = KHTH
         vm%khth_slope_cff = CFF
         vm%visbeck_l_scale = LSCALE
         call build_uniform_static(vm, metrics, grid, 1.0e-4_wp)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)

         ! slope_y = 0 ⇒ S² = slope_x² ⇒ SN_u = sqrt(SLOPE²·N2).  (The
         ! orthogonal slope is folded into S² in varmix_sn_u, not a separate
         ! combine — the slope_y/=0 case is covered by test_sn_orthogonal.)
         sn_expect = sqrt(SLOPE*SLOPE*N2)
         kh_expect = KHTH + CFF*LSCALE*LSCALE*sn_expect
         sn_got = vm%sn_u(3, 3)
         kh_got = vm%khth_u(3, 3)
         call destroy_all_static(grid, ms, metrics, sl, ws, vm)

         call check(error, abs(sn_got - sn_expect) < 1.0e-12_wp, "SN = sqrt(S^2 N^2)")
         if (allocated(error)) exit checks
         call check(error, abs(kh_got - kh_expect) < 1.0e-9_wp*kh_expect, &
                    "KhTh base = KHTH + cff*L^2*SN")
      end block checks
   end subroutine test_visbeck_scaling

   ! Test 6: SN with BOTH slopes nonzero — the orthogonal slope_y is folded
   ! into S^2 (MOM6 calc_Visbeck_coeffs_old), so SN_u = sqrt((Sx^2+Sy^2)*N2).
   ! This is the regression gate for the double-count bug: the spurious
   ! 4-corner SN_v combine would give sqrt(2*(Sx^2+Sy^2))*N — caught here.
   subroutine test_sn_orthogonal(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 6, NY = 5, NZ = 6
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: SX = 0.005_wp, SY = 0.003_wp, N2 = 1.0e-5_wp
      real(wp), parameter :: LSCALE = 50000.0_wp, CFF = 0.1_wp
      real(wp) :: sn_expect, sn_got, sn_combine_bug
      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         ! Uniform slope_x AND slope_y ⇒ S^2 = Sx^2 + Sy^2 (h-weighted
         ! orthogonal average of a uniform slope_y^2 is slope_y^2).
         sl%slope_x = SX
         sl%slope_y = SY
         sl%n2_u = N2
         sl%n2_v = N2
         call setup_ws(ws, grid, cg1_val=2.0_wp, f_val=1.0e-4_wp)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%use_visbeck = .true.
         vm%resoln_scaled_khth = .false.
         vm%khth = KHTH
         vm%khth_slope_cff = CFF
         vm%visbeck_l_scale = LSCALE
         call build_uniform_static(vm, metrics, grid, 1.0e-4_wp)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)

         sn_expect = sqrt((SX*SX + SY*SY)*N2)              ! orthogonal in S^2
         sn_combine_bug = sqrt(2.0_wp*(SX*SX + SY*SY)*N2)  ! the rejected combine
         sn_got = vm%sn_u(3, 3)
         call destroy_all_static(grid, ms, metrics, sl, ws, vm)

         call check(error, abs(sn_got - sn_expect) < 1.0e-12_wp, &
                    "SN_u = sqrt((Sx^2+Sy^2) N^2): orthogonal folded into S^2")
         if (allocated(error)) exit checks
         ! Guard: the value must NOT be the double-counted combine.
         call check(error, abs(sn_got - sn_combine_bug) > 1.0e-9_wp, &
                    "SN_u must not double-count the orthogonal slope (no SN_v combine)")
      end block checks
   end subroutine test_sn_orthogonal

   ! ------------------------------------------------------------------
   ! Test 7: SN_u 4-corner slope_y^2 weighting uses MOM6's h4_v PRODUCT
   ! (the 4 thicknesses straddling each corner's slope_y v-point, at the
   ! two layers the interface separates), NOT a single cell-centre h.
   ! Variable h ⇒ the four corner weights differ, so the choice is
   ! observable.  Oracle is hand-computed from the same h4 products; the
   ! OLD single-thickness weighting is asserted to give a DIFFERENT value.
   ! NZ=2 ⇒ single interior interface (k=2), so the hgeom column-reduction
   ! cancels and SN_u = sqrt(sy2 * N2) with sy2 the h4-weighted average.
   ! ------------------------------------------------------------------
   subroutine test_sn_h4_weighted(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 6, NY = 5, NZ = 2
      integer, parameter :: IF_U = 3, JR = 3   ! u-face index + interior cell-row
      integer, parameter :: IW = IF_U - 1       ! west cell of the u-face
      real(wp), parameter :: DX = 10000.0_wp
      real(wp), parameter :: N2 = 1.0e-5_wp
      ! Distinct corner slope_y values (squared in S^2).
      real(wp), parameter :: SY_NW = 0.001_wp, SY_NE = 0.004_wp
      real(wp), parameter :: SY_SW = 0.002_wp, SY_SE = 0.006_wp
      real(wp), parameter :: LSCALE = 50000.0_wp, CFF = 0.1_wp
      real(wp) :: wnw, wne, wsw, wse, sy2, sn_expect, sn_got
      real(wp) :: w_old_nw, w_old_ne, w_old_sw, w_old_se, sy2_old, sn_old
      integer :: i, j, k
      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         ! NZ=2: setup_ms sets a uniform h; we then overwrite with a
         ! distinct value per (cell, layer) so the 4 h4 weights differ.
         call setup_ms(ms, grid, NZ, 50.0_wp)
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  ! Smoothly varying, strictly positive, distinct per cell+layer.
                  ms%h_layer(i, j, k) = 40.0_wp + 3.0_wp*real(i, wp) &
                                        + 5.0_wp*real(j, wp) + 7.0_wp*real(k, wp)
               end do
            end do
         end do

         call setup_slopes(sl, grid, NZ)
         sl%slope_x = 0.0_wp              ! isolate the orthogonal slope_y term
         sl%slope_y = 0.0_wp
         ! The 4 v-faces around u-face IF_U at cell-row JR (interface k=2):
         !   NW = slope_y(IW,   JR+1, 2),  NE = slope_y(IF_U, JR+1, 2)
         !   SW = slope_y(IW,   JR,   2),  SE = slope_y(IF_U, JR,   2)
         sl%slope_y(IW, JR + 1, 2) = SY_NW
         sl%slope_y(IF_U, JR + 1, 2) = SY_NE
         sl%slope_y(IW, JR, 2) = SY_SW
         sl%slope_y(IF_U, JR, 2) = SY_SE
         sl%n2_u = N2
         sl%n2_v = N2

         call setup_ws(ws, grid, cg1_val=2.0_wp, f_val=1.0e-4_wp)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%use_visbeck = .true.
         vm%resoln_scaled_khth = .false.
         vm%visbeck_max_slope = 0.0_wp    ! no S2 limiting
         vm%khth = KHTH
         vm%khth_slope_cff = CFF
         vm%visbeck_l_scale = LSCALE
         call build_uniform_static(vm, metrics, grid, 1.0e-4_wp)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)
         sn_got = vm%sn_u(IF_U, JR)

         ! Oracle: MOM6 h4_v product co-located with each corner's slope_y
         ! v-point — the two cell-rows the v-face separates, layers k=2 & 1.
         !   NW (cell IW,  v-face JR+1): rows JR, JR+1
         !   NE (cell IF_U,v-face JR+1): rows JR, JR+1
         !   SW (cell IW,  v-face JR  ): rows JR-1, JR
         !   SE (cell IF_U,v-face JR  ): rows JR-1, JR
         wnw = (ms%h_layer(IW, JR, 2)*ms%h_layer(IW, JR + 1, 2)) &
               *(ms%h_layer(IW, JR, 1)*ms%h_layer(IW, JR + 1, 1))
         wne = (ms%h_layer(IF_U, JR, 2)*ms%h_layer(IF_U, JR + 1, 2)) &
               *(ms%h_layer(IF_U, JR, 1)*ms%h_layer(IF_U, JR + 1, 1))
         wsw = (ms%h_layer(IW, JR - 1, 2)*ms%h_layer(IW, JR, 2)) &
               *(ms%h_layer(IW, JR - 1, 1)*ms%h_layer(IW, JR, 1))
         wse = (ms%h_layer(IF_U, JR - 1, 2)*ms%h_layer(IF_U, JR, 2)) &
               *(ms%h_layer(IF_U, JR - 1, 1)*ms%h_layer(IF_U, JR, 1))
         sy2 = ((wnw*SY_NW*SY_NW + wse*SY_SE*SY_SE) &
                + (wne*SY_NE*SY_NE + wsw*SY_SW*SY_SW))/(wnw + wne + wsw + wse)
         sn_expect = sqrt(sy2*N2)

         ! The REJECTED single-thickness weighting (old code): the layer-k=2
         ! cell-centre h.  Asserted DIFFERENT so the test guards the fix.
         w_old_nw = ms%h_layer(IW, JR, 2)
         w_old_ne = ms%h_layer(IF_U, JR, 2)
         w_old_sw = ms%h_layer(IW, JR, 2)
         w_old_se = ms%h_layer(IF_U, JR, 2)
         sy2_old = ((w_old_nw*SY_NW*SY_NW + w_old_se*SY_SE*SY_SE) &
                    + (w_old_ne*SY_NE*SY_NE + w_old_sw*SY_SW*SY_SW)) &
                   /(w_old_nw + w_old_ne + w_old_sw + w_old_se)
         sn_old = sqrt(sy2_old*N2)

         call destroy_all_static(grid, ms, metrics, sl, ws, vm)

         call check(error, abs(sn_got - sn_expect) < 1.0e-12_wp, &
                    "SN_u 4-corner slope_y^2 weighted by the MOM6 h4_v product")
         if (allocated(error)) exit checks
         ! Guard: the h4 weighting must differ from the old single-h weighting.
         call check(error, abs(sn_expect - sn_old) > 1.0e-9_wp, &
                    "h4 weighting must differ from single-thickness (guards the fix)")
      end block checks
   end subroutine test_sn_h4_weighted

   ! ------------------------------------------------------------------
   ! Test 4: assembly order — Res_fn applied BEFORE the clamp.
   ! ------------------------------------------------------------------
   subroutine test_assembly_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 6, NY = 5, NZ = 4
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp, FMID = 1.0e-4_wp
      real(wp), allocatable :: f_centre(:, :)
      real(wp) :: r, kh_min, kh_max, kh_got, kh_expect, kh_clamp_first

      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         ! Coarse grid ⇒ Res_fn ~ 1 but not exactly; we read it back below.
         call setup_ws(ws, grid, cg1_val=2.0_wp, f_val=FMID)
         ! kh_max < KHTH so the clamp BINDS: Res_fn-before-clamp gives
         ! clamp(KHTH*r) = kh_max (since KHTH*r > kh_max), while clamp-first
         ! gives clamp(KHTH)*r = kh_max*r — the two orders differ, so the
         ! test genuinely discriminates the ordering.
         kh_min = 50.0_wp
         kh_max = 900.0_wp
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%use_visbeck = .false.
         vm%resoln_scaled_khth = .true.
         vm%khth = KHTH
         vm%khth_min = kh_min
         vm%khth_max = kh_max
         allocate (f_centre(grid%nx_total, grid%ny_total), source=FMID)
         call vm%build_static(metrics, f_centre)

         call map_in(ms, sl, ws, vm)
         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         call map_out(ms, sl, ws, vm)
         r = vm%res_fn_u(3, 3)
         kh_got = vm%khth_u(3, 3)
         call destroy_all(grid, ms, metrics, sl, ws, vm, f_centre)

         ! Correct (Res_fn BEFORE clamp): clamp(KHTH * r).
         kh_expect = max(kh_min, min(KHTH*r, kh_max))
         ! WRONG order (clamp first): clamp(KHTH) * r — must DIFFER here so
         ! the test actually discriminates the ordering.
         kh_clamp_first = max(kh_min, min(KHTH, kh_max))*r
         call check(error, abs(kh_got - kh_expect) < 1.0e-9_wp*max(kh_expect, 1.0_wp), &
                    "KhTh = clamp(KHTH * Res_fn) (Res_fn before clamp)")
         if (allocated(error)) exit checks
         call check(error, abs(kh_expect - kh_clamp_first) > 1.0e-6_wp, &
                    "the two orderings must differ so the test discriminates")
      end block checks
   end subroutine test_assembly_order

   ! ------------------------------------------------------------------
   ! Test 5: VarMix feeds GM — GM khth_u = min(CFL, varmix base); off ⇒
   !         min(CFL, const) byte-identical.
   ! ------------------------------------------------------------------
   subroutine test_feeds_gm(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_wave_speed_t) :: ws
      type(ocean_varmix_t) :: vm
      type(ocean_gm_t) :: gm_on, gm_off
      integer, parameter :: NX = 6, NY = 5, NZ = 4
      real(wp), parameter :: DX = 10000.0_wp, DZ = 50.0_wp, FMID = 1.0e-4_wp
      real(wp), parameter :: DTT = 1800.0_wp
      real(wp), allocatable :: f_centre(:, :)
      real(wp) :: base, kh_on, kh_off, mx

      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         call setup_ms(ms, grid, NZ, DZ)
         call setup_slopes(sl, grid, NZ)
         call setup_ws(ws, grid, cg1_val=0.2_wp, f_val=FMID)  ! coarse ⇒ Res~1
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%resoln_scaled_khth = .true.
         vm%khth = KHTH
         ! Cap the VarMix base WELL BELOW khth so the base is unambiguously
         ! distinct from the constant khth the off-path uses.
         vm%khth_max = 300.0_wp
         allocate (f_centre(grid%nx_total, grid%ny_total), source=FMID)
         call vm%build_static(metrics, f_centre)
         ! Both GM slots share khth = KHTH; no CFL cap (so the seam alone
         ! sets the difference).
         call setup_gm(gm_on, grid, NZ)
         call setup_gm(gm_off, grid, NZ)

         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(sl)
         call sl%enter_data()
         !$acc enter data copyin(ws)
         call ws%enter_data()
         !$acc enter data copyin(vm)
         call vm%enter_data()
         !$acc enter data copyin(gm_on)
         call gm_on%enter_data()
         !$acc enter data copyin(gm_off)
         call gm_off%enter_data()

         call varmix_compute(grid, metrics, vm, sl, ws, ms)
         ! GM with the VarMix external base.
         call gm_compute_transports(grid, metrics, gm_on, sl, ms, DTT, &
                                    khth_ext_u=vm%khth_u, khth_ext_v=vm%khth_v)
         ! GM with the constant khth (VarMix off path).
         call gm_compute_transports(grid, metrics, gm_off, sl, ms, DTT)

         !$acc update self(vm%khth_u, gm_on%khth_u, gm_off%khth_u)
         call gm_on%exit_data()
         !$acc exit data delete(gm_on)
         call gm_off%exit_data()
         !$acc exit data delete(gm_off)
         call vm%exit_data()
         !$acc exit data delete(vm)
         call ws%exit_data()
         !$acc exit data delete(ws)
         call sl%exit_data()
         !$acc exit data delete(sl)
         call ms%exit_data()
         !$acc exit data delete(ms)

         ! VarMix base capped at 300; gm_max_cfl is huge so no CFL bind:
         !   GM-on  face = min(CFL, base)  = base    (~300, the VarMix cap)
         !   GM-off face = min(CFL, khth)  = khth    (=1000, the constant)
         base = vm%khth_u(3, 3)
         kh_on = gm_on%khth_u(3, 3)
         kh_off = gm_off%khth_u(3, 3)
         call destroy_gm(gm_on)
         call destroy_gm(gm_off)
         call ms%destroy()
         call sl%destroy()
         call ws%destroy()
         call vm%destroy()
         call destroy_cartesian_metrics(metrics)
         deallocate (f_centre)

         mx = 1.0e-9_wp*max(base, 1.0_wp)
         ! GM-on must take the VarMix base (no CFL bind ⇒ exactly the base).
         call check(error, abs(kh_on - base) < mx, "GM-on face = VarMix base")
         if (allocated(error)) exit checks
         ! GM-off must take the constant khth (byte-identical to pre-VarMix).
         call check(error, abs(kh_off - KHTH) < 1.0e-9_wp*KHTH, &
                    "GM-off face = constant khth (bit-identical fallback)")
         if (allocated(error)) exit checks
         ! And the two paths genuinely differ (the seam is exercised).
         call check(error, abs(kh_on - kh_off) > 1.0_wp, &
                    "VarMix base differs from the constant ⇒ seam active")
      end block checks
   end subroutine test_feeds_gm

   ! ------------------------------------------------------------------
   ! Shared GM + teardown helpers
   ! ------------------------------------------------------------------

   subroutine setup_gm(gm, grid, nz)
      type(ocean_gm_t), intent(inout) :: gm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%khth = KHTH
      gm%khth_max_cfl = 1.0e6_wp   ! effectively no CFL clamp
      gm%khth_slope_max = 0.01_wp
      gm%rho0 = RHO0
   end subroutine setup_gm

   subroutine destroy_gm(gm)
      type(ocean_gm_t), intent(inout) :: gm
      call gm%destroy()
   end subroutine destroy_gm

   subroutine build_uniform_static(vm, metrics, grid, f_val)
      type(ocean_varmix_t), intent(inout) :: vm
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: f_val
      real(wp), allocatable :: f_centre(:, :)
      allocate (f_centre(grid%nx_total, grid%ny_total), source=f_val)
      call vm%build_static(metrics, f_centre)
      deallocate (f_centre)
   end subroutine build_uniform_static

   subroutine destroy_all(grid, ms, metrics, sl, ws, vm, f_centre)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_wave_speed_t), intent(inout) :: ws
      type(ocean_varmix_t), intent(inout) :: vm
      real(wp), allocatable, intent(inout) :: f_centre(:, :)
      call vm%destroy()
      call ws%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      if (allocated(f_centre)) deallocate (f_centre)
   end subroutine destroy_all

   subroutine destroy_all_static(grid, ms, metrics, sl, ws, vm)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_slopes_t), intent(inout) :: sl
      type(ocean_wave_speed_t), intent(inout) :: ws
      type(ocean_varmix_t), intent(inout) :: vm
      call vm%destroy()
      call ws%destroy()
      call sl%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine destroy_all_static

end module test_ocean_varmix
