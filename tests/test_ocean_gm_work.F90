!! Analytical + device tests for MEKE's GM energy source from the work of the
!! FINAL bolus transports, `&ocean_gm_nml gm_src_form = "work"`
!! (`rdb_ocean_gm :: gm_work_source`, the inverse of MOM6 `MEKE_GM_SRC_ALT`).
!!
!! Bottom-up convention: layer k lies between interfaces k (bed side) and
!! k+1; `U_K = Sum_{k>=K} uhD_k` is the bolus transport ABOVE interface K.
!!
!! The two-column fixture (nx=2, nz=2, one face at i=2) is shared by tests 1-3:
!! a tilted interface (e_2 = -20 m left, -16 m right), flat bed, level
!! surface, linear EOS, a prescribed bolus transport +F in the top layer and
!! -F in the bed layer (U_2 = F).
!!
!!  1. gm_work_two_column     — gm_src equals the hand-evaluated
!!                              W = g [F dRho_2 ebar_2 - F dRho_1 ebar_1
!!                                     + F (D_R - D_L)] to 1e-12 (relative),
!!                              the same value in both cells (0.5 W / A).
!!  2. gm_work_pe_budget      — W equals -dPE/dt of an independent
!!                              conservative advance of the same two columns
!!                              (central difference of the PE) to 3e-3 (the
!!                              linear EOS' (rho/rho0)^2 factor in drho/dT is
!!                              ~2e-4 here; a wrong sign or a missing term is
!!                              O(1)); and reversing the transport reverses
!!                              the sign of W exactly.
!!  3. gm_work_zero_transport — uhD = 0 gives gm_src = 0 exactly.
!!  4. gm_work_sign_front     — real GM (`gm_compute_transports`) on a tilted
!!                              stable front: interior gm_src > 0.
!!  5. gm_work_vs_local_small — uniform N^2, small slope, uniform KH, no
!!                              limiter: "work" = 2 x "local" to 3 %.  The
!!                              local form's 1/4 averages the four faces'
!!                              rho0 KH S^2 N^2 h, i.e. HALF the full release
!!                              rho0 KH (Sx^2+Sy^2) N^2 h that test 2 verifies
!!                              against dPE/dt; the residual 0.3 % is the
!!                              stencil difference (one-sided interface vs
!!                              centred layer, 6 layers, 20 m).
!!  6. gm_work_fgnv_sensitive — with use_fgnv the "work" source changes
!!                              (> 1e-3 relative) while the "local" one is
!!                              bit-identical to the non-FGNV run.
!!  7. gm_work_config         — gm_src_form accepts local/work and refuses "work" without GM (an unknown
!!                              spelling stops at the schema parse).
module test_ocean_gm_work
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t, ocean_slopes_compute
   use rdb_ocean_gm, only: ocean_gm_t, gm_compute_transports, gm_work_source, gm_work_source_impl
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_gm_work_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: AREA = 1.0e6_wp
   real(wp), parameter :: FLUX = 100.0_wp

   ! Two-column fixture: (column, layer) = (L/R, bed/top).
   real(wp), parameter :: H_L(2) = [30.0_wp, 20.0_wp]
   real(wp), parameter :: H_R(2) = [34.0_wp, 16.0_wp]
   real(wp), parameter :: T_L(2) = [10.1_wp, 10.5_wp]
   real(wp), parameter :: T_R(2) = [10.3_wp, 10.8_wp]
   real(wp), parameter :: S_L(2) = [35.00_wp, 35.04_wp]
   real(wp), parameter :: S_R(2) = [34.98_wp, 35.02_wp]
   real(wp), parameter :: BED = 50.0_wp

contains

   subroutine collect_ocean_gm_work_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("gm_work_two_column", test_two_column), &
                  new_unittest("gm_work_pe_budget", test_pe_budget), &
                  new_unittest("gm_work_zero_transport", test_zero_transport), &
                  new_unittest("gm_work_sign_front", test_sign_front), &
                  new_unittest("gm_work_vs_local_small", test_vs_local_small), &
                  new_unittest("gm_work_fgnv_sensitive", test_fgnv_sensitive), &
                  new_unittest("gm_work_config", test_config) &
                  ]
   end subroutine collect_ocean_gm_work_tests

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine make_eos

   pure function rho_lin(t, s) result(r)
      real(wp), intent(in) :: t, s
      real(wp) :: r
      r = RHO0 - ALPHA_T*(t - T0) + BETA_S*(s - S0)
   end function rho_lin

   subroutine two_column_src(flux_top, src)
      !! Run the production flat kernel on the two-column fixture with
      !! uhD(2,1,2) = flux_top, uhD(2,1,1) = -flux_top; returns gm_src(2,1).
      real(wp), intent(in) :: flux_top
      real(wp), intent(out) :: src(2, 1)
      type(eos_t) :: eos
      real(wp) :: areaT(2, 1), h(2, 1, 2), bathy(2, 1), t(2, 1, 2), s(2, 1, 2)
      real(wp) :: uhD(3, 1, 2), vhD(2, 2, 2), ou(3, 1, 2), ov(2, 2, 2)
      integer :: k
      call make_eos(eos)
      areaT = AREA
      bathy = BED
      do k = 1, 2
         h(1, 1, k) = H_L(k); h(2, 1, k) = H_R(k)
         t(1, 1, k) = T_L(k); t(2, 1, k) = T_R(k)
         s(1, 1, k) = S_L(k); s(2, 1, k) = S_R(k)
      end do
      uhD = 0.0_wp
      vhD = 0.0_wp
      ou = 1.0_wp
      ov = 1.0_wp
      uhD(2, 1, 2) = flux_top
      uhD(2, 1, 1) = -flux_top
      src = 0.0_wp
      !$acc enter data copyin(areaT, h, bathy, t, s, uhD, vhD, ou, ov, src)
      call gm_work_source_impl(2, 1, 2, eos, RHO0, .false., areaT, h, bathy, t, s, &
                               uhD, vhD, ou, ov, src)
      !$acc update self(src)
      !$acc exit data delete(areaT, h, bathy, t, s, uhD, vhD, ou, ov, src)
   end subroutine two_column_src

   ! ------------------------------------------------------------------
   ! Test 1: hand-evaluated two-column release
   ! ------------------------------------------------------------------
   subroutine test_two_column(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: src(2, 1), w, f, tm, sm, rm
      real(wp) :: d2, d1, eb2, eb1, drl, drr, dt_c, ds_c
      real(wp), parameter :: EL(3) = [-50.0_wp, -20.0_wp, 0.0_wp]
      real(wp), parameter :: ER(3) = [-50.0_wp, -16.0_wp, 0.0_wp]
      checks: block
         call two_column_src(FLUX, src)
         ! Layer 2 (top): factor f = (rho/rho0)^2 at the face-mean T/S, then
         ! drho/dT = -f alpha, drho/dS = +f beta (locally referenced slope).
         tm = 0.5_wp*(T_L(2) + T_R(2)); sm = 0.5_wp*(S_L(2) + S_R(2))
         rm = rho_lin(tm, sm); f = (rm/RHO0)**2
         d2 = -f*ALPHA_T*(T_R(2) - T_L(2)) + f*BETA_S*(S_R(2) - S_L(2))
         eb2 = 0.25_wp*((EL(2) + EL(3)) + (ER(2) + ER(3)))
         ! Layer 1 (bed).
         tm = 0.5_wp*(T_L(1) + T_R(1)); sm = 0.5_wp*(S_L(1) + S_R(1))
         rm = rho_lin(tm, sm); f = (rm/RHO0)**2
         d1 = -f*ALPHA_T*(T_R(1) - T_L(1)) + f*BETA_S*(S_R(1) - S_L(1))
         eb1 = 0.25_wp*((EL(1) + EL(2)) + (ER(1) + ER(2)))
         ! Interface 2: four-point T/S; vertical (below - above) differences.
         tm = 0.25_wp*(T_L(1) + T_L(2) + T_R(1) + T_R(2))
         sm = 0.25_wp*(S_L(1) + S_L(2) + S_R(1) + S_R(2))
         rm = rho_lin(tm, sm); f = (rm/RHO0)**2
         dt_c = T_L(1) - T_L(2); ds_c = S_L(1) - S_L(2)
         drl = -f*ALPHA_T*dt_c + f*BETA_S*ds_c
         dt_c = T_R(1) - T_R(2); ds_c = S_R(1) - S_R(2)
         drr = -f*ALPHA_T*dt_c + f*BETA_S*ds_c
         w = GRAVITY*FLUX*(d2*eb2 - d1*eb1 + (drr*ER(2) - drl*EL(2)))
         call check(error, abs(w) > 1.0e-3_wp, "two_column: hand value must be non-trivial")
         if (allocated(error)) exit checks
         call check(error, abs(src(1, 1) - 0.5_wp*w/AREA) <= 1.0e-12_wp*abs(0.5_wp*w/AREA), &
                    "two_column: left cell gm_src = 0.5 W / A")
         if (allocated(error)) exit checks
         call check(error, abs(src(2, 1) - 0.5_wp*w/AREA) <= 1.0e-12_wp*abs(0.5_wp*w/AREA), &
                    "two_column: right cell gm_src = 0.5 W / A")
      end block checks
   end subroutine test_two_column

   ! ------------------------------------------------------------------
   ! Test 2: W = -dPE/dt of an independent conservative advance
   ! ------------------------------------------------------------------
   function two_col_pe(eps) result(pe)
      !! PE (J) after advancing the two columns by `eps` seconds under the
      !! bolus fluxes: layer mass moves by uhD, T/S ride the face-mean value.
      real(wp), intent(in) :: eps
      real(wp) :: pe
      real(wp) :: h(2, 2), t(2, 2), s(2, 2), fl(2), e(3), tb, sb
      integer :: k, c
      fl = [-FLUX, FLUX]                 ! bed layer, top layer (L -> R)
      do k = 1, 2
         h(1, k) = H_L(k); h(2, k) = H_R(k)
         t(1, k) = T_L(k); t(2, k) = T_R(k)
         s(1, k) = S_L(k); s(2, k) = S_R(k)
         tb = 0.5_wp*(t(1, k) + t(2, k))
         sb = 0.5_wp*(s(1, k) + s(2, k))
         ! (h T) conservative update, then divide by the new h.
         t(1, k) = (h(1, k)*t(1, k) - eps*fl(k)/AREA*tb)/(h(1, k) - eps*fl(k)/AREA)
         t(2, k) = (h(2, k)*T_R(k) + eps*fl(k)/AREA*tb)/(h(2, k) + eps*fl(k)/AREA)
         s(1, k) = (h(1, k)*S_L(k) - eps*fl(k)/AREA*sb)/(h(1, k) - eps*fl(k)/AREA)
         s(2, k) = (h(2, k)*S_R(k) + eps*fl(k)/AREA*sb)/(h(2, k) + eps*fl(k)/AREA)
         h(1, k) = h(1, k) - eps*fl(k)/AREA
         h(2, k) = h(2, k) + eps*fl(k)/AREA
      end do
      pe = 0.0_wp
      do c = 1, 2
         e(1) = -BED
         e(2) = e(1) + h(c, 1)
         e(3) = e(2) + h(c, 2)
         do k = 1, 2
            pe = pe + AREA*GRAVITY*rho_lin(t(c, k), s(c, k))*0.5_wp*(e(k + 1)**2 - e(k)**2)
         end do
      end do
   end function two_col_pe

   subroutine test_pe_budget(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: src(2, 1), srcn(2, 1), w, dpe
      real(wp), parameter :: EPS = 10.0_wp
      checks: block
         call two_column_src(FLUX, src)
         call two_column_src(-FLUX, srcn)
         w = 2.0_wp*src(1, 1)*AREA                     ! face work (W)
         dpe = (two_col_pe(-EPS) - two_col_pe(EPS))/(2.0_wp*EPS)   ! -dPE/dt
         call check(error, abs(w) > 1.0_wp, "pe_budget: non-trivial release")
         if (allocated(error)) exit checks
         call check(error, abs(w - dpe) <= 3.0e-3_wp*abs(dpe), &
                    "pe_budget: W must equal -dPE/dt of the conservative advance")
         if (allocated(error)) exit checks
         call check(error, srcn(1, 1) == -src(1, 1) .and. srcn(2, 1) == -src(2, 1), &
                    "pe_budget: reversing the transport must reverse W exactly")
      end block checks
   end subroutine test_pe_budget

   subroutine test_zero_transport(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: src(2, 1)
      call two_column_src(0.0_wp, src)
      call check(error, src(1, 1) == 0.0_wp .and. src(2, 1) == 0.0_wp, &
                 "zero_transport: no transport, no work")
   end subroutine test_zero_transport

   ! ------------------------------------------------------------------
   ! Real-GM harness (8x6 cells, tilted stable front)
   ! ------------------------------------------------------------------
   subroutine run_gm(nz, gx, gy, khth, use_fgnv, cg1_val, src_work, src, nxt, nyt)
      integer, intent(in) :: nz
      real(wp), intent(in) :: gx, gy, khth, cg1_val
      logical, intent(in) :: use_fgnv, src_work
      real(wp), allocatable, intent(out) :: src(:, :)
      integer, intent(out) :: nxt, nyt
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_slopes_t) :: sl
      type(ocean_gm_t) :: gm
      type(eos_t) :: eos
      real(wp), allocatable :: cg1(:, :)
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp, DZ = 20.0_wp, GZ = 0.02_wp
      integer :: i, j, k
      call grid%init(8, 6, NGHOST, DX, DX)
      nxt = grid%nx_total
      nyt = grid%ny_total
      ms%nz_ml = nz
      call ms%init(grid)
      call make_cartesian_metrics(metrics, grid)
      call sl%init(grid, nz_ml=nz)
      sl%enable = .true.
      sl%rho0 = RHO0
      sl%kd_smooth = 0.0_wp
      sl%min_dz_for_n2 = 1.0_wp
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%khth = khth
      gm%khth_max_cfl = 1.0e6_wp
      gm%khth_slope_max = 0.01_wp
      gm%rho0 = RHO0
      gm%use_fgnv = use_fgnv
      gm%fgnv_filter_scale = 0.1_wp
      gm%src_work = src_work
      call make_eos(eos)
      ms%h_layer = DZ
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         do j = 1, nyt
            do i = 1, nxt
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = DZ*(T0 + GZ*(real(k, wp) - 0.5_wp)*DZ &
                                                                 + gx*real(i, wp)*DX + gy*real(j, wp)*DX)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
            end do
         end do
      end do
      allocate (cg1(nxt, nyt), source=cg1_val)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call sl%set_bathymetry(sum(ms%h_layer, dim=3))
      !$acc enter data copyin(sl)
      call sl%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
      !$acc enter data copyin(cg1)

      call ocean_slopes_compute(grid, metrics, eos, sl, ms, DT)
      if (use_fgnv) then
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT, cg1=cg1)
      else
         call gm_compute_transports(grid, metrics, gm, sl, ms, DT)
      end if
      call gm_work_source(grid, metrics, eos, gm, sl, ms)
      !$acc update self(gm%gm_src)
      src = gm%gm_src

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
   end subroutine run_gm

   ! ------------------------------------------------------------------
   ! Test 4: GM flattening a front releases PE
   ! ------------------------------------------------------------------
   subroutine test_sign_front(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: src(:, :)
      integer :: nxt, nyt
      checks: block
         call run_gm(4, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .false., 0.0_wp, .true., src, nxt, nyt)
         call check(error, minval(src(3:nxt - 2, 3:nyt - 2)) > 0.0_wp, &
                    "sign_front: GM release must give gm_src > 0 in every interior cell")
         if (allocated(error)) exit checks
         call check(error, minval(src) >= -1.0e-12_wp*maxval(src), &
                    "sign_front: no cell may take energy from the front")
      end block checks
   end subroutine test_sign_front

   ! ------------------------------------------------------------------
   ! Test 5: agreement with the local form for small slopes
   ! ------------------------------------------------------------------
   subroutine test_vs_local_small(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: sw(:, :), sl(:, :)
      integer :: nxt, nyt
      real(wp) :: ratio
      checks: block
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .false., 0.0_wp, .true., sw, nxt, nyt)
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .false., 0.0_wp, .false., sl, nxt, nyt)
         ratio = sum(sw(3:nxt - 2, 3:nyt - 2))/sum(sl(3:nxt - 2, 3:nyt - 2))
         write (*, "(a,es12.5,a,es12.5,a,f8.5)") "  [gm_work] work=", sum(sw(3:nxt - 2, 3:nyt - 2)), &
            " local=", sum(sl(3:nxt - 2, 3:nyt - 2)), " ratio=", ratio
         call check(error, abs(0.5_wp*ratio - 1.0_wp) <= 0.03_wp, &
                    "vs_local: work = 2 x local (full vs quarter-weighted release) within 3%")
      end block checks
   end subroutine test_vs_local_small

   ! ------------------------------------------------------------------
   ! Test 6: FGNV changes "work" but not "local"
   ! ------------------------------------------------------------------
   subroutine test_fgnv_sensitive(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w0(:, :), w1(:, :), l0(:, :), l1(:, :)
      integer :: nxt, nyt
      real(wp) :: rel
      checks: block
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .false., 0.4_wp, .true., w0, nxt, nyt)
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .true., 0.4_wp, .true., w1, nxt, nyt)
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .false., 0.4_wp, .false., l0, nxt, nyt)
         call run_gm(6, 2.0e-5_wp, 1.0e-5_wp, 100.0_wp, .true., 0.4_wp, .false., l1, nxt, nyt)
         rel = maxval(abs(w1 - w0))/maxval(abs(w0))
         call check(error, rel > 1.0e-3_wp, "fgnv: the work source must respond to FGNV")
         if (allocated(error)) exit checks
         call check(error, all(l1 == l0), "fgnv: the local source must not see FGNV")
      end block checks
   end subroutine test_fgnv_sensitive

   ! ------------------------------------------------------------------
   ! Test 7: config
   ! ------------------------------------------------------------------
   function validates(extra) result(ok)
      character(len=*), intent(in) :: extra
      logical :: ok
      type(config_t) :: cfg
      integer :: ierr
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
      call validate_config(cfg, ierr=ierr)
      ok = (ierr == OCEAN_STATUS_OK)
   end function validates

   subroutine test_config(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: SLO = "&ocean_slopes_nml enable = .true. /"
      character(len=*), parameter :: GM = "&ocean_gm_nml enable = .true., khth = 500.0, "
      checks: block
         call check(error, validates(SLO//GM//'gm_src_form = "work" /'), &
                    "config: gm_src_form='work' with GM must pass")
         if (allocated(error)) exit checks
         call check(error, validates(SLO//GM//'gm_src_form = "local" /'), &
                    "config: gm_src_form='local' must pass")
         if (allocated(error)) exit checks
         ! An unknown spelling is rejected by the schema's enum parse
         ! (ERROR STOP at read time, not a recoverable validate_config error).
         call check(error,.not. validates('&ocean_gm_nml enable = .false., gm_src_form = "work" /'), &
                    "config: gm_src_form='work' without GM must be refused")
      end block checks
   end subroutine test_config

end module test_ocean_gm_work
