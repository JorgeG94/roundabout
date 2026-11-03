!! Tests for the isolated balanced-vortex IC (`ic_config = "balanced_vortex"`).
!!
!! Static (no time stepping), against the analytic construction:
!!   * `phi_k = A_k G(r)`, `G = exp(-r^2/R^2)`, `A_k = (g alpha dT/rho0)(S_all - S_k)`
!!     (hydrostatic anomaly, zero bottom pressure); the face velocity is the
!!     geostrophic `u = -phi_y/f`, `v = phi_x/f` with the LOCAL f of the face.
!!     The scheme (4-cell corner average + centred difference) is second
!!     order, so the error is checked to be <= TOL_CONST (dx/R)^2 of the peak
!!     speed AND to fall by >= 3.5x when dx is halved.
!!   * thermal wind between neighbouring layers, from the stored T;
!!   * zero bottom pressure anomaly from the stored SSH + T;
!!   * the swirl sense (warm core => anticyclone for f > 0, cold core =>
!!     cyclone);
!!   * the gradient-wind quadratic `v^2/r + f v = d(phi)/dr` and its fail-loud
!!     on a too-intense anticyclone.
!!
!! Dynamic (the full split-explicit core, CPU or GPU through the C API):
!!   * f-plane steadiness: the balanced vortex keeps its centre (< 1 cell),
!!     surface KE (a few %) and axisymmetry; the same density with zero
!!     velocity (`vortex_balance = "none"`) does NOT (its KE change is
!!     reported -- this is what makes the first test non-vacuous);
!!   * beta drift (McWilliams & Flierl 1979): northern-hemisphere cyclones
!!     move north-west, anticyclones south-west.
module test_ocean_balanced_vortex
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, GRAVITY
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, &
                            rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_balanced_vortex_tests

   integer, parameter :: NGHOST = 4
   integer, parameter :: NZ = 10
   real(wp), parameter :: F0 = 1.0e-4_wp
   real(wp), parameter :: ALPHA = 0.17_wp
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: HDEPTH = 2000.0_wp
   real(wp), parameter :: HV = 500.0_wp
   real(wp), parameter :: TOL_CONST = 1.75_wp
      !! Error budget.  Per Fourier mode k on a square grid the 4-cell corner
      !! average has symbol cos(k dx/2)^2 ~ 1 - k^2 dx^2/4 (both axes) and the
      !! centred difference sinc(k dx/2) ~ 1 - k^2 dx^2/24, so the relative error
      !! of the gradient is O((7/24)(k dx)^2) = O((dx/R)^2) with k^2 R^2 <~ 6 on
      !! the Gaussian's peak-speed ring: <= 1.75 (dx/R)^2.  Measured 0.99 (dx/R)^2,
      !! and a 4.0x drop on halving dx (second order).
   ! dynamic cases
   integer, parameter :: DNX = 32
   real(wp), parameter :: DDX = 10000.0_wp, DR = 40000.0_wp
   real(wp), parameter :: DDT = 1200.0_wp

contains

   subroutine collect_ocean_balanced_vortex_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vortex_geostrophic_velocity_second_order", test_velocity_order), &
                  new_unittest("vortex_thermal_wind_and_zero_bottom_pressure", test_thermal_wind), &
                  new_unittest("vortex_swirl_sense", test_swirl_sense), &
                  new_unittest("vortex_gradient_wind_quadratic_and_fail_loud", test_gradient_wind), &
                  new_unittest("vortex_refuses_statically_unstable_cold_core", test_unstable_refused), &
                  new_unittest("vortex_fplane_balanced_steady", test_fplane_steady), &
                  new_unittest("vortex_beta_drift_sign", test_beta_drift) &
                  ]
   end subroutine collect_ocean_balanced_vortex_tests

   ! ------------------------------------------------------------------
   ! static IC helpers
   ! ------------------------------------------------------------------

   subroutine setup_ic(cfg, grid, state, nxp, dx, dT, balance, beta, ierr, dtdz)
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state
      integer, intent(in) :: nxp
      real(wp), intent(in) :: dx, dT, beta
      character(len=*), intent(in) :: balance
      integer, intent(out), optional :: ierr
      real(wp), intent(in), optional :: dtdz

      cfg%sim_type = "ocean"
      cfg%ocean%topo%topo_config = "flat"
      cfg%ocean%topo%max_depth = HDEPTH
      cfg%ocean%ic%ic_config = "balanced_vortex"
      cfg%coriolis_f = F0
      cfg%ocean%topo%coriolis_beta = beta
      cfg%ocean%topo%coriolis_y_ref = 0.5_wp*real(nxp, wp)*dx
      cfg%ocean%ic%alpha_T = ALPHA
      cfg%ocean%ic%rho_0 = RHO0
      cfg%ocean%ic%vortex_radius = 8.0_wp*5000.0_wp
      cfg%ocean%ic%vortex_hv = HV
      cfg%ocean%ic%vortex_dT = dT
      cfg%ocean%ic%vortex_dT_dz = 0.03_wp   ! statically stable for every |dT| <= 12 K here
      if (present(dtdz)) cfg%ocean%ic%vortex_dT_dz = dtdz
      cfg%ocean%ic%vortex_balance = balance
      cfg%initial_salinity = 35.0_wp
      cfg%initial_temperature = 10.0_wp
      cfg%T_init_surface = 0.0_wp
      cfg%T_init_bottom = 0.0_wp

      call grid%init(nxp, nxp, NGHOST, dx, dx)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      state%eos%alpha_T = cfg%ocean%ic%alpha_T
      state%eos%rho0 = cfg%ocean%ic%rho_0
      call ocean_state_seed_from_cfg(state, grid, cfg, ierr=ierr)
   end subroutine setup_ic

   pure function layer_amplitude(k, dT) result(a)
      !! A_k = (g alpha dT / rho0)(S_all - S_k), S_k = sum_{m>k} E_m dz + E_k dz/2.
      integer, intent(in) :: k
      real(wp), intent(in) :: dT
      real(wp) :: a
      real(wp) :: dz, s_all, s_k, d
      integer :: m
      dz = HDEPTH/real(NZ, wp)
      s_all = 0.0_wp
      s_k = 0.0_wp
      do m = 1, NZ
         d = (real(NZ - m, wp) + 0.5_wp)*dz
         s_all = s_all + exp(-(d/HV)**2)*dz
         if (m > k) s_k = s_k + exp(-(d/HV)**2)*dz
         if (m == k) s_k = s_k + 0.5_wp*exp(-(d/HV)**2)*dz
      end do
      a = GRAVITY*ALPHA*dT/RHO0*(s_all - s_k)
   end function layer_amplitude

   subroutine velocity_error(nxp, dx, dT, umax, err_u, err_v)
      !! Max |u - u_analytic| / max|u_analytic| over the interior faces and
      !! all layers (same for v), at grid spacing dx with R fixed.
      integer, intent(in) :: nxp
      real(wp), intent(in) :: dx, dT
      real(wp), intent(out) :: umax, err_u, err_v
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: x0, y0, R, f, xx, yy, ga, ua, va, du, dv
      integer :: i, j, k, ng

      call setup_ic(cfg, grid, state, nxp, dx, dT, "geostrophic", 0.0_wp)
      ng = NGHOST
      R = cfg%ocean%ic%vortex_radius
      x0 = 0.5_wp*real(nxp, wp)*dx
      y0 = x0
      f = F0
      umax = 0.0_wp
      du = 0.0_wp
      dv = 0.0_wp
      do k = 1, NZ
         do j = ng + 2, ng + nxp - 1
            do i = ng + 2, ng + nxp - 1
               xx = real(i - 1 - ng, wp)*dx
               yy = (real(j - ng, wp) - 0.5_wp)*dx
               ga = exp(-((xx - x0)**2 + (yy - y0)**2)/R**2)
               ua = 2.0_wp*layer_amplitude(k, dT)*(yy - y0)*ga/(f*R**2)
               umax = max(umax, abs(ua))
               du = max(du, abs(state%multilayer%u_face_x_layer(i, j, k) - ua))
               xx = (real(i - ng, wp) - 0.5_wp)*dx
               yy = real(j - 1 - ng, wp)*dx
               ga = exp(-((xx - x0)**2 + (yy - y0)**2)/R**2)
               va = -2.0_wp*layer_amplitude(k, dT)*(xx - x0)*ga/(f*R**2)
               dv = max(dv, abs(state%multilayer%v_face_y_layer(i, j, k) - va))
            end do
         end do
      end do
      err_u = du/umax
      err_v = dv/umax
      call state%destroy()
   end subroutine velocity_error

   subroutine test_velocity_order(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: umax1, eu1, ev1, umax2, eu2, ev2, ratio, tol
      character(len=200) :: msg

      ! Same physical vortex (R = 40 km, domain 480 km) at dx = 5 km, 2.5 km.
      call velocity_error(96, 5000.0_wp, 1.0_wp, umax1, eu1, ev1)
      call velocity_error(192, 2500.0_wp, 1.0_wp, umax2, eu2, ev2)
      write (*, '(a,es10.3,a,es10.3,a,es10.3)') &
         "  [vortex] peak |u| (m/s) ", umax1, "; rel err (u, v) at dx=5km: ", eu1, ", ", ev1
      write (*, '(a,es10.3,a,es10.3)') "  [vortex] rel err (u, v) at dx=2.5km: ", eu2, ", ", ev2
      tol = TOL_CONST*(5000.0_wp/40000.0_wp)**2
      write (msg, '(a,es10.3,a,es10.3)') "u error ", eu1, " exceeds TOL_CONST (dx/R)^2 = ", tol
      call check(error, eu1 <= tol .and. ev1 <= tol, trim(msg))
      if (allocated(error)) return
      ratio = max(eu1, ev1)/max(eu2, ev2)
      write (msg, '(a,f6.2)') "halving dx must cut the error ~4x (second order); got ratio ", ratio
      call check(error, ratio > 3.5_wp, trim(msg))
      if (allocated(error)) return
      call check(error, umax1 > 0.05_wp .and. umax1 < 1.0_wp, "peak speed is a sane mesoscale value")
   end subroutine test_velocity_order

   subroutine test_thermal_wind(error)
      !! du/dz = (g/(f rho0)) d(rho)/dy between layers, from the analytic T
      !! profile; and the bottom pressure anomaly from the STORED eta + T.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer, parameter :: NXP = 96
      real(wp), parameter :: DX = 5000.0_wp
      real(wp) :: R, x0, y0, xx, yy, ga, dz, ek, ek1, shear_an, shear, emax, smax
      real(wp) :: pbot, pscale, T_bg, rhop, hk, eta
      integer :: i, j, k, ng, idx_T
      character(len=200) :: msg

      call setup_ic(cfg, grid, state, NXP, DX, 1.0_wp, "geostrophic", 0.0_wp)
      ng = NGHOST
      R = cfg%ocean%ic%vortex_radius
      x0 = 0.5_wp*real(NXP, wp)*DX
      y0 = x0
      dz = HDEPTH/real(NZ, wp)
      idx_T = state%multilayer%idx_temperature

      emax = 0.0_wp
      smax = 0.0_wp
      do k = 1, NZ - 1
         ek = exp(-((real(NZ - k, wp) + 0.5_wp)*dz/HV)**2)
         ek1 = exp(-((real(NZ - k - 1, wp) + 0.5_wp)*dz/HV)**2)
         do j = ng + 2, ng + NXP - 1
            do i = ng + 2, ng + NXP - 1
               xx = real(i - 1 - ng, wp)*DX
               yy = (real(j - ng, wp) - 0.5_wp)*DX
               ga = exp(-((xx - x0)**2 + (yy - y0)**2)/R**2)
               shear_an = GRAVITY*ALPHA*1.0_wp/(F0*RHO0)*0.5_wp*(ek + ek1)*2.0_wp*(yy - y0)*ga/R**2
               shear = (state%multilayer%u_face_x_layer(i, j, k + 1) &
                        - state%multilayer%u_face_x_layer(i, j, k))/dz
               emax = max(emax, abs(shear - shear_an))
               smax = max(smax, abs(shear_an))
            end do
         end do
      end do
      write (msg, '(a,es10.3,a,es10.3)') "thermal-wind shear error ", emax/smax, &
         " exceeds TOL_CONST (dx/R)^2 = ", TOL_CONST*(DX/R)**2
      call check(error, emax/smax <= TOL_CONST*(DX/R)**2, trim(msg))
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! bottom pressure anomaly at the vortex centre column
      i = ng + NXP/2
      j = ng + NXP/2
      eta = sum(state%multilayer%h_layer(i, j, :)) - state%barotropic%b(i, j)
      pbot = GRAVITY*eta
      pscale = abs(pbot)
      do k = 1, NZ
         T_bg = cfg%ocean%ic%vortex_T_surf - cfg%ocean%ic%vortex_dT_dz*(real(NZ - k, wp) + 0.5_wp)*dz
         hk = state%multilayer%h_layer(i, j, k)
         rhop = -ALPHA*(state%multilayer%tracers(idx_T)%hTr(i, j, k)/hk - T_bg)
         pbot = pbot + GRAVITY*rhop*hk/RHO0
      end do
      write (msg, '(a,es10.3,a,es10.3)') "bottom pressure anomaly ", pbot, " vs g*eta ", GRAVITY*eta
      call check(error, eta > 1.0e-3_wp .and. abs(pbot) < 1.0e-3_wp*pscale, trim(msg))
      call state%destroy()
   end subroutine test_thermal_wind

   subroutine test_swirl_sense(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer, parameter :: NXP = 48
      real(wp), parameter :: DX = 5000.0_wp
      integer :: ic, jc, ng, s
      real(wp) :: zeta
      real(wp), parameter :: SGN(2) = [1.0_wp, -1.0_wp]

      ng = NGHOST
      do s = 1, 2
         ! s=1: warm core (dT>0) => anticyclone (zeta < 0); s=2: cold => cyclone
         call setup_ic(cfg, grid, state, NXP, DX, SGN(s), "geostrophic", 0.0_wp)
         ic = ng + NXP/2 + 1  ! corner at the vortex centre (x0 = 24 dx)
         jc = ng + NXP/2 + 1
         zeta = (state%multilayer%v_face_y_layer(ic, jc, NZ) &
                 - state%multilayer%v_face_y_layer(ic - 1, jc, NZ))/DX &
                - (state%multilayer%u_face_x_layer(ic, jc, NZ) &
                   - state%multilayer%u_face_x_layer(ic, jc - 1, NZ))/DX
         call check(error, SGN(s)*zeta < 0.0_wp, &
                    "warm core must be anticyclonic (zeta<0), cold core cyclonic, for f>0")
         if (allocated(error)) then
            call state%destroy()
            return
         end if
         ! warm core: mean flow at the surface is clockwise => v<0 east of centre
         call check(error, SGN(s)*state%multilayer%v_face_y_layer(ic + 3, jc, NZ) < 0.0_wp .and. &
                    SGN(s)*state%multilayer%u_face_x_layer(ic, jc + 3, NZ) > 0.0_wp, &
                    "warm core: v<0 east of centre, u>0 north of it (clockwise)")
         call state%destroy()
         if (allocated(error)) return
      end do
   end subroutine test_swirl_sense

   subroutine test_gradient_wind(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer, parameter :: NXP = 96
      real(wp), parameter :: DX = 5000.0_wp
      real(wp) :: R, x0, y0, xx, yy, ga, gr, res, rmax, vgeo, vgrad, rad, grpk
      integer :: i, j, ng, s, ierr, m
      real(wp), parameter :: SGN(2) = [-1.0_wp, 1.0_wp]
      character(len=200) :: msg

      ng = NGHOST
      do s = 1, 2
         call setup_ic(cfg, grid, state, NXP, DX, 2.0_wp*SGN(s), "gradient", 0.0_wp)
         R = cfg%ocean%ic%vortex_radius
         x0 = 0.5_wp*real(NXP, wp)*DX
         y0 = x0
         j = ng + NXP/2 + 1   ! v-face row at y = y0 exactly
         rmax = 0.0_wp
         grpk = 0.0_wp
         do m = 1, 8     ! rad = (m - 1/2) dx up to ~R: the vortex core
            i = ng + NXP/2 + m   ! v-face centred (m - 1/2) dx east of x0
            xx = (real(i - ng, wp) - 0.5_wp)*DX
            yy = real(j - 1 - ng, wp)*DX
            rad = xx - x0
            ga = exp(-((xx - x0)**2 + (yy - y0)**2)/R**2)
            gr = layer_amplitude(NZ, 2.0_wp*SGN(s))*ga*(-2.0_wp*rad/R**2)
            vgrad = state%multilayer%v_face_y_layer(i, j, NZ)
            res = abs(vgrad*vgrad/rad + F0*vgrad - gr)
            rmax = max(rmax, res)
            grpk = max(grpk, abs(gr))
         end do
         rmax = rmax/grpk
         write (msg, '(a,i0,a,es10.3)') "gradient-wind residual (dT sign ", nint(SGN(s)), &
            ") relative to d(phi)/dr: ", rmax
         call check(error, rmax < 2.0_wp*TOL_CONST*(DX/R)**2, trim(msg))
         if (allocated(error)) then
            call state%destroy()
            return
         end if
         ! gradient-wind speed vs geostrophic: a cyclone is slower, an
         ! anticyclone faster, than the geostrophic estimate at the same point.
         i = ng + NXP/2 + 6
         xx = (real(i - ng, wp) - 0.5_wp)*DX
         rad = xx - x0
         ga = exp(-(rad*rad)/R**2)
         vgeo = -2.0_wp*layer_amplitude(NZ, 2.0_wp*SGN(s))*rad*ga/(F0*R**2)
         vgrad = state%multilayer%v_face_y_layer(i, j, NZ)
         if (SGN(s) > 0.0_wp) then
            call check(error, abs(vgrad) > abs(vgeo), "anticyclone: gradient-wind speed exceeds geostrophic")
         else
            call check(error, abs(vgrad) < abs(vgeo), "cyclone: gradient-wind speed is below geostrophic")
         end if
         call state%destroy()
         if (allocated(error)) return
      end do

      ! fail loud: a 12 degC warm core has 4 G_r/(f^2 r) < -1 near r ~ R/sqrt(2)... no real root.
      call setup_ic(cfg, grid, state, NXP, DX, 12.0_wp, "gradient", 0.0_wp, ierr)
      call check(error, ierr /= 0, "gradient wind with no real root must fail loud (ierr /= 0)")
      call state%destroy()
      if (allocated(error)) return
      ! the same strength as a cyclone, and as a geostrophic anticyclone, is fine
      call setup_ic(cfg, grid, state, NXP, DX, -12.0_wp, "gradient", 0.0_wp, ierr)
      call check(error, ierr == 0, "a strong cyclone always has a gradient-wind root")
      call state%destroy()
      if (allocated(error)) return
      call setup_ic(cfg, grid, state, NXP, DX, 12.0_wp, "geostrophic", 0.0_wp, ierr)
      call check(error, ierr == 0, "geostrophic balance has no root condition")
      call state%destroy()
   end subroutine test_gradient_wind

   subroutine test_unstable_refused(error)
      !! A cold core warms with depth by up to sqrt(2) exp(-1/2) |dT|/Hv
      !! (K/m, at d = Hv/sqrt(2)); a weaker background dT/dz gives N^2 < 0.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr
      real(wp) :: dtdz_min

      dtdz_min = sqrt(2.0_wp)*exp(-0.5_wp)*2.0_wp/HV   ! dT = -2 K, Hv = 500 m: 3.4e-3
      call setup_ic(cfg, grid, state, 48, 5000.0_wp, -2.0_wp, "geostrophic", 0.0_wp, ierr, 0.9_wp*dtdz_min)
      call check(error, ierr /= 0, "cold core with dT_dz below 0.858|dT|/Hv must be refused")
      call state%destroy()
      if (allocated(error)) return
      call setup_ic(cfg, grid, state, 48, 5000.0_wp, -2.0_wp, "geostrophic", 0.0_wp, ierr, 1.3_wp*dtdz_min)
      call check(error, ierr == 0, "the same core with margin above the minimum is accepted")
      call state%destroy()
      if (allocated(error)) return
      call setup_ic(cfg, grid, state, 48, 5000.0_wp, 2.0_wp, "geostrophic", 0.0_wp, ierr, 1.0e-4_wp)
      call check(error, ierr == 0, "a warm core only adds stratification: weak background is fine")
      call state%destroy()
   end subroutine test_unstable_refused

   ! ------------------------------------------------------------------
   ! dynamic cases
   ! ------------------------------------------------------------------

   function vortex_nml(dT, beta, balance, dtdz) result(nml)
      real(wp), intent(in) :: dT, beta, dtdz
      character(len=*), intent(in) :: balance
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      character(len=32) :: s_dt, s_beta, s_yref, s_dtdz
      write (s_dtdz, '(es14.6)') dtdz
      write (s_dt, '(es14.6)') dT
      write (s_beta, '(es14.6)') beta
      write (s_yref, '(es14.6)') 0.5_wp*real(DNX, wp)*DDX
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 32, ny = 32, dx = 10000.0, dy = 10000.0, nghost = 4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 6 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 1200.0 /"//nl// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//nl// &
            "&tracer_nml initial_salinity = 35.0 /"//nl// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 2000.0, "// &
            "wind_config = 'constant', coriolis_beta = "//trim(adjustl(s_beta))// &
            ", coriolis_y_ref = "//trim(adjustl(s_yref))//" /"//nl// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//nl// &
            "&ocean_eos_nml eos = 'linear' /"//nl// &
            "&ocean_ic_nml ic_config = 'balanced_vortex', alpha_T = 0.17, rho_0 = 1035.0, "// &
     "vortex_radius = 40000.0, vortex_hv = 500.0, vortex_dT_dz = "//trim(adjustl(s_dtdz))//", vortex_dT = "//trim(adjustl(s_dt))// &
            ", vortex_balance = '"//balance//"' /"//nl// &
            "&ocean_hvisc_nml nu_h = 350.0 /"//nl// &
            "&ocean_vmix_nml use_closure = .false. /"//nl// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .true. /"//nl// &
            "&vcoord_nml vcoord_type = 'sigma' /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function vortex_nml

   subroutine surface_diag(handle, sgn, xbar, ybar, aniso, ke, finite)
      !! Surface-layer diagnostics: KE sum (m^2/s^2, per-cell sum), and the
      !! same-sign vorticity-weighted centroid (xbar, ybar; m) and
      !! anisotropy |Ixx - Iyy| / (Ixx + Iyy) of the weight distribution,
      !! over a window of 3R about the initial centre.  `sgn` = +1 selects
      !! cyclonic (zeta > 0) weights, -1 anticyclonic.
      type(c_ptr), intent(in) :: handle
      real(wp), intent(in) :: sgn
      real(wp), intent(out) :: xbar, ybar, aniso, ke
      logical, intent(out) :: finite
      type(c_ptr) :: ptr
      integer(c_int) :: status, nxp, nyp, nzp, ng, nx, ny, nz, gen
      real(wp), pointer :: u(:, :, :), v(:, :, :)
      real(wp) :: zeta, w, sw, sx, sy, xx, yy, uc, vc, ixx, iyy, x0, y0
      integer :: i, j

      status = rdb_ocean_refresh_host(handle)
      status = rdb_ocean_get_grid_info(handle, nxp, nyp, nzp, ng)
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, u, [nx, ny, nz])
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, v, [nx, ny, nz])
      x0 = 0.5_wp*real(nxp, wp)*DDX
      y0 = x0
      ke = 0.0_wp
      sw = 0.0_wp
      sx = 0.0_wp
      sy = 0.0_wp
      do j = ng + 1, ng + nyp
         do i = ng + 1, ng + nxp
            uc = 0.5_wp*(u(i, j, nzp) + u(i + 1, j, nzp))
            vc = 0.5_wp*(v(i, j, nzp) + v(i, j + 1, nzp))
            ke = ke + 0.5_wp*(uc*uc + vc*vc)
         end do
      end do
      do j = ng + 2, ng + nyp
         do i = ng + 2, ng + nxp
            xx = real(i - 1 - ng, wp)*DDX
            yy = real(j - 1 - ng, wp)*DDX
            if ((xx - x0)**2 + (yy - y0)**2 > (3.0_wp*DR)**2) cycle
            zeta = (v(i, j, nzp) - v(i - 1, j, nzp))/DDX - (u(i, j, nzp) - u(i, j - 1, nzp))/DDX
            w = max(0.0_wp, sgn*zeta)
            sw = sw + w
            sx = sx + w*xx
            sy = sy + w*yy
         end do
      end do
      finite = ieee_is_finite(ke) .and. ieee_is_finite(sw)
      if (.not. finite .or. sw <= 0.0_wp) then   ! no flow (unbalanced start): no centroid
         xbar = 0.0_wp
         ybar = 0.0_wp
         aniso = 0.0_wp
         return
      end if
      xbar = sx/sw
      ybar = sy/sw
      ixx = 0.0_wp
      iyy = 0.0_wp
      do j = ng + 2, ng + nyp
         do i = ng + 2, ng + nxp
            xx = real(i - 1 - ng, wp)*DDX
            yy = real(j - 1 - ng, wp)*DDX
            if ((xx - x0)**2 + (yy - y0)**2 > (3.0_wp*DR)**2) cycle
            zeta = (v(i, j, nzp) - v(i - 1, j, nzp))/DDX - (u(i, j, nzp) - u(i, j - 1, nzp))/DDX
            w = max(0.0_wp, sgn*zeta)
            ixx = ixx + w*(xx - xbar)**2
            iyy = iyy + w*(yy - ybar)**2
         end do
      end do
      aniso = abs(ixx - iyy)/(ixx + iyy)
   end subroutine surface_diag

   subroutine run_vortex(dT, beta, balance, n_steps, xbar0, ybar0, xbar1, ybar1, &
                         aniso0, aniso1, ke0, ke1, ok, dtdz)
      real(wp), intent(in) :: dT, beta, dtdz
      character(len=*), intent(in) :: balance
      integer, intent(in) :: n_steps
      real(wp), intent(out) :: xbar0, ybar0, xbar1, ybar1, aniso0, aniso1, ke0, ke1
      logical, intent(out) :: ok
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: nml
      real(wp) :: sgn
      logical :: fin0, fin1

      ok = .false.
      xbar0 = 0.0_wp; ybar0 = 0.0_wp; xbar1 = 0.0_wp; ybar1 = 0.0_wp
      aniso0 = 0.0_wp; aniso1 = 0.0_wp; ke0 = 0.0_wp; ke1 = 0.0_wp
      sgn = merge(-1.0_wp, 1.0_wp, dT > 0.0_wp)   ! warm core => anticyclone
      nml = vortex_nml(dT, beta, balance, dtdz)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return
      call surface_diag(handle, sgn, xbar0, ybar0, aniso0, ke0, fin0)
      status = rdb_ocean_step(handle, int(n_steps, c_int))
      if (status == OCEAN_STATUS_OK) call surface_diag(handle, sgn, xbar1, ybar1, aniso1, ke1, fin1)
      ok = (status == OCEAN_STATUS_OK) .and. fin0 .and. fin1
      status = rdb_ocean_destroy(handle)
   end subroutine run_vortex

   subroutine test_fplane_steady(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x0, y0, x1, y1, a0, a1, ke0, ke1, d, dke, dke_rest
      real(wp) :: ru0, rv0, ra0, ra1, rke0, rke1
      integer, parameter :: NSTEP = 144   ! 2 days at dt = 1200 s
      logical :: ok
      character(len=240) :: msg

      call run_vortex(1.0_wp, 0.0_wp, "geostrophic", NSTEP, x0, y0, x1, y1, a0, a1, ke0, ke1, ok, 2.0e-3_wp)
      call check(error, ok .and. x0 > 0.0_wp .and. x1 > 0.0_wp, &
                 "balanced f-plane vortex run must stay finite and keep a vorticity centroid")
      if (allocated(error)) return
      d = sqrt((x1 - x0)**2 + (y1 - y0)**2)/DDX
      dke = (ke1 - ke0)/ke0
      write (*, '(a,f7.3,a,f7.3,a,f8.5,a,f8.5)') &
         "  [vortex] f-plane 2 d, balanced: centre moved ", d, " cells; surface KE change ", &
         100.0_wp*dke, " %; aniso ", a0, " -> ", a1
      write (msg, '(a,f6.3,a)') "balanced vortex centre moved ", d, " cells (limit 1)"
      call check(error, d < 1.0_wp, trim(msg))
      if (allocated(error)) return
      write (msg, '(a,f7.3,a)') "balanced vortex surface KE changed ", 100.0_wp*dke, " % (limit 3 %)"
      call check(error, abs(dke) < 0.03_wp, trim(msg))
      if (allocated(error)) return
      write (msg, '(a,es10.3)') "balanced vortex lost axisymmetry: aniso = ", a1
      call check(error, a1 < 0.05_wp, trim(msg))
      if (allocated(error)) return

      ! Non-vacuity: identical density + SSH but zero velocity must NOT be steady.
      call run_vortex(1.0_wp, 0.0_wp, "none", NSTEP, ru0, rv0, x1, y1, ra0, ra1, rke0, rke1, ok, 2.0e-3_wp)
      call check(error, ok, "unbalanced run must execute and stay finite")
      if (allocated(error)) return
      dke_rest = (rke1 - ke0)/ke0
      write (*, '(a,es10.3,a,es10.3,a,f9.2,a)') &
         "  [vortex] unbalanced (zero-velocity) start: surface KE ", rke0, " -> ", rke1, &
         " (", 100.0_wp*dke_rest, " % of the balanced KE)"
      write (msg, '(a,f9.2,a)') "unbalanced start changed KE by only ", 100.0_wp*dke_rest, &
         " % -- the steadiness test would pass without balance"
      call check(error, abs(dke_rest) > 10.0_wp*abs(dke) .and. abs(dke_rest) > 0.2_wp, trim(msg))
   end subroutine test_fplane_steady

   subroutine test_beta_drift(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x0, y0, x1, y1, a0, a1, ke0, ke1
      integer, parameter :: NSTEP = 360   ! 5 days
      logical :: ok
      integer :: s
      real(wp), parameter :: DTS(2) = [-2.0_wp, 2.0_wp]   ! cold (cyclone), warm (anticyclone)
      character(len=240) :: msg

      do s = 1, 2
         call run_vortex(DTS(s), 1.2e-10_wp, "geostrophic", NSTEP, x0, y0, x1, y1, a0, a1, ke0, ke1, ok, 5.0e-3_wp)
         call check(error, ok .and. x0 > 0.0_wp .and. x1 > 0.0_wp, &
                    "beta-plane vortex run must stay finite and keep a vorticity centroid")
         if (allocated(error)) return
         write (*, '(a,a,a,f8.2,a,f8.2,a)') "  [vortex] beta 1.2e-10, 5 d, ", &
            merge("cyclone    ", "anticyclone", DTS(s) < 0.0_wp), ": displacement (dx, dy) = (", &
            (x1 - x0)/1000.0_wp, ", ", (y1 - y0)/1000.0_wp, ") km"
         write (msg, '(a,f8.2,a)') "vortex must drift WEST; dx = ", (x1 - x0)/1000.0_wp, " km"
         call check(error, x1 - x0 < -500.0_wp, trim(msg))
         if (allocated(error)) return
         if (DTS(s) < 0.0_wp) then
            write (msg, '(a,f8.2,a)') "cyclone must drift NORTH (NW); dy = ", (y1 - y0)/1000.0_wp, " km"
            call check(error, y1 - y0 > 1000.0_wp, trim(msg))
         else
            write (msg, '(a,f8.2,a)') "anticyclone must drift SOUTH (SW); dy = ", (y1 - y0)/1000.0_wp, " km"
            call check(error, y1 - y0 < -1000.0_wp, trim(msg))
         end if
         if (allocated(error)) return
      end do
   end subroutine test_beta_drift

end module test_ocean_balanced_vortex
