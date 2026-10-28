!! Biharmonic inside the stress tensor (`&ocean_hvisc_nml biharm_stress`).
!!
!! (a) Dispersion: u = sin(kx) on a uniform grid, uniform h, gives
!!     du = -(Kh*lam + Ah*lam^2)*u*h/(h+H_VANISHED), lam = (2-2cos(k dx))/dx^2,
!!     to 1e-12 (Ah alone, and Kh+Ah together).
!! (b) Solid-body rotation (power-of-two numbers) gives EXACTLY zero tendency
!!     in the interior, with a spatially varying Ah field.
!! (c) Interior momentum conservation of the biharmonic part with strongly
!!     varying h: sum(du*(h_u+eps)) = 0 for the stress form; the velocity-form
!!     add-on FAILS it (positive control).
!! (d) The Ah formula, |S|*(C_b*D^4 + D^6*f*|S|/(5V^2)) max'ed against
!!     max(AH, c_A*D^3), with |S| recomputed here (square-then-average), on a
!!     stretched dx/=dy grid with a varying f.
!! (e) The joint bound: Kh at its limit gives Ah = 0; Kh = Kh_max/2 gives
!!     Ah = Ah_max/2; Kh = 0 gives Ah_max = c_b D^4/(64 dt); stencil 1/D4 =
!!     D^4/32 on a uniform square grid.
!! (f) Knob off: no statics allocated and the stress + nu_4 tendency equals the Laplacian-only stress tendency plus the
!!     velocity-form biharmonic tendency.
module test_ocean_hvisc_biharm_stress
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_biharm_stress_tests

   real(wp), parameter :: TWO_PI = 6.28318530717958647692528676655900577_wp

contains

   subroutine collect_ocean_hvisc_biharm_stress_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sinusoid_dispersion_ah_and_kh", test_dispersion), &
                  new_unittest("solid_body_rotation_exactly_zero", test_solid_body), &
                  new_unittest("momentum_conservation_vs_velocity_form", test_conservation), &
                  new_unittest("ah_formula_smag_coriolis_floor", test_ah_formula), &
                  new_unittest("joint_bound_cases", test_joint_bound), &
                  new_unittest("knob_off_bit_identity", test_knob_off) &
                  ]
   end subroutine collect_ocean_hvisc_biharm_stress_tests

   subroutine make_state(grid, ms, nx, ny, dx, dy, h0)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dx, dy, h0
      call grid%init(nx, ny, 2, dx, dy)
      ms%nz_ml = 1
      call ms%init(grid)
      ms%h_layer = h0
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine make_state

   subroutine run_hv(grid, ms, hv, dt, want_bh, fc, ahv, cb, smag, cor, vcor, vary)
      !! Metrics, optional biharmonic statics (host, before the device map),
      !! device map, one tendency evaluation, readback, unmap.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      real(wp), intent(in) :: dt
      logical, intent(in) :: want_bh
      real(wp), intent(in) :: fc(:, :)
      real(wp), intent(in) :: ahv, cb, vcor
      logical, intent(in) :: smag, cor
      logical, intent(in), optional :: vary
      type(ocean_metrics_t) :: metrics
      integer :: i, j

      call make_cartesian_metrics(metrics, grid)
      if (want_bh) then
         call hv%build_biharm_stress(metrics, fc, ahv, cb, smag, cor, vcor)
         if (present(vary)) then
            if (vary) then
               do j = 1, size(hv%bh_bg_t, 2)
                  do i = 1, size(hv%bh_bg_t, 1)
                     hv%bh_bg_t(i, j) = hv%bh_bg_t(i, j)*(1.0_wp + 0.5_wp*sin(real(i + 2*j, wp)))
                  end do
               end do
               do j = 1, size(hv%bh_bg_q, 2)
                  do i = 1, size(hv%bh_bg_q, 1)
                     hv%bh_bg_q(i, j) = hv%bh_bg_q(i, j)*(1.0_wp + 0.5_wp*cos(real(2*i + j, wp)))
                  end do
               end do
            end if
         end if
      end if
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=dt)
      !$acc update self(hv%du_visc%data, hv%dv_visc%data, hv%ah4_t%data, hv%ah4_q%data)
      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_hv

   subroutine test_dispersion(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DX = 1024.0_wp, H0 = 100.0_wp, AH = 1.0e9_wp, KH = 1.0e3_wp, DT = 1.0_wp
      integer, parameter :: NX = 40, NY = 16
      integer :: i, j, nxt, nyt, pass
      real(wp) :: kw, lam, expct, obs, worst
      real(wp), allocatable :: fc(:, :)

      do pass = 1, 2
         call make_state(grid, ms, NX, NY, DX, DX, H0)
         nxt = grid%nx_total
         nyt = grid%ny_total
         kw = TWO_PI*3.0_wp/(real(nxt, wp)*DX)
         lam = (2.0_wp - 2.0_wp*cos(kw*DX))/(DX*DX)
         do j = 1, nyt
            do i = 1, nxt + 1
               ms%u_face_x_layer(i, j, 1) = sin(kw*real(i - 1, wp)*DX)
            end do
         end do
         allocate (fc(nxt + 1, nyt + 1), source=0.0_wp)
         call hv%init(grid, nz_ml=1)
         hv%stress_tensor = .true.
         hv%nu_4 = AH
         hv%nu_h = merge(KH, 0.0_wp, pass == 2)
         call run_hv(grid, ms, hv, DT, .true., fc, 0.0_wp, 0.0_wp, .false., .false., 6.0_wp)
         worst = 0.0_wp
         do j = 4, nyt - 3
            do i = 7, nxt - 5
               expct = -(hv%nu_h*lam + AH*lam*lam)*ms%u_face_x_layer(i, j, 1)*H0/(H0 + H_VANISHED)
               obs = hv%du_visc%data(i, j, 1)
               worst = max(worst, abs(obs - expct)/(AH*lam*lam))
            end do
         end do
         call check(error, worst <= 1.0e-12_wp, "sinusoid: tendency must equal -(Kh lam + Ah lam^2) u to 1e-12")
         call hv%destroy()
         call ms%destroy()
         deallocate (fc)
         if (allocated(error)) return
      end do
   end subroutine test_dispersion

   subroutine test_solid_body(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DX = 1024.0_wp, OM = 6.103515625e-05_wp   ! 2^-14
      integer, parameter :: NX = 30, NY = 30
      integer :: i, j, nxt, nyt
      real(wp) :: worst
      real(wp), allocatable :: fc(:, :)

      call make_state(grid, ms, NX, NY, DX, DX, 100.0_wp)
      nxt = grid%nx_total
      nyt = grid%ny_total
      do j = 1, nyt
         do i = 1, nxt + 1
            ms%u_face_x_layer(i, j, 1) = -OM*512.0_wp*real(2*j - 1, wp)
         end do
      end do
      do j = 1, nyt + 1
         do i = 1, nxt
            ms%v_face_y_layer(i, j, 1) = OM*512.0_wp*real(2*i - 1, wp)
         end do
      end do
      allocate (fc(nxt + 1, nyt + 1), source=1.0e-4_wp)
      call hv%init(grid, nz_ml=1)
      hv%stress_tensor = .true.
      hv%nu_4 = 1.0e8_wp
      hv%nu_h = 1.0e3_wp
      call run_hv(grid, ms, hv, 1.0_wp, .true., fc, 0.0_wp, 0.0_wp, .false., .false., 6.0_wp, vary=.true.)
      worst = max(maxval(abs(hv%du_visc%data(7:nxt - 5, 7:nyt - 6, :))), &
                  maxval(abs(hv%dv_visc%data(7:nxt - 6, 7:nyt - 5, :))))
      call check(error, worst == 0.0_wp, "solid-body rotation must give an exactly zero tendency (varying Ah)")
      call hv%destroy()
      call ms%destroy()
   end subroutine test_solid_body

   subroutine conserved_residual(hh, du, dv, nxt, nyt, ru, rv)
      !! sum(du*(h_u+eps)) / sum(|du|*(h_u+eps)) and the v twin, interior faces.
      integer, intent(in) :: nxt, nyt
      real(wp), intent(in) :: hh(nxt, nyt), du(nxt + 1, nyt), dv(nxt, nyt + 1)
      real(wp), intent(out) :: ru, rv
      real(wp) :: su, sa, sv, sb, w
      integer :: i, j
      su = 0.0_wp
      sa = 0.0_wp
      sv = 0.0_wp
      sb = 0.0_wp
      do j = 2, nyt - 1
         do i = 2, nxt
            w = 0.5_wp*(hh(i - 1, j) + hh(i, j)) + H_VANISHED
            su = su + du(i, j)*w
            sa = sa + abs(du(i, j))*w
         end do
      end do
      do j = 2, nyt
         do i = 2, nxt - 1
            w = 0.5_wp*(hh(i, j - 1) + hh(i, j)) + H_VANISHED
            sv = sv + dv(i, j)*w
            sb = sb + abs(dv(i, j))*w
         end do
      end do
      ru = abs(su)/max(sa, tiny(1.0_wp))
      rv = abs(sv)/max(sb, tiny(1.0_wp))
   end subroutine conserved_residual

   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DX = 1024.0_wp, AH = 1.0e8_wp, DT = 1.0_wp
      integer, parameter :: NX = 40, NY = 40
      integer :: i, j, nxt, nyt, pass
      real(wp) :: r2, ru, rv, ru_stress, rv_stress
      real(wp), allocatable :: fc(:, :), hh(:, :)

      do pass = 1, 2
         call make_state(grid, ms, NX, NY, DX, DX, 100.0_wp)
         nxt = grid%nx_total
         nyt = grid%ny_total
         allocate (hh(nxt, nyt))
         do j = 1, nyt
            do i = 1, nxt
               hh(i, j) = 100.0_wp*(1.0_wp + 0.5_wp*sin(0.37_wp*real(i, wp) + 0.21_wp*real(j, wp)) &
                                    *cos(0.13_wp*real(i - j, wp)))
               ms%h_layer(i, j, 1) = hh(i, j)
            end do
         end do
         do j = 1, nyt
            do i = 1, nxt + 1
               r2 = (real(i, wp) - 22.0_wp)**2 + (real(j, wp) - 21.5_wp)**2
               if (r2 < 100.0_wp) ms%u_face_x_layer(i, j, 1) = exp(-r2/18.0_wp)*(1.0_wp + 0.3_wp*sin(real(i, wp)))
            end do
         end do
         do j = 1, nyt + 1
            do i = 1, nxt
               r2 = (real(i, wp) - 21.5_wp)**2 + (real(j, wp) - 22.0_wp)**2
               if (r2 < 100.0_wp) ms%v_face_y_layer(i, j, 1) = -0.7_wp*exp(-r2/14.0_wp)*(1.0_wp + 0.2_wp*cos(real(j, wp)))
            end do
         end do
         allocate (fc(nxt + 1, nyt + 1), source=0.0_wp)
         call hv%init(grid, nz_ml=1)
         hv%stress_tensor = .true.
         hv%nu_4 = AH
         hv%nu_h = 0.0_wp
         ! pass 1: stress-form biharmonic; pass 2: velocity-form add-on.
         call run_hv(grid, ms, hv, DT, pass == 1, fc, 0.0_wp, 0.0_wp, .false., .false., 6.0_wp)
         call conserved_residual(hh, hv%du_visc%data(:, :, 1), hv%dv_visc%data(:, :, 1), nxt, nyt, ru, rv)
         if (pass == 1) then
            ru_stress = ru
            rv_stress = rv
            call check(error, maxval(abs(hv%du_visc%data)) > 0.0_wp, "bump must be damped (nonzero tendency)")
         end if
         call hv%destroy()
         call ms%destroy()
         deallocate (fc, hh)
         if (allocated(error)) return
      end do
      call check(error, ru_stress <= 1.0e-12_wp .and. rv_stress <= 1.0e-12_wp, &
                 "stress-form biharmonic must conserve interior momentum sum(h du)")
      if (allocated(error)) return
      call check(error, ru > 1.0e-4_wp .and. rv > 1.0e-4_wp, &
                 "positive control: the velocity-form biharmonic must FAIL the momentum-conservation test")
   end subroutine test_conservation

   subroutine test_ah_formula(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DX = 40000.0_wp, DY = 20000.0_wp, CB = 0.06_wp, CA = 0.08_wp, V = 6.0_wp
      real(wp), parameter :: AH0 = 1.0e3_wp
      integer, parameter :: NX = 24, NY = 24
      integer :: i, j, ii, jj, nxt, nyt, nwin_s
      real(wp) :: d2, d3, e, s(0:1, 0:1), sm, fm, expct, worst, ee(0:1, 0:1), sq, fq
      real(wp), allocatable :: fc(:, :)

      call make_state(grid, ms, NX, NY, DX, DY, 500.0_wp)
      nxt = grid%nx_total
      nyt = grid%ny_total
      do j = 1, nyt
         do i = 1, nxt + 1
            ms%u_face_x_layer(i, j, 1) = 3.0_wp*sin(0.5_wp*real(i, wp) + 0.3_wp*real(j, wp))
         end do
      end do
      do j = 1, nyt + 1
         do i = 1, nxt
            ms%v_face_y_layer(i, j, 1) = 2.5_wp*cos(0.4_wp*real(i, wp) - 0.6_wp*real(j, wp))
         end do
      end do
      allocate (fc(nxt + 1, nyt + 1))
      do j = 1, nyt + 1
         do i = 1, nxt + 1
            fc(i, j) = 1.0e-4_wp + 1.0e-6_wp*real(j, wp) - 2.0e-7_wp*real(i, wp)
         end do
      end do
      call hv%init(grid, nz_ml=1)
      hv%stress_tensor = .true.
      hv%nu_4 = AH0
      call run_hv(grid, ms, hv, 1.0_wp, .true., fc, CA, CB, .true., .true., V)

      d2 = 2.0_wp*DX*DX*DY*DY/(DX*DX + DY*DY)
      d3 = d2*sqrt(d2)
      worst = 0.0_wp
      nwin_s = 0
      ! T cells 4..nxt-3: Ah_T = max(bg, |S|(Cb D^4 + D^6 fmax |S|/(5V^2)))
      do j = 4, nyt - 3
         do i = 4, nxt - 3
            e = tension(i, j)
            sm = 0.0_wp
            do jj = 0, 1
               do ii = 0, 1
                  s(ii, jj) = shear(i + ii, j + jj)
               end do
            end do
            sm = sqrt(e*e + 0.25_wp*(sum(s*s)))
            fm = max(abs(fc(i, j)), abs(fc(i + 1, j)), abs(fc(i, j + 1)), abs(fc(i + 1, j + 1)))
            expct = max(max(AH0, CA*d3), sm*(CB*d2*d2 + d2*d2*d2*fm*sm/(5.0_wp*V*V)))
            if (sm*(CB*d2*d2 + d2*d2*d2*fm*sm/(5.0_wp*V*V)) > max(AH0, CA*d3)) nwin_s = nwin_s + 1
            worst = max(worst, abs(hv%ah4_t%data(i, j, 1) - expct)/expct)
         end do
      end do
      ! corners 5..nxt-3
      do j = 5, nyt - 3
         do i = 5, nxt - 3
            sq = shear(i, j)
            do jj = 0, 1
               do ii = 0, 1
                  ee(ii, jj) = tension(i - 1 + ii, j - 1 + jj)
               end do
            end do
            sm = sqrt(sq*sq + 0.25_wp*(sum(ee*ee)))
            fq = abs(fc(i, j))
            expct = max(max(AH0, CA*d3), sm*(CB*d2*d2 + d2*d2*d2*fq*sm/(5.0_wp*V*V)))
            worst = max(worst, abs(hv%ah4_q%data(i, j, 1) - expct)/expct)
         end do
      end do
      call check(error, worst <= 1.0e-12_wp, "Ah at T and Bu must equal max(bg, |S|(Cb D^4 + D^6 f |S|/(5V^2)))")
      if (allocated(error)) return
      call check(error, nwin_s > 20 .and. nwin_s < (nxt - 6)*(nyt - 6), &
                 "both the strain term and the floor must win somewhere (branches exercised)")
      call hv%destroy()
      call ms%destroy()

   contains

      function tension(a, b) result(t)
         integer, intent(in) :: a, b
         real(wp) :: t
         t = (ms%u_face_x_layer(a + 1, b, 1) - ms%u_face_x_layer(a, b, 1))/DX &
             - (ms%v_face_y_layer(a, b + 1, 1) - ms%v_face_y_layer(a, b, 1))/DY
      end function tension

      function shear(a, b) result(t)
         integer, intent(in) :: a, b
         real(wp) :: t
         t = (ms%v_face_y_layer(a, b, 1) - ms%v_face_y_layer(a - 1, b, 1))/DX &
             + (ms%u_face_x_layer(a, b, 1) - ms%u_face_x_layer(a, b - 1, 1))/DY
      end function shear

   end subroutine test_ah_formula

   subroutine joint_case(nu_h_in, ah_in, dt, ah_t, mx_t, mx_q)
      real(wp), intent(in) :: nu_h_in, ah_in, dt
      real(wp), allocatable, intent(out) :: ah_t(:, :), mx_t(:, :), mx_q(:, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), allocatable :: fc(:, :)
      call make_state(grid, ms, 20, 20, 1000.0_wp, 1000.0_wp, 100.0_wp)
      allocate (fc(grid%nx_total + 1, grid%ny_total + 1), source=0.0_wp)
      call hv%init(grid, nz_ml=1)
      hv%stress_tensor = .true.
      hv%nu_h = nu_h_in
      hv%nu_4 = ah_in
      call run_hv(grid, ms, hv, dt, .true., fc, 0.0_wp, 0.0_wp, .false., .false., 6.0_wp)
      allocate (ah_t, source=hv%ah4_t%data(:, :, 1))
      allocate (mx_t, source=hv%bh_mx_t)
      allocate (mx_q, source=hv%bh_mx_q)
      call hv%destroy()
      call ms%destroy()
   end subroutine joint_case

   subroutine test_joint_bound(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: D = 1000.0_wp, DT = 2.0_wp, CBND = 0.8_wp
      real(wp), allocatable :: a0(:, :), a1(:, :), a2(:, :), a3(:, :), mx_t(:, :), mx_q(:, :), m2(:, :), m3(:, :)
      real(wp) :: ah_max, kh_max
      integer :: n

      ah_max = CBND*D**4/(64.0_wp*DT)
      kh_max = CBND*D**2/(8.0_wp*DT)
      ! Kh = 0: Ah = min(huge floor, Ah_max) = Ah_max, and 1/D4 = D^4/32.
      call joint_case(0.0_wp, 100.0_wp*ah_max, DT, a0, mx_t, mx_q)
      n = size(a0, 1)
      call check(error, abs(mx_t(n/2, n/2) - D**4/32.0_wp) <= 1.0e-14_wp*D**4/32.0_wp .and. &
                 abs(mx_q(n/2, n/2) - D**4/32.0_wp) <= 1.0e-14_wp*D**4/32.0_wp, &
                 "uniform square grid: stencil 1/D4 = D^4/32 at T and Bu")
      if (allocated(error)) return
      call check(error, abs(a0(n/2, n/2) - ah_max) <= 1.0e-13_wp*ah_max, &
                 "Kh=0: Ah is capped at Ah_max = c_b D^4/(64 dt)")
      if (allocated(error)) return
      ! Kh at its limit: Ah = 0.
      call joint_case(1.0e3_wp*kh_max, 100.0_wp*ah_max, DT, a1, m2, m3)
      call check(error, maxval(abs(a1(5:n - 4, 5:n - 4))) <= 1.0e-12_wp*ah_max, &
                 "Kh = Kh_max: the Laplacian spent the budget, Ah = 0")
      if (allocated(error)) return
      ! Kh = Kh_max / 2: Ah = Ah_max / 2.
      call joint_case(0.5_wp*kh_max, 100.0_wp*ah_max, DT, a2, m2, m3)
      call check(error, abs(a2(n/2, n/2) - 0.5_wp*ah_max) <= 1.0e-13_wp*ah_max, &
                 "Kh = Kh_max/2: Ah = Ah_max/2")
      if (allocated(error)) return
      ! Small floor under the bound is untouched.
      call joint_case(0.5_wp*kh_max, 1.0e-3_wp*ah_max, DT, a3, m2, m3)
      call check(error, abs(a3(n/2, n/2) - 1.0e-3_wp*ah_max) <= 1.0e-13_wp*ah_max, &
                 "Ah below its joint bound is not modified")
   end subroutine test_joint_bound

   subroutine test_knob_off(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), allocatable :: du_a(:, :, :), du_b(:, :, :), du_c(:, :, :), fc(:, :)
      integer :: pass, i, j, nxt, nyt

      do pass = 1, 3
         call make_state(grid, ms, 24, 24, 5000.0_wp, 5000.0_wp, 80.0_wp)
         nxt = grid%nx_total
         nyt = grid%ny_total
         do j = 1, nyt
            do i = 1, nxt + 1
               ms%u_face_x_layer(i, j, 1) = 0.3_wp*sin(0.9_wp*real(i, wp) + 1.3_wp*real(j, wp))
            end do
         end do
         do j = 1, nyt + 1
            do i = 1, nxt
               ms%v_face_y_layer(i, j, 1) = 0.3_wp*cos(1.1_wp*real(i, wp) - 0.7_wp*real(j, wp))
            end do
         end do
         allocate (fc(nxt + 1, nyt + 1), source=0.0_wp)
         call hv%init(grid, nz_ml=1)
         hv%stress_tensor = pass /= 3
         hv%nu_h = merge(500.0_wp, 0.0_wp, pass /= 3)
         hv%nu_4 = merge(0.0_wp, 1.0e8_wp, pass == 1)
         call run_hv(grid, ms, hv, 10.0_wp, .false., fc, 0.0_wp, 0.0_wp, .false., .false., 6.0_wp)
         if (.not. allocated(error)) then
            call check(error,.not. allocated(hv%bh_bg_t) .and. .not. hv%biharm_stress, &
                       "knob off: no biharmonic statics")
         end if
         select case (pass)
         case (1)
            allocate (du_a, source=hv%du_visc%data)
         case (2)
            allocate (du_b, source=hv%du_visc%data)
         case (3)
            allocate (du_c, source=hv%du_visc%data)
         end select
         call hv%destroy()
         call ms%destroy()
         deallocate (fc)
         if (allocated(error)) return
      end do
      ! pass 1: Laplacian only (stress, nu_h); pass 2: + velocity biharmonic;
      ! pass 3: velocity biharmonic only (nu_h = 0 scalar Laplacian = 0).
      ! Not `==`: the velocity-form subtraction may contract to an FMA when it
      ! starts from a nonzero Laplacian tendency.  Exact knob-off bit-identity
      ! against the previous build is the job of the gates.sh cases.
      call check(error, maxval(abs(du_b - (du_a + du_c))) <= 1.0e-13_wp*maxval(abs(du_c)) .and. &
                 maxval(abs(du_c)) > 0.0_wp, &
                 "knob off: stress + nu_4 equals Laplacian-only + velocity-form biharmonic")
   end subroutine test_knob_off

end module test_ocean_hvisc_biharm_stress
