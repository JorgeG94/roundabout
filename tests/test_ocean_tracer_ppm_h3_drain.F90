!! Unit tests for `ppm_h3` inside the windowed tracer-advection drain
!! (`&ocean_vmix_nml tracer_advection_scheme = "ppm_h3"`,
!! `dt_tracer_advect_ratio > 1`).  The drain moves per-pass VOLUMES
!! (`uhh`), not rates, so the tests drive `drain_swept_flux_x_h3` with
!! volumes and update `hTr`/`hprev` exactly as `continuity_tracer_drain`
!! does (`X -= (F(i+1) - F(i))*iareaT`):
!!
!!   * a uniform tracer stays uniform under a divergent (hprev-changing)
!!     windowed flow, to round-off;
!!   * a 1-D sine translated one revolution converges at better than 1.7
!!     order and beats the legacy drain reconstruction (`drain_parabola_x`
!!     + `drain_swept_flux_x`) on L1 (both reported);
!!   * a step profile under a divergent flow creates no new extrema.
module test_ocean_tracer_ppm_h3_drain
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: output_unit
   use rdb_constants, only: wp
   use rdb_continuity, only: drain_swept_flux_x_h3, drain_swept_flux_x, drain_parabola_x, &
                             drain_fill_conc
   implicit none
   private

   public :: collect_ocean_tracer_ppm_h3_drain_tests

   integer, parameter :: NG = 3
   integer, parameter :: NY = 2
   real(wp), parameter :: PI = 3.14159265358979323846_wp

contains

   subroutine collect_ocean_tracer_ppm_h3_drain_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("drain_uniform_stays_uniform", test_uniform), &
                  new_unittest("drain_sine_translation_order", test_sine), &
                  new_unittest("drain_step_no_new_extrema", test_step) &
                  ]
   end subroutine collect_ocean_tracer_ppm_h3_drain_tests

   subroutine wrap_c(np, a)
      integer, intent(in) :: np
      real(wp), intent(inout) :: a(np + 2*NG, NY)
      a(1:NG, :) = a(np + 1:np + NG, :)
      a(NG + np + 1:np + 2*NG, :) = a(NG + 1:2*NG, :)
   end subroutine wrap_c

   subroutine one_pass(np, uhh, hprev, hTr, h3)
      !! One drain pass (zonal): flux from the chosen reconstruction, then
      !! the production update of hTr and hprev.  Periodic in x.
      integer, intent(in) :: np
      real(wp), intent(in) :: uhh(np + 2*NG + 1, NY, 1)
      real(wp), intent(inout) :: hprev(np + 2*NG, NY, 1), hTr(np + 2*NG, NY, 1)
      logical, intent(in) :: h3
      integer :: nx, i
      real(wp), allocatable :: areaT(:, :), iareaT(:, :), wet(:, :), tr(:, :, :), F(:, :, :)
      real(wp), allocatable :: pal(:, :, :), par(:, :, :), pa6(:, :, :)
      nx = np + 2*NG
      allocate (areaT(nx, NY), iareaT(nx, NY), wet(nx, NY), tr(nx, NY, 1), F(nx + 1, NY, 1))
      allocate (pal(nx, NY, 1), par(nx, NY, 1), pa6(nx, NY, 1))
      areaT = 1.0_wp
      iareaT = 1.0_wp
      wet = 1.0_wp
      tr = 0.0_wp
      F = 0.0_wp
      pal = 0.0_wp
      par = 0.0_wp
      pa6 = 0.0_wp
      call wrap_c(np, hprev(:, :, 1))
      call wrap_c(np, hTr(:, :, 1))
      if (h3) then
         !$acc data copyin(areaT, uhh, hprev, hTr, wet) create(tr) copyout(F)
         call drain_fill_conc(nx, NY, 1, hTr, hprev, tr)
         !$acc update self(tr)
         call wrap_c(np, tr(:, :, 1))
         !$acc update device(tr)
         call drain_swept_flux_x_h3(nx, NY, 1, areaT, uhh, hprev, tr, wet, F)
         !$acc end data
      else
         !$acc data copyin(areaT, uhh, hprev, hTr, wet) create(tr, pal, par, pa6) copyout(F)
         call drain_parabola_x(nx, NY, 1, wet, hTr, hprev, tr, pal, par, pa6)
         !$acc update self(pal, par, pa6)
         call wrap_c(np, pal(:, :, 1))
         call wrap_c(np, par(:, :, 1))
         call wrap_c(np, pa6(:, :, 1))
         !$acc update device(pal, par, pa6)
         call drain_swept_flux_x(nx, NY, 1, areaT, uhh, hprev, pal, par, pa6, F)
         !$acc end data
      end if
      do i = 1, nx
         hTr(i, :, 1) = hTr(i, :, 1) - (F(i + 1, :, 1) - F(i, :, 1))*iareaT(i, :)
         hprev(i, :, 1) = hprev(i, :, 1) - (uhh(i + 1, :, 1) - uhh(i, :, 1))*iareaT(i, :)
      end do
      call wrap_c(np, hprev(:, :, 1))
      call wrap_c(np, hTr(:, :, 1))
   end subroutine one_pass

   subroutine div_flow(np, amp, uhh)
      !! Periodic divergent face volumes: |uhh| < 0.45, hprev changes by
      !! O(amp) each pass.
      integer, intent(in) :: np
      real(wp), intent(in) :: amp
      real(wp), intent(out) :: uhh(np + 2*NG + 1, NY, 1)
      integer :: i
      do i = 1, np + 2*NG + 1
         uhh(i, :, 1) = 0.25_wp + amp*sin(2.0_wp*PI*real(modulo(i - NG - 1, np), wp)/real(np, wp))
      end do
   end subroutine div_flow

   subroutine test_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: np = 48
      real(wp) :: uhh(np + 2*NG + 1, NY, 1), hprev(np + 2*NG, NY, 1), hTr(np + 2*NG, NY, 1)
      real(wp), parameter :: C0 = 7.5_wp
      real(wp) :: err, exc
      integer :: n
      exc = 0.0_wp
      hprev = 1.0_wp
      hTr = C0
      call div_flow(np, 0.15_wp, uhh)
      do n = 1, 12
         ! alternate the sign so hprev oscillates instead of drifting
         if (mod(n, 2) == 0) uhh = -uhh
         call one_pass(np, uhh, hprev, hTr, .true.)
         exc = max(exc, maxval(abs(hprev - 1.0_wp)))
      end do
      err = maxval(abs(hTr(NG + 1:NG + np, :, 1)/hprev(NG + 1:NG + np, :, 1) - C0))/C0
      write (output_unit, '(a,es10.3,a,es10.3)') '  [drain h3] uniform rel err = ', err, &
         '  max hprev excursion = ', exc
      call check(error, exc > 0.01_wp, "flow must actually be divergent")
      if (allocated(error)) return
      call check(error, err < 1.0e-13_wp, "uniform tracer must stay uniform")
   end subroutine test_uniform

   function sine_l1(np, h3) result(l1)
      !! L1 error after one revolution of sin(2 pi x) at Courant 0.4.
      integer, intent(in) :: np
      logical, intent(in) :: h3
      real(wp) :: l1
      real(wp), parameter :: CU = 0.4_wp
      real(wp), allocatable :: uhh(:, :, :), hprev(:, :, :), hTr(:, :, :), ex(:)
      integer :: n, nsteps, i
      allocate (uhh(np + 2*NG + 1, NY, 1), hprev(np + 2*NG, NY, 1), hTr(np + 2*NG, NY, 1), ex(np))
      uhh = CU
      hprev = 1.0_wp
      do i = 1, np + 2*NG
         hTr(i, :, 1) = 1.0_wp + 0.5_wp*sin(2.0_wp*PI*(real(modulo(i - NG - 1, np), wp) + 0.5_wp)/real(np, wp))
      end do
      nsteps = nint(real(np, wp)/CU)
      do n = 1, nsteps
         call one_pass(np, uhh, hprev, hTr, h3)
      end do
      ! nsteps*CU = np exactly only when np/CU is an integer: np = 32,64,128
      ! and CU = 0.4 give 80, 160, 320.
      do i = 1, np
         ex(i) = 1.0_wp + 0.5_wp*sin(2.0_wp*PI*(real(i - 1, wp) + 0.5_wp)/real(np, wp))
      end do
      l1 = sum(abs(hTr(NG + 1:NG + np, 1, 1)/hprev(NG + 1:NG + np, 1, 1) - ex))/real(np, wp)
   end function sine_l1

   subroutine test_sine(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h32, h64, h128, l64, o1, o2
      h32 = sine_l1(32, .true.)
      h64 = sine_l1(64, .true.)
      h128 = sine_l1(128, .true.)
      l64 = sine_l1(64, .false.)
      o1 = log(h32/h64)/log(2.0_wp)
      o2 = log(h64/h128)/log(2.0_wp)
      write (output_unit, '(a,3es11.3,a,2f7.2)') '  [drain h3] sine L1 n=32/64/128 = ', h32, h64, h128, &
         '  order = ', o1, o2
      write (output_unit, '(a,es11.3,a,es11.3)') '  [drain] n=64 L1: legacy = ', l64, '  h3 = ', h64
      call check(error, min(o1, o2) > 1.7_wp, "h3 drain must converge at > 1.7 order")
      if (allocated(error)) return
      call check(error, h64 < l64, "h3 drain must beat the legacy drain scheme on L1")
   end subroutine test_sine

   subroutine test_step(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: np = 64
      real(wp) :: uhh(np + 2*NG + 1, NY, 1), hprev(np + 2*NG, NY, 1), hTr(np + 2*NG, NY, 1)
      real(wp) :: c, cmin, cmax
      integer :: i, n
      hprev = 1.0_wp
      do i = 1, np + 2*NG
         c = 0.0_wp
         if (modulo(i - NG - 1, np) >= np/4 .and. modulo(i - NG - 1, np) < np/2) c = 1.0_wp
         hTr(i, :, 1) = c
      end do
      call div_flow(np, 0.2_wp, uhh)
      cmin = 0.0_wp
      cmax = 1.0_wp
      do n = 1, 6
         if (mod(n, 2) == 0) uhh = -uhh
         call one_pass(np, uhh, hprev, hTr, .true.)
         cmin = min(cmin, minval(hTr(NG + 1:NG + np, :, 1)/hprev(NG + 1:NG + np, :, 1)))
         cmax = max(cmax, maxval(hTr(NG + 1:NG + np, :, 1)/hprev(NG + 1:NG + np, :, 1)))
      end do
      write (output_unit, '(a,2es11.3)') '  [drain h3] step min/max excursion = ', cmin, cmax - 1.0_wp
      call check(error, cmin > -1.0e-12_wp .and. cmax < 1.0_wp + 1.0e-12_wp, "no new extrema")
   end subroutine test_step

end module test_ocean_tracer_ppm_h3_drain
