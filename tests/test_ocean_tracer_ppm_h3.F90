!! Unit tests for the `ppm_h3` horizontal tracer-advection scheme
!! (`&ocean_vmix_nml tracer_advection_scheme = "ppm_h3"`): Huynh (1997)
!! 3rd-order edges, Colella & Woodward (1984, J. Comput. Phys. 54) limiter,
!! swept-volume flux.  All tests drive `tracer_advect_zonal_one_impl`, the
!! production kernel, on explicit arrays (device-mapped on GPU builds):
!!
!!   * linear profile: the flux equals the EXACT swept mean of the line,
!!     both flow directions (wrong edge weights or a sign slip fail this);
!!   * local extremum: the flux is the plain donor value (PCM fallback);
!!   * smooth bump, CFL 0.5, one revolution: L1 error is reported for `ppm`
!!     and `ppm_h3`; `ppm_h3` must beat `ppm` by at least 2x, and the
!!     reversed flow must reproduce the mirror-image error;
!!   * step function: no new extrema;
!!   * conservation of sum(h*T) over the periodic revolution to 1e-13.
module test_ocean_tracer_ppm_h3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: output_unit
   use rdb_constants, only: wp
   use rdb_continuity, only: tracer_advect_zonal_one_impl
   implicit none
   private

   public :: collect_ocean_tracer_ppm_h3_tests

   integer, parameter :: NP = 64          !! physical cells
   integer, parameter :: NG = 3           !! ghost columns each side
   integer, parameter :: NX = NP + 2*NG
   integer, parameter :: NY = 3
   integer, parameter :: NZ = 1
   real(wp), parameter :: DT = 2.0_wp   !! step length; a unit step would hide a dt slip in the flux

contains

   subroutine collect_ocean_tracer_ppm_h3_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("linear_profile_exact", test_linear), &
                  new_unittest("extremum_is_pcm", test_extremum), &
                  new_unittest("bump_revolution_l1", test_bump), &
                  new_unittest("step_no_new_extrema", test_step), &
                  new_unittest("revolution_conserves", test_conserve) &
                  ]
   end subroutine collect_ocean_tracer_ppm_h3_tests

   subroutine advect_steps(nsteps, mass, T, h3, tot0, tot1)
      !! Advance the periodic row T(NP) `nsteps` times with the production
      !! zonal kernel.  `mass` is the per-step volume swept in cell units (= the
      !! Courant number), so the face mass flux handed to the kernel is `mass/DT`.
      integer, intent(in) :: nsteps
      real(wp), intent(in) :: mass
      real(wp), intent(inout) :: T(NP)
      logical, intent(in) :: h3
      real(wp), intent(out), optional :: tot0, tot1
      real(wp) :: iareaT(NX, NY), wet(NX, NY), h(NX, NY, NZ), hTr(NX, NY, NZ)
      real(wp) :: mfx(NX + 1, NY, NZ), fxl(NX + 1, NY, NZ), fxr(NX + 1, NY, NZ)
      integer :: n, j

      iareaT = 1.0_wp
      wet = 1.0_wp
      h = 1.0_wp
      mfx = mass/DT
      fxl = 0.0_wp
      fxr = 0.0_wp
      do j = 1, NY
         hTr(NG + 1:NG + NP, j, 1) = T
      end do
      if (present(tot0)) tot0 = sum(hTr(NG + 1:NG + NP, :, :))
      !$acc enter data copyin(iareaT, wet, h, hTr, mfx) create(fxl, fxr)
      do n = 1, nsteps
         call wrap_ghosts(hTr)
         !$acc update device(hTr)
         call tracer_advect_zonal_one_impl(NX, NY, NZ, DT, iareaT, wet, h, hTr, &
                                           mfx, fxl, fxr, h3=h3)
         !$acc update self(hTr)
      end do
      !$acc exit data delete(iareaT, wet, h, hTr, mfx, fxl, fxr)
      T = hTr(NG + 1:NG + NP, 1, 1)
      if (present(tot1)) tot1 = sum(hTr(NG + 1:NG + NP, :, :))
   end subroutine advect_steps

   subroutine wrap_ghosts(f)
      real(wp), intent(inout) :: f(NX, NY, NZ)
      f(1:NG, :, :) = f(NP + 1:NP + NG, :, :)
      f(NG + NP + 1:NX, :, :) = f(NG + 1:2*NG, :, :)
   end subroutine wrap_ghosts

   subroutine bump(T)
      real(wp), intent(out) :: T(NP)
      integer :: i
      do i = 1, NP
         T(i) = exp(-((real(i, wp) - 0.5_wp*NP)/5.0_wp)**2)
      end do
   end subroutine bump

   subroutine test_linear(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: iareaT(NX, NY), wet(NX, NY), h(NX, NY, NZ), hTr(NX, NY, NZ)
      real(wp) :: mfx(NX + 1, NY, NZ), fxl(NX + 1, NY, NZ), fxr(NX + 1, NY, NZ)
      real(wp), parameter :: A = 3.0_wp, B = 0.7_wp, CFL = 0.4_wp
      real(wp) :: err, xm, exact
      integer :: i, sgn

      err = 0.0_wp
      do sgn = 1, -1, -2
         iareaT = 1.0_wp
         wet = 1.0_wp
         h = 1.0_wp
         mfx = real(sgn, wp)*CFL/DT
         fxl = 0.0_wp
         fxr = 0.0_wp
         do i = 1, NX
            hTr(i, :, 1) = A + B*real(i, wp)
         end do
         !$acc enter data copyin(iareaT, wet, h, hTr, mfx) create(fxl, fxr)
         call tracer_advect_zonal_one_impl(NX, NY, NZ, DT, iareaT, wet, h, hTr, &
                                           mfx, fxl, fxr, h3=.true.)
         !$acc update self(fxl)
         !$acc exit data delete(iareaT, wet, h, hTr, mfx, fxl, fxr)
         do i = 4, NX - 2
            ! face i sits at x = i - 1/2; the flow sweeps CFL cells upstream.
            xm = real(i, wp) - 0.5_wp - real(sgn, wp)*0.5_wp*CFL
            exact = real(sgn, wp)*CFL/DT*(A + B*xm)
            err = max(err, abs(fxl(i, 2, 1) - exact))
         end do
      end do
      call check(error, err < 1.0e-13_wp, "linear flux != exact swept mean, err="//fmt(err))
   end subroutine test_linear

   subroutine test_extremum(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: iareaT(NX, NY), wet(NX, NY), h(NX, NY, NZ), hTr(NX, NY, NZ)
      real(wp) :: mfx(NX + 1, NY, NZ), fxl(NX + 1, NY, NZ), fxr(NX + 1, NY, NZ)
      real(wp), parameter :: CFL = 0.45_wp
      integer, parameter :: IP = 20     !! peak cell
      real(wp) :: e_pos, e_neg, e_other

      iareaT = 1.0_wp
      wet = 1.0_wp
      h = 1.0_wp
      fxl = 0.0_wp
      fxr = 0.0_wp
      hTr = 0.0_wp
      hTr(IP, :, 1) = 1.0_wp
      mfx = CFL/DT
      !$acc enter data copyin(iareaT, wet, h, hTr, mfx) create(fxl, fxr)
      call tracer_advect_zonal_one_impl(NX, NY, NZ, DT, iareaT, wet, h, hTr, &
                                        mfx, fxl, fxr, h3=.true.)
      !$acc update self(fxl)
      !$acc exit data delete(iareaT, wet, h, hTr, mfx, fxl, fxr)
      ! Right face of the peak cell (face IP+1), flow to the right: PCM => U*Tc.
      e_pos = abs(fxl(IP + 1, 2, 1) - CFL/DT)
      ! Flow to the left through the peak's left face (face IP): redo reversed.
      hTr = 0.0_wp
      hTr(IP, :, 1) = 1.0_wp
      mfx = -CFL/DT
      fxl = 0.0_wp
      fxr = 0.0_wp
      !$acc enter data copyin(iareaT, wet, h, hTr, mfx) create(fxl, fxr)
      call tracer_advect_zonal_one_impl(NX, NY, NZ, DT, iareaT, wet, h, hTr, &
                                        mfx, fxl, fxr, h3=.true.)
      !$acc update self(fxl)
      !$acc exit data delete(iareaT, wet, h, hTr, mfx, fxl, fxr)
      e_neg = abs(fxl(IP, 2, 1) + CFL/DT)
      ! A flat zero region far from the peak must carry exactly nothing.
      e_other = abs(fxl(40, 2, 1))
      call check(error, max(e_pos, e_neg, e_other) < 1.0e-15_wp, &
                 "extremum not PCM: "//fmt(e_pos)//" "//fmt(e_neg)//" "//fmt(e_other))
   end subroutine test_extremum

   subroutine test_bump(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: T0(NP), Tp(NP), Th(NP), Tn(NP)
      real(wp) :: l1_ppm, l1_h3, l1_neg
      integer, parameter :: NSTEP = 128   !! CFL 0.5 x 128 = one revolution of 64 cells

      call bump(T0)
      Tp = T0
      Th = T0
      Tn = T0
      call advect_steps(NSTEP, 0.5_wp, Tp, .false.)
      call advect_steps(NSTEP, 0.5_wp, Th, .true.)
      call advect_steps(NSTEP, -0.5_wp, Tn, .true.)
      l1_ppm = sum(abs(Tp - T0))/real(NP, wp)
      l1_h3 = sum(abs(Th - T0))/real(NP, wp)
      l1_neg = sum(abs(Tn - T0))/real(NP, wp)
      write (output_unit, '(a,es11.4,a,es11.4,a,es11.4)') &
         "ppm_h3 bump L1 (CFL 0.5, 1 rev): ppm=", l1_ppm, " ppm_h3=", l1_h3, &
         " ppm_h3 reversed=", l1_neg
      call check(error, l1_h3 < 0.5_wp*l1_ppm, &
                 "ppm_h3 L1 "//fmt(l1_h3)//" not < half of ppm L1 "//fmt(l1_ppm))
      if (allocated(error)) return
      call check(error, abs(l1_neg - l1_h3) < 1.0e-10_wp, &
                 "reversed flow is not the mirror image: "//fmt(l1_neg)//" vs "//fmt(l1_h3))
   end subroutine test_bump

   subroutine test_step(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: T(NP)
      integer :: i

      T = 0.0_wp
      do i = 20, 40
         T(i) = 1.0_wp
      end do
      call advect_steps(128, 0.5_wp, T, .true.)
      call check(error, minval(T) >= -1.0e-14_wp .and. maxval(T) <= 1.0_wp + 1.0e-14_wp, &
                 "new extrema: min="//fmt(minval(T))//" max="//fmt(maxval(T)))
   end subroutine test_step

   subroutine test_conserve(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: T(NP), tot0, tot1

      call bump(T)
      call advect_steps(128, 0.5_wp, T, .true., tot0, tot1)
      call check(error, abs(tot1 - tot0) <= 1.0e-13_wp*abs(tot0), &
                 "sum(h*T) drift "//fmt(tot1 - tot0))
   end subroutine test_conserve

   function fmt(x) result(s)
      real(wp), intent(in) :: x
      character(len=24) :: s
      write (s, '(es12.4)') x
   end function fmt

end module test_ocean_tracer_ppm_h3
