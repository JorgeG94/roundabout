!> Golden-vector tests for the Redi continuous neutral-surface sweep (R1).
!!
!! Reproduces the MOM6 `ndiff_unit_tests_continuous` fixtures
!! (`MOM_neutral_diffusion.F90` lines ~2657-2854) to ~1e-12:
!!   * `redi_interface_scalar` PPM edge reconstruction (line ~2657);
!!   * `redi_interpolate_position` (lines ~2662-2682);
!!   * `redi_neutral_positions_continuous` — identical / slightly-cooler /
!!     no-overlap column pairs (lines ~2684-2779).
!!
!! INDEXING: R1 is transcribed TOP-DOWN (faithful to MOM6, k=1 = surface), so
!! the golden vectors are used AS-IS with no k-flip.  The bottom-up flip is
!! deferred to R2's flux-kernel caller (see rdb_ocean_redi module header).  The
!! transcribed `find_neutral_surface_positions_continuous` fixtures pass T/S as
!! INTERFACE arrays (nk+1) directly; PoL/PoR are fractional, KoL/KoR are layer
!! indices, hEff has 2*nk+1 entries.
module test_ocean_redi
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_redi, only: redi_interface_scalar, redi_interpolate_position, &
                             redi_neutral_positions_continuous, &
                             ocean_redi_t, redi_calc_coeffs, redi_apply_flux
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_redi_tests

   real(wp), parameter :: TOL = 1.0e-12_wp
   integer, parameter :: NGHOST = 1
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: T0 = 10.0_wp, S0 = 35.0_wp
   real(wp), parameter :: ALPHA_T = 0.2_wp     ! kg/m^3/degC
   real(wp), parameter :: BETA_S = 0.78_wp     ! kg/m^3/PSU
   real(wp), parameter :: DT = 1800.0_wp

contains

   subroutine collect_ocean_redi_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("redi_interface_scalar_ppm", test_interface_scalar), &
                  new_unittest("redi_interpolate_position", test_interpolate_position), &
                  new_unittest("redi_positions_identical", test_positions_identical), &
                  new_unittest("redi_positions_cooler", test_positions_cooler), &
                  new_unittest("redi_positions_no_overlap", test_positions_no_overlap), &
                  new_unittest("redi_flux_conserves", test_flux_conserves), &
                  new_unittest("redi_along_isopycnal_zero", test_along_isopycnal_zero), &
                  new_unittest("redi_cross_isopycnal_downgradient", test_cross_downgradient), &
                  new_unittest("redi_varmix_khtr_consumed", test_varmix_khtr), &
                  new_unittest("redi_wall_no_leak", test_wall_no_leak) &
                  ]
   end subroutine collect_ocean_redi_tests

   ! --------------------------------------------------------------------------
   ! No along-isopycnal flux crosses a no-normal-flow WALL: a cross-gradient
   ! tracer in a closed basin (no `bc` => all edges WALL) must conserve its
   ! PHYSICAL-domain content to round-off.  Before the wall closure the flux
   ! crossed the physical-boundary faces into the ghost halo — the full-array
   ! sum stayed closed but the physical domain drifted (the benchmark_ALE
   ! full-5 heat-leak signature).
   ! --------------------------------------------------------------------------
   subroutine test_wall_no_leak(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8, NSTEP = 40
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp, GZ = 0.02_wp
      real(wp), parameter :: GX = 5.0e-4_wp   ! x-tilt of T => tilted isopycnals
      real(wp) :: z_c, x_c, t_val, pt0, pt1, pchange
      real(wp), allocatable :: t_before(:, :, :)
      integer :: i, j, k, ni, nj, i0, i1, j0, j1, step
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_redi(rd, grid, NZ, khtr=100.0_wp)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Stratified T tilted in x => sloping isopycnals => STRONG
         ! along-isopycnal Redi flux (incl. across the physical-wall faces).
         ! S uniform (its content is the closed-wall invariant we check).
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
               end do
            end do
         end do
         ! PHYSICAL interior (exclude ghosts).
         i0 = grid%nghost + 1
         i1 = grid%nghost + grid%nx_phys
         j0 = grid%nghost + 1
         j1 = grid%nghost + grid%ny_phys
         pt0 = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))
         allocate (t_before, source=ms%tracers(ms%idx_temperature)%hTr)

         call map_in(ms, metrics, grid, rd, eos)
         ! Accumulate the leak over many steps (no bc => closed walls).
         do step = 1, NSTEP
            call run_redi(grid, metrics, eos, rd, ms)
         end do
         call map_out(ms, metrics, rd)

         pt1 = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))
         ! Non-triviality: the along-isopycnal flux genuinely redistributed T
         ! within the physical domain (else conservation is vacuous).
         pchange = maxval(abs(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :) &
                              - t_before(i0:i1, j0:j1, :)))
         deallocate (t_before)
         call check(error, abs(pt1 - pt0) < 1.0e-11_wp*abs(pt0), &
                    "Redi must conserve PHYSICAL-domain Sum(h*T) over "// &
                    "many steps (no flux through closed walls)")
         if (allocated(error)) exit checks
         call check(error, pchange > 1.0e-6_wp, "Redi flux must be active (non-trivial)")
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_wall_no_leak

   ! --------------------------------------------------------------------------
   ! MOM6 line ~2657:
   !   call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 2, h_neglect)
   !   test_data1d(v,5, Tio, (/24.,22.,15.,8.,6./), 'Linear profile, PPM interface temperatures')
   ! --------------------------------------------------------------------------
   subroutine test_interface_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NK = 4
      real(wp) :: h(NK), tr(NK), edge(NK + 1)
      real(wp) :: expect(NK + 1)

      h = [10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp]
      tr = [24.0_wp, 18.0_wp, 12.0_wp, 6.0_wp]
      expect = [24.0_wp, 22.0_wp, 15.0_wp, 8.0_wp, 6.0_wp]

      call redi_interface_scalar(NK, h, tr, edge)
      call check(error, maxval(abs(edge - expect)) < TOL, &
                 "PPM interface scalar reconstruction mismatch")
   end subroutine test_interface_scalar

   ! --------------------------------------------------------------------------
   ! MOM6 lines ~2662-2682 (test_ifndp(v, dRhoNeg, Pneg, dRhoPos, Ppos, Ptrue)):
   !   (-1, 0, 1, 1) -> 0.5   'mid-point'
   !   ( 0, 0, 1, 1) -> 0.0   'bottom'
   !   ( 0.1,0,1.1,1)-> 0.0   'below'
   !   (-1, 0, 0, 1) -> 1.0   'top'
   !   (-1, 0,-0.1,1)-> 1.0   'above'
   !   (-1, 0, 3, 1) -> 0.25  '1/4'
   !   (-3, 0, 1, 1) -> 0.75  '3/4'
   !   ( 1, 0, 1, 1) -> 0.0   'dRho=0 below'
   !   (-1, 0,-1, 1) -> 1.0   'dRho=0 above'
   !   ( 0, 0, 0, 1) -> 0.5   'dRho=0 mid'
   !   (-2,.5, 5,.5) -> 0.5   'dP=0'
   ! --------------------------------------------------------------------------
   subroutine test_interpolate_position(error)
      type(error_type), allocatable, intent(out) :: error

      call chk(error, redi_interpolate_position(-1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp), 0.5_wp, "mid-point")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp), 0.0_wp, "bottom")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(0.1_wp, 0.0_wp, 1.1_wp, 1.0_wp), 0.0_wp, "below")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp), 1.0_wp, "top")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-1.0_wp, 0.0_wp, -0.1_wp, 1.0_wp), 1.0_wp, "above")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-1.0_wp, 0.0_wp, 3.0_wp, 1.0_wp), 0.25_wp, "1/4")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-3.0_wp, 0.0_wp, 1.0_wp, 1.0_wp), 0.75_wp, "3/4")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp), 0.0_wp, "dRho=0 below")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-1.0_wp, 0.0_wp, -1.0_wp, 1.0_wp), 1.0_wp, "dRho=0 above")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp), 0.5_wp, "dRho=0 mid")
      if (allocated(error)) return
      call chk(error, redi_interpolate_position(-2.0_wp, 0.5_wp, 5.0_wp, 0.5_wp), 0.5_wp, "dP=0")
   end subroutine test_interpolate_position

   ! --------------------------------------------------------------------------
   ! MOM6 lines ~2684-2697 'Identical columns':
   !   nk=3, Pl=Pr=(0,10,20,30), Tl=Tr=(22,18,14,10), Sl=Sr=0,
   !   dRdT=-1, dRdS=1 (both columns)
   !   KoL=(1,1,2,2,3,3,3,3) KoR=(1,1,2,2,3,3,3,3)
   !   pL=pR=(0,0,0,0,0,0,1,1)  hEff=(0,10,0,10,0,10,0)
   ! --------------------------------------------------------------------------
   subroutine test_positions_identical(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NK = 3
      real(wp) :: P(NK + 1), T(NK + 1), S(NK + 1), dRdT(NK + 1), dRdS(NK + 1)
      real(wp) :: PoL(2*NK + 2), PoR(2*NK + 2), hEff(2*NK + 1)
      integer :: KoL(2*NK + 2), KoR(2*NK + 2)

      P = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      T = [22.0_wp, 18.0_wp, 14.0_wp, 10.0_wp]
      S = 0.0_wp
      dRdT = -1.0_wp
      dRdS = 1.0_wp

      call redi_neutral_positions_continuous(NK, P, T, S, dRdT, dRdS, &
                                             P, T, S, dRdT, dRdS, &
                                             PoL, PoR, KoL, KoR, hEff)

      call check(error, all(KoL == [1, 1, 2, 2, 3, 3, 3, 3]), "identical KoL")
      if (allocated(error)) return
      call check(error, all(KoR == [1, 1, 2, 2, 3, 3, 3, 3]), "identical KoR")
      if (allocated(error)) return
      call check(error, maxval(abs(PoL - [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp])) < TOL, &
                 "identical PoL")
      if (allocated(error)) return
      call check(error, maxval(abs(PoR - [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp])) < TOL, &
                 "identical PoR")
      if (allocated(error)) return
      call check(error, maxval(abs(hEff - [0.0_wp, 10.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, 10.0_wp, 0.0_wp])) < TOL, &
                 "identical hEff")
   end subroutine test_positions_identical

   ! --------------------------------------------------------------------------
   ! MOM6 lines ~2715-2728 'Right column slightly cooler':
   !   nk=3, P=(0,10,20,30); Tl=(22,18,14,10), Tr=(20,16,12,8); S=0;
   !   dRdT=-1, dRdS=1
   !   KoL=(1,1,2,2,3,3,3,3) KoR=(1,1,1,2,2,3,3,3)
   !   pL=(0,.5,0,.5,0,.5,1,1)  pR=(0,0,.5,0,.5,0,.5,1)  hEff=(0,5,5,5,5,5,0)
   ! --------------------------------------------------------------------------
   subroutine test_positions_cooler(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NK = 3
      real(wp) :: P(NK + 1), Tl(NK + 1), Tr(NK + 1), S(NK + 1), dRdT(NK + 1), dRdS(NK + 1)
      real(wp) :: PoL(2*NK + 2), PoR(2*NK + 2), hEff(2*NK + 1)
      integer :: KoL(2*NK + 2), KoR(2*NK + 2)

      P = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      Tl = [22.0_wp, 18.0_wp, 14.0_wp, 10.0_wp]
      Tr = [20.0_wp, 16.0_wp, 12.0_wp, 8.0_wp]
      S = 0.0_wp
      dRdT = -1.0_wp
      dRdS = 1.0_wp

      call redi_neutral_positions_continuous(NK, P, Tl, S, dRdT, dRdS, &
                                             P, Tr, S, dRdT, dRdS, &
                                             PoL, PoR, KoL, KoR, hEff)

      call check(error, all(KoL == [1, 1, 2, 2, 3, 3, 3, 3]), "cooler KoL")
      if (allocated(error)) return
      call check(error, all(KoR == [1, 1, 1, 2, 2, 3, 3, 3]), "cooler KoR")
      if (allocated(error)) return
      call check(error, maxval(abs(PoL - [0.0_wp, 0.5_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.5_wp, 1.0_wp, 1.0_wp])) < TOL, &
                 "cooler PoL")
      if (allocated(error)) return
      call check(error, maxval(abs(PoR - [0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.5_wp, 1.0_wp])) < TOL, &
                 "cooler PoR")
      if (allocated(error)) return
      call check(error, maxval(abs(hEff - [0.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 0.0_wp])) < TOL, &
                 "cooler hEff")
   end subroutine test_positions_cooler

   ! --------------------------------------------------------------------------
   ! MOM6 lines ~2766-2779 'Right column much cooler than left with no overlap':
   !   nk=3, P=(0,10,20,30); Tl=(22,18,14,10), Tr=(9,7,5,3); S=0; dRdT=-1, dRdS=1
   !   KoL=(1,2,3,3,3,3,3,3) KoR=(1,1,1,1,1,2,3,3)
   !   pL=(0,0,0,1,1,1,1,1)  pR=(0,0,0,0,0,0,0,1)  hEff=(0,0,0,0,0,0,0)
   ! --------------------------------------------------------------------------
   subroutine test_positions_no_overlap(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NK = 3
      real(wp) :: P(NK + 1), Tl(NK + 1), Tr(NK + 1), S(NK + 1), dRdT(NK + 1), dRdS(NK + 1)
      real(wp) :: PoL(2*NK + 2), PoR(2*NK + 2), hEff(2*NK + 1)
      integer :: KoL(2*NK + 2), KoR(2*NK + 2)

      P = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      Tl = [22.0_wp, 18.0_wp, 14.0_wp, 10.0_wp]
      Tr = [9.0_wp, 7.0_wp, 5.0_wp, 3.0_wp]
      S = 0.0_wp
      dRdT = -1.0_wp
      dRdS = 1.0_wp

      call redi_neutral_positions_continuous(NK, P, Tl, S, dRdT, dRdS, &
                                             P, Tr, S, dRdT, dRdS, &
                                             PoL, PoR, KoL, KoR, hEff)

      call check(error, all(KoL == [1, 2, 3, 3, 3, 3, 3, 3]), "no-overlap KoL")
      if (allocated(error)) return
      call check(error, all(KoR == [1, 1, 1, 1, 1, 2, 3, 3]), "no-overlap KoR")
      if (allocated(error)) return
      call check(error, maxval(abs(PoL - [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp])) < TOL, &
                 "no-overlap PoL")
      if (allocated(error)) return
      call check(error, maxval(abs(PoR - [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp])) < TOL, &
                 "no-overlap PoR")
      if (allocated(error)) return
      call check(error, maxval(abs(hEff)) < TOL, "no-overlap hEff all zero")
   end subroutine test_positions_no_overlap

   !> Scalar check helper for redi_interpolate_position cases.
   subroutine chk(error, got, expect, label)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: got, expect
      character(*), intent(in) :: label

      call check(error, abs(got - expect) < TOL, &
                 "interpolate_position '"//label//"' mismatch")
   end subroutine chk

   ! ==================================================================
   ! R2 device tests — drive the production calc_coeffs + apply_flux on
   ! the NVHPC GPU build with a linear EOS (rho monotone in -T at fixed S).
   ! Bottom-up: k=1 bed, k=nz surface; exercises the k-flip end-to-end.
   ! ==================================================================

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%T_ref = T0
      eos%S_ref = S0
   end subroutine make_eos

   subroutine setup_redi(rd, grid, nz, khtr)
      type(ocean_redi_t), intent(inout) :: rd
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: khtr
      call rd%init(grid, nz_ml=nz)
      rd%enable = .true.
      rd%continuous = .true.
      rd%khtr = khtr
   end subroutine setup_redi

   subroutine map_in(ms, metrics, grid, rd, eos)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(ocean_redi_t), intent(inout) :: rd
      type(eos_t), intent(in) :: eos
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(rd)
      call rd%enter_data()
   end subroutine map_in

   subroutine run_redi(grid, metrics, eos, rd, ms)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(eos_t), intent(in) :: eos
      type(ocean_redi_t), intent(inout) :: rd
      type(multilayer_state_t), intent(inout) :: ms
      call redi_calc_coeffs(grid, metrics, eos, rd, ms)
      call redi_apply_flux(grid, metrics, rd, ms, DT)
   end subroutine run_redi

   subroutine map_out(ms, metrics, rd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_redi_t), intent(inout) :: rd
      ! Pull the updated tracers back to host.
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      call rd%exit_data()
      !$acc exit data delete(rd)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! Test: closed-domain conservation.  Sum(hTr) over the interior is
   ! unchanged after a Redi step (flux-form triad scatter telescopes).
   ! ------------------------------------------------------------------
   subroutine test_flux_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: GZ = 0.02_wp, GX = 5.0e-5_wp   ! gentle tilt (CFL-stable)
      ! Salinity carries a CROSS-isopycnal x-gradient (small in density terms vs
      ! the T stratification), so the neutral flux is genuinely NON-ZERO and the
      ! telescoping conservation is actually exercised — not trivially satisfied
      ! by an on-isopycnal tracer (T alone defines the isopycnals here).
      real(wp), parameter :: GSX = 1.0e-5_wp
      real(wp) :: z_c, x_c, t_val, s_val
      real(wp) :: tsum0, tsum1, ssum0, ssum1, s_change
      real(wp), allocatable :: s0f(:, :, :)
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_redi(rd, grid, NZ, khtr=100.0_wp)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c   ! tilted isopycnals (T-dominated)
                  s_val = S0 + GSX*x_c           ! cross-isopycnal S perturbation
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = s_val*DZ
               end do
            end do
         end do
         ! FULL wall-closed domain (outer faces i=1/ni+1, j=1/nj+1 are walls =>
         ! zero flux), so the full-array sum is the closed-system invariant.
         tsum0 = sum(ms%tracers(ms%idx_temperature)%hTr(1:ni, 1:nj, :))
         ssum0 = sum(ms%tracers(ms%idx_salinity)%hTr(1:ni, 1:nj, :))
         s0f = ms%tracers(ms%idx_salinity)%hTr   ! full snapshot for the change check

         call map_in(ms, metrics, grid, rd, eos)
         call run_redi(grid, metrics, eos, rd, ms)
         call map_out(ms, metrics, rd)

         tsum1 = sum(ms%tracers(ms%idx_temperature)%hTr(1:ni, 1:nj, :))
         ssum1 = sum(ms%tracers(ms%idx_salinity)%hTr(1:ni, 1:nj, :))
         ! Max change over the FULL domain.  A LINEAR S gradient has zero
         ! diffusive divergence in the interior (nabla^2 of a linear field = 0);
         ! the neutral flux only changes the boundary-adjacent cells — so the
         ! non-triviality guard must scan the whole array, not an interior cell.
         s_change = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - s0f))
         call check(error, abs(tsum1 - tsum0) < 1.0e-9_wp*abs(tsum0), &
                    "Redi must conserve full-domain Sum(h*T) to round-off")
         if (allocated(error)) exit checks
         call check(error, abs(ssum1 - ssum0) < 1.0e-9_wp*abs(ssum0), &
                    "Redi must conserve full-domain Sum(h*S) under NON-ZERO neutral flux (telescoping)")
         if (allocated(error)) exit checks
         ! Guard the conservation test is not trivial: the cross-isopycnal flux
         ! must actually move S (else Sum(h*S) is conserved for the wrong reason).
         call check(error, s_change > 1.0e-10_wp, &
                    "cross-isopycnal S flux must be non-zero (conservation must be exercised)")
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_flux_conserves

   ! ------------------------------------------------------------------
   ! Test: a tracer that is a FUNCTION OF DENSITY (constant on isopycnals)
   ! => near-zero interior neutral flux (the defining property; the
   ! interior column T is unchanged to a tight tolerance).
   ! ------------------------------------------------------------------
   subroutine test_along_isopycnal_zero(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      ! Gentle horizontal tilt: GX*DX = 0.05 degC is a small fraction of the
      ! per-layer vertical contrast GZ*DZ = 1.0 degC, so the neutral surfaces
      ! stay close to horizontal and the densest/lightest isopycnals only
      ! marginally outcrop at the bed/surface (a finite-depth tilted tracer
      ! ALWAYS outcrops there — the up-slope column's densest water has no
      ! neutral counterpart down-slope; MOM6 produces the same boundary term).
      ! Paired with a CFL-stable khtr (khtr*dt/dx^2 = 100*1800/1e6 = 0.18 < 0.5)
      ! the residual interior change stays at round-off.
      real(wp), parameter :: GZ = 0.02_wp, GX = 5.0e-5_wp
      real(wp), allocatable :: t_before(:, :, :)
      real(wp) :: z_c, x_c, t_val, mx
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call setup_redi(rd, grid, NZ, khtr=100.0_wp)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Linear EOS at fixed S => rho monotone in -T => "constant on
         ! isopycnals" == "constant T".  Set T = f(rho) by making T itself
         ! the tracer with a sloped field; along a neutral surface (same
         ! rho => same T) the matched tracer difference is ~0.
         allocate (t_before(ni, nj, NZ))
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            do j = 1, nj
               do i = 1, ni
                  x_c = real(i, wp)*DX
                  t_val = T0 + GZ*z_c + GX*x_c
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*DZ
                  t_before(i, j, k) = t_val
               end do
            end do
         end do

         call map_in(ms, metrics, grid, rd, eos)
         call run_redi(grid, metrics, eos, rd, ms)
         call map_out(ms, metrics, rd)

         ! Interior cells: T (=hTr/h) must be ~unchanged — the only tracer
         ! is T itself, which is constant along the neutral surfaces it
         ! defines, so the along-neutral difference (and flux) vanishes.
         mx = 0.0_wp
         do k = 1, NZ
            do j = 2, nj - 1
               do i = 2, ni - 1
                  mx = max(mx, abs(ms%tracers(ms%idx_temperature)%hTr(i, j, k)/DZ &
                                   - t_before(i, j, k)))
               end do
            end do
         end do
         call check(error, mx < 1.0e-6_wp, &
                    "tracer constant on isopycnals => ~zero interior neutral flux")
         deallocate (t_before)
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_along_isopycnal_zero

   ! ------------------------------------------------------------------
   ! Test: a cross-isopycnal tracer gradient (a PASSIVE tracer flat in z
   ! but varying in x, NOT aligned with density) => down-gradient flux
   ! only (sign guard holds: the high-x side loses, the low-x side gains;
   ! no new extrema).  Here T/S set the (flat) isopycnals and we inspect
   ! the salinity tracer carrying the cross-gradient.
   ! ------------------------------------------------------------------
   subroutine test_cross_downgradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp
      real(wp), parameter :: GZ = 0.02_wp
      real(wp), allocatable :: s_before(:, :, :)
      real(wp) :: z_c, t_val, s_val, mn, mx, lo_gain, hi_loss
      integer :: i, j, k, ni, nj
      checks: block
         call make_grid(grid, NX, NY, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         ! CFL-stable diffusion number: khtr*dt/dx^2 = 100*1800/1e6 = 0.18 < 0.5.
         call setup_redi(rd, grid, NZ, khtr=100.0_wp)
         call make_eos(eos)
         ni = grid%nx_total
         nj = grid%ny_total
         ms%h_layer = DZ
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! FLAT isopycnals: density set by T(z) only (no x-tilt) at fixed S
         ! contribution from temperature; salinity carries a horizontal
         ! (cross-isopycnal) step so the neutral surfaces are horizontal and
         ! the S-difference along them is purely down-gradient.
         !
         ! The salinity step must be SMALL: salinity feeds density too
         ! (dR/dS = beta_S = 0.78 kg/m^3/PSU), so a large dS would tilt the
         ! neutral surfaces and connect different depths (the surfaces would no
         ! longer be horizontal — defeating the "cross-isopycnal" framing and
         ! letting the along-neutral redistribution create apparent new
         ! extrema).  dS = 0.01 PSU/cell => 0.0078 kg/m^3, ~25x below the
         ! per-layer T-stratification (alpha_T*GZ*DZ = 0.2 kg/m^3), so the
         ! isopycnals stay ~horizontal and the S transport is purely
         ! down-gradient.
         allocate (s_before(ni, nj, NZ))
         do k = 1, NZ
            z_c = (real(k, wp) - 0.5_wp)*DZ
            t_val = T0 + GZ*z_c
            do j = 1, nj
               do i = 1, ni
                  s_val = S0 + 0.01_wp*real(i, wp)   ! gentle cross-isopycnal step
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = s_val*DZ
                  s_before(i, j, k) = s_val
               end do
            end do
         end do

         call map_in(ms, metrics, grid, rd, eos)
         call run_redi(grid, metrics, eos, rd, ms)
         call map_out(ms, metrics, rd)

         ! Down-gradient: no new extrema (S stays within the original range)
         ! and the low-x interior column gains while the high-x loses.
         mn = minval(s_before)
         mx = maxval(s_before)
         lo_gain = ms%tracers(ms%idx_salinity)%hTr(2, 3, NZ/2)/DZ - s_before(2, 3, NZ/2)
         hi_loss = ms%tracers(ms%idx_salinity)%hTr(ni - 1, 3, NZ/2)/DZ - s_before(ni - 1, 3, NZ/2)
         block
            real(wp) :: smin_after, smax_after
            integer :: ii, jj, kk
            smin_after = huge(1.0_wp)
            smax_after = -huge(1.0_wp)
            do kk = 1, NZ
               do jj = 2, nj - 1
                  do ii = 2, ni - 1
                     smin_after = min(smin_after, ms%tracers(ms%idx_salinity)%hTr(ii, jj, kk)/DZ)
                     smax_after = max(smax_after, ms%tracers(ms%idx_salinity)%hTr(ii, jj, kk)/DZ)
                  end do
               end do
            end do
            call check(error, smin_after >= mn - 1.0e-9_wp .and. smax_after <= mx + 1.0e-9_wp, &
                       "cross-isopycnal: down-gradient flux creates no new extrema")
         end block
         if (allocated(error)) exit checks
         call check(error, lo_gain >= -1.0e-12_wp .and. hi_loss <= 1.0e-12_wp, &
                    "cross-isopycnal: low-x column gains, high-x loses (down-gradient)")
         deallocate (s_before)
      end block checks
      call rd%destroy()
      call ms%destroy()
   end subroutine test_cross_downgradient

   ! --------------------------------------------------------------------------
   ! VarMix seam: a UNIFORM VarMix KhTr field (khtr_u/v = K) fed through the
   ! `khtr_u_ext`/`khtr_v_ext` path must reproduce the SCALAR-khtr=K result
   ! exactly — AND it must drive the flux even when the scalar `khtr` is 0
   ! (the VarMix field overrides the scalar no-flux gate).  Run-twice-compare:
   !   A: scalar khtr=K, no ext.   B: scalar khtr=0, uniform ext=K.
   ! S_after(A) == S_after(B) to round-off.
   ! --------------------------------------------------------------------------
   subroutine cross_iso_salinity(khtr_scalar, use_ext, khtr_ext_val, sout)
      real(wp), intent(in) :: khtr_scalar, khtr_ext_val
      logical, intent(in) :: use_ext
      real(wp), allocatable, intent(out) :: sout(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_redi_t) :: rd
      type(eos_t) :: eos
      integer, parameter :: NX = 6, NY = 5, NZ = 8
      real(wp), parameter :: DX = 1000.0_wp, DZ = 50.0_wp, GZ = 0.02_wp
      real(wp), allocatable :: khu(:, :), khv(:, :)
      real(wp) :: z_c, t_val, s_val
      integer :: i, j, k, ni, nj
      call make_grid(grid, NX, NY, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call setup_redi(rd, grid, NZ, khtr=khtr_scalar)
      call make_eos(eos)
      ni = grid%nx_total
      nj = grid%ny_total
      ms%h_layer = DZ
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         z_c = (real(k, wp) - 0.5_wp)*DZ
         t_val = T0 + GZ*z_c
         do j = 1, nj
            do i = 1, ni
               s_val = S0 + 0.01_wp*real(i, wp)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_val*DZ
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = s_val*DZ
            end do
         end do
      end do
      call map_in(ms, metrics, grid, rd, eos)
      if (use_ext) then
         allocate (khu(ni + 1, nj), source=khtr_ext_val)
         allocate (khv(ni, nj + 1), source=khtr_ext_val)
         !$acc enter data copyin(khu, khv)
         call redi_calc_coeffs(grid, metrics, eos, rd, ms)
         call redi_apply_flux(grid, metrics, rd, ms, DT, khtr_u_ext=khu, khtr_v_ext=khv)
         !$acc exit data delete(khu, khv)
         deallocate (khu, khv)
      else
         call run_redi(grid, metrics, eos, rd, ms)
      end if
      call map_out(ms, metrics, rd)
      allocate (sout(ni, nj, NZ))
      sout = ms%tracers(ms%idx_salinity)%hTr/DZ
      call rd%destroy()
      call ms%destroy()
   end subroutine cross_iso_salinity

   subroutine test_varmix_khtr(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: sa(:, :, :), sb(:, :, :), sc(:, :, :)
      real(wp) :: dmax
      checks: block
         call cross_iso_salinity(100.0_wp, .false., 0.0_wp, sa)   ! scalar khtr=100 (flux on)
         call cross_iso_salinity(0.0_wp, .true., 100.0_wp, sb)    ! scalar 0 + uniform ext=100
         call cross_iso_salinity(0.0_wp, .false., 0.0_wp, sc)     ! scalar 0, no ext => gate no-op (= IC)
         ! Non-triviality: the flux genuinely moved salinity (sa differs from
         ! the no-flux IC), so the equivalence below can't pass vacuously.
         call check(error, maxval(abs(sa - sc)) > 1.0e-6_wp, &
                    "redi flux must move salinity vs the no-flux gate (non-trivial)")
         if (allocated(error)) exit checks
         ! Wiring: uniform VarMix ext reproduces the scalar path exactly AND
         ! drives flux despite scalar khtr=0 (overrides the no-flux gate).
         dmax = maxval(abs(sa - sb))
         call check(error, dmax < 1.0e-12_wp, &
                    "uniform VarMix KhTr must reproduce scalar khtr exactly (and override khtr=0 gate)")
         deallocate (sa, sb, sc)
      end block checks
   end subroutine test_varmix_khtr

end module test_ocean_redi
