!! Analytic tests for the FV_MOM6 constant-by-layer (PCM) pressure gradient
!! at the IN-SITU density (`&ocean_pgf_nml insitu_density`, MOM6
!! `int_density_dz_generic_pcm` parity).
!!
!! 1. thermobaric_face_matches_insitu_oracle — flat, aligned z-level
!!    columns at rest (eta = 0) whose T/S vary in x such that the SURFACE
!!    potential density is horizontally uniform.  The in-situ density is
!!    not: the Wright expansion coefficients depend on pressure, so the
!!    compensation breaks with depth and the true pressure gradient grows
!!    downward.  Oracle: the layer-mean in-situ pressure difference from an
!!    independent fine trapezoid integration of `EOS(T, S, -g*rho0*z)`.
!!    The in-situ path must match it to 1e-6; the legacy potential-density
!!    path (`insitu_density = .false.`, `rho_layer` at p_ref = 0) sees a
!!    uniform density and returns ZERO — this is the fails-before half, the
!!    defect that held the global 1-degree Drake Passage transport at about
!!    half of MOM6's.
!! 2. rest_partial_steps_zero_pgf — horizontally uniform (per layer) T/S
!!    stratification, z-level layers with partial bottom cells of
!!    different depth: the compressible in-situ density varies strongly in
!!    z, and the tilted bed edges must not turn that into a spurious
!!    acceleration (the reason the in-situ path integrates across the face
!!    by Boole's rule instead of the trapezoid).  With and without the
!!    near-bottom mass weighting.
!! 3. linear_eos_bit_identical — the knob is inert for the linear EOS
!!    (in-situ and potential density coincide): on and off agree to the
!!    bit.
module test_ocean_pgf_insitu
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_insitu_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 6, NYP = 4
   real(wp), parameter :: DX = 5.0e4_wp
   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: RHO_REF = 1035.0_wp
   real(wp), parameter :: H_FILL = 1.0e-4_wp
      !! Inert filler thickness below a partial bottom cell (< H_VANISHED).

contains

   subroutine collect_ocean_pgf_insitu_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("thermobaric_face_matches_insitu_oracle", test_thermobaric), &
                  new_unittest("rest_partial_steps_zero_pgf", test_rest_partial_steps), &
                  new_unittest("linear_eos_bit_identical", test_linear_bitident) &
                  ]
   end subroutine collect_ocean_pgf_insitu_tests

   subroutine make_wright(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%rho0 = RHO0
      eos%is_init = .true.
   end subroutine make_wright

   subroutine run_pgf(grid, ms, pgf, eos)
      !! One compute pass on the device (or host), results back on the host.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(eos_t), intent(in) :: eos
      type(ocean_metrics_t) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   subroutine make_pgf(pgf, grid, nz, b, insitu, mass_weight)
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: b(:, :)
      logical, intent(in) :: insitu, mass_weight
      call pgf%init(grid, nz_ml=nz)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO0
      pgf%rho_ref = RHO_REF
      pgf%insitu_density = insitu
      pgf%mass_weight = mass_weight
      call pgf%set_bathymetry(b)
   end subroutine make_pgf

   pure function salinity_for_sigma0(eos, t, rho_target) result(s)
      !! The salinity at which `EOS(t, s, 0) = rho_target` (Newton).
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: t, rho_target
      real(wp) :: s, f, df
      integer :: it
      s = 34.7_wp
      do it = 1, 30
         f = eos_density_point(eos, t, s, 0.0_wp) - rho_target
         df = (eos_density_point(eos, t, s + 1.0e-3_wp, 0.0_wp) - &
               eos_density_point(eos, t, s - 1.0e-3_wp, 0.0_wp))/2.0e-3_wp
         s = s - f/df
      end do
   end function salinity_for_sigma0

   pure function layer_mean_pressure(eos, t, s, depth_top, dz) result(pbar)
      !! Mean over [depth_top, depth_top + dz] of the in-situ pressure
      !! anomaly p(d) = g * int_0^d (rho(t, s, g*rho0*d') - rho_ref) dd'
      !! for a column of UNIFORM t/s, by a fine trapezoid march — the
      !! oracle, independent of the kernel's Boole rules.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: t, s, depth_top, dz
      real(wp) :: pbar
      integer, parameter :: NSUB = 4000
      real(wp) :: d, dd, p, r_prev, r_now, acc, p_prev
      integer :: n, n_top
      ! March from the surface to depth_top, then average through the layer.
      n_top = nint(depth_top/dz*real(NSUB, wp))
      dd = dz/real(NSUB, wp)
      p = 0.0_wp
      d = 0.0_wp
      r_prev = eos_density_point(eos, t, s, 0.0_wp) - RHO_REF
      do n = 1, n_top
         d = real(n, wp)*dd
         r_now = eos_density_point(eos, t, s, GRAVITY*RHO0*d) - RHO_REF
         p = p + GRAVITY*0.5_wp*(r_prev + r_now)*dd
         r_prev = r_now
      end do
      acc = 0.0_wp
      do n = 1, NSUB
         p_prev = p
         d = depth_top + real(n, wp)*dd
         r_now = eos_density_point(eos, t, s, GRAVITY*RHO0*d) - RHO_REF
         p = p + GRAVITY*0.5_wp*(r_prev + r_now)*dd
         r_prev = r_now
         acc = acc + 0.5_wp*(p_prev + p)*dd
      end do
      pbar = acc/dz
   end function layer_mean_pressure

   subroutine test_thermobaric(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 20
      real(wp), parameter :: HLAY = 200.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_new, pgf_old
      type(eos_t) :: eos
      real(wp), allocatable :: b(:, :), tcol(:), scol(:)
      real(wp) :: rho_target, oracle, err_new, err_old, omax, vmax, depth_top
      integer :: i, j, k, nx, ny, jc
      character(len=160) :: msg
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_wright(eos)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (b(nx, ny), source=real(NZ, wp)*HLAY)
         allocate (tcol(nx), scol(nx))
         ! Cold/fresh to warm/salty in x, compensated at the surface.
         rho_target = eos_density_point(eos, 1.0_wp, 34.7_wp, 0.0_wp)
         do i = 1, nx
            tcol(i) = 1.0_wp + 1.5_wp*real(i - 1, wp)
            scol(i) = salinity_for_sigma0(eos, tcol(i), rho_target)
         end do
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = HLAY
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = tcol(i)*HLAY
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = scol(i)*HLAY
                  ! Legacy consumer: potential density at p_ref = 0.
                  ms%rho_layer(i, j, k) = eos_density_point(eos, tcol(i), scol(i), 0.0_wp)
               end do
            end do
         end do

         call make_pgf(pgf_new, grid, NZ, b, .true., .false.)
         call run_pgf(grid, ms, pgf_new, eos)
         call make_pgf(pgf_old, grid, NZ, b, .false., .false.)
         call run_pgf(grid, ms, pgf_old, eos)

         jc = NGHOST + NYP/2
         err_new = 0.0_wp
         err_old = 0.0_wp
         omax = 0.0_wp
         do k = 1, NZ
            ! k = 1 is the bed layer; its top sits (NZ - k) layers down.
            depth_top = real(NZ - k, wp)*HLAY
            do i = NGHOST + 2, NGHOST + NXP
               oracle = -(layer_mean_pressure(eos, tcol(i), scol(i), depth_top, HLAY) - &
                          layer_mean_pressure(eos, tcol(i - 1), scol(i - 1), depth_top, HLAY)) &
                        /(RHO0*DX)
               omax = max(omax, abs(oracle))
               err_new = max(err_new, abs(pgf_new%dpdx_face%data(i, jc, k) - oracle))
               err_old = max(err_old, abs(pgf_old%dpdx_face%data(i, jc, k) - oracle))
            end do
         end do
         vmax = maxval(abs(pgf_new%dpdy_face%data(NGHOST + 2:NGHOST + NXP, &
                                                  NGHOST + 2:NGHOST + NYP, :)))
         write (msg, '(a,es10.3,a,es10.3,a,es10.3,a,es10.3)') "max|oracle| ", omax, &
            " err in-situ ", err_new, " err legacy ", err_old, " max|PFv| ", vmax
         ! The thermobaric signal is real and sizeable at 4000 m.
         call check(error, omax > 1.0e-7_wp, "thermobaric: oracle PGF unexpectedly small: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, err_new <= 1.0e-6_wp*omax, &
                    "thermobaric: in-situ PCM PGF does not match the in-situ oracle: "//trim(msg))
         if (allocated(error)) exit checks
         ! Fails-before: the potential-density integral sees a uniform
         ! density and misses the whole signal.
         call check(error, err_old >= 0.5_wp*omax, &
                    "thermobaric: legacy potential-density path unexpectedly matches: "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, vmax <= 1.0e-12_wp, "thermobaric: spurious meridional PGF: "//trim(msg))
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf_new%destroy()
      call pgf_old%destroy()
      call ms%destroy()
   end subroutine test_thermobaric

   subroutine seed_zlevel_steps(ms, eos, nz, hlay, depth, b)
      !! z-level layers of thickness `hlay` down to each column's depth,
      !! a PARTIAL bottom cell, inert fillers below (h = H_FILL).  T/S are
      !! a function of the layer INDEX only (horizontally uniform per
      !! layer), stratified in the vertical.
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      integer, intent(in) :: nz
      real(wp), intent(in) :: hlay
      real(wp), intent(in) :: depth(:)
      real(wp), intent(inout) :: b(:, :)
      integer :: i, j, k, kk
      real(wp) :: z_top, tk, sk, hk
      do j = 1, size(ms%h_layer, 2)
         do i = 1, size(ms%h_layer, 1)
            b(i, j) = 0.0_wp
            do k = nz, 1, -1
               kk = nz - k                        ! layers above this one
               z_top = real(kk, wp)*hlay
               hk = min(hlay, max(depth(i) - z_top, 0.0_wp))
               if (hk <= 0.0_wp) hk = H_FILL
               tk = 2.0_wp + 16.0_wp*exp(-(z_top + 0.5_wp*hlay)/700.0_wp)
               sk = 34.6_wp + 0.5_wp*exp(-(z_top + 0.5_wp*hlay)/900.0_wp)
               ms%h_layer(i, j, k) = hk
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = tk*hk
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = sk*hk
               ms%rho_layer(i, j, k) = eos_density_point(eos, tk, sk, 0.0_wp)
               b(i, j) = b(i, j) + hk
            end do
         end do
      end do
   end subroutine seed_zlevel_steps

   subroutine test_rest_partial_steps(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 16
      real(wp), parameter :: HLAY = 250.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(eos_t) :: eos
      real(wp), allocatable :: b(:, :), depth(:)
      real(wp) :: worst(2), pmax
      integer :: i, j, k, nx, ny, imw
      logical :: mw
      character(len=96) :: msg
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call make_wright(eos)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (b(nx, ny), source=0.0_wp)
         allocate (depth(nx))
         ! Column depths with partial bottom cells of every size.
         do i = 1, nx
            depth(i) = 2600.0_wp + 137.0_wp*real(mod(7*i, 10), wp)
         end do
         call seed_zlevel_steps(ms, eos, NZ, HLAY, depth, b)
         do imw = 1, 2
            mw = imw == 2
            call make_pgf(pgf, grid, NZ, b, .true., mw)
            call run_pgf(grid, ms, pgf, eos)
            pmax = 0.0_wp
            do k = 1, NZ
               do j = NGHOST + 1, NGHOST + NYP
                  do i = NGHOST + 2, NGHOST + NXP
                     ! Open faces only: a filler on either side is a
                     ! closed (masked) face in the dynamics.
                     if (ms%h_layer(i - 1, j, k) > H_FILL .and. &
                         ms%h_layer(i, j, k) > H_FILL) then
                        pmax = max(pmax, abs(pgf%dpdx_face%data(i, j, k)))
                     end if
                  end do
               end do
            end do
            worst(imw) = pmax
            call pgf%destroy()
         end do
         write (msg, '(a,es10.3,a,es10.3)') "max|PFu| no-mw ", worst(1), " mw ", worst(2)
         call check(error, worst(1) <= 1.0e-10_wp, &
                    "rest steps: spurious in-situ PGF at partial cells (no mass weight): "//trim(msg))
         if (allocated(error)) exit checks
         call check(error, worst(2) <= 1.0e-10_wp, &
                    "rest steps: spurious in-situ PGF at partial cells (mass weight): "//trim(msg))
      end block checks
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine test_rest_partial_steps

   subroutine test_linear_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 16
      real(wp), parameter :: HLAY = 250.0_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_on, pgf_off
      type(eos_t) :: eos
      real(wp), allocatable :: b(:, :), depth(:)
      real(wp) :: dmax
      integer :: i, nx, ny
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         eos%variant = EOS_VARIANT_LINEAR
         eos%rho0 = RHO0
         eos%alpha_T = 0.2_wp
         eos%beta_S = 0.78_wp
         eos%is_init = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (b(nx, ny), source=0.0_wp)
         allocate (depth(nx))
         do i = 1, nx
            depth(i) = 2600.0_wp + 137.0_wp*real(mod(7*i, 10), wp)
         end do
         call seed_zlevel_steps(ms, eos, NZ, HLAY, depth, b)
         ! Break the horizontal uniformity so the PGF is not zero.
         do i = 1, nx
            ms%tracers(ms%idx_temperature)%hTr(i, :, :) = &
               ms%tracers(ms%idx_temperature)%hTr(i, :, :) + 0.1_wp*real(i, wp)*ms%h_layer(i, :, :)
            ms%rho_layer(i, :, :) = ms%rho_layer(i, :, :) - 0.02_wp*real(i, wp)
         end do
         call make_pgf(pgf_on, grid, NZ, b, .true., .false.)
         call run_pgf(grid, ms, pgf_on, eos)
         call make_pgf(pgf_off, grid, NZ, b, .false., .false.)
         call run_pgf(grid, ms, pgf_off, eos)
         dmax = max(maxval(abs(pgf_on%dpdx_face%data - pgf_off%dpdx_face%data)), &
                    maxval(abs(pgf_on%dpdy_face%data - pgf_off%dpdy_face%data)))
         call check(error, dmax == 0.0_wp, &
                    "linear EOS: insitu_density on/off not bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(pgf_on%dpdx_face%data)) > 0.0_wp, &
                    "linear EOS: the comparison state has no PGF at all")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf_on%destroy()
      call pgf_off%destroy()
      call ms%destroy()
   end subroutine test_linear_bitident

end module test_ocean_pgf_insitu
