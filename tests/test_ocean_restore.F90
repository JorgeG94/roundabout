!! Analytical unit tests for surface buoyancy restoring (MOM6
!! `RESTOREBUOY`): `ocean_surface_restore_apply_tracers` in
!! `rdb_ocean_surface_flux`.
!!
!! The kernel relaxes the top-layer (`k = nz`) T / S toward a scalar
!! target with a piston velocity `p` [m/s].  Per thermo step, in
!! hTr-space (no other forcing):
!!   T_{n+1} = T_n + (dt·p/h)·(T_target − T_n)
!! i.e. discrete forward-Euler relaxation at rate `lambda = p/h`.
!! Closed-form oracle after N steps:
!!   T_N = T_target + (T_0 − T_target)·(1 − lambda·dt)^N
!! All tests run the kernel ON DEVICE.
!!
!! Bottom-up convention: surface = `k = nz`, bed = `k = 1`.
!!
!! Cases:
!!   * RELAX_T: a single column with no other forcing relaxes SST toward
!!     T_target as the forward-Euler exponential (closed-form oracle);
!!     monotone approach + ~3 e-folds reached.
!!   * SIGN: T_0 < T_target warms; T_0 > T_target cools.
!!   * RELAX_S: the salt branch mirrors with `lambda = p_S/h`.
!!   * BUDGET: Σ_steps heat_budget_surface(nz) equals the accumulated
!!     restoring source.
!!   * OFF: enable_restore_* = .false. ⇒ the kernel makes no change
!!     (byte-for-byte) — the default-off bit-identity guard.
!!   * LAND: a wet_mask = 0 column with a target ≠ its T must NOT change.
module test_ocean_restore
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_restore_apply_tracers
   implicit none
   private

   public :: collect_ocean_restore_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H_LAYER = 10.0_wp    ! m per layer (h_top = 10 m)
   real(wp), parameter :: DT = 3600.0_wp       ! s
   ! Piston velocities (m/day → m/s in the setter).  Pick lambda·dt = O(0.1):
   ! p = 0.24 m/day = 2.7778e-6 m/s; lambda = p/h_top = 2.7778e-7 1/s;
   ! lambda·dt = 1.0e-3 per step (monotone regime).
   real(wp), parameter :: PISTON_T_DAY = 0.24_wp
   real(wp), parameter :: PISTON_S_DAY = 0.24_wp
   integer, parameter :: NSTEPS = 3000          ! ~3 e-folds (lambda·dt·N = 3)

contains

   subroutine collect_ocean_restore_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("restore_T_exponential", test_relax_T), &
                  new_unittest("restore_sign", test_sign), &
                  new_unittest("restore_S_exponential", test_relax_S), &
                  new_unittest("restore_budget_contributor", test_budget), &
                  new_unittest("restore_off_bit_identity", test_off), &
                  new_unittest("restore_land_no_op", test_land) &
                  ]
   end subroutine collect_ocean_restore_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   !> Set up a uniform column: surface tracers seeded to T0 / S0 (degC,
   !> PSU), all layers thickness H_LAYER, optional wet_mask override.
   subroutine setup(grid, ms, eos, sf, T0, S0)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: T0, S0
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      call sf%init(grid)
      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T0*H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S0*H_LAYER
   end subroutine setup

   !> Run the restore kernel NSTEPS times on device.  Returns the device
   !> state pulled back to the host on exit.
   subroutine run_restore(grid, sf, ms, n)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: n
      integer :: s
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      do s = 1, n
         call ocean_surface_restore_apply_tracers(grid, sf, ms, DT)
      end do
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 hbud => ms%heat_budget_surface, &
                 sbud => ms%salt_budget_surface)
         !$acc update self(hT, hS, hbud, sbud)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_restore

   !> Forward-Euler discrete oracle: X_N = X* + (X0 − X*)·(1 − lam·dt)^N.
   pure function euler_oracle(X0, Xstar, piston_day, n) result(xn)
      real(wp), intent(in) :: X0, Xstar, piston_day
      integer, intent(in) :: n
      real(wp) :: xn, lam_dt
      lam_dt = (piston_day/86400.0_wp)/H_LAYER*DT   ! (p/h)·dt
      xn = Xstar + (X0 - Xstar)*(1.0_wp - lam_dt)**n
   end function euler_oracle

   ! -----------------------------------------------------------------

   subroutine test_relax_T(error)
      !! SST relaxes toward T_target as the forward-Euler exponential.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: T0 = 10.0_wp, TSTAR = 20.0_wp
      real(wp) :: sst, oracle, efold_band
      integer :: i, j
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, T0, 35.0_wp)
         call sf%set_restore(.true., .false., PISTON_T_DAY, 0.0_wp, TSTAR, 0.0_wp)
         call run_restore(grid, sf, ms, NSTEPS)

         i = grid%nx_total/2
         j = grid%ny_total/2
         sst = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H_LAYER
         oracle = euler_oracle(T0, TSTAR, PISTON_T_DAY, NSTEPS)

         ! Same arithmetic as the kernel ⇒ tight tol (pins rate + sign).
         call check(error, abs(sst - oracle) < 1.0e-10_wp*abs(TSTAR) + 1.0e-12_wp, &
                    "SST must match the forward-Euler relaxation oracle")
         if (allocated(error)) exit checks
         ! Monotone approach: sst is between T0 and T* and closer to T*.
         call check(error, sst > T0 .and. sst < TSTAR, &
                    "SST must lie strictly between T0 and T_target")
         if (allocated(error)) exit checks
         ! ~3 e-folds: residual within (T0−T*)·exp(−3) band of T*.
         efold_band = abs(T0 - TSTAR)*exp(-3.0_wp)
         call check(error, abs(sst - TSTAR) < efold_band, &
                    "SST must relax ~3 e-folds toward T_target")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_relax_T

   subroutine test_sign(error)
      !! T0 < T* warms; T0 > T* cools — guards the (target − SST) sign.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: sst_warm, sst_cool
      integer :: i, j, nsmall
      checks: block
         nsmall = 50
         ! Warming case: T0 = 5 < T* = 15.
         call make_grid(grid)
         call setup(grid, ms, eos, sf, 5.0_wp, 35.0_wp)
         call sf%set_restore(.true., .false., PISTON_T_DAY, 0.0_wp, 15.0_wp, 0.0_wp)
         call run_restore(grid, sf, ms, nsmall)
         i = grid%nx_total/2; j = grid%ny_total/2
         sst_warm = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H_LAYER
         call sf%destroy(); call eos%destroy(); call ms%destroy()

         ! Cooling case: T0 = 25 > T* = 15.
         call setup(grid, ms, eos, sf, 25.0_wp, 35.0_wp)
         call sf%set_restore(.true., .false., PISTON_T_DAY, 0.0_wp, 15.0_wp, 0.0_wp)
         call run_restore(grid, sf, ms, nsmall)
         sst_cool = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H_LAYER

         call check(error, sst_warm > 5.0_wp, "T0 < T_target must warm the surface")
         if (allocated(error)) exit checks
         call check(error, sst_cool < 25.0_wp, "T0 > T_target must cool the surface")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_sign

   subroutine test_relax_S(error)
      !! SSS relaxes toward S_target with lambda = p_S / h.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: S0 = 35.0_wp, SSTAR = 36.0_wp
      real(wp) :: sss, oracle
      integer :: i, j
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, 15.0_wp, S0)
         call sf%set_restore(.false., .true., 0.0_wp, PISTON_S_DAY, 0.0_wp, SSTAR)
         call run_restore(grid, sf, ms, NSTEPS)

         i = grid%nx_total/2; j = grid%ny_total/2
         sss = ms%tracers(ms%idx_salinity)%hTr(i, j, NZ)/H_LAYER
         oracle = euler_oracle(S0, SSTAR, PISTON_S_DAY, NSTEPS)

         call check(error, abs(sss - oracle) < 1.0e-10_wp*abs(SSTAR) + 1.0e-12_wp, &
                    "SSS must match the forward-Euler relaxation oracle")
         if (allocated(error)) exit checks
         ! Temperature untouched (salt-only restoring).
         call check(error, &
                    abs(ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H_LAYER - 15.0_wp) &
                    == 0.0_wp, &
                    "salt-only restoring must leave SST unchanged")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_relax_S

   subroutine test_budget(error)
      !! Σ_steps heat_budget_surface(nz) equals the accumulated restoring
      !! source (final − initial hTr at the surface, since restoring is
      !! the only source and budget mirrors the same increment).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: T0 = 10.0_wp, TSTAR = 20.0_wp
      real(wp) :: bud, applied
      integer :: i, j
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, T0, 35.0_wp)
         call sf%set_restore(.true., .false., PISTON_T_DAY, 0.0_wp, TSTAR, 0.0_wp)
         call run_restore(grid, sf, ms, NSTEPS)

         i = grid%nx_total/2; j = grid%ny_total/2
         bud = ms%heat_budget_surface(i, j, NZ)
         ! Total applied increment to hTr = final − initial.
         applied = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ) - T0*H_LAYER
         call check(error, abs(bud - applied) < 1.0e-12_wp*abs(applied) + 1.0e-14_wp, &
                    "heat_budget_surface must record the restoring source")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_budget

   subroutine test_off(error)
      !! enable_restore_* = .false. ⇒ has_restore_* = .false. ⇒ the
      !! kernel early-returns; T / S byte-for-byte unchanged on device.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, 10.0_wp, 35.0_wp)
         ! Switches OFF but pistons + targets non-zero — must still no-op.
         call sf%set_restore(.false., .false., PISTON_T_DAY, PISTON_S_DAY, &
                             20.0_wp, 36.0_wp)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         call run_restore(grid, sf, ms, NSTEPS)
         call check(error, &
                    maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp &
                    .and. &
                    maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) == 0.0_wp, &
                    "restore OFF must leave T/S byte-for-byte unchanged")
      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(hS_ic)) deallocate (hS_ic)
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_off

   subroutine test_land(error)
      !! A wet_mask = 0 column with a target ≠ its T must NOT change —
      !! guards the mask multiply (no restoring on land).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: sst_land, sst_wet
      integer :: i, j
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, 10.0_wp, 35.0_wp)
         call sf%set_restore(.true., .false., PISTON_T_DAY, 0.0_wp, 20.0_wp, 0.0_wp)
         i = grid%nx_total/2; j = grid%ny_total/2
         ms%wet_mask = 1.0_wp
         ms%wet_mask(i, j) = 0.0_wp   ! make the probe column land
         call run_restore(grid, sf, ms, NSTEPS)

         sst_land = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ)/H_LAYER
         sst_wet = ms%tracers(ms%idx_temperature)%hTr(i + 1, j, NZ)/H_LAYER

         call check(error, sst_land == 10.0_wp, &
                    "land (wet_mask=0) column must not be restored")
         if (allocated(error)) exit checks
         call check(error, sst_wet > 10.0_wp, &
                    "neighbouring wet column must still be restored")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_land

end module test_ocean_restore
