!! Analytical tests for file-based SSS restoring (`rdb_ocean_sss_restore`):
!! per-cell clamp, global net-zero adjustment (scaling and area-mean),
!! the file -> reader -> target -> deposit chain, salt conservation, the
!! SP -> SA conversion and the configure refusals.
!!
!! Every oracle is written from the physical definition, not from the
!! kernel: with piston velocity `w` and clamp `c`, a wet cell with surface
!! salinity `S` and target `T` gains, per thermo step, `dt * w * dS`
!! of `S * h` (`dS = clamp(T - S, +-c)`), times the net-adjustment factor.
!! In the model's `Q_salt` unit convention that is the flux
!! `F = rho0 * w * dS` (SI salt mass `1e-3 * F`) deposited as
!! `(dt/rho0) * F`.
!!
!!   * clamp        : dS = +-0.8 with cap 0.5 gives +-0.5, 0.3 unclamped.
!!   * scale_pos/neg: 4 cells with known areas; sum A F = 0, signs and
!!                    zero cells unchanged, only the larger side scaled by
!!                    exactly min(P,-N)/max(P,-N); P = 0 unchanged.
!!   * mean         : area-mean subtraction nets to round-off.
!!   * file_blend   : a 2-record cyclic NetCDF, queried mid-way and on the
!!                    wrap seam, gives the linear blend as the target.
!!   * conservation : closed basin; scaling keeps total salt, no scaling
!!                    changes it by exactly the integrated flux.
!!   * practical    : `sss_convention = practical` multiplies by the
!!                    Reference-Salinity factor 35.16504/35.
!!   * refusals     : `validate_config` and the reader refuse the bad
!!                    combinations.
!! NetCDF I/O: runs with OMP_NUM_THREADS=1.  GPU: every array the kernels
!! touch is mapped (state, surface flux incl. the target, metrics).
module test_ocean_sss_restore
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_config, only: config_t, read_config_from_string, validate_config, &
                         ocean_restore_config_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_ocean_data_forcing, only: ocean_data_forcing_t, ocean_data_forcing_register_sss, &
                                     ocean_data_forcing_apply
   use rdb_ocean_sss_restore, only: ocean_sss_restore_apply, sss_net_factors, &
                                    SSS_NET_NONE, SSS_NET_MEAN, SSS_NET_SCALE
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            rdb_def_var_1d, nc_enddef, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_sss_restore_tests

   integer, parameter :: NG = 2
   integer, parameter :: NXP = 6, NYP = 4
   integer, parameter :: NZ = 2
   real(wp), parameter :: H0 = 10.0_wp
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: DT = 1800.0_wp
   real(wp), parameter :: PISTON_DAY = 8.64_wp
   real(wp), parameter :: W = PISTON_DAY/86400.0_wp   ! 1e-4 m/s
   real(wp), parameter :: SR_FACTOR = 35.16504_wp/35.0_wp

contains

   subroutine collect_ocean_sss_restore_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sss_restore_clamp", test_clamp), &
                  new_unittest("sss_restore_scaling", test_scaling), &
                  new_unittest("sss_restore_mean_subtract", test_mean), &
                  new_unittest("sss_restore_file_blend", test_file_blend), &
                  new_unittest("sss_restore_conservation", test_conservation), &
                  new_unittest("sss_restore_practical_convention", test_practical), &
                  new_unittest("sss_restore_refusals", test_refusals) &
                  ]
   end subroutine collect_ocean_sss_restore_tests

   ! -----------------------------------------------------------------
   ! Harness
   ! -----------------------------------------------------------------

   subroutine build(grid, metrics, ms, eos, sf, net_mode, max_delta, conv)
      !! Plane with uniform `S = S0`, all wet, unit-area metrics, the
      !! file-restoring state configured.  Nothing is device-mapped yet
      !! except the metrics (the caller edits host arrays first).
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: net_mode
      real(wp), intent(in) :: max_delta
      logical, intent(in) :: conv
      call comm_once()
      call grid%init(NXP, NYP, NG, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      call sf%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*H0
      ms%tracers(ms%idx_salinity)%hTr = S0*H0
      call sf%set_restore(.false., .true., 0.0_wp, PISTON_DAY, 0.0_wp, 0.0_wp)
      call sf%set_restore_file(grid, PISTON_DAY, max_delta, net_mode, conv)
      sf%sss_target = S0
   end subroutine build

   subroutine map_all(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
   end subroutine map_all

   subroutine pull_unmap(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      associate (hS => ms%tracers(ms%idx_salinity)%hTr, sbud => ms%salt_budget_surface)
         !$acc update self(hS, sbud)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine pull_unmap

   subroutine set_area(metrics, a)
      type(ocean_metrics_t), intent(inout) :: metrics
      real(wp), intent(in) :: a(:, :)
      metrics%areaT = a
      !$acc update device(metrics%areaT)
   end subroutine set_area

   subroutine teardown(metrics, ms, sf)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      call sf%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown

   pure function inc_at(ms, i, j) result(inc)
      !! Deposit of one cell, read from the salt budget contributor (which
      !! starts at zero, so there is no cancellation against `S*h = 350`).
      !! `test_clamp` separately pins that the tracer moved by the same.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: i, j
      real(wp) :: inc
      inc = ms%salt_budget_surface(i, j, NZ)
   end function inc_at

   ! -----------------------------------------------------------------
   ! 1. Clamp
   ! -----------------------------------------------------------------

   subroutine test_clamp(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: CAP = 0.5_wp
      real(wp) :: tol
      tol = 1.0e-13_wp*DT*W*CAP

      call build(grid, metrics, ms, eos, sf, SSS_NET_NONE, CAP, .false.)
      sf%sss_target(NG + 1, NG + 1) = S0 + 0.8_wp    ! clamped up
      sf%sss_target(NG + 2, NG + 1) = S0 - 0.8_wp    ! clamped down
      sf%sss_target(NG + 3, NG + 1) = S0 + 0.3_wp    ! below the cap
      sf%sss_target(NG + 4, NG + 1) = S0             ! at target
      call map_all(ms, sf)
      call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
      call pull_unmap(ms, sf)

      checks: block
         call check(error, abs(inc_at(ms, NG + 1, NG + 1) - DT*W*CAP) < tol, &
                    "dS=+0.8 must clamp to +cap")
         if (allocated(error)) exit checks
         call check(error, abs(inc_at(ms, NG + 2, NG + 1) + DT*W*CAP) < tol, &
                    "dS=-0.8 must clamp to -cap")
         if (allocated(error)) exit checks
         call check(error, abs(inc_at(ms, NG + 3, NG + 1) - DT*W*0.3_wp) < tol, &
                    "dS=0.3 must pass unclamped")
         if (allocated(error)) exit checks
         call check(error, inc_at(ms, NG + 4, NG + 1) == 0.0_wp, &
                    "a cell at target must not change")
         if (allocated(error)) exit checks
         ! The budget contributor mirrors the deposit.
         call check(error, abs((ms%tracers(ms%idx_salinity)%hTr(NG + 1, NG + 1, NZ) - S0*H0) &
                               - inc_at(ms, NG + 1, NG + 1)) < 1.0e-12_wp, &
                    "the tracer must move by exactly what salt_budget_surface records")
      end block checks
      call teardown(metrics, ms, sf)
   end subroutine test_clamp

   ! -----------------------------------------------------------------
   ! 2. Net-zero by scaling
   ! -----------------------------------------------------------------

   subroutine test_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: fp, fn
      integer :: pass

      ! The factor predicate itself.
      call sss_net_factors(2.0_wp, -1.0_wp, fp, fn)
      call check(error, fp == 0.5_wp .and. fn == 1.0_wp, "P=2,N=-1 scales positives by 0.5")
      if (allocated(error)) return
      call sss_net_factors(1.0_wp, -4.0_wp, fp, fn)
      call check(error, fp == 1.0_wp .and. fn == 0.25_wp, "P=1,N=-4 scales negatives by 0.25")
      if (allocated(error)) return
      call sss_net_factors(0.0_wp, -3.0_wp, fp, fn)
      call check(error, fp == 1.0_wp .and. fn == 1.0_wp, "P=0 must leave F unchanged")
      if (allocated(error)) return
      call sss_net_factors(3.0_wp, 0.0_wp, fp, fn)
      call check(error, fp == 1.0_wp .and. fn == 1.0_wp, "N=0 must leave F unchanged")
      if (allocated(error)) return

      do pass = 1, 3
         call scaling_case(pass, error)
         if (allocated(error)) return
      end do
   end subroutine test_scaling

   subroutine scaling_case(pass, error)
      !! pass 1: positives dominate; 2: negatives dominate; 3: P = 0.
      !! Four wet cells at (3..6, 3) with areas 2, 1, 3, 5.
      integer, intent(in) :: pass
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: AR(4) = [2.0_wp, 1.0_wp, 3.0_wp, 5.0_wp]
      real(wp) :: ds(4), fac(4), a(NXP + 2*NG, NYP + 2*NG), exp_inc, net, mag
      integer :: c, i, j
      real(wp) :: fpos, fneg

      select case (pass)
      case (1)
         ds = [0.4_wp, 0.2_wp, -0.1_wp, 0.0_wp]
      case (2)
         ds = [-0.4_wp, -0.2_wp, 0.1_wp, 0.0_wp]
      case default
         ds = [-0.4_wp, -0.2_wp, -0.1_wp, 0.0_wp]
      end select
      ! Hand values.  pass 1: P = 2*0.4+1*0.2 = 1.0, N = -3*0.1 = -0.3 =>
      ! positives x 0.3.  pass 2: mirrored => negatives x 0.3.  pass 3: P=0.
      fpos = 1.0_wp
      fneg = 1.0_wp
      if (pass == 1) fpos = 0.3_wp
      if (pass == 2) fneg = 0.3_wp
      do c = 1, 4
         fac(c) = merge(fpos, merge(fneg, 1.0_wp, ds(c) < 0.0_wp), ds(c) > 0.0_wp)
      end do

      call build(grid, metrics, ms, eos, sf, SSS_NET_SCALE, 999.0_wp, .false.)
      ms%wet_mask = 0.0_wp
      a = 1.0_wp
      do c = 1, 4
         ms%wet_mask(NG + c, NG + 1) = 1.0_wp
         a(NG + c, NG + 1) = AR(c)
         sf%sss_target(NG + c, NG + 1) = S0 + ds(c)
      end do
      call set_area(metrics, a)
      call map_all(ms, sf)
      call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
      call pull_unmap(ms, sf)

      checks: block
         net = 0.0_wp
         mag = 0.0_wp
         do c = 1, 4
            exp_inc = DT*W*ds(c)*fac(c)
            call check(error, abs(inc_at(ms, NG + c, NG + 1) - exp_inc) <= 1.0e-13_wp*DT*W*0.4_wp, &
                       "cell deposit must be dt*w*dS*factor")
            if (allocated(error)) exit checks
            net = net + AR(c)*inc_at(ms, NG + c, NG + 1)
            mag = mag + AR(c)*abs(inc_at(ms, NG + c, NG + 1))
         end do
         if (pass < 3) then
            call check(error, abs(net) <= 1.0e-15_wp*mag, &
                       "sum A F must vanish after scaling (relative 1e-15)")
            if (allocated(error)) exit checks
            ! Only the larger-magnitude side is scaled, by exactly the ratio.
         else
            call check(error, fpos == 1.0_wp .and. fneg == 1.0_wp, "P=0: unchanged")
            if (allocated(error)) exit checks
         end if
         ! Signs preserved, the zero cell and every dry cell untouched.
         call check(error, inc_at(ms, NG + 4, NG + 1) == 0.0_wp, "zero cell stays zero")
         if (allocated(error)) exit checks
         do j = NG + 1, NG + NYP
            do i = NG + 1, NG + NXP
               if (j == NG + 1 .and. i <= NG + 4) cycle
               call check(error, inc_at(ms, i, j) == 0.0_wp, "dry cell must not change")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
      call teardown(metrics, ms, sf)
   end subroutine scaling_case

   ! -----------------------------------------------------------------
   ! 3. Net-zero by subtracting the area mean
   ! -----------------------------------------------------------------

   subroutine test_mean(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: AR(4) = [2.0_wp, 1.0_wp, 3.0_wp, 5.0_wp]
      real(wp), parameter :: DS(4) = [0.4_wp, 0.2_wp, -0.1_wp, 0.0_wp]
      real(wp) :: a(NXP + 2*NG, NYP + 2*NG), mean_ds, net, mag
      integer :: c

      call build(grid, metrics, ms, eos, sf, SSS_NET_MEAN, 999.0_wp, .false.)
      ms%wet_mask = 0.0_wp
      a = 1.0_wp
      do c = 1, 4
         ms%wet_mask(NG + c, NG + 1) = 1.0_wp
         a(NG + c, NG + 1) = AR(c)
         sf%sss_target(NG + c, NG + 1) = S0 + DS(c)
      end do
      call set_area(metrics, a)
      call map_all(ms, sf)
      call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
      call pull_unmap(ms, sf)

      mean_ds = sum(AR*DS)/sum(AR)      ! 0.9/11 : wet cells only, zero cell included
      net = 0.0_wp
      mag = 0.0_wp
      checks: block
         do c = 1, 4
            call check(error, abs(inc_at(ms, NG + c, NG + 1) - DT*W*(DS(c) - mean_ds)) &
                       <= 1.0e-13_wp*DT*W*0.4_wp, "deposit must be dt*w*(dS - mean dS)")
            if (allocated(error)) exit checks
            net = net + AR(c)*inc_at(ms, NG + c, NG + 1)
            mag = mag + AR(c)*abs(inc_at(ms, NG + c, NG + 1))
         end do
         call check(error, abs(net) <= 1.0e-14_wp*mag, "net must vanish to round-off")
         if (allocated(error)) exit checks
         call check(error, inc_at(ms, NG + 5, NG + 1) == 0.0_wp, &
                    "a dry cell must not receive the mean shift")
      end block checks
      call teardown(metrics, ms, sf)
   end subroutine test_mean

   ! -----------------------------------------------------------------
   ! 4. File -> reader -> target -> deposit
   ! -----------------------------------------------------------------

   subroutine write_two_records(filename, v0, v1, t0, t1, varname)
      character(len=*), intent(in) :: filename, varname
      real(wp), intent(in) :: v0(NXP, NYP), v1(NXP, NYP), t0, t1
      integer :: ncid, dx, dy, dtm, vf, vt, ierr
      real(wp) :: f(NXP, NYP, 2)
      f(:, :, 1) = v0
      f(:, :, 2) = v1
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", NXP, dx)
      call nc_def_dim(ncid, "y", NYP, dy)
      call nc_def_dim(ncid, "time", 2, dtm)
      call nc_def_var_3d(ncid, varname, [dx, dy, dtm], vf)
      call rdb_def_var_1d(ncid, "time", dtm, vt)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vf, f)
      call rdb_put_var_1d(ncid, vt, [t0, t1])
      call nc_close(ncid)
   end subroutine write_two_records

   subroutine comm_once()
      !! Comms must be up before any apply: the net-adjustment path does a
      !! global EFP allreduce.
      logical, save :: comm_inited = .false.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine comm_once

   subroutine init_halo()
      type(decomp_t) :: decomp
      call comm_once()
      call decomp_init(decomp, NXP, NYP, 1, 1, 0)
      call ocean_halo_init(decomp, NG, .false., .false.)
   end subroutine init_halo

   subroutine test_file_blend(error)
      !! Records at t = 0 and 100 s of a 200 s cycle.  t = 50 is the
      !! linear midpoint; t = 125 is on the (last -> first) wrap seam, a
      !! quarter of the way back: 0.75 v1 + 0.25 v0.  After each
      !! restoring step the surface salt change must be dt*w*(T - S).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      type(ocean_surface_stress_t) :: ss
      type(ocean_bc_state_t) :: bc
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_restore_config_t) :: rc
      real(wp) :: v0(NXP, NYP), v1(NXP, NYP), tq(2), wt(2), tgt, s_before, expect
      integer :: ip, jp, i, j, k, ierr
      character(len=*), parameter :: FN = "test_sss_restore_two_records.nc"

      do jp = 1, NYP
         do ip = 1, NXP
            v0(ip, jp) = 34.0_wp + 0.05_wp*real(ip, wp) + 0.03_wp*real(jp, wp)
            v1(ip, jp) = v0(ip, jp) + 0.4_wp
         end do
      end do
      call write_two_records(FN, v0, v1, 0.0_wp, 100.0_wp, "salt")
      tq = [50.0_wp, 125.0_wp]
      wt = [0.5_wp, 0.75_wp]       ! weight of record 2

      call build(grid, metrics, ms, eos, sf, SSS_NET_NONE, 999.0_wp, .false.)
      call init_halo()
      call ss%init(grid)
      call reader%init()
      rc%sss_file = FN
      rc%sss_var = "salt"
      rc%sss_time_cyclic = .true.
      rc%sss_cycle_period = 200.0_wp
      call ocean_data_forcing_register_sss(rc, reader, grid, sf, df, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK .and. df%active .and. df%id_sss > 0, &
                 "registration of a valid file must succeed")
      if (allocated(error)) return

      !$acc enter data copyin(reader)
      call reader%enter_data()
      call map_all(ms, sf)
      checks: do k = 1, 2
         call ocean_data_input_update_all(reader, tq(k))
         call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, tq(k))
         call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
         associate (hS => ms%tracers(ms%idx_salinity)%hTr)
            !$acc update self(hS)
            do j = NG + 1, NG + NYP
               do i = NG + 1, NG + NXP
                  ip = i - NG
                  jp = j - NG
                  tgt = (1.0_wp - wt(k))*v0(ip, jp) + wt(k)*v1(ip, jp)
                  if (k == 1) then
                     s_before = S0
                  else
                     s_before = S0 + (DT*W*(0.5_wp*(v0(ip, jp) + v1(ip, jp)) - S0))/H0
                  end if
                  expect = (s_before*H0 + DT*W*(tgt - s_before))/H0
                  call check(error, abs(hS(i, j, NZ)/H0 - expect) <= 1.0e-12_wp*S0, &
                             "surface salinity must follow the piston law toward "// &
                             "the time-blended file target")
                  if (allocated(error)) exit checks
               end do
            end do
         end associate
      end do checks
      call pull_unmap(ms, sf)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call reader%destroy()
      call ocean_halo_destroy()
      call ss%destroy()
      call teardown(metrics, ms, sf)
   end subroutine test_file_blend

   ! -----------------------------------------------------------------
   ! 5. Conservation in a closed basin
   ! -----------------------------------------------------------------

   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dt_scaled, dt_none, integ_none, scale_ref

      call conservation_run(SSS_NET_SCALE, dt_scaled, integ_none, scale_ref)
      call conservation_run(SSS_NET_NONE, dt_none, integ_none, scale_ref)
      call check(error, abs(dt_scaled) < 1.0e-12_wp*scale_ref, &
                 "scaling: total salt must change by < 1e-12 relative over 10 steps")
      if (allocated(error)) return
      call check(error, abs(dt_none) > 1.0e-6_wp*scale_ref, &
                 "control: without the adjustment the total must actually move")
      if (allocated(error)) return
      call check(error, abs(dt_none - integ_none) < 1.0e-12_wp*abs(integ_none), &
                 "no adjustment: total salt changes by exactly the integrated flux")
   end subroutine test_conservation

   subroutine conservation_run(net_mode, dtotal, integ, scale_ref)
      !! Total `sum A hS` change over 10 steps (`dtotal`), the integral of
      !! the budget contributor (`integ`), and the initial total.
      integer, intent(in) :: net_mode
      real(wp), intent(out) :: dtotal, integ, scale_ref
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: a(NXP + 2*NG, NYP + 2*NG), t_before, t_after
      integer :: i, j, s

      call build(grid, metrics, ms, eos, sf, net_mode, 999.0_wp, .false.)
      do j = 1, NYP + 2*NG
         do i = 1, NXP + 2*NG
            a(i, j) = 1.0_wp + 0.1_wp*real(i, wp) + 0.05_wp*real(j, wp)
            sf%sss_target(i, j) = S0 + 0.4_wp*sin(1.3_wp*real(i, wp) + 0.7_wp*real(j, wp)) + 0.1_wp
         end do
      end do
      call set_area(metrics, a)
      t_before = 0.0_wp
      do j = NG + 1, NG + NYP
         do i = NG + 1, NG + NXP
            t_before = t_before + a(i, j)*ms%tracers(ms%idx_salinity)%hTr(i, j, NZ)
         end do
      end do
      scale_ref = t_before
      call map_all(ms, sf)
      do s = 1, 10
         call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
      end do
      call pull_unmap(ms, sf)
      t_after = 0.0_wp
      integ = 0.0_wp
      do j = NG + 1, NG + NYP
         do i = NG + 1, NG + NXP
            t_after = t_after + a(i, j)*ms%tracers(ms%idx_salinity)%hTr(i, j, NZ)
            integ = integ + a(i, j)*ms%salt_budget_surface(i, j, NZ)
         end do
      end do
      dtotal = t_after - t_before
      call teardown(metrics, ms, sf)
   end subroutine conservation_run

   ! -----------------------------------------------------------------
   ! 6. Practical-salinity file under cons_abs
   ! -----------------------------------------------------------------

   subroutine test_practical(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf

      call build(grid, metrics, ms, eos, sf, SSS_NET_NONE, 999.0_wp, .true.)
      sf%sss_target = 35.0_wp                      ! SP; model S is SA = 35.0 g/kg
      sf%sss_target(NG + 1, NG + 1) = 34.0_wp      ! SP 34 -> SA 34.1359...
      call map_all(ms, sf)
      call ocean_sss_restore_apply(grid, metrics, sf, ms, DT)
      call pull_unmap(ms, sf)
      call check(error, abs(inc_at(ms, NG + 2, NG + 1) - DT*W*(35.0_wp*SR_FACTOR - S0)) &
                 < 1.0e-12_wp, "SP=35 must restore toward SA = 35*35.16504/35")
      if (allocated(error)) return
      call check(error, abs(inc_at(ms, NG + 1, NG + 1) - DT*W*(34.0_wp*SR_FACTOR - S0)) &
                 < 1.0e-12_wp, "SP=34 must restore toward SA = 34*35.16504/35")
      call teardown(metrics, ms, sf)
   end subroutine test_practical

   ! -----------------------------------------------------------------
   ! 7. Refusals
   ! -----------------------------------------------------------------

   function validate_restore(extra) result(ierr)
      !! `validate_config` status of a minimal ocean namelist with the
      !! given `&ocean_restore_nml` body.
      character(len=*), intent(in) :: extra
      integer :: ierr
      type(config_t) :: cfg
      ierr = -999
      call read_config_from_string( &
         "&sim_nml sim_type = 'ocean' /"//new_line('a')// &
         "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"// &
         new_line('a')//"&ocean_restore_nml "//extra//" /"//new_line('a'), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr=ierr)
   end function validate_restore

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: OK = &
                                     "enable_restore_salt = .true., piston_s = 0.11, sss_file = 'x.nc', "// &
                                     "sss_var = 'asalt', sss_as_sflux = .true., "

      call check(error, validate_restore("sss_file = ''") == OCEAN_STATUS_OK, &
                 "defaults must validate")
      if (allocated(error)) return
      call check(error, validate_restore(OK//"sss_adjust_net_zero = .true., "// &
                                         "sss_adjust_by_scaling = .true., "// &
                                         "max_delta_sss = 0.5") == OCEAN_STATUS_OK, &
                 "the OM3 combination must validate")
      if (allocated(error)) return
      call check(error, validate_restore("enable_restore_salt = .false., piston_s = 0.11, "// &
                                         "sss_file = 'x.nc', sss_as_sflux = .true.") &
                 /= OCEAN_STATUS_OK, "sss_file without enable_restore_salt must be refused")
      if (allocated(error)) return
      call check(error, validate_restore("enable_restore_salt = .true., piston_s = 0.11, "// &
                                         "sss_file = 'x.nc', sss_as_sflux = .false.") &
                 /= OCEAN_STATUS_OK, "sss_file with sss_as_sflux = .false. must be refused")
      if (allocated(error)) return
      call check(error, validate_restore("enable_restore_salt = .true., piston_s = 0.0, "// &
                                         "sss_file = 'x.nc', sss_as_sflux = .true.") &
                 /= OCEAN_STATUS_OK, "sss_file with piston_s = 0 must be refused")
      if (allocated(error)) return
      call check(error, validate_restore(OK//"sss_adjust_by_scaling = .true.") &
                 /= OCEAN_STATUS_OK, "by_scaling without net_zero must be refused")
      if (allocated(error)) return
      call check(error, validate_restore(OK//"sss_convention = 'wrong'") &
                 /= OCEAN_STATUS_OK, "an unknown sss_convention must be refused")
      if (allocated(error)) return
      call check(error, validate_restore("max_delta_sss = 0.5") /= OCEAN_STATUS_OK, &
                 "max_delta_sss without sss_file must be refused")
      if (allocated(error)) return

      call reader_refusals(error)
   end subroutine test_refusals

   subroutine reader_refusals(error)
      !! A missing file and a missing variable are fatal at registration.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_restore_config_t) :: rc
      real(wp) :: v(NXP, NYP)
      integer :: ierr
      character(len=*), parameter :: FN = "test_sss_restore_refusal.nc"

      v = 35.0_wp
      call write_two_records(FN, v, v, 0.0_wp, 100.0_wp, "salt")
      call build(grid, metrics, ms, eos, sf, SSS_NET_NONE, 999.0_wp, .false.)
      call reader%init()
      rc%sss_time_cyclic = .true.
      rc%sss_cycle_period = 200.0_wp

      rc%sss_file = "test_sss_restore_no_such_file.nc"
      rc%sss_var = "salt"
      call ocean_data_forcing_register_sss(rc, reader, grid, sf, df, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, "a missing file must be refused")
      if (allocated(error)) then
         call teardown(metrics, ms, sf)
         return
      end if

      rc%sss_file = FN
      rc%sss_var = "no_such_variable"
      call ocean_data_forcing_register_sss(rc, reader, grid, sf, df, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, "a missing variable must be refused")
      call reader%destroy()
      call teardown(metrics, ms, sf)
   end subroutine reader_refusals

end module test_ocean_sss_restore
