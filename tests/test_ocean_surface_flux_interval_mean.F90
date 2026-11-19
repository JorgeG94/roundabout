!! Interval-mean surface fluxes at `&ocean_vmix_nml dt_therm_ratio = N > 1`.
!!
!! The thermodynamic step spends its surface flux over `N*dt`, so the assembled `Q_heat` / `Q_salt` /
!! `fw_net` must be built from the TIME MEAN of the atmospheric / data-forcing components received during
!! the interval, not from the single dynamic step that happened to precede the assembly; and `sw_net`
!! (what the penetration, KPP and EPBL read) must be the shortwave that is INSIDE `Q_heat`.  Wind stress
!! is not averaged.
!!
!! Gates:
!!   1. `ifm_assembly_of_means` (unit level, "virtual" and "mass"): N = 4 steps whose components differ
!!      on every step and in every cell; accumulate, apply the mean, assemble.  `Q_heat`, `Q_salt`,
!!      `fw_net`, `Q_salt_tracer` and `sw_net` equal those of a slot that was handed the arithmetic
!!      MEANS directly, and differ from the assembly of the last sample alone (non-vacuous).
!!   2. `ifm_energy_closure`: the heat the ocean receives over the interval, `Q_heat * N * dt`, equals
!!      `dt * sum over the steps` of the per-step open-water heat (bands, enthalpy of the mass fluxes,
!!      `cp*SST*evap`), to round-off.
!!   3. `ifm_accumulators_empty_and_ratio1_inert`: after the assembly the sums and the count are zero;
!!      with the treatment off (ratio 1) nothing is allocated and the accumulate / mean / clear calls
!!      leave every component bitwise untouched.
!!   4. `ifm_engine_wiring`: through the real engine at ratio 4, for `diabatic_first` true and false:
!!      the assembly fires on exactly the expected steps (after steps `N-1, 2N-1, ...`, resp.
!!      `N-2, 2N-2, ...`), the count after every other step is the number of dynamic steps since the
!!      last assembly, every mean is over `N` samples (`N-1` in the first deferred interval), and
!!      `Q_heat`, `fw_net`, `Q_salt`, `sw_net` are those means with the sums empty afterwards and
!!      `sw_net` held in between.  Two further cases drive a REAL per-step filler (time-varying
!!      `&ocean_dataovr_nml` heat and lprec ramps, nothing hand-set), so an accumulate placed before
!!      the fillers (a one-step shift of a strictly increasing ramp) fails.
!!   5. `ifm_restart_registered`: the sums, the count and `sw_net` are restart state; ratio-1 restart
!!      files do not carry them.
!!   6. `ifm_midinterval_restart`: ratio 4, checkpoint after step 5 (two samples in the sums), resume:
!!      bitwise equal to the straight run at the next two assemblies; a checkpoint without
!!      `sf_sw_net` resumes with `sw_net` formed once from the first resumed step.
module test_ocean_surface_flux_interval_mean
   use, intrinsic :: iso_fortran_env, only: output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_assemble, &
                                     ocean_surface_flux_fw_virtual_step, &
                                     ocean_surface_flux_refresh_sw_net, &
                                     ocean_surface_flux_accumulate, &
                                     ocean_surface_flux_apply_interval_mean, &
                                     ocean_surface_flux_clear_accum, &
                                     SURF_FW_NONE, SURF_FW_VIRTUAL, SURF_FW_MASS, SF_NACC
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_restart_write_drop_field
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, nc_enddef, &
                            rdb_def_var_1d, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_surface_flux_interval_mean_tests

   integer, parameter :: NZ = 4
   integer, parameter :: NXP = 6
   integer, parameter :: NYP = 4
   integer, parameter :: NGHOST = 3
   integer, parameter :: NSTEP = 4
   real(wp), parameter :: DT = 300.0_wp
   character(len=*), parameter :: HEAT_FILE = "test_ifm_heat.nc"
   character(len=*), parameter :: LPREC_FILE = "test_ifm_lprec.nc"
   character(len=*), parameter :: NL = new_line("a")
   logical, save :: comm_inited = .false.

contains

   subroutine collect_ocean_surface_flux_interval_mean_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ifm_assembly_of_means", test_assembly_of_means), &
                  new_unittest("ifm_energy_closure", test_energy_closure), &
                  new_unittest("ifm_accumulators_empty_and_ratio1_inert", test_empty_and_inert), &
                  new_unittest("ifm_engine_wiring", test_engine_wiring), &
                  new_unittest("ifm_restart_registered", test_restart_registered), &
                  new_unittest("ifm_midinterval_restart", test_midinterval_restart) &
                  ]
   end subroutine collect_ocean_surface_flux_interval_mean_tests

   ! -----------------------------------------------------------------
   ! unit-level scaffolding
   ! -----------------------------------------------------------------

   subroutine setup_slot(grid, ms, sf, mode, accum)
      !! A components-on slot over a 6x4 tile, stratified T and S, `sw_source = "q_sw"`.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: mode
      logical, intent(in) :: accum
      integer :: k

      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call sf%set_interval_mean(grid, accum)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
      call sf%set_sw_penetration(0.6_wp, 1.0_wp, 0.35_wp, 23.0_wp, sw_source="q_sw")
      call sf%set_freshwater(mode)
      ms%h_layer = 10.0_wp
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%c(:, :, k) = 8.0_wp + 0.5_wp*real(k, wp)
         ms%tracers(ms%idx_salinity)%c(:, :, k) = 35.0_wp - 0.1_wp*real(k, wp)
      end do
      sf%has_heat = .true.
      sf%sw_ice = 7.0_wp
   end subroutine setup_slot

   subroutine map_all(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
   end subroutine map_all

   subroutine unmap_all(ms, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine unmap_all

   pure function samp(step, i, j, base, slope) result(v)
      !! A component sample: differs on every step and in every cell.
      integer, intent(in) :: step, i, j
      real(wp), intent(in) :: base, slope
      real(wp) :: v
      v = base*(1.0_wp + 0.35_wp*real(step, wp)*real(step, wp) - 0.2_wp*real(step, wp)) &
          + slope*real(i + 2*j, wp)*(1.0_wp + 0.1_wp*real(step, wp))
   end function samp

   subroutine set_sample(sf, step)
      !! Fill every accumulated component (host) with step `step`'s sample.
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: step
      integer :: i, j

      do j = 1, size(sf%q_sw, 2)
         do i = 1, size(sf%q_sw, 1)
            sf%q_sw(i, j) = samp(step, i, j, 150.0_wp, 3.0_wp)
            sf%q_lw(i, j) = samp(step, i, j, -45.0_wp, -0.7_wp)
            sf%q_lat(i, j) = samp(step, i, j, -60.0_wp, -1.1_wp)
            sf%q_sens(i, j) = samp(step, i, j, -12.0_wp, -0.4_wp)
            sf%q_fus(i, j) = samp(step, i, j, -3.0_wp, -0.05_wp)
            sf%q_fus_rof(i, j) = samp(step, i, j, -0.5_wp, -0.01_wp)
            sf%heat_added(i, j) = samp(step, i, j, 2.0_wp, 0.3_wp)
            sf%evap(i, j) = -abs(samp(step, i, j, 2.5e-5_wp, 1.0e-6_wp))
            sf%lprec(i, j) = abs(samp(step, i, j, 1.5e-5_wp, 2.0e-6_wp))
            sf%fprec(i, j) = abs(samp(step, i, j, 4.0e-6_wp, 5.0e-7_wp))
            sf%lrunoff(i, j) = abs(samp(step, i, j, 2.0e-6_wp, 3.0e-7_wp))
            sf%frunoff(i, j) = abs(samp(step, i, j, 1.0e-6_wp, 1.0e-7_wp))
            sf%salt_flux(i, j) = samp(step, i, j, 2.0e-6_wp, 1.0e-7_wp)
            sf%heat_content_lprec(i, j) = samp(step, i, j, 5.0e-2_wp, 1.0e-3_wp)
            sf%heat_content_fprec(i, j) = samp(step, i, j, 1.0e-2_wp, 1.0e-3_wp)
         end do
      end do
   end subroutine set_sample

   subroutine push_components(sf)
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$acc update device(sf%q_sw, sf%q_lw, sf%q_lat, sf%q_sens, sf%q_fus, sf%q_fus_rof, sf%heat_added)
      !$acc update device(sf%evap, sf%lprec, sf%fprec, sf%lrunoff, sf%frunoff, sf%salt_flux)
      !$acc update device(sf%heat_content_lprec, sf%heat_content_fprec, sf%sw_ice)
   end subroutine push_components

   subroutine pull_assembled(sf)
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$acc update self(sf%Q_heat, sf%Q_salt, sf%fw_net, sf%Q_salt_tracer, sf%sw_net, sf%salt_fw)
      !$acc update self(sf%q_sw, sf%evap, sf%lprec, sf%accum)
   end subroutine pull_assembled

   subroutine assemble_all(grid, ms, sf)
      !! The assembling finalize's order, with the interval mean applied only when the slot has it.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_flux_t), intent(inout) :: sf

      if (sf%accum_on) then
         call ocean_surface_flux_apply_interval_mean(sf)
         call ocean_surface_flux_refresh_sw_net(sf, at_assembly=.true.)
      else
         call ocean_surface_flux_refresh_sw_net(sf)
      end if
      call ocean_surface_flux_fw_virtual_step(grid, sf, ms)
      call ocean_surface_flux_assemble(grid, sf, ms)
      if (sf%accum_on) call ocean_surface_flux_clear_accum(sf)
   end subroutine assemble_all

   function max_rel_diff(a, b) result(d)
      real(wp), intent(in) :: a(:, :), b(:, :)
      real(wp) :: d, scale

      scale = max(maxval(abs(b)), tiny(1.0_wp))
      d = maxval(abs(a - b))/scale
   end function max_rel_diff

   ! -----------------------------------------------------------------
   ! 1. assembly of the means
   ! -----------------------------------------------------------------
   subroutine assembly_case(error, mode, label)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: mode
      character(len=*), intent(in) :: label
      type(hgrid_t) :: grid, grid_r, grid_l
      type(multilayer_state_t) :: ms, ms_r, ms_l
      type(ocean_surface_flux_t) :: sf, sf_r, sf_l
      real(wp), allocatable :: qh(:, :), qs(:, :), fw(:, :), qst(:, :), swn(:, :)
      integer :: s, n

      call setup_slot(grid, ms, sf, mode, .true.)
      call setup_slot(grid_r, ms_r, sf_r, mode, .false.)
      call setup_slot(grid_l, ms_l, sf_l, mode, .false.)
      call map_all(ms, sf)
      call map_all(ms_r, sf_r)
      call map_all(ms_l, sf_l)

      ! the slot under test: four steps, each one's sample put in place and accumulated
      do s = 0, NSTEP - 1
         call set_sample(sf, s)
         call push_components(sf)
         call ocean_surface_flux_accumulate(sf)
      end do
      ! the last sample is what a sampling run would have assembled
      call set_sample(sf_l, NSTEP - 1)
      call push_components(sf_l)
      ! the reference: the arithmetic means, handed to a slot with no treatment
      call set_mean_sample(sf_r)
      call push_components(sf_r)

      call assemble_all(grid, ms, sf)
      call assemble_all(grid_r, ms_r, sf_r)
      call assemble_all(grid_l, ms_l, sf_l)
      call pull_assembled(sf)
      call pull_assembled(sf_r)
      call pull_assembled(sf_l)

      qh = sf%Q_heat
      qs = sf%Q_salt
      fw = sf%fw_net
      qst = sf%Q_salt_tracer
      swn = sf%sw_net
      call check(error, max_rel_diff(qh, sf_r%Q_heat) < 1.0e-14_wp, label//": Q_heat is the assembly of the means")
      if (allocated(error)) return
      call check(error, max_rel_diff(qs, sf_r%Q_salt) < 1.0e-14_wp, label//": Q_salt is the assembly of the means")
      if (allocated(error)) return
      call check(error, max_rel_diff(swn, sf_r%sw_net) < 1.0e-14_wp, label//": sw_net is the mean q_sw + sw_ice")
      if (allocated(error)) return
      if (mode == SURF_FW_MASS) then
         call check(error, max_rel_diff(fw, sf_r%fw_net) < 1.0e-14_wp, label//": fw_net is the assembly of the means")
         if (allocated(error)) return
         call check(error, max_rel_diff(qst, sf_r%Q_salt_tracer) < 1.0e-14_wp, &
                    label//": Q_salt_tracer is the assembly of the means")
         if (allocated(error)) return
         call check(error, maxval(abs(fw)) > 0.0_wp, label//": non-vacuous: fw_net is non-zero")
         if (allocated(error)) return
      end if
      ! sw_net really is q_sw_mean + sw_ice, formed independently here
      s = NXP/2 + NGHOST
      n = NYP/2 + NGHOST
      call check(error, abs(swn(s, n) - (mean_of(150.0_wp, 3.0_wp, s, n) + 7.0_wp)) < 1.0e-12_wp*abs(swn(s, n)), &
                 label//": sw_net equals the hand-formed mean(q_sw) + sw_ice in a cell")
      if (allocated(error)) return
      ! non-vacuity: the sampled (last-step) assembly is far from the mean one
      call check(error, max_rel_diff(sf_l%Q_heat, qh) > 1.0e-2_wp, &
                 label//": non-vacuous: the last sample alone gives a different Q_heat")
      if (allocated(error)) return
      call check(error, max_rel_diff(sf_l%sw_net, swn) > 1.0e-2_wp, &
                 label//": non-vacuous: the last sample alone gives a different sw_net")
      if (allocated(error)) return

      call unmap_all(ms_l, sf_l)
      call unmap_all(ms_r, sf_r)
      call unmap_all(ms, sf)
      call sf%destroy()
      call ms%destroy()
      call sf_r%destroy()
      call ms_r%destroy()
      call sf_l%destroy()
      call ms_l%destroy()
   end subroutine assembly_case

   function mean_of(base, slope, i, j) result(m)
      !! Host mean over steps 0..NSTEP-1 of `samp` (same summation order as the accumulator).
      real(wp), intent(in) :: base, slope
      integer, intent(in) :: i, j
      real(wp) :: m, acc
      integer :: s

      acc = 0.0_wp
      do s = 0, NSTEP - 1
         acc = acc + samp(s, i, j, base, slope)
      end do
      m = acc/real(NSTEP, wp)
   end function mean_of

   subroutine set_mean_sample(sf)
      !! Host-form the arithmetic means of the samples of steps 0..NSTEP-1, same order as the sums.
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer :: i, j, s
      real(wp) :: a(SF_NACC)

      do j = 1, size(sf%q_sw, 2)
         do i = 1, size(sf%q_sw, 1)
            a = 0.0_wp
            do s = 0, NSTEP - 1
               a(1) = a(1) + samp(s, i, j, 150.0_wp, 3.0_wp)
               a(2) = a(2) + samp(s, i, j, -45.0_wp, -0.7_wp)
               a(3) = a(3) + samp(s, i, j, -60.0_wp, -1.1_wp)
               a(4) = a(4) + samp(s, i, j, -12.0_wp, -0.4_wp)
               a(5) = a(5) + samp(s, i, j, -3.0_wp, -0.05_wp)
               a(6) = a(6) + samp(s, i, j, -0.5_wp, -0.01_wp)
               a(7) = a(7) + samp(s, i, j, 2.0_wp, 0.3_wp)
               a(8) = a(8) - abs(samp(s, i, j, 2.5e-5_wp, 1.0e-6_wp))
               a(9) = a(9) + abs(samp(s, i, j, 1.5e-5_wp, 2.0e-6_wp))
               a(10) = a(10) + abs(samp(s, i, j, 4.0e-6_wp, 5.0e-7_wp))
               a(11) = a(11) + abs(samp(s, i, j, 2.0e-6_wp, 3.0e-7_wp))
               a(12) = a(12) + abs(samp(s, i, j, 1.0e-6_wp, 1.0e-7_wp))
               a(13) = a(13) + samp(s, i, j, 2.0e-6_wp, 1.0e-7_wp)
               a(14) = a(14) + samp(s, i, j, 5.0e-2_wp, 1.0e-3_wp)
               a(15) = a(15) + samp(s, i, j, 1.0e-2_wp, 1.0e-3_wp)
            end do
            a = a/real(NSTEP, wp)
            sf%q_sw(i, j) = a(1)
            sf%q_lw(i, j) = a(2)
            sf%q_lat(i, j) = a(3)
            sf%q_sens(i, j) = a(4)
            sf%q_fus(i, j) = a(5)
            sf%q_fus_rof(i, j) = a(6)
            sf%heat_added(i, j) = a(7)
            sf%evap(i, j) = a(8)
            sf%lprec(i, j) = a(9)
            sf%fprec(i, j) = a(10)
            sf%lrunoff(i, j) = a(11)
            sf%frunoff(i, j) = a(12)
            sf%salt_flux(i, j) = a(13)
            sf%heat_content_lprec(i, j) = a(14)
            sf%heat_content_fprec(i, j) = a(15)
         end do
      end do
   end subroutine set_mean_sample

   subroutine test_assembly_of_means(error)
      type(error_type), allocatable, intent(out) :: error

      call assembly_case(error, SURF_FW_VIRTUAL, "virtual")
      if (allocated(error)) return
      call assembly_case(error, SURF_FW_MASS, "mass")
   end subroutine test_assembly_of_means

   ! -----------------------------------------------------------------
   ! 2. energy closure
   ! -----------------------------------------------------------------
   subroutine test_energy_closure(error)
      !! `freshwater = "none"`: `Q_heat` carries the six mass-flux enthalpies and `cp*SST*evap` as well,
      !! so the whole open-water heat budget is in it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      real(wp) :: expect, got, worst, sst
      integer :: s, i, j

      call setup_slot(grid, ms, sf, SURF_FW_NONE, .true.)
      sf%sw_ice = 0.0_wp
      call map_all(ms, sf)
      call push_components(sf)
      do s = 0, NSTEP - 1
         call set_sample(sf, s)
         call push_components(sf)
         call ocean_surface_flux_accumulate(sf)
      end do
      call assemble_all(grid, ms, sf)
      call pull_assembled(sf)
      !$acc update self(sf%Q_heat)

      worst = 0.0_wp
      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            sst = ms%tracers(ms%idx_temperature)%c(i, j, NZ)
            expect = 0.0_wp
            do s = 0, NSTEP - 1
               expect = expect + DT*( &
                        samp(s, i, j, 150.0_wp, 3.0_wp) + samp(s, i, j, -45.0_wp, -0.7_wp) &
                        + samp(s, i, j, -60.0_wp, -1.1_wp) + samp(s, i, j, -12.0_wp, -0.4_wp) &
                        + samp(s, i, j, -3.0_wp, -0.05_wp) + samp(s, i, j, -0.5_wp, -0.01_wp) &
                        + samp(s, i, j, 2.0_wp, 0.3_wp) &
                        + samp(s, i, j, 5.0e-2_wp, 1.0e-3_wp) + samp(s, i, j, 1.0e-2_wp, 1.0e-3_wp) &
                        + sf%cp*sst*(-abs(samp(s, i, j, 2.5e-5_wp, 1.0e-6_wp))))
            end do
            got = sf%Q_heat(i, j)*real(NSTEP, wp)*DT
            worst = max(worst, abs(got - expect)/max(abs(expect), 1.0_wp))
         end do
      end do
      call check(error, worst < 1.0e-12_wp, &
                 "Q_heat * N*dt equals dt * sum over the steps of the per-step open-water heat (to round-off)")
      call unmap_all(ms, sf)
      call sf%destroy()
      call ms%destroy()
   end subroutine test_energy_closure

   ! -----------------------------------------------------------------
   ! 3. empty after assembly, ratio 1 inert
   ! -----------------------------------------------------------------
   subroutine test_empty_and_inert(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid, grid1
      type(multilayer_state_t) :: ms, ms1
      type(ocean_surface_flux_t) :: sf, sf1
      real(wp), allocatable :: q_sw0(:, :), evap0(:, :), lprec0(:, :)
      integer :: s

      call setup_slot(grid, ms, sf, SURF_FW_VIRTUAL, .true.)
      call map_all(ms, sf)
      do s = 0, NSTEP - 1
         call set_sample(sf, s)
         call push_components(sf)
         call ocean_surface_flux_accumulate(sf)
      end do
      !$acc update self(sf%accum)
      call check(error, sf%accum_n == real(NSTEP, wp), "the count is the number of accumulated steps")
      if (allocated(error)) return
      call check(error, maxval(abs(sf%accum)) > 0.0_wp, "non-vacuous: the sums are non-zero before the assembly")
      if (allocated(error)) return
      call assemble_all(grid, ms, sf)
      !$acc update self(sf%accum)
      call check(error, sf%accum_n == 0.0_wp, "the count is zero after the assembly")
      if (allocated(error)) return
      call check(error, all(sf%accum == 0.0_wp), "the sums are zero after the assembly")
      if (allocated(error)) return
      call unmap_all(ms, sf)
      call sf%destroy()
      call ms%destroy()

      ! ratio 1: no treatment, nothing allocated, every call a bitwise no-op
      call setup_slot(grid1, ms1, sf1, SURF_FW_VIRTUAL, .false.)
      call check(error,.not. sf1%accum_on .and. .not. allocated(sf1%accum), &
                 "ratio 1: the accumulator is not allocated")
      if (allocated(error)) return
      call set_sample(sf1, 2)
      allocate (q_sw0, source=sf1%q_sw)
      allocate (evap0, source=sf1%evap)
      allocate (lprec0, source=sf1%lprec)
      call ocean_surface_flux_accumulate(sf1)
      call ocean_surface_flux_apply_interval_mean(sf1)
      call ocean_surface_flux_clear_accum(sf1)
      call check(error, all(sf1%q_sw == q_sw0) .and. all(sf1%evap == evap0) .and. all(sf1%lprec == lprec0), &
                 "ratio 1: the components are untouched, bit for bit")
      if (allocated(error)) return
      call sf1%set_interval_mean(grid1, .true.)
      call check(error, allocated(sf1%accum), "set_interval_mean(.true.) allocates")
      if (allocated(error)) return
      call sf1%set_interval_mean(grid1, .false.)
      call check(error,.not. allocated(sf1%accum) .and. .not. sf1%accum_on, "set_interval_mean(.false.) frees")
      call sf1%destroy()
      call ms1%destroy()
   end subroutine test_empty_and_inert

   ! -----------------------------------------------------------------
   ! engine
   ! -----------------------------------------------------------------
   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function engine_nml(ratio, deferred, filler) result(nml)
      !! A closed 12 x 8 x 4 basin, components on, real fresh-water mass, `sw_source = "q_sw"`.
      !! `filler`: a time-varying data-forcing file for `heat` (-> `heat_added`) and `lprec`, so a real
      !! per-step filler writes components inside `engine_step`.
      integer, intent(in) :: ratio
      logical, intent(in) :: deferred
      logical, intent(in), optional :: filler
      character(len=:), allocatable :: nml, vmix, ovr
      character(len=16) :: r_s

      vmix = ""
      if (ratio > 1) then
         write (r_s, '(I0)') ratio
         vmix = "&ocean_vmix_nml dt_therm_ratio = "//trim(adjustl(r_s))
         if (deferred) vmix = vmix//", diabatic_first = .false., dt_tracer_advect_ratio = "//trim(adjustl(r_s))
         vmix = vmix//" /"//NL
      end if
      ovr = ""
      if (present(filler)) then
         if (filler) ovr = "&ocean_dataovr_nml enable = .true., time_mode = 'linear', heat_file = '"//HEAT_FILE// &
                           "', heat_var = 'field', lprec_file = '"//LPREC_FILE//"', lprec_var = 'field' /"//NL
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 12, ny = 8, nghost = 3, dx = 4000.0, dy = 4000.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 400.0, wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_bottom = 6.0, T_init_surface = 12.0 /"//NL// &
            "&ocean_eos_nml eos = 'linear' /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .true. /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true., sw_pen_frac = 0.6, sw_source = 'q_sw' /"//NL// &
            "&ocean_forcing_nml enable_components = .true., freshwater = 'mass' /"//NL// &
            "&ocean_hvisc_nml nu_h = 100.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL// &
            "&vcoord_nml vcoord_type = 'zstar' /"//NL//ovr//vmix
   end function engine_nml

   subroutine start_engine(nml, engine, cfg, ok)
      character(len=*), intent(in) :: nml
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      integer :: ierr

      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= 0) return
      call validate_config(cfg, ierr)
      if (ierr /= 0) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine start_engine

   subroutine start_engine_restart(nml, engine, cfg, ok, restart_file, t0, step0)
      character(len=*), intent(in) :: nml, restart_file
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      real(wp), intent(out) :: t0
      integer, intent(out) :: step0
      integer :: ierr

      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= 0) return
      call validate_config(cfg, ierr)
      if (ierr /= 0) return
      call engine_setup(engine, cfg, ierr, restart_file=restart_file, t_restart=t0, step_restart=step0)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine start_engine_restart

   subroutine write_ramp_file(filename, amp0, amp1)
      !! `field(x, y, t)`, uniform in space, linear in time from `amp0` at t = 0 to `amp1` at t = 1e6 s.
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: amp0, amp1
      integer, parameter :: NX = 12, NY = 8
      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr
      real(wp) :: f(NX, NY, 2)

      f(:, :, 1) = amp0
      f(:, :, 2) = amp1
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", NX, dim_x)
      call nc_def_dim(ncid, "y", NY, dim_y)
      call nc_def_dim(ncid, "time", 2, dim_t)
      call nc_def_var_3d(ncid, "field", [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, [0.0_wp, 1.0e6_wp])
      call nc_close(ncid)
   end subroutine write_ramp_file

   pure function heat_at(t) result(q)
      !! The data-forcing heat at model time `t` (linear ramp, 100 -> 1100 W m-2 over 1e6 s).
      real(wp), intent(in) :: t
      real(wp) :: q
      q = 100.0_wp + 1.0e-3_wp*t
   end function heat_at

   pure function lprec_at(t) result(p)
      !! The data-forcing liquid precipitation at `t` (1e-5 -> 3e-5 kg m-2 s-1 over 1e6 s).
      real(wp), intent(in) :: t
      real(wp) :: p
      p = 1.0e-5_wp + 2.0e-11_wp*t
   end function lprec_at

   pure function assembles_after(n, ratio, deferred) result(yes)
      !! Does the finalize of step index `n` (0-based) assemble?  The assembly fires when the
      !! post-increment counter `n + 1` is a thermo step: `diabatic_first` thermo steps are
      !! `mod(c, N) == 0`, the deferred ones `mod(c + 1, N) == 0`, so the assemblies follow steps
      !! `N-1, 2N-1, ...` and `N-2, 2N-2, ...` respectively.
      integer, intent(in) :: n, ratio
      logical, intent(in) :: deferred
      logical :: yes

      if (deferred) then
         yes = (mod(n + 2, ratio) == 0)
      else
         yes = (mod(n + 1, ratio) == 0)
      end if
   end function assembles_after

   subroutine pull_ms(ms)
      type(multilayer_state_t), intent(inout) :: ms

      !$acc update self(ms%wet_mask)
      !$acc update self(ms%tracers(ms%idx_salinity)%c)
   end subroutine pull_ms

   subroutine engine_case(error, ratio, deferred, label)
      !! Hand-set components per step: the means, the NUMBER of samples (N; N-1 in the first
      !! deferred interval), the step indices the assembly fires on, and the count after every step.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ratio
      logical, intent(in) :: deferred
      character(len=*), intent(in) :: label
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok, fires
      integer, parameter :: NSTEPS_RUN = 14
      integer :: n, ierr, i, j, s0, nas, ng, nxp, nyp, k, n_expect
      real(wp) :: t, hq, fq, sq, wq, worst_h, worst_f, worst_s, worst_w
      real(wp), allocatable :: sw_seen(:, :)

      call start_engine(engine_nml(ratio, deferred), engine, cfg, ok)
      call check(error, ok, label//": the engine starts")
      if (allocated(error)) return
      gate: block
         call check(error, engine%state%surface_flux%accum_on, label//": the interval mean is on at ratio > 1")
         if (allocated(error)) exit gate
         engine%state%surface_flux%has_mass_flux = .true.
         ng = engine%grid%nghost
         nxp = engine%grid%nx_phys
         nyp = engine%grid%ny_phys
         k = engine%state%multilayer%nz_ml
         s0 = 0
         nas = 0
         t = 0.0_wp
         do n = 0, NSTEPS_RUN - 1
            call set_sample(engine%state%surface_flux, n)
            call push_components(engine%state%surface_flux)
            call engine_step(engine, 300.0_wp, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": engine_step")
            if (allocated(error)) exit gate
            call engine_step_ice(engine, cfg, 300.0_wp, t, ierr=ierr)
            call engine_step_finalize(engine, 300.0_wp, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": engine_step_finalize")
            if (allocated(error)) exit gate
            t = t + 300.0_wp
            fires = assembles_after(n, ratio, deferred)
            ! the assembly fires on exactly the expected steps, and after every other step the count is
            ! the number of dynamic steps since the last one (a double accumulate, a skipped one or an
            ! assembly on the wrong step shows here)
            call check(error, (engine%state%surface_flux%accum_n == 0.0_wp) .eqv. fires, &
                       label//": the assembly fires on exactly the expected steps")
            if (allocated(error)) exit gate
            if (.not. fires) then
               call check(error, engine%state%surface_flux%accum_n == real(n + 1 - s0, wp), &
                          label//": the count is one per dynamic step since the last assembly")
               if (allocated(error)) exit gate
            end if
            if (fires) then
               nas = nas + 1
               n_expect = ratio
               if (deferred .and. nas == 1) n_expect = ratio - 1
               call check(error, n - s0 + 1 == n_expect, &
                          label//": the mean is over N samples (N-1 in the first deferred interval)")
               if (allocated(error)) exit gate
               call pull_assembled(engine%state%surface_flux)
               call pull_ms(engine%state%multilayer)
               worst_h = 0.0_wp
               worst_f = 0.0_wp
               worst_s = 0.0_wp
               worst_w = 0.0_wp
               do j = ng + 1, ng + nyp
                  do i = ng + 1, ng + nxp
                     hq = mean_range(s0, n, i, j, 150.0_wp, 3.0_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, -60.0_wp, -1.1_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, -45.0_wp, -0.7_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, -12.0_wp, -0.4_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, -3.0_wp, -0.05_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, -0.5_wp, -0.01_wp, 1.0_wp) &
                          + mean_range(s0, n, i, j, 2.0_wp, 0.3_wp, 1.0_wp)
                     worst_h = max(worst_h, abs(engine%state%surface_flux%Q_heat(i, j) - hq)/max(abs(hq), 1.0_wp))
                     fq = mean_range(s0, n, i, j, 1.5e-5_wp, 2.0e-6_wp, 2.0_wp) &
                          + mean_range(s0, n, i, j, 4.0e-6_wp, 5.0e-7_wp, 2.0_wp) &
                          + mean_range(s0, n, i, j, 2.5e-5_wp, 1.0e-6_wp, -2.0_wp) &
                          + mean_range(s0, n, i, j, 2.0e-6_wp, 3.0e-7_wp, 2.0_wp) &
                          + mean_range(s0, n, i, j, 1.0e-6_wp, 1.0e-7_wp, 2.0_wp)
                     worst_f = max(worst_f, abs(engine%state%surface_flux%fw_net(i, j) &
                                                - engine%state%multilayer%wet_mask(i, j)*fq)/1.0e-5_wp)
                     sq = -engine%state%multilayer%wet_mask(i, j)*engine%state%multilayer%tracers( &
                          engine%state%multilayer%idx_salinity)%c(i, j, k)*fq &
                          + mean_range(s0, n, i, j, 2.0e-6_wp, 1.0e-7_wp, 1.0_wp)
                     worst_s = max(worst_s, abs(engine%state%surface_flux%Q_salt(i, j) - sq)/max(abs(sq), 1.0e-6_wp))
                     wq = mean_range(s0, n, i, j, 150.0_wp, 3.0_wp, 1.0_wp)
                     worst_w = max(worst_w, abs(engine%state%surface_flux%sw_net(i, j) - wq)/max(abs(wq), 1.0_wp))
                  end do
               end do
               call check(error, worst_h < 1.0e-12_wp, label//": Q_heat is the mean of the interval's samples")
               if (allocated(error)) exit gate
               call check(error, worst_f < 1.0e-9_wp, label//": fw_net is the mean of the interval's samples")
               if (allocated(error)) exit gate
               call check(error, worst_s < 1.0e-9_wp, label//": Q_salt is built from the mean fresh water")
               if (allocated(error)) exit gate
               call check(error, worst_w < 1.0e-12_wp, label//": sw_net is the mean q_sw (+ sw_ice = 0)")
               if (allocated(error)) exit gate
               call check(error, all(engine%state%surface_flux%accum == 0.0_wp), &
                          label//": the sums are empty after the assembly")
               if (allocated(error)) exit gate
               sw_seen = engine%state%surface_flux%sw_net
               s0 = n + 1
            else if (nas > 0) then
               call pull_assembled(engine%state%surface_flux)
               call check(error, all(engine%state%surface_flux%sw_net == sw_seen), &
                          label//": sw_net is held between assemblies")
               if (allocated(error)) exit gate
            end if
         end do
         call check(error, nas >= 3, label//": at least three assemblies were observed")
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine engine_case

   subroutine filler_case(error, ratio, deferred, label)
      !! A REAL per-step filler (`&ocean_dataovr_nml` heat and lprec files, time-varying) writes the
      !! components inside `engine_step`; nothing is hand-set.  Both ramps are strictly increasing, so a
      !! sample taken BEFORE the filler ran (the previous step's value) shifts every mean by a step and
      !! fails the tight tolerance.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ratio
      logical, intent(in) :: deferred
      character(len=*), intent(in) :: label
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer, parameter :: NSTEPS_RUN = 14
      integer :: n, ierr, s0, nas, ng, i, j, k, m
      real(wp) :: t, hm, pm, worst_h, worst_f, worst_s

      call write_ramp_file(HEAT_FILE, 100.0_wp, 1100.0_wp)
      call write_ramp_file(LPREC_FILE, 1.0e-5_wp, 3.0e-5_wp)
      call start_engine(engine_nml(ratio, deferred, filler=.true.), engine, cfg, ok)
      call check(error, ok, label//": the engine starts")
      if (allocated(error)) return
      gate: block
         ng = engine%grid%nghost
         k = engine%state%multilayer%nz_ml
         s0 = 0
         nas = 0
         t = 0.0_wp
         do n = 0, NSTEPS_RUN - 1
            call engine_step(engine, 300.0_wp, t, ierr=ierr)
            call engine_step_ice(engine, cfg, 300.0_wp, t, ierr=ierr)
            call engine_step_finalize(engine, 300.0_wp, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": the step succeeds")
            if (allocated(error)) exit gate
            t = t + 300.0_wp
            if (assembles_after(n, ratio, deferred)) then
               nas = nas + 1
               call pull_assembled(engine%state%surface_flux)
               call pull_ms(engine%state%multilayer)
               ! the filler of step m writes the value at the START of the step, t_m = 300*m
               hm = 0.0_wp
               pm = 0.0_wp
               do m = s0, n
                  hm = hm + heat_at(300.0_wp*real(m, wp))
                  pm = pm + lprec_at(300.0_wp*real(m, wp))
               end do
               hm = hm/real(n - s0 + 1, wp)
               pm = pm/real(n - s0 + 1, wp)
               worst_h = 0.0_wp
               worst_f = 0.0_wp
               worst_s = 0.0_wp
               do j = ng + 1, ng + engine%grid%ny_phys
                  do i = ng + 1, ng + engine%grid%nx_phys
                     worst_h = max(worst_h, abs(engine%state%surface_flux%Q_heat(i, j) - hm)/hm)
                     worst_f = max(worst_f, abs(engine%state%surface_flux%fw_net(i, j) &
                                                - engine%state%multilayer%wet_mask(i, j)*pm)/pm)
                     worst_s = max(worst_s, abs(engine%state%surface_flux%Q_salt(i, j) &
                                                + engine%state%multilayer%wet_mask(i, j)*pm &
                                                *engine%state%multilayer%tracers( &
                                                engine%state%multilayer%idx_salinity)%c(i, j, k))/(pm*30.0_wp))
                  end do
               end do
               call check(error, worst_h < 1.0e-11_wp, &
                          label//": Q_heat is the mean of the filler's per-step values (accumulate runs after the filler)")
               if (allocated(error)) exit gate
               call check(error, worst_f < 1.0e-11_wp, label//": fw_net is the mean of the filler's lprec")
               if (allocated(error)) exit gate
               call check(error, worst_s < 1.0e-11_wp, label//": Q_salt is the virtual salt of the mean lprec")
               if (allocated(error)) exit gate
               s0 = n + 1
            end if
         end do
         call check(error, nas >= 3, label//": at least three assemblies were observed")
      end block gate
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine filler_case

   pure function mean_range(s_lo, s_hi, i, j, base, slope, sgn) result(m)
      !! Mean over steps `s_lo..s_hi` of `samp` (sgn = 1), `|samp|` (sgn = 2) or `-|samp|` (sgn = -2).
      integer, intent(in) :: s_lo, s_hi, i, j
      real(wp), intent(in) :: base, slope, sgn
      real(wp) :: m, acc, v
      integer :: s

      acc = 0.0_wp
      do s = s_lo, s_hi
         v = samp(s, i, j, base, slope)
         if (sgn > 1.5_wp) v = abs(v)
         if (sgn < -1.5_wp) v = -abs(v)
         acc = acc + v
      end do
      m = acc/real(s_hi - s_lo + 1, wp)
   end function mean_range

   subroutine test_engine_wiring(error)
      type(error_type), allocatable, intent(out) :: error

      call engine_case(error, 4, .false., "ratio 4")
      if (allocated(error)) return
      call engine_case(error, 4, .true., "ratio 4 deferred")
      if (allocated(error)) return
      call filler_case(error, 4, .false., "filler, ratio 4")
      if (allocated(error)) return
      call filler_case(error, 4, .true., "filler, ratio 4 deferred")
   end subroutine test_engine_wiring

   ! -----------------------------------------------------------------
   ! mid-interval restart
   ! -----------------------------------------------------------------
   subroutine set_qsw(sf, step)
      !! Only the shortwave is hand-set (the fillers write heat and lprec), so `sw_net` is non-zero.
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: step
      integer :: i, j

      do j = 1, size(sf%q_sw, 2)
         do i = 1, size(sf%q_sw, 1)
            sf%q_sw(i, j) = samp(step, i, j, 150.0_wp, 3.0_wp)
         end do
      end do
      !$acc update device(sf%q_sw)
   end subroutine set_qsw

   subroutine run_steps(engine, cfg, t, n_lo, n_hi, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: n_lo, n_hi
      logical, intent(out) :: ok
      integer :: n, ierr

      ok = .false.
      do n = n_lo, n_hi
         call set_qsw(engine%state%surface_flux, n)
         call engine_step(engine, 300.0_wp, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_ice(engine, cfg, 300.0_wp, t, ierr=ierr)
         call engine_step_finalize(engine, 300.0_wp, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + 300.0_wp
      end do
      ok = .true.
   end subroutine run_steps

   subroutine test_midinterval_restart(error)
      !! Ratio 4, checkpoint after step 5 (mid-interval: two samples in the sums).  A resumed run is
      !! bitwise the straight run at the next two assemblies (`Q_heat`, `Q_salt`, `fw_net`,
      !! `sw_net`); a checkpoint WITHOUT `sf_sw_net` resumes with `sw_net` formed once from the first
      !! resumed step's `q_sw + sw_ice`, never zero against the restored non-zero `Q_heat`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: FN = "test_ifm_restart.nc"
      type(ocean_engine_t) :: ea, eb, ec
      type(config_t) :: cfg_a, cfg_b, cfg_c
      real(wp) :: t_a, t_b, t_c
      real(wp), allocatable :: qh7(:, :), qs7(:, :), fw7(:, :), sw7(:, :)
      real(wp), allocatable :: qh11(:, :), qs11(:, :), fw11(:, :), sw11(:, :)
      integer :: step_b, step_c, ierr
      logical :: ok

      call write_ramp_file(HEAT_FILE, 100.0_wp, 1100.0_wp)
      call write_ramp_file(LPREC_FILE, 1.0e-5_wp, 3.0e-5_wp)
      ! straight run to step 11, recording the assemblies after steps 7 and 11
      call start_engine(engine_nml(4, .false., filler=.true.), ea, cfg_a, ok)
      call check(error, ok, "straight engine starts")
      if (allocated(error)) return
      t_a = 0.0_wp
      call run_steps(ea, cfg_a, t_a, 0, 7, ok)
      call check(error, ok, "straight run to step 7")
      if (allocated(error)) return
      call pull_assembled(ea%state%surface_flux)
      qh7 = ea%state%surface_flux%Q_heat
      qs7 = ea%state%surface_flux%Q_salt
      fw7 = ea%state%surface_flux%fw_net
      sw7 = ea%state%surface_flux%sw_net
      call run_steps(ea, cfg_a, t_a, 8, 11, ok)
      call check(error, ok, "straight run to step 11")
      if (allocated(error)) return
      call pull_assembled(ea%state%surface_flux)
      qh11 = ea%state%surface_flux%Q_heat
      qs11 = ea%state%surface_flux%Q_salt
      fw11 = ea%state%surface_flux%fw_net
      sw11 = ea%state%surface_flux%sw_net
      call check(error, maxval(abs(sw7)) > 0.0_wp .and. maxval(abs(qh7)) > 0.0_wp, &
                 "non-vacuous: sw_net and Q_heat are non-zero")
      if (allocated(error)) return
      call engine_exit_data(ea)
      call engine_teardown(ea)

      ! a second run to the checkpoint (after step 5), written in full
      call start_engine(engine_nml(4, .false., filler=.true.), ea, cfg_a, ok)
      t_a = 0.0_wp
      call run_steps(ea, cfg_a, t_a, 0, 5, ok)
      call check(error, ok .and. ea%state%surface_flux%accum_n == 2.0_wp, &
                 "two samples are in the sums at the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, 6, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write")
      if (allocated(error)) return
      call engine_exit_data(ea)
      call engine_teardown(ea)

      call start_engine_restart(engine_nml(4, .false., filler=.true.), eb, cfg_b, ok, FN, t_b, step_b)
      call check(error, ok .and. step_b == 6, "resumed engine starts")
      if (allocated(error)) return
      call check(error, eb%state%surface_flux%accum_n == 2.0_wp .and. .not. eb%state%surface_flux%sw_net_stale, &
                 "the resumed run carries the two samples and a valid sw_net")
      if (allocated(error)) return
      call run_steps(eb, cfg_b, t_b, 6, 7, ok)
      call check(error, ok, "resumed run to step 7")
      if (allocated(error)) return
      call pull_assembled(eb%state%surface_flux)
      call check(error, all(eb%state%surface_flux%Q_heat == qh7) .and. all(eb%state%surface_flux%Q_salt == qs7) &
                 .and. all(eb%state%surface_flux%fw_net == fw7) .and. all(eb%state%surface_flux%sw_net == sw7), &
                 "resumed == straight, bitwise, at the next assembly (Q_heat, Q_salt, fw_net, sw_net)")
      if (allocated(error)) return
      call run_steps(eb, cfg_b, t_b, 8, 11, ok)
      call pull_assembled(eb%state%surface_flux)
      call check(error, ok .and. all(eb%state%surface_flux%Q_heat == qh11) &
                 .and. all(eb%state%surface_flux%Q_salt == qs11) &
                 .and. all(eb%state%surface_flux%fw_net == fw11) .and. all(eb%state%surface_flux%sw_net == sw11), &
                 "resumed == straight, bitwise, at the following assembly")
      if (allocated(error)) return
      call engine_exit_data(eb)
      call engine_teardown(eb)

      ! the same checkpoint with sf_sw_net dropped (a ratio-1 / pre-fix file)
      call start_engine(engine_nml(4, .false., filler=.true.), ea, cfg_a, ok)
      t_a = 0.0_wp
      call run_steps(ea, cfg_a, t_a, 0, 5, ok)
      call ocean_state_restart_write_drop_field(ea%state, ea%grid, ea%decomp, FN, t_a, 6, "sf_sw_net", ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write (sf_sw_net dropped)")
      if (allocated(error)) return
      call engine_exit_data(ea)
      call engine_teardown(ea)
      call start_engine_restart(engine_nml(4, .false., filler=.true.), ec, cfg_c, ok, FN, t_c, step_c)
      call check(error, ok, "resumed engine (no sf_sw_net) starts")
      if (allocated(error)) return
      call check(error, ec%state%surface_flux%sw_net_stale, "the missing sf_sw_net is latched")
      if (allocated(error)) return
      call run_steps(ec, cfg_c, t_c, 6, 6, ok)
      call pull_assembled(ec%state%surface_flux)
      call check(error, ok .and. .not. ec%state%surface_flux%sw_net_stale .and. &
                 abs(ec%state%surface_flux%sw_net(4, 4) - samp(6, 4, 4, 150.0_wp, 3.0_wp)) < 1.0e-12_wp, &
                 "sw_net is formed once from the first resumed step's q_sw + sw_ice, not left at zero")
      call engine_exit_data(ec)
      call engine_teardown(ec)
   end subroutine test_midinterval_restart

   ! -----------------------------------------------------------------
   ! restart registration
   ! -----------------------------------------------------------------
   subroutine test_restart_registered(error)
      !! The accumulators, the count and `sw_net` are registered at ratio > 1 (a checkpoint can fall inside
      !! an interval) and absent at ratio 1 (restart files unchanged).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok

      call start_engine(engine_nml(4, .false.), engine, cfg, ok)
      call check(error, ok, "ratio 4: the engine starts")
      if (allocated(error)) return
      call check(error, registry_has(engine, "sf_accum") .and. registry_has(engine, "sf_accum_n") &
                 .and. registry_has(engine, "sf_sw_net"), "ratio 4: sf_accum, sf_accum_n and sf_sw_net are registered")
      call engine_exit_data(engine)
      call engine_teardown(engine)
      if (allocated(error)) return

      call start_engine(engine_nml(1, .false.), engine, cfg, ok)
      call check(error, ok, "ratio 1: the engine starts")
      if (allocated(error)) return
      call check(error,.not. (registry_has(engine, "sf_accum") .or. registry_has(engine, "sf_accum_n") &
                              .or. registry_has(engine, "sf_sw_net")), "ratio 1: the restart file is unchanged")
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine test_restart_registered

   function registry_has(engine, tag) result(found)
      use rdb_ocean_restart, only: restart_registry_t
      use rdb_ocean_state, only: ocean_state_build_restart_registry
      type(ocean_engine_t), intent(inout) :: engine
      character(len=*), intent(in) :: tag
      logical :: found
      type(restart_registry_t) :: reg

      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      found = any(reg%entries(1:reg%n)%tag == tag)
   end function registry_has

end module test_ocean_surface_flux_interval_mean
