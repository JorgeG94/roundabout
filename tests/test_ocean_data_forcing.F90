!! Tests for file-backed surface forcing (`rdb_ocean_data_forcing`,
!! PR-15) — the `(file, variable) -> forcing slot` consumer of the PR-14
!! reader.
!!
!! Tests (NetCDF I/O — must run with OMP_NUM_THREADS=1):
!!   data_forcing_disabled_is_noop  — T0: `enable = .false.` registers
!!                                    nothing and `apply` leaves every
!!                                    slot at its configure-time value.
!!                                    The bit-identity guarantee.
!!   data_forcing_tau_half_amp      — T1: a 2-record tau file (0 at t=0,
!!                                    TAU at t=T) queried at T/2 puts
!!                                    exactly TAU/2 on every physical
!!                                    face — the half-amplitude wind that
!!                                    drives the half-amplitude Ekman
!!                                    response.  Device-resident: the
!!                                    whole session runs inside one
!!                                    mapped region and the result is
!!                                    read back with `update self`.
!!   data_forcing_ghosts_extended   — T2: ghost faces carry the
!!                                    zero-gradient extension of their
!!                                    nearest physical face, not the zero
!!                                    the reader leaves behind.
!!   data_forcing_stress_mag_fresh  — T3: `stress_mag` tracks the BLENDED
!!                                    wind across a bracket advance.  The
!!                                    regression guard for the trap that
!!                                    KPP/EPBL would otherwise mix on the
!!                                    configure-time wind forever.
!!   data_forcing_heat_to_q_heat    — T4: with the component set off, the
!!                                    heat tag lands in `Q_heat` and sets
!!                                    `has_heat`.
!!   data_forcing_cyclic_wrap       — T5: a climatology queried past the
!!                                    end of its period wraps to the
!!                                    (nt,1) seam bracket rather than
!!                                    aborting.
!!   data_forcing_components_route  — T6: with the component set ON, the
!!                                    heat/salt tags route to
!!                                    `heat_added`/`salt_flux` (NOT
!!                                    `Q_heat`/`Q_salt`, which the
!!                                    assembler would overwrite) and the
!!                                    freshwater tags reach `evap`/
!!                                    `lprec`.  Also asserts `Q_heat`
!!                                    stays untouched, which is the
!!                                    whole point of the routing rule.
!!   data_forcing_ekman_matches_const
!!                                  — T7: THE wiring test.  A
!!                                    time-constant forcing file driven
!!                                    through the full RK2 dyn step
!!                                    reproduces the equivalent
!!                                    `set_wind_stress_const` run to
!!                                    round-off, AND shows the physical
!!                                    Ekman signature (surface u grows,
!!                                    Coriolis turns it southward).  The
!!                                    earlier tests prove the file
!!                                    reaches `tau_x`; this one proves
!!                                    `tau_x` reaches the flow.
!!   data_forcing_ekman_half_amp    — T8: the half-amplitude claim,
!!                                    stated exactly.  One step taken at
!!                                    the midpoint of a 0->TAU bracket
!!                                    produces half the surface velocity
!!                                    of one step under the full TAU.
!!
!!   data_forcing_config_predicates — T9: both branches of every
!!                                    `&ocean_dataovr_nml` validity
!!                                    predicate.  These are what
!!                                    `validate_config` turns into an
!!                                    abort, so this is the coverage for
!!                                    the fail-loud paths — reached
!!                                    through the pure predicate rather
!!                                    than an in-process death test
!!                                    (`error stop` would take the test
!!                                    binary down with it; test-drive's
!!                                    `should_fail` does NOT help,
!!                                    it expects a RETURNED error).
!!
!! Why the predicates live in `rdb_config` and run in `validate_config`
!! rather than at registration: a broken config must not run at all.
!! `ocean_data_forcing_configure` executes deep in the driver setup,
!! after the grid, bathymetry and IC are built; a namelist typo caught
!! there has already cost the user all of that. Everything checkable
!! from `cfg` alone is therefore checked before any work happens. The
!! registration-time guards remain for the programmatic (non-namelist)
!! caller.
module test_ocean_data_forcing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: ocean_dataovr_config_t, dataovr_entry_config_t, &
                         dataovr_entry_is_valid, dataovr_time_is_valid, &
                         dataovr_freshwater_needs_components, dataovr_any_tag_set
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_ocean_data_forcing, only: ocean_data_forcing_t, &
                                     ocean_data_forcing_configure, &
                                     ocean_data_forcing_apply
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            rdb_def_var_1d, nc_enddef, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_data_forcing_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 4, NYP = 3
   real(wp), parameter :: TOL = 1.0e-12_wp

   ! --- dyn-step (Ekman) test parameters ---
   integer, parameter :: NZ = 4
   integer, parameter :: EK_NX = 16, EK_NY = 12
   real(wp), parameter :: EK_DX = 1.0e3_wp, EK_DY = 1.0e3_wp
   real(wp), parameter :: EK_H_LAYER = 25.0_wp
   real(wp), parameter :: EK_DT = 5.0_wp
   real(wp), parameter :: EK_F0 = 1.0e-4_wp
   real(wp), parameter :: EK_TAU = 0.1_wp
   integer, parameter :: EK_STEPS = 20

contains

   subroutine collect_ocean_data_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("data_forcing_disabled_is_noop", test_disabled_is_noop), &
                  new_unittest("data_forcing_tau_half_amp", test_tau_half_amp), &
                  new_unittest("data_forcing_ghosts_extended", test_ghosts_extended), &
                  new_unittest("data_forcing_stress_mag_fresh", test_stress_mag_fresh), &
                  new_unittest("data_forcing_heat_to_q_heat", test_heat_to_q_heat), &
                  new_unittest("data_forcing_cyclic_wrap", test_cyclic_wrap), &
                  new_unittest("data_forcing_components_route", test_components_route), &
                  new_unittest("data_forcing_ekman_matches_const", test_ekman_matches_const), &
                  new_unittest("data_forcing_ekman_half_amp", test_ekman_half_amp), &
                  new_unittest("data_forcing_config_predicates", test_config_predicates) &
                  ]
   end subroutine collect_ocean_data_forcing_tests

   ! =================================================================
   ! Test support
   ! =================================================================

   subroutine write_ramp_file(filename, nx, ny, nt, t_axis, amp)
      !! `field(x, y, time)` whose value is `amp(k)` everywhere at record
      !! `k` — spatially uniform on purpose, so an assertion failure
      !! reads as a time-interpolation or destination-offset bug rather
      !! than a horizontal indexing one.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nt
      real(wp), intent(in) :: t_axis(nt), amp(nt)

      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr, i, j, k
      real(wp) :: f(nx, ny, nt)

      do k = 1, nt
         do j = 1, ny
            do i = 1, nx
               f(i, j, k) = amp(k)
            end do
         end do
      end do

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_3d(ncid, "field", [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_ramp_file

   subroutine write_gradient_file(filename, nx, ny, nt, t_axis)
      !! `field(i, j, k) = 100*i + 10*j` for every record — a horizontally
      !! varying, time-constant field.  Distinguishes a correct
      !! `(dest_i0, dest_j0)` landing from an off-by-`nghost` one, and
      !! gives the ghost-extension test something with a gradient to
      !! extend.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny, nt
      real(wp), intent(in) :: t_axis(nt)

      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr, i, j, k
      real(wp) :: f(nx, ny, nt)

      do k = 1, nt
         do j = 1, ny
            do i = 1, nx
               f(i, j, k) = 100.0_wp*real(i, wp) + 10.0_wp*real(j, wp)
            end do
         end do
      end do

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", nt, dim_t)
      call nc_def_var_3d(ncid, "field", [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_gradient_file

   function make_grid() result(g)
      type(hgrid_t) :: g
      call g%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
   end function make_grid

   subroutine make_slots(grid, ss, sf)
      !! A minimal stress + flux pair, both zeroed, standing in for the
      !! configure-time seed.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(out) :: ss
      type(ocean_surface_flux_t), intent(out) :: sf
      call ss%init(grid)
      call sf%init(grid)
   end subroutine make_slots

   subroutine tau_x_config(cfg, file, var)
      !! `&ocean_dataovr_nml` with only the tau_x tag driven.
      type(ocean_dataovr_config_t), intent(out) :: cfg
      character(len=*), intent(in) :: file, var
      cfg%enable = .true.
      cfg%time_mode = "linear"
      cfg%tau_x%file = file
      cfg%tau_x%var = var
   end subroutine tau_x_config

   ! =================================================================
   ! T0 — disabled is a no-op (the bit-identity guarantee).
   ! =================================================================

   subroutine test_disabled_is_noop(error)
      type(error_type), allocatable, intent(out) :: error

      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      real(wp), parameter :: SEED = 0.077_wp

      grid = make_grid()
      call make_slots(grid, ss, sf)
      call reader%init()

      ! A file IS named — the master switch alone must suppress it.
      call tau_x_config(cfg, "/tmp/test_data_forcing_never_read.nc", "field")
      cfg%enable = .false.

      ss%tau_x = SEED
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      call check(error,.not. df%active, "T0: disabled config must not activate")
      if (allocated(error)) return
      call check(error, reader%nfields == 0, "T0: disabled config must register no fields")
      if (allocated(error)) return

      ! `apply` on an inactive slot must not touch the seeded wind (and
      ! must not trip the reader's freshness check either).
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 0.0_wp)
      call check(error, all(abs(ss%tau_x - SEED) < TOL), &
                 "T0: apply must leave the configure-time wind untouched")
      if (allocated(error)) return

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_disabled_is_noop

   ! =================================================================
   ! T1 — half-amplitude wind at the midpoint of a 2-record bracket.
   ! =================================================================

   subroutine test_tau_half_amp(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: TAU = 0.2_wp, TEND = 1000.0_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), amp(NT)
      real(wp) :: tq
      integer :: i, j

      t_axis = [0.0_wp, TEND]
      amp = [0.0_wp, TAU]
      fname = "/tmp/test_data_forcing_tau_half.nc"
      call write_ramp_file(trim(fname), NXP + 1, NYP, NT, t_axis, amp)

      grid = make_grid()
      call make_slots(grid, ss, sf)
      call reader%init()
      call tau_x_config(cfg, trim(fname), "field")
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      call check(error, df%active .and. df%id_tau_x > 0, "T1: tau_x must register")
      if (allocated(error)) return

      ! Full device-resident session: map the reader's brackets and the
      ! destination, blend on-device, pull the result back.  On a host
      ! build every directive is inert and the same assertions hold.
      tq = 0.5_wp*TEND
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ss%tau_x, ss%tau_y, ss%stress_mag)
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, tq)
      !$acc update self(ss%tau_x, ss%stress_mag)
      !$acc exit data delete(ss%tau_x, ss%tau_y, ss%stress_mag)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(ss%tau_x(i, j) - 0.5_wp*TAU) < TOL, &
                       "T1: midpoint blend must be half amplitude")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_tau_half_amp

   ! =================================================================
   ! T2 — ghosts are NOT extrapolated; periodic ghosts come from the wrap.
   ! =================================================================

   subroutine test_ghosts_extended(error)
      !! The inverse of what an earlier revision asserted.
      !!
      !! Stress ghosts are exchange-owned: `ocean_seam_refresh_surface_stress`
      !! fills them by MPI exchange, periodic wrap, or fold — never by
      !! copying this rank's edge value outward.  Two halves:
      !!
      !!   (a) NON-PERIODIC single rank: every ghost is BC-owned, so the
      !!       refresh must leave them exactly as it found them.  A
      !!       zero-gradient extension (the old behaviour) would fail this.
      !!   (b) PERIODIC single rank: the ghosts must equal their wrapped
      !!       interior partner, which is a genuinely different value from
      !!       the adjacent edge cell — so this also distinguishes a correct
      !!       wrap from an edge-copy.
      !!
      !! Both halves also confirm the physical block landed at
      !! `(nghost+1, nghost+1)`; an off-by-nghost landing would otherwise
      !! satisfy (a) trivially.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: SENTINEL = -777.0_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT)
      integer :: i, j, ilo, ihi, jlo, jhi

      t_axis = [0.0_wp, 1000.0_wp]
      fname = "/tmp/test_data_forcing_ghosts.nc"
      call write_gradient_file(trim(fname), NXP + 1, NYP, NT, t_axis)

      grid = make_grid()
      ilo = NGHOST + 1
      ihi = NGHOST + NXP + 1
      jlo = NGHOST + 1
      jhi = NGHOST + NYP

      ! ---- (a) non-periodic: ghosts must be left untouched ----
      call make_slots(grid, ss, sf)
      call reader%init()
      call tau_x_config(cfg, trim(fname), "field")
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      ss%tau_x = SENTINEL
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ss%tau_x, ss%tau_y, ss%stress_mag)
      !$acc update device(ss%tau_x)
      call ocean_data_input_update_all(reader, 0.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 0.0_wp)
      !$acc update self(ss%tau_x)
      !$acc exit data delete(ss%tau_x, ss%tau_y, ss%stress_mag)
      call reader%exit_data()
      !$acc exit data delete(reader)

      call check(error, abs(ss%tau_x(ilo, jlo) - 110.0_wp) < TOL, &
                 "T2a: file element (1,1) must land at (nghost+1, nghost+1)")
      if (allocated(error)) return
      do j = 1, size(ss%tau_x, 2)
         do i = 1, size(ss%tau_x, 1)
            if (i >= ilo .and. i <= ihi .and. j >= jlo .and. j <= jhi) cycle
            call check(error, abs(ss%tau_x(i, j) - SENTINEL) < TOL, &
                       "T2a: a non-periodic ghost must NOT be written "// &
                       "(no extrapolation from the edge)")
            if (allocated(error)) return
         end do
      end do
      call reader%destroy()
      call ss%destroy()
      call sf%destroy()

      ! ---- (b) periodic-x: ghosts must equal the wrapped partner ----
      call make_slots(grid, ss, sf)
      call reader%init()
      call tau_x_config(cfg, trim(fname), "field")
      bc%periodic_x = .true.
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      ss%tau_x = SENTINEL
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ss%tau_x, ss%tau_y, ss%stress_mag)
      !$acc update device(ss%tau_x)
      call ocean_data_input_update_all(reader, 0.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 0.0_wp)
      !$acc update self(ss%tau_x)
      !$acc exit data delete(ss%tau_x, ss%tau_y, ss%stress_mag)
      call reader%exit_data()
      !$acc exit data delete(reader)

      ! West ghosts wrap from the east physical block, and vice versa.
      ! The gradient file varies in i, so the wrapped value differs from
      ! the adjacent edge value — an edge-copy would fail here.
      do j = jlo, jhi
         do i = 1, NGHOST
            call check(error, abs(ss%tau_x(i, j) - ss%tau_x(i + NXP, j)) < TOL, &
                       "T2b: west ghost must equal its periodic partner")
            if (allocated(error)) return
         end do
      end do
      call check(error, abs(ss%tau_x(NGHOST, jlo) - ss%tau_x(ilo, jlo)) > 1.0_wp, &
                 "T2b: TEETH — the wrapped ghost must DIFFER from the adjacent "// &
                 "edge cell, else this test cannot tell a wrap from an edge-copy")
      if (allocated(error)) return

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_ghosts_extended

   ! =================================================================
   ! T3 — stress_mag tracks the blended wind across a bracket advance.
   ! =================================================================

   subroutine test_stress_mag_fresh(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 3
      real(wp), parameter :: TAU = 0.4_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), amp(NT)
      integer :: i, j

      ! Records at 0 / 100 / 200 s with amplitudes 0 / TAU / 2*TAU: the
      ! two queries below sit in DIFFERENT brackets, so the second only
      ! passes if the slab advance reached the device AND stress_mag was
      ! re-derived afterwards.
      t_axis = [0.0_wp, 100.0_wp, 200.0_wp]
      amp = [0.0_wp, TAU, 2.0_wp*TAU]
      fname = "/tmp/test_data_forcing_stress_mag.nc"
      call write_ramp_file(trim(fname), NXP + 1, NYP, NT, t_axis, amp)

      grid = make_grid()
      call make_slots(grid, ss, sf)
      call reader%init()
      call tau_x_config(cfg, trim(fname), "field")
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ss%tau_x, ss%tau_y, ss%stress_mag)

      call ocean_data_input_update_all(reader, 50.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 50.0_wp)
      !$acc update self(ss%stress_mag)
      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(ss%stress_mag(i, j) - 0.5_wp*TAU) < TOL, &
                       "T3: stress_mag must equal |tau| in the first bracket")
            if (allocated(error)) return
         end do
      end do

      call ocean_data_input_update_all(reader, 150.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 150.0_wp)
      !$acc update self(ss%stress_mag)
      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(ss%stress_mag(i, j) - 1.5_wp*TAU) < TOL, &
                       "T3: stress_mag must follow the wind across a bracket advance")
            if (allocated(error)) return
         end do
      end do

      !$acc exit data delete(ss%tau_x, ss%tau_y, ss%stress_mag)
      call reader%exit_data()
      !$acc exit data delete(reader)

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_stress_mag_fresh

   ! =================================================================
   ! T4 — heat lands in Q_heat when the component set is off.
   ! =================================================================

   subroutine test_heat_to_q_heat(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: Q = 50.0_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), amp(NT)
      integer :: i, j

      t_axis = [0.0_wp, 1000.0_wp]
      amp = [Q, Q]
      fname = "/tmp/test_data_forcing_heat.nc"
      call write_ramp_file(trim(fname), NXP + 1, NYP, NT, t_axis, amp)

      grid = make_grid()
      call make_slots(grid, ss, sf)
      call reader%init()

      cfg%enable = .true.
      cfg%time_mode = "linear"
      cfg%heat%file = trim(fname)
      cfg%heat%var = "field"

      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      call check(error,.not. df%heat_to_component, &
                 "T4: with components off the heat tag must target Q_heat")
      if (allocated(error)) return
      call check(error, sf%has_heat, "T4: registering a heat file must latch has_heat")
      if (allocated(error)) return

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(sf%Q_heat)
      call ocean_data_input_update_all(reader, 500.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 500.0_wp)
      !$acc update self(sf%Q_heat)
      !$acc exit data delete(sf%Q_heat)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(sf%Q_heat(i, j) - Q) < TOL, &
                       "T4: Q_heat must carry the file value")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_heat_to_q_heat

   ! =================================================================
   ! T5 — a climatology wraps instead of aborting.
   ! =================================================================

   subroutine test_cyclic_wrap(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: TAU = 0.3_wp, PERIOD = 400.0_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      character(len=256) :: fname
      real(wp) :: t_axis(NT), amp(NT)
      real(wp) :: tq, expect
      integer :: i, j

      ! Records at 0 and 200 s, period 400 s.  A query at 300 s is on the
      ! (nt,1) seam bracket, halfway from the last record back round to
      ! the first: expect the mean of the two amplitudes.
      t_axis = [0.0_wp, 200.0_wp]
      amp = [0.0_wp, TAU]
      fname = "/tmp/test_data_forcing_cyclic.nc"
      call write_ramp_file(trim(fname), NXP + 1, NYP, NT, t_axis, amp)

      grid = make_grid()
      call make_slots(grid, ss, sf)
      call reader%init()
      call tau_x_config(cfg, trim(fname), "field")
      cfg%time_mode = "cyclic"
      cfg%cycle_period = PERIOD
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      tq = 300.0_wp
      expect = 0.5_wp*(TAU + 0.0_wp)
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(ss%tau_x, ss%tau_y, ss%stress_mag)
      call ocean_data_input_update_all(reader, tq)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, tq)
      !$acc update self(ss%tau_x)
      !$acc exit data delete(ss%tau_x, ss%tau_y, ss%stress_mag)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(ss%tau_x(i, j) - expect) < TOL, &
                       "T5: seam bracket must blend last->first, not abort")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_cyclic_wrap

   ! =================================================================
   ! T6 — component-set routing.
   ! =================================================================

   subroutine test_components_route(error)
      !! With `use_components` on, `ocean_surface_flux_assemble` rebuilds
      !! `Q_heat`/`Q_salt` from the component fields every thermo step.
      !! A file writer must therefore feed the COMPONENTS, or its values
      !! would be overwritten before anything read them.  This asserts
      !! the routing goes where it should AND that `Q_heat` is left alone
      !! — the second half is the part that would have silently rotted.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: Q = 40.0_wp, EV = -3.0e-5_wp
      real(wp), parameter :: LP = 7.0e-5_wp, SL = 1.5e-6_wp
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(hgrid_t) :: grid
      real(wp) :: t_axis(NT)
      integer :: i, j

      t_axis = [0.0_wp, 1000.0_wp]
      call write_ramp_file("/tmp/test_data_forcing_c_heat.nc", NXP, NYP, NT, t_axis, [Q, Q])
      call write_ramp_file("/tmp/test_data_forcing_c_evap.nc", NXP, NYP, NT, t_axis, [EV, EV])
      call write_ramp_file("/tmp/test_data_forcing_c_lprec.nc", NXP, NYP, NT, t_axis, [LP, LP])
      call write_ramp_file("/tmp/test_data_forcing_c_salt.nc", NXP, NYP, NT, t_axis, [SL, SL])

      grid = make_grid()
      call ss%init(grid)
      call sf%init(grid)
      ! The gate under test: allocate the component set BEFORE configure,
      ! exactly as the driver does (set_components precedes the forcing
      ! registration, which is what makes the destination resolvable).
      call sf%set_components(grid, .true.)
      call reader%init()

      cfg%enable = .true.
      cfg%time_mode = "linear"
      cfg%heat%file = "/tmp/test_data_forcing_c_heat.nc"
      cfg%heat%var = "field"
      cfg%evap%file = "/tmp/test_data_forcing_c_evap.nc"
      cfg%evap%var = "field"
      cfg%lprec%file = "/tmp/test_data_forcing_c_lprec.nc"
      cfg%lprec%var = "field"
      cfg%salt%file = "/tmp/test_data_forcing_c_salt.nc"
      cfg%salt%var = "field"

      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)

      call check(error, df%heat_to_component .and. df%salt_to_component, &
                 "T6: components on must route heat/salt to the component fields")
      if (allocated(error)) return
      call check(error, sf%has_heat .and. sf%has_salt .and. sf%has_mass_flux, &
                 "T6: all three flux latches must be set")
      if (allocated(error)) return

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(sf%heat_added, sf%evap, sf%lprec, sf%salt_flux, sf%Q_heat)
      call ocean_data_input_update_all(reader, 500.0_wp)
      call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, 500.0_wp)
      !$acc update self(sf%heat_added, sf%evap, sf%lprec, sf%salt_flux, sf%Q_heat)
      !$acc exit data delete(sf%heat_added, sf%evap, sf%lprec, sf%salt_flux, sf%Q_heat)
      call reader%exit_data()
      !$acc exit data delete(reader)

      do j = NGHOST + 1, NGHOST + NYP
         do i = NGHOST + 1, NGHOST + NXP
            call check(error, abs(sf%heat_added(i, j) - Q) < TOL, &
                       "T6: heat must land in heat_added")
            if (allocated(error)) return
            call check(error, abs(sf%evap(i, j) - EV) < TOL, &
                       "T6: evap must land in evap (negative-definite convention preserved)")
            if (allocated(error)) return
            call check(error, abs(sf%lprec(i, j) - LP) < TOL, &
                       "T6: lprec must land in lprec")
            if (allocated(error)) return
            call check(error, abs(sf%salt_flux(i, j) - SL) < TOL, &
                       "T6: salt must land in salt_flux")
            if (allocated(error)) return
            ! The routing rule's whole purpose: Q_heat is the assembler's
            ! output, never a file destination under components.
            call check(error, abs(sf%Q_heat(i, j)) < TOL, &
                       "T6: Q_heat must be left for the assembler, not written directly")
            if (allocated(error)) return
         end do
      end do

      call reader%destroy()
      call ss%destroy()
      call sf%destroy()
   end subroutine test_components_route

   ! =================================================================
   ! Dyn-step scaffold for T7/T8.
   !
   ! Deliberately a local copy of `test_ocean_validation`'s init/map/
   ! destroy helpers rather than a shared module: that file carries no
   ! NetCDF dependency, and making it link the reader to share four
   ! boilerplate subroutines would drag NetCDF into a test that has no
   ! business needing it.
   ! =================================================================

   subroutine ek_init(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      integer :: k

      call grid%init(EK_NX, EK_NY, NGHOST, EK_DX, EK_DY)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = EK_F0
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call sf%init(grid)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid)

      ! Uniform resting column, uniform (T, S) — so the ONLY thing that
      ! can move the fluid is the surface stress.
      ms%h_layer = EK_H_LAYER
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*EK_H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*EK_H_LAYER
         ms%rho_layer(:, :, k) = eos%rho0
      end do
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      ! vmix off: keep the discriminator Coriolis-only, matching
      ! `wind_ekman_spinup_f_plane`.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
   end subroutine ek_init

   subroutine ek_map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
   end subroutine ek_map_in

   subroutine ek_map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      call destroy_cartesian_metrics(metrics)
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call sf%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine ek_map_out

   subroutine ek_destroy(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call sf%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine ek_destroy

   subroutine ek_run_const(tau, n_steps, u_top, v_top, mass_drift)
      !! Reference run: spatially-uniform constant wind via the existing
      !! `set_wind_stress_const` setter (seeded pre-map, as production
      !! does), no reader involved at all.
      real(wp), intent(in) :: tau
      integer, intent(in) :: n_steps
      real(wp), intent(out) :: u_top, v_top, mass_drift

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      real(wp) :: mass_ic
      integer :: step

      call ek_init(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
      call ss%set_wind_stress_const(tau, 0.0_wp)
      mass_ic = sum(ms%h_layer)

      call ek_map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      do step = 1, n_steps
         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, &
                             vmix, ms, EK_DT, sf=sf)
      end do
      call ek_map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)

      u_top = ms%u_face_x_layer(grid%nx_total/2, grid%ny_total/2, NZ)
      v_top = ms%v_face_y_layer(grid%nx_total/2, grid%ny_total/2, NZ)
      mass_drift = abs(sum(ms%h_layer) - mass_ic)/abs(mass_ic)

      call ek_destroy(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine ek_run_const

   subroutine ek_run_file(file, n_steps, t_start, u_top, v_top, mass_drift)
      !! File-driven run: identical setup, but the wind arrives every
      !! step through `update_all` + `ocean_data_forcing_apply`, exactly
      !! as the driver sequences it.  `t_start` seeds the model clock so
      !! the caller can land the query anywhere in the file's bracket.
      character(len=*), intent(in) :: file
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: t_start
      real(wp), intent(out) :: u_top, v_top, mass_drift

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_dataovr_config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(ocean_data_forcing_t) :: df
      type(ocean_bc_state_t) :: bc
      real(wp) :: mass_ic, t
      integer :: step

      call ek_init(grid, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)

      call reader%init()
      call tau_x_config(cfg, file, "field")
      call ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, df)
      mass_ic = sum(ms%h_layer)

      call ek_map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)
      !$acc enter data copyin(reader)
      call reader%enter_data()

      t = t_start
      do step = 1, n_steps
         ! The driver's order, verbatim: refresh brackets, blend into the
         ! slots, then step.
         call ocean_data_input_update_all(reader, t)
         call ocean_data_forcing_apply(df, reader, grid, ss, sf, bc, t)
         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, &
                             vmix, ms, EK_DT, sf=sf)
         t = t + EK_DT
      end do

      call reader%exit_data()
      !$acc exit data delete(reader)
      call ek_map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix)

      u_top = ms%u_face_x_layer(grid%nx_total/2, grid%ny_total/2, NZ)
      v_top = ms%v_face_y_layer(grid%nx_total/2, grid%ny_total/2, NZ)
      mass_drift = abs(sum(ms%h_layer) - mass_ic)/abs(mass_ic)

      call reader%destroy()
      call ek_destroy(ms, ct, cor, pgf, hv, bd, ss, sf, va, hd, vd, vmix, eos, dyn)
   end subroutine ek_run_file

   ! =================================================================
   ! T7 — file-driven wind reproduces the constant-wind Ekman spin-up.
   ! =================================================================

   subroutine test_ekman_matches_const(error)
      !! The wiring test.  Everything before this proves the file value
      !! reaches `surface_stress%tau_x`; this proves `tau_x` reaches the
      !! momentum equation, by driving the full RK2 dyn step with a
      !! TIME-CONSTANT forcing file and demanding round-off agreement
      !! with the equivalent `set_wind_stress_const` run.
      !!
      !! The physical-signature checks are not redundant with the
      !! comparison: if both runs were broken the same way (e.g. the
      !! stress never reached the top layer at all), the difference would
      !! still be zero.  Asserting the Ekman signature makes the
      !! agreement meaningful.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp) :: t_axis(NT)
      real(wp) :: u_ref, v_ref, drift_ref
      real(wp) :: u_file, v_file, drift_file
      character(len=256) :: fname

      ! Records spanning well past the integration so the run stays
      ! inside one bracket, both records equal so the blend is constant
      ! in time for ANY weight — this isolates the wiring from the
      ! time-interpolation logic (T1/T8 cover the weight).
      t_axis = [0.0_wp, 1.0e6_wp]
      fname = "/tmp/test_data_forcing_ekman_const.nc"
      call write_ramp_file(trim(fname), EK_NX + 1, EK_NY, NT, t_axis, [EK_TAU, EK_TAU])

      call ek_run_const(EK_TAU, EK_STEPS, u_ref, v_ref, drift_ref)
      call ek_run_file(trim(fname), EK_STEPS, 0.0_wp, u_file, v_file, drift_file)

      ! --- the physical Ekman signature must actually be present ---
      call check(error, u_ref > 1.0e-5_wp, &
                 "T7: reference run produced no surface u — stress not reaching the layer")
      if (allocated(error)) return
      call check(error, v_ref < 0.0_wp, &
                 "T7: reference run shows no southward Ekman turning")
      if (allocated(error)) return

      ! --- and the file-driven run must reproduce it to round-off ---
      call check(error, abs(u_file - u_ref) <= 1.0e-13_wp*max(abs(u_ref), 1.0e-30_wp), &
                 "T7: file-driven surface u differs from the constant-wind run")
      if (allocated(error)) return
      call check(error, abs(v_file - v_ref) <= 1.0e-13_wp*max(abs(v_ref), 1.0e-30_wp), &
                 "T7: file-driven surface v differs from the constant-wind run")
      if (allocated(error)) return
      call check(error, drift_file < 1.0e-12_wp, &
                 "T7: file-driven run broke closed-domain mass conservation")
   end subroutine test_ekman_matches_const

   ! =================================================================
   ! T8 — half amplitude at the bracket midpoint, through the solver.
   ! =================================================================

   subroutine test_ekman_half_amp(error)
      !! The half-amplitude claim from the plan, stated exactly: one step
      !! taken at the midpoint of a 0 -> TAU bracket must produce half
      !! the surface velocity of one step under the full TAU.
      !!
      !! ONE step, from rest, on purpose.  The response is then linear in
      !! tau to within the momentum-advection term (quadratic in
      !! velocity, ~1e-4 relative here), which is why the tolerance is
      !! 1e-3 relative and not round-off.  Over many steps the file wind
      !! would keep ramping while the reference stayed constant, so a
      !! multi-step form of this comparison would not be a clean
      !! statement about the blend weight.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NT = 2
      real(wp), parameter :: TEND = 1000.0_wp
      real(wp) :: t_axis(NT)
      real(wp) :: u_full, v_full, drift_full
      real(wp) :: u_half, v_half, drift_half
      character(len=256) :: fname

      ! Record 1 = 0, record 2 = EK_TAU: the query at TEND/2 blends to
      ! exactly EK_TAU/2.
      t_axis = [0.0_wp, TEND]
      fname = "/tmp/test_data_forcing_ekman_ramp.nc"
      call write_ramp_file(trim(fname), EK_NX + 1, EK_NY, NT, t_axis, [0.0_wp, EK_TAU])

      call ek_run_const(EK_TAU, 1, u_full, v_full, drift_full)
      call ek_run_file(trim(fname), 1, 0.5_wp*TEND, u_half, v_half, drift_half)

      call check(error, u_full > 1.0e-9_wp, &
                 "T8: full-amplitude single step produced no surface u")
      if (allocated(error)) return
      call check(error, abs(u_half - 0.5_wp*u_full) < 1.0e-3_wp*abs(u_full), &
                 "T8: midpoint-blend response is not half the full-amplitude response")
      if (allocated(error)) return
      call check(error, drift_half < 1.0e-12_wp, &
                 "T8: half-amplitude run broke closed-domain mass conservation")
   end subroutine test_ekman_half_amp

   ! =================================================================
   ! T9 — the config predicates behind the fail-loud paths.
   ! =================================================================

   subroutine test_config_predicates(error)
      !! Each predicate is exercised on BOTH branches.  A one-sided test
      !! (`ok` for a good config) would pass against a predicate hard-
      !! wired to `.true.`, which is precisely the regression that would
      !! let a broken namelist through.
      type(error_type), allocatable, intent(out) :: error

      type(ocean_dataovr_config_t) :: cfg
      type(dataovr_entry_config_t) :: e

      ! --- dataovr_entry_is_valid ---
      e%file = ""
      e%var = ""
      call check(error, dataovr_entry_is_valid(e), "T9: blank entry is valid (tag not driven)")
      if (allocated(error)) return

      e%file = "wind.nc"
      e%var = "taux"
      call check(error, dataovr_entry_is_valid(e), "T9: file+var is valid")
      if (allocated(error)) return

      e%file = "wind.nc"
      e%var = ""
      call check(error,.not. dataovr_entry_is_valid(e), &
                 "T9: file WITHOUT var must be rejected")
      if (allocated(error)) return

      ! A var with no file is harmless — the tag simply never registers.
      e%file = ""
      e%var = "taux"
      call check(error, dataovr_entry_is_valid(e), "T9: var without file is harmless")
      if (allocated(error)) return

      ! --- dataovr_time_is_valid ---
      cfg%time_mode = "linear"
      cfg%cycle_period = 0.0_wp
      call check(error, dataovr_time_is_valid(cfg), &
                 "T9: linear mode must not require a cycle_period")
      if (allocated(error)) return

      cfg%time_mode = "cyclic"
      cfg%cycle_period = 0.0_wp
      call check(error,.not. dataovr_time_is_valid(cfg), &
                 "T9: cyclic without cycle_period must be rejected")
      if (allocated(error)) return

      cfg%time_mode = "cyclic"
      cfg%cycle_period = 31536000.0_wp
      call check(error, dataovr_time_is_valid(cfg), "T9: cyclic with a period is valid")
      if (allocated(error)) return

      ! --- dataovr_freshwater_needs_components (reports the PROBLEM) ---
      cfg%evap%file = ""
      cfg%lprec%file = ""
      call check(error,.not. dataovr_freshwater_needs_components(cfg, .false.), &
                 "T9: no freshwater tag => no component requirement")
      if (allocated(error)) return

      cfg%evap%file = "evap.nc"
      call check(error, dataovr_freshwater_needs_components(cfg, .false.), &
                 "T9: evap without components must be flagged")
      if (allocated(error)) return
      call check(error,.not. dataovr_freshwater_needs_components(cfg, .true.), &
                 "T9: evap WITH components is fine")
      if (allocated(error)) return

      cfg%evap%file = ""
      cfg%lprec%file = "prec.nc"
      call check(error, dataovr_freshwater_needs_components(cfg, .false.), &
                 "T9: lprec without components must be flagged too")
      if (allocated(error)) return

      ! --- dataovr_any_tag_set ---
      block
         type(ocean_dataovr_config_t) :: empty
         call check(error,.not. dataovr_any_tag_set(empty), &
                    "T9: a pristine group drives nothing")
         if (allocated(error)) return
      end block

      block
         type(ocean_dataovr_config_t) :: one
         one%salt%file = "salt.nc"
         call check(error, dataovr_any_tag_set(one), &
                    "T9: any single tag counts (salt is the last one checked)")
         if (allocated(error)) return
      end block
   end subroutine test_config_predicates

end module test_ocean_data_forcing
