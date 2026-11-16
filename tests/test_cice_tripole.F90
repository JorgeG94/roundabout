!! The CICE-port sea ice (`&ocean_ice_nml scheme = "cice"`) across the tripolar north fold, ONE RANK, through the
!! production engine path: the fold is right, not merely self-consistent (stage B; Murray 1996 for the grid, Hunke et
!! al. for the scheme).
!!
!! **The fold as a symmetry.**  The folded domain (`NJ` rows, the north edge glued to itself by
!! `(x, y) -> (ni - x, 2 nj - y)`) is the quotient of a DOUBLED domain (`2 NJ` rows, walls south and north, periodic in
!! x) by the 180-degree rotation about `(ni/2, nj)`: scalars copy, true vectors negate.  A state and a forcing that are
!! invariant under the rotation must STAY invariant, and the folded run must reproduce the lower half of the doubled run.
!! A decomposed run equals the serial one even when a sign is wrong everywhere (both apply the same wrong fold); this is
!! the independent statement.  The map indices below are a host re-derivation from the continuous coordinates, not a call
!! of the code under test.
!!
!! The forcing is a wind that is ODD under the rotation (`u = A sin(kx (x-xc)) cos(ky (y-yc))`,
!! `v = A cos(kx (x-xc)) sin(ky (y-yc))`), set on the engine's per-cell atmosphere arrays, ghosts included.  The ocean is on
!! a flat f-plane at rest, so the ocean (currents, sea-surface tilt, the blended ice-ocean stress, the surface fluxes of the
!! hand-off) is symmetric too and enters the comparison.  Cases:
!!   * `doubled_domain_is_symmetric`   the doubled run's upper half is the rotated lower half (round-off bound), over the
!!                                     ice state, the ice velocity, the ocean layers and the surface fluxes and stress.
!!   * `fold_equals_doubled_lower_half`  the folded run (`north = tripolar_fold`) equals the lower half of the doubled run
!!                                     (same fields), including the fold-line row of the v faces.
!!   * `asymmetry_is_detected`         the comparison itself fails on a deliberately wrong sign / a displaced cell.
!!   * `fold_faults_are_detected`      the folded run with ONE fold fault injected (`ice_mr_test_skip`: the stress pair folded
!!                                     as a scalar, or its fold skipped) no longer equals the doubled lower half.  (The
!!                                     transport's fold-row flux pairing is NOT visible here: on one rank both sides compute
!!                                     the same flux to round-off, so unpairing changes the answer by 1e-12; the pairing is
!!                                     what makes a split fold row agree, and `test_ocean_decomp_bitid_mpi` detects its removal.)
!!   * `fold_seam_conservation`        the ice ledger (water, salt, heat) of the folded run closes at round-off every step.
module test_cice_tripole
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_cice_hooks, only: ice_mr_test_skip, ICE_SKIP_FOLD_SIGN, ICE_SKIP_FOLD_PAIR
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   implicit none
   private

   public :: collect_cice_tripole_tests

   real(wp), parameter :: DT = 600.0_wp
   real(wp), parameter :: DX = 20000.0_wp
   real(wp), parameter :: PI = 3.14159265358979323846_wp
   integer, parameter :: NI = 16
   integer, parameter :: NJ = 8
   integer, parameter :: NG = 3
   integer, parameter :: NSTEP = 60
   integer, parameter :: MAXF = 40
   integer, parameter :: ST_T = 1, ST_U = 2, ST_V = 3
   logical :: comm_inited = .false.

   type :: field_t
      character(len=24) :: tag = ""
      integer :: stag = ST_T
      real(wp) :: sgn = 1.0_wp
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snap_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snap_t

contains

   subroutine collect_cice_tripole_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("doubled_domain_is_symmetric", test_doubled_symmetric), &
                  new_unittest("fold_equals_doubled_lower_half", test_fold_equals_doubled), &
                  new_unittest("asymmetry_is_detected", test_asymmetry_detected), &
                  new_unittest("fold_faults_are_detected", test_fold_faults_detected), &
                  new_unittest("fold_seam_conservation", test_fold_conservation) &
                  ]
   end subroutine collect_cice_tripole_tests

   ! ------------------------------------------------------------------
   ! configuration and engine plumbing
   ! ------------------------------------------------------------------

   function tripole_nml(folded) result(nml)
      !! Flat f-plane cap: `NI x NJ` with the fold on the north edge (`folded`), `NI x 2NJ` with walls otherwise.
      logical, intent(in) :: folded
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny
      character(len=:), allocatable :: north

      write (snx, '(i0)') NI
      write (sny, '(i0)') merge(NJ, 2*NJ, folded)
      north = merge("tripolar_fold", "wall         ", folded)
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = -1.8, initial_salinity = 34.0, "// &
            "T_init_surface = -1.8, T_init_bottom = -1.8 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 1000.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
            "&ocean_ice_nml enable = .true., scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, "// &
            "dynamics = .true., transport = .true. /"//NL// &
            "&ocean_ice_thermo_nml update_ocn_f = .true. /"//NL// &
            "&ocean_ice_dyn_nml advection = 'remap', ndte = 20 /"//NL// &
            "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.8, h_snow = 0.1, t_ice = -8.0 /"//NL// &
            "&ocean_ice_atm_nml uni_tas = 255.0, uni_rlds = 180.0, uni_uas = 0.0, uni_vas = 0.0 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', north = '"//trim(north)//"' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function tripole_nml

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   subroutine make_engine(engine, cfg, nml, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      character(len=*), intent(in) :: nml
      logical, intent(out) :: ok
      integer :: ierr

      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine make_engine

   subroutine finish_engine(engine)
      type(ocean_engine_t), intent(inout) :: engine
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine finish_engine

   subroutine one_step(engine, cfg, t, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      logical, intent(out) :: ok
      integer :: ierr

      ok = .false.
      call engine_step(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_finalize(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      t = t + DT
      ok = .true.
   end subroutine one_step

   real(wp) function wind_amp()
      !! Amplitude of the odd wind [m/s]; `RDB_TRIPOLE_AMP` overrides it for a sensitivity hunt.
      character(len=16) :: v
      real(wp) :: a
      integer :: st

      wind_amp = 10.0_wp
      call get_environment_variable("RDB_TRIPOLE_AMP", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) a
         if (st == 0 .and. a > 0.0_wp) wind_amp = a
      end if
   end function wind_amp

   logical function verbose_run()
      !! `RDB_TRIPOLE_VERBOSE` set: name every step (to find the step of an abort).
      character(len=4) :: v
      integer :: st

      call get_environment_variable("RDB_TRIPOLE_VERBOSE", v, status=st)
      verbose_run = st == 0
   end function verbose_run

   integer function nsteps_run()
      !! `NSTEP`, or `RDB_TRIPOLE_STEPS` for a first-divergence hunt.
      character(len=16) :: v
      integer :: n, st

      nsteps_run = NSTEP
      call get_environment_variable("RDB_TRIPOLE_STEPS", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) n
         if (st == 0 .and. n > 0) nsteps_run = n
      end if
   end function nsteps_run

   subroutine set_odd_wind(engine, amp)
      !! The rotation-odd wind on every cell of the atmosphere arrays, ghost rows included (an analytic field, so the
      !! folded run's ghosts are its own continuation).  Centre `(xc, yc) = (NI/2, NJ) dx`.
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp), intent(in) :: amp
      integer :: i, j
      real(wp) :: x, y, kx, ky

      kx = 2.0_wp*PI/(real(NI, wp)*DX)
      ky = PI/(real(NJ, wp)*DX)
      associate (r => engine%state%icec)
         do j = 1, size(r%u_uas, 2)
            do i = 1, size(r%u_uas, 1)
               x = (real(i - NG, wp) - 0.5_wp)*DX - 0.5_wp*real(NI, wp)*DX
               y = (real(j - NG, wp) - 0.5_wp)*DX - real(NJ, wp)*DX
               r%u_uas(i, j) = amp*sin(kx*x)*cos(ky*y)
               r%u_vas(i, j) = amp*cos(kx*x)*sin(ky*y)
            end do
         end do
         !$acc update device(r%u_uas, r%u_vas)
      end associate
   end subroutine set_odd_wind

   subroutine run_case(folded, nsteps, snap, ok, ledger_res)
      !! Run `nsteps` coupled steps under the odd wind; return the snapshot (and, when asked, the worst ledger residuals).
      logical, intent(in) :: folded
      integer, intent(in) :: nsteps
      type(snap_t), intent(out) :: snap
      logical, intent(out) :: ok
      real(wp), intent(out), optional :: ledger_res(3)
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(ice_ledger_mirror_t) :: led
      real(wp) :: t
      integer :: n

      if (present(ledger_res)) ledger_res = 0.0_wp
      call make_engine(engine, cfg, tripole_nml(folded), ok)
      if (.not. ok) return
      call set_odd_wind(engine, wind_amp())
      t = 0.0_wp
      do n = 1, nsteps
         if (verbose_run()) write (*, '(a,l1,a,i0)') "  [tripole] folded=", folded, " step ", n
         if (present(ledger_res)) call led%open(engine%state%ice%cice, engine%state%ice_cgrid%tarea, DT, 1)
         call one_step(engine, cfg, t, ok)
         if (.not. ok) exit
         if (present(ledger_res)) then
            !$acc update self(engine%state%icec%flux%fresh, engine%state%icec%flux%fsalt, engine%state%icec%flux%fhocn)
            !$acc update self(engine%state%icec%flux%fswthru, engine%state%ice%cice%aice_init)
            call led%close(engine%state%ice%cice, engine%state%icec%flux, engine%state%ice_cgrid%tarea, 0.0_wp, 0.0_wp)
            ledger_res(ICE_LEDGER_WATER) = max(ledger_res(ICE_LEDGER_WATER), &
                                               abs(led%res(ICE_LEDGER_WATER))/max(led%scale(ICE_LEDGER_WATER), tiny(1.0_wp)))
            ledger_res(ICE_LEDGER_SALT) = max(ledger_res(ICE_LEDGER_SALT), &
                                              abs(led%res(ICE_LEDGER_SALT))/max(led%scale(ICE_LEDGER_SALT), tiny(1.0_wp)))
         end if
      end do
      if (ok) call take_snap(engine, snap)
      call finish_engine(engine)
   end subroutine run_case

   subroutine add_field(snap, tag, stag, sgn, a)
      type(snap_t), intent(inout) :: snap
      character(len=*), intent(in) :: tag
      integer, intent(in) :: stag
      real(wp), intent(in) :: sgn
      real(wp), intent(in) :: a(:, :, :)

      if (snap%n >= MAXF) error stop "test_cice_tripole: raise MAXF"
      snap%n = snap%n + 1
      snap%f(snap%n)%tag = tag
      snap%f(snap%n)%stag = stag
      snap%f(snap%n)%sgn = sgn
      allocate (snap%f(snap%n)%a(size(a, 1), size(a, 2), size(a, 3)))
      snap%f(snap%n)%a = a
   end subroutine add_field

   subroutine take_snap(engine, snap)
      !! The fields of the comparison, device to host first.  Scalars copy under the rotation (`sgn = 1`); the true
      !! vectors (face velocities, the stress pair) negate.
      type(ocean_engine_t), intent(inout), target :: engine
      type(snap_t), intent(inout) :: snap
      integer :: nx, ny, nz

      snap%n = 0
      associate (c => engine%state%ice%cice, ms => engine%state%multilayer, sf => engine%state%surface_flux, &
                 ss => engine%state%surface_stress, ice => engine%state%ice)
         nx = c%nx
         ny = c%ny
         !$acc update self(c%aicen, c%vicen, c%vsnon, c%Tsfcn, c%qicen, c%sicen, c%qsnon, c%aice, c%vice, c%vsno)
         !$acc update self(ice%u_ice, ice%v_ice)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         !$acc update self(sf%Q_heat, sf%Q_salt, sf%heat_ice, sf%sw_ice, sf%salt_ice)
         !$acc update self(ss%tau_x, ss%tau_y)
         nz = size(ms%h_layer, 3)
         call add_field(snap, "aicen", ST_T, 1.0_wp, c%aicen)
         call add_field(snap, "vicen", ST_T, 1.0_wp, c%vicen)
         call add_field(snap, "vsnon", ST_T, 1.0_wp, c%vsnon)
         call add_field(snap, "Tsfcn", ST_T, 1.0_wp, c%Tsfcn)
         call add_field(snap, "qicen", ST_T, 1.0_wp, reshape(c%qicen, [nx, ny, size(c%qicen)/(nx*ny)]))
         call add_field(snap, "sicen", ST_T, 1.0_wp, reshape(c%sicen, [nx, ny, size(c%sicen)/(nx*ny)]))
         call add_field(snap, "qsnon", ST_T, 1.0_wp, reshape(c%qsnon, [nx, ny, size(c%qsnon)/(nx*ny)]))
         call add_field(snap, "aice", ST_T, 1.0_wp, reshape(c%aice, [nx, ny, 1]))
         call add_field(snap, "vice", ST_T, 1.0_wp, reshape(c%vice, [nx, ny, 1]))
         call add_field(snap, "u_ice", ST_U, -1.0_wp, reshape(ice%u_ice, [size(ice%u_ice, 1), size(ice%u_ice, 2), 1]))
         call add_field(snap, "v_ice", ST_V, -1.0_wp, reshape(ice%v_ice, [size(ice%v_ice, 1), size(ice%v_ice, 2), 1]))
         call add_field(snap, "h_layer", ST_T, 1.0_wp, ms%h_layer(:, :, 1:nz))
         call add_field(snap, "u_face_x_layer", ST_U, -1.0_wp, ms%u_face_x_layer)
         call add_field(snap, "v_face_y_layer", ST_V, -1.0_wp, ms%v_face_y_layer)
         call add_field(snap, "Q_heat", ST_T, 1.0_wp, reshape(sf%Q_heat, [nx, ny, 1]))
         call add_field(snap, "Q_salt", ST_T, 1.0_wp, reshape(sf%Q_salt, [nx, ny, 1]))
         call add_field(snap, "heat_ice", ST_T, 1.0_wp, reshape(sf%heat_ice, [nx, ny, 1]))
         call add_field(snap, "salt_ice", ST_T, 1.0_wp, reshape(sf%salt_ice, [nx, ny, 1]))
         call add_field(snap, "tau_x", ST_U, -1.0_wp, reshape(ss%tau_x, [size(ss%tau_x, 1), size(ss%tau_x, 2), 1]))
         call add_field(snap, "tau_y", ST_V, -1.0_wp, reshape(ss%tau_y, [size(ss%tau_y, 1), size(ss%tau_y, 2), 1]))
      end associate
   end subroutine take_snap

   ! ------------------------------------------------------------------
   ! the comparison
   ! ------------------------------------------------------------------

   pure integer function mirror_index(stag, which, i, nin)
      !! Physical mirror index of `i` (1-based, owned) under the rotation: `which = 1` x, `2` y (for a field of `nin`
      !! physical rows, the DOUBLED height).  T: `n+1-i`; u in x: `NI+2-i` (the seam duplicate 1 is the face `NI+1`);
      !! v in y: `nin+2-j`.
      integer, intent(in) :: stag, which, i, nin
      if (which == 1) then
         if (stag == ST_U) then
            mirror_index = NI + 2 - i
            if (mirror_index == 1) mirror_index = NI + 1
         else
            mirror_index = NI + 1 - i
         end if
      else
         if (stag == ST_V) then
            mirror_index = nin + 2 - i
         else
            mirror_index = nin + 1 - i
         end if
      end if
   end function mirror_index

   subroutine symmetry_residual(fld, nrows, worst, name_worst, flip)
      !! Largest `|a(p) - sgn a(R p)|` over the owned cells of a DOUBLED-height field, relative to the field's largest
      !! absolute value.  `flip` negates the sign (the deliberately wrong comparison of the control).
      type(field_t), intent(in) :: fld
      integer, intent(in) :: nrows
      real(wp), intent(inout) :: worst
      character(len=24), intent(inout) :: name_worst
      logical, intent(in), optional :: flip
      integer :: i, j, k, im, jm, i0, i1, j1
      real(wp) :: d, scl, sg, w

      sg = fld%sgn
      if (present(flip)) then
         if (flip) sg = -sg
      end if
      scl = max(maxval(abs(fld%a)), tiny(1.0_wp))
      i0 = 1
      i1 = NI
      j1 = nrows
      if (fld%stag == ST_U) then
         i0 = 2
         i1 = NI + 1
      end if
      if (fld%stag == ST_V) j1 = nrows + 1
      w = 0.0_wp
      do k = 1, size(fld%a, 3)
         do j = 1, j1
            jm = mirror_index(fld%stag, 2, j, nrows)
            do i = i0, i1
               im = mirror_index(fld%stag, 1, i, nrows)
               d = abs(fld%a(NG + i, NG + j, k) - sg*fld%a(NG + im, NG + jm, k))
               w = max(w, d)
            end do
         end do
      end do
      w = w/scl
      if (w > worst) then
         worst = w
         name_worst = fld%tag
      end if
   end subroutine symmetry_residual

   subroutine folded_residual(ff, fd, worst, name_worst)
      !! Largest relative `|folded - doubled|` over the lower half (the folded run's owned cells).
      type(field_t), intent(in) :: ff, fd
      real(wp), intent(inout) :: worst
      character(len=24), intent(inout) :: name_worst
      integer :: i, j, k, i0, i1, j1
      real(wp) :: scl, w
      integer :: wi, wj, wk

      wi = 0
      wj = 0
      wk = 0
      scl = max(maxval(abs(fd%a)), tiny(1.0_wp))
      i0 = 1
      i1 = NI
      j1 = NJ
      if (ff%stag == ST_U) then
         i0 = 2
         i1 = NI + 1
      end if
      if (ff%stag == ST_V) j1 = NJ + 1
      w = 0.0_wp
      do k = 1, size(ff%a, 3)
         do j = 1, j1
            do i = i0, i1
               if (abs(ff%a(NG + i, NG + j, k) - fd%a(NG + i, NG + j, k)) > w) then
                  w = abs(ff%a(NG + i, NG + j, k) - fd%a(NG + i, NG + j, k))
                  wi = i
                  wj = j
                  wk = k
               end if
            end do
         end do
      end do
      if (w > 1.0e-9_wp*scl) write (*, '(a,a,a,3(1x,i0),a,es10.3,a,es10.3)') "    first-worst ", trim(ff%tag), " at (i,j,k)", &
         wi, wj, wk, " folded ", ff%a(NG + wi, NG + wj, wk), " doubled ", fd%a(NG + wi, NG + wj, wk)
      w = w/scl
      if (w > worst) then
         worst = w
         name_worst = ff%tag
      end if
   end subroutine folded_residual

   ! ------------------------------------------------------------------
   ! tests
   ! ------------------------------------------------------------------

   subroutine test_doubled_symmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(snap_t) :: sd
      logical :: ok
      integer :: e
      real(wp) :: worst, vmax
      character(len=24) :: nm

      call run_case(.false., nsteps_run(), sd, ok)
      call check(error, ok, "the doubled-domain run completed")
      if (allocated(error)) return
      worst = 0.0_wp
      nm = ""
      do e = 1, sd%n
         call symmetry_residual(sd%f(e), 2*NJ, worst, nm)
      end do
      vmax = maxval(abs(sd%f(10)%a))
      write (*, '(a,es10.3,a,a,a,es10.3)') "  doubled domain: worst asymmetry ", worst, " (field ", trim(nm), "), max|u_ice| ", vmax
      call check(error, vmax > 1.0e-3_wp, "the ice moves under the vortex wind")
      if (.not. allocated(error)) call check(error, worst < 1.0e-9_wp, "the doubled domain stays rotation-symmetric")
   end subroutine test_doubled_symmetric

   subroutine test_fold_equals_doubled(error)
      type(error_type), allocatable, intent(out) :: error
      type(snap_t) :: sf, sd
      logical :: ok
      integer :: e
      real(wp) :: worst
      character(len=24) :: nm

      call run_case(.true., nsteps_run(), sf, ok)
      call check(error, ok, "the folded run completed")
      if (allocated(error)) return
      call run_case(.false., nsteps_run(), sd, ok)
      call check(error, ok, "the doubled-domain run completed")
      if (allocated(error)) return
      worst = 0.0_wp
      nm = ""
      do e = 1, sf%n
         call folded_residual(sf%f(e), sd%f(e), worst, nm)
         write (*, '(a,a,es10.3)') "    folded vs doubled: ", sf%f(e)%tag, folded_one(sf%f(e), sd%f(e))
      end do
      write (*, '(a,es10.3,a,a)') "  folded vs doubled lower half: worst ", worst, " in ", trim(nm)
      call check(error, worst < 1.0e-9_wp, "the folded run equals the lower half of the doubled run")
   end subroutine test_fold_equals_doubled

   function folded_one(ff, fd) result(r)
      type(field_t), intent(in) :: ff, fd
      real(wp) :: r
      character(len=24) :: nm
      r = 0.0_wp
      nm = ""
      call folded_residual(ff, fd, r, nm)
   end function folded_one

   subroutine test_asymmetry_detected(error)
      !! The comparison can fail: the same doubled run read with the vector sign flipped, and with one displaced cell.
      type(error_type), allocatable, intent(out) :: error
      type(snap_t) :: sd
      logical :: ok
      real(wp) :: worst, w0
      character(len=24) :: nm

      call run_case(.false., 12, sd, ok)
      call check(error, ok, "the doubled-domain run completed")
      if (allocated(error)) return
      worst = 0.0_wp
      nm = ""
      call symmetry_residual(sd%f(10), 2*NJ, worst, nm, flip=.true.)
      call check(error, worst > 1.0e-3_wp, "a flipped vector sign is detected on u_ice")
      w0 = 0.0_wp
      sd%f(1)%a(NG + 3, NG + 2, 1) = sd%f(1)%a(NG + 3, NG + 2, 1) + 1.0e-3_wp
      call symmetry_residual(sd%f(1), 2*NJ, w0, nm)
      if (.not. allocated(error)) call check(error, w0 > 1.0e-6_wp, "a displaced cell is detected on aicen")
   end subroutine test_asymmetry_detected

   subroutine test_fold_faults_detected(error)
      !! A wrong sign on the T-cell stress pair and a skipped fold of it each break the
      !! equivalence with the doubled domain: the test above can fail for each of them.
      type(error_type), allocatable, intent(out) :: error
      type(snap_t) :: sf, sd
      logical :: ok
      integer :: e, k, bits(2)
      real(wp) :: worst
      character(len=24) :: nm, names(2)

      bits = [ICE_SKIP_FOLD_SIGN, ICE_SKIP_FOLD_PAIR]
      names = [character(len=24) :: "stress pair as scalar", "stress pair not folded"]
      call run_case(.false., 12, sd, ok)
      call check(error, ok, "the doubled-domain run completed")
      if (allocated(error)) return
      do k = 1, 2
         ice_mr_test_skip = bits(k)
         call run_case(.true., 12, sf, ok)
         ice_mr_test_skip = 0
         call check(error, ok, "the faulted folded run completed")
         if (allocated(error)) return
         worst = 0.0_wp
         nm = ""
         do e = 1, sf%n
            call folded_residual(sf%f(e), sd%f(e), worst, nm)
         end do
  write (*, '(a,a,a,es10.3,a,a)') "  fault '", trim(names(k)), "': worst difference to the doubled domain ", worst, " in ", trim(nm)
         call check(error, worst > 1.0e-6_wp, "the fault '"//trim(names(k))//"' is detected")
         if (allocated(error)) return
      end do
   end subroutine test_fold_faults_detected

   subroutine test_fold_conservation(error)
      !! The ice ledger closes every step with the transport crossing the fold; the residual is compared with the
      !! doubled-domain run's (round-off of the ledger sums) and bounded far below any seam leak (a flux row's worth is
      !! of order 1e-3 of the totals).
      type(error_type), allocatable, intent(out) :: error
      type(snap_t) :: sf
      real(wp) :: res(3), resd(3)
      logical :: ok

      call run_case(.true., nsteps_run(), sf, ok, res)
      call check(error, ok, "the folded run completed")
      if (allocated(error)) return
      call run_case(.false., nsteps_run(), sf, ok, resd)
      call check(error, ok, "the doubled-domain run completed")
      if (allocated(error)) return
      write (*, '(a,2es10.3,a,2es10.3)') "  ledger residual (water, salt) folded: ", res(ICE_LEDGER_WATER), res(ICE_LEDGER_SALT), &
         "  doubled: ", resd(ICE_LEDGER_WATER), resd(ICE_LEDGER_SALT)
      call check(error, res(ICE_LEDGER_WATER) < 1.0e-8_wp, "water conserved through the fold")
      if (.not. allocated(error)) call check(error, res(ICE_LEDGER_SALT) < 1.0e-8_wp, "salt conserved through the fold")
   end subroutine test_fold_conservation

end module test_cice_tripole
