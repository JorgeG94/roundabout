!! Engine-driven (namelist) dynamic wet/dry case on ONE rank: an emergent
!! Gaussian shoal whose wetting/drying front sweeps across the lines
!! `i = nx/2 | nx/2+1` and `j = ny/2 | ny/2+1` -- the seams of the 2x1, 1x2 and
!! 2x2 decompositions of `tests/mpi/test_ocean_decomp_bitid_mpi.F90`, whose
!! grid size (26 x 18), `dt` (900 s) and step count (48) this case reuses so
!! it drops into that gate unchanged.
!!
!! Every other wet/dry test builds its state by hand; this is the first that
!! goes `read_config_from_string` -> `validate_config` -> `engine_setup` ->
!! `engine_step`.  The test asserts
!!   (a) every prognostic stays finite,
!!   (b) the layer volume is conserved to round-off and no column goes
!!       negative (the gates of the sibling `test_ocean_wetdry*` suites),
!!   (c) the hysteresis mask `wd_wet_dyn` FLIPS in each seam band, on each
!!       side of each centre line and in the 2x2 corner block, in BOTH
!!       directions (wet->dry and dry->wet): the non-vacuity evidence a
!!       later multi-rank bit-identity case relies on.  The counts are printed.
!! Flips are sampled once per outer step (substep flips that revert within a
!! step are not seen, so the counts are a lower bound).
module test_ocean_wetdry_engine
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_wetdry_engine_tests

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: N_STEPS = 48
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: VOL_RTOL = 1.0e-12_wp
      !! Relative layer-volume drift: the round-off gate of the sibling
      !! `test_ocean_wetdry` / `test_ocean_wetdry_driver` conservation
      !! checks (closed walls, no sources; measured drift here ~2e-16).
   real(wp), parameter :: COL_FLOOR = -1.0e-9_wp
      !! Most negative column depth tolerated (sibling positivity gate).
   integer, parameter :: NBANDS = 9
   character(len=*), parameter :: BAND_NAME(NBANDS) = [character(len=14) :: &
                                                       "i-line west", "i-line east", &
                                                       "j-line south", "j-line north", &
                                                       "corner 4x4", &
                                                       "i west lower", "i west upper", &
                                                       "i east lower", "i east upper"]

   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_wetdry_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wetdry_shoal_engine_front_crosses_seams", test_shoal) &
                  ]
   end subroutine collect_ocean_wetdry_engine_tests

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function shoal_nml() result(nml)
      !! The case text (shared verbatim with the decomp gate).
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny

      write (snx, '(i0)') NX_G
      write (sny, '(i0)') NY_G
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'ssp_rk2' /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
            "dx = 4000.0, dy = 4000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.05, wind_stress_y = 0.05 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 2.0, "// &
            "edge_depth = 0.0, slope_scale = 28000.0, wind_config = 'constant' /"//NL// &
            "&ocean_continuity_nml ppm_limit_pos = .true. /"//NL// &
            "&ocean_wetdry_nml enable = .true. /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function shoal_nml

   subroutine test_shoal(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      integer :: ierr, n, i, j, k, ng, nx, ny, nz, b
      integer :: ilo(NBANDS), ihi(NBANDS), jlo(NBANDS), jhi(NBANDS)
      integer :: w2d(NBANDS), d2w(NBANDS)
      real(wp), allocatable :: prev(:, :)
      real(wp) :: t, vol, vol0, mind, mincol, col, emin, emax, dxy, drift
      logical :: finite_ok

      call ensure_comm()
      call read_config_from_string(shoal_nml(), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "read_config_from_string failed")
      if (allocated(error)) return
      call validate_config(cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "validate_config refused the shoal case")
      if (allocated(error)) return
      call engine_setup(engine, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "engine_setup failed")
      if (allocated(error)) return
      call engine_enter_data(engine, cfg)

      ng = engine%grid%nghost
      nx = NX_G
      ny = NY_G
      nz = size(engine%state%multilayer%h_layer, 3)
      dxy = engine%grid%dx*engine%grid%dy
      ! Bands in GLOBAL physical indices (serial: local = global + ng).
      ilo = 1; ihi = nx; jlo = 1; jhi = ny
      ilo(1) = nx/2 - 1; ihi(1) = nx/2          ! west of the i seam
      ilo(2) = nx/2 + 1; ihi(2) = nx/2 + 2      ! east of the i seam
      jlo(3) = ny/2 - 1; jhi(3) = ny/2          ! south of the j seam
      jlo(4) = ny/2 + 1; jhi(4) = ny/2 + 2      ! north of the j seam
      ilo(5) = nx/2 - 1; ihi(5) = nx/2 + 2; jlo(5) = ny/2 - 1; jhi(5) = ny/2 + 2
      ! i-seam bands split by the j seam (so each quadrant is seen).
      ilo(6) = nx/2 - 1; ihi(6) = nx/2; jhi(6) = ny/2
      ilo(7) = nx/2 - 1; ihi(7) = nx/2; jlo(7) = ny/2 + 1
      ilo(8) = nx/2 + 1; ihi(8) = nx/2 + 2; jhi(8) = ny/2
      ilo(9) = nx/2 + 1; ihi(9) = nx/2 + 2; jlo(9) = ny/2 + 1
      w2d = 0
      d2w = 0

      allocate (prev(size(engine%state%dyn%bt_work%wd_wet_dyn, 1), &
                     size(engine%state%dyn%bt_work%wd_wet_dyn, 2)))
      prev = engine%state%dyn%bt_work%wd_wet_dyn
      vol0 = volume()
      mincol = huge(1.0_wp)
      emin = huge(1.0_wp)
      emax = -huge(1.0_wp)
      finite_ok = .true.
      t = 0.0_wp
      ierr = OCEAN_STATUS_OK
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
         !$omp target update from(engine%state%multilayer%h_layer)
         !$omp target update from(engine%state%dyn%bt_work%wd_wet_dyn)
         !$omp target update from(engine%state%multilayer%u_face_x_layer)
         !$omp target update from(engine%state%multilayer%v_face_y_layer)
         if (.not. all(ieee_is_finite(engine%state%multilayer%h_layer)) .or. &
             .not. all(ieee_is_finite(engine%state%multilayer%u_face_x_layer)) .or. &
             .not. all(ieee_is_finite(engine%state%multilayer%v_face_y_layer))) then
            finite_ok = .false.
            exit
         end if
         do j = ng + 1, ng + ny
            do i = ng + 1, ng + nx
               col = 0.0_wp
               do k = 1, nz
                  col = col + engine%state%multilayer%h_layer(i, j, k)
               end do
               mincol = min(mincol, col)
               emin = min(emin, col - engine%state%barotropic%b(i, j))
               emax = max(emax, col - engine%state%barotropic%b(i, j))
            end do
         end do
         do b = 1, NBANDS
            do j = jlo(b), jhi(b)
               do i = ilo(b), ihi(b)
                  associate (now => engine%state%dyn%bt_work%wd_wet_dyn(i + ng, j + ng), &
                             was => prev(i + ng, j + ng))
                     if (was > 0.5_wp .and. now < 0.5_wp) w2d(b) = w2d(b) + 1
                     if (was < 0.5_wp .and. now > 0.5_wp) d2w(b) = d2w(b) + 1
                  end associate
               end do
            end do
         end do
         prev = engine%state%dyn%bt_work%wd_wet_dyn
      end do
      call check(error, ierr == OCEAN_STATUS_OK .and. finite_ok, &
                 "shoal run failed or produced non-finite prognostics")
      if (.not. allocated(error)) then
         vol = volume()
         drift = abs(vol - vol0)/abs(vol0)
         mind = mincol
         print '(a,es10.3,a,es10.3,a,es10.3)', "wetdry_shoal: vol drift ", drift, &
            "  min column depth ", mind, "  eta range [col-b] ", emin
         print '(a,es10.3)', "wetdry_shoal: max (col-b) ", emax
         do b = 1, NBANDS
            print '(a,a14,a,i5,a,i5)', "wetdry_shoal: band ", BAND_NAME(b), &
               "  wet->dry ", w2d(b), "  dry->wet ", d2w(b)
         end do
         call check(error, drift <= VOL_RTOL, "layer volume not conserved")
         if (.not. allocated(error)) call check(error, mind >= COL_FLOOR, "negative column depth")
         do b = 1, NBANDS
            if (allocated(error)) exit
            call check(error, w2d(b) > 0 .and. d2w(b) > 0, &
                       "seam band '"//trim(BAND_NAME(b))//"' saw no wet<->dry flip in one direction")
         end do
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)

   contains

      function volume() result(v)
         real(wp) :: v
         integer :: ii, jj, kk
         v = 0.0_wp
         do kk = 1, nz
            do jj = ng + 1, ng + ny
               do ii = ng + 1, ng + nx
                  v = v + engine%state%multilayer%h_layer(ii, jj, kk)*dxy
               end do
            end do
         end do
      end function volume

   end subroutine test_shoal

end module test_ocean_wetdry_engine
