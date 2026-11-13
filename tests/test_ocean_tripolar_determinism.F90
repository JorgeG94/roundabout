!! Analytic tripolar grid through the PRODUCTION engine: two identical
!! short runs must stay physical and end BITWISE identical, ghosts
!! included, on every restart-registry field.
!!
!! Why this exists (2026-10-04 compatibility matrix, GPU leg): every
!! tripolar cell either stopped in steps 1-3 or ended 1e-6..1e-3 away
!! from gfortran on the nvfortran cc70 build.  The cause was not a device
!! race but the analytic cap generator: the column of supergrid nodes
!! through a cap pole maps onto the pole, so its along-j face (`dyCu`)
!! and its corners are geometrically zero -- but the bipolar map reached
!! that pole through finite transcendental round-off.  Under nvfortran's
!! `-fast` libm the east periodic seam column (pseudo-longitude 360, the
!! image of the lon_pole = 0 pole) landed ~1e-9 m off the pole where
!! gfortran happened to land exactly on it.  A 1e-9 m face passes every
!! `/= 0` reciprocal guard, so `idyCu`, `iareaCu`, `dx_dyBu` came out
!! 1e3..1e14 at the pole-adjacent seam faces, the Smagorinsky / biharmonic
!! viscosity there hit its CFL cap, and `hvisc_du_visc` reached 1e16 in
!! OWNED faces in the second stage of step 1.  The fix snaps every cap
!! node of a pole column exactly onto the pole
!! (`rdb_ocean_metrics::tripolar_node_latlon`), so the face is exactly
!! zero on every compiler.
!!
!! The gate is the matrix's own tripolar geometry (24 x 16 cells of
!! 15 x 1 degree from 59N, cap above 70N, lon_pole = 0, nghost = 4) with
!! the lateral closure that turned the sliver into a blow-up (Smagorinsky
!! Kh + Ah).  Run A and run B are configured and stepped independently;
!! each must stay finite and physically bounded (a 0.1 Pa wind cannot
!! drive 1 m/s in N_STEPS * DT = 1.5 h), and every restart-registry field
!! -- the full local array, ghost ring and fold rows included -- must
!! match bitwise between them.  The bound is what fails on the unfixed
!! GPU build (|u| ~ 1e4 m/s at step 2); the bitwise comparison is the
!! run-to-run determinism gate for the fold / cap path on the device.
module test_ocean_tripolar_determinism
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_tripolar_determinism_tests

   integer, parameter :: N_STEPS = 6
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: U_BOUND = 1.0_wp
      !! m/s.  The wind-driven flow after 1.5 h is O(1e-3) m/s; the
      !! sliver-face blow-up is O(1e4).  Any bound between is the same gate.

   logical :: comm_inited = .false.

   type :: field_snap_t
      !! Host copy of one registry entry (full local extent, ghosts included).
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:)
   end type field_snap_t

contains

   subroutine collect_ocean_tripolar_determinism_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("tripolar_engine_two_runs_bitwise_and_bounded", &
                               test_two_runs_bitwise) &
                  ]
   end subroutine collect_ocean_tripolar_determinism_tests

   subroutine ensure_comm()
      !! See `test_ocean_budget_periodic_sponge_serial::ensure_comm`.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function case_nml() result(nml)
      !! The compatibility matrix's `geometry = tripolar` row on its base
      !! closures (`tests/regression/compat_matrix.py`), with an analytic
      !! flat bed and a linear T profile in place of its input files.
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 4, dx = 15.0, dy = 1.0 /"//NL// &
            "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, lat_south = 59.0, "// &
            "phi_join = 70.0, lon_pole = 0.0, coriolis_scheme = 'planetary' /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'tripolar_fold' /"//NL// &
            "&nonhydrostatic_nml nz_layers = 5 /"//NL// &
            "&time_nml t_end = 5400.0, dt_fixed = 900.0 /"//NL// &
            "&tracer_nml initial_temperature = 9.0, initial_salinity = 34.8, "// &
            "T_init_surface = 16.0, T_init_bottom = 3.0 /"//NL// &
            "&ocean_topo_nml max_depth = 2000.0 /"//NL// &
            "&physics_nml wind_stress_x = 0.1 /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "c_smag = 0.15, smag_ah = .true. /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine run_case(snaps, max_u, n_nonfinite, ok)
      !! Configure, map, step N_STEPS, pull every registry field host-ward.
      type(field_snap_t), allocatable, intent(out) :: snaps(:)
      real(wp), intent(out) :: max_u
         !! max |u|, |v| over the layer velocities (whole storage).
      integer, intent(out) :: n_nonfinite
         !! Non-finite values over every registry field.
      logical, intent(out) :: ok
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      type(restart_registry_t) :: reg
      real(wp) :: t
      integer :: ierr, n, e

      ok = .false.
      max_u = huge(1.0_wp)
      n_nonfinite = 0
      call ensure_comm()
      call read_config_from_string(case_nml(), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr /= OCEAN_STATUS_OK) then
         ! Unmap and tear down on the failure path too, so a failed
         ! step leaves no device mappings behind for the next run.
         call engine_exit_data(engine)
         call engine_teardown(engine)
         return
      end if

      ! The restart registry is the checkpoint's field set: pull each
      ! device-mapped entry down by its COMPONENT pointer (never an
      ! aggregate derived type -- see CLAUDE.md), exactly as
      ! `ocean_state_restart_write` does, and keep the full local array.
      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      allocate (snaps(reg%n))
      do e = 1, reg%n
         associate (en => reg%entries(e))
            snaps(e)%tag = en%tag
            select case (en%rank)
            case (0)
               snaps(e)%a = [en%p0]
            case (2)
               if (en%device_mapped) then
                  !$omp target update from(en%p2)
               end if
               snaps(e)%a = reshape(en%p2, [size(en%p2)])
            case (3)
               if (en%device_mapped) then
                  !$omp target update from(en%p3)
               end if
               snaps(e)%a = reshape(en%p3, [size(en%p3)])
            end select
            n_nonfinite = n_nonfinite + count(.not. ieee_is_finite(snaps(e)%a))
         end associate
      end do
      ! Explicit pull: the registry loop above already refreshes these
      ! (`ml_u_face_x_layer` / `ml_v_face_y_layer` alias them), but the
      ! bound must not depend on that aliasing on the mem:separate build.
      !$omp target update from(engine%state%multilayer%u_face_x_layer, &
      !$omp&            engine%state%multilayer%v_face_y_layer)
      max_u = max(maxval(abs(engine%state%multilayer%u_face_x_layer)), &
                  maxval(abs(engine%state%multilayer%v_face_y_layer)))
      call reg%clear()
      call engine_exit_data(engine)
      call engine_teardown(engine)
      ok = .true.
   end subroutine run_case

   subroutine test_two_runs_bitwise(error)
      type(error_type), allocatable, intent(out) :: error
      type(field_snap_t), allocatable :: sa(:), sb(:)
      real(wp) :: umax_a, umax_b
      integer :: nf_a, nf_b, e, n_diff
      logical :: ok
      character(len=256) :: msg

      call run_case(sa, umax_a, nf_a, ok)
      call check(error, ok, "run A failed (setup or step)")
      if (allocated(error)) return
      call run_case(sb, umax_b, nf_b, ok)
      call check(error, ok, "run B failed (setup or step)")
      if (allocated(error)) return

      write (msg, '(a,i0,a,i0)') "non-finite registry values: run A ", nf_a, ", run B ", nf_b
      call check(error, nf_a == 0 .and. nf_b == 0, trim(msg))
      if (allocated(error)) return
      write (msg, '(a,es10.3,a,es10.3,a,es8.1,a)') "max |u|,|v| run A ", umax_a, &
         ", run B ", umax_b, " m/s > ", U_BOUND, &
         " (degenerate cap-pole face? see the module header)"
      call check(error, umax_a < U_BOUND .and. umax_b < U_BOUND, trim(msg))
      if (allocated(error)) return

      call check(error, size(sa) == size(sb) .and. size(sa) > 0, &
                 "restart registries differ in size between the two runs")
      if (allocated(error)) return
      do e = 1, size(sa)
         call check(error, sa(e)%tag == sb(e)%tag .and. size(sa(e)%a) == size(sb(e)%a), &
                    "registry entry mismatch at "//trim(sa(e)%tag))
         if (allocated(error)) return
         ! Bitwise: compare the IEEE bit patterns, so -0.0 /= +0.0 and a
         ! NaN payload counts (a value compare would hide both).
         n_diff = count(transfer(sa(e)%a, 0_int64, size(sa(e)%a)) /= &
                        transfer(sb(e)%a, 0_int64, size(sb(e)%a)))
         write (msg, '(a,a,a,i0,a,i0,a)') "registry field ", trim(sa(e)%tag), &
            " differs between two identical runs in ", n_diff, " of ", size(sa(e)%a), &
            " values (ghosts included)"
         call check(error, n_diff == 0, trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_two_runs_bitwise

end module test_ocean_tripolar_determinism
