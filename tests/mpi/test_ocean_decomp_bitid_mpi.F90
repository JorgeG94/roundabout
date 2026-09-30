!! Dycore bit-identity across domain decompositions, through the
!! production engine.
!!
!! "MPI works" means a decomposed run IS the serial run: every owned value
!! of every prognostic field bit-for-bit, not an integral that agrees to
!! round-off.  A round-off-level difference is where every seam bug
!! starts (a one-sided stencil at a seam face, an allreduce feeding back
!! into the state, a seam treated as a wall, a global-index formula
!! evaluated on local indices); integrals hide it for hundreds of steps.
!!
!! On EACH launch (ctest runs np = 1, 2, 4) the binary, for every case
!! below:
!!   1. runs the SERIAL reference (px = py = 1, compute_size = 1) on every
!!      rank — deterministic, no collective that could differ — through
!!      `engine_setup` -> `engine_enter_data` -> `engine_step` x N_STEPS;
!!   2. runs every px x py factorisation of the launched rank count
!!      (np = 2: 2x1, 1x2; np = 4: 4x1, 1x4, 2x2) on the same namelist;
!!   3. compares, BITWISE, each rank's OWNED window (ghosts excluded; a
!!      staggered face array includes both of its tile's edge faces, so a
!!      seam face is checked on both ranks that compute it) of every field
!!      the restart registry checkpoints — h, u, v, every tracer hTr, the
!!      barotropic prognostics, the KPP/EPBL/kappa-shear persistent state,
!!      rho, the OBC reservoirs — plus the BT end-of-step eta, against the
!!      matching window of the reference.
!!
!! Cases (each under `pred_corr` AND `ssp_rk2`, 48 steps):
!!   * island_basin — closed Cartesian basin with an interior land block
!!     (static land mask), beta plane, `2gyre` wind, sigma;
!!   * periodic_channel_zstar — re-entrant channel over a seamount on z*
!!     (ALE remap every step), with porous barriers;
!!   * open_obc — tidal west edge (eta target) + Flather east edge,
!!     zstar_sigma, over a seamount;
!!   * spherical — lon-lat sector, planetary Coriolis, spoon basin, Wright
!!     EOS, `2gyre` wind;
!!   * obc_radiation_sponge — Orlanski-radiating open south edge with tracer
!!     reservoirs, clamped north inflow, legacy relaxing sponge band west;
!!   * closures — the spherical case with EPBL (instead of KPP), Fox-Kemper
!!     MLE, GM + MEKE, Redi, kappa-shear, tidal mixing, convective
!!     adjustment, geothermal heating and tracer hdiff;
!!   * file_readers — the per-rank windowed readers (a periodic 360-degree
!!     MOM6 supergrid, a C-order bathymetry file with land, a z-level T/S
!!     IC), all written by rank 0 first, with the global-1-degree physics
!!     set (z_fixed + closed faces, fv_mom6, energy Coriolis, Wright).
!! All are stratified with a boundary-layer scheme on, so the tiles exchange real
!! flow and real tracer structure.  26 x 18 cells, nghost = 3: every
!! factorisation above is uneven somewhere.
!!
!! The barotropic march-in (`&ocean_bt_nml bt_halo > 0`) is NOT covered: it
!! is opt-in precisely because it is not bit-identical to the serial run over
!! variable bathymetry or with open boundaries (see `resolve_bt_halo`).
!!
!! Each bug this test found, and the fix that made it pass, is named in the
!! commit that fixed it; a regression anywhere prints the first mismatching
!! cell of every mismatching field on every rank.
#ifdef RDB_ENABLE_MPI
program test_ocean_decomp_bitid_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, engine_exit_data, &
                               engine_teardown
   use rdb_ocean_state, only: ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use rdb_console_stats, only: console_stats_t, conservation_budget_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, ocean_budget_stage_weight
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
#ifndef RDB_NO_NETCDF
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_def_var_3d, nc_enddef, nc_put_var_2d, rdb_def_var_1d, &
                            rdb_put_var_1d
   use netcdf, only: nf90_put_var
#endif
   implicit none

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: NG = 3
   integer, parameter :: N_STEPS = 48
   real(wp), parameter :: DT = 900.0_wp
   integer, parameter :: MAXF = 64
   character(len=*), parameter :: SG_FILE = "bitid_supergrid.nc"
   character(len=*), parameter :: BATHY_FILE = "bitid_bathy.nc"
   character(len=*), parameter :: ZINIT_FILE = "bitid_zinit.nc"

   type :: field_t
      !! One compared field: a host copy of the whole local array (2-D
      !! fields are stored with a unit third extent) + its tag.
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snap_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
      integer :: io = 0, jo = 0, nxl = 0, nyl = 0
      type(conservation_budget_t) :: bud
         !! Closed-budget mass/salt/heat out+src totals at the final step —
         !! EFP-reduced (`reproducing_sums`, default on), so these must be
         !! bit-identical between the serial reference and every
         !! decomposition on a periodic configuration (see `compare_budget`).
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail, ic
   type(comm_t) :: comm
   character(len=16), parameter :: SCHEMES(2) = [character(len=16) :: "pred_corr", "ssp_rk2"]
   character(len=24), parameter :: CASES(8) = [character(len=24) :: &
                                               "island_basin", "periodic_channel_zstar", &
                                               "periodic_sponge", &
                                               "open_obc", "spherical", "obc_radiation_sponge", &
                                               "closures", "file_readers"]

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0
#ifndef RDB_NO_NETCDF
   if (rank == 0) call write_input_files()
   call comm%barrier()
#endif

   do ic = 1, size(CASES)
#ifdef RDB_NO_NETCDF
      if (trim(CASES(ic)) == "file_readers") cycle
#endif
      call run_case(trim(CASES(ic)), trim(SCHEMES(1)))
      call run_case(trim(CASES(ic)), trim(SCHEMES(2)))
   end do

   if (nprocs >= 2) call check_single_rank_fences()

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_decomp_bitid_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_decomp_bitid_mpi: all checks PASSED (nprocs=", &
            nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   function case_nml(label, scheme, px, py) result(nml)
      !! The namelist of one case on a px x py process grid.
      character(len=*), intent(in) :: label, scheme
      integer, intent(in) :: px, py
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: spx, spy, snx, sny
      character(len=:), allocatable :: common

      write (spx, '(i0)') px
      write (spy, '(i0)') py
      write (snx, '(i0)') NX_G
      write (sny, '(i0)') NY_G
      common = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               ""

      select case (label)
      case ("island_basin")
         ! Closed Cartesian basin, interior land block (static land mask),
         ! beta plane, eastward wind onto the island.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("periodic_channel_zstar")
         ! Re-entrant channel over a seamount, z* (ALE remap every step).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = -1.0e-4, wind_stress_x = 0.1 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_porous_nml enable = .true. /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("periodic_sponge")
         ! Re-entrant channel, periodic west/east, with a relaxing sponge
         ! band on the closed north edge -- the shape of the real Southern
         ! Ocean 1-degree configuration the periodic-seam / sponge-seam-ghost
         ! bug (this test's `bud` capture) was found on: a periodic axis
         ! carrying a physical-span-only relaxation on the orthogonal one.
         ! Legacy per-edge band sponge (not map-driven): simplest namelist
         ! that reaches `ocean_sponge_apply`/`ocean_sponge_apply_tracers`
         ! and the seam-ghost refresh gated on `sponge_seam` in
         ! `rdb_ocean_dyn.F90`.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.2 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'sponge', sponge_width = 3, sponge_strength = 1.0e-2, "// &
               "sponge_relax_tracers = .true. /"//NL
      case ("open_obc")
         ! Tidal west edge (eta target) + Flather open east edge.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.05, wind_stress_y = 0.02 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar_sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 200.0, slope_scale = 40000.0 /"//NL// &
               "&ocean_bc_nml west = 'tidal', east = 'open', south = 'wall', north = 'wall', "// &
               "west_n_tidal = 1, west_tidal_amp = 0.5, west_tidal_phase = 0.0, "// &
               "west_tidal_omega = 1.4051890e-4 /"//NL
      case ("spherical", "closures")
         ! Lon-lat sector, spoon basin, Wright EOS.  "closures" adds the
         ! lateral and vertical parameterisation set on top (EPBL instead of
         ! KPP, Fox-Kemper MLE, GM + MEKE, Redi, kappa-shear, tidal mixing,
         ! convective adjustment, geothermal heating, tracer hdiff).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 1.0, dy = 1.0 /"//NL// &
               "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, lat_south = 20.0, "// &
               "rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_topo_nml topo_config = 'spoon', max_depth = 3000.0, "// &
               "edge_depth = 300.0, slope_scale = 300000.0, wind_config = '2gyre', "// &
               "taux_magnitude = 0.1 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
         if (label == "closures") then
            nml = nml// &
                  "&ocean_vmix_nml use_kpp = .false. /"//NL// &
                  "&ocean_epbl_nml enable = .true. /"//NL// &
                  "&ocean_foxkemper_nml enable = .true. /"//NL// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_gm_nml enable = .true. /"//NL// &
                  "&ocean_meke_nml enable = .true. /"//NL// &
                  "&ocean_redi_nml enable = .true. /"//NL// &
                  "&ocean_kappa_shear_nml enable = .true. /"//NL// &
                  "&ocean_tidal_mixing_nml enable = .true., e_uniform = 1.0e-3 /"//NL// &
                  "&ocean_conv_nml enable = .true. /"//NL// &
                  "&ocean_geothermal_nml enable = .true. /"//NL// &
                  "&ocean_hdiff_nml kappa_h = 100.0 /"//NL
         end if
      case ("obc_radiation_sponge")
         ! Orlanski-radiating open south edge with tracer reservoirs, a
         ! clamped (inflow) north edge, and a relaxing sponge band west.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 500.0, slope_scale = 40000.0, wind_config = '2gyre', "// &
               "taux_magnitude = 0.05 /"//NL// &
               "&ocean_bc_nml west = 'sponge', east = 'wall', south = 'open', "// &
               "north = 'clamped', north_clamped_v = -0.01, sponge_width = 3, "// &
               "sponge_strength = 1.0e-4, sponge_relax_tracers = .true., "// &
               "radiation_scheme = 'orlanski', res_lscale_out = 20000.0, "// &
               "res_lscale_in = 20000.0 /"//NL
      case ("file_readers")
         ! The three per-rank windowed readers (supergrid, bathymetry,
         ! z-level T/S IC) on files the test writes, with the global-1-degree
         ! physics set: z_fixed + closed partial-step faces, fv_mom6 PGF,
         ! energy Coriolis, Wright EOS, periodic in x.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3 /"//NL// &
               "&ocean_grid_nml grid_config = 'supergrid', supergrid_file = '"//SG_FILE// &
               "', coriolis_scheme = 'planetary', rad_earth = 6.371e6 /"//NL// &
               "&physics_nml wind_stress_x = 0.08, wind_stress_y = 0.0 /"//NL// &
               "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true., "// &
               "check_vanished_content = .true. /"//NL// &
               "&ocean_topo_nml topo_config = 'file', max_depth = 3000.0 /"//NL// &
               "&output_nml bathymetry_file = '"//BATHY_FILE//"', output_to_file = .false. /"//NL// &
               "&ocean_zinit_nml enable = .true., source = 'file', file = '"//ZINIT_FILE//"' /"//NL// &
               "&ocean_pgf_nml form = 'fv_mom6' /"//NL// &
               "&ocean_coriolis_nml form = 'sadourny_energy' /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3, hbbl = 10.0, bg_vel = 0.1 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case default
         error stop "test_ocean_decomp_bitid_mpi: unknown case"
      end select
      if (label /= "file_readers") nml = nml//"&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine run_one(nml, csize, crank, snap, ok)
      !! Configure, step N_STEPS, snapshot every registry field to host.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: csize, crank
      type(snap_t), intent(out) :: snap
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(console_stats_t) :: cstats
      integer :: ierr, n
      real(wp) :: t

      ok = .false.
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr, compute_rank=crank, compute_size=csize)
      if (ierr /= OCEAN_STATUS_OK) return

      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      ierr = OCEAN_STATUS_OK
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         ! Closed-budget mass/salt/heat out+src at the final step, EFP-reduced
         ! (reproducing_sums default on) so it is order-invariant across rank
         ! counts -- the periodic-seam / tripolar-fold regression this test
         ! case exists for (see `compare_budget`).  A fresh `cstats` per call
         ! is fine: `budget_out` is computed independently of the t=0
         ! reference `cstats` latches internally.
         !
         ! HAZARD (found in review): `ocean_console_stats_report` takes no
         ! `compute_size` -- it is COLLECTIVE over whatever the REAL
         ! communicator is (`comm_env_setup_roles`), regardless of the
         ! logical `csize`/`crank` the caller is emulating.  The serial
         ! REFERENCE call below always passes `csize=1, crank=0` to
         ! `engine_setup` so every physical rank independently computes a
         ! full, non-decomposed copy of the SAME problem -- calling the
         ! collective report from inside that call would allreduce across
         ! every physical rank's independent copy (an `nprocs`-fold bogus
         ! sum for the mass/salt/heat totals feeding `bud`), not the single
         ! logical rank the reference is standing in for.  It is safe ONLY
         ! when the logical size this call is emulating equals the REAL
         ! physical size (`csize == nprocs`): always true for the decomposed
         ! `px x py` runs (`csize` is literally the real launch size there),
         ! and true for the reference only when the whole job is launched
         ! at `nprocs == 1` (ctest's np=1 leg), where the reference's own
         ! `csize=1` collective is trivially correct.  `snap%bud` stays its
         ! default (all-zero, `*_active=.false.`) when this is skipped;
         ! `run_case` only compares budgets it knows were computed safely.
         if (csize == nprocs) then
            call ocean_console_stats_report(cstats, engine%grid, engine%state%metrics, &
                                            engine%state%multilayer, t, DT, N_STEPS, &
                                            compute_rank=crank, reproducing_sums=.true., &
                                            budget_stage_weight=ocean_budget_stage_weight( &
                                            engine%state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR), &
                                            budget_out=snap%bud)
         end if
         call take_snapshot(engine, snap)
         ok = .true.
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_one

   subroutine take_snapshot(engine, s)
      !! Host copy of every restart-registry field (device -> host first)
      !! plus the BT end-of-step eta.
      type(ocean_engine_t), intent(inout), target :: engine
      type(snap_t), intent(inout) :: s
      type(restart_registry_t) :: reg
      integer :: e

      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      s%n = 0
      do e = 1, reg%n
         associate (en => reg%entries(e))
            if (en%rank == 0) cycle     ! rank-local host scalars (Chapman corners)
            if (en%device_mapped) then
               if (en%rank == 2) then
                  !$acc update self(en%p2)
               else
                  !$acc update self(en%p3)
               end if
            end if
            s%n = s%n + 1
            if (s%n > MAXF) error stop "test_ocean_decomp_bitid_mpi: raise MAXF"
            s%f(s%n)%tag = en%tag
            if (en%rank == 2) then
               allocate (s%f(s%n)%a(size(en%p2, 1), size(en%p2, 2), 1))
               s%f(s%n)%a(:, :, 1) = en%p2
            else
               s%f(s%n)%a = en%p3
            end if
         end associate
      end do
      !$acc update self(engine%state%dyn%bt_work%bt_eta_end)
      s%n = s%n + 1
      s%f(s%n)%tag = "bt_eta_end"
      allocate (s%f(s%n)%a(size(engine%state%dyn%bt_work%bt_eta_end, 1), &
                           size(engine%state%dyn%bt_work%bt_eta_end, 2), 1))
      s%f(s%n)%a(:, :, 1) = engine%state%dyn%bt_work%bt_eta_end
      s%io = engine%grid%i_offset_global
      s%jo = engine%grid%j_offset_global
      s%nxl = engine%grid%nx_phys
      s%nyl = engine%grid%ny_phys
   end subroutine take_snapshot

   subroutine compare(dec, ref, nbad, nfin, report)
      !! Bitwise comparison of each field's owned window.  A face array is
      !! one wider than the cell array along its stagger; its owned window
      !! then includes both edge faces of the tile.
      type(snap_t), intent(in) :: dec, ref
      integer, intent(out) :: nbad, nfin
      character(len=*), intent(in) :: report
      integer :: e, i, j, k, ex, ey, nb
      character(len=160) :: first

      nbad = 0
      nfin = 0
      if (dec%n /= ref%n) then
         write (*, '(a,i0,a,i0,a,i0)') "  rank ", rank, ": field count differs ", dec%n, " vs ", ref%n
         nbad = 1
         return
      end if
      do e = 1, dec%n
         if (dec%f(e)%tag /= ref%f(e)%tag) then
            write (*, '(a,i0,4a)') "  rank ", rank, ": field order differs ", &
               trim(dec%f(e)%tag), " vs ", trim(ref%f(e)%tag)
            nbad = nbad + 1
            cycle
         end if
         associate (a => dec%f(e)%a, b => ref%f(e)%a)
            ex = size(a, 1) - (dec%nxl + 2*NG)
            ey = size(a, 2) - (dec%nyl + 2*NG)
            nb = 0
            first = ""
            do k = 1, size(a, 3)
               do j = NG + 1, NG + dec%nyl + ey
                  do i = NG + 1, NG + dec%nxl + ex
                     if (.not. ieee_is_finite(b(i + dec%io, j + dec%jo, k))) nfin = nfin + 1
                     if (transfer(a(i, j, k), 0_int64) /= &
                         transfer(b(i + dec%io, j + dec%jo, k), 0_int64)) then
                        if (nb == 0) write (first, '(a,3(i0,1x),a,es24.16,a,es24.16)') &
                           "first at local (i,j,k)=", i, j, k, ": ", a(i, j, k), " vs ", &
                           b(i + dec%io, j + dec%jo, k)
                        nb = nb + 1
                     end if
                  end do
               end do
            end do
         end associate
         if (nb > 0) then
            write (*, '(a,i0,5a,i0,2a)') "  rank ", rank, " ", report, " ", &
               trim(dec%f(e)%tag), ": mismatches=", nb, "  ", trim(first)
         end if
         nbad = nbad + nb
      end do
   end subroutine compare

   subroutine run_case(label, scheme)
      character(len=*), intent(in) :: label, scheme
      type(snap_t) :: ref, dec
      logical :: ok_ref, ok_dec
      integer :: px, nbad, nfin, glob(3), nbud
      character(len=96) :: tag
      type(conservation_budget_t) :: buds(8)

      call run_one(case_nml(label, scheme, 1, 1), 1, 0, ref, ok_ref)
      if (.not. ok_ref) then
         write (*, '(5a,i0)') "FAIL ", label, "/", scheme, ": serial reference failed on rank ", rank
         n_fail = n_fail + 1
         return
      end if
      ! Real single-rank launch (ctest np=1): `ref`'s budget was computed
      ! safely (csize == nprocs == 1, see the hazard note in `run_one`) and
      ! IS the historical bug's exact scenario -- a genuinely undecomposed
      ! run, no MPI seam anywhere, so a periodic axis's local wrap is the
      ! ONLY thing that can keep the sponge-relaxed seam consistent.  This
      ! is the one place in THIS binary that can re-detect the original
      ! bug (see `compare_budget`'s docstring for why the decomposed
      ! factorizations below cannot).
      if (nprocs == 1) call check_periodic_sponge_serial_out(label, scheme, ref%bud)

      nbud = 0
      do px = nprocs, 1, -1
         if (mod(nprocs, px) /= 0) cycle
         write (tag, '(4a,i0,a,i0)') label, "/", scheme, " ", px, "x", nprocs/px
         call run_one(case_nml(label, scheme, px, nprocs/px), nprocs, rank, dec, ok_dec)
         if (.not. ok_dec) then
            write (*, '(3a,i0)') "FAIL ", trim(tag), ": decomposed run failed on rank ", rank
            glob = [1, 0, 1]
         else
            call compare(dec, ref, nbad, nfin, trim(tag))
            glob = [nbad, nfin, 0]
            if (nbud < size(buds)) then
               nbud = nbud + 1
               buds(nbud) = dec%bud
            end if
         end if
         call allreduce(comm, glob, op=MPI_SUM)
         if (rank == 0) then
            if (glob(1) == 0 .and. glob(3) == 0 .and. glob(2) == 0) then
               write (*, '(3a,i0,a)') "case ", trim(tag), ": IDENTICAL (", ref%n, " fields)"
            else
               write (*, '(3a,i0,a,i0,a,i0)') "FAIL ", trim(tag), ": mismatches=", glob(1), &
                  " nonfinite_ref=", glob(2), " run_failures=", glob(3)
            end if
         end if
         if (sum(glob) > 0) n_fail = n_fail + 1
      end do
      ! Every collected `dec%bud` above has `csize == nprocs` (always true
      ! for the decomposed loop), so its collective was safe regardless of
      ! `nprocs`; comparing them to each other is meaningful once nprocs > 1
      ! (>= 2 distinct factorizations).
      if (nbud >= 2) call compare_budget(buds(1:nbud), label, scheme)
   end subroutine run_case

   subroutine compare_budget(buds, label, scheme)
      !! Bit-identity check across every decomposed-factorization budget of
      !! one case/scheme collected by `run_case` (2x1 vs 1x2 for nprocs=2;
      !! 4x1 vs 2x2 vs 1x4 for nprocs=4) -- every entry has `csize ==
      !! nprocs` so its `ocean_console_stats_report` collective was taken
      !! over the REAL communicator correctly (see the hazard note in
      !! `run_one`); no unsafe reference-vs-decomposed comparison here.
      !! The EFP reproducing-sums path (default on) makes the closed
      !! mass/salt/heat out+src totals order-invariant across rank counts,
      !! so any two factorizations of the SAME real launch must land on the
      !! identical bit pattern.
      !!
      !! NOTE what this does and does not catch: every valid `px x py`
      !! factorization of nprocs >= 2 has `px > 1 .or. py > 1`, so the OLD
      !! buggy gate (`ocean_halo_is_decomposed_x() .or. ..._y()`) already
      !! fired the seam-ghost refresh for EVERY entry compared here, both
      !! before and after the fix -- a mismatch here is a regression in
      !! the general refresh mechanism, not a re-detection of the original
      !! single-rank bug.  Only a genuinely undecomposed run (nprocs == 1,
      !! `check_periodic_sponge_serial_out` below) exercises the code path
      !! the fix changed.
      type(conservation_budget_t), intent(in) :: buds(:)
      character(len=*), intent(in) :: label, scheme
      integer :: k, nbad

      nbad = 0
      do k = 2, size(buds)
         if (transfer(buds(k)%mass_out, 0_int64) /= transfer(buds(1)%mass_out, 0_int64) .or. &
             transfer(buds(k)%salt_out, 0_int64) /= transfer(buds(1)%salt_out, 0_int64) .or. &
             transfer(buds(k)%heat_out, 0_int64) /= transfer(buds(1)%heat_out, 0_int64) .or. &
             transfer(buds(k)%mass_src, 0_int64) /= transfer(buds(1)%mass_src, 0_int64) .or. &
             transfer(buds(k)%salt_src, 0_int64) /= transfer(buds(1)%salt_src, 0_int64) .or. &
             transfer(buds(k)%heat_src, 0_int64) /= transfer(buds(1)%heat_src, 0_int64)) then
            nbad = nbad + 1
         end if
      end do
      if (nbad > 0) then
         n_fail = n_fail + nbad
         if (rank == 0) write (*, '(5a,i0,a)') "FAIL budget ", trim(label), "/", trim(scheme), &
            ": ", nbad, " decomposed factorization(s) disagree on the closed mass/salt/heat budget"
      else if (rank == 0) then
         write (*, '(4a)') "case ", trim(label), "/", trim(scheme), &
            " budget: IDENTICAL across factorizations"
      end if
   end subroutine compare_budget

   subroutine check_periodic_sponge_serial_out(label, scheme, bud)
      !! The real single-rank regression for the periodic-seam / sponge bug:
      !! on "periodic_sponge" (periodic west/east, a relaxing sponge band
      !! on the closed north edge) at a genuinely undecomposed nprocs == 1
      !! launch, the closed-budget `mass_out` must stay near a WALL-EDGE
      !! CONTROL's own -- the same namelist with west/east = "wall" instead
      !! of "periodic", so there is no periodic seam at all and the old
      !! gate (`ocean_halo_is_decomposed_x() .or. ..._y()`, both false on a
      !! real single rank) was never wrong there.  The bug left the
      !! periodic-seam ghost columns un-relaxed after the sponge touched
      !! only the tile's physical cells, so the corrector's advection saw a
      !! spurious mass/salt jump at the seam that the control cannot have
      !! (the real Southern Ocean 1-degree case: Salt out 1.6e8 on 1 rank
      !! vs -81 on 2 ranks -- this toy grid/duration reaches nothing like
      !! that magnitude, but `mass_out` still separates monotonically:
      !! broken is consistently LARGER than fixed on both split schemes,
      !! see the commit message).  `mass_out` is used rather than
      !! `salt_out` because the latter is noisy at this scale under
      !! `pred_corr` (its own resting-state round-off floor, see CLAUDE.md,
      !! dominates before the seam term does); `mass_out` is monotonic on
      !! both schemes. `sponge_width` is deliberately narrow (3 cells,
      !! `<=` the smallest per-rank tile at any tested `px x py`): a wider
      !! band that crosses an uneven y-decomposition tile boundary (found
      !! at `sponge_width = 6`, `1x4`: the legacy band sponge does not
      !! handle a band split across ranks) is a SEPARATE, pre-existing
      !! defect this case must not also exercise.
      !!
      !! Thresholds are PER-SCHEME (`ssp_rk2`'s own noise floor sits an
      !! order above `pred_corr`'s, CLAUDE.md's split-scheme section) and
      !! are pinned directly to measured before/after `mass_out` on this
      !! exact namelist (gfortran, default flags):
      !!   pred_corr: fixed 9.2255E-05, broken 1.5045E-04 (1.6x)
      !!   ssp_rk2:   fixed 8.3477E-06, broken 1.6591E-04 (19.9x)
      !! Each threshold sits above its scheme's fixed value and below its
      !! broken value.
      character(len=*), intent(in) :: label, scheme
      type(conservation_budget_t), intent(in) :: bud
      real(wp) :: tol

      if (trim(label) /= "periodic_sponge") return
      if (.not. bud%mass_active) return

      select case (trim(scheme))
      case ("pred_corr")
         tol = 1.2e-4_wp
      case ("ssp_rk2")
         tol = 4.0e-5_wp
      case default
         return
      end select

      if (abs(bud%mass_out) > tol) then
         n_fail = n_fail + 1
         write (*, '(5a,es14.6,a,es14.6)') "FAIL serial-out ", trim(label), "/", trim(scheme), &
            ": |mass_out|=", abs(bud%mass_out), " exceeds ", tol
      else if (rank == 0) then
         write (*, '(4a)') "case ", trim(label), "/", trim(scheme), &
            " serial (1 rank) out: below the fixed/broken separation bound (ok)"
      end if
   end subroutine check_periodic_sponge_serial_out

   subroutine check_single_rank_fences()
      !! The single-rank features must be REFUSED at configure on more than
      !! one rank -- with the process grid left unset (px = py = 1, which
      !! the engine auto-factors and which `validate_config`'s px*py fences
      !! cannot see), so the engine-side gate is the one exercised.
      character(len=*), parameter :: NL = new_line("a")
      character(len=48), parameter :: KNOBS(4) = [character(len=48) :: &
                                                  "&ocean_wetdry_nml enable = .true. /", &
                                                  "&ocean_ice_nml enable = .true. /", &
                                                  "&ocean_cavity_dyn_nml enable = .true. /", &
                                                  "&ocean_bc_nml east = 'chapman' /"]
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ik, ierr

      do ik = 1, size(KNOBS)
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&grid_nml nx = 26, ny = 18, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'ssp_rk2' /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&output_nml output_to_file = .false. /"//NL// &
               trim(KNOBS(ik))//NL
         call read_config_from_string(nml, cfg, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call validate_config(cfg, ierr)
         if (ierr == OCEAN_STATUS_OK) then
            call engine_setup(engine, cfg, ierr, compute_rank=rank, compute_size=nprocs)
            call engine_teardown(engine)
         end if
         if (ierr == OCEAN_STATUS_OK) then
            write (*, '(3a,i0)') "FAIL fence: '", trim(KNOBS(ik)), &
               "' was ACCEPTED on multi-rank, rank ", rank
            n_fail = n_fail + 1
         else if (rank == 0) then
            write (*, '(3a)') "case fence: '", trim(KNOBS(ik)), "' refused on multi-rank (ok)"
         end if
      end do
   end subroutine check_single_rank_fences

#ifndef RDB_NO_NETCDF
   subroutine write_input_files()
      !! Rank 0 writes the three whole-grid input files the "file_readers"
      !! case reads through the per-rank windowed readers: a periodic lon-lat
      !! MOM6 supergrid (360 deg in x, 20-56 N), a bathymetry with a land
      !! block stored C-ORDER (Fortran dims (y, x) — the reader's transpose
      !! path) and a z-level T/S initial condition (x, y, z).  Every field
      !! varies in both directions, so a tile that read the wrong window
      !! cannot match the serial run.
      integer, parameter :: NZS = 8
      real(wp), parameter :: DEG = 3.14159265358979323846_wp/180.0_wp
      real(wp), parameter :: RE = 6.371e6_wp, LAT0 = 20.0_wp, DLAT = 2.0_wp
      real(wp) :: dlon
      integer :: ncid, dnxp, dnyp, dnx, dny, dx_, dy_, dz_, v1, v2, v3, v4, v5, ierr
      integer :: m, n, i, j, k
      real(wp), allocatable :: sx(:, :), sy(:, :), sdx(:, :), sdy(:, :), sar(:, :)
      real(wp), allocatable :: byx(:, :), tt(:, :, :), ss(:, :, :)
      real(wp) :: zs(NZS), lat, x, y

      dlon = 360.0_wp/real(NX_G, wp)
      allocate (sx(2*NX_G + 1, 2*NY_G + 1), sy(2*NX_G + 1, 2*NY_G + 1))
      allocate (sdx(2*NX_G, 2*NY_G + 1), sdy(2*NX_G + 1, 2*NY_G), sar(2*NX_G, 2*NY_G))
      do n = 1, 2*NY_G + 1
         lat = LAT0 + real(n - 1, wp)*0.5_wp*DLAT
         do m = 1, 2*NX_G + 1
            sx(m, n) = real(m - 1, wp)*0.5_wp*dlon
            sy(m, n) = lat
         end do
         do m = 1, 2*NX_G
            sdx(m, n) = RE*cos(lat*DEG)*0.5_wp*dlon*DEG
         end do
      end do
      sdy = RE*0.5_wp*DLAT*DEG
      do n = 1, 2*NY_G
         do m = 1, 2*NX_G
            sar(m, n) = sdx(m, n)*sdy(m, n)
         end do
      end do
      call nc_create_file(SG_FILE, ncid)
      call nc_def_dim(ncid, "nxp", 2*NX_G + 1, dnxp)
      call nc_def_dim(ncid, "nyp", 2*NY_G + 1, dnyp)
      call nc_def_dim(ncid, "nx", 2*NX_G, dnx)
      call nc_def_dim(ncid, "ny", 2*NY_G, dny)
      call nc_def_var_2d(ncid, "x", [dnxp, dnyp], v1)
      call nc_def_var_2d(ncid, "y", [dnxp, dnyp], v2)
      call nc_def_var_2d(ncid, "dx", [dnx, dnyp], v3)
      call nc_def_var_2d(ncid, "dy", [dnxp, dny], v4)
      call nc_def_var_2d(ncid, "area", [dnx, dny], v5)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, v1, sx)
      call nc_put_var_2d(ncid, v2, sy)
      call nc_put_var_2d(ncid, v3, sdx)
      call nc_put_var_2d(ncid, v4, sdy)
      call nc_put_var_2d(ncid, v5, sar)
      call nc_close(ncid)

      ! Bathymetry, positive-down, C-order: byx(j, i).  A shelf rising to
      ! the north, a Gaussian ridge, and a small land block.
      allocate (byx(NY_G, NX_G))
      do i = 1, NX_G
         do j = 1, NY_G
            x = real(i, wp)/real(NX_G, wp)
            y = real(j, wp)/real(NY_G, wp)
            byx(j, i) = 3000.0_wp - 1800.0_wp*y**2 - &
                        900.0_wp*exp(-((x - 0.3_wp)**2 + (y - 0.5_wp)**2)/0.02_wp)
            if (i >= 17 .and. i <= 19 .and. j >= 7 .and. j <= 10) byx(j, i) = 0.0_wp
         end do
      end do
      call nc_create_file(BATHY_FILE, ncid)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_var_2d(ncid, "depth", [dy_, dx_], v1)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, v1, byx)
      call nc_close(ncid)

      ! z-level T/S (x, y, z), positive-down source depths.
      allocate (tt(NX_G, NY_G, NZS), ss(NX_G, NY_G, NZS))
      do k = 1, NZS
         zs(k) = real(k - 1, wp)*450.0_wp
      end do
      do k = 1, NZS
         do j = 1, NY_G
            do i = 1, NX_G
               x = real(i, wp)/real(NX_G, wp)
               y = real(j, wp)/real(NY_G, wp)
               tt(i, j, k) = 22.0_wp - 18.0_wp*y - zs(k)/250.0_wp + 1.5_wp*sin(6.2831853_wp*x)
               tt(i, j, k) = max(tt(i, j, k), -1.5_wp)
               ss(i, j, k) = 34.5_wp + 0.6_wp*y - 0.2_wp*cos(6.2831853_wp*x) + zs(k)*1.0e-4_wp
            end do
         end do
      end do
      call nc_create_file(ZINIT_FILE, ncid)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "z", NZS, dz_)
      call nc_def_var_3d(ncid, "temp", [dx_, dy_, dz_], v1)
      call nc_def_var_3d(ncid, "salt", [dx_, dy_, dz_], v2)
      call rdb_def_var_1d(ncid, "z_src", dz_, v3)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, v1, tt)
      ierr = nf90_put_var(ncid, v2, ss)
      call rdb_put_var_1d(ncid, v3, zs)
      call nc_close(ncid)
   end subroutine write_input_files
#endif

end program test_ocean_decomp_bitid_mpi
#else
program test_ocean_decomp_bitid_mpi
   implicit none
   write (*, '(a)') "test_ocean_decomp_bitid_mpi: skipped (RDB_ENABLE_MPI=OFF)"
end program test_ocean_decomp_bitid_mpi
#endif
