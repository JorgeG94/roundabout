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
!!     (ALE remap every step);
!!   * open_obc — tidal west edge (eta target) + Flather east edge,
!!     zstar_sigma, over a seamount;
!!   * spherical — lon-lat sector, planetary Coriolis, spoon basin, Wright
!!     EOS, `2gyre` wind;
!!   * obc_radiation_sponge — Orlanski-radiating open south edge with tracer
!!     reservoirs, clamped north inflow, legacy relaxing sponge band west.
!! All are stratified with KPP on (the default), so the tiles exchange real
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
   implicit none

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: NG = 3
   integer, parameter :: N_STEPS = 48
   real(wp), parameter :: DT = 900.0_wp
   integer, parameter :: MAXF = 64

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
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail, ic
   type(comm_t) :: comm
   character(len=16), parameter :: SCHEMES(2) = [character(len=16) :: "pred_corr", "ssp_rk2"]
   character(len=24), parameter :: CASES(5) = [character(len=24) :: &
                                               "island_basin", "periodic_channel_zstar", &
                                               "open_obc", "spherical", "obc_radiation_sponge"]

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0

   do ic = 1, size(CASES)
      call run_case(trim(CASES(ic)), trim(SCHEMES(1)))
      call run_case(trim(CASES(ic)), trim(SCHEMES(2)))
   end do

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
               "&output_nml output_to_file = .false. /"//NL

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
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
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
      case ("spherical")
         ! Lon-lat sector, spoon basin, Wright EOS.
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
      case default
         error stop "test_ocean_decomp_bitid_mpi: unknown case"
      end select
   end function case_nml

   subroutine run_one(nml, csize, crank, snap, ok)
      !! Configure, step N_STEPS, snapshot every registry field to host.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: csize, crank
      type(snap_t), intent(out) :: snap
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
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
      integer :: px, nbad, nfin, glob(3)
      character(len=96) :: tag

      call run_one(case_nml(label, scheme, 1, 1), 1, 0, ref, ok_ref)
      if (.not. ok_ref) then
         write (*, '(5a,i0)') "FAIL ", label, "/", scheme, ": serial reference failed on rank ", rank
         n_fail = n_fail + 1
         return
      end if

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
   end subroutine run_case

end program test_ocean_decomp_bitid_mpi
#else
program test_ocean_decomp_bitid_mpi
   implicit none
   write (*, '(a)') "test_ocean_decomp_bitid_mpi: skipped (RDB_ENABLE_MPI=OFF)"
end program test_ocean_decomp_bitid_mpi
#endif
