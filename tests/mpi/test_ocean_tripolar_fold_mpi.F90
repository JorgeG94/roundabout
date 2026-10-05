!! Tripolar north fold under MPI: 1-rank vs decomposed bit-identity through
!! the production engine.
!!
!! The fold is applied by the ranks that own the north edge (`bc%north_fold`
!! is rank-local).  Before that gate existed, every rank of a px = 1, py > 1
!! run folded its OWN north ghosts — overwriting, with a mirror of its own
!! tile, the rows the MPI exchange had just filled — and the analytic
!! tripolar generator built every tile as a whole globe of the tile's size.
!! Neither failed loud; the answers were simply wrong.
!!
!! On EACH launch (ctest runs np = 1, 2, 3, 4) the binary:
!!   1. runs the SERIAL reference (px = py = 1, compute_size = 1) on every
!!      rank — deterministic, no collective that could differ — through
!!      `engine_setup` -> seeded structure -> `engine_step` x N_STEPS;
!!   2. runs every DECOMPOSED case of the launched rank count on the same
!!      namelist: the north-south split (px = 1, py = nprocs) and the
!!      east-west splits through the distributed fold exchange (2x1; 3x1;
!!      4x1 and 2x2), the latter via `engine_setup`'s TEST-ONLY
!!      `allow_distributed_fold` until the px > 1 refusal is lifted;
!!   3. compares, BITWISE, each rank's whole storage window (ghost rows and
!!      columns included) of h, u, v, every tracer hTr (S, T) and the BT
!!      end-of-step eta against the matching window of the reference, plus
!!      the tripolar metrics at configure time (the generator gate);
!!   4. on np >= 2, asserts that px = 2 WITHOUT the bypass is still REFUSED
!!      at configure time.
!!
!! Two grids:
!!   * "spanning": ny = 24 from 59N at 1 deg, phi_join = 74 — the bipolar
!!     cap (rows 16-24) crosses the 4-rank seam at row 18 (an even split:
!!     12/12 on 2 ranks, 6/6/6/6 on 4);
!!   * "north_tile": ny = 23 from 47N at 1.5 deg, phi_join = 76 — the cap
!!     lies inside the north rank's tile, and the split is uneven (the
!!     remainder goes to the southern ranks, so the north tile is the
!!     SHORTEST: 12/11, 6/6/6/5).
!! Both caps are kept at >= 73N on purpose: a larger analytic cap on this
!! 32-column grid goes non-finite within two steps even single-rank (on the
!! V100 build from ~72N, on gfortran from ~66N) — a property of the analytic
!! generator's pole meridian, not of the decomposition — and a NaN
!! trajectory is trivially "bit-identical", so the reference is also
!! required to stay finite.
!!
!! The seeded state (after `engine_setup`, before the device map) puts a
!! structured u, v and T/S anomaly across the whole grid — the fold row
!! included — then refills the ghosts the way `engine_setup` does (halo
!! exchange, then the fold), so the first step starts from a
!! fold-consistent state on every decomposition.
#ifdef RDB_ENABLE_MPI
program test_ocean_tripolar_fold_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, engine_exit_data, &
                               engine_teardown
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_state
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   integer, parameter :: NX_G = 32
   integer, parameter :: NZ = 3
   integer, parameter :: N_STEPS = 24
   real(wp), parameter :: DT = 1800.0_wp

   type :: snap_t
      !! Host copy of one run's fields (whole local storage).
      real(wp), allocatable :: h(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: tr(:, :, :, :)
      real(wp), allocatable :: eta(:, :)
      real(wp), allocatable :: dxt(:, :), areat(:, :), dycu(:, :), dxcv(:, :)
      real(wp), allocatable :: dxbu(:, :), areabu(:, :), geolatbu(:, :), angle(:, :)
      real(wp), allocatable :: fcorner(:, :)
      integer :: io = 0
         !! Global i offset of the tile.
      integer :: jo = 0
         !! Global j offset of the tile.
      logical :: has_south = .true.
      logical :: has_north = .true.
      integer :: px = 1
         !! Tiles along x of the run that produced the snapshot.
      integer :: ng = 0
         !! `nghost` of the run that produced the snapshot.
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail
   logical :: ok_ref_g
   type(comm_t) :: comm

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0

   call run_case("spanning", 24, 59.0_wp, 1.0_wp, 74.0_wp)
   call run_case("north_tile", 23, 47.0_wp, 1.5_wp, 76.0_wp)
   if (nprocs >= 2) call check_px_refused()

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_tripolar_fold_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_tripolar_fold_mpi: all checks PASSED (nprocs=", &
            nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   function case_nml(ny, lat_south, dlat, phi_join, px, py) result(nml)
      !! Small analytic tripolar ocean: periodic x, wall south, fold north,
      !! flat 1000 m bed, 3 layers, planetary Coriolis, the production split.
      integer, intent(in) :: ny, px, py
      real(wp), intent(in) :: lat_south, dlat, phi_join
      character(len=:), allocatable :: nml
      character(len=32) :: sny, sphi, spx, spy, slat, sdlat
      character(len=*), parameter :: NL = new_line("a")
      write (sny, '(i0)') ny
      write (sphi, '(f8.3)') phi_join
      write (slat, '(f8.3)') lat_south
      write (sdlat, '(f8.3)') dlat
      write (spx, '(i0)') px
      write (spy, '(i0)') py
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 32, ny = "//trim(sny)//", nghost = 3, "// &
            "dx = 11.25, dy = "//trim(adjustl(sdlat))//" /"//NL// &
            "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
            "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, "// &
            "lat_south = "//trim(adjustl(slat))//", phi_join = "//trim(adjustl(sphi))//", lon_pole = 0.0, "// &
            "rad_earth = 6.378e6, coriolis_scheme = 'planetary' /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 1800.0 /"//NL// &
            "&ocean_topo_nml max_depth = 1000.0 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'tripolar_fold' /"//NL// &
            "&ocean_hvisc_nml nu_h = 2000.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false., reproducing_sums = .true. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine setup_engine(engine, nml, csize, crank, ierr)
      !! Parse + validate + `engine_setup` (the production configure path).
      type(ocean_engine_t), intent(inout) :: engine
      character(len=*), intent(in) :: nml
      integer, intent(in) :: csize, crank
      integer, intent(out) :: ierr
      type(config_t) :: cfg
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr, compute_rank=crank, compute_size=csize)
   end subroutine setup_engine

   subroutine run_one(nml, csize, crank, snap_cfg, snap_end, ok, allow_dfold)
      !! Configure, seed the structure, step N_STEPS, snapshot to host.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: csize, crank
      logical, intent(in), optional :: allow_dfold
         !! Pass engine_setup's test-only `allow_distributed_fold` (px > 1).
      type(snap_t), intent(out) :: snap_cfg, snap_end
      logical, intent(out) :: ok
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      integer :: ierr, n
      real(wp) :: t

      ok = .false.
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr, compute_rank=crank, compute_size=csize, &
                        allow_distributed_fold=allow_dfold)
      if (ierr /= OCEAN_STATUS_OK) return

      call seed_structure(engine)
      call take_snapshot(engine, snap_cfg, .false.)

      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         call take_snapshot(engine, snap_end, .true.)
         ok = .true.
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_one

   subroutine seed_structure(engine)
      !! A structured, decomposition-invariant anomaly (functions of the
      !! GLOBAL indices only) in u, v and every tracer, reaching the fold row;
      !! then the ghosts are refilled the way `engine_setup` does it.  Exact
      !! rational values, NOT transcendentals: a vectorised libm `sin`/`cos`
      !! (gfortran -O3 -march=native -> libmvec) differs from the scalar one
      !! in the last bit, and which cells take which path depends on the
      !! tile's row length — a seed artifact that differs by decomposition.
      !! u uses the REDUCED face index, so face NX_G+1 equals face 1 bitwise
      !! (the periodic-seam invariant every decomposition relies on).
      type(ocean_engine_t), intent(inout) :: engine
      integer :: i, j, k, it, ng, nxl, nyl, io, jo, ig, jg

      ng = engine%grid%nghost
      nxl = engine%grid%nx_phys
      nyl = engine%grid%ny_phys
      io = engine%grid%i_offset_global
      jo = engine%grid%j_offset_global
      associate (ms => engine%state%multilayer)
         do k = 1, NZ
            do j = ng + 1, ng + nyl
               do i = ng + 1, ng + nxl + 1
                  ig = modulo(i - ng + io - 1, NX_G) + 1
                  jg = j - ng + jo
                  ms%u_face_x_layer(i, j, k) = real(modulo(5*ig + 3*jg + k, 13) - 6, wp)/128.0_wp
               end do
            end do
            ! v: every local south face plus the tile's north face — on the
            ! north rank that last row is the fold line itself.
            do j = ng + 1, ng + nyl + 1
               do i = ng + 1, ng + nxl
                  ig = i - ng + io
                  jg = j - ng + jo
                  ms%v_face_y_layer(i, j, k) = real(modulo(7*ig + 2*jg + 3*k, 11) - 5, wp)/256.0_wp
               end do
            end do
            do it = 1, size(ms%tracers)
               do j = ng + 1, ng + nyl
                  do i = ng + 1, ng + nxl
                     ig = i - ng + io
                     jg = j - ng + jo
                     ms%tracers(it)%hTr(i, j, k) = ms%tracers(it)%hTr(i, j, k)* &
                                                   (1.0_wp + real(modulo(3*ig + jg + it, 7) - 3, wp)/1024.0_wp)
                  end do
               end do
            end do
         end do
         call ocean_halo_exchange_ml_state(ms, device_resident=.false.)
         ! Host-only state (before `engine_enter_data`): the distributed
         ! path's pack/unpack must run on the host copies.
         call ocean_fold_wrap_state(engine%grid, engine%state%bc, ms, device_resident=.false.)
      end associate
   end subroutine seed_structure

   subroutine take_snapshot(engine, s, after_steps)
      !! Host copy of the compared fields (device -> host first when mapped).
      type(ocean_engine_t), intent(inout) :: engine
      type(snap_t), intent(out) :: s
      logical, intent(in) :: after_steps
      integer :: it, nt

      s%io = engine%grid%i_offset_global
      s%jo = engine%grid%j_offset_global
      s%has_south = engine%decomp%has_south
      s%has_north = engine%decomp%has_north
      s%px = engine%decomp%px
      s%ng = engine%grid%nghost
      associate (ms => engine%state%multilayer, mt => engine%state%metrics)
         if (after_steps) then
            !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
            do it = 1, size(ms%tracers)
               !$acc update self(ms%tracers(it)%hTr)
            end do
            !$acc update self(engine%state%dyn%bt_work%bt_eta_end)
         end if
         s%h = ms%h_layer
         s%u = ms%u_face_x_layer
         s%v = ms%v_face_y_layer
         nt = size(ms%tracers)
         allocate (s%tr(size(ms%h_layer, 1), size(ms%h_layer, 2), NZ, nt))
         do it = 1, nt
            s%tr(:, :, :, it) = ms%tracers(it)%hTr
         end do
         if (after_steps) s%eta = engine%state%dyn%bt_work%bt_eta_end
         s%dxt = mt%dxT
         s%areat = mt%areaT
         s%dycu = mt%dyCu
         s%dxcv = mt%dxCv
         s%dxbu = mt%dxBu
         s%areabu = mt%areaBu
         s%geolatbu = mt%geolatBu
         s%angle = mt%angle_dx
         s%fcorner = engine%state%coriolis_adv%f_corner
      end associate
   end subroutine take_snapshot

   subroutine compare_2d(name, tile, whole, d, is_vface, nbad_phys, nbad_ghost, is_uface)
      !! Bitwise tile-vs-reference-window comparison.  The window is
      !! `whole(i + io, j + jo)`; a mismatch is counted as PHYSICAL when it
      !! lies in the tile's own rows (for a v array, its south faces plus its
      !! north face — on the north tile, the fold row), otherwise as a GHOST
      !! mismatch.  BOTH fail the test.  One row is exempt: the OUTERMOST v
      !! row beyond an MPI y seam (`j = 1` above a south seam,
      !! `j = ny_total+1` below a north seam).  The face-y halo exchanges
      !! `nghost` rows, and a south-face v array has `nghost + 1` rows past
      !! the last owned face on the north side and `nghost` + the array-edge
      !! row on the south side: that one row is outside every stencil and is
      !! never refreshed on ANY decomposed run (fold or not), so it holds
      !! whatever the previous stage left there.  Its x twin is exempt too:
      !! the OUTERMOST u column (`i = 1`, `i = nx_total+1`) of an x-split run
      !! (px > 1; the periodic wrap link is an MPI seam as well).  The ALE
      !! face remap (`remap_x_face_velocity`, `I = 1:nx+1`) gives an
      !! array-edge face the adjacent cell's thickness, one-sided, AFTER the
      !! stage-end exchange; on the tile that column is the neighbour's
      !! interior face, on the serial run it is an interior face remapped
      !! two-sided, and it is refreshed at the next stage entry before any
      !! owned value reads it.
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: tile(:, :), whole(:, :)
      type(snap_t), intent(in) :: d
         !! The decomposed snapshot (tile offsets and edge flags).
      logical, intent(in) :: is_vface
      integer, intent(inout) :: nbad_phys, nbad_ghost
      logical, intent(in), optional :: is_uface
         !! A u (x-face) array: apply the outermost-column exemption.
      integer :: i, j, jhi, nb_p, nb_g, ng, nyl
      logical :: phys, uface

      uface = .false.
      if (present(is_uface)) uface = is_uface
      ng = d%ng
      nyl = size(tile, 2) - 2*ng
      if (is_vface) nyl = nyl - 1
      jhi = ng + nyl
      if (is_vface) jhi = jhi + 1
      nb_p = 0
      nb_g = 0
      do j = 1, size(tile, 2)
         phys = (j >= ng + 1 .and. j <= jhi)
         if (is_vface .and. j == 1 .and. .not. d%has_south) cycle
         if (is_vface .and. j == size(tile, 2) .and. .not. d%has_north) cycle
         do i = 1, size(tile, 1)
            if (uface .and. d%px > 1 .and. (i == 1 .or. i == size(tile, 1))) cycle
            if (transfer(tile(i, j), 0_int64) /= transfer(whole(i + d%io, j + d%jo), 0_int64)) then
               if (phys) then
                  nb_p = nb_p + 1
               else
                  nb_g = nb_g + 1
               end if
            end if
         end do
      end do
      if (nb_p + nb_g > 0) then
         write (*, '(a,i0,2a,a,i0,a,i0)') "  rank ", rank, " ", name, &
            ": physical mismatches=", nb_p, " ghost mismatches=", nb_g
      end if
      nbad_phys = nbad_phys + nb_p
      nbad_ghost = nbad_ghost + nb_g
   end subroutine compare_2d

   subroutine compare_3d(name, tile, whole, d, is_vface, nbad_phys, nbad_ghost, is_uface)
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: tile(:, :, :), whole(:, :, :)
      type(snap_t), intent(in) :: d
      logical, intent(in) :: is_vface
      integer, intent(inout) :: nbad_phys, nbad_ghost
      logical, intent(in), optional :: is_uface
      integer :: k
      character(len=64) :: lbl
      do k = 1, size(tile, 3)
         write (lbl, '(a,a,i0)') name, " k=", k
         call compare_2d(trim(lbl), tile(:, :, k), whole(:, :, k), d, is_vface, &
                         nbad_phys, nbad_ghost, is_uface)
      end do
   end subroutine compare_3d

   subroutine run_case(label, ny, lat_south, dlat, phi_join)
      !! The serial reference once, then every decomposition of the
      !! launched rank count: the north-south split (px = 1, py = nprocs)
      !! and the east-west splits through the distributed fold (engine_setup's
      !! test-only `allow_distributed_fold`): np 2 -> 2x1; np 3 -> 3x1;
      !! np 4 -> 4x1, 2x2.
      character(len=*), intent(in) :: label
      integer, intent(in) :: ny
      real(wp), intent(in) :: lat_south, dlat, phi_join
      type(snap_t) :: ref_cfg, ref_end

      call run_one(case_nml(ny, lat_south, dlat, phi_join, 1, 1), 1, 0, ref_cfg, ref_end, ok_ref_g)
      if (.not. ok_ref_g) then
         write (*, '(3a,i0)') "FAIL ", label, ": the serial reference failed on rank ", rank
         n_fail = n_fail + 1
         return
      end if
      ! Teeth: a non-finite reference would compare equal to a non-finite
      ! decomposed run.  Require the serial trajectory to stay finite.
      if (.not. (all(ieee_is_finite(ref_end%h)) .and. all(ieee_is_finite(ref_end%v)) .and. &
                 all(ieee_is_finite(ref_end%tr)))) then
         write (*, '(3a,i0)') "FAIL ", label, ": the serial reference went non-finite on rank ", rank
         n_fail = n_fail + 1
      end if
      if (rank == 0) then
         write (*, '(a,a,a,es23.15,a,es23.15)') "case ", label, ": ref sum(h)=", &
            sum(ref_end%h(4:3 + NX_G, 4:3 + ny, :)), "  sum(v^2)=", &
            sum(ref_end%v(4:3 + NX_G, 4:3 + ny + 1, :)**2)
      end if

      if (nprocs == 1) then
         call run_decomposed(label, ny, lat_south, dlat, phi_join, 1, 1, ref_cfg, ref_end)
         return
      end if
      call run_decomposed(label, ny, lat_south, dlat, phi_join, 1, nprocs, ref_cfg, ref_end)
      call run_decomposed(label, ny, lat_south, dlat, phi_join, nprocs, 1, ref_cfg, ref_end)
      if (nprocs == 4) call run_decomposed(label, ny, lat_south, dlat, phi_join, 2, 2, &
                                           ref_cfg, ref_end)
   end subroutine run_case

   subroutine run_decomposed(label, ny, lat_south, dlat, phi_join, px, py, ref_cfg, ref_end)
      character(len=*), intent(in) :: label
      integer, intent(in) :: ny, px, py
      real(wp), intent(in) :: lat_south, dlat, phi_join
      type(snap_t), intent(in) :: ref_cfg, ref_end
      type(snap_t) :: dec_cfg, dec_end
      logical :: ok_dec
      integer :: it, bad_cfg_p, bad_cfg_g, bad_p, bad_g, glob(4), nok

      call run_one(case_nml(ny, lat_south, dlat, phi_join, px, py), nprocs, rank, dec_cfg, &
                   dec_end, ok_dec, allow_dfold=(px > 1))
      nok = merge(0, 1, ok_dec)
      call allreduce(comm, nok, op=MPI_SUM)
      if (nok > 0) then
         if (.not. ok_dec) write (*, '(3a,i0,a,i0,a,i0)') "FAIL ", label, ": rank ", rank, &
            " decomposed setup/step failed at px=", px, " py=", py
         n_fail = n_fail + 1
         return
      end if

      ! ---- configure-time generator gate (metrics + seeded state) ----
      bad_cfg_p = 0
      bad_cfg_g = 0
      call compare_2d("dxT", dec_cfg%dxt, ref_cfg%dxt, dec_cfg, .false., bad_cfg_p, bad_cfg_g)
      call compare_2d("areaT", dec_cfg%areat, ref_cfg%areat, dec_cfg, .false., bad_cfg_p, bad_cfg_g)
      call compare_2d("dyCu", dec_cfg%dycu, ref_cfg%dycu, dec_cfg, .false., bad_cfg_p, bad_cfg_g)
      call compare_2d("dxCv", dec_cfg%dxcv, ref_cfg%dxcv, dec_cfg, .true., bad_cfg_p, bad_cfg_g)
      call compare_2d("dxBu", dec_cfg%dxbu, ref_cfg%dxbu, dec_cfg, .true., bad_cfg_p, bad_cfg_g)
      call compare_2d("areaBu", dec_cfg%areabu, ref_cfg%areabu, dec_cfg, .true., bad_cfg_p, bad_cfg_g)
      call compare_2d("geolatBu", dec_cfg%geolatbu, ref_cfg%geolatbu, dec_cfg, .true., &
                      bad_cfg_p, bad_cfg_g)
      call compare_2d("angle_dx", dec_cfg%angle, ref_cfg%angle, dec_cfg, .false., bad_cfg_p, bad_cfg_g)
      call compare_2d("f_corner", dec_cfg%fcorner, ref_cfg%fcorner, dec_cfg, .true., &
                      bad_cfg_p, bad_cfg_g)
      call compare_3d("seed h", dec_cfg%h, ref_cfg%h, dec_cfg, .false., bad_cfg_p, bad_cfg_g)
      call compare_3d("seed u", dec_cfg%u, ref_cfg%u, dec_cfg, .false., bad_cfg_p, bad_cfg_g, &
                      is_uface=.true.)
      call compare_3d("seed v", dec_cfg%v, ref_cfg%v, dec_cfg, .true., bad_cfg_p, bad_cfg_g)

      ! ---- after N_STEPS: the prognostic fields ----
      bad_p = 0
      bad_g = 0
      call compare_3d("h", dec_end%h, ref_end%h, dec_end, .false., bad_p, bad_g)
      call compare_3d("u", dec_end%u, ref_end%u, dec_end, .false., bad_p, bad_g, is_uface=.true.)
      call compare_3d("v", dec_end%v, ref_end%v, dec_end, .true., bad_p, bad_g)
      do it = 1, size(ref_end%tr, 4)
         call compare_3d("hTr", dec_end%tr(:, :, :, it), ref_end%tr(:, :, :, it), dec_end, &
                         .false., bad_p, bad_g)
      end do
      call compare_2d("eta", dec_end%eta, ref_end%eta, dec_end, .false., bad_p, bad_g)

      glob = [bad_cfg_p, bad_cfg_g, bad_p, bad_g]
      call allreduce(comm, glob, op=MPI_SUM)
      if (rank == 0) then
         write (*, '(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "case ", label, ": ny=", ny, &
            " px=", px, " py=", py, "  configure mismatches phys/ghost=", glob(1), "/", glob(2), &
            "  after ", N_STEPS, " steps phys/ghost=", glob(3), "/", glob(4)
      end if
      if (sum(glob) > 0) then
         if (rank == 0) write (*, '(3a,i0,a,i0)') "FAIL ", label, &
            ": not bit-identical to the serial run at px=", px, " py=", py
         n_fail = n_fail + 1
      end if
   end subroutine run_decomposed

   subroutine check_px_refused()
      !! An east-west split of the fold row must fail loud at configure.
      type(ocean_engine_t) :: engine
      integer :: ierr
      if (mod(nprocs, 2) /= 0) return
      call setup_engine(engine, case_nml(24, 59.0_wp, 1.0_wp, 74.0_wp, 2, nprocs/2), nprocs, rank, ierr)
      if (ierr == OCEAN_STATUS_OK) then
         write (*, '(a,i0,a)') "FAIL px_refused: rank ", rank, &
            " px = 2 with the tripolar fold was ACCEPTED (must be refused)"
         n_fail = n_fail + 1
      else if (rank == 0) then
         write (*, '(a)') "case px_refused: px = 2 tripolar fold refused at configure (ok)"
      end if
      call engine_teardown(engine)
   end subroutine check_px_refused

end program test_ocean_tripolar_fold_mpi
#else
program test_ocean_tripolar_fold_mpi
   implicit none
   write (*, '(a)') "test_ocean_tripolar_fold_mpi: skipped (RDB_ENABLE_MPI=OFF)"
end program test_ocean_tripolar_fold_mpi
#endif
