!! Distributed tripolar fold (`px > 1`) at the baroclinic / continuity /
!! stress seam sites, on a real decomposed engine state, vs the serial run.
!!
!! A COMPONENT test, through the production entry points, of every
!! baroclinic / continuity / stress fold site in isolation, each on a
!! field whose fold-line row is deliberately NOT projected, so a site that
!! skipped its fold or its projection fails here even if a stepped run
!! would mask it (stepped bit-identity: test_ocean_tripolar_fold_mpi and
!! test_ocean_decomp_bitid_mpi).  On each launch (np 2, 3, 4) and each
!! factorisation with px > 1 (2x1; 3x1; 4x1 = 8/8/7/7 and 2x2), the binary
!! runs the same sequence twice — SERIAL reference (px = py = 1, every
!! rank) and DECOMPOSED —
!! snapshots every compared field after every operation, and compares each
!! rank's WHOLE storage window (ghost rows and columns) of the decomposed
!! run bitwise against the matching window of the reference:
!!
!!   E   `engine_setup` (cold start): b, bt_H_ref, h, u, v, every hTr —
!!       the init-time folds, deferred after the host halo on px > 1.
!!   S   a structured seed + host halo + `ocean_fold_wrap_state`
!!       (`device_resident=.false.`): the host path of D1/D2/D4.
!!   then, on the DEVICE-MAPPED state (after `engine_enter_data`), every
!!   field filled from an index-encoded global function (north ghosts NaN,
!!   fold-line row unprojected) and pushed to the device before each site:
!!   D   `ocean_fold_wrap_state`            h, u, v, hTr  (D1, D2, D4)
!!   C   `ocean_fold_wrap_centre_3d_state`  h, hTr        (C1, D5)
!!   M   `ocean_fold_wrap_time_means`       u_av, v_av, h_av (D3)
!!   F   `ocean_fold_north_v_face`          mass_flux_y_layer (C2)
!!   H   `ocean_seam_refresh_surface_stress` tau_x, tau_y, stress_mag (H1)
!! On the GPU build the D..H legs run the device pack/unpack and hand MPI
!! device buffers, exactly as at step time.
#ifdef RDB_ENABLE_MPI
program test_ocean_fold_sites_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state, &
                                   ocean_seam_refresh_surface_stress
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_state, ocean_fold_wrap_centre_3d_state, &
                                   ocean_fold_wrap_time_means
   use rdb_ocean_fold_exchange, only: ocean_fold_north_v_face, ocean_fold_is_distributed
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   integer, parameter :: NX_G = 30
      !! Global columns (4x1 gives the uneven 8/8/7/7 split).
   integer, parameter :: NY_G = 24
      !! Global rows (the "spanning" cap of test_ocean_tripolar_fold_mpi).
   integer, parameter :: NZ = 3
   integer, parameter :: NG = 3

   type :: fld_t
      !! One snapshotted field (whole local storage, as 3D).
      character(len=32) :: name = ""
      real(wp), allocatable :: a(:, :, :)
      logical :: yface = .false.
   end type fld_t

   type :: run_t
      !! Every snapshot of one run, in operation order.
      type(fld_t), allocatable :: f(:)
      integer :: n = 0
      integer :: io = 0
      integer :: jo = 0
      logical :: has_south = .true.
      logical :: has_north = .true.
      logical :: distributed = .false.
   end type run_t

   integer :: rank, nprocs, n_fail, total_fail
   type(comm_t) :: comm
   real(wp) :: nan

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0
   nan = ieee_value(1.0_wp, ieee_quiet_nan)

   select case (nprocs)
   case (2)
      call run_case(2, 1)
   case (3)
      call run_case(3, 1)
   case (4)
      call run_case(4, 1)
      call run_case(2, 2)
   case default
      if (rank == 0) write (*, '(a,i0)') "SKIP: test_ocean_fold_sites_mpi needs 2-4 ranks, got ", &
         nprocs
   end select

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_fold_sites_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_fold_sites_mpi: all checks PASSED (nprocs=", &
            nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   function case_nml(px, py) result(nml)
      !! Small analytic tripolar ocean (test_ocean_tripolar_fold_mpi's
      !! "spanning" grid at nx = 30, dx = 12 deg — 360/30 is exact).
      integer, intent(in) :: px, py
      character(len=:), allocatable :: nml
      character(len=16) :: spx, spy
      character(len=*), parameter :: NL = new_line("a")
      write (spx, '(i0)') px
      write (spy, '(i0)') py
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 30, ny = 24, nghost = 3, dx = 12.0, dy = 1.0 /"//NL// &
            "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
            "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, "// &
            "lat_south = 59.0, phi_join = 74.0, lon_pole = 0.0, "// &
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

   ! ----------------------------------------------------------------
   ! Snapshots
   ! ----------------------------------------------------------------

   subroutine snap3(r, name, a, yface)
      type(run_t), intent(inout) :: r
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: a(:, :, :)
      logical, intent(in) :: yface
      type(fld_t), allocatable :: tmp(:)
      if (.not. allocated(r%f)) allocate (r%f(64))
      if (r%n == size(r%f)) then
         allocate (tmp(2*size(r%f)))
         tmp(1:r%n) = r%f(1:r%n)
         call move_alloc(tmp, r%f)
      end if
      r%n = r%n + 1
      r%f(r%n)%name = name
      r%f(r%n)%a = a
      r%f(r%n)%yface = yface
   end subroutine snap3

   subroutine snap2(r, name, a, yface)
      type(run_t), intent(inout) :: r
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: a(:, :)
      logical, intent(in) :: yface
      call snap3(r, name, reshape(a, [size(a, 1), size(a, 2), 1]), yface)
   end subroutine snap2

   subroutine snap_ml(r, tag, e)
      type(run_t), intent(inout) :: r
      character(len=*), intent(in) :: tag
      type(ocean_engine_t), intent(in) :: e
      integer :: it
      character(len=32) :: nm
      call snap3(r, tag//" h", e%state%multilayer%h_layer, .false.)
      call snap3(r, tag//" u", e%state%multilayer%u_face_x_layer, .false.)
      call snap3(r, tag//" v", e%state%multilayer%v_face_y_layer, .true.)
      do it = 1, size(e%state%multilayer%tracers)
         write (nm, '(a,a,i0)') tag, " hTr", it
         call snap3(r, trim(nm), e%state%multilayer%tracers(it)%hTr, .false.)
      end do
   end subroutine snap_ml

   ! ----------------------------------------------------------------
   ! Index-encoded fills (global storage coordinates; x-periodic)
   ! ----------------------------------------------------------------

   subroutine fill3(a, e, fid, yext)
      !! Fill a whole local storage array from the global function: every
      !! row up to the last pre-fold row (north ghosts NaN), the fold-line
      !! row deliberately NOT antisymmetric (the site must project it).
      real(wp), intent(inout) :: a(:, :, :)
      type(ocean_engine_t), intent(in) :: e
      integer, intent(in) :: fid, yext
      integer :: i, j, k, c, jg
      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            jg = j + e%grid%j_offset_global
            do i = 1, size(a, 1)
               c = modulo(i + e%grid%i_offset_global - NG - 1, NX_G) + 1
               if (jg > NG + NY_G + yext) then
                  a(i, j, k) = nan
               else
                  a(i, j, k) = real(fid, wp)*1.0e6_wp + real(c, wp)*1.0e3_wp + &
                               real(jg, wp)*10.0_wp + real(k, wp) + 0.25_wp
               end if
            end do
         end do
      end do
   end subroutine fill3

   subroutine fill2(a, e, fid, yext)
      real(wp), intent(inout) :: a(:, :)
      type(ocean_engine_t), intent(in) :: e
      integer, intent(in) :: fid, yext
      real(wp), allocatable :: t(:, :, :)
      allocate (t(size(a, 1), size(a, 2), 1))
      call fill3(t, e, fid, yext)
      a = t(:, :, 1)
   end subroutine fill2

   subroutine seed_structure(e)
      !! Decomposition-invariant structure (functions of GLOBAL indices) in
      !! u, v (fold row included) and every tracer, then the ghosts refilled
      !! as `engine_setup` does: host halo, then the fold.  Exact rational
      !! values, NOT transcendentals: a vectorised libm `sin`/`cos` (gfortran
      !! -O3 -march=native → libmvec) differs from the scalar one in the last
      !! bit, and which cells take which path depends on the tile's row
      !! length — a seed artifact that would differ by decomposition.
      type(ocean_engine_t), intent(inout) :: e
      integer :: i, j, k, it, nxl, nyl, io, jo, ig, jg

      nxl = e%grid%nx_phys
      nyl = e%grid%ny_phys
      io = e%grid%i_offset_global
      jo = e%grid%j_offset_global
      do k = 1, NZ
         do j = NG + 1, NG + nyl
            do i = NG + 1, NG + nxl + 1
               ! Reduced face index: face NX_G+1 IS face 1 (bit-equal
               ! copies, the periodic-seam invariant).
               ig = modulo(i - NG + io - 1, NX_G) + 1
               jg = j - NG + jo
               e%state%multilayer%u_face_x_layer(i, j, k) = &
                  real(modulo(5*ig + 3*jg + k, 13) - 6, wp)/128.0_wp
            end do
         end do
         do j = NG + 1, NG + nyl + 1
            do i = NG + 1, NG + nxl
               ig = i - NG + io
               jg = j - NG + jo
               e%state%multilayer%v_face_y_layer(i, j, k) = &
                  real(modulo(7*ig + 2*jg + 3*k, 11) - 5, wp)/256.0_wp
            end do
         end do
         do it = 1, size(e%state%multilayer%tracers)
            do j = NG + 1, NG + nyl
               do i = NG + 1, NG + nxl
                  ig = i - NG + io
                  jg = j - NG + jo
                  e%state%multilayer%tracers(it)%hTr(i, j, k) = &
                     e%state%multilayer%tracers(it)%hTr(i, j, k)* &
                     (1.0_wp + real(modulo(3*ig + jg + it, 7) - 3, wp)/1024.0_wp)
               end do
            end do
         end do
      end do
      call ocean_halo_exchange_ml_state(e%state%multilayer, device_resident=.false.)
      call ocean_fold_wrap_state(e%grid, e%state%bc, e%state%multilayer, device_resident=.false.)
   end subroutine seed_structure

   ! ----------------------------------------------------------------
   ! One run: the whole operation sequence, snapshotted
   ! ----------------------------------------------------------------

   subroutine run_one(px, py, csize, crank, r, ok)
      integer, intent(in) :: px, py, csize, crank
      type(run_t), intent(out) :: r
      logical, intent(out) :: ok
      type(ocean_engine_t) :: e
      type(config_t) :: cfg
      integer :: ierr, it, nt

      ok = .false.
      call read_config_from_string(case_nml(px, py), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(e, cfg, ierr, compute_rank=crank, compute_size=csize)
      if (ierr /= OCEAN_STATUS_OK) return
      r%io = e%grid%i_offset_global
      r%jo = e%grid%j_offset_global
      r%has_south = e%decomp%has_south
      r%has_north = e%decomp%has_north
      r%distributed = ocean_fold_is_distributed()

      associate (ms => e%state%multilayer, ss => e%state%surface_stress)
         nt = size(ms%tracers)
         ! E: engine_setup's init-time folds (cold start).
         call snap2(r, "E b", e%state%barotropic%b, .false.)
         call snap2(r, "E bt_H_ref", e%state%dyn%bt_work%bt_H_ref, .false.)
         call snap_ml(r, "E", e)

         ! S: host-side dispatcher on a structured state.
         call seed_structure(e)
         call snap_ml(r, "S", e)

         call engine_enter_data(e, cfg)

         ! D: the stage-entry / post-continuity / stage-end group.
         call fill3(ms%h_layer, e, 1, 0)
         call fill3(ms%u_face_x_layer, e, 2, 0)
         call fill3(ms%v_face_y_layer, e, 3, 1)
         do it = 1, nt
            call fill3(ms%tracers(it)%hTr, e, 10 + it, 0)
         end do
         !$acc update device(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         do it = 1, nt
            !$acc update device(ms%tracers(it)%hTr)
         end do
         call ocean_fold_wrap_state(e%grid, e%state%bc, ms)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         do it = 1, nt
            !$acc update self(ms%tracers(it)%hTr)
         end do
         call snap_ml(r, "D", e)

         ! C: the centre-field group (mid Lie split, tracer refresh).
         call fill3(ms%h_layer, e, 4, 0)
         do it = 1, nt
            call fill3(ms%tracers(it)%hTr, e, 20 + it, 0)
         end do
         !$acc update device(ms%h_layer)
         do it = 1, nt
            !$acc update device(ms%tracers(it)%hTr)
         end do
         call ocean_fold_wrap_centre_3d_state(e%grid, e%state%bc, ms)
         !$acc update self(ms%h_layer)
         do it = 1, nt
            !$acc update self(ms%tracers(it)%hTr)
         end do
         call snap3(r, "C h", ms%h_layer, .false.)
         do it = 1, nt
            call snap3(r, "C hTr", ms%tracers(it)%hTr, .false.)
         end do

         ! M: the pred_corr time-means (allocated under the default scheme).
         if (allocated(ms%u_av_layer)) then
            call fill3(ms%u_av_layer, e, 5, 0)
            call fill3(ms%v_av_layer, e, 6, 1)
            call fill3(ms%h_av_layer, e, 7, 0)
            !$acc update device(ms%u_av_layer, ms%v_av_layer, ms%h_av_layer)
            call ocean_fold_wrap_time_means(e%grid, e%state%bc, ms)
            !$acc update self(ms%u_av_layer, ms%v_av_layer, ms%h_av_layer)
            call snap3(r, "M u_av", ms%u_av_layer, .false.)
            call snap3(r, "M v_av", ms%v_av_layer, .true.)
            call snap3(r, "M h_av", ms%h_av_layer, .false.)
         end if

         ! F: the final meridional mass flux projection (no halo before it).
         call fill3(ms%mass_flux_y_layer, e, 8, 1)
         !$acc update device(ms%mass_flux_y_layer)
         if (e%state%bc%north_fold) then
            call ocean_fold_north_v_face(ms%mass_flux_y_layer, e%grid%nx_total, &
                                         e%grid%ny_total + 1, NZ, e%grid%nx_phys, &
                                         e%grid%ny_phys, NG)
         end if
         !$acc update self(ms%mass_flux_y_layer)
         call snap3(r, "F mass_flux_y", ms%mass_flux_y_layer, .true.)

         ! H: the surface-stress seam refresh (halo + wrap + fold + |tau|).
         call fill2(ss%tau_x, e, 30, 0)
         call fill2(ss%tau_y, e, 31, 1)
         !$acc update device(ss%tau_x, ss%tau_y)
         call ocean_seam_refresh_surface_stress(ss, e%grid, e%state%bc)
         !$acc update self(ss%tau_x, ss%tau_y, ss%stress_mag)
         call snap2(r, "H tau_x", ss%tau_x, .false.)
         call snap2(r, "H tau_y", ss%tau_y, .true.)
         call snap2(r, "H stress_mag", ss%stress_mag, .false.)
      end associate

      call engine_exit_data(e)
      call engine_teardown(e)
      ok = .true.
   end subroutine run_one

   subroutine run_case(px, py)
      integer, intent(in) :: px, py
      type(run_t) :: ref, dec
      logical :: ok_ref, ok_dec
      integer :: n, nbad, nbad_case, glob(2)
      character(len=24) :: lbl

      write (lbl, '(i0,a,i0)') px, "x", py
      call run_one(1, 1, 1, 0, ref, ok_ref)
      call run_one(px, py, nprocs, rank, dec, ok_dec)
      if (.not. (ok_ref .and. ok_dec)) then
         write (*, '(3a,i0,a,l1,a,l1)') "FAIL ", trim(lbl), ": rank ", rank, &
            " setup failed: ref ok=", ok_ref, " decomposed ok=", ok_dec
         n_fail = n_fail + 1
         return
      end if
      if (dec%has_north .and. .not. dec%distributed) then
         write (*, '(3a,i0)') "FAIL ", trim(lbl), ": the distributed fold is not active on rank ", rank
         n_fail = n_fail + 1
      end if
      if (ref%n /= dec%n) then
         write (*, '(3a,i0)') "FAIL ", trim(lbl), ": snapshot count differs on rank ", rank
         n_fail = n_fail + 1
         return
      end if

      nbad_case = 0
      do n = 1, ref%n
         call compare(ref%f(n), dec%f(n), dec, nbad)
         if (nbad > 0) write (*, '(a,i0,4a,i0)') "  rank ", rank, " ", trim(lbl), " ", &
            trim(dec%f(n)%name)//": mismatches=", nbad
         nbad_case = nbad_case + nbad
         if (.not. all(ieee_is_finite(ref%f(n)%a))) then
            write (*, '(5a)') "FAIL ", trim(lbl), ": reference ", trim(ref%f(n)%name), &
               " is not finite"
            n_fail = n_fail + 1
         end if
      end do
      glob = [nbad_case, ref%n]
      call allreduce(comm, glob, op=MPI_SUM)
      if (rank == 0) write (*, '(3a,i0,a,i0)') "case ", trim(lbl), ": fields compared (sum over ranks)=", &
         glob(2), "  mismatching slots=", glob(1)
      if (nbad_case > 0) n_fail = n_fail + 1
   end subroutine run_case

   subroutine compare(fr, fd, dec, nbad)
      !! Bitwise: the decomposed tile window vs the reference window.  One
      !! row is exempt, as in test_ocean_tripolar_fold_mpi: the outermost
      !! y-face row beyond an MPI y seam, which no halo pass refreshes.
      type(fld_t), intent(in) :: fr, fd
      type(run_t), intent(in) :: dec
      integer, intent(out) :: nbad
      integer :: i, j, k
      nbad = 0
      do k = 1, size(fd%a, 3)
         do j = 1, size(fd%a, 2)
            if (fd%yface .and. j == 1 .and. .not. dec%has_south) cycle
            if (fd%yface .and. j == size(fd%a, 2) .and. .not. dec%has_north) cycle
            do i = 1, size(fd%a, 1)
               if (transfer(fd%a(i, j, k), 0_int64) /= &
                   transfer(fr%a(i + dec%io, j + dec%jo, k), 0_int64)) nbad = nbad + 1
            end do
         end do
      end do
   end subroutine compare

end program test_ocean_fold_sites_mpi
#endif
