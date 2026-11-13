!! Decomposition bit-identity of the bulk-forcing wind-stress sequence under
!! the tripolar north fold (`bulk_forcing_apply`: scalar halo of the
!! cell-centred pair, vector fold of the pair into the north ghost rows, face
!! average, face-pair seam refresh incl. the distributed fold at px > 1,
!! `stress_mag`).
!!
!! The atmosphere is analytic in the GLOBAL index and obeys the fold symmetry
!! of a true vector (`uas`, `vas` odd about the fold's mirror column), so the
!! stress is a continuous vector field across the fold and the fold-line
!! `tau_y` face must be the two-sided average.  The 1-rank run (every rank
!! builds the whole sector on one logical tile) is the reference; the 2x2,
!! 4x1 (px = 4: distributed fold) and 1x4 splits must reproduce its owned
!! `tau_x`, `tau_y` and `stress_mag` bitwise.  Teeth: the reference fold-line
!! `tau_y` is checked against the two-sided average and found nonzero.
#ifdef RDB_ENABLE_MPI
program test_ocean_bulk_stress_fold_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_ocean_fold_exchange, only: ocean_fold_exchange_init, ocean_fold_exchange_destroy
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_atm_state, only: atm_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_bulk_forcing, only: bulk_forcing_t, bulk_forcing_apply
   implicit none

   integer, parameter :: NI = 16, NJ = 16, NG = 3
   real(wp), parameter :: SST_DEGC = 15.0_wp, H0 = 10.0_wp

   type :: snap_t
      real(wp), allocatable :: tau_x(:, :), tau_y(:, :), smag(:, :), tc_last(:)
      integer :: io = 0, jo = 0, nxl = 0, nyl = 0
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail, nbad
   type(comm_t) :: comm
   type(snap_t) :: ref, dec

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0
   if (nprocs /= 1 .and. nprocs /= 4) error stop "test_ocean_bulk_stress_fold_mpi: 1 or 4 ranks"

   call run_one(1, 1, 0, ref)
   call check_reference(ref)

   if (nprocs == 4) then
      call run_one(2, 2, rank, dec)
      call compare_all(dec, ref, "2x2")
      call run_one(4, 1, rank, dec)
      call compare_all(dec, ref, "4x1")
      call run_one(1, 4, rank, dec)
      call compare_all(dec, ref, "1x4")
   end if

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_bulk_stress_fold_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_bulk_stress_fold_mpi: all checks PASSED (nprocs=", nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   subroutine check_reference(r)
      !! The reference fold-line tau_y (row NG+NJ+1) is the two-sided
      !! average of tauy_c, hence odd in i and nonzero; and the whole
      !! stress is live and finite.
      type(snap_t), intent(in) :: r
      integer :: i, jf
      real(wp) :: odd_err, two_sided_err

      jf = NG + NJ + 1
      odd_err = 0.0_wp
      two_sided_err = 0.0_wp
      do i = 1, NI
         two_sided_err = max(two_sided_err, abs(r%tau_y(NG + i, jf) - 0.5_wp*(r%tc_last(i) - r%tc_last(NI + 1 - i))))
         odd_err = max(odd_err, abs(r%tau_y(NG + i, jf) + r%tau_y(NG + NI + 1 - i, jf)))
      end do
      if (.not. all(ieee_is_finite(r%tau_y)) .or. maxval(abs(r%tau_y(NG + 1:NG + NI, jf))) <= 1.0e-3_wp &
          .or. odd_err /= 0.0_wp .or. two_sided_err > 1.0e-14_wp .or. maxval(abs(r%smag)) <= 0.0_wp) then
         if (rank == 0) write (*, '(a)') "reference fold-line tau_y vacuous / not antisymmetric"
         n_fail = n_fail + 1
      end if
   end subroutine check_reference

   subroutine compare_all(d, r, tag)
      type(snap_t), intent(in) :: d, r
      character(len=*), intent(in) :: tag
      call compare2d(d%tau_x, d%io, d%jo, d%nxl, d%nyl, r%tau_x, tag//" tau_x", nbad)
      n_fail = n_fail + nbad
      call compare2d(d%tau_y, d%io, d%jo, d%nxl, d%nyl, r%tau_y, tag//" tau_y", nbad)
      n_fail = n_fail + nbad
      call compare2d(d%smag, d%io, d%jo, d%nxl, d%nyl, r%smag, tag//" stress_mag", nbad)
      n_fail = n_fail + nbad
   end subroutine compare_all

   subroutine run_one(px, py, crank, snap)
      integer, intent(in) :: px, py, crank
      type(snap_t), intent(out) :: snap

      type(decomp_t) :: decomp
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer :: i, j, nx, ny
      real(wp) :: d

      call decomp_init(decomp, NI, NJ, px, py, crank)
      call grid%init(decomp%nx_local, decomp%ny_local, NG, 1.0_wp, 1.0_wp)
      grid%i_offset_global = decomp%i_start - 1
      grid%j_offset_global = decomp%j_start - 1
      grid%nx_global = decomp%nx_global
      grid%ny_global = decomp%ny_global
      nx = grid%nx_total
      ny = grid%ny_total

      call ocean_halo_destroy()
      call ocean_halo_init(decomp, NG, .true., .false.)
      bc%periodic_x = .true.
      bc%north_fold = decomp%has_north
      call ocean_fold_exchange_init(decomp, NG, bc%north_fold)

      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = 1
      call ms%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%c = SST_DEGC
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      allocate (atm%tas(nx, ny), source=290.0_wp)
      allocate (atm%huss(nx, ny), source=0.008_wp)
      allocate (atm%uas(nx, ny), atm%vas(nx, ny))
      allocate (atm%psl(nx, ny), source=101000.0_wp)
      allocate (atm%rlds(nx, ny), source=300.0_wp)
      allocate (atm%rsds(nx, ny), source=200.0_wp)
      allocate (atm%prra(nx, ny), source=2.0e-5_wp)
      allocate (atm%prsn(nx, ny), source=1.0e-6_wp)
      do j = 1, ny
         do i = 1, nx
            ! Global column about the fold's mirror point: odd in it.
            d = real(i - NG + grid%i_offset_global, wp) - 0.5_wp*real(NI + 1, wp)
            atm%uas(i, j) = 0.4_wp*d + 0.01_wp*d**3
            atm%vas(i, j) = 0.6_wp*d - 0.02_wp*d**3 + 0.05_wp*d*real(j + grid%j_offset_global, wp)
         end do
      end do
      atm%enable = .true.
      atm%nx_total = nx
      atm%ny_total = ny
      atm%nx_phys = grid%nx_phys
      atm%ny_phys = grid%ny_phys
      atm%nghost = NG

      bf%enable = .true.
      bf%nx_total = nx
      bf%ny_total = ny
      bf%cfg_atm%wind_relative = .true.
      bf%cfg_atm%ocean_albedo = 0.066_wp
      bf%cfg_bulk%n_iter = 5
      bf%cfg_bulk%z_wind = 10.0_wp
      bf%cfg_bulk%z_ta = 2.0_wp
      bf%cfg_bulk%emissivity = 0.97_wp
      allocate (bf%taux_c(nx, ny), source=0.0_wp)
      allocate (bf%tauy_c(nx, ny), source=0.0_wp)
      allocate (bf%open_water(nx, ny), source=1.0_wp)

      call ss%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)

      ! mem:separate: map everything the apply touches (component arrays
      ! only are read back; directives are inert on CPU builds).
      !$acc enter data copyin(atm, ms, bf, ss, sf)
      call atm%enter_data()
      call ms%enter_data()
      call bf%enter_data()
      call ss%enter_data()
      call sf%enter_data()

      call bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)
      !$acc update self(ss%tau_x, ss%tau_y, ss%stress_mag, bf%tauy_c)

      call sf%exit_data()
      call ss%exit_data()
      call bf%exit_data()
      call ms%exit_data()
      call atm%exit_data()
      !$acc exit data delete(atm, ms, bf, ss, sf)

      allocate (snap%tau_x(size(ss%tau_x, 1), size(ss%tau_x, 2)), source=ss%tau_x)
      allocate (snap%tau_y(size(ss%tau_y, 1), size(ss%tau_y, 2)), source=ss%tau_y)
      allocate (snap%smag(size(ss%stress_mag, 1), size(ss%stress_mag, 2)), source=ss%stress_mag)
      allocate (snap%tc_last(NI), source=0.0_wp)
      if (px == 1 .and. py == 1) snap%tc_last = bf%tauy_c(NG + 1:NG + NI, NG + NJ)
      snap%io = grid%i_offset_global
      snap%jo = grid%j_offset_global
      snap%nxl = grid%nx_phys
      snap%nyl = grid%ny_phys

      call ocean_fold_exchange_destroy()
      call ocean_halo_destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_one

   subroutine compare2d(dec_a, io, jo, nxl, nyl, ref_a, tag, nbad_out)
      !! Bitwise, owned window (a face array's extra high-side row/column
      !! included); `ref_a` is the global-sized 1-rank array.
      real(wp), intent(in) :: dec_a(:, :), ref_a(:, :)
      integer, intent(in) :: io, jo, nxl, nyl
      character(len=*), intent(in) :: tag
      integer, intent(out) :: nbad_out
      integer :: i, j, ex, ey

      ex = size(dec_a, 1) - (nxl + 2*NG)
      ey = size(dec_a, 2) - (nyl + 2*NG)
      nbad_out = 0
      do j = NG + 1, NG + nyl + ey
         do i = NG + 1, NG + nxl + ex
            if (transfer(dec_a(i, j), 0_int64) /= transfer(ref_a(i + io, j + jo), 0_int64)) then
               if (nbad_out == 0) write (*, '(a,i0,2a,2i4,2es24.16)') "  rank ", rank, " ", tag, i, j, &
                  dec_a(i, j), ref_a(i + io, j + jo)
               nbad_out = nbad_out + 1
            end if
         end do
      end do
   end subroutine compare2d

end program test_ocean_bulk_stress_fold_mpi
#endif
