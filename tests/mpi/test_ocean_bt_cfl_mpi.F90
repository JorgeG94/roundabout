!! Regression gate: barotropic `n_inner` must be GLOBALLY consistent across
!! ranks on a spherical grid under a meridional (py>1) decomposition.
!!
!! THE BUG (feat/ocean-mpi-marchin): `configure_ocean_bt_split` derived the
!! external-gravity-wave CFL substep count `n_inner` from the LOCAL rank's
!! `metrics_bt_cfl_length` (min over the rank's own physical cells) while only
!! the wave speed `c_ext = sqrt(g*H_max)` was globally reduced.  On a spherical
!! lon-lat grid `dxT = R*cos(lat)*dlon` shrinks poleward, so a N/S (py>1) split
!! gives each rank a DIFFERENT local min length -> a DIFFERENT `n_inner` -> the
!! ranks run a different number of barotropic substeps -> their per-substep
!! grouped halo exchanges desync and MPI aborts with MPI_ERR_TRUNCATE at the
!! first N/S exchange.  Cartesian (uniform dxT) and a zonal (px>1) spherical
!! split give every rank the same local length, which is why only spherical
!! py>1 exposed it.  Fix: reduce the CFL length with `halo_allreduce_min`.
!!
!! This test rebuilds ONLY the CFL-length -> n_inner derivation (not a full
!! solver step) on each rank's spherical subdomain and asserts:
!!   TEETH  (nprocs>1): the LOCAL-length n_inner DIFFERS across ranks (proves
!!                      the latitude band genuinely diverges — without this the
!!                      test could pass on a degenerate grid).
!!   FIX    (all nprocs): the GLOBAL-min-length n_inner is IDENTICAL across
!!                        ranks (the desync can no longer happen).
!! On nprocs=1 the reduction is an identity stub, TEETH is vacuous, FIX trivial.
#ifdef RDB_ENABLE_MPI
program test_ocean_bt_cfl_mpi
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
   use rdb_ocean_setup, only: metrics_bt_cfl_length, bt_auto_n_inner
   use rdb_halo, only: halo_allreduce_min
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM, MPI_MIN, MPI_MAX
   implicit none

   ! --- Global spherical sector (divisible by 1, 2 and 4 for the py legs) ---
   integer, parameter :: NX_G = 16
   integer, parameter :: NY_G = 32
   integer, parameter :: NG = 3
   real(wp), parameter :: DLON = 1.0_wp     !! deg
   real(wp), parameter :: DLAT = 1.0_wp     !! deg
   real(wp), parameter :: LON_W = 0.0_wp
   real(wp), parameter :: LAT_S = -45.0_wp  !! high-|lat| south edge => strong cos(lat) spread
   real(wp), parameter :: RAD_EARTH = 6.371e6_wp
   real(wp), parameter :: H_MAX = 5500.0_wp
   real(wp), parameter :: DT_OUTER = 1800.0_wp
   real(wp), parameter :: CFL_SAFETY = 0.65_wp

   integer :: rank, nprocs, n_fail, total_fail
   integer :: j_off, nyl
   integer :: n_local, n_global, n_loc_min, n_loc_max, n_glob_min, n_glob_max
   real(wp) :: l_cfl_local, l_cfl_global, c_ext
   type(comm_t) :: comm
   type(decomp_t) :: decomp
   type(hgrid_t) :: grid
   type(ocean_metrics_t) :: metrics

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0

   ! --- Build this rank's spherical subdomain (py = nprocs meridional split) ---
   call decomp_init(decomp, NX_G, NY_G, 1, nprocs, rank)
   nyl = decomp%ny_local
   j_off = decomp%j_start - 1
   call grid%init(NX_G/1, nyl, NG, DLON, DLAT)
   grid%i_offset_global = 0
   grid%j_offset_global = j_off
   call make_spherical_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

   ! --- Derive n_inner two ways: from the LOCAL length (the bug) and from the
   !     GLOBAL-min length (the fix). ---
   c_ext = sqrt(GRAVITY*H_MAX)
   l_cfl_local = metrics_bt_cfl_length(metrics, grid)
   call halo_allreduce_min(l_cfl_local, l_cfl_global)
   n_local = bt_auto_n_inner(DT_OUTER, CFL_SAFETY, c_ext, l_cfl_local)
   n_global = bt_auto_n_inner(DT_OUTER, CFL_SAFETY, c_ext, l_cfl_global)

   ! Spread of each across ranks (min/max reductions; equal <=> identical).
   n_loc_min = n_local; call allreduce(comm, n_loc_min, MPI_MIN)
   n_loc_max = n_local; call allreduce(comm, n_loc_max, MPI_MAX)
   n_glob_min = n_global; call allreduce(comm, n_glob_min, MPI_MIN)
   n_glob_max = n_global; call allreduce(comm, n_glob_max, MPI_MAX)

   if (rank == 0) then
      write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') &
         "nprocs=", nprocs, "  local n_inner spread=[", n_loc_min, ",", n_loc_max, &
         "]  global n_inner spread=[", n_glob_min, ",", n_glob_max
      write (*, '(a,es14.6,a,es14.6)') &
         "  l_cfl_local(rank0)=", l_cfl_local, "  l_cfl_global=", l_cfl_global
   end if

   ! --- FIX: the global-min-length n_inner is IDENTICAL on every rank. ---
   if (n_glob_min /= n_glob_max) then
      n_fail = n_fail + 1
      write (*, '(a,i0,a,i0,a,i0)') &
         "FAIL fix: rank ", rank, " global n_inner not consistent across ranks: min=", &
         n_glob_min, " max=", n_glob_max
   end if

   ! --- TEETH: on a real meridional split the LOCAL-length n_inner DIFFERS
   !     across ranks (this is exactly the desync that truncated the halo). ---
   if (nprocs > 1) then
      if (n_loc_min == n_loc_max) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,i0,a)') &
            "FAIL teeth: rank ", rank, " LOCAL n_inner did NOT diverge across the ", &
            nprocs, "-way meridional spherical split — test has no teeth (check grid)."
      end if
      ! And the fix must pick the GLOBAL MIN length => the LARGER n_inner
      ! (finest-cell / most-substeps rank), so every rank is at least as
      ! finely subcycled as the tightest cell anywhere in the domain.
      if (n_glob_min < n_loc_max) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,i0,a,i0)') &
            "FAIL fix-bound: rank ", rank, " global n_inner=", n_glob_min, &
            " is below the domain-max local requirement=", n_loc_max
      end if
   end if

   call destroy_cartesian_metrics(metrics)

   ! --- Report ---
   call comm%barrier()
   if (n_fail > 0) then
      write (*, '(a,i0,a,i0,a,i0,a)') &
         "Rank ", rank, ": ", n_fail, " check(s) FAILED (nprocs=", nprocs, ")"
   else
      write (*, '(a,i0,a,i0,a)') &
         "Rank ", rank, ": all checks PASSED (nprocs=", nprocs, ")"
   end if
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

end program test_ocean_bt_cfl_mpi
#else
program test_ocean_bt_cfl_mpi
   implicit none
   ! Non-MPI build: nothing to exercise (the reduction is a single-rank identity).
   write (*, '(a)') "test_ocean_bt_cfl_mpi: skipped (RDB_ENABLE_MPI=OFF)"
end program test_ocean_bt_cfl_mpi
#endif
