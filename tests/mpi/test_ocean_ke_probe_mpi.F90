!! MPI test of the KE-attribution probe's reductions (`debug_ke_attr`).
!! Run with 1, 2 or 4 ranks: mpirun -np N ./test_ocean_ke_probe_mpi
!!
!! A non-uniform, non-symmetric analytic C-grid state (a closed form of the
!! GLOBAL cell index, ghosts included) is cut into tiles on every px x py
!! factorisation of the rank count (uneven splits).  Each rank calls the
!! probe's rank-local EFP routines, the partials are combined with
!! `halo_allreduce_efp_list`, and the result is compared bitwise against an
!! independent serial EFP sum over the global field with plainly written
!! region masks.  Control: the pre-fix semantics (tile-local indices, tile
!! extents as the "global" extent) must DIFFER on more than one rank, so the
!! test has teeth.
#ifdef RDB_ENABLE_MPI
program test_ocean_ke_probe_mpi
   use, intrinsic :: iso_fortran_env, only: real64
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_from_real, efp_plus, efp_to_real
   use rdb_halo, only: halo_allreduce_efp_list
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size
   use rdb_ocean_ke_probe, only: ke_regions_efp_local, coradv_split_efp_local, N_KE, N_SPLIT
   implicit none

   integer, parameter :: NXG = 12, NYG = 10, NZ = 3, NG = 2
   integer, parameter :: RIM = 3, EAST = 4
   real(wp), parameter :: DT = 0.37_wp
   integer :: rank, nprocs, px, py, nfail, ntry

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   nfail = 0
   ntry = 0

   do px = 1, nprocs
      if (mod(nprocs, px) /= 0) cycle
      py = nprocs/px
      call run_case(px, py)
   end do

   if (rank == 0) write (*, "(a,i0,a,i0,a)") "ke_probe_mpi: ", ntry, " checks, ", nfail, " failures"
   call comm_env_finalize()
   if (nfail > 0) error stop 1

contains

   pure function hf(ig, jg, k) result(r)
      integer, intent(in) :: ig, jg, k
      real(wp) :: r
      r = 10.0_wp + 3.0_wp*sin(0.7_wp*ig + 1.3_wp*jg + k) + 0.05_wp*ig*jg - 0.2_wp*k
   end function hf
   pure function uf(ig, jg, k) result(r)
      integer, intent(in) :: ig, jg, k
      real(wp) :: r
      r = 0.4_wp*cos(0.9_wp*ig - 0.4_wp*jg + 2.0_wp*k) + 0.01_wp*ig
   end function uf
   pure function vf(ig, jg, k) result(r)
      integer, intent(in) :: ig, jg, k
      real(wp) :: r
      r = 0.3_wp*sin(0.5_wp*ig + 1.1_wp*jg*jg*0.1_wp - k) - 0.02_wp*jg
   end function vf
   pure function ff(ig, jg, k, m) result(r)
      integer, intent(in) :: ig, jg, k, m
      real(wp) :: r
      r = 1.0e-3_wp*sin(0.3_wp*m*ig + 0.8_wp*jg + k)
   end function ff
   pure function kf(ig, jg, k) result(r)
      integer, intent(in) :: ig, jg, k
      real(wp) :: r
      r = 0.2_wp*cos(0.6_wp*ig + 0.2_wp*jg) + 0.1_wp*k
   end function kf

   subroutine split(n, p, r, off, cnt)
      integer, intent(in) :: n, p, r
      integer, intent(out) :: off, cnt
      cnt = n/p + merge(1, 0, r < mod(n, p))
      off = r*(n/p) + min(r, mod(n, p))
   end subroutine split

   subroutine local_partials(nxp, nyp, ioff, joff, nxg_use, nyg_use, kt, ct)
      !! Fill a tile (ghosts from the closed form) and return the probe's
      !! rank-local partials.  `nxg_use/nyg_use` are what the routine is told
      !! the global extent is.
      integer, intent(in) :: nxp, nyp, ioff, joff, nxg_use, nyg_use
      type(efp_t), intent(out) :: kt(N_KE), ct(N_SPLIT)
      real(wp), allocatable :: h(:, :, :), u(:, :, :), v(:, :, :), fx(:, :, :), fy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: iu(:, :), iv(:, :)
      integer :: nx, ny, i, j, k, ig, jg
      nx = nxp + 2*NG
      ny = nyp + 2*NG
      allocate (h(nx, ny, NZ), u(nx + 1, ny, NZ), v(nx, ny + 1, NZ), fx(nx + 1, ny, NZ), &
                fy(nx, ny + 1, NZ), kec(nx, ny, NZ), iu(nx + 1, ny), iv(nx, ny + 1))
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx + 1
               ig = i - NG + ioff
               jg = j - NG + joff
               if (i <= nx .and. j <= ny) then
                  h(i, j, k) = hf(ig, jg, k)
                  kec(i, j, k) = kf(ig, jg, k)
               end if
               if (j <= ny) then
                  u(i, j, k) = uf(ig, jg, k)
                  fx(i, j, k) = ff(ig, jg, k, 1)
               end if
               if (i <= nx) then
                  v(i, j, k) = vf(ig, jg, k)
                  fy(i, j, k) = ff(ig, jg, k, 2)
               end if
            end do
         end do
      end do
      do j = 1, ny
         do i = 1, nx + 1
            iu(i, j) = 1.0_wp + 0.01_wp*(i - NG + ioff)
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            iv(i, j) = 1.0_wp + 0.02_wp*(j - NG + joff)
         end do
      end do
      !$acc enter data copyin(h, u, v, fx, fy, kec, iu, iv)
      call ke_regions_efp_local(h, u, v, nx, ny, NZ, NG, nxp, nyp, ioff, joff, &
                                nxg_use, nyg_use, RIM, EAST, kt)
      call coradv_split_efp_local(h, u, v, fx, fy, kec, iu, iv, nx, ny, NZ, NG, nxp, nyp, &
                                  ioff, joff, nxg_use, nyg_use, RIM, EAST, DT, ct)
      !$acc exit data delete(h, u, v, fx, fy, kec, iu, iv)
   end subroutine local_partials

   subroutine serial_ref(kt, ct)
      !! Independent serial EFP sums over the GLOBAL field, masks written out.
      type(efp_t), intent(out) :: kt(N_KE), ct(N_SPLIT)
      integer :: i, j, k, r
      real(wp) :: hw, he, hs, hn, ke, hfc, ft, fk, fp, um, w
      logical :: inr
      kt = efp_t()
      ct = efp_t()
      do k = 1, NZ
         do j = 1, NYG
            do i = 1, NXG
               hw = 0.5_wp*(hf(i - 1, j, k) + hf(i, j, k))
               he = 0.5_wp*(hf(i, j, k) + hf(i + 1, j, k))
               hs = 0.5_wp*(hf(i, j - 1, k) + hf(i, j, k))
               hn = 0.5_wp*(hf(i, j, k) + hf(i, j + 1, k))
               ke = 0.25_wp*(hw*uf(i, j, k)*uf(i, j, k) + he*uf(i + 1, j, k)*uf(i + 1, j, k)) + &
                    0.25_wp*(hs*vf(i, j, k)*vf(i, j, k) + hn*vf(i, j + 1, k)*vf(i, j + 1, k))
               kt(1) = efp_plus(kt(1), efp_from_real(real(ke, real64)))
               if (i - 1 < RIM .or. NXG - i < RIM .or. j - 1 < RIM .or. NYG - j < RIM) &
                  kt(2) = efp_plus(kt(2), efp_from_real(real(ke, real64)))
               if (NXG - i < EAST) kt(3) = efp_plus(kt(3), efp_from_real(real(ke, real64)))
               ! u-face i (west face), region from west cell i-1
               if (i > 1) then
                  hfc = 0.5_wp*(hf(i - 1, j, k) + hf(i, j, k))
                  ft = ff(i, j, k, 1)
                  fk = -(kf(i, j, k) - kf(i - 1, j, k))*(1.0_wp + 0.01_wp*i)
                  fp = ft - fk
                  um = uf(i, j, k) - 0.5_wp*DT*ft
                  do r = 1, 3
                     inr = r == 1
                     if (r == 2) inr = (i - 2 < RIM .or. NXG - (i - 1) < RIM .or. j - 1 < RIM .or. NYG - j < RIM)
                     if (r == 3) inr = (NXG - (i - 1) < EAST)
                     if (inr) then
                        w = hfc*um*DT*fp
                        ct(r) = efp_plus(ct(r), efp_from_real(real(w, real64)))
                        w = hfc*um*DT*fk
                        ct(3 + r) = efp_plus(ct(3 + r), efp_from_real(real(w, real64)))
                     end if
                  end do
               end if
               if (j > 1) then
                  hfc = 0.5_wp*(hf(i, j - 1, k) + hf(i, j, k))
                  ft = ff(i, j, k, 2)
                  fk = -(kf(i, j, k) - kf(i, j - 1, k))*(1.0_wp + 0.02_wp*j)
                  fp = ft - fk
                  um = vf(i, j, k) - 0.5_wp*DT*ft
                  do r = 1, 3
                     inr = r == 1
                     if (r == 2) inr = (i - 1 < RIM .or. NXG - i < RIM .or. j - 2 < RIM .or. NYG - (j - 1) < RIM)
                     if (r == 3) inr = (NXG - i < EAST)
                     if (inr) then
                        w = hfc*um*DT*fp
                        ct(r) = efp_plus(ct(r), efp_from_real(real(w, real64)))
                        w = hfc*um*DT*fk
                        ct(3 + r) = efp_plus(ct(3 + r), efp_from_real(real(w, real64)))
                     end if
                  end do
               end if
            end do
         end do
      end do
   end subroutine serial_ref

   subroutine run_case(px, py)
      integer, intent(in) :: px, py
      type(efp_t) :: kt(N_KE), ct(N_SPLIT), kg(N_KE), cg(N_SPLIT), kr(N_KE), cr(N_SPLIT), k1(N_KE), c1(N_SPLIT), ko(N_KE), co(N_SPLIT), kog(N_KE), cog(N_SPLIT)
      integer :: rx, ry, ioff, joff, nxp, nyp, m
      logical :: ok_new, any_diff

      rx = mod(rank, px)
      ry = rank/px
      call split(NXG, px, rx, ioff, nxp)
      call split(NYG, py, ry, joff, nyp)
      call serial_ref(kr, cr)

      ! Fixed behaviour: global offsets + global extents.
      call local_partials(nxp, nyp, ioff, joff, NXG, NYG, kt, ct)
      call halo_allreduce_efp_list(kt, kg, N_KE)
      call halo_allreduce_efp_list(ct, cg, N_SPLIT)
      ! Bitwise against the routine run over the whole domain as ONE tile
      ! (what the 1-rank run computes), and to round-off against the
      ! independent serial formula (masks / face ownership: an error there is
      ! O(1), so 1e-12 relative cannot hide one; the last digits differ only
      ! through compiler FMA contraction of the same expressions).
      call local_partials(NXG, NYG, 0, 0, NXG, NYG, k1, c1)
      ok_new = .true.
      do m = 1, N_KE
         if (efp_to_real(kg(m)) /= efp_to_real(k1(m))) ok_new = .false.
         if (abs(efp_to_real(kg(m)) - efp_to_real(kr(m))) > 1.0e-12_real64*abs(efp_to_real(kr(m)))) ok_new = .false.
      end do
      do m = 1, N_SPLIT
         if (efp_to_real(cg(m)) /= efp_to_real(c1(m))) ok_new = .false.
         if (abs(efp_to_real(cg(m)) - efp_to_real(cr(m))) > 1.0e-12_real64*abs(efp_to_real(cr(m)))) ok_new = .false.
      end do
      ntry = ntry + 1
      if (.not. ok_new) then
         nfail = nfail + 1
         if (rank == 0) write (*, "(a,i0,a,i0,a)") "FAIL: global result differs on ", px, "x", py, " tiles"
      end if

      ! Control: the pre-fix semantics (tile-local indices and extents).
      call local_partials(nxp, nyp, 0, 0, nxp, nyp, ko, co)
      call halo_allreduce_efp_list(ko, kog, N_KE)
      call halo_allreduce_efp_list(co, cog, N_SPLIT)
      any_diff = .false.
      do m = 1, N_KE
         if (efp_to_real(kog(m)) /= efp_to_real(k1(m))) any_diff = .true.
      end do
      do m = 1, N_SPLIT
         if (efp_to_real(cog(m)) /= efp_to_real(c1(m))) any_diff = .true.
      end do
      if (px*py > 1) then
         ntry = ntry + 1
         if (.not. any_diff) then
            nfail = nfail + 1
            if (rank == 0) write (*, "(a)") "FAIL: control (tile-local masks) matched; test has no teeth"
         end if
      end if
      if (rank == 0) write (*, "(a,i0,a,i0,a,l1,a,l1)") "  ", px, "x", py, " global_match=", ok_new, &
         " tilelocal_differs=", any_diff
   end subroutine run_case

end program test_ocean_ke_probe_mpi
#endif
