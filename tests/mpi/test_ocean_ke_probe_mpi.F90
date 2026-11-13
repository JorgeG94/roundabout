!! MPI test of the KE-attribution probe's reductions (`&ocean_debug_nml ke_attr`).
!! Run with 1, 2 or 4 ranks: mpirun -np N ./test_ocean_ke_probe_mpi
!!
!! A non-uniform, non-symmetric analytic C-grid state (a closed form of the
!! GLOBAL cell index, ghosts included) is cut into tiles on every px x py
!! factorisation of the rank count (uneven splits).  Each rank calls the
!! probe's rank-local EFP routines, the partials are combined with
!! `halo_allreduce_efp_list`, and the result is compared bitwise against an
!! independent serial EFP sum over the global field with plainly written
!! region masks.  Control: the field is filled with the TRUE global offsets,
!! but the probe is told tile-local offsets and extents (the pre-fix masks);
!! that must DIFFER on more than one rank, so the test has teeth.
!!
!! Three boundary modes: closed, periodic in x, doubly periodic.  On a
!! periodic axis the state is the wrapped closed form (h, kec ghosts are wrap
!! copies of the true neighbours), and EVERY non-owned u/v/fx/fy value is
!! poisoned -- including the east (north) slot nx_g+1, which on a real
!! periodic domain is a stale wrap copy of face 1.  The serial reference uses
!! periodic neighbours (face 1 is the east face of cell nx_g).
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
   integer :: rank, nprocs, px, py, nfail, ntry, mode
   logical :: per_x = .false., per_y = .false.

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   nfail = 0
   ntry = 0

   do px = 1, nprocs
      if (mod(nprocs, px) /= 0) cycle
      py = nprocs/px
      do mode = 0, 2
         per_x = mode >= 1
         per_y = mode >= 2
         call run_case(px, py)
      end do
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

   pure function wr(i, n, per) result(r)
      !! Wrap a global index into 1..n on a periodic axis (identity otherwise).
      integer, intent(in) :: i, n
      logical, intent(in) :: per
      integer :: r
      r = merge(modulo(i - 1, n) + 1, i, per)
   end function wr

   subroutine split(n, p, r, off, cnt)
      integer, intent(in) :: n, p, r
      integer, intent(out) :: off, cnt
      cnt = n/p + merge(1, 0, r < mod(n, p))
      off = r*(n/p) + min(r, mod(n, p))
   end subroutine split

   subroutine local_partials(nxp, nyp, ioff, joff, pnx, pny, poff, qoff, kt, ct)
      !! Fill a tile (TRUE global offsets `ioff/joff`, ghosts from the
      !! closed form, wrapped on periodic axes) and return the probe's
      !! rank-local partials.  `pnx/pny/poff/qoff` are the global extents
      !! and tile offsets the PROBE is told (the true ones, except in the
      !! control).
      integer, intent(in) :: nxp, nyp, ioff, joff, pnx, pny, poff, qoff
      type(efp_t), intent(out) :: kt(N_KE), ct(N_SPLIT)
      real(wp), allocatable :: h(:, :, :), u(:, :, :), v(:, :, :), fx(:, :, :), fy(:, :, :), kec(:, :, :)
      real(wp), allocatable :: iu(:, :), iv(:, :)
      integer :: nx, ny, i, j, k, ig, jg
      logical :: own_u, own_v
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
                  h(i, j, k) = hf(wr(ig, NXG, per_x), wr(jg, NYG, per_y), k)
                  kec(i, j, k) = kf(wr(ig, NXG, per_x), wr(jg, NYG, per_y), k)
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
      ! POISON every u/v/fx/fy value this rank does not own (ghost faces are
      ! stale at the production stage-close sample): owned = west faces of the
      ! owned cells, plus the global east / north boundary face on a
      ! NON-periodic axis (periodic: that slot is a stale wrap copy).  h and kec
      ! ghosts stay exact (the probe documents that it relies on them).
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx + 1
               own_u = (j >= NG + 1 .and. j <= NG + nyp) .and. (i >= NG + 1 .and. &
                                                 (i <= NG + nxp .or. (i == NG + nxp + 1 .and. ioff + nxp == NXG .and. .not. per_x)))
               own_v = (i >= NG + 1 .and. i <= NG + nxp) .and. (j >= NG + 1 .and. &
                                                 (j <= NG + nyp .or. (j == NG + nyp + 1 .and. joff + nyp == NYG .and. .not. per_y)))
               if (.not. own_u .and. j <= ny) then
                  u(i, j, k) = 1.0e3_wp*(rank + 1)*(i + j + k)
                  fx(i, j, k) = 7.0e2_wp*(rank + 1)*(i - j + k)
               end if
               if (.not. own_v .and. i <= nx) then
                  v(i, j, k) = -1.0e3_wp*(rank + 1)*(i + j + k)
                  fy(i, j, k) = 7.0e2_wp*(rank + 1)*(i - j + k)
               end if
            end do
         end do
      end do
      !$acc enter data copyin(h, u, v, fx, fy, kec, iu, iv)
      ! (copyin happens after every host fill, poison included)
      call ke_regions_efp_local(h, u, v, nx, ny, NZ, NG, nxp, nyp, poff, qoff, &
                                pnx, pny, RIM, EAST, per_x, per_y, kt)
      call coradv_split_efp_local(h, u, v, fx, fy, kec, iu, iv, nx, ny, NZ, NG, nxp, nyp, &
                                  poff, qoff, pnx, pny, RIM, EAST, DT, per_x, per_y, ct)
      !$acc exit data delete(h, u, v, fx, fy, kec, iu, iv)
   end subroutine local_partials

   subroutine serial_ref(kt, ct)
      !! Independent serial EFP sums over the GLOBAL field, masks written out.
      type(efp_t), intent(out) :: kt(N_KE), ct(N_SPLIT)
      integer :: i, j, k, r, iw, ie, js, jn
      real(wp) :: hw, he, hs, hn, ke, hfc, ft, fk, fp, um, w
      logical :: inr
      kt = efp_t()
      ct = efp_t()
      do k = 1, NZ
         do j = 1, NYG
            do i = 1, NXG
               iw = wr(i - 1, NXG, per_x)
               ie = wr(i + 1, NXG, per_x)
               js = wr(j - 1, NYG, per_y)
               jn = wr(j + 1, NYG, per_y)
               hw = 0.5_wp*(hf(iw, j, k) + hf(i, j, k))
               he = 0.5_wp*(hf(i, j, k) + hf(ie, j, k))
               hs = 0.5_wp*(hf(i, js, k) + hf(i, j, k))
               hn = 0.5_wp*(hf(i, j, k) + hf(i, jn, k))
               ke = 0.25_wp*(hw*uf(i, j, k)*uf(i, j, k) + he*uf(ie, j, k)*uf(ie, j, k)) + &
                    0.25_wp*(hs*vf(i, j, k)*vf(i, j, k) + hn*vf(i, jn, k)*vf(i, jn, k))
               kt(1) = efp_plus(kt(1), efp_from_real(real(ke, real64)))
               if (i - 1 < RIM .or. NXG - i < RIM .or. j - 1 < RIM .or. NYG - j < RIM) &
                  kt(2) = efp_plus(kt(2), efp_from_real(real(ke, real64)))
               if (NXG - i < EAST) kt(3) = efp_plus(kt(3), efp_from_real(real(ke, real64)))
               ! u-face i (west face), region from west cell i-1
               if (i > 1 .or. per_x) then
                  hfc = 0.5_wp*(hf(iw, j, k) + hf(i, j, k))
                  ft = ff(i, j, k, 1)
                  fk = -(kf(i, j, k) - kf(iw, j, k))*(1.0_wp + 0.01_wp*i)
                  fp = ft - fk
                  um = uf(i, j, k) - 0.5_wp*DT*ft
                  do r = 1, 3
                     inr = r == 1
                     if (r == 2) inr = (iw - 1 < RIM .or. NXG - iw < RIM .or. j - 1 < RIM .or. NYG - j < RIM)
                     if (r == 3) inr = (NXG - iw < EAST)
                     if (inr) then
                        w = hfc*um*DT*fp
                        ct(r) = efp_plus(ct(r), efp_from_real(real(w, real64)))
                        w = hfc*um*DT*fk
                        ct(3 + r) = efp_plus(ct(3 + r), efp_from_real(real(w, real64)))
                     end if
                  end do
               end if
               if (j > 1 .or. per_y) then
                  hfc = 0.5_wp*(hf(i, js, k) + hf(i, j, k))
                  ft = ff(i, j, k, 2)
                  fk = -(kf(i, j, k) - kf(i, js, k))*(1.0_wp + 0.02_wp*j)
                  fp = ft - fk
                  um = vf(i, j, k) - 0.5_wp*DT*ft
                  do r = 1, 3
                     inr = r == 1
                     if (r == 2) inr = (i - 1 < RIM .or. NXG - i < RIM .or. js - 1 < RIM .or. NYG - js < RIM)
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
      call local_partials(nxp, nyp, ioff, joff, NXG, NYG, ioff, joff, kt, ct)
      call halo_allreduce_efp_list(kt, kg, N_KE)
      call halo_allreduce_efp_list(ct, cg, N_SPLIT)
      ! Bitwise against the routine run over the whole domain as ONE tile
      ! (what the 1-rank run computes), and to round-off against the
      ! independent serial formula (masks / face ownership: an error there is
      ! O(1), so 1e-12 relative cannot hide one; the last digits differ only
      ! through compiler FMA contraction of the same expressions).
      call local_partials(NXG, NYG, 0, 0, NXG, NYG, 0, 0, k1, c1)
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
        if (rank == 0) write (*, "(a,i0,a,i0,a,l1,l1)") "FAIL: global result differs on ", px, "x", py, " tiles, per=", per_x, per_y
      end if

      ! Control: the field holds the TRUE global state; the probe is told
      ! tile-local offsets and extents (the pre-fix masks).
      call local_partials(nxp, nyp, ioff, joff, nxp, nyp, 0, 0, ko, co)
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
     if (rank == 0) write (*, "(a,l1,l1,a,i0,a,i0,a,l1,a,l1)") "  per=", per_x, per_y, " ", px, "x", py, " global_match=", ok_new, &
         " tilelocal_differs=", any_diff
   end subroutine run_case

end program test_ocean_ke_probe_mpi
#endif
