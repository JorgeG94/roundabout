!! `ice_cice_state_t` (CICE-port ice state slot): lifecycle, device residency,
!! the aggregate sweep, the uniform seed, restart registration + round trip and
!! the error field.
!!
!! Device tests map the slot with `enter_data`, run `do concurrent` kernels that
!! read and write every array (canary values) and compare CELL BY CELL with a
!! host reference on states whose cells all differ; component arrays only are
!! pulled back with `!$acc update self`.  On host builds the directives are
!! inert and the same code checks the kernels against the references.
!! `ice_update_state` is a plain `do concurrent` body with the category sum
!! inline (no kernel call, so no loop index is passed to a routine).
module test_cice_state
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data, &
                              ocean_state_seed_from_cfg, ocean_state_build_restart_registry
   use rdb_state, only: register_default_tracers
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, &
                                   ocean_restart_metadata_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_constants_t, ICE_CONSTANTS_ICEPACK
   use rdb_ice_itd_shift, only: itd_aggregate_area
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_NERR, ice_err_summary_t, ice_err_reduce, &
                             ICE_ERR_FATAL_MASK
   use rdb_ice_seed, only: ICE_SEED_OK, ICE_SEED_BAD_ARG
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_update_state, ice_state_seed_uniform, &
                                 ICE_STATE_VERSION, ICE_STATE_INIT_BAD_PARAMS
   implicit none
   private

   public :: collect_cice_state_tests

   integer, parameter :: NXP = 7, NYP = 5, NG = 2, NCAT = 5, NILYR = 7
   integer, parameter :: NX = NXP + 2*NG, NY = NYP + 2*NG

contains

   subroutine collect_cice_state_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_destroy_bytes", test_init_destroy_bytes), &
                  new_unittest("init_refuses_bad_params", test_init_refuses), &
                  new_unittest("device_canary_every_array", test_device_canary), &
                  new_unittest("update_state_random_cell_by_cell", test_update_state), &
                  new_unittest("update_state_nonfinite_flags_cell", test_update_state_nan), &
                  new_unittest("seed_uniform_closed_forms", test_seed), &
                  new_unittest("seed_uniform_on_device", test_seed_device), &
                  new_unittest("seed_uniform_refusals", test_seed_refusals), &
                  new_unittest("restart_registration", test_restart_registration), &
                  new_unittest("restart_round_trip", test_restart_round_trip), &
                  new_unittest("sis2_default_untouched", test_sis2_untouched) &
                  ]
   end subroutine collect_cice_state_tests

   subroutine make_state(st, ierr, ncat_in, nilyr_in)
      type(ice_cice_state_t), intent(inout) :: st
      integer, intent(out) :: ierr
      integer, intent(in), optional :: ncat_in, nilyr_in
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      integer :: nc, nl

      nc = NCAT
      nl = NILYR
      if (present(ncat_in)) nc = ncat_in
      if (present(nilyr_in)) nl = nilyr_in
      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, nc, nl, 1)
      call st%init(g, p, ierr)
   end subroutine make_state

   pure function tol_for(scale) result(t)
      !! Device vs host: FMA contraction may differ in the last bits.
      real(wp), intent(in) :: scale
      real(wp) :: t
      t = 1.0e-13_wp*max(abs(scale), 1.0_wp)
   end function tol_for

   subroutine same(error, got, want, what)
      !! Finite-checked comparison.
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: got, want
      character(len=*), intent(in) :: what
      if (allocated(error)) return
      if (.not. ieee_is_finite(got)) then
         call test_failed(error, what//": non-finite value")
         return
      end if
      if (abs(got - want) > tol_for(want)) call test_failed(error, what//": mismatch")
   end subroutine same

   ! ------------------------------------------------------------------
   subroutine test_init_destroy_bytes(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      integer :: ierr, nxy, n
      integer(int64) :: want

      call make_state(st, ierr)
      call check(error, ierr, 0, "init ok")
      if (allocated(error)) return
      call check(error, st%is_init, "is_init")
      call check(error, size(st%aicen, 1) == NX .and. size(st%aicen, 2) == NY .and. size(st%aicen, 3) == NCAT, &
                 "aicen shape (ghosted, ncat)")
      call check(error, size(st%qicen, 4) == NILYR .and. size(st%sicen, 4) == NILYR, "layer extent")
      call check(error, size(st%qsnon, 4) == 1, "snow layer extent")
      call check(error, lbound(st%hin_max, 1) == 0 .and. ubound(st%hin_max, 1) == NCAT, "hin_max(0:ncat)")
      if (allocated(error)) return
      call same(error, st%hin_max(0), 0.0_wp, "hin_max(0)")
      call same(error, st%hin_max(1), 0.6_wp, "hin_max(1) = d1 = 3/ncat")
      call same(error, st%hin_max(NCAT), 999.9_wp, "hin_max(ncat) = top bound")
      call check(error, all(st%aice0 == 1.0_wp) .and. all(st%aicen == 0.0_wp) .and. all(st%ice_err == 0_int32), &
                 "inited as open water, no errors")
      nxy = NX*NY
      n = NCAT
      ! per cell: aicen,vicen,vsnon,Tsfcn + 4 start-of-step copies = 8n; qicen,sicen = 2 n nilyr;
      ! qsnon = n; aice0,aice,vice,vsno,aice_init = 5 reals; ice_err = one int32.
      want = int(nxy, int64)*8_int64*(8_int64*n + 2_int64*n*NILYR + n + 5_int64) &
             + 8_int64*(n + 1_int64) + int(nxy, int64)*4_int64
      call check(error, st%bytes() == want, "bytes() counts every array (closed form)")
      call st%destroy()
      call check(error,.not. st%is_init .and. .not. allocated(st%aicen) .and. st%bytes() == 0_int64, &
                 "destroy releases everything")
   end subroutine test_init_destroy_bytes

   subroutine test_init_refuses(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      integer :: ierr

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 9, 4, 1)   ! ncat > ICE_NCAT_MAX
      call st%init(g, p, ierr)
      call check(error, ierr, ICE_STATE_INIT_BAD_PARAMS, "ncat beyond the compile-time bound is refused")
      call check(error,.not. st%is_init .and. .not. allocated(st%aicen), "nothing allocated on refusal")
   end subroutine test_init_refuses

   ! ------------------------------------------------------------------
   subroutine canary_3d(nx_, ny_, nc, c, hin_max, aicen, vicen, vsnon, Tsfcn, ai, vi, vs, ts)
      integer, intent(in) :: nx_, ny_, nc
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: hin_max(0:nc)
      real(wp), intent(inout) :: aicen(nx_, ny_, nc)
      real(wp), intent(inout) :: vicen(nx_, ny_, nc), vsnon(nx_, ny_, nc), Tsfcn(nx_, ny_, nc)
      real(wp), intent(inout) :: ai(nx_, ny_, nc), vi(nx_, ny_, nc), vs(nx_, ny_, nc), ts(nx_, ny_, nc)
      integer :: i, j, n

      !$acc data present(c, hin_max, aicen, vicen, vsnon, Tsfcn, ai, vi, vs, ts)
      do concurrent(n=1:nc, j=1:ny_, i=1:nx_)
         ai(i, j, n) = aicen(i, j, n)
         vi(i, j, n) = 2.0_wp*aicen(i, j, n)
         vs(i, j, n) = vi(i, j, n) + 1.0_wp
         ts(i, j, n) = 3.0_wp*vs(i, j, n)
         vicen(i, j, n) = ts(i, j, n)
         vsnon(i, j, n) = ai(i, j, n) + c%rhoi + hin_max(n)
         Tsfcn(i, j, n) = Tsfcn(i, j, n) + 1.0_wp
         aicen(i, j, n) = aicen(i, j, n) + 4.0_wp
      end do
      !$acc end data
   end subroutine canary_3d

   subroutine canary_4d(nx_, ny_, nc, nl, ns, qicen, sicen, qsnon)
      integer, intent(in) :: nx_, ny_, nc, nl, ns
      real(wp), intent(inout) :: qicen(nx_, ny_, nc, nl), sicen(nx_, ny_, nc, nl)
      real(wp), intent(inout) :: qsnon(nx_, ny_, nc, ns)
      integer :: i, j, n, k

      !$acc data present(qicen, sicen, qsnon)
      do concurrent(k=1:nl, n=1:nc, j=1:ny_, i=1:nx_)
         qicen(i, j, n, k) = qicen(i, j, n, k) + 2.0_wp*sicen(i, j, n, k)
         sicen(i, j, n, k) = sicen(i, j, n, k) + qsnon(i, j, n, 1)
      end do
      do concurrent(n=1:nc, j=1:ny_, i=1:nx_)
         qsnon(i, j, n, 1) = qsnon(i, j, n, 1) + 5.0_wp
      end do
      !$acc end data
   end subroutine canary_4d

   subroutine canary_2d(nx_, ny_, aice0, aice, vice, vsno, a_init, ice_err)
      integer, intent(in) :: nx_, ny_
      real(wp), intent(inout) :: aice0(nx_, ny_), aice(nx_, ny_), vice(nx_, ny_), vsno(nx_, ny_)
      real(wp), intent(inout) :: a_init(nx_, ny_)
      integer(int32), intent(inout) :: ice_err(nx_, ny_)
      integer :: i, j

      !$acc data present(aice0, aice, vice, vsno, a_init, ice_err)
      do concurrent(j=1:ny_, i=1:nx_)
         a_init(i, j) = aice0(i, j) + 1.0_wp
         aice(i, j) = aice0(i, j) + 10.0_wp
         vice(i, j) = vice(i, j) + real(ice_err(i, j), wp)
         vsno(i, j) = vsno(i, j)*2.0_wp
         ice_err(i, j) = ice_err(i, j) + 1_int32
         aice0(i, j) = aice0(i, j)*0.5_wp
      end do
      !$acc end data
   end subroutine canary_2d

   subroutine test_device_canary(error)
      !! Every array is read and written on the device with values that differ per
      !! cell; the persistent ones come back through component `update self`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      real(wp), allocatable :: a0(:, :, :), q0(:, :, :, :), s0(:, :, :, :), qs0(:, :, :, :)
      real(wp), allocatable :: t0(:, :, :), v0(:, :), ai0(:, :), vs0(:, :), ao0(:, :)
      integer(int32), allocatable :: e0(:, :)
      integer :: ierr, i, j, n, k
      real(wp) :: tsn, vi_ref, vs_ref

      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               st%aicen(i, j, n) = 0.5_wp*real(i + 3*j + 7*n, wp)
               st%Tsfcn(i, j, n) = 0.25_wp*real(2*i - j + n, wp)
               do k = 1, NILYR
                  st%qicen(i, j, n, k) = real(i + j*n + k, wp)
                  st%sicen(i, j, n, k) = 0.5_wp*real(i - j + n*k, wp)
               end do
               st%qsnon(i, j, n, 1) = real(i*j + n, wp)
            end do
         end do
      end do
      do j = 1, NY
         do i = 1, NX
            st%aice0(i, j) = 0.125_wp*real(i + 2*j, wp)
            st%vice(i, j) = real(i - j, wp)
            st%vsno(i, j) = 0.5_wp*real(i*j, wp)
            st%ice_err(i, j) = int(i + 11*j, int32)
         end do
      end do
      a0 = st%aicen
      q0 = st%qicen
      s0 = st%sicen
      qs0 = st%qsnon
      t0 = st%Tsfcn
      v0 = st%vice
      vs0 = st%vsno
      ai0 = st%aice0
      e0 = st%ice_err

      call st%enter_data()
      call canary_3d(NX, NY, NCAT, st%p%c, st%hin_max, st%aicen, st%vicen, st%vsnon, st%Tsfcn, &
                     st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init)
      call canary_4d(NX, NY, NCAT, NILYR, 1, st%qicen, st%sicen, st%qsnon)
      call canary_2d(NX, NY, st%aice0, st%aice, st%vice, st%vsno, st%aice_init, st%ice_err)
      ! Scratch and `ice_err` are deleted by exit_data: pull them first.  The
      ! persistent state and the aggregates must come back through exit_data's
      ! copyout alone (no `update self` of them here).
      !$acc update self(st%ice_err)
      !$acc update self(st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init, st%aice_init)
      call st%exit_data()
      call check(error,.not. st%mapped, "exit_data unmaps")

      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               vi_ref = 3.0_wp*(2.0_wp*a0(i, j, n) + 1.0_wp)
               vs_ref = a0(i, j, n) + st%p%c%rhoi + st%hin_max(n)
               call same(error, st%aicen(i, j, n), a0(i, j, n) + 4.0_wp, "aicen (device-written, copied out)")
               call same(error, st%aicen_init(i, j, n), a0(i, j, n), "aicen_init")
               call same(error, st%vicen_init(i, j, n), 2.0_wp*a0(i, j, n), "vicen_init")
               call same(error, st%vsnon_init(i, j, n), 2.0_wp*a0(i, j, n) + 1.0_wp, "vsnon_init")
               call same(error, st%Tsfcn_init(i, j, n), vi_ref, "Tsfcn_init")
               call same(error, st%vicen(i, j, n), vi_ref, "vicen")
               call same(error, st%vsnon(i, j, n), vs_ref, "vsnon (reads device p and hin_max)")
               tsn = t0(i, j, n) + 1.0_wp
               call same(error, st%Tsfcn(i, j, n), tsn, "Tsfcn")
               do k = 1, NILYR
                  call same(error, st%qicen(i, j, n, k), q0(i, j, n, k) + 2.0_wp*s0(i, j, n, k), "qicen")
                  call same(error, st%sicen(i, j, n, k), s0(i, j, n, k) + qs0(i, j, n, 1), "sicen")
               end do
               call same(error, st%qsnon(i, j, n, 1), qs0(i, j, n, 1) + 5.0_wp, "qsnon (device-written, copied out)")
            end do
         end do
      end do
      do j = 1, NY
         do i = 1, NX
            call same(error, st%aice_init(i, j), ai0(i, j) + 1.0_wp, "aice_init")
            call same(error, st%aice(i, j), ai0(i, j) + 10.0_wp, "aice")
            call same(error, st%vice(i, j), v0(i, j) + real(e0(i, j), wp), "vice")
            call same(error, st%vsno(i, j), 2.0_wp*vs0(i, j), "vsno")
            call same(error, st%aice0(i, j), 0.5_wp*ai0(i, j), "aice0")
            if (st%ice_err(i, j) /= e0(i, j) + 1_int32) call test_failed(error, "ice_err canary")
            if (allocated(error)) return
         end do
      end do
      call st%exit_data()
      call check(error,.not. st%mapped, "exit_data unmaps")
      call st%destroy()
   end subroutine test_device_canary

   ! ------------------------------------------------------------------
   subroutine fill_random(st, seed, empty_frac)
      !! Per-cell different area/volume state; a fraction of cells is empty.
      type(ice_cice_state_t), intent(inout) :: st
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: empty_frac
      integer(int64) :: s
      integer :: i, j, n
      real(wp) :: u

      s = seed
      do j = 1, NY
         do i = 1, NX
            do n = 1, NCAT
               s = lcg_next(s)
               u = lcg_unit(s)
               st%aicen(i, j, n) = 0.0_wp
               st%vicen(i, j, n) = 0.0_wp
               st%vsnon(i, j, n) = 0.0_wp
               if (u > empty_frac) then
                  st%aicen(i, j, n) = 0.2_wp*u
                  st%vicen(i, j, n) = 3.0_wp*u*u
                  st%vsnon(i, j, n) = 0.4_wp*(1.0_wp - u)
               end if
            end do
         end do
      end do
   end subroutine fill_random

   subroutine test_update_state(error)
      !! `ice_update_state` against a brute-force host sum, per cell, over the owned
      !! range (ghosts must stay untouched), then again with ghosts included.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      integer :: ierr, i, j, n, is, ie, js, je
      real(wp) :: a, v, s, a_chk, a0_chk, differ_min, differ_max
      logical :: owned

      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_random(st, 123456789_int64, 0.3_wp)
      call st%enter_data()
      call ice_update_state(st, .false.)
      !$acc update self(st%aice, st%vice, st%vsno, st%aice0, st%ice_err)
      is = NG + 1
      ie = NG + NXP
      js = NG + 1
      je = NG + NYP
      differ_min = huge(1.0_wp)
      differ_max = -huge(1.0_wp)
      do j = 1, NY
         do i = 1, NX
            owned = i >= is .and. i <= ie .and. j >= js .and. j <= je
            a = 0.0_wp
            v = 0.0_wp
            s = 0.0_wp
            do n = 1, NCAT
               a = a + st%aicen(i, j, n)
               v = v + st%vicen(i, j, n)
               s = s + st%vsnon(i, j, n)
            end do
            if (owned) then
               call same(error, st%aice(i, j), a, "aice owned")
               call same(error, st%vice(i, j), v, "vice owned")
               call same(error, st%vsno(i, j), s, "vsno owned")
               call same(error, st%aice0(i, j), max(1.0_wp - a, 0.0_wp), "aice0 owned")
               call itd_aggregate_area(NCAT, st%aicen(i, j, :), a_chk, a0_chk)
               call same(error, st%aice(i, j), a_chk, "aice vs itd_aggregate_area")
               call same(error, st%aice0(i, j), a0_chk, "aice0 vs itd_aggregate_area")
               differ_min = min(differ_min, a)
               differ_max = max(differ_max, a)
            else
               call same(error, st%aice(i, j), 0.0_wp, "ghost aice untouched")
               call same(error, st%aice0(i, j), 1.0_wp, "ghost aice0 untouched")
            end if
            if (st%ice_err(i, j) /= 0_int32) call test_failed(error, "no error bits on a finite state")
            if (allocated(error)) return
         end do
      end do
      call check(error, differ_max - differ_min > 0.05_wp, "cells genuinely differ (test is not degenerate)")
      if (allocated(error)) return
      call ice_update_state(st, .true.)
      !$acc update self(st%aice, st%vice, st%vsno, st%aice0)
      do j = 1, NY
         do i = 1, NX
            a = 0.0_wp
            v = 0.0_wp
            do n = 1, NCAT
               a = a + st%aicen(i, j, n)
               v = v + st%vicen(i, j, n)
            end do
            call same(error, st%aice(i, j), a, "aice with ghosts")
            call same(error, st%vice(i, j), v, "vice with ghosts")
            call same(error, st%aice0(i, j), max(1.0_wp - a, 0.0_wp), "aice0 with ghosts")
            if (allocated(error)) return
         end do
      end do
      call st%exit_data()
      call st%destroy()
   end subroutine test_update_state

   subroutine test_update_state_nan(error)
      !! A non-finite category area poisons its own cell (no clamp launders it),
      !! sets the error bit there only, and `ice_err_reduce` over the OWNED range
      !! finds it (and ignores a poisoned ghost).
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(ice_err_summary_t) :: sm
      integer :: ierr, i, j, ib, jb, nbad
      real(wp) :: nan

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call fill_random(st, 99_int64, 0.2_wp)
      ib = NG + 2
      jb = NG + 3
      st%aicen(ib, jb, 2) = nan
      st%aicen(1, 1, 1) = nan      ! ghost: outside the owned sweep and the reduction
      call st%enter_data()
      call ice_update_state(st, .false.)
      !$acc update self(st%aice, st%aice0, st%ice_err)
      nbad = 0
      do j = 1, NY
         do i = 1, NX
            if (st%ice_err(i, j) /= 0_int32) then
               nbad = nbad + 1
               call check(error, i == ib .and. j == jb, "error bit only on the poisoned owned cell")
            end if
         end do
      end do
      call check(error, nbad, 1, "exactly one flagged cell")
      call check(error,.not. ieee_is_finite(st%aice(ib, jb)) .and. .not. ieee_is_finite(st%aice0(ib, jb)), &
                 "aice and aice0 stay non-finite")
      if (allocated(error)) return
      call ice_err_reduce(st%ice_err, NX, NY, NG + 1, NG + NXP, NG + 1, NG + NYP, &
                          st%i_offset_global, st%j_offset_global, st%nx_global, ICE_ERR_FATAL_MASK, sm)
      call check(error, iand(sm%mask, ICE_ERR_NONFINITE) /= 0_int32, "reduction sees the bit")
      call check(error, int(sm%count(5)), 1, "one cell carries bit 5 (NONFINITE)")
      call check(error, sm%first_idx == int(jb - NG - 1, int64)*NXP + int(ib - NG, int64), &
                 "first offender global index")
      call st%exit_data()
      call st%destroy()
   end subroutine test_update_state_nan

   ! ------------------------------------------------------------------
   subroutine seed_check(error, st, aice_in, hi, hs, Tsfc, Tf, ncell)
      !! Closed forms (Icepack constants written here, not taken from the module).
      type(error_type), allocatable, intent(inout) :: error
      type(ice_cice_state_t), intent(in) :: st
      real(wp), intent(in) :: aice_in, hi, hs, Tsfc, Tf
      integer, intent(in) :: ncell
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: RHOI = 917.0_wp, RHOS = 330.0_wp, CPI = 2106.0_wp, CPO = 4218.0_wp
      real(wp), parameter :: LF = 2.835e6_wp - 2.501e6_wp, DT = 0.054_wp
      integer :: n, k, i, j
      real(wp) :: zn, sal, Tmlt, T, q, Ts, vol

      vol = 0.0_wp
      do n = 1, NCAT
         if (n == ncell) then
            Ts = min(Tsfc, 0.0_wp)
         else
            Ts = min(Tf, 0.0_wp)
         end if
         do j = 1, NY, 3
            do i = 1, NX, 4
               if (n == ncell) then
                  call same(error, st%aicen(i, j, n), aice_in, "aicen slab")
                  call same(error, st%vicen(i, j, n), aice_in*hi, "vicen slab")
                  call same(error, st%vsnon(i, j, n), aice_in*hs, "vsnon slab")
               else
                  call same(error, st%aicen(i, j, n) + 1.0_wp, 1.0_wp, "aicen empty")
                  call same(error, st%vicen(i, j, n) + 1.0_wp, 1.0_wp, "vicen empty")
               end if
               call same(error, st%Tsfcn(i, j, n), Ts, "Tsfcn")
               call same(error, st%qsnon(i, j, n, 1), -RHOS*(LF - CPI*min(0.0_wp, Ts)), "qsnon")
               do k = 1, NILYR
                  zn = (real(k, wp) - 0.5_wp)/real(NILYR, wp)
                  sal = max(0.5_wp*3.2_wp*(1.0_wp - cos(PI*zn**(0.407_wp/(0.573_wp + zn)))), 0.1_wp)
                  Tmlt = -DT*sal
                  if (n == ncell) then
                     T = min(Ts + (Tf - Ts)*zn, Tmlt)
                  else
                     T = min(Ts, Tmlt)
                  end if
                  q = -RHOI*(CPI*(Tmlt - T) + LF*(1.0_wp - Tmlt/T) - CPO*Tmlt)
                  call same(error, st%sicen(i, j, n, k), sal, "sicen profile")
                  call same(error, st%qicen(i, j, n, k), q, "qicen profile")
               end do
               if (allocated(error)) return
            end do
         end do
         vol = vol + st%vicen(1, 1, n)
      end do
      call same(error, vol, aice_in*hi, "requested volume conserved")
      call same(error, st%aice(3, 3), aice_in, "aice")
      call same(error, st%vice(3, 3), aice_in*hi, "vice")
      call same(error, st%aice0(3, 3), 1.0_wp - aice_in, "aice0")
   end subroutine seed_check

   subroutine test_seed(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      integer :: ierr

      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      ! kcatbound = 1, ncat = 5: bounds n*(0.6 + 0.1*(n-1)) = 0.6, 1.4, 2.4, 3.6, ...; 1.5 m -> category 3.
      call ice_state_seed_uniform(st, 0.8_wp, 1.5_wp, 0.2_wp, -10.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_OK, "seed ok")
      if (allocated(error)) return
      call seed_check(error, st, 0.8_wp, 1.5_wp, 0.2_wp, -10.0_wp, -1.8_wp, 3)
      if (allocated(error)) return
      ! A thin slab on the bound's lower side, melting-point surface, no snow: category 1.
      call ice_state_seed_uniform(st, 1.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_OK, "seed 2 ok")
      if (allocated(error)) return
      call seed_check(error, st, 1.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, -1.8_wp, 1)
      if (allocated(error)) return
      ! Empty: every category the isothermal base column, aggregates empty.
      call ice_state_seed_uniform(st, 0.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_OK, "seed empty ok")
      if (allocated(error)) return
      call seed_check(error, st, 0.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, -1.8_wp, 0)
      call st%destroy()
   end subroutine test_seed

   subroutine probe_device(nx_, ny_, nc, nl, aicen, vicen, vsnon, qicen, sicen, o1, o2, o3, o4)
      integer, intent(in) :: nx_, ny_, nc, nl
      real(wp), intent(in) :: aicen(nx_, ny_, nc), vicen(nx_, ny_, nc), vsnon(nx_, ny_, nc)
      real(wp), intent(in) :: qicen(nx_, ny_, nc, nl), sicen(nx_, ny_, nc, nl)
      real(wp), intent(inout) :: o1(nx_, ny_, nc), o2(nx_, ny_, nc), o3(nx_, ny_, nc), o4(nx_, ny_, nc)
      integer :: i, j, n

      !$acc data present(aicen, vicen, vsnon, qicen, sicen, o1, o2, o3, o4)
      do concurrent(n=1:nc, j=1:ny_, i=1:nx_)
         o1(i, j, n) = aicen(i, j, n)
         o2(i, j, n) = vicen(i, j, n)
         o3(i, j, n) = vsnon(i, j, n)
         o4(i, j, n) = qicen(i, j, n, nl) + sicen(i, j, n, 1)
      end do
      !$acc end data
   end subroutine probe_device

   subroutine probe_device2(nx_, ny_, nc, Tsfcn, qsnon, aice0, aice, vice, vsno, ice_err, o1, o2, o3, o4)
      integer, intent(in) :: nx_, ny_, nc
      real(wp), intent(in) :: Tsfcn(nx_, ny_, nc), qsnon(nx_, ny_, nc, 1)
      real(wp), intent(in) :: aice0(nx_, ny_), aice(nx_, ny_), vice(nx_, ny_), vsno(nx_, ny_)
      integer(int32), intent(in) :: ice_err(nx_, ny_)
      real(wp), intent(inout) :: o1(nx_, ny_, nc), o2(nx_, ny_, nc), o3(nx_, ny_, nc), o4(nx_, ny_, nc)
      integer :: i, j, n

      !$acc data present(Tsfcn, qsnon, aice0, aice, vice, vsno, ice_err, o1, o2, o3, o4)
      do concurrent(n=1:nc, j=1:ny_, i=1:nx_)
         o1(i, j, n) = Tsfcn(i, j, n)
         o2(i, j, n) = qsnon(i, j, n, 1)
         o3(i, j, n) = aice0(i, j) + 10.0_wp*aice(i, j)
         o4(i, j, n) = vice(i, j) + 100.0_wp*vsno(i, j) + real(ice_err(i, j), wp)
      end do
      !$acc end data
   end subroutine probe_device2

   subroutine test_seed_device(error)
      !! Seeding a MAPPED slot reaches the device (host clobbered, then pulled back).
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      integer :: ierr

      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call st%enter_data()
      st%ice_err = 7_int32      ! stale device error words: the seed must clear them on the device
      !$acc update device(st%ice_err)
      call ice_state_seed_uniform(st, 0.5_wp, 2.0_wp, 0.1_wp, -20.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_OK, "seed ok")
      if (allocated(error)) return
      ! Copy what the DEVICE holds into the scratch arrays on the device, then pull
      ! those back: stale device data (a missing `update device`) shows as zeros.
      call probe_device(NX, NY, NCAT, NILYR, st%aicen, st%vicen, st%vsnon, st%qicen, st%sicen, &
                        st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init)
      !$acc update self(st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init)
      call check(error, maxval(abs(st%aicen_init - st%aicen)) == 0.0_wp, "device aicen = seeded host aicen")
      call check(error, maxval(abs(st%vicen_init - st%vicen)) == 0.0_wp, "device vicen = seeded host vicen")
      call check(error, maxval(abs(st%vsnon_init - st%vsnon)) == 0.0_wp, "device vsnon = seeded host vsnon")
      call check(error, maxval(abs(st%Tsfcn_init - (st%qicen(:, :, :, NILYR) + st%sicen(:, :, :, 1)))) &
                 <= 1.0e-6_wp*maxval(abs(st%qicen)), "device qicen/sicen = seeded host values")
      call check(error, any(st%aicen_init /= 0.0_wp), "the seeded slab is non-zero on the device")
      call probe_device2(NX, NY, NCAT, st%Tsfcn, st%qsnon, st%aice0, st%aice, st%vice, st%vsno, st%ice_err, &
                         st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init)
      !$acc update self(st%aicen_init, st%vicen_init, st%vsnon_init, st%Tsfcn_init)
      call check(error, maxval(abs(st%aicen_init - st%Tsfcn)) == 0.0_wp, "device Tsfcn = seeded host Tsfcn")
      call check(error, maxval(abs(st%vicen_init - st%qsnon(:, :, :, 1))) == 0.0_wp, "device qsnon = seeded host qsnon")
      call same(error, st%vsnon_init(3, 3, 1), 0.5_wp + 10.0_wp*0.5_wp, "device aice0 + 10 aice")
      call same(error, st%Tsfcn_init(3, 3, 1), 1.0_wp + 100.0_wp*0.05_wp + 0.0_wp, "device vice + 100 vsno + ice_err (cleared)")
      if (allocated(error)) return
      call seed_check(error, st, 0.5_wp, 2.0_wp, 0.1_wp, -20.0_wp, -1.8_wp, 3)
      call st%exit_data()
      call st%destroy()
   end subroutine test_seed_device

   subroutine test_seed_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      integer :: ierr
      real(wp) :: nan

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call ice_state_seed_uniform(st, 1.5_wp, 1.0_wp, 0.0_wp, -5.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_BAD_ARG, "aice > 1")
      call ice_state_seed_uniform(st, 0.5_wp, 0.0_wp, 0.0_wp, -5.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_BAD_ARG, "occupied with zero thickness")
      call ice_state_seed_uniform(st, 0.5_wp, 1000.0_wp, 0.0_wp, -5.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_BAD_ARG, "thickness beyond the top bound")
      call ice_state_seed_uniform(st, 0.5_wp, 1.0_wp, nan, -5.0_wp, -1.8_wp, ierr)
      call check(error, ierr, ICE_SEED_BAD_ARG, "non-finite snow depth")
      call check(error, all(st%aicen == 0.0_wp) .and. all(st%qicen == 0.0_wp), "arrays untouched on refusal")
      call st%destroy()
   end subroutine test_seed_refusals

   ! ------------------------------------------------------------------
   subroutine test_restart_registration(error)
      !! Tags, ranks, shapes, optional flags and the version scalar.
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: st
      type(restart_registry_t) :: reg
      integer :: ierr, e, k
      character(len=64) :: want

      call make_state(st, ierr)
      call check(error, ierr, 0, "init")
      if (allocated(error)) return
      call st%register_restart(reg)
 call check(error, reg%n, 1 + 4 + 2*NILYR + 1 + 3, "entry count: version + 4 + qicen/sicen per layer + qsnon + 3 structure scalars")
      call check(error, trim(reg%entries(1)%tag) == "ice_state_version" .and. reg%entries(1)%rank == 0, &
                 "version scalar first")
      call same(error, st%state_version, ICE_STATE_VERSION, "state_version = 2")
      call same(error, ICE_STATE_VERSION, 2.0_wp, "ICE_STATE_VERSION is 2")
      if (allocated(error)) return
      do e = 1, reg%n
         call check(error, reg%entries(e)%optional, "every ice entry is optional")
         if (reg%entries(e)%rank == 3) then
            call check(error, reg%entries(e)%nx_phys == NX .and. reg%entries(e)%ny_phys == NY .and. &
                       reg%entries(e)%ng == 0, "full ghosted local array")
            call check(error, reg%entries(e)%nk == NCAT, "third extent is ncat")
         end if
      end do
      call check(error, trim(reg%entries(2)%tag) == "ice_aicen" .and. trim(reg%entries(5)%tag) == "ice_Tsfcn", &
                 "state tags")
      do k = 1, NILYR
         write (want, "(a,i0)") "ice_qicen_k", k
         call check(error, trim(reg%entries(4 + 2*k)%tag) == trim(want), trim(want))
         write (want, "(a,i0)") "ice_sicen_k", k
         call check(error, trim(reg%entries(5 + 2*k)%tag) == trim(want), trim(want))
      end do
      call check(error, trim(reg%entries(reg%n - 3)%tag) == "ice_qsnon_k1" .and. trim(reg%entries(reg%n)%tag) == "ice_nslyr", &
                 "snow layer tag, then the structure scalars")
      call st%destroy()
   end subroutine test_restart_registration

   subroutine test_restart_round_trip(error)
      !! Write the registered fields with the registry's NetCDF writer into a fresh
      !! file and read them into a second slot: bitwise identical, ghosts included.
      !! (The engine-level path, `scheme = "cice"` through the model, is owed by
      !! the PR that makes it runnable.)
      type(error_type), allocatable, intent(out) :: error
      type(ice_cice_state_t) :: a, b
      type(restart_registry_t) :: ra, rb
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      character(len=*), parameter :: FN = "test_cice_state_restart.nc"
      real(wp) :: t
      integer :: ierr, step, osc, i, j, n, k, lun
      integer(int64) :: s

      call make_state(a, ierr)
      call check(error, ierr, 0, "init a")
      if (allocated(error)) return
      call make_state(b, ierr)
      call check(error, ierr, 0, "init b")
      if (allocated(error)) return
      call fill_random(a, 777_int64, 0.2_wp)
      s = 12345_int64
      do k = 1, NILYR
         do n = 1, NCAT
            do j = 1, NY
               do i = 1, NX
                  s = lcg_next(s)
                  a%qicen(i, j, n, k) = -3.0e8_wp*lcg_unit(s)
                  a%sicen(i, j, n, k) = 3.0_wp*real(i + j + n + k, wp)/17.0_wp
               end do
            end do
         end do
      end do
      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               a%Tsfcn(i, j, n) = -0.37_wp*real(i + 13*j + 101*n, wp)/17.0_wp
               a%qsnon(i, j, n, 1) = -1.1e8_wp*real(3*i + 5*j + 7*n, wp)/29.0_wp
            end do
         end do
      end do
      call a%register_restart(ra)
      call b%register_restart(rb)
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = NG
      open (newunit=lun, file=FN, status="replace")
      close (lun, status="delete")
      call ocean_restart_write_local(FN, ra, d, meta, 3.5_wp, 7, 2, ierr)
      call check(error, ierr, 0, "write ok")
      if (allocated(error)) return
      call ocean_restart_read_local(FN, rb, d, meta, t, step, osc, ierr)
      call check(error, ierr, 0, "read ok")
      if (allocated(error)) return
      call check(error, rb%entries(1)%found .and. rb%entries(2)%found, "version scalar and arrays found in the file")
      call same(error, b%state_version, 2.0_wp, "version read back")
      call check(error, all(a%aicen == b%aicen) .and. all(a%vicen == b%vicen) .and. all(a%vsnon == b%vsnon), &
                 "aicen/vicen/vsnon bitwise")
      call check(error, all(a%Tsfcn == b%Tsfcn), "Tsfcn bitwise")
      call check(error, all(a%qicen == b%qicen) .and. all(a%sicen == b%sicen), "qicen/sicen bitwise (every layer)")
      call check(error, all(a%qsnon == b%qsnon), "qsnon bitwise")
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      call a%destroy()
      call b%destroy()
   end subroutine test_restart_round_trip

   subroutine test_sis2_untouched(error)
      !! `scheme = "sis2"` (the default): the cice slot is never initialised, mapped,
      !! counted or registered by the parent state.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(restart_registry_t) :: reg
      integer :: e
      integer(int64) :: with_slot

      cfg%sim_type = "ocean"
      cfg%nx = NXP
      cfg%ny = NYP
      cfg%dx = 10.0e3_wp
      cfg%dy = 10.0e3_wp
      cfg%nz_layers = 2
      cfg%nghost = NG
      cfg%ocean%topo%max_depth = 1000.0_wp
      cfg%initial_temperature = 10.0_wp
      cfg%initial_salinity = 34.0_wp
      call grid%init(NXP, NYP, NG, 10.0e3_wp, 10.0e3_wp)
      state%multilayer%nz_ml = 2
      state%ice%enable = .true.
      state%ice%ncat = 2
      state%ice%nk_ice = 2
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call ocean_state_seed_from_cfg(state, grid, cfg)
      call register_default_tracers( &
         state%multilayer%tracers(state%multilayer%idx_salinity), &
         state%multilayer%tracers(state%multilayer%idx_temperature), cfg)
      call ocean_state_enter_data(state)

      call check(error, state%ice%is_init, "the SIS2 slot is live (the test would be vacuous otherwise)")
      call check(error,.not. state%ice%cice%is_init .and. .not. state%ice%cice%mapped, &
                 "cice slot neither initialised nor mapped by the parent enter_data")
      call check(error,.not. allocated(state%ice%cice%aicen), "cice arrays not allocated")
      call check(error, state%ice%cice%bytes() == 0_int64, "cice slot contributes 0 bytes")
      with_slot = state%ice%bytes()
      call check(error, with_slot > 0_int64, "SIS2 slot bytes are counted")

      call ocean_state_build_restart_registry(state, grid, reg)
      do e = 1, reg%n
         call check(error, index(reg%entries(e)%tag, "ice_state_version") == 0 .and. &
                    index(reg%entries(e)%tag, "ice_aicen") == 0 .and. index(reg%entries(e)%tag, "ice_vicen") == 0 &
                    .and. index(reg%entries(e)%tag, "ice_vsnon") == 0 .and. index(reg%entries(e)%tag, "ice_Tsfcn") == 0 &
                    .and. index(reg%entries(e)%tag, "ice_qicen") == 0 .and. index(reg%entries(e)%tag, "ice_sicen") == 0 &
                    .and. index(reg%entries(e)%tag, "ice_qsnon") == 0, &
                    "no v2 ice tag registered: "//trim(reg%entries(e)%tag))
         if (allocated(error)) exit
      end do
      call ocean_state_exit_data(state)
      call state%destroy()
   end subroutine test_sis2_untouched

end module test_cice_state
