!! Sea-ice transport (D10.1): tracer dependency table, pack/unpack between the
!! category state and the transport fields, and the six transport checks.
module test_cice_transport_pack
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Test structure derived from the behaviour list of CICE's transport
   ! pack/unpack unit test (CICE Consortium, https://github.com/CICE-Consortium);
   ! written from the stated behaviours, not copied.  Modified by the roundabout
   ! developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   !! **What is pinned.**  (1) the table builder for the shipped tracer set at
   !! nilyr 1, 4, 7, 8 and each of its refusals; (2) the six pack/unpack
   !! behaviours (pack zeroes `trm` where the area is not above `puny`; round
   !! trip on `aicen > puny`; `0 < aicen <= puny` keeps the area and zeroes
   !! volumes and tracers; ice removed; ice added to an ice-free cell;
   !! ice-free before and after) plus the empty block; (3) the round trip on
   !! 100 000 random category cells, with the measured worst; (4) the conserved
   !! quantities of the packed fields against the state; (5) each of the six
   !! checks triggered and not triggered, NaN raising its bit, the EFP sums
   !! against plain sums and permutation-invariant bitwise, and the
   !! monotonicity bounds against an independent two-pass oracle; (6) the
   !! device against an independent host reference on EVERY field of EVERY
   !! cell, with the ring outside the range untouched; (7) throughput.
   !!
   !! The host reference (`ref_pack`, `ref_unpack`) spells the tracer layout
   !! out (`Tsfc | qice(k) | qsno | sice(k)`) instead of calling the library
   !! accessors, so a layout disagreement is caught.
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_decompose, efp_to_real, EFP_DIGITS
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, &
                             ICE_ERR_TRANSPORT_NEG_AREA, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK
   use rdb_ice_transport_pack, only: transport_tracer_table_t, transport_tracer_table_from_depend, &
                                     transport_tracer_table_build, transport_tracer_table_message, &
                                     transport_workspace_t, transport_pack, transport_unpack, &
                                     TT_OK, TT_ERR_COUNT, TT_ERR_UNKNOWN_DEP, TT_ERR_ORDER, TT_ERR_TYPE3
   use rdb_ice_transport_checks, only: transport_check_workspace_t, transport_nquantities, &
                                       transport_check_finite, transport_check_neg_area, transport_check_cfl, &
                                       transport_conserved_sums, transport_conservation_compare, &
                                       transport_bounds, transport_check_monotone, transport_efp_decompose_seq, &
                                       transport_check_slivers
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_transport_pack_tests

   type :: st_t
      !! A host state with the shapes of `ice_cice_state_t` (explicit components so
      !! the tests do not depend on the ocean-linked state module).
      integer :: nx = 0, ny = 0, nc = 0, nil = 0, nsl = 1
      real(wp), allocatable :: aice0(:, :), aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
   end type st_t

   real(wp), parameter :: PUNY = 1.0e-11_wp

contains

   subroutine collect_cice_transport_pack_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("table_shipped_sets", test_table_shipped), &
                  new_unittest("table_builder_refusals", test_table_refusals), &
                  new_unittest("behaviour_1_pack_values", test_b1_pack), &
                  new_unittest("behaviour_2_roundtrip_above_puny", test_b2_roundtrip), &
                  new_unittest("behaviour_3_tiny_area_keeps_area_only", test_b3_tiny), &
                  new_unittest("behaviour_4_5_6_remove_add_empty", test_b456), &
                  new_unittest("behaviour_7_empty_block", test_b7_empty), &
                  new_unittest("identity_100k_random_cells", test_identity_100k), &
                  new_unittest("conserved_quantities_pack", test_conserved_pack), &
                  new_unittest("check_nonfinite_entry", test_chk_finite), &
                  new_unittest("check_negative_area", test_chk_neg), &
                  new_unittest("check_cfl_both_directions", test_chk_cfl), &
                  new_unittest("check_sliver_precondition", test_chk_sliver), &
                  new_unittest("check_workspace_shape_guards", test_chk_guards), &
                  new_unittest("efp_decompose_pinned", test_efp_pinned), &
                  new_unittest("check_conservation_area_and_tracers", test_chk_conservation), &
                  new_unittest("check_monotonicity_bounds_oracle", test_chk_monotone), &
                  new_unittest("device_vs_host_every_cell", test_device_every_cell), &
                  new_unittest("nan_in_raises_bit", test_nan_bit), &
                  new_unittest("throughput_512x512", test_throughput) &
                  ]
   end subroutine collect_cice_transport_pack_tests

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   subroutine st_alloc(st, nx, ny, nc, nil)
      type(st_t), intent(inout) :: st
      integer, intent(in) :: nx, ny, nc, nil

      st%nx = nx; st%ny = ny; st%nc = nc; st%nil = nil; st%nsl = 1
      allocate (st%aice0(nx, ny), st%aicen(nx, ny, nc), st%vicen(nx, ny, nc), st%vsnon(nx, ny, nc))
      allocate (st%Tsfcn(nx, ny, nc), st%qicen(nx, ny, nc, nil), st%sicen(nx, ny, nc, nil), st%qsnon(nx, ny, nc, 1))
      st%aice0 = 1.0_wp
      st%aicen = 0.0_wp; st%vicen = 0.0_wp; st%vsnon = 0.0_wp; st%Tsfcn = 0.0_wp
      st%qicen = 0.0_wp; st%sicen = 0.0_wp; st%qsnon = 0.0_wp
   end subroutine st_alloc

   subroutine st_copy(src, dst)
      type(st_t), intent(in) :: src
      type(st_t), intent(inout) :: dst
      dst%aice0 = src%aice0; dst%aicen = src%aicen; dst%vicen = src%vicen; dst%vsnon = src%vsnon
      dst%Tsfcn = src%Tsfcn; dst%qicen = src%qicen; dst%sicen = src%sicen; dst%qsnon = src%qsnon
   end subroutine st_copy

   function tr_get(st, it, i, j, n) result(v)
      !! Tracer `it` of the packed order, spelled out (independent of the library accessor).
      type(st_t), intent(in) :: st
      integer, intent(in) :: it, i, j, n
      real(wp) :: v
      if (it == 1) then
         v = st%Tsfcn(i, j, n)
      else if (it <= 1 + st%nil) then
         v = st%qicen(i, j, n, it - 1)
      else if (it == 2 + st%nil) then
         v = st%qsnon(i, j, n, 1)
      else
         v = st%sicen(i, j, n, it - 2 - st%nil)
      end if
   end function tr_get

   subroutine tr_set(st, it, i, j, n, v)
      type(st_t), intent(inout) :: st
      integer, intent(in) :: it, i, j, n
      real(wp), intent(in) :: v
      if (it == 1) then
         st%Tsfcn(i, j, n) = v
      else if (it <= 1 + st%nil) then
         st%qicen(i, j, n, it - 1) = v
      else if (it == 2 + st%nil) then
         st%qsnon(i, j, n, 1) = v
      else
         st%sicen(i, j, n, it - 2 - st%nil) = v
      end if
   end subroutine tr_set

   subroutine fill_random(st, seed, tiny_frac, empty_frac, aice_max)
      !! Random state: a fraction `empty_frac` of the category cells empty (zero area and volumes,
      !! junk tracers), `tiny_frac` with `0 < aicen <= puny`, the rest regular.
      type(st_t), intent(inout) :: st
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: tiny_frac, empty_frac, aice_max
      integer :: i, j, n, k
      integer(int64) :: s
      real(wp) :: u, a, h, rq
      real(wp), parameter :: SHIFT = 330.0_wp*3.337e5_wp

      s = lcg_next(seed)
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               s = lcg_next(s); u = lcg_unit(s)
               if (u < empty_frac) then
                  a = 0.0_wp
               else if (u < empty_frac + tiny_frac) then
                  s = lcg_next(s)
                  a = PUNY*(0.01_wp + 0.99_wp*lcg_unit(s))
                  if (mod(i + j, 7) == 0) a = PUNY            ! exactly the threshold
               else
                  s = lcg_next(s)
                  a = 1.0e-3_wp + (aice_max - 1.0e-3_wp)*lcg_unit(s)
               end if
               s = lcg_next(s); h = 0.05_wp + 5.0_wp*lcg_unit(s)
               st%aicen(i, j, n) = a
               st%vicen(i, j, n) = a*h
               s = lcg_next(s)
               if (lcg_unit(s) < 0.3_wp) then
                  st%vsnon(i, j, n) = 0.0_wp
               else
                  s = lcg_next(s); st%vsnon(i, j, n) = a*0.5_wp*lcg_unit(s)
               end if
               s = lcg_next(s); st%Tsfcn(i, j, n) = -30.0_wp*lcg_unit(s)
               do k = 1, st%nil
                  s = lcg_next(s); st%qicen(i, j, n, k) = -(1.0e8_wp + 2.0e8_wp*lcg_unit(s))
                  s = lcg_next(s); st%sicen(i, j, n, k) = 10.0_wp*lcg_unit(s)
               end do
               s = lcg_next(s); rq = lcg_unit(s)
               st%qsnon(i, j, n, 1) = -SHIFT*(1.0_wp - 0.1_wp*rq)
            end do
         end do
      end do
      do j = 1, st%ny
         do i = 1, st%nx
            st%aice0(i, j) = 1.0_wp - sum(st%aicen(i, j, :))
         end do
      end do
   end subroutine fill_random

   subroutine fill_transportchk(st, c)
      !! The deterministic state of the pack/unpack behaviour list: category cells of
      !! zero area, half the threshold, exactly the threshold, and regular ones.
      type(st_t), intent(inout) :: st
      type(ice_constants_t), intent(in) :: c
      integer :: i, j, n, it, ntr

      ntr = 1 + 2*st%nil + 1
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               select case (mod(i + 3*j + 5*n, 8))
               case (0, 3)
                  st%aicen(i, j, n) = 0.0_wp
               case (1)
                  st%aicen(i, j, n) = 0.5_wp*c%puny
               case (2)
                  st%aicen(i, j, n) = c%puny
               case default
                  st%aicen(i, j, n) = 0.01_wp + 0.1_wp*real(mod(i*j + n, 9), wp)
               end select
               st%vicen(i, j, n) = st%aicen(i, j, n)*(0.3_wp + 0.7_wp*real(n, wp) + 0.01_wp*real(i, wp))
               if (mod(i + j, 3) == 0) then
                  st%vsnon(i, j, n) = 0.0_wp
               else
                  st%vsnon(i, j, n) = st%aicen(i, j, n)*0.03_wp*real(mod(j, 7) + 1, wp)
               end if
               do it = 1, ntr
                  if (it == 2 + st%nil) then
                     call tr_set(st, it, i, j, n, -c%rhos*c%Lfresh - 1.0e5_wp*real(it + n, wp))
                  else
                     call tr_set(st, it, i, j, n, -1.7_wp + 0.37_wp*real(it, wp) + 0.013_wp*real(i - j, wp) &
                                 + 0.11_wp*real(n, wp))
                  end if
               end do
            end do
         end do
      end do
      do j = 1, st%ny
         do i = 1, st%nx
            st%aice0(i, j) = 1.0_wp - sum(st%aicen(i, j, :))
         end do
      end do
   end subroutine fill_transportchk

   subroutine ref_pack(st, c, is, ie, js, je, aim, trm)
      !! Independent host reference of the pack (plain loops, scalar shift).
      type(st_t), intent(in) :: st
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(inout) :: aim(:, :, 0:), trm(:, :, :, :)
      integer :: i, j, n, it, ntr
      real(wp) :: w1, shift

      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      do j = js, je
         do i = is, ie
            aim(i, j, 0) = st%aice0(i, j)
            do n = 1, st%nc
               aim(i, j, n) = st%aicen(i, j, n)
               if (st%aicen(i, j, n) > c%puny) then
                  w1 = 1.0_wp/st%aicen(i, j, n)
                  trm(i, j, 1, n) = st%vicen(i, j, n)*w1
                  trm(i, j, 2, n) = st%vsnon(i, j, n)*w1
                  do it = 1, ntr
                     if (it == 2 + st%nil) then
                        trm(i, j, 2 + it, n) = tr_get(st, it, i, j, n) + shift
                     else
                        trm(i, j, 2 + it, n) = tr_get(st, it, i, j, n)
                     end if
                  end do
               else
                  trm(i, j, :, n) = 0.0_wp
               end if
            end do
         end do
      end do
   end subroutine ref_pack

   subroutine ref_unpack(st, c, is, ie, js, je, aim, trm)
      type(st_t), intent(inout) :: st
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(in) :: aim(:, :, 0:), trm(:, :, :, :)
      integer :: i, j, n, it, ntr
      real(wp) :: shift

      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      do j = js, je
         do i = is, ie
            st%aice0(i, j) = aim(i, j, 0)
            do n = 1, st%nc
               if (aim(i, j, n) > 0.0_wp) then
                  st%aicen(i, j, n) = aim(i, j, n)
                  st%vicen(i, j, n) = aim(i, j, n)*trm(i, j, 1, n)
                  st%vsnon(i, j, n) = aim(i, j, n)*trm(i, j, 2, n)
                  do it = 1, ntr
                     if (it == 2 + st%nil) then
                        call tr_set(st, it, i, j, n, trm(i, j, 2 + it, n) - shift)
                     else
                        call tr_set(st, it, i, j, n, trm(i, j, 2 + it, n))
                     end if
                  end do
               else
                  st%aicen(i, j, n) = 0.0_wp
                  st%vicen(i, j, n) = 0.0_wp
                  st%vsnon(i, j, n) = 0.0_wp
               end if
            end do
         end do
      end do
   end subroutine ref_unpack

   subroutine map_in(ws, st)
      !! Device-map the state components and the workspace (component arrays only).
      type(transport_workspace_t), intent(inout) :: ws
      type(st_t), intent(inout) :: st
      !$acc enter data copyin(st%aice0, st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon)
      call ws%enter_data()
   end subroutine map_in

   subroutine map_out(ws, st)
      !! Pull the workspace and the state back and unmap.
      type(transport_workspace_t), intent(inout) :: ws
      type(st_t), intent(inout) :: st
      !$acc update self(ws%aim, ws%trm)
      !$acc update self(st%aice0, st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon)
      !$acc exit data delete(st%aice0, st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon)
      call ws%exit_data()
   end subroutine map_out

   subroutine run_pack(ws, c, st, is, ie, js, je, ierr)
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(in) :: c
      type(st_t), intent(inout) :: st
      integer, intent(in) :: is, ie, js, je
      integer(int32), intent(out) :: ierr
      call map_in(ws, st)
      call transport_pack(ws, c, st%nil, 1, is, ie, js, je, st%aice0, st%aicen, st%vicen, st%vsnon, &
                          st%Tsfcn, st%qicen, st%sicen, st%qsnon, ierr)
      call map_out(ws, st)
   end subroutine run_pack

   subroutine run_unpack(ws, c, st, is, ie, js, je, ierr)
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(in) :: c
      type(st_t), intent(inout) :: st
      integer, intent(in) :: is, ie, js, je
      integer(int32), intent(out) :: ierr
      call map_in(ws, st)
      call transport_unpack(ws, c, st%nil, 1, is, ie, js, je, st%aice0, st%aicen, st%vicen, st%vsnon, &
                            st%Tsfcn, st%qicen, st%sicen, st%qsnon, ierr)
      call map_out(ws, st)
   end subroutine run_unpack

   subroutine make_ws(ws, st, error)
      type(transport_workspace_t), intent(inout) :: ws
      type(st_t), intent(in) :: st
      type(error_type), allocatable, intent(out) :: error
      type(transport_tracer_table_t) :: tab
      integer :: ierr

      call transport_tracer_table_build(st%nil, 1, tab)
      call ws%init(st%nx, st%ny, st%nc, tab, ierr)
      call check(error, ierr, 0, "workspace init")
   end subroutine make_ws

   pure function relerr(a, b) result(e)
      !! Relative difference; the absolute one where the reference is exactly zero.
      real(wp), intent(in) :: a, b
      real(wp) :: e
      if (b == 0.0_wp) then
         e = abs(a)
      else
         e = abs(a - b)/abs(b)
      end if
   end function relerr

   function finite_all(a) result(ok)
      real(wp), intent(in) :: a(:)
      logical :: ok
      integer :: i
      ok = .true.
      do i = 1, size(a)
         if (.not. ieee_is_finite(a(i))) ok = .false.
      end do
   end function finite_all

   ! ------------------------------------------------------------------
   ! 1: tables
   ! ------------------------------------------------------------------

   subroutine test_table_shipped(error)
      type(error_type), allocatable, intent(out) :: error
      type(transport_tracer_table_t) :: tab
      integer, parameter :: NILS(4) = [1, 4, 7, 8]
      integer :: is, nil, nt, ntrace, et, ed, eb
      logical :: eh

      do is = 1, size(NILS)
         nil = NILS(is)
         call transport_tracer_table_build(nil, 1, tab)
         call check(error, tab%ierr, TT_OK, "shipped set builds")
         if (allocated(error)) return
         ntrace = 2 + 1 + 2*nil + 1
         call check(error, tab%ntrace, ntrace, "ntrace = 2 + ntrcr")
         if (allocated(error)) return
         call check(error, tab%nt_qsno, 2 + nil, "snow-enthalpy tracer index")
         if (allocated(error)) return
         do nt = 1, ntrace
            ! closed form of the dependency structure
            if (nt <= 3) then
               et = 1; ed = 0; eb = 0
            else if (nt <= 3 + nil) then
               et = 2; ed = 1; eb = 1             ! qice: depends on the ice thickness
            else if (nt == 4 + nil) then
               et = 2; ed = 2; eb = 2             ! qsno: depends on the snow thickness
            else
               et = 2; ed = 1; eb = 1             ! sice
            end if
            eh = nt <= 2
            call check(error, tab%tracer_type(nt), et, "tracer type")
            if (allocated(error)) return
            call check(error, tab%depend(nt), ed, "tracer depend")
            if (allocated(error)) return
            call check(error, tab%base_field(nt), eb, "base field")
            if (allocated(error)) return
            call check(error, tab%has_dependents(nt), eh, "has_dependents")
            if (allocated(error)) return
         end do
         call check(error, all(tab%tracer_type(ntrace + 1:) == 0), "entries beyond ntrace are zero")
         if (allocated(error)) return
      end do
      ! one layer above the fixed maximum: refused
      call transport_tracer_table_build(9, 1, tab)
      call check(error, tab%ierr, TT_ERR_COUNT, "count above the maximum refused")
   end subroutine test_table_shipped

   subroutine test_table_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(transport_tracer_table_t) :: tab
      type(transport_workspace_t) :: ws
      integer :: d3(3), dmax(ICE_NTRCR_MAX + 1), ierr

      ! a dependent tracer on a tracer (type 2 on a type-1 parent): valid
      d3 = [0, 3, 1]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_OK, "chain on a tracer accepted")
      if (allocated(error)) return
      call check(error, tab%depend(4), 3, "dependent names its parent by trm index")
      if (allocated(error)) return
      call check(error, tab%has_dependents(3), .true., "parent has a dependent")
      if (allocated(error)) return
      call check(error, tab%tracer_type(4), 2, "dependent is type 2")
      if (allocated(error)) return

      ! parent listed after its dependent
      d3 = [4, 0, 0]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_ORDER, "parent after dependent refused")
      if (allocated(error)) return
      call check(error, tab%bad_tracer, 1, "names the tracer")
      if (allocated(error)) return
      call check(error, all(tab%tracer_type == 0) .and. all(tab%depend == 0), "failed table is empty")
      if (allocated(error)) return
      call check(error, len_trim(transport_tracer_table_message(tab)) > 0, "refusal has a message")
      if (allocated(error)) return
      ! self dependence
      d3 = [0, 4, 0]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_ORDER, "self dependence refused")
      if (allocated(error)) return
      ! unknown dependencies
      d3 = [0, -1, 0]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_UNKNOWN_DEP, "negative dependency refused")
      if (allocated(error)) return
      d3 = [0, 0, 6]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_UNKNOWN_DEP, "dependency beyond the list refused")
      if (allocated(error)) return
      ! type 3: the parent is itself dependent
      d3 = [1, 3, 0]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_TYPE3, "type 3 refused")
      if (allocated(error)) return
      call check(error, tab%bad_tracer, 2, "names the type-3 tracer")
      if (allocated(error)) return
      ! counts
      call transport_tracer_table_from_depend(0, d3, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_COUNT, "zero tracers refused")
      if (allocated(error)) return
      dmax = 0
      call transport_tracer_table_from_depend(ICE_NTRCR_MAX + 1, dmax, 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_COUNT, "count above the fixed maximum refused")
      if (allocated(error)) return
      call transport_tracer_table_from_depend(3, dmax(1:2), 0, 0, tab)
      call check(error, tab%ierr, TT_ERR_COUNT, "dependency list shorter than the count refused")
      if (allocated(error)) return
      call transport_tracer_table_from_depend(ICE_NTRCR_MAX, dmax, 0, 0, tab)
      call check(error, tab%ierr, TT_OK, "exactly the maximum accepted")
      if (allocated(error)) return
      ! a failed table cannot seed a workspace
      d3 = [4, 0, 0]
      call transport_tracer_table_from_depend(3, d3, 0, 0, tab)
      call ws%init(8, 8, 3, tab, ierr)
      call check(error, ierr, TT_ERR_ORDER, "workspace refuses a failed table")
      if (allocated(error)) return
      call check(error, ws%is_init, .false., "nothing allocated")
   end subroutine test_table_refusals

   ! ------------------------------------------------------------------
   ! 2: the pack/unpack behaviours
   ! ------------------------------------------------------------------

   subroutine setup_chk_state(st, ws, c, error)
      type(st_t), intent(inout) :: st
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(out) :: c
      type(error_type), allocatable, intent(out) :: error

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, 24, 17, 5, 4)
      call fill_transportchk(st, c)
      call make_ws(ws, st, error)
   end subroutine setup_chk_state

   subroutine test_b1_pack(error)
      !! Behaviour 1: `aim` equals the area; `trm` holds thickness and the shifted tracers where
      !! `aicen > puny` and is exactly zero elsewhere.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, it, ntr, nreg, nzero
      real(wp) :: err_h, shift, expect

      call setup_chk_state(st, ws, c, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      call check(error, ierr, ICE_ERR_NONE, "clean state, no bit")
      if (allocated(error)) return
      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      err_h = 0.0_wp
      nreg = 0; nzero = 0
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (ws%aim(i, j, n) /= st%aicen(i, j, n)) call test_failed(error, "aim /= aicen")
               if (allocated(error)) return
               if (st%aicen(i, j, n) > c%puny) then
                  nreg = nreg + 1
                  err_h = max(err_h, abs(ws%trm(i, j, 1, n) - st%vicen(i, j, n)/st%aicen(i, j, n)) &
                              /(st%vicen(i, j, n)/st%aicen(i, j, n)))
                  do it = 1, ntr
                     expect = tr_get(st, it, i, j, n)
                     if (it == 2 + st%nil) expect = expect + shift
                     if (ws%trm(i, j, 2 + it, n) /= expect) then
                        call test_failed(error, "packed tracer is the state tracer (+ snow shift)")
                        return
                     end if
                  end do
               else
                  nzero = nzero + 1
                  if (any(ws%trm(i, j, :, n) /= 0.0_wp)) then
                     call test_failed(error, "trm is zero where aicen <= puny")
                     return
                  end if
               end if
            end do
         end do
      end do
      if (any(ws%aim(:, :, 0) /= st%aice0)) call test_failed(error, "aim(0) is aice0")
      if (allocated(error)) return
      write (error_unit, "(a,es10.3,2(a,i0))") "  [transport_pack] b1 max rel err hi = ", err_h, &
         "  cells above puny = ", nreg, "  at/below = ", nzero
      call check(error, nreg > 0 .and. nzero > 0, "both classes present")
      if (allocated(error)) return
      call check(error, err_h <= 4.0_wp*epsilon(1.0_wp), "hi = vicen*(1/aicen) within 4 eps of the quotient")
   end subroutine test_b1_pack

   subroutine test_b2_roundtrip(error)
      !! Behaviour 2: pack then unpack where `aicen > puny`: the area exactly, the rest to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, it, ntr, nreg
      real(wp) :: errmax, e, shift

      call setup_chk_state(st, ws, c, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      call st_alloc(out, st%nx, st%ny, st%nc, st%nil)
      out%aice0 = -1.0e36_wp; out%aicen = -1.0e36_wp; out%vicen = -1.0e36_wp; out%vsnon = -1.0e36_wp
      out%Tsfcn = -1.0e36_wp; out%qicen = -1.0e36_wp; out%sicen = -1.0e36_wp; out%qsnon = -1.0e36_wp
      call run_unpack(ws, c, out, 1, st%nx, 1, st%ny, ierr)
      call check(error, ierr, ICE_ERR_NONE, "no bit")
      if (allocated(error)) return
      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      errmax = 0.0_wp; nreg = 0
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > c%puny) then
                  nreg = nreg + 1
                  if (out%aicen(i, j, n) /= st%aicen(i, j, n)) call test_failed(error, "area exact")
                  if (allocated(error)) return
                  e = abs(out%vicen(i, j, n) - st%vicen(i, j, n))/st%vicen(i, j, n)
                  if (.not. ieee_is_finite(e)) call test_failed(error, "finite vicen")
                  if (allocated(error)) return
                  errmax = max(errmax, e)
                  if (st%vsnon(i, j, n) > 0.0_wp) then
                     e = abs(out%vsnon(i, j, n) - st%vsnon(i, j, n))/st%vsnon(i, j, n)
                     if (.not. ieee_is_finite(e)) call test_failed(error, "finite vsnon")
                     if (allocated(error)) return
                     errmax = max(errmax, e)
                  else if (out%vsnon(i, j, n) /= 0.0_wp) then
                     call test_failed(error, "zero snow stays zero")
                     return
                  end if
                  do it = 1, ntr
                     e = relerr(tr_get(out, it, i, j, n), tr_get(st, it, i, j, n))
                     if (.not. ieee_is_finite(e)) call test_failed(error, "finite tracer")
                     if (allocated(error)) return
                     errmax = max(errmax, e)
                  end do
               end if
            end do
         end do
      end do
      write (error_unit, "(a,es10.3,a,i0)") "  [transport_pack] b2 worst rel err = ", errmax, "  cells = ", nreg
      call check(error, nreg > 0, "regular cells present")
      if (allocated(error)) return
      call check(error, errmax <= 1.0e-14_wp, "round trip within 1e-14 relative")
   end subroutine test_b2_roundtrip

   subroutine test_b3_tiny(error)
      !! Behaviour 3: `0 < aicen <= puny` comes back with its area, zero volumes, zero tracers
      !! (the snow enthalpy minus the shift); `aicen == 0` comes back with zero volumes.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, it, ntr, ntiny, nzero
      real(wp) :: shift, expect

      call setup_chk_state(st, ws, c, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      call st_alloc(out, st%nx, st%ny, st%nc, st%nil)
      call st_copy(st, out)
      call run_unpack(ws, c, out, 1, st%nx, 1, st%ny, ierr)
      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      ntiny = 0; nzero = 0
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > 0.0_wp .and. st%aicen(i, j, n) <= c%puny) then
                  ntiny = ntiny + 1
                  if (out%aicen(i, j, n) /= st%aicen(i, j, n)) call test_failed(error, "tiny area kept")
                  if (allocated(error)) return
                  if (out%vicen(i, j, n) /= 0.0_wp .or. out%vsnon(i, j, n) /= 0.0_wp) then
                     call test_failed(error, "tiny volumes zeroed")
                     return
                  end if
                  do it = 1, ntr
                     expect = 0.0_wp
                     if (it == 2 + st%nil) expect = 0.0_wp - shift
                     if (tr_get(out, it, i, j, n) /= expect) then
                        call test_failed(error, "tiny tracers are zero (snow enthalpy minus the shift)")
                        return
                     end if
                  end do
               else if (st%aicen(i, j, n) == 0.0_wp) then
                  nzero = nzero + 1
                  if (out%aicen(i, j, n) /= 0.0_wp .or. out%vicen(i, j, n) /= 0.0_wp .or. &
                      out%vsnon(i, j, n) /= 0.0_wp) then
                     call test_failed(error, "empty stays empty")
                     return
                  end if
               end if
            end do
         end do
      end do
      call check(error, ntiny > 0 .and. nzero > 0, "tiny and empty cells present")
   end subroutine test_b3_tiny

   subroutine test_b456(error)
      !! Behaviours 4-6: ice removed (volumes zeroed, tracers untouched), ice added to an
      !! ice-free cell (area, volumes and tracers from the transport fields), ice-free before
      !! and after (nothing changes).
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, it, ntr, nrem, nadd, nfree
      real(wp) :: shift, expect
      real(wp), parameter :: AIMADD = 0.2_wp

      call setup_chk_state(st, ws, c, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > c%puny .and. mod(i + 2*j + n, 5) == 0) then
                  ws%aim(i, j, n) = 0.0_wp
               else if (st%aicen(i, j, n) == 0.0_wp .and. mod(i + j, 2) == 0) then
                  ws%aim(i, j, n) = AIMADD
                  ws%trm(i, j, 1, n) = 1.1_wp
                  ws%trm(i, j, 2, n) = 0.05_wp
                  do it = 1, ntr
                     ws%trm(i, j, 2 + it, n) = -0.5_wp*real(it, wp)
                     if (it == 2 + st%nil) ws%trm(i, j, 2 + it, n) = ws%trm(i, j, 2 + it, n) + shift
                  end do
               end if
            end do
         end do
      end do
      call st_alloc(out, st%nx, st%ny, st%nc, st%nil)
      call st_copy(st, out)
      call run_unpack(ws, c, out, 1, st%nx, 1, st%ny, ierr)
      nrem = 0; nadd = 0; nfree = 0
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > 0.0_wp .and. ws%aim(i, j, n) == 0.0_wp) then
                  nrem = nrem + 1
                  if (out%aicen(i, j, n) /= 0.0_wp .or. out%vicen(i, j, n) /= 0.0_wp .or. &
                      out%vsnon(i, j, n) /= 0.0_wp) then
                     call test_failed(error, "removed ice: state zeroed")
                     return
                  end if
                  do it = 1, ntr
                     if (tr_get(out, it, i, j, n) /= tr_get(st, it, i, j, n)) then
                        call test_failed(error, "removed ice: tracers untouched")
                        return
                     end if
                  end do
               else if (st%aicen(i, j, n) == 0.0_wp .and. ws%aim(i, j, n) > 0.0_wp) then
                  nadd = nadd + 1
                  if (out%aicen(i, j, n) /= AIMADD) call test_failed(error, "added area")
                  if (allocated(error)) return
                  if (out%vicen(i, j, n) /= AIMADD*1.1_wp) call test_failed(error, "added ice volume")
                  if (allocated(error)) return
                  if (out%vsnon(i, j, n) /= AIMADD*0.05_wp) call test_failed(error, "added snow volume")
                  if (allocated(error)) return
                  do it = 1, ntr
                     expect = ws%trm(i, j, 2 + it, n)
                     if (it == 2 + st%nil) expect = expect - shift
                     if (tr_get(out, it, i, j, n) /= expect) then
                        call test_failed(error, "added ice: tracers from the transport fields")
                        return
                     end if
                  end do
               else if (st%aicen(i, j, n) == 0.0_wp .and. ws%aim(i, j, n) == 0.0_wp) then
                  nfree = nfree + 1
                  if (out%aicen(i, j, n) /= 0.0_wp .or. out%vicen(i, j, n) /= 0.0_wp .or. &
                      out%vsnon(i, j, n) /= 0.0_wp) then
                     call test_failed(error, "ice-free stays ice-free")
                     return
                  end if
                  do it = 1, ntr
                     if (tr_get(out, it, i, j, n) /= tr_get(st, it, i, j, n)) then
                        call test_failed(error, "ice-free: tracers untouched")
                        return
                     end if
                  end do
               end if
            end do
         end do
      end do
      if (any(out%aice0 /= ws%aim(:, :, 0))) call test_failed(error, "aice0 from aim(0)")
      if (allocated(error)) return
      call check(error, nrem > 0 .and. nadd > 0 .and. nfree > 0, "all three classes present")
   end subroutine test_b456

   subroutine test_b7_empty(error)
      !! An all-open-water block: pack gives `aim(0) = 1`, zero everything else; the round trip
      !! returns open water and leaves the tracers alone.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr

      call setup_chk_state(st, ws, c, error)
      if (allocated(error)) return
      st%aicen = 0.0_wp; st%vicen = 0.0_wp; st%vsnon = 0.0_wp; st%aice0 = 1.0_wp
      ws%aim = -1.0e36_wp; ws%trm = -1.0e36_wp
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      call check(error, all(ws%aim(:, :, 0) == 1.0_wp) .and. all(ws%aim(:, :, 1:) == 0.0_wp), "aim")
      if (allocated(error)) return
      call check(error, all(ws%trm == 0.0_wp), "trm all zero")
      if (allocated(error)) return
      call st_alloc(out, st%nx, st%ny, st%nc, st%nil)
      call st_copy(st, out)
      out%aice0 = -1.0e36_wp; out%aicen = -1.0e36_wp
      call run_unpack(ws, c, out, 1, st%nx, 1, st%ny, ierr)
      call check(error, all(out%aice0 == 1.0_wp) .and. all(out%aicen == 0.0_wp) .and. &
                 all(out%vicen == 0.0_wp) .and. all(out%vsnon == 0.0_wp), "open water back")
      if (allocated(error)) return
      call check(error, all(out%qicen == st%qicen) .and. all(out%Tsfcn == st%Tsfcn), "tracers untouched")
   end subroutine test_b7_empty

   ! ------------------------------------------------------------------
   ! 3: identity and conservation on random states
   ! ------------------------------------------------------------------

   subroutine test_identity_100k(error)
      !! 200 x 100 cells x 5 categories = 100 000 category cells, nilyr = 7, with empty, tiny
      !! (`0 < a <= puny`) and regular categories.  Reports the measured worst per field class.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, it, ntr, nreg, ntiny, nempty
      real(wp) :: worst_area, worst_v, worst_t, worst_q, shift, e, spc, expect
      logical :: ok

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, 200, 100, 5, 7)
      call fill_random(st, 20261113_int64, 0.1_wp, 0.15_wp, 0.2_wp)
      call make_ws(ws, st, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      call st_alloc(out, st%nx, st%ny, st%nc, st%nil)
      call st_copy(st, out)
      call run_unpack(ws, c, out, 1, st%nx, 1, st%ny, ierr)
      call check(error, ierr, ICE_ERR_NONE, "no bit on a clean state")
      if (allocated(error)) return
      ntr = 1 + 2*st%nil + 1
      shift = c%rhos*c%Lfresh
      spc = spacing(shift)
      worst_area = 0.0_wp; worst_v = 0.0_wp; worst_t = 0.0_wp; worst_q = 0.0_wp
      nreg = 0; ntiny = 0; nempty = 0
      ok = .true.
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > c%puny) then
                  nreg = nreg + 1
                  worst_area = max(worst_area, abs(out%aicen(i, j, n) - st%aicen(i, j, n)))
                  if (st%vicen(i, j, n) > 0.0_wp) then
                     e = abs(out%vicen(i, j, n) - st%vicen(i, j, n))/st%vicen(i, j, n)
                     if (.not. ieee_is_finite(e)) ok = .false.
                     worst_v = max(worst_v, e)
                  end if
                  if (st%vsnon(i, j, n) > 0.0_wp) then
                     e = abs(out%vsnon(i, j, n) - st%vsnon(i, j, n))/st%vsnon(i, j, n)
                     if (.not. ieee_is_finite(e)) ok = .false.
                     worst_v = max(worst_v, e)
                  else if (out%vsnon(i, j, n) /= 0.0_wp) then
                     ok = .false.
                  end if
                  do it = 1, ntr
                     e = abs(tr_get(out, it, i, j, n) - tr_get(st, it, i, j, n))
                     if (.not. ieee_is_finite(e)) ok = .false.
                     if (it == 2 + st%nil) then
                        worst_q = max(worst_q, e/spc)
                     else
                        worst_t = max(worst_t, e)
                     end if
                  end do
               else if (st%aicen(i, j, n) > 0.0_wp) then
                  ntiny = ntiny + 1
                  if (out%aicen(i, j, n) /= st%aicen(i, j, n) .or. out%vicen(i, j, n) /= 0.0_wp) ok = .false.
                  do it = 1, ntr
                     expect = 0.0_wp
                     if (it == 2 + st%nil) expect = -shift
                     if (tr_get(out, it, i, j, n) /= expect) ok = .false.
                  end do
               else
                  nempty = nempty + 1
                  if (out%aicen(i, j, n) /= 0.0_wp .or. out%vicen(i, j, n) /= 0.0_wp) ok = .false.
                  do it = 1, ntr
                     if (tr_get(out, it, i, j, n) /= tr_get(st, it, i, j, n)) ok = .false.
                  end do
               end if
            end do
         end do
      end do
      write (error_unit, "(a,3(i0,1x),a,4es10.3)") "  [transport_pack] identity: regular/tiny/empty = ", nreg, ntiny, nempty, &
         " worst area, vol(rel), tracer(abs), qsno(ulp of shift) = ", worst_area, worst_v, worst_t, worst_q
      call check(error, ok, "every compared value finite and every class behaves")
      if (allocated(error)) return
      call check(error, nreg + ntiny + nempty == 100000, "100000 category cells")
      if (allocated(error)) return
      call check(error, nreg > 1000 .and. ntiny > 1000 .and. nempty > 1000, "classes populated")
      if (allocated(error)) return
      call check(error, worst_area == 0.0_wp, "area exact")
      if (allocated(error)) return
      call check(error, worst_t == 0.0_wp, "tracers other than the snow enthalpy exact")
      if (allocated(error)) return
      call check(error, worst_v <= 4.0_wp*epsilon(1.0_wp), "volumes within 4 eps")
      if (allocated(error)) return
      call check(error, worst_q <= 2.0_wp, "snow enthalpy within 2 ulp of the shift")
   end subroutine test_identity_100k

   subroutine test_conserved_pack(error)
      !! Conserved quantities of the packed fields equal those of the state: area, ice and snow
      !! volume, ice and snow energy, salt (area-weighted sums), EXCEPT the volumes, energies and
      !! salt of sliver categories (`0 < aicen <= puny`), which the pack drops by design: the test
      !! asserts that documented loss explicitly (area is kept) and that the sliver check reports
      !! exactly the dropped ice and snow volume.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n, k, nil
      real(wp), allocatable :: area(:, :)
      real(wp) :: s(6), p(6), sabs(6), lost(6), rel, worst, hi, hs, a, t, vi_d, vs_d, vi_t, vs_t, big
      integer(int32), allocatable :: ie(:, :)
      integer :: nsl, nsl_t

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, 120, 80, 5, 7)
      call fill_random(st, 777_int64, 0.1_wp, 0.2_wp, 0.2_wp)
      call make_ws(ws, st, error)
      if (allocated(error)) return
      allocate (ie(st%nx, st%ny)); ie = 0_int32
      !$acc enter data copyin(st%aicen, st%vicen, st%vsnon, ie)
      call transport_check_slivers(st%nx, st%ny, st%nc, 1, st%nx, 1, st%ny, c%puny, st%aicen, st%vicen, &
                                   st%vsnon, ie, nsl, vi_d, vs_d)
      !$acc update self(ie)
      !$acc exit data delete(st%aicen, st%vicen, st%vsnon, ie)
      call run_pack(ws, c, st, 1, st%nx, 1, st%ny, ierr)
      allocate (area(st%nx, st%ny))
      do j = 1, st%ny
         do i = 1, st%nx
            area(i, j) = 1.0e8_wp*(1.0_wp + 0.1_wp*real(mod(3*i + 5*j, 11), wp))
         end do
      end do
      nil = st%nil
      s = 0.0_wp; p = 0.0_wp; sabs = 0.0_wp; lost = 0.0_wp
      vi_t = 0.0_wp; vs_t = 0.0_wp; nsl_t = 0
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               ! state side.  Sliver categories (0 < aicen <= puny) are accounted in `lost`, not in `s`
               ! (adding 1e-11-sized terms into a 1e20 total only adds accumulation noise).
               s(1) = s(1) + area(i, j)*st%aicen(i, j, n)
               if (st%aicen(i, j, n) > 0.0_wp .and. st%aicen(i, j, n) <= c%puny) then
                  nsl_t = nsl_t + 1
                  vi_t = vi_t + st%vicen(i, j, n)
                  vs_t = vs_t + st%vsnon(i, j, n)
                  lost(2) = lost(2) + area(i, j)*st%vicen(i, j, n)
                  lost(3) = lost(3) + area(i, j)*st%vsnon(i, j, n)
               else
                  s(2) = s(2) + area(i, j)*st%vicen(i, j, n)
                  s(3) = s(3) + area(i, j)*st%vsnon(i, j, n)
                  do k = 1, nil
                     s(4) = s(4) + area(i, j)*st%vicen(i, j, n)*st%qicen(i, j, n, k)/real(nil, wp)
                     s(6) = s(6) + area(i, j)*st%vicen(i, j, n)*st%sicen(i, j, n, k)/real(nil, wp)
                  end do
                  s(5) = s(5) + area(i, j)*st%vsnon(i, j, n)*st%qsnon(i, j, n, 1)
               end if
               ! packed side: thickness times area, tracers in the packed order (qsno shifted back)
               a = ws%aim(i, j, n)
               hi = ws%trm(i, j, 1, n); hs = ws%trm(i, j, 2, n)
               p(1) = p(1) + area(i, j)*a
               p(2) = p(2) + area(i, j)*a*hi
               p(3) = p(3) + area(i, j)*a*hs
               do k = 1, nil
                  p(4) = p(4) + area(i, j)*a*hi*ws%trm(i, j, 3 + k, n)/real(nil, wp)
                  p(6) = p(6) + area(i, j)*a*hi*ws%trm(i, j, 4 + nil + k, n)/real(nil, wp)
               end do
               t = ws%trm(i, j, 4 + nil, n) - c%rhos*c%Lfresh
               p(5) = p(5) + area(i, j)*a*hs*t
               sabs(1) = sabs(1) + abs(area(i, j)*st%aicen(i, j, n))
               sabs(2) = sabs(2) + abs(area(i, j)*st%vicen(i, j, n))
               sabs(3) = sabs(3) + abs(area(i, j)*st%vsnon(i, j, n))
               sabs(4) = sabs(4) + abs(area(i, j)*st%vicen(i, j, n)*st%qicen(i, j, n, 1))
               sabs(5) = sabs(5) + abs(area(i, j)*st%vsnon(i, j, n)*st%qsnon(i, j, n, 1))
               sabs(6) = sabs(6) + abs(area(i, j)*st%vicen(i, j, n)*st%sicen(i, j, n, 1))
            end do
         end do
      end do
      worst = 0.0_wp
      do k = 1, 6
         if (.not. (ieee_is_finite(s(k)) .and. ieee_is_finite(p(k)))) then
            call test_failed(error, "non-finite total")
            return
         end if
         rel = abs(p(k) - s(k))/sabs(k)      ! s excludes the sliver categories (documented loss)
         worst = max(worst, rel)
      end do
      ! the loss is real: without the correction the volumes disagree well above round-off
      big = lost(2)/s(2)
      call check(error, nsl_t > 100, "slivers present")
      if (allocated(error)) return
      call check(error, big > 1.0e-13_wp, "the dropped sliver ice volume is a measurable fraction of the total")
      if (allocated(error)) return
      call check(error, nsl, nsl_t, "sliver check counts the sliver categories")
      if (allocated(error)) return
      call check(error, abs(vi_d - vi_t) <= 1.0e-12_wp*vi_t .and. abs(vs_d - vs_t) <= 1.0e-12_wp*vs_t, &
                 "sliver check reports the dropped ice and snow volume")
      if (allocated(error)) return
      call check(error, count(ie /= 0_int32) > 0 .and. all(iand(ie, ICE_ERR_TRANSPORT_CONSERV) == ie), &
                 "sliver cells carry the conservation bit")
      if (allocated(error)) return
      write (error_unit, "(a,es10.3)") "  [transport_pack] conserved totals worst rel (to the sum of |terms|) = ", worst
      call check(error, worst <= 1.0e-14_wp, "area, volumes, energies and salt conserved by the pack")
   end subroutine test_conserved_pack

   ! ------------------------------------------------------------------
   ! 4: the six checks
   ! ------------------------------------------------------------------

   subroutine packed_random(st, ws, c, error, seed, nx, ny, nc, nil)
      !! A packed random state (regular and empty categories, tiny excluded) in `ws`, host side.
      type(st_t), intent(inout) :: st
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(out) :: c
      type(error_type), allocatable, intent(out) :: error
      integer(int64), intent(in) :: seed
      integer, intent(in) :: nx, ny, nc, nil

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, nx, ny, nc, nil)
      call fill_random(st, seed, 0.0_wp, 0.2_wp, 0.2_wp)
      call make_ws(ws, st, error)
      if (allocated(error)) return
      call ref_pack(st, c, 1, nx, 1, ny, ws%aim, ws%trm)
   end subroutine packed_random

   subroutine test_chk_finite(error)
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32), allocatable :: ie(:, :)
      integer :: nviol, i, j
      real(wp) :: nan

      nan = ieee_value(nan, ieee_quiet_nan)
      call packed_random(st, ws, c, error, 11_int64, 20, 14, 5, 4)
      if (allocated(error)) return
      allocate (ie(st%nx, st%ny)); ie = 0_int32
      call ws%enter_data()
      !$acc enter data copyin(ie)
      call transport_check_finite(st%nx, st%ny, st%nc, ws%ntrace, 1, st%nx, 1, st%ny, ws%aim, ws%trm, ie, nviol)
      !$acc update self(ie)
      call check(error, nviol, 0, "clean: no cell flagged")
      if (allocated(error)) return
      call check(error, all(ie == 0_int32), "clean: no bit")
      if (allocated(error)) return
      !$acc exit data delete(ie)
      call ws%exit_data()
      ws%trm(5, 7, 3, 2) = nan
      ws%aim(9, 9, 0) = ieee_value(nan, ieee_positive_inf)
      call ws%enter_data()
      !$acc enter data copyin(ie)
      call transport_check_finite(st%nx, st%ny, st%nc, ws%ntrace, 1, st%nx, 1, st%ny, ws%aim, ws%trm, ie, nviol)
      !$acc update self(ie)
      !$acc exit data delete(ie)
      call ws%exit_data()
      call check(error, nviol, 2, "two cells flagged (NaN tracer, Inf open water)")
      if (allocated(error)) return
      do j = 1, st%ny
         do i = 1, st%nx
            if ((i == 5 .and. j == 7) .or. (i == 9 .and. j == 9)) then
               if (ie(i, j) /= ICE_ERR_NONFINITE) call test_failed(error, "flagged cell carries the bit")
            else
               if (ie(i, j) /= 0_int32) call test_failed(error, "no other cell flagged")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_chk_finite

   subroutine test_chk_neg(error)
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32), allocatable :: ie(:, :)
      integer :: nviol, i, j

      call packed_random(st, ws, c, error, 12_int64, 20, 14, 5, 4)
      if (allocated(error)) return
      allocate (ie(st%nx, st%ny)); ie = 0_int32
      ws%aim = abs(ws%aim)
      call ws%enter_data()
      !$acc enter data copyin(ie)
      call transport_check_neg_area(st%nx, st%ny, st%nc, 1, st%nx, 1, st%ny, c%puny, ws%aim, ie, nviol)
      !$acc update self(ie)
      !$acc exit data delete(ie)
      call ws%exit_data()
      call check(error, nviol, 0, "non-negative: nothing flagged")
      if (allocated(error)) return
      ws%aim(3, 3, 1) = -1.0e-3_wp        ! below -puny: flagged
      ws%aim(4, 4, 2) = -1.0e-20_wp       ! round-off negative: tolerated
      ws%aim(5, 5, 0) = -2.0_wp*c%puny    ! open water counts
      ws%aim(6, 6, 3) = -c%puny           ! exactly -puny: tolerated (the test is strict)
      ie = 0_int32
      call ws%enter_data()
      !$acc enter data copyin(ie)
      call transport_check_neg_area(st%nx, st%ny, st%nc, 1, st%nx, 1, st%ny, c%puny, ws%aim, ie, nviol)
      !$acc update self(ie)
      !$acc exit data delete(ie)
      call ws%exit_data()
      call check(error, nviol, 2, "two cells flagged")
      if (allocated(error)) return
      do j = 1, st%ny
         do i = 1, st%nx
            if ((i == 3 .and. j == 3) .or. (i == 5 .and. j == 5)) then
               if (ie(i, j) /= ICE_ERR_TRANSPORT_NEG_AREA) call test_failed(error, "bit on the flagged cell")
            else
               if (ie(i, j) /= 0_int32) call test_failed(error, "tolerated negatives are clean")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_chk_neg

   subroutine test_chk_cfl(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 20, NY = 12
      real(wp) :: uE(NX, NY), vN(NX, NY), dxT(NX, NY), dyT(NX, NY), dt, nan
      integer(int32) :: ie(NX, NY)
      integer :: nviol, i, j

      nan = ieee_value(nan, ieee_quiet_nan)
      dt = 1000.0_wp
      dxT = 1.0e4_wp; dyT = 1.0e4_wp
      uE = 5.0_wp; vN = -4.0_wp                 ! 5 km and 4 km displacements: inside
      ie = 0_int32
      !$acc enter data copyin(uE, vN, dxT, dyT, ie)
      call transport_check_cfl(NX, NY, 2, NX - 1, 2, NY - 1, dt, uE, vN, dxT, dyT, ie, nviol)
      !$acc update self(ie)
      !$acc exit data delete(uE, vN, dxT, dyT, ie)
      call check(error, nviol, 0, "inside the CFL limit: clean")
      if (allocated(error)) return
      uE(5, 5) = -15.0_wp           ! westward, 15 km > 10 km: the one-sided reference test misses this
      vN(6, 6) = 11.0_wp            ! northward, 11 km
      uE(7, 7) = 10.0_wp            ! exactly the cell width: allowed
      uE(8, 8) = nan                ! non-finite velocity
      dxT(12, 3) = 2000.0_wp        ! narrow neighbour of the east face of cell (11, 3)
      uE(11, 3) = 5.0_wp
      ie = 0_int32
      !$acc enter data copyin(uE, vN, dxT, dyT, ie)
      call transport_check_cfl(NX, NY, 2, NX - 1, 2, NY - 1, dt, uE, vN, dxT, dyT, ie, nviol)
      !$acc update self(ie)
      !$acc exit data delete(uE, vN, dxT, dyT, ie)
      call check(error, nviol, 5, "five cells flagged")
      if (allocated(error)) return
      do j = 1, NY
         do i = 1, NX
            if ((i == 5 .and. j == 5) .or. (i == 6 .and. j == 6) .or. (i == 11 .and. j == 3) .or. (i == 12 .and. j == 3)) then
               if (ie(i, j) /= ICE_ERR_TRANSPORT_CFL) call test_failed(error, "CFL bit")
            else if (i == 8 .and. j == 8) then
               if (ie(i, j) /= ICE_ERR_NONFINITE) call test_failed(error, "non-finite velocity: its own bit")
            else
               if (ie(i, j) /= 0_int32) call test_failed(error, "no other cell flagged")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_chk_cfl

   subroutine test_chk_sliver(error)
      !! A state with slivers (`0 < aicen <= puny`) trips the precondition check; a zapped state
      !! (slivers set to zero area and volumes) does not.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32), allocatable :: ie(:, :)
      integer :: nsl
      real(wp) :: vi, vs
      integer :: i, j, n

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, 20, 14, 4, 4)
      call fill_random(st, 5_int64, 0.1_wp, 0.2_wp, 0.2_wp)
      allocate (ie(st%nx, st%ny)); ie = 0_int32
      !$acc enter data copyin(st%aicen, st%vicen, st%vsnon, ie)
      call transport_check_slivers(st%nx, st%ny, st%nc, 1, st%nx, 1, st%ny, c%puny, st%aicen, st%vicen, &
                                   st%vsnon, ie, nsl, vi, vs)
      !$acc update self(ie)
      !$acc exit data delete(st%aicen, st%vicen, st%vsnon, ie)
      call check(error, nsl > 0 .and. vi > 0.0_wp, "slivers detected")
      if (allocated(error)) return
      call check(error, any(iand(ie, ICE_ERR_TRANSPORT_CONSERV) /= 0_int32), "bit set on sliver cells")
      if (allocated(error)) return
      do n = 1, st%nc
         do j = 1, st%ny
            do i = 1, st%nx
               if (st%aicen(i, j, n) > 0.0_wp .and. st%aicen(i, j, n) <= c%puny) then
                  st%aicen(i, j, n) = 0.0_wp; st%vicen(i, j, n) = 0.0_wp; st%vsnon(i, j, n) = 0.0_wp
               end if
            end do
         end do
      end do
      ie = 0_int32
      !$acc enter data copyin(st%aicen, st%vicen, st%vsnon, ie)
      call transport_check_slivers(st%nx, st%ny, st%nc, 1, st%nx, 1, st%ny, c%puny, st%aicen, st%vicen, &
                                   st%vsnon, ie, nsl, vi, vs)
      !$acc update self(ie)
      !$acc exit data delete(st%aicen, st%vicen, st%vsnon, ie)
      call check(error, nsl, 0, "zapped state: no slivers")
      if (allocated(error)) return
      call check(error, all(ie == 0_int32) .and. vi == 0.0_wp .and. vs == 0.0_wp, "zapped state: clean")
   end subroutine test_chk_sliver

   subroutine test_chk_guards(error)
      !! The checks refuse (ok = .false., nothing run) on an uninitialised workspace or one whose
      !! extents or quantity count do not match the call.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(transport_check_workspace_t) :: cw, cw0
      type(ice_constants_t) :: c
      type(efp_t), allocatable :: s(:)
      real(wp), allocatable :: area(:, :)
      integer(int32), allocatable :: ie(:, :)
      logical :: ok
      integer :: nviol, ntrace

      call packed_random(st, ws, c, error, 21_int64, 12, 10, 3, 4)
      if (allocated(error)) return
      ntrace = ws%ntrace
      allocate (area(12, 10), s(transport_nquantities(3, ntrace)), ie(12, 10))
      area = 1.0_wp; ie = 0_int32
      ! uninitialised
      call transport_conserved_sums(cw0, 12, 10, 3, ntrace, 3, 10, 3, 8, ws%tab%depend(1:ntrace), area, ws%aim, ws%trm, s, ok)
      call check(error, ok, .false., "uninitialised workspace refused (sums)")
      if (allocated(error)) return
      call transport_bounds(cw0, 12, 10, 3, ntrace, 3, 10, 3, 8, c%puny, ws%tab%tracer_type(1:ntrace), &
                            ws%tab%depend(1:ntrace), ws%aim, ws%trm, ok)
      call check(error, ok, .false., "uninitialised workspace refused (bounds)")
      if (allocated(error)) return
      ! initialised for a different shape / category count
      call cw%init(12, 10, 2, ntrace, .true.)
      call transport_conserved_sums(cw, 12, 10, 3, ntrace, 3, 10, 3, 8, ws%tab%depend(1:ntrace), area, ws%aim, ws%trm, s, ok)
      call check(error, ok, .false., "quantity-count mismatch refused (sums)")
      if (allocated(error)) return
      call cw%destroy()
      call cw%init(11, 10, 3, ntrace, .true.)
      call transport_conserved_sums(cw, 12, 10, 3, ntrace, 3, 10, 3, 8, ws%tab%depend(1:ntrace), area, ws%aim, ws%trm, s, ok)
      call check(error, ok, .false., "row-extent mismatch refused (sums)")
      if (allocated(error)) return
      call transport_bounds(cw, 12, 10, 3, ntrace, 3, 10, 3, 8, c%puny, ws%tab%tracer_type(1:ntrace), &
                            ws%tab%depend(1:ntrace), ws%aim, ws%trm, ok)
      call check(error, ok, .false., "bound-extent mismatch refused (bounds)")
      if (allocated(error)) return
      call transport_check_monotone(cw, 12, 10, 3, ntrace, 3, 10, 3, 8, c%puny, ws%tab%tracer_type(1:ntrace), &
                                    ws%tab%depend(1:ntrace), ws%aim, ws%trm, ie, nviol, ok)
      call check(error, ok .or. nviol /= 0, .false., "bound-extent mismatch refused (monotone)")
      if (allocated(error)) return
      call cw%destroy()
      ! bounds not requested at init
      call cw%init(12, 10, 3, ntrace, .false.)
      call transport_bounds(cw, 12, 10, 3, ntrace, 3, 10, 3, 8, c%puny, ws%tab%tracer_type(1:ntrace), &
                            ws%tab%depend(1:ntrace), ws%aim, ws%trm, ok)
      call check(error, ok, .false., "no bounds allocated: refused")
      call cw%destroy()
   end subroutine test_chk_guards

   subroutine test_efp_pinned(error)
      !! The device decompose equals `rdb_efp::efp_decompose` bin for bin over a magnitude table,
      !! and poisons exactly where the canonical one flags NaN/Inf/overflow.
      type(error_type), allocatable, intent(out) :: error
      real(real64) :: vals(14), r
      integer(int64) :: e(EFP_DIGITS), d1, d2, d3, d4, d5, d6, dp
      logical :: is_nan, is_ovf
      integer :: k

      vals = [0.0_real64, 1.0_real64, -1.0_real64, 1.0e-3_real64, -7.3e-12_real64, 123456.789_real64, 3.0e15_real64, &
              -2.5e-20_real64, 1.0e8_real64/3.0_real64, ieee_value(r, ieee_quiet_nan), 1.0e50_real64, &
              -1.0e50_real64, huge(1.0_real64), tiny(1.0_real64)]
      do k = 1, size(vals)
         r = vals(k)
         call efp_decompose(r, e, is_nan, is_ovf)
         call transport_efp_decompose_seq(r, d1, d2, d3, d4, d5, d6, dp)
         if (is_nan .or. is_ovf) then
            if (dp /= 1_int64) call test_failed(error, "poison where the canonical flags")
         else
            if (dp /= 0_int64) call test_failed(error, "no poison where the canonical does not")
            if (any([d1, d2, d3, d4, d5, d6] /= e)) call test_failed(error, "bins equal the canonical bins")
         end if
         if (allocated(error)) return
      end do
      r = ieee_value(r, ieee_positive_inf)
      call transport_efp_decompose_seq(r, d1, d2, d3, d4, d5, d6, dp)
      call check(error, dp, 1_int64, "+Inf poisons")
   end subroutine test_efp_pinned

   subroutine test_chk_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 30, NY = 22, NC = 5, NIL = 4
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(transport_check_workspace_t) :: cw
      type(ice_constants_t) :: c
      real(wp), allocatable :: area(:, :), aim2(:, :, :), trm2(:, :, :, :), ar2(:, :)
      type(efp_t), allocatable :: s0(:), s1(:), sp(:)
      logical :: ok
      integer(int32) :: ierr
      integer :: nq, q, i, j, n, nt, ntrace, is, ie, js, je, worst_q
      real(wp) :: worst, ref, plain, sc

      call packed_random(st, ws, c, error, 13_int64, NX, NY, NC, NIL)
      if (allocated(error)) return
      ntrace = ws%ntrace
      nq = transport_nquantities(NC, ntrace)
      call check(error, nq, (NC + 1) + NC*ntrace, "quantity count")
      if (allocated(error)) return
      allocate (area(NX, NY), s0(nq), s1(nq), sp(nq))
      do j = 1, NY
         do i = 1, NX
            area(i, j) = 1.0e8_wp*(1.0_wp + 0.1_wp*real(mod(3*i + 5*j, 11), wp))
         end do
      end do
      is = 3; ie = NX - 2; js = 3; je = NY - 2         ! owned range; the ring holds garbage
      ws%aim(1, :, :) = 1.0e30_wp; ws%trm(2, :, :, :) = 1.0e30_wp
      call cw%init(NX, NY, NC, ntrace, .false.)
      call cw%enter_data()
      call ws%enter_data()
      !$acc enter data copyin(area)
      call transport_conserved_sums(cw, NX, NY, NC, ntrace, is, ie, js, je, ws%tab%depend(1:ntrace), area, &
                                    ws%aim, ws%trm, s0, ok)
      !$acc exit data delete(area)
      call ws%exit_data()
      call check(error, ok, "sums computed")
      if (allocated(error)) return
      ! against plain host sums (same summands), relative to the sum of magnitudes
      worst = 0.0_wp
      do q = 1, nq
         plain = 0.0_wp; sc = 0.0_wp
         do j = js, je
            do i = is, ie
               if (q <= NC + 1) then
                  n = q - 1
                  ref = area(i, j)*ws%aim(i, j, n)
               else
                  n = (q - NC - 2)/ntrace + 1
                  nt = q - NC - 1 - (n - 1)*ntrace
                  if (ws%tab%depend(nt) == 0) then
                     ref = area(i, j)*ws%aim(i, j, n)*ws%trm(i, j, nt, n)
                  else
                     ref = area(i, j)*ws%aim(i, j, n)*(ws%trm(i, j, nt, n)*ws%trm(i, j, ws%tab%depend(nt), n))
                  end if
               end if
               plain = plain + ref
               sc = sc + abs(ref)
            end do
         end do
         if (sc > 0.0_wp) worst = max(worst, abs(efp_to_real(s0(q)) - plain)/sc)
         if (.not. ieee_is_finite(efp_to_real(s0(q)))) call test_failed(error, "finite EFP total")
         if (allocated(error)) return
      end do
      write (error_unit, "(a,es10.3)") "  [transport_checks] EFP vs plain sums worst rel (to sum |terms|) = ", worst
      call check(error, worst <= 1.0e-14_wp, "EFP totals equal the plain sums")
      if (allocated(error)) return
      ! bitwise independence of the summation order: reverse the owned range in i and j
      allocate (aim2(NX, NY, 0:NC), trm2(NX, NY, ntrace, NC), ar2(NX, NY))
      aim2 = ws%aim; trm2 = ws%trm; ar2 = area
      do j = js, je
         do i = is, ie
            aim2(i, j, :) = ws%aim(ie + is - i, je + js - j, :)
            trm2(i, j, :, :) = ws%trm(ie + is - i, je + js - j, :, :)
            ar2(i, j) = area(ie + is - i, je + js - j)
         end do
      end do
      call ws%enter_data()
      !$acc enter data copyin(ar2, aim2, trm2)
      call transport_conserved_sums(cw, NX, NY, NC, ntrace, is, ie, js, je, ws%tab%depend(1:ntrace), ar2, &
                                    aim2, trm2, sp, ok)
      !$acc exit data delete(ar2, aim2, trm2)
      call ws%exit_data()
      do q = 1, nq
         if (efp_to_real(sp(q)) /= efp_to_real(s0(q))) call test_failed(error, "totals bitwise order independent")
         if (allocated(error)) return
      end do
      ! unchanged: no violation
      call transport_conservation_compare(NC, ntrace, s0, s0, c%puny, ierr, worst, worst_q)
      call check(error, ierr, ICE_ERR_NONE, "identical totals pass")
      if (allocated(error)) return
      call check(error, worst, 0.0_wp, "no change")
      if (allocated(error)) return
      ! a clear open-water area error (check 4; open water appears in no tracer product)
      aim2 = ws%aim; aim2(10, 10, 0) = aim2(10, 10, 0) + 0.3_wp
      call sums_of(aim2, ws%trm, s1)
      call transport_conservation_compare(NC, ntrace, s0, s1, c%puny, ierr, worst, worst_q)
      call check(error, ierr, ICE_ERR_TRANSPORT_CONSERV, "area change flagged")
      if (allocated(error)) return
      call check(error, worst_q, 1, "worst quantity is the open-water area (check 4)")
      if (allocated(error)) return
      ! a clear tracer error: Tsfc of category 2 (check 5)
      trm2 = ws%trm; trm2(12, 12, 3, 2) = trm2(12, 12, 3, 2) - 3.0_wp
      call sums_of(ws%aim, trm2, s1)
      call transport_conservation_compare(NC, ntrace, s0, s1, c%puny, ierr, worst, worst_q)
      call check(error, ierr, ICE_ERR_TRANSPORT_CONSERV, "tracer change flagged")
      if (allocated(error)) return
      call check(error, worst_q, (NC + 1) + (2 - 1)*ntrace + 3, "worst quantity is Tsfc of category 2")
      if (allocated(error)) return
      ! a change below the tolerance: not flagged
      aim2 = ws%aim; aim2(10, 10, 1) = aim2(10, 10, 1)*(1.0_wp + 1.0e-9_wp)
      call sums_of(aim2, ws%trm, s1)
      call transport_conservation_compare(NC, ntrace, s0, s1, c%puny, ierr, worst, worst_q)
      call check(error, ierr, ICE_ERR_NONE, "a 1e-9 change of one cell is below the 1e-11 tolerance")
      if (allocated(error)) return
      ! a NaN anywhere poisons the total: flagged, never a plausible number
      trm2 = ws%trm; trm2(14, 9, 5, 3) = ieee_value(worst, ieee_quiet_nan)
      call sums_of(ws%aim, trm2, s1)
      call transport_conservation_compare(NC, ntrace, s0, s1, c%puny, ierr, worst, worst_q)
      call check(error, ierr, ICE_ERR_TRANSPORT_CONSERV, "NaN poisons the sums")
      if (allocated(error)) return
      call check(error, ieee_is_finite(worst), .false., "worst is NaN")
      call cw%destroy()
   contains
      subroutine sums_of(a_, t_, out_)
         real(wp), intent(in) :: a_(:, :, 0:), t_(:, :, :, :)
         type(efp_t), intent(out) :: out_(:)
         real(wp), allocatable :: aa(:, :, :), tt(:, :, :, :), ww(:, :)
         allocate (aa(NX, NY, 0:NC), tt(NX, NY, ntrace, NC), ww(NX, NY))
         aa = a_; tt = t_; ww = area
         call ws%enter_data()
         !$acc enter data copyin(aa, tt, ww)
         call transport_conserved_sums(cw, NX, NY, NC, ntrace, is, ie, js, je, ws%tab%depend(1:ntrace), ww, &
                                       aa, tt, out_, ok)
         !$acc exit data delete(aa, tt, ww)
         call ws%exit_data()
      end subroutine sums_of
   end subroutine test_chk_conservation

   subroutine test_chk_monotone(error)
      !! The quasi-local bounds equal an independent two-pass oracle (3x3 local extrema with
      !! home-cell substitution, then the 3x3 of those) on random fields with empty cells; the
      !! check passes on unchanged fields and flags an excursion, a NaN, and nothing outside the
      !! defined region.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 26, NY = 20, NC = 3, NIL = 3
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(transport_check_workspace_t) :: cw
      type(ice_constants_t) :: c
      integer(int32), allocatable :: ie(:, :)
      real(wp), allocatable :: lmin(:, :, :, :), lmax(:, :, :, :), qmin(:, :, :, :), qmax(:, :, :, :), trm2(:, :, :, :)
      integer :: i, j, n, nt, ntrace, nviol, is, iee, js, je, di, dj, par, j2
      logical :: ok, okb
      real(wp) :: x, w

      call packed_random(st, ws, c, error, 14_int64, NX, NY, NC, NIL)
      if (allocated(error)) return
      ntrace = ws%ntrace
      is = 3; iee = NX - 2; js = 3; je = NY - 2
      call cw%init(NX, NY, NC, ntrace, .true.)
      allocate (lmin(NX, NY, ntrace, NC), lmax(NX, NY, ntrace, NC), qmin(NX, NY, ntrace, NC), qmax(NX, NY, ntrace, NC))
      lmin = 0.0_wp; lmax = 0.0_wp; qmin = 0.0_wp; qmax = 0.0_wp
      ! oracle pass 1: local 3x3 over the ring-1 range
      do n = 1, NC
         do nt = 1, ntrace
            par = ws%tab%depend(nt)
            do j = js - 1, je + 1
               do i = is - 1, iee + 1
                  lmin(i, j, nt, n) = ws%trm(i, j, nt, n)
                  lmax(i, j, nt, n) = ws%trm(i, j, nt, n)
                  do dj = -1, 1
                     do di = -1, 1
                        if (ws%tab%tracer_type(nt) == 1) then
                           ok = ws%aim(i + di, j + dj, n) > c%puny
                        else
                           ok = abs(ws%trm(i + di, j + dj, par, n)) > c%puny
                        end if
                        x = ws%trm(i, j, nt, n)
                        if (ok) x = ws%trm(i + di, j + dj, nt, n)
                        lmin(i, j, nt, n) = min(lmin(i, j, nt, n), x)
                        lmax(i, j, nt, n) = max(lmax(i, j, nt, n), x)
                     end do
                  end do
               end do
            end do
            ! pass 2: extend by one cell
            do j = js, je
               do i = is, iee
                  qmin(i, j, nt, n) = minval(lmin(i - 1:i + 1, j - 1:j + 1, nt, n))
                  qmax(i, j, nt, n) = maxval(lmax(i - 1:i + 1, j - 1:j + 1, nt, n))
               end do
            end do
         end do
      end do
      allocate (ie(NX, NY)); ie = 0_int32
      call ws%enter_data()
      call cw%enter_data()
      !$acc enter data copyin(ie)
      call transport_bounds(cw, NX, NY, NC, ntrace, is, iee, js, je, c%puny, ws%tab%tracer_type(1:ntrace), &
                            ws%tab%depend(1:ntrace), ws%aim, ws%trm, okb)
      call check(error, okb, "bounds run")
      if (allocated(error)) return
      !$acc update self(cw%tmin, cw%tmax)
      do n = 1, NC
         do nt = 1, ntrace
            do j = js, je
               do i = is, iee
                  if (cw%tmin(i, j, nt, n) /= qmin(i, j, nt, n) .or. cw%tmax(i, j, nt, n) /= qmax(i, j, nt, n)) then
                     call test_failed(error, "bounds equal the two-pass oracle exactly")
                     return
                  end if
               end do
            end do
         end do
      end do
      ! unchanged field: inside its own bounds
      call transport_check_monotone(cw, NX, NY, NC, ntrace, is, iee, js, je, c%puny, ws%tab%tracer_type(1:ntrace), &
                                    ws%tab%depend(1:ntrace), ws%aim, ws%trm, ie, nviol, okb)
      !$acc update self(ie)
      call check(error, nviol, 0, "unchanged field is monotone")
      if (allocated(error)) return
      !$acc exit data delete(ie)
      ! excursions on a COPY of the packed field
      allocate (trm2(NX, NY, ntrace, NC))
      trm2 = ws%trm
      ! find cells with ice for categories 1 and 2 near the middle
      do i = 8, 12
         if (ws%aim(i, 10, 1) > c%puny) exit
      end do
      call check(error, ws%aim(i, 10, 1) > c%puny, "a cell with ice to perturb")
      if (allocated(error)) return
      w = max(1.0_wp, abs(qmax(i, 10, 4, 1)))
      trm2(i, 10, 4, 1) = qmax(i, 10, 4, 1) + 0.5_wp*w*c%puny    ! inside the round-off allowance
      ie = 0_int32
      !$acc enter data copyin(ie, trm2)
      call transport_check_monotone(cw, NX, NY, NC, ntrace, is, iee, js, je, c%puny, ws%tab%tracer_type(1:ntrace), &
                                    ws%tab%depend(1:ntrace), ws%aim, trm2, ie, nviol, okb)
      !$acc update self(ie)
      !$acc exit data delete(ie, trm2)
      call check(error, nviol, 0, "within the allowance: not flagged")
      if (allocated(error)) return
      trm2(i, 10, 4, 1) = qmax(i, 10, 4, 1) + 2.0_wp*w*c%puny + 1.0_wp  ! beyond
      do j2 = 11, 18
         if (ws%aim(i, j2, 1) > c%puny) exit
      end do
      call check(error, ws%aim(i, j2, 1) > c%puny, "a second cell with ice")
      if (allocated(error)) return
      trm2(i, j2, 3, 1) = ieee_value(x, ieee_quiet_nan)                  ! NaN where ice is present
      ie = 0_int32
      !$acc enter data copyin(ie, trm2)
      call transport_check_monotone(cw, NX, NY, NC, ntrace, is, iee, js, je, c%puny, ws%tab%tracer_type(1:ntrace), &
                                    ws%tab%depend(1:ntrace), ws%aim, trm2, ie, nviol, okb)
      !$acc update self(ie)
      !$acc exit data delete(ie, trm2)
      call check(error, ie(i, 10), ICE_ERR_TRANSPORT_CONSERV, "excursion flagged with the conservation bit")
      if (allocated(error)) return
      call check(error, nviol >= 1, "counted")
      if (allocated(error)) return
      call check(error, ie(i, j2), ICE_ERR_TRANSPORT_CONSERV, "the NaN excursion flags its own cell")
      if (allocated(error)) return
      ! an excursion where no ice is present is not checked
      ie = 0_int32
      trm2 = ws%trm
      do j = js, je
         do i = is, iee
            if (ws%aim(i, j, 2) <= c%puny) then
               trm2(i, j, 3, 2) = 1.0e6_wp
               exit
            end if
         end do
      end do
      !$acc enter data copyin(ie, trm2)
      call transport_check_monotone(cw, NX, NY, NC, ntrace, is, iee, js, je, c%puny, ws%tab%tracer_type(1:ntrace), &
                                    ws%tab%depend(1:ntrace), ws%aim, trm2, ie, nviol, okb)
      !$acc update self(ie)
      !$acc exit data delete(ie, trm2)
      call check(error, nviol, 0, "no ice, no check")
      call ws%exit_data()
      call cw%destroy()
   end subroutine test_chk_monotone

   ! ------------------------------------------------------------------
   ! 5: device vs host, NaN, throughput
   ! ------------------------------------------------------------------

   subroutine test_device_every_cell(error)
      !! Pack and unpack through the library (device on the GPU build) against the independent
      !! host reference: EVERY field of EVERY cell, ring included.  The range is the interior; the
      !! workspace and the state outside it must stay exactly as set.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 37, NY = 29, NC = 5, NIL = 7
      type(st_t) :: st, dv, hr
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      real(wp), allocatable :: aimh(:, :, :), trmh(:, :, :, :)
      integer(int32) :: ierr
      integer :: is, ie, js, je, nfld
      real(wp) :: tol, worst

#ifdef RDB_GPU_OFFLOAD
      tol = 1.0e-13_wp
#else
      tol = 0.0_wp
#endif
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, NX, NY, NC, NIL)
      call fill_random(st, 31415_int64, 0.1_wp, 0.15_wp, 0.2_wp)
      call make_ws(ws, st, error)
      if (allocated(error)) return
      is = 3; ie = NX - 2; js = 3; je = NY - 2
      ws%aim = 7.0_wp; ws%trm = -5.0_wp                        ! sentinels outside the range
      allocate (aimh(NX, NY, 0:NC), trmh(NX, NY, ws%ntrace, NC))
      aimh = ws%aim; trmh = ws%trm
      call ref_pack(st, c, is, ie, js, je, aimh, trmh)
      call run_pack(ws, c, st, is, ie, js, je, ierr)
      call check(error, ierr, ICE_ERR_NONE, "no bit")
      if (allocated(error)) return
      worst = 0.0_wp
      call cmp3(ws%aim, aimh, worst)
      call cmp4(ws%trm, trmh, worst)
      if (.not. (worst <= tol)) then
         write (error_unit, "(a,es10.3)") "  [transport_pack] pack device vs host worst = ", worst
         call test_failed(error, "packed fields differ from the host reference in some cell")
         return
      end if
      if (any(ws%aim(1:2, :, :) /= 7.0_wp) .or. any(ws%trm(:, NY - 1:NY, :, :) /= -5.0_wp)) then
         call test_failed(error, "ring outside the range untouched")
         return
      end if
      ! unpack: perturb the packed fields so the state really changes, then compare every state field
      ws%aim(is:ie, js:je, 1:NC) = ws%aim(is:ie, js:je, 1:NC)*0.97_wp
      ws%trm(is:ie, js:je, 3, 1:NC) = ws%trm(is:ie, js:je, 3, 1:NC) + 0.5_wp
      aimh = ws%aim; trmh = ws%trm
      call st_alloc(dv, NX, NY, NC, NIL)
      call st_copy(st, dv)
      call st_alloc(hr, NX, NY, NC, NIL)
      call st_copy(st, hr)
      call ref_unpack(hr, c, is, ie, js, je, aimh, trmh)
      call run_unpack(ws, c, dv, is, ie, js, je, ierr)
      call check(error, ierr, ICE_ERR_NONE, "no bit")
      if (allocated(error)) return
      worst = 0.0_wp
      nfld = 0
      call cmp2(dv%aice0, hr%aice0, worst)
      call cmp3(dv%aicen, hr%aicen, worst)
      call cmp3(dv%vicen, hr%vicen, worst)
      call cmp3(dv%vsnon, hr%vsnon, worst)
      call cmp3(dv%Tsfcn, hr%Tsfcn, worst)
      call cmp4(dv%qicen, hr%qicen, worst)
      call cmp4(dv%sicen, hr%sicen, worst)
      call cmp4(dv%qsnon, hr%qsnon, worst)
      if (.not. (worst <= tol)) then
         write (error_unit, "(a,es10.3)") "  [transport_pack] unpack device vs host worst = ", worst
         call test_failed(error, "unpacked state differs from the host reference in some cell")
         return
      end if
      ! the state outside the range is exactly the original
    if (any(dv%aicen(1:2, :, :) /= st%aicen(1:2, :, :)) .or. any(dv%qicen(:, NY - 1:NY, :, :) /= st%qicen(:, NY - 1:NY, :, :))) then
         call test_failed(error, "state outside the range untouched")
         return
      end if
      write (error_unit, "(a,es10.3)") "  [transport_pack] device vs host, every field of every cell: worst = ", worst
   contains
      subroutine upd(d, h, w)
         real(wp), intent(in) :: d, h
         real(wp), intent(inout) :: w
         real(wp) :: x
         if (.not. (ieee_is_finite(d) .and. ieee_is_finite(h))) then
            w = huge(1.0_wp)
            return
         end if
         x = abs(d - h)/max(abs(h), 1.0e-30_wp)
         if (d == h) x = 0.0_wp
         if (x > w) w = x
      end subroutine upd
      subroutine cmp2(d, h, w)
         real(wp), intent(in) :: d(:, :), h(:, :)
         real(wp), intent(inout) :: w
         integer :: a, b
         do b = 1, size(d, 2)
            do a = 1, size(d, 1)
               call upd(d(a, b), h(a, b), w)
            end do
         end do
      end subroutine cmp2
      subroutine cmp3(d, h, w)
         real(wp), intent(in) :: d(:, :, :), h(:, :, :)
         real(wp), intent(inout) :: w
         integer :: a, b, k
         do k = 1, size(d, 3)
            do b = 1, size(d, 2)
               do a = 1, size(d, 1)
                  call upd(d(a, b, k), h(a, b, k), w)
               end do
            end do
         end do
      end subroutine cmp3
      subroutine cmp4(d, h, w)
         real(wp), intent(in) :: d(:, :, :, :), h(:, :, :, :)
         real(wp), intent(inout) :: w
         integer :: a, b, k, m
         do m = 1, size(d, 4)
            do k = 1, size(d, 3)
               do b = 1, size(d, 2)
                  do a = 1, size(d, 1)
                     call upd(d(a, b, k, m), h(a, b, k, m), w)
                  end do
               end do
            end do
         end do
      end subroutine cmp4
   end subroutine test_device_every_cell

   subroutine test_nan_bit(error)
      !! A non-finite value that pack or unpack READS raises `ICE_ERR_NONFINITE`; it never
      !! disappears into the `> puny` / `> 0` branches.  A non-finite tracer of an EMPTY category
      !! is not read and is not an error.
      type(error_type), allocatable, intent(out) :: error
      type(st_t) :: st, out
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer :: i, j, n
      logical :: found
      real(wp) :: nan

      nan = ieee_value(nan, ieee_quiet_nan)
      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, 16, 12, 3, 4)
      call fill_random(st, 99_int64, 0.0_wp, 0.3_wp, 0.2_wp)
      call make_ws(ws, st, error)
      if (allocated(error)) return
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call check(error, ierr, ICE_ERR_NONE, "clean")
      if (allocated(error)) return
      ! pack: NaN volume on a cell with ice
      n = 1
      found = .false.
      do j = 2, 11
         do i = 2, 15
            if (st%aicen(i, j, n) > c%puny) then
               found = .true.
               exit
            end if
         end do
         if (found) exit
      end do
      call check(error, found, "a cell with ice")
      if (allocated(error)) return
      st%vicen(i, j, n) = nan
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0, "NaN volume raises the bit (pack)")
      if (allocated(error)) return
      st%vicen(i, j, n) = 1.0_wp
      ! NaN area: `a > puny` is false for NaN; it must still raise the bit
      st%aicen(i, j, n) = nan
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0, "NaN area raises the bit (pack)")
      if (allocated(error)) return
      st%aicen(i, j, n) = 0.1_wp
      ! NaN tracer on a cell with ice
      st%qicen(i, j, n, 2) = nan
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0, "NaN tracer raises the bit (pack)")
      if (allocated(error)) return
      st%qicen(i, j, n, 2) = -1.0e8_wp
      ! NaN tracer of an EMPTY category: not read, no bit
      st%aicen(i, j, n) = 0.0_wp; st%vicen(i, j, n) = 0.0_wp; st%vsnon(i, j, n) = 0.0_wp
      st%Tsfcn(i, j, n) = nan
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call check(error, ierr, ICE_ERR_NONE, "empty category tracers are not read")
      if (allocated(error)) return
      ! unpack: NaN area, NaN transport thickness
      call fill_random(st, 99_int64, 0.0_wp, 0.3_wp, 0.2_wp)
      call run_pack(ws, c, st, 1, 16, 1, 12, ierr)
      call st_alloc(out, 16, 12, 3, 4)
      call st_copy(st, out)
      ws%aim(i, j, n) = nan
      call run_unpack(ws, c, out, 1, 16, 1, 12, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0, "NaN aim raises the bit (unpack)")
      if (allocated(error)) return
      ws%aim(i, j, n) = st%aicen(i, j, n)
      ws%trm(i, j, 1, n) = nan
      call run_unpack(ws, c, out, 1, 16, 1, 12, ierr)
      call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0, "NaN hi raises the bit (unpack)")
   end subroutine test_nan_bit

   subroutine test_throughput(error)
      !! Pack and unpack on 512 x 512 x ncat 5 x nilyr 7: category-cells per second (GPU build:
      !! device; host build: the CPU reference number).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 512, NY = 512, NC = 5, NIL = 7, NREP = 10
      type(st_t) :: st
      type(transport_workspace_t) :: ws
      type(ice_constants_t) :: c
      integer(int32) :: ierr
      integer(int64) :: t0, t1, t2, rate
      integer :: rep, i, j, n, k
      real(wp) :: tp, tu, cells

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call st_alloc(st, NX, NY, NC, NIL)
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               st%aicen(i, j, n) = merge(0.15_wp, 0.0_wp, mod(i + 2*j + n, 5) /= 0)
               st%vicen(i, j, n) = st%aicen(i, j, n)*(0.5_wp + 0.3_wp*real(n, wp))
               st%vsnon(i, j, n) = 0.1_wp*st%aicen(i, j, n)
               st%Tsfcn(i, j, n) = -10.0_wp
               st%qsnon(i, j, n, 1) = -1.1e8_wp
               do k = 1, NIL
                  st%qicen(i, j, n, k) = -2.0e8_wp
                  st%sicen(i, j, n, k) = 5.0_wp
               end do
            end do
         end do
      end do
      call make_ws(ws, st, error)
      if (allocated(error)) return
      call map_in(ws, st)
      call transport_pack(ws, c, NIL, 1, 1, NX, 1, NY, st%aice0, st%aicen, st%vicen, st%vsnon, &
                          st%Tsfcn, st%qicen, st%sicen, st%qsnon, ierr)
      call system_clock(t0, rate)
      do rep = 1, NREP
         call transport_pack(ws, c, NIL, 1, 1, NX, 1, NY, st%aice0, st%aicen, st%vicen, st%vsnon, &
                             st%Tsfcn, st%qicen, st%sicen, st%qsnon, ierr)
      end do
      call system_clock(t1)
      do rep = 1, NREP
         call transport_unpack(ws, c, NIL, 1, 1, NX, 1, NY, st%aice0, st%aicen, st%vicen, st%vsnon, &
                               st%Tsfcn, st%qicen, st%sicen, st%qsnon, ierr)
      end do
      call system_clock(t2)
      call map_out(ws, st)
      tp = real(t1 - t0, wp)/real(rate, wp)/real(NREP, wp)
      tu = real(t2 - t1, wp)/real(rate, wp)/real(NREP, wp)
      cells = real(NX, wp)*real(NY, wp)*real(NC, wp)
      write (error_unit, "(a,2(a,es10.3),a)") "  [transport_pack] 512x512x5x7: ", "pack ", cells/tp, " cat-cells/s, unpack ", &
         cells/tu, " cat-cells/s"
      call check(error, ierr, ICE_ERR_NONE, "clean")
      if (allocated(error)) return
      call check(error, finite_all([cells/tp, cells/tu]), "timings finite")
   end subroutine test_throughput

end module test_cice_transport_pack
