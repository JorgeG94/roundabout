!! Sea-ice upwind transport (D10.2): the donor-cell kernel, the transport driver
!! (checks, pack, advection, conservation, unpack, post step) and its acceptance set.
module test_cice_upwind
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Tests written for the roundabout port of the CICE donor-cell transport (CICE
   ! Consortium, https://github.com/CICE-Consortium); the independent host oracle
   ! follows the published scheme, not the port.  Modified by the roundabout
   ! developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   !! **What is pinned.**  (1) uniform translation of a smooth periodic field
   !! conserves area, ice and snow volume, ice and snow energy, salt and the
   !! surface-temperature content: independent EFP totals of the state (order
   !! independent, so the only error is the arithmetic of the scheme), 150
   !! steps, relative change <= 1e-13, together with the driver's own
   !! conservation check and monotonicity bound; (2) a spatially constant
   !! tracer under a divergent flow with a varying area stays constant to
   !! round-off (the hierarchy creates no extremum); (3) CFL: a face just
   !! below the limit is accepted, one just above is counted (one violating
   !! cell for one violating face, either sign) and the step is refused with
   !! the state untouched; (4) mirror symmetry in x and in y; (5) first-order
   !! convergence of a smooth profile (error ratio on grid doubling ~ 2);
   !! (6) a field pushed into a wall and a land block neither leaks nor
   !! creates ice in the closed cells; (7) open water plus ice is 1 after the
   !! clean-up; (8) convergent flow exceeds 1 (documented: the ridging's job),
   !! divergent flow dilutes; (9) an empty domain and a zero velocity are
   !! no-ops (packed fields bitwise); (10) NaN velocity or field raises the
   !! bit, never produces a finite state; (11) entry slivers are refused,
   !! a negative new area is counted; (12) the library against an independent
   !! host implementation of the scheme (flux arrays, the published recovery
   !! of tracers) on EVERY field of EVERY cell, ring untouched; (13)
   !! configuration refusals; (14) throughput.
   use, intrinsic :: iso_fortran_env, only: int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_from_real, efp_plus, efp_real_diff, efp_to_real
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_NEG_AREA, &
                             ICE_ERR_TRANSPORT_CONSERV, ICE_ERR_AICE_RANGE
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_ADVECT_UPWIND, ICE_ADVECT_REMAP
   use rdb_ice_itd_params, only: itd_category_bounds
   use rdb_ice_transport_cice, only: transport_cice_ws_t, transport_result_t, ice_step_transport, &
                                     transport_conservation_compare_floor, &
                                     ICE_TRANSPORT_POST_NONE, ICE_TRANSPORT_POST_AGGREGATE, ICE_TRANSPORT_POST_CLEANUP, &
                                     TRANSPORT_CFG_OK, TRANSPORT_CFG_PARAMS, TRANSPORT_CFG_RANGE, TRANSPORT_CFG_WORKSPACE
   use rdb_ice_upwind, only: upwind_advect, UPWIND_THICK_FLOOR
   use rdb_ice_transport_pack, only: transport_pack
   use rdb_ice_transport_checks, only: transport_conservation_compare
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_upwind_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TWOPI = 6.283185307179586476925286766559_wp
   real(wp), parameter :: SHIFT = 330.0_wp*3.337e5_wp
   real(wp), parameter :: HMID(5) = [0.15_wp, 0.5_wp, 0.95_wp, 1.6_wp, 2.5_wp]
      !! Mid-bin thicknesses of the default 5-category limits (0.3, 0.7, 1.2, 2.0 m).

   type :: sim_t
      !! A host-owned C-grid ice state, its metrics and the driver workspace.  Every array is
      !! mapped once; after each `sim_step` the host copy holds the result.
      integer :: nx = 0, ny = 0, ncat = 0, nil = 0, ng = 2
      integer :: nxp = 0, nyp = 0, is = 0, ie = 0, js = 0, je = 0
      logical :: per_x = .true., per_y = .true., mapped = .false.
      real(wp) :: dx = 1.0e4_wp, dy = 1.0e4_wp
      type(ice_params_t) :: p
      type(transport_cice_ws_t) :: w
      real(wp), allocatable :: u(:, :), v(:, :), HTE(:, :), HTN(:, :), tarea(:, :), dxT(:, :), dyT(:, :)
      real(wp), allocatable :: epm(:, :), npm(:, :), hm(:, :), Tf(:, :), hin_max(:)
      real(wp), allocatable :: aice0(:, :), aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), Tsfc_c(:, :), qice_c(:, :, :), sice_c(:, :, :)
      real(wp), allocatable :: qsno_c(:, :, :), fresh(:, :), fsalt(:, :), fhocn(:, :)
      integer(int32), allocatable :: ice_err(:, :)
   end type sim_t

contains

   subroutine collect_cice_upwind_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("translation_conserves_everything", test_translation), &
                  new_unittest("translation_ncat3_nilyr4", test_translation_small), &
                  new_unittest("ice_edge_translation_tails_keep_volume", test_edge_translation), &
                  new_unittest("category_empties_in_one_step", test_one_cell_empties), &
                  new_unittest("pack_thick_floor_and_default", test_pack_floor), &
                  new_unittest("conservation_compare_floor", test_compare_floor), &
                  new_unittest("cfl_checked_with_checks_off", test_cfl_unconditional), &
                  new_unittest("constant_tracer_divergent_flow", test_constant_tracer), &
                  new_unittest("cfl_counted_and_refused", test_cfl), &
                  new_unittest("mirror_symmetry_x", test_mirror_x), &
                  new_unittest("mirror_symmetry_y", test_mirror_y), &
                  new_unittest("first_order_convergence", test_convergence), &
                  new_unittest("wall_and_land_do_not_leak", test_wall_land), &
                  new_unittest("open_water_plus_ice_is_one_after_cleanup", test_cleanup_sums_to_one), &
                  new_unittest("convergent_and_divergent_flow", test_conv_div), &
                  new_unittest("empty_domain_and_zero_velocity", test_noops), &
                  new_unittest("nan_raises_bit_state_untouched", test_nan), &
                  new_unittest("slivers_refused_negative_area_counted", test_sliver_neg), &
                  new_unittest("oracle_every_field_every_cell", test_oracle), &
                  new_unittest("oracle_doubly_periodic", test_oracle_periodic), &
                  new_unittest("configuration_refusals", test_cfg), &
                  new_unittest("throughput", test_throughput) &
                  ]
   end subroutine collect_cice_upwind_tests

   ! ------------------------------------------------------------------
   ! simulation helper
   ! ------------------------------------------------------------------

   subroutine sim_init(s, nxp, nyp, ncat, nil, per_x, per_y, check_on, want_bounds)
      !! Uniform Cartesian grid `dx = dy = 1e4 m`, a 2-cell ghost ring, all faces open (periodic) or
      !! closed at the domain edge (the ghost band is land).
      type(sim_t), intent(out) :: s
      integer, intent(in) :: nxp, nyp, ncat, nil
      logical, intent(in) :: per_x, per_y, check_on, want_bounds
      integer :: ierr, status

      s%nxp = nxp; s%nyp = nyp; s%ncat = ncat; s%nil = nil
      s%per_x = per_x; s%per_y = per_y
      s%nx = nxp + 2*s%ng; s%ny = nyp + 2*s%ng
      s%is = s%ng + 1; s%ie = s%ng + nxp; s%js = s%ng + 1; s%je = s%ng + nyp
      s%p = ice_params_default(ICE_CONSTANTS_ICEPACK, ncat, nil, 1)
      s%p%transport%advection = ICE_ADVECT_UPWIND
      s%p%transport%check = check_on
      call s%w%init(s%nx, s%ny, ncat, nil, 1, want_bounds, ierr)
      allocate (s%u(s%nx + 1, s%ny), s%v(s%nx, s%ny + 1))
      allocate (s%HTE(s%nx, s%ny), s%HTN(s%nx, s%ny), s%tarea(s%nx, s%ny), s%dxT(s%nx, s%ny), s%dyT(s%nx, s%ny))
      allocate (s%epm(s%nx, s%ny), s%npm(s%nx, s%ny), s%hm(s%nx, s%ny), s%Tf(s%nx, s%ny), s%hin_max(0:ncat))
      allocate (s%aice0(s%nx, s%ny), s%aicen(s%nx, s%ny, ncat), s%vicen(s%nx, s%ny, ncat), s%vsnon(s%nx, s%ny, ncat))
      allocate (s%Tsfcn(s%nx, s%ny, ncat), s%qicen(s%nx, s%ny, ncat, nil), s%sicen(s%nx, s%ny, ncat, nil))
      allocate (s%qsnon(s%nx, s%ny, ncat, 1))
      allocate (s%aice(s%nx, s%ny), s%vice(s%nx, s%ny), s%vsno(s%nx, s%ny), s%Tsfc_c(s%nx, s%ny))
      allocate (s%qice_c(s%nx, s%ny, nil), s%sice_c(s%nx, s%ny, nil), s%qsno_c(s%nx, s%ny, 1))
      allocate (s%fresh(s%nx, s%ny), s%fsalt(s%nx, s%ny), s%fhocn(s%nx, s%ny), s%ice_err(s%nx, s%ny))
      s%u = 0.0_wp; s%v = 0.0_wp
      s%HTE = s%dy; s%HTN = s%dx; s%tarea = s%dx*s%dy; s%dxT = s%dx; s%dyT = s%dy
      s%Tf = -1.8_wp
      s%aice0 = 1.0_wp
      s%aicen = 0.0_wp; s%vicen = 0.0_wp; s%vsnon = 0.0_wp; s%Tsfcn = -1.8_wp
      s%qicen = 0.0_wp; s%sicen = 0.0_wp; s%qsnon = 0.0_wp
      s%aice = 0.0_wp; s%vice = 0.0_wp; s%vsno = 0.0_wp; s%Tsfc_c = 0.0_wp
      s%qice_c = 0.0_wp; s%sice_c = 0.0_wp; s%qsno_c = 0.0_wp
      s%fresh = 0.0_wp; s%fsalt = 0.0_wp; s%fhocn = 0.0_wp; s%ice_err = 0_int32
      s%hin_max = 0.0_wp
      if (ncat == 5) call itd_category_bounds(s%p%itd%kcatbound, s%p%itd%kitd, ncat, s%p%thermo%hi_min, s%hin_max, status)
      call sim_masks(s)
   end subroutine sim_init

   subroutine sim_masks(s)
      !! `hm` from the periodicity (ghost band open when periodic, land otherwise) then the face
      !! masks from it (a face is open when both cells it joins are).  Land blocks are put in `hm`
      !! by the caller BETWEEN `sim_init` and this routine's second call.
      type(sim_t), intent(inout) :: s
      integer :: i, j

      do j = 1, s%ny
         do i = 1, s%nx
            s%hm(i, j) = 1.0_wp
            if (.not. s%per_x .and. (i < s%is .or. i > s%ie)) s%hm(i, j) = 0.0_wp
            if (.not. s%per_y .and. (j < s%js .or. j > s%je)) s%hm(i, j) = 0.0_wp
         end do
      end do
      call sim_face_masks(s)
   end subroutine sim_masks

   subroutine sim_face_masks(s)
      type(sim_t), intent(inout) :: s
      integer :: i, j

      s%epm = 0.0_wp
      s%npm = 0.0_wp
      do j = 1, s%ny
         do i = 1, s%nx - 1
            s%epm(i, j) = s%hm(i, j)*s%hm(i + 1, j)
         end do
      end do
      do j = 1, s%ny - 1
         do i = 1, s%nx
            s%npm(i, j) = s%hm(i, j)*s%hm(i, j + 1)
         end do
      end do
   end subroutine sim_face_masks

   subroutine sim_map(s)
      type(sim_t), intent(inout) :: s
      if (s%mapped) return
      !$acc enter data copyin(s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max)
      !$acc enter data copyin(s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon)
      !$acc enter data copyin(s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c)
      !$acc enter data copyin(s%fresh, s%fsalt, s%fhocn, s%ice_err)
      call s%w%enter_data()
      s%mapped = .true.
   end subroutine sim_map

   subroutine sim_unmap(s)
      type(sim_t), intent(inout) :: s
      if (.not. s%mapped) return
      call s%w%exit_data()
      !$acc exit data delete(s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max)
      !$acc exit data delete(s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon)
      !$acc exit data delete(s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c)
      !$acc exit data delete(s%fresh, s%fsalt, s%fhocn, s%ice_err)
      s%mapped = .false.
   end subroutine sim_unmap

   subroutine sim_push(s)
      !! Host to device: every array the host may have edited.
      type(sim_t), intent(inout) :: s
      !$acc update device(s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max)
      !$acc update device(s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon)
      !$acc update device(s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c)
      !$acc update device(s%fresh, s%fsalt, s%fhocn, s%ice_err)
   end subroutine sim_push

   subroutine sim_pull(s)
      type(sim_t), intent(inout) :: s
      !$acc update self(s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon)
      !$acc update self(s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c)
      !$acc update self(s%fresh, s%fsalt, s%fhocn, s%ice_err)
   end subroutine sim_pull

   subroutine sim_step(s, dt, post, res, wrap)
      !! Wrap the ghost ring on the host (periodic axes; the host holds the latest state), push,
      !! run the driver, pull.  `wrap = .false.` skips the wrap (the caller prepared the ring).
      type(sim_t), intent(inout) :: s
      real(wp), intent(in) :: dt
      integer, intent(in) :: post
      type(transport_result_t), intent(out) :: res
      logical, intent(in), optional :: wrap

      if (.not. s%mapped) call sim_map(s)
      if (.not. present(wrap)) then
         call sim_wrap(s)
      else if (wrap) then
         call sim_wrap(s)
      end if
      call sim_push(s)
      call ice_step_transport(s%w, s%p, dt, post, s%is, s%ie, s%js, s%je, s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, &
                              s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, s%aicen, s%vicen, s%vsnon, &
                              s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, &
                              s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res)
      call sim_pull(s)
   end subroutine sim_step

   subroutine chk_ierr(error, res, msg)
      !! `res%ierr == 0`, printing the whole result when not.
      type(error_type), allocatable, intent(inout) :: error
      type(transport_result_t), intent(in) :: res
      character(len=*), intent(in) :: msg

      if (res%ierr /= ICE_ERR_NONE .or. res%cfg /= TRANSPORT_CFG_OK) then
        write (error_unit, "(a,a,a,z0,a,i0,a,5(a,i0),a,es10.3,a,i0)") "  [upwind] ", msg, ": ierr 0x", res%ierr, " cfg ", res%cfg, &
            ";", " n_cfl ", res%n_cfl, " n_slivers ", res%n_slivers, " n_nonfinite ", res%n_nonfinite, " n_neg ", res%n_neg, &
            " n_mono ", res%n_mono, " worst_rel ", res%worst_rel, " worst_q ", res%worst_q
         call test_failed(error, msg)
      end if
   end subroutine chk_ierr

   subroutine sim_destroy(s)
      type(sim_t), intent(inout) :: s
      call sim_unmap(s)
      call s%w%destroy()
   end subroutine sim_destroy

   pure function wrapmap_c(i, is, nxp, per) result(m)
      !! Owned cell a ghost index copies (identity when the axis is not periodic).
      integer, intent(in) :: i, is, nxp
      logical, intent(in) :: per
      integer :: m
      if (per) then
         m = is + modulo(i - is, nxp)
      else
         m = i
      end if
   end function wrapmap_c

   subroutine sim_wrap(s)
      !! Periodic ghost fill on the host: every cell array, and the face arrays (face `ie+1` is
      !! face `is`).
      type(sim_t), intent(inout) :: s
      integer :: i, j, k, m, mi, mj

      do j = 1, s%ny
         mj = wrapmap_c(j, s%js, s%nyp, s%per_y)
         do i = 1, s%nx
            mi = wrapmap_c(i, s%is, s%nxp, s%per_x)
            if (mi == i .and. mj == j) cycle
            s%aice0(i, j) = s%aice0(mi, mj)
            do k = 1, s%ncat
               s%aicen(i, j, k) = s%aicen(mi, mj, k)
               s%vicen(i, j, k) = s%vicen(mi, mj, k)
               s%vsnon(i, j, k) = s%vsnon(mi, mj, k)
               s%Tsfcn(i, j, k) = s%Tsfcn(mi, mj, k)
               s%qsnon(i, j, k, 1) = s%qsnon(mi, mj, k, 1)
               do m = 1, s%nil
                  s%qicen(i, j, k, m) = s%qicen(mi, mj, k, m)
                  s%sicen(i, j, k, m) = s%sicen(mi, mj, k, m)
               end do
            end do
         end do
      end do
      if (s%per_x .or. s%per_y) then
         do j = 1, s%ny
            mj = wrapmap_c(j, s%js, s%nyp, s%per_y)
            do i = 1, s%nx + 1
               mi = wrapmap_c(i, s%is, s%nxp, s%per_x)
               s%u(i, j) = s%u(mi, mj)
            end do
         end do
         do j = 1, s%ny + 1
            mj = wrapmap_c(j, s%js, s%nyp, s%per_y)
            do i = 1, s%nx
               mi = wrapmap_c(i, s%is, s%nxp, s%per_x)
               s%v(i, j) = s%v(mi, mj)
            end do
         end do
      end if
   end subroutine sim_wrap

   subroutine fill_smooth(s, phase, snow_free_band)
      !! Smooth periodic fields with a floor well above `puny` in every category: area
      !! `0.05 + 0.03 (1 + sin cos)`, thickness, snow (zero in a band when asked), temperature,
      !! enthalpies, salinity.  `aice0 = 1 - sum`.
      type(sim_t), intent(inout) :: s
      real(wp), intent(in) :: phase
      logical, intent(in) :: snow_free_band
      integer :: i, j, n, k
      real(wp) :: x, y, c1, c2, a, fk

      do j = 1, s%ny
         y = (real(j - s%js, wp) + 0.5_wp)/real(s%nyp, wp)
         do i = 1, s%nx
            x = (real(i - s%is, wp) + 0.5_wp)/real(s%nxp, wp)
            do n = 1, s%ncat
               c1 = sin(TWOPI*x + phase + real(n, wp))
               c2 = cos(TWOPI*y - 0.7_wp*real(n, wp))
               a = 0.05_wp + 0.03_wp*(1.0_wp + c1*c2)
               s%aicen(i, j, n) = a
               if (s%ncat == 5) then
                  s%vicen(i, j, n) = a*HMID(n)*(1.0_wp + 0.05_wp*c1)     ! inside its own thickness bin: the rebin has no work
               else
                  s%vicen(i, j, n) = a*(0.4_wp + 0.35_wp*real(n, wp) + 0.2_wp*c1)
               end if
               if (snow_free_band .and. x > 0.3_wp .and. x < 0.55_wp) then
                  s%vsnon(i, j, n) = 0.0_wp
               else
                  s%vsnon(i, j, n) = a*0.06_wp*(1.0_wp + 0.5_wp*c2)
               end if
               s%Tsfcn(i, j, n) = -6.0_wp - 8.0_wp*(1.0_wp + c1)*0.5_wp - 0.3_wp*real(n, wp)
               do k = 1, s%nil
                  fk = real(k, wp)
                  s%qicen(i, j, n, k) = -(1.8e8_wp + 6.0e6_wp*fk + 1.0e7_wp*c2*c1 + 2.0e6_wp*real(n, wp))
                  s%sicen(i, j, n, k) = 4.0_wp + 0.4_wp*fk + 1.5_wp*c1
               end do
               s%qsnon(i, j, n, 1) = -SHIFT*(1.1_wp + 0.05_wp*c1 + 0.03_wp*c2)
            end do
            s%aice0(i, j) = 1.0_wp - sum(s%aicen(i, j, :))
         end do
      end do
   end subroutine fill_smooth

   subroutine set_uniform_velocity(s, u0, v0)
      type(sim_t), intent(inout) :: s
      real(wp), intent(in) :: u0, v0
      s%u = u0
      s%v = v0
   end subroutine set_uniform_velocity

   ! ------------------------------------------------------------------
   ! measurements
   ! ------------------------------------------------------------------

   subroutine totals(s, q)
      !! Independent EFP totals of the OWNED cells (area-weighted): 1 area (ice + open water),
      !! 2 ice volume, 3 snow volume, 4 ice enthalpy `sum_k vicen*qice_k`, 5 snow enthalpy,
      !! 6 salt `sum_k vicen*sice_k`, 7 `aicen*Tsfc`.  Order independent: a pure function of the
      !! cell values.
      type(sim_t), intent(in) :: s
      type(efp_t), intent(out) :: q(7)
      integer :: i, j, n, k
      real(wp) :: ar

      do n = 1, 7
         q(n) = efp_from_real(0.0_wp)
      end do
      do j = s%js, s%je
         do i = s%is, s%ie
            ar = s%tarea(i, j)
            q(1) = efp_plus(q(1), efp_from_real(ar*s%aice0(i, j)))
            do n = 1, s%ncat
               q(1) = efp_plus(q(1), efp_from_real(ar*s%aicen(i, j, n)))
               q(2) = efp_plus(q(2), efp_from_real(ar*s%vicen(i, j, n)))
               q(3) = efp_plus(q(3), efp_from_real(ar*s%vsnon(i, j, n)))
               q(5) = efp_plus(q(5), efp_from_real(ar*s%vsnon(i, j, n)*s%qsnon(i, j, n, 1)))
               q(7) = efp_plus(q(7), efp_from_real(ar*s%aicen(i, j, n)*s%Tsfcn(i, j, n)))
               do k = 1, s%nil
                  q(4) = efp_plus(q(4), efp_from_real(ar*s%vicen(i, j, n)*s%qicen(i, j, n, k)))
                  q(6) = efp_plus(q(6), efp_from_real(ar*s%vicen(i, j, n)*s%sicen(i, j, n, k)))
               end do
            end do
         end do
      end do
   end subroutine totals

   function worst_change(q0, q1, nq) result(w)
      !! Largest relative change among the seven totals (NaN-safe: a non-finite total is `huge`).
      type(efp_t), intent(in) :: q0(7), q1(7)
      integer, intent(in), optional :: nq
      real(wp) :: w
      integer :: n, nn
      real(wp) :: b, d, r

      w = 0.0_wp
      nn = 7
      if (present(nq)) nn = nq
      do n = 1, nn
         b = efp_to_real(q0(n))
         d = efp_real_diff(q1(n), q0(n))
         if (.not. (ieee_is_finite(b) .and. ieee_is_finite(d))) then
            w = huge(1.0_wp)
            return
         end if
         r = abs(d)/max(abs(b), 1.0e-300_wp)
         if (r > w) w = r
      end do
   end function worst_change

   pure function rd(a, b) result(r)
      !! Relative difference, `huge` when either value is not finite.
      real(wp), intent(in) :: a, b
      real(wp) :: r
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) then
         r = huge(1.0_wp)
      else if (a == b) then
         r = 0.0_wp
      else
         r = abs(a - b)/max(abs(a), abs(b), 1.0e-30_wp)
      end if
   end function rd

   subroutine state_diff(a, b, worst, only_owned)
      !! Largest relative difference over every state field of two sims (ring included unless
      !! `only_owned`).
      type(sim_t), intent(in) :: a, b
      real(wp), intent(out) :: worst
      logical, intent(in) :: only_owned
      integer :: i, j, n, k, i0, i1, j0, j1

      i0 = 1; i1 = a%nx; j0 = 1; j1 = a%ny
      if (only_owned) then
         i0 = a%is; i1 = a%ie; j0 = a%js; j1 = a%je
      end if
      worst = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            worst = max(worst, rd(a%aice0(i, j), b%aice0(i, j)))
            do n = 1, a%ncat
               worst = max(worst, rd(a%aicen(i, j, n), b%aicen(i, j, n)), rd(a%vicen(i, j, n), b%vicen(i, j, n)), &
                           rd(a%vsnon(i, j, n), b%vsnon(i, j, n)), rd(a%Tsfcn(i, j, n), b%Tsfcn(i, j, n)), &
                           rd(a%qsnon(i, j, n, 1), b%qsnon(i, j, n, 1)))
               do k = 1, a%nil
                  worst = max(worst, rd(a%qicen(i, j, n, k), b%qicen(i, j, n, k)), &
                              rd(a%sicen(i, j, n, k), b%sicen(i, j, n, k)))
               end do
            end do
         end do
      end do
   end subroutine state_diff

   subroutine copy_state(a, b)
      !! `b` state := `a` state (all state arrays, ring included).
      type(sim_t), intent(in) :: a
      type(sim_t), intent(inout) :: b
      b%aice0 = a%aice0; b%aicen = a%aicen; b%vicen = a%vicen; b%vsnon = a%vsnon; b%Tsfcn = a%Tsfcn
      b%qicen = a%qicen; b%sicen = a%sicen; b%qsnon = a%qsnon
   end subroutine copy_state

   function tol_const() result(t)
      !! Round-off bound of a constant tracer through 40 divergent steps: measured worst 1.4e-15 on
      !! gfortran, 1.7e-15 on the GPU build (a few ulp per step), bound 1.5e-14 on both.
      real(wp) :: t
      t = 1.5e-14_wp
   end function tol_const

   function tol_arith() result(t)
      !! Tolerance for arithmetic-only differences (FMA contraction is on for the GPU build).
      real(wp) :: t
#ifdef RDB_GPU_OFFLOAD
      t = 1.5e-14_wp            ! measured worst 1.2e-15 (FMA on)
#else
      t = 0.0_wp                ! measured: bitwise on gfortran
#endif
   end function tol_arith

   ! ------------------------------------------------------------------
   ! (1) uniform translation conserves everything
   ! ------------------------------------------------------------------

   subroutine test_translation(error)
      type(error_type), allocatable, intent(out) :: error
      call run_translation(error, 5, 7, 150)
   end subroutine test_translation

   subroutine test_translation_small(error)
      !! The same with `ncat = 3`, `nilyr = 4` (nothing is hard-wired to 5 and 7).
      type(error_type), allocatable, intent(out) :: error
      call run_translation(error, 3, 4, 60)
   end subroutine test_translation_small

   subroutine run_translation(error, ncat, nil, nstep)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ncat, nil, nstep
      type(sim_t) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      real(wp) :: worst, w, dt, lib_worst
      integer :: it

      call sim_init(s, 40, 28, ncat, nil, .true., .true., .true., .true.)
      call fill_smooth(s, 0.3_wp, .true.)
      call set_uniform_velocity(s, 0.35_wp, -0.22_wp)
      dt = 8000.0_wp
      call totals(s, q0)
      worst = 0.0_wp
      lib_worst = 0.0_wp
      do it = 1, nstep
         call sim_step(s, dt, ICE_TRANSPORT_POST_AGGREGATE, res)
         call chk_ierr(error, res, "no error bit")
         if (allocated(error)) exit
         if (.not. (res%stepped .and. res%conserv_done)) then
            call test_failed(error, "step taken and the conservation comparison ran")
            exit
         end if
         lib_worst = max(lib_worst, res%worst_rel)
         call totals(s, q1)
         w = worst_change(q0, q1)
         worst = max(worst, w)
      end do
write (error_unit, "(a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,i0,a)") "  [upwind] translation ncat ", ncat, " nilyr ", nil, ", ", nstep, &
         " steps: independent totals worst rel change ", worst, ", driver's own comparison worst ", lib_worst, &
         ", launches/step ", res%n_launch, " (check on)"
      if (.not. allocated(error)) then
         if (.not. (worst <= 1.0e-13_wp)) call test_failed(error, "conserved totals drift beyond 1e-13")
      end if
      if (.not. allocated(error)) then
         if (.not. (lib_worst <= 1.0e-13_wp)) call test_failed(error, "driver conservation worst beyond 1e-13")
      end if
      call sim_destroy(s)
   end subroutine run_translation

   subroutine test_edge_translation(error)
      !! An ice EDGE (ice in half the domain, none in the rest) translated with the checks ON: the
      !! diffusive tails reach areas far below `puny` within a few steps; their volume must be carried
      !! (no refusal, volume, snow, enthalpies, salt and area conserved), only the surface-temperature
      !! content of cells below `puny` is dropped by design.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: i, j, n, it, nsl
      real(wp) :: worst6, worst7, dt, amin

      call sim_init(s, 64, 8, 5, 7, .true., .true., .true., .true.)
      call fill_smooth(s, 0.6_wp, .false.)
      do j = 1, s%ny
         do i = 1, s%nx
            if (wrapmap_c(i, s%is, s%nxp, .true.) - s%is >= s%nxp/2) then
               s%aicen(i, j, :) = 0.0_wp
               s%vicen(i, j, :) = 0.0_wp
               s%vsnon(i, j, :) = 0.0_wp
               s%Tsfcn(i, j, :) = s%Tf(i, j)
               s%qicen(i, j, :, :) = 0.0_wp
               s%sicen(i, j, :, :) = 0.0_wp
               s%qsnon(i, j, :, :) = 0.0_wp
               s%aice0(i, j) = 1.0_wp
            end if
         end do
      end do
      call set_uniform_velocity(s, 0.35_wp, 0.0_wp)
      dt = 8000.0_wp
      call totals(s, q0)
      worst6 = 0.0_wp
      worst7 = 0.0_wp
      nsl = 0
      do it = 1, 40
         call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
         call chk_ierr(error, res, "edge translation runs without refusal")
         if (allocated(error)) exit
         call totals(s, q1)
         worst6 = max(worst6, worst_change(q0, q1, 6))
         worst7 = max(worst7, worst_change(q0, q1, 7))
      end do
      if (.not. allocated(error)) then
         amin = huge(1.0_wp)
         do n = 1, 5
            do j = s%js, s%je
               do i = s%is, s%ie
                  if (s%aicen(i, j, n) > 0.0_wp .and. s%aicen(i, j, n) <= PUNY) nsl = nsl + 1
                  amin = min(amin, s%aicen(i, j, n))
               end do
            end do
         end do
         write (error_unit, "(a,es10.3,a,es10.3,a,i0,a)") "  [upwind] ice edge, 40 steps, check on: area/volumes/enthalpies/salt worst change ", &
            worst6, "; Tsfc content (sub-puny cells dropped by design) ", worst7, "; sliver categories at the end ", nsl, ""
         if (.not. (worst6 <= 1.0e-13_wp)) call test_failed(error, "volumes lost through the tails")
      end if
      if (.not. allocated(error)) then
         if (nsl == 0) call test_failed(error, "the test must produce slivers to mean anything")
      end if
      call sim_destroy(s)
   end subroutine test_edge_translation

   subroutine test_one_cell_empties(error)
      !! A category that empties in one step: a single occupied cell shifted by exactly one cell
      !! (`|u| dt = dx`, the largest displacement the CFL check accepts); the source cell comes back
      !! empty (Tsfc = Tf), the destination carries the tracers, volume conserved.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: k
      real(wp) :: worst

      call sim_init(s, 16, 12, 3, 4, .true., .true., .true., .true.)
      s%aicen(8, 6, 2) = 0.3_wp
      s%vicen(8, 6, 2) = 0.3_wp
      s%vsnon(8, 6, 2) = 0.03_wp
      s%Tsfcn(8, 6, 2) = -9.5_wp
      do k = 1, 4
         s%qicen(8, 6, 2, k) = -2.5e8_wp - 1.0e7_wp*real(k, wp)
         s%sicen(8, 6, 2, k) = 4.0_wp + 0.3_wp*real(k, wp)
      end do
      s%qsnon(8, 6, 2, 1) = -SHIFT*1.2_wp
      s%aice0(8, 6) = 0.7_wp
      call set_uniform_velocity(s, 1.0_wp, 0.0_wp)
      call totals(s, q0)
      call sim_step(s, 1.0e4_wp, ICE_TRANSPORT_POST_NONE, res)
      call chk_ierr(error, res, "one-cell shift")
      if (allocated(error)) return
      call totals(s, q1)
      worst = max(rd(s%aicen(9, 6, 2), 0.3_wp), rd(s%vicen(9, 6, 2), 0.3_wp), rd(s%vsnon(9, 6, 2), 0.03_wp), &
                  rd(s%Tsfcn(9, 6, 2), -9.5_wp), rd(s%qsnon(9, 6, 2, 1), -SHIFT*1.2_wp))
      do k = 1, 4
 worst = max(worst, rd(s%qicen(9, 6, 2, k), -2.5e8_wp - 1.0e7_wp*real(k, wp)), rd(s%sicen(9, 6, 2, k), 4.0_wp + 0.3_wp*real(k, wp)))
      end do
write (error_unit, "(a,es10.3,a,es10.3,a,es10.3)") "  [upwind] category emptied in one step: destination vs source worst ", worst, &
         ", source area ", s%aicen(8, 6, 2), ", totals change ", worst_change(q0, q1, 6)
      if (.not. (worst <= 1.0e-14_wp)) call test_failed(error, "the destination must carry the source cell")
      if (.not. allocated(error)) then
         if (.not. (s%aicen(8, 6, 2) == 0.0_wp .and. s%vicen(8, 6, 2) == 0.0_wp .and. s%Tsfcn(8, 6, 2) == s%Tf(8, 6))) &
            call test_failed(error, "the source category is empty with Tsfc = Tf")
      end if
      if (.not. allocated(error)) then
         if (.not. (worst_change(q0, q1, 6) <= 1.0e-13_wp)) call test_failed(error, "totals conserved")
      end if
      call sim_destroy(s)
   end subroutine test_one_cell_empties

   ! ------------------------------------------------------------------
   ! (2) constant tracer under a divergent flow
   ! ------------------------------------------------------------------

   subroutine test_constant_tracer(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(transport_result_t) :: res
      integer :: i, j, n, k, it
      real(wp) :: x, worst, dt, amax, amin
      real(wp), parameter :: C_T = -7.3_wp, C_Q = -2.2e8_wp, C_S = 5.1_wp, C_QS = -1.1e8_wp

      call sim_init(s, 36, 24, 5, 7, .true., .true., .true., .true.)
      call fill_smooth(s, 1.1_wp, .false.)
      s%Tsfcn = C_T
      s%qicen = C_Q
      s%sicen = C_S
      s%qsnon = C_QS
      do j = 1, s%ny
         do i = 1, s%nx + 1
            x = (real(i - s%is, wp))/real(s%nxp, wp)
            s%u(i, j) = 0.4_wp*sin(TWOPI*x)
         end do
      end do
      do j = 1, s%ny + 1
         do i = 1, s%nx
            x = (real(j - s%js, wp))/real(s%nyp, wp)
            s%v(i, j) = 0.25_wp*sin(TWOPI*x + 0.5_wp)
         end do
      end do
      dt = 6000.0_wp
      worst = 0.0_wp
      amax = 0.0_wp
      amin = huge(1.0_wp)
      do it = 1, 40
         call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
         call chk_ierr(error, res, "no error bit")
         if (allocated(error)) exit
         do n = 1, 5
            do j = s%js, s%je
               do i = s%is, s%ie
                  amax = max(amax, s%aicen(i, j, n))
                  amin = min(amin, s%aicen(i, j, n))
                  worst = max(worst, rd(s%Tsfcn(i, j, n), C_T), rd(s%qsnon(i, j, n, 1), C_QS))
                  do k = 1, 7
                     worst = max(worst, rd(s%qicen(i, j, n, k), C_Q), rd(s%sicen(i, j, n, k), C_S))
                  end do
               end do
            end do
         end do
      end do
     write (error_unit, "(a,es10.3,a,2es10.3)") "  [upwind] constant tracer, divergent flow, 40 steps: worst relative deviation ", &
         worst, "; area range ", amin, amax
      if (.not. allocated(error)) then
         if (.not. (amin > 1.0e-6_wp)) then
            call test_failed(error, "the area field must stay well above puny for this test to mean anything")
         else if (.not. (worst <= 1.0e-13_wp)) then
            call test_failed(error, "a constant tracer drifted beyond 1e-13")
         end if
      end if
      call sim_destroy(s)
   end subroutine test_constant_tracer

   ! ------------------------------------------------------------------
   ! (3) CFL
   ! ------------------------------------------------------------------

   subroutine test_cfl(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: dt, lim, worst
      integer :: sgn

      call sim_init(s, 20, 16, 5, 7, .true., .true., .true., .false.)
      call sim_init(ref, 20, 16, 5, 7, .true., .true., .true., .false.)
      call fill_smooth(s, 0.2_wp, .false.)
      call copy_state(s, ref)
      dt = 1000.0_wp
      lim = s%dx/dt                              ! face speed with a displacement of exactly one cell
      ! all velocities a hair below the limit: accepted (the sum condition for positivity is the
      ! scheme's, the check is the per-face displacement)
      s%u = 0.0_wp
      s%v = 0.0_wp
      s%u(10, 8) = 0.999_wp*lim
      call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
      call check(error, res%n_cfl, 0, "just below the limit: no violation")
      if (allocated(error)) return
      call chk_ierr(error, res, "just below the limit: no bit")
      if (allocated(error)) return
      ! one face just above, both signs, u and v: exactly one violating cell each
      do sgn = -1, 1, 2
         call copy_state(ref, s)
         s%u = 0.0_wp
         s%v = 0.0_wp
         s%u(10, 8) = real(sgn, wp)*1.001_wp*lim
         call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
         call check(error, res%n_cfl, 1, "one violating u face counted once")
         if (allocated(error)) return
         call check(error, iand(res%ierr, ICE_ERR_TRANSPORT_CFL) /= 0, "CFL bit raised")
         if (allocated(error)) return
         call check(error, res%stepped, .false., "the step is refused, not clamped")
         if (allocated(error)) return
         call state_diff(s, ref, worst, .true.)
         if (.not. (worst == 0.0_wp)) then
            call test_failed(error, "state untouched after a refused step")
            return
         end if
         call check(error, iand(s%ice_err(9, 8), ICE_ERR_TRANSPORT_CFL) /= 0 .or. &
                    iand(s%ice_err(10, 8), ICE_ERR_TRANSPORT_CFL) /= 0, "per-cell bit on the cell of the face")
         if (allocated(error)) return
         s%ice_err = 0_int32
      end do
      ! v face
      call copy_state(ref, s)
      s%u = 0.0_wp
      s%v = 0.0_wp
      s%v(7, 9) = -1.001_wp*s%dy/dt
      call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
      call check(error, res%n_cfl, 1, "one violating v face counted once")
      if (allocated(error)) return
      write (error_unit, "(a)") "  [upwind] CFL: 0.999 accepted, 1.001 counted (u and v, both signs), step refused, state untouched"
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_cfl

   ! ------------------------------------------------------------------
   ! (4) mirror symmetry
   ! ------------------------------------------------------------------

   subroutine mirror_state(a, b, axis)
      !! `b` owned state := the mirror image of `a` along `axis` (1 x, 2 y); velocities mirrored with
      !! their sign flip on the mirrored axis.
      type(sim_t), intent(in) :: a
      type(sim_t), intent(inout) :: b
      integer, intent(in) :: axis
      integer :: i, j, n, k, mi, mj, f, g

      do j = a%js, a%je
         mj = j
         if (axis == 2) mj = a%js + a%je - j
         do i = a%is, a%ie
            mi = i
            if (axis == 1) mi = a%is + a%ie - i
            b%aice0(i, j) = a%aice0(mi, mj)
            do n = 1, a%ncat
               b%aicen(i, j, n) = a%aicen(mi, mj, n)
               b%vicen(i, j, n) = a%vicen(mi, mj, n)
               b%vsnon(i, j, n) = a%vsnon(mi, mj, n)
               b%Tsfcn(i, j, n) = a%Tsfcn(mi, mj, n)
               b%qsnon(i, j, n, 1) = a%qsnon(mi, mj, n, 1)
               do k = 1, a%nil
                  b%qicen(i, j, n, k) = a%qicen(mi, mj, n, k)
                  b%sicen(i, j, n, k) = a%sicen(mi, mj, n, k)
               end do
            end do
         end do
      end do
      do j = a%js, a%je
         mj = j
         if (axis == 2) mj = a%js + a%je - j
         do f = a%is, a%ie
            if (axis == 1) then
               b%u(f, j) = -a%u(a%is + a%ie + 1 - f, mj)
            else
               b%u(f, j) = a%u(f, mj)
            end if
         end do
      end do
      do g = a%js, a%je
         do i = a%is, a%ie
            mi = i
            if (axis == 1) mi = a%is + a%ie - i
            if (axis == 2) then
               b%v(i, g) = -a%v(mi, a%js + a%je + 1 - g)
            else
               b%v(i, g) = a%v(mi, g)
            end if
         end do
      end do
   end subroutine mirror_state

   subroutine run_mirror(error, axis)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: axis
      type(sim_t) :: a, b
      type(transport_result_t) :: res
      integer :: i, j, it, mi, mj, n, k
      real(wp) :: worst, dt
      integer(int64) :: st

      call sim_init(a, 24, 20, 5, 7, .true., .true., .false., .false.)
      call sim_init(b, 24, 20, 5, 7, .true., .true., .false., .false.)
      call fill_smooth(a, 0.9_wp, .true.)
      ! asymmetric random perturbation of the fields so no accidental symmetry exists
      st = lcg_next(777_int64)
      do n = 1, 5
         do j = a%js, a%je
            do i = a%is, a%ie
               st = lcg_next(st)
               a%aicen(i, j, n) = a%aicen(i, j, n)*(1.0_wp + 0.3_wp*lcg_unit(st))
               st = lcg_next(st)
               a%vicen(i, j, n) = a%aicen(i, j, n)*(0.5_wp + 2.0_wp*lcg_unit(st))
               st = lcg_next(st)
               a%Tsfcn(i, j, n) = -2.0_wp - 20.0_wp*lcg_unit(st)
               do k = 1, 7
                  st = lcg_next(st)
                  a%qicen(i, j, n, k) = -1.0e8_wp - 2.0e8_wp*lcg_unit(st)
               end do
            end do
         end do
      end do
      do j = a%js, a%je
         do i = a%is, a%ie
            a%aice0(i, j) = 1.0_wp - sum(a%aicen(i, j, :))
         end do
      end do
      do j = a%js, a%je
         do i = a%is, a%ie + 1
            st = lcg_next(st)
            a%u(i, j) = 0.5_wp*(lcg_unit(st) - 0.5_wp)
         end do
      end do
      do j = a%js, a%je + 1
         do i = a%is, a%ie
            st = lcg_next(st)
            a%v(i, j) = 0.5_wp*(lcg_unit(st) - 0.5_wp)
         end do
      end do
      call sim_wrap(a)
      call mirror_state(a, b, axis)
      call sim_wrap(b)
      dt = 6000.0_wp
      do it = 1, 5
         call sim_step(a, dt, ICE_TRANSPORT_POST_NONE, res)
         call chk_ierr(error, res, "A: no error bit")
         if (allocated(error)) return
         call sim_step(b, dt, ICE_TRANSPORT_POST_NONE, res)
         call chk_ierr(error, res, "B: no error bit")
         if (allocated(error)) return
      end do
      ! compare B with the mirror of A on the owned cells
      worst = 0.0_wp
      do j = a%js, a%je
         mj = j
         if (axis == 2) mj = a%js + a%je - j
         do i = a%is, a%ie
            mi = i
            if (axis == 1) mi = a%is + a%ie - i
            worst = max(worst, rd(b%aice0(i, j), a%aice0(mi, mj)))
            do n = 1, 5
               worst = max(worst, rd(b%aicen(i, j, n), a%aicen(mi, mj, n)), rd(b%vicen(i, j, n), a%vicen(mi, mj, n)), &
                           rd(b%vsnon(i, j, n), a%vsnon(mi, mj, n)), rd(b%Tsfcn(i, j, n), a%Tsfcn(mi, mj, n)), &
                           rd(b%qsnon(i, j, n, 1), a%qsnon(mi, mj, n, 1)))
               do k = 1, 7
                  worst = max(worst, rd(b%qicen(i, j, n, k), a%qicen(mi, mj, n, k)), &
                              rd(b%sicen(i, j, n, k), a%sicen(mi, mj, n, k)))
               end do
            end do
         end do
      end do
      write (error_unit, "(a,i0,a,es10.3)") "  [upwind] mirror symmetry axis ", axis, ", 5 steps, worst relative difference ", worst
      if (.not. (worst <= tol_arith())) call test_failed(error, "the mirrored run is not the mirror of the run")
      call sim_destroy(a)
      call sim_destroy(b)
   end subroutine run_mirror

   subroutine test_mirror_x(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mirror(error, 1)
   end subroutine test_mirror_x

   subroutine test_mirror_y(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mirror(error, 2)
   end subroutine test_mirror_y

   ! ------------------------------------------------------------------
   ! (5) first-order convergence
   ! ------------------------------------------------------------------

   subroutine conv_run(n, err_a, err_t)
      !! One category, translation by a quarter of the domain at fixed CFL 0.4: L1 errors of the
      !! area and of the surface temperature against the exact translated profile.
      integer, intent(in) :: n
      real(wp), intent(out) :: err_a, err_t
      type(sim_t) :: s
      type(transport_result_t) :: res
      integer :: i, j, it, nstep
      real(wp) :: dt, u0, xl, x, ea, et, tsum, xs

      call sim_init(s, n, 4, 1, 1, .true., .true., .false., .false.)
      xl = real(n, wp)*s%dx
      u0 = 0.5_wp
      dt = 0.4_wp*s%dx/u0
      nstep = (5*n)/8
      do j = 1, s%ny
         do i = 1, s%nx
            x = (real(i - s%is, wp) + 0.5_wp)*s%dx
            s%aicen(i, j, 1) = 0.3_wp + 0.2_wp*sin(TWOPI*x/xl)
            s%vicen(i, j, 1) = s%aicen(i, j, 1)
            s%Tsfcn(i, j, 1) = -5.0_wp + 3.0_wp*cos(TWOPI*x/xl)
            s%aice0(i, j) = 1.0_wp - s%aicen(i, j, 1)
         end do
      end do
      call set_uniform_velocity(s, u0, 0.0_wp)
      do it = 1, nstep
         call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res)
      end do
      xs = u0*dt*real(nstep, wp)
      ea = 0.0_wp
      et = 0.0_wp
      tsum = 0.0_wp
      j = s%js
      do i = s%is, s%ie
         x = (real(i - s%is, wp) + 0.5_wp)*s%dx - xs
         ea = ea + abs(s%aicen(i, j, 1) - (0.3_wp + 0.2_wp*sin(TWOPI*x/xl)))
         et = et + abs(s%Tsfcn(i, j, 1) - (-5.0_wp + 3.0_wp*cos(TWOPI*x/xl)))
         tsum = tsum + 1.0_wp
      end do
      err_a = ea/tsum
      err_t = et/tsum
      call sim_destroy(s)
   end subroutine conv_run

   subroutine test_convergence(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS(4) = [32, 64, 128, 256]
      real(wp) :: ea(4), et(4), ra, rt
      integer :: k

      do k = 1, 4
         call conv_run(NS(k), ea(k), et(k))
      end do
      do k = 2, 4
         ra = ea(k - 1)/ea(k)
         rt = et(k - 1)/et(k)
         write (error_unit, "(a,i0,a,f6.3,a,f6.3,a,2es10.3)") "  [upwind] convergence N ", NS(k - 1), "->", ra, " (area), ", rt, &
            " (Tsfc); errors ", ea(k), et(k)
         if (.not. (ra > 1.8_wp .and. ra < 2.2_wp .and. rt > 1.8_wp .and. rt < 2.2_wp)) then
            call test_failed(error, "error ratio on grid doubling is not ~2 (first order)")
            return
         end if
      end do
   end subroutine test_convergence

   ! ------------------------------------------------------------------
   ! (6) wall and land
   ! ------------------------------------------------------------------

   subroutine test_wall_land(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: i, j, n, it, nleak
      real(wp) :: worst, w, dt

      call sim_init(s, 24, 16, 5, 7, .false., .false., .true., .true.)
      call fill_smooth(s, 0.0_wp, .false.)
      ! land block, empty and zero-tracer; the ghost band is land too
      do j = s%js + 5, s%js + 9
         do i = s%is + 10, s%is + 12
            s%hm(i, j) = 0.0_wp
         end do
      end do
      call sim_face_masks(s)
      do j = 1, s%ny
         do i = 1, s%nx
            if (s%hm(i, j) <= 0.5_wp) then
               s%aice0(i, j) = 0.0_wp
               s%aicen(i, j, :) = 0.0_wp
               s%vicen(i, j, :) = 0.0_wp
               s%vsnon(i, j, :) = 0.0_wp
               s%Tsfcn(i, j, :) = 0.0_wp
               s%qicen(i, j, :, :) = 0.0_wp
               s%sicen(i, j, :, :) = 0.0_wp
               s%qsnon(i, j, :, :) = 0.0_wp
            end if
         end do
      end do
      call set_uniform_velocity(s, 0.4_wp, 0.15_wp)    ! pushes into the east wall, the north wall and the block
      dt = 8000.0_wp
      call totals(s, q0)
      worst = 0.0_wp
      nleak = 0
      do it = 1, 30
         call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
         call chk_ierr(error, res, "no error bit")
         if (allocated(error)) return
         call totals(s, q1)
         w = worst_change(q0, q1)
         worst = max(worst, w)
         do j = 1, s%ny
            do i = 1, s%nx
               if (s%hm(i, j) <= 0.5_wp) then
                  if (s%aice0(i, j) /= 0.0_wp) nleak = nleak + 1
                  do n = 1, 5
                     if (s%aicen(i, j, n) /= 0.0_wp .or. s%vicen(i, j, n) /= 0.0_wp .or. s%vsnon(i, j, n) /= 0.0_wp &
                         .or. s%Tsfcn(i, j, n) /= 0.0_wp) nleak = nleak + 1
                  end do
               end if
            end do
         end do
      end do
      write (error_unit, "(a,es10.3,a,i0)") "  [upwind] wall + land block, 30 steps pushing into them: worst total change ", &
         worst, ", closed cells touched: ", nleak
      if (.not. (worst <= 1.0e-13_wp)) call test_failed(error, "ice leaked through a closed face")
      if (.not. allocated(error)) call check(error, nleak, 0, "closed cells stay exactly empty (and Tsfc 0 on land)")
      call sim_destroy(s)
   end subroutine test_wall_land

   ! ------------------------------------------------------------------
   ! (7) open water + ice = 1 after the clean-up
   ! ------------------------------------------------------------------

   subroutine test_cleanup_sums_to_one(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(transport_result_t) :: res
      integer :: i, j, it
      real(wp) :: worst

      call sim_init(s, 30, 22, 5, 7, .true., .true., .true., .false.)
      call fill_smooth(s, 0.4_wp, .true.)
      call set_uniform_velocity(s, 0.3_wp, 0.2_wp)
      worst = 0.0_wp
      do it = 1, 20
         call sim_step(s, 8000.0_wp, ICE_TRANSPORT_POST_CLEANUP, res)
         call chk_ierr(error, res, "no transport bit")
         if (allocated(error)) return
         do j = s%js, s%je
            do i = s%is, s%ie
               worst = max(worst, abs(s%aice(i, j) + s%aice0(i, j) - 1.0_wp))
               if (abs(s%aice(i, j) - sum(s%aicen(i, j, :))) > 1.0e-15_wp) worst = max(worst, 1.0_wp)
            end do
         end do
         if (any(s%ice_err(s%is:s%ie, s%js:s%je) /= 0)) then
            call test_failed(error, "the clean-up raised a bit on a smooth field")
            return
         end if
      end do
      write (error_unit, "(a,es10.3)") "  [upwind] after transport + clean-up: |aice + aice0 - 1| worst ", worst
      if (.not. (worst <= 1.0e-15_wp)) call test_failed(error, "ice + open water is not 1")
      call sim_destroy(s)
   end subroutine test_cleanup_sums_to_one

   ! ------------------------------------------------------------------
   ! (8) convergent and divergent flow
   ! ------------------------------------------------------------------

   subroutine test_conv_div(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, c
      type(transport_result_t) :: res
      integer :: i, j, n, it, ic, id
      real(wp) :: x, amax, amin, a0
      logical :: any_range

      call sim_init(s, 40, 6, 5, 7, .true., .true., .false., .false.)
      call sim_init(c, 40, 6, 5, 7, .true., .true., .false., .false.)
      do n = 1, 5
         do j = 1, s%ny
            do i = 1, s%nx
               s%aicen(i, j, n) = 0.19_wp
               s%vicen(i, j, n) = 0.19_wp*(0.5_wp + 0.4_wp*real(n, wp))
               s%vsnon(i, j, n) = 0.19_wp*0.05_wp
               s%Tsfcn(i, j, n) = -8.0_wp
               s%qicen(i, j, n, :) = -2.0e8_wp
               s%sicen(i, j, n, :) = 5.0_wp
               s%qsnon(i, j, n, 1) = -SHIFT*1.1_wp
            end do
         end do
      end do
      s%aice0 = 0.05_wp
      do j = 1, s%ny
         do i = 1, s%nx + 1
            x = real(i - s%is, wp)/real(s%nxp, wp)
            s%u(i, j) = 0.5_wp*sin(TWOPI*x)           ! + then -: converges at x = 1/2, diverges at x = 0
         end do
      end do
      call copy_state(s, c)
      c%u = s%u
      ic = s%is + s%nxp/2        ! convergence cell (the face between ic-1 and ic is at x = 1/2)
      id = s%is                  ! divergence cell
      a0 = sum(s%aicen(ic, s%js, :))
      any_range = .false.
      do it = 1, 50
         call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_AGGREGATE, res)
         call chk_ierr(error, res, "aggregate: no bit (the advection does not limit)")
         if (allocated(error)) return
         call sim_step(c, 6000.0_wp, ICE_TRANSPORT_POST_CLEANUP, res)
         if (any(c%ice_err(c%is:c%ie, c%js:c%je) /= 0)) any_range = .true.
         if (any_range) exit
      end do
      amax = maxval(s%aice(s%is:s%ie, s%js:s%je))
      amin = minval(s%aice(s%is:s%ie, s%js:s%je))
      write (error_unit, "(a,f7.4,a,f7.4,a,f7.4)") "  [upwind] converging/diverging flow: initial aice ", a0, ", max ", amax, &
         " (aggregate, unlimited), min ", amin
     if (.not. (amax > 1.0_wp + 1.0e-6_wp)) call test_failed(error, "convergent flow should exceed 1 (no limiter in the advection)")
      if (.not. allocated(error)) then
         if (.not. (amin < a0 - 1.0e-3_wp)) call test_failed(error, "divergent flow should dilute")
      end if
      if (.not. allocated(error)) then
         call check(error, any_range, "the clean-up reports the over-full cells (ICE_ERR_AICE_RANGE)")
      end if
      if (.not. allocated(error)) then
         call check(error, iand(c%ice_err(ic, c%js), ICE_ERR_AICE_RANGE) /= 0 .or. &
                    iand(c%ice_err(ic - 1, c%js), ICE_ERR_AICE_RANGE) /= 0, "bit sits at the convergence")
      end if
      if (.not. allocated(error)) then
         call check(error, iand(c%ice_err(id, c%js), ICE_ERR_AICE_RANGE), 0, "no bit at the divergence")
      end if
      call sim_destroy(s)
      call sim_destroy(c)
   end subroutine test_conv_div

   ! ------------------------------------------------------------------
   ! (9) no-ops
   ! ------------------------------------------------------------------

   subroutine test_noops(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst, wp_pack
      integer :: i, j, n, nt

      ! (a) empty domain, nonzero velocity: every state array bitwise unchanged
      call sim_init(s, 16, 12, 5, 7, .true., .true., .true., .false.)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .true., .false.)
      s%qicen = -3.0e8_wp                    ! junk in the empty categories must not be touched
      s%sicen = 7.0_wp
      call copy_state(s, ref)
      call set_uniform_velocity(s, 0.4_wp, -0.3_wp)
      call sim_step(s, 8000.0_wp, ICE_TRANSPORT_POST_NONE, res)
      call chk_ierr(error, res, "empty: no bit")
      if (allocated(error)) return
      call check(error, res%stepped, "empty: stepped")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      write (error_unit, "(a,es10.3)") "  [upwind] empty domain: worst difference ", worst
      if (.not. (worst == 0.0_wp)) then
         call test_failed(error, "empty domain must be a bitwise no-op")
         return
      end if
      call sim_destroy(s)
      call sim_destroy(ref)
      ! (b) zero velocity on a full smooth state: packed fields bitwise, state to round-off
      call sim_init(s, 20, 16, 5, 7, .true., .true., .false., .false.)
      call sim_init(ref, 20, 16, 5, 7, .true., .true., .false., .false.)
      call fill_smooth(s, 0.5_wp, .true.)
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, 8000.0_wp, ICE_TRANSPORT_POST_NONE, res)
      call check(error, res%stepped, "zero velocity: stepped")
      if (allocated(error)) return
      !$acc update self(s%w%a%aim, s%w%a%trm, s%w%b%aim, s%w%b%trm)
      wp_pack = 0.0_wp
      do j = s%js, s%je
         do i = s%is, s%ie
            wp_pack = max(wp_pack, rd(s%w%b%aim(i, j, 0), s%w%a%aim(i, j, 0)))
            do n = 1, 5
               wp_pack = max(wp_pack, rd(s%w%b%aim(i, j, n), s%w%a%aim(i, j, n)))
               do nt = 1, s%w%ntrace
                  wp_pack = max(wp_pack, rd(s%w%b%trm(i, j, nt, n), s%w%a%trm(i, j, nt, n)))
               end do
            end do
         end do
      end do
      call state_diff(s, ref, worst, .false.)
      write (error_unit, "(a,es10.3,a,es10.3)") "  [upwind] zero velocity: packed fields worst difference ", wp_pack, &
         ", state worst difference (pack round trip) ", worst
      if (.not. (wp_pack == 0.0_wp)) call test_failed(error, "zero velocity must leave the packed fields bitwise unchanged")
      if (.not. allocated(error)) then
         if (.not. (worst <= 1.0e-15_wp)) call test_failed(error, "state after the pack round trip beyond 1e-15")
      end if
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_noops

   ! ------------------------------------------------------------------
   ! (10) NaN
   ! ------------------------------------------------------------------

   subroutine test_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst, nan
      integer :: pass, i
      logical :: chk

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      do pass = 1, 4
         chk = (pass >= 3)
         call sim_init(s, 16, 12, 5, 7, .true., .true., chk, .false.)
         call sim_init(ref, 16, 12, 5, 7, .true., .true., chk, .false.)
         call fill_smooth(s, 0.1_wp, .false.)
         call set_uniform_velocity(s, 0.3_wp, 0.1_wp)
         if (mod(pass, 2) == 1) then
            s%u(9, 7) = nan                              ! NaN face velocity
         else
            s%vicen(8, 6, 2) = nan                       ! NaN field value in an owned cell
         end if
         call sim_wrap(s)
         call copy_state(s, ref)
         call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_AGGREGATE, res, wrap=.false.)
         if (iand(res%ierr, ICE_ERR_NONFINITE) == 0) then
            call test_failed(error, "NaN did not raise ICE_ERR_NONFINITE")
            return
         end if
         if (res%stepped) then
            call test_failed(error, "a state with NaN was advanced")
            return
         end if
         ! the state is exactly what went in (the NaN stays a NaN, nothing finite was made up)
         worst = 0.0_wp
         call state_diff_nan(s, ref, worst)
         if (.not. (worst == 0.0_wp)) then
            call test_failed(error, "state changed by a refused step")
            return
         end if
         if (chk) then
            i = merge(1, 0, any(iand(s%ice_err, ICE_ERR_NONFINITE) /= 0))
            call check(error, i, 1, "per-cell bit set under transport_check")
            if (allocated(error)) return
         end if
         write (error_unit, "(a,i0,a,l1,a)") "  [upwind] NaN case ", pass, " check=", chk, ": bit raised, state untouched"
         call sim_destroy(s)
         call sim_destroy(ref)
      end do
   end subroutine test_nan

   subroutine state_diff_nan(a, b, worst)
      !! Bitwise equality of every state array where NaN equals NaN; `worst` becomes 1 when different.
      type(sim_t), intent(in) :: a, b
      real(wp), intent(inout) :: worst
      if (.not. same1(reshape(a%aice0, [size(a%aice0)]), reshape(b%aice0, [size(b%aice0)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%aicen, [size(a%aicen)]), reshape(b%aicen, [size(b%aicen)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%vicen, [size(a%vicen)]), reshape(b%vicen, [size(b%vicen)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%vsnon, [size(a%vsnon)]), reshape(b%vsnon, [size(b%vsnon)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%Tsfcn, [size(a%Tsfcn)]), reshape(b%Tsfcn, [size(b%Tsfcn)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%qsnon, [size(a%qsnon)]), reshape(b%qsnon, [size(b%qsnon)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%qicen, [size(a%qicen)]), reshape(b%qicen, [size(b%qicen)]))) worst = 1.0_wp
      if (.not. same1(reshape(a%sicen, [size(a%sicen)]), reshape(b%sicen, [size(b%sicen)]))) worst = 1.0_wp
   contains
      pure function same1(x, y) result(ok)
         real(wp), intent(in) :: x(:), y(:)
         logical :: ok
         ok = all((x == y) .or. (.not. ieee_is_finite(x) .and. .not. ieee_is_finite(y)))
      end function same1
   end subroutine state_diff_nan

   ! ------------------------------------------------------------------
   ! (11) slivers and negative area
   ! ------------------------------------------------------------------

   subroutine test_sliver_neg(error)
      !! (a) a sliver keeps its volume; (b) the (-puny, 0] clip; (c) a negative new area from a flow
      !! within the per-face CFL (the two axes together exceed 1), checks OFF: counted, per-cell bit,
      !! step refused, state untouched.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      real(wp) :: worst, dt
      integer :: i, j

      ! (a) sliver category with volume, checks ON, one step: accepted, volume conserved
      call sim_init(s, 16, 12, 5, 7, .true., .true., .true., .true.)
      call fill_smooth(s, 0.1_wp, .false.)
      s%aicen(8, 6, 2) = 0.5_wp*PUNY
      s%vicen(8, 6, 2) = 0.5_wp*PUNY*1.2_wp
      s%vsnon(8, 6, 2) = 0.5_wp*PUNY*0.07_wp
      call set_uniform_velocity(s, 0.3_wp, 0.1_wp)
      call sim_wrap(s)
      call totals(s, q0)
      call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call chk_ierr(error, res, "a sliver is carried, not refused")
      if (allocated(error)) return
      call totals(s, q1)
  write (error_unit, "(a,es10.3)") "  [upwind] sliver in the entry state: volumes/enthalpies conserved to ", worst_change(q0, q1, 6)
      if (.not. (worst_change(q0, q1, 6) <= 1.0e-13_wp)) then
         call test_failed(error, "sliver volume lost")
         return
      end if
      call sim_destroy(s)
      ! (b) new area in (-puny, 0]: unpacked as 0
      call sim_init(s, 16, 12, 5, 7, .true., .true., .true., .false.)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .true., .false.)
      call fill_smooth(s, 0.1_wp, .false.)
      s%aicen(8, 6, 2) = -0.5_wp*PUNY
      s%vicen(8, 6, 2) = 0.0_wp
      s%vsnon(8, 6, 2) = 0.0_wp
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call chk_ierr(error, res, "round-off negative area is tolerated")
      if (allocated(error)) return
      ref%aicen(8, 6, 2) = 0.0_wp
      ref%Tsfcn(8, 6, 2) = ref%Tf(8, 6)        ! an emptied category gets Tf
      call state_diff(s, ref, worst, .true.)
      if (.not. (worst <= 1.0e-15_wp)) then
         call test_failed(error, "(-puny, 0] area is unpacked as 0 and the rest is unchanged")
         return
      end if
      call sim_destroy(s)
      call sim_destroy(ref)
      ! (c) negative area from the two axes together, check OFF
      call sim_init(s, 24, 12, 5, 7, .true., .true., .false., .false.)
      call sim_init(ref, 24, 12, 5, 7, .true., .true., .false., .false.)
      call fill_smooth(s, 0.1_wp, .false.)
      dt = 1000.0_wp
      do j = 1, s%ny
         do i = 1, s%nx
            s%aicen(i, j, :) = merge(0.15_wp, 0.05_wp, mod(i + j, 2) == 0)
            s%vicen(i, j, :) = s%aicen(i, j, :)
            s%aice0(i, j) = 1.0_wp - sum(s%aicen(i, j, :))
         end do
      end do
      call set_uniform_velocity(s, 0.9_wp*s%dx/dt, 0.9_wp*s%dy/dt)
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      write (error_unit, "(a,i0,a)") "  [upwind] per-face CFL 0.9 on both axes, checks off: ", res%n_neg, " negative-area categories counted"
      call check(error, iand(res%ierr, ICE_ERR_TRANSPORT_NEG_AREA) /= 0, "negative area bit")
      if (allocated(error)) return
      call check(error, any(iand(s%ice_err, ICE_ERR_TRANSPORT_NEG_AREA) /= 0), "the cell error word carries it with the checks off")
      if (allocated(error)) return
      call check(error, res%stepped, .false., "refused")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      if (.not. (worst == 0.0_wp)) call test_failed(error, "state untouched")
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_sliver_neg

   subroutine test_cfl_unconditional(error)
      !! Checks OFF: a face above the CFL displacement is still an error, in the cell error word.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst

      call sim_init(s, 16, 12, 5, 7, .true., .true., .false., .false.)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .false., .false.)
      call fill_smooth(s, 0.1_wp, .false.)
      call set_uniform_velocity(s, 0.3_wp, 0.0_wp)
      s%u(9, 6) = 1.5_wp*s%dx/6000.0_wp
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call check(error, res%n_cfl, 1, "one cell above the CFL displacement")
      if (allocated(error)) return
      call check(error, any(iand(s%ice_err, ICE_ERR_TRANSPORT_CFL) /= 0), "cell bit with the checks off")
      if (allocated(error)) return
      call check(error, res%stepped, .false., "refused")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      if (.not. (worst == 0.0_wp)) call test_failed(error, "state untouched")
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_cfl_unconditional

   subroutine test_pack_floor(error)
      !! The pack with `thick_floor` against the default: bitwise identical wherever `aim > puny`;
      !! in a sliver category the default packs zeros and the floor packs the thicknesses and the
      !! volume tracers (not the surface temperature); `thick_floor = puny` is the default exactly.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      integer(int32) :: ierr
      integer :: i, j, n, nt, nbit, nsl
      real(wp) :: tr0(2, 20)

      call sim_init(s, 12, 10, 3, 4, .true., .true., .false., .false.)
      call fill_smooth(s, 0.2_wp, .false.)
      s%aicen(6, 5, 2) = 0.4_wp*PUNY
      s%vicen(6, 5, 2) = 0.4_wp*PUNY*1.5_wp
      s%vsnon(6, 5, 2) = 0.4_wp*PUNY*0.1_wp
      s%aicen(7, 5, 2) = 0.0_wp
      s%vicen(7, 5, 2) = 0.0_wp
      s%vsnon(7, 5, 2) = 0.0_wp
      call sim_map(s)
      call sim_push(s)
      call transport_pack(s%w%a, s%p%c, 4, 1, 1, s%nx, 1, s%ny, s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, &
                          s%sicen, s%qsnon, ierr)
      call transport_pack(s%w%b, s%p%c, 4, 1, 1, s%nx, 1, s%ny, s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, &
                          s%sicen, s%qsnon, ierr, thick_floor=PUNY)
      !$acc update self(s%w%a%aim, s%w%a%trm, s%w%b%aim, s%w%b%trm)
      nbit = 0
      do n = 1, 3
         do j = 1, s%ny
            do i = 1, s%nx
               if (s%w%a%aim(i, j, n) /= s%w%b%aim(i, j, n)) nbit = nbit + 1
               do nt = 1, s%w%ntrace
                  if (s%w%a%trm(i, j, nt, n) /= s%w%b%trm(i, j, nt, n)) nbit = nbit + 1
               end do
            end do
         end do
      end do
      call check(error, nbit, 0, "thick_floor = puny is the default, bitwise")
      if (allocated(error)) return
      tr0(1, 1:s%w%ntrace) = s%w%a%trm(6, 5, 1:s%w%ntrace, 2)      ! default: all zero in the sliver
      if (any(tr0(1, 1:s%w%ntrace) /= 0.0_wp)) then
         call test_failed(error, "the default packs zeros in a sliver")
         return
      end if
      call transport_pack(s%w%b, s%p%c, 4, 1, 1, s%nx, 1, s%ny, s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, &
                          s%sicen, s%qsnon, ierr, thick_floor=UPWIND_THICK_FLOOR)
      !$acc update self(s%w%b%aim, s%w%b%trm)
      nsl = 0
      do nt = 1, s%w%ntrace
         if (s%w%b%trm(6, 5, nt, 2) /= 0.0_wp) nsl = nsl + 1
      end do
      if (rd(s%w%b%trm(6, 5, 1, 2), 1.5_wp) > 1.0e-15_wp .or. rd(s%w%b%trm(6, 5, 2, 2), 0.1_wp) > 1.0e-15_wp) then
         call test_failed(error, "the floor packs the thicknesses of a sliver")
         return
      end if
      if (s%w%b%trm(6, 5, 3, 2) /= 0.0_wp) then
         call test_failed(error, "the surface temperature of a sliver is still not packed")
         return
      end if
      if (nsl < 2 + 4 + 4 + 1) then
         call test_failed(error, "the volume tracers of a sliver are packed")
         return
      end if
      if (any(s%w%b%trm(7, 5, 1:s%w%ntrace, 2) /= 0.0_wp)) then
         call test_failed(error, "an empty category packs zeros")
         return
      end if
      call sim_destroy(s)
   end subroutine test_pack_floor

   subroutine test_compare_floor(error)
      !! The conservation comparison with a floor: the negligible content change of a nearly empty
      !! category passes (the reference comparison trips on it), real losses still trip.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NC = 2, NT = 3, NQ = 9
      type(efp_t) :: b(NQ), a(NQ)
      real(wp) :: v(NQ), vv(NQ), wrel
      integer(int32) :: ierr
      integer :: q, wq

      v = [1.0_wp, 0.5_wp, 1.0e-9_wp, -1.5_wp, -2.0_wp, -0.7_wp, -1.0e-8_wp, -2.0e-9_wp, -7.0e-10_wp]
      do q = 1, NQ
         b(q) = efp_from_real(v(q))
      end do
      ! (i) negligible: 3e-12 of a 1e-8 content in a category of area 1e-9
      vv = v
      vv(7) = v(7) + 3.0e-12_wp
      do q = 1, NQ
         a(q) = efp_from_real(vv(q))
      end do
      call transport_conservation_compare(NC, NT, b, a, PUNY, ierr, wrel, wq)
      call check(error, ierr /= 0, "the reference comparison trips on a nearly empty category")
      if (allocated(error)) return
      call transport_conservation_compare_floor(NC, NT, b, a, PUNY, ierr, wrel, wq)
      call check(error, ierr, ICE_ERR_NONE, "the floor accepts the negligible change")
      if (allocated(error)) return
      ! (ii) a real loss in a populated category
      vv = v
      vv(4) = v(4)*(1.0_wp - 1.0e-6_wp)
      do q = 1, NQ
         a(q) = efp_from_real(vv(q))
      end do
      call transport_conservation_compare_floor(NC, NT, b, a, PUNY, ierr, wrel, wq)
      call check(error, ierr /= 0, "a genuine loss trips")
      if (allocated(error)) return
      call check(error, wq, 4, "and names its quantity")
      if (allocated(error)) return
      ! (iii) the whole content of the nearly empty category lost: far above the floor
      vv = v
      vv(7) = 0.0_wp
      do q = 1, NQ
         a(q) = efp_from_real(vv(q))
      end do
      call transport_conservation_compare_floor(NC, NT, b, a, PUNY, ierr, wrel, wq)
      call check(error, ierr /= 0, "losing a whole (small) category content trips")
      if (allocated(error)) return
      ! (iv) an area lost
      vv = v
      vv(2) = v(2) - 1.0e-6_wp
      do q = 1, NQ
         a(q) = efp_from_real(vv(q))
      end do
      call transport_conservation_compare_floor(NC, NT, b, a, PUNY, ierr, wrel, wq)
      call check(error, ierr /= 0, "an area loss trips")
   end subroutine test_compare_floor

   ! ------------------------------------------------------------------
   ! (12) oracle: independent host implementation of the published scheme
   ! ------------------------------------------------------------------

   pure function ref_flux(y1, y2, a, h, dt) result(f)
      real(wp), intent(in) :: y1, y2, a, h, dt
      real(wp) :: f
      f = 0.5_wp*dt*h*((a + abs(a))*y1 + (a - abs(a))*y2)
   end function ref_flux

   subroutine oracle_step(s, dt, r)
      !! The scheme on the STATE with flux arrays, as published: advect `[aice0 | aicen, vicen, vsnon,
      !! tracer products]` with the face-normal velocities times the face masks, then recover the
      !! tracers (area tracer: `P/a` if `a > puny` else `Tf` (0 on land); volume tracer: `P/V` if
      !! `V > 0` else 0).  `r` receives the new owned state (ring and everything else copied from `s`).
      type(sim_t), intent(in) :: s
      real(wp), intent(in) :: dt
      type(sim_t), intent(inout) :: r
      integer :: narr, nf, i, j, n, k, it, nt
      real(wp), allocatable :: works(:, :, :), fa(:, :), fb(:, :)
      real(wp) :: uee, vnn

      call copy_state(s, r)
      nt = 1 + 2*s%nil + 1
      narr = 1 + s%ncat*(3 + nt)
      allocate (works(s%nx, s%ny, narr), fa(s%nx, s%ny), fb(s%nx, s%ny))
      works = 0.0_wp
      works(:, :, 1) = s%aice0
      nf = 1
      do n = 1, s%ncat
         works(:, :, nf + 1) = s%aicen(:, :, n)
         works(:, :, nf + 2) = s%vicen(:, :, n)
         works(:, :, nf + 3) = s%vsnon(:, :, n)
         nf = nf + 3
         works(:, :, nf + 1) = s%aicen(:, :, n)*s%Tsfcn(:, :, n)
         do k = 1, s%nil
            works(:, :, nf + 1 + k) = s%vicen(:, :, n)*s%qicen(:, :, n, k)
            works(:, :, nf + 1 + s%nil + 1 + k) = s%vicen(:, :, n)*s%sicen(:, :, n, k)
         end do
         works(:, :, nf + 1 + s%nil + 1) = s%vsnon(:, :, n)*s%qsnon(:, :, n, 1)
         nf = nf + nt
      end do
      do it = 1, narr
         fa = 0.0_wp
         fb = 0.0_wp
         do j = s%js - 1, s%je
            do i = s%is - 1, s%ie
               uee = s%u(i + 1, j)*s%epm(i, j)
               vnn = s%v(i, j + 1)*s%npm(i, j)
               fa(i, j) = ref_flux(works(i, j, it), works(i + 1, j, it), uee, s%HTE(i, j), dt)
               fb(i, j) = ref_flux(works(i, j, it), works(i, j + 1, it), vnn, s%HTN(i, j), dt)
            end do
         end do
         do j = s%js, s%je
            do i = s%is, s%ie
               works(i, j, it) = works(i, j, it) - (fa(i, j) - fa(i - 1, j) + fb(i, j) - fb(i, j - 1))/s%tarea(i, j)
            end do
         end do
      end do
      do j = s%js, s%je
         do i = s%is, s%ie
            r%aice0(i, j) = works(i, j, 1)
            nf = 1
            do n = 1, s%ncat
               r%aicen(i, j, n) = works(i, j, nf + 1)
               r%vicen(i, j, n) = works(i, j, nf + 2)
               r%vsnon(i, j, n) = works(i, j, nf + 3)
               nf = nf + 3
               if (r%aicen(i, j, n) > PUNY) then
                  r%Tsfcn(i, j, n) = works(i, j, nf + 1)/r%aicen(i, j, n)
               else
                  r%Tsfcn(i, j, n) = s%Tf(i, j)
               end if
               do k = 1, s%nil
                  if (r%vicen(i, j, n) > 0.0_wp) then
                     r%qicen(i, j, n, k) = works(i, j, nf + 1 + k)/r%vicen(i, j, n)
                     r%sicen(i, j, n, k) = works(i, j, nf + 1 + s%nil + 1 + k)/r%vicen(i, j, n)
                  else
                     r%qicen(i, j, n, k) = 0.0_wp
                     r%sicen(i, j, n, k) = 0.0_wp
                  end if
               end do
               if (r%vsnon(i, j, n) > 0.0_wp) then
                  r%qsnon(i, j, n, 1) = works(i, j, nf + 1 + s%nil + 1)/r%vsnon(i, j, n)
               else
                  r%qsnon(i, j, n, 1) = 0.0_wp
               end if
               if (s%hm(i, j) <= 0.5_wp) r%Tsfcn(i, j, n) = 0.0_wp
               nf = nf + nt
            end do
         end do
      end do
   end subroutine oracle_step

   subroutine test_oracle(error)
      type(error_type), allocatable, intent(out) :: error
      call run_oracle(error, .false.)
   end subroutine test_oracle

   subroutine test_oracle_periodic(error)
      type(error_type), allocatable, intent(out) :: error
      call run_oracle(error, .true.)
   end subroutine test_oracle_periodic

   subroutine run_oracle(error, per)
      !! Random masked state: the library (device on the GPU build) against the oracle, every field of
      !! every cell; compared where the library defines the value (tracers of a category with no ice
      !! are left untouched by the unpack, the reference zeroes them).  Ring and untouched state exactly.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: per
      type(sim_t) :: s, o, s0
      type(transport_result_t) :: res
      integer :: i, j, n, k
      integer(int64) :: st
      real(wp) :: worst, dt, ncmp, tol
      logical :: land
      character(len=12) :: wname
      integer :: wi, wj, wn

      call sim_init(s, 30, 22, 5, 7, per, per, .true., .true.)
      call sim_init(o, 30, 22, 5, 7, per, per, .true., .true.)
      call sim_init(s0, 30, 22, 5, 7, per, per, .true., .true.)
      st = lcg_next(20240611_int64)
      do j = s%js, s%je
         do i = s%is, s%ie
            st = lcg_next(st)
            land = lcg_unit(st) < 0.08_wp .and. .not. per
            if (land) s%hm(i, j) = 0.0_wp
         end do
      end do
      call sim_face_masks(s)
      do j = 1, s%ny
         do i = 1, s%nx
            s%aice0(i, j) = 0.0_wp
            s%Tsfcn(i, j, :) = 0.0_wp
            if (s%hm(i, j) <= 0.5_wp) cycle
            do n = 1, 5
               st = lcg_next(st)
               if (lcg_unit(st) < 0.15_wp) then
                  s%aicen(i, j, n) = 0.0_wp
                  s%vicen(i, j, n) = 0.0_wp
                  s%vsnon(i, j, n) = 0.0_wp
                  s%Tsfcn(i, j, n) = s%Tf(i, j)
                  s%qicen(i, j, n, :) = 0.0_wp
                  s%sicen(i, j, n, :) = 0.0_wp
                  s%qsnon(i, j, n, 1) = 0.0_wp
                  cycle
               end if
               st = lcg_next(st); s%aicen(i, j, n) = 0.02_wp + 0.16_wp*lcg_unit(st)
               st = lcg_next(st); s%vicen(i, j, n) = s%aicen(i, j, n)*(0.1_wp + 4.0_wp*lcg_unit(st))
               st = lcg_next(st)
               if (lcg_unit(st) < 0.3_wp) then
                  s%vsnon(i, j, n) = 0.0_wp
               else
                  st = lcg_next(st); s%vsnon(i, j, n) = s%aicen(i, j, n)*0.4_wp*lcg_unit(st)
               end if
               st = lcg_next(st); s%Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
               do k = 1, 7
                  st = lcg_next(st); s%qicen(i, j, n, k) = -(1.0e8_wp + 2.0e8_wp*lcg_unit(st))
                  st = lcg_next(st); s%sicen(i, j, n, k) = 10.0_wp*lcg_unit(st)
               end do
               st = lcg_next(st); s%qsnon(i, j, n, 1) = -SHIFT*(1.05_wp + 0.1_wp*lcg_unit(st))
            end do
            s%aice0(i, j) = max(1.0_wp - sum(s%aicen(i, j, :)), 0.02_wp)
         end do
      end do
      do j = 1, s%ny
         do i = 1, s%nx + 1
            st = lcg_next(st)
            s%u(i, j) = (lcg_unit(st) - 0.5_wp)*0.9_wp
            if (abs(s%u(i, j)) < 0.05_wp) s%u(i, j) = 0.05_wp
         end do
      end do
      do j = 1, s%ny + 1
         do i = 1, s%nx
            st = lcg_next(st)
            s%v(i, j) = (lcg_unit(st) - 0.5_wp)*0.9_wp
            if (abs(s%v(i, j)) < 0.05_wp) s%v(i, j) = 0.05_wp
         end do
      end do
      if (per) call sim_wrap(s)
      call copy_state(s, s0)
      dt = 5000.0_wp                               ! |u| dt/dx <= 0.23 per axis: positive-definite
      call oracle_step(s, dt, o)
      call sim_step(s, dt, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call chk_ierr(error, res, "no bit on the random state")
      if (allocated(error)) return
      call check(error, res%stepped, "stepped")
      if (allocated(error)) return
      worst = 0.0_wp
      ncmp = 0.0_wp
      wname = "none"
      wi = 0; wj = 0; wn = 0
      do j = 1, s%ny
         do i = 1, s%nx
            if (i < s%is .or. i > s%ie .or. j < s%js .or. j > s%je) then
               ! ring: untouched
               call upd("ring aice0", s%aice0(i, j), s0%aice0(i, j), i, j, 0)
               do n = 1, 5
                  call upd("ring aicen", s%aicen(i, j, n), s0%aicen(i, j, n), i, j, n)
                  call upd("ring vicen", s%vicen(i, j, n), s0%vicen(i, j, n), i, j, n)
                  call upd("ring Tsfcn", s%Tsfcn(i, j, n), s0%Tsfcn(i, j, n), i, j, n)
               end do
               cycle
            end if
            call upd("aice0", s%aice0(i, j), o%aice0(i, j), i, j, 0)
            do n = 1, 5
               call upd("aicen", s%aicen(i, j, n), o%aicen(i, j, n), i, j, n)
               call upd("vicen", s%vicen(i, j, n), o%vicen(i, j, n), i, j, n)
               call upd("vsnon", s%vsnon(i, j, n), o%vsnon(i, j, n), i, j, n)
               call upd("Tsfcn", s%Tsfcn(i, j, n), o%Tsfcn(i, j, n), i, j, n)
               ncmp = ncmp + 4.0_wp
               if (o%aicen(i, j, n) > 0.0_wp) then
                  if (o%vsnon(i, j, n) > 0.0_wp) then
                     ! a snow-free layer comes back at the melting-point value of the shifted frame
                     ! (inert), the reference writes 0: not compared
                     call upd("qsnon", s%qsnon(i, j, n, 1), o%qsnon(i, j, n, 1), i, j, n)
                     ncmp = ncmp + 1.0_wp
                  end if
                  do k = 1, 7
                     call upd("qicen", s%qicen(i, j, n, k), o%qicen(i, j, n, k), i, j, n)
                     call upd("sicen", s%sicen(i, j, n, k), o%sicen(i, j, n, k), i, j, n)
                     ncmp = ncmp + 2.0_wp
                  end do
               end if
            end do
         end do
      end do
      if (worst > 1.0e-12_wp) write (error_unit, "(a,a,3(a,i0))") "  [upwind] oracle worst field ", trim(wname), " at i ", wi, " j ", wj, " n ", wn
      tol = 1.5e-14_wp          ! measured worst 1.05e-15 on gfortran, 1.2e-15 on the GPU build (FMA on)
      write (error_unit, "(a,es10.3,a,f10.0,a)") "  [upwind] library vs independent host scheme: worst relative difference ", &
         worst, " over ", ncmp, " field values (ring exact)"
      if (.not. (worst <= tol)) call test_failed(error, "library differs from the independent host scheme")
      call sim_destroy(s)
      call sim_destroy(o)
      call sim_destroy(s0)
   contains
      subroutine upd(name, a, b, i, j, n)
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: a, b
         integer, intent(in) :: i, j, n
         real(wp) :: x
         x = rd(a, b)
         if (x > worst) then
            worst = x
            wname = name
            wi = i; wj = j; wn = n
         end if
      end subroutine upd
   end subroutine run_oracle

   ! ------------------------------------------------------------------
   ! (13) configuration refusals
   ! ------------------------------------------------------------------

   subroutine test_cfg(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst

      call sim_init(s, 16, 12, 5, 7, .true., .true., .false., .false.)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .false., .false.)
      call fill_smooth(s, 0.1_wp, .false.)
      call set_uniform_velocity(s, 0.3_wp, 0.1_wp)
      call sim_wrap(s)
      call copy_state(s, ref)
      s%p%transport%advection = ICE_ADVECT_REMAP
      call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call check(error, res%cfg, TRANSPORT_CFG_PARAMS, "remap is refused, never run as upwind")
      if (allocated(error)) return
      s%p%transport%advection = ICE_ADVECT_UPWIND
      s%p%ncat = 4
      call sim_step(s, 6000.0_wp, ICE_TRANSPORT_POST_NONE, res, wrap=.false.)
      call check(error, res%cfg, TRANSPORT_CFG_PARAMS, "ncat mismatch refused")
      if (allocated(error)) return
      s%p%ncat = 5
      ! range without room for the ghost ring
      call ice_step_transport(s%w, s%p, 6000.0_wp, ICE_TRANSPORT_POST_NONE, 1, s%ie, s%js, s%je, s%u, s%v, s%HTE, s%HTN, &
                              s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, s%aicen, s%vicen, &
                              s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, &
                              s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res)
      call check(error, res%cfg, TRANSPORT_CFG_RANGE, "range without a ghost ring refused")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      if (.not. (worst == 0.0_wp)) call test_failed(error, "refused calls leave the state alone")
      call check(error, TRANSPORT_CFG_OK == 0 .and. TRANSPORT_CFG_WORKSPACE /= 0, "codes")
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_cfg

   ! ------------------------------------------------------------------
   ! (14) throughput
   ! ------------------------------------------------------------------

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      type(transport_result_t) :: res
      integer :: nxp, nyp, i, j, n, k, rep, nbad, nneg, cfg
      integer, parameter :: NREP = 5
      integer(int64) :: t0, t1, t2, rate
      real(wp) :: tstep, tker, cells

#ifdef RDB_GPU_OFFLOAD
      nxp = 1440; nyp = 1152
#else
      nxp = 512; nyp = 512
#endif
      call sim_init(s, nxp, nyp, 5, 7, .true., .true., .false., .false.)
      do n = 1, 5
         do j = 1, s%ny
            do i = 1, s%nx
               s%aicen(i, j, n) = 0.04_wp + 0.01_wp*real(mod(i + 2*j + n, 7), wp)
               s%vicen(i, j, n) = s%aicen(i, j, n)*(0.5_wp + 0.3_wp*real(n, wp))
               s%vsnon(i, j, n) = 0.1_wp*s%aicen(i, j, n)
               s%Tsfcn(i, j, n) = -10.0_wp
               s%qsnon(i, j, n, 1) = -1.1e8_wp
               do k = 1, 7
                  s%qicen(i, j, n, k) = -2.0e8_wp
                  s%sicen(i, j, n, k) = 5.0_wp
               end do
            end do
         end do
      end do
      s%aice0 = 0.5_wp
      call set_uniform_velocity(s, 0.3_wp, -0.2_wp)
      call sim_map(s)
      call sim_push(s)
      call run_one(s, res)
      call system_clock(t0, rate)
      do rep = 1, NREP
         call run_one(s, res)
      end do
      call system_clock(t1)
      do rep = 1, NREP
         call upwind_advect(s%w%a, s%w%b, s%is, s%ie, s%js, s%je, 6000.0_wp, PUNY, s%u, s%v, s%HTE, s%HTN, s%tarea, &
                            s%epm, s%npm, nbad, nneg, cfg)
      end do
      call system_clock(t2)
      call sim_pull(s)
      tstep = real(t1 - t0, wp)/real(rate, wp)/real(NREP, wp)
      tker = real(t2 - t1, wp)/real(rate, wp)/real(NREP, wp)
      cells = real(s%nx, wp)*real(s%ny, wp)*5.0_wp
      write (error_unit, "(a,i0,a,i0,a,es10.3,a,es10.3,a,i0,a,es10.3,a)") "  [upwind] ", nxp, "x", nyp, "x5x7: full step (check off, aggregate) ", &
         cells/tstep, " cat-cells/s (", tstep*1000.0_wp, " ms, ", res%n_launch, " launches); advection kernel alone ", &
         cells/tker, " cat-cells/s"
      call chk_ierr(error, res, "clean")
      if (.not. allocated(error)) call check(error, ieee_is_finite(cells/tstep) .and. ieee_is_finite(cells/tker), "finite")
      call sim_destroy(s)
   end subroutine test_throughput

   subroutine run_one(s, res)
      !! One driver call on the device-resident state (no host wrap, no transfers).
      type(sim_t), intent(inout) :: s
      type(transport_result_t), intent(out) :: res
      call ice_step_transport(s%w, s%p, 6000.0_wp, ICE_TRANSPORT_POST_AGGREGATE, s%is, s%ie, s%js, s%je, s%u, s%v, &
                              s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, &
                              s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, &
                              s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res)
   end subroutine run_one

end module test_cice_upwind
