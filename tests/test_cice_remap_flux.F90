!! Sea-ice incremental-remapping transport (D10.5): the flux integrals over the departure triangles,
!! the conservative update, the batched buffers and the remap branch of the transport driver.
module test_cice_remap_flux
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Tests written for the roundabout port of the CICE incremental-remapping transport (CICE
   ! Consortium, https://github.com/CICE-Consortium); the host oracle integrates the
   ! reconstructed polynomials over each triangle by exact monomial formulas in barycentric
   ! coordinates (independent of the kernel's quadrature ONLY: it re-reads the kernel's triangles,
   ! reconstruction and gradient mask) and re-implements the published update.  Modified by the roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   !! **What is pinned.**  (1) every edge flux of every category and tracer equals a host
   !! integration of the reconstruction over the stored triangles by exact barycentric monomial
   !! means (independent of the quadrature only, not of the geometry or the reconstruction), on a random flow that runs the geometry's cancelling-pair branches, and the whole
   !! new state (area, every tracer, every cell, every category) equals the host re-implementation
   !! of the update, the ring of the output untouched; (2) area, ice and snow volume, ice and snow
   !! energy, salt and the surface-temperature content are conserved by EFP sums over many steps
   !! of uniform translation, solid-body rotation and a deformation field; (3) no new extremum of
   !! any tracer beyond the independent 5x5 bound on sharp noisy fronts under a random shearing
   !! flow, and the driver's own check agrees; (4) one-cell translation reproduces the shifted
   !! field to round-off; (5) half-cell translation is second order where upwind is first order
   !! (same case, error ratios under refinement); (6) a cone after one revolution keeps its peak
   !! (remap against upwind); (7) in a convergent / divergent flow a uniform field changes by the
   !! flux convergence and a constant tracer stays constant to round-off; (8) closed faces and land
   !! stay untouched; (9) empty domain and zero velocity are bitwise no-ops; (10) NaN velocity or
   !! field raises the bit and leaves the state; (11) CFL violations are refused; (12) the category
   !! batch does not change one bit; (13) configuration refusals and the halo requirement; (14)
   !! throughput, launches per step and memory.
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_from_real, efp_plus, efp_real_diff, efp_to_real
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_ADVECT_UPWIND, ICE_ADVECT_REMAP
   use rdb_ice_itd_params, only: itd_category_bounds
   use rdb_ice_remap_geom, only: remap_src_di, remap_src_dj, REMAP_NSLOT, remap_locate_east, remap_locate_edge, remap_geometry
   use rdb_ice_remap_recon_batch, only: remap_recon_batch, remap_input_check
   use rdb_ice_remap_flux, only: remap_flux_integrals, remap_update
   use rdb_ice_transport_cice, only: transport_cice_ws_t, transport_result_t, ice_step_transport, &
                                     ICE_TRANSPORT_POST_NONE, TRANSPORT_CFG_OK, TRANSPORT_CFG_PARAMS, TRANSPORT_CFG_RANGE, &
                                     TRANSPORT_CFG_REMAP
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_remap_flux_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TWOPI = 6.283185307179586476925286766559_wp
   real(wp), parameter :: SHIFT = 330.0_wp*3.337e5_wp
   real(wp), parameter :: HMID(5) = [0.15_wp, 0.5_wp, 0.95_wp, 1.6_wp, 2.5_wp]
   integer, parameter :: NG = 3
      !! Ghost cells: the remap's minimum.
   real(wp), parameter :: TOL_T7 = 1.0e-8_wp
      !! The surface-temperature content (an area tracer) of the cells whose area is `<= puny` is dropped by
      !! design (at most `puny` times the temperature per event); upwind loses 2.5e-10 on the same cone.
   integer, parameter :: FLOW_UNIFORM = 1, FLOW_ROTATION = 2, FLOW_CELLULAR = 3, FLOW_CONVERGENT = 4, FLOW_RANDOM = 5

   type :: sim_t
      !! A host-owned C-grid ice state, its metrics and the driver workspace.  Every array is
      !! mapped once; after each `sim_step` the host copy holds the result.
      integer :: nx = 0, ny = 0, ncat = 0, nil = 0
      integer :: nxp = 0, nyp = 0, is = 0, ie = 0, js = 0, je = 0
      logical :: per_x = .true., per_y = .true., mapped = .false.
      real(wp) :: dx = 1.0e4_wp, dy = 1.0e4_wp
      type(ice_params_t) :: p
      type(transport_cice_ws_t) :: w
      real(wp), allocatable :: u(:, :), v(:, :), HTE(:, :), HTN(:, :), tarea(:, :), dxT(:, :), dyT(:, :)
      real(wp), allocatable :: dxU(:, :), dyU(:, :), earea(:, :), narea(:, :), uvm(:, :)
      real(wp), allocatable :: epm(:, :), npm(:, :), hm(:, :), Tf(:, :), hin_max(:)
      real(wp), allocatable :: aice0(:, :), aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), Tsfc_c(:, :), qice_c(:, :, :), sice_c(:, :, :)
      real(wp), allocatable :: qsno_c(:, :, :), fresh(:, :), fsalt(:, :), fhocn(:, :)
      integer(int32), allocatable :: ice_err(:, :)
   end type sim_t

contains

   subroutine collect_cice_remap_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fluxes_and_update_vs_exact_host_oracle", test_oracle), &
                  new_unittest("conservation_translation_rotation_deformation", test_conservation), &
                  new_unittest("monotone_noisy_fronts_cancelling_pairs", test_monotone), &
                  new_unittest("one_cell_translation_reproduces_shift", test_one_cell), &
                  new_unittest("half_cell_second_order_vs_upwind", test_half_cell), &
                  new_unittest("rotation_peak_retention_vs_upwind", test_rotation), &
                  new_unittest("divergence_consistency_constant_tracer", test_divergence), &
                  new_unittest("walls_and_land_stay_closed", test_land), &
                  new_unittest("empty_domain_and_zero_velocity_noops", test_noops), &
                  new_unittest("nan_raises_bit_state_untouched", test_nan), &
                  new_unittest("cfl_refused", test_cfl), &
                  new_unittest("category_batch_is_bit_identical", test_batch), &
                  new_unittest("tiny_departure_closes_the_flux_adjustment", test_tiny), &
                  new_unittest("nonuniform_metrics_conserve_and_stay_monotone", test_nonuniform), &
                  new_unittest("configuration_refusals_and_halo", test_cfg), &
                  new_unittest("throughput_launches_memory", test_throughput) &
                  ]
   end subroutine collect_cice_remap_flux_tests

   ! ------------------------------------------------------------------
   ! simulation helper
   ! ------------------------------------------------------------------

   subroutine sim_init(s, nxp, nyp, ncat, nil, per_x, per_y, check_on, want_bounds, dx, nbatch, bytes_max)
      !! Uniform Cartesian grid `dx = dy` with an `NG`-cell ghost ring, all faces open (periodic) or
      !! closed at the domain edge (the ghost band is land).  `nbatch`: category indices per flux batch.
      type(sim_t), intent(out) :: s
      integer, intent(in) :: nxp, nyp, ncat, nil
      logical, intent(in) :: per_x, per_y, check_on, want_bounds
      real(wp), intent(in), optional :: dx
      integer, intent(in), optional :: nbatch
      integer(int64), intent(in), optional :: bytes_max
      integer :: ierr, status

      s%nxp = nxp; s%nyp = nyp; s%ncat = ncat; s%nil = nil
      s%per_x = per_x; s%per_y = per_y
      if (present(dx)) then
         s%dx = dx
         s%dy = dx
      end if
      s%nx = nxp + 2*NG; s%ny = nyp + 2*NG
      s%is = NG + 1; s%ie = NG + nxp; s%js = NG + 1; s%je = NG + nyp
      s%p = ice_params_default(ICE_CONSTANTS_ICEPACK, ncat, nil, 1)
      s%p%transport%advection = ICE_ADVECT_REMAP
      s%p%transport%check = check_on
      if (present(nbatch)) then
         call s%w%init(s%nx, s%ny, ncat, nil, 1, want_bounds, ierr, remap=.true., remap_nbatch=nbatch)
      else if (present(bytes_max)) then
         call s%w%init(s%nx, s%ny, ncat, nil, 1, want_bounds, ierr, remap=.true., remap_flux_bytes_max=bytes_max)
      else
         call s%w%init(s%nx, s%ny, ncat, nil, 1, want_bounds, ierr, remap=.true.)
      end if
      allocate (s%u(s%nx + 1, s%ny), s%v(s%nx, s%ny + 1))
      allocate (s%HTE(s%nx, s%ny), s%HTN(s%nx, s%ny), s%tarea(s%nx, s%ny), s%dxT(s%nx, s%ny), s%dyT(s%nx, s%ny))
      allocate (s%dxU(s%nx, s%ny), s%dyU(s%nx, s%ny), s%earea(s%nx, s%ny), s%narea(s%nx, s%ny), s%uvm(s%nx, s%ny))
      allocate (s%epm(s%nx, s%ny), s%npm(s%nx, s%ny), s%hm(s%nx, s%ny), s%Tf(s%nx, s%ny), s%hin_max(0:ncat))
      allocate (s%aice0(s%nx, s%ny), s%aicen(s%nx, s%ny, ncat), s%vicen(s%nx, s%ny, ncat), s%vsnon(s%nx, s%ny, ncat))
      allocate (s%Tsfcn(s%nx, s%ny, ncat), s%qicen(s%nx, s%ny, ncat, nil), s%sicen(s%nx, s%ny, ncat, nil))
      allocate (s%qsnon(s%nx, s%ny, ncat, 1))
      allocate (s%aice(s%nx, s%ny), s%vice(s%nx, s%ny), s%vsno(s%nx, s%ny), s%Tsfc_c(s%nx, s%ny))
      allocate (s%qice_c(s%nx, s%ny, nil), s%sice_c(s%nx, s%ny, nil), s%qsno_c(s%nx, s%ny, 1))
      allocate (s%fresh(s%nx, s%ny), s%fsalt(s%nx, s%ny), s%fhocn(s%nx, s%ny), s%ice_err(s%nx, s%ny))
      s%u = 0.0_wp; s%v = 0.0_wp
      s%HTE = s%dy; s%HTN = s%dx; s%tarea = s%dx*s%dy; s%dxT = s%dx; s%dyT = s%dy
      s%dxU = s%dx; s%dyU = s%dy; s%earea = s%dx*s%dy; s%narea = s%dx*s%dy
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
      !! masks from it.  Land blocks are put in `hm` by the caller, then `sim_face_masks` is called again.
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
      !! `epm`, `npm` (both cells wet) and `uvm` (all four cells wet) from `hm`.
      type(sim_t), intent(inout) :: s
      integer :: i, j

      s%epm = 0.0_wp
      s%npm = 0.0_wp
      s%uvm = 0.0_wp
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
      do j = 1, s%ny - 1
         do i = 1, s%nx - 1
            s%uvm(i, j) = s%hm(i, j)*s%hm(i + 1, j)*s%hm(i, j + 1)*s%hm(i + 1, j + 1)
         end do
      end do
   end subroutine sim_face_masks

   subroutine sim_map(s)
      type(sim_t), intent(inout) :: s
      if (s%mapped) return
      !$acc enter data copyin(s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max)
      !$acc enter data copyin(s%dxU, s%dyU, s%earea, s%narea, s%uvm)
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
      !$acc exit data delete(s%dxU, s%dyU, s%earea, s%narea, s%uvm)
      !$acc exit data delete(s%aice0, s%aicen, s%vicen, s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon)
      !$acc exit data delete(s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, s%sice_c, s%qsno_c)
      !$acc exit data delete(s%fresh, s%fsalt, s%fhocn, s%ice_err)
      s%mapped = .false.
   end subroutine sim_unmap

   subroutine sim_push(s)
      !! Host to device: every array the host may have edited.
      type(sim_t), intent(inout) :: s
      !$acc update device(s%u, s%v, s%HTE, s%HTN, s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max)
      !$acc update device(s%dxU, s%dyU, s%earea, s%narea, s%uvm)
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

   subroutine sim_step(s, dt, res, wrap)
      !! Wrap the ghost ring on the host (periodic axes; the host holds the latest state), push,
      !! run the driver, pull.  `wrap = .false.` skips the wrap (the caller prepared the ring).
      type(sim_t), intent(inout) :: s
      real(wp), intent(in) :: dt
      type(transport_result_t), intent(out) :: res
      logical, intent(in), optional :: wrap

      if (.not. s%mapped) call sim_map(s)
      if (.not. present(wrap)) then
         call sim_wrap(s)
      else if (wrap) then
         call sim_wrap(s)
      end if
      call sim_push(s)
      call ice_step_transport(s%w, s%p, dt, ICE_TRANSPORT_POST_NONE, s%is, s%ie, s%js, s%je, s%u, s%v, s%HTE, s%HTN, &
                              s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, s%aicen, s%vicen, &
                              s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, &
                              s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res, &
                              dxU=s%dxU, dyU=s%dyU, earea=s%earea, narea=s%narea, uvm=s%uvm)
      call sim_pull(s)
   end subroutine sim_step

   subroutine chk_ierr(error, res, msg)
      !! `res%ierr == 0` and the call accepted, printing the whole result when not.
      type(error_type), allocatable, intent(inout) :: error
      type(transport_result_t), intent(in) :: res
      character(len=*), intent(in) :: msg

      if (res%ierr /= ICE_ERR_NONE .or. res%cfg /= TRANSPORT_CFG_OK) then
         write (error_unit, "(a,a,a,z0,a,i0,a,5(a,i0),a,es10.3,a,i0)") "  [remap] ", msg, ": ierr 0x", res%ierr, " cfg ", &
            res%cfg, ";", " n_cfl ", res%n_cfl, " n_geom ", res%n_geom, " n_nonfinite ", res%n_nonfinite, " n_neg ", res%n_neg, &
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

   ! ------------------------------------------------------------------
   ! flows and fields
   ! ------------------------------------------------------------------

   pure function cps(s, dt, c) result(v)
      !! Velocity [m/s] that moves a field `c` cells in one step of `dt`.
      type(sim_t), intent(in) :: s
      real(wp), intent(in) :: dt, c
      real(wp) :: v
      v = c*s%dx/dt
   end function cps

   subroutine set_velocity(s, kind, a, b)
      !! The face velocities (ghost ring included).  `FLOW_UNIFORM`: `u = a`, `v = b` [m/s].
      !! `FLOW_ROTATION`: solid-body rotation of angular rate `a` [rad/s] about the domain centre, tapered
      !! to zero outside 0.35 of the half-width (so it is periodic).  `FLOW_CELLULAR`: the discretely
      !! non-divergent cellular flow `psi = P sin(kx x) sin(ky y)` with peak speed `a`.  `FLOW_CONVERGENT`:
      !! `u = a sin(2 pi x/L)`, `v = b sin(2 pi y/L)`, zero on closed faces.  `FLOW_RANDOM`: the cellular
      !! flow of peak speed `a` plus a random periodic stream function of face speeds up to `b` (both
      !! discretely non-divergent: `u = (psi_N - psi_S)/dy`, `v = -(psi_E - psi_W)/dx` on corner values),
      !! so the flow shears sharply without a divergence the scheme's positivity does not allow.
      type(sim_t), intent(inout) :: s
      integer, intent(in) :: kind
      real(wp), intent(in) :: a, b
      integer :: i, j
      integer(int64) :: st
      real(wp) :: x, y, xc, yc, lx, ly, r, kx, ky, pk
      real(wp), allocatable :: rnd(:, :)

      lx = real(s%nxp, wp)*s%dx
      ly = real(s%nyp, wp)*s%dy
      xc = 0.5_wp*lx
      yc = 0.5_wp*ly
      kx = TWOPI/lx
      ky = TWOPI/ly
      if (kind == FLOW_RANDOM) then
         allocate (rnd(s%nxp, s%nyp))
         st = lcg_next(4242_int64)
         do j = 1, s%nyp
            do i = 1, s%nxp
               st = lcg_next(st)
               rnd(i, j) = 2.0_wp*lcg_unit(st) - 1.0_wp
            end do
         end do
      end if
      do j = 1, s%ny
         do i = 1, s%nx + 1
            x = real(i - s%is, wp)*s%dx
            y = (real(j - s%js, wp) + 0.5_wp)*s%dy
            select case (kind)
            case (FLOW_UNIFORM)
               s%u(i, j) = a
            case (FLOW_ROTATION)
               r = sqrt((x - xc)**2 + (y - yc)**2)
               s%u(i, j) = -a*(y - yc)*taper(r, 0.35_wp*min(lx, ly))
            case (FLOW_CELLULAR, FLOW_RANDOM)
               pk = a/ky
               s%u(i, j) = (pk*sin(kx*x)*sin(ky*(y + 0.5_wp*s%dy)) - pk*sin(kx*x)*sin(ky*(y - 0.5_wp*s%dy)))/s%dy
               if (kind == FLOW_RANDOM) s%u(i, j) = s%u(i, j) + 0.5_wp*b*(psi_r(i - 1, j) - psi_r(i - 1, j - 1))
            case (FLOW_CONVERGENT)
               s%u(i, j) = a*sin(TWOPI*x/lx)
            case default
               s%u(i, j) = 0.0_wp
            end select
            if (kind == FLOW_CONVERGENT .and. (i <= s%is .or. i > s%ie)) s%u(i, j) = 0.0_wp
         end do
      end do
      do j = 1, s%ny + 1
         do i = 1, s%nx
            x = (real(i - s%is, wp) + 0.5_wp)*s%dx
            y = real(j - s%js, wp)*s%dy
            select case (kind)
            case (FLOW_UNIFORM)
               s%v(i, j) = b
            case (FLOW_ROTATION)
               r = sqrt((x - xc)**2 + (y - yc)**2)
               s%v(i, j) = a*(x - xc)*taper(r, 0.35_wp*min(lx, ly))
            case (FLOW_CELLULAR, FLOW_RANDOM)
               pk = a/kx
               s%v(i, j) = -(pk*sin(kx*(x + 0.5_wp*s%dx))*sin(ky*y) - pk*sin(kx*(x - 0.5_wp*s%dx))*sin(ky*y))/s%dx
               if (kind == FLOW_RANDOM) s%v(i, j) = s%v(i, j) - 0.5_wp*b*(psi_r(i, j - 1) - psi_r(i - 1, j - 1))
            case (FLOW_CONVERGENT)
               s%v(i, j) = b*sin(TWOPI*y/ly)
            case default
               s%v(i, j) = 0.0_wp
            end select
            if (kind == FLOW_CONVERGENT .and. (j <= s%js .or. j > s%je)) s%v(i, j) = 0.0_wp
         end do
      end do
      if (.not. (s%per_x .and. s%per_y)) call close_faces(s)
      if (allocated(rnd)) deallocate (rnd)
   contains
      pure function psi_r(ic, jc) result(p)
         !! Random periodic corner stream function in units of `dx*dx/dt`-free face speed (|u| <= b).
         integer, intent(in) :: ic, jc
         real(wp) :: p
         p = rnd(1 + modulo(ic - s%is, s%nxp), 1 + modulo(jc - s%js, s%nyp))
      end function psi_r
      pure function taper(rad, r1) result(f)
         real(wp), intent(in) :: rad, r1
         real(wp) :: f
         if (rad <= r1) then
            f = 1.0_wp
         else if (rad >= 1.4_wp*r1) then
            f = 0.0_wp
         else
            f = cos(0.25_wp*TWOPI*(rad - r1)/(0.4_wp*r1))**2
         end if
      end function taper
   end subroutine set_velocity

   subroutine close_faces(s)
      !! Zero the normal velocity on every closed face (the dynamics hands over zeros there).
      type(sim_t), intent(inout) :: s
      integer :: i, j
      do j = 1, s%ny
         do i = 1, s%nx - 1
            if (s%epm(i, j) == 0.0_wp) s%u(i + 1, j) = 0.0_wp
         end do
         s%u(1, j) = 0.0_wp
         s%u(s%nx + 1, j) = 0.0_wp
      end do
      do j = 1, s%ny - 1
         do i = 1, s%nx
            if (s%npm(i, j) == 0.0_wp) s%v(i, j + 1) = 0.0_wp
         end do
      end do
      s%v(:, 1) = 0.0_wp
      s%v(:, s%ny + 1) = 0.0_wp
   end subroutine close_faces

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
                  s%vicen(i, j, n) = a*HMID(n)*(1.0_wp + 0.05_wp*c1)
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

   subroutine fill_front(s, seed)
      !! Sharp, noisy fields: every category holds a rectangular patch of area 0.15 (the patches overlap
      !! and differ per category) on an ice-free background, each cell with its own random thickness,
      !! temperature, enthalpies and salinity, and snow absent in about a third of the cells.
      type(sim_t), intent(inout) :: s
      integer(int64), intent(in) :: seed
      integer :: i, j, n, k, i0, i1, j0, j1
      integer(int64) :: st
      real(wp) :: a, u1, hh, hs

      st = lcg_next(seed)
      s%aicen = 0.0_wp
      s%vicen = 0.0_wp
      s%vsnon = 0.0_wp
      s%Tsfcn = -1.8_wp
      s%qicen = 0.0_wp
      s%sicen = 0.0_wp
      s%qsnon = 0.0_wp
      do n = 1, s%ncat
         i0 = s%is + 2 + (n - 1)*(s%nxp/12)
         i1 = min(s%ie - 2, i0 + s%nxp/2)
         j0 = s%js + 2 + mod(n*3, 5)
         j1 = min(s%je - 2, j0 + s%nyp/2)
         do j = j0, j1
            do i = i0, i1
               st = lcg_next(st)
               u1 = lcg_unit(st)
               a = 0.15_wp
               hh = HMID(min(n, 5))*(0.7_wp + 0.6_wp*u1)
               s%aicen(i, j, n) = a
               s%vicen(i, j, n) = a*hh
               st = lcg_next(st)
               if (lcg_unit(st) > 0.33_wp) then
                  st = lcg_next(st)
                  hs = 0.02_wp + 0.1_wp*lcg_unit(st)
                  s%vsnon(i, j, n) = a*hs
                  st = lcg_next(st)
                  s%qsnon(i, j, n, 1) = -SHIFT*(1.0_wp + 0.2_wp*lcg_unit(st))
               end if
               st = lcg_next(st)
               s%Tsfcn(i, j, n) = -2.0_wp - 20.0_wp*lcg_unit(st)
               do k = 1, s%nil
                  st = lcg_next(st)
                  s%qicen(i, j, n, k) = -(1.5e8_wp + 1.0e8_wp*lcg_unit(st))
                  st = lcg_next(st)
                  s%sicen(i, j, n, k) = 2.0_wp + 6.0_wp*lcg_unit(st)
               end do
            end do
         end do
      end do
      do j = 1, s%ny
         do i = 1, s%nx
            s%aice0(i, j) = 1.0_wp - sum(s%aicen(i, j, :))
         end do
      end do
   end subroutine fill_front

   subroutine fill_cone(s, n0, cx, cy, rad, peak)
      !! A cone of area `peak (1 - r/rad)` (in cell units) in category `n0`, every tracer uniform.
      type(sim_t), intent(inout) :: s
      integer, intent(in) :: n0
      real(wp), intent(in) :: cx, cy, rad, peak
      integer :: i, j, k
      real(wp) :: r, a

      s%aicen = 0.0_wp
      s%vicen = 0.0_wp
      s%vsnon = 0.0_wp
      s%Tsfcn = -1.8_wp
      s%qicen = 0.0_wp
      s%sicen = 0.0_wp
      s%qsnon = 0.0_wp
      do j = 1, s%ny
         do i = 1, s%nx
            r = sqrt((real(i - s%is, wp) + 0.5_wp - cx)**2 + (real(j - s%js, wp) + 0.5_wp - cy)**2)
            a = peak*max(0.0_wp, 1.0_wp - r/rad)
            s%aicen(i, j, n0) = a
            if (a > 0.0_wp) then
               s%vicen(i, j, n0) = a*HMID(min(n0, 5))
               s%vsnon(i, j, n0) = 0.1_wp*a*HMID(min(n0, 5))
               s%Tsfcn(i, j, n0) = -8.0_wp
               s%qsnon(i, j, n0, 1) = -SHIFT*1.1_wp
               do k = 1, s%nil
                  s%qicen(i, j, n0, k) = -2.0e8_wp
                  s%sicen(i, j, n0, k) = 5.0_wp
               end do
            end if
            s%aice0(i, j) = 1.0_wp - sum(s%aicen(i, j, :))
         end do
      end do
   end subroutine fill_cone

   ! ------------------------------------------------------------------
   ! independent host oracle: exact triangle integrals and the update
   ! ------------------------------------------------------------------

   pure function fact(n) result(f)
      integer, intent(in) :: n
      real(wp) :: f
      integer :: k
      f = 1.0_wp
      do k = 2, n
         f = f*real(k, wp)
      end do
   end function fact

   pure function tri_mean(a, b, x, y) result(m)
      !! EXACT mean over the triangle `(x(k), y(k))` of the monomial `x^a y^b`: expand each factor in
      !! barycentric coordinates, `x = sum_k l_k x_k`, and use `mean(l1^c1 l2^c2 l3^c3) =
      !! 2 c1! c2! c3!/(c1 + c2 + c3 + 2)!` over every choice of a vertex per factor.
      integer, intent(in) :: a, b
      real(wp), intent(in) :: x(3), y(3)
      real(wp) :: m
      integer :: d, t, tt, k, v, c(3), ntup
      real(wp) :: prod

      d = a + b
      ntup = 3**d
      m = 0.0_wp
      do t = 0, ntup - 1
         tt = t
         c = 0
         prod = 1.0_wp
         do k = 1, d
            v = mod(tt, 3) + 1
            tt = tt/3
            c(v) = c(v) + 1
            if (k <= a) then
               prod = prod*x(v)
            else
               prod = prod*y(v)
            end if
         end do
         m = m + prod*2.0_wp*fact(c(1))*fact(c(2))*fact(c(3))/fact(d + 2)
      end do
   end function tri_mean

   pure function plin(c0, cx, cy) result(p)
      real(wp), intent(in) :: c0, cx, cy
      real(wp) :: p(0:3, 0:3)
      p = 0.0_wp
      p(0, 0) = c0
      p(1, 0) = cx
      p(0, 1) = cy
   end function plin

   pure function pmul(p, q) result(r)
      !! Product of two polynomials in `x, y` (coefficients `(a, b)` of `x^a y^b`), total degree `<= 3`.
      real(wp), intent(in) :: p(0:3, 0:3), q(0:3, 0:3)
      real(wp) :: r(0:3, 0:3)
      integer :: a1, b1, a2, b2
      r = 0.0_wp
      do b1 = 0, 3
         do a1 = 0, 3
            if (a1 + b1 > 3 .or. p(a1, b1) == 0.0_wp) cycle
            do b2 = 0, 3
               do a2 = 0, 3
                  if (a2 + b2 > 3 .or. a1 + a2 + b1 + b2 > 3) cycle
                  r(a1 + a2, b1 + b2) = r(a1 + a2, b1 + b2) + p(a1, b1)*q(a2, b2)
               end do
            end do
         end do
      end do
   end function pmul

   pure function pint(p, mom) result(v)
      !! Mean of the polynomial given the monomial means `mom(a, b)`.
      real(wp), intent(in) :: p(0:3, 0:3), mom(0:3, 0:3)
      real(wp) :: v
      integer :: a, b
      v = 0.0_wp
      do b = 0, 3
         do a = 0, 3
            if (a + b <= 3) v = v + p(a, b)*mom(a, b)
         end do
      end do
   end function pint

   pure function pintabs(p, mom) result(v)
      !! Sum of the absolute terms of `pint` (a bound on its round-off scale).
      real(wp), intent(in) :: p(0:3, 0:3), mom(0:3, 0:3)
      real(wp) :: v
      integer :: a, b
      v = 0.0_wp
      do b = 0, 3
         do a = 0, 3
            if (a + b <= 3) v = v + abs(p(a, b))*abs(mom(a, b))
         end do
      end do
   end function pintabs

   subroutine host_face_flux(nx, ny, ncat, ntrace, i, j, n, src, xp, yp, area, aim, mx, my, tc, tx, ty, ttype, depend, &
                             fm, ft, sc)
      !! Flux of category `n` through the edge `(i, j)` by exact means of the polynomials over
      !! each stored triangle (independent of the kernel's quadrature only: the triangles, the
      !! reconstruction and the gradient mask are read from the kernel's own output).  `sc` receives the sum of the absolute contributions (a scale for the
      !! tolerance) per field (`sc(0)` mass, `sc(nt)` tracers).
      integer, intent(in) :: nx, ny, ncat, ntrace, i, j, n
      integer(int8), intent(in) :: src(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), mx(nx, ny, 0:ncat), my(nx, ny, 0:ncat)
      real(wp), intent(in) :: tc(nx, ny, ntrace, 0:ncat), tx(nx, ny, ntrace, 0:ncat), ty(nx, ny, ntrace, 0:ncat)
      integer, intent(in) :: ttype(ntrace), depend(ntrace)
      real(wp), intent(out) :: fm, ft(ntrace), sc(0:ntrace)
      integer :: s, nt, a, b, i2, j2, p
      real(wp) :: x(3), y(3), mom(0:3, 0:3), pm(0:3, 0:3), pt(0:3, 0:3), pp(0:3, 0:3), contrib

      fm = 0.0_wp
      ft = 0.0_wp
      sc = 0.0_wp
      do s = 1, REMAP_NSLOT
         if (src(i, j, s) == 0_int8) cycle
         i2 = i + remap_src_di(src(i, j, s))
         j2 = j + remap_src_dj(src(i, j, s))
         x = xp(i, j, 1:3, s)
         y = yp(i, j, 1:3, s)
         do b = 0, 3
            do a = 0, 3
               mom(a, b) = 0.0_wp
               if (a + b <= 3) mom(a, b) = tri_mean(a, b, x, y)
            end do
         end do
         if (aim(i2, j2, n) > PUNY) then
            pm = plin(aim(i2, j2, n), mx(i2, j2, n), my(i2, j2, n))
         else
            pm = plin(0.0_wp, 0.0_wp, 0.0_wp)    ! the gradients are read through mc (module docstring)
         end if
         contrib = area(i, j, s)*pint(pm, mom)
         fm = fm + contrib
         sc(0) = sc(0) + abs(area(i, j, s))*pintabs(pm, mom)
         if (n >= 1) then
            do nt = 1, ntrace
               pt = plin(tc(i2, j2, nt, n), tx(i2, j2, nt, n), ty(i2, j2, nt, n))
               if (ttype(nt) == 1) then
                  pp = pmul(pm, pt)
               else
                  p = depend(nt)
                  pp = pmul(pmul(pm, plin(tc(i2, j2, p, n), tx(i2, j2, p, n), ty(i2, j2, p, n))), pt)
               end if
               contrib = area(i, j, s)*pint(pp, mom)
               ft(nt) = ft(nt) + contrib
               sc(nt) = sc(nt) + abs(area(i, j, s))*pintabs(pp, mom)
            end do
         end if
      end do
   end subroutine host_face_flux

   ! ------------------------------------------------------------------
   ! (1) fluxes and update against the exact host oracle
   ! ------------------------------------------------------------------

   function tol_flux() result(t)
      !! Flux vs exact integral, relative to the sum of the absolute contributions.
      real(wp) :: t
#ifdef RDB_GPU_OFFLOAD
      t = 2.0e-13_wp          ! FMA contraction is on for the device build
#else
      t = 1.0e-14_wp          ! measured worst 1.1e-15 (gfortran, no contraction)
#endif
   end function tol_flux

   function tol_update() result(t)
      !! New area (absolute, areas are fractions) and new tracer (relative to the cancellation scale of its quotient).
      real(wp) :: t(2)
#ifdef RDB_GPU_OFFLOAD
      t = [1.0e-13_wp, 1.0e-11_wp]
#else
      t = [2.0e-15_wp, 3.0e-12_wp]    ! measured worst 2.2e-16 and 2.2e-13
#endif
   end function tol_update

   subroutine pull_remap_work(s)
      !! Device to host: the packed fields, the geometry, the reconstruction and the fluxes.
      type(sim_t), intent(inout) :: s
      !$acc update self(s%w%a%aim, s%w%a%trm, s%w%b%aim, s%w%b%trm)
      !$acc update self(s%w%gw%src_e, s%w%gw%xp_e, s%w%gw%yp_e, s%w%gw%area_e, s%w%gw%icase_e)
      !$acc update self(s%w%gw%src_n, s%w%gw%xp_n, s%w%gw%yp_n, s%w%gw%area_n, s%w%gw%icase_n)
      !$acc update self(s%w%rb%mx, s%w%rb%my, s%w%rb%tc, s%w%rb%tx, s%w%rb%ty)
      !$acc update self(s%w%fw%mflx_e, s%w%fw%mflx_n, s%w%fw%mtflx_e, s%w%fw%mtflx_n)
   end subroutine pull_remap_work

   subroutine count_geometry(s, ncancel, nout, ncentre24, nedge)
      !! Edges with a pair of opposite-sign triangles in one source cell (cancelling pairs), triangle
      !! vertices lying outside their source cell, edges that took centre branch 2 or 4 (the "rare"
      !! ones), and edges with any triangle.  Reads the pulled geometry.
      type(sim_t), intent(in) :: s
      integer, intent(out) :: ncancel, nout, ncentre24, nedge
      integer :: i, j, s1, s2, v, dir, centre, ic
      logical :: any_tri, canc

      ncancel = 0
      nout = 0
      ncentre24 = 0
      nedge = 0
      do dir = 1, 2
         do j = s%js - 1, s%je
            do i = s%is - 1, s%ie
               any_tri = .false.
               canc = .false.
               do s1 = 1, REMAP_NSLOT
                  if (merge(s%w%gw%src_e(i, j, s1), s%w%gw%src_n(i, j, s1), dir == 1) == 0_int8) cycle
                  any_tri = .true.
                  do v = 1, 3
                     if (dir == 1) then
                        if (abs(s%w%gw%xp_e(i, j, v, s1)) > 0.5_wp + 1.0e-12_wp .or. &
                            abs(s%w%gw%yp_e(i, j, v, s1)) > 0.5_wp + 1.0e-12_wp) nout = nout + 1
                     else
                        if (abs(s%w%gw%xp_n(i, j, v, s1)) > 0.5_wp + 1.0e-12_wp .or. &
                            abs(s%w%gw%yp_n(i, j, v, s1)) > 0.5_wp + 1.0e-12_wp) nout = nout + 1
                     end if
                  end do
                  do s2 = s1 + 1, REMAP_NSLOT
                     if (dir == 1) then
  if (s%w%gw%src_e(i, j, s2) == s%w%gw%src_e(i, j, s1) .and. s%w%gw%area_e(i, j, s1)*s%w%gw%area_e(i, j, s2) < 0.0_wp) canc = .true.
                     else
  if (s%w%gw%src_n(i, j, s2) == s%w%gw%src_n(i, j, s1) .and. s%w%gw%area_n(i, j, s1)*s%w%gw%area_n(i, j, s2) < 0.0_wp) canc = .true.
                     end if
                  end do
               end do
               if (any_tri) nedge = nedge + 1
               if (canc) ncancel = ncancel + 1
               ic = merge(s%w%gw%icase_e(i, j), s%w%gw%icase_n(i, j), dir == 1)
               centre = mod(ic, 325)/25
               if (centre == 2 .or. centre == 4) ncentre24 = ncentre24 + 1
            end do
         end do
      end do
   end subroutine count_geometry

   subroutine test_oracle(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call oracle_body(error, s)
      call sim_destroy(s)
   end subroutine test_oracle

   subroutine oracle_body(error, s)
      !! One step of the whole remap on a sharp noisy state and a random shearing flow, then: every edge
      !! flux against the exact host integration, the new state against the host update, the output ring
      !! untouched.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      type(transport_result_t) :: res
      integer :: pass, i, j, n, nt, ncat, ntr, nx, ny, ncancel, nout, nc24, nedge, nmis, ntot_c, ntot_o, ntot_24
      real(wp) :: worst_f, worst_s, worst_u, fm, e, sc0, fe, tol
      real(wp), allocatable :: ft(:), sc(:)
      real(wp), allocatable :: fem(:, :, :), fnm(:, :, :), fet(:, :, :, :), fnt(:, :, :, :)
      real(wp), allocatable :: gm(:, :, :), gt(:, :, :, :)
      real(wp) :: dm, wt, rt, mm0, mnew, tm0, tp0, tpn, tnew, den, scl, dt
      integer :: p
      real(wp) :: tu(2)

      dt = 1000.0_wp
      tol = tol_flux()
      tu = tol_update()
      worst_f = 0.0_wp
      worst_s = 0.0_wp
      worst_u = 0.0_wp
      ntot_c = 0
      ntot_o = 0
      ntot_24 = 0
      nmis = 0
      do pass = 1, 3
         call sim_init(s, 24, 20, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
         if (pass == 1) then
            call fill_front(s, 11_int64)
            call set_velocity(s, FLOW_RANDOM, cps(s, dt, 0.40_wp), cps(s, dt, 0.10_wp))
         else if (pass == 2) then
            call fill_front(s, 77_int64)
            call set_velocity(s, FLOW_RANDOM, cps(s, dt, 0.20_wp), cps(s, dt, 0.12_wp))
         else
            call fill_smooth(s, 0.3_wp, .true.)
            call set_velocity(s, FLOW_CELLULAR, cps(s, dt, 0.45_wp), 0.0_wp)
         end if
         call sim_step(s, dt, res)
         call chk_ierr(error, res, "oracle step")
         if (allocated(error)) return
         call check(error, s%w%fw%nbatch == s%ncat + 1, "one batch holds every category")
         if (allocated(error)) return
         call pull_remap_work(s)
         nx = s%nx
         ny = s%ny
         ncat = s%ncat
         ntr = s%w%ntrace
         allocate (ft(ntr), sc(0:ntr))
         allocate (fem(nx, ny, 0:ncat), fnm(nx, ny, 0:ncat), fet(nx, ny, ntr, ncat), fnt(nx, ny, ntr, ncat))
         allocate (gm(nx, ny, 0:ncat), gt(nx, ny, ntr, ncat))
         fem = 0.0_wp; fnm = 0.0_wp; fet = 0.0_wp; fnt = 0.0_wp
         ! host fluxes of every edge of the owned window, compared with the stored ones
         do n = 0, ncat
            do j = s%js - 1, s%je
               do i = s%is - 1, s%ie
                  if (j >= s%js) then
                     call host_face_flux(nx, ny, ncat, ntr, i, j, n, s%w%gw%src_e, s%w%gw%xp_e, s%w%gw%yp_e, s%w%gw%area_e, &
                                         s%w%a%aim, s%w%rb%mx, s%w%rb%my, s%w%rb%tc, s%w%rb%tx, s%w%rb%ty, &
                                         s%w%a%tab%tracer_type(1:ntr), s%w%a%tab%depend(1:ntr), fm, ft, sc)
                     fem(i, j, n) = fm
                     e = abs(fm - s%w%fw%mflx_e(i, j, n + 1))/(sc(0) + 1.0e-3_wp*s%dx*s%dy*PUNY)
                     if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
                     worst_f = max(worst_f, e)
                     do nt = 1, ntr
                        if (n == 0) exit
                        fet(i, j, nt, n) = ft(nt)
                        e = abs(ft(nt) - s%w%fw%mtflx_e(i, j, nt, n + 1))/(sc(nt) + 1.0e-3_wp*s%dx*s%dy*PUNY)
                        if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
                        worst_f = max(worst_f, e)
                     end do
                  end if
                  if (i >= s%is) then
                     call host_face_flux(nx, ny, ncat, ntr, i, j, n, s%w%gw%src_n, s%w%gw%xp_n, s%w%gw%yp_n, s%w%gw%area_n, &
                                         s%w%a%aim, s%w%rb%mx, s%w%rb%my, s%w%rb%tc, s%w%rb%tx, s%w%rb%ty, &
                                         s%w%a%tab%tracer_type(1:ntr), s%w%a%tab%depend(1:ntr), fm, ft, sc)
                     fnm(i, j, n) = fm
                     e = abs(fm - s%w%fw%mflx_n(i, j, n + 1))/(sc(0) + 1.0e-3_wp*s%dx*s%dy*PUNY)
                     if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
                     worst_f = max(worst_f, e)
                     do nt = 1, ntr
                        if (n == 0) exit
                        fnt(i, j, nt, n) = ft(nt)
                        e = abs(ft(nt) - s%w%fw%mtflx_n(i, j, nt, n + 1))/(sc(nt) + 1.0e-3_wp*s%dx*s%dy*PUNY)
                        if (.not. ieee_is_finite(e)) e = huge(1.0_wp)
                        worst_f = max(worst_f, e)
                     end do
                  end if
               end do
            end do
         end do
         ! host update (the published rule) from the host fluxes
         gm = 0.0_wp
         gt = 0.0_wp
         do n = 0, ncat
            do j = s%js, s%je
               do i = s%is, s%ie
                  rt = 1.0_wp/s%tarea(i, j)
                  mm0 = s%w%a%aim(i, j, n)
                  dm = ((fem(i, j, n) - fem(i - 1, j, n)) + fnm(i, j, n)) - fnm(i, j - 1, n)
                  mnew = mm0 - dm*rt
                  if (mnew < 0.0_wp .and. mnew >= -PUNY) mnew = 0.0_wp
                  gm(i, j, n) = mnew
                  e = abs(mnew - s%w%b%aim(i, j, n))
                  worst_s = max(worst_s, e)
                  if (.not. ieee_is_finite(s%w%b%aim(i, j, n))) worst_s = huge(1.0_wp)
                  if (n >= 1) then
                     do nt = 1, ntr
                        tm0 = s%w%a%trm(i, j, nt, n)
                        wt = ((fet(i, j, nt, n) - fet(i - 1, j, nt, n)) + fnt(i, j, nt, n)) - fnt(i, j - 1, nt, n)
                        tnew = 0.0_wp
                        den = 1.0_wp
                        scl = abs(mm0*tm0) + abs(wt*rt)
                        if (mnew > 0.0_wp) then
                           if (s%w%a%tab%tracer_type(nt) == 1) then
                              tnew = (mm0*tm0 - wt*rt)/mnew
                              den = mnew
                           else
                              p = s%w%a%tab%depend(nt)
                              tp0 = s%w%a%trm(i, j, p, n)
                              tpn = gt(i, j, p, n)
                              if (abs(tpn) > 0.0_wp) then
                                 tnew = ((mm0*tp0)*tm0 - wt*rt)/(mnew*tpn)
                                 den = mnew*abs(tpn)
                                 scl = abs((mm0*tp0)*tm0) + abs(wt*rt)
                              end if
                           end if
                        end if
                        gt(i, j, nt, n) = tnew
                        e = abs(tnew - s%w%b%trm(i, j, nt, n))/(scl/abs(den) + 1.0e-30_wp)
                        if (.not. ieee_is_finite(s%w%b%trm(i, j, nt, n))) e = huge(1.0_wp)
                        worst_u = max(worst_u, e)
                     end do
                  end if
                  if (((mnew > 0.0_wp) .neqv. (s%w%b%aim(i, j, n) > 0.0_wp)) .and. abs(mnew) > 1.0e-9_wp) nmis = nmis + 1
               end do
            end do
         end do
         ! the ring of the output is untouched (the workspace starts at zero)
         do n = 0, ncat
            do j = 1, ny
               do i = 1, nx
                  if (i >= s%is .and. i <= s%ie .and. j >= s%js .and. j <= s%je) cycle
                  if (s%w%b%aim(i, j, n) /= 0.0_wp) nmis = nmis + 1000
                  if (n >= 1) then
                     do nt = 1, ntr
                        if (s%w%b%trm(i, j, nt, n) /= 0.0_wp) nmis = nmis + 1000
                     end do
                  end if
               end do
            end do
         end do
         call count_geometry(s, ncancel, nout, nc24, nedge)
         ntot_c = ntot_c + ncancel
         ntot_o = ntot_o + nout
         ntot_24 = ntot_24 + nc24
         write (error_unit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "  [remap] oracle case ", pass, ": edges with triangles ", nedge, &
            ", cancelling pairs ", ncancel, ", vertices outside their cell ", nout, ", centre branch 2/4 edges ", nc24
         deallocate (ft, sc, fem, fnm, fet, fnt, gm, gt)
         if (pass < 3) call sim_destroy(s)
      end do
     write (error_unit, "(a,es10.3,a,es10.3,a,es10.3,a,i0)") "  [remap] oracle: worst flux error / absolute-term scale ", worst_f, &
         ", worst area error ", worst_s, ", worst tracer error / scale ", worst_u, ", decision / ring mismatches ", nmis
      if (.not. (worst_f <= tol)) call test_failed(error, "an edge flux differs from the exact integral")
      if (.not. allocated(error)) then
         if (.not. (worst_s <= tu(1))) call test_failed(error, "new area differs from the host update")
      end if
      if (.not. allocated(error)) then
         if (.not. (worst_u <= tu(2))) call test_failed(error, "new tracer differs from the host update")
      end if
      if (.not. allocated(error)) call check(error, nmis, 0, "decisions and the untouched ring agree")
      if (.not. allocated(error)) call check(error, ntot_c > 0, "the flow ran cancelling-pair triangles")
      if (.not. allocated(error)) call check(error, ntot_24 > 0, "the flow ran the rare centre branches 2 and 4")
   end subroutine oracle_body

   ! ------------------------------------------------------------------
   ! (2) conservation
   ! ------------------------------------------------------------------

   function tol_cons() result(t)
      !! Relative change of a conserved total over a run (round-off of one division by the cell area per
      !! step plus the quadrature rounding).
      real(wp) :: t
#ifdef RDB_GPU_OFFLOAD
      t = 1.0e-13_wp          ! device build: FMA contraction in the quadrature
#else
      t = 1.0e-14_wp          ! measured worst 7.7e-17 over 80-100 steps
#endif
   end function tol_cons

   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call conservation_body(error, s)
      call sim_destroy(s)
   end subroutine test_conservation

   subroutine conservation_body(error, s)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: kind, it, nsteps
      real(wp) :: worst, dt, worst7, cum(6), d(5), b(5), acc_err
      integer :: nclip_tot, k
      character(len=14) :: nm(3) = [character(len=14) :: "translation", "rotation", "deformation"]

      dt = 1000.0_wp
      do kind = 1, 3
         call sim_init(s, 40, 32, 5, 7, .true., .true., .true., .true., 1.0e4_wp)
         select case (kind)
         case (1)
            call fill_smooth(s, 0.2_wp, .true.)
            call set_velocity(s, FLOW_UNIFORM, cps(s, dt, 0.37_wp), cps(s, dt, 0.21_wp))
            nsteps = 100
         case (2)
            call fill_front(s, 5_int64)
            call set_velocity(s, FLOW_ROTATION, 0.03_wp/dt, 0.0_wp)
            nsteps = 80
         case default
            call fill_front(s, 9_int64)
            call set_velocity(s, FLOW_CELLULAR, cps(s, dt, 0.45_wp), 0.0_wp)
            nsteps = 80
         end select
         call sim_wrap(s)
         call totals(s, q0)
         worst = 0.0_wp
         worst7 = 0.0_wp
         cum = 0.0_wp
         nclip_tot = 0
         do it = 1, nsteps
            call sim_step(s, dt, res, wrap=.true.)
            call chk_ierr(error, res, "conservation step (driver conservation + monotonicity checks on)")
            if (allocated(error)) return
            nclip_tot = nclip_tot + res%ledger%n_clip
            cum = cum + [res%ledger%clip_area, res%ledger%drop_vice, res%ledger%drop_vsno, res%ledger%drop_eice, &
                         res%ledger%drop_esno - SHIFT*res%ledger%drop_vsno, 0.0_wp]
            call totals(s, q1)
            worst = max(worst, worst_change(q0, q1, 1, 6))
            worst7 = max(worst7, worst_change(q0, q1, 7, 7))
         end do
         write (error_unit, "(a,a,a,i0,a,es10.3,a,es10.3,a,es10.3)") "  [remap] conservation ", trim(nm(kind)), ", steps ", &
            nsteps, ": worst relative change of area, ice and snow volume, ice and snow energy, salt ", worst, &
            "; surface-temperature content ", worst7, "; driver worst_rel (last step) ", res%worst_rel
         if (.not. (worst <= tol_cons())) then
            call test_failed(error, "a conserved total moved beyond round-off")
            return
         end if
         ! the ledger of the guards closes the books: change of the total = + area created by the clips,
         ! - content dropped where the new area (or the parent) is 0, with NO absolute floor
         do k = 1, 5
            b(k) = efp_to_real(q0(k))
            d(k) = efp_real_diff(q1(k), q0(k))
         end do
         ! q1 area, q2 ice volume, q3 snow volume, q4 ice enthalpy, q5 snow enthalpy (unshifted)
         d(1) = d(1) - cum(1)
         d(2) = d(2) + cum(2)
         d(3) = d(3) + cum(3)
         d(4) = d(4) + cum(4)
         d(5) = d(5) + cum(5)
         acc_err = 0.0_wp
         do k = 1, 5
            acc_err = max(acc_err, abs(d(k))/max(abs(b(k)), 1.0e-300_wp))
         end do
         write (error_unit, "(a,i0,a,5es10.2,a,es10.2)") "  [remap]   guards: clipped cells ", nclip_tot, ", created/dropped totals (area, vice, vsno, eice, esno) ", &
            cum(1:5), "; books closed to ", acc_err
         if (.not. (acc_err <= tol_cons())) then
            call test_failed(error, "the totals do not close with the clip/drop ledger")
            return
         end if
         ! the area-tracer content of a cell whose area falls to puny or below is dropped by design
         if (.not. (worst7 <= TOL_T7)) then
            call test_failed(error, "the surface-temperature content moved more than the sliver bound")
            return
         end if
         call sim_destroy(s)
      end do
   end subroutine conservation_body

   subroutine totals(s, q)
      !! Independent EFP totals of the OWNED cells (area-weighted): 1 area (ice + open water),
      !! 2 ice volume, 3 snow volume, 4 ice enthalpy, 5 snow enthalpy, 6 salt, 7 `aicen*Tsfc`.
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

   function worst_change(q0, q1, n0, n1) result(w)
      !! Largest relative change among the totals `n0..n1` (NaN-safe: a non-finite total is `huge`).
      type(efp_t), intent(in) :: q0(7), q1(7)
      integer, intent(in), optional :: n0, n1
      real(wp) :: w
      integer :: n, lo, hi
      real(wp) :: b, d, r

      w = 0.0_wp
      lo = 1
      hi = 7
      if (present(n0)) lo = n0
      if (present(n1)) hi = n1
      do n = lo, hi
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

   pure function rdf(a, b, fl) result(r)
      !! `|a - b|/max(|a|, |b|, fl)`, `huge` when either value is not finite (a floor for a field that is zero in places).
      real(wp), intent(in) :: a, b, fl
      real(wp) :: r
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) then
         r = huge(1.0_wp)
      else
         r = abs(a - b)/max(abs(a), abs(b), fl)
      end if
   end function rdf

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
      !! Largest relative difference over every state field of two sims (ring included unless `only_owned`).
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
   ! (3) monotonicity
   ! ------------------------------------------------------------------

   function tracer_value(s, n, t, i, j, defined) result(v)
      !! The six checked tracers in their physical frame: 1 `hi`, 2 `hs`, 3 `Tsfc` (tracers of area),
      !! 4 `qice(1)`, 5 `sice(3)` (of the ice volume), 6 snow enthalpy (of the snow volume);
      !! `defined` is where the tracer has a meaning (the area, or the parent, above `puny`).
      type(sim_t), intent(in) :: s
      integer, intent(in) :: n, t, i, j
      logical, intent(out) :: defined
      real(wp) :: v
      real(wp) :: a, hi, hs

      a = s%aicen(i, j, n)
      v = 0.0_wp
      defined = .false.
      if (.not. (a > PUNY)) return
      hi = s%vicen(i, j, n)/a
      hs = s%vsnon(i, j, n)/a
      select case (t)
      case (1)
         v = hi
         defined = .true.
      case (2)
         v = hs
         defined = .true.
      case (3)
         v = s%Tsfcn(i, j, n)
         defined = .true.
      case (4)
         v = s%qicen(i, j, n, 1)
         defined = abs(hi) > PUNY
      case (5)
         v = s%sicen(i, j, n, min(3, s%nil))
         defined = abs(hi) > PUNY
      case default
         v = s%qsnon(i, j, n, 1)
         defined = abs(hs) > PUNY
      end select
   end function tracer_value

   function packed_home(s, n, t, i, j) result(v)
      !! The value of tracer `t` the pack hands over for a cell that is NOT `defined` (area at or below
      !! `puny`): the thicknesses and the volume tracers of a sliver are carried, the surface temperature
      !! is 0.  A sliver's area stays in its cell, so this value is part of the mixture that forms the
      !! new tracer.
      type(sim_t), intent(in) :: s
      integer, intent(in) :: n, t, i, j
      real(wp) :: v
      real(wp) :: a

      a = s%aicen(i, j, n)
      v = 0.0_wp
      if (a > 1.0e-100_wp) then
         select case (t)
         case (1)
            v = s%vicen(i, j, n)/a
         case (2)
            v = s%vsnon(i, j, n)/a
         case (4)
            v = s%qicen(i, j, n, 1)
         case (5)
            v = s%sicen(i, j, n, min(3, s%nil))
         case (6)
            v = s%qsnon(i, j, n, 1)
         case default
            v = 0.0_wp                 ! the surface temperature of a sliver is not carried
         end select
      end if
   end function packed_home

   subroutine mono_check(old, new, nviol, worst)
      !! Independent monotonicity check: each tracer of the NEW state, where defined, lies within the
      !! min/max of the tracer over the 5x5 neighbourhood of the OLD state (cells where it is defined),
      !! up to `max(1, |bound|) puny`.  `worst` is the largest excursion beyond the bound relative to
      !! `max(1, |bound|)`.
      type(sim_t), intent(in) :: old, new
      integer, intent(out) :: nviol
      real(wp), intent(out) :: worst
      integer :: i, j, n, t, di, dj
      logical :: dnew, dold, any_def
      real(wp) :: vnew, vold, lo, hi, w, ex

      nviol = 0
      worst = 0.0_wp
      do n = 1, old%ncat
         do t = 1, 6
            do j = old%js, old%je
               do i = old%is, old%ie
                  vnew = tracer_value(new, n, t, i, j, dnew)
                  if (.not. dnew) cycle
                  lo = huge(1.0_wp)
                  hi = -huge(1.0_wp)
                  any_def = .false.
                  ! the cell's own old content stays (a sliver is not a defined value but its area is kept)
                  vold = tracer_value(old, n, t, i, j, dold)
                  if (.not. dold) vold = packed_home(old, n, t, i, j)
                  if (old%aicen(i, j, n) > 1.0e-100_wp) then
                     lo = vold
                     hi = vold
                     any_def = .true.
                  end if
                  do dj = -2, 2
                     do di = -2, 2
                        vold = tracer_value(old, n, t, i + di, j + dj, dold)
                        if (.not. dold) cycle
                        lo = min(lo, vold)
                        hi = max(hi, vold)
                        any_def = .true.
                     end do
                  end do
                  if (.not. ieee_is_finite(vnew)) then
                     nviol = nviol + 1
                     worst = huge(1.0_wp)
                     cycle
                  end if
                  if (.not. any_def) then
                     ! a tracer appeared with no defined neighbour: an unbounded extremum
                     nviol = nviol + 1
                     worst = huge(1.0_wp)
                     cycle
                  end if
                  w = max(1.0_wp, abs(lo), abs(hi))
                  ex = max(lo - vnew, vnew - hi, 0.0_wp)/w
                  worst = max(worst, ex)
                  if (ex > PUNY) nviol = nviol + 1
               end do
            end do
         end do
      end do
   end subroutine mono_check

   subroutine test_monotone(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, old
      call monotone_body(error, s, old)
      call sim_destroy(s)
      call sim_destroy(old)
   end subroutine test_monotone

   subroutine monotone_body(error, s, old)
      !! Sharp noisy fronts under a random shearing flow with a large noise part: no new extremum of any
      !! tracer beyond the 5x5 bound (independent check and the driver's own), and the geometry's
      !! cancelling-pair triangles and vertices outside their cell ran.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, old
      type(transport_result_t) :: res
      integer :: it, nv, nvt, ncancel, nout, nc24, nedge, tc, to, t24, pass
      real(wp) :: worst, w, dt

      dt = 1000.0_wp
      nvt = 0
      worst = 0.0_wp
      tc = 0
      to = 0
      t24 = 0
      do pass = 1, 2
         call sim_init(s, 32, 26, 5, 7, .true., .true., .true., .true., 1.0e4_wp)
         call sim_init(old, 32, 26, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
         call fill_front(s, 31_int64 + 100_int64*pass)
         if (pass == 1) then
            call set_velocity(s, FLOW_RANDOM, cps(s, dt, 0.30_wp), cps(s, dt, 0.10_wp))
         else
            call set_velocity(s, FLOW_RANDOM, cps(s, dt, 0.10_wp), cps(s, dt, 0.15_wp))
         end if
         do it = 1, 30
            call sim_wrap(s)
            call copy_state(s, old)
            call sim_step(s, dt, res, wrap=.false.)
            call chk_ierr(error, res, "monotone step (driver monotonicity check on)")
            if (allocated(error)) return
            call mono_check(old, s, nv, w)
            nvt = nvt + nv
            worst = max(worst, w)
            call pull_remap_work(s)
            call count_geometry(s, ncancel, nout, nc24, nedge)
            tc = tc + ncancel
            to = to + nout
            t24 = t24 + nc24
         end do
         call sim_destroy(s)
         call sim_destroy(old)
      end do
     write (error_unit, "(a,i0,a,es10.3,a,i0,a,i0,a,i0)") "  [remap] monotonicity: violations of the independent 5x5 bound ", nvt, &
         ", worst excursion ", worst, "; edge-steps with cancelling pairs ", tc, ", vertices outside their cell ", to, &
         ", centre branch 2/4 ", t24
      call check(error, nvt, 0, "no tracer left the 5x5 bound")
      if (.not. allocated(error)) call check(error, tc > 0, "cancelling-pair triangles were exercised")
      if (.not. allocated(error)) call check(error, t24 > 0, "the rare centre branches 2 and 4 were exercised")
   end subroutine monotone_body

   ! ------------------------------------------------------------------
   ! (4) one-cell translation
   ! ------------------------------------------------------------------

   subroutine test_one_cell(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, old
      call one_cell_body(error, s, old)
      call sim_destroy(s)
      call sim_destroy(old)
   end subroutine test_one_cell

   subroutine one_cell_body(error, s, old)
      !! Translation by exactly one cell per step reproduces the shifted field to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, old
      type(transport_result_t) :: res
      integer :: c, i, j, n, k, di, dj, it
      real(wp) :: worst, w, dt, cu(3), cv(3), one

      dt = 1000.0_wp
      one = 1.0_wp - 2.0_wp**(-47)
      cu = [one, 0.0_wp, -one]
      cv = [0.0_wp, one, one]
      worst = 0.0_wp
      do c = 1, 3
         call sim_init(s, 28, 20, 5, 7, .true., .true., .true., .true., 1.0e4_wp)
         call sim_init(old, 28, 20, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
         call fill_smooth(s, 0.4_wp*real(c, wp), .true.)
         call set_velocity(s, FLOW_UNIFORM, cu(c)*s%dx/dt, cv(c)*s%dy/dt)
         di = nint(cu(c))
         dj = nint(cv(c))
         ! (the geometry refuses a displacement longer than the corner length; round-off in the midpoint
         ! interpolation can push exactly 1 cell over, so the step moves 1 - 2^-47 cells)
         do it = 1, 3
            call sim_wrap(s)
            call copy_state(s, old)
            call sim_step(s, dt, res, wrap=.false.)
            call chk_ierr(error, res, "one-cell step")
            if (allocated(error)) return
            w = 0.0_wp
            do j = s%js, s%je
               do i = s%is, s%ie
                  w = max(w, rd(s%aice0(i, j), old%aice0(i - di, j - dj)))
                  do n = 1, s%ncat
                     w = max(w, rd(s%aicen(i, j, n), old%aicen(i - di, j - dj, n)), &
                             rd(s%vicen(i, j, n), old%vicen(i - di, j - dj, n)), &
                             rdf(s%vsnon(i, j, n), old%vsnon(i - di, j - dj, n), 1.0e-3_wp), &
                             rd(s%Tsfcn(i, j, n), old%Tsfcn(i - di, j - dj, n)))
                     if (old%vsnon(i - di, j - dj, n) > 1.0e-6_wp) &      ! the enthalpy of a negligible snow volume is not a value
                        w = max(w, rd(s%qsnon(i, j, n, 1), old%qsnon(i - di, j - dj, n, 1)))
                     do k = 1, s%nil
                        w = max(w, rd(s%qicen(i, j, n, k), old%qicen(i - di, j - dj, n, k)), &
                                rd(s%sicen(i, j, n, k), old%sicen(i - di, j - dj, n, k)))
                     end do
                  end do
               end do
            end do
            worst = max(worst, w)
         end do
         call sim_destroy(s)
         call sim_destroy(old)
      end do
      write (error_unit, "(a,es10.3)") "  [remap] one-cell translation (x, y, diagonal; 3 steps each, 1 - 2^-47 cells): worst relative difference from the shifted field ", worst
      if (.not. (worst <= 1.0e-12_wp)) call test_failed(error, "one-cell translation does not reproduce the shifted field")
   end subroutine one_cell_body

   ! ------------------------------------------------------------------
   ! (5) half-cell translation: order of accuracy
   ! ------------------------------------------------------------------

   subroutine test_half_cell(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call half_cell_body(error, s)
      call sim_destroy(s)
   end subroutine test_half_cell

   subroutine half_cell_error(s, nxp, upwind, e2)
      !! L2 error of the category-2 area of a smooth sine moved half a cell per step over half the domain.
      type(sim_t), intent(inout) :: s
      integer, intent(in) :: nxp
      logical, intent(in) :: upwind
      real(wp), intent(out) :: e2
      type(transport_result_t) :: res
      integer :: i, j, it, nn
      real(wp) :: lx, x, dt, ex

      lx = 3.2e5_wp
      call sim_init(s, nxp, 4, 5, 7, .true., .true., .false., .false., lx/real(nxp, wp))
      if (upwind) s%p%transport%advection = ICE_ADVECT_UPWIND
      dt = 1000.0_wp
      do j = 1, s%ny
         do i = 1, s%nx
            x = (real(i - s%is, wp) + 0.5_wp)*s%dx
            s%aicen(i, j, 2) = 0.3_wp + 0.2_wp*sin(TWOPI*x/lx)
            s%vicen(i, j, 2) = s%aicen(i, j, 2)*HMID(2)
            s%Tsfcn(i, j, 2) = -5.0_wp
            s%qsnon(i, j, 2, 1) = -SHIFT
            s%qicen(i, j, 2, :) = -2.0e8_wp
            s%sicen(i, j, 2, :) = 5.0_wp
            s%aice0(i, j) = 1.0_wp - s%aicen(i, j, 2)
         end do
      end do
      call set_velocity(s, FLOW_UNIFORM, 0.5_wp*s%dx/dt, 0.0_wp)
      nn = nxp
      do it = 1, nn
         call sim_step(s, dt, res)
         if (res%ierr /= 0 .or. res%cfg /= 0) then
            e2 = huge(1.0_wp)
            return
         end if
      end do
      e2 = 0.0_wp
      do j = s%js, s%je
         do i = s%is, s%ie
            x = (real(i - s%is, wp) + 0.5_wp)*s%dx
            ex = 0.3_wp + 0.2_wp*sin(TWOPI*(x - 0.5_wp*real(nn, wp)*s%dx)/lx)
            e2 = e2 + (s%aicen(i, j, 2) - ex)**2
         end do
      end do
      e2 = sqrt(e2/real(s%nxp*s%nyp, wp))
   end subroutine half_cell_error

   subroutine half_cell_body(error, s)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      integer :: k, nxp(3)
      real(wp) :: er(3), eu(3), r1, r2, u1, u2

      nxp = [32, 64, 128]
      do k = 1, 3
         call half_cell_error(s, nxp(k), .false., er(k))
         call sim_destroy(s)
         call half_cell_error(s, nxp(k), .true., eu(k))
         call sim_destroy(s)
      end do
      r1 = er(1)/er(2)
      r2 = er(2)/er(3)
      u1 = eu(1)/eu(2)
      u2 = eu(2)/eu(3)
      write (error_unit, "(a,3es11.3,a,2f7.3)") "  [remap] half-cell translation, L2 error at 32/64/128 cells (remap): ", er, &
         "; error ratios ", r1, r2
      write (error_unit, "(a,3es11.3,a,2f7.3)") "  [remap] half-cell translation, L2 error at 32/64/128 cells (upwind): ", eu, &
         "; error ratios ", u1, u2
      if (.not. (ieee_is_finite(r1) .and. ieee_is_finite(r2) .and. ieee_is_finite(u1) .and. ieee_is_finite(u2))) then
         call test_failed(error, "non-finite error")
         return
      end if
      if (.not. (r1 > 3.3_wp .and. r2 > 3.3_wp)) call test_failed(error, "remap is not second order on a half-cell translation")
      if (.not. allocated(error)) then
         if (.not. (u1 > 1.6_wp .and. u1 < 2.4_wp .and. u2 > 1.6_wp .and. u2 < 2.4_wp)) &
            call test_failed(error, "upwind is not first order on the same case")
      end if
      if (.not. allocated(error)) then
         if (.not. (er(3) < 0.05_wp*eu(3))) call test_failed(error, "remap is not more accurate than upwind")
      end if
   end subroutine half_cell_body

   ! ------------------------------------------------------------------
   ! (6) rotation of a cone
   ! ------------------------------------------------------------------

   subroutine test_rotation(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call rotation_body(error, s)
      call sim_destroy(s)
   end subroutine test_rotation

   subroutine cone_run(s, upwind, retention, mass_err, error)
      type(sim_t), intent(inout) :: s
      logical, intent(in) :: upwind
      real(wp), intent(out) :: retention, mass_err
      type(error_type), allocatable, intent(inout) :: error
      type(transport_result_t) :: res
      type(efp_t) :: q0(7), q1(7)
      integer :: it, nsteps
      real(wp) :: dt, peak0

      nsteps = 288
      dt = 1000.0_wp
      call sim_init(s, 64, 64, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      if (upwind) s%p%transport%advection = ICE_ADVECT_UPWIND
      call fill_cone(s, 2, 46.0_wp, 32.0_wp, 8.0_wp, 0.8_wp)
      call set_velocity(s, FLOW_ROTATION, TWOPI/real(nsteps, wp)/dt, 0.0_wp)
      peak0 = maxval(s%aicen(s%is:s%ie, s%js:s%je, 2))
      call sim_wrap(s)
      call totals(s, q0)
      retention = -1.0_wp
      mass_err = huge(1.0_wp)
      do it = 1, nsteps
         call sim_step(s, dt, res)
         call chk_ierr(error, res, "cone step")
         if (allocated(error)) return
      end do
      call totals(s, q1)
      mass_err = worst_change(q0, q1, 1, 6)
      retention = maxval(s%aicen(s%is:s%ie, s%js:s%je, 2))/peak0
   end subroutine cone_run

   subroutine rotation_body(error, s)
      !! One revolution (288 steps, at most 0.48 cell per step) of a cone of radius 8 cells in solid-body
      !! rotation: the peak kept by remap against upwind.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      real(wp) :: rr, ru, mr, mu

      call cone_run(s, .false., rr, mr, error)
      if (allocated(error)) return
      call sim_destroy(s)
      call cone_run(s, .true., ru, mu, error)
      if (allocated(error)) return
      write (error_unit, "(a,f7.4,a,f7.4,a,es9.2,a,es9.2)") "  [remap] cone, one revolution: peak retained, remap ", rr, ", upwind ", ru, &
         "; area change remap ", mr, ", upwind ", mu
      if (.not. (rr > 0.70_wp)) call test_failed(error, "remap lost the peak of the cone")
      if (.not. allocated(error)) then
         if (.not. (rr > 2.0_wp*ru)) call test_failed(error, "remap does not retain the peak better than upwind")
      end if
      if (.not. allocated(error)) then
         if (.not. (mr <= tol_cons())) call test_failed(error, "the area of the cone is not conserved")
      end if
   end subroutine rotation_body

   ! ------------------------------------------------------------------
   ! (7) divergence consistency
   ! ------------------------------------------------------------------

   subroutine test_divergence(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call divergence_body(error, s)
      call sim_destroy(s)
   end subroutine test_divergence

   subroutine divergence_body(error, s)
      !! A uniform field in a convergent / divergent flow in a closed basin changes by exactly the flux
      !! convergence of the face velocities, and a constant tracer stays constant to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      type(transport_result_t) :: res
      integer :: i, j, k, it, n
      real(wp) :: dt, c, ex, e1, e0, tc, worst_tr, h0, v0
      real(wp) :: hvals(3)

      dt = 1000.0_wp
      call sim_init(s, 32, 24, 5, 7, .false., .false., .true., .true., 1.0e4_wp)
      s%aicen = 0.0_wp
      s%vicen = 0.0_wp
      s%vsnon = 0.0_wp
      do j = s%js, s%je
         do i = s%is, s%ie
            s%aicen(i, j, 2) = 0.4_wp
            s%vicen(i, j, 2) = 0.4_wp*1.3_wp
            s%vsnon(i, j, 2) = 0.4_wp*0.15_wp
            s%Tsfcn(i, j, 2) = -7.0_wp
            s%qsnon(i, j, 2, 1) = -SHIFT*1.2_wp
            do k = 1, s%nil
               s%qicen(i, j, 2, k) = -2.0e8_wp - 1.0e6_wp*real(k, wp)
               s%sicen(i, j, 2, k) = 4.0_wp + 0.5_wp*real(k, wp)
            end do
            s%aice0(i, j) = 0.6_wp
         end do
      end do
      s%aice0(1:s%is - 1, :) = 0.0_wp
      s%aice0(s%ie + 1:s%nx, :) = 0.0_wp
      s%aice0(:, 1:s%js - 1) = 0.0_wp
      s%aice0(:, s%je + 1:s%ny) = 0.0_wp
      call set_velocity(s, FLOW_CONVERGENT, cps(s, dt, 0.30_wp), cps(s, dt, 0.25_wp))
      ! one step: the exact new area from the face velocities
      call sim_step(s, dt, res, wrap=.false.)
      call chk_ierr(error, res, "divergence step")
      if (allocated(error)) return
      e1 = 0.0_wp
      e0 = 0.0_wp
      do j = s%js, s%je
         do i = s%is, s%ie
            c = dt*((s%u(i + 1, j) - s%u(i, j))*s%HTE(i, j) + (s%v(i, j + 1) - s%v(i, j))*s%HTN(i, j))/s%tarea(i, j)
            ex = 0.4_wp*(1.0_wp - c)
            e1 = max(e1, abs(s%aicen(i, j, 2) - ex))
            ex = 0.6_wp*(1.0_wp - c)
            e0 = max(e0, abs(s%aice0(i, j) - ex))
         end do
      end do
      write (error_unit, "(a,es10.3,a,es10.3)") "  [remap] convergent flow, one step: area error vs the flux convergence, category ", e1, &
         ", open water ", e0
      if (.not. (e1 <= 1.0e-14_wp .and. e0 <= 1.0e-14_wp)) call test_failed(error, "area did not change by the flux convergence")
      if (allocated(error)) return
      ! more steps, then the constant tracers
      do it = 2, 12
         call sim_step(s, dt, res, wrap=.false.)
         call chk_ierr(error, res, "divergence step")
         if (allocated(error)) return
      end do
      worst_tr = 0.0_wp
      do j = s%js, s%je
         do i = s%is, s%ie
            n = 2
            hvals = [s%vicen(i, j, n)/s%aicen(i, j, n)/1.3_wp - 1.0_wp, s%vsnon(i, j, n)/s%aicen(i, j, n)/0.15_wp - 1.0_wp, &
                     s%Tsfcn(i, j, n)/(-7.0_wp) - 1.0_wp]
            worst_tr = max(worst_tr, abs(hvals(1)), abs(hvals(2)), abs(hvals(3)), abs(s%qsnon(i, j, n, 1)/(-SHIFT*1.2_wp) - 1.0_wp))
            do k = 1, s%nil
               worst_tr = max(worst_tr, abs(s%qicen(i, j, n, k)/(-2.0e8_wp - 1.0e6_wp*real(k, wp)) - 1.0_wp), &
                              abs(s%sicen(i, j, n, k)/(4.0_wp + 0.5_wp*real(k, wp)) - 1.0_wp))
            end do
         end do
      end do
      tc = maxval(abs(s%aicen(s%is:s%ie, s%js:s%je, 2) - 0.4_wp))
      h0 = 0.0_wp
      v0 = 0.0_wp
      write (error_unit, "(a,es10.3,a,es10.3)") "  [remap] 12 steps of the convergent flow: area ranges over ", tc, &
         "; constant tracers, worst relative departure ", worst_tr
      if (.not. (tc > 0.02_wp)) call test_failed(error, "the flow did not change the area")
      if (.not. allocated(error)) then
         if (.not. (worst_tr <= 2.0e-14_wp)) call test_failed(error, "a constant tracer did not stay constant")
      end if
   end subroutine divergence_body

   ! ------------------------------------------------------------------
   ! (8) closed faces and land
   ! ------------------------------------------------------------------

   subroutine test_land(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call land_body(error, s)
      call sim_destroy(s)
   end subroutine test_land

   subroutine land_body(error, s)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: i, j, n, it, nleak
      real(wp) :: worst, dt, worst7

      dt = 8000.0_wp
      call sim_init(s, 24, 16, 5, 7, .false., .false., .true., .true., 1.0e4_wp)
      call fill_smooth(s, 0.0_wp, .false.)
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
      call set_velocity(s, FLOW_UNIFORM, cps(s, dt, 0.4_wp), cps(s, dt, 0.15_wp))
      call close_faces(s)
      call totals(s, q0)
      worst = 0.0_wp
      worst7 = 0.0_wp
      nleak = 0
      do it = 1, 30
         call sim_step(s, dt, res, wrap=.false.)
         call chk_ierr(error, res, "land step")
         if (allocated(error)) return
         call totals(s, q1)
         worst = max(worst, worst_change(q0, q1, 1, 6))
         worst7 = max(worst7, worst_change(q0, q1, 7, 7))
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
    write (error_unit, "(a,es10.3,a,es10.3,a,i0)") "  [remap] wall + land block, 30 steps pushing into them: worst total change ", &
         worst, " (surface-temperature content ", worst7, "), closed cells touched: ", nleak
      if (.not. (worst <= tol_cons() .and. worst7 <= TOL_T7)) call test_failed(error, "ice leaked through a closed face")
      if (.not. allocated(error)) call check(error, nleak, 0, "closed cells stay exactly empty (and Tsfc 0 on land)")
   end subroutine land_body

   ! ------------------------------------------------------------------
   ! (9) no-ops
   ! ------------------------------------------------------------------

   subroutine test_noops(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      call noops_body(error, s, ref)
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_noops

   subroutine noops_body(error, s, ref)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst, wp_pack
      integer :: i, j, n, nt

      ! (a) empty domain, nonzero velocity: every state array bitwise unchanged
      call sim_init(s, 16, 12, 5, 7, .true., .true., .true., .false., 1.0e4_wp)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .true., .false., 1.0e4_wp)
      s%qicen = -3.0e8_wp                    ! junk in the empty categories must not be touched
      s%sicen = 7.0_wp
      call copy_state(s, ref)
      call set_velocity(s, FLOW_UNIFORM, 4.0_wp, -3.0_wp)
      call sim_step(s, 1000.0_wp, res)
      call chk_ierr(error, res, "empty: no bit")
      if (allocated(error)) return
      call check(error, res%stepped, "empty: stepped")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      write (error_unit, "(a,es10.3)") "  [remap] empty domain: worst difference ", worst
      if (.not. (worst == 0.0_wp)) then
         call test_failed(error, "empty domain must be a bitwise no-op")
         return
      end if
      call sim_destroy(s)
      call sim_destroy(ref)
      ! (b) zero velocity on a full smooth state: packed fields bitwise
      call sim_init(s, 20, 16, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call sim_init(ref, 20, 16, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call fill_smooth(s, 0.5_wp, .true.)
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, 1000.0_wp, res)
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
      write (error_unit, "(a,es10.3,a,es10.3)") "  [remap] zero velocity: packed fields worst difference ", wp_pack, &
         ", state worst difference (pack round trip) ", worst
      if (.not. (wp_pack == 0.0_wp)) call test_failed(error, "zero velocity must leave the packed fields bitwise unchanged")
      if (.not. allocated(error)) then
         if (.not. (worst <= 1.0e-15_wp)) call test_failed(error, "state after the pack round trip beyond 1e-15")
      end if
   end subroutine noops_body

   ! ------------------------------------------------------------------
   ! (10) NaN
   ! ------------------------------------------------------------------

   subroutine test_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      call nan_body(error, s, ref)
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_nan

   subroutine nan_body(error, s, ref)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst, nan
      integer :: pass
      logical :: chk

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      do pass = 1, 5
         chk = (pass >= 4)
         call sim_init(s, 16, 12, 5, 7, .true., .true., chk, .false., 1.0e4_wp)
         call sim_init(ref, 16, 12, 5, 7, .true., .true., chk, .false., 1.0e4_wp)
         call fill_smooth(s, 0.1_wp, .false.)
         call set_velocity(s, FLOW_UNIFORM, 3.0_wp, 1.0_wp)
         if (mod(pass, 3) == 2) s%vicen(8, 6, 2) = nan  ! NaN field value in an owned cell
         call sim_wrap(s)
         select case (mod(pass, 3))
         case (1)
            s%u(9, 7) = nan                              ! NaN face velocity in the window
         case (0)
            s%u(s%is - 1, 7) = nan                       ! NaN face velocity at the outer edge of the geometry's window
         case default
         end select
         call copy_state(s, ref)
         call sim_step(s, 1000.0_wp, res, wrap=.false.)
         if (iand(res%ierr, ICE_ERR_NONFINITE) == 0) then
            call test_failed(error, "NaN did not raise ICE_ERR_NONFINITE")
            return
         end if
         if (res%stepped) then
            call test_failed(error, "a state with NaN was advanced")
            return
         end if
         worst = 0.0_wp
         call state_diff_nan(s, ref, worst)
         if (.not. (worst == 0.0_wp)) then
            call test_failed(error, "state changed by a refused step")
            return
         end if
         if (.not. any(iand(s%ice_err, ICE_ERR_NONFINITE) /= 0)) then
            call test_failed(error, "no cell carries the non-finite bit")
            return
         end if
         write (error_unit, "(a,i0,a,l1,a)") "  [remap] NaN case ", pass, " check=", chk, ": bit raised, state untouched"
         call sim_destroy(s)
         call sim_destroy(ref)
      end do
   end subroutine nan_body

   ! ------------------------------------------------------------------
   ! (11) CFL
   ! ------------------------------------------------------------------

   subroutine test_cfl(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      call cfl_body(error, s, ref)
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_cfl

   subroutine cfl_body(error, s, ref)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst, dt

      dt = 1000.0_wp
      call sim_init(s, 16, 12, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call fill_smooth(s, 0.1_wp, .false.)
      call set_velocity(s, FLOW_UNIFORM, 0.5_wp*s%dx/dt, 0.0_wp)
      s%u(9, 7) = 0.999_wp*s%dx/dt                      ! just inside
      call sim_wrap(s)
      call sim_step(s, dt, res, wrap=.false.)
      call chk_ierr(error, res, "a face just inside the limit is accepted")
      if (allocated(error)) return
      call set_velocity(s, FLOW_UNIFORM, 0.5_wp*s%dx/dt, 0.0_wp)
      s%u(9, 7) = 1.001_wp*s%dx/dt                      ! just outside
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_step(s, dt, res, wrap=.false.)
 write (error_unit, "(a,i0,a,z0)") "  [remap] one face just above the CFL limit: violating cells ", res%n_cfl, ", ierr 0x", res%ierr
      call check(error, res%n_cfl, 1, "one violating cell for one violating face")
      if (allocated(error)) return
      call check(error, iand(res%ierr, ICE_ERR_TRANSPORT_CFL) /= 0, "CFL bit raised")
      if (allocated(error)) return
      call check(error, res%stepped, .false., "a CFL violation is refused")
      if (allocated(error)) return
      call state_diff(s, ref, worst, .false.)
      if (.not. (worst == 0.0_wp)) call test_failed(error, "the refused step left the state alone")
   end subroutine cfl_body

   ! ------------------------------------------------------------------
   ! (12) category batches
   ! ------------------------------------------------------------------

   subroutine test_batch(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      call batch_body(error, s, ref)
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_batch

   subroutine batch_body(error, s, ref)
      !! The same five steps with the flux buffers holding 6, 4, 2 and 1 category indices: every state
      !! array identical bit for bit.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, ref
      type(transport_result_t) :: res
      integer :: k, it, nb(4)
      real(wp) :: worst, dt

      dt = 1000.0_wp
      nb = [6, 4, 2, 1]
      call sim_init(ref, 22, 18, 5, 7, .true., .true., .false., .false., 1.0e4_wp, nbatch=6)
      call fill_front(ref, 61_int64)
      call set_velocity(ref, FLOW_RANDOM, cps(ref, dt, 0.30_wp), cps(ref, dt, 0.10_wp))
      do it = 1, 5
         call sim_step(ref, dt, res)
         call chk_ierr(error, res, "batch reference step")
         if (allocated(error)) return
      end do
      do k = 2, 4
         call sim_init(s, 22, 18, 5, 7, .true., .true., .false., .false., 1.0e4_wp, nbatch=nb(k))
         call fill_front(s, 61_int64)
         call set_velocity(s, FLOW_RANDOM, cps(s, dt, 0.30_wp), cps(s, dt, 0.10_wp))
         do it = 1, 5
            call sim_step(s, dt, res)
            call chk_ierr(error, res, "batch step")
            if (allocated(error)) return
         end do
         call state_diff(s, ref, worst, .false.)
  write (error_unit, "(a,i0,a,i0,a,es10.3,a,i0,a)") "  [remap] flux and reconstruction buffers of ", nb(k), " category indices (", &
            (s%ncat + nb(k))/nb(k), " batches): worst difference from the single batch ", worst, "; flux bytes ", &
            s%w%fw%bytes(), ""
         if (.not. (worst == 0.0_wp)) then
            call test_failed(error, "batching changed the answer")
            return
         end if
         call sim_destroy(s)
      end do
   end subroutine batch_body

   ! ------------------------------------------------------------------
   ! (12b) a departure point a hair off a cell side
   ! ------------------------------------------------------------------

   subroutine test_tiny(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call tiny_body(error, s)
      call sim_destroy(s)
   end subroutine test_tiny

   subroutine tiny_sweep(error)
      !! Every component of the departure points, perturbed by `+-10^k` cells, `k = -20..-12`, around three
      !! base configurations (a corner at zero displacement on either side, both at zero), on both edge
      !! orientations: the flux adjustment closes at every value (no error bit, `|sum - edgearea| <= 1e-13` of
      !! the cell area).  The failing input was a departure point a hair off the cell side, whose axis
      !! crossing rounded ONTO the corner and missed the (then open) pinned-crossing interval.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: base(5, 3), d(4), px(3, 6), py(3, 6), ar(6), asum, fc, ea, delta, worst
      integer(int8) :: sc(6)
      integer(int32) :: ic, bits
      integer :: ib, comp, k, sg, ie, nfail, ncall

      fc = 1.0e8_wp
      base(:, 1) = [0.0_wp, 0.0_wp, -2.00978e-3_wp, 2.79308e-2_wp, -1.51604e-4_wp]
      base(:, 2) = [2.00978e-3_wp, -2.79308e-2_wp, 0.0_wp, 0.0_wp, 1.51604e-4_wp]
      base(:, 3) = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0e-4_wp]
      nfail = 0
      ncall = 0
      worst = 0.0_wp
      do ie = 1, 2
         do ib = 1, 3
            do comp = 1, 4
               do k = 12, 20
                  do sg = -1, 1, 2
                     delta = real(sg, wp)*10.0_wp**(-real(k, wp))
                     d = base(1:4, ib)
                     d(comp) = d(comp) + delta
                     ea = base(5, ib)*fc
                     call remap_locate_edge(ie, 1, PUNY, d(1), d(2), d(3), d(4), fc, fc, fc, fc, fc, fc, fc, ea, px, py, ar, sc, &
                                            ic, asum, bits)
                     ncall = ncall + 1
                     if (bits /= 0_int32) then
                        nfail = nfail + 1
                     else
                        worst = max(worst, abs(asum - ea)/fc)
                        if (abs(asum - ea) > 1.0e-13_wp*fc) nfail = nfail + 1
                     end if
                  end do
               end do
            end do
         end do
      end do
      write (error_unit, "(a,i0,a,i0,a,es10.3)") "  [remap] tiny departure sweep 1e-20..1e-12, both signs, 4 components, 3 bases, 2 edges: ", &
         ncall, " edges, failures ", nfail, ", worst |sum - edgearea| / cell area ", worst
      call check(error, nfail, 0, "the flux adjustment closes at every tiny departure offset")
   end subroutine tiny_sweep

   subroutine tiny_body(error, s)
      !! Regression: the departure point DL of an east edge a displacement of `-1e-20` cells off a cell side
      !! (a corner velocity of 1e-21 m/s, as the tapered edge of a rotation or ice at rest gives) sent the
      !! case table of the departure geometry down a branch that left the flux adjustment `1e-2` of a cell
      !! area short of the prescribed edge area.  The displacements below 1e-16 of a corner length are now zero.
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      integer :: i, j, k, nbad
      real(wp) :: sgn, y, ea, fc, err

      call tiny_sweep(error)
      if (allocated(error)) return
      call sim_init(s, 12, 10, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call sim_map(s)
      call sim_push(s)
      i = s%is + 3
      j = s%js + 3
      fc = s%dx*s%dy
      nbad = 0
      err = 0.0_wp
      do k = 1, 6
         sgn = merge(-1.0_wp, 1.0_wp, mod(k, 2) == 1)
         y = sgn*10.0_wp**(-17.0_wp - real(k, wp))     ! 1e-18 ... 1e-23 cells, both signs
         s%w%gw%dpx = 0.0_wp
         s%w%gw%dpy = 0.0_wp
         s%w%gw%dpx(i, j) = y*s%dx                         ! DL: NE corner, along the edge axis
         s%w%gw%dpy(i, j - 1) = 2.00978e-3_wp*s%dy         ! DR: SE corner
         s%w%gw%dpx(i, j - 1) = 2.79308e-2_wp*s%dx
         s%w%gw%edgearea_e = 0.0_wp
         ea = -1.51604e-4_wp*fc
         s%w%gw%edgearea_e(i, j) = ea
         s%ice_err = 0_int32
         !$acc update device(s%w%gw%dpx, s%w%gw%dpy, s%w%gw%edgearea_e, s%ice_err)
         call remap_locate_east(s%nx, s%ny, s%is, s%ie, s%js, s%je, 1, PUNY, s%w%gw%dpx, s%w%gw%dpy, s%dxU, s%dyU, s%earea, &
                                s%narea, s%w%gw%edgearea_e, s%w%gw%xp_e, s%w%gw%yp_e, s%w%gw%area_e, s%w%gw%src_e, &
                                s%w%gw%icase_e, s%ice_err)
         !$acc update self(s%w%gw%area_e, s%ice_err)
         if (any(s%ice_err /= 0_int32)) nbad = nbad + 1
         err = max(err, abs(sum(s%w%gw%area_e(i, j, :)) - ea)/fc)
      end do
  write (error_unit, "(a,i0,a,es10.3)") "  [remap] tiny departure offsets (1e-18..1e-23 cells, both signs): edges refused ", nbad, &
         ", worst |sum(triarea) - edgearea| / cell area ", err
      call check(error, nbad, 0, "no edge with a tiny departure offset is refused")
      if (.not. allocated(error)) then
         if (.not. (err <= 1.0e-13_wp)) call test_failed(error, "the flux adjustment does not close for a tiny departure offset")
      end if
   end subroutine tiny_body

   ! ------------------------------------------------------------------
   ! (12c) non-uniform metrics
   ! ------------------------------------------------------------------

   subroutine set_metrics_nonuniform(s)
      !! A stretched, converging grid: cell height `H(i)` stretched by +-25 % along x, cell width `W(j)`
      !! shrinking by 45 % toward the north (the convergence of meridians), the face and corner lengths
      !! and the face-cell areas from them.
      type(sim_t), intent(inout) :: s
      integer :: i, j, ie, jn
      real(wp), allocatable :: w(:), h(:)

      allocate (w(s%ny), h(s%nx))
      do j = 1, s%ny
         w(j) = s%dx*(1.0_wp - 0.45_wp*real(j - 1, wp)/real(s%ny - 1, wp))
      end do
      do i = 1, s%nx
         h(i) = s%dy*(1.0_wp + 0.25_wp*sin(TWOPI*real(i, wp)/real(s%nx, wp)))
      end do
      do j = 1, s%ny
         do i = 1, s%nx
            ie = min(i + 1, s%nx)
            jn = min(j + 1, s%ny)
            s%dxT(i, j) = w(j)
            s%dyT(i, j) = h(i)
            s%tarea(i, j) = w(j)*h(i)
            s%HTE(i, j) = 0.5_wp*(h(i) + h(ie))
            s%HTN(i, j) = 0.5_wp*(w(j) + w(jn))
            s%dxU(i, j) = 0.5_wp*(w(j) + w(jn))
            s%dyU(i, j) = 0.5_wp*(h(i) + h(ie))
            s%earea(i, j) = 0.5_wp*(w(j)*h(i) + w(j)*h(ie))
            s%narea(i, j) = 0.5_wp*(w(j)*h(i) + w(jn)*h(i))
         end do
      end do
      deallocate (w, h)
   end subroutine set_metrics_nonuniform

   subroutine set_flux_flow(s, c, seed)
      !! Random flow whose volume FLUXES are discretely non-divergent on any metrics: `u HTE` and
      !! `v HTN` are differences of a corner stream function (zero in a ring of three cells at the
      !! boundary, so the basin is closed), peak speed `c` m/s at the reference lengths.
      type(sim_t), intent(inout) :: s
      real(wp), intent(in) :: c
      integer(int64), intent(in) :: seed
      real(wp), allocatable :: psi(:, :)
      integer :: i, j
      integer(int64) :: st

      allocate (psi(0:s%nx + 1, 0:s%ny + 1))
      psi = 0.0_wp
      st = lcg_next(seed)
      do j = 4, s%ny - 3
         do i = 4, s%nx - 3
            st = lcg_next(st)
            psi(i, j) = 2.0_wp*lcg_unit(st) - 1.0_wp
         end do
      end do
      s%u = 0.0_wp
      s%v = 0.0_wp
      do j = 2, s%ny
         do i = 2, s%nx
            s%u(i, j) = 0.5_wp*c*s%dy*(psi(i - 1, j) - psi(i - 1, j - 1))/s%HTE(i - 1, j)
            s%v(i, j) = -0.5_wp*c*s%dx*(psi(i, j - 1) - psi(i - 1, j - 1))/s%HTN(i, j - 1)
         end do
      end do
      call close_faces(s)
      deallocate (psi)
   end subroutine set_flux_flow

   subroutine test_nonuniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, old
      call nonuniform_body(error, s, old)
      call sim_destroy(s)
      call sim_destroy(old)
   end subroutine test_nonuniform

   subroutine nonuniform_body(error, s, old)
      !! Stretched / converging metrics, closed basin, non-divergent fluxes: every conserved total is
      !! conserved, a constant tracer stays constant, no tracer leaves the 5x5 bound (driver check on).
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, old
      type(efp_t) :: q0(7), q1(7)
      type(transport_result_t) :: res
      integer :: it, nv, i, j, k, n
      real(wp) :: dt, worst, worst7, w, wmono, wconst, vv
      integer :: nvt

      dt = 1000.0_wp
      ! (a) fronts: conservation and monotonicity
      call sim_init(s, 30, 26, 5, 7, .false., .false., .true., .true., 1.0e4_wp)
      call sim_init(old, 30, 26, 5, 7, .false., .false., .false., .false., 1.0e4_wp)
      call set_metrics_nonuniform(s)
      call fill_front(s, 17_int64)
      call set_flux_flow(s, cps(s, dt, 0.20_wp), 5_int64)
      call totals(s, q0)
      worst = 0.0_wp
      worst7 = 0.0_wp
      wmono = 0.0_wp
      nvt = 0
      do it = 1, 40
         call copy_state(s, old)
         call sim_step(s, dt, res, wrap=.false.)
         call chk_ierr(error, res, "non-uniform metrics step (driver checks on)")
         if (allocated(error)) return
         call totals(s, q1)
         worst = max(worst, worst_change(q0, q1, 1, 6))
         worst7 = max(worst7, worst_change(q0, q1, 7, 7))
         call mono_check(old, s, nv, w)
         nvt = nvt + nv
         wmono = max(wmono, w)
      end do
      write (error_unit, "(a,es10.3,a,es10.3,a,i0,a,es10.3)") "  [remap] non-uniform metrics, 40 steps: conserved totals change by ", worst, &
         " (surface-temperature content ", worst7, "), 5x5 violations ", nvt, ", worst excursion ", wmono
      if (.not. (worst <= tol_cons() .and. worst7 <= TOL_T7)) then
         call test_failed(error, "a total moved on non-uniform metrics")
         return
      end if
      call check(error, nvt, 0, "no tracer left the 5x5 bound on non-uniform metrics")
      if (allocated(error)) return
      call sim_destroy(s)
      call sim_destroy(old)
      ! (b) a constant tracer under a varying area
      call sim_init(s, 30, 26, 5, 7, .false., .false., .true., .false., 1.0e4_wp)
      call set_metrics_nonuniform(s)
      call fill_smooth(s, 0.3_wp, .false.)
      s%aicen(:, :, 2) = 0.0_wp
      do n = 1, s%ncat
         do j = 1, s%ny
            do i = 1, s%nx
               if (s%hm(i, j) > 0.5_wp) then
                  s%vicen(i, j, n) = s%aicen(i, j, n)*HMID(n)
                  s%vsnon(i, j, n) = s%aicen(i, j, n)*0.1_wp
                  s%Tsfcn(i, j, n) = -7.0_wp
                  s%qsnon(i, j, n, 1) = -SHIFT*1.2_wp
                  do k = 1, s%nil
                     s%qicen(i, j, n, k) = -2.0e8_wp - 1.0e6_wp*real(k, wp)
                     s%sicen(i, j, n, k) = 4.0_wp + 0.5_wp*real(k, wp)
                  end do
               else
                  s%aicen(i, j, n) = 0.0_wp
                  s%vicen(i, j, n) = 0.0_wp
                  s%vsnon(i, j, n) = 0.0_wp
               end if
            end do
         end do
      end do
      do j = 1, s%ny
         do i = 1, s%nx
            s%aice0(i, j) = merge(1.0_wp - sum(s%aicen(i, j, :)), 0.0_wp, s%hm(i, j) > 0.5_wp)
         end do
      end do
      call set_flux_flow(s, cps(s, dt, 0.20_wp), 9_int64)
      do it = 1, 30
         call sim_step(s, dt, res, wrap=.false.)
         call chk_ierr(error, res, "non-uniform constant-tracer step")
         if (allocated(error)) return
      end do
      wconst = 0.0_wp
      do n = 1, s%ncat
         do j = s%js, s%je
            do i = s%is, s%ie
               if (.not. (s%aicen(i, j, n) > PUNY)) cycle
               vv = abs(s%vicen(i, j, n)/s%aicen(i, j, n)/HMID(n) - 1.0_wp)
               wconst = max(wconst, vv, abs(s%Tsfcn(i, j, n)/(-7.0_wp) - 1.0_wp), &
                            abs(s%vsnon(i, j, n)/s%aicen(i, j, n)/0.1_wp - 1.0_wp))
               do k = 1, s%nil
                  wconst = max(wconst, abs(s%qicen(i, j, n, k)/(-2.0e8_wp - 1.0e6_wp*real(k, wp)) - 1.0_wp), &
                               abs(s%sicen(i, j, n, k)/(4.0_wp + 0.5_wp*real(k, wp)) - 1.0_wp))
               end do
            end do
         end do
      end do
write (error_unit, "(a,es10.3)") "  [remap] non-uniform metrics, constant tracers after 30 steps: worst relative departure ", wconst
      if (.not. (wconst <= 1.0e-13_wp)) call test_failed(error, "a constant tracer did not stay constant on non-uniform metrics")
   end subroutine nonuniform_body

   ! ------------------------------------------------------------------
   ! (13) configuration refusals and the halo
   ! ------------------------------------------------------------------

   subroutine test_cfg(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s, ref
      call cfg_body(error, s, ref)
      call sim_destroy(s)
      call sim_destroy(ref)
   end subroutine test_cfg

   subroutine call_driver(s, is, ie, js, je, res, with_metrics)
      type(sim_t), intent(inout) :: s
      integer, intent(in) :: is, ie, js, je
      type(transport_result_t), intent(out) :: res
      logical, intent(in) :: with_metrics
      if (with_metrics) then
         call ice_step_transport(s%w, s%p, 1000.0_wp, ICE_TRANSPORT_POST_NONE, is, ie, js, je, s%u, s%v, s%HTE, s%HTN, &
                                 s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, s%aicen, s%vicen, &
                                 s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, &
                                 s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res, &
                                 dxU=s%dxU, dyU=s%dyU, earea=s%earea, narea=s%narea, uvm=s%uvm)
      else
         call ice_step_transport(s%w, s%p, 1000.0_wp, ICE_TRANSPORT_POST_NONE, is, ie, js, je, s%u, s%v, s%HTE, s%HTN, &
                                 s%tarea, s%dxT, s%dyT, s%epm, s%npm, s%hm, s%Tf, s%hin_max, s%aice0, s%aicen, s%vicen, &
                                 s%vsnon, s%Tsfcn, s%qicen, s%sicen, s%qsnon, s%aice, s%vice, s%vsno, s%Tsfc_c, s%qice_c, &
                                 s%sice_c, s%qsno_c, s%fresh, s%fsalt, s%fhocn, s%ice_err, res)
      end if
   end subroutine call_driver

   subroutine cfg_body(error, s, ref)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s, ref
      type(transport_result_t) :: res
      real(wp) :: worst

      call sim_init(s, 16, 12, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call sim_init(ref, 16, 12, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call fill_smooth(s, 0.1_wp, .false.)
      call set_velocity(s, FLOW_UNIFORM, 3.0_wp, 1.0_wp)
      call sim_wrap(s)
      call copy_state(s, ref)
      call sim_map(s)
      call sim_push(s)
      call call_driver(s, s%is, s%ie, s%js, s%je, res, .false.)
      call check(error, res%cfg, TRANSPORT_CFG_PARAMS, "remap without the corner metrics is refused")
      if (allocated(error)) return
      call call_driver(s, s%is - 1, s%ie, s%js, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_RANGE, "a window with two ghost cells is refused")
      if (allocated(error)) return
      call call_driver(s, s%is, s%ie + 1, s%js, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_RANGE, "a window with two ghost cells (east) is refused")
      if (allocated(error)) return
      call call_driver(s, s%is, s%ie, s%js - 1, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_RANGE, "a window with two ghost cells (south) is refused")
      if (allocated(error)) return
      s%p%ncat = 4
      call call_driver(s, s%is, s%ie, s%js, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_PARAMS, "ncat mismatch refused")
      if (allocated(error)) return
      s%p%ncat = 5
      ! a land mask that is not 0/1
      s%hm(8, 6) = 0.5_wp
      !$acc update device(s%hm)
      call call_driver(s, s%is, s%ie, s%js, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_REMAP, "a land mask other than 0/1 is refused")
      if (allocated(error)) return
      s%hm(8, 6) = 1.0_wp
      !$acc update device(s%hm)
      call sim_pull(s)
      call state_diff(s, ref, worst, .false.)
      if (.not. (worst == 0.0_wp)) call test_failed(error, "refused calls leave the state alone")
      if (allocated(error)) return
      ! a workspace built without the remap buffers
      call sim_destroy(s)
      call sim_init_noremap(s)
      call sim_map(s)
      call sim_push(s)
      call call_driver(s, s%is, s%ie, s%js, s%je, res, .true.)
      call check(error, res%cfg, TRANSPORT_CFG_REMAP, "remap on a workspace built without the remap buffers is refused")
   end subroutine cfg_body

   subroutine sim_init_noremap(s)
      type(sim_t), intent(out) :: s
      integer :: ierr
      call sim_init(s, 16, 12, 5, 7, .true., .true., .false., .false., 1.0e4_wp)
      call s%w%destroy()
      call s%w%init(s%nx, s%ny, 5, 7, 1, .false., ierr)
   end subroutine sim_init_noremap

   ! ------------------------------------------------------------------
   ! (14) throughput, launches, memory
   ! ------------------------------------------------------------------

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t) :: s
      call throughput_body(error, s)
      call sim_destroy(s)
   end subroutine test_throughput

   subroutine stage_breakdown(s)
      !! Time of each stage of the remap on the device (category 2 for the per-category stages; one
      !! step costs the geometry once and the per-category stages six times), `!$acc wait` around each.
      type(sim_t), intent(inout) :: s
      integer, parameter :: NREP = 5
      integer :: rep, cfg, nbad, nneg
      integer(int64) :: t0, t1, rate
      real(wp) :: tg, tc, tf, tu, ti
      integer(int32), allocatable :: ierr(:, :)

      allocate (ierr(s%nx, s%ny))
      ierr = 0_int32
      !$acc enter data copyin(ierr)
      call system_clock(t0, rate)
      do rep = 1, NREP
         call remap_geometry(s%w%gw, s%is, s%ie, s%js, s%je, 1000.0_wp, 1, PUNY, s%u, s%v, s%dxU, s%dyU, s%HTE, s%HTN, s%earea, &
                             s%narea, s%epm, s%npm, s%uvm, ierr)
         !$acc wait
      end do
      call system_clock(t1)
      tg = real(t1 - t0, wp)/real(NREP*rate, wp)*1000.0_wp
      call system_clock(t0)
      do rep = 1, NREP
         call remap_input_check(s%w%a, s%is, s%ie, s%js, s%je, s%hm, nbad, cfg)
         !$acc wait
      end do
      call system_clock(t1)
      ti = real(t1 - t0, wp)/real(NREP*rate, wp)*1000.0_wp
      call system_clock(t0)
      do rep = 1, NREP
         call remap_recon_batch(s%w%a, s%w%rb, s%is - 1, s%ie + 1, s%js - 1, s%je + 1, 2, 2, PUNY, s%hm, cfg)
         !$acc wait
      end do
      call system_clock(t1)
      tc = real(t1 - t0, wp)/real(NREP*rate, wp)*1000.0_wp
      call system_clock(t0)
      do rep = 1, NREP
         call remap_flux_integrals(s%w%a, s%w%gw, s%w%rb, s%w%fw, s%is, s%ie, s%js, s%je, 2, 2, PUNY)
         !$acc wait
      end do
      call system_clock(t1)
      tf = real(t1 - t0, wp)/real(NREP*rate, wp)*1000.0_wp
      call system_clock(t0)
      do rep = 1, NREP
         call remap_update(s%w%a, s%w%b, s%w%fw, s%is, s%ie, s%js, s%je, 2, 2, PUNY, s%tarea, nbad, nneg)
         !$acc wait
      end do
      call system_clock(t1)
      tu = real(t1 - t0, wp)/real(NREP*rate, wp)*1000.0_wp
      !$acc exit data delete(ierr)
      write (error_unit, "(a,5(a,f8.2))") "  [remap] stage times [ms]:", " geometry ", tg, ", input check ", ti, &
         ", reconstruction (1 category) ", tc, ", fluxes e+n (1 category) ", tf, ", update (1 category) ", tu
   end subroutine stage_breakdown

   subroutine throughput_body(error, s)
      type(error_type), allocatable, intent(out) :: error
      type(sim_t), intent(inout) :: s
      type(transport_result_t) :: res
      integer :: nxp, nyp, i, j, n, k, rep, nxg, nyg
      integer, parameter :: NREP = 5
      integer(int64) :: t0, t1, rate
      real(wp) :: tstep, cells, bpc, full_cells, dt
      real(wp) :: gb
      character(len=32) :: envv
      integer :: ist, mib

#ifdef RDB_GPU_OFFLOAD
      nxp = 1440; nyp = 1152
#else
      nxp = 384; nyp = 384
#endif
      dt = 1000.0_wp
      ! RDB_REMAP_TEST_MIB: the flux + reconstruction budget [MiB] of one category batch (default 1024)
      call get_environment_variable("RDB_REMAP_TEST_MIB", envv, status=ist)
      mib = 1024
      if (ist == 0) read (envv, *, iostat=ist) mib
      call sim_init(s, nxp, nyp, 5, 7, .true., .true., .false., .false., 1.0e4_wp, bytes_max=int(mib, int64)*1048576_int64)
      do n = 1, 5
         do j = 1, s%ny
            do i = 1, s%nx
               s%aicen(i, j, n) = 0.04_wp + 0.01_wp*real(mod(i + 2*j + n, 7), wp)
               s%vicen(i, j, n) = s%aicen(i, j, n)*(0.5_wp + 0.3_wp*real(n, wp))
               s%vsnon(i, j, n) = 0.1_wp*s%aicen(i, j, n)
               s%Tsfcn(i, j, n) = -10.0_wp + 0.1_wp*real(mod(i + j, 5), wp)
               s%qsnon(i, j, n, 1) = -1.1e8_wp
               do k = 1, 7
                  s%qicen(i, j, n, k) = -2.0e8_wp - 1.0e5_wp*real(mod(i + 3*j, 4), wp)
                  s%sicen(i, j, n, k) = 5.0_wp
               end do
            end do
         end do
      end do
      s%aice0 = 0.5_wp
      call set_velocity(s, FLOW_UNIFORM, cps(s, dt, 0.3_wp), cps(s, dt, -0.2_wp))
      call sim_wrap(s)
      call sim_map(s)
      call sim_push(s)
      call run_one(s, res)
      call chk_ierr(error, res, "throughput warm-up")
      if (allocated(error)) return
      call system_clock(t0, rate)
      do rep = 1, NREP
         call run_one(s, res)
      end do
      call system_clock(t1)
      call sim_pull(s)
      tstep = real(t1 - t0, wp)/real(NREP*rate, wp)
      nxg = s%nx
      nyg = s%ny
      cells = real(nxg, wp)*real(nyg, wp)*5.0_wp
      bpc = real(s%w%bytes(), wp)/(real(nxg, wp)*real(nyg, wp))
      full_cells = 1446.0_wp*1158.0_wp
      gb = bpc*full_cells/1.0e9_wp
    write (error_unit, "(a,i0,a,i0,a,es10.3,a,f9.2,a,i0,a)") "  [remap] ", nxp, "x", nyp, "x5x7: full step (check off, no post) ", &
         cells/tstep, " cat-cells/s (", tstep*1000.0_wp, " ms, ", res%n_launch, " launches; flux batch of the default budget)"
      write (error_unit, "(a,i0,a,f8.1,a,f8.2,a)") "  [remap] workspace bytes ", s%w%bytes(), " = ", bpc, " per cell; at 1440x1152x5 (", gb, " GB)"
      write (error_unit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "  [remap] of which geometry ", s%w%gw%bytes(), ", reconstruction ", s%w%rb%bytes(), &
         ", flux ", s%w%fw%bytes(), ", packed 2x ", s%w%a%bytes() + s%w%b%bytes(), ", flux batch (category indices) ", s%w%fw%nbatch
      call check(error, ieee_is_finite(cells/tstep), "finite throughput")
      call stage_breakdown(s)
   contains
      subroutine run_one(sim, r)
         type(sim_t), intent(inout) :: sim
         type(transport_result_t), intent(out) :: r
         call ice_step_transport(sim%w, sim%p, dt, ICE_TRANSPORT_POST_NONE, sim%is, sim%ie, sim%js, sim%je, sim%u, sim%v, &
                                 sim%HTE, sim%HTN, sim%tarea, sim%dxT, sim%dyT, sim%epm, sim%npm, sim%hm, sim%Tf, &
                                 sim%hin_max, sim%aice0, sim%aicen, sim%vicen, sim%vsnon, sim%Tsfcn, sim%qicen, &
                                 sim%sicen, sim%qsnon, sim%aice, sim%vice, sim%vsno, sim%Tsfc_c, sim%qice_c, sim%sice_c, &
                                 sim%qsno_c, sim%fresh, sim%fsalt, sim%fhocn, sim%ice_err, r, &
                                 dxU=sim%dxU, dyU=sim%dyU, earea=sim%earea, narea=sim%narea, uvm=sim%uvm)
      end subroutine run_one
   end subroutine throughput_body

end module test_cice_remap_flux
