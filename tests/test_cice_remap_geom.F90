!! Incremental-remapping geometry (D10.3): departure points, edge areas and the departure-region
!! triangles of the C-grid transport, and their acceptance set.
module test_cice_remap_geom
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Tests written for the roundabout port of the CICE incremental-remapping geometry (CICE
   ! Consortium, https://github.com/CICE-Consortium); the independent oracles (a polygon clipper, a
   ! solver of the edge-flux-adjusted midpoint, closed forms, a stdlib Python transcript) follow the
   ! published scheme, not the port.  Modified by the roundabout developers; see NOTICE for the
   ! full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   !! **What is pinned.**  (1) With the edge flux adjustment OFF the triangles of an edge are exactly the
   !! departure quadrilateral: against an INDEPENDENT Sutherland-Hodgman clipper (plain host code,
   !! sharing nothing with the kernel) the area attributed to every source cell equals the clipped
   !! polygon area times the cell's area factor, to 1e-13 of the cell area, over 300000 random
   !! departure-point pairs, together with the translation of the vertices into the source cell
   !! frame; (2) with it ON, the triangle areas sum to the prescribed edge area to round-off and the
   !! area of EVERY source cell (side cells from the clipper, the centre column from an independent
   !! solve of the shifted midpoint) equals the oracle; (3) a branch coverage counter proves every
   !! branch of the table (4 left, 4 right, 12 centre, 3 flux-adjustment) was run, and lists them;
   !! (4) on a whole grid, with a uniform metric: the area leaving each cell through its four edges
   !! equals the cell area minus the area of its back-traced image (no gap, no double count at the
   !! corners), each edge equals its swept quadrilateral, and a single moved corner gives the closed
   !! forms `+-(a+-b)/2`; (5) uniform translation, solid-body rotation and a deformation field on a
   !! non-uniform metric (the area factors of every triangle are one of the documented neighbours');
   !! (6) the one-cell channel, zero velocity (everything exactly zero, stale data cleared), CFL
   !! (both signs, both components, crossing points), NaN, mirror symmetry in x and y; (7) a stdlib
   !! Python transcript at 1e-12; (8) every output of every edge against a host reference, the ring
   !! untouched; (9) throughput.
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_remap_geom, only: remap_geom_workspace_t, remap_corner_velocity, remap_departure_points, &
                                 remap_edge_areas, remap_locate_east, remap_locate_north, remap_geometry, &
                                 remap_locate_edge, remap_src_di, remap_src_dj, REMAP_NSLOT, REMAP_EDGE_NORTH, &
                                 REMAP_EDGE_EAST
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_remap_geom_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: SENT = 7.25e3_wp
      !! Sentinel the output buffers are filled with: proves what the kernels do not write.
   integer, parameter :: NG = 3
   integer, parameter :: NUU = 5
   integer, parameter :: NLAB = 6, NWORK = 40, NPENT = 5, NQ = 8, NCC = 12
      !! Cell labels, polygon work size, pentagon vertices, test queries, centre branches.

   type :: cov_t
      !! Branch coverage: counts of each branch taken (index 0 = none).
      integer :: left(0:4) = 0, right(0:4) = 0, centre(0:NCC) = 0, efa(0:4) = 0
      integer :: nsamp = 0, nflag_cfl = 0, nflag_other = 0, ncross = 0
      real(wp) :: worst_area = 0.0_wp, worst_sum = 0.0_wp, worst_vert = 0.0_wp, worst_mom = 0.0_wp
   end type cov_t

   type :: geo_t
      !! A host-owned C-grid geometry problem: metrics, face velocities, the workspace, the error word.
      integer :: nx = 0, ny = 0, is = 0, ie = 0, js = 0, je = 0
      logical :: mapped = .false.
      real(wp) :: dt = 1.0_wp
      real(wp), allocatable :: u(:, :), v(:, :), dxU(:, :), dyU(:, :), HTE(:, :), HTN(:, :), earea(:, :), narea(:, :)
      real(wp), allocatable :: epm(:, :), npm(:, :), uvm(:, :), hm(:, :)
         !! Open masks of the E / N faces and NE corners, and the T-cell mask (all 1 unless a test sets land).
      integer(int32), allocatable :: ierr(:, :)
      type(remap_geom_workspace_t) :: ws
   end type geo_t

#include "test_cice_remap_geom_golden.inc"

contains

   subroutine collect_cice_remap_geom_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("source_code_roundtrip", test_src_code), &
                  new_unittest("every_branch_vs_clipper_no_efa", test_clipper_noefa), &
                  new_unittest("efa_sum_and_cells_vs_oracle", test_efa_oracle), &
                  new_unittest("case_table_coverage", test_coverage), &
                  new_unittest("python_transcript", test_python_transcript), &
                  new_unittest("corner_consistency_uniform_metric", test_corner_identity), &
                  new_unittest("single_moved_corner_closed_form", test_single_corner), &
                  new_unittest("uniform_translation", test_translation), &
                  new_unittest("solid_body_rotation", test_rotation), &
                  new_unittest("deformation_nonuniform_metric_efa", test_deformation), &
                  new_unittest("one_cell_channel", test_channel), &
                  new_unittest("zero_velocity_clears_everything", test_zero), &
                  new_unittest("cfl_bits_both_signs_and_crossing", test_cfl), &
                  new_unittest("nan_in_raises_bit", test_nan), &
                  new_unittest("mirror_symmetry_x_and_y", test_mirror), &
                  new_unittest("coast_and_island_no_flux_through_walls", test_coast), &
                  new_unittest("stated_input_halo_is_enough", test_halo_window), &
                  new_unittest("kernels_vs_host_reference_every_output", test_vs_host), &
                  new_unittest("throughput", test_throughput) &
                  ]
   end subroutine collect_cice_remap_geom_tests

   ! ------------------------------------------------------------------
   ! independent oracles: polygon clipper, cell tables
   ! ------------------------------------------------------------------

   subroutine clip_keep(n, x, y, a, b, c)
      !! Sutherland-Hodgman against the half-plane `a x + b y + c >= 0`, in place.
      integer, intent(inout) :: n
      real(wp), intent(inout) :: x(:), y(:)
      real(wp), intent(in) :: a, b, c
      real(wp) :: sp, sq, t
      real(wp) :: xo(size(x)), yo(size(x))
      integer :: m, k, kn

      m = 0
      do k = 1, n
         kn = mod(k, n) + 1
         sp = a*x(k) + b*y(k) + c
         sq = a*x(kn) + b*y(kn) + c
         if (sp >= 0.0_wp) then
            m = m + 1
            xo(m) = x(k)
            yo(m) = y(k)
            if (sq < 0.0_wp) then
               t = sp/(sp - sq)
               m = m + 1
               xo(m) = x(k) + t*(x(kn) - x(k))
               yo(m) = y(k) + t*(y(kn) - y(k))
            end if
         else if (sq >= 0.0_wp) then
            t = sp/(sp - sq)
            m = m + 1
            xo(m) = x(k) + t*(x(kn) - x(k))
            yo(m) = y(k) + t*(y(kn) - y(k))
         end if
      end do
      n = m
      x(1:m) = xo(1:m)
      y(1:m) = yo(1:m)
   end subroutine clip_keep

   function poly_area(n, x, y) result(a)
      !! Signed shoelace area (counter-clockwise positive).
      integer, intent(in) :: n
      real(wp), intent(in) :: x(:), y(:)
      real(wp) :: a
      integer :: k, kn
      a = 0.0_wp
      do k = 1, n
         kn = mod(k, n) + 1
         a = a + x(k)*y(kn) - x(kn)*y(k)
      end do
      a = 0.5_wp*a
   end function poly_area

   function clipped_area(n0, x0, y0, label) result(a)
      !! Signed area of the polygon inside the cell `label` (1 TL, 2 TR, 3 BL, 4 BR, 5 TC, 6 BC) of the
      !! edge frame (top = y > 0; left, centre, right = x < -1/2, between, > 1/2); labels 7 and 8 are the
      !! half-planes y >= 0 and y <= 0 whatever x (the table puts a shifted midpoint beyond the column wholly in TC or BC).
      integer, intent(in) :: n0, label
      real(wp), intent(in) :: x0(:), y0(:)
      real(wp) :: a
      real(wp) :: xlo, xhi, ylo, yhi
      real(wp) :: x(NWORK), y(NWORK)
      integer :: n

      select case (label)
      case (1, 3)
         xlo = -9.0_wp
         xhi = -0.5_wp
      case (2, 4)
         xlo = 0.5_wp
         xhi = 9.0_wp
      case (7, 8)
         xlo = -9.0_wp
         xhi = 9.0_wp
      case default
         xlo = -0.5_wp
         xhi = 0.5_wp
      end select
      if (label == 1 .or. label == 2 .or. label == 5 .or. label == 7) then
         ylo = 0.0_wp
         yhi = 9.0_wp
      else
         ylo = -9.0_wp
         yhi = 0.0_wp
      end if
      n = n0
      x(1:n) = x0(1:n)
      y(1:n) = y0(1:n)
      call clip_keep(n, x, y, 1.0_wp, 0.0_wp, -xlo)
      call clip_keep(n, x, y, -1.0_wp, 0.0_wp, xhi)
      call clip_keep(n, x, y, 0.0_wp, 1.0_wp, -ylo)
      call clip_keep(n, x, y, 0.0_wp, -1.0_wp, yhi)
      a = poly_area(n, x, y)
   end function clipped_area

   subroutine clipped_moments(n0, x0, y0, label, a, mx, my)
      !! Signed area and first moments `int x dA`, `int y dA` of the polygon inside the cell `label`
      !! (as `clipped_area`), edge frame.
      integer, intent(in) :: n0, label
      real(wp), intent(in) :: x0(:), y0(:)
      real(wp), intent(out) :: a, mx, my
      real(wp) :: xlo, xhi, ylo, yhi, x(NWORK), y(NWORK), cr
      integer :: n, k, kn

      select case (label)
      case (1, 3)
         xlo = -9.0_wp
         xhi = -0.5_wp
      case (2, 4)
         xlo = 0.5_wp
         xhi = 9.0_wp
      case (7, 8)
         xlo = -9.0_wp
         xhi = 9.0_wp
      case default
         xlo = -0.5_wp
         xhi = 0.5_wp
      end select
      if (label == 1 .or. label == 2 .or. label == 5 .or. label == 7) then
         ylo = 0.0_wp
         yhi = 9.0_wp
      else
         ylo = -9.0_wp
         yhi = 0.0_wp
      end if
      n = n0
      x(1:n) = x0(1:n)
      y(1:n) = y0(1:n)
      call clip_keep(n, x, y, 1.0_wp, 0.0_wp, -xlo)
      call clip_keep(n, x, y, -1.0_wp, 0.0_wp, xhi)
      call clip_keep(n, x, y, 0.0_wp, 1.0_wp, -ylo)
      call clip_keep(n, x, y, 0.0_wp, -1.0_wp, yhi)
      a = 0.0_wp
      mx = 0.0_wp
      my = 0.0_wp
      do k = 1, n
         kn = mod(k, n) + 1
         cr = x(k)*y(kn) - x(kn)*y(k)
         a = a + cr
         mx = mx + (x(k) + x(kn))*cr
         my = my + (y(k) + y(kn))*cr
      end do
      a = 0.5_wp*a
      mx = mx/6.0_wp
      my = my/6.0_wp
   end subroutine clipped_moments

   subroutine label_offset(iedge, label, di, dj)
      !! Offset of cell `label` from the edge's own cell (own tables).
      integer, intent(in) :: iedge, label
      integer, intent(out) :: di, dj
      integer :: tn_i(NLAB), tn_j(NLAB), te_i(NLAB), te_j(NLAB)
      tn_i = [-1, 1, -1, 1, 0, 0]
      tn_j = [1, 1, 0, 0, 1, 0]
      te_i = [1, 1, 0, 0, 1, 0]
      te_j = [1, -1, 1, -1, 0, 0]
      if (iedge == REMAP_EDGE_NORTH) then
         di = tn_i(label)
         dj = tn_j(label)
      else
         di = te_i(label)
         dj = te_j(label)
      end if
   end subroutine label_offset

   subroutine to_source_frame(iedge, label, a, mxe, mye, mxs, mys)
      !! First moments of a region of area `a` taken in the edge frame, expressed in the frame of
      !! source cell `label` (north: translate; east: translate and turn a quarter back).
      integer, intent(in) :: iedge, label
      real(wp), intent(in) :: a, mxe, mye
      real(wp), intent(out) :: mxs, mys
      integer :: di, dj
      call label_offset(iedge, label, di, dj)
      if (iedge == REMAP_EDGE_NORTH) then
         mxs = mxe - real(di, wp)*a
         mys = mye + (0.5_wp - real(dj, wp))*a
      else
         mxs = mye + (0.5_wp - real(di, wp))*a
         mys = -mxe - real(dj, wp)*a
      end if
   end subroutine to_source_frame

   subroutine kernel_cell_moments(iedge, sc, ar, px, py, kmx, kmy)
      !! `sum_s area_s * centroid_s` per source-cell label from the stored vertices (source frame).
      integer, intent(in) :: iedge
      integer(int8), intent(in) :: sc(:)
      real(wp), intent(in) :: ar(:), px(:, :), py(:, :)
      real(wp), intent(out) :: kmx(NLAB), kmy(NLAB)
      integer :: g, l
      kmx = 0.0_wp
      kmy = 0.0_wp
      do g = 1, REMAP_NSLOT
         if (sc(g) == 0_int8) cycle
         l = label_of(iedge, remap_src_di(sc(g)), remap_src_dj(sc(g)))
         if (l < 1) cycle
         kmx(l) = kmx(l) + ar(g)*(px(1, g) + px(2, g) + px(3, g))/3.0_wp
         kmy(l) = kmy(l) + ar(g)*(py(1, g) + py(2, g) + py(3, g))/3.0_wp
      end do
   end subroutine kernel_cell_moments

   function label_of(iedge, di, dj) result(l)
      !! Cell label of the offset `(di, dj)` from the edge's own cell (own tables, 0 = not a neighbour).
      integer, intent(in) :: iedge, di, dj
      integer :: l
      l = 0
      if (iedge == REMAP_EDGE_NORTH) then
         if (di == -1 .and. dj == 1) l = 1
         if (di == 1 .and. dj == 1) l = 2
         if (di == -1 .and. dj == 0) l = 3
         if (di == 1 .and. dj == 0) l = 4
         if (di == 0 .and. dj == 1) l = 5
         if (di == 0 .and. dj == 0) l = 6
      else
         if (di == 1 .and. dj == 1) l = 1
         if (di == 1 .and. dj == -1) l = 2
         if (di == 0 .and. dj == 1) l = 3
         if (di == 0 .and. dj == -1) l = 4
         if (di == 1 .and. dj == 0) l = 5
         if (di == 0 .and. dj == 0) l = 6
      end if
   end function label_of

   subroutine cell_bounds(label, xlo, xhi, ylo, yhi)
      integer, intent(in) :: label
      real(wp), intent(out) :: xlo, xhi, ylo, yhi
      select case (label)
      case (1, 3)
         xlo = -1.5_wp
         xhi = -0.5_wp
      case (2, 4)
         xlo = 0.5_wp
         xhi = 1.5_wp
      case default
         xlo = -0.5_wp
         xhi = 0.5_wp
      end select
      if (label == 1 .or. label == 2 .or. label == 5) then
         ylo = 0.0_wp
         yhi = 1.0_wp
      else
         ylo = -1.0_wp
         yhi = 0.0_wp
      end if
   end subroutine cell_bounds

   subroutine efa_oracle(dlx, dly, drx, dry, ea, fc, osum, omx, omy, ok)
      !! Independent expectation of the area in each source cell with the edge flux adjustment, for a
      !! UNIFORM area factor `fc`: side cells from the clipper on the unadjusted quadrilateral; the
      !! centre column from the polygon whose midpoint is solved (linear in the shift) for the area
      !! left over after the sides, with the lone triangle taken off first when the departure
      !! points lie on opposite sides.
      real(wp), intent(in) :: dlx, dly, drx, dry, ea, fc
      real(wp), intent(out) :: osum(NLAB), omx(NLAB), omy(NLAB)
      !! Area and first moments (edge frame) per source cell.
      logical, intent(out) :: ok
      real(wp) :: a0, a1, tt, target, side, slope
      real(wp) :: qx(4), qy(4), px(NPENT), py(NPENT), lx(3), ly(3)
      real(wp) :: xdl, ydl, xdr, ydr, xl, yl, xr, yr, xic, nx_, ny_, xm, ym
      integer :: k

      ok = .true.
      osum = 0.0_wp
      omx = 0.0_wp
      omy = 0.0_wp
      xdl = -0.5_wp + dlx
      ydl = dly
      xdr = 0.5_wp + drx
      ydr = dry
      qx = [-0.5_wp, xdl, xdr, 0.5_wp]
      qy = [0.0_wp, ydl, ydr, 0.0_wp]
      do k = 1, 4
         call add(k, 4, qx, qy, k)
      end do
      side = sum(osum(1:4))
      target = (ea - side)/fc
      slope = (ydr - ydl)/(xdr - xdl)
      xl = xdl
      yl = ydl
      if (xdl < -0.5_wp) then
         xl = -0.5_wp
         yl = ydl + slope*(-0.5_wp - xdl)
      end if
      xr = xdr
      yr = ydr
      if (xdr > 0.5_wp) then
         xr = 0.5_wp
         yr = ydl + slope*(0.5_wp - xdl)
      end if
      if (yl*yr >= 0.0_wp) then
         nx_ = yr - yl
         ny_ = -(xr - xl)
         ! the reference shifts the ORIGINAL midpoint of the departure points (not that of the redefined ones)
         xm = 0.5_wp*(xdl + xdr)
         ym = 0.5_wp*(ydl + ydr)
         px = [-0.5_wp, xl, xm, xr, 0.5_wp]
         py = [0.0_wp, yl, ym, yr, 0.0_wp]
         a0 = poly_area(5, px, py)
         px(3) = xm + nx_
         py(3) = ym + ny_
         a1 = poly_area(5, px, py)
         tt = (target - a0)/(a1 - a0)
         px(3) = xm + tt*nx_
         py(3) = ym + tt*ny_
         call add(5, 5, px, py, 7)
         call add(6, 5, px, py, 8)
      else
         xic = xl - yl*(xr - xl)/(yr - yl)
         if (xic < 0.0_wp .and. xic > -0.5_wp) then
            lx = [-0.5_wp, xl, xic]
            ly = [0.0_wp, yl, 0.0_wp]
            call add(5, 3, lx, ly, 7)
            call add(6, 3, lx, ly, 8)
            nx_ = yr
            ny_ = -(xr - xic)
            xm = 0.5_wp*(xr + xic)
            ym = 0.5_wp*yr
            px(1:4) = [xic, xm, xr, 0.5_wp]
            py(1:4) = [0.0_wp, ym, yr, 0.0_wp]
            a0 = poly_area(4, px, py) + poly_area(3, lx, ly)
            px(2) = xm + nx_
            py(2) = ym + ny_
            a1 = poly_area(4, px, py) + poly_area(3, lx, ly)
            tt = (target - a0)/(a1 - a0)
            px(2) = xm + tt*nx_
            py(2) = ym + tt*ny_
            call add(5, 4, px, py, 7)
            call add(6, 4, px, py, 8)
         else if (xic >= 0.0_wp .and. xic < 0.5_wp) then
            lx = [xic, xr, 0.5_wp]
            ly = [0.0_wp, yr, 0.0_wp]
            call add(5, 3, lx, ly, 7)
            call add(6, 3, lx, ly, 8)
            nx_ = -yl
            ny_ = -(xic - xl)
            xm = 0.5_wp*(xic + xl)
            ym = 0.5_wp*yl
            px(1:4) = [-0.5_wp, xl, xm, xic]
            py(1:4) = [0.0_wp, yl, ym, 0.0_wp]
            a0 = poly_area(4, px, py) + poly_area(3, lx, ly)
            px(3) = xm + nx_
            py(3) = ym + ny_
            a1 = poly_area(4, px, py) + poly_area(3, lx, ly)
            tt = (target - a0)/(a1 - a0)
            px(3) = xm + tt*nx_
            py(3) = ym + tt*ny_
            call add(5, 4, px, py, 7)
            call add(6, 4, px, py, 8)
         else
            ok = .false.
         end if
      end if
   contains
      subroutine add(lout, n, x, y, lab)
         integer, intent(in) :: lout, n, lab
         real(wp), intent(in) :: x(:), y(:)
         real(wp) :: a, mx, my
         call clipped_moments(n, x, y, lab, a, mx, my)
         osum(lout) = osum(lout) + fc*a
         omx(lout) = omx(lout) + fc*mx
         omy(lout) = omy(lout) + fc*my
      end subroutine add
   end subroutine efa_oracle

   subroutine kernel_cell_sums(iedge, sc, ar, ksum)
      !! Physical area per source-cell label of one edge's output.
      integer, intent(in) :: iedge
      integer(int8), intent(in) :: sc(:)
      real(wp), intent(in) :: ar(:)
      real(wp), intent(out) :: ksum(NLAB)
      integer :: g, l
      ksum = 0.0_wp
      do g = 1, REMAP_NSLOT
         if (sc(g) == 0_int8) cycle
         l = label_of(iedge, remap_src_di(sc(g)), remap_src_dj(sc(g)))
         if (l >= 1) ksum(l) = ksum(l) + ar(g)
      end do
   end subroutine kernel_cell_sums

   ! ------------------------------------------------------------------
   ! sampling sweeps
   ! ------------------------------------------------------------------

   subroutine draw_departure(state, dlx, dly, drx, dry)
      integer(int64), intent(inout) :: state
      real(wp), intent(out) :: dlx, dly, drx, dry
      real(wp) :: s, u0
      state = lcg_next(state)
      u0 = lcg_unit(state)
      if (u0 < 0.3_wp) then
         s = 0.15_wp
      else if (u0 < 0.6_wp) then
         s = 0.5_wp
      else
         s = 1.0_wp
      end if
      state = lcg_next(state)
      dlx = s*(2.0_wp*lcg_unit(state) - 1.0_wp)
      state = lcg_next(state)
      dly = s*(2.0_wp*lcg_unit(state) - 1.0_wp)
      state = lcg_next(state)
      drx = s*(2.0_wp*lcg_unit(state) - 1.0_wp)
      state = lcg_next(state)
      dry = s*(2.0_wp*lcg_unit(state) - 1.0_wp)
   end subroutine draw_departure

   subroutine tally(cov, icase)
      type(cov_t), intent(inout) :: cov
      integer(int32), intent(in) :: icase
      integer :: il, ir, ic, ie
      il = mod(int(icase), 5)
      ir = mod(int(icase)/5, 5)
      ic = mod(int(icase)/25, 13)
      ie = int(icase)/325
      cov%left(il) = cov%left(il) + 1
      cov%right(ir) = cov%right(ir) + 1
      cov%centre(ic) = cov%centre(ic) + 1
      cov%efa(ie) = cov%efa(ie) + 1
   end subroutine tally

   subroutine sweep_noefa(error, cov, nsamp)
      !! Random departure points, edge flux adjustment off, against the clipper.
      type(error_type), allocatable, intent(out) :: error
      type(cov_t), intent(inout) :: cov
      integer, intent(in) :: nsamp
      real(wp), parameter :: FC = 1.0e8_wp
      real(wp), parameter :: TOLA = 1.0e-13_wp*FC
      real(wp) :: dlx, dly, drx, dry, asum
      real(wp) :: f(NLAB), px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      real(wp) :: want, f1, f2, err, xlo, xhi, ylo, yhi, vx, vy, cross, rat, quad
      real(wp) :: qx(4), qy(4), ksum(NLAB)
      real(wp) :: kmx(NLAB), kmy(NLAB), aq, mxe, mye, mxs, mys, ff
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: ic, bits
      integer(int64) :: st
      integer :: n, iedge, l, g, k, di, dj
      character(len=160) :: msg

      st = 12345_int64
      do n = 1, nsamp
         call draw_departure(st, dlx, dly, drx, dry)
         do k = 1, 6
            st = lcg_next(st)
            f(k) = FC*(0.8_wp + 0.4_wp*lcg_unit(st))
         end do
         st = lcg_next(st)
         iedge = merge(REMAP_EDGE_NORTH, REMAP_EDGE_EAST, lcg_unit(st) < 0.5_wp)
         cov%nsamp = cov%nsamp + 1
         call remap_locate_edge(iedge, 0, PUNY, dlx, dly, drx, dry, FC, f(1), f(2), f(3), f(4), f(5), f(6), &
                                0.0_wp, px, py, ar, sc, ic, asum, bits)
         if (1.0_wp + drx - dlx <= 0.0_wp) then
            cov%ncross = cov%ncross + 1
            if (bits /= ICE_ERR_TRANSPORT_CFL) then
               call test_failed(error, "crossing departure points must raise the CFL bit")
               return
            end if
            cycle
         end if
         if (bits /= 0_int32) then
            write (msg, "(a,i0,4(1x,es12.4))") "unexpected bits ", bits, dlx, dly, drx, dry
            call test_failed(error, trim(msg))
            return
         end if
         call tally(cov, ic)
         qx = [-0.5_wp, -0.5_wp + dlx, 0.5_wp + drx, 0.5_wp]
         qy = [0.0_wp, dly, dry, 0.0_wp]
         call kernel_cell_sums(iedge, sc, ar, ksum)
         call kernel_cell_moments(iedge, sc, ar, px, py, kmx, kmy)
         quad = 0.0_wp
         do l = 1, 6
            select case (l)
            case (1)
               f1 = f(1)
               f2 = f(5)
            case (2)
               f1 = f(3)
               f2 = f(6)
            case (3)
               f1 = f(2)
               f2 = f(5)
            case (4)
               f1 = f(4)
               f2 = f(6)
            case default
               f1 = FC
               f2 = FC
            end select
            want = clipped_area(4, qx, qy, l)
            call clipped_moments(4, qx, qy, l, aq, mxe, mye)
            call to_source_frame(iedge, l, aq, mxe, mye, mxs, mys)
            ff = merge(f1, f2, abs(ksum(l) - f1*want) <= abs(ksum(l) - f2*want))
            err = max(abs(kmx(l) - ff*mxs), abs(kmy(l) - ff*mys))
            cov%worst_mom = max(cov%worst_mom, err/FC)
            if (.not. (err <= TOLA)) then
               write (msg, "(a,i0,a,i0,a,4(1x,es12.4))") "cell ", l, " edge ", iedge, ": first moments vs clipper ", &
                  kmx(l), ff*mxs, kmy(l), ff*mys
               call test_failed(error, trim(msg))
               return
            end if
            err = min(abs(ksum(l) - f1*want), abs(ksum(l) - f2*want))
            cov%worst_area = max(cov%worst_area, err/FC)
            if (.not. (err <= TOLA)) then
               write (msg, "(a,i0,a,i0,a,3(1x,es12.4))") "cell ", l, " edge ", iedge, ": area vs clipper ", ksum(l), &
                  f1*want, f2*want
               call test_failed(error, trim(msg))
               return
            end if
            quad = quad + want
         end do
         ! total
         err = abs(asum - sum(ksum))
         if (.not. (err <= TOLA)) then
            call test_failed(error, "areasum is not the sum of the slots")
            return
         end if
         ! vertices: in the source cell, area consistent with the vertices
         do g = 1, REMAP_NSLOT
            if (sc(g) == 0_int8) then
               if (ar(g) /= 0.0_wp) then
                  call test_failed(error, "unused slot has an area")
                  return
               end if
               cycle
            end if
            di = remap_src_di(sc(g))
            dj = remap_src_dj(sc(g))
            l = label_of(iedge, di, dj)
            if (l < 1) then
               call test_failed(error, "source code is not a neighbour of the edge")
               return
            end if
            call cell_bounds(l, xlo, xhi, ylo, yhi)
            do k = 1, 3
               ! back to the edge frame: north: x = xp + di, y = yp - 1/2 + dj; east: x = -(yp + dj), y = xp - 1/2 + di
               if (iedge == REMAP_EDGE_NORTH) then
                  vx = px(k, g) + real(di, wp)
                  vy = py(k, g) - 0.5_wp + real(dj, wp)
               else
                  vx = -(py(k, g) + real(dj, wp))
                  vy = px(k, g) - 0.5_wp + real(di, wp)
               end if
               err = max(xlo - vx, vx - xhi, ylo - vy, vy - yhi)
               cov%worst_vert = max(cov%worst_vert, err)
               if (err > 1.0e-12_wp) cov%nflag_other = cov%nflag_other + 1
               if (.not. (err <= 1.0_wp)) then
                  write (msg, "(a,i0,a,i0,a,es10.2)") "vertex of slot ", g, " outside source cell ", l, " by ", err
                  call test_failed(error, trim(msg))
                  return
               end if
            end do
            cross = 0.5_wp*((px(2, g) - px(1, g))*(py(3, g) - py(1, g)) - (py(2, g) - py(1, g))*(px(3, g) - px(1, g)))
            rat = abs(ar(g))
            err = min(abs(rat - f(min(l, 6))*abs(cross)), abs(rat - FC*abs(cross)))
            do k = 1, 6
               err = min(err, abs(rat - f(k)*abs(cross)))
            end do
            if (.not. (err <= TOLA)) then
               call test_failed(error, "area is not (an area factor) x (area of the vertices)")
               return
            end if
         end do
      end do
   end subroutine sweep_noefa

   subroutine sweep_efa(error, cov, nsamp)
      !! Random departure points and edge areas, flux adjustment on: sum identity and the independent
      !! per-cell oracle (uniform factors).
      type(error_type), allocatable, intent(out) :: error
      type(cov_t), intent(inout) :: cov
      integer, intent(in) :: nsamp
      real(wp), parameter :: FC = 1.0e8_wp
      real(wp), parameter :: TOLA = 1.0e-13_wp*FC
      real(wp) :: dlx, dly, drx, dry, asum, ea
      real(wp) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      real(wp) :: err, u0
      real(wp) :: qx(4), qy(4), ksum(NLAB), osum(NLAB), omx(NLAB), omy(NLAB), kmx(NLAB), kmy(NLAB), mxs, mys
      real(wp) :: fr(NLAB), ea2
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: ic, bits
      integer(int64) :: st
      integer :: n, iedge, l, k
      logical :: ok
      character(len=260) :: msg

      st = 987654321_int64
      do n = 1, nsamp
         call draw_departure(st, dlx, dly, drx, dry)
         st = lcg_next(st)
         iedge = merge(REMAP_EDGE_NORTH, REMAP_EDGE_EAST, lcg_unit(st) < 0.5_wp)
         st = lcg_next(st)
         u0 = lcg_unit(st)
         qx = [-0.5_wp, -0.5_wp + dlx, 0.5_wp + drx, 0.5_wp]
         qy = [0.0_wp, dly, dry, 0.0_wp]
         st = lcg_next(st)
         if (u0 < 0.7_wp) then
            ea = FC*poly_area(4, qx, qy)*(1.0_wp + 0.3_wp*(2.0_wp*lcg_unit(st) - 1.0_wp))
         else
            ea = FC*0.6_wp*(2.0_wp*lcg_unit(st) - 1.0_wp)
         end if
         cov%nsamp = cov%nsamp + 1
         call remap_locate_edge(iedge, 1, PUNY, dlx, dly, drx, dry, FC, FC, FC, FC, FC, FC, FC, ea, px, py, ar, sc, &
                                ic, asum, bits)
         if (1.0_wp + drx - dlx <= 0.0_wp) then
            cov%ncross = cov%ncross + 1
            if (bits /= ICE_ERR_TRANSPORT_CFL) then
               call test_failed(error, "crossing departure points must raise the CFL bit (efa)")
               return
            end if
            cycle
         end if
         if (bits /= 0_int32) then
            ! a refusal must say CFL (the shifted midpoint left the two cells of the edge); a sum that does not close
            ! WITHOUT it, or a non-finite result, is a failure
            if (iand(bits, ICE_ERR_TRANSPORT_CFL) == 0_int32 .or. iand(bits, ICE_ERR_NONFINITE) /= 0_int32) then
               write (msg, "(a,i0,5(1x,es12.4))") "efa raised bits ", bits, dlx, dly, drx, dry, ea/FC
               call test_failed(error, trim(msg))
               return
            end if
            cov%nflag_cfl = cov%nflag_cfl + 1
            cycle
         end if
         call tally(cov, ic)
         cov%worst_sum = max(cov%worst_sum, abs(asum - ea)/FC)
         if (.not. (abs(asum - ea) <= 1.0e-13_wp*FC)) then
            write (msg, "(a,3(1x,es12.4))") "sum(triarea) /= edgearea: ", asum, ea, asum - ea
            call test_failed(error, trim(msg))
            return
         end if
         call efa_oracle(dlx, dly, drx, dry, ea, FC, osum, omx, omy, ok)
         if (.not. ok) then
            call test_failed(error, "oracle: axis crossing outside the centre column")
            return
         end if
         call kernel_cell_sums(iedge, sc, ar, ksum)
         call kernel_cell_moments(iedge, sc, ar, px, py, kmx, kmy)
         do l = 1, 6
            call to_source_frame(iedge, l, osum(l), omx(l), omy(l), mxs, mys)
            err = max(abs(kmx(l) - mxs), abs(kmy(l) - mys))
            cov%worst_mom = max(cov%worst_mom, err/FC)
            if (.not. (err <= TOLA)) then
               write (msg, "(a,i0,a,i0,a,4(1x,es12.4),a,i0)") "cell ", l, " edge ", iedge, " efa first moments vs oracle ", &
                  kmx(l), mxs, kmy(l), mys, " case ", ic
               call test_failed(error, trim(msg))
               return
            end if
            err = abs(ksum(l) - osum(l))
            cov%worst_area = max(cov%worst_area, err/FC)
            if (.not. (err <= TOLA)) then
               write (msg, "(a,i0,a,i0,a,3(1x,es12.4),a,i0,a,5es11.3)") "cell ", l, " edge ", iedge, " efa area vs oracle ", &
                  ksum(l), osum(l), err/FC, " case ", ic, " in ", dlx, dly, drx, dry, ea/FC
               call test_failed(error, trim(msg))
               return
            end if
         end do
         ! non-uniform area factors: the adjustment still closes the sum (or refuses with the CFL bit)
         do k = 1, NLAB
            st = lcg_next(st)
            fr(k) = FC*(0.8_wp + 0.4_wp*lcg_unit(st))
         end do
         call remap_locate_edge(iedge, 1, PUNY, dlx, dly, drx, dry, fr(1), fr(2), fr(3), fr(4), fr(5), fr(6), fr(3), ea, &
                                px, py, ar, sc, ic, asum, bits)
         ea2 = ea
         if (bits == 0_int32) then
            if (.not. (abs(asum - ea2) <= 1.0e-13_wp*FC)) then
               call test_failed(error, "non-uniform factors: sum(triarea) /= edgearea")
               return
            end if
         else if (iand(bits, ICE_ERR_TRANSPORT_CFL) == 0_int32 .or. iand(bits, ICE_ERR_NONFINITE) /= 0_int32) then
            call test_failed(error, "non-uniform factors: unexpected error bits")
            return
         end if
      end do
   end subroutine sweep_efa

   subroutine test_clipper_noefa(error)
      type(error_type), allocatable, intent(out) :: error
      type(cov_t) :: cov
      call sweep_noefa(error, cov, 300000)
      if (allocated(error)) return
      write (error_unit, "(a,i0,a,es9.2,a,es9.2,a,i0,a,es9.2)") "  [remap_geom] clipper (no efa): samples ", cov%nsamp, &
         ", worst cell-area error / cell area ", cov%worst_area, ", worst vertex excursion outside its cell ", cov%worst_vert, &
         ", vertices outside their cell ", cov%nflag_other, ", worst first-moment difference / cell area ", cov%worst_mom
   end subroutine test_clipper_noefa

   subroutine test_efa_oracle(error)
      type(error_type), allocatable, intent(out) :: error
      type(cov_t) :: cov
      call sweep_efa(error, cov, 300000)
      if (allocated(error)) return
    write (error_unit, "(a,i0,a,i0,a,es9.2,a,es9.2,a,es9.2)") "  [remap_geom] efa: samples ", cov%nsamp, ", refused (cell exit) ", &
         cov%nflag_cfl, ", worst sum error / cell area ", cov%worst_sum, ", worst cell-area error ", cov%worst_area, &
         ", worst first-moment difference ", cov%worst_mom
      call check(error, cov%nflag_cfl < cov%nsamp/2, "most flux-adjusted samples stay inside their cells")
   end subroutine test_efa_oracle

   subroutine test_coverage(error)
      !! Every branch of the case table was run (union of both sweeps); no flux-adjustment fall-through.
      type(error_type), allocatable, intent(out) :: error
      type(cov_t) :: cov, cn, ce
      integer :: k
      character(len=200) :: msg

      call sweep_noefa(error, cn, 300000)
      if (allocated(error)) return
      call sweep_efa(error, ce, 300000)
      if (allocated(error)) return
      cov%left = cn%left + ce%left
      cov%right = cn%right + ce%right
      cov%centre = cn%centre + ce%centre
      cov%efa = ce%efa
      write (error_unit, "(a)") "  [remap_geom] branch coverage (no-efa + efa samples):"
      write (error_unit, "(a,4(1x,i8))") "     left  TL TL.. : 1 TL, 2 BL, 3 TL1+BL1, 4 TL2+BL2 ->", cov%left(1:4)
      write (error_unit, "(a,4(1x,i8))") "     right TR..    : 1 TR, 2 BR, 3 TR1+BR1, 4 TR2+BR2 ->", cov%right(1:4)
      write (error_unit, "(a,12(1x,i7))") "     centre 1..12 ->", cov%centre(1:12)
      write (error_unit, "(a,4(1x,i8))") "     efa 1 same side, 2 pin ICL, 3 pin ICR, 4 fall-through ->", cov%efa(1:4)
      write (error_unit, "(a,12(1x,i7))") "     centre (no-efa only) ->", cn%centre(1:12)
      do k = 1, 4
         if (cov%left(k) == 0 .or. cov%right(k) == 0) then
            write (msg, "(a,i0)") "side branch never run: ", k
            call test_failed(error, trim(msg))
            return
         end if
      end do
      do k = 1, 12
         if (cov%centre(k) == 0) then
            write (msg, "(a,i0)") "centre branch never run: ", k
            call test_failed(error, trim(msg))
            return
         end if
      end do
      do k = 1, 3
         if (cov%efa(k) == 0) then
            write (msg, "(a,i0)") "flux-adjustment branch never run: ", k
            call test_failed(error, trim(msg))
            return
         end if
      end do
      call check(error, cov%efa(4) == 0, "the flux adjustment never falls through its three branches")
   end subroutine test_coverage

   subroutine test_src_code(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: di, dj
      integer(int8) :: c
      do dj = -1, 1
         do di = -1, 1
            c = int(1 + (di + 1) + 3*(dj + 1), int8)
            call check(error, remap_src_di(c) == di .and. remap_src_dj(c) == dj, "source code round trip")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_src_code

   ! ------------------------------------------------------------------
   ! stdlib Python transcript
   ! ------------------------------------------------------------------

   subroutine test_python_transcript(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 5
      real(wp) :: uc(N, N), vc(N, N), dxU(N, N), dyU(N, N), dpx(N, N), dpy(N, N), uc2(N, N), vc2(N, N)
      integer(int32) :: ierr(N, N)
      real(wp) :: asum, w, err
      real(wp) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT), in(PY_NIN)
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: ic, bits
      integer :: i, j, e, g, k, nused, iedge

      do j = 1, N
         do i = 1, N
            uc(i, j) = 0.30_wp + 0.05_wp*i - 0.02_wp*j + 0.01_wp*i*j
            vc(i, j) = -0.20_wp + 0.03_wp*i + 0.04_wp*j - 0.01_wp*i*j
            dxU(i, j) = 1.0e4_wp*(1.0_wp + 0.05_wp*i)
            dyU(i, j) = 1.0e4_wp*(1.0_wp + 0.04_wp*j)
         end do
      end do
      ierr = 0_int32
      dpx = 0.0_wp
      dpy = 0.0_wp
      uc2 = -uc
      vc2 = -0.7_wp*vc
      !$acc enter data copyin(uc, vc, uc2, vc2, dxU, dyU, dpx, dpy, ierr)
      call remap_departure_points(N, N, 3, 3, 3, 3, 9000.0_wp, uc, vc, dxU, dyU, dpx, dpy, ierr)
      !$acc update self(dpx, dpy, ierr)
      call check_dp(error, dpx(2, 2), dpy(2, 2), 1)
      if (.not. allocated(error)) call check_dp(error, dpx(3, 2), dpy(3, 2), 2)
      if (.not. allocated(error)) call check_dp(error, dpx(2, 3), dpy(2, 3), 3)
      if (.not. allocated(error)) call check_dp(error, dpx(3, 3), dpy(3, 3), 4)
      call remap_departure_points(N, N, 3, 3, 3, 3, 9000.0_wp, uc2, vc2, dxU, dyU, dpx, dpy, ierr)
      !$acc update self(dpx, dpy, ierr)
      if (.not. allocated(error)) call check_dp(error, dpx(2, 2), dpy(2, 2), 5)
      vc2 = -vc2
      !$acc update device(vc2)
      call remap_departure_points(N, N, 3, 3, 3, 3, 9000.0_wp, uc2, vc2, dxU, dyU, dpx, dpy, ierr)
      !$acc update self(dpx, dpy, ierr)
      !$acc exit data delete(uc, vc, uc2, vc2, dxU, dyU, dpx, dpy, ierr)
      if (.not. allocated(error)) call check_dp(error, dpx(3, 2), dpy(3, 2), 6)
      if (allocated(error)) return
      call check(error, all(ierr == 0_int32), "departure points of the transcript clean")
      if (allocated(error)) return

      err = 0.0_wp
      do e = 1, PY_NEDGE
         in = PY_EDGE_IN(:, e)
         iedge = nint(in(2))
         call remap_locate_edge(iedge, nint(in(1)), PUNY, in(3), in(4), in(5), in(6), in(7), in(8), in(9), in(10), &
                                in(11), in(12), in(13), in(14), px, py, ar, sc, ic, asum, bits)
         call check(error, bits == 0_int32, "transcript edge accepted")
         if (allocated(error)) return
         nused = count(sc /= 0_int8)
         call check(error, nused == PY_NUSED(e), "transcript: number of triangles")
         if (allocated(error)) return
         do g = 1, REMAP_NSLOT
            call check(error, (sc(g) /= 0_int8) .eqv. (PY_SLOT(1, g, e) /= 0.0_wp), "transcript: slot in use")
            if (allocated(error)) return
            if (sc(g) == 0_int8) cycle
            call check(error, remap_src_di(sc(g)) == nint(PY_SLOT(2, g, e)) .and. &
                       remap_src_dj(sc(g)) == nint(PY_SLOT(3, g, e)), "transcript: source cell")
            if (allocated(error)) return
            w = abs(ar(g) - PY_SLOT(4, g, e))/in(7)
            err = max(err, w)
            if (.not. (w <= 1.0e-12_wp)) then
               call test_failed(error, "transcript: triangle area")
               return
            end if
            do k = 1, 3
               w = max(abs(px(k, g) - PY_SLOT(3 + 2*k, g, e)), abs(py(k, g) - PY_SLOT(4 + 2*k, g, e)))
               err = max(err, w)
               if (.not. (w <= 1.0e-12_wp)) then
                  call test_failed(error, "transcript: vertex")
                  return
               end if
            end do
         end do
      end do
      write (error_unit, "(a,es9.2)") "  [remap_geom] python transcript: worst edge deviation ", err
   contains
      subroutine check_dp(err_, a, b, q)
         type(error_type), allocatable, intent(inout) :: err_
         real(wp), intent(in) :: a, b
         integer, intent(in) :: q
         real(wp) :: d
         d = max(abs(a - PY_DP(1, q))/abs(PY_DP(1, q)), abs(b - PY_DP(2, q))/abs(PY_DP(2, q)))
         if (.not. (d <= 1.0e-12_wp)) call test_failed(err_, "transcript: departure point")
      end subroutine check_dp
   end subroutine test_python_transcript

   ! ------------------------------------------------------------------
   ! grid helpers
   ! ------------------------------------------------------------------

   subroutine geo_init(g, nxp, nyp, nonuniform)
      type(geo_t), intent(out) :: g
      integer, intent(in) :: nxp, nyp
      logical, intent(in) :: nonuniform
      integer :: i, j
      real(wp), parameter :: DX = 1.0e4_wp, DY = 1.2e4_wp

      g%nx = nxp + 2*NG
      g%ny = nyp + 2*NG
      g%is = NG + 1
      g%ie = NG + nxp
      g%js = NG + 1
      g%je = NG + nyp
      allocate (g%u(g%nx + 1, g%ny), g%v(g%nx, g%ny + 1))
      allocate (g%dxU(g%nx, g%ny), g%dyU(g%nx, g%ny), g%HTE(g%nx, g%ny), g%HTN(g%nx, g%ny))
      allocate (g%earea(g%nx, g%ny), g%narea(g%nx, g%ny), g%ierr(g%nx, g%ny))
      allocate (g%epm(g%nx, g%ny), g%npm(g%nx, g%ny), g%uvm(g%nx, g%ny), g%hm(g%nx, g%ny))
      g%epm = 1.0_wp
      g%npm = 1.0_wp
      g%uvm = 1.0_wp
      g%hm = 1.0_wp
      g%u = 0.0_wp
      g%v = 0.0_wp
      g%ierr = 0_int32
      do j = 1, g%ny
         do i = 1, g%nx
            if (nonuniform) then
               g%dxU(i, j) = DX*(1.0_wp + 0.10_wp*sin(0.37_wp*i + 0.21_wp*j))
               g%dyU(i, j) = DY*(1.0_wp + 0.08_wp*cos(0.29_wp*i - 0.17_wp*j))
               g%HTE(i, j) = DY*(1.0_wp + 0.07_wp*sin(0.31_wp*i + 0.13_wp*j + 1.0_wp))
               g%HTN(i, j) = DX*(1.0_wp + 0.09_wp*cos(0.23_wp*i + 0.27_wp*j + 2.0_wp))
               g%earea(i, j) = DX*DY*(1.0_wp + 0.10_wp*sin(0.19_wp*i + 0.33_wp*j + 0.5_wp))
               g%narea(i, j) = DX*DY*(1.0_wp + 0.10_wp*cos(0.41_wp*i + 0.11_wp*j + 1.5_wp))
            else
               g%dxU(i, j) = DX
               g%dyU(i, j) = DY
               g%HTE(i, j) = DY
               g%HTN(i, j) = DX
               g%earea(i, j) = DX*DY
               g%narea(i, j) = DX*DY
            end if
         end do
      end do
      call g%ws%init(g%nx, g%ny)
   end subroutine geo_init

   subroutine geo_sentinel(g)
      !! Fill every workspace output with the sentinel (host copy).
      type(geo_t), intent(inout) :: g
      g%ws%uc = SENT
      g%ws%vc = SENT
      g%ws%dpx = SENT
      g%ws%dpy = SENT
      g%ws%edgearea_e = SENT
      g%ws%edgearea_n = SENT
      g%ws%xp_e = SENT
      g%ws%yp_e = SENT
      g%ws%area_e = SENT
      g%ws%xp_n = SENT
      g%ws%yp_n = SENT
      g%ws%area_n = SENT
      g%ws%src_e = 77_int8
      g%ws%src_n = 77_int8
      g%ws%icase_e = 7777_int32
      g%ws%icase_n = 7777_int32
   end subroutine geo_sentinel

   subroutine geo_map(g)
      type(geo_t), intent(inout) :: g
      if (g%mapped) return
      !$acc enter data copyin(g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, g%narea, g%epm, g%npm, g%uvm, g%ierr)
      call g%ws%enter_data()
      g%mapped = .true.
   end subroutine geo_map

   subroutine geo_unmap(g)
      type(geo_t), intent(inout) :: g
      if (.not. g%mapped) return
      call g%ws%exit_data()
      !$acc exit data delete(g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, g%narea, g%epm, g%npm, g%uvm, g%ierr)
      g%mapped = .false.
   end subroutine geo_unmap

   subroutine geo_push(g)
      !! Host to device: the inputs and every workspace array the host may have edited.
      type(geo_t), intent(inout) :: g
      !$acc update device(g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, g%narea, g%epm, g%npm, g%uvm, g%ierr)
      !$acc update device(g%ws%uc, g%ws%vc, g%ws%dpx, g%ws%dpy, g%ws%edgearea_e, g%ws%edgearea_n)
      !$acc update device(g%ws%xp_e, g%ws%yp_e, g%ws%area_e, g%ws%src_e, g%ws%icase_e)
      !$acc update device(g%ws%xp_n, g%ws%yp_n, g%ws%area_n, g%ws%src_n, g%ws%icase_n)
   end subroutine geo_push

   subroutine geo_pull(g)
      type(geo_t), intent(inout) :: g
      !$acc update self(g%ierr)
      !$acc update self(g%ws%uc, g%ws%vc, g%ws%dpx, g%ws%dpy, g%ws%edgearea_e, g%ws%edgearea_n)
      !$acc update self(g%ws%xp_e, g%ws%yp_e, g%ws%area_e, g%ws%src_e, g%ws%icase_e)
      !$acc update self(g%ws%xp_n, g%ws%yp_n, g%ws%area_n, g%ws%src_n, g%ws%icase_n)
   end subroutine geo_pull

   subroutine geo_destroy(g)
      type(geo_t), intent(inout) :: g
      call geo_unmap(g)
      call g%ws%destroy()
   end subroutine geo_destroy

   subroutine geo_run(g, efa)
      !! The whole geometry on the (mapped) problem, the host copy refreshed after.
      type(geo_t), intent(inout) :: g
      integer, intent(in) :: efa
      call geo_map(g)
      call geo_push(g)
      call remap_geometry(g%ws, g%is, g%ie, g%js, g%je, g%dt, efa, PUNY, g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, &
                          g%narea, g%epm, g%npm, g%uvm, g%ierr)
      call geo_pull(g)
   end subroutine geo_run

   subroutine geo_locate(g, efa)
      !! Only the triangle stage, from the host-set departure points and edge areas.
      type(geo_t), intent(inout) :: g
      integer, intent(in) :: efa
      call geo_map(g)
      call geo_push(g)
      call remap_locate_east(g%nx, g%ny, g%is, g%ie, g%js, g%je, efa, PUNY, g%ws%dpx, g%ws%dpy, g%dxU, g%dyU, g%earea, &
                             g%narea, g%ws%edgearea_e, g%ws%xp_e, g%ws%yp_e, g%ws%area_e, g%ws%src_e, g%ws%icase_e, g%ierr)
      call remap_locate_north(g%nx, g%ny, g%is, g%ie, g%js, g%je, efa, PUNY, g%ws%dpx, g%ws%dpy, g%dxU, g%dyU, g%earea, &
                              g%narea, g%ws%edgearea_n, g%ws%xp_n, g%ws%yp_n, g%ws%area_n, g%ws%src_n, g%ws%icase_n, g%ierr)
      call geo_pull(g)
   end subroutine geo_locate

   function edge_total(g, iedge, i, j) result(s)
      !! Sum of the slot areas of an edge.
      type(geo_t), intent(in) :: g
      integer, intent(in) :: iedge, i, j
      real(wp) :: s
      if (iedge == REMAP_EDGE_EAST) then
         s = sum(g%ws%area_e(i, j, :))
      else
         s = sum(g%ws%area_n(i, j, :))
      end if
   end function edge_total

   subroutine edge_cell_sums(g, iedge, i, j, sums)
      !! Physical area per source offset `(di, dj)` of an edge.
      type(geo_t), intent(in) :: g
      integer, intent(in) :: iedge, i, j
      real(wp), intent(out) :: sums(-1:1, -1:1)
      integer :: s, di, dj
      integer(int8) :: c
      sums = 0.0_wp
      do s = 1, REMAP_NSLOT
         if (iedge == REMAP_EDGE_EAST) then
            c = g%ws%src_e(i, j, s)
         else
            c = g%ws%src_n(i, j, s)
         end if
         if (c == 0_int8) cycle
         di = remap_src_di(c)
         dj = remap_src_dj(c)
         if (iedge == REMAP_EDGE_EAST) then
            sums(di, dj) = sums(di, dj) + g%ws%area_e(i, j, s)
         else
            sums(di, dj) = sums(di, dj) + g%ws%area_n(i, j, s)
         end if
      end do
   end subroutine edge_cell_sums

   subroutine set_faces_from_fn(g, kind, p1, p2)
      !! Face velocities of a prescribed field on the uniform grid (geometry only: positions
      !! `x = (f-1) dx` at the x-faces, cell centres in y; the field is evaluated at the face).
      type(geo_t), intent(inout) :: g
      integer, intent(in) :: kind
      real(wp), intent(in) :: p1, p2
      integer :: i, j
      real(wp) :: x, y, xc, yc, dx, dy

      dx = g%dxU(1, 1)
      dy = g%dyU(1, 1)
      xc = 0.5_wp*real(g%nx, wp)*dx
      yc = 0.5_wp*real(g%ny, wp)*dy
      do j = 1, g%ny
         do i = 1, g%nx + 1
            x = real(i - 1, wp)*dx
            y = (real(j, wp) - 0.5_wp)*dy
            select case (kind)
            case (1)    ! uniform translation (p1, p2)
               g%u(i, j) = p1
            case (2)    ! solid-body rotation, angular rate p1
               g%u(i, j) = -p1*(y - yc)
            case (3)    ! smooth deformation
               g%u(i, j) = p1*sin(6.2831853_wp*(x/(real(g%nx, wp)*dx) + 0.3_wp)) &
                           *cos(3.1415927_wp*y/(real(g%ny, wp)*dy)) + p2*cos(0.002_wp*x/dx*dy)
            case default
               g%u(i, j) = 0.0_wp
            end select
         end do
      end do
      do j = 1, g%ny + 1
         do i = 1, g%nx
            x = (real(i, wp) - 0.5_wp)*dx
            y = real(j - 1, wp)*dy
            select case (kind)
            case (1)
               g%v(i, j) = p2
            case (2)
               g%v(i, j) = p1*(x - xc)
            case (3)
               g%v(i, j) = -p1*cos(6.2831853_wp*(x/(real(g%nx, wp)*dx) + 0.3_wp)) &
                           *sin(3.1415927_wp*y/(real(g%ny, wp)*dy)) + 0.5_wp*p2*sin(0.003_wp*y/dy*dx)
            case default
               g%v(i, j) = 0.0_wp
            end select
         end do
      end do
   end subroutine set_faces_from_fn

   subroutine set_random_faces(g, amp, seed)
      !! Rough random face velocities in `[-amp, amp]` (every edge, every case).
      type(geo_t), intent(inout) :: g
      real(wp), intent(in) :: amp
      integer(int64), intent(in) :: seed
      integer(int64) :: st
      integer :: i, j
      st = seed
      do j = 1, g%ny
         do i = 1, g%nx + 1
            st = lcg_next(st)
            g%u(i, j) = amp*(2.0_wp*lcg_unit(st) - 1.0_wp)
         end do
      end do
      do j = 1, g%ny + 1
         do i = 1, g%nx
            st = lcg_next(st)
            g%v(i, j) = amp*(2.0_wp*lcg_unit(st) - 1.0_wp)
         end do
      end do
   end subroutine set_random_faces

   function quad_area_phys(x, y) result(a)
      real(wp), intent(in) :: x(4), y(4)
      real(wp) :: a
      a = poly_area(4, x, y)
   end function quad_area_phys

   ! ------------------------------------------------------------------
   ! (4) corner consistency on a uniform metric
   ! ------------------------------------------------------------------

   subroutine test_corner_identity(error)
      !! No-adjustment geometry on a uniform metric: each edge equals its swept quadrilateral, and
      !! the area leaving a cell through its four edges equals the cell area minus the area of its
      !! back-traced image (no gap, no double count where the edges meet at a corner).
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: i, j, mode
      real(wp) :: dx, dy, a, wanted, worst_e, worst_c, out_c, img, tol
      real(wp) :: qx(4), qy(4)
      character(len=160) :: msg

      call geo_init(g, 14, 12, .false.)
      dx = g%dxU(1, 1)
      dy = g%dyU(1, 1)
      worst_e = 0.0_wp
      worst_c = 0.0_wp
      do mode = 1, 3
         g%dt = 1.0_wp
         if (mode == 1) call set_random_faces(g, 0.45_wp*dx, 11_int64)
         if (mode == 2) call set_faces_from_fn(g, 2, 0.05_wp, 0.0_wp)
         if (mode == 3) call set_faces_from_fn(g, 3, 0.35_wp*dx, 0.1_wp*dx)
         call geo_run(g, 0)
         if (any(g%ierr /= 0_int32)) then
            call test_failed(error, "corner identity: an error bit was raised")
            call geo_destroy(g)
            return
         end if
         do j = g%js, g%je
            do i = g%is, g%ie
               ! east edge E(i,j): swept quad U(i,j-1), U(i,j), U(i,j)+dp, U(i,j-1)+dp
               qx = [0.0_wp, 0.0_wp, g%ws%dpx(i, j), g%ws%dpx(i, j - 1)]
               qy = [-dy, 0.0_wp, g%ws%dpy(i, j), -dy + g%ws%dpy(i, j - 1)]
               wanted = quad_area_phys(qx, qy)
               a = edge_total(g, REMAP_EDGE_EAST, i, j)
               worst_e = max(worst_e, abs(a - wanted)/(dx*dy))
               ! north edge N(i,j): U(i,j), U(i-1,j), U(i-1,j)+dp, U(i,j)+dp
               qx = [0.0_wp, -dx, -dx + g%ws%dpx(i - 1, j), g%ws%dpx(i, j)]
               qy = [0.0_wp, 0.0_wp, g%ws%dpy(i - 1, j), g%ws%dpy(i, j)]
               wanted = quad_area_phys(qx, qy)
               a = edge_total(g, REMAP_EDGE_NORTH, i, j)
               worst_e = max(worst_e, abs(a - wanted)/(dx*dy))
               ! cell: outflow = area - area of the image
               out_c = edge_total(g, REMAP_EDGE_EAST, i, j) - edge_total(g, REMAP_EDGE_EAST, i - 1, j) &
                       + edge_total(g, REMAP_EDGE_NORTH, i, j) - edge_total(g, REMAP_EDGE_NORTH, i, j - 1)
               ! image corners: SW U(i-1,j-1), SE U(i,j-1), NE U(i,j), NW U(i-1,j); cell [-dx,0]x[-dy,0] about U(i,j)
               qx = [-dx + g%ws%dpx(i - 1, j - 1), g%ws%dpx(i, j - 1), g%ws%dpx(i, j), -dx + g%ws%dpx(i - 1, j)]
               qy = [-dy + g%ws%dpy(i - 1, j - 1), -dy + g%ws%dpy(i, j - 1), g%ws%dpy(i, j), g%ws%dpy(i - 1, j)]
               img = poly_area(4, qx, qy)
               worst_c = max(worst_c, abs(out_c - (dx*dy - img))/(dx*dy))
            end do
         end do
      end do
      write (error_unit, "(a,es9.2,a,es9.2)") "  [remap_geom] uniform metric: worst edge vs swept quad ", worst_e, &
         ", worst cell outflow vs (area - image) ", worst_c
      tol = 1.0e-13_wp
      if (.not. (worst_e <= tol .and. worst_c <= tol)) then
         write (msg, "(a,2es10.2)") "corner identity violated: ", worst_e, worst_c
         call test_failed(error, trim(msg))
      end if
      call geo_destroy(g)
   end subroutine test_corner_identity

   subroutine test_single_corner(error)
      !! One corner displaced by `(a, b)` cell units (all four quadrants): the four cells around it
      !! have outflow `-(a+b)/2`, `(a-b)/2`, `(a+b)/2`, `(b-a)/2` cell areas (SW, SE, NE, NW), all
      !! other cells and all other edges exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: i, j, ic, jc, q, k
      integer :: ci(4), cj(4)
      real(wp) :: a, b, dx, dy, out_c, worst
      real(wp) :: want(4), ab(2, NQ)
      character(len=160) :: msg

      call geo_init(g, 10, 10, .false.)
      dx = g%dxU(1, 1)
      dy = g%dyU(1, 1)
      ab = reshape([0.4_wp, 0.7_wp, -0.4_wp, 0.7_wp, -0.6_wp, -0.3_wp, 0.5_wp, -0.9_wp, 0.9_wp, 0.0_wp, 0.0_wp, -0.8_wp, &
                    0.2_wp, 0.2_wp, -0.95_wp, 0.6_wp], [2, 8])
      ic = 8
      jc = 7
      worst = 0.0_wp
      do q = 1, 8
         a = ab(1, q)
         b = ab(2, q)
         g%ws%dpx = 0.0_wp
         g%ws%dpy = 0.0_wp
         g%ws%edgearea_e = 0.0_wp
         g%ws%edgearea_n = 0.0_wp
         g%ws%dpx(ic, jc) = a*dx
         g%ws%dpy(ic, jc) = b*dy
         call geo_locate(g, 0)
         if (any(g%ierr /= 0_int32)) then
            call test_failed(error, "single corner: error bit")
            call geo_destroy(g)
            return
         end if
         ! cells (SW, SE, NE, NW) of the corner U(ic, jc)
         ci = [ic, ic + 1, ic + 1, ic]
         cj = [jc, jc, jc + 1, jc + 1]
         want = [-(a + b)/2, (a - b)/2, (a + b)/2, (b - a)/2]
         do k = 1, 4
            i = ci(k)
            j = cj(k)
            out_c = edge_total(g, REMAP_EDGE_EAST, i, j) - edge_total(g, REMAP_EDGE_EAST, i - 1, j) &
                    + edge_total(g, REMAP_EDGE_NORTH, i, j) - edge_total(g, REMAP_EDGE_NORTH, i, j - 1)
            worst = max(worst, abs(out_c/(dx*dy) - want(k)))
         end do
         ! everything else is zero: total of |edge areas| equals the four incident edges only
         do j = g%js - 1, g%je
            do i = g%is - 1, g%ie
               if (i >= g%is - 1 .and. j >= g%js .and. j <= g%je) then
                  if (.not. ((i == ic .and. (j == jc .or. j == jc + 1)))) then
                     if (edge_total(g, REMAP_EDGE_EAST, i, j) /= 0.0_wp) then
                        write (msg, "(a,2i4)") "single corner: stray east area at ", i, j
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                  end if
               end if
               if (i >= g%is .and. i <= g%ie .and. j >= g%js - 1) then
                  if (.not. (j == jc .and. (i == ic .or. i == ic + 1))) then
                     if (edge_total(g, REMAP_EDGE_NORTH, i, j) /= 0.0_wp) then
                        write (msg, "(a,2i4)") "single corner: stray north area at ", i, j
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                  end if
               end if
            end do
         end do
      end do
      write (error_unit, "(a,es9.2)") "  [remap_geom] single moved corner: worst deviation from +-(a+-b)/2 ", worst
      if (.not. (worst <= 1.0e-13_wp)) call test_failed(error, "single moved corner: closed form")
      call geo_destroy(g)
   end subroutine test_single_corner

   ! ------------------------------------------------------------------
   ! (5) translation, rotation, deformation
   ! ------------------------------------------------------------------

   subroutine test_translation(error)
      !! Uniform velocity: every edge's area is `u HTE dt` (resp. `v HTN dt`) both with and without
      !! the adjustment, a pure x translation lies in ONE source cell (the cell it leaves for
      !! eastward flow, the one it enters for westward), and any translation touches at most two.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: efa, i, j, q, nsrc, di, dj
      real(wp) :: a, want, worst
      real(wp) :: uu(2, NUU), sums(-1:1, -1:1)
      character(len=160) :: msg

      call geo_init(g, 12, 10, .false.)
      g%dt = 1.0_wp
      uu = reshape([0.3_wp, 0.0_wp, -0.45_wp, 0.0_wp, 0.0_wp, 0.4_wp, 0.35_wp, -0.25_wp, -0.2_wp, 0.15_wp], [2, 5])
      worst = 0.0_wp
      do efa = 0, 1
         do q = 1, 5
            call set_faces_from_fn(g, 1, uu(1, q)*g%dxU(1, 1), uu(2, q)*g%dyU(1, 1))
            call geo_run(g, efa)
            if (any(g%ierr /= 0_int32)) then
               call test_failed(error, "translation: error bit")
               call geo_destroy(g)
               return
            end if
            do j = g%js, g%je
               do i = g%is - 1, g%ie
                  want = uu(1, q)*g%dxU(1, 1)*g%HTE(i, j)
                  a = edge_total(g, REMAP_EDGE_EAST, i, j)
                  worst = max(worst, abs(a - want)/(g%dxU(1, 1)*g%dyU(1, 1)))
                  call edge_cell_sums(g, REMAP_EDGE_EAST, i, j, sums)
                  nsrc = count(sums /= 0.0_wp)
                  if (nsrc > 2) then
                     call test_failed(error, "translation: east edge draws on more than two cells")
                     call geo_destroy(g)
                     return
                  end if
                  if (uu(2, q) == 0.0_wp .and. uu(1, q) /= 0.0_wp) then
                     di = merge(0, 1, uu(1, q) > 0.0_wp)
                     dj = 0
                     if (nsrc /= 1 .or. sums(di, dj) == 0.0_wp) then
                        write (msg, "(a,2i4)") "translation: x-translation not confined to its source cell, edge ", i, j
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                  end if
               end do
            end do
            do j = g%js - 1, g%je
               do i = g%is, g%ie
                  want = uu(2, q)*g%dyU(1, 1)*g%HTN(i, j)
                  a = edge_total(g, REMAP_EDGE_NORTH, i, j)
                  worst = max(worst, abs(a - want)/(g%dxU(1, 1)*g%dyU(1, 1)))
                  call edge_cell_sums(g, REMAP_EDGE_NORTH, i, j, sums)
                  if (count(sums /= 0.0_wp) > 2) then
                     call test_failed(error, "translation: north edge draws on more than two cells")
                     call geo_destroy(g)
                     return
                  end if
               end do
            end do
         end do
      end do
      write (error_unit, "(a,es9.2)") "  [remap_geom] translation: worst edge area error / cell area ", worst
      if (.not. (worst <= 1.0e-14_wp)) call test_failed(error, "translation: edge area is not u HTE dt")
      call geo_destroy(g)
   end subroutine test_translation

   subroutine test_rotation(error)
      !! Solid-body rotation by `omega dt = 0.05` rad: the corner identity holds exactly, the net outflow
      !! of every cell (the image of a rotated cell has the cell's area) is below the third-order
      !! bound of the midpoint trajectories, and with the adjustment on every edge sums to its face flux.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: i, j
      real(wp) :: dx, dy, out_c, worst, worst_sum, wsum
      character(len=160) :: msg

      call geo_init(g, 16, 14, .false.)
      dx = g%dxU(1, 1)
      dy = g%dyU(1, 1)
      g%dt = 1.0_wp
      call set_faces_from_fn(g, 2, 0.05_wp, 0.0_wp)
      ! u = -w (y - yc), v = w (x - xc) with w dt = 0.05: the displacement w dt |r| stays below one cell
      call geo_run(g, 0)
      if (any(g%ierr /= 0_int32)) then
         call test_failed(error, "rotation: error bit")
         call geo_destroy(g)
         return
      end if
      worst = 0.0_wp
      do j = g%js, g%je
         do i = g%is, g%ie
            out_c = edge_total(g, REMAP_EDGE_EAST, i, j) - edge_total(g, REMAP_EDGE_EAST, i - 1, j) &
                    + edge_total(g, REMAP_EDGE_NORTH, i, j) - edge_total(g, REMAP_EDGE_NORTH, i, j - 1)
            worst = max(worst, abs(out_c)/(dx*dy))
         end do
      end do
      call geo_run(g, 1)
      worst_sum = 0.0_wp
      do j = g%js, g%je
         do i = g%is - 1, g%ie
            wsum = max(0.0_wp, abs(edge_total(g, REMAP_EDGE_EAST, i, j) - g%ws%edgearea_e(i, j))/(dx*dy))
            worst_sum = max(worst_sum, wsum)
         end do
      end do
      do j = g%js - 1, g%je
         do i = g%is, g%ie
            wsum = abs(edge_total(g, REMAP_EDGE_NORTH, i, j) - g%ws%edgearea_n(i, j))/(dx*dy)
            worst_sum = max(worst_sum, wsum)
         end do
      end do
      write (error_unit, "(a,es9.2,a,es9.2)") "  [remap_geom] rotation: worst |net outflow|/area ", worst, &
         ", efa worst |sum - edgearea|/area ", worst_sum
      if (any(g%ierr /= 0_int32)) call test_failed(error, "rotation (efa): error bit")
      if (allocated(error)) then
         call geo_destroy(g)
         return
      end if
      if (.not. (worst <= 2.0e-3_wp)) then
         write (msg, "(a,es10.2)") "rotation: net outflow of a rigidly rotated cell too large ", worst
         call test_failed(error, trim(msg))
      else if (.not. (worst_sum <= 1.0e-13_wp)) then
         call test_failed(error, "rotation (efa): sum(triarea) /= edgearea")
      end if
      call geo_destroy(g)
   end subroutine test_rotation

   subroutine test_deformation(error)
      !! A smooth deformation field on a NON-uniform metric, both modes: no error bit; with the
      !! adjustment each edge sums to its face flux (round-off); and every triangle's area is
      !! (an area factor the table documents for that source cell) x (area of its vertices in the
      !! source frame): the corner pieces carry the neighbouring edge's cell area.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: efa, i, j, s, di, dj, l, iedge, k
      real(wp) :: fc, cross, rat, err, worst_sum, worst_fac, ea, ar, cx, cy
      real(wp) :: cand(3)
      integer(int8) :: c
      logical :: found
      character(len=160) :: msg

      call geo_init(g, 14, 12, .true.)
      g%dt = 1.0_wp
      worst_sum = 0.0_wp
      worst_fac = 0.0_wp
      do efa = 0, 1
         call set_faces_from_fn(g, 3, 0.04_wp*g%dxU(1, 1), 0.02_wp*g%dxU(1, 1))
         call geo_run(g, efa)
         if (any(g%ierr /= 0_int32)) then
            write (msg, "(a,i0,a,i0)") "deformation: error bit ", maxval(g%ierr), " efa ", efa
            call test_failed(error, trim(msg))
            call geo_destroy(g)
            return
         end if
         do iedge = REMAP_EDGE_NORTH, REMAP_EDGE_EAST
            do j = merge(g%js - 1, g%js, iedge == REMAP_EDGE_NORTH), g%je
               do i = merge(g%is, g%is - 1, iedge == REMAP_EDGE_NORTH), g%ie
                  if (iedge == REMAP_EDGE_NORTH) then
                     fc = g%narea(i, j)
                     ea = g%ws%edgearea_n(i, j)
                  else
                     fc = g%earea(i, j)
                     ea = g%ws%edgearea_e(i, j)
                  end if
                  if (efa == 1) then
                     worst_sum = max(worst_sum, abs(edge_total(g, iedge, i, j) - ea)/fc)
                  end if
                  do s = 1, REMAP_NSLOT
                     if (iedge == REMAP_EDGE_NORTH) then
                        c = g%ws%src_n(i, j, s)
                        ar = g%ws%area_n(i, j, s)
                     else
                        c = g%ws%src_e(i, j, s)
                        ar = g%ws%area_e(i, j, s)
                     end if
                     if (c == 0_int8) cycle
                     di = remap_src_di(c)
                     dj = remap_src_dj(c)
                     l = label_of(iedge, di, dj)
                     if (iedge == REMAP_EDGE_NORTH) then
                        cross = vcross(g%ws%xp_n(i, j, :, s), g%ws%yp_n(i, j, :, s))
                        select case (l)
                        case (1)
                           cand = [g%earea(i - 1, j + 1), g%narea(i - 1, j), fc]
                        case (2)
                           cand = [g%earea(i, j + 1), g%narea(i + 1, j), fc]
                        case (3)
                           cand = [g%earea(i - 1, j), g%narea(i - 1, j), fc]
                        case (4)
                           cand = [g%earea(i, j), g%narea(i + 1, j), fc]
                        case default
                           cand = [fc, fc, fc]
                        end select
                     else
                        cross = vcross(g%ws%xp_e(i, j, :, s), g%ws%yp_e(i, j, :, s))
                        select case (l)
                        case (1)
                           cand = [g%narea(i + 1, j), g%earea(i, j + 1), fc]
                        case (2)
                           cand = [g%narea(i + 1, j - 1), g%earea(i, j - 1), fc]
                        case (3)
                           cand = [g%narea(i, j), g%earea(i, j + 1), fc]
                        case (4)
                           cand = [g%narea(i, j - 1), g%earea(i, j - 1), fc]
                        case default
                           cand = [fc, fc, fc]
                        end select
                     end if
                     rat = abs(ar)
                     err = huge(1.0_wp)
                     found = .false.
                     do k = 1, 3
                        err = min(err, abs(rat - cand(k)*abs(cross)))
                     end do
                     worst_fac = max(worst_fac, err/fc)
                     if (.not. (err <= 1.0e-13_wp*fc)) then
                        write (msg, "(a,i0,a,2i4,a,es10.2)") "deformation: area factor of slot ", s, " at ", i, j, &
                           " is not a documented neighbour's, off by ", err/fc
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                     cx = 0.0_wp
                     cy = cx
                  end do
               end do
            end do
         end do
      end do
   write (error_unit, "(a,es9.2,a,es9.2)") "  [remap_geom] deformation (non-uniform metric): efa worst |sum - edgearea|/areafac ", &
         worst_sum, ", worst area-factor deviation ", worst_fac
      if (.not. (worst_sum <= 1.0e-13_wp)) call test_failed(error, "deformation: sum(triarea) /= edgearea")
      call geo_destroy(g)
   end subroutine test_deformation

   function vcross(x, y) result(c)
      real(wp), intent(in) :: x(3), y(3)
      real(wp) :: c
      c = 0.5_wp*((x(2) - x(1))*(y(3) - y(1)) - (y(2) - y(1))*(x(3) - x(1)))
   end function vcross

   ! ------------------------------------------------------------------
   ! (6) channel, zero, CFL, NaN, mirror
   ! ------------------------------------------------------------------

   subroutine test_channel(error)
      !! One-cell-wide channel: departure points zero, edge area non-zero: ONE triangle of exactly
      !! the prescribed area, in the cell the fluid leaves (`BC`) for eastward flow, in the one it
      !! enters (`TC`) for westward; zero edge area and zero displacement is inactive.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: q, i, j, nused, s
      real(wp) :: fc
      real(wp) :: ea(3), sums(-1:1, -1:1)
      integer, parameter :: ci = 7, cj = 6

      call geo_init(g, 10, 10, .true.)
      ea = [0.3e8_wp, -0.25e8_wp, 0.0_wp]
      do q = 1, 3
         g%ws%dpx = 0.0_wp
         g%ws%dpy = 0.0_wp
         g%ws%edgearea_e = 0.0_wp
         g%ws%edgearea_n = 0.0_wp
         g%ws%edgearea_e(ci, cj) = ea(q)
         g%ws%edgearea_n(ci, cj) = -ea(q)
         call geo_locate(g, 1)
         if (any(g%ierr /= 0_int32)) then
            call test_failed(error, "channel: error bit")
            call geo_destroy(g)
            return
         end if
         do j = g%js, g%je
            do i = g%is - 1, g%ie
               nused = count(g%ws%src_e(i, j, :) /= 0_int8)
               if (i == ci .and. j == cj .and. ea(q) /= 0.0_wp) then
                  fc = g%earea(i, j)
                  call edge_cell_sums(g, REMAP_EDGE_EAST, i, j, sums)
                  if (nused < 1 .or. .not. (abs(edge_total(g, REMAP_EDGE_EAST, i, j) - ea(q)) <= 1.0e-13_wp*fc)) then
                     call test_failed(error, "channel: east area is not the prescribed one")
                     call geo_destroy(g)
                     return
                  end if
                  if (ea(q) > 0.0_wp .and. sums(0, 0) == 0.0_wp) then
                     call test_failed(error, "channel: eastward flow must come from the edge's own cell")
                     call geo_destroy(g)
                     return
                  end if
                  if (ea(q) < 0.0_wp .and. sums(1, 0) == 0.0_wp) then
                     call test_failed(error, "channel: westward flow must come from the cell across the edge")
                     call geo_destroy(g)
                     return
                  end if
                  if (count(sums /= 0.0_wp) /= 1) then
                     call test_failed(error, "channel: more than one source cell")
                     call geo_destroy(g)
                     return
                  end if
                  do s = 1, REMAP_NSLOT
                     if (g%ws%src_e(i, j, s) /= 0_int8 .and. g%ws%icase_e(i, j) <= 0) then
                        call test_failed(error, "channel: case code")
                        call geo_destroy(g)
                        return
                     end if
                  end do
               else
                  if (nused /= 0 .or. g%ws%icase_e(i, j) /= 0_int32) then
                     call test_failed(error, "channel: a stray triangle on an inactive east edge")
                     call geo_destroy(g)
                     return
                  end if
               end if
            end do
         end do
         if (ea(q) /= 0.0_wp) then
            call edge_cell_sums(g, REMAP_EDGE_NORTH, ci, cj, sums)
            if (.not. (abs(edge_total(g, REMAP_EDGE_NORTH, ci, cj) + ea(q)) <= 1.0e-13_wp*g%narea(ci, cj))) then
               call test_failed(error, "channel: north area is not the prescribed one")
               call geo_destroy(g)
               return
            end if
         end if
      end do
      call geo_destroy(g)
   end subroutine test_channel

   subroutine test_zero(error)
      !! Zero velocity: every area is exactly zero, no slot in use, no case code, zero displacement,
      !! no error bit, with and without the adjustment, also after a non-zero step (stale data cleared).
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: efa, i, j, s

      call geo_init(g, 12, 10, .true.)
      g%dt = 1.0_wp
      do efa = 0, 1
         call set_random_faces(g, 0.1_wp*g%dxU(1, 1), 5_int64)
         call geo_run(g, efa)
         if (count(g%ws%src_e(g%is:g%ie, g%js:g%je, :) /= 0_int8) == 0) then
            call test_failed(error, "zero: the non-zero step produced no triangles")
            call geo_destroy(g)
            return
         end if
         g%u = 0.0_wp
         g%v = 0.0_wp
         call geo_run(g, efa)
         if (any(g%ierr /= 0_int32)) then
            call test_failed(error, "zero: error bit")
            call geo_destroy(g)
            return
         end if
         do j = g%js - 1, g%je
            do i = g%is - 1, g%ie
               if (j >= g%js) then
                  if (g%ws%icase_e(i, j) /= 0_int32 .or. any(g%ws%src_e(i, j, :) /= 0_int8) .or. &
                      any(g%ws%area_e(i, j, :) /= 0.0_wp) .or. g%ws%edgearea_e(i, j) /= 0.0_wp) then
                     call test_failed(error, "zero: east edge not exactly empty")
                     call geo_destroy(g)
                     return
                  end if
               end if
               if (i >= g%is) then
                  if (g%ws%icase_n(i, j) /= 0_int32 .or. any(g%ws%src_n(i, j, :) /= 0_int8) .or. &
                      any(g%ws%area_n(i, j, :) /= 0.0_wp) .or. g%ws%edgearea_n(i, j) /= 0.0_wp) then
                     call test_failed(error, "zero: north edge not exactly empty")
                     call geo_destroy(g)
                     return
                  end if
               end if
               if (g%ws%dpx(i, j) /= 0.0_wp .or. g%ws%dpy(i, j) /= 0.0_wp) then
                  call test_failed(error, "zero: displacement not exactly zero")
                  call geo_destroy(g)
                  return
               end if
            end do
         end do
         s = 0
      end do
      call geo_destroy(g)
   end subroutine test_zero

   subroutine test_cfl(error)
      !! A corner displaced by more than a cell raises the bit at that corner (+-x, +-y, both
      !! components); one just below passes; departure points that cross raise it on the edge.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: q, ic, jc, i, j, n
      real(wp) :: asum
      real(wp) :: fac(2, NQ), px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: icase, bits
      logical :: want_bit

      call geo_init(g, 10, 10, .true.)
      g%dt = 1.0_wp
      fac = reshape([1.02_wp, 0.0_wp, -1.02_wp, 0.0_wp, 0.0_wp, 1.02_wp, 0.0_wp, -1.02_wp, 0.98_wp, 0.0_wp, -0.98_wp, 0.0_wp, &
                     0.0_wp, 0.98_wp, 0.0_wp, -0.98_wp], [2, 8])
      ic = 8
      jc = 7
      do q = 1, 8
         ! corner velocities set directly: the departure kernel is what is under test
         g%ws%uc = 0.0_wp
         g%ws%vc = 0.0_wp
         g%ierr = 0_int32
         g%ws%uc(ic, jc) = -fac(1, q)*g%dxU(ic, jc)
         g%ws%vc(ic, jc) = -fac(2, q)*g%dyU(ic, jc)
         call geo_map(g)
         call geo_push(g)
         call remap_departure_points(g%nx, g%ny, g%is, g%ie, g%js, g%je, g%dt, g%ws%uc, g%ws%vc, g%dxU, g%dyU, &
                                     g%ws%dpx, g%ws%dpy, g%ierr)
         call geo_pull(g)
         want_bit = maxval(abs(fac(:, q))) > 1.0_wp
         n = count(g%ierr /= 0_int32)
         if (want_bit) then
            if (n /= 1 .or. iand(g%ierr(ic, jc), ICE_ERR_TRANSPORT_CFL) == 0_int32) then
               call test_failed(error, "cfl: a displacement beyond one cell must raise the bit at its corner only")
               call geo_destroy(g)
               return
            end if
         else
            if (n /= 0) then
               call test_failed(error, "cfl: a displacement just below one cell must pass")
               call geo_destroy(g)
               return
            end if
         end if
      end do
      do j = 1, g%ny
         do i = 1, g%nx
            if (iand(g%ierr(i, j), not(ICE_ERR_TRANSPORT_CFL)) /= 0_int32) then
               call test_failed(error, "cfl: only the CFL bit may be raised")
               call geo_destroy(g)
               return
            end if
         end do
      end do
      ! crossing departure points and a displacement beyond one cell at the edge routine
      call remap_locate_edge(REMAP_EDGE_NORTH, 1, PUNY, 0.6_wp, 0.1_wp, -0.6_wp, 0.1_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, &
                             1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 0.0_wp, px, py, ar, sc, icase, asum, bits)
      if (bits /= ICE_ERR_TRANSPORT_CFL) call test_failed(error, "cfl: crossing departure points")
      call remap_locate_edge(REMAP_EDGE_EAST, 0, PUNY, 0.2_wp, 1.2_wp, 0.1_wp, 0.1_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, &
                             1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 0.0_wp, px, py, ar, sc, icase, asum, bits)
      if (bits /= ICE_ERR_TRANSPORT_CFL) call test_failed(error, "cfl: displacement beyond one cell at the edge")
      call geo_destroy(g)
   end subroutine test_cfl

   subroutine test_nan(error)
      !! A NaN face velocity: the corners that average it raise the bit with a zero (finite)
      !! displacement; an edge with a NaN area is refused with the bit and no triangle; no output is NaN.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: i, j, s
      real(wp) :: nan, asum
      real(wp) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: icase, bits
      integer, parameter :: fi = 8, fj = 7

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call geo_init(g, 12, 10, .false.)
      g%dt = 1.0_wp
      call set_random_faces(g, 0.1_wp*g%dxU(1, 1), 3_int64)
      g%u(fi, fj) = nan
      call geo_run(g, 1)
      ! u(fi, fj) is the east face of cell (fi-1, fj): corners U(fi-1, fj) and U(fi-1, fj-1) average it
    if (iand(g%ierr(fi - 1, fj), ICE_ERR_NONFINITE) == 0_int32 .or. iand(g%ierr(fi - 1, fj - 1), ICE_ERR_NONFINITE) == 0_int32) then
         call test_failed(error, "nan: the corners that average the face must raise the bit")
         call geo_destroy(g)
         return
      end if
      if (g%ws%dpx(fi - 1, fj) /= 0.0_wp .or. g%ws%dpy(fi - 1, fj) /= 0.0_wp) then
         call test_failed(error, "nan: displacement of a refused corner must be zero")
         call geo_destroy(g)
         return
      end if
      ! the east edge of that face has a NaN prescribed area: refused
      if (iand(g%ierr(fi - 1, fj), ICE_ERR_NONFINITE) == 0_int32 .or. any(g%ws%src_e(fi - 1, fj, :) /= 0_int8)) then
         call test_failed(error, "nan: an edge with a NaN prescribed area carries no triangle")
         call geo_destroy(g)
         return
      end if
      do j = g%js - 1, g%je
         do i = g%is - 1, g%ie
            if (.not. (ieee_is_finite(g%ws%dpx(i, j)) .and. ieee_is_finite(g%ws%dpy(i, j)))) then
               call test_failed(error, "nan: a displacement is non-finite")
               call geo_destroy(g)
               return
            end if
            do s = 1, REMAP_NSLOT
               if (.not. ieee_is_finite(g%ws%area_e(i, j, s)) .or. .not. ieee_is_finite(g%ws%area_n(i, j, s))) then
                  call test_failed(error, "nan: a triangle area is non-finite")
                  call geo_destroy(g)
                  return
               end if
            end do
         end do
      end do
      call remap_locate_edge(REMAP_EDGE_NORTH, 0, PUNY, nan, 0.1_wp, 0.1_wp, 0.1_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 1.0e8_wp, &
                             1.0e8_wp, 1.0e8_wp, 1.0e8_wp, 0.0_wp, px, py, ar, sc, icase, asum, bits)
      call check(error, bits == ICE_ERR_NONFINITE, "nan: a NaN displacement at the edge routine")
      call geo_destroy(g)
   end subroutine test_nan

   subroutine test_mirror(error)
      !! Reflecting the velocity field in x (resp. y) reflects every edge's area per source cell,
      !! with the sign the flux convention gives (east/north edges: see the comments), to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g, h
      integer :: efa, i, j, ip, jp, di, dj, nxx, nyy, dir
      real(wp) :: worst, dx, dy, sgn, d
      real(wp) :: sa(-1:1, -1:1), sb(-1:1, -1:1)
      character(len=160) :: msg

      call geo_init(g, 12, 10, .false.)
      call geo_init(h, 12, 10, .false.)
      g%dt = 1.0_wp
      h%dt = 1.0_wp
      nxx = g%nx
      nyy = g%ny
      dx = g%dxU(1, 1)
      dy = g%dyU(1, 1)
      worst = 0.0_wp
      do efa = 0, 1
         call set_random_faces(g, 0.4_wp*dx, 99_int64)
         call geo_run(g, efa)
         do dir = 1, 2
            if (dir == 1) then
               do j = 1, nyy
                  do i = 1, nxx + 1
                     h%u(i, j) = -g%u(nxx + 2 - i, j)
                  end do
               end do
               do j = 1, nyy + 1
                  do i = 1, nxx
                     h%v(i, j) = g%v(nxx + 1 - i, j)
                  end do
               end do
            else
               do j = 1, nyy
                  do i = 1, nxx + 1
                     h%u(i, j) = g%u(i, nyy + 1 - j)
                  end do
               end do
               do j = 1, nyy + 1
                  do i = 1, nxx
                     h%v(i, j) = -g%v(i, nyy + 2 - j)
                  end do
               end do
            end if
            call geo_run(h, efa)
            if (any(g%ierr /= 0_int32) .or. any(h%ierr /= 0_int32)) then
               call test_failed(error, "mirror: error bit")
               call geo_destroy(g)
               call geo_destroy(h)
               return
            end if
            do j = g%js, g%je
               do i = g%is - 1, g%ie
                  ! east edge E(i,j)
                  if (dir == 1) then
                     ip = nxx - i
                     jp = j
                  else
                     ip = i
                     jp = nyy + 1 - j
                  end if
                  call edge_cell_sums(g, REMAP_EDGE_EAST, i, j, sa)
                  call edge_cell_sums(h, REMAP_EDGE_EAST, ip, jp, sb)
                  do dj = -1, 1
                     do di = -1, 1
                        if (dir == 1) then
                           d = abs(sb(1 - di, dj) + sa(di, dj))
                        else
                           d = abs(sb(di, -dj) - sa(di, dj))
                        end if
                        if (di < 0 .or. dj > 1) cycle
                        worst = max(worst, d/(dx*dy))
                     end do
                  end do
               end do
            end do
            do j = g%js - 1, g%je
               do i = g%is, g%ie
                  if (dir == 1) then
                     ip = nxx + 1 - i
                     jp = j
                  else
                     ip = i
                     jp = nyy - j
                  end if
                  call edge_cell_sums(g, REMAP_EDGE_NORTH, i, j, sa)
                  call edge_cell_sums(h, REMAP_EDGE_NORTH, ip, jp, sb)
                  do dj = -1, 1
                     do di = -1, 1
                        if (dir == 1) then
                           d = abs(sb(-di, dj) - sa(di, dj))
                        else
                           d = abs(sb(di, 1 - dj) + sa(di, dj))
                        end if
                        if (dj < 0 .or. dj > 1) cycle
                        worst = max(worst, d/(dx*dy))
                     end do
                  end do
               end do
            end do
         end do
      end do
      write (error_unit, "(a,es9.2)") "  [remap_geom] mirror symmetry x and y: worst per-source-cell deviation / cell area ", worst
      sgn = 1.0e-13_wp
      if (.not. (worst <= sgn)) then
         write (msg, "(a,es10.2)") "mirror symmetry broken: ", worst
         call test_failed(error, trim(msg))
      end if
      call geo_destroy(g)
      call geo_destroy(h)
   end subroutine test_mirror

   subroutine poison_outside(a, i0, i1, j0, j1)
      !! NaN everywhere in a 2-D array outside the index box `[i0:i1, j0:j1]`.
      real(wp), intent(inout) :: a(:, :)
      integer, intent(in) :: i0, i1, j0, j1
      integer :: i, j
      real(wp) :: nan
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (i < i0 .or. i > i1 .or. j < j0 .or. j > j1) a(i, j) = nan
         end do
      end do
   end subroutine poison_outside

   subroutine test_halo_window(error)
      !! The stated input windows are enough: with NaN everywhere outside them (faces `u, v` on `(is-2:ie+2,
      !! js-2:je+2)`, masks on `(is-2:ie+1, js-2:je+1)`, `dxU, dyU` on `is-1:ie, js-1:je`, `HTE` on `is-1:ie, js:je`,
      !! `HTN` on `is:ie, js-1:je`, `earea` on `is-3:ie+1, js-2:je+2`, `narea` on `is-2:ie+2, js-3:je+1`) both modes give outputs
      !! bit-identical to the run on clean data and raise no bit.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g, h
      integer :: efa, i, j, s
      logical :: same

      call geo_init(g, 12, 10, .true.)
      call geo_init(h, 12, 10, .true.)
      g%dt = 1.0_wp
      h%dt = 1.0_wp
      do efa = 0, 1
         call set_random_faces(g, 0.4_wp*g%dxU(1, 1), 31_int64)
         call geo_run(g, efa)
         h%u = g%u
         h%v = g%v
         h%dxU = g%dxU
         h%dyU = g%dyU
         h%HTE = g%HTE
         h%HTN = g%HTN
         h%earea = g%earea
         h%narea = g%narea
         call poison_outside(h%u, g%is - 2, g%ie + 2, g%js - 2, g%je + 2)
         call poison_outside(h%v, g%is - 2, g%ie + 2, g%js - 2, g%je + 2)
         h%epm = g%epm
         h%npm = g%npm
         h%uvm = g%uvm
         call poison_outside(h%epm, g%is - 2, g%ie + 1, g%js - 2, g%je + 1)
         call poison_outside(h%npm, g%is - 2, g%ie + 1, g%js - 2, g%je + 1)
         call poison_outside(h%uvm, g%is - 2, g%ie + 1, g%js - 2, g%je + 1)
         call poison_outside(h%dxU, g%is - 1, g%ie, g%js - 1, g%je)
         call poison_outside(h%dyU, g%is - 1, g%ie, g%js - 1, g%je)
         call poison_outside(h%HTE, g%is - 1, g%ie, g%js, g%je)
         call poison_outside(h%HTN, g%is, g%ie, g%js - 1, g%je)
         call poison_outside(h%earea, g%is - 3, g%ie + 1, g%js - 2, g%je + 2)
         call poison_outside(h%narea, g%is - 2, g%ie + 2, g%js - 3, g%je + 1)
         h%ierr = 0_int32
         call geo_sentinel(h)
         call geo_run(h, efa)
         if (any(h%ierr /= 0_int32) .or. any(g%ierr /= 0_int32)) then
            call test_failed(error, "halo window: a NaN outside the stated windows was read")
            call geo_destroy(g)
            call geo_destroy(h)
            return
         end if
         same = .true.
         do j = g%js - 1, g%je
            do i = g%is - 1, g%ie
               if (j >= g%js) then
                  same = same .and. all(g%ws%src_e(i, j, :) == h%ws%src_e(i, j, :)) .and. &
                         all(g%ws%area_e(i, j, :) == h%ws%area_e(i, j, :)) .and. g%ws%icase_e(i, j) == h%ws%icase_e(i, j)
               end if
               if (i >= g%is) then
                  same = same .and. all(g%ws%src_n(i, j, :) == h%ws%src_n(i, j, :)) .and. &
                         all(g%ws%area_n(i, j, :) == h%ws%area_n(i, j, :)) .and. g%ws%icase_n(i, j) == h%ws%icase_n(i, j)
               end if
               same = same .and. g%ws%dpx(i, j) == h%ws%dpx(i, j) .and. g%ws%dpy(i, j) == h%ws%dpy(i, j)
            end do
         end do
         if (.not. same) then
            call test_failed(error, "halo window: results depend on data outside the stated windows")
            call geo_destroy(g)
            call geo_destroy(h)
            return
         end if
         s = 0
      end do
      call geo_destroy(g)
      call geo_destroy(h)
   end subroutine test_halo_window

   subroutine set_land(g)
      !! West coast (`i <= is`) and a 2 x 2 island: `hm = 0` there, the face masks `epm, npm` and the corner
      !! mask `uvm` (all four T cells wet) as `rdb_ice_cgrid_metrics` defines them.
      type(geo_t), intent(inout) :: g
      integer :: i, j
      g%hm = 1.0_wp
      do j = 1, g%ny
         do i = 1, g%nx
            if (i <= g%is) g%hm(i, j) = 0.0_wp
            if (i >= g%is + 6 .and. i <= g%is + 7 .and. j >= g%js + 4 .and. j <= g%js + 5) g%hm(i, j) = 0.0_wp
         end do
      end do
      g%epm = 0.0_wp
      g%npm = 0.0_wp
      g%uvm = 0.0_wp
      do j = 1, g%ny - 1
         do i = 1, g%nx - 1
            g%epm(i, j) = g%hm(i, j)*g%hm(i + 1, j)
            g%npm(i, j) = g%hm(i, j)*g%hm(i, j + 1)
            g%uvm(i, j) = g%hm(i, j)*g%hm(i + 1, j)*g%hm(i, j + 1)*g%hm(i + 1, j + 1)
         end do
      end do
   end subroutine set_land

   subroutine test_coast(error)
      !! A coast and an island in a uniform flow (open faces carry the flow, closed faces zero): with
      !! the dynamics' masked corner velocities a closed face has NO triangle and zero area, no edge
      !! with a triangle touches a land cell, no source cell is land; and the same run with the corner
      !! mask switched off DOES put triangles on closed faces (the test has teeth); corner velocities
      !! supplied by the caller give the identical result.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      real(wp), allocatable :: ucin(:, :), vcin(:, :), area_e0(:, :, :), area_n0(:, :, :)
      integer :: efa, pass, i, j, s, di, dj, nopen, nbad
      integer(int8) :: c

      call geo_init(g, 14, 12, .true.)
      g%dt = 1.0_wp
      call set_land(g)
      allocate (ucin(g%nx, g%ny), vcin(g%nx, g%ny), area_e0(g%nx, g%ny, REMAP_NSLOT), area_n0(g%nx, g%ny, REMAP_NSLOT))
      do efa = 0, 1
         do j = 1, g%ny
            do i = 1, g%nx + 1
               g%u(i, j) = 0.0_wp
               if (i >= 2 .and. i <= g%nx) g%u(i, j) = 0.30_wp*g%dxU(1, 1)*g%epm(i - 1, j)
            end do
         end do
         do j = 1, g%ny + 1
            do i = 1, g%nx
               g%v(i, j) = 0.0_wp
               if (j >= 2 .and. j <= g%ny) g%v(i, j) = 0.20_wp*g%dxU(1, 1)*g%npm(i, j - 1)
            end do
         end do
         do pass = 1, 3
            if (pass == 2) g%uvm = 1.0_wp
            if (pass == 3) call set_land(g)
            call geo_run(g, efa)
            nopen = 0
            nbad = 0
            do j = g%js - 1, g%je
               do i = g%is - 1, g%ie
                  if (j >= g%js) then
                     call check_edge(REMAP_EDGE_EAST, i, j, g%epm(i, j), nopen, nbad)
                  end if
                  if (i >= g%is) then
                     call check_edge(REMAP_EDGE_NORTH, i, j, g%npm(i, j), nopen, nbad)
                  end if
               end do
            end do
            if (pass == 2) then
           call check(error, nbad > 0, "coast: with the corner mask off a closed face or land cell gets triangles (test has teeth)")
            else
               call check(error, nbad == 0 .and. nopen > 0 .and. all(g%ierr == 0_int32), &
                          "coast: closed faces carry no triangle, land is never a source or a destination")
            end if
            if (allocated(error)) then
               call geo_destroy(g)
               return
            end if
            if (pass == 1) then
               ucin = g%ws%uc
               vcin = g%ws%vc
               area_e0 = g%ws%area_e
               area_n0 = g%ws%area_n
            end if
         end do
         ! caller-supplied corner velocities (the dynamics' own uvelU/vvelU) give the identical geometry
         call geo_map(g)
         !$acc enter data copyin(ucin, vcin)
         call geo_push(g)
         call remap_geometry(g%ws, g%is, g%ie, g%js, g%je, g%dt, efa, PUNY, g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, &
                             g%narea, g%epm, g%npm, g%uvm, g%ierr, ucin, vcin)
         !$acc exit data delete(ucin, vcin)
         call geo_pull(g)
         do s = 1, REMAP_NSLOT
            if (any(g%ws%area_e(g%is - 1:g%ie, g%js:g%je, s) /= area_e0(g%is - 1:g%ie, g%js:g%je, s)) .or. &
                any(g%ws%area_n(g%is:g%ie, g%js - 1:g%je, s) /= area_n0(g%is:g%ie, g%js - 1:g%je, s))) then
               call test_failed(error, "coast: supplied corner velocities give a different geometry")
               call geo_destroy(g)
               return
            end if
         end do
      end do
      call geo_destroy(g)
   contains
      subroutine check_edge(iedge, i, j, open_mask, nopen, nbad)
         integer, intent(in) :: iedge, i, j
         real(wp), intent(in) :: open_mask
         integer, intent(inout) :: nopen, nbad
         integer :: s, di, dj, ic, jc
         logical :: any_used
         any_used = .false.
         do s = 1, REMAP_NSLOT
            if (iedge == REMAP_EDGE_EAST) then
               c = g%ws%src_e(i, j, s)
            else
               c = g%ws%src_n(i, j, s)
            end if
            if (c == 0_int8) cycle
            any_used = .true.
            di = remap_src_di(c)
            dj = remap_src_dj(c)
            if (g%hm(i + di, j + dj) == 0.0_wp) nbad = nbad + 1
         end do
         if (any_used) then
            ic = merge(i + 1, i, iedge == REMAP_EDGE_EAST)
            jc = merge(j, j + 1, iedge == REMAP_EDGE_EAST)
            if (open_mask == 0.0_wp .or. g%hm(i, j) == 0.0_wp .or. g%hm(ic, jc) == 0.0_wp) nbad = nbad + 1
            nopen = nopen + 1
         end if
      end subroutine check_edge
   end subroutine test_coast

   ! ------------------------------------------------------------------
   ! (8) kernels against a host reference
   ! ------------------------------------------------------------------

   subroutine test_vs_host(error)
      !! Library kernels (device on the offload build) against a plain host reference on a random
      !! field over a non-uniform metric: every output of every edge, source codes and case codes
      !! EXACT, areas and coordinates to a measured tolerance, the ring and the unused vertices untouched.
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: efa, i, j, s, k, i2, j2, di, dj
      real(wp) :: dx, dy, uu, vv, ex, ey, sx, sy, mpx, mpy, mx, my, ump, vmp, worst_dp, worst_a, worst_v, w, tol_dp, tol_a
      real(wp) :: asum, dlx, dly, drx, dry, fc, ea
      real(wp) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      real(wp), allocatable :: rdpx(:, :), rdpy(:, :), rea_e(:, :), rea_n(:, :)
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: ic, bits
      integer :: iedge, ii, jj
      character(len=200) :: msg
      real(wp), allocatable :: puc(:, :), pvc(:, :)

      call geo_init(g, 16, 13, .true.)
      g%dt = 1.0_wp
      worst_dp = 0.0_wp
      worst_a = 0.0_wp
      worst_v = 0.0_wp
#ifdef RDB_GPU_OFFLOAD
      tol_dp = 1.0e-10_wp
      tol_a = 1.0e-14_wp
      ! measured on V100: corner stages bitwise, areas 4e-17, vertices 1.1e-16 (FMA contraction on)
#else
      tol_dp = 0.0_wp
      tol_a = 0.0_wp
#endif
      allocate (rdpx(g%nx, g%ny), rdpy(g%nx, g%ny), rea_e(g%nx, g%ny), rea_n(g%nx, g%ny), puc(g%nx, g%ny), pvc(g%nx, g%ny))
      do efa = 0, 1
         call set_random_faces(g, 0.45_wp*g%dxU(1, 1), 4242_int64 + int(efa, int64))
         call geo_sentinel(g)
         call geo_map(g)
         call geo_push(g)
         call remap_geometry(g%ws, g%is, g%ie, g%js, g%je, g%dt, efa, PUNY, g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, &
                             g%narea, g%epm, g%npm, g%uvm, g%ierr)
         call geo_pull(g)
         ! ---- host reference ----
         puc = SENT
         pvc = SENT
         rdpx = SENT
         rdpy = SENT
         rea_e = SENT
         rea_n = SENT
         do j = g%js - 2, g%je + 1
            do i = g%is - 2, g%ie + 1
               ! the dynamics' average: area-weighted, times the corner mask
               puc(i, j) = (g%u(i + 1, j)*g%earea(i, j) + g%u(i + 1, j + 1)*g%earea(i, j + 1)) &
                           /(g%earea(i, j) + g%earea(i, j + 1))*g%uvm(i, j)
               pvc(i, j) = (g%v(i, j + 1)*g%narea(i, j) + g%v(i + 1, j + 1)*g%narea(i + 1, j)) &
                           /(g%narea(i, j) + g%narea(i + 1, j))*g%uvm(i, j)
            end do
         end do
         do j = g%js - 1, g%je
            do i = g%is - 1, g%ie
               uu = puc(i, j)
               vv = pvc(i, j)
               rdpx(i, j) = 0.0_wp
               rdpy(i, j) = 0.0_wp
               if (uu /= 0.0_wp .or. vv /= 0.0_wp) then
                  ex = -g%dt*uu
                  ey = -g%dt*vv
                  sx = ex/g%dxU(i, j)
                  sy = ey/g%dyU(i, j)
                  mpx = 0.5_wp*sx
                  mpy = 0.5_wp*sy
                  if (mpx >= 0.0_wp .and. mpy >= 0.0_wp) then
                     i2 = i + 1
                     j2 = j + 1
                     mx = mpx - 0.5_wp
                     my = mpy - 0.5_wp
                  else if (mpx < 0.0_wp .and. mpy < 0.0_wp) then
                     i2 = i
                     j2 = j
                     mx = mpx + 0.5_wp
                     my = mpy + 0.5_wp
                  else if (mpx >= 0.0_wp) then
                     i2 = i + 1
                     j2 = j
                     mx = mpx - 0.5_wp
                     my = mpy + 0.5_wp
                  else
                     i2 = i
                     j2 = j + 1
                     mx = mpx + 0.5_wp
                     my = mpy - 0.5_wp
                  end if
                  ump = puc(i2 - 1, j2 - 1)*(mx - 0.5_wp)*(my - 0.5_wp) - puc(i2, j2 - 1)*(mx + 0.5_wp)*(my - 0.5_wp) &
                        + puc(i2, j2)*(mx + 0.5_wp)*(my + 0.5_wp) - puc(i2 - 1, j2)*(mx - 0.5_wp)*(my + 0.5_wp)
                  vmp = pvc(i2 - 1, j2 - 1)*(mx - 0.5_wp)*(my - 0.5_wp) - pvc(i2, j2 - 1)*(mx + 0.5_wp)*(my - 0.5_wp) &
                        + pvc(i2, j2)*(mx + 0.5_wp)*(my + 0.5_wp) - pvc(i2 - 1, j2)*(mx - 0.5_wp)*(my + 0.5_wp)
                  rdpx(i, j) = -g%dt*ump
                  rdpy(i, j) = -g%dt*vmp
               end if
            end do
         end do
         do j = g%js, g%je
            do i = g%is - 1, g%ie
               rea_e(i, j) = g%u(i + 1, j)*g%HTE(i, j)*g%dt
            end do
         end do
         do j = g%js - 1, g%je
            do i = g%is, g%ie
               rea_n(i, j) = g%v(i, j + 1)*g%HTN(i, j)*g%dt
            end do
         end do
         ! ---- compare the corner stages (whole arrays: the ring must still hold the sentinel) ----
         do j = 1, g%ny
            do i = 1, g%nx
               worst_dp = max(worst_dp, abs(g%ws%dpx(i, j) - rdpx(i, j)), abs(g%ws%dpy(i, j) - rdpy(i, j)), &
                              abs(g%ws%uc(i, j) - puc(i, j)), abs(g%ws%vc(i, j) - pvc(i, j)))
            end do
         end do
         if (.not. (worst_dp <= tol_dp)) then
            write (msg, "(a,es10.3)") "host reference: corner stages differ by ", worst_dp
            call test_failed(error, trim(msg))
            call geo_destroy(g)
            return
         end if
         ! ---- edges: gather scalars, call the table, compare slot by slot ----
         do iedge = REMAP_EDGE_NORTH, REMAP_EDGE_EAST
            do j = merge(g%js - 1, g%js, iedge == REMAP_EDGE_NORTH), g%je
               do i = merge(g%is, g%is - 1, iedge == REMAP_EDGE_NORTH), g%ie
                  call gather_edge(g, iedge, i, j, rdpx, rdpy, dlx, dly, drx, dry, fc, px, py, ar)
                  ! reuse px.. as scratch: gather returns the scalars only
                  call call_table(g, iedge, i, j, efa, rdpx, rdpy, rea_e, rea_n, px, py, ar, sc, ic, asum, bits, ea)
                  if (bits /= 0_int32) cycle
                  if (iedge == REMAP_EDGE_NORTH) then
                     w = abs(g%ws%edgearea_n(i, j) - merge(rea_n(i, j), asum, efa /= 0))/fc
                  else
                     w = abs(g%ws%edgearea_e(i, j) - merge(rea_e(i, j), asum, efa /= 0))/fc
                  end if
                  worst_a = max(worst_a, w)
                  do s = 1, REMAP_NSLOT
                     if (iedge == REMAP_EDGE_NORTH) then
                        if (g%ws%src_n(i, j, s) /= sc(s)) then
                           call test_failed(error, "host reference: source code differs")
                           call geo_destroy(g)
                           return
                        end if
                        w = abs(g%ws%area_n(i, j, s) - ar(s))/fc
                     else
                        if (g%ws%src_e(i, j, s) /= sc(s)) then
                           call test_failed(error, "host reference: source code differs")
                           call geo_destroy(g)
                           return
                        end if
                        w = abs(g%ws%area_e(i, j, s) - ar(s))/fc
                     end if
                     worst_a = max(worst_a, w)
                     if (sc(s) /= 0_int8) then
                        do k = 1, 3
                           if (iedge == REMAP_EDGE_NORTH) then
                              w = max(abs(g%ws%xp_n(i, j, k, s) - px(k, s)), abs(g%ws%yp_n(i, j, k, s) - py(k, s)))
                           else
                              w = max(abs(g%ws%xp_e(i, j, k, s) - px(k, s)), abs(g%ws%yp_e(i, j, k, s) - py(k, s)))
                           end if
                           worst_v = max(worst_v, w)
                        end do
                     else
                        do k = 1, 3
                           if (iedge == REMAP_EDGE_NORTH) then
                              if (g%ws%xp_n(i, j, k, s) /= SENT .or. g%ws%yp_n(i, j, k, s) /= SENT) then
                                 call test_failed(error, "host reference: vertices of an unused slot were written")
                                 call geo_destroy(g)
                                 return
                              end if
                           else
                              if (g%ws%xp_e(i, j, k, s) /= SENT .or. g%ws%yp_e(i, j, k, s) /= SENT) then
                                 call test_failed(error, "host reference: vertices of an unused slot were written")
                                 call geo_destroy(g)
                                 return
                              end if
                           end if
                        end do
                     end if
                  end do
                  if (iedge == REMAP_EDGE_NORTH) then
                     if (g%ws%icase_n(i, j) /= ic) then
                        write (msg, "(a,2i4,2i8)") "host reference: case code differs at north edge ", i, j, &
                           g%ws%icase_n(i, j), ic
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                  else
                     if (g%ws%icase_e(i, j) /= ic) then
                        write (msg, "(a,2i4,2i8)") "host reference: case code differs at east edge ", i, j, &
                           g%ws%icase_e(i, j), ic
                        call test_failed(error, trim(msg))
                        call geo_destroy(g)
                        return
                     end if
                  end if
               end do
            end do
         end do
         ! ---- ring untouched: outside the edge windows the sentinels survive ----
         do j = 1, g%ny
            do i = 1, g%nx
               ii = i
               jj = j
               if (ii >= g%is - 1 .and. ii <= g%ie .and. jj >= g%js .and. jj <= g%je) cycle
               if (g%ws%icase_e(i, j) /= 7777_int32 .or. any(g%ws%src_e(i, j, :) /= 77_int8) .or. &
                   any(g%ws%area_e(i, j, :) /= SENT) .or. any(g%ws%xp_e(i, j, :, :) /= SENT)) then
                  call test_failed(error, "ring: an east-edge output outside its window was written")
                  call geo_destroy(g)
                  return
               end if
            end do
         end do
         do j = 1, g%ny
            do i = 1, g%nx
               if (i >= g%is .and. i <= g%ie .and. j >= g%js - 1 .and. j <= g%je) cycle
               if (g%ws%icase_n(i, j) /= 7777_int32 .or. any(g%ws%src_n(i, j, :) /= 77_int8) .or. &
                   any(g%ws%area_n(i, j, :) /= SENT) .or. any(g%ws%yp_n(i, j, :, :) /= SENT)) then
                  call test_failed(error, "ring: a north-edge output outside its window was written")
                  call geo_destroy(g)
                  return
               end if
            end do
         end do
         di = 0
         dj = 0
      end do
      write (error_unit, "(a,es9.2,a,es9.2,a,es9.2)") "  [remap_geom] vs host reference: worst corner-stage deviation ", worst_dp, &
         ", area/areafac ", worst_a, ", vertex ", worst_v
      if (.not. (worst_a <= tol_a .and. worst_v <= tol_a)) then
         call test_failed(error, "host reference: areas or vertices differ beyond the measured tolerance")
      end if
      call geo_destroy(g)
      deallocate (rdpx, rdpy, rea_e, rea_n, puc, pvc)
   end subroutine test_vs_host

   subroutine gather_edge(g, iedge, i, j, dpx, dpy, dlx, dly, drx, dry, fc, px, py, ar)
      !! The scalars of one edge from the host departure points (own indexing, written from the
      !! documented conventions).
      type(geo_t), intent(in) :: g
      integer, intent(in) :: iedge, i, j
      real(wp), intent(in) :: dpx(:, :), dpy(:, :)
      real(wp), intent(out) :: dlx, dly, drx, dry, fc
      real(wp), intent(out) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      px = 0.0_wp
      py = 0.0_wp
      ar = 0.0_wp
      if (iedge == REMAP_EDGE_NORTH) then
         dlx = dpx(i - 1, j)/g%dxU(i - 1, j)
         dly = dpy(i - 1, j)/g%dyU(i - 1, j)
         drx = dpx(i, j)/g%dxU(i, j)
         dry = dpy(i, j)/g%dyU(i, j)
         fc = g%narea(i, j)
      else
         dlx = -dpy(i, j)/g%dyU(i, j)
         dly = dpx(i, j)/g%dxU(i, j)
         drx = -dpy(i, j - 1)/g%dyU(i, j - 1)
         dry = dpx(i, j - 1)/g%dxU(i, j - 1)
         fc = g%earea(i, j)
      end if
   end subroutine gather_edge

   subroutine call_table(g, iedge, i, j, efa, dpx, dpy, rea_e, rea_n, px, py, ar, sc, ic, asum, bits, ea)
      !! The table on the host reference scalars (the neighbours' area factors from the documented positions).
      type(geo_t), intent(in) :: g
      integer, intent(in) :: iedge, i, j, efa
      real(wp), intent(in) :: dpx(:, :), dpy(:, :), rea_e(:, :), rea_n(:, :)
      real(wp), intent(out) :: asum, ea
      real(wp), intent(out) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      integer(int8), intent(out) :: sc(REMAP_NSLOT)
      integer(int32), intent(out) :: ic, bits
      real(wp) :: dlx, dly, drx, dry, fc
      real(wp), parameter :: ZERO = 0.0_wp
      logical :: active

      call gather_edge(g, iedge, i, j, dpx, dpy, dlx, dly, drx, dry, fc, px, py, ar)
      if (iedge == REMAP_EDGE_NORTH) then
         ea = merge(rea_n(i, j), ZERO, efa /= 0)
         active = dpx(i - 1, j) /= ZERO .or. dpy(i - 1, j) /= ZERO .or. dpx(i, j) /= ZERO .or. dpy(i, j) /= ZERO
         if (efa /= 0) active = active .or. rea_n(i, j) /= ZERO
         sc = 0_int8
         ic = 0_int32
         asum = ZERO
         bits = 0_int32
         if (.not. active) return
         call remap_locate_edge(iedge, efa, PUNY, dlx, dly, drx, dry, fc, g%earea(i - 1, j + 1), g%earea(i - 1, j), &
                                g%earea(i, j + 1), g%earea(i, j), g%narea(i - 1, j), g%narea(i + 1, j), ea, px, py, ar, sc, &
                                ic, asum, bits)
      else
         ea = merge(rea_e(i, j), ZERO, efa /= 0)
         active = dpx(i, j - 1) /= ZERO .or. dpy(i, j - 1) /= ZERO .or. dpx(i, j) /= ZERO .or. dpy(i, j) /= ZERO
         if (efa /= 0) active = active .or. rea_e(i, j) /= ZERO
         sc = 0_int8
         ic = 0_int32
         asum = ZERO
         bits = 0_int32
         if (.not. active) return
         call remap_locate_edge(iedge, efa, PUNY, dlx, dly, drx, dry, fc, g%narea(i + 1, j), g%narea(i, j), &
                                g%narea(i + 1, j - 1), g%narea(i, j - 1), g%earea(i, j + 1), g%earea(i, j - 1), ea, px, py, &
                                ar, sc, ic, asum, bits)
      end if
   end subroutine call_table

   ! ------------------------------------------------------------------
   ! (9) throughput
   ! ------------------------------------------------------------------

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      type(geo_t) :: g
      integer :: nxp, nyp, rep
      integer, parameter :: NREP = 5
      integer(int64) :: t0, t1, rate
      real(wp) :: tg, edges

#ifdef RDB_GPU_OFFLOAD
      nxp = 1440
      nyp = 1152
#else
      nxp = 512
      nyp = 512
#endif
      call geo_init(g, nxp, nyp, .true.)
      g%dt = 1.0_wp
      call set_faces_from_fn(g, 3, 0.2_wp*g%dxU(1, 1), 0.05_wp*g%dxU(1, 1))
      call geo_map(g)
      call geo_push(g)
      call remap_geometry(g%ws, g%is, g%ie, g%js, g%je, g%dt, 1, PUNY, g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, &
                          g%narea, g%epm, g%npm, g%uvm, g%ierr)
      call system_clock(t0, rate)
      do rep = 1, NREP
         call remap_geometry(g%ws, g%is, g%ie, g%js, g%je, g%dt, 1, PUNY, g%u, g%v, g%dxU, g%dyU, g%HTE, g%HTN, g%earea, &
                             g%narea, g%epm, g%npm, g%uvm, g%ierr)
      end do
      call system_clock(t1)
      call geo_pull(g)
      tg = real(t1 - t0, wp)/real(rate, wp)/real(NREP, wp)
      edges = 2.0_wp*real(nxp, wp)*real(nyp, wp)
      write (error_unit, "(a,i0,a,i0,a,es10.3,a,f8.2,a,es10.3,a,f8.1,a)") "  [remap_geom] ", nxp, "x", nyp, ": geometry (5 launches, efa) ", &
         edges/tg, " edges/s (", tg*1000.0_wp, " ms), ", real(nxp, wp)*real(nyp, wp)/tg, " cells/s; workspace ", &
         real(g%ws%bytes(), wp)/1.0e6_wp, " MB"
      call check(error, ieee_is_finite(edges/tg) .and. all(g%ierr(g%is:g%ie, g%js:g%je) == 0_int32), "clean and finite")
      call geo_destroy(g)
   end subroutine test_throughput

end module test_cice_remap_geom
