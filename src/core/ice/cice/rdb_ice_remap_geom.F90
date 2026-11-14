!! Incremental-remapping geometry of the sea-ice transport on the C grid: the
!! departure points, the edge areas and the departure-region triangles
!! (Dukowicz and Baumgardner 2000; Lipscomb and Hunke 2004).
module rdb_ice_remap_geom
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_remap (departure points, the departure-region triangle table, the edge flux adjustment).
   ! Translated to roundabout's do-concurrent/GPU form and modified by the
   ! roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   ! Original authors (from the upstream header): William H. Lipscomb (LANL),
   ! John R. Baumgardner (LANL), Elizabeth C. Hunke (LANL).
   ! Physics: the incremental-remapping transport of Dukowicz and Baumgardner
   ! (2000, J. Comput. Phys. 160 318) and Lipscomb and Hunke (2004, Mon. Wea.
   ! Rev. 132 1341); the edge flux adjustment is the idea of M. Bentsen
   ! discussed there.
   !
   !! **What this module is.**  The geometric half of the remapping: for every
   !! cell edge, the region the fluid crossing the edge in one step came from
   !! (the DEPARTURE REGION) cut into at most six triangles, each lying in one
   !! source cell, with its vertices in that cell's own frame and its signed
   !! area.  The flux integrals (a later PR) integrate the reconstructed
   !! fields over those triangles.  Nothing here depends on a category or a
   !! tracer, so the geometry is built ONCE per step and reused by every
   !! `ncat x ntrace` integral.
   !!
   !! **Stages (each ONE `do concurrent` launch over its range).**
   !!  1. `remap_corner_velocity`: the face velocities to the NE corners with
   !!     the dynamics' OWN averaging (`evp_average_impl`, the same expression):
   !!     `uc(i,j)` the `earea`-weighted mean of the two east faces `E(i,j)`,
   !!     `E(i,j+1)`, `vc(i,j)` the `narea`-weighted mean of the two north faces,
   !!     both MULTIPLIED BY THE CORNER MASK `uvm` (zero at any corner touching
   !!     land).  A caller that already holds the dynamics' `uvelU/vvelU` passes
   !!     them to `remap_geometry` and this stage is skipped.  The mask is what
   !!     keeps a closed face (zero face velocity, land across it) from getting
   !!     triangles: both its corners touch land, so both displacements are zero,
   !!     and its prescribed area is zero, so the edge is inactive.
   !!  2. `remap_departure_points`: at each corner the back trajectory
   !!     `-dt (u, v)` improved with the velocity at the trajectory MIDPOINT,
   !!     interpolated bilinearly from the four corner velocities around it
   !!     (Dukowicz and Baumgardner 2000), giving `dpx, dpy` (m).
   !!  3. `remap_edge_areas`: the area prescribed for each edge,
   !!     `edgearea_e(i,j) = u(i+1,j) HTE(i,j) dt`, `edgearea_n(i,j) =
   !!     v(i,j+1) HTN(i,j) dt` (m^2), POSITIVE for eastward / northward
   !!     flow, i.e. positive when fluid leaves cell `(i,j)` across the edge.
   !!     It is exactly the continuity flux of the same face velocity.
   !!  4. `remap_locate_east` / `remap_locate_north`: the triangles.
   !! `remap_geometry` runs all of them (six launches: the corner average
   !! computes four arrays in one launch, then departure points, the two edge
   !! areas, east and north triangles).
   !!
   !! **Index and frame conventions** (those of `rdb_ice_cgrid_metrics`).
   !! Arrays are `(nx, ny)` with ghost rings, addressed by the cell they
   !! belong to: the EAST edge `E(i,j)` is the east face of cell `(i,j)`, the
   !! NORTH edge `N(i,j)` its north face, the corner `U(i,j)` its NE corner;
   !! `u(nx+1, ny)` has the east face of `(i,j)` at `u(i+1,j)`, `v(nx, ny+1)`
   !! the north face at `v(i,j+1)`.  `dxU, dyU` are the corner lengths,
   !! `HTE, HTN` the east / north face lengths, `earea, narea` the areas of
   !! the cells centred on the east / north faces (the area factors).
   !!
   !! **Departure region.**  In a frame with its origin at the edge midpoint,
   !! the edge runs from `CL = (-1/2, 0)` to `CR = (1/2, 0)` in units of the
   !! corner lengths; for a NORTH edge `x` is east, `y` north (`CL` is the NW
   !! corner `U(i-1,j)`, `CR` the NE corner `U(i,j)`); for an EAST edge the
   !! frame is turned by 90 degrees (`x` points south, `y` east; `CL` is the
   !! NE corner `U(i,j)`, `CR` the SE corner `U(i,j-1)`), so one routine
   !! serves both.  `DL`, `DR` are the departure points of `CL`, `CR`
   !! (`CL + dp/dxU`, scaled by the corner lengths), `DM = (DL + DR)/2`.  The
   !! region is `CL-CR-DR-DL`; it overlaps up to six cells: `TL TC TR` on the
   !! far side of the edge, `BL BC BR` on the near side (`BC` is the cell
   !! `(i,j)` itself, `TC` the one the edge leads into).  A back trajectory
   !! must stay within one cell: a corner whose scaled displacement exceeds 1
   !! in either component, or departure points that cross (`DR` not east of
   !! `DL`), raise `ICE_ERR_TRANSPORT_CFL` (the reference checks only one
   !! sign, against the wrong length); outputs are then unspecified.
   !!
   !! **The case table.**  The side cells are decided by where the lines
   !! `DL-DM`, `DM-DR` cut the sides `x = -+1/2` and the line `DL-DR` cuts the
   !! axis: left `{TL, BL, TL1+BL1, TL2+BL2}`, right `{TR, BR, TR1+BR1,
   !! TR2+BR2}`, and the centre column (two departure points inside it) has
   !! twelve branches on the signs of `ydl, ydr, ydm` and `xic`.  They are the
   !! reference's, one for one; `icase` returns the branch taken per edge as
   !! `left + 5 right + 25 centre + 325 efa` (left, right in 0..4, centre in
   !! 0..12, efa in 0..4; zero = branch not taken), which is how the tests
   !! prove every branch was run.  Twelve centre branches: 1 `TC1a/2a/3a`,
   !! 2 the rare `TC1b/2b + BC3b` (both points above, `DM` below), 3
   !! `BC1a/2a/3a`, 4 the rare mirror, 5-8 `DL` above and `DR` below with
   !! `xic >=< 0` and `ydm >=< 0`, 9-12 the mirror (`DL` below, `DR` above).
   !!
   !! **Edge flux adjustment (EFA, `efa = 1`: the C-grid mode).**  The
   !! departure region is tweaked, following Bentsen, so that its area is the
   !! prescribed `edgearea`: the side triangles are kept, and the midpoint
   !! `DM` of the central part is shifted along the normal of `DL-DR` until
   !! the centre polygon has the area left over (`edgearea` minus the side
   !! areas); when the two departure points lie in different cells the point
   !! where `DL-DR` meets the axis is pinned and the lone triangle's area is
   !! taken off first.  This is what makes the area transported across an
   !! edge identical to the continuity flux (`sum(triarea) = edgearea`), and
   !! it is why the full table is live.  The sum equals `edgearea` TO ROUND-OFF
   !! (a few ulp of the cell area, not bitwise: the reference computes the
   !! side areas, the shifted centre and the final triangles by different
   !! expressions, and drops a triangle smaller than `1e-16` cell areas); the
   !! check `|sum - edgearea| <= 1e-13 areafac` raises
   !! `ICE_ERR_TRANSPORT_CONSERV`.  A one-cell-wide channel (both
   !! departure points at the corners, `edgearea /= 0`) is an active edge:
   !! the shifted midpoint alone carries the area.  `efa = 0` is the
   !! unadjusted (B-grid) region: `edgearea` is then an OUTPUT, the sum of the
   !! triangle areas.
   !!
   !! **Output layout (for the flux integration).**  Per edge direction
   !! (suffix `_e` east, `_n` north), `s` = triangle slot 1..6 (the
   !! reference's groups; at most five are non-empty at once):
   !!  - `src(i,j,s)` `integer(int8)`: 0 = slot UNUSED (also when the edge is
   !!    inactive or the triangle was dropped as smaller than `1e-16` cell
   !!    areas); else `1 + (di+1) + 3 (dj+1)` with `(di, dj)` in `-1:1` the
   !!    offset of the SOURCE cell from the edge's own cell `(i,j)`
   !!    (decode with `remap_src_di`, `remap_src_dj`).
   !!  - `area(i,j,s)`: the signed physical area (m^2) = (positive area factor)
   !!    x (area of the triangle in cell units, taken with the sign of the
   !!    departure region's orientation).  POSITIVE when the fluid it holds
   !!    moves from `BC` (the cell `(i,j)`) toward `TC`: a triangle in `TC`
   !!    carries a NEGATIVE area (inflow).  Their sum over the slots is the
   !!    net area leaving cell `(i,j)` through the edge.  Zero for an unused
   !!    slot.
   !!  - `xp(i,j,v,s)`, `yp(i,j,v,s)`, `v = 1..3`: the vertices in the SOURCE
   !!    cell's own frame (origin at its centre, `x` east, `y` north, cell
   !!    sides at `+-1/2`), after the vertex translation (and, for an
   !!    east edge, the quarter turn back).  Defined only where `src /= 0`.
   !!    Most triangles lie inside their cell, but a pair of triangles of
   !!    OPPOSITE sign in one source cell (the "rare" centre branches, where the
   !!    line `DL-DR` meets the axis outside the edge) may reach beyond it, by
   !!    at most one cell: the signed areas cancel exactly, the part outside
   !!    the cell being counted once positive and once negative, so the flux
   !!    integral must evaluate the (linear) reconstruction there by
   !!    extrapolation, not clip it.  A vertex beyond `|x|,|y| > 3/2`
   !!    raises `ICE_ERR_TRANSPORT_CFL` (a shifted midpoint that left the two
   !!    cells of the edge).
   !!  - `icase(i,j)` `integer(int32)`: the branch code above; 0 for an
   !!    inactive edge.
   !!  - `edgearea_e/_n(i,j)`: input (EFA) or output (`efa = 0`).
   !! An edge is ACTIVE when either of its two corners has a non-zero
   !! departure displacement or (with `efa`) `edgearea /= 0`; an inactive
   !! edge gets `src = 0`, `area = 0`, `icase = 0`.
   !!
   !! **Rule for the flux integration (matches the reference).**  With
   !! `<R_c, T>` the integral of the reconstruction of cell `c` over triangle
   !! `T` evaluated at the STORED vertices in THAT cell's frame, UNCLIPPED:
   !! `F_e(i,j) = sum_s area_e(i,j,s) * <R_{(i+di, j+dj)}, triangle s>` (same for
   !! `F_n`, with the source offset of its own slot), and
   !! `delta(mass*area)(i,j) = -(F_e(i,j) - F_e(i-1,j) + F_n(i,j) - F_n(i,j-1))`.
   !! Never take `|area|`: the sign carries the inflow and the cancelling pairs.
   !! Vertices may lie up to one cell outside their source cell (see below), so
   !! each triangle integral must be exact for the (linear) reconstruction
   !! there, and the monotonicity tests of the flux PR must hit those branches
   !! (centre branches 2 and 4 and the pinned-axis ones).
   !!
   !! **Ranges and halos.**  The caller gives the owned window `is:ie, js:je`
   !! and `nghost >= 3` of valid input around it.  East edges `E(is-1:ie,
   !! js:je)` and north edges `N(is:ie, js-1:je)` are produced (the cell
   !! `(i,j)` gathers `E(i-1,j), E(i,j), N(i,j-1), N(i,j)`); the departure
   !! points on `U(is-1:ie, js-1:je)`; the corner velocities on
   !! `U(is-2:ie+1, js-2:je+1)`.  Inputs, when the corner velocities are
   !! computed here: faces `u, v` on `(is-2:ie+2, js-2:je+2)`; masks `epm,
   !! npm, uvm` on `(is-2:ie+1, js-2:je+1)`; `earea` on `(is-3:ie+1,
   !! js-2:je+2)`, `narea` on `(is-2:ie+2, js-3:je+1)`; (with caller-supplied
   !! corner velocities the faces are needed only for the edge areas, `u(is:ie+1,
   !! js:je)`, `v(is-1:ie, js:je+1)` of the face index, and the corner arrays on
   !! the window above.)  Always: `dxU, dyU` on the departure window, `HTE` on
   !! `E(is-1:ie, js:je)`, `HTN` on `N(is:ie, js-1:je)`, `earea` on
   !! `(is-1:ie, js-1:je+1)` and `narea` on `(is-1:ie+1, js-1:je)` (a corner
   !! triangle takes the area factor of the neighbouring edge's cell, so the
   !! same corner region seen from its two edges has the same physical area).
   !! Nothing outside the output windows is written; the input windows are
   !! pinned by `test_cice_remap_geom` (NaN outside them changes nothing).
   !!
   !! **GPU shape.**  Every stage is one `do concurrent` over its edges or
   !! corners (loop ranges, no edge guards), kernels in the outer-shim +
   !! flat-impl form with explicit-shape dummies.  The case table is branch
   !! code in ONE `pure` `!$acc routine seq` per-edge routine
   !! (`remap_edge_one`, indices `integer, value`) over fixed-size locals (six
   !! slots of 3 x 2 coordinates, area and source), writing the global arrays
   !! at `(i,j)` only: no scatter, no atomics, no allocation, no per-call
   !! bundle.  `remap_locate_edge` is the same table on scalars for
   !! host-side use.  The persistent workspace maps component-wise.
   !!
   !! **Errors.**  Bits OR-ed into the cell `ierr(i,j)` of the corner or edge:
   !! `ICE_ERR_TRANSPORT_CFL` (displacement > 1 cell, crossing departure
   !! points, a vertex more than one cell outside its source cell),
   !! `ICE_ERR_NONFINITE` (non-finite velocity, displacement, area or metric where one is
   !! needed, and unusable geometry input: a non-positive corner length under a non-zero
   !! velocity, a non-positive area factor under the adjustment; the registry has no separate
   !! "bad geometry" bit),
   !! `ICE_ERR_TRANSPORT_CONSERV` (the EFA sum off `edgearea`).  Nothing is
   !! clamped; the driver does not commit a step with `ierr /= 0`.
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_evp_c, only: evp_average_impl
   implicit none
   private

   public :: remap_geom_workspace_t
   public :: remap_corner_velocity, remap_departure_points, remap_edge_areas
   public :: remap_locate_east, remap_locate_north, remap_geometry
   public :: remap_locate_edge, remap_src_di, remap_src_dj
   public :: REMAP_NSLOT, REMAP_NCASE_LEFT, REMAP_NCASE_RIGHT, REMAP_NCASE_CENTRE, REMAP_NCASE_EFA
   public :: REMAP_EDGE_NORTH, REMAP_EDGE_EAST, REMAP_VERTEX_MAX

   integer, parameter :: REMAP_NSLOT = 6
      !! Triangle slots per edge (the reference's groups).
   integer, parameter :: REMAP_NCASE_LEFT = 4
   integer, parameter :: REMAP_NCASE_RIGHT = 4
   integer, parameter :: REMAP_NCASE_CENTRE = 12
   integer, parameter :: REMAP_NCASE_EFA = 4
      !! Number of branches per table (1..N; 0 = not taken).
   integer, parameter :: REMAP_EDGE_NORTH = 1
   integer, parameter :: REMAP_EDGE_EAST = 2
   real(wp), parameter :: REMAP_VERTEX_MAX = 1.5_wp
      !! Largest `|x|, |y|` a triangle vertex may have in its source cell's frame (cell units).

   real(wp), parameter :: EPS16 = 1.0e-16_wp
   real(wp), parameter :: EPS13 = 1.0e-13_wp

   ! Cell labels: 1 TL, 2 TR, 3 BL, 4 BR, 5 TC, 6 BC; offsets of the cell from the edge's own cell.
   integer, parameter :: LTL = 1, LTR = 2, LBL = 3, LBR = 4, LTC = 5, LBC = 6
   integer(int8), parameter :: SRC_CODE(-1:1, -1:1) = reshape([1_int8, 2_int8, 3_int8, 4_int8, 5_int8, 6_int8, 7_int8, 8_int8, 9_int8], [3, 3])
      !! `1 + (di+1) + 3 (dj+1)` of the source offset (a table: no run-time narrowing in device code).
   integer, parameter :: NLABEL = 6, NEDGE_KIND = 2
   integer, parameter :: SH_I(NLABEL, NEDGE_KIND) = reshape([-1, 1, -1, 1, 0, 0, 1, 1, 0, 0, 1, 0], [NLABEL, NEDGE_KIND])
   integer, parameter :: SH_J(NLABEL, NEDGE_KIND) = reshape([1, 1, 0, 0, 1, 0, 1, -1, 1, -1, 0, 0], [NLABEL, NEDGE_KIND])

   type :: remap_geom_workspace_t
      !! Persistent geometry buffers (see the module docstring for the
      !! layout).  `init` allocates, `enter_data` maps the component arrays
      !! ONCE (never the aggregate), `exit_data` unmaps, `destroy` frees.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      real(wp), allocatable :: uc(:, :), vc(:, :)
         !! Corner velocities (m/s), `U(i,j)`.
      real(wp), allocatable :: uN(:, :), vE(:, :)
         !! Scratch of the dynamics' averaging (tangential components at N and E; not an output).
      real(wp), allocatable :: dpx(:, :), dpy(:, :)
         !! Departure displacements (m) at the corners.
      real(wp), allocatable :: edgearea_e(:, :), edgearea_n(:, :)
         !! Edge areas (m^2), `> 0` for eastward / northward flow.
      real(wp), allocatable :: xp_e(:, :, :, :), yp_e(:, :, :, :), area_e(:, :, :)
      real(wp), allocatable :: xp_n(:, :, :, :), yp_n(:, :, :, :), area_n(:, :, :)
      integer(int8), allocatable :: src_e(:, :, :), src_n(:, :, :)
      integer(int32), allocatable :: icase_e(:, :), icase_n(:, :)
   contains
      procedure :: init => workspace_init
      procedure :: enter_data => workspace_enter_data
      procedure :: exit_data => workspace_exit_data
      procedure :: destroy => workspace_destroy
      procedure :: bytes => remap_ws_bytes
   end type remap_geom_workspace_t

contains

   ! ------------------------------------------------------------------
   ! workspace
   ! ------------------------------------------------------------------

   subroutine workspace_init(this, nx, ny)
      class(remap_geom_workspace_t), intent(inout) :: this
      integer, intent(in) :: nx, ny

      call this%destroy()
      this%nx = nx
      this%ny = ny
      allocate (this%uc(nx, ny), this%vc(nx, ny), this%uN(nx, ny), this%vE(nx, ny), this%dpx(nx, ny), this%dpy(nx, ny))
      allocate (this%edgearea_e(nx, ny), this%edgearea_n(nx, ny))
      allocate (this%xp_e(nx, ny, 3, REMAP_NSLOT), this%yp_e(nx, ny, 3, REMAP_NSLOT), this%area_e(nx, ny, REMAP_NSLOT))
      allocate (this%xp_n(nx, ny, 3, REMAP_NSLOT), this%yp_n(nx, ny, 3, REMAP_NSLOT), this%area_n(nx, ny, REMAP_NSLOT))
      allocate (this%src_e(nx, ny, REMAP_NSLOT), this%src_n(nx, ny, REMAP_NSLOT))
      allocate (this%icase_e(nx, ny), this%icase_n(nx, ny))
      this%uc = 0.0_wp
      this%vc = 0.0_wp
      this%uN = 0.0_wp
      this%vE = 0.0_wp
      this%dpx = 0.0_wp
      this%dpy = 0.0_wp
      this%edgearea_e = 0.0_wp
      this%edgearea_n = 0.0_wp
      this%xp_e = 0.0_wp
      this%yp_e = 0.0_wp
      this%area_e = 0.0_wp
      this%xp_n = 0.0_wp
      this%yp_n = 0.0_wp
      this%area_n = 0.0_wp
      this%src_e = 0_int8
      this%src_n = 0_int8
      this%icase_e = 0_int32
      this%icase_n = 0_int32
      this%is_init = .true.
   end subroutine workspace_init

   subroutine workspace_enter_data(this)
      !! Map the buffers once (`copyin`: the host values).
      class(remap_geom_workspace_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%uN, this%vE, this%uc, this%vc, this%dpx, this%dpy, this%edgearea_e, this%edgearea_n)
      !$acc enter data copyin(this%xp_e, this%yp_e, this%area_e, this%src_e, this%icase_e)
      !$acc enter data copyin(this%xp_n, this%yp_n, this%area_n, this%src_n, this%icase_n)
      this%mapped = .true.
   end subroutine workspace_enter_data

   subroutine workspace_exit_data(this)
      class(remap_geom_workspace_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%uN, this%vE, this%uc, this%vc, this%dpx, this%dpy, this%edgearea_e, this%edgearea_n)
      !$acc exit data delete(this%xp_e, this%yp_e, this%area_e, this%src_e, this%icase_e)
      !$acc exit data delete(this%xp_n, this%yp_n, this%area_n, this%src_n, this%icase_n)
      this%mapped = .false.
   end subroutine workspace_exit_data

   pure function arr_bytes_r2(a) result(nb)
      real(wp), allocatable, intent(in) :: a(:, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes_r2

   pure function arr_bytes_r3(a) result(nb)
      real(wp), allocatable, intent(in) :: a(:, :, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes_r3

   pure function arr_bytes_r4(a) result(nb)
      real(wp), allocatable, intent(in) :: a(:, :, :, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes_r4

   pure function arr_bytes_i8(a) result(nb)
      integer(int8), allocatable, intent(in) :: a(:, :, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes_i8

   pure function arr_bytes_i32(a) result(nb)
      integer(int32), allocatable, intent(in) :: a(:, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes_i32

   pure function remap_ws_bytes(this) result(b)
      !! Bytes of the buffers (host copy; the device mirror is the same size).
      class(remap_geom_workspace_t), intent(in) :: this
      integer(int64) :: b
      b = arr_bytes_r2(this%uc) + arr_bytes_r2(this%vc) + arr_bytes_r2(this%uN) + arr_bytes_r2(this%vE) + arr_bytes_r2(this%dpx) + arr_bytes_r2(this%dpy) &
          + arr_bytes_r2(this%edgearea_e) + arr_bytes_r2(this%edgearea_n) &
          + arr_bytes_r4(this%xp_e) + arr_bytes_r4(this%yp_e) + arr_bytes_r3(this%area_e) &
          + arr_bytes_r4(this%xp_n) + arr_bytes_r4(this%yp_n) + arr_bytes_r3(this%area_n) &
          + arr_bytes_i8(this%src_e) + arr_bytes_i8(this%src_n) &
          + arr_bytes_i32(this%icase_e) + arr_bytes_i32(this%icase_n)
   end function remap_ws_bytes

   subroutine workspace_destroy(this)
      class(remap_geom_workspace_t), intent(inout) :: this
      call this%exit_data()
      if (allocated(this%uc)) deallocate (this%uc, this%vc, this%uN, this%vE, this%dpx, this%dpy, this%edgearea_e, this%edgearea_n)
      if (allocated(this%xp_e)) deallocate (this%xp_e, this%yp_e, this%area_e, this%xp_n, this%yp_n, this%area_n)
      if (allocated(this%src_e)) deallocate (this%src_e, this%src_n, this%icase_e, this%icase_n)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
   end subroutine workspace_destroy

   ! ------------------------------------------------------------------
   ! source-cell code
   ! ------------------------------------------------------------------

   pure elemental function remap_src_di(code) result(di)
      !! `i` offset (-1:1) of the source cell encoded in `src` (code /= 0).
      !$acc routine seq
      integer(int8), intent(in) :: code
      integer :: di
      di = mod(int(code) - 1, 3) - 1
   end function remap_src_di

   pure elemental function remap_src_dj(code) result(dj)
      !! `j` offset (-1:1) of the source cell encoded in `src` (code /= 0).
      !$acc routine seq
      integer(int8), intent(in) :: code
      integer :: dj
      dj = (int(code) - 1)/3 - 1
   end function remap_src_dj

   ! ------------------------------------------------------------------
   ! stage 1: corner velocities
   ! ------------------------------------------------------------------

   subroutine remap_corner_velocity(nx, ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm, uN, vE, uc, vc)
      !! Face velocities to the NE corners over `U(is-2:ie+1, js-2:je+1)` with the
      !! dynamics' own averaging (`evp_average_impl`, called on the window
      !! `is:ie-1, js:je-1` whose extension by two is exactly this one): area-weighted,
      !! times the corner mask `uvm`.  `uN, vE` are scratch (the tangential
      !! components the same kernel forms at the N and E points).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), epm(nx, ny), npm(nx, ny), uvm(nx, ny)
      real(wp), intent(inout) :: uN(nx, ny), vE(nx, ny), uc(nx, ny), vc(nx, ny)

      call evp_average_impl(nx, ny, is, ie - 1, js, je - 1, u, v, earea, narea, epm, npm, uvm, uN, vE, uc, vc)
   end subroutine remap_corner_velocity

   ! ------------------------------------------------------------------
   ! stage 2: departure points
   ! ------------------------------------------------------------------

   subroutine remap_departure_points(nx, ny, is, ie, js, je, dt, uc, vc, dxU, dyU, dpx, dpy, ierr)
      !! Back trajectories from the corners `U(is-1:ie, js-1:je)` with the
      !! midpoint velocity correction.  The first estimate is `-dt (u, v)`; its
      !! midpoint lies in the quadrant cell the displacement points into, the
      !! velocity there is the bilinear interpolation of the four corner
      !! velocities around it (in the frame stretched to unit cells), and the
      !! displacement is recomputed from it.  A zero corner velocity gives an
      !! exactly zero displacement.  CFL: a first-estimate or final displacement
      !! longer than the corner length (`dxU` in x, `dyU` in y) raises
      !! `ICE_ERR_TRANSPORT_CFL` (symmetric in sign); a non-finite velocity or
      !! result raises `ICE_ERR_NONFINITE` and the displacement is 0.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: uc(nx, ny), vc(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(inout) :: dpx(nx, ny), dpy(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j, i2, j2
      real(wp) :: uu, vv, ex, ey, sx, sy, mpx, mpy, mpxt, mpyt, ump, vmp, rx, ry
      integer(int32) :: bits

      !$acc data present(uc, vc, dxU, dyU, dpx, dpy, ierr)
      do concurrent(j=js - 1:je, i=is - 1:ie) local(uu, vv, ex, ey, sx, sy, mpx, mpy, mpxt, mpyt, ump, vmp, rx, ry, &
                                                    i2, j2, bits)
         uu = uc(i, j)
         vv = vc(i, j)
         bits = 0_int32
         dpx(i, j) = 0.0_wp
         dpy(i, j) = 0.0_wp
         if (.not. (ieee_is_finite(uu) .and. ieee_is_finite(vv))) then
            bits = ICE_ERR_NONFINITE
         else if (uu /= 0.0_wp .or. vv /= 0.0_wp) then
            ex = -dt*uu
            ey = -dt*vv
            if (.not. (dxU(i, j) > 0.0_wp .and. dyU(i, j) > 0.0_wp)) then
               ! nonzero velocity on a corner without lengths: unusable geometry input
               bits = ICE_ERR_NONFINITE
            else
               sx = ex/dxU(i, j)
               sy = ey/dyU(i, j)
               if (abs(sx) > 1.0_wp .or. abs(sy) > 1.0_wp) bits = ICE_ERR_TRANSPORT_CFL
               mpx = 0.5_wp*sx
               mpy = 0.5_wp*sy
               if (mpx >= 0.0_wp .and. mpy >= 0.0_wp) then
                  i2 = i + 1
                  j2 = j + 1
                  mpxt = mpx - 0.5_wp
                  mpyt = mpy - 0.5_wp
               else if (mpx < 0.0_wp .and. mpy < 0.0_wp) then
                  i2 = i
                  j2 = j
                  mpxt = mpx + 0.5_wp
                  mpyt = mpy + 0.5_wp
               else if (mpx >= 0.0_wp .and. mpy < 0.0_wp) then
                  i2 = i + 1
                  j2 = j
                  mpxt = mpx - 0.5_wp
                  mpyt = mpy + 0.5_wp
               else
                  i2 = i
                  j2 = j + 1
                  mpxt = mpx + 0.5_wp
                  mpyt = mpy - 0.5_wp
               end if
               ump = uc(i2 - 1, j2 - 1)*(mpxt - 0.5_wp)*(mpyt - 0.5_wp) &
                     - uc(i2, j2 - 1)*(mpxt + 0.5_wp)*(mpyt - 0.5_wp) &
                     + uc(i2, j2)*(mpxt + 0.5_wp)*(mpyt + 0.5_wp) &
                     - uc(i2 - 1, j2)*(mpxt - 0.5_wp)*(mpyt + 0.5_wp)
               vmp = vc(i2 - 1, j2 - 1)*(mpxt - 0.5_wp)*(mpyt - 0.5_wp) &
                     - vc(i2, j2 - 1)*(mpxt + 0.5_wp)*(mpyt - 0.5_wp) &
                     + vc(i2, j2)*(mpxt + 0.5_wp)*(mpyt + 0.5_wp) &
                     - vc(i2 - 1, j2)*(mpxt - 0.5_wp)*(mpyt + 0.5_wp)
               rx = -dt*ump
               ry = -dt*vmp
               if (.not. (ieee_is_finite(rx) .and. ieee_is_finite(ry))) then
                  bits = ior(bits, ICE_ERR_NONFINITE)
               else
                  dpx(i, j) = rx
                  dpy(i, j) = ry
                  if (abs(rx) > dxU(i, j) .or. abs(ry) > dyU(i, j)) bits = ior(bits, ICE_ERR_TRANSPORT_CFL)
               end if
            end if
         end if
         if (bits /= 0_int32) ierr(i, j) = ior(ierr(i, j), bits)
      end do
      !$acc end data
   end subroutine remap_departure_points

   ! ------------------------------------------------------------------
   ! stage 3: prescribed edge areas
   ! ------------------------------------------------------------------

   subroutine remap_edge_areas(nx, ny, is, ie, js, je, dt, u, v, HTE, HTN, edgearea_e, edgearea_n)
      !! `edgearea_e(i,j) = u(i+1,j) HTE(i,j) dt` over `E(is-1:ie, js:je)` and
      !! `edgearea_n(i,j) = v(i,j+1) HTN(i,j) dt` over `N(is:ie, js-1:je)`: the
      !! area crossing the face in one step, positive eastward / northward.  The
      !! velocities are used as given (a closed face carries zero).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1), HTE(nx, ny), HTN(nx, ny)
      real(wp), intent(inout) :: edgearea_e(nx, ny), edgearea_n(nx, ny)
      integer :: i, j

      !$acc data present(u, v, HTE, HTN, edgearea_e, edgearea_n)
      do concurrent(j=js:je, i=is - 1:ie)
         edgearea_e(i, j) = u(i + 1, j)*HTE(i, j)*dt
      end do
      do concurrent(j=js - 1:je, i=is:ie)
         edgearea_n(i, j) = v(i, j + 1)*HTN(i, j)*dt
      end do
      !$acc end data
   end subroutine remap_edge_areas

   ! ------------------------------------------------------------------
   ! stage 4: the triangles of one edge (the case table)
   ! ------------------------------------------------------------------

   pure subroutine put_tri(px, py, af, lab, g, x1, y1, x2, y2, x3, y3, f, label)
      !! Store one triangle in slot `g` (vertices, positive area factor with
      !! its sign, source-cell label).
      !$acc routine seq
      real(wp), intent(inout) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), af(REMAP_NSLOT)
      integer, intent(inout) :: lab(REMAP_NSLOT)
      integer, intent(in) :: g, label
      real(wp), intent(in) :: x1, y1, x2, y2, x3, y3, f

      px(1, g) = x1
      py(1, g) = y1
      px(2, g) = x2
      py(2, g) = y2
      px(3, g) = x3
      py(3, g) = y3
      af(g) = f
      lab(g) = label
   end subroutine put_tri

   pure function tri_area(px, py, af, g) result(a)
      !! `0.5 cross((v2-v1), (v3-v1)) * af` of slot `g` (the reference's expression).
      !$acc routine seq
      real(wp), intent(in) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), af(REMAP_NSLOT)
      integer, intent(in) :: g
      real(wp) :: a

      a = 0.5_wp*((px(2, g) - px(1, g))*(py(3, g) - py(1, g)) - (py(2, g) - py(1, g))*(px(3, g) - px(1, g)))*af(g)
   end function tri_area

   pure subroutine remap_locate_edge(iedge, efa, puny, dlx, dly, drx, dry, fc, f_tl, f_bl, f_tr, f_br, f_l, f_r, &
                                     edgearea, px, py, ar, sc, icase, areasum, ibits)
      !! The case table for ONE edge, on scalars.  `iedge` is
      !! `REMAP_EDGE_NORTH` or `REMAP_EDGE_EAST`, `efa` 1 (adjusted to
      !! `edgearea`) or 0.  `(dlx, dly)`, `(drx, dry)` are the SCALED
      !! displacements of the left and right corners in the EDGE frame
      !! (north: `dp/dxU`, `dp/dyU` of `U(i-1,j)`, `U(i,j)`; east: `(-dp_y/dyU,
      !! dp_x/dxU)` of `U(i,j)`, `U(i,j-1)`).  Area factors: `fc` the edge
      !! cell's own (`narea` / `earea`), `f_tl, f_bl, f_tr, f_br` the
      !! other-edge cell areas of the four side cells (`earea` of
      !! `(i-1,j+1), (i-1,j), (i,j+1), (i,j)` for a north edge; `narea` of
      !! `(i+1,j), (i,j), (i+1,j-1), (i,j-1)` for an east edge), `f_l, f_r` the
      !! same-edge areas of the left / right neighbour (`narea(i-1,j)`,
      !! `narea(i+1,j)`; `earea(i,j+1)`, `earea(i,j-1)`).  Outputs: vertices
      !! `px, py(3, 6)` in the SOURCE cell frame, signed areas `ar(6)`, source
      !! codes `sc(6)` (0 = unused), the branch code `icase`, the area sum
      !! `areasum` and the error bits `ibits`.  On `ibits /= 0` the other
      !! outputs are unspecified.
      !$acc routine seq
      integer, value :: iedge, efa
      real(wp), intent(in) :: puny, dlx, dly, drx, dry, fc, f_tl, f_bl, f_tr, f_br, f_l, f_r, edgearea
      real(wp), intent(out) :: areasum
      real(wp), intent(out) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      integer(int8), intent(out) :: sc(REMAP_NSLOT)
      integer(int32), intent(out) :: icase, ibits
      real(wp) :: af(REMAP_NSLOT)
      integer :: lab(REMAP_NSLOT)
      real(wp) :: xcl, ycl, xcr, ycr, xdl, ydl, xdr, ydr, xdm, ydm, xil, yil, xir, yir, md, xic, yic
      real(wp) :: xicl, yicl, xicr, yicr, mdl, mdr, area1, area2, area3, area4, area_c, w1, w2, wx, wy, a
      integer :: g, v, k, ic_l, ic_r, ic_c, ic_e, di, dj
      integer(int32) :: bits

      ibits = 0_int32
      icase = 0_int32
      areasum = 0.0_wp
      do g = 1, REMAP_NSLOT
         ar(g) = 0.0_wp
         sc(g) = 0_int8
         af(g) = 0.0_wp
         lab(g) = 0
         do v = 1, 3
            px(v, g) = 0.0_wp
            py(v, g) = 0.0_wp
         end do
      end do

      if (.not. (ieee_is_finite(dlx) .and. ieee_is_finite(dly) .and. ieee_is_finite(drx) .and. ieee_is_finite(dry) &
                 .and. ieee_is_finite(edgearea) .and. ieee_is_finite(fc) .and. ieee_is_finite(f_tl) &
                 .and. ieee_is_finite(f_bl) .and. ieee_is_finite(f_tr) .and. ieee_is_finite(f_br) &
                 .and. ieee_is_finite(f_l) .and. ieee_is_finite(f_r))) then
         ibits = ICE_ERR_NONFINITE
         return
      end if
      if (abs(dlx) > 1.0_wp .or. abs(dly) > 1.0_wp .or. abs(drx) > 1.0_wp .or. abs(dry) > 1.0_wp) then
         ibits = ICE_ERR_TRANSPORT_CFL
         return
      end if

      xcl = -0.5_wp
      ycl = 0.0_wp
      xcr = 0.5_wp
      ycr = 0.0_wp
      xdl = xcl + dlx
      ydl = ycl + dly
      xdr = xcr + drx
      ydr = ycr + dry
      if (.not. (xdr > xdl)) then
         ! departure points crossed (or coincident): the region has no orientation
         ibits = ICE_ERR_TRANSPORT_CFL
         return
      end if
      if (efa /= 0 .and. .not. (fc > 0.0_wp)) then
         ibits = ICE_ERR_NONFINITE
         return
      end if

      xdm = 0.5_wp*(xdr + xdl)
      ydm = 0.5_wp*(ydr + ydl)

      xil = xcl
      yil = (xcl*(ydm - ydl) + xdm*ydl - xdl*ydm)/(xdm - xdl)
      xir = xcr
      yir = (xcr*(ydr - ydm) - xdm*ydr + xdr*ydm)/(xdr - xdm)
      md = (ydr - ydl)/(xdr - xdl)
      if (abs(md) > puny) then
         xic = xdl - ydl/md
      else
         xic = 0.0_wp
      end if
      yic = 0.0_wp
      xicl = xic
      yicl = yic
      xicr = xic
      yicr = yic

      ic_l = 0
      ic_r = 0
      ic_c = 0
      ic_e = 0

      ! ---- left side: TL, BL, TL1+BL1, TL2+BL2 ----
      if (yil > 0.0_wp .and. xdl < xcl .and. ydl >= 0.0_wp) then
         ic_l = 1
         call put_tri(px, py, af, lab, 1, xcl, ycl, xil, yil, xdl, ydl, -f_tl, LTL)
      else if (yil < 0.0_wp .and. xdl < xcl .and. ydl < 0.0_wp) then
         ic_l = 2
         call put_tri(px, py, af, lab, 1, xcl, ycl, xdl, ydl, xil, yil, f_bl, LBL)
      else if (yil < 0.0_wp .and. xdl < xcl .and. ydl >= 0.0_wp) then
         ic_l = 3
         call put_tri(px, py, af, lab, 1, xcl, ycl, xdl, ydl, xic, yic, f_l, LTL)
         call put_tri(px, py, af, lab, 3, xcl, ycl, xic, yic, xil, yil, f_bl, LBL)
      else if (yil > 0.0_wp .and. xdl < xcl .and. ydl < 0.0_wp) then
         ic_l = 4
         call put_tri(px, py, af, lab, 3, xcl, ycl, xil, yil, xic, yic, -f_tl, LTL)
         call put_tri(px, py, af, lab, 1, xcl, ycl, xic, yic, xdl, ydl, -f_l, LBL)
      end if

      ! ---- right side: TR, BR, TR1+BR1, TR2+BR2 ----
      if (yir > 0.0_wp .and. xdr >= xcr .and. ydr >= 0.0_wp) then
         ic_r = 1
         call put_tri(px, py, af, lab, 2, xcr, ycr, xdr, ydr, xir, yir, -f_tr, LTR)
      else if (yir < 0.0_wp .and. xdr >= xcr .and. ydr < 0.0_wp) then
         ic_r = 2
         call put_tri(px, py, af, lab, 2, xcr, ycr, xir, yir, xdr, ydr, f_br, LBR)
      else if (yir < 0.0_wp .and. xdr >= xcr .and. ydr >= 0.0_wp) then
         ic_r = 3
         call put_tri(px, py, af, lab, 2, xcr, ycr, xic, yic, xdr, ydr, f_r, LTR)
         call put_tri(px, py, af, lab, 3, xcr, ycr, xir, yir, xic, yic, f_br, LBR)
      else if (yir > 0.0_wp .and. xdr >= xcr .and. ydr < 0.0_wp) then
         ic_r = 4
         call put_tri(px, py, af, lab, 3, xcr, ycr, xic, yic, xir, yir, -f_tr, LTR)
         call put_tri(px, py, af, lab, 2, xcr, ycr, xdr, ydr, xic, yic, -f_r, LBR)
      end if

      ! ---- departure points outside the centre column move to the side lines ----
      if (xdl < xcl) then
         xdl = xil
         ydl = yil
      end if
      if (xdr > xcr) then
         xdr = xir
         ydr = yir
      end if

      ! ---- edge flux adjustment: shift the midpoint so the region has the prescribed area ----
      if (efa /= 0) then
         area1 = tri_area(px, py, af, 1)
         area2 = tri_area(px, py, af, 2)
         area3 = tri_area(px, py, af, 3)
         if (ydl*ydr >= 0.0_wp) then
            ! both departure points on one side of the axis
            ic_e = 1
            area_c = edgearea - area1 - area2 - area3
            w1 = 2.0_wp*area_c/fc + (xdr - xcl)*ydl + (xcr - xdl)*ydr
            w2 = (xdr - xdl)**2 + (ydr - ydl)**2
            w1 = w1/w2
            xdm = xdm + (ydr - ydl)*w1
            ydm = ydm - (xdr - xdl)*w1
            mdl = (ydm - ydl)/(xdm - xdl)
            mdr = (ydr - ydm)/(xdr - xdm)
            if (abs(mdl) > puny) then
               xicl = xdl - ydl/mdl
            else
               xicl = 0.0_wp
            end if
            yicl = 0.0_wp
            if (abs(mdr) > puny) then
               xicr = xdr - ydr/mdr
            else
               xicr = 0.0_wp
            end if
            yicr = 0.0_wp
         else if (xic < 0.0_wp .and. xic > xcl) then
            ! different sides, axis crossing left of centre: pin ICL = IC
            ic_e = 2
            xicl = xic
            yicl = yic
            xdm = 0.5_wp*(xdr + xicl)
            ydm = 0.5_wp*ydr
            area4 = 0.5_wp*(xcl - xic)*ydl*fc
            area_c = edgearea - area1 - area2 - area3 - area4
            w1 = 2.0_wp*area_c/fc + (xcr - xic)*ydr
            w2 = (xdr - xic)**2 + ydr**2
            w1 = w1/w2
            xdm = xdm + ydr*w1
            ydm = ydm - (xdr - xic)*w1
            mdr = (ydr - ydm)/(xdr - xdm)
            if (abs(mdr) > puny) then
               xicr = xdr - ydr/mdr
            else
               xicr = 0.0_wp
            end if
            yicr = 0.0_wp
         else if (xic >= 0.0_wp .and. xic < xcr) then
            ! different sides, axis crossing right of centre: pin ICR = IC
            ic_e = 3
            xicr = xic
            yicr = yic
            xdm = 0.5_wp*(xicr + xdl)
            ydm = 0.5_wp*ydl
            area4 = 0.5_wp*(xic - xcr)*ydr*fc
            area_c = edgearea - area1 - area2 - area3 - area4
            w1 = 2.0_wp*area_c/fc + (xic - xcl)*ydl
            w2 = (xic - xdl)**2 + ydl**2
            w1 = w1/w2
            xdm = xdm - ydl*w1
            ydm = ydm - (xic - xdl)*w1
            mdl = (ydm - ydl)/(xdm - xdl)
            if (abs(mdl) > puny) then
               xicl = xdl - ydl/mdl
            else
               xicl = 0.0_wp
            end if
            yicl = 0.0_wp
         else
            ic_e = 4
         end if
      end if

      ! ---- centre column: TC and BC ----
      if (ydl >= 0.0_wp .and. ydr >= 0.0_wp .and. ydm >= 0.0_wp) then
         ic_c = 1
         call put_tri(px, py, af, lab, 4, xcl, ycl, xcr, ycr, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xdl, ydl, xdr, ydr, xdm, ydm, -fc, LTC)
      else if (ydl >= 0.0_wp .and. ydr >= 0.0_wp .and. ydm < 0.0_wp) then
         ic_c = 2
         call put_tri(px, py, af, lab, 4, xcl, ycl, xicl, yicl, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xicr, yicr, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xicr, yicr, xicl, yicl, xdm, ydm, fc, LBC)
      else if (ydl <= 0.0_wp .and. ydr <= 0.0_wp .and. ydm <= 0.0_wp) then
         ic_c = 3
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xcr, ycr, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdl, ydl, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xdl, ydl, xdm, ydm, xdr, ydr, fc, LBC)
      else if (ydl <= 0.0_wp .and. ydr <= 0.0_wp .and. ydm > 0.0_wp) then
         ic_c = 4
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xicl, yicl, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xicr, yicr, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xicl, yicl, xicr, yicr, xdm, ydm, -fc, LTC)
      else if (ydl > 0.0_wp .and. ydr < 0.0_wp .and. xic >= 0.0_wp .and. ydm >= 0.0_wp) then
         ic_c = 5
         call put_tri(px, py, af, lab, 4, xcl, ycl, xicr, yicr, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xicr, yicr, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xdl, ydl, xicr, yicr, xdm, ydm, -fc, LTC)
      else if (ydl > 0.0_wp .and. ydr < 0.0_wp .and. xic >= 0.0_wp .and. ydm < 0.0_wp) then
         ic_c = 6
         call put_tri(px, py, af, lab, 4, xcl, ycl, xicl, yicl, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xicr, yicr, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xicr, yicr, xicl, yicl, xdm, ydm, fc, LBC)
      else if (ydl > 0.0_wp .and. ydr < 0.0_wp .and. xic < 0.0_wp .and. ydm < 0.0_wp) then
         ic_c = 7
         call put_tri(px, py, af, lab, 4, xcl, ycl, xicl, yicl, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xicl, yicl, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xdr, ydr, xicl, yicl, xdm, ydm, fc, LBC)
      else if (ydl > 0.0_wp .and. ydr < 0.0_wp .and. xic < 0.0_wp .and. ydm >= 0.0_wp) then
         ic_c = 8
         call put_tri(px, py, af, lab, 4, xcl, ycl, xicl, yicl, xdl, ydl, -fc, LTC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xicr, yicr, xdr, ydr, fc, LBC)
         call put_tri(px, py, af, lab, 6, xicl, yicl, xicr, yicr, xdm, ydm, -fc, LTC)
      else if (ydl < 0.0_wp .and. ydr > 0.0_wp .and. xic < 0.0_wp .and. ydm >= 0.0_wp) then
         ic_c = 9
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xicl, yicl, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xicl, yicl, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xicl, yicl, xdr, ydr, xdm, ydm, -fc, LTC)
      else if (ydl < 0.0_wp .and. ydr > 0.0_wp .and. xic < 0.0_wp .and. ydm < 0.0_wp) then
         ic_c = 10
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xicl, yicl, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xicr, yicr, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xicr, yicr, xicl, yicl, xdm, ydm, fc, LBC)
      else if (ydl < 0.0_wp .and. ydr > 0.0_wp .and. xic >= 0.0_wp .and. ydm < 0.0_wp) then
         ic_c = 11
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xicr, yicr, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xicr, yicr, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xicr, yicr, xdl, ydl, xdm, ydm, fc, LBC)
      else if (ydl < 0.0_wp .and. ydr > 0.0_wp .and. xic >= 0.0_wp .and. ydm >= 0.0_wp) then
         ic_c = 12
         call put_tri(px, py, af, lab, 4, xcl, ycl, xdl, ydl, xicl, yicl, fc, LBC)
         call put_tri(px, py, af, lab, 5, xcr, ycr, xdr, ydr, xicr, yicr, -fc, LTC)
         call put_tri(px, py, af, lab, 6, xicl, yicl, xicr, yicr, xdm, ydm, -fc, LTC)
      end if

      icase = int(ic_l + 5*ic_r + 25*ic_c + 325*ic_e, int32)

      ! ---- areas: sign by the area factor, drop what is smaller than 1e-16 cell areas ----
      bits = 0_int32
      do g = 1, REMAP_NSLOT
         if (lab(g) == 0) cycle
         a = tri_area(px, py, af, g)
         if (abs(a) < EPS16*fc) then
            lab(g) = 0
         else
            ar(g) = a
            areasum = areasum + a
         end if
      end do
      if (efa /= 0) then
         if (abs(areasum - edgearea) > EPS13*fc) bits = ior(bits, ICE_ERR_TRANSPORT_CONSERV)
      end if

      ! ---- vertices into the source cell's own frame ----
      do g = 1, REMAP_NSLOT
         if (lab(g) == 0) then
            ar(g) = 0.0_wp
            cycle
         end if
         di = SH_I(lab(g), iedge)
         dj = SH_J(lab(g), iedge)
         sc(g) = SRC_CODE(di, dj)
         do k = 1, 3
            if (iedge == REMAP_EDGE_NORTH) then
               wx = px(k, g) - real(di, wp)
               wy = py(k, g) + 0.5_wp - real(dj, wp)
            else
               ! east edge: quarter turn back
               wx = py(k, g) + 0.5_wp - real(di, wp)
               wy = -px(k, g) - real(dj, wp)
            end if
            px(k, g) = wx
            py(k, g) = wy
            if (.not. (ieee_is_finite(wx) .and. ieee_is_finite(wy))) then
               bits = ior(bits, ICE_ERR_NONFINITE)
            else if (abs(wx) > REMAP_VERTEX_MAX .or. abs(wy) > REMAP_VERTEX_MAX) then
               bits = ior(bits, ICE_ERR_TRANSPORT_CFL)
            end if
         end do
      end do
      if (.not. ieee_is_finite(areasum)) bits = ior(bits, ICE_ERR_NONFINITE)
      ibits = bits
   end subroutine remap_locate_edge

   pure function scaled_disp(dp, d) result(s)
      !! `dp/d` for a corner length `d > 0`; a non-zero displacement on a corner
      !! without a length maps to 2 (beyond one cell: the edge routine flags it).
      !$acc routine seq
      real(wp), intent(in) :: dp, d
      real(wp) :: s
      if (dp == 0.0_wp) then
         s = 0.0_wp
      else if (d > 0.0_wp) then
         s = dp/d
      else
         s = 2.0_wp
      end if
   end function scaled_disp

   pure subroutine remap_edge_one(i, j, nx, ny, iedge, efa, puny, dlx, dly, drx, dry, fc, f_tl, f_bl, f_tr, f_br, &
                                  f_l, f_r, edgearea, xp, yp, area, src, icase, ierr)
      !! One edge: `remap_locate_edge` on fixed-size locals, then the stores at
      !! `(i, j)` (vertices only of the slots in use).  `edgearea(i,j)` is read
      !! with `efa /= 0` and written otherwise.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny
      integer, value :: iedge, efa
      real(wp), intent(in) :: puny, dlx, dly, drx, dry, fc, f_tl, f_bl, f_tr, f_br, f_l, f_r
      real(wp), intent(inout) :: edgearea(nx, ny)
      real(wp), intent(inout) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      integer(int8), intent(inout) :: src(nx, ny, REMAP_NSLOT)
      integer(int32), intent(inout) :: icase(nx, ny), ierr(nx, ny)
      real(wp) :: asum, ea
      real(wp) :: px(3, REMAP_NSLOT), py(3, REMAP_NSLOT), ar(REMAP_NSLOT)
      integer(int8) :: sc(REMAP_NSLOT)
      integer(int32) :: ic, bits
      integer :: g, v

      ea = 0.0_wp
      if (efa /= 0) ea = edgearea(i, j)
      call remap_locate_edge(iedge, efa, puny, dlx, dly, drx, dry, fc, f_tl, f_bl, f_tr, f_br, f_l, f_r, ea, &
                             px, py, ar, sc, ic, asum, bits)
      if (bits /= 0_int32) then
         do g = 1, REMAP_NSLOT
            src(i, j, g) = 0_int8
            area(i, j, g) = 0.0_wp
         end do
         icase(i, j) = 0_int32
         ierr(i, j) = ior(ierr(i, j), bits)
         return
      end if
      do g = 1, REMAP_NSLOT
         src(i, j, g) = sc(g)
         area(i, j, g) = ar(g)
         if (sc(g) /= 0_int8) then
            do v = 1, 3
               xp(i, j, v, g) = px(v, g)
               yp(i, j, v, g) = py(v, g)
            end do
         end if
      end do
      icase(i, j) = ic
      if (efa == 0) edgearea(i, j) = asum
   end subroutine remap_edge_one

   subroutine remap_locate_north(nx, ny, is, ie, js, je, efa, puny, dpx, dpy, dxU, dyU, earea, narea, edgearea, &
                                 xp, yp, area, src, icase, ierr)
      !! The triangles of the north edges `N(is:ie, js-1:je)`.
      integer, intent(in) :: nx, ny, is, ie, js, je, efa
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: dpx(nx, ny), dpy(nx, ny), dxU(nx, ny), dyU(nx, ny), earea(nx, ny), narea(nx, ny)
      real(wp), intent(inout) :: edgearea(nx, ny)
      real(wp), intent(inout) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      integer(int8), intent(inout) :: src(nx, ny, REMAP_NSLOT)
      integer(int32), intent(inout) :: icase(nx, ny), ierr(nx, ny)
      integer :: i, j, g
      real(wp) :: dlx, dly, drx, dry
      logical :: active

      !$acc data present(dpx, dpy, dxU, dyU, earea, narea, edgearea, xp, yp, area, src, icase, ierr)
      do concurrent(j=js - 1:je, i=is:ie) local(dlx, dly, drx, dry, active, g)
         active = dpx(i - 1, j) /= 0.0_wp .or. dpy(i - 1, j) /= 0.0_wp .or. dpx(i, j) /= 0.0_wp .or. dpy(i, j) /= 0.0_wp
         if (efa /= 0) active = active .or. edgearea(i, j) /= 0.0_wp
         if (active) then
            dlx = scaled_disp(dpx(i - 1, j), dxU(i - 1, j))
            dly = scaled_disp(dpy(i - 1, j), dyU(i - 1, j))
            drx = scaled_disp(dpx(i, j), dxU(i, j))
            dry = scaled_disp(dpy(i, j), dyU(i, j))
            call remap_edge_one(i, j, nx, ny, REMAP_EDGE_NORTH, efa, puny, dlx, dly, drx, dry, narea(i, j), &
                                earea(i - 1, j + 1), earea(i - 1, j), earea(i, j + 1), earea(i, j), &
                                narea(i - 1, j), narea(i + 1, j), edgearea, xp, yp, area, src, icase, ierr)
         else
            do g = 1, REMAP_NSLOT
               src(i, j, g) = 0_int8
               area(i, j, g) = 0.0_wp
            end do
            icase(i, j) = 0_int32
            if (efa == 0) edgearea(i, j) = 0.0_wp
         end if
      end do
      !$acc end data
   end subroutine remap_locate_north

   subroutine remap_locate_east(nx, ny, is, ie, js, je, efa, puny, dpx, dpy, dxU, dyU, earea, narea, edgearea, &
                                xp, yp, area, src, icase, ierr)
      !! The triangles of the east edges `E(is-1:ie, js:je)`: the same table in
      !! the frame turned by 90 degrees (`CL` the NE corner, `CR` the SE one).
      integer, intent(in) :: nx, ny, is, ie, js, je, efa
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: dpx(nx, ny), dpy(nx, ny), dxU(nx, ny), dyU(nx, ny), earea(nx, ny), narea(nx, ny)
      real(wp), intent(inout) :: edgearea(nx, ny)
      real(wp), intent(inout) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      integer(int8), intent(inout) :: src(nx, ny, REMAP_NSLOT)
      integer(int32), intent(inout) :: icase(nx, ny), ierr(nx, ny)
      integer :: i, j, g
      real(wp) :: dlx, dly, drx, dry
      logical :: active

      !$acc data present(dpx, dpy, dxU, dyU, earea, narea, edgearea, xp, yp, area, src, icase, ierr)
      do concurrent(j=js:je, i=is - 1:ie) local(dlx, dly, drx, dry, active, g)
         active = dpx(i, j - 1) /= 0.0_wp .or. dpy(i, j - 1) /= 0.0_wp .or. dpx(i, j) /= 0.0_wp .or. dpy(i, j) /= 0.0_wp
         if (efa /= 0) active = active .or. edgearea(i, j) /= 0.0_wp
         if (active) then
            ! left = NE corner U(i,j), right = SE corner U(i,j-1); frame turned: (x, y) = (-dy, dx)
            dlx = -scaled_disp(dpy(i, j), dyU(i, j))
            dly = scaled_disp(dpx(i, j), dxU(i, j))
            drx = -scaled_disp(dpy(i, j - 1), dyU(i, j - 1))
            dry = scaled_disp(dpx(i, j - 1), dxU(i, j - 1))
            call remap_edge_one(i, j, nx, ny, REMAP_EDGE_EAST, efa, puny, dlx, dly, drx, dry, earea(i, j), &
                                narea(i + 1, j), narea(i, j), narea(i + 1, j - 1), narea(i, j - 1), &
                                earea(i, j + 1), earea(i, j - 1), edgearea, xp, yp, area, src, icase, ierr)
         else
            do g = 1, REMAP_NSLOT
               src(i, j, g) = 0_int8
               area(i, j, g) = 0.0_wp
            end do
            icase(i, j) = 0_int32
            if (efa == 0) edgearea(i, j) = 0.0_wp
         end if
      end do
      !$acc end data
   end subroutine remap_locate_east

   ! ------------------------------------------------------------------
   ! the whole geometry
   ! ------------------------------------------------------------------

   subroutine remap_geometry(ws, is, ie, js, je, dt, efa, puny, u, v, dxU, dyU, HTE, HTN, earea, narea, epm, npm, uvm, ierr, &
                             uc_in, vc_in)
      !! Corner velocities, departure points, edge areas and both triangle
      !! sets from the face velocities (six launches).  `efa /= 0`: the
      !! C-grid mode (the triangles are adjusted to the face flux).  The corner
      !! velocities are the dynamics' masked, area-weighted average unless the
      !! caller passes its own (`uc_in, vc_in`, e.g. the dynamics' `uvelU/vvelU`;
      !! both or neither), in which case that stage is skipped.  Everything
      !! device-resident on an offload build.  See the module docstring for
      !! the ranges and the outputs (`ws%xp_e` ... `ws%icase_n`).
      type(remap_geom_workspace_t), intent(inout) :: ws
      integer, intent(in) :: is, ie, js, je, efa
      real(wp), intent(in) :: dt, puny
      real(wp), intent(in) :: u(ws%nx + 1, ws%ny), v(ws%nx, ws%ny + 1)
      real(wp), intent(in) :: dxU(ws%nx, ws%ny), dyU(ws%nx, ws%ny), HTE(ws%nx, ws%ny), HTN(ws%nx, ws%ny)
      real(wp), intent(in) :: earea(ws%nx, ws%ny), narea(ws%nx, ws%ny)
      real(wp), intent(in) :: epm(ws%nx, ws%ny), npm(ws%nx, ws%ny), uvm(ws%nx, ws%ny)
      integer(int32), intent(inout) :: ierr(ws%nx, ws%ny)
      real(wp), intent(in), optional :: uc_in(ws%nx, ws%ny), vc_in(ws%nx, ws%ny)

      if (present(uc_in) .and. present(vc_in)) then
         call remap_departure_points(ws%nx, ws%ny, is, ie, js, je, dt, uc_in, vc_in, dxU, dyU, ws%dpx, ws%dpy, ierr)
      else
         call remap_corner_velocity(ws%nx, ws%ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm, ws%uN, ws%vE, ws%uc, ws%vc)
         call remap_departure_points(ws%nx, ws%ny, is, ie, js, je, dt, ws%uc, ws%vc, dxU, dyU, ws%dpx, ws%dpy, ierr)
      end if
      call remap_edge_areas(ws%nx, ws%ny, is, ie, js, je, dt, u, v, HTE, HTN, ws%edgearea_e, ws%edgearea_n)
      call remap_locate_east(ws%nx, ws%ny, is, ie, js, je, efa, puny, ws%dpx, ws%dpy, dxU, dyU, earea, narea, &
                             ws%edgearea_e, ws%xp_e, ws%yp_e, ws%area_e, ws%src_e, ws%icase_e, ierr)
      call remap_locate_north(ws%nx, ws%ny, is, ie, js, je, efa, puny, ws%dpx, ws%dpy, dxU, dyU, earea, narea, &
                              ws%edgearea_n, ws%xp_n, ws%yp_n, ws%area_n, ws%src_n, ws%icase_n, ierr)
   end subroutine remap_geometry

end module rdb_ice_remap_geom
