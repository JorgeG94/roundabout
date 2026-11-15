!! Sea-ice incremental remapping, step 3: the flux integrals over the
!! departure triangles, the conservative field update and the advection
!! driver that strings geometry, reconstruction, fluxes and update together
!! (Dukowicz and Baumgardner 2000, J. Comput. Phys. 160, 318; Lipscomb and
!! Hunke 2004, Mon. Wea. Rev. 132, 1341).
module rdb_ice_remap_flux
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_remap (transport_integrals, triangle_coordinates, update_fields).
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
   ! Elizabeth C. Hunke (LANL).
   ! Physics: Lipscomb and Hunke (2004, Mon. Wea. Rev. 132 1341); Dukowicz and
   ! Baumgardner (2000, J. Comput. Phys. 160 318); the triangle quadrature is
   ! the four-point cubic formula of Stroud (1971, Approximate Calculation of
   ! Multiple Integrals, section 8.8, formula 3.1).
   !
   !! **What this module does.**  Given the departure geometry
   !! (`rdb_ice_remap_geom`: per edge and triangle slot a signed area, a source
   !! cell and three vertices in that cell's frame) and the limited linear
   !! reconstruction of the area and of the tracer hierarchy
   !! (`rdb_ice_remap_recon`), it
   !!  1. integrates the reconstructed fields over every departure triangle and
   !!     accumulates them into one flux per edge, category and field
   !!     (`remap_flux_east_impl`, `remap_flux_north_impl`);
   !!  2. updates the area and the tracers from the flux divergence
   !!     (`remap_update_impl`);
   !!  3. (`remap_advect`) runs geometry, reconstruction and both of those, a
   !!     category batch at a time.
   !!
   !! **The flux rule** (the geometry module's header, restated).  With
   !! `<R_c, T>` the MEAN of the reconstruction of cell `c` over triangle `T`
   !! at the stored vertices, in `c`'s own frame, UNCLIPPED:
   !! `F_e(i,j) = sum_s area_e(i,j,s) <R_(i+di, j+dj), triangle s>`, the same
   !! for the north edge, and for each field
   !! `delta(field*area)(i,j) = -(F_e(i,j) - F_e(i-1,j) + F_n(i,j) - F_n(i,j-1))`.
   !! The signed area is never replaced by its modulus: a triangle in the cell
   !! the edge leads into carries a negative area (inflow), and a pair of
   !! opposite-sign triangles in one source cell (the "rare" branches of the
   !! case table) may reach a cell beyond its own, where the linear
   !! reconstruction is evaluated by extrapolation; the two contributions
   !! cancel to round-off wherever the extrapolation is the same polynomial.
   !! Each edge's flux is computed ONCE and stored; both cells it joins read
   !! the same number, which is what makes the update conservative.
   !!
   !! **The integrals and the quadrature rule.**  Per triangle with vertices
   !! `P1, P2, P3` and centroid `C`, the rule evaluates at the four points `C`
   !! (weight `-9/16`) and `Q_k = C + (2/5)(P_k - C)` (weight `25/48` each; the
   !! weights sum to 1, so the rule gives the MEAN over the triangle and the
   !! signed triangle area is applied outside).  It is exact for polynomials
   !! up to degree THREE (Stroud 1971; the reference's `integral_order = 3`).
   !! The reconstruction is `m = mc + mx x + my y` and `T = tc + tx x + ty y`,
   !! so with `<.>` the mean over the triangle:
   !!
   !! | product | polynomial degree | integrals formed | exact? |
   !! |---|---|---|---|
   !! | mass `m` | 1 | the six moments `<m>`, `<m x>`, `<m y>`, `<m x^2>`, `<m x y>`, `<m y^2>` | exact (`m` times degree `<= 2`, so degree 3 at most) |
   !! | `m T` (type 1: `depend = 0`) | 2 | `<m T> = tc <m> + tx <m x> + ty <m y>` | quadratic: exact |
   !! | `m P T` (type 2, parent `P`) | 3 | `<m P T> = tc <m P> + tx <m P x> + ty <m P y>` | cubic: exact |
   !!
   !! with, for the parent `P = pc + px x + py y`: `<m P> = pc <m> + px <m x> + py <m y>`,
   !! `<m P x> = pc <m x> + px <m x^2> + py <m x y>` and the `y` analogue.
   !!
   !! The tracer fluxes are formed in hierarchy order (a type-2 tracer reads its
   !! parent's `tc, tx, ty`, never a neighbour's value), so no `ntrace`-sized
   !! per-thread array exists: the six mass moments of each of the (at most six)
   !! slots are the only per-thread storage (six scalars, one slot at a time).  The tests
   !! compare every flux against an independent host integration by exact
   !! monomial formulas in barycentric coordinates.
   !!
   !! **The update** (`remap_update_impl`, the reference's guards).  Per
   !! category `n = 0..ncat`: `m_new = m - (F_e(i,j) - F_e(i-1,j) + F_n(i,j) -
   !! F_n(i,j-1))/A`; a `m_new` in `[-puny, 0)` is set to 0 (round-off), one
   !! below `-puny` is COUNTED (`nneg`; the value stays as computed and the
   !! driver refuses the step; the reference aborts).  A tracer is recovered
   !! where `m_new > 0`: type 1 `T_new = (m T - div F_T/A)/m_new`; type 2
   !! `T_new = (m P T - div F_T/A)/(m_new P_new)` where `|P_new| > 0` (else 0);
   !! everywhere else 0 (the reference's `tm = 0`, which in the packed frame is
   !! the melting point for the snow enthalpy).  Where the net flux of the field
   !! and of its base are exactly zero the OLD value is copied, so an empty
   !! domain and a zero velocity are bitwise no-ops of the packed fields (the
   !! quotient `(m T)/m` is not `T` in the last bit; for a type-2 tracer also where its
   !! parent is zero, so the packed snow enthalpy of a snow-free category is not reset).  Non-finite new values are
   !! counted (`nbad`).  Note the area-tracer of a cell with `0 < m_new <=
   !! puny` is NOT zeroed (as the reference; the upwind kernel zeroes it): such
   !! a sliver is removed by the ridging step's clean-up.
   !!
   !! **Conservation** is to round-off, not bitwise: one division by the cell
   !! area per update and the quadrature rounding in the flux, which both
   !! neighbours share.  **Monotonicity** is the property of the limited
   !! reconstruction (`rdb_ice_remap_recon`), checked by the driver under
   !! `transport_check` with the 5x5 bound of `rdb_ice_transport_checks`.
   !!
   !! **GPU shape (category-batched buffers).**  The geometry is
   !! category-independent and built once per step.  The reconstruction
   !! (`rdb_ice_remap_recon_batch`, no stored masks) and the fluxes of a BATCH of
   !! `nb` category indices (`nlo..nhi`, `0` = open water) live in
   !! `mflx_{e,n}(nx, ny, nb)` and `mtflx_{e,n}(nx, ny, ntrace, nb)`; per batch
   !! there are FIVE launches: the reconstruction, east edges `do concurrent (b, j, i)` over
   !! `E(is-1:ie, js:je)`, north edges over `N(is:ie, js-1:je)` (two launches,
   !! not one over the union, so no edge guard sits in the loop body), then the
   !! update and the ledger (what the guards create or drop) over the owned cells.  The face routine is `pure` +
   !! `!$acc routine seq` with every index `integer, value`.  With `nb = ncat+1`
   !! (the default while it fits) one batch does everything; a smaller `nb`
   !! gives bit-identical results with a flux buffer `nb/(ncat+1)` the size
   !! (the plan's category-batched fallback for a 16 GB card at 1440 x 1152 x
   !! 5 categories: see `bytes` and the driver's `remap_flux_bytes_max`); the budget counts the reconstruction too.
   !!
   !! **Ranges and halos.**  Owned window `is:ie, js:je`; the packed input
   !! valid on `is-2:ie+2, js-2:je+2` (two rings), the geometry's inputs on
   !! its documented windows (three ghost cells), the reconstruction on
   !! `is-1:ie+1, js-1:je+1`, the fluxes on the edges of the owned window, the
   !! output on the owned cells only.  `is >= 4`, `ie <= nx - 3` (likewise in
   !! `j`) is therefore the minimum (`REMAP_CFG_RANGE`).
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, ICE_ERR_TRANSPORT_CONSERV
   use rdb_ice_remap_geom, only: remap_geom_workspace_t, remap_geometry, remap_src_di, remap_src_dj, REMAP_NSLOT
   use rdb_ice_remap_recon, only: recon_table_ok
   use rdb_ice_remap_recon_batch, only: remap_recon_batch_t, remap_recon_batch, remap_input_check, RECON_CFG_OK, &
                                        RECON_CFG_TABLE, RECON_CFG_RANGE, RECON_CFG_MASK
   use rdb_ice_transport_pack, only: transport_workspace_t
   use rdb_ice_transport_seam, only: transport_seam_t, seam_sum_int, seam_pair_fold_flux
   implicit none
   private

   public :: remap_flux_workspace_t
   public :: remap_advect, remap_flux_integrals, remap_update
   public :: remap_flux_east_impl, remap_flux_north_impl, remap_update_impl, remap_face_flux
   public :: remap_geometry_errors_impl
   public :: remap_flux_bytes_per_index, remap_default_nbatch
   public :: REMAP_CFG_OK, REMAP_CFG_WORKSPACE, REMAP_CFG_TABLE, REMAP_CFG_RANGE, REMAP_CFG_MASK
   public :: REMAP_QUAD_W0, REMAP_QUAD_W1, REMAP_QUAD_FRAC
   public :: remap_ledger_t, remap_ledger_impl, REMAP_TGROUP

   integer, parameter :: NGV = 3
      !! Counts reduced together over the ranks at one decision point.
   integer, parameter :: REMAP_TGROUP = 9
      !! Tracers accumulated together in registers by the face flux kernel.

   integer, parameter :: REMAP_CFG_OK = 0
      !! Configuration accepted.
   integer, parameter :: REMAP_CFG_WORKSPACE = 1
      !! A workspace is not initialised, or the extents / category counts / batch differ.
   integer, parameter :: REMAP_CFG_TABLE = 2
      !! The tracer table has a tracer the reconstruction does not carry (not type 1, or type 2 under a type-1 parent with dependents).
   integer, parameter :: REMAP_CFG_RANGE = 3
      !! The window leaves no room for the three ghost cells.
   integer, parameter :: REMAP_CFG_MASK = 4
      !! The land mask holds a value other than 0 and 1.

   real(wp), parameter :: REMAP_QUAD_W0 = -9.0_wp/16.0_wp
      !! Weight of the centroid in the four-point cubic rule (Stroud 1971).
   real(wp), parameter :: REMAP_QUAD_W1 = 25.0_wp/48.0_wp
      !! Weight of each of the three outer points.
   real(wp), parameter :: REMAP_QUAD_FRAC = 0.4_wp
      !! The outer points lie this fraction of the way from the centroid to the vertices.

   type :: remap_ledger_t
      !! What the update's guards create or drop, for the caller's conservation ledger (sums over the owned
      !! cells and the categories of one call; the tracer contents are in the PACKED frame, times the cell
      !! area, so `drop_esno` includes the `rhos Lfresh` shift of the snow enthalpy).
      integer :: n_clip = 0
         !! Cells whose new area fell in `[-puny, 0)` and was set to 0.
      real(wp) :: clip_area = 0.0_wp
         !! Area CREATED by those clips [m^2] (`sum(-m_new) A`, `>= 0`).
      real(wp) :: drop_vice = 0.0_wp
         !! Ice volume content `m hi` left in cells whose new area is 0 (or whose parent is 0) [m^3].
      real(wp) :: drop_vsno = 0.0_wp
         !! Snow volume likewise [m^3].
      real(wp) :: drop_eice = 0.0_wp
         !! Ice enthalpy content `m hi q` dropped likewise [J].
      real(wp) :: drop_esno = 0.0_wp
         !! Snow enthalpy content (packed frame) dropped likewise [J].
   end type remap_ledger_t

   type :: remap_flux_workspace_t
      !! Persistent flux buffers of one category batch (layout: module
      !! docstring) and the geometry error word.  OWNERSHIP: the transport
      !! driver slot owns one; `init` allocates (zeroed), `enter_data` maps
      !! the component arrays ONCE (never the aggregate), `exit_data` unmaps,
      !! `destroy` frees.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      integer :: ncat = 0
      integer :: ntrace = 0
      integer :: nbatch = 0
         !! Category indices per batch (`1..ncat+1`).
      real(wp), allocatable :: mflx_e(:, :, :), mflx_n(:, :, :)
         !! `(nx, ny, nbatch)`: area flux through the east / north edge of the cell.
      real(wp), allocatable :: mtflx_e(:, :, :, :), mtflx_n(:, :, :, :)
         !! `(nx, ny, ntrace, nbatch)`: tracer-product flux.
      integer(int32), allocatable :: gerr(:, :)
         !! `(nx, ny)` error word the geometry ORs into; zero between calls.
   contains
      procedure :: init => flux_ws_init
      procedure :: enter_data => flux_ws_enter_data
      procedure :: exit_data => flux_ws_exit_data
      procedure :: destroy => flux_ws_destroy
      procedure :: bytes => flux_ws_bytes
   end type remap_flux_workspace_t

contains

   ! ------------------------------------------------------------------
   ! workspace
   ! ------------------------------------------------------------------

   pure function remap_flux_bytes_per_index(nx, ny, ntrace) result(b)
      !! Bytes of the flux buffers AND the reconstruction buffers per category index of a batch.
      integer, intent(in) :: nx, ny, ntrace
      integer(int64) :: b
      b = int(4 + 5*ntrace, int64)*int(nx, int64)*int(ny, int64)*8_int64
   end function remap_flux_bytes_per_index

   pure function remap_default_nbatch(nx, ny, ncat, ntrace, bytes_max) result(nb)
      !! Largest batch (in category indices, `1..ncat+1`) whose flux and reconstruction buffers
      !! stay within `bytes_max`; at least 1.
      integer, intent(in) :: nx, ny, ncat, ntrace
      integer(int64), intent(in) :: bytes_max
      integer :: nb
      integer(int64) :: per

      per = remap_flux_bytes_per_index(nx, ny, ntrace)
      nb = int(min(int(ncat + 1, int64), max(1_int64, bytes_max/max(per, 1_int64))))
   end function remap_default_nbatch

   subroutine flux_ws_init(this, nx, ny, ncat, ntrace, nbatch, ierr)
      !! Allocate (zeroed).  `ierr /= 0`: invalid extents or batch, nothing allocated.
      class(remap_flux_workspace_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, ntrace, nbatch
      integer, intent(out) :: ierr

      ierr = 0
      if (nx < 3 .or. ny < 3 .or. ncat < 1 .or. ntrace < 1 .or. nbatch < 1 .or. nbatch > ncat + 1) then
         ierr = REMAP_CFG_WORKSPACE
         return
      end if
      call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%ntrace = ntrace
      this%nbatch = nbatch
      allocate (this%mflx_e(nx, ny, nbatch), this%mflx_n(nx, ny, nbatch))
      allocate (this%mtflx_e(nx, ny, ntrace, nbatch), this%mtflx_n(nx, ny, ntrace, nbatch))
      allocate (this%gerr(nx, ny))
      this%mflx_e = 0.0_wp
      this%mflx_n = 0.0_wp
      this%mtflx_e = 0.0_wp
      this%mtflx_n = 0.0_wp
      this%gerr = 0_int32
      this%is_init = .true.
   end subroutine flux_ws_init

   subroutine flux_ws_enter_data(this)
      !! Map the component arrays once (`copyin`: the host zeros).
      class(remap_flux_workspace_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%mflx_e, this%mflx_n)
      !$acc enter data copyin(this%mtflx_e, this%mtflx_n)
      !$acc enter data copyin(this%gerr)
      this%mapped = .true.
   end subroutine flux_ws_enter_data

   subroutine flux_ws_exit_data(this)
      class(remap_flux_workspace_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%mflx_e, this%mflx_n)
      !$acc exit data delete(this%mtflx_e, this%mtflx_n)
      !$acc exit data delete(this%gerr)
      this%mapped = .false.
   end subroutine flux_ws_exit_data

   subroutine flux_ws_destroy(this)
      class(remap_flux_workspace_t), intent(inout) :: this
      call this%exit_data()
      if (allocated(this%mflx_e)) deallocate (this%mflx_e)
      if (allocated(this%mflx_n)) deallocate (this%mflx_n)
      if (allocated(this%mtflx_e)) deallocate (this%mtflx_e)
      if (allocated(this%mtflx_n)) deallocate (this%mtflx_n)
      if (allocated(this%gerr)) deallocate (this%gerr)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%ntrace = 0
      this%nbatch = 0
   end subroutine flux_ws_destroy

   pure function flux_ws_bytes(this) result(b)
      !! Bytes of the buffers (host copy; the device mirror is the same size).
      class(remap_flux_workspace_t), intent(in) :: this
      integer(int64) :: b
      b = 0_int64
      if (allocated(this%mflx_e)) b = b + int(size(this%mflx_e), int64)*8_int64
      if (allocated(this%mflx_n)) b = b + int(size(this%mflx_n), int64)*8_int64
      if (allocated(this%mtflx_e)) b = b + int(size(this%mtflx_e), int64)*8_int64
      if (allocated(this%mtflx_n)) b = b + int(size(this%mtflx_n), int64)*8_int64
      if (allocated(this%gerr)) b = b + int(size(this%gerr), int64)*4_int64
   end function flux_ws_bytes

   ! ------------------------------------------------------------------
   ! face kernel
   ! ------------------------------------------------------------------

   pure subroutine remap_slot_moments(x1, x2, x3, y1, y2, y3, mcv, mxv, myv, m1s, mxs, mys, mxx, mxy, myy)
      !! The six mass moments `<m>, <m x>, <m y>, <m x^2>, <m x y>, <m y^2>` of the triangle `(x, y)` for the
      !! mass profile `mcv + mxv x + myv y` by the four-point rule (weights `-9/16`, `25/48` x 3 at the centroid
      !! and 2/5 of the way to the vertices), as MEANS over the triangle.
      !$acc routine seq
      real(wp), intent(in) :: x1, x2, x3, y1, y2, y3, mcv, mxv, myv
      real(wp), intent(out) :: m1s, mxs, mys, mxx, mxy, myy
      real(wp) :: cx, cy, qx1, qx2, qx3, qy1, qy2, qy3, m0, m1, m2, m3, w0, w1, w2, w3

      cx = (x1 + x2 + x3)/3.0_wp
      cy = (y1 + y2 + y3)/3.0_wp
      qx1 = REMAP_QUAD_FRAC*x1 + (1.0_wp - REMAP_QUAD_FRAC)*cx
      qx2 = REMAP_QUAD_FRAC*x2 + (1.0_wp - REMAP_QUAD_FRAC)*cx
      qx3 = REMAP_QUAD_FRAC*x3 + (1.0_wp - REMAP_QUAD_FRAC)*cx
      qy1 = REMAP_QUAD_FRAC*y1 + (1.0_wp - REMAP_QUAD_FRAC)*cy
      qy2 = REMAP_QUAD_FRAC*y2 + (1.0_wp - REMAP_QUAD_FRAC)*cy
      qy3 = REMAP_QUAD_FRAC*y3 + (1.0_wp - REMAP_QUAD_FRAC)*cy
      m0 = REMAP_QUAD_W0*(mcv + cx*mxv + cy*myv)
      m1 = REMAP_QUAD_W1*(mcv + qx1*mxv + qy1*myv)
      m2 = REMAP_QUAD_W1*(mcv + qx2*mxv + qy2*myv)
      m3 = REMAP_QUAD_W1*(mcv + qx3*mxv + qy3*myv)
      m1s = ((m0 + m1) + m2) + m3
      w0 = m0*cx
      w1 = m1*qx1
      w2 = m2*qx2
      w3 = m3*qx3
      mxs = ((w0 + w1) + w2) + w3
      mxx = ((w0*cx + w1*qx1) + w2*qx2) + w3*qx3
      mxy = ((w0*cy + w1*qy1) + w2*qy2) + w3*qy3
      w0 = m0*cy
      w1 = m1*qy1
      w2 = m2*qy2
      w3 = m3*qy3
      mys = ((w0 + w1) + w2) + w3
      myy = ((w0*cy + w1*qy1) + w2*qy2) + w3*qy3
   end subroutine remap_slot_moments

   pure subroutine remap_face_flux(i, j, n, b, nx, ny, ncat, ntrace, nb, puny, tracer_type, depend, src, xp, yp, area, &
                                   aim, mx, my, tc, tx, ty, mflx, mtflx)
      !! Flux of category `n` (batch slot `b`) through the edge `(i, j)` of ONE
      !! direction: the integrals over its triangle slots (slot outer loop, every
      !! tracer inside it), accumulated in slot order into `mflx(i,j,b)` and (for `n >= 1`) `mtflx(i,j,nt,b)`.  Writes
      !! nothing else.  Indices are passed BY VALUE (nvfortran `do concurrent`
      !! index by-reference bug).
      !$acc routine seq
      integer, value :: i, j, n, b
      integer, intent(in) :: nx, ny, ncat, ntrace, nb
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      integer(int8), intent(in) :: src(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), mx(nx, ny, nb), my(nx, ny, nb)
      real(wp), intent(in) :: tc(nx, ny, ntrace, nb), tx(nx, ny, ntrace, nb), ty(nx, ny, ntrace, nb)
      real(wp), intent(inout) :: mflx(nx, ny, nb), mtflx(nx, ny, ntrace, nb)
      integer :: s, nt, p, i2, j2, ig, ngrp, nt0, g
      real(wp) :: ar, m1s, mxs, mys, mxx, mxy, myy, mcv, mxv, myv
      real(wp) :: tcv, txv, tyv, tcp, txp, typ, mtp, mtxp, mtyp, mt, fm
      real(wp) :: acc(REMAP_TGROUP)

      ! TRACER GROUPS of `REMAP_TGROUP`: inside a group the slots are the outer loop (the six mass moments of a
      ! slot live only inside its pass, the group's fluxes accumulate in `REMAP_TGROUP` scalars-in-registers
      ! from 0 in slot order 1..6, the same sum order as a tracer-outer loop), and the group is stored once.
      ! The mass flux is formed in the first group.  (A slot-outer loop over ALL tracers read-modified-wrote
      ! `mtflx` per slot and per tracer and was 45 % slower; a tracer-outer loop needs the moments of all
      ! slots at once, 255 registers.)
      fm = 0.0_wp
      ngrp = 1
      if (n >= 1) ngrp = (ntrace + REMAP_TGROUP - 1)/REMAP_TGROUP
      !$acc loop seq
      do ig = 1, ngrp
         nt0 = (ig - 1)*REMAP_TGROUP
         acc = 0.0_wp
         !$acc loop seq
         do s = 1, REMAP_NSLOT
            if (src(i, j, s) == 0_int8) cycle
            i2 = i + remap_src_di(src(i, j, s))
            j2 = j + remap_src_dj(src(i, j, s))
            ar = area(i, j, s)
            mcv = aim(i2, j2, n)
            mxv = mx(i2, j2, b)
            myv = my(i2, j2, b)
            ! the centre value is the area where it exceeds `puny` and 0 otherwise (`mc`), and the gradients are
            ! read through it: a cell with `mc = 0` has none.  The reconstruction buffers hold the batch's
            ! categories only: batch slot `b`
            if (.not. (mcv > puny)) then
               mcv = 0.0_wp
               mxv = 0.0_wp
               myv = 0.0_wp
            end if
            call remap_slot_moments(xp(i, j, 1, s), xp(i, j, 2, s), xp(i, j, 3, s), yp(i, j, 1, s), yp(i, j, 2, s), &
                                    yp(i, j, 3, s), mcv, mxv, myv, m1s, mxs, mys, mxx, mxy, myy)
            if (ig == 1) fm = fm + ar*m1s
            if (n >= 1) then
               !$acc loop seq
               do g = 1, REMAP_TGROUP
                  nt = nt0 + g
                  if (nt > ntrace) exit
                  tcv = tc(i2, j2, nt, b)
                  txv = tx(i2, j2, nt, b)
                  tyv = ty(i2, j2, nt, b)
                  if (tracer_type(nt) == 1) then
                     mt = m1s*tcv + mxs*txv + mys*tyv
                  else
                     p = depend(nt)
                     tcp = tc(i2, j2, p, b)
                     txp = tx(i2, j2, p, b)
                     typ = ty(i2, j2, p, b)
                     mtp = m1s*tcp + mxs*txp + mys*typ
                     mtxp = mxs*tcp + mxx*txp + mxy*typ
                     mtyp = mys*tcp + mxy*txp + myy*typ
                     mt = mtp*tcv + mtxp*txv + mtyp*tyv
                  end if
                  acc(g) = acc(g) + ar*mt
               end do
            end if
         end do
         if (n >= 1) then
            !$acc loop seq
            do g = 1, REMAP_TGROUP
               nt = nt0 + g
               if (nt > ntrace) exit
               mtflx(i, j, nt, b) = acc(g)
            end do
         end if
      end do
      mflx(i, j, b) = fm
   end subroutine remap_face_flux

   subroutine remap_flux_east_impl(nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, puny, tracer_type, depend, src, xp, yp, &
                                   area, aim, mx, my, tc, tx, ty, mflx, mtflx)
      !! East-edge fluxes of the categories `nlo..nhi` over `E(is-1:ie, js:je)`:
      !! ONE `do concurrent` over `(b, j, i)`, the triangle and tracer loops inside.
      integer, intent(in) :: nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      integer(int8), intent(in) :: src(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), mx(nx, ny, nb), my(nx, ny, nb)
      real(wp), intent(in) :: tc(nx, ny, ntrace, nb), tx(nx, ny, ntrace, nb), ty(nx, ny, ntrace, nb)
      real(wp), intent(inout) :: mflx(nx, ny, nb), mtflx(nx, ny, ntrace, nb)
      integer :: i, j, b

      !$acc data present(tracer_type, depend, src, xp, yp, area, aim, mx, my, tc, tx, ty, mflx, mtflx)
      do concurrent(b=1:nhi - nlo + 1, j=js:je, i=is - 1:ie)
         call remap_face_flux(i, j, nlo + b - 1, b, nx, ny, ncat, ntrace, nb, puny, tracer_type, depend, src, xp, yp, area, &
                              aim, mx, my, tc, tx, ty, mflx, mtflx)
      end do
      !$acc end data
   end subroutine remap_flux_east_impl

   subroutine remap_flux_north_impl(nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, puny, tracer_type, depend, src, xp, yp, &
                                    area, aim, mx, my, tc, tx, ty, mflx, mtflx)
      !! North-edge fluxes over `N(is:ie, js-1:je)` (see `remap_flux_east_impl`).
      integer, intent(in) :: nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      integer(int8), intent(in) :: src(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: xp(nx, ny, 3, REMAP_NSLOT), yp(nx, ny, 3, REMAP_NSLOT), area(nx, ny, REMAP_NSLOT)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), mx(nx, ny, nb), my(nx, ny, nb)
      real(wp), intent(in) :: tc(nx, ny, ntrace, nb), tx(nx, ny, ntrace, nb), ty(nx, ny, ntrace, nb)
      real(wp), intent(inout) :: mflx(nx, ny, nb), mtflx(nx, ny, ntrace, nb)
      integer :: i, j, b

      !$acc data present(tracer_type, depend, src, xp, yp, area, aim, mx, my, tc, tx, ty, mflx, mtflx)
      do concurrent(b=1:nhi - nlo + 1, j=js - 1:je, i=is:ie)
         call remap_face_flux(i, j, nlo + b - 1, b, nx, ny, ncat, ntrace, nb, puny, tracer_type, depend, src, xp, yp, area, &
                              aim, mx, my, tc, tx, ty, mflx, mtflx)
      end do
      !$acc end data
   end subroutine remap_flux_north_impl

   ! ------------------------------------------------------------------
   ! update kernel
   ! ------------------------------------------------------------------

   subroutine remap_update_impl(nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, puny, tracer_type, depend, tarea, &
                                aim, trm, mflx_e, mflx_n, mtflx_e, mtflx_n, aimo, trmo, nbad, nneg)
      !! New area and tracers of the categories `nlo..nhi` on the owned cells
      !! (one `do concurrent` over `(b, j, i)`; the tracer loop is sequential
      !! inside, in hierarchy order: a dependent reads its parent's NEW value
      !! from the output it has just written).  See the module docstring.
      integer, intent(in) :: nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: tarea(nx, ny)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(in) :: mflx_e(nx, ny, nb), mflx_n(nx, ny, nb), mtflx_e(nx, ny, ntrace, nb), mtflx_n(nx, ny, ntrace, nb)
      real(wp), intent(inout) :: aimo(nx, ny, 0:ncat), trmo(nx, ny, ntrace, ncat)
      integer, intent(out) :: nbad, nneg
      integer :: i, j, b, n, nt, p
      real(wp) :: mm0, dm, tar, rt, mnew, tm0, wt, tnew, tp0, tpn

      nbad = 0
      nneg = 0
      !$acc data present(tracer_type, depend, tarea, aim, trm, mflx_e, mflx_n, mtflx_e, mtflx_n, aimo, trmo)
      do concurrent(b=1:nhi - nlo + 1, j=js:je, i=is:ie) local(n, nt, p, mm0, dm, tar, rt, mnew, tm0, wt, tnew, tp0, tpn) &
         reduce(+:nbad, nneg)
         n = nlo + b - 1
         mm0 = aim(i, j, n)
         dm = ((mflx_e(i, j, b) - mflx_e(i - 1, j, b)) + mflx_n(i, j, b)) - mflx_n(i, j - 1, b)
         tar = tarea(i, j)
         rt = 0.0_wp
         if (tar > 0.0_wp) rt = 1.0_wp/tar
         mnew = mm0 - dm*rt
         if (.not. ieee_is_finite(mnew)) nbad = nbad + 1
         if (mnew < -puny) then
            nneg = nneg + 1
         else if (mnew < 0.0_wp) then
            mnew = 0.0_wp
         end if
         aimo(i, j, n) = mnew
         if (n >= 1) then
            !$acc loop seq
            do nt = 1, ntrace
               tm0 = trm(i, j, nt, n)
               wt = ((mtflx_e(i, j, nt, b) - mtflx_e(i - 1, j, nt, b)) + mtflx_n(i, j, nt, b)) - mtflx_n(i, j - 1, nt, b)
               tnew = 0.0_wp
               if (mnew > 0.0_wp) then
                  if (tracer_type(nt) == 1) then
                     if (wt == 0.0_wp .and. dm == 0.0_wp) then
                        tnew = tm0
                     else
                        tnew = (mm0*tm0 - wt*rt)/mnew
                     end if
                  else
                     p = depend(nt)
                     tp0 = trm(i, j, p, n)
                     tpn = trmo(i, j, p, n)
                     if (wt == 0.0_wp .and. dm == 0.0_wp .and. tpn == tp0) then
                        tnew = tm0
                     else if (abs(tpn) > 0.0_wp) then
                        tnew = ((mm0*tp0)*tm0 - wt*rt)/(mnew*tpn)
                     end if
                  end if
               end if
               trmo(i, j, nt, n) = tnew
               if (.not. ieee_is_finite(tnew)) nbad = nbad + 1
            end do
         end if
      end do
      !$acc end data
   end subroutine remap_update_impl

   subroutine remap_ledger_impl(nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, puny, tracer_type, depend, &
                                nq1, nq2, nqs, tarea, aim, trm, mflx_e, mflx_n, mtflx_e, mtflx_n, &
                                nclip, clip_area, d_vi, d_vs, d_ei, d_es)
      !! What the guards of `remap_update_impl` create or drop on the owned cells of the categories
      !! `nlo..nhi`: the cells clipped from `[-puny, 0)` to 0 (count and area created), and the content
      !! of the tracer products left behind where the new area is 0 or a type-2 tracer's parent is 0
      !! (ice and snow volume, ice enthalpy `nq1..nq2`, snow enthalpy `nqs`: `trm` indices), times the
      !! cell area.  One `do concurrent` with six `reduce(+)` accumulators; it reads the fluxes the
      !! update read and recomputes the new area (the same expression).
      integer, intent(in) :: nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, nq1, nq2, nqs
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: tarea(nx, ny)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(in) :: mflx_e(nx, ny, nb), mflx_n(nx, ny, nb), mtflx_e(nx, ny, ntrace, nb), mtflx_n(nx, ny, ntrace, nb)
      integer, intent(out) :: nclip
      real(wp), intent(out) :: clip_area, d_vi, d_vs, d_ei, d_es
      integer :: i, j, b, n, nt, p
      real(wp) :: mm0, dm, tar, rt, mnew, wt, wtp, res, tpn, tp0

      nclip = 0
      clip_area = 0.0_wp
      d_vi = 0.0_wp
      d_vs = 0.0_wp
      d_ei = 0.0_wp
      d_es = 0.0_wp
      !$acc data present(tracer_type, depend, tarea, aim, trm, mflx_e, mflx_n, mtflx_e, mtflx_n)
      do concurrent(b=1:nhi - nlo + 1, j=js:je, i=is:ie) local(n, nt, p, mm0, dm, tar, rt, mnew, wt, wtp, res, tpn, tp0) &
         reduce(+:nclip, clip_area, d_vi, d_vs, d_ei, d_es)
         n = nlo + b - 1
         mm0 = aim(i, j, n)
         dm = ((mflx_e(i, j, b) - mflx_e(i - 1, j, b)) + mflx_n(i, j, b)) - mflx_n(i, j - 1, b)
         tar = tarea(i, j)
         rt = 0.0_wp
         if (tar > 0.0_wp) rt = 1.0_wp/tar
         mnew = mm0 - dm*rt
         if (mnew < 0.0_wp .and. mnew >= -puny) then
            nclip = nclip + 1
            clip_area = clip_area - mnew*tar
         end if
         if (n >= 1) then
            !$acc loop seq
            do nt = 1, ntrace
               if (.not. (nt == 1 .or. nt == 2 .or. (nt >= nq1 .and. nt <= nq2) .or. nt == nqs)) cycle
               wt = ((mtflx_e(i, j, nt, b) - mtflx_e(i - 1, j, nt, b)) + mtflx_n(i, j, nt, b)) - mtflx_n(i, j - 1, nt, b)
               res = 0.0_wp
               if (.not. (mnew > 0.0_wp)) then
                  if (tracer_type(nt) == 1) then
                     res = mm0*trm(i, j, nt, n) - wt*rt
                  else
                     p = depend(nt)
                     res = (mm0*trm(i, j, p, n))*trm(i, j, nt, n) - wt*rt
                  end if
               else if (tracer_type(nt) /= 1) then
                  p = depend(nt)
                  tp0 = trm(i, j, p, n)
                  wtp = ((mtflx_e(i, j, p, b) - mtflx_e(i - 1, j, p, b)) + mtflx_n(i, j, p, b)) - mtflx_n(i, j - 1, p, b)
                  tpn = (mm0*tp0 - wtp*rt)/mnew
                  if (tpn == 0.0_wp .and. tp0 /= 0.0_wp) then
                     res = (mm0*tp0)*trm(i, j, nt, n) - wt*rt
                  end if
               end if
               if (nt == 1) then
                  d_vi = d_vi + res*tar
               else if (nt == 2) then
                  d_vs = d_vs + res*tar
               else if (nt == nqs) then
                  d_es = d_es + res*tar
               else
                  d_ei = d_ei + res*tar
               end if
            end do
         end if
      end do
      !$acc end data
   end subroutine remap_ledger_impl

   ! ------------------------------------------------------------------
   ! geometry error word
   ! ------------------------------------------------------------------

   subroutine remap_geometry_errors_impl(nx, ny, is, ie, js, je, gerr, ice_err, ncfl, nnonfinite, nconserv)
      !! Fold the geometry's error word `gerr` into the cell error word
      !! `ice_err` over `is:ie, js:je`, count the cells carrying each fatal
      !! bit (`ICE_ERR_TRANSPORT_CFL`, `ICE_ERR_NONFINITE`,
      !! `ICE_ERR_TRANSPORT_CONSERV`) and reset `gerr` to zero (the invariant
      !! the geometry relies on).  ONE `do concurrent`, three `reduce(+)`
      !! accumulators.
      integer, intent(in) :: nx, ny, is, ie, js, je
      integer(int32), intent(inout) :: gerr(nx, ny), ice_err(nx, ny)
      integer, intent(out) :: ncfl, nnonfinite, nconserv
      integer :: i, j
      integer(int32) :: e

      ncfl = 0
      nnonfinite = 0
      nconserv = 0
      !$acc data present(gerr, ice_err)
      do concurrent(j=js:je, i=is:ie) local(e) reduce(+:ncfl, nnonfinite, nconserv)
         e = gerr(i, j)
         if (e /= 0_int32) then
            if (iand(e, ICE_ERR_TRANSPORT_CFL) /= 0_int32) ncfl = ncfl + 1
            if (iand(e, ICE_ERR_NONFINITE) /= 0_int32) nnonfinite = nnonfinite + 1
            if (iand(e, ICE_ERR_TRANSPORT_CONSERV) /= 0_int32) nconserv = nconserv + 1
            ice_err(i, j) = ior(ice_err(i, j), e)
            gerr(i, j) = 0_int32
         end if
      end do
      !$acc end data
   end subroutine remap_geometry_errors_impl

   ! ------------------------------------------------------------------
   ! host shims
   ! ------------------------------------------------------------------

   subroutine remap_flux_integrals(ws, gw, rb, fw, is, ie, js, je, nlo, nhi, puny)
      !! Fluxes of the categories `nlo..nhi` (`nhi - nlo + 1 <= fw%nbatch`) into
      !! `fw` (two launches).  `ws` supplies the tracer table; geometry and
      !! reconstruction must be current.  Everything device-present.
      type(transport_workspace_t), intent(in) :: ws
      type(remap_geom_workspace_t), intent(in) :: gw
      type(remap_recon_batch_t), intent(in) :: rb
      type(remap_flux_workspace_t), intent(inout) :: fw
      integer, intent(in) :: is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny

      call remap_flux_east_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, fw%nbatch, is, ie, js, je, nlo, nhi, puny, &
                                ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), gw%src_e, gw%xp_e, gw%yp_e, &
                                gw%area_e, ws%aim, rb%mx, rb%my, rb%tc, rb%tx, rb%ty, fw%mflx_e, fw%mtflx_e)
      call remap_flux_north_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, fw%nbatch, is, ie, js, je, nlo, nhi, puny, &
                                 ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), gw%src_n, gw%xp_n, gw%yp_n, &
                                 gw%area_n, ws%aim, rb%mx, rb%my, rb%tc, rb%tx, rb%ty, fw%mflx_n, fw%mtflx_n)
   end subroutine remap_flux_integrals

   subroutine remap_update(ws, wso, fw, is, ie, js, je, nlo, nhi, puny, tarea, nbad, nneg)
      !! New fields of the categories `nlo..nhi` from the fluxes in `fw` into
      !! `wso` over the owned cells (one launch).
      type(transport_workspace_t), intent(in) :: ws
      type(transport_workspace_t), intent(inout) :: wso
      type(remap_flux_workspace_t), intent(in) :: fw
      integer, intent(in) :: is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: tarea(ws%nx, ws%ny)
      integer, intent(out) :: nbad, nneg

      call remap_update_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, fw%nbatch, is, ie, js, je, nlo, nhi, puny, &
                             ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), tarea, ws%aim, ws%trm, &
                             fw%mflx_e, fw%mflx_n, fw%mtflx_e, fw%mtflx_n, wso%aim, wso%trm, nbad, nneg)
   end subroutine remap_update

   subroutine remap_advect(ws, wso, gw, rb, fw, is, ie, js, je, dt, puny, u, v, dxU, dyU, HTE, HTN, tarea, earea, narea, &
                           epm, npm, uvm, hm, ice_err, nbad, nneg, ngeo, nrecon, nlaunch, ledger, cfg, seam)
      !! Advect the packed fields of `ws` into `wso` over the owned cells
      !! `is:ie, js:je`: geometry (six launches plus the error fold), the input check (one), then
      !! per category batch the fluxes (two), the update (one) and the ledger (one).  `wso` is
      !! written on the owned cells only.  `cfg /= REMAP_CFG_OK`: nothing was
      !! run.  `ngeo(1:3)` counts the cells whose geometry raised
      !! `ICE_ERR_TRANSPORT_CFL`, `ICE_ERR_NONFINITE`, `ICE_ERR_TRANSPORT_CONSERV`
      !! on the windows the geometry writes, `is-1:ie, js-1:je` (the bits are also ORed into `ice_err`); when any is nonzero the run
      !! stops after the geometry (nothing else is formed, `wso` untouched).
      !! `nrecon` is the non-finite count of the reconstruction's read ring (the
      !! caller turns it into `ICE_ERR_NONFINITE`; `wso` untouched; the ring has no cell of its
      !! own, so the driver marks the FIRST OWNED cell, not the offender).  `ledger` is what the update's
      !! guards create or drop (`remap_ledger_t`).  `nbad`,
      !! `nneg` as the upwind kernel.  `nlaunch` counts the launches issued.
      !! Every array must be device-present on an offload build.
      type(transport_workspace_t), intent(in) :: ws
      type(transport_workspace_t), intent(inout) :: wso
      type(remap_geom_workspace_t), intent(inout) :: gw
      type(remap_recon_batch_t), intent(inout) :: rb
      type(remap_flux_workspace_t), intent(inout) :: fw
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(in) :: dt, puny
      real(wp), intent(in) :: u(ws%nx + 1, ws%ny), v(ws%nx, ws%ny + 1)
      real(wp), intent(in) :: dxU(ws%nx, ws%ny), dyU(ws%nx, ws%ny), HTE(ws%nx, ws%ny), HTN(ws%nx, ws%ny)
      real(wp), intent(in) :: tarea(ws%nx, ws%ny), earea(ws%nx, ws%ny), narea(ws%nx, ws%ny)
      real(wp), intent(in) :: epm(ws%nx, ws%ny), npm(ws%nx, ws%ny), uvm(ws%nx, ws%ny), hm(ws%nx, ws%ny)
      integer(int32), intent(inout) :: ice_err(ws%nx, ws%ny)
      integer, intent(out) :: nbad, nneg, nrecon, nlaunch, cfg
      integer, intent(out) :: ngeo(3)
      type(remap_ledger_t), intent(out) :: ledger
      type(transport_seam_t), intent(in), optional :: seam
         !! The decomposition (module `rdb_ice_transport_seam`; default one rank, no fold).  With it
         !! `ngeo`, `nrecon`, `nbad`, `nneg` and `cfg` are GLOBAL (the same on every rank, so a stop
         !! is a stop everywhere) and the fold-row fluxes of every batch are paired before the
         !! update; `ledger` stays the sum over this rank's owned cells.
      type(transport_seam_t) :: sm
      integer :: gv(NGV)
      integer :: lnc
      real(wp) :: lca, lvi, lvs, lei, les
      integer :: nlo, nhi, nbb, nnn, rcfg, i0, i1, j0, j1

      if (present(seam)) sm = seam
      nbad = 0
      nneg = 0
      ngeo = 0
      ledger = remap_ledger_t()
      nrecon = 0
      nlaunch = 0
      cfg = REMAP_CFG_OK
      if (.not. (ws%is_init .and. wso%is_init .and. gw%is_init .and. rb%is_init .and. fw%is_init)) then
         cfg = REMAP_CFG_WORKSPACE
         return
      end if
      if (ws%nx /= wso%nx .or. ws%ny /= wso%ny .or. ws%ncat /= wso%ncat .or. ws%ntrace /= wso%ntrace &
          .or. gw%nx /= ws%nx .or. gw%ny /= ws%ny .or. rb%nx /= ws%nx .or. rb%ny /= ws%ny &
          .or. rb%ncat /= ws%ncat .or. rb%ntrace /= ws%ntrace .or. rb%nbatch /= fw%nbatch .or. fw%nx /= ws%nx &
          .or. fw%ny /= ws%ny .or. fw%ncat /= ws%ncat .or. fw%ntrace /= ws%ntrace) then
         cfg = REMAP_CFG_WORKSPACE
         return
      end if
      if (.not. recon_table_ok(ws)) then
         cfg = REMAP_CFG_TABLE
         return
      end if
      if (is < 4 .or. js < 4 .or. ie > ws%nx - 3 .or. je > ws%ny - 3 .or. ie < is .or. je < js) then
         cfg = REMAP_CFG_RANGE
         return
      end if

      ! geometry once; the error word is folded into the cell word and reset
      call remap_geometry(gw, is, ie, js, je, dt, 1, puny, u, v, dxU, dyU, HTE, HTN, earea, narea, epm, npm, uvm, fw%gerr)
      nlaunch = nlaunch + 6
      ! the windows the geometry writes: corners U(is-1:ie, js-1:je), east edges E(is-1:ie, js:je), north N(is:ie, js-1:je)
      i0 = is - 1
      i1 = ie
      j0 = js - 1
      j1 = je
      call remap_geometry_errors_impl(ws%nx, ws%ny, i0, i1, j0, j1, fw%gerr, ice_err, ngeo(1), ngeo(2), ngeo(3))
      nlaunch = nlaunch + 1
      call seam_sum_int(sm, 3, ngeo)
      if (any(ngeo > 0)) return

      ! the read ring of the reconstruction: finite values, a 0/1 land mask, the table and the range
      call remap_input_check(ws, is, ie, js, je, hm, nrecon, rcfg)
      nlaunch = nlaunch + 1
      ! decided on the global counts: a land mask or a table refused by one tile refuses the step everywhere
      gv(1) = nrecon
      gv(2) = merge(1, 0, rcfg == RECON_CFG_MASK)
      gv(3) = merge(1, 0, rcfg /= RECON_CFG_OK .and. rcfg /= RECON_CFG_MASK)
      call seam_sum_int(sm, 3, gv(1:3))
      nrecon = gv(1)
      if (gv(3) > 0) then
         select case (rcfg)
         case (RECON_CFG_TABLE)
            cfg = REMAP_CFG_TABLE
         case (RECON_CFG_RANGE)
            cfg = REMAP_CFG_RANGE
         case default
            cfg = REMAP_CFG_WORKSPACE
         end select
         return
      end if
      if (gv(2) > 0) then
         cfg = REMAP_CFG_MASK
         return
      end if
      if (nrecon > 0) return

      ! fluxes and update, a category batch at a time
      nlo = 0
      do while (nlo <= ws%ncat)
         nhi = min(ws%ncat, nlo + fw%nbatch - 1)
         ! the batch's reconstruction on the owned cells plus one ring (one launch)
         call remap_recon_batch(ws, rb, is - 1, ie + 1, js - 1, je + 1, nlo, nhi, puny, hm, rcfg)
         if (rcfg /= RECON_CFG_OK) then
            cfg = REMAP_CFG_WORKSPACE
            return
         end if
         call remap_flux_integrals(ws, gw, rb, fw, is, ie, js, je, nlo, nhi, puny)
         ! the fold seam: both sides of a duplicated fold-row face take ONE flux (no-op off the fold)
         call seam_pair_fold_flux(sm, ws%nx, ws%ny, is, ie, je, nhi - nlo + 1, ws%ntrace, fw%mflx_n, fw%mtflx_n)
         call remap_update(ws, wso, fw, is, ie, js, je, nlo, nhi, puny, tarea, nbb, nnn)
         nbad = nbad + nbb
         nneg = nneg + nnn
         call remap_ledger_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, fw%nbatch, is, ie, js, je, nlo, nhi, puny, &
                                ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), 4, ws%tab%nt_qsno + 1, &
                                ws%tab%nt_qsno + 2, tarea, ws%aim, ws%trm, fw%mflx_e, fw%mflx_n, fw%mtflx_e, fw%mtflx_n, &
                                lnc, lca, lvi, lvs, lei, les)
         ledger%n_clip = ledger%n_clip + lnc
         ledger%clip_area = ledger%clip_area + lca
         ledger%drop_vice = ledger%drop_vice + lvi
         ledger%drop_vsno = ledger%drop_vsno + lvs
         ledger%drop_eice = ledger%drop_eice + lei
         ledger%drop_esno = ledger%drop_esno + les
         nlaunch = nlaunch + 5
         nlo = nhi + 1
      end do
      gv(1) = nbad
      gv(2) = nneg
      call seam_sum_int(sm, 2, gv(1:2))
      nbad = gv(1)
      nneg = gv(2)
   end subroutine remap_advect

end module rdb_ice_remap_flux
