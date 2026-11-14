!! Sea-ice transport checks (the bounds follow the monotonicity
!! argument of incremental remapping, Lipscomb and Hunke 2004).
module rdb_ice_transport_checks
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_driver (global_conservation, local_max_min, quasilocal_max_min, check_monotonicity).
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
   ! Baumgardner (2000, J. Comput. Phys. 160 318).
   !
   !! Six runtime checks around the horizontal transport, switched together by
   !! `&ocean_ice_dyn_nml transport_check`.  Each is a device kernel; the
   !! per-cell ones set a bit in the per-cell word `ice_err(nx, ny)` of
   !! `rdb_ice_errors` (a thread owns its cell, so no race and no atomics) and
   !! return the number of flagged cells through a ONE-accumulator
   !! `reduce(+)`:
   !!
   !! | # | check | kernel | bit |
   !! |---|---|---|---|
   !! | 1 | non-finite area / tracer / face velocity | `transport_check_finite`, `transport_check_cfl` | `ICE_ERR_NONFINITE` |
   !! | 2 | negative area below `-puny` (a finite round-off negative is tolerated) | `transport_check_neg_area` | `ICE_ERR_TRANSPORT_NEG_AREA` |
   !! | 2b | slivers `0 < aicen <= puny` entering the pack (its volumes would be dropped): precondition | `transport_check_slivers` | `ICE_ERR_TRANSPORT_CONSERV` |
   !! | 3 | CFL: a face displacement `|u| dt` beyond the narrower of the two cells it joins (checked in BOTH directions) | `transport_check_cfl` | `ICE_ERR_TRANSPORT_CFL` |
   !! | 4 | conservation of area, per category incl. open water | `transport_conserved_sums` + `transport_conservation_compare` | `ICE_ERR_TRANSPORT_CONSERV` (global) |
   !! | 5 | conservation of area-weighted tracer content per category | same | same |
   !! | 6 | monotonicity: new tracer inside the quasi-local (5x5) bound | `transport_bounds`, `transport_check_monotone` | `ICE_ERR_TRANSPORT_CONSERV` (per cell) |
   !!
   !! The registry has no separate monotonicity bit; it shares the conservation
   !! bit (both mean "an invariant of the remapping failed").  Upwind is
   !! monotone only inside its CFL limit, so check 6 is a remap-PR tool and an
   !! upwind regression guard.
   !!
   !! **Conservation sums are EFP reproducing sums** (`rdb_efp`): each summand
   !! `area*aim` or `area*aim*tracer[*parent]` is decomposed into six `int64`
   !! fixed-point bins on the device, so the totals are bitwise independent of
   !! summation order, backend and (with a bin-wise `MPI_SUM`) rank count.  A
   !! non-finite summand poisons the total (`efp_to_real` returns NaN).  The
   !! kernel shape is thread-per-`(quantity, i)` with a sequential `j` loop and
   !! seven `int64` locals, writing one row of bins per thread: no reduction
   !! clause is needed (integer addition is exact and order-free), so there is
   !! no accumulator count to keep small, `i` is last (coalesced), and one
   !! launch sums all `(ncat + 1) + ncat*ntrace` quantities; a second tiny
   !! launch sums the rows.  Conservation compares the totals before and after
   !! with the reference tolerance (relative `puny`, only where the initial
   !! total exceeds `puny` in magnitude).
   !!
   !! **Monotonicity bounds.**  `transport_bounds` stores, per `(i, j, nt, n)`,
   !! the min/max over the 5x5 neighbourhood of the 3x3 local extrema, where a
   !! neighbour that holds no ice (area test for tracers of area; parent test
   !! for dependent tracers) is replaced by the home cell, so non-physical
   !! values never widen the bound.  It reads `i +- 2`, `j +- 2`: the caller's
   !! range must keep two valid rings (`nghost >= 2` beyond the owned range or
   !! the range shrunk by two).  The scratch is `2*ntrace*ncat` arrays of the
   !! grid size: a DEBUG cost, allocated only when `transport_check` is on.
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, EFP_DIGITS, EFP_PREC_WIDTH, EFP_MAX_SUMMANDS, efp_carry, efp_to_real
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL, &
                             ICE_ERR_TRANSPORT_NEG_AREA, ICE_ERR_TRANSPORT_CONSERV
   implicit none
   private

   public :: transport_check_workspace_t
   public :: transport_nquantities
   public :: transport_check_finite, transport_check_neg_area, transport_check_cfl
   public :: transport_check_slivers, transport_check_slivers_impl
   public :: transport_conserved_sums, transport_conservation_compare
   public :: transport_bounds, transport_check_monotone
   public :: transport_check_finite_impl, transport_check_neg_area_impl, transport_check_cfl_impl
   public :: transport_efp_rows_impl, transport_efp_total_impl
   public :: transport_bounds_impl, transport_check_monotone_impl
   public :: transport_efp_decompose_seq

   integer, parameter :: NBIN = EFP_DIGITS + 1
      !! Six fixed-point bins plus the poison counter.

   type :: transport_check_workspace_t
      !! Persistent scratch of the checks.  `rows`/`tot` serve the EFP sums
      !! (always allocated by `init`); `tmin`/`tmax` serve the monotonicity
      !! bounds (allocated only with `want_bounds`).  `enter_data` maps the
      !! component arrays ONCE.
      logical :: is_init = .false.
      logical :: mapped = .false.
      logical :: has_bounds = .false.
      integer :: nq = 0
      integer(int64), allocatable :: rows(:, :, :)
         !! `(NBIN, nq, nx)`: per-column bin sums.
      integer(int64), allocatable :: tot(:, :)
         !! `(NBIN, nq)`: bin totals (host copy refreshed by each sum).
      real(wp), allocatable :: tmin(:, :, :, :)
         !! `(nx, ny, ntrace, ncat)` lower bound.
      real(wp), allocatable :: tmax(:, :, :, :)
         !! `(nx, ny, ntrace, ncat)` upper bound.
   contains
      procedure :: init => check_ws_init
      procedure :: enter_data => check_ws_enter_data
      procedure :: exit_data => check_ws_exit_data
      procedure :: destroy => check_ws_destroy
      procedure :: bytes => check_ws_bytes
   end type transport_check_workspace_t

contains

   pure function ws_shape_ok(ws, nx, ny, ncat, ntrace, want_bounds) result(ok)
      !! The workspace is initialised and its arrays match the extents the
      !! explicit-shape kernels are about to be given (bounds: also `tmin`,
      !! `tmax`).
      type(transport_check_workspace_t), intent(in) :: ws
      integer, intent(in) :: nx, ny, ncat, ntrace
      logical, intent(in) :: want_bounds
      logical :: ok

      ok = ws%is_init
      if (.not. ok) return
      if (want_bounds) then
         ok = ws%has_bounds .and. allocated(ws%tmin) .and. allocated(ws%tmax)
         if (.not. ok) return
         ok = all(shape(ws%tmin) == [nx, ny, ntrace, ncat]) .and. all(shape(ws%tmax) == [nx, ny, ntrace, ncat])
      else
         ok = allocated(ws%rows) .and. allocated(ws%tot)
         if (.not. ok) return
         ok = ws%nq == transport_nquantities(ncat, ntrace) .and. all(shape(ws%rows) == [NBIN, ws%nq, nx]) &
              .and. all(shape(ws%tot) == [NBIN, ws%nq])
      end if
   end function ws_shape_ok

   pure function transport_nquantities(ncat, ntrace) result(nq)
      !! Number of conserved quantities: area of categories `0..ncat` plus
      !! `ncat*ntrace` area-weighted tracer contents.  Order: `q = n + 1` for
      !! the area of category `n`; `q = ncat + 1 + (n - 1)*ntrace + nt` for
      !! tracer `nt` of category `n`.
      integer, intent(in) :: ncat, ntrace
      integer :: nq
      nq = (ncat + 1) + ncat*ntrace
   end function transport_nquantities

   ! ------------------------------------------------------------------
   ! workspace
   ! ------------------------------------------------------------------

   subroutine check_ws_init(this, nx, ny, ncat, ntrace, want_bounds)
      class(transport_check_workspace_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, ntrace
      logical, intent(in) :: want_bounds

      call this%destroy()
      this%nq = transport_nquantities(ncat, ntrace)
      allocate (this%rows(NBIN, this%nq, nx), this%tot(NBIN, this%nq))
      this%rows = 0_int64
      this%tot = 0_int64
      this%has_bounds = want_bounds
      if (want_bounds) then
         allocate (this%tmin(nx, ny, ntrace, ncat), this%tmax(nx, ny, ntrace, ncat))
         this%tmin = 0.0_wp
         this%tmax = 0.0_wp
      end if
      this%is_init = .true.
   end subroutine check_ws_init

   subroutine check_ws_enter_data(this)
      class(transport_check_workspace_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%rows, this%tot)
      if (this%has_bounds) then
         !$acc enter data copyin(this%tmin, this%tmax)
      end if
      this%mapped = .true.
   end subroutine check_ws_enter_data

   subroutine check_ws_exit_data(this)
      class(transport_check_workspace_t), intent(inout) :: this
      if (.not. this%mapped) return
      if (this%has_bounds) then
         !$acc exit data delete(this%tmin, this%tmax)
      end if
      !$acc exit data delete(this%rows, this%tot)
      this%mapped = .false.
   end subroutine check_ws_exit_data

   pure function check_ws_bytes(this) result(b)
      !! Bytes of the scratch (host copy; the device mirror is the same size).
      class(transport_check_workspace_t), intent(in) :: this
      integer(int64) :: b
      b = 0_int64
      if (allocated(this%rows)) b = b + int(size(this%rows), int64)*8_int64
      if (allocated(this%tot)) b = b + int(size(this%tot), int64)*8_int64
      if (allocated(this%tmin)) b = b + int(size(this%tmin), int64)*int(storage_size(this%tmin)/8, int64)
      if (allocated(this%tmax)) b = b + int(size(this%tmax), int64)*int(storage_size(this%tmax)/8, int64)
   end function check_ws_bytes

   subroutine check_ws_destroy(this)
      class(transport_check_workspace_t), intent(inout) :: this
      call this%exit_data()
      if (allocated(this%rows)) deallocate (this%rows)
      if (allocated(this%tot)) deallocate (this%tot)
      if (allocated(this%tmin)) deallocate (this%tmin)
      if (allocated(this%tmax)) deallocate (this%tmax)
      this%is_init = .false.
      this%has_bounds = .false.
      this%nq = 0
   end subroutine check_ws_destroy

   ! ------------------------------------------------------------------
   ! 1: non-finite entry
   ! ------------------------------------------------------------------

   subroutine transport_check_finite(nx, ny, ncat, ntrace, is, ie, js, je, aim, trm, ice_err, nviol)
      !! Check 1a: every packed `aim`, `trm` over the range is finite; a cell
      !! with any non-finite value gets `ICE_ERR_NONFINITE`.  `nviol` is the
      !! number of such cells.
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      call transport_check_finite_impl(nx, ny, ncat, ntrace, is, ie, js, je, aim, trm, ice_err, nviol)
   end subroutine transport_check_finite

   pure subroutine transport_check_finite_impl(nx, ny, ncat, ntrace, is, ie, js, je, aim, trm, ice_err, nviol)
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      integer :: i, j, n, nt
      logical :: bad

      nviol = 0
      !$acc data present(aim, trm, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, nt, bad) reduce(+:nviol)
         bad = .not. ieee_is_finite(aim(i, j, 0))
         !$acc loop seq
         do n = 1, ncat
            if (.not. ieee_is_finite(aim(i, j, n))) bad = .true.
            !$acc loop seq
            do nt = 1, ntrace
               if (.not. ieee_is_finite(trm(i, j, nt, n))) bad = .true.
            end do
         end do
         if (bad) then
            ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
            nviol = nviol + 1
         end if
      end do
      !$acc end data
   end subroutine transport_check_finite_impl

   ! ------------------------------------------------------------------
   ! 2: negative area
   ! ------------------------------------------------------------------

   subroutine transport_check_neg_area(nx, ny, ncat, is, ie, js, je, puny, aim, ice_err, nviol)
      !! Check 2: `aim < -puny` in any category (open water included) sets
      !! `ICE_ERR_TRANSPORT_NEG_AREA`.  A finite round-off negative inside
      !! `[-puny, 0)` is tolerated.  Non-finite values are NOT seen here
      !! (`x < -puny` is false for NaN): they are check 1's.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: aim(nx, ny, 0:ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      call transport_check_neg_area_impl(nx, ny, ncat, is, ie, js, je, puny, aim, ice_err, nviol)
   end subroutine transport_check_neg_area

   pure subroutine transport_check_neg_area_impl(nx, ny, ncat, is, ie, js, je, puny, aim, ice_err, nviol)
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: aim(nx, ny, 0:ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      integer :: i, j, n
      logical :: bad

      nviol = 0
      !$acc data present(aim, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, bad) reduce(+:nviol)
         bad = .false.
         !$acc loop seq
         do n = 0, ncat
            if (aim(i, j, n) < -puny) bad = .true.
         end do
         if (bad) then
            ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_TRANSPORT_NEG_AREA)
            nviol = nviol + 1
         end if
      end do
      !$acc end data
   end subroutine transport_check_neg_area_impl

   ! ------------------------------------------------------------------
   ! 2b: slivers (precondition of the pack)
   ! ------------------------------------------------------------------

   subroutine transport_check_slivers(nx, ny, ncat, is, ie, js, je, puny, aicen, vicen, vsnon, ice_err, &
                                      nsliver, vi_drop, vs_drop)
      !! Precondition check of the PACK, run on the state BEFORE `transport_pack`:
      !! a category with `0 < aicen <= puny` (a sliver) loses its ice and snow
      !! volume and its tracers in the pack (the reference behaviour; it is
      !! safe there only because the small-area cleanup has already zapped
      !! slivers with the fluxes booked).  Each cell holding a sliver gets
      !! `ICE_ERR_TRANSPORT_CONSERV`; `nsliver` counts sliver CATEGORIES and
      !! `vi_drop`, `vs_drop` are the ice and snow volume per cell area that
      !! would be dropped, summed over the range (unweighted by cell area).
      !! Without this check a conservation "before" sum taken after the pack
      !! cannot see the loss.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nsliver
      real(wp), intent(out) :: vi_drop, vs_drop
      call transport_check_slivers_impl(nx, ny, ncat, is, ie, js, je, puny, aicen, vicen, vsnon, ice_err, &
                                        nsliver, vi_drop, vs_drop)
   end subroutine transport_check_slivers

   pure subroutine transport_check_slivers_impl(nx, ny, ncat, is, ie, js, je, puny, aicen, vicen, vsnon, ice_err, &
                                                nsliver, vi_drop, vs_drop)
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nsliver
      real(wp), intent(out) :: vi_drop, vs_drop
      integer :: i, j, n
      integer :: k
      real(wp) :: vi, vs

      nsliver = 0
      vi_drop = 0.0_wp
      vs_drop = 0.0_wp
      !$acc data present(aicen, vicen, vsnon, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, k, vi, vs) reduce(+:nsliver, vi_drop, vs_drop)
         k = 0
         vi = 0.0_wp
         vs = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            if (aicen(i, j, n) > 0.0_wp .and. aicen(i, j, n) <= puny) then
               k = k + 1
               vi = vi + vicen(i, j, n)
               vs = vs + vsnon(i, j, n)
            end if
         end do
         if (k > 0) then
            ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_TRANSPORT_CONSERV)
            nsliver = nsliver + k
            vi_drop = vi_drop + vi
            vs_drop = vs_drop + vs
         end if
      end do
      !$acc end data
   end subroutine transport_check_slivers_impl

   ! ------------------------------------------------------------------
   ! 3: CFL
   ! ------------------------------------------------------------------

   subroutine transport_check_cfl(nx, ny, is, ie, js, je, dt, uE, vN, dxT, dyT, ice_err, nviol)
      !! Check 3 (with 1b): face velocities `uE(i,j)` on the EAST face of cell
      !! `(i,j)` and `vN(i,j)` on its NORTH face.  The displacement `|u| dt`
      !! must not exceed the width of either cell the face joins
      !! (`min(dxT(i,j), dxT(i+1,j))`, and `dyT` for `v`): the remapping
      !! departure region then stays inside the four neighbouring cells.  The
      !! test is SYMMETRIC in the sign of the velocity (the reference tests
      !! only one direction).  A non-finite velocity sets `ICE_ERR_NONFINITE`,
      !! an excess `ICE_ERR_TRANSPORT_CFL`.  The kernel checks exactly the faces
      !! `(i, j)`, `i = is:ie`, `j = js:je`: the east/north face of each cell in
      !! the range, so to cover every face bounding the owned cells
      !! `ilo:ihi, jlo:jhi` the caller passes `is = ilo - 1`, `ie = ihi`,
      !! `js = jlo - 1`, `je = jhi` (one west/south face more).  It also reads
      !! `dxT(ie+1, :)` and `dyT(:, je+1)`.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: uE(nx, ny), vN(nx, ny), dxT(nx, ny), dyT(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      call transport_check_cfl_impl(nx, ny, is, ie, js, je, dt, uE, vN, dxT, dyT, ice_err, nviol)
   end subroutine transport_check_cfl

   pure subroutine transport_check_cfl_impl(nx, ny, is, ie, js, je, dt, uE, vN, dxT, dyT, ice_err, nviol)
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: uE(nx, ny), vN(nx, ny), dxT(nx, ny), dyT(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      integer :: i, j
      integer(int32) :: w
      real(wp) :: dx, dy, du, dv

      nviol = 0
      !$acc data present(uE, vN, dxT, dyT, ice_err)
      do concurrent(j=js:je, i=is:ie) local(w, dx, dy, du, dv) reduce(+:nviol)
         w = ICE_ERR_NONE
         dx = min(dxT(i, j), dxT(i + 1, j))
         dy = min(dyT(i, j), dyT(i, j + 1))
         if (.not. (ieee_is_finite(uE(i, j)) .and. ieee_is_finite(vN(i, j)))) then
            w = ior(w, ICE_ERR_NONFINITE)
         else
            du = abs(uE(i, j))*dt
            dv = abs(vN(i, j))*dt
            if (du > dx .or. dv > dy) w = ior(w, ICE_ERR_TRANSPORT_CFL)
         end if
         if (w /= ICE_ERR_NONE) then
            ice_err(i, j) = ior(ice_err(i, j), w)
            nviol = nviol + 1
         end if
      end do
      !$acc end data
   end subroutine transport_check_cfl_impl

   ! ------------------------------------------------------------------
   ! 4, 5: conservation (EFP sums)
   ! ------------------------------------------------------------------

   pure subroutine transport_efp_decompose_seq(r, e1, e2, e3, e4, e5, e6, epoison)
      !! Device twin of `rdb_efp::efp_decompose` (same arithmetic, scalar
      !! outputs so the callers accumulate in registers; `epoison = 1` for a
      !! non-finite or bin-1-overflowing input, bins zero then).  Pinned
      !! bin-for-bin against the canonical procedure by the tests.
      !$acc routine seq
      real(real64), intent(in) :: r
      integer(int64), intent(out) :: e1, e2, e3, e4, e5, e6, epoison
      real(real64) :: rs, s
      real(real64), parameter :: PR1 = 2.0_real64**(2*EFP_PREC_WIDTH)
      real(real64), parameter :: PR2 = 2.0_real64**(1*EFP_PREC_WIDTH)
      real(real64), parameter :: PR3 = 1.0_real64
      real(real64), parameter :: PR4 = 2.0_real64**(-1*EFP_PREC_WIDTH)
      real(real64), parameter :: PR5 = 2.0_real64**(-2*EFP_PREC_WIDTH)
      real(real64), parameter :: PR6 = 2.0_real64**(-3*EFP_PREC_WIDTH)
      real(real64), parameter :: IPR1 = 1.0_real64/PR1
      real(real64), parameter :: IPR2 = 1.0_real64/PR2
      real(real64), parameter :: IPR3 = 1.0_real64/PR3
      real(real64), parameter :: IPR4 = 1.0_real64/PR4
      real(real64), parameter :: IPR5 = 1.0_real64/PR5
      real(real64), parameter :: IPR6 = 1.0_real64/PR6
      real(real64), parameter :: MAX_E1 = real(huge(0_int64), real64)

      e1 = 0_int64
      e2 = 0_int64
      e3 = 0_int64
      e4 = 0_int64
      e5 = 0_int64
      e6 = 0_int64
      epoison = 0_int64
      if (.not. ieee_is_finite(r)) then
         epoison = 1_int64
         return
      end if
      s = 1.0_real64
      rs = r
      if (rs < 0.0_real64) then
         s = -1.0_real64
         rs = -rs
      end if
      if (rs*IPR1 >= MAX_E1) then
         epoison = 1_int64
         return
      end if
      e1 = int(s*aint(rs*IPR1), int64)
      rs = rs - real(abs(e1), real64)*PR1
      e2 = int(s*aint(rs*IPR2), int64)
      rs = rs - real(abs(e2), real64)*PR2
      e3 = int(s*aint(rs*IPR3), int64)
      rs = rs - real(abs(e3), real64)*PR3
      e4 = int(s*aint(rs*IPR4), int64)
      rs = rs - real(abs(e4), real64)*PR4
      e5 = int(s*aint(rs*IPR5), int64)
      rs = rs - real(abs(e5), real64)*PR5
      e6 = int(s*aint(rs*IPR6), int64)
   end subroutine transport_efp_decompose_seq

   subroutine transport_conserved_sums(ws, nx, ny, ncat, ntrace, is, ie, js, je, depend, area, aim, trm, &
                                       sums, ok)
      !! Global (this tile) EFP sums of all `transport_nquantities` conserved
      !! quantities over the OWNED range `is:ie, js:je` (never ghosts):
      !! `area*aim(n)` and `area*aim(n)*trm(nt)[*trm(depend(nt))]` per
      !! category.  `depend` is `tab%depend`.  `ok = .false.` (and `sums`
      !! zero) if the range holds more summands than `EFP_MAX_SUMMANDS`.  A
      !! cross-rank total is a bin-wise `MPI_SUM` of `sums(q)%v` and `%poison`
      !! (`efp_to_transport`), added by the multi-rank PR.
      type(transport_check_workspace_t), intent(inout) :: ws
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      integer, intent(in) :: depend(ntrace)
      real(wp), intent(in) :: area(nx, ny)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      type(efp_t), intent(out) :: sums(:)
      logical, intent(out) :: ok
      integer :: q, nq

      nq = ws%nq
      ok = .false.
      do q = 1, size(sums)
         sums(q)%v = 0_int64
         sums(q)%poison = 0_int64
      end do
      if (.not. ws_shape_ok(ws, nx, ny, ncat, ntrace, .false.)) return
      if (size(sums) < nq) return
      if (int(ie - is + 1, int64)*int(je - js + 1, int64) > EFP_MAX_SUMMANDS) return
      call transport_efp_rows_impl(nx, ny, ncat, ntrace, nq, is, ie, js, je, depend, area, aim, trm, ws%rows)
      call transport_efp_total_impl(nq, is, ie, nx, ws%rows, ws%tot)
      !$acc update self(ws%tot)
      do q = 1, nq
         sums(q)%v = ws%tot(1:EFP_DIGITS, q)
         sums(q)%poison = ws%tot(NBIN, q)
         call efp_carry(sums(q)%v)
      end do
      ok = .true.
   end subroutine transport_conserved_sums

   pure subroutine transport_efp_rows_impl(nx, ny, ncat, ntrace, nq, is, ie, js, je, depend, area, aim, trm, rows)
      !! Launch 1 of the sums: thread `(q, i)` sums its column `j = js:je`
      !! into seven `int64` locals and stores them as `rows(:, q, i)`.
      integer, intent(in) :: nx, ny, ncat, ntrace, nq, is, ie, js, je
      integer, intent(in) :: depend(ntrace)
      real(wp), intent(in) :: area(nx, ny), aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int64), intent(inout) :: rows(NBIN, nq, nx)
      integer :: q, i, j, n, nt, qq, par
      integer(int64) :: s1, s2, s3, s4, s5, s6, sp, d1, d2, d3, d4, d5, d6, dp
      real(real64) :: v

      !$acc data present(depend, area, aim, trm, rows)
      do concurrent(q=1:nq, i=is:ie) local(j, n, nt, qq, par, s1, s2, s3, s4, s5, s6, sp, &
                                           d1, d2, d3, d4, d5, d6, dp, v)
         s1 = 0_int64
         s2 = 0_int64
         s3 = 0_int64
         s4 = 0_int64
         s5 = 0_int64
         s6 = 0_int64
         sp = 0_int64
         if (q <= ncat + 1) then
            n = q - 1
            nt = 0
         else
            qq = q - ncat - 1
            n = (qq - 1)/ntrace + 1
            nt = qq - (n - 1)*ntrace
         end if
         !$acc loop seq
         do j = js, je
            if (nt == 0) then
               v = area(i, j)*aim(i, j, n)
            else
               par = depend(nt)
               if (par == 0) then
                  v = area(i, j)*aim(i, j, n)*trm(i, j, nt, n)
               else
                  v = area(i, j)*aim(i, j, n)*(trm(i, j, nt, n)*trm(i, j, par, n))
               end if
            end if
            call transport_efp_decompose_seq(v, d1, d2, d3, d4, d5, d6, dp)
            s1 = s1 + d1
            s2 = s2 + d2
            s3 = s3 + d3
            s4 = s4 + d4
            s5 = s5 + d5
            s6 = s6 + d6
            sp = sp + dp
         end do
         rows(1, q, i) = s1
         rows(2, q, i) = s2
         rows(3, q, i) = s3
         rows(4, q, i) = s4
         rows(5, q, i) = s5
         rows(6, q, i) = s6
         rows(7, q, i) = sp
      end do
      !$acc end data
   end subroutine transport_efp_rows_impl

   pure subroutine transport_efp_total_impl(nq, is, ie, nx, rows, tot)
      !! Launch 2: `tot(b, q) = sum_i rows(b, q, i)` (exact integer sums).
      integer, intent(in) :: nq, is, ie, nx
      integer(int64), intent(in) :: rows(NBIN, nq, nx)
      integer(int64), intent(inout) :: tot(NBIN, nq)
      integer :: b, q, i
      integer(int64) :: s

      !$acc data present(rows, tot)
      do concurrent(q=1:nq, b=1:NBIN) local(i, s)
         s = 0_int64
         !$acc loop seq
         do i = is, ie
            s = s + rows(b, q, i)
         end do
         tot(b, q) = s
      end do
      !$acc end data
   end subroutine transport_efp_total_impl

   subroutine transport_conservation_compare(ncat, ntrace, before, after, puny, ierr, worst_rel, worst_q)
      !! Checks 4 and 5 (host, a few scalars): per quantity, where
      !! `|before| > puny` the relative change `|after - before|/|before|`
      !! must not exceed `puny`.  A poisoned (non-finite) total is a
      !! violation.  `ierr = ICE_ERR_TRANSPORT_CONSERV` on any violation;
      !! `worst_rel` is the largest relative change among the compared
      !! quantities (NaN if any total is poisoned) and `worst_q` its index
      !! (the first poisoned one).  The first `ncat + 1` quantities are the
      !! areas (check 4), the rest the tracer contents (check 5).
      integer, intent(in) :: ncat, ntrace
      type(efp_t), intent(in) :: before(:), after(:)
      real(wp), intent(in) :: puny
      integer(int32), intent(out) :: ierr
      real(wp), intent(out) :: worst_rel
      integer, intent(out) :: worst_q
      integer :: q, nq
      real(real64) :: b, a, rel

      nq = transport_nquantities(ncat, ntrace)
      ierr = ICE_ERR_NONE
      worst_rel = 0.0_wp
      worst_q = 0
      if (size(before) < nq .or. size(after) < nq) then
         ierr = ICE_ERR_TRANSPORT_CONSERV
         worst_rel = ieee_value(worst_rel, ieee_quiet_nan)
         return
      end if
      do q = 1, nq
         b = efp_to_real(before(q))
         a = efp_to_real(after(q))
         if (.not. (ieee_is_finite(b) .and. ieee_is_finite(a))) then
            ierr = ICE_ERR_TRANSPORT_CONSERV
            worst_rel = ieee_value(worst_rel, ieee_quiet_nan)
            worst_q = q
            return
         end if
         if (abs(b) > puny) then
            rel = abs((a - b)/b)
            if (rel > worst_rel) then
               worst_rel = rel
               worst_q = q
            end if
            if (rel > puny) ierr = ICE_ERR_TRANSPORT_CONSERV
         end if
      end do
   end subroutine transport_conservation_compare

   ! ------------------------------------------------------------------
   ! 6: monotonicity
   ! ------------------------------------------------------------------

   subroutine transport_bounds(ws, nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, aim, trm, ok)
      !! Fill `ws%tmin`, `ws%tmax` (quasi-local bounds of the PRE-transport
      !! `aim`, `trm`) over `is:ie, js:je`; reads two cells beyond.  `ok =
      !! .false.` (nothing run) if the workspace is not initialised with
      !! bounds or its extents do not match `nx, ny, ncat, ntrace`.
      !!
      !! The bound is a clean 5x5 (3x3 of the 3x3 local extrema, each with its
      !! own masks), computed from the pre-transport field.  The reference
      !! extends the local bound in place with a row sweep, which can reach
      !! beyond 5x5 in one direction, so OUR bound is never wider and is
      !! STRICTLY TIGHTER: the check can flag a value the reference passes.
      type(transport_check_workspace_t), intent(inout) :: ws
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      logical, intent(out) :: ok
      ok = ws_shape_ok(ws, nx, ny, ncat, ntrace, .true.)
      if (.not. ok) return
      call transport_bounds_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, aim, trm, &
                                 ws%tmin, ws%tmax)
   end subroutine transport_bounds

   pure subroutine transport_bounds_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, &
                                         aim, trm, tmin, tmax)
      !! ONE launch over `(n, nt, j, i)`.  Per thread: the nine points of the
      !! 3x3 block, each contributing the extrema of its own masked 3x3
      !! (home-cell substitution for a masked-out neighbour), accumulated in
      !! two scalars.  The reference computes the same bound in two passes
      !! over stored 3x3 extrema; recomputing them keeps the scratch at the
      !! two bound arrays.
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(inout) :: tmin(nx, ny, ntrace, ncat), tmax(nx, ny, ntrace, ncat)
      integer :: i, j, n, nt, di, dj, ei, ej, ii, jj, par
      real(wp) :: lo, hi, home, v
      logical :: m

      !$acc data present(tracer_type, depend, aim, trm, tmin, tmax)
      do concurrent(n=1:ncat, nt=1:ntrace, j=js:je, i=is:ie) local(di, dj, ei, ej, ii, jj, par, lo, hi, home, v, m)
         par = depend(nt)
         lo = trm(i, j, nt, n)
         hi = lo
         !$acc loop seq
         do dj = -1, 1
            !$acc loop seq
            do di = -1, 1
               ! local 3x3 of the point (i+di, j+dj), home value there
               if (i + di < 1 .or. i + di > nx .or. j + dj < 1 .or. j + dj > ny) cycle
               home = trm(i + di, j + dj, nt, n)
               lo = min(lo, home)
               hi = max(hi, home)
               !$acc loop seq
               do ej = -1, 1
                  !$acc loop seq
                  do ei = -1, 1
                     ii = i + di + ei
                     jj = j + dj + ej
                     if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) cycle
                     if (tracer_type(nt) == 1) then
                        m = aim(ii, jj, n) > puny
                     else
                        m = abs(trm(ii, jj, par, n)) > puny
                     end if
                     if (m) then
                        v = trm(ii, jj, nt, n)
                     else
                        v = home
                     end if
                     lo = min(lo, v)
                     hi = max(hi, v)
                  end do
               end do
            end do
         end do
         tmin(i, j, nt, n) = lo
         tmax(i, j, nt, n) = hi
      end do
      !$acc end data
   end subroutine transport_bounds_impl

   subroutine transport_check_monotone(ws, nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, &
                                       aim, trm, ice_err, nviol, ok)
      !! Check 6: after transport, every tracer where it is defined (area
      !! above `puny`, or the parent magnitude above `puny`) lies inside
      !! `[tmin - w*puny, tmax + w*puny]`, `w = max(1, |bound|)` (round-off
      !! allowance).  A violating or non-finite value sets
      !! `ICE_ERR_TRANSPORT_CONSERV` on the cell; `nviol` counts cells.  `ok =
      !! .false.` (nothing run) on a workspace/extent mismatch.
      type(transport_check_workspace_t), intent(in) :: ws
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      logical, intent(out) :: ok
      nviol = 0
      ok = ws_shape_ok(ws, nx, ny, ncat, ntrace, .true.)
      if (.not. ok) return
      call transport_check_monotone_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, &
                                         aim, trm, ws%tmin, ws%tmax, ice_err, nviol)
   end subroutine transport_check_monotone

   pure subroutine transport_check_monotone_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, &
                                                 depend, aim, trm, tmin, tmax, ice_err, nviol)
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(in) :: tmin(nx, ny, ntrace, ncat), tmax(nx, ny, ntrace, ncat)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer, intent(out) :: nviol
      integer :: i, j, n, nt, par
      logical :: bad, chk
      real(wp) :: w1, w2, x

      nviol = 0
      !$acc data present(tracer_type, depend, aim, trm, tmin, tmax, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, nt, par, bad, chk, w1, w2, x) reduce(+:nviol)
         bad = .false.
         !$acc loop seq
         do n = 1, ncat
            !$acc loop seq
            do nt = 1, ntrace
               par = depend(nt)
               if (tracer_type(nt) == 1) then
                  chk = aim(i, j, n) > puny
               else
                  chk = abs(trm(i, j, par, n)) > puny
               end if
               if (chk) then
                  x = trm(i, j, nt, n)
                  if (.not. ieee_is_finite(x)) then
                     bad = .true.
                  else
                     w1 = max(1.0_wp, abs(tmin(i, j, nt, n)))
                     w2 = max(1.0_wp, abs(tmax(i, j, nt, n)))
                     if (x < tmin(i, j, nt, n) - w1*puny) bad = .true.
                     if (x > tmax(i, j, nt, n) + w2*puny) bad = .true.
                  end if
               end if
            end do
         end do
         if (bad) then
            ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_TRANSPORT_CONSERV)
            nviol = nviol + 1
         end if
      end do
      !$acc end data
   end subroutine transport_check_monotone_impl

end module rdb_ice_transport_checks
