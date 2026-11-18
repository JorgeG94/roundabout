!! Sea-ice incremental remapping, step 1: the masks and the limited linear
!! reconstruction of the area and of the tracer hierarchy inside each cell
!! (Dukowicz and Baumgardner 2000, J. Comput. Phys. 160, 318; Lipscomb and
!! Hunke 2004, Mon. Wea. Rev. 132, 1341).
module rdb_ice_remap_recon
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_remap (make_masks, construct_fields, limited_gradient).
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
   !! **What is reconstructed.**  Per category `n = 0..ncat` (`0` is open
   !! water, mass only) the packed area `aim` (the "mass" field `m`) and per
   !! category `n >= 1` the intensive tracers `trm(:,:,nt,n)` of
   !! `rdb_ice_transport_pack` are given a LINEAR profile inside each cell,
   !! in SCALED coordinates: the cell is the unit square `x, y in [-1/2, 1/2]`
   !! about its geometric centre, so the gradients are per cell width (index
   !! space) and NO metric factor enters this module: the cell lengths and
   !! areas belong to the departure-point and flux-integration steps.  A
   !! consequence is that the reconstruction is exact for fields linear in the
   !! INDEX coordinates, and on a stretched grid a field linear in physical
   !! distance is reproduced only to the order of the stretching (second
   !! order), as in the reference.
   !!
   !! ```
   !!    m(x, y)  = mc + mx x + my y                    (mass)
   !!    T(x, y)  = tc + tx x + ty y                    (tracer)
   !! ```
   !!
   !! **Masks** (`recon_masks_impl`).  `mmask(i,j,n) = 1` where
   !! `aim(i,j,n) > puny` (open water included), else `0`.  For a tracer that
   !! HAS DEPENDENTS (`hi`, `hs` in the shipped set) `tmask(i,j,nt,n) = 1`
   !! where `mmask = 1` and `|trm| > puny`; every other `tmask` is `0`.  The
   !! masks are `integer(int8)` 0/1 (the reference carries them as real 0/1
   !! and blends arithmetically; the selection used here is the same value for
   !! every finite neighbour).
   !!
   !! **Limited gradient** (`recon_limited_gradient`, Dukowicz and Baumgardner
   !! 2000 section 4, Lipscomb and Hunke 2004 section 2).  For a field `phi`
   !! with validity mask `M`, in a cell where `M` holds: the eight neighbours
   !! with `M = 0` are REPLACED by the cell's own value (so a masked neighbour
   !! contributes a zero difference and the stencil never sees a meaningless
   !! value); centred differences `gx = (phi_e - phi_w)/2`,
   !! `gy = (phi_n - phi_s)/2`; the extreme deviations of the neighbourhood
   !! `pmn = min(3x3) - phi`, `pmx = max(3x3) - phi`; the deviations of the
   !! linear profile at the four cell corners from the cell VALUE, with the
   !! value located at `(cnx, cny)` (the centroid of the weight the value is
   !! the mean of): `w = (+-1/2 - cnx) gx + (+-1/2 - cny) gy`, `qmn`, `qmx` their
   !! extremes; and the single limiter
   !! `alpha = min(a1, a2)`, `a1 = max(0, pmn/qmn)` if `|qmn| > |pmn|` else
   !! `1`, `a2 = max(0, pmx/qmx)` if `|qmx| > |pmx|` else `1` (a quotient is
   !! formed only when its denominator is strictly the larger magnitude, so
   !! it never divides by zero).  Result `alpha*gx`, `alpha*gy`; zero
   !! where `M = 0`.  Consequence (the monotonicity bound): the profile at
   !! every point of the cell lies within `[min, max]` of the nine cell values,
   !! to round-off; a local extremum gets a zero gradient.
   !!
   !! **The reconstruction (`construct_fields`).**  Mass: `phi = aim`, mask =
   !! the land mask `hm` (an ice-free ocean neighbour counts as value zero:
   !! its mass IS zero), `cn = 0`; `mc = aim` where `aim > puny` else `0`.
   !! The centre of mass of the cell is `mxav = mx/(12 mc)`, `myav = my/(12 mc)`
   !! (`xxav = yyav = 1/12` is the second moment of the unit square).
   !! Tracers, in index order (parents first):
   !!
   !! * type 1 (`depend = 0`, a density of area): mask `mmask`, `cn = (mxav,
   !!   myav)`, `tc = tm - tx mxav - ty myav`; when the tracer has dependents
   !!   also the centroid of area times tracer
   !!   `mtxav = (mc tx + mx tc) xxav/(mm tm)`, `mtyav` likewise (zero where
   !!   `tmask = 0`);
   !! * type 2 (`depend = parent`, a density of area times the parent): mask
   !!   `tmask(parent)`, `cn = (mtxav, mtyav)` OF THE PARENT in the same cell,
   !!   `tc = tm - tx mtxav(parent) - ty mtyav(parent)`;
   !! * type 3 is not carried: the tracer table builder refuses it.
   !!
   !! Every value is nonzero only where `aim > puny` (otherwise all of
   !! `mc, tc, tx, ty` are written `0`); `mx, my` follow the mass mask `hm`
   !! alone, exactly as the reference (they vanish where `aim = 0`, and the
   !! flux integration reads them through `mc`).
   !!
   !! **Identities (the mean-preservation contract the flux integration relies
   !! on).**  With `<.>` the mean over the unit cell, in every cell with
   !! `aim > puny`: `<m> = mc = aim`; for a type-1 tracer
   !! `<m T> = aim tm`; for a type-2 tracer whose parent has `tmask = 1`,
   !! `<m T_parent T> = aim tm_parent tm`.  (Products of linear profiles:
   !! `<x^2> = <y^2> = 1/12`, `<x> = <y> = <xy> = 0`.)  They hold for any
   !! gradients, limited or not; the tests integrate with a quadrature that
   !! shares no code with this module.
   !!
   !! **Layout of the outputs** (`remap_recon_t`; every array `(nx, ny, ...)`
   !! with the packed state's extents, horizontal indices first):
   !!
   !! | array | shape | meaning |
   !! |---|---|---|
   !! | `mmask` | `(nx, ny, 0:ncat)` int8 | area present |
   !! | `tmask` | `(nx, ny, ntrace, ncat)` int8 | parent tracer present (has-dependents slots only) |
   !! | `mc`, `mx`, `my` | `(nx, ny, 0:ncat)` | area at the cell centre, gradients per cell width |
   !! | `tc`, `tx`, `ty` | `(nx, ny, ntrace, ncat)` | tracer likewise (`ntrace` index = the `trm` index; category last) |
   !!
   !! **Halo contract.**  The outputs are written on the cells `is:ie, js:je`
   !! (every entry of every output there, zeros included) and nothing else is
   !! touched.  The masks are written on `is-1:ie+1, js-1:je+1`.  The INPUTS
   !! `aim`, `trm` and `hm` must be valid on `is-1:ie+1, js-1:je+1`; the range
   !! must satisfy `is >= 2`, `ie <= nx-1`, `js >= 2`, `je <= ny-1`.  So the
   !! reconstruction is valid on the range `is:ie, js:je`, which the caller sets to the
   !! owned cells plus one ring (`is:ie` = owned `+- 1`), because the flux
   !! integration of an edge reads the cells on both of its sides (and the
   !! diagonal ones for a departure region crossing a corner).  The packed
   !! fields (and `hm`) must therefore be valid on owned `+- 2` rings: THE
   !! MINIMUM VALID GHOST WIDTH OF THE PACKED FIELDS IS 2 RINGS.
   !!
   !! **Contracts the caller owns.**  (1) The land mask `hm` must be EXACTLY
   !! 0 or 1: the reference blends `hm*phi + (1 - hm)*phi_centre`, this port
   !! selects on `hm > puny`, and the two agree only for 0/1.  The masks
   !! kernel counts finite `hm` values other than 0 and 1 on the whole read
   !! ring and `remap_reconstruct` refuses (`RECON_CFG_MASK`, no
   !! reconstruction launched).  (2) The non-finite count `nbad` covers the
   !! whole READ ring, land cells included.  A NaN area fails every `>` test
   !! and would come out as an ice-free cell (zero outputs), so the caller
   !! (the flux-integration driver) MUST turn `nbad > 0` into the fatal
   !! `ICE_ERR_NONFINITE` bit; nothing else flags it.
   !!
   !! **GPU shape.**  Two launches, both `do concurrent` over
   !! `(n = 0:ncat, j, i)` on persistent device buffers: (1) the masks (they
   !! need only the packed fields, and the next stage reads NEIGHBOURS' masks,
   !! so they cannot be formed in the same launch); (2) the reconstruction,
   !! which reads neighbours' VALUES and masks and never a neighbour's
   !! reconstruction.  The tracer loop is sequential INSIDE the thread
   !! (hierarchy order): a child needs only its parent's quantities in the
   !! SAME cell, which the thread re-derives from the parent's already-written
   !! `tc, tx, ty, mc, mx, my` (the identical expression, so no centroid
   !! array exists); one launch per hierarchy level would cost a launch and a
   !! full re-read per level for no extra parallelism, the tracer count being
   !! far smaller than the cell count.  No stencil array is formed: the nine
   !! values go through scalars into `recon_limited_gradient`
   !! (`!$acc routine seq`).  Non-finite packed input is counted (`nbad`,
   !! `ICE_ERR_NONFINITE`) because a NaN fails every `>` test and would pass
   !! as an ice-free cell.
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE
   use rdb_ice_transport_pack, only: transport_workspace_t
   implicit none
   private

   public :: RECON_XXAV, RECON_CFG_OK, RECON_CFG_TABLE, RECON_CFG_WORKSPACE, RECON_CFG_RANGE, RECON_CFG_MASK
   public :: remap_recon_t
   public :: remap_reconstruct, recon_masks_impl, recon_fields_impl
   public :: recon_limited_gradient, recon_table_ok

   real(wp), parameter :: RECON_XXAV = 1.0_wp/12.0_wp
      !! Mean of `x^2` (and `y^2`) over the unit cell.

   integer, parameter :: RECON_CFG_OK = 0
      !! Configuration accepted.
   integer, parameter :: RECON_CFG_TABLE = 1
      !! The tracer table has a tracer that is not type 1 or type 2, or a type-2 tracer whose parent is not an earlier type-1 tracer with dependents.
   integer, parameter :: RECON_CFG_WORKSPACE = 2
      !! A workspace is not initialised, or the packed workspace and the reconstruction workspace differ in extents.
   integer, parameter :: RECON_CFG_MASK = 4
      !! The land mask holds a finite value that is neither 0 nor 1 (checked on the whole read ring; nothing reconstructed).
   integer, parameter :: RECON_CFG_RANGE = 3
      !! The range leaves no one-cell input ring inside the arrays.

   type :: remap_recon_t
      !! Persistent masks and reconstruction (layout: module docstring).
      !! OWNERSHIP: the transport driver slot owns one; `init` allocates
      !! (zeroed), `enter_data` maps the component arrays ONCE (never the
      !! aggregate), `exit_data` unmaps, `destroy` frees.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      integer :: ncat = 0
      integer :: ntrace = 0
      integer(int8), allocatable :: mmask(:, :, :)
         !! `(nx, ny, 0:ncat)`
      integer(int8), allocatable :: tmask(:, :, :, :)
         !! `(nx, ny, ntrace, ncat)`
      real(wp), allocatable :: mc(:, :, :)
         !! `(nx, ny, 0:ncat)`
      real(wp), allocatable :: mx(:, :, :)
      real(wp), allocatable :: my(:, :, :)
      real(wp), allocatable :: tc(:, :, :, :)
         !! `(nx, ny, ntrace, ncat)`
      real(wp), allocatable :: tx(:, :, :, :)
      real(wp), allocatable :: ty(:, :, :, :)
   contains
      procedure :: init => recon_init
      procedure :: enter_data => recon_enter_data
      procedure :: exit_data => recon_exit_data
      procedure :: destroy => recon_destroy
      procedure :: bytes => recon_bytes
   end type remap_recon_t

contains

   ! ------------------------------------------------------------------
   ! Workspace
   ! ------------------------------------------------------------------

   subroutine recon_init(this, nx, ny, ncat, ntrace, ierr)
      !! Allocate (zeroed) for a packed state of extent `(nx, ny, ncat)` with
      !! `ntrace` tracers.  `ierr /= 0`: invalid extent, nothing allocated.
      class(remap_recon_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, ntrace
      integer, intent(out) :: ierr

      ierr = 0
      if (nx < 3 .or. ny < 3 .or. ncat < 1 .or. ntrace < 1) then
         ierr = RECON_CFG_WORKSPACE
         return
      end if
      call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%ntrace = ntrace
      allocate (this%mmask(nx, ny, 0:ncat), this%tmask(nx, ny, ntrace, ncat))
      allocate (this%mc(nx, ny, 0:ncat), this%mx(nx, ny, 0:ncat), this%my(nx, ny, 0:ncat))
      allocate (this%tc(nx, ny, ntrace, ncat), this%tx(nx, ny, ntrace, ncat), this%ty(nx, ny, ntrace, ncat))
      this%mmask = 0_int8
      this%tmask = 0_int8
      this%mc = 0.0_wp
      this%mx = 0.0_wp
      this%my = 0.0_wp
      this%tc = 0.0_wp
      this%tx = 0.0_wp
      this%ty = 0.0_wp
      this%is_init = .true.
   end subroutine recon_init

   subroutine recon_enter_data(this)
      !! Map the component arrays once (`copyin`: the host zeros).
      class(remap_recon_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%mmask, this%tmask)
      !$acc enter data copyin(this%mc, this%mx, this%my)
      !$acc enter data copyin(this%tc, this%tx, this%ty)
      this%mapped = .true.
   end subroutine recon_enter_data

   subroutine recon_exit_data(this)
      class(remap_recon_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%mmask, this%tmask)
      !$acc exit data delete(this%mc, this%mx, this%my)
      !$acc exit data delete(this%tc, this%tx, this%ty)
      this%mapped = .false.
   end subroutine recon_exit_data

   subroutine recon_destroy(this)
      class(remap_recon_t), intent(inout) :: this
      call this%exit_data()
      if (allocated(this%mmask)) deallocate (this%mmask)
      if (allocated(this%tmask)) deallocate (this%tmask)
      if (allocated(this%mc)) deallocate (this%mc)
      if (allocated(this%mx)) deallocate (this%mx)
      if (allocated(this%my)) deallocate (this%my)
      if (allocated(this%tc)) deallocate (this%tc)
      if (allocated(this%tx)) deallocate (this%tx)
      if (allocated(this%ty)) deallocate (this%ty)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%ntrace = 0
   end subroutine recon_destroy

   pure function recon_bytes(this) result(b)
      !! Bytes of the buffers (host copy; the device mirror is the same size).
      class(remap_recon_t), intent(in) :: this
      integer(int64) :: b
      b = 0_int64
      if (allocated(this%mmask)) b = b + int(size(this%mmask), int64)*int(storage_size(this%mmask)/8, int64)
      if (allocated(this%tmask)) b = b + int(size(this%tmask), int64)*int(storage_size(this%tmask)/8, int64)
      if (allocated(this%mc)) b = b + int(size(this%mc), int64)*int(storage_size(this%mc)/8, int64)
      if (allocated(this%mx)) b = b + int(size(this%mx), int64)*int(storage_size(this%mx)/8, int64)
      if (allocated(this%my)) b = b + int(size(this%my), int64)*int(storage_size(this%my)/8, int64)
      if (allocated(this%tc)) b = b + int(size(this%tc), int64)*int(storage_size(this%tc)/8, int64)
      if (allocated(this%tx)) b = b + int(size(this%tx), int64)*int(storage_size(this%tx)/8, int64)
      if (allocated(this%ty)) b = b + int(size(this%ty), int64)*int(storage_size(this%ty)/8, int64)
   end function recon_bytes

   ! ------------------------------------------------------------------
   ! Host shim
   ! ------------------------------------------------------------------

   pure function recon_table_ok(ws) result(ok)
      !! `.true.` when every tracer is type 1, or type 2 with a parent that is an
      !! EARLIER type-1 tracer with dependents (the hierarchy the kernel carries).
      type(transport_workspace_t), intent(in) :: ws
      logical :: ok
      integer :: nt, p

      ok = ws%is_init
      if (.not. ok) return
      do nt = 1, ws%ntrace
         select case (ws%tab%tracer_type(nt))
         case (1)
            if (ws%tab%depend(nt) /= 0) ok = .false.
         case (2)
            p = ws%tab%depend(nt)
            if (p < 1 .or. p >= nt) then
               ok = .false.
            else if (ws%tab%tracer_type(p) /= 1 .or. .not. ws%tab%has_dependents(p)) then
               ok = .false.
            end if
         case default
            ok = .false.
         end select
      end do
   end function recon_table_ok

   subroutine remap_reconstruct(ws, rw, is, ie, js, je, puny, hm, nbad, cfg)
      !! Masks and reconstruction of the packed fields of `ws` into `rw` over
      !! the cells `is:ie, js:je` (two launches; see the module docstring for
      !! the input halo).  `hm` is the 0/1 land mask (`> puny` is ocean).
      !! `nbad` counts non-finite inputs read on the whole read ring: the
      !! caller MUST turn `nbad > 0` into `ICE_ERR_NONFINITE` (see the module
      !! docstring).  `hm` must be exactly 0/1 (`RECON_CFG_MASK` otherwise).
      !! `cfg /= RECON_CFG_OK`: the reconstruction was not run.  Every array
      !! must be device-present on an offload build.
      type(transport_workspace_t), intent(in) :: ws
      type(remap_recon_t), intent(inout) :: rw
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: hm(ws%nx, ws%ny)
      integer, intent(out) :: nbad, cfg
      integer :: nhm

      nbad = 0
      cfg = RECON_CFG_OK
      if (.not. (ws%is_init .and. rw%is_init)) then
         cfg = RECON_CFG_WORKSPACE
         return
      end if
      if (ws%nx /= rw%nx .or. ws%ny /= rw%ny .or. ws%ncat /= rw%ncat .or. ws%ntrace /= rw%ntrace) then
         cfg = RECON_CFG_WORKSPACE
         return
      end if
      if (.not. recon_table_ok(ws)) then
         cfg = RECON_CFG_TABLE
         return
      end if
      if (is < 2 .or. ie > ws%nx - 1 .or. js < 2 .or. je > ws%ny - 1 .or. is > ie .or. js > je) then
         cfg = RECON_CFG_RANGE
         return
      end if
      call recon_masks_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, is - 1, ie + 1, js - 1, je + 1, puny, &
                            ws%tab%has_dependents(1:ws%ntrace), hm, ws%aim, ws%trm, rw%mmask, rw%tmask, nbad, nhm)
      if (nhm > 0) then
         cfg = RECON_CFG_MASK
         return
      end if
      call recon_fields_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, is, ie, js, je, puny, &
                             ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), hm, ws%aim, ws%trm, &
                             rw%mmask, rw%tmask, rw%mc, rw%mx, rw%my, rw%tc, rw%tx, rw%ty)
   end subroutine remap_reconstruct

   ! ------------------------------------------------------------------
   ! Kernels
   ! ------------------------------------------------------------------

   pure subroutine recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, cnx, cny, gx, gy)
      !! Limited gradient of one cell from its value `pc` and the eight
      !! neighbour values (already substituted by `pc` where masked), with the
      !! value located at `(cnx, cny)` in the unit cell.  Scalars only.
      !$acc routine seq
      real(wp), intent(in) :: pc, pnw, pn, pne, pw, pe, psw, ps, pse, cnx, cny
      real(wp), intent(out) :: gx, gy
      real(wp) :: gxt, gyt, pmn, pmx, qmn, qmx, w1, w2, w3, w4

      gxt = (pe - pw)*0.5_wp
      gyt = (pn - ps)*0.5_wp
      pmn = min(pnw, pn, pne, pw, pc, pe, psw, ps, pse)
      pmx = max(pnw, pn, pne, pw, pc, pe, psw, ps, pse)
      pmn = pmn - pc
      pmx = pmx - pc
      w1 = (0.5_wp - cnx)*gxt + (0.5_wp - cny)*gyt
      w2 = (0.5_wp - cnx)*gxt - (0.5_wp + cny)*gyt
      w3 = -(0.5_wp + cnx)*gxt - (0.5_wp + cny)*gyt
      w4 = (0.5_wp - cny)*gyt - (0.5_wp + cnx)*gxt
      qmn = min(w1, w2, w3, w4)
      qmx = max(w1, w2, w3, w4)
      if (abs(qmn) > abs(pmn)) then
         w1 = max(0.0_wp, pmn/qmn)
      else
         w1 = 1.0_wp
      end if
      if (abs(qmx) > abs(pmx)) then
         w2 = max(0.0_wp, pmx/qmx)
      else
         w2 = 1.0_wp
      end if
      w1 = min(w1, w2)
      gx = w1*gxt
      gy = w1*gyt
   end subroutine recon_limited_gradient

   pure subroutine recon_masks_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, has_dep, hm, aim, trm, mmask, tmask, nbad, nhm)
      !! Flat kernel of the masks: ONE `do concurrent` over `(n = 0:ncat, j, i)`
      !! on `is:ie, js:je`; the tracer loop is inside the thread.  `nbad`
      !! counts non-finite `aim`, `trm` (categories `>= 1`) and `hm` read;
      !! `nhm` the finite `hm` values that are neither 0 nor 1.
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      logical, intent(in) :: has_dep(ntrace)
      real(wp), intent(in) :: hm(nx, ny), aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int8), intent(inout) :: mmask(nx, ny, 0:ncat), tmask(nx, ny, ntrace, ncat)
      integer, intent(out) :: nbad, nhm
      integer :: i, j, n, nt
      real(wp) :: a, tv

      nbad = 0
      nhm = 0
      !$acc data present(has_dep, hm, aim, trm, mmask, tmask)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(a, tv, nt) reduce(+:nbad, nhm)
         a = aim(i, j, n)
         if (.not. ieee_is_finite(a)) nbad = nbad + 1
         if (n == 0) then
            if (.not. ieee_is_finite(hm(i, j))) then
               nbad = nbad + 1
            else if (hm(i, j) /= 0.0_wp .and. hm(i, j) /= 1.0_wp) then
               nhm = nhm + 1
            end if
         end if
         if (a > puny) then
            mmask(i, j, n) = 1_int8
         else
            mmask(i, j, n) = 0_int8
         end if
         if (n >= 1) then
            !$acc loop seq
            do nt = 1, ntrace
               tv = trm(i, j, nt, n)
               if (.not. ieee_is_finite(tv)) nbad = nbad + 1
               tmask(i, j, nt, n) = 0_int8
               if (has_dep(nt) .and. a > puny) then
                  if (abs(tv) > puny) tmask(i, j, nt, n) = 1_int8
               end if
            end do
         end if
      end do
      !$acc end data
   end subroutine recon_masks_impl

   pure subroutine recon_fields_impl(nx, ny, ncat, ntrace, is, ie, js, je, puny, tracer_type, depend, hm, aim, trm, &
                                     mmask, tmask, mc, mx, my, tc, tx, ty)
      !! Flat kernel of the reconstruction: ONE `do concurrent` over
      !! `(n = 0:ncat, j, i)` on `is:ie, js:je`, reading the masks on the
      !! one-cell-wider range.  See the module docstring.
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: hm(nx, ny), aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer(int8), intent(in) :: mmask(nx, ny, 0:ncat), tmask(nx, ny, ntrace, ncat)
      real(wp), intent(inout) :: mc(nx, ny, 0:ncat), mx(nx, ny, 0:ncat), my(nx, ny, 0:ncat)
      real(wp), intent(inout) :: tc(nx, ny, ntrace, ncat), tx(nx, ny, ntrace, ncat), ty(nx, ny, ntrace, ncat)
      integer :: i, j, n, nt, p
      real(wp) :: mm, mcv, gx, gy, mxav, myav, tm, tgx, tgy, tcv, mtx, mty, w2, w3, w7
      real(wp) :: pc, pnw, pn, pne, pw, pe, psw, ps, pse

      !$acc data present(tracer_type, depend, hm, aim, trm, mmask, tmask, mc, mx, my, tc, tx, ty)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(mm, mcv, gx, gy, mxav, myav, tm, tgx, tgy, tcv, mtx, mty, &
                                                      w2, w3, w7, pc, pnw, pn, pne, pw, pe, psw, ps, pse, nt, p)
         mm = aim(i, j, n)
         ! mass: validity mask = land mask
         gx = 0.0_wp
         gy = 0.0_wp
         if (hm(i, j) > puny) then
            pc = mm
            pnw = merge(aim(i - 1, j + 1, n), pc, hm(i - 1, j + 1) > puny)
            pn = merge(aim(i, j + 1, n), pc, hm(i, j + 1) > puny)
            pne = merge(aim(i + 1, j + 1, n), pc, hm(i + 1, j + 1) > puny)
            pw = merge(aim(i - 1, j, n), pc, hm(i - 1, j) > puny)
            pe = merge(aim(i + 1, j, n), pc, hm(i + 1, j) > puny)
            psw = merge(aim(i - 1, j - 1, n), pc, hm(i - 1, j - 1) > puny)
            ps = merge(aim(i, j - 1, n), pc, hm(i, j - 1) > puny)
            pse = merge(aim(i + 1, j - 1, n), pc, hm(i + 1, j - 1) > puny)
            call recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, 0.0_wp, 0.0_wp, gx, gy)
         end if
         mx(i, j, n) = gx
         my(i, j, n) = gy
         mcv = 0.0_wp
         if (mm > puny) mcv = mm
         mc(i, j, n) = mcv
         if (n >= 1) then
            mxav = 0.0_wp
            myav = 0.0_wp
            if (mm > puny) then
               mxav = gx*RECON_XXAV/mm
               myav = gy*RECON_XXAV/mm
            end if
            !$acc loop seq
            do nt = 1, ntrace
               tm = trm(i, j, nt, n)
               tgx = 0.0_wp
               tgy = 0.0_wp
               tcv = 0.0_wp
               if (mm > puny) then
                  if (tracer_type(nt) == 1) then
                     pc = tm
                     pnw = merge(trm(i - 1, j + 1, nt, n), pc, mmask(i - 1, j + 1, n) /= 0_int8)
                     pn = merge(trm(i, j + 1, nt, n), pc, mmask(i, j + 1, n) /= 0_int8)
                     pne = merge(trm(i + 1, j + 1, nt, n), pc, mmask(i + 1, j + 1, n) /= 0_int8)
                     pw = merge(trm(i - 1, j, nt, n), pc, mmask(i - 1, j, n) /= 0_int8)
                     pe = merge(trm(i + 1, j, nt, n), pc, mmask(i + 1, j, n) /= 0_int8)
                     psw = merge(trm(i - 1, j - 1, nt, n), pc, mmask(i - 1, j - 1, n) /= 0_int8)
                     ps = merge(trm(i, j - 1, nt, n), pc, mmask(i, j - 1, n) /= 0_int8)
                     pse = merge(trm(i + 1, j - 1, nt, n), pc, mmask(i + 1, j - 1, n) /= 0_int8)
                     call recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, mxav, myav, tgx, tgy)
                     tcv = tm - tgx*mxav - tgy*myav
                  else
                     ! type 2: the parent's centroid of (area x tracer) in this cell
                     p = depend(nt)
                     mtx = 0.0_wp
                     mty = 0.0_wp
                     if (tmask(i, j, p, n) /= 0_int8) then
                        w2 = mcv*tx(i, j, p, n) + gx*tc(i, j, p, n)
                        w3 = mcv*ty(i, j, p, n) + gy*tc(i, j, p, n)
                        w7 = 1.0_wp/(mm*trm(i, j, p, n))
                        mtx = w2*RECON_XXAV*w7
                        mty = w3*RECON_XXAV*w7
                        pc = tm
                        pnw = merge(trm(i - 1, j + 1, nt, n), pc, tmask(i - 1, j + 1, p, n) /= 0_int8)
                        pn = merge(trm(i, j + 1, nt, n), pc, tmask(i, j + 1, p, n) /= 0_int8)
                        pne = merge(trm(i + 1, j + 1, nt, n), pc, tmask(i + 1, j + 1, p, n) /= 0_int8)
                        pw = merge(trm(i - 1, j, nt, n), pc, tmask(i - 1, j, p, n) /= 0_int8)
                        pe = merge(trm(i + 1, j, nt, n), pc, tmask(i + 1, j, p, n) /= 0_int8)
                        psw = merge(trm(i - 1, j - 1, nt, n), pc, tmask(i - 1, j - 1, p, n) /= 0_int8)
                        ps = merge(trm(i, j - 1, nt, n), pc, tmask(i, j - 1, p, n) /= 0_int8)
                        pse = merge(trm(i + 1, j - 1, nt, n), pc, tmask(i + 1, j - 1, p, n) /= 0_int8)
                        call recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, mtx, mty, tgx, tgy)
                     end if
                     tcv = tm - tgx*mtx - tgy*mty
                  end if
               end if
               tx(i, j, nt, n) = tgx
               ty(i, j, nt, n) = tgy
               tc(i, j, nt, n) = tcv
            end do
         end if
      end do
      !$acc end data
   end subroutine recon_fields_impl

end module rdb_ice_remap_recon
