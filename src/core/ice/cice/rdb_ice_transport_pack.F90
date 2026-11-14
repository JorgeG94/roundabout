!! Sea-ice transport: tracer dependency table and the pack/unpack between the
!! category state and the transport fields (packing convention of incremental
!! remapping, Lipscomb and Hunke 2004; tracer hierarchy of Dukowicz and
!! Baumgardner 2000).
module rdb_ice_transport_pack
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_driver (init_transport tracer table, state_to_tracers, tracers_to_state).
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
   !! **What is transported.**  The horizontal transport advects, per category
   !! `n = 1..ncat`, the category area `aim(:,:,n)` and the INTENSIVE tracers
   !! `trm(:,:,nt,n)`; open water `aim(:,:,0)` is a mass-only category.  The
   !! tracer list is `nt = 1` ice thickness `hi = vicen/aicen`, `nt = 2` snow
   !! thickness `hs = vsnon/aicen`, then the `ntrcr` state tracers
   !! (`Tsfc | qice(1:nilyr) | qsno | sice(1:nilyr)`, the packed order of
   !! `rdb_ice_params`), so `ntrace = 2 + ntrcr`.
   !!
   !! **The dependency hierarchy** (`transport_tracer_table_t`): a tracer is
   !! a density of the field it depends on (`depend`): `0` = area (type 1:
   !! `aim*T` is conserved: `hi`, `hs`, `Tsfc`), `1` = ice volume (type 2:
   !! `aim*hi*T`: ice enthalpy, ice salinity), `2` = snow volume (type 2:
   !! `aim*hs*T`: snow enthalpy); a value `> 2` names the PARENT TRACER by its
   !! `trm` index (`2 + parent tracer`).  A parent must precede its dependent
   !! (the update divides by the already-updated parent), a parent may not
   !! itself depend on a tracer (type 3 is not carried), and counts above the
   !! fixed maxima are refused: `transport_tracer_table_from_depend` returns a
   !! nonzero `ierr` (`TT_ERR_*`) instead of building a half-valid table.
   !!
   !! **Layout of the packed buffers** (`transport_workspace_t`, mapped once):
   !!
   !! | array | shape | role |
   !! |---|---|---|
   !! | `aim` | `(nx, ny, 0:ncat)` | area; `0` is open water |
   !! | `trm` | `(nx, ny, ntrace, ncat)` | tracers; category LAST, so `trm(:,:,:,n)` is contiguous (a category-batched consumer passes it to an explicit-shape dummy, no copy) |
   !!
   !! Horizontal indices first (coalesced in `i`), the full local extent
   !! `(nx_total, ny_total)` of the state, so packed and state arrays share
   !! the SAME `(i, j)`.  The workspace is allocated and device-mapped ONCE
   !! (`init` + `enter_data`), never per call.  Pack and unpack work on the
   !! index range the caller gives (`is:ie, js:je`: the owned cells plus
   !! whatever ghost ring the advection scheme needs; entries outside it are
   !! not touched).
   !!
   !! **Pack** (per `(n, j, i)`, ONE `do concurrent` launch over
   !! `n = 0:ncat`): `aim = aicen`; where `aim > puny`: `hi = vicen*(1/aim)`,
   !! `hs = vsnon*(1/aim)`, `trm(2+it) = trcrn(it)` EXCEPT the snow enthalpy,
   !! which is shifted by `+rhos*Lfresh` ([J/m3], so it is `<= 0` and `0` at
   !! the melting point: the monotone limiter needs a sign-definite field);
   !! elsewhere `trm = 0`.  The reciprocal multiply is the reference
   !! arithmetic, so the packed values agree with it bit for bit.
   !!
   !! **Precondition: the entry state is already zapped.**  A category with
   !! `0 < aicen <= puny` (a sliver) loses its volumes and tracers in the pack,
   !! so ice mass and energy vanish from the transported fields.  The
   !! reference is safe only because the small-area cleanup kernel has already
   !! removed slivers, with the removed mass booked in the fluxes.  The driver
   !! MUST run that cleanup before the pack; under `transport_check` it is
   !! verified by `rdb_ice_transport_checks::transport_check_slivers` (run on
   !! the state before the pack: a conservation sum taken after the pack cannot
   !! see the loss).  Only `nslyr = 1` is supported (one snow-enthalpy tracer).
   !!
   !! **Unpack** (same shape): where `aim > 0`: `aicen = aim`,
   !! `vicen = aim*hi`, `vsnon = aim*hs`, `trcrn = trm` (snow enthalpy shifted
   !! back); where `aim <= 0`: `aicen = vicen = vsnon = 0` and the tracers are
   !! LEFT UNTOUCHED (stale values of an empty category are never read, and
   !! writing them would cost bandwidth).  So a category with
   !! `0 < aim <= puny` comes back with its area but zero volumes and zero
   !! tracers (the pack zeroed them): the reference behaviour, which only
   !! matters for sub-`puny` slivers.
   !!
   !! **Identity.**  `unpack(pack(state))` reproduces the state exactly for the
   !! area, for every tracer except the snow enthalpy and for any category
   !! with `aicen == 0` (and zero volumes); it reproduces `vicen`, `vsnon` to
   !! round-off (`aim*(v*(1/aim))` is not `v`: at most ~2 ulp) and the snow
   !! enthalpy to round-off of the shift (`(q + c) - c`, error of order
   !! `ulp(rhos*Lfresh)`, about 1e-16 relative to `q`).  The tests state the
   !! measured worst.
   !!
   !! **Signs and units.**  `aim` is a fraction of the cell area; `hi`, `hs`
   !! are m per unit CATEGORY area; `Tsfc` degC; `qice` J/m3 of ice;
   !! `qsno + rhos*Lfresh` J/m3 of snow; `sice` ppt.  Ice layers are top-down.
   !!
   !! **Errors.**  Non-finite input never disappears into the `> puny` /
   !! `> 0` branch: pack and unpack count every non-finite value they read
   !! (`nbad`, a one-accumulator `reduce(+)`) and return
   !! `ICE_ERR_NONFINITE` in `ierr`.  The per-cell attribution is the job of
   !! `rdb_ice_transport_checks`.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONE, ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ice_nt_qsno, ice_ntrcr, ice_trcr_get, ice_trcr_set, &
                             ice_trcr_depend
   implicit none
   private

   public :: ICE_TRANSPORT_NTRACE_MAX
   public :: TT_OK, TT_ERR_COUNT, TT_ERR_UNKNOWN_DEP, TT_ERR_ORDER, TT_ERR_TYPE3
   public :: transport_tracer_table_t
   public :: transport_tracer_table_from_depend, transport_tracer_table_build
   public :: transport_tracer_table_message
   public :: transport_workspace_t
   public :: transport_pack, transport_unpack
   public :: transport_pack_impl, transport_unpack_impl

   integer, parameter :: ICE_TRANSPORT_NTRACE_MAX = 2 + ICE_NTRCR_MAX
      !! Maximum number of transported tracers (`hi`, `hs` and the state tracers).

   integer, parameter :: TT_OK = 0
      !! Table built.
   integer, parameter :: TT_ERR_COUNT = 1
      !! Tracer count outside `1..ICE_NTRCR_MAX` (or a dependency list shorter than the count).
   integer, parameter :: TT_ERR_UNKNOWN_DEP = 2
      !! A dependency code that names no field and no tracer (negative, or beyond the list).
   integer, parameter :: TT_ERR_ORDER = 3
      !! A dependent listed at or before its parent.
   integer, parameter :: TT_ERR_TYPE3 = 4
      !! A parent that itself depends on a tracer (type 3: not carried).

   type :: transport_tracer_table_t
      !! Dependency table of the transported tracers (plain integers: copyable
      !! to the device, passed to kernels as explicit-shape arrays).  Indices
      !! are `trm` indices: `1` = `hi`, `2` = `hs`, `2 + it` = state tracer `it`.
      !! No default initialisers; `transport_tracer_table_from_depend` sets
      !! every component (entries beyond `ntrace` are zero).
      integer :: ierr                 !! `TT_OK` or the first `TT_ERR_*`
      integer :: bad_tracer           !! state-tracer index the error refers to (0 if none)
      integer :: ntrcr                !! state tracers
      integer :: ntrace               !! `2 + ntrcr`
      integer :: nt_qsno              !! state-tracer index of the first snow enthalpy (0: none)
      integer :: nslyr                !! snow layers (snow enthalpy tracers `nt_qsno..nt_qsno+nslyr-1`)
      integer :: tracer_type(ICE_TRANSPORT_NTRACE_MAX)    !! 1 = depends on area, 2 = depends on one tracer/field
      integer :: depend(ICE_TRANSPORT_NTRACE_MAX)         !! `trm` index of the parent (`1` hi, `2` hs), 0 = area
      integer :: base_field(ICE_TRANSPORT_NTRACE_MAX)     !! 0 area, 1 ice volume, 2 snow volume
      logical :: has_dependents(ICE_TRANSPORT_NTRACE_MAX)  !! some tracer lists this one as its parent
   end type transport_tracer_table_t

   type :: transport_workspace_t
      !! Persistent packed-transport buffers plus the table (see the module
      !! docstring for the layout).  OWNERSHIP: the transport driver slot owns
      !! one; `init` allocates, `enter_data` maps ONCE (component arrays and
      !! the table, never the aggregate), `exit_data` unmaps, `destroy` frees.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
         !! Local extent in x, ghosts included (the state's `nx_total`).
      integer :: ny = 0
      integer :: ncat = 0
      integer :: ntrace = 0
      type(transport_tracer_table_t) :: tab
      real(wp), allocatable :: aim(:, :, :)
         !! `(nx, ny, 0:ncat)`
      real(wp), allocatable :: trm(:, :, :, :)
         !! `(nx, ny, ntrace, ncat)`
   contains
      procedure :: init => workspace_init
      procedure :: enter_data => workspace_enter_data
      procedure :: exit_data => workspace_exit_data
      procedure :: destroy => workspace_destroy
      procedure :: bytes => workspace_bytes
   end type transport_workspace_t

contains

   ! ------------------------------------------------------------------
   ! Tracer table
   ! ------------------------------------------------------------------

   pure subroutine transport_tracer_table_from_depend(ntrcr, trcr_depend, nt_qsno, nslyr, tab)
      !! Build the table from a dependency list `trcr_depend(1:ntrcr)` (`0`
      !! area, `1` ice volume, `2` snow volume, `> 2` the `trm` index of the
      !! parent tracer) and the position of the snow enthalpy.  The builder
      !! checks, in this order per tracer: count; code unknown; parent not
      !! before its dependent; parent itself dependent (type 3).  On failure
      !! `tab%ierr /= TT_OK`, `tab%bad_tracer` names the tracer, and the rest
      !! of the table is zero: a failed table is never usable by accident.
      integer, intent(in) :: ntrcr
      integer, intent(in) :: trcr_depend(:)
      integer, intent(in) :: nt_qsno, nslyr
      type(transport_tracer_table_t), intent(out) :: tab
      integer :: it, d, parent, nt

      tab%ierr = TT_OK
      tab%bad_tracer = 0
      tab%ntrcr = 0
      tab%ntrace = 0
      tab%nt_qsno = 0
      tab%nslyr = 0
      tab%tracer_type = 0
      tab%depend = 0
      tab%base_field = 0
      tab%has_dependents = .false.

      if (ntrcr < 1 .or. ntrcr > ICE_NTRCR_MAX .or. size(trcr_depend) < ntrcr) then
         tab%ierr = TT_ERR_COUNT
         return
      end if

      do it = 1, ntrcr
         d = trcr_depend(it)
         if (d < 0 .or. d > 2 + ntrcr) then
            call fail(tab, TT_ERR_UNKNOWN_DEP, it)
            return
         end if
         if (d > 2) then
            if (d >= 2 + it) then
               call fail(tab, TT_ERR_ORDER, it)
               return
            end if
            if (d > 2 .and. trcr_depend(d - 2) /= 0) then
               call fail(tab, TT_ERR_TYPE3, it)
               return
            end if
         end if
      end do

      tab%ntrcr = ntrcr
      tab%ntrace = 2 + ntrcr
      tab%nt_qsno = nt_qsno
      tab%nslyr = nslyr
      tab%tracer_type(1:2) = 1
      tab%depend(1:2) = 0
      tab%base_field(1:2) = 0
      do it = 1, ntrcr
         nt = 2 + it
         d = trcr_depend(it)
         tab%depend(nt) = d
         if (d == 0) then
            tab%tracer_type(nt) = 1
            tab%base_field(nt) = 0
         else
            tab%tracer_type(nt) = 2
            if (d <= 2) then
               tab%base_field(nt) = d
            else
               parent = d
               tab%base_field(nt) = tab%base_field(parent)
            end if
         end if
      end do
      do nt = 1, tab%ntrace
         if (tab%depend(nt) > 0) tab%has_dependents(tab%depend(nt)) = .true.
      end do

   contains

      pure subroutine fail(t, code, which)
         type(transport_tracer_table_t), intent(inout) :: t
         integer, intent(in) :: code, which
         t%ierr = code
         t%bad_tracer = which
         t%tracer_type = 0
         t%depend = 0
         t%base_field = 0
         t%has_dependents = .false.
      end subroutine fail

   end subroutine transport_tracer_table_from_depend

   pure subroutine transport_tracer_table_build(nilyr, nslyr, tab)
      !! The table of the shipped tracer set at `nilyr`/`nslyr` (dependency
      !! codes from `ice_trcr_depend`: `Tsfc` area, `qice`/`sice` ice volume,
      !! `qsno` snow volume).
      integer, intent(in) :: nilyr, nslyr
      type(transport_tracer_table_t), intent(out) :: tab
      integer :: ntrcr, it
      integer :: dep(ICE_NTRCR_MAX)

      dep = 0
      if (nilyr < 1 .or. nslyr < 0 .or. ice_ntrcr(nilyr, nslyr) > ICE_NTRCR_MAX) then
         call transport_tracer_table_from_depend(0, dep, 0, nslyr, tab)
         return
      end if
      ntrcr = ice_ntrcr(nilyr, nslyr)
      do it = 1, ntrcr
         dep(it) = ice_trcr_depend(it, nilyr, nslyr)
      end do
      call transport_tracer_table_from_depend(ntrcr, dep, ice_nt_qsno(nilyr), nslyr, tab)
   end subroutine transport_tracer_table_build

   function transport_tracer_table_message(tab) result(msg)
      !! Host-side refusal text for a failed table (empty when it is valid).
      type(transport_tracer_table_t), intent(in) :: tab
      character(len=160) :: msg
      character(len=12) :: s

      write (s, "(i0)") tab%bad_tracer
      select case (tab%ierr)
      case (TT_OK)
         msg = ""
      case (TT_ERR_COUNT)
         msg = "transport tracer table: tracer count outside 1..ICE_NTRCR_MAX"
      case (TT_ERR_UNKNOWN_DEP)
         msg = "transport tracer table: tracer "//trim(s)//" depends on an unknown field or tracer"
      case (TT_ERR_ORDER)
         msg = "transport tracer table: tracer "//trim(s)//" is listed at or before its parent"
      case (TT_ERR_TYPE3)
         msg = "transport tracer table: tracer "//trim(s)//" depends on a tracer that is itself "// &
               "dependent (type 3 is not carried)"
      case default
         msg = "transport tracer table: unknown error"
      end select
   end function transport_tracer_table_message

   ! ------------------------------------------------------------------
   ! Workspace
   ! ------------------------------------------------------------------

   subroutine workspace_init(this, nx, ny, ncat, tab, ierr)
      !! Allocate the packed buffers (zeroed) for a state of extent
      !! `(nx, ny, ncat)` and keep `tab`.  `ierr /= 0` for an invalid table
      !! or extent (nothing is allocated).
      class(transport_workspace_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat
      type(transport_tracer_table_t), intent(in) :: tab
      integer, intent(out) :: ierr

      ierr = 0
      if (tab%ierr /= TT_OK .or. nx < 1 .or. ny < 1 .or. ncat < 1) then
         ierr = TT_ERR_COUNT
         if (tab%ierr /= TT_OK) ierr = tab%ierr
         return
      end if
      call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%tab = tab
      this%ntrace = tab%ntrace
      allocate (this%aim(nx, ny, 0:ncat), this%trm(nx, ny, tab%ntrace, ncat))
      this%aim = 0.0_wp
      this%trm = 0.0_wp
      this%is_init = .true.
   end subroutine workspace_init

   subroutine workspace_enter_data(this)
      !! Map the buffers and the table once (`copyin`: host zeros).  The
      !! table is all static integers/logicals, so mapping it is safe; the
      !! aggregate (with its allocatable descriptors) is never mapped.
      class(transport_workspace_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%tab)
      !$acc enter data copyin(this%aim, this%trm)
      this%mapped = .true.
   end subroutine workspace_enter_data

   subroutine workspace_exit_data(this)
      class(transport_workspace_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%aim, this%trm)
      !$acc exit data delete(this%tab)
      this%mapped = .false.
   end subroutine workspace_exit_data

   pure function workspace_bytes(this) result(b)
      !! Bytes of the packed buffers (host copy; the device mirror is the same size).
      class(transport_workspace_t), intent(in) :: this
      integer(int64) :: b
      b = 0_int64
      if (allocated(this%aim)) b = b + int(size(this%aim), int64)*int(storage_size(this%aim)/8, int64)
      if (allocated(this%trm)) b = b + int(size(this%trm), int64)*int(storage_size(this%trm)/8, int64)
   end function workspace_bytes

   subroutine workspace_destroy(this)
      class(transport_workspace_t), intent(inout) :: this
      call this%exit_data()
      if (allocated(this%aim)) deallocate (this%aim)
      if (allocated(this%trm)) deallocate (this%trm)
      this%is_init = .false.
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%ntrace = 0
   end subroutine workspace_destroy

   ! ------------------------------------------------------------------
   ! Pack / unpack
   ! ------------------------------------------------------------------

   subroutine transport_pack(ws, c, nilyr, nslyr, is, ie, js, je, aice0, aicen, vicen, vsnon, &
                             Tsfcn, qicen, sicen, qsnon, ierr)
      !! Pack the state into `ws%aim`, `ws%trm` over `is:ie, js:je` (one
      !! launch).  `c` supplies `puny`, `rhos`, `Lfresh`.  `ierr` is
      !! `ICE_ERR_NONFINITE` when any read value was non-finite.  All arrays
      !! must be device-present on an offload build.
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nilyr, nslyr, is, ie, js, je
      real(wp), intent(in) :: aice0(ws%nx, ws%ny), aicen(ws%nx, ws%ny, ws%ncat), vicen(ws%nx, ws%ny, ws%ncat)
      real(wp), intent(in) :: vsnon(ws%nx, ws%ny, ws%ncat), Tsfcn(ws%nx, ws%ny, ws%ncat)
      real(wp), intent(in) :: qicen(ws%nx, ws%ny, ws%ncat, nilyr), sicen(ws%nx, ws%ny, ws%ncat, nilyr)
      real(wp), intent(in) :: qsnon(ws%nx, ws%ny, ws%ncat, nslyr)
      integer(int32), intent(out) :: ierr
      integer :: nbad
      real(wp) :: puny, qshift

      puny = c%puny
      qshift = c%rhos*c%Lfresh
      call transport_pack_impl(ws%nx, ws%ny, ws%ncat, nilyr, nslyr, ws%tab%ntrcr, ws%tab%nt_qsno, &
                               is, ie, js, je, puny, qshift, aice0, aicen, vicen, vsnon, &
                               Tsfcn, qicen, sicen, qsnon, ws%aim, ws%trm, nbad)
      ierr = ICE_ERR_NONE
      if (nbad > 0) ierr = ICE_ERR_NONFINITE
   end subroutine transport_pack

   subroutine transport_unpack(ws, c, nilyr, nslyr, is, ie, js, je, aice0, aicen, vicen, vsnon, &
                               Tsfcn, qicen, sicen, qsnon, ierr)
      !! Unpack `ws%aim`, `ws%trm` into the state over `is:ie, js:je` (one
      !! launch); see the module docstring for the empty-category rules.
      type(transport_workspace_t), intent(inout) :: ws
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nilyr, nslyr, is, ie, js, je
      real(wp), intent(inout) :: aice0(ws%nx, ws%ny), aicen(ws%nx, ws%ny, ws%ncat), vicen(ws%nx, ws%ny, ws%ncat)
      real(wp), intent(inout) :: vsnon(ws%nx, ws%ny, ws%ncat), Tsfcn(ws%nx, ws%ny, ws%ncat)
      real(wp), intent(inout) :: qicen(ws%nx, ws%ny, ws%ncat, nilyr), sicen(ws%nx, ws%ny, ws%ncat, nilyr)
      real(wp), intent(inout) :: qsnon(ws%nx, ws%ny, ws%ncat, nslyr)
      integer(int32), intent(out) :: ierr
      integer :: nbad
      real(wp) :: qshift

      qshift = c%rhos*c%Lfresh
      call transport_unpack_impl(ws%nx, ws%ny, ws%ncat, nilyr, nslyr, ws%tab%ntrcr, ws%tab%nt_qsno, &
                                 is, ie, js, je, qshift, ws%aim, ws%trm, aice0, aicen, vicen, vsnon, &
                                 Tsfcn, qicen, sicen, qsnon, nbad)
      ierr = ICE_ERR_NONE
      if (nbad > 0) ierr = ICE_ERR_NONFINITE
   end subroutine transport_unpack

   pure subroutine transport_pack_impl(nx, ny, ncat, nilyr, nslyr, ntrcr, nt_qsno, is, ie, js, je, &
                                       puny, qshift, aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, &
                                       qsnon, aim, trm, nbad)
      !! Flat kernel behind `transport_pack`: ONE `do concurrent` over
      !! `(n = 0:ncat, j, i)`; `n = 0` copies the open-water area.  The
      !! tracer loop runs inside the thread (its trip count is loop-uniform,
      !! so no divergence); the state tracers are read through the shared
      !! accessor `ice_trcr_get`, so the packed order has one source of truth.
      !! `qshift = rhos*Lfresh`.  `nbad` counts non-finite values read.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, ntrcr, nt_qsno, is, ie, js, je
      real(wp), intent(in) :: puny, qshift
      real(wp), intent(in) :: aice0(nx, ny), aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aim(nx, ny, 0:ncat), trm(nx, ny, 2 + ntrcr, ncat)
      integer, intent(out) :: nbad
      integer :: i, j, n, it
      real(wp) :: a, w1, tv

      nbad = 0
      !$acc data present(aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aim, trm)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(a, w1, tv, it) reduce(+:nbad)
         if (n == 0) then
            aim(i, j, 0) = aice0(i, j)
            if (.not. ieee_is_finite(aice0(i, j))) nbad = nbad + 1
         else
            a = aicen(i, j, n)
            aim(i, j, n) = a
            if (.not. ieee_is_finite(a)) nbad = nbad + 1
            if (.not. ieee_is_finite(vicen(i, j, n))) nbad = nbad + 1
            if (.not. ieee_is_finite(vsnon(i, j, n))) nbad = nbad + 1
            if (a > puny) then
               w1 = 1.0_wp/a
               trm(i, j, 1, n) = vicen(i, j, n)*w1
               trm(i, j, 2, n) = vsnon(i, j, n)*w1
               !$acc loop seq
               do it = 1, ntrcr
                  tv = ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon)
                  if (.not. ieee_is_finite(tv)) nbad = nbad + 1
                  if (it >= nt_qsno .and. it < nt_qsno + nslyr) tv = tv + qshift
                  trm(i, j, 2 + it, n) = tv
               end do
            else
               !$acc loop seq
               do it = 1, 2 + ntrcr
                  trm(i, j, it, n) = 0.0_wp
               end do
            end if
         end if
      end do
      !$acc end data
   end subroutine transport_pack_impl

   pure subroutine transport_unpack_impl(nx, ny, ncat, nilyr, nslyr, ntrcr, nt_qsno, is, ie, js, je, &
                                         qshift, aim, trm, aice0, aicen, vicen, vsnon, Tsfcn, qicen, &
                                         sicen, qsnon, nbad)
      !! Flat kernel behind `transport_unpack` (one launch over
      !! `(n = 0:ncat, j, i)`).  `nbad` counts non-finite `aim`, and
      !! non-finite `trm` where `aim > 0`.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, ntrcr, nt_qsno, is, ie, js, je
      real(wp), intent(in) :: qshift
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, 2 + ntrcr, ncat)
      real(wp), intent(inout) :: aice0(nx, ny), aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: qsnon(nx, ny, ncat, nslyr)
      integer, intent(out) :: nbad
      integer :: i, j, n, it
      real(wp) :: a, tv

      nbad = 0
      !$acc data present(aim, trm, aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(a, tv, it) reduce(+:nbad)
         if (n == 0) then
            aice0(i, j) = aim(i, j, 0)
            if (.not. ieee_is_finite(aim(i, j, 0))) nbad = nbad + 1
         else
            a = aim(i, j, n)
            if (.not. ieee_is_finite(a)) nbad = nbad + 1
            if (a > 0.0_wp) then
               aicen(i, j, n) = a
               vicen(i, j, n) = a*trm(i, j, 1, n)
               vsnon(i, j, n) = a*trm(i, j, 2, n)
               if (.not. ieee_is_finite(trm(i, j, 1, n))) nbad = nbad + 1
               if (.not. ieee_is_finite(trm(i, j, 2, n))) nbad = nbad + 1
               !$acc loop seq
               do it = 1, ntrcr
                  tv = trm(i, j, 2 + it, n)
                  if (.not. ieee_is_finite(tv)) nbad = nbad + 1
                  if (it >= nt_qsno .and. it < nt_qsno + nslyr) tv = tv - qshift
                  call ice_trcr_set(it, i, j, n, nx, ny, ncat, nilyr, nslyr, tv, Tsfcn, qicen, sicen, qsnon)
               end do
            else
               aicen(i, j, n) = 0.0_wp
               vicen(i, j, n) = 0.0_wp
               vsnon(i, j, n) = 0.0_wp
            end if
         end if
      end do
      !$acc end data
   end subroutine transport_unpack_impl

end module rdb_ice_transport_pack
