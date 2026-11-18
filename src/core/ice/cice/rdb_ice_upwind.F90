!! First-order donor-cell (upwind) advection of the packed sea-ice transport
!! fields on the C grid (Lipscomb and Hunke 2004 for the tracer hierarchy).
module rdb_ice_upwind
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: CICE 2d271f5, module ice_transport_driver (the upwind option: edge-velocity flux and cell update, tracer recovery).
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
   ! Physics: donor-cell scheme; Lipscomb and Hunke (2004, Mon. Wea. Rev. 132 1341)
   ! and Dukowicz and Baumgardner (2000, J. Comput. Phys. 160 318) for the
   ! tracer hierarchy.
   !
   !! **What is advected.**  The EXTENSIVE products of the packed fields
   !! (`rdb_ice_transport_pack`): per category the area `aim`, the ice volume
   !! `aim*hi` and snow volume `aim*hs` (the first two tracers), and for each
   !! state tracer the product with the field it is a density of: `aim*T`
   !! (`depend = 0`), `(aim*hi)*T` (`1`), `(aim*hs)*T` (`2`); open water
   !! `aim(:,:,0)` is a mass-only field.  Each product `P` obeys
   !! `P_new = P - ((Fe - Fw) + (Fn - Fs))/tarea` (the reference sums
   !! `((Fe - Fw) + Fn) - Fs`; the pairing here makes the arithmetic exactly
   !! mirror-symmetric in both axes, and differs from it by round-off) with the donor-cell face flux
   !! `F = (dt/2) L ((u+|u|) P_donor_west + (u-|u|) P_donor_east)`, `L` the face
   !! length and `u` the face-normal velocity times the face mask (a land face
   !! carries no flux).  The new intensive tracers are recovered from the new
   !! products, the new base first: `hi = vice/aim`, `hs = vsno/aim` (any
   !! `aim > 0`, so the volume of a sliver category is kept for the clean-up to
   !! book; only a vanishing `aim <= 1e-100` is not divided), an area tracer
   !! `T = P/aim` where `aim > puny` (else 0), a volume tracer `T = P/V` where
   !! `V > 0` (else 0).  Where a cell's net product flux AND its base flux are
   !! exactly zero the OLD intensive value is copied, so an empty domain and
   !! zero velocity are bitwise no-ops of the packed fields.  A volume tracer
   !! on a zero base volume is 0 in the PACKED frame, which for the snow
   !! enthalpy (stored shifted by `+rhos*Lfresh`) comes back from the unpack as
   !! `-rhos*Lfresh` (the melting point) where the reference writes 0: inert
   !! (no snow), and the physical one of the two.
   !!
   !! **Velocities.**  The C-grid face velocities are used DIRECTLY (no
   !! corner-to-edge averaging): `u(nx+1, ny)` with the EAST face of cell
   !! `(i,j)` at `u(i+1,j)` (its WEST face at `u(i,j)`), `v(nx, ny+1)` with
   !! the NORTH face at `v(i,j+1)`; `HTE(i,j)`, `HTN(i,j)` the east / north
   !! face lengths of cell `(i,j)`, `epm(i,j)`, `npm(i,j)` the 0/1 masks of
   !! those faces (`rdb_ice_cgrid_metrics`).
   !!
   !! **Ranges and halos.**  The kernel updates the cells `is:ie, js:je` and
   !! reads the packed fields one cell beyond (`is-1:ie+1, js-1:je+1`), the
   !! face velocities at `u(is:ie+1, js:je)`, `v(is:ie, js:je+1)` and the
   !! lengths / masks at `i = is-1:ie`, `j = js-1:je`; nothing outside
   !! `is:ie, js:je` of the OUTPUT is written (the ghost ring of the output is
   !! the caller's to refresh).  The packed input must therefore have been
   !! packed over the one-cell-wider range from a state whose ghost ring is
   !! valid.
   !!
   !! **GPU shape.**  ONE launch for the whole advection: `do concurrent`
   !! over `(n = 0:ncat, j, i)`, the tracer loop sequential inside the thread
   !! (it carries the hierarchy: a dependent tracer reads the new products of
   !! its base, which the same thread has just formed).  No flux array exists:
   !! each face flux is a pure function of the two cells it joins, evaluated
   !! by both neighbours from the OLD input (the output goes to a second
   !! workspace), so the two evaluations are the same expression and a
   !! neighbour's new value is never read.  The remaining per-thread state is
   !! scalars.  Only `depend <= 2` (area, ice volume, snow volume) is carried;
   !! a table with a dependent of a parent tracer (`depend > 2`, not in the
   !! shipped set) is refused (`cfg = UPWIND_CFG_DEPEND`).
   !!
   !! **Volumes of negligible areas.**  The products are the reference's
   !! extensive fields, so a category of tiny area keeps its volume and
   !! enthalpy through the advection (the driver packs with `thick_floor =
   !! UPWIND_THICK_FLOOR`); only the area tracer (the surface temperature) of
   !! a cell whose new area is `<= puny` is set to 0 (then `Tf`).  A new area
   !! in `(-puny, 0]` is unpacked as 0, with its volume: the content lost is
   !! at most `puny` times the thickness, round-off of the scheme and not
   !! booked.  Conservation is TO ROUND-OFF: one division by the cell area, the
   !! packed-frame re-multiply and, on the device, FMA contraction (which can
   !! make the two evaluations of a face flux differ in the last bit).
   !!
   !! **Counters.**  `nbad` counts non-finite new products (a NaN velocity or
   !! field reaches every product it touches, so it is never silently
   !! absorbed by the `> puny` branches), `nneg` the categories (open water
   !! included) whose new area is below `-puny`.  Positivity holds only inside
   !! `dt (|u|/dx + |v|/dy) <= 1` per cell; outside it `nneg` counts the
   !! violation and the new area stays as computed.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_transport_pack, only: transport_workspace_t
   implicit none
   private

   public :: UPWIND_CFG_OK, UPWIND_CFG_DEPEND, UPWIND_CFG_WORKSPACE, UPWIND_THICK_FLOOR, upwind_table_ok
   public :: upwind_advect, upwind_advect_impl
   public :: upwind_empty_tsfc, upwind_empty_tsfc_impl

   integer, parameter :: UPWIND_CFG_OK = 0
      !! Configuration accepted.
   integer, parameter :: UPWIND_CFG_DEPEND = 1
      !! The tracer table has a dependent of a parent tracer (`depend > 2`).
   integer, parameter :: UPWIND_CFG_WORKSPACE = 2
      !! A workspace is not initialised, or the two do not have the same extents/table.

   real(wp), parameter :: UPWIND_THICK_FLOOR = 1.0e-100_wp
      !! Smallest area the thickness quotient `P/aim` is formed for, and the
      !! `thick_floor` the driver gives the pack so the volume of a category of
      !! negligible area survives: it only keeps the quotient of a numerically
      !! cancelled cell (`aim` ~ 1e-17 of round-off) from overflowing; every
      !! physical sliver is far above it.

contains

   pure function upw(y1, y2, a, h, dt) result(f)
      !! Donor-cell flux through a face of length `h` with normal velocity
      !! `a` between a cell holding `y1` (behind the face) and one holding
      !! `y2` (ahead of it).  Scalars only: it inlines into the kernel.
      !$acc routine seq
      real(wp), intent(in) :: y1, y2, a, h, dt
      real(wp) :: f

      f = 0.5_wp*dt*h*((a + abs(a))*y1 + (a - abs(a))*y2)
   end function upw

   pure function upwind_table_ok(ws) result(ok)
      !! `.true.` when the tracer table of `ws` is one the kernel carries
      !! (`depend <= 2`).  Host check, cheap: run it before any launch.
      type(transport_workspace_t), intent(in) :: ws
      logical :: ok
      ok = ws%is_init
      if (ok) ok = .not. any(ws%tab%depend(1:ws%ntrace) > 2)
   end function upwind_table_ok

   subroutine upwind_advect(ws, wso, is, ie, js, je, dt, puny, u, v, HTE, HTN, tarea, epm, npm, nbad, nneg, cfg)
      !! Advect the packed fields of `ws` into `wso` over the cells
      !! `is:ie, js:je` (one launch).  `wso` must be initialised with the same
      !! extents and table; it is written only on the range.  Everything must
      !! be device-present on an offload build.  `cfg /= UPWIND_CFG_OK`:
      !! nothing was run.
      type(transport_workspace_t), intent(in) :: ws
      type(transport_workspace_t), intent(inout) :: wso
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(in) :: dt, puny
      real(wp), intent(in) :: u(ws%nx + 1, ws%ny), v(ws%nx, ws%ny + 1)
      real(wp), intent(in) :: HTE(ws%nx, ws%ny), HTN(ws%nx, ws%ny), tarea(ws%nx, ws%ny)
      real(wp), intent(in) :: epm(ws%nx, ws%ny), npm(ws%nx, ws%ny)
      integer, intent(out) :: nbad, nneg, cfg

      nbad = 0
      nneg = 0
      cfg = UPWIND_CFG_OK
      if (.not. (ws%is_init .and. wso%is_init)) then
         cfg = UPWIND_CFG_WORKSPACE
         return
      end if
      if (ws%nx /= wso%nx .or. ws%ny /= wso%ny .or. ws%ncat /= wso%ncat .or. ws%ntrace /= wso%ntrace) then
         cfg = UPWIND_CFG_WORKSPACE
         return
      end if
      if (any(ws%tab%depend(1:ws%ntrace) > 2)) then
         cfg = UPWIND_CFG_DEPEND
         return
      end if
      call upwind_advect_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, is, ie, js, je, dt, puny, &
                              ws%tab%depend(1:ws%ntrace), u, v, HTE, HTN, tarea, epm, npm, &
                              ws%aim, ws%trm, wso%aim, wso%trm, nbad, nneg)
   end subroutine upwind_advect

   pure subroutine upwind_advect_impl(nx, ny, ncat, ntrace, is, ie, js, je, dt, puny, depend, u, v, HTE, HTN, tarea, &
                                      epm, npm, aim, trm, aimo, trmo, nbad, nneg)
      !! Flat kernel behind `upwind_advect` (see the module docstring).
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: dt, puny
      integer, intent(in) :: depend(ntrace)
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: HTE(nx, ny), HTN(nx, ny), tarea(nx, ny), epm(nx, ny), npm(nx, ny)
      real(wp), intent(in) :: aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(inout) :: aimo(nx, ny, 0:ncat), trmo(nx, ny, ntrace, ncat)
      integer, intent(out) :: nbad, nneg
      integer :: i, j, n, nt, d
      real(wp) :: ue, uw, vn, vs, he, hw, hn, hs_, tar
      real(wp) :: ac, ae, aw, an, as_, nfa, anew
      real(wp) :: pc, pe, pw, pn, ps, nf, pnew, bnew, bnf, tnew, nf1, nf2, pn1, pn2

      nbad = 0
      nneg = 0
      !$acc data present(depend, u, v, HTE, HTN, tarea, epm, npm, aim, trm, aimo, trmo)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(ue, uw, vn, vs, he, hw, hn, hs_, tar, ac, ae, aw, an, as_, &
                                                      nfa, anew, nt, d, pc, pe, pw, pn, ps, nf, pnew, bnew, bnf, &
                                                      tnew, nf1, nf2, pn1, pn2) reduce(+:nbad, nneg)
         ue = u(i + 1, j)*epm(i, j)
         uw = u(i, j)*epm(i - 1, j)
         vn = v(i, j + 1)*npm(i, j)
         vs = v(i, j)*npm(i, j - 1)
         he = HTE(i, j)
         hw = HTE(i - 1, j)
         hn = HTN(i, j)
         hs_ = HTN(i, j - 1)
         tar = tarea(i, j)
         ac = aim(i, j, n)
         ae = aim(i + 1, j, n)
         aw = aim(i - 1, j, n)
         an = aim(i, j + 1, n)
         as_ = aim(i, j - 1, n)
         nfa = (upw(ac, ae, ue, he, dt) - upw(aw, ac, uw, hw, dt)) + (upw(ac, an, vn, hn, dt) - upw(as_, ac, vs, hs_, dt))
         anew = ac - nfa/tar
         aimo(i, j, n) = anew
         if (.not. ieee_is_finite(anew)) nbad = nbad + 1
         if (anew < -puny) nneg = nneg + 1
         if (n >= 1) then
            nf1 = 0.0_wp
            nf2 = 0.0_wp
            pn1 = 0.0_wp
            pn2 = 0.0_wp
            !$acc loop seq
            do nt = 1, ntrace
               d = depend(nt)
               if (d == 0) then
                  pc = ac*trm(i, j, nt, n)
                  pe = ae*trm(i + 1, j, nt, n)
                  pw = aw*trm(i - 1, j, nt, n)
                  pn = an*trm(i, j + 1, nt, n)
                  ps = as_*trm(i, j - 1, nt, n)
               else
                  pc = (ac*trm(i, j, d, n))*trm(i, j, nt, n)
                  pe = (ae*trm(i + 1, j, d, n))*trm(i + 1, j, nt, n)
                  pw = (aw*trm(i - 1, j, d, n))*trm(i - 1, j, nt, n)
                  pn = (an*trm(i, j + 1, d, n))*trm(i, j + 1, nt, n)
                  ps = (as_*trm(i, j - 1, d, n))*trm(i, j - 1, nt, n)
               end if
               nf = (upw(pc, pe, ue, he, dt) - upw(pw, pc, uw, hw, dt)) + (upw(pc, pn, vn, hn, dt) - upw(ps, pc, vs, hs_, dt))
               pnew = pc - nf/tar
               if (.not. ieee_is_finite(pnew)) nbad = nbad + 1
               if (d == 0) then
                  bnew = anew
                  bnf = nfa
               else if (d == 1) then
                  bnew = pn1
                  bnf = nf1
               else
                  bnew = pn2
                  bnf = nf2
               end if
               if (nt == 1) then
                  pn1 = pnew
                  nf1 = nf
               else if (nt == 2) then
                  pn2 = pnew
                  nf2 = nf
               end if
               if (nf == 0.0_wp .and. bnf == 0.0_wp) then
                  tnew = trm(i, j, nt, n)
               else if (nt <= 2) then
                  if (bnew > UPWIND_THICK_FLOOR) then
                     tnew = pnew/bnew
                  else
                     tnew = 0.0_wp
                  end if
               else if (d == 0) then
                  if (bnew > puny) then
                     tnew = pnew/bnew
                  else
                     tnew = 0.0_wp
                  end if
               else
                  if (bnew > 0.0_wp) then
                     tnew = pnew/bnew
                  else
                     tnew = 0.0_wp
                  end if
               end if
               trmo(i, j, nt, n) = tnew
            end do
         end if
      end do
      !$acc end data
   end subroutine upwind_advect_impl

   subroutine upwind_empty_tsfc(nx, ny, ncat, is, ie, js, je, puny, Tf, hm, aicen, Tsfcn)
      !! Surface temperature of the categories that hold no ice after the
      !! unpack: `Tsfc = Tf` where `aicen <= puny` (the unpack leaves the
      !! tracers of an empty category untouched), and `Tsfc = 0` on land
      !! (`hm <= 0.5`), over `is:ie, js:je` (one launch).
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: Tf(nx, ny), hm(nx, ny), aicen(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat)

      call upwind_empty_tsfc_impl(nx, ny, ncat, is, ie, js, je, puny, Tf, hm, aicen, Tsfcn)
   end subroutine upwind_empty_tsfc

   pure subroutine upwind_empty_tsfc_impl(nx, ny, ncat, is, ie, js, je, puny, Tf, hm, aicen, Tsfcn)
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: Tf(nx, ny), hm(nx, ny), aicen(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat)
      integer :: i, j, n

      !$acc data present(Tf, hm, aicen, Tsfcn)
      do concurrent(n=1:ncat, j=js:je, i=is:ie)
         if (aicen(i, j, n) <= puny) Tsfcn(i, j, n) = Tf(i, j)
         if (hm(i, j) <= 0.5_wp) Tsfcn(i, j, n) = 0.0_wp
      end do
      !$acc end data
   end subroutine upwind_empty_tsfc_impl

end module rdb_ice_upwind
