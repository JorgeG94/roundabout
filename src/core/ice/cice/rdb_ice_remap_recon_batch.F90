!! Sea-ice incremental remapping, step 1 for ONE CATEGORY BATCH: the limited linear
!! reconstruction of the area and of the tracer hierarchy (Dukowicz and Baumgardner 2000,
!! J. Comput. Phys. 160, 318; Lipscomb and Hunke 2004, Mon. Wea. Rev. 132, 1341) without
!! stored masks.
module rdb_ice_remap_recon_batch
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
   !
   !! **What this is.**  The same reconstruction as `rdb_ice_remap_recon` (read that header
   !! for the profile, the limiter and the identities), reorganised for memory: (1) the masks
   !! are not stored: `mmask = aim > puny`, `tmask(parent) = aim > puny .and. |trm(parent)| > puny`
   !! (only parents carry one) and `mc = aim where aim > puny, else 0` are re-evaluated from the packed fields
   !! wherever a neighbour's mask or a cell's centre value is read, so the masks launch is fused
   !! into the fields launch; (2) only the categories of ONE BATCH are reconstructed, into
   !! buffers `mx, my(nx, ny, nb)` and `tc, tx, ty(nx, ny, ntrace, nb)` with the batch slot
   !! `b = n - nlo + 1` (`n = 0` is open water: mass gradients only).  The values are the
   !! stored-mask module's, bit for bit.  The centre value `mc` is not stored: the flux
   !! integral reads `aim` (see `rdb_ice_remap_flux`).
   !!
   !! **Halo.**  Outputs on the cells `is:ie, js:je` (the caller gives owned `+- 1`), inputs `aim`,
   !! `trm`, `hm` valid one cell wider; `hm` exactly 0 or 1.  `remap_input_check` counts, over the
   !! whole read ring and every category, non-finite packed values and `hm` values that are
   !! neither 0 nor 1 (the caller turns a non-finite count into `ICE_ERR_NONFINITE`).
   !!
   !! **GPU shape.**  ONE `do concurrent (b, j, i)` launch per batch, the tracer loop sequential
   !! inside the thread in hierarchy order, the nine stencil values passed as scalars to
   !! `recon_limited_gradient`.  `remap_input_check` is one launch per step, two `reduce(+)`
   !! accumulators.
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_remap_recon, only: recon_limited_gradient, recon_table_ok, RECON_XXAV, RECON_CFG_OK, RECON_CFG_TABLE, &
                                  RECON_CFG_WORKSPACE, RECON_CFG_RANGE, RECON_CFG_MASK
   use rdb_ice_transport_pack, only: transport_workspace_t
   implicit none
   private

   public :: remap_recon_batch_t, remap_recon_batch, remap_input_check, remap_recon_batch_impl, remap_input_check_impl
   public :: RECON_CFG_OK, RECON_CFG_TABLE, RECON_CFG_WORKSPACE, RECON_CFG_RANGE, RECON_CFG_MASK

   type :: remap_recon_batch_t
      !! Persistent reconstruction buffers of one category batch.  `init` allocates (zeroed),
      !! `enter_data` maps the component arrays ONCE, `exit_data` unmaps, `destroy` frees.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      integer :: ncat = 0
      integer :: ntrace = 0
      integer :: nbatch = 0
      real(wp), allocatable :: mx(:, :, :), my(:, :, :)
         !! `(nx, ny, nbatch)` mass gradients per cell width.
      real(wp), allocatable :: tc(:, :, :, :), tx(:, :, :, :), ty(:, :, :, :)
         !! `(nx, ny, ntrace, nbatch)` tracer centre value and gradients.
   contains
      procedure :: init => rb_init
      procedure :: enter_data => rb_enter_data
      procedure :: exit_data => rb_exit_data
      procedure :: destroy => rb_destroy
      procedure :: bytes => rb_bytes
   end type remap_recon_batch_t

contains

   subroutine rb_init(this, nx, ny, ncat, ntrace, nbatch, ierr)
      !! Allocate (zeroed).  `ierr /= 0`: invalid extents or batch, nothing allocated.
      class(remap_recon_batch_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, ntrace, nbatch
      integer, intent(out) :: ierr

      ierr = 0
      if (nx < 3 .or. ny < 3 .or. ncat < 1 .or. ntrace < 1 .or. nbatch < 1 .or. nbatch > ncat + 1) then
         ierr = RECON_CFG_WORKSPACE
         return
      end if
      call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%ntrace = ntrace
      this%nbatch = nbatch
      allocate (this%mx(nx, ny, nbatch), this%my(nx, ny, nbatch))
      allocate (this%tc(nx, ny, ntrace, nbatch), this%tx(nx, ny, ntrace, nbatch), this%ty(nx, ny, ntrace, nbatch))
      this%mx = 0.0_wp
      this%my = 0.0_wp
      this%tc = 0.0_wp
      this%tx = 0.0_wp
      this%ty = 0.0_wp
      this%is_init = .true.
   end subroutine rb_init

   subroutine rb_enter_data(this)
      class(remap_recon_batch_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%mx, this%my)
      !$acc enter data copyin(this%tc, this%tx, this%ty)
      this%mapped = .true.
   end subroutine rb_enter_data

   subroutine rb_exit_data(this)
      class(remap_recon_batch_t), intent(inout) :: this
      if (.not. this%mapped) return
      !$acc exit data delete(this%mx, this%my)
      !$acc exit data delete(this%tc, this%tx, this%ty)
      this%mapped = .false.
   end subroutine rb_exit_data

   subroutine rb_destroy(this)
      class(remap_recon_batch_t), intent(inout) :: this
      call this%exit_data()
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
      this%nbatch = 0
   end subroutine rb_destroy

   pure function rb_bytes(this) result(b)
      !! Bytes of the buffers (host copy; the device mirror is the same size).
      class(remap_recon_batch_t), intent(in) :: this
      integer(int64) :: b
      b = 0_int64
      if (allocated(this%mx)) b = b + int(size(this%mx), int64)*8_int64
      if (allocated(this%my)) b = b + int(size(this%my), int64)*8_int64
      if (allocated(this%tc)) b = b + int(size(this%tc), int64)*8_int64
      if (allocated(this%tx)) b = b + int(size(this%tx), int64)*8_int64
      if (allocated(this%ty)) b = b + int(size(this%ty), int64)*8_int64
   end function rb_bytes

   subroutine remap_input_check(ws, is, ie, js, je, hm, nbad, cfg)
      !! Count, over the read ring `is-2:ie+2, js-2:je+2` and every category, the non-finite packed
      !! values and `hm` (`nbad`), and the finite `hm` values other than 0 and 1 (`cfg =
      !! RECON_CFG_MASK`).  Also checks the table and the range.  One launch.
      type(transport_workspace_t), intent(in) :: ws
      integer, intent(in) :: is, ie, js, je
      real(wp), intent(in) :: hm(ws%nx, ws%ny)
      integer, intent(out) :: nbad, cfg
      integer :: nhm

      nbad = 0
      cfg = RECON_CFG_OK
      if (.not. ws%is_init) then
         cfg = RECON_CFG_WORKSPACE
         return
      end if
      if (.not. recon_table_ok(ws)) then
         cfg = RECON_CFG_TABLE
         return
      end if
      if (is < 3 .or. ie > ws%nx - 2 .or. js < 3 .or. je > ws%ny - 2 .or. is > ie .or. js > je) then
         cfg = RECON_CFG_RANGE
         return
      end if
      call remap_input_check_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, is - 2, ie + 2, js - 2, je + 2, hm, ws%aim, ws%trm, nbad, nhm)
      if (nhm > 0) cfg = RECON_CFG_MASK
   end subroutine remap_input_check

   subroutine remap_input_check_impl(nx, ny, ncat, ntrace, is, ie, js, je, hm, aim, trm, nbad, nhm)
      integer, intent(in) :: nx, ny, ncat, ntrace, is, ie, js, je
      real(wp), intent(in) :: hm(nx, ny), aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      integer, intent(out) :: nbad, nhm
      integer :: i, j, n, nt

      nbad = 0
      nhm = 0
      !$acc data present(hm, aim, trm)
      do concurrent(n=0:ncat, j=js:je, i=is:ie) local(nt) reduce(+:nbad, nhm)
         if (.not. ieee_is_finite(aim(i, j, n))) nbad = nbad + 1
         if (n == 0) then
            if (.not. ieee_is_finite(hm(i, j))) then
               nbad = nbad + 1
            else if (hm(i, j) /= 0.0_wp .and. hm(i, j) /= 1.0_wp) then
               nhm = nhm + 1
            end if
         else
            !$acc loop seq
            do nt = 1, ntrace
               if (.not. ieee_is_finite(trm(i, j, nt, n))) nbad = nbad + 1
            end do
         end if
      end do
      !$acc end data
   end subroutine remap_input_check_impl

   subroutine remap_recon_batch(ws, rb, is, ie, js, je, nlo, nhi, puny, hm, cfg)
      !! Reconstruct the categories `nlo..nhi` (`nhi - nlo + 1 <= rb%nbatch`) of the packed fields
      !! of `ws` into `rb` over `is:ie, js:je` (one launch).  The range and the inputs are checked
      !! by `remap_input_check` (once per step).  Everything device-present.
      type(transport_workspace_t), intent(in) :: ws
      type(remap_recon_batch_t), intent(inout) :: rb
      integer, intent(in) :: is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: hm(ws%nx, ws%ny)
      integer, intent(out) :: cfg

      cfg = RECON_CFG_OK
      if (.not. (ws%is_init .and. rb%is_init)) then
         cfg = RECON_CFG_WORKSPACE
         return
      end if
      if (ws%nx /= rb%nx .or. ws%ny /= rb%ny .or. ws%ncat /= rb%ncat .or. ws%ntrace /= rb%ntrace &
          .or. nhi - nlo + 1 > rb%nbatch .or. nlo < 0 .or. nhi > ws%ncat) then
         cfg = RECON_CFG_WORKSPACE
         return
      end if
      call remap_recon_batch_impl(ws%nx, ws%ny, ws%ncat, ws%ntrace, rb%nbatch, is, ie, js, je, nlo, nhi, puny, &
                                  ws%tab%tracer_type(1:ws%ntrace), ws%tab%depend(1:ws%ntrace), hm, ws%aim, ws%trm, &
                                  rb%mx, rb%my, rb%tc, rb%tx, rb%ty)
   end subroutine remap_recon_batch

   subroutine remap_recon_batch_impl(nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi, puny, tracer_type, depend, hm, &
                                     aim, trm, mx, my, tc, tx, ty)
      !! Flat kernel (see the module docstring).
      integer, intent(in) :: nx, ny, ncat, ntrace, nb, is, ie, js, je, nlo, nhi
      real(wp), intent(in) :: puny
      integer, intent(in) :: tracer_type(ntrace), depend(ntrace)
      real(wp), intent(in) :: hm(nx, ny), aim(nx, ny, 0:ncat), trm(nx, ny, ntrace, ncat)
      real(wp), intent(inout) :: mx(nx, ny, nb), my(nx, ny, nb)
      real(wp), intent(inout) :: tc(nx, ny, ntrace, nb), tx(nx, ny, ntrace, nb), ty(nx, ny, ntrace, nb)
      integer :: i, j, b, n, nt, p
      real(wp) :: mm, mcv, gx, gy, mxav, myav, tm, tgx, tgy, tcv, mtx, mty, w2, w3, w7
      real(wp) :: pc, pnw, pn, pne, pw, pe, psw, ps, pse
      logical :: tpar

      !$acc data present(tracer_type, depend, hm, aim, trm, mx, my, tc, tx, ty)
      do concurrent(b=1:nhi - nlo + 1, j=js:je, i=is:ie) local(n, mm, mcv, gx, gy, mxav, myav, tm, tgx, tgy, tcv, mtx, mty, &
                                                               w2, w3, w7, pc, pnw, pn, pne, pw, pe, psw, ps, pse, nt, p, tpar)
         n = nlo + b - 1
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
         mx(i, j, b) = gx
         my(i, j, b) = gy
         mcv = 0.0_wp
         if (mm > puny) mcv = mm
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
                     ! neighbours with area above puny count (the stored mask `mmask`)
                     pc = tm
                     pnw = merge(trm(i - 1, j + 1, nt, n), pc, aim(i - 1, j + 1, n) > puny)
                     pn = merge(trm(i, j + 1, nt, n), pc, aim(i, j + 1, n) > puny)
                     pne = merge(trm(i + 1, j + 1, nt, n), pc, aim(i + 1, j + 1, n) > puny)
                     pw = merge(trm(i - 1, j, nt, n), pc, aim(i - 1, j, n) > puny)
                     pe = merge(trm(i + 1, j, nt, n), pc, aim(i + 1, j, n) > puny)
                     psw = merge(trm(i - 1, j - 1, nt, n), pc, aim(i - 1, j - 1, n) > puny)
                     ps = merge(trm(i, j - 1, nt, n), pc, aim(i, j - 1, n) > puny)
                     pse = merge(trm(i + 1, j - 1, nt, n), pc, aim(i + 1, j - 1, n) > puny)
                     call recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, mxav, myav, tgx, tgy)
                     tcv = tm - tgx*mxav - tgy*myav
                  else
                     ! type 2: the parent's centroid of (area x tracer) in this cell
                     p = depend(nt)
                     mtx = 0.0_wp
                     mty = 0.0_wp
                     tpar = abs(trm(i, j, p, n)) > puny
                     if (tpar) then
                        w2 = mcv*tx(i, j, p, b) + gx*tc(i, j, p, b)
                        w3 = mcv*ty(i, j, p, b) + gy*tc(i, j, p, b)
                        w7 = 1.0_wp/(mm*trm(i, j, p, n))
                        mtx = w2*RECON_XXAV*w7
                        mty = w3*RECON_XXAV*w7
                        pc = tm
                        ! neighbours whose parent tracer is present (the stored mask `tmask(parent)`)
                    pnw = merge(trm(i - 1, j + 1, nt, n), pc, aim(i - 1, j + 1, n) > puny .and. abs(trm(i - 1, j + 1, p, n)) > puny)
                        pn = merge(trm(i, j + 1, nt, n), pc, aim(i, j + 1, n) > puny .and. abs(trm(i, j + 1, p, n)) > puny)
                    pne = merge(trm(i + 1, j + 1, nt, n), pc, aim(i + 1, j + 1, n) > puny .and. abs(trm(i + 1, j + 1, p, n)) > puny)
                        pw = merge(trm(i - 1, j, nt, n), pc, aim(i - 1, j, n) > puny .and. abs(trm(i - 1, j, p, n)) > puny)
                        pe = merge(trm(i + 1, j, nt, n), pc, aim(i + 1, j, n) > puny .and. abs(trm(i + 1, j, p, n)) > puny)
                    psw = merge(trm(i - 1, j - 1, nt, n), pc, aim(i - 1, j - 1, n) > puny .and. abs(trm(i - 1, j - 1, p, n)) > puny)
                        ps = merge(trm(i, j - 1, nt, n), pc, aim(i, j - 1, n) > puny .and. abs(trm(i, j - 1, p, n)) > puny)
                    pse = merge(trm(i + 1, j - 1, nt, n), pc, aim(i + 1, j - 1, n) > puny .and. abs(trm(i + 1, j - 1, p, n)) > puny)
                        call recon_limited_gradient(pc, pnw, pn, pne, pw, pe, psw, ps, pse, mtx, mty, tgx, tgy)
                     end if
                     tcv = tm - tgx*mtx - tgy*mty
                  end if
               end if
               tx(i, j, nt, b) = tgx
               ty(i, j, nt, b) = tgy
               tc(i, j, nt, b) = tcv
            end do
         end if
      end do
      !$acc end data
   end subroutine remap_recon_batch_impl

end module rdb_ice_remap_recon_batch
