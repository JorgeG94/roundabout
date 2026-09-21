! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_vanished_layer.inc — THE vanished-layer rule, as inlinable bodies.
!
! `#include` this file in the `contains` section of any module that needs to
! know whether a layer is live, what the concentration of a layer is, or that
! a column satisfies invariant I1.  You get a module-LOCAL, `private` copy of
! each procedure, which NVHPC's device codegen CAN inline into a
! `do concurrent` kernel — it does not inline a `pure !$acc routine seq`
! helper across a module boundary (CLAUDE.md, "Cross-TU helper inlining").
! One source of truth, no hand-duplicated `*_impl` twins to drift.
!
! Contract (see `src/core/ocean/README.md`, "The vanished-layer content rule"
! and `src/shared_module_utilities/README.md`):
!
!   I1  `h <= H_VANISHED  ⇒  hTr = 0`, for every registered tracer.
!
! The consuming module must already `use rdb_constants, only: wp, H_VANISHED`
! (and `NZ_STACK_MAX` if it includes the column merge).  Everything here is
! `pure` and device-callable.
!
! Do NOT add state, logging, or anything impure to this file.

   pure logical function rdb_vl_is_live(h)
      !$acc routine seq
      !! THE live-layer predicate.  A layer at or below `H_VANISHED` is a
      !! vanished / inert filler: it carries no meaningful concentration and,
      !! by I1, no content.  The test is a STRICT `>` — a layer sitting
      !! exactly ON the marker reads as vanished, which is what the geometric
      !! vcoord families rely on (`rdb_vcoord :: vcoord_h_min_role`).
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      rdb_vl_is_live = h > H_VANISHED
   end function rdb_vl_is_live

   pure real(wp) function rdb_vl_conc(hTr, h)
      !$acc routine seq
      !! THE concentration of layer `k`: `hTr/h` on a live layer, `0` on a
      !! vanished one.
      !!
      !! Zero, not a floored divide and not a reference value — under I1 a
      !! vanished layer holds no content, so `0` is its honest intensive
      !! value, and it can never be a near-zero divisor feeding a limiter
      !! stencil.  A consumer that needs something else on a vanished layer
      !! (the EOS wants reference T/S so `rho = rho_0`; the sponge wants the
      !! nearest massive layer; the diagnostics want NaN) is making a
      !! DIFFERENT, documented substitution and must not be routed here —
      !! see the survey table in the commit that introduced this file.
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`), e.g. `tracers(t)%hTr(i,j,k)`.
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      if (h > H_VANISHED) then
         rdb_vl_conc = hTr/h
      else
         rdb_vl_conc = 0.0_wp
      end if
   end function rdb_vl_conc

   pure subroutine rdb_vl_merge_content(nz, h_col, q_col)
      !$acc routine seq
      !! Establish I1 on one column: `h <= H_VANISHED ⇒ q = 0`, with the
      !! column sum `Σ q` preserved to round-off.
      !!
      !! `q_col` is CONTENT (`hTr = h·c`, or any extensive per-layer
      !! quantity).  Content held by a sub-threshold layer is handed to the
      !! NEAREST LIVE layer — the next live layer above, or, for a run of
      !! fillers that reaches the top of the column, the topmost live layer
      !! below them — and the filler is zeroed.  A D4 skip/MERGE, not a
      !! clamp: nothing is created, nothing is destroyed, nothing leaves the
      !! column, so a budget accumulator fed `q_new - q_old` still telescopes
      !! to zero per column and stays a valid leak detector.
      !!
      !! **Bit-identity:** every mutation sits behind `orphan /= 0` or the
      !! sub-threshold branch, so a column with no layer at or below
      !! `H_VANISHED` leaves `q_col` textually untouched.
      !!
      !! **Column with no live layer at all** (land / fully grounded): every
      !! layer is a filler, the content is zero by construction, and `q_col`
      !! is zeroed.
      integer, intent(in) :: nz
         !! Number of layers in the column; `k = 1` is the bed.
      real(wp), intent(in) :: h_col(NZ_STACK_MAX)
         !! Layer thicknesses (m).
      real(wp), intent(inout) :: q_col(NZ_STACK_MAX)
         !! Layer content, rewritten in place.
      integer :: k_vl
      real(wp) :: orphan_vl

      orphan_vl = 0.0_wp
      do k_vl = 1, nz
         if (h_col(k_vl) > H_VANISHED) then
            if (orphan_vl /= 0.0_wp) then
               q_col(k_vl) = q_col(k_vl) + orphan_vl
               orphan_vl = 0.0_wp
            end if
         else
            orphan_vl = orphan_vl + q_col(k_vl)
            q_col(k_vl) = 0.0_wp
         end if
      end do
      if (orphan_vl /= 0.0_wp) then
         ! Fillers above the last live layer: hand it back downward.
         do k_vl = nz, 1, -1
            if (h_col(k_vl) > H_VANISHED) then
               q_col(k_vl) = q_col(k_vl) + orphan_vl
               orphan_vl = 0.0_wp
               exit
            end if
         end do
      end if
   end subroutine rdb_vl_merge_content
