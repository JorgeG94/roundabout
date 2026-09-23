! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_vanished_layer.inc — THE vanished-layer rule, as inlinable bodies.
!
! `#include` this file in the `contains` section of any module that needs to
! know whether a layer is live, what the concentration of a layer is, or that
! a column satisfies invariant I1′.  You get a module-LOCAL, `private` copy of
! each procedure, which NVHPC's device codegen CAN inline into a
! `do concurrent` kernel — it does not inline a `pure !$acc routine seq`
! helper across a module boundary (CLAUDE.md, "Cross-TU helper inlining").
! One source of truth, no hand-duplicated `*_impl` twins to drift.
!
! Contract (see `src/core/ocean/README.md`, "The vanished-layer content rule"
! and `src/shared_module_utilities/README.md`):
!
!   I1′  `h <= H_VANISHED  ⇒  hTr = h·c_live`, for every registered tracer,
!        where `c_live` is the concentration of the filler's DONOR: the
!        nearest live layer ABOVE it or, for a run of fillers that reaches
!        the top of the column (the fillers inside an ice shelf), the
!        topmost live layer (`k_top`).  A column with no live layer at all
!        (land, fully grounded) holds `hTr = 0`.
!
! The consuming module must already `use rdb_constants, only: wp, H_VANISHED`
! (and `NZ_STACK_MAX` if it includes the column procedures).  Everything here
! is `pure` and device-callable.
!
! Do NOT add state, logging, or anything impure to this file.

   pure logical function rdb_vl_is_live(h)
      !$acc routine seq
      !! THE live-layer predicate.  A layer at or below `H_VANISHED` is a
      !! vanished / inert filler: it has no concentration of its OWN and, by
      !! I1′, carries its donor's.  The test is a STRICT `>` — a layer
      !! sitting exactly ON the marker reads as vanished, which is what the
      !! geometric vcoord families rely on (`rdb_vcoord :: vcoord_h_min_role`).
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      rdb_vl_is_live = h > H_VANISHED
   end function rdb_vl_is_live

   pure real(wp) function rdb_vl_conc(hTr, h)
      !$acc routine seq
      !! THE concentration of ONE layer, from that layer alone: `hTr/h`.
      !!
      !! On a live layer that is the definition.  On a vanished layer with
      !! `h > 0` it is, by I1′, the donor's concentration `c_live` to
      !! round-off — exactly so on a state the enforcement point
      !! (`multilayer_state_t%enforce_vanished_content`) has just swept.  A
      !! zero-thickness layer has no recoverable concentration and reads
      !! `0`.  A caller that needs `c_live` on a vanished layer ROBUSTLY —
      !! whatever its thickness, and between two sweeps — uses the column
      !! form, `rdb_vl_column_conc`, which reads it off the donor.
      !!
      !! A consumer that wants something ELSE on a vanished layer (the EOS
      !! wants reference T/S so `rho = rho_0`; the diagnostics want NaN) is
      !! making a DIFFERENT, documented substitution and must not be routed
      !! here — see the consumer table in `src/core/ocean/README.md`.
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`), e.g. `tracers(t)%hTr(i,j,k)`.
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      if (h > 0.0_wp) then
         rdb_vl_conc = hTr/h
      else
         rdb_vl_conc = 0.0_wp
      end if
   end function rdb_vl_conc

   pure logical function rdb_vl_holds_live_conc(hTr, h, c_live)
      !$acc routine seq
      !! THE I1′ test on one vanished layer: does it hold `h·c_live` to
      !! round-off?  Relative tolerance `1e-12` on `|h·c_live|`.  The pooled
      !! merge (`rdb_vl_merge_content`) leaves a few ulp times the number of
      !! layers in the pool, so `1e-12` is loose enough for any
      !! `nz <= NZ_STACK_MAX` and tight enough to catch any real violation
      !! (a 0 °C filler under 15 °C water is off by 100 %).  `c_live = 0` (a
      !! column with no live layer, or a genuinely zero tracer) demands
      !! `hTr == 0`.
      real(wp), intent(in) :: hTr
         !! The vanished layer's content.
      real(wp), intent(in) :: h
         !! The vanished layer's thickness (m).
      real(wp), intent(in) :: c_live
         !! Its donor's concentration (`0` for a column with no live layer).
      rdb_vl_holds_live_conc = abs(hTr - h*c_live) <= 1.0e-12_wp*abs(h*c_live)
   end function rdb_vl_holds_live_conc

   pure subroutine rdb_vl_column_conc(nz, h_col, q_col, c_col)
      !$acc routine seq
      !! THE concentration of every layer of a column: `q/h` on a live
      !! layer; on a vanished layer, the concentration of its DONOR (the
      !! nearest live layer above it; the topmost live layer for a run of
      !! fillers that reaches the top of the column) — `c_live` itself, read
      !! off the donor, never recovered from a near-zero divisor.  A column
      !! with no live layer reads `0` throughout.
      !!
      !! **Bit-identity:** on a live layer this is `q/h`, textually the
      !! expression `rdb_vl_conc` evaluates.
      integer, intent(in) :: nz
         !! Number of layers in the column; `k = 1` is the bed.
      real(wp), intent(in) :: h_col(NZ_STACK_MAX)
         !! Layer thicknesses (m).
      real(wp), intent(in) :: q_col(NZ_STACK_MAX)
         !! Layer content.
      real(wp), intent(out) :: c_col(NZ_STACK_MAX)
         !! Layer concentration.
      integer :: k_vl, k_top_vl, k_don_vl

      k_top_vl = 0
      do k_vl = nz, 1, -1
         if (h_col(k_vl) > H_VANISHED) then
            k_top_vl = k_vl
            exit
         end if
      end do
      if (k_top_vl == 0) then
         do k_vl = 1, nz
            c_col(k_vl) = 0.0_wp
         end do
         return
      end if
      ! Walk DOWN: `k_don_vl` is the nearest live layer above `k_vl` — or the
      ! topmost live layer while the walk is still in the run above it.
      k_don_vl = k_top_vl
      do k_vl = nz, 1, -1
         if (h_col(k_vl) > H_VANISHED) then
            c_col(k_vl) = q_col(k_vl)/h_col(k_vl)
            k_don_vl = k_vl
         else
            c_col(k_vl) = q_col(k_don_vl)/h_col(k_don_vl)
         end if
      end do
   end subroutine rdb_vl_column_conc

   pure subroutine rdb_vl_merge_content(nz, h_col, q_col)
      !$acc routine seq
      !! Establish I1′ on one column — `h <= H_VANISHED ⇒ q = h·c_live` —
      !! with the column sum `Σ q` preserved to round-off.
      !!
      !! `q_col` is CONTENT (`hTr = h·c`, or any extensive per-layer
      !! quantity).  Each live layer and the fillers it is DONOR to form one
      !! contiguous POOL: the run of fillers directly below it (down to the
      !! previous live layer, or the bed) plus, for the topmost live layer
      !! only, the run of fillers above it up to the top of the column (the
      !! fillers inside an ice shelf).  The pool is mixed to one
      !! concentration `c = Σq/Σh`: every filler gets `h·c`, and the live
      !! layer gives up exactly what the fillers gained.  So the pool —
      !! hence the column — keeps its content to round-off, and afterwards
      !! the filler
      !! holds exactly the concentration its donor reads
      !! (`rdb_vl_holds_live_conc`).  A D4 skip/MERGE, not a clamp: nothing
      !! is created, nothing destroyed, nothing leaves the column, so a
      !! budget accumulator fed `q_new − q_old` still telescopes to zero per
      !! column and stays a valid leak detector.
      !!
      !! Why a POOL and not "copy the donor's concentration": a filler whose
      !! content has drifted off `h·c_live` (a remap target, a deposit) must
      !! hand the difference back to its donor, which moves the donor's
      !! concentration by a tiny amount — mixing the pool is the one
      !! assignment after which filler and donor agree exactly.  A pool that
      !! already satisfies I1′ (a uniform tracer) mixes back to the same
      !! concentration to round-off, which is what keeps a uniform tracer
      !! uniform.  The previous rule, `h <= H_VANISHED ⇒ q = 0`, did not: the
      !! thickness the continuity step moved OUT of an emptied filler
      !! arrived in the live layers as fresh, 0 °C water.
      !!
      !! **Bit-identity:** every mutation sits behind "this pool contains a
      !! filler", so a column with no layer at or below `H_VANISHED` leaves
      !! `q_col` textually untouched.
      !!
      !! **Column with no live layer at all** (land / fully grounded): there
      !! is no `c_live`, the content is zero by construction (the land-state
      !! contract), and `q_col` is zeroed.
      integer, intent(in) :: nz
         !! Number of layers in the column; `k = 1` is the bed.
      real(wp), intent(in) :: h_col(NZ_STACK_MAX)
         !! Layer thicknesses (m).
      real(wp), intent(inout) :: q_col(NZ_STACK_MAX)
         !! Layer content, rewritten in place.
      integer :: k_vl, kk_vl, k_top_vl, k_lo_vl, k_hi_vl
      real(wp) :: q_sum_vl, h_sum_vl, c_vl, q_new_vl, dq_vl

      k_top_vl = 0
      do k_vl = nz, 1, -1
         if (h_col(k_vl) > H_VANISHED) then
            k_top_vl = k_vl
            exit
         end if
      end do
      if (k_top_vl == 0) then
         do k_vl = 1, nz
            q_col(k_vl) = 0.0_wp
         end do
         return
      end if

      ! `k_lo_vl` = first layer of the filler run waiting below the next
      ! live layer.
      k_lo_vl = 1
      do k_vl = 1, k_top_vl
         if (h_col(k_vl) > H_VANISHED) then
            k_hi_vl = k_vl
            if (k_vl == k_top_vl) k_hi_vl = nz
            if (k_lo_vl < k_vl .or. k_hi_vl > k_vl) then
               ! Pool [k_lo_vl, k_hi_vl]: exactly one live layer, `k_vl`.
               q_sum_vl = 0.0_wp
               h_sum_vl = 0.0_wp
               do kk_vl = k_lo_vl, k_hi_vl
                  q_sum_vl = q_sum_vl + q_col(kk_vl)
                  h_sum_vl = h_sum_vl + h_col(kk_vl)
               end do
               c_vl = q_sum_vl/h_sum_vl
               ! Move only the DIFFERENCE: the donor gives up exactly what
               ! the fillers gain (`dq_vl`), in one subtraction, rather
               ! than being rewritten as `Σq − Σ_fillers h·c` (which
               ! re-rounds the donor through the pool sum on every call).
               ! A pool that already satisfies I1′ then has `dq_vl` at
               ! round-off of the FILLER content, below one ulp of the
               ! donor, and the donor is left bit-for-bit untouched.
               dq_vl = 0.0_wp
               do kk_vl = k_lo_vl, k_hi_vl
                  if (kk_vl /= k_vl) then
                     q_new_vl = h_col(kk_vl)*c_vl
                     dq_vl = dq_vl + (q_new_vl - q_col(kk_vl))
                     q_col(kk_vl) = q_new_vl
                  end if
               end do
               q_col(k_vl) = q_col(k_vl) - dq_vl
            end if
            k_lo_vl = k_vl + 1
         end if
      end do
   end subroutine rdb_vl_merge_content
