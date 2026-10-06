! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_frhat_face.inc — the barotropic face-layer-weight ("frhat") RULE, as an
! inlinable body: given the two columns abutting a u/v face, the per-layer
! face thickness `h_face(k)` that every barotropic depth mean in
! `rdb_barotropic_coupling` (`derive_bt_from_layers`, `face_depth_mean_u/v`,
! `face_depth_mean_rem_u/v`, `apply_bt_correction`'s open/visc_rem folds,
! `compute_bt_rem_from_visc_rem`'s `av_rem`) must use IDENTICALLY, or the
! `dt*F_bt` subtracted back out of the layers is not the quantity the fast
! loop actually integrated (see `face_depth_mean_u`'s own docstring on the
! ISOMIP+ Ocean0 x22.7 bug from exactly one inconsistent call site).
!
! Two schemes (`&ocean_bt_nml frhat_scheme`, `rdb_constants::FRHAT_*`):
!
!   FRHAT_ARITHMETIC (default) — today's plain two-abutting-cell mean,
!     h_face(k) = 0.5*(h_L(k) + h_R(k)).
!
!   FRHAT_HYBRID — MOM6's HVEL_SCHEME=HYBRID (the MOM6 default), ported from
!     `btcalc` (`MOM_barotropic.F90:4546-4790`).  MOM6 is top-down (k=1
!     surface, k=nz bed) and sweeps k=nz downto 1 (bed to surface), building
!     an elevation `e_u` from a shared mean-bed datum upward using the
!     arithmetic-mean thickness of each layer in turn, and comparing it
!     against `D_shallow_u`, the shallower of the two columns' own bed
!     elevations:
!
!       e_u(nz+1) = -0.5*(D_L + D_R)            ! mean-bed datum
!       D_shallow = -min(D_L, D_R)              ! shallower column's bed
!       do k = nz, 1, -1                        ! bed -> surface
!          h_arith = 0.5*(h_L(k) + h_R(k))
!          e_u(k)  = e_u(k+1) + h_arith         ! top of layer k
!          if (e_u(k+1) >= D_shallow) then
!             h_face(k) = h_arith               ! layer wholly above the shelf
!          else
!             h_harm = h_L(k)*h_R(k) / (h_arith + eps)
!             if (e_u(k) <= D_shallow) then
!                h_face(k) = h_harm             ! layer wholly below the shelf
!             else
!                wt = (e_u(k) - D_shallow) / (h_arith + eps)
!                h_face(k) = wt*h_arith + (1-wt)*h_harm   ! straddles it
!             end if
!          end if
!       end do
!
!     roundabout is bottom-up (k=1 bed, k=nz surface), so the SAME sweep runs
!     k = 1, 2, ..., nz (also bed to surface — the physical direction is
!     unchanged, only the label is) and is isomorphic term for term: nothing
!     here depends on the absolute value of k, only on sweep order.  `eps` is
!     `rdb_constants::H_DIV_EPS` (pure 1/0 armour, MOM6's `h_neglect`/
!     `H_subroundoff` role) — NOT a physical floor.
!
!     Net effect: above the shelf (both columns have real water at this
!     layer) the face thickness is the arithmetic mean, same as today; below
!     it, the harmonic mean `h_L*h_R/(h_L+h_R)` VANISHES with the thinner
!     side instead of being dragged up toward the thicker one — the
!     suppression `docs/visc_rem_bt_rem_plan.md` Section 7 identifies as
!     missing.  A face where both columns are identical (a wall, where the
!     caller passes the same column twice) degenerates to `h_arith == h_harm
!     == h_L(k)`, so `FRHAT_HYBRID` is byte-identical to `FRHAT_ARITHMETIC`
!     at a wall for either scheme, matching the pre-port wall fallback.
!
! `frhat_h_face_step` is a per-k STEP of the sweep above: the caller carries
! the running elevation `e_prev` as a `local()` scalar across its own
! `do k = 1, nz` loop (exactly the loop shape every consumer already has) and
! initialises it to `-0.5*(href_L + href_R)` before k=1.  `href_L`/`href_R`
! are the two columns' reference bathymetric depths (roundabout
! `barotropic_workstate_t%bt_H_ref`, MOM6 `G%bathyT`, positive down).  Under
! `FRHAT_ARITHMETIC` `e_prev`/`href_L`/`href_R` are read but never advance the
! sweep state in a way that changes `h_face` — the branch is unconditional.
!
! The consuming module must already `use rdb_constants, only: wp, H_DIV_EPS,
! FRHAT_HYBRID`.  Everything here is `pure` and device-callable
! (`!$acc routine seq`).
!
! Do NOT add state, logging, or anything impure to this file.

   pure subroutine frhat_h_face_step(h_l, h_r, href_l, href_r, scheme, e_prev, h_face)
      !$acc routine seq
      !! One k-step of the MOM6 `btcalc` HYBRID face-thickness sweep (see
      !! this file's header). Call for k = 1, 2, ..., nz IN ORDER with
      !! `e_prev` initialised to `-0.5*(href_l + href_r)` before k=1 — it is
      !! updated in place for the next call. `scheme /= FRHAT_HYBRID` ⇒ the
      !! plain arithmetic mean, `e_prev` left untouched.
      real(wp), intent(in) :: h_l, h_r
         !! This layer's thickness (m) at the left/west (or south) and
         !! right/east (or north) abutting column. Equal at a wall, where
         !! the caller passes the single available column for both.
      real(wp), intent(in) :: href_l, href_r
         !! Reference column depth (m, positive down) at the same two
         !! columns — `barotropic_workstate_t%bt_H_ref`. Unused under
         !! `FRHAT_ARITHMETIC`.
      integer, intent(in) :: scheme
         !! A `FRHAT_*` constant (`rdb_constants`).
      real(wp), intent(inout) :: e_prev
         !! Running elevation (m, mean-bed datum, increasing upward) of the
         !! bottom of the layer about to be processed. Advanced to the top
         !! of this layer on return. Meaningless (and untouched) under
         !! `FRHAT_ARITHMETIC`.
      real(wp), intent(out) :: h_face
         !! This layer's face thickness (m) under the selected scheme.
      real(wp) :: h_arith, h_harm, e_cur, d_shallow, wt_arith

      h_arith = 0.5_wp*(h_l + h_r)
      if (scheme == FRHAT_HYBRID) then
         d_shallow = -min(href_l, href_r)
         e_cur = e_prev + h_arith
         if (e_prev >= d_shallow) then
            ! This layer lies wholly above the shelf formed by the
            ! shallower column's bed: both columns have real water here.
            h_face = h_arith
         else
            h_harm = (h_l*h_r)/(h_arith + H_DIV_EPS)
            if (e_cur <= d_shallow) then
               ! Wholly below the shelf: suppress the thin side instead of
               ! averaging it up toward the thick one.
               h_face = h_harm
            else
               ! Straddles the shelf: linear blend, wt_arith in (0, 1) by
               ! construction (e_prev < d_shallow <= e_cur).
               wt_arith = (e_cur - d_shallow)/(h_arith + H_DIV_EPS)
               h_face = wt_arith*h_arith + (1.0_wp - wt_arith)*h_harm
            end if
         end if
         e_prev = e_cur
      else
         h_face = h_arith
      end if
   end subroutine frhat_h_face_step
