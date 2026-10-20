!! Ice-thickness-distribution (ITD) restore: adjust_ice_categories (PR 4a).
module rdb_ice_itd
   !! Port of SIS2 `adjust_ice_categories` (`SIS_transport.F90:611-891`,
   !! Apache-2.0). SIS2's shipped algorithm is NOT the Lipscomb (2001)
   !! linear-profile remap — that is a `SIS_transport.F90:737` TODO
   !! comment, never code. What SIS2 actually runs (and what this module
   !! ports) is the *"for now move all of it"* WHOLE-CATEGORY shift: when
   !! a category's mean ice mass crosses a bin boundary, its ENTIRE
   !! area+mass moves to the adjacent category, merged into that
   !! category's existing area+mass via an area-weighted destination-
   !! thickness average (mass AND area conserving, SIS2:730-741) with
   !! upwind mass-weighted tracer mixing (equivalent in effect to SIS2's
   !! deferred `advect_tracers_thicker`, `SIS_tracer_advect.F90:1791-1850`
   !! — each category boundary is visited once per pass in a fixed order
   !! with running masses, so a single closed-form mass-weighted mix at
   !! the transfer site reproduces the same result; verified by the
   !! validated Python prototype, `tmp_local_artifacts/proto_itd_adjust.py`).
   !!
   !! Ported with `ice_cover_discard = -1` (off, SIS2 default),
   !! `ocean_part_min = 0`, `inconsistent_cover_bug = .false.` (i.e. DO
   !! the open-water resum).
   !!
   !! Algorithm per cell (all steps operate on the SAME category array —
   !! order matters, see the docstrings of each impl step below):
   !!   1. Massless cleanup (SIS2:664-685): any category with
   !!      `m_ice <= 0` has its area returned to open water (flagged for
   !!      the final resum, not applied immediately).
   !!   2. Upward pass, `c = 1..ncat-1` ascending (SIS2:714-767): a
   !!      category whose ice mass exceeds ITS OWN upper bin edge
   !!      (`mh_lim(c+1)`) moves ALL its area+mass into category `c+1`.
   !!   3. Downward pass, `c = ncat..2` descending (SIS2:786-839): mirror
   !!      of step 2 — a category whose ice mass falls below its OWN
   !!      lower bin edge (`mh_lim(c)`) moves ALL its area+mass into
   !!      category `c-1`.
   !!   4. Cat-1 minimum-thickness compress (SIS2:850-864): category 1 has
   !!      no category below it to demote into (`mh_lim(1) > 0` is the
   !!      physical floor), so instead of a transfer it AREA-SHRINKS:
   !!      `part(1) *= m_ice(1)/mh_lim(1)`, pinning `m_ice(1)` at the
   !!      floor (mass- and per-remaining-area-tracer-conserving; the
   !!      lost area returns to open water via the resum).
   !!   5. Open-water resum (SIS2:874-889): only if step 1 or step 4
   !!      fired (the up/down transfers of steps 2-3 preserve `Σ part`
   !!      exactly, add/subtract the same `part_trans` on both sides, so
   !!      they never need a resum) — `part(0) = max(1 - Σ_{c=1..ncat}
   !!      part(c), 0)`.
   !!
   !! Deliberately NOT ported (SIS2 features orthogonal to the "for now
   !! move all of it" shift, or infra this branch does not have yet):
   !! melt ponds (no pond state on `ocean_sea_ice_t`), `t_surf` category
   !! remap (`tsurf_out` here is per-thermo-window scratch, recomputed
   !! from scratch every step — nothing to remap), `ice_cover_discard`
   !! (SIS2 default off), the FATAL input-consistency checks (a GPU
   !! `do concurrent` kernel cannot abort mid-loop — negative mass /
   !! snow-on-no-ice is treated as "massless", the same taxonomy as step
   !! 1), and SIS2's `do_j` per-row early-exit optimisation (the per-cell
   !! `if` gate below is the DC-kernel equivalent). Ridging is a SEPARATE,
   !! still-unported gap — NOT closed by `compress_ice` (PR 4b,
   !! `rdb_ice_transport`): `compress_ice` is SIS2's own
   !! `DO_RIDGING=.false.` fallback (area compaction, in-category, zero
   !! energetic cost), not the participation/redistribution ridging
   !! scheme SIS2 ships as an Icepack wrapper (`ice_ridge.F90`, default
   !! off). See `docs/CAPABILITIES_AND_LIMITATIONS.md` "Sea ice" for the
   !! physical consequence.
   !!
   !! **ncat==1 short-circuit (bit-identity contract, module docstring of
   !! `rdb_ice_state`)**: the outer shim returns immediately when
   !! `ice%ncat == 1` — the legacy lumped mode never runs any part of this
   !! module, by construction (not by an exactness argument). SIS2 itself
   !! would still run the cat-1 compress at nCat=1 (`mh_lim(1) = 1e-10 m`
   !! makes it physically unreachable in practice), but skipping it
   !! outright is the simplest statement of "ncat=1 is untouched" and is
   !! the explicit contract this branch commits to.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private

   public :: ice_adjust_categories

contains

   pure subroutine ice_adjust_categories(grid, ms, ice)
      !! Outer shim (outer-shim + flat-impl pattern). No-op when the ice
      !! slot is not live (`is_init` gate) or when `ice%ncat == 1` (the
      !! legacy-mode bit-identity contract — see module docstring). No
      !! `dt`/`eos` needed: this is a pure area/mass reshuffle, no
      !! thermodynamics.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
         !! READ-ONLY: only `wet_mask` is read (the gate).
      type(ocean_sea_ice_t), intent(inout) :: ice

      if (.not. ice%is_init) return
      if (ice%ncat == 1) return

      call ice_adjust_categories_impl(ms%wet_mask, ice%part_size, ice%m_ice, ice%m_snow, &
                                      ice%enth_ice, ice%enth_snow, ice%sal_ice, ice%mh_lim, &
                                      grid%nghost, ice%ncat, ice%nk_ice, &
                                      grid%nx_total, grid%ny_total)
   end subroutine ice_adjust_categories

   pure subroutine ice_adjust_categories_impl(wet_mask, part_size, m_ice, m_snow, &
                                              enth_ice, enth_snow, sal_ice, mh_lim, &
                                              nghost, ncat, nk, nx, ny)
      !! Device kernel: ONE `do concurrent(j, i)` over PHYSICAL cells
      !! (ghosts excluded), `wet_mask > 0.5` inner gate (never a masked DC
      !! header), serial `cat` loops inside — each cell touches only its
      !! own `(i, j, :)` slice, race-free. No per-cell gather into local
      !! category-sized arrays (avoids register pressure and any
      !! compile-time `ncat` cap, unlike the `ICE_NK_MAX`-capped
      !! PER-LAYER arrays elsewhere in the ice model) — scalar
      !! temporaries only, all declared in `local(...)`.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. No locals named after intrinsics.
      integer, intent(in) :: nghost, ncat, nk, nx, ny
         !! Grid + category + layer extents (declared first — decl-order).
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: enth_snow(nx, ny, ncat, 1)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
      real(wp), intent(in) :: mh_lim(ncat + 1)

      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      integer :: cat, lay
      logical :: resum
      real(wp) :: part_trans, part_sum
      real(wp) :: mca_ice_src, mca_ice_dst, mca_snow_src, mca_snow_dst, mnew, msnew

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(cat, lay, resum, part_trans, part_sum, &
                                                    mca_ice_src, mca_ice_dst, mca_snow_src, &
                                                    mca_snow_dst, mnew, msnew)
         if (wet_mask(i, j) > 0.5_wp) then
            resum = .false.

            ! ---- 1. Massless cleanup (SIS2:664-685) ----
            do cat = 1, ncat
               if (m_ice(i, j, cat) <= 0.0_wp) then
                  if (part_size(i, j, cat) > 0.0_wp) resum = .true.
                  part_size(i, j, cat) = 0.0_wp
               end if
            end do

            ! ---- 2. Upward pass, c = 1..ncat-1 ascending (SIS2:714-767) ----
            do cat = 1, ncat - 1
               if (part_size(i, j, cat)*m_ice(i, j, cat) > 0.0_wp .and. &
                   m_ice(i, j, cat) > mh_lim(cat + 1)) then
                  part_trans = part_size(i, j, cat)
                  mca_snow_src = part_size(i, j, cat)*m_snow(i, j, cat)
                  mca_snow_dst = part_size(i, j, cat + 1)*m_snow(i, j, cat + 1)
                  mca_ice_src = part_trans*m_ice(i, j, cat)
                  mca_ice_dst = part_size(i, j, cat + 1)*m_ice(i, j, cat + 1)

                  m_ice(i, j, cat + 1) = (part_trans*m_ice(i, j, cat) &
                                          + part_size(i, j, cat + 1)*m_ice(i, j, cat + 1)) &
                                         /(part_trans + part_size(i, j, cat + 1))
                  m_ice(i, j, cat) = mh_lim(cat + 1)
                  part_size(i, j, cat + 1) = part_size(i, j, cat + 1) + part_trans
                  part_size(i, j, cat) = part_size(i, j, cat) - part_trans

                  m_snow(i, j, cat) = 0.0_wp
                  if (part_size(i, j, cat + 1) > 0.0_wp) then
                     m_snow(i, j, cat + 1) = (mca_snow_src + mca_snow_dst)/part_size(i, j, cat + 1)
                  else
                     m_snow(i, j, cat + 1) = 0.0_wp
                  end if

                  mnew = mca_ice_src + mca_ice_dst
                  if (mca_ice_src > 0.0_wp .and. mnew > 0.0_wp) then
                     do lay = 1, nk
                        enth_ice(i, j, cat + 1, lay) = (mca_ice_src*enth_ice(i, j, cat, lay) &
                                                        + mca_ice_dst*enth_ice(i, j, cat + 1, lay)) &
                                                       /mnew
                        sal_ice(i, j, cat + 1, lay) = (mca_ice_src*sal_ice(i, j, cat, lay) &
                                                       + mca_ice_dst*sal_ice(i, j, cat + 1, lay)) &
                                                      /mnew
                     end do
                  end if
                  msnew = mca_snow_src + mca_snow_dst
                  if (mca_snow_src > 0.0_wp .and. msnew > 0.0_wp) then
                     enth_snow(i, j, cat + 1, 1) = (mca_snow_src*enth_snow(i, j, cat, 1) &
                                                    + mca_snow_dst*enth_snow(i, j, cat + 1, 1)) &
                                                   /msnew
                  end if
               end if
            end do

            ! ---- 3. Downward pass, c = ncat..2 descending (SIS2:786-839) ----
            do cat = ncat, 2, -1
               if (part_size(i, j, cat)*m_ice(i, j, cat) > 0.0_wp .and. &
                   m_ice(i, j, cat) < mh_lim(cat)) then
                  part_trans = part_size(i, j, cat)
                  mca_snow_src = part_size(i, j, cat)*m_snow(i, j, cat)
                  mca_snow_dst = part_size(i, j, cat - 1)*m_snow(i, j, cat - 1)
                  mca_ice_src = part_trans*m_ice(i, j, cat)
                  mca_ice_dst = part_size(i, j, cat - 1)*m_ice(i, j, cat - 1)

                  m_ice(i, j, cat - 1) = (part_trans*m_ice(i, j, cat) &
                                          + part_size(i, j, cat - 1)*m_ice(i, j, cat - 1)) &
                                         /(part_trans + part_size(i, j, cat - 1))
                  m_ice(i, j, cat) = mh_lim(cat)
                  part_size(i, j, cat - 1) = part_size(i, j, cat - 1) + part_trans
                  part_size(i, j, cat) = part_size(i, j, cat) - part_trans

                  m_snow(i, j, cat) = 0.0_wp
                  if (part_size(i, j, cat - 1) > 0.0_wp) then
                     m_snow(i, j, cat - 1) = (mca_snow_src + mca_snow_dst)/part_size(i, j, cat - 1)
                  else
                     m_snow(i, j, cat - 1) = 0.0_wp
                  end if

                  mnew = mca_ice_src + mca_ice_dst
                  if (mca_ice_src > 0.0_wp .and. mnew > 0.0_wp) then
                     do lay = 1, nk
                        enth_ice(i, j, cat - 1, lay) = (mca_ice_src*enth_ice(i, j, cat, lay) &
                                                        + mca_ice_dst*enth_ice(i, j, cat - 1, lay)) &
                                                       /mnew
                        sal_ice(i, j, cat - 1, lay) = (mca_ice_src*sal_ice(i, j, cat, lay) &
                                                       + mca_ice_dst*sal_ice(i, j, cat - 1, lay)) &
                                                      /mnew
                     end do
                  end if
                  msnew = mca_snow_src + mca_snow_dst
                  if (mca_snow_src > 0.0_wp .and. msnew > 0.0_wp) then
                     enth_snow(i, j, cat - 1, 1) = (mca_snow_src*enth_snow(i, j, cat, 1) &
                                                    + mca_snow_dst*enth_snow(i, j, cat - 1, 1)) &
                                                   /msnew
                  end if
               end if
            end do

            ! ---- 4. Cat-1 minimum-thickness compress (SIS2:850-864) ----
            if (mh_lim(1) > 0.0_wp) then
               if (m_ice(i, j, 1)*part_size(i, j, 1) > 0.0_wp .and. m_ice(i, j, 1) < mh_lim(1)) then
                  part_size(i, j, 1) = part_size(i, j, 1)*(m_ice(i, j, 1)/mh_lim(1))
                  if (m_ice(i, j, 1) > 0.0_wp) then
                     m_snow(i, j, 1) = m_snow(i, j, 1)*(mh_lim(1)/m_ice(i, j, 1))
                  end if
                  m_ice(i, j, 1) = mh_lim(1)
                  resum = .true.
               end if
            end if

            ! ---- 5. Open-water resum (SIS2:874-889) ----
            if (resum) then
               part_sum = 0.0_wp
               do cat = 1, ncat
                  part_sum = part_sum + part_size(i, j, cat)
               end do
               part_size(i, j, 0) = max(1.0_wp - part_sum, 0.0_wp)
            end if
         end if
      end do
   end subroutine ice_adjust_categories_impl

end module rdb_ice_itd
