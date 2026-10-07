! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_bed_live_face.inc — the "is this layer live at this face" RULE used to
! find `kb_live`, the first LIVE layer at/above a face's static bed index
! (`ms%k_bot_u/v`, filled once at configure from the `eta=0` target — see
! "Vertical Layer Convention" in CLAUDE.md).  Under a `z_fixed`-family
! vcoord a free-surface drawdown can vanish the nominal bed sliver `kb`
! points at (its layer thickness falls to/below `H_VANISHED` on one or both
! sides of the face) WITHOUT moving `kb` itself, which is static for the
! run.  Everything that applies bed drag or BBL physics to "the bed row" —
! `rdb_ocean_bottom_drag`'s explicit/implicit-fold drag and its HBBL band
! walk, and `rdb_ocean_vdiff`'s bed-BC row placement in
! `diffuse_velocity_columns_impl` — must walk up from `kb` to the same
! first-live `kb_live`, or the drag law samples one row and the solver
! glues the bed BC to another (#178: the two silently disagreed, the drag
! piston decoupled from the live column, and a bed jet ran away with
! `visc_rem == 1` exactly).
!
! Two regimes, matching `diffuse_velocity_columns_impl`'s own `hvel(k)`
! construction for the ONE case every consumer must agree on bit for bit:
!
!   zlevel_faces (`&vcoord_nml zfixed_closed_faces`) — the face column is
!     the OVERLAP of its two abutting cell columns, so a layer vanished on
!     EITHER side decouples the whole face (MIN of the two thicknesses);
!     `diffuse_velocity_columns_impl` applies this unconditionally,
!     overriding whatever blend built `hvel(k)` beforehand.
!
!   otherwise — every row is live, so `kb_live = kb`.  Off closed faces
!     nothing decouples the static bed row from the live column (the
!     vdiff coupling to a thin row is strong, not zeroed), so #178 cannot
!     occur and the bed row stays where it was: bit-identical for every
!     non-`zfixed_closed_faces` configuration, including `zstar_full`'s
!     dynamic bed pinch (whose HBBL band walk is unchanged).
!
! Needs from the consuming module: `wp` (`rdb_constants`) AND the vanished-
! layer rule itself, `#include "rdb_vanished_layer.inc"` — the liveness
! VERDICT is `rdb_vl_is_live`, not a second hand-rolled `> H_VANISHED`
! (the `vanished-layer` pre-commit hook enforces this tree-wide).
pure function rdb_blf_is_live(hl, hr, zlevel_faces) result(is_live)
   !$acc routine seq
   real(wp), intent(in) :: hl, hr
      !! The two abutting cells' layer thickness at this face/layer.
   logical, intent(in) :: zlevel_faces
   logical :: is_live

   if (zlevel_faces) then
      is_live = rdb_vl_is_live(min(hl, hr))
   else
      ! Off closed faces nothing decouples the static bed row from the
      ! live column, so it stays the bed row: `kb_live = kb` ⇒
      ! bit-identical everywhere #178 cannot occur.  `hl`/`hr` are
      ! deliberately unread here.
      is_live = .true.
   end if
end function rdb_blf_is_live
