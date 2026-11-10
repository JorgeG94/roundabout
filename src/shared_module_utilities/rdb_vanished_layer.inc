! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_vanished_layer.inc — THE vanished-layer predicate and donor read, as
! inlinable bodies.
!
! `#include` this file in the `contains` section of any module that needs to
! know whether a layer is live, or what a vanished layer should read as when a
! consumer integrates its concentration over a depth the layer does not
! occupy.  You get a module-LOCAL, `private` copy of each procedure, which
! NVHPC's device codegen CAN inline into a `do concurrent` kernel — it does
! not inline a `pure !$acc routine seq` helper across a module boundary
! (CLAUDE.md, "Cross-TU helper inlining").  One source of truth, no
! hand-duplicated `*_impl` twins to drift.
!
! Contract (see also `src/shared_module_utilities/README.md`):
!
!   Tracers are stored as CONCENTRATION `c`; the conserved content is
!   `h_layer·c`, formed only where a budget or console total needs it.
!   Several vertical coordinates (`z_fixed`, `zstar`, `zstar_full`, wet/dry)
!   place inert FILLER layers with `h <= H_VANISHED`.  A filler keeps the
!   last concentration it held: transport weights it by its own thickness, so
!   it carries no content of its own to maintain, and nothing is pooled,
!   merged or policed.  Only a consumer that integrates a filler's
!   concentration over a depth the filler does not occupy (the FV-MOM6
!   cross-face Boole integral, the BBL stratification, Redi) reads through
!   the donor substitution below.
!
!   Substituting something OTHER than `c_live` on a vanished layer is
!   legal when the consumer answers a different question.  Do NOT route
!   these through `rdb_vl_column_conc`; document any new one in the table:
!
!     consumer                          substitution            why
!     rdb_eos (eos_*_impl)              reference T/S, rho=rho0 a filler must not perturb the PGF density column
!     rdb_ocean_sponge                  nearest massive layer   relaxation target must be a physical water mass
!     diag fills / rho_layer / MLD      NaN / never marks       a plot shows a gap; a filler's rho0 is no pycnocline
!     rdb_ocean_pseudo_salt deviation   NaN                     0 is the perfect score, a filler must not read as one
!     diag output-vcoord remap          zero weight             a NaN would be smeared over the target cells
!     rdb_ocean_cavity_flux sampler     `cycle`                 a filler carries no water to melt against
!     rdb_ocean_kappa_shear             massless merge          the solve needs a well-conditioned grid
!
!   A filler's own stored concentration is only safe where the value is
!   multiplied by the filler's OWN thickness.  The FV-MOM6 PGF's cross-face
!   Boole integral reads a filler's T/S over tens of metres of live water, so
!   it reads `c_live` off the donor
!   (`test_ocean_pgf_insitu :: open_step_filler_faces_*`).
!   The FV-lite / FV-Wright `rho_face*dz_centre` correction has the same
!   hazard and is not yet ported.
!
! The consuming module must already `use rdb_constants, only: wp, H_VANISHED`
! (and `NZ_STACK_MAX` if it includes the column procedures).  Everything here
! is `pure` and device-callable.
!
! Do NOT add state, logging, or anything impure to this file.

   pure logical function rdb_vl_is_live(h)
      !$acc routine seq
      !! THE live-layer predicate.  A layer at or below `H_VANISHED` is a
      !! vanished / inert filler: it holds no water of its own.  The test is a STRICT `>` — a layer
      !! sitting exactly ON the marker reads as vanished, which is what the
      !! geometric vcoord families rely on (`rdb_vcoord :: vcoord_h_min_role`).
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      rdb_vl_is_live = h > H_VANISHED
   end function rdb_vl_is_live

   pure subroutine rdb_vl_column_conc(nz, h_col, c_in, c_col)
      !$acc routine seq
      !! DONOR-SUBSTITUTED concentration of every layer of a column: the
      !! stored `c` on a live layer; on a vanished layer, the concentration of
      !! its DONOR (the nearest live layer above it; the topmost live layer
      !! for a run of fillers that reaches the top of the column).  A column
      !! with no live layer reads `0` throughout.
      !!
      !! The tracer is stored as concentration, so a filler holds a finite
      !! value of its own (the last water it held).  That value is right for
      !! transport (it is weighted by the filler's own thickness) and WRONG
      !! for a consumer that integrates it over a depth the filler does not
      !! occupy (the FV-MOM6 cross-face Boole integral, the BBL stratification):
      !! those read through this substitution.
      integer, intent(in) :: nz
         !! Number of layers in the column; `k = 1` is the bed.
      real(wp), intent(in) :: h_col(NZ_STACK_MAX)
         !! Layer thicknesses (m).
      real(wp), intent(in) :: c_in(NZ_STACK_MAX)
         !! Stored layer concentration.
      real(wp), intent(out) :: c_col(NZ_STACK_MAX)
         !! Layer concentration with fillers replaced by their donor's.
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
            c_col(k_vl) = c_in(k_vl)
            k_don_vl = k_vl
         else
            c_col(k_vl) = c_in(k_don_vl)
         end if
      end do
   end subroutine rdb_vl_column_conc
