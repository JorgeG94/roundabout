! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_tracer_conc.inc -- tracer-concentration source deposits.
!
! `#include` this file in the `contains` section of any module whose kernels
! deposit a tracer SOURCE into a layer.  Like `rdb_vanished_layer.inc` you get
! a module-LOCAL, `private` copy of each body, which NVHPC's device codegen
! inlines into a `do concurrent` kernel.
!
! The tracer prognostic is the CONCENTRATION `c` (Griffies, Adcroft &
! Hallberg 2020), so a reader loads `c(i,j,k)` directly and no accessor is
! needed.  What stays here are the two deposits whose divide must be written
! once and written right: an extensive source `src` (content, `h*dc`) added to
! a layer that either keeps its thickness or changes it in the SAME operator.
!
!   rdb_conc_add_source(c, h, src)               c + src/(h + h_neglect)
!   rdb_conc_add_source_mass(c, hb, ha, src)     (c*hb + src)/ha
!
! `h_neglect` is `H_DIV_EPS` (1e-20 m): far below any thickness that carries
! content we care about and above `tiny()`.  The consuming module must already
! `use rdb_constants, only: wp, H_DIV_EPS`.  Everything here is `pure` and
! device-callable.  Do NOT add state, logging or anything impure to this file.

   pure real(wp) function rdb_conc_add_source(c, h, src)
      !$acc routine seq
      !! Concentration after adding a CONTENT source `src` (`h·Δc`, units of
      !! `c·m`) to a layer of thickness `h` that the operator does not move:
      !! `c + src/(h + h_neglect)`.
      real(wp), intent(in) :: c
         !! Layer concentration before the source.
      real(wp), intent(in) :: h
         !! Layer thickness (m), unchanged by the operator.
      real(wp), intent(in) :: src
         !! Content source over the step (`dt·tendency`).
      rdb_conc_add_source = c + src/(h + H_DIV_EPS)
   end function rdb_conc_add_source

   pure real(wp) function rdb_conc_add_source_mass(c, h_before, h_after, src)
      !$acc routine seq
      !! Concentration after adding a CONTENT source `src` to a layer whose
      !! thickness the SAME operator changes, `h_before -> h_after` (melt
      !! with mass): `(c*h_before + src)/h_after`.  A layer the operator leaves
      !! empty (`h_after <= h_neglect`) keeps its `c`.
      real(wp), intent(in) :: c
         !! Layer concentration before the operator.
      real(wp), intent(in) :: h_before
         !! Layer thickness (m) before the operator.
      real(wp), intent(in) :: h_after
         !! Layer thickness (m) the operator writes.
      real(wp), intent(in) :: src
         !! Content source over the step.
      if (h_after > H_DIV_EPS) then
         rdb_conc_add_source_mass = (c*h_before + src)/h_after
      else
         rdb_conc_add_source_mass = c
      end if
   end function rdb_conc_add_source_mass
