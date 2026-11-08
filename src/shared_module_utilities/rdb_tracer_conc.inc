! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_tracer_conc.inc -- THE tracer-concentration accessor.
!
! `#include` this file in the `contains` section of any module whose kernels
! need the concentration `c` of one layer.  Like `rdb_vanished_layer.inc` you
! get a module-LOCAL, `private` copy of each body, which NVHPC's device codegen
! inlines into a `do concurrent` kernel (README.md in this directory).
!
! Why this file exists.  Tracers are stored as CONTENT, `hTr = h*c`, so every
! reader recovers `c = hTr/h`, and each reader historically spelled the divide
! (and its guard against a thin or empty layer) in its own way.  The tracer
! prognostic is going to become the concentration itself; when it does, every
! accessor below collapses to a plain load of `c(i,j,k)` and the thickness
! argument disappears.  Routing every reader through ONE named accessor now
! makes that change local: one call site per read, no divide to find and
! reason about.  Today each accessor is a thin wrapper whose expression TEXT is
! the one the call site used before, so answers are bit-identical.
!
! The four accessors differ ONLY in the regularisation of the divide, which is
! a property of the consumer and is kept exactly as it was:
!
!   rdb_conc_of(hTr, h)          h > 0 ? hTr/h : 0     (the invariant-I1' read;
!                                  identical to `rdb_vl_conc`)
!   rdb_conc_of_raw(hTr, h)      hTr/h                 (transport kernels: the
!                                  caller guarantees h > 0 on every cell read)
!   rdb_conc_of_floor(hTr,h,hf)  hTr/max(h, hf)        (solvers / deposits that
!                                  floor the thickness)
!   rdb_conc_of_eps(hTr,h,eps)   hTr/(h + eps)         (1/0 armour)
!
! Write-side twins (disappear, or collapse to a plain store, with the storage
! change):
!
!   rdb_content_of_floor(c,h,hf) c*max(h, hf)          (solver reconstitute)
!   rdb_content_of(c, h)         c*h                   (seed / reset / OBC ghost)
!   rdb_conc_add_source(hTr,h,src)             hTr + src   (deposit, h unchanged)
!
! The consuming module must already `use rdb_constants, only: wp`.  Everything
! here is `pure` and device-callable.  Do NOT add state, logging or anything
! impure to this file.

   pure real(wp) function rdb_conc_of(hTr, h)
      !$acc routine seq
      !! Concentration of ONE layer from its content and thickness, `0` on a
      !! zero-thickness layer (the I1′-aware read; see `rdb_vl_conc`).
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`).
      real(wp), intent(in) :: h
         !! Layer thickness (m) the content is paired with.
      if (h > 0.0_wp) then
         rdb_conc_of = hTr/h  ! vanished-ok: the accessor's own definition
      else
         rdb_conc_of = 0.0_wp
      end if
   end function rdb_conc_of

   pure real(wp) function rdb_conc_of_raw(hTr, h)
      !$acc routine seq
      !! Concentration `hTr/h` with NO guard: the caller guarantees `h > 0`
      !! (or accepts the non-finite result, as the kernel always did).
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`).
      real(wp), intent(in) :: h
         !! Layer thickness (m) the content is paired with.
      rdb_conc_of_raw = hTr/h  ! vanished-ok: the accessor's own definition
   end function rdb_conc_of_raw

   pure real(wp) function rdb_conc_of_floor(hTr, h, h_floor)
      !$acc routine seq
      !! Concentration `hTr/max(h, h_floor)`.
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`).
      real(wp), intent(in) :: h
         !! Layer thickness (m) the content is paired with.
      real(wp), intent(in) :: h_floor
         !! Thickness floor (m) of the consumer.
      rdb_conc_of_floor = hTr/max(h, h_floor)  ! vanished-ok: the accessor's own definition
   end function rdb_conc_of_floor

   pure real(wp) function rdb_conc_of_eps(hTr, h, eps)
      !$acc routine seq
      !! Concentration `hTr/(h + eps)`.
      real(wp), intent(in) :: hTr
         !! Layer content (`h·c`).
      real(wp), intent(in) :: h
         !! Layer thickness (m) the content is paired with.
      real(wp), intent(in) :: eps
         !! Additive 1/0 armour (m).
      rdb_conc_of_eps = hTr/(h + eps)  ! vanished-ok: the accessor's own definition
   end function rdb_conc_of_eps

   pure real(wp) function rdb_content_of_floor(c, h, h_floor)
      !$acc routine seq
      !! The WRITE-side twin of `rdb_conc_of_floor`: content `c*max(h, h_floor)`
      !! from a concentration and the SAME floored thickness the read used, so
      !! a solve that reads `c` through the floor reconstitutes through it too.
      !! Disappears with the storage change (the solve then stores `c`).
      real(wp), intent(in) :: c
         !! Concentration.
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      real(wp), intent(in) :: h_floor
         !! Thickness floor (m) of the consumer; must equal the read's.
      rdb_content_of_floor = c*max(h, h_floor)
   end function rdb_content_of_floor

   pure real(wp) function rdb_content_of(c, h)
      !$acc routine seq
      !! Content `c*h` of a layer seeded at concentration `c` (the write-side
      !! twin of `rdb_conc_of_raw`; no guard, no floor).  Disappears with the
      !! storage change (the seed then stores `c`).
      real(wp), intent(in) :: c
         !! Concentration.
      real(wp), intent(in) :: h
         !! Layer thickness (m).
      rdb_content_of = c*h
   end function rdb_content_of

   pure real(wp) function rdb_conc_add_source(hTr, h, src)
      !$acc routine seq
      !! Tracer state after adding a CONTENT source `src` (`h·Δc`, units of
      !! `hTr`) to a layer of thickness `h` that the operator does not move:
      !! `hTr + src`.  The thickness argument is the one the concentration
      !! form divides by (`c + src/(h + h_neglect)`); it is unused while the
      !! tracer is stored as content.
      real(wp), intent(in) :: hTr
         !! Layer content before the source.
      real(wp), intent(in) :: h
         !! Layer thickness (m), unchanged by the operator.
      real(wp), intent(in) :: src
         !! Content source over the step (`dt·tendency`).
      rdb_conc_add_source = hTr + src
   end function rdb_conc_add_source
