!! Compile-time bounds for the sea-ice column code (roundabout-original).
module rdb_ice_dims
   !! Fixed upper bounds used to size stack-local column work arrays inside
   !! `!$acc routine seq` procedures, where an automatic or allocatable array
   !! would stall the device.  The RUNTIME extents `ncat`, `nilyr`, `nslyr`
   !! (namelist, validated against these bounds) are always `<=` the bounds.
   !!
   !! **Layer order.**  The sea-ice column indexes its layers TOP-DOWN: `k = 1`
   !! is the layer under the surface and `k = nilyr` is the layer at the ice
   !! base (Bitz and Lipscomb 1999 vertical grid).  This is a documented
   !! exception to the ocean's bottom-up rule (`k = nz` is the surface): the
   !! ice never exchanges layer-by-layer with the ocean, so the two orders
   !! never meet.  Every ice module states this where it indexes layers.
   !!
   !! **Tracer count.**  The packed tracer column of the Bitz-Lipscomb
   !! thermodynamics holds one surface temperature, `nilyr` ice enthalpies,
   !! `nslyr` snow enthalpies and `nilyr` ice salinities, so
   !! `ICE_NTRCR_MAX = 1 + 2*ICE_NILYR_MAX + ICE_NSLYR_MAX`.
   implicit none
   private

   public :: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX, ICE_NTRCR_MAX
   public :: ICE_NTHX, ICE_TH_EVAPI, ICE_TH_EVAPS, ICE_TH_MELTS

   integer, parameter :: ICE_NCAT_MAX = 8
      !! Maximum number of thickness categories.
   integer, parameter :: ICE_NILYR_MAX = 8
      !! Maximum number of ice layers per category.
   integer, parameter :: ICE_NSLYR_MAX = 1
      !! Maximum number of snow layers per category.
   integer, parameter :: ICE_NTRCR_MAX = 1 + 2*ICE_NILYR_MAX + ICE_NSLYR_MAX
      !! Maximum length of the packed tracer column (= 18).

   integer, parameter :: ICE_NTHX = 3
      !! Number of extra per-category thermodynamic history outputs (the last extent of `thxn`, the
      !! first two extents of the cell array `thx`).
   integer, parameter :: ICE_TH_EVAPI = 1
      !! Slot: evaporation (< 0 sublimation) of the ICE part of the category, [kg m-2 s-1] per category area.
   integer, parameter :: ICE_TH_EVAPS = 2
      !! Slot: evaporation of the SNOW part, [kg m-2 s-1] per category area.
   integer, parameter :: ICE_TH_MELTS = 3
      !! Slot: snow melt, [m] of snow per category area per step (+ melt).

end module rdb_ice_dims
