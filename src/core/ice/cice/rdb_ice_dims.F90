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

   integer, parameter :: ICE_NCAT_MAX = 8
      !! Maximum number of thickness categories.
   integer, parameter :: ICE_NILYR_MAX = 8
      !! Maximum number of ice layers per category.
   integer, parameter :: ICE_NSLYR_MAX = 1
      !! Maximum number of snow layers per category.
   integer, parameter :: ICE_NTRCR_MAX = 1 + 2*ICE_NILYR_MAX + ICE_NSLYR_MAX
      !! Maximum length of the packed tracer column (= 18).

end module rdb_ice_dims
