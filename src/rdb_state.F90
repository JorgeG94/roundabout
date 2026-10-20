!! Shared tracer-registry defaults for the ocean dyn-core.
module rdb_state
   !! Overlay of the scalar `&tracer_nml` configuration onto the salinity /
   !! temperature slots of a `tracer_t` registry.
   !!
   !! Historically this module also carried the A-grid coastal `state_t`
   !! god-object; the ocean path composes `ocean_state_t` instead
   !! (`src/core/ocean/state/rdb_ocean_state.F90`), so all that survives here
   !! is the registry overlay both the driver and the benches call at setup.
   use rdb_config, only: config_t
   use rdb_tracer, only: tracer_t
   implicit none
   private

   public :: register_default_tracers

contains

   subroutine register_default_tracers(tr_S, tr_T, cfg)
      !! Overlay scalar config onto the already-allocated salinity /
      !! temperature tracer slots (by reference, so every caller shares the
      !! overlay).  Does not touch hTr / hTr0 — those are populated by the
      !! path's IC step once h_layer is known.
      type(tracer_t), intent(inout) :: tr_S
      type(tracer_t), intent(inout) :: tr_T
      type(config_t), intent(in)    :: cfg

      ! Salinity
      tr_S%name = "salinity"
      tr_S%long_name = "Sea water salinity"
      tr_S%units = "PSU"
      tr_S%standard_name = "sea_water_salinity"
      tr_S%tr_init = cfg%initial_salinity
      tr_S%tr_inflow = cfg%inflow_salinity
      tr_S%tr_min = cfg%S_min
      tr_S%tr_max = cfg%S_max
      tr_S%kappa_bg = cfg%kappa_S_bg
      tr_S%hdiff_kappa = cfg%hdiff_kappa
      tr_S%eos_coeff = cfg%beta_S
      tr_S%eos_ref = cfg%S_ref

      ! Temperature
      tr_T%name = "temperature"
      tr_T%long_name = "Sea water potential temperature"
      tr_T%units = "degC"
      tr_T%standard_name = "sea_water_potential_temperature"
      tr_T%tr_init = cfg%initial_temperature
      tr_T%tr_inflow = cfg%inflow_temperature
      tr_T%tr_min = cfg%T_min
      tr_T%tr_max = cfg%T_max
      tr_T%kappa_bg = cfg%kappa_T_bg
      tr_T%hdiff_kappa = cfg%hdiff_kappa
      tr_T%eos_coeff = -cfg%alpha_T
      tr_T%eos_ref = cfg%T_ref

   end subroutine register_default_tracers

end module rdb_state
