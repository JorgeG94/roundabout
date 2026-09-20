!! Ocean-side EOS shim over the C-grid multilayer state.
module rdb_ocean_eos_compute
   !! Hosts `ocean_eos_compute`, the outer shim that pulls the
   !! registered S, T tracer arrays off `multilayer_state_t` and
   !! forwards them as bare 3D arrays to the regime-agnostic
   !! `eos_compute_arrays` dispatch in `rdb_eos`.
   !!
   !! This shim is the ONLY EOS code that touches the ocean C-grid
   !! state, and it lives here — not in `rdb_eos` — so the shared
   !! EOS module carries zero ocean dependencies: the coastal opt-in
   !! nonlinear EOS (`&tracer_nml eos_form`, via
   !! `rdb_ml_eos::ml_equation_of_state_shared`) reaches
   !! `eos_compute_arrays` without transitively dragging in
   !! `rdb_multilayer_state`.
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_compute_arrays
   implicit none
   private

   public :: ocean_eos_compute

contains

   subroutine ocean_eos_compute(eos, ms, active)
      !! Compute `ms%rho_layer` from the registered (T, S) tracers.
      !! Outer-shim pattern: pulls the tracer hTr arrays off the
      !! registry on the host and forwards them as bare 3D arrays
      !! to a flat-impl chosen by `eos%variant`.
      !!
      !! Supported variants:
      !!   EOS_VARIANT_LINEAR    — `eos_linear_impl`
      !!   EOS_VARIANT_WRIGHT_97 — `eos_wright_impl` at `eos%p_ref`
      !!
      !! `ms%rho_layer` is a POTENTIAL density at the single, horizontally
      !! uniform `eos%p_ref` (`&ocean_eos_nml p_ref`).  It is deliberately
      !! NOT offset by the surface load `ms%p_top`: its consumers
      !! difference it along a layer and vertically, so a per-column
      !! reference pressure would manufacture a spurious along-layer
      !! density gradient — see the contract on `eos_compute_arrays`.
      !! The N² builders in `rdb_ocean_vmix` inherit that uniform
      !! reference and stay consistent with it.
      !!
      !! The in-situ-pressure Wright branch lives in the FV-PGF
      !! column sweep (`eos_wright_pgf_column_sweep_impl`), which
      !! owns its own pressure/density column scratch on the PGF
      !! state — not an EOS field.  THAT is where `ms%p_top` lands,
      !! because there the pressure is a true per-layer hydrostatic one.
      !!
      !! Optional `active` lets the dyn-step driver call this
      !! unconditionally — when present and false (thermodynamics
      !! disabled or non-thermo substep) the kernel is a no-op.
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      logical, intent(in), optional :: active

      integer :: nx, ny, nz

      if (present(active)) then
         if (.not. active) return
      end if

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = ms%nz_ml

      if (ms%idx_salinity <= 0 .or. ms%idx_temperature <= 0) return

      call eos_compute_arrays(eos, ms%h_layer, &
                              ms%tracers(ms%idx_salinity)%hTr, &
                              ms%tracers(ms%idx_temperature)%hTr, &
                              ms%rho_layer, nx, ny, nz)
   end subroutine ocean_eos_compute

end module rdb_ocean_eos_compute
