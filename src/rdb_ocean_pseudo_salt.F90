!! Pseudo-salt verification tracer for the ocean path.
module rdb_ocean_pseudo_salt
   !! Pseudo-salt (Shao 2016): a passive tracer seeded to salinity's
   !! initial condition and given exactly the operators salinity
   !! receives that the registry does NOT deliver automatically — the
   !! surface salt flux and the KPP/EPBL non-local counter-gradient
   !! `gamma_s` (both mirrored in-place in the kernels that own them,
   !! `rdb_ocean_surface_flux` / `rdb_ocean_vmix` — this module adds no
   !! new dyn-step call site).  Every other operator (advection, ALE
   !! remap, vertical/horizontal diffusion, vertical exchange, halo,
   !! sponge, OBC) already rides the generic tracer registry, so
   !! pseudo-salt receives it "for free" the moment it is registered.
   !!
   !! The deviation `D = pseudo_salt - S` is a direct, measured proxy
   !! for how far the passive-tracer transport path has drifted from
   !! the active-tracer (salinity) path — an architectural claim turned
   !! into a number.  Reference: Shao, A. (2016), MOM6
   !! `pseudo_salt_tracer.F90`; the idea is standard practice in
   !! offline/online transport verification and carries no separate
   !! citation in MOM6 either.
   !!
   !! Excluded by design (fail-loud at configure, `rdb_config`
   !! `validate_config`) because Roundabout has no mirror for them and an
   !! un-mirrored salinity source would make `D` measure that source
   !! instead of the passive-path error: SSS piston restoring
   !! (`&ocean_restore_nml enable_restore_salt`) and sea-ice frazil /
   !! basal salt exchange (`&ocean_ice_nml enable`).
   use rdb_constants, only: wp, H_VANISHED, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use pic_logger, only: global_logger
   use rdb_error_ring, only: error_ring_push
   implicit none
   private

   public :: ocean_pseudo_salt_register
   public :: ocean_pseudo_salt_seed
   public :: ocean_pseudo_salt_deviation
   public :: pseudo_salt_conflicts_restore
   public :: pseudo_salt_conflicts_ice
   public :: pseudo_salt_needs_thermo_warning

contains

   subroutine ocean_pseudo_salt_register(ms, grid)
      !! Setup-time entry point: registers the "pseudo_salt" passive
      !! tracer (name/units/long_name copied verbatim from Shao 2016)
      !! and records its slot on `ms%idx_pseudo_salt`.  MUST be called
      !! after `ms%init` and before `enter_data` — and, on the ocean
      !! path, before `ocean_bc_state_init` sizes `bc%n_tracers` (the
      !! `register_passive_tracer` contract).  `error stop`s on
      !! refusal (mirrors `rdb_state.F90`'s treatment of a fatal setup
      !! misconfiguration — a caller that reaches here has already
      !! decided to register, so a silent no-op would be worse than a
      !! loud abort).
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid

      integer :: idx

      call ms%register_passive_tracer(grid, "pseudo_salt", "psu", &
                                      "Pseudo salt passive tracer", idx)
      if (idx <= 0) then
         call error_ring_push("ocean_pseudo_salt_register: "// &
                              "register_passive_tracer refused (registry "// &
                              "not initialised, or already locked by enter_data)")
         call global_logger%error("ocean_pseudo_salt_register: "// &
                                  "register_passive_tracer refused (registry "// &
                                  "not initialised, or already locked by enter_data)")
         error stop "ocean_pseudo_salt_register: registration refused"
      end if
      ms%idx_pseudo_salt = idx
   end subroutine ocean_pseudo_salt_register

   pure subroutine ocean_pseudo_salt_seed(ms)
      !! Seed `hTr(idx_pseudo_salt) = hTr(idx_salinity)` over the FULL
      !! array shape (nx_total, ny_total, nz_ml) — ghosts included,
      !! matching MOM6's `isd:ied, jsd:jed` seed (§4 of the plan). A
      !! halved ghost band would make the deviation diagnostic non-zero
      !! at the very first halo exchange, before any real transport has
      !! run. Self-gates on either index being unregistered. Must run
      !! AFTER every write to salinity's initial condition (analytical
      !! IC, then any z-file overlay) and BEFORE `ocean_state_enter_data`.
      type(multilayer_state_t), intent(inout) :: ms

      if (ms%idx_pseudo_salt <= 0 .or. ms%idx_salinity <= 0) return
      ms%tracers(ms%idx_pseudo_salt)%hTr = ms%tracers(ms%idx_salinity)%hTr
   end subroutine ocean_pseudo_salt_seed

   pure subroutine ocean_pseudo_salt_deviation(h_layer, hTr_ps, hTr_s, buf, nx, ny, nz, &
                                               missing)
      !! Diagnostic helper: `D(i,j,k) = hTr_ps/h - hTr_s/h`, i.e. the
      !! pseudo-salt concentration minus the salinity concentration.
      !!
      !! A vanished layer (`rdb_vl_is_live` false — the ONE predicate, see
      !! `src/core/ocean/README.md`, "The vanished-layer content rule")
      !! writes `missing`, which the caller passes as the IEEE NaN sentinel
      !! `fill_tracer_impl` uses for the two concentrations this is the
      !! difference of.  It used to write `0` below a private `1e-6 m`
      !! floor: a deviation of zero is the one value this diagnostic
      !! exists to report as "the two transport paths agree", so a filler
      !! read as a perfect score, and layers between `1e-6 m` and
      !! `H_VANISHED` were divided through at all.
      integer, intent(in)     :: nx, ny, nz
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: hTr_ps(nx, ny, nz)
      real(wp), intent(in)    :: hTr_s(nx, ny, nz)
      real(wp), intent(inout) :: buf(nx, ny, nz)
      real(wp), intent(in)    :: missing
         !! Value written on a vanished layer (the diagnostics' NaN).
      integer :: i, j, k
      real(wp) :: h

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         h = h_layer(i, j, k)
         if (rdb_vl_is_live(h)) then
            buf(i, j, k) = rdb_vl_conc(hTr_ps(i, j, k), h) - rdb_vl_conc(hTr_s(i, j, k), h)
         else
            buf(i, j, k) = missing
         end if
      end do
   end subroutine ocean_pseudo_salt_deviation

   pure function pseudo_salt_conflicts_restore(enable_ps, enable_restore_salt) result(conflict)
      !! `.true.` iff pseudo-salt is enabled alongside SSS piston
      !! restoring — an un-mirrored salinity source
      !! (`rdb_ocean_surface_flux.F90` restore branch) that would make
      !! the deviation diagnostic measure the restoring term instead of
      !! the passive-transport-path error. Drives a configure-time
      !! fail-loud abort.
      logical, intent(in) :: enable_ps
      logical, intent(in) :: enable_restore_salt
      logical :: conflict
      conflict = enable_ps .and. enable_restore_salt
   end function pseudo_salt_conflicts_restore

   pure function pseudo_salt_conflicts_ice(enable_ps, enable_ice) result(conflict)
      !! `.true.` iff pseudo-salt is enabled alongside the sea-ice
      !! model — frazil / basal salt exchange are un-mirrored salinity
      !! sources (`rdb_ice_frazil.F90`, `rdb_ice_frazil_uptake.F90`,
      !! `rdb_ice_basal_flux.F90`). Drives a configure-time fail-loud
      !! abort.
      logical, intent(in) :: enable_ps
      logical, intent(in) :: enable_ice
      logical :: conflict
      conflict = enable_ps .and. enable_ice
   end function pseudo_salt_conflicts_ice

   pure function pseudo_salt_needs_thermo_warning(enable_ps, enable_thermodynamics) &
      result(warn)
      !! `.true.` iff pseudo-salt is enabled with thermodynamics off.
      !! NOT a hard error — pseudo-salt still measures pure transport,
      !! a legitimate use — but the caller should warn: without
      !! thermodynamics, the tracer never receives the surface salt
      !! flux / KPP nonlocal mirrors, so `D` degenerates to a pure
      !! advection/diffusion probe.
      logical, intent(in) :: enable_ps
      logical, intent(in) :: enable_thermodynamics
      logical :: warn
      warn = enable_ps .and. .not. enable_thermodynamics
   end function pseudo_salt_needs_thermo_warning

#include "rdb_vanished_layer.inc"

end module rdb_ocean_pseudo_salt
