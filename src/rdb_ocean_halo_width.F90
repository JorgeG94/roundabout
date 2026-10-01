!! One `required_halo` query folding every scattered minimum-`nghost` rule
!! (Python runtime API plan, P2.5).
module rdb_ocean_halo_width
   !! Roundabout has (at the time of writing) six independent minimum-`nghost`
   !! rules spread across five files, each enforced at its own call site in
   !! `rdb_config.F90`/`rdb_ocean_boundary_types.F90`/`rdb_ocean_halo.F90`/
   !! `rdb_ocean_setup.F90`. A Python caller building a grid in memory
   !! (`docs/ocean_python_api_plan.md` S5b) needs to be able to ask "how
   !! wide must my halo be" WITHOUT reproducing that scatter — the model
   !! Oceananigans uses is a per-scheme trait folded by `max()` over every
   !! tendency term (`inflate_halo_size`, `automatic_halo_sizing.jl:73-81`;
   !! `06_python_surface_design.md` B2.4). `required_halo` below is that
   !! fold, reusing the EXISTING per-rule functions (`pv_adv_required_nghost`,
   !! `tracer_recon_required_nghost`) rather than re-deriving their numbers,
   !! so this module and the enforcement sites it does not replace can never
   !! silently disagree.
   !!
   !! Two deliberate departures from Oceananigans, both already Roundabout's
   !! existing policy and kept here on purpose (`06_python_surface_design.md`
   !! B2.4): fail loud, never silently inflate the halo or downgrade a
   !! scheme's order to fit a small grid — `required_halo` only ADVISES the
   !! minimum; the fail-loud enforcement stays exactly where it already is
   !! (`rdb_config.F90`'s `validate_config`, `ocean_bc_validate_periodic`,
   !! `ocean_bc_validate_fold`, `ocean_halo_init`, `configure_ocean_kappa_shear`
   !! or equivalent) — this module does not remove or duplicate any of that
   !! enforcement.
   !!
   !! The six rules folded here:
   !!   1. global floor                        1   (`rdb_config.F90` nghost >= 1)
   !!   2. PV-advection weno5/weno7             3/4 (`rdb_coriolis_adv::pv_adv_required_nghost`)
   !!   3. tracer reconstruction weno5/7/9      3/4/5 (`rdb_recon_weno::tracer_recon_required_nghost`)
   !!   4. any PERIODIC topology                3   (`rdb_ocean_boundary_types::ocean_bc_validate_periodic`)
   !!   5. tripolar north fold                  3   (`rdb_ocean_boundary_types::ocean_bc_validate_fold`)
   !!   6. MPI-decomposed run                   3   (`rdb_ocean_halo` halo-exchange stencil depth)
   !!   (+) kappa-shear at_vertex               2   (`rdb_ocean_setup::configure_ocean_kappa_shear`)
   use rdb_coriolis_adv, only: pv_adv_required_nghost, parse_pv_adv_scheme
   use rdb_recon_weno, only: tracer_recon_required_nghost, parse_tracer_recon
   implicit none
   private

   public :: required_halo

   integer, parameter :: REQUIRED_HALO_FLOOR = 1
      !! Global minimum `nghost` regardless of scheme (rule 1).
   integer, parameter :: REQUIRED_HALO_PERIODIC = 3
      !! Minimum `nghost` for any periodic axis (rule 4) — PPM 5-point +
      !! biharmonic stencil depth.
   integer, parameter :: REQUIRED_HALO_TRIPOLAR_FOLD = 3
      !! Minimum `nghost` for the tripolar north fold (rule 5) — same
      !! stencil-depth requirement as periodic (the fold implies periodic
      !! west/east).
   integer, parameter :: REQUIRED_HALO_DECOMPOSED = 3
      !! Minimum `nghost` for an MPI-decomposed (multi-rank) run (rule 6).
   integer, parameter :: REQUIRED_HALO_KAPPA_SHEAR_VERTEX = 2
      !! Minimum `nghost` for kappa-shear's opt-in vertex form
      !! (`&ocean_kappa_shear_nml at_vertex = .true.`).

contains

   pure function required_halo(pv_adv_scheme, tracer_recon, periodic, &
                               tripolar_fold, decomposed, kappa_shear_at_vertex) result(ng)
      !! `max()` over every rule that applies to the given scheme
      !! selection — the minimum `nghost` a grid must carry. Every argument
      !! is optional; an absent selector contributes its BASELINE (2, the
      !! centered/PPM/non-periodic/single-rank floor), an absent logical
      !! contributes nothing (as if `.false.`).
      character(len=*), intent(in), optional :: pv_adv_scheme
         !! `&ocean_coriolis_nml pv_adv_scheme` string ("centered", "weno3",
         !! "weno5", "weno7"). Unrecognised/absent -> baseline (2).
      character(len=*), intent(in), optional :: tracer_recon
         !! `&ocean_tracers_nml recon`-style string ("ppm", "weno5",
         !! "weno7", "weno9"). Unrecognised/absent -> baseline (2).
      logical, intent(in), optional :: periodic
         !! `.true.` if EITHER grid axis is periodic (rule 4).
      logical, intent(in), optional :: tripolar_fold
         !! `.true.` for a tripolar north-fold grid (rule 5).
      logical, intent(in), optional :: decomposed
         !! `.true.` for a multi-rank (MPI-decomposed) run (rule 6).
      logical, intent(in), optional :: kappa_shear_at_vertex
         !! `.true.` when `&ocean_kappa_shear_nml at_vertex = .true.`.
      integer :: ng

      character(len=:), allocatable :: pv_s, tr_s

      pv_s = ""
      if (present(pv_adv_scheme)) pv_s = trim(pv_adv_scheme)
      tr_s = ""
      if (present(tracer_recon)) tr_s = trim(tracer_recon)

      ng = REQUIRED_HALO_FLOOR
      ng = max(ng, pv_adv_required_nghost(parse_pv_adv_scheme(pv_s)))
      ng = max(ng, tracer_recon_required_nghost(parse_tracer_recon(tr_s)))
      if (present(periodic)) then
         if (periodic) ng = max(ng, REQUIRED_HALO_PERIODIC)
      end if
      if (present(tripolar_fold)) then
         if (tripolar_fold) ng = max(ng, REQUIRED_HALO_TRIPOLAR_FOLD)
      end if
      if (present(decomposed)) then
         if (decomposed) ng = max(ng, REQUIRED_HALO_DECOMPOSED)
      end if
      if (present(kappa_shear_at_vertex)) then
         if (kappa_shear_at_vertex) ng = max(ng, REQUIRED_HALO_KAPPA_SHEAR_VERTEX)
      end if
   end function required_halo

end module rdb_ocean_halo_width
