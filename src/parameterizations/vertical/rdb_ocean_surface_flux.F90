!! Surface heat + salt fluxes for the ocean dynamical core.  Sister
!! module to `rdb_ocean_surface_stress` — same 2D-field shape, applies
!! a tracer flux at the surface layer (`k = nz` under the ROMS-style
!! k=1-bed convention).
!!
!! Heat: `dT/dt|_surface = Q_heat / (rho_0 * cp * h_top)` with
!! Q_heat in W/m^2 (positive downward into the ocean).  Salt:
!! `dS/dt|_surface = +Q_salt / (rho_0 * h_top)` with `Q_salt` the net
!! virtual salt flux (PSU.kg m-2 s-1, POSITIVE SALINIFIES — evaporation
!! excess, or sea-ice brine rejection via `rdb_ice_ocean_coupler`).
!!
!! The 2D fields `Q_heat(:,:)` / `Q_salt(:,:)` are seeded at configure
!! time from the `&ocean_thermo_nml q_heat / q_salt` scalars via
!! `set_surface_flux_const`.  They live on the device between
!! `enter_data` and `exit_data`; Area-A3 data-override and Area-A4
!! restoring write directly into these fields after seeding.
!!
!! `Q_heat_const` / `Q_salt_const` are kept as the fill-path scalars
!! (mirroring `set_wind_stress_const` on the stress side) so
!! `set_surface_flux_const` produces identical arithmetic to the old
!! scalar broadcast — existing nmls and tests are bit-identical.
!!
!! Phase 2 of the KPP / EPBL build-out also reads the heat/salt fluxes
!! to compute B_0.  Those shimmed paths now pick up `Q_heat(i,j)` /
!! `Q_salt(i,j)` per column instead of the broadcast scalar.
!!
!! Restart: the fields are configure-time-filled static forcing; they
!! are NOT registered in the restart registry (re-seeded from the
!! namelist scalar on every resume — see `rdb_driver.F90` for the
!! seeding call-site).
!!
!! **Component set (PR-12, `&ocean_forcing_nml enable_components`).**
!! `Q_heat` / `Q_salt` are the derived views every downstream kernel
!! keeps reading; when `use_components` is on, `ocean_surface_flux_assemble`
!! rebuilds them every thermo step from `Q_heat_const`/`Q_salt_const`
!! plus the per-component fields below.  Fill contract for any future
!! filler (a reader, the ice coupler, a rivers PR, ...):
!!   (a) write your OWN component(s) — NEVER `Q_heat` / `Q_salt` directly
!!       (those are assembler-owned; a second writer is a review-rejectable
!!       error, see `ocean_surface_flux_assemble`'s docstring);
!!   (b) `!$acc update device(...)` your write (or write on-device) before
!!       the assembler's next call — components are device-resident between
!!       `enter_data`/`exit_data`;
!!   (c) set `has_heat`/`has_salt` (existing contract) and, as needed,
!!       `has_mass_flux` (after any mass-flux write) / `has_q_sw` (after a
!!       `q_sw` write) host-side — never from a device reduction;
!!   (d) fill `heat_content_<flux>` for every mass flux you fill (the
!!       source owns the enthalpy of the mass it injects; the ocean can
!!       only compute the enthalpy of mass it loses, `heat_content_massout`,
!!       which the assembler derives from SST — never write it yourself);
!! **Ice-shelf cover (`&ocean_cavity_dyn_nml`).**  Under a shelf there
!! is no atmosphere, so every ATMOSPHERIC contribution must be zero in a
!! covered cell and unchanged everywhere else, while the cavity's own
!! `heat_cavity` / `salt_cavity` must NOT be masked.  The mask is an
!! OPTIONAL `cover_frac` argument (absent ⇒ the original kernel,
!! byte-identical) and it is applied in the **assembler**, not at apply
!! time.  That choice is load-bearing, for two reasons:
!!
!!   1. `Q_heat` / `Q_salt` are not only the apply kernel's input — they
!!      are what KPP and EPBL read to build `B_0`.  Masking at apply
!!      time would leave both boundary-layer schemes forced by an
!!      atmosphere that is not there, with the tracer deposit correct
!!      and the mixing wrong: a plausible, publishable, wrong answer.
!!   2. The assembler is the ONE place where the atmospheric bands and
!!      the cavity bands are still distinguishable.  After it, `Q_heat`
!!      is a sum, and any factor applied to the sum would also scale
!!      `heat_cavity` — i.e. mask away the melt flux the cover is
!!      supposed to admit.
!!
!! The two surface kernels that do NOT route through `Q_heat`/`Q_salt`
!! carry the same optional argument and mask themselves:
!! `ocean_surface_flux_apply_sw_penetration` (it reads `q_sw`, a
!! pristine INPUT component, on the `sw_source="q_sw"` branch — and even
!! on the `net_heat` branch it must not *move* heat that the masked
!! deposit never added) and `ocean_surface_restore_apply_tracers` (it
!! forms its flux in-kernel from the live SST/SSS).  With the component
!! set OFF there is no assembler, so the static `Q_heat`/`Q_salt` fill is
!! masked once at configure by `ocean_surface_flux_apply_cover_const`.
!!
!!   (e) register your component in the restart registry if it is
!!       time-varying (`registry_register_2d`, `optional=.true.`,
!!       `device_mapped=.true.`) — the assembler's own outputs
!!       (`Q_heat`/`Q_salt`, `heat_content_massin/massout`) are never
!!       registered (derived-field rule, `rdb_ocean_state.F90`).
module rdb_ocean_surface_flux
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_surface_flux_t
   public :: ocean_surface_flux_apply_tracers
   public :: ocean_surface_flux_apply_sw_penetration
   public :: ocean_surface_restore_apply_tracers
   public :: ocean_surface_flux_assemble
   public :: ocean_surface_flux_apply_cover_const
   public :: sw_transmission
   public :: sw_pe_cost_shape
   public :: sw_source_is_implemented

   ! Seconds per day — converts the MOM6 `FLUXCONST` piston velocity
   ! (specified in m/day) to MKS m/s at the seed call.
   real(wp), parameter :: SECONDS_PER_DAY = 86400.0_wp

   ! Specific heat capacity of seawater (J/kg/K).  Wright (1997)
   ! Table A1 reference value; same as MOM6's `CP_SW` default.
   real(wp), parameter, public :: SEAWATER_CP = 3992.0_wp

   type :: ocean_surface_flux_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: has_heat = .false.
         !! True when `Q_heat` carries a non-zero fill (set by
         !! `set_surface_flux_const` when `q_heat_val /= 0`).
         !! Any future field-fill path (A3 data-override, A4 restoring)
         !! MUST set `has_heat = .true.` after writing into `Q_heat`
         !! so the apply-tracers kernel fires.  Host-side flag only
         !! (early-return guard in `ocean_surface_flux_apply_tracers`).
      logical :: has_salt = .false.
         !! True when `Q_salt` carries a non-zero fill (set by
         !! `set_surface_flux_const` when `q_salt_val /= 0`).
         !! Same contract as `has_heat` for any field-fill path.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3) — the `dt/(rho0*cp)` heat
         !! and `dt/rho0` salt divisors applied to EVERY surface tracer
         !! source, including whatever the sea-ice coupler writes into
         !! `Q_heat`/`Q_salt`.
         !!
         !! ASSIGNED FROM CONFIG by `configure_ocean_reference_density`,
         !! which copies the single rho0 of record (`&ocean_ic_nml rho_0`
         !! -> `eos%rho0`).  The literal here is only the pre-configure
         !! type default; do not read it as the value a run uses.  Host
         !! scalar: the divisor is folded into the `inv_scale` argument on
         !! the host, so the assignment owes no `!$acc update device`.
      real(wp) :: cp = SEAWATER_CP
         !! Specific heat capacity (J/kg/K).
      real(wp) :: h_min = 1.0e-3_wp
         !! Floor on the surface-layer thickness in the `1/h_top`
         !! division — keeps the kernel finite when the top layer
         !! pinches out.
      real(wp) :: Q_heat_const = 0.0_wp
         !! Scalar fill source for `Q_heat(:,:)`.  Seeded from
         !! `&ocean_thermo_nml q_heat` by `set_surface_flux_const`.
         !! Kept for diagnostic logging; kernels read `Q_heat` directly.
      real(wp) :: Q_salt_const = 0.0_wp
         !! Scalar fill source for `Q_salt(:,:)`.  Seeded from
         !! `&ocean_thermo_nml q_salt` by `set_surface_flux_const`.

      logical :: has_sw = .false.
         !! True when shortwave penetration is active (set by
         !! `set_sw_penetration` when `sw_pen_frac /= 0`).  Host-side
         !! gate only — `ocean_surface_flux_apply_sw_penetration`
         !! early-returns unless this is set, so the default-off path
         !! leaves the surface-flux deposition byte-for-byte unchanged.
      logical :: sw_from_qsw = .false.
         !! Selects the irradiance source for shortwave penetration and
         !! the boundary-layer SW coupling.  `.false.` (default): the
         !! source is the NET heat flux `Q_heat` (`sw_source="net_heat"`,
         !! legacy, bit-identical).  `.true.` (`sw_source="q_sw"`): the
         !! source is the dedicated `q_sw` component (>= 0), which
         !! removes the night-time negative-`I0` hazard where
         !! `sw_pen_frac*Q_heat < 0` drives unphysical negative
         !! irradiance down the two-band profile.  Set by
         !! `set_sw_penetration`; host-side gate only (never a device
         !! reduction) — the shim selects the source array on the host so
         !! the conditionally-allocated `q_sw` is never dereferenced in a
         !! device kernel.
      real(wp) :: sw_pen_frac = 0.0_wp
         !! Penetrating fraction of `Q_heat` carried below the surface
         !! layer as a two-band exponential (Paulson & Simpson 1977).
         !! 0 = off (all of `Q_heat` lands at `k = nz`, legacy path).
      real(wp) :: sw_band_ratio = 0.58_wp
         !! Band-1 weight `R` of the two-band irradiance decay.  Jerlov
         !! type I (clear open ocean) default.
      real(wp) :: sw_zeta1 = 0.35_wp
         !! Band-1 e-folding depth (m) — the rapidly-absorbed
         !! red/near-IR band.
      real(wp) :: sw_zeta2 = 23.0_wp
         !! Band-2 e-folding depth (m) — the slowly-absorbed
         !! blue/green band.
      logical :: has_restore_T = .false.
         !! True when SST restoring is active (`enable_restore_temp .and.
         !! restore_piston_T /= 0`).  Host-side gate only —
         !! `ocean_surface_restore_apply_tracers` early-returns unless
         !! this or `has_restore_S` is set, so the default-off path is
         !! byte-for-byte unchanged.
      logical :: has_restore_S = .false.
         !! True when SSS restoring is active (`enable_restore_salt .and.
         !! restore_piston_S /= 0`).
      real(wp) :: restore_piston_T = 0.0_wp
         !! SST piston velocity (m/s), seeded from `&ocean_restore_nml
         !! piston_t` (m/day) via `/86400`.  The surface relaxation rate
         !! for a top layer of thickness `h_top` is
         !! `lambda = restore_piston_T / h_top` [1/s].
      real(wp) :: restore_piston_S = 0.0_wp
         !! SSS piston velocity (m/s).
      real(wp) :: restore_T_target = 0.0_wp
         !! Scalar target SST (degC).  Read by-value into the device
         !! `_impl` kernel — no per-cell field (so no extra device
         !! array, the `enter_data` orchestrator is untouched).  A
         !! 2D-field target is the documented A4-v2 follow-up.
      real(wp) :: restore_S_target = 0.0_wp
         !! Scalar target SSS (PSU).
      real(wp), allocatable :: Q_heat(:, :)
         !! 2D net surface heat flux (W/m^2, positive downward),
         !! shape `(nx, ny)`.  Fill via `set_surface_flux_const` for
         !! spatially-uniform forcing (the default); Area-A3 override
         !! or Area-A4 restoring writes the field directly.
      real(wp), allocatable :: Q_salt(:, :)
         !! 2D net surface salt flux (kg salt/m^2/s, positive salinifies),
         !! shape `(nx, ny)`.

      ! ---- PR-12 component set — allocated iff `use_components` ----
      logical :: use_components = .false.
         !! Master gate (`&ocean_forcing_nml enable_components`).  `.false.`
         !! (default): none of the arrays below are allocated,
         !! `ocean_surface_flux_assemble` is a no-op, and `Q_heat`/`Q_salt`
         !! are filled exactly as today — bit-identical.  `.true.`:
         !! allocates the component set (`set_components`) and the
         !! assembler rebuilds `Q_heat`/`Q_salt` every thermo step.
      logical :: has_mass_flux = .false.
         !! Host-side latch — set by a filler after writing ANY of
         !! `evap`/`lprec`/`fprec`/`vprec`/`lrunoff`/`frunoff`/
         !! `seaice_melt`.  Cheap early-return gate for a future
         !! freshwater kernel (real-mass PR).  Same contract as
         !! `has_heat` (`:56-62`) — never set from a device reduction.
      logical :: has_q_sw = .false.
         !! Host-side latch — set by a filler after writing `q_sw`
         !! (e.g. an ice `sw_thru` coupler).  **Not** the same as
         !! `has_sw` below (that gates shortwave *penetration*, an
         !! unrelated pre-existing switch) — do not conflate the two.

      real(wp), allocatable :: q_sw(:, :)
         !! Shortwave into the ocean (W/m^2, **>= 0**, positive down).
      real(wp), allocatable :: q_lw(:, :)
         !! Net longwave (W/m^2, typically **< 0**, positive down).
      real(wp), allocatable :: q_lat(:, :)
         !! Latent heat flux (W/m^2, typically **< 0**, positive down).
      real(wp), allocatable :: q_sens(:, :)
         !! Sensible heat flux (W/m^2, typically **< 0**, positive down).
      real(wp), allocatable :: heat_added(:, :)
         !! Restoring / flux-adjustment / "other" net heat term not
         !! decomposed into the radiative/turbulent bands above (W/m^2,
         !! either sign; MOM6 `heat_added`).  The v1 ice coupler's
         !! `heat_flux_diag` lands here (§5.4 of the PR-12 plan) —
         !! it is already a net W/m^2, not further decomposable.
      real(wp), allocatable :: heat_cavity(:, :)
         !! **Ice-shelf cavity basal-melt heat component** (W/m^2, same
         !! positive-DOWN-into-the-ocean convention as every other heat
         !! band; `&ocean_cavity_melt_nml`).  OWNED by
         !! `rdb_ocean_cavity_flux`; written `heat_cavity = -q_ocean`,
         !! where `q_ocean = rho_w*c_w*gamma_t*(T_w - T_b) > 0` is the
         !! kernel's turbulent heat flux OCEAN -> INTERFACE, so warm
         !! water under a shelf COOLS the top of the column.  It is a
         !! SEPARATE field from `heat_added` precisely because the
         !! sea-ice coupler full-overwrites `heat_added` — two writers
         !! on one slot clobber silently (cavity x sea ice is refused
         !! today, but the ownership rule must not depend on that).
         !! Zero unless a cavity melt step ran.

      real(wp), allocatable :: evap(:, :)
         !! Evaporative mass flux (kg/m^2/s, **<= 0** — MOM6 convention,
         !! `(-1)*flux out of the ocean`).  v1: enthalpy + salt
         !! bookkeeping only — does NOT change column mass (real
         !! freshwater is a named follow-up, see the module docstring).
      real(wp), allocatable :: lprec(:, :)
         !! Liquid precipitation (kg/m^2/s, **>= 0** into the ocean).
      real(wp), allocatable :: fprec(:, :)
         !! Frozen precipitation / snowfall (kg/m^2/s, **>= 0**).
      real(wp), allocatable :: vprec(:, :)
         !! Virtual precipitation (kg/m^2/s, either sign — SSS-restoring
         !! convention; NOT wired to the restoring kernel in v1, see
         !! `set_restore` below and the module docstring's follow-up note).
      real(wp), allocatable :: lrunoff(:, :)
         !! Liquid river runoff (kg/m^2/s, **>= 0**).
      real(wp), allocatable :: frunoff(:, :)
         !! Frozen (calving/ice) runoff (kg/m^2/s, **>= 0**).
      real(wp), allocatable :: seaice_melt(:, :)
         !! Sea-ice melt-water mass flux (kg/m^2/s, **>0** = melt into
         !! the ocean, **<0** = formation / freezing withdraws mass).

      real(wp), allocatable :: heat_content_lprec(:, :)
         !! Enthalpy carried by `lprec` (W/m^2).  A filler that writes
         !! `lprec` MUST fill this — the v1 convenience is
         !! `SEAWATER_CP * T_source * lprec`.  The source (not the
         !! ocean) owns this enthalpy — see the module docstring §(d).
      real(wp), allocatable :: heat_content_fprec(:, :)
         !! Enthalpy carried by `fprec` (W/m^2).
      real(wp), allocatable :: heat_content_vprec(:, :)
         !! Enthalpy carried by `vprec` (W/m^2).
      real(wp), allocatable :: heat_content_lrunoff(:, :)
         !! Enthalpy carried by `lrunoff` (W/m^2).
      real(wp), allocatable :: heat_content_frunoff(:, :)
         !! Enthalpy carried by `frunoff` (W/m^2).
      real(wp), allocatable :: heat_content_seaice_melt(:, :)
         !! Enthalpy carried by `seaice_melt` (W/m^2).
      real(wp), allocatable :: heat_content_massin(:, :)
         !! **Assembler output — do NOT write.**  Sum of the six
         !! `heat_content_<flux>` companions above (W/m^2, >= 0 for warm
         !! inflow).  Filled by `ocean_surface_flux_assemble`.
      real(wp), allocatable :: heat_content_massout(:, :)
         !! **Assembler output — do NOT write.**  `SEAWATER_CP * T_sst *
         !! evap` (W/m^2, <= 0 since `evap <= 0`) — the enthalpy the ocean
         !! loses with evaporating mass, computed from the ocean's own
         !! surface temperature (the ocean, not a filler, owns this
         !! number).  There is deliberately NO `heat_content_evap` field
         !! — see the PR-12 plan §11.6.  Filled by
         !! `ocean_surface_flux_assemble`.

      real(wp), allocatable :: salt_flux(:, :)
         !! Net surface salt-flux COMPONENT (kg salt/m^2/s, **positive
         !! salinifies**) — a filler writes this (e.g. the ice brine
         !! coupler); the assembler adds `Q_salt_const` to produce
         !! `Q_salt`.  Virtual in v1 (no column-mass change).
      real(wp), allocatable :: salt_cavity(:, :)
         !! **Ice-shelf cavity basal-melt salt component**, same units
         !! and sign as `salt_flux` (positive salinifies;
         !! `&ocean_cavity_melt_nml`).  OWNED by `rdb_ocean_cavity_flux`
         !! and never written by the ice coupler, which full-overwrites
         !! `salt_flux`.
         !!
         !! **VIRTUAL salt flux.**  Melting adds freshwater MASS the
         !! Boussinesq column does not yet carry (Phase 3), so the
         !! dilution is emulated by removing salt:
         !!
         !!   `salt_cavity = -m_mass*(S_far - s_ice)`
         !!
         !! which is the exact fixed-mass equivalent of adding mass
         !! `m_mass` at salinity `s_ice` — see the derivation in
         !! `rdb_ocean_cavity_flux`'s module docstring.  Melting
         !! (`m_mass > 0`, `S_far > s_ice`) therefore FRESHENS.

      real(wp), allocatable :: p_surf_atm(:, :)
         !! **Input component.**  Atmospheric surface-pressure load
         !! (Pa, >= 0).  Filled by an external reader / configure-time
         !! scalar seed — the sea-ice path never writes this field.
         !! Ships zeroed with no consumer in this PR (the inverse-
         !! barometer PGF fold is a same-release-cycle follow-up).
      real(wp), allocatable :: p_surf(:, :)
         !! **Assembled total** (Pa, >= 0) — `p_surf_atm` plus any ice
         !! mass-loading term, **full overwrite, never `+=`** (a `+=`
         !! ratchets the load across outer steps with no bound).  No
         !! consumer in this PR; ships zeroed alongside `p_surf_atm` so
         !! the follow-up PGF fold needs no further plumbing.
   contains
      procedure, non_overridable :: init => ocean_surfflux_init
      procedure, non_overridable :: destroy => ocean_surfflux_destroy
      procedure, non_overridable :: enter_data => ocean_surfflux_enter_data
      procedure, non_overridable :: exit_data => ocean_surfflux_exit_data
      procedure, non_overridable :: set_surface_flux_const => ocean_surfflux_set_const
      procedure, non_overridable :: set_sw_penetration => ocean_surfflux_set_sw
      procedure, non_overridable :: set_restore => ocean_surfflux_set_restore
      procedure, non_overridable :: set_components => ocean_surfflux_set_components
      procedure, non_overridable :: set_p_surf_const => ocean_surfflux_set_p_surf_const
      procedure, non_overridable :: bytes => ocean_surface_flux_bytes
   end type ocean_surface_flux_t

contains

   subroutine ocean_surfflux_init(this, grid)
      class(ocean_surface_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (this%Q_heat(nx, ny), source=0.0_wp)
      allocate (this%Q_salt(nx, ny), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_surfflux_init

   subroutine ocean_surfflux_destroy(this)
      class(ocean_surface_flux_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%Q_heat)) deallocate (this%Q_heat)
      if (allocated(this%Q_salt)) deallocate (this%Q_salt)
      call ocean_surfflux_dealloc_components(this)
   end subroutine ocean_surfflux_destroy

   subroutine ocean_surfflux_dealloc_components(this)
      !! Deallocate the component set (no-op on an already-unallocated
      !! slot — every `deallocate` is `if (allocated(...))`-guarded).
      !! Shared by `destroy` and by `set_components` re-entry.
      class(ocean_surface_flux_t), intent(inout) :: this
      if (allocated(this%q_sw)) deallocate (this%q_sw)
      if (allocated(this%q_lw)) deallocate (this%q_lw)
      if (allocated(this%q_lat)) deallocate (this%q_lat)
      if (allocated(this%q_sens)) deallocate (this%q_sens)
      if (allocated(this%heat_added)) deallocate (this%heat_added)
      if (allocated(this%heat_cavity)) deallocate (this%heat_cavity)
      if (allocated(this%evap)) deallocate (this%evap)
      if (allocated(this%lprec)) deallocate (this%lprec)
      if (allocated(this%fprec)) deallocate (this%fprec)
      if (allocated(this%vprec)) deallocate (this%vprec)
      if (allocated(this%lrunoff)) deallocate (this%lrunoff)
      if (allocated(this%frunoff)) deallocate (this%frunoff)
      if (allocated(this%seaice_melt)) deallocate (this%seaice_melt)
      if (allocated(this%heat_content_lprec)) deallocate (this%heat_content_lprec)
      if (allocated(this%heat_content_fprec)) deallocate (this%heat_content_fprec)
      if (allocated(this%heat_content_vprec)) deallocate (this%heat_content_vprec)
      if (allocated(this%heat_content_lrunoff)) deallocate (this%heat_content_lrunoff)
      if (allocated(this%heat_content_frunoff)) deallocate (this%heat_content_frunoff)
      if (allocated(this%heat_content_seaice_melt)) deallocate (this%heat_content_seaice_melt)
      if (allocated(this%heat_content_massin)) deallocate (this%heat_content_massin)
      if (allocated(this%heat_content_massout)) deallocate (this%heat_content_massout)
      if (allocated(this%salt_flux)) deallocate (this%salt_flux)
      if (allocated(this%salt_cavity)) deallocate (this%salt_cavity)
      if (allocated(this%p_surf_atm)) deallocate (this%p_surf_atm)
      if (allocated(this%p_surf)) deallocate (this%p_surf)
   end subroutine ocean_surfflux_dealloc_components

   subroutine ocean_surfflux_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so the
      !! device-attach map base is the heap object, not a polymorphic stack
      !! box (AMD libomptarget cross-slot-overlap fix).
      class(ocean_surface_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_surface_flux_t)
         call ocean_surfflux_enter_data_impl(this)
      end select
   end subroutine ocean_surfflux_enter_data

   subroutine ocean_surfflux_enter_data_impl(this)
      type(ocean_surface_flux_t), intent(inout) :: this
      !$acc enter data copyin(this%Q_heat, this%Q_salt)
      !$acc update device(this%Q_heat, this%Q_salt)
      if (this%use_components) then
         !$acc enter data copyin(this%q_sw, this%q_lw, this%q_lat, this%q_sens, &
         !$acc&                  this%heat_added, this%heat_cavity, this%evap, &
         !$acc&                  this%lprec, this%fprec, &
         !$acc&                  this%vprec, this%lrunoff, this%frunoff, this%seaice_melt, &
         !$acc&                  this%heat_content_lprec, this%heat_content_fprec, &
         !$acc&                  this%heat_content_vprec, this%heat_content_lrunoff, &
         !$acc&                  this%heat_content_frunoff, this%heat_content_seaice_melt, &
         !$acc&                  this%heat_content_massin, this%heat_content_massout, &
         !$acc&                  this%salt_flux, this%salt_cavity, &
         !$acc&                  this%p_surf_atm, this%p_surf)
         !$acc update device(this%q_sw, this%q_lw, this%q_lat, this%q_sens, &
         !$acc&               this%heat_added, this%heat_cavity, this%evap, &
         !$acc&               this%lprec, this%fprec, &
         !$acc&               this%vprec, this%lrunoff, this%frunoff, this%seaice_melt, &
         !$acc&               this%heat_content_lprec, this%heat_content_fprec, &
         !$acc&               this%heat_content_vprec, this%heat_content_lrunoff, &
         !$acc&               this%heat_content_frunoff, this%heat_content_seaice_melt, &
         !$acc&               this%heat_content_massin, this%heat_content_massout, &
         !$acc&               this%salt_flux, this%salt_cavity, &
         !$acc&               this%p_surf_atm, this%p_surf)
      end if
   end subroutine ocean_surfflux_enter_data_impl

   subroutine ocean_surfflux_exit_data(this)
      class(ocean_surface_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_surface_flux_t)
         call ocean_surfflux_exit_data_impl(this)
      end select
   end subroutine ocean_surfflux_exit_data

   subroutine ocean_surfflux_exit_data_impl(this)
      type(ocean_surface_flux_t), intent(inout) :: this
      if (this%use_components) then
         !$acc exit data delete(this%q_sw, this%q_lw, this%q_lat, this%q_sens, &
         !$acc&                 this%heat_added, this%heat_cavity, this%evap, &
         !$acc&                 this%lprec, this%fprec, &
         !$acc&                 this%vprec, this%lrunoff, this%frunoff, this%seaice_melt, &
         !$acc&                 this%heat_content_lprec, this%heat_content_fprec, &
         !$acc&                 this%heat_content_vprec, this%heat_content_lrunoff, &
         !$acc&                 this%heat_content_frunoff, this%heat_content_seaice_melt, &
         !$acc&                 this%heat_content_massin, this%heat_content_massout, &
         !$acc&                 this%salt_flux, this%salt_cavity, &
         !$acc&                 this%p_surf_atm, this%p_surf)
      end if
      !$acc exit data delete(this%Q_heat, this%Q_salt)
   end subroutine ocean_surfflux_exit_data_impl

   subroutine ocean_surfflux_set_const(this, q_heat_val, q_salt_val)
      !! Fill `Q_heat` / `Q_salt` uniformly from scalar values and set
      !! the `has_heat` / `has_salt` flags so the apply-tracers kernel
      !! fires.  Mirrors `set_wind_stress_const` on the stress side.
      !! Host only — call `enter_data` afterwards (or `!$acc update
      !! device` if already mapped) to sync to the GPU.
      class(ocean_surface_flux_t), intent(inout) :: this
      real(wp), intent(in) :: q_heat_val, q_salt_val
      this%Q_heat_const = q_heat_val
      this%Q_salt_const = q_salt_val
      this%Q_heat = q_heat_val
      this%Q_salt = q_salt_val
      this%has_heat = (q_heat_val /= 0.0_wp)
      this%has_salt = (q_salt_val /= 0.0_wp)
   end subroutine ocean_surfflux_set_const

   subroutine ocean_surfflux_set_sw(this, sw_pen_frac, sw_band_ratio, &
                                    sw_zeta1, sw_zeta2, sw_source)
      !! Seed the shortwave-penetration band parameters and set the
      !! `has_sw` gate (`sw_pen_frac /= 0`).  Sibling to
      !! `set_surface_flux_const` — kept separate so existing callers of
      !! the heat/salt setter are unchanged.  Host only; the scalars are
      !! read host-side by the apply kernel (they parameterise the
      !! by-value arguments passed into the device `_impl`), so no extra
      !! device sync is needed beyond the existing `copyin(this)`.
      !!
      !! `sw_source` (optional; default `"net_heat"`) selects the
      !! irradiance source: `"net_heat"` ⇒ `I0 = sw_pen_frac*Q_heat`
      !! (legacy, bit-identical); `"q_sw"` ⇒ `I0 = sw_pen_frac*q_sw`
      !! (requires the PR-12 component set — the caller / `validate_config`
      !! guards allocation).  Sets the host-side `sw_from_qsw` gate.
      class(ocean_surface_flux_t), intent(inout) :: this
      real(wp), intent(in) :: sw_pen_frac, sw_band_ratio, sw_zeta1, sw_zeta2
      character(len=*), intent(in), optional :: sw_source
      this%sw_pen_frac = sw_pen_frac
      this%sw_band_ratio = sw_band_ratio
      this%sw_zeta1 = sw_zeta1
      this%sw_zeta2 = sw_zeta2
      this%has_sw = (sw_pen_frac /= 0.0_wp)
      this%sw_from_qsw = .false.
      if (present(sw_source)) this%sw_from_qsw = (trim(sw_source) == "q_sw")
   end subroutine ocean_surfflux_set_sw

   pure function sw_source_is_implemented(name) result(ok)
      !! Fail-loud predicate for the `&ocean_thermo_nml sw_source`
      !! selector — `validate_config` aborts on any string this rejects.
      !! The two recognised sources are the net-heat legacy path and the
      !! PR-12 `q_sw` component.
      character(len=*), intent(in) :: name
      logical :: ok
      ok = (trim(name) == "net_heat" .or. trim(name) == "q_sw")
   end function sw_source_is_implemented

   subroutine ocean_surfflux_set_restore(this, enable_T, enable_S, &
                                         piston_t_day, piston_s_day, &
                                         T_target, S_target)
      !! Seed the surface buoyancy restoring (MOM6 `RESTOREBUOY`)
      !! parameters and set the `has_restore_T` / `has_restore_S` gates.
      !! Sibling to `set_surface_flux_const` / `set_sw_penetration` —
      !! kept separate so existing callers are unchanged.  The piston
      !! velocities arrive in **m/day** (the MOM6 `FLUXCONST_*` unit) and
      !! are converted to MKS m/s here.  Effective-enable guard:
      !! `has_restore_* = enable_* .and. piston /= 0`, so an enabled
      !! switch with a zero piston is a silent no-op (rather than
      !! restoring everything toward the 0-degC / 0-PSU default target).
      !! Host only; all knobs are read host-side as by-value arguments to
      !! the device `_impl`, so no extra device sync beyond `copyin(this)`.
      class(ocean_surface_flux_t), intent(inout) :: this
      logical, intent(in) :: enable_T, enable_S
      real(wp), intent(in) :: piston_t_day, piston_s_day
      real(wp), intent(in) :: T_target, S_target
      this%restore_piston_T = piston_t_day/SECONDS_PER_DAY
      this%restore_piston_S = piston_s_day/SECONDS_PER_DAY
      this%restore_T_target = T_target
      this%restore_S_target = S_target
      this%has_restore_T = (enable_T .and. this%restore_piston_T /= 0.0_wp)
      this%has_restore_S = (enable_S .and. this%restore_piston_S /= 0.0_wp)
   end subroutine ocean_surfflux_set_restore

   subroutine ocean_surfflux_set_p_surf_const(this, p_surf_val)
      !! Seed the atmospheric surface-pressure INPUT component
      !! `p_surf_atm` (Pa) uniformly from a scalar namelist value (PR-17
      !! `&ocean_psurf_nml p_surf_const`).  Full overwrite of the pristine
      !! atmospheric base; the assembled total `p_surf` is built from it
      !! once per outer step in `p_surf_update_seam`.  No-op when the
      !! component set is not allocated (`use_components=.false.`) — the
      !! `&ocean_psurf_nml enable` guard in `validate_config` already
      !! requires `enable_components=.true.`, so a live consumer never hits
      !! the no-op.  Host only — call `enter_data` afterwards (or `!$acc
      !! update device` if already mapped) to sync to the GPU.
      !!
      !! Also seeds the assembled total `p_surf` to the same value — with no
      !! ice mass-loading `p_surf == p_surf_atm` and the atmospheric base is
      !! static, so this configure-time seed IS the assembly (PR-17 reads
      !! `p_surf` read-only in the dyn step).  PR-18's ice path overwrites
      !! `p_surf` per outer step via its own `inout` access.
      class(ocean_surface_flux_t), intent(inout) :: this
      real(wp), intent(in) :: p_surf_val
      if (.not. allocated(this%p_surf_atm)) return
      this%p_surf_atm = p_surf_val
      this%p_surf = p_surf_val
   end subroutine ocean_surfflux_set_p_surf_const

   subroutine ocean_surfflux_set_components(this, grid, enable)
      !! Configure-time gate for the PR-12 component set
      !! (`&ocean_forcing_nml enable_components`).  `init` runs before the
      !! namelist gate is known, so allocation happens HERE rather than in
      !! `init`: `enable = .false.` (default) leaves `use_components`
      !! false and allocates nothing — bit-identical, zero extra device
      !! memory.  `enable = .true.` allocates the full component set
      !! (`source=0.0_wp`) and flips the gate so
      !! `ocean_surface_flux_assemble` stops early-returning.  Must be
      !! called BEFORE `enter_data` (`rdb_ocean_state.F90`'s orchestrator)
      !! so the freshly-allocated arrays get mapped.  Re-entrant: calling
      !! again with a different `enable` deallocates first.
      class(ocean_surface_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: enable
      call ocean_surfflux_dealloc_components(this)
      this%use_components = enable
      if (enable) call ocean_surfflux_alloc_components(this, grid)
   end subroutine ocean_surfflux_set_components

   subroutine ocean_surfflux_alloc_components(this, grid)
      !! Allocate the 22-field component set + the two `p_surf*` fields,
      !! all `source=0.0_wp`, shape `(nx_total, ny_total)`.  Private —
      !! called only from `set_components`.
      class(ocean_surface_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (this%q_sw(nx, ny), source=0.0_wp)
      allocate (this%q_lw(nx, ny), source=0.0_wp)
      allocate (this%q_lat(nx, ny), source=0.0_wp)
      allocate (this%q_sens(nx, ny), source=0.0_wp)
      allocate (this%heat_added(nx, ny), source=0.0_wp)
      allocate (this%heat_cavity(nx, ny), source=0.0_wp)
      allocate (this%evap(nx, ny), source=0.0_wp)
      allocate (this%lprec(nx, ny), source=0.0_wp)
      allocate (this%fprec(nx, ny), source=0.0_wp)
      allocate (this%vprec(nx, ny), source=0.0_wp)
      allocate (this%lrunoff(nx, ny), source=0.0_wp)
      allocate (this%frunoff(nx, ny), source=0.0_wp)
      allocate (this%seaice_melt(nx, ny), source=0.0_wp)
      allocate (this%heat_content_lprec(nx, ny), source=0.0_wp)
      allocate (this%heat_content_fprec(nx, ny), source=0.0_wp)
      allocate (this%heat_content_vprec(nx, ny), source=0.0_wp)
      allocate (this%heat_content_lrunoff(nx, ny), source=0.0_wp)
      allocate (this%heat_content_frunoff(nx, ny), source=0.0_wp)
      allocate (this%heat_content_seaice_melt(nx, ny), source=0.0_wp)
      allocate (this%heat_content_massin(nx, ny), source=0.0_wp)
      allocate (this%heat_content_massout(nx, ny), source=0.0_wp)
      allocate (this%salt_flux(nx, ny), source=0.0_wp)
      allocate (this%salt_cavity(nx, ny), source=0.0_wp)
      allocate (this%p_surf_atm(nx, ny), source=0.0_wp)
      allocate (this%p_surf(nx, ny), source=0.0_wp)
   end subroutine ocean_surfflux_alloc_components

   subroutine ocean_surface_flux_apply_tracers(grid, sf, ms, dt, active, wet_dyn)
      !! Add the surface heat + salt fluxes directly to the top
      !! tracer layer.  Operates in `hTr` space (concentration·
      !! thickness): for temperature
      !!   d(hT_top)/dt = Q_heat(i,j) / (rho_0 · cp)
      !! For salinity
      !!   d(hS_top)/dt = Q_salt(i,j) / rho_0
      !! (Both expressed in units that match the `hTr` convention:
      !! `hTr = T·h` so the forcing has units of T·h/s = K·m/s.
      !! Q_heat / (rho_0·cp) has units (W/m^2)/(kg/m^3·J/kg/K) =
      !! K·m/s ✓.)
      !!
      !! Reads the 2D `Q_heat(:,:)` / `Q_salt(:,:)` fields per column
      !! (seeded uniformly from the scalar knobs by default; overwritten
      !! pointwise by Area-A3 data-override or Area-A4 restoring).
      !! No-op when `sf%has_heat` / `sf%has_salt` are false (set by
      !! `set_surface_flux_const` when the fill value is non-zero; any
      !! field-override path must set them before calling), or when no
      !! temperature / salinity tracer is registered.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in), optional :: sf
         !! Optional — when absent the kernel is a no-op (no surface
         !! forcing configured).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional gate (thermo cadence).  Absent ⇒ kernel runs;
         !! present-and-false ⇒ early return.
      real(wp), intent(in), optional :: wet_dyn(:, :)
         !! Optional DYNAMIC cell wet mask (wet/dry,
         !! docs/ocean_wetdry_plan.md §4.4) composed multiplicatively
         !! with the static `ms%wet_mask` — surface fluxes must not
         !! enter a dynamically dry column (heating a mm-scale residual
         !! sliver blows its temperature up).  Absent ⇒ the original
         !! static-mask path, byte-identical.

      integer :: nx, ny, nz, idx_T, idx_S, idx_ps

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. present(sf)) return
      ! Guard: has_heat / has_salt flags unset → skip (set by
      ! set_surface_flux_const or any field-override path).
      if (.not. sf%has_heat .and. .not. sf%has_salt) return
      if (.not. allocated(ms%tracers)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      idx_ps = ms%idx_pseudo_salt

      ! Cell-centred wet mask: heat / salt flux only enters ocean cells.
      ! Land cells (mask = 0) accumulate nothing — including in the
      ! budget contributor, so the conservation residual stays clean.
      ! The shim+_impl split keeps `tracers(idx)%hTr` deref on the host
      ! (array-of-DT registry indirection blocks NVHPC device codegen).
      ! The division Q(i,j)/(rho0*cp) [or Q(i,j)/rho0] is applied as
      ! the inv_scale multiplier inside the kernel.
      if (present(wet_dyn)) then
         if (idx_T > 0 .and. sf%has_heat) then
            call apply_surface_src_2d_dyn_impl(ms%tracers(idx_T)%hTr, &
                                               ms%heat_budget_surface, &
                                               ms%wet_mask, wet_dyn, sf%Q_heat, &
                                               dt/(sf%rho0*sf%cp), nz, nx, ny)
         end if
         if (idx_S > 0 .and. sf%has_salt) then
            call apply_surface_src_2d_dyn_impl(ms%tracers(idx_S)%hTr, &
                                               ms%salt_budget_surface, &
                                               ms%wet_mask, wet_dyn, sf%Q_salt, &
                                               dt/sf%rho0, nz, nx, ny)
         end if
         ! Pseudo-salt mirror: exactly salinity's surface salt flux,
         ! but through the NOBUDGET twin — budget_id = NONE so it must
         ! not add into salt_budget_surface (§5.6 / test
         ! pseudo_salt_no_budget_contribution).
         if (idx_ps > 0 .and. sf%has_salt) then
            call apply_surface_src_2d_dyn_nobudget_impl(ms%tracers(idx_ps)%hTr, &
                                                        ms%wet_mask, wet_dyn, sf%Q_salt, &
                                                        dt/sf%rho0, nz, nx, ny)
         end if
         return
      end if
      if (idx_T > 0 .and. sf%has_heat) then
         call apply_surface_src_2d_impl(ms%tracers(idx_T)%hTr, &
                                        ms%heat_budget_surface, &
                                        ms%wet_mask, sf%Q_heat, &
                                        dt/(sf%rho0*sf%cp), nz, nx, ny)
      end if
      if (idx_S > 0 .and. sf%has_salt) then
         call apply_surface_src_2d_impl(ms%tracers(idx_S)%hTr, &
                                        ms%salt_budget_surface, &
                                        ms%wet_mask, sf%Q_salt, &
                                        dt/sf%rho0, nz, nx, ny)
      end if
      ! Pseudo-salt mirror (plain, static-mask path) — see note above.
      if (idx_ps > 0 .and. sf%has_salt) then
         call apply_surface_src_2d_nobudget_impl(ms%tracers(idx_ps)%hTr, &
                                                 ms%wet_mask, sf%Q_salt, &
                                                 dt/sf%rho0, nz, nx, ny)
      end if
   end subroutine ocean_surface_flux_apply_tracers

   pure subroutine apply_surface_src_2d_impl(hTr, budget, wet_mask, Q_field, &
                                             inv_scale, nz, nx, ny)
      !! Stamp `inv_scale · Q_field(i,j) · wet_mask(i,j)` onto the top
      !! layer (k = nz) of a tracer's hTr array, mirror into the matching
      !! budget contributor.  Explicit-shape dummies so NVHPC stdpar can
      !! compile device kernels against static bounds.
      !!
      !! `inv_scale` = dt/(rho0·cp) for heat, dt/rho0 for salt — a
      !! column-invariant multiplier that the caller derives from `sf`.
      !! `Q_field` carries any (i,j) spatial variation; for the default
      !! constant-fill case it is uniform, giving arithmetic identical to
      !! the old scalar-broadcast path.
      !!
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: Q_field(nx, ny)
      real(wp), intent(in)    :: inv_scale
      integer :: i, j
      real(wp) :: cell
      do concurrent(j=1:ny, i=1:nx) local(cell)
         cell = inv_scale*Q_field(i, j)*wet_mask(i, j)
         hTr(i, j, nz) = hTr(i, j, nz) + cell
         budget(i, j, nz) = budget(i, j, nz) + cell
      end do
   end subroutine apply_surface_src_2d_impl

   pure subroutine apply_surface_src_2d_dyn_impl(hTr, budget, wet_mask, wet_dyn, &
                                                 Q_field, inv_scale, nz, nx, ny)
      !! Wet/dry variant of `apply_surface_src_2d_impl`: the DYNAMIC cell
      !! wet mask composes multiplicatively with the static one, so a
      !! dynamically dry column (total depth below `&ocean_wetdry_nml
      !! dry_depth`) receives NO surface flux — heating a mm-scale
      !! residual sliver would blow its temperature up
      !! (docs/ocean_wetdry_plan.md §4.4).  Separate _impl (not an
      !! in-loop optional test): the knob-off path keeps the original
      !! kernel untouched, byte-identical.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: wet_dyn(nx, ny)
      real(wp), intent(in)    :: Q_field(nx, ny)
      real(wp), intent(in)    :: inv_scale
      integer :: i, j
      real(wp) :: cell
      do concurrent(j=1:ny, i=1:nx) local(cell)
         cell = inv_scale*Q_field(i, j)*wet_mask(i, j)*wet_dyn(i, j)
         hTr(i, j, nz) = hTr(i, j, nz) + cell
         budget(i, j, nz) = budget(i, j, nz) + cell
      end do
   end subroutine apply_surface_src_2d_dyn_impl

   pure subroutine apply_surface_src_2d_nobudget_impl(hTr, wet_mask, Q_field, &
                                                      inv_scale, nz, nx, ny)
      !! Byte-for-byte copy of `apply_surface_src_2d_impl` with the
      !! `budget` dummy and its accumulation line removed — the
      !! pseudo-salt mirror of salinity's surface flux, which by
      !! contract (`budget_id = TRACER_BUDGET_NONE`) must not touch
      !! `salt_budget_surface`.  Separate `_impl`, not an in-loop
      !! `present(budget)` test (house idiom, see
      !! `apply_surface_src_2d_dyn_impl`'s docstring) — this keeps the
      !! production S/T impl untouched and the increment `hTr` receives
      !! bit-identical to salinity's.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: Q_field(nx, ny)
      real(wp), intent(in)    :: inv_scale
      integer :: i, j
      real(wp) :: cell
      do concurrent(j=1:ny, i=1:nx) local(cell)
         cell = inv_scale*Q_field(i, j)*wet_mask(i, j)
         hTr(i, j, nz) = hTr(i, j, nz) + cell
      end do
   end subroutine apply_surface_src_2d_nobudget_impl

   pure subroutine apply_surface_src_2d_dyn_nobudget_impl(hTr, wet_mask, wet_dyn, &
                                                          Q_field, inv_scale, nz, nx, ny)
      !! Wet/dry NOBUDGET twin — see `apply_surface_src_2d_nobudget_impl`
      !! and `apply_surface_src_2d_dyn_impl`.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: wet_dyn(nx, ny)
      real(wp), intent(in)    :: Q_field(nx, ny)
      real(wp), intent(in)    :: inv_scale
      integer :: i, j
      real(wp) :: cell
      do concurrent(j=1:ny, i=1:nx) local(cell)
         cell = inv_scale*Q_field(i, j)*wet_mask(i, j)*wet_dyn(i, j)
         hTr(i, j, nz) = hTr(i, j, nz) + cell
      end do
   end subroutine apply_surface_src_2d_dyn_nobudget_impl

   subroutine ocean_surface_flux_apply_sw_penetration(grid, sf, ms, dt, active, cover_frac)
      !! Additive correction that redistributes the penetrating
      !! shortwave fraction of `Q_heat` through the upper water column
      !! as a two-band exponential (Paulson & Simpson 1977; Jerlov
      !! types), instead of leaving all of it deposited at the surface
      !! layer by `ocean_surface_flux_apply_tracers`.
      !!
      !! Penetrating irradiance at downward depth `d`:
      !!   I(d) = I0 · [ R·exp(-d/zeta1) + (1-R)·exp(-d/zeta2) ]
      !! with `I0 = sw_pen_frac · Q_heat(i,j)`.  Per-layer absorbed SW =
      !! `I(d_top) - I(d_bot)`; the bed (`k = 1`) is treated as opaque
      !! (`I_bot := 0`) so the column absorbs all of `I0` and energy is
      !! conserved exactly (Σ_k absorbed_k = I0).
      !!
      !! Additive-correction structure: the surface kernel already
      !! deposited the full `I0` lump at `k = nz`; this kernel removes it
      !! there (`-inv_scale · I0`) and adds the distributed profile, so
      !! the net column heat change versus the legacy all-at-`nz`
      !! deposition is zero — shortwave only MOVES heat in depth.  Both
      !! the `hTr` and `heat_budget_surface` increments are mirrored, so
      !! the total surface heat budget is unchanged (just depth-spread).
      !!
      !! Bottom-up convention (`k = nz` surface, `k = 1` bed).  No-op
      !! unless `sf%has_sw .and. sf%has_heat`, a temperature tracer is
      !! registered, and (optionally) `active` is true.  Default-off
      !! (`sw_pen_frac = 0` ⇒ `has_sw = .false.`) leaves the path
      !! byte-for-byte unchanged.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in), optional :: sf
         !! Optional — absent ⇒ no-op (no forcing configured).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional thermo-cadence gate.  Present-and-false ⇒ early
         !! return; absent ⇒ kernel runs.
      real(wp), intent(in), optional :: cover_frac(:, :)
         !! Optional ice-shelf cover fraction (`metrics%cover_frac`,
         !! v1 binary).  Present ⇒ the penetrating irradiance is scaled
         !! by `1 - cover_frac`, so a covered column absorbs NOTHING —
         !! no sunlight reaches the ocean through several hundred metres
         !! of ice.  Needed even on the `sw_source="net_heat"` branch,
         !! where `Q_heat` is already assembler-masked: this kernel's
         !! job is to MOVE a surface lump down the column, and on a
         !! masked column the lump it would remove was never deposited
         !! (the same argument `&ocean_wetdry_nml` makes for a dry
         !! column).  Absent ⇒ the original kernel, byte-identical.
      ! assumed-shape-ok: thermo-cadence shim, forwarded to an
      ! explicit-shape `_impl` before the device loop.

      integer :: nx, ny, nz, idx_T
      logical :: masked

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. present(sf)) return
      if (.not. (sf%has_sw .and. sf%has_heat)) return
      if (.not. allocated(ms%tracers)) return

      idx_T = ms%idx_temperature
      if (idx_T <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Shim+_impl split: keep the `tracers(idx_T)%hTr` registry deref on
      ! the host (array-of-DT indirection blocks NVHPC device codegen).
      ! The irradiance source is selected HOST-SIDE (`sw_from_qsw`): the
      ! PR-12 `q_sw` component is allocated only under `use_components`,
      ! so passing it as an actual argument is only legal on the branch
      ! guarded by the host flag (validate_config forces
      ! enable_components when sw_source="q_sw", making this total).  The
      ! additive-correction identity `-I0 + Σ_k I0·(T_top - T_bot) = 0`
      ! holds for ANY I0, so the source swap cannot break conservation.
      masked = .false.
      if (present(cover_frac)) then
         masked = (size(cover_frac, 1) == nx .and. size(cover_frac, 2) == ny)
      end if
      if (masked) then
         if (sf%sw_from_qsw) then
            call apply_sw_penetration_cover_impl(ms%tracers(idx_T)%hTr, &
                                                 ms%heat_budget_surface, &
                                                 ms%h_layer, ms%wet_mask, cover_frac, sf%q_sw, &
                                                 dt/(sf%rho0*sf%cp), sf%sw_pen_frac, &
                                                 sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                                 nz, nx, ny)
         else
            call apply_sw_penetration_cover_impl(ms%tracers(idx_T)%hTr, &
                                                 ms%heat_budget_surface, &
                                                 ms%h_layer, ms%wet_mask, cover_frac, sf%Q_heat, &
                                                 dt/(sf%rho0*sf%cp), sf%sw_pen_frac, &
                                                 sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                                 nz, nx, ny)
         end if
         return
      end if
      if (sf%sw_from_qsw) then
         call apply_sw_penetration_impl(ms%tracers(idx_T)%hTr, &
                                        ms%heat_budget_surface, &
                                        ms%h_layer, ms%wet_mask, sf%q_sw, &
                                        dt/(sf%rho0*sf%cp), sf%sw_pen_frac, &
                                        sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                        nz, nx, ny)
      else
         call apply_sw_penetration_impl(ms%tracers(idx_T)%hTr, &
                                        ms%heat_budget_surface, &
                                        ms%h_layer, ms%wet_mask, sf%Q_heat, &
                                        dt/(sf%rho0*sf%cp), sf%sw_pen_frac, &
                                        sf%sw_band_ratio, sf%sw_zeta1, sf%sw_zeta2, &
                                        nz, nx, ny)
      end if
   end subroutine ocean_surface_flux_apply_sw_penetration

   pure function sw_transmission(d, R, zeta1, zeta2) result(trans)
      !! Two-band (Paulson & Simpson 1977) normalised downward
      !! irradiance transmission at depth `d` below the free surface,
      !! `T(d) = R·exp(-d/zeta1) + (1-R)·exp(-d/zeta2)`, `T(0) = 1`.
      !! This is THE single shared definition — the SW deposition
      !! kernel, the KPP `MXL_SW`/`LV1_SW` boundary-layer correction, and
      !! the EPBL in-layer PE-cost ledger all consume it, so the three
      !! consumers cannot disagree about where the sunlight went.  Marked
      !! `!$acc routine seq` so it inlines into same-module device kernels
      !! and is callable from cross-module `do concurrent` kernels.
      !$acc routine seq
      real(wp), intent(in) :: d, R, zeta1, zeta2
      real(wp) :: trans
      trans = R*exp(-d/zeta1) + (1.0_wp - R)*exp(-d/zeta2)
   end function sw_transmission

   pure function sw_pe_cost_shape(tau) result(phi)
      !! In-layer potential-energy-cost shape function `Phi(tau)` for the
      !! EPBL TKE ledger, `tau = h/zeta` the in-layer optical depth of a
      !! single band.  It is the fraction of the pure-skin PE cost that
      !! homogenising an EXPONENTIALLY distributed in-layer heating
      !! actually incurs (Paulson & Simpson 1977 profile; the EPBL
      !! energetics of Reichl & Hallberg 2018):
      !!
      !!   Phi(tau) = [ tau·(1+e^-tau) - 2·(1-e^-tau) ] / [ tau·(1-e^-tau) ]
      !!
      !! Limits: `Phi(0) = 0` (heating already uniform through the layer
      !! ⇒ homogenising it costs nothing) and `Phi(∞) = 1` (all heating
      !! at the layer top ⇒ full skin cost, recovering Roundabout's existing
      !! `ctke_sfc` skin form exactly).  The closed form is `0/0` as
      !! `tau -> 0` and cancels catastrophically for small `tau`; a thin
      !! layer (`h << zeta2 = 23 m`) is the COMMON case, so the Taylor
      !! branch `Phi ≈ (tau/6)·(1 - tau²/60)` is mandatory below the
      !! `tau = 1e-2` seam.  `!$acc routine seq` for cross-module device
      !! calls (the EPBL prep sweep).
      !$acc routine seq
      real(wp), intent(in) :: tau
      real(wp) :: phi
      real(wp) :: em1
      real(wp), parameter :: TAU_TAYLOR = 1.0e-2_wp
      real(wp), parameter :: C1_6 = 1.0_wp/6.0_wp
      real(wp), parameter :: C1_60 = 1.0_wp/60.0_wp
      if (tau <= TAU_TAYLOR) then
         phi = C1_6*tau*(1.0_wp - C1_60*tau*tau)
      else
         em1 = 1.0_wp - exp(-tau)
         phi = (tau*(1.0_wp + exp(-tau)) - 2.0_wp*em1)/(tau*em1)
      end if
   end function sw_pe_cost_shape

   pure subroutine apply_sw_penetration_impl(hTr, budget, h_layer, wet_mask, &
                                             sw_src, inv_scale, sw_pen_frac, R, &
                                             zeta1, zeta2, nz, nx, ny)
      !! Per-column two-band shortwave redistribution.  Explicit-shape
      !! dummies so NVHPC stdpar compiles device kernels against static
      !! bounds.  Difference form (`I(d_top) - I(d_bot)`) — no division
      !! by `h`, so a vanishing layer (`h → 0 ⇒ d_top == d_bot`) absorbs
      !! zero automatically with no guard.  The transmission `T(d)` is
      !! the shared `sw_transmission` (same-module ⇒ inlined by NVHPC),
      !! so the deposition and the boundary-layer coupling cannot diverge.
      !!
      !! `inv_scale` = dt/(rho0·cp); `R`/`zeta1`/`zeta2` are the two-band
      !! parameters; `sw_pen_frac` scales `sw_src` to the penetrating
      !! irradiance `I0`.  `sw_src` is the caller-selected source
      !! (`Q_heat` for the legacy net-heat path, `q_sw` for the PR-12
      !! component path); the additive-correction conservation identity
      !! `-I0 + Σ_k I0·(T_top - T_bot) = 0` holds for ANY `sw_src`, so the
      !! source swap cannot break conservation.  At `k = nz` the legacy
      !! `I0` lump is subtracted before adding the surface band, keeping
      !! the net column change against the all-at-`nz` baseline at zero.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: sw_src(nx, ny)
      real(wp), intent(in)    :: inv_scale, sw_pen_frac, R, zeta1, zeta2
      integer :: i, j, k
      real(wp) :: i0col, d_top, d_bot, trans_top, trans_bot, absorbed, add

      do concurrent(j=1:ny, i=1:nx) local(i0col, d_top, d_bot, trans_top, &
                                          trans_bot, absorbed, add, k)
         i0col = sw_pen_frac*sw_src(i, j)*wet_mask(i, j)
         d_top = 0.0_wp
         do k = nz, 1, -1
            d_bot = d_top + h_layer(i, j, k)
            trans_top = sw_transmission(d_top, R, zeta1, zeta2)
            if (k > 1) then
               trans_bot = sw_transmission(d_bot, R, zeta1, zeta2)
            else
               trans_bot = 0.0_wp   ! bed opaque: column absorbs all of I0
            end if
            absorbed = i0col*(trans_top - trans_bot)
            add = inv_scale*absorbed
            if (k == nz) add = add - inv_scale*i0col   ! remove the surface lump
            hTr(i, j, k) = hTr(i, j, k) + add
            budget(i, j, k) = budget(i, j, k) + add
            d_top = d_bot
         end do
      end do
   end subroutine apply_sw_penetration_impl

   pure subroutine apply_sw_penetration_cover_impl(hTr, budget, h_layer, wet_mask, &
                                                   cover_frac, sw_src, inv_scale, &
                                                   sw_pen_frac, R, zeta1, zeta2, nz, nx, ny)
      !! Ice-shelf-cover twin of `apply_sw_penetration_impl`: the
      !! open-water factor `1 - cover_frac` composes multiplicatively
      !! with `wet_mask` into the column irradiance `I0`, so a fully
      !! covered column neither removes the surface lump nor deposits a
      !! profile — it is left EXACTLY untouched.  Separate `_impl`, not
      !! an in-loop `present()` test (house idiom, see
      !! `apply_surface_src_2d_dyn_impl`) — the cover-off path keeps the
      !! original kernel byte-identical.
      !!
      !! The additive-correction conservation identity
      !! `-I0 + Σ_k I0·(T_top - T_bot) = 0` holds for ANY `I0`, and
      !! `I0 = 0` is the degenerate case of it, so scaling the source
      !! cannot break column heat conservation.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: cover_frac(nx, ny)
      real(wp), intent(in)    :: sw_src(nx, ny)
      real(wp), intent(in)    :: inv_scale, sw_pen_frac, R, zeta1, zeta2
      integer :: i, j, k
      real(wp) :: i0col, d_top, d_bot, trans_top, trans_bot, absorbed, add

      do concurrent(j=1:ny, i=1:nx) local(i0col, d_top, d_bot, trans_top, &
                                          trans_bot, absorbed, add, k)
         i0col = sw_pen_frac*sw_src(i, j)*wet_mask(i, j)*(1.0_wp - cover_frac(i, j))
         d_top = 0.0_wp
         do k = nz, 1, -1
            d_bot = d_top + h_layer(i, j, k)
            trans_top = sw_transmission(d_top, R, zeta1, zeta2)
            if (k > 1) then
               trans_bot = sw_transmission(d_bot, R, zeta1, zeta2)
            else
               trans_bot = 0.0_wp   ! bed opaque: column absorbs all of I0
            end if
            absorbed = i0col*(trans_top - trans_bot)
            add = inv_scale*absorbed
            if (k == nz) add = add - inv_scale*i0col   ! remove the surface lump
            hTr(i, j, k) = hTr(i, j, k) + add
            budget(i, j, k) = budget(i, j, k) + add
            d_top = d_bot
         end do
      end do
   end subroutine apply_sw_penetration_cover_impl

   subroutine ocean_surface_restore_apply_tracers(grid, sf, ms, dt, active, cover_frac)
      !! Surface buoyancy restoring (MOM6 `RESTOREBUOY`): relax the
      !! top-layer (`k = nz`) temperature / salinity toward scalar
      !! targets with a piston velocity `p` [m/s].  Unlike
      !! `ocean_surface_flux_apply_tracers` (which reads a pre-filled
      !! static `Q_*` field), the restoring flux is DYNAMIC — it depends
      !! on the live SST / SSS each thermo step — so it is computed
      !! in-kernel from `(target - surface_concentration)` rather than a
      !! stored field.  This keeps the const-flux path byte-for-byte
      !! unchanged and adds no second device sync of `Q_*`.
      !!
      !! Bottom-up convention: surface = `k = nz`, bed = `k = 1`.  The
      !! surface concentration is `hTr(i,j,nz) / max(h_layer(i,j,nz),
      !! h_min)`.  Per thermo step, in hTr-space:
      !!   d(hT_top) = dt · p_T · (T_target - SST) · wet_mask   [K·m]
      !!   d(hS_top) = dt · p_S · (S_target - SSS) · wet_mask   [PSU·m]
      !! (The `rho0·cp` of the equivalent W/m^2 restoring heat flux
      !! cancels against the `dt/(rho0·cp)` apply scaling — see the spec.)
      !! The same increment is mirrored into `heat_budget_surface(:,:,nz)`
      !! / `salt_budget_surface(:,:,nz)` so the surface-budget diagnostics
      !! see restoring as an explicit (non-conservative) source.
      !!
      !! Restoring is a relaxation forcing, NOT a conservative process:
      !! it deliberately injects / removes heat + salt to nudge the
      !! surface.  The budget contributors account for it so the
      !! conservation diagnostics do not flag it as a leak.
      !!
      !! MOM6-fidelity note (salt path): the temperature branch is an
      !! exact analogue of MOM6 `RESTOREBUOY` `heat_added` (the `rho0*cp`
      !! cancellation above makes the tendency identical for a given
      !! `FLUXCONST_T`).  MOM6's SALT branch is instead a virtual
      !! freshwater flux `vprec = -rho0*p_S*(S*-SSS)/(0.5*(SSS+S*))` that
      !! changes the column mass and dilutes a CONSERVED salt content.
      !! We use a linearised salt-CONTENT injection (no thickness change),
      !! which reproduces MOM6's surface-salinity tendency TO FIRST ORDER
      !! — the dilution factor `SSS/(0.5*(SSS+S*)) ~ 1` for realistic
      !! anomalies — but is non-conservative and drops that second-order
      !! factor.  Exact `vprec` parity is a documented v2 follow-up.
      !!
      !! No-op unless `sf%has_restore_T .or. sf%has_restore_S` (set by
      !! `set_restore` when the matching switch is on AND the piston is
      !! non-zero), a matching tracer is registered, and (optionally)
      !! `active` is true.  Default-off leaves the surface-flux path
      !! byte-for-byte unchanged.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in), optional :: sf
         !! Optional — absent ⇒ no-op (no forcing configured).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional thermo-cadence gate.  Present-and-false ⇒ early
         !! return; absent ⇒ kernel runs.
      real(wp), intent(in), optional :: cover_frac(:, :)
         !! Optional ice-shelf cover fraction (`metrics%cover_frac`,
         !! v1 binary).  Present ⇒ the restoring increment is scaled by
         !! `1 - cover_frac`, so a covered column is NOT relaxed toward
         !! an atmospheric target — under a shelf the surface is a
         !! melting ice interface, and restoring there would overwhelm
         !! the melt signal with a number the atmosphere never set.
         !! Absent ⇒ the original kernel, byte-identical.
      ! assumed-shape-ok: thermo-cadence shim, forwarded to an
      ! explicit-shape `_impl` before the device loop.

      integer :: nx, ny, nz, idx_T, idx_S
      logical :: masked

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. present(sf)) return
      if (.not. sf%has_restore_T .and. .not. sf%has_restore_S) return
      if (.not. allocated(ms%tracers)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity

      masked = .false.
      if (present(cover_frac)) then
         masked = (size(cover_frac, 1) == nx .and. size(cover_frac, 2) == ny)
      end if
      if (masked) then
         if (idx_T > 0 .and. sf%has_restore_T) then
            call apply_surface_restore_2d_cover_impl(ms%tracers(idx_T)%hTr, &
                                                     ms%heat_budget_surface, &
                                                     ms%h_layer, ms%wet_mask, cover_frac, &
                                                     dt*sf%restore_piston_T, &
                                                     sf%restore_T_target, sf%h_min, &
                                                     nz, nx, ny)
         end if
         if (idx_S > 0 .and. sf%has_restore_S) then
            call apply_surface_restore_2d_cover_impl(ms%tracers(idx_S)%hTr, &
                                                     ms%salt_budget_surface, &
                                                     ms%h_layer, ms%wet_mask, cover_frac, &
                                                     dt*sf%restore_piston_S, &
                                                     sf%restore_S_target, sf%h_min, &
                                                     nz, nx, ny)
         end if
         return
      end if

      ! Shim+_impl split: keep the `tracers(idx)%hTr` registry deref on
      ! the host (array-of-DT indirection blocks NVHPC device codegen).
      ! Scalar targets / piston / h_min pass by value into the kernel —
      ! no per-cell field, so no extra device array.
      if (idx_T > 0 .and. sf%has_restore_T) then
         call apply_surface_restore_2d_impl(ms%tracers(idx_T)%hTr, &
                                            ms%heat_budget_surface, &
                                            ms%h_layer, ms%wet_mask, &
                                            dt*sf%restore_piston_T, &
                                            sf%restore_T_target, sf%h_min, &
                                            nz, nx, ny)
      end if
      if (idx_S > 0 .and. sf%has_restore_S) then
         call apply_surface_restore_2d_impl(ms%tracers(idx_S)%hTr, &
                                            ms%salt_budget_surface, &
                                            ms%h_layer, ms%wet_mask, &
                                            dt*sf%restore_piston_S, &
                                            sf%restore_S_target, sf%h_min, &
                                            nz, nx, ny)
      end if
   end subroutine ocean_surface_restore_apply_tracers

   pure subroutine apply_surface_restore_2d_impl(hTr, budget, h_layer, &
                                                 wet_mask, dt_piston, tgt, &
                                                 h_min, nz, nx, ny)
      !! Stamp the per-step restoring increment `dt·p·(target - surf)·
      !! wet_mask` onto the top layer (k = nz) of a tracer's hTr array
      !! and mirror it into the matching budget contributor.
      !! Explicit-shape dummies so NVHPC stdpar compiles device kernels
      !! against static bounds.
      !!
      !! `dt_piston` = dt·piston [m] is the column-invariant multiplier
      !! the caller derives from `sf` (the `rho0·cp` cancels — see the
      !! caller doc).  The surface concentration is recovered as
      !! `hTr(nz) / max(h_layer(nz), h_min)`; the `h_min` floor keeps the
      !! relaxation finite under a pinched top layer (rate `p/h`
      !! saturates rather than diverges).
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: dt_piston, tgt, h_min
      integer :: i, j
      real(wp) :: surf, inc
      do concurrent(j=1:ny, i=1:nx) local(surf, inc)
         surf = hTr(i, j, nz)/max(h_layer(i, j, nz), h_min)
         inc = dt_piston*(tgt - surf)*wet_mask(i, j)
         hTr(i, j, nz) = hTr(i, j, nz) + inc
         budget(i, j, nz) = budget(i, j, nz) + inc
      end do
   end subroutine apply_surface_restore_2d_impl

   pure subroutine apply_surface_restore_2d_cover_impl(hTr, budget, h_layer, &
                                                       wet_mask, cover_frac, dt_piston, &
                                                       tgt, h_min, nz, nx, ny)
      !! Ice-shelf-cover twin of `apply_surface_restore_2d_impl`: the
      !! open-water factor `1 - cover_frac` composes multiplicatively
      !! with `wet_mask`, so a covered column receives EXACTLY zero
      !! restoring — and, because the same factor multiplies the budget
      !! mirror, exactly zero restoring shows up in the heat/salt
      !! surface budget there too.  Separate `_impl`, not an in-loop
      !! `present()` test (house idiom).
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in)    :: h_layer(nx, ny, nz)
      real(wp), intent(in)    :: wet_mask(nx, ny)
      real(wp), intent(in)    :: cover_frac(nx, ny)
      real(wp), intent(in)    :: dt_piston, tgt, h_min
      integer :: i, j
      real(wp) :: surf, inc
      do concurrent(j=1:ny, i=1:nx) local(surf, inc)
         surf = hTr(i, j, nz)/max(h_layer(i, j, nz), h_min)
         inc = dt_piston*(tgt - surf)*wet_mask(i, j)*(1.0_wp - cover_frac(i, j))
         hTr(i, j, nz) = hTr(i, j, nz) + inc
         budget(i, j, nz) = budget(i, j, nz) + inc
      end do
   end subroutine apply_surface_restore_2d_cover_impl

   pure function ocean_surface_flux_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the surface flux slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_surface_flux_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%Q_heat) &
               + arr_bytes(this%Q_salt) &
               + arr_bytes(this%q_sw) + arr_bytes(this%q_lw) &
               + arr_bytes(this%q_lat) + arr_bytes(this%q_sens) &
               + arr_bytes(this%heat_added) &
               + arr_bytes(this%heat_cavity) &
               + arr_bytes(this%evap) + arr_bytes(this%lprec) &
               + arr_bytes(this%fprec) + arr_bytes(this%vprec) &
               + arr_bytes(this%lrunoff) + arr_bytes(this%frunoff) &
               + arr_bytes(this%seaice_melt) &
               + arr_bytes(this%heat_content_lprec) &
               + arr_bytes(this%heat_content_fprec) &
               + arr_bytes(this%heat_content_vprec) &
               + arr_bytes(this%heat_content_lrunoff) &
               + arr_bytes(this%heat_content_frunoff) &
               + arr_bytes(this%heat_content_seaice_melt) &
               + arr_bytes(this%heat_content_massin) &
               + arr_bytes(this%heat_content_massout) &
               + arr_bytes(this%salt_flux) &
               + arr_bytes(this%salt_cavity) &
               + arr_bytes(this%p_surf_atm) + arr_bytes(this%p_surf)
   end function ocean_surface_flux_bytes

   pure subroutine ocean_surface_flux_apply_cover_const(sf, cover_frac)
      !! Mask the STATIC scalar `Q_heat` / `Q_salt` fill with the
      !! ice-shelf cover, for the `use_components = .false.` path only.
      !!
      !! With the component set on, `Q_heat`/`Q_salt` are assembler
      !! outputs and the mask belongs there (`ocean_surface_flux_assemble`'s
      !! `cover_frac`); writing them here would be a second writer on an
      !! assembler-owned slot, which the fill contract forbids.  With the
      !! component set OFF there is no assembler, `Q_heat`/`Q_salt` are
      !! the configure-time `set_surface_flux_const` fill and nothing
      !! rewrites them per step — so masking them once, at configure
      !! after the cover is built, is the whole job.  Without this, a
      !! geometry-only cavity run (`&ocean_cavity_dyn_nml` with no melt)
      !! would still push a uniform `&ocean_thermo_nml q_heat` through
      !! the ice.
      !!
      !! No-op when the component set is on, when the fields are
      !! unallocated, or when `cover_frac` is the `(1,1)` placeholder.
      !! Idempotent (multiplies by 0 or 1).  Host-side at configure —
      !! call BEFORE `enter_data` (or follow with an `!$acc update
      !! device`).
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: cover_frac(:, :)
         !! Ice-cover fraction at cell centres (`metrics%cover_frac`).
      ! assumed-shape-ok: configure-time, one call per run.
      integer :: nx, ny

      if (sf%use_components) return
      if (.not. allocated(sf%Q_heat) .or. .not. allocated(sf%Q_salt)) return
      nx = size(sf%Q_heat, 1)
      ny = size(sf%Q_heat, 2)
      if (size(cover_frac, 1) /= nx .or. size(cover_frac, 2) /= ny) return
      call ocean_surfflux_cover_const_impl(sf%Q_heat, sf%Q_salt, cover_frac, nx, ny)
   end subroutine ocean_surface_flux_apply_cover_const

   pure subroutine ocean_surfflux_cover_const_impl(Q_heat, Q_salt, cover_frac, nx, ny)
      !! Flat `do concurrent` kernel behind
      !! `ocean_surface_flux_apply_cover_const` — explicit-shape dummies,
      !! integer dims first (decl-order, ifx #8586).
      integer, intent(in)    :: nx, ny
      real(wp), intent(inout) :: Q_heat(nx, ny), Q_salt(nx, ny)
      real(wp), intent(in)    :: cover_frac(nx, ny)
      integer :: i, j
      real(wp) :: open_f
      do concurrent(j=1:ny, i=1:nx) local(open_f)
         open_f = 1.0_wp - cover_frac(i, j)
         Q_heat(i, j) = Q_heat(i, j)*open_f
         Q_salt(i, j) = Q_salt(i, j)*open_f
      end do
   end subroutine ocean_surfflux_cover_const_impl

   pure subroutine ocean_surface_flux_assemble(grid, sf, ms, active, cover_frac)
      !! **The single gate** that derives `Q_heat`/`Q_salt` from the
      !! component set (§3.1/§3.3 of the PR-12 plan) — the exact analogue
      !! of `vmix_assemble`: fillers contribute components, this routine
      !! alone derives the net fields every downstream kernel reads.  A
      !! no-op unless `sf%use_components` — with components off, `Q_heat`
      !! / `Q_salt` are exactly what `set_surface_flux_const` (or a
      !! field-override path) left them, byte-for-byte.
      !!
      !! Net surface heat into the ocean:
      !!   Q_heat = Q_heat_const + q_sw + q_lw + q_lat + q_sens + heat_added
      !!          + heat_cavity + heat_content_massin + heat_content_massout
      !!   heat_content_massin  = Σ heat_content_{lprec,fprec,vprec,
      !!                            lrunoff,frunoff,seaice_melt}
      !!   heat_content_massout = SEAWATER_CP * T_sst * evap
      !! Net surface salt flux:
      !!   Q_salt = Q_salt_const + salt_flux + salt_cavity
      !! All four outputs multiplied by `ms%wet_mask` (land carries
      !! exactly zero; interior loop bounds are NOT restricted — see
      !! CLAUDE.md's "nghost and the assembler loop bounds" gotcha).
      !!
      !! `has_heat`/`has_salt` are deliberately NOT touched here — a
      !! components-on run with no live filler must reproduce EXACTLY
      !! what `set_surface_flux_const` left them (the "+0.0 bit-identity
      !! trap": forcing them true would change `apply_surface_src_2d_impl`
      !! from an early-return to a `+0.0` stamp, altering a budget sum's
      !! operand COUNT even though the value is unchanged).  Setting
      !! `has_heat`/`has_salt` is the FILLER's job (see the module
      !! docstring's fill contract).
      !!
      !! Outer-shim + flat-impl (`ocean_surfflux_assemble_impl`): the SST
      !! read needs `ms%tracers(idx_temperature)%hTr`, an array-of-DT
      !! registry deref that must happen on the host before the `do
      !! concurrent` (NVHPC device codegen constraint — see
      !! `ocean_surface_flux_apply_tracers` for the identical pattern).
      !! No-op when no temperature tracer is registered (SST is
      !! undefined without one).
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(multilayer_state_t), intent(in) :: ms
      logical, intent(in), optional :: active
         !! Optional thermo-cadence gate.  Present-and-false ⇒ early
         !! return; absent ⇒ kernel runs (matches
         !! `ocean_surface_flux_apply_tracers`'s convention).
      real(wp), intent(in), optional :: cover_frac(:, :)
         !! Optional ice-shelf cover fraction (`metrics%cover_frac`,
         !! v1 binary).  Present ⇒ every ATMOSPHERIC contribution
         !! (`Q_heat_const`, `q_sw`, `q_lw`, `q_lat`, `q_sens`,
         !! `heat_added`, both mass-enthalpy terms, `Q_salt_const`,
         !! `salt_flux`) is scaled by `1 - cover_frac`, while the
         !! cavity's OWN `heat_cavity` / `salt_cavity` pass through
         !! unmasked.  This is the single place those two groups are
         !! still distinguishable — see the module docstring for why the
         !! mask lives here and not at apply time.  Absent ⇒ the
         !! original kernel, byte-identical.
      ! assumed-shape-ok: thermo-cadence shim, forwarded to an
      ! explicit-shape `_impl` before the device loop.

      integer :: nx, ny, nz, idx_T
      logical :: masked

      if (.not. sf%use_components) return
      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      idx_T = ms%idx_temperature
      if (idx_T <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      masked = .false.
      if (present(cover_frac)) then
         masked = (size(cover_frac, 1) == nx .and. size(cover_frac, 2) == ny)
      end if
      if (masked) then
         call ocean_surfflux_assemble_cover_impl( &
            sf%heat_content_massin, sf%heat_content_massout, sf%Q_heat, sf%Q_salt, &
            sf%q_sw, sf%q_lw, sf%q_lat, sf%q_sens, sf%heat_added, sf%heat_cavity, &
            sf%heat_content_lprec, sf%heat_content_fprec, sf%heat_content_vprec, &
            sf%heat_content_lrunoff, sf%heat_content_frunoff, sf%heat_content_seaice_melt, &
            sf%evap, sf%salt_flux, sf%salt_cavity, &
            ms%tracers(idx_T)%hTr, ms%h_layer, ms%wet_mask, cover_frac, &
            sf%Q_heat_const, sf%Q_salt_const, sf%cp, sf%h_min, nz, nx, ny)
         return
      end if

      call ocean_surfflux_assemble_impl( &
         sf%heat_content_massin, sf%heat_content_massout, sf%Q_heat, sf%Q_salt, &
         sf%q_sw, sf%q_lw, sf%q_lat, sf%q_sens, sf%heat_added, sf%heat_cavity, &
         sf%heat_content_lprec, sf%heat_content_fprec, sf%heat_content_vprec, &
         sf%heat_content_lrunoff, sf%heat_content_frunoff, sf%heat_content_seaice_melt, &
         sf%evap, sf%salt_flux, sf%salt_cavity, &
         ms%tracers(idx_T)%hTr, ms%h_layer, ms%wet_mask, &
         sf%Q_heat_const, sf%Q_salt_const, sf%cp, sf%h_min, nz, nx, ny)
   end subroutine ocean_surface_flux_assemble

   pure subroutine ocean_surfflux_assemble_impl(heat_content_massin, heat_content_massout, &
                                                Q_heat, Q_salt, &
                                                q_sw, q_lw, q_lat, q_sens, heat_added, &
                                                heat_cavity, &
                                                heat_content_lprec, heat_content_fprec, &
                                                heat_content_vprec, heat_content_lrunoff, &
                                                heat_content_frunoff, heat_content_seaice_melt, &
                                                evap, salt_flux, salt_cavity, &
                                                hTr_T, h_layer, wet_mask, &
                                                Q_heat_const, Q_salt_const, cp, h_min, &
                                                nz, nx, ny)
      !! Flat `do concurrent` kernel — explicit-shape dummies, integer
      !! dims declared first (decl-order, ifx #8586).  See
      !! `ocean_surface_flux_assemble` for the physics; this is the
      !! arithmetic verbatim.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: heat_content_massin(nx, ny), heat_content_massout(nx, ny)
      real(wp), intent(inout) :: Q_heat(nx, ny), Q_salt(nx, ny)
      real(wp), intent(in)    :: q_sw(nx, ny), q_lw(nx, ny), q_lat(nx, ny), q_sens(nx, ny)
      real(wp), intent(in)    :: heat_added(nx, ny), heat_cavity(nx, ny)
      real(wp), intent(in)    :: heat_content_lprec(nx, ny), heat_content_fprec(nx, ny)
      real(wp), intent(in)    :: heat_content_vprec(nx, ny), heat_content_lrunoff(nx, ny)
      real(wp), intent(in)    :: heat_content_frunoff(nx, ny), heat_content_seaice_melt(nx, ny)
      real(wp), intent(in)    :: evap(nx, ny), salt_flux(nx, ny), salt_cavity(nx, ny)
      real(wp), intent(in)    :: hTr_T(nx, ny, nz), h_layer(nx, ny, nz), wet_mask(nx, ny)
      real(wp), intent(in)    :: Q_heat_const, Q_salt_const, cp, h_min
      integer :: i, j
      real(wp) :: sst, massin, massout

      do concurrent(j=1:ny, i=1:nx) local(sst, massin, massout)
         massin = heat_content_lprec(i, j) + heat_content_fprec(i, j) &
                  + heat_content_vprec(i, j) + heat_content_lrunoff(i, j) &
                  + heat_content_frunoff(i, j) + heat_content_seaice_melt(i, j)
         sst = hTr_T(i, j, nz)/max(h_layer(i, j, nz), h_min)
         massout = cp*sst*evap(i, j)
         heat_content_massin(i, j) = wet_mask(i, j)*massin
         heat_content_massout(i, j) = wet_mask(i, j)*massout
         Q_heat(i, j) = wet_mask(i, j)* &
                        (Q_heat_const + q_sw(i, j) + q_lw(i, j) + q_lat(i, j) + &
                         q_sens(i, j) + heat_added(i, j) + heat_cavity(i, j) + &
                         massin + massout)
         Q_salt(i, j) = wet_mask(i, j)*(Q_salt_const + salt_flux(i, j) + &
                                        salt_cavity(i, j))
      end do
   end subroutine ocean_surfflux_assemble_impl

   pure subroutine ocean_surfflux_assemble_cover_impl(heat_content_massin, &
                                                      heat_content_massout, &
                                                      Q_heat, Q_salt, &
                                                      q_sw, q_lw, q_lat, q_sens, heat_added, &
                                                      heat_cavity, &
                                                      heat_content_lprec, heat_content_fprec, &
                                                      heat_content_vprec, heat_content_lrunoff, &
                                                      heat_content_frunoff, heat_content_seaice_melt, &
                                                      evap, salt_flux, salt_cavity, &
                                                      hTr_T, h_layer, wet_mask, cover_frac, &
                                                      Q_heat_const, Q_salt_const, cp, h_min, &
                                                      nz, nx, ny)
      !! Ice-shelf-cover twin of `ocean_surfflux_assemble_impl`: the
      !! open-water factor `open_f = 1 - cover_frac` multiplies the
      !! ATMOSPHERIC group and NOT the cavity group.  Grouping, spelt
      !! out because it is the whole point of this kernel:
      !!
      !!   masked   — `Q_heat_const`, `q_sw`, `q_lw`, `q_lat`, `q_sens`,
      !!              `heat_added`, `heat_content_massin` (the six
      !!              mass-enthalpy companions), `heat_content_massout`
      !!              (`cp·SST·evap`), `Q_salt_const`, `salt_flux`;
      !!   UNmasked — `heat_cavity`, `salt_cavity`.
      !!
      !! The two `heat_content_mass*` OUTPUTS carry the factor too: they
      !! are diagnostics of atmospheric mass exchange, which under a
      !! shelf is zero, and reporting the unmasked value next to a
      !! masked `Q_heat` would make the ledger not add up.
      !!
      !! Separate `_impl`, not an in-loop `present()` test (house
      !! idiom) — the cover-off path keeps the production assembler
      !! byte-identical.
      integer, intent(in)    :: nz, nx, ny
      real(wp), intent(inout) :: heat_content_massin(nx, ny), heat_content_massout(nx, ny)
      real(wp), intent(inout) :: Q_heat(nx, ny), Q_salt(nx, ny)
      real(wp), intent(in)    :: q_sw(nx, ny), q_lw(nx, ny), q_lat(nx, ny), q_sens(nx, ny)
      real(wp), intent(in)    :: heat_added(nx, ny), heat_cavity(nx, ny)
      real(wp), intent(in)    :: heat_content_lprec(nx, ny), heat_content_fprec(nx, ny)
      real(wp), intent(in)    :: heat_content_vprec(nx, ny), heat_content_lrunoff(nx, ny)
      real(wp), intent(in)    :: heat_content_frunoff(nx, ny), heat_content_seaice_melt(nx, ny)
      real(wp), intent(in)    :: evap(nx, ny), salt_flux(nx, ny), salt_cavity(nx, ny)
      real(wp), intent(in)    :: hTr_T(nx, ny, nz), h_layer(nx, ny, nz), wet_mask(nx, ny)
      real(wp), intent(in)    :: cover_frac(nx, ny)
      real(wp), intent(in)    :: Q_heat_const, Q_salt_const, cp, h_min
      integer :: i, j
      real(wp) :: sst, massin, massout, open_f

      do concurrent(j=1:ny, i=1:nx) local(sst, massin, massout, open_f)
         open_f = 1.0_wp - cover_frac(i, j)
         massin = open_f*(heat_content_lprec(i, j) + heat_content_fprec(i, j) &
                          + heat_content_vprec(i, j) + heat_content_lrunoff(i, j) &
                          + heat_content_frunoff(i, j) + heat_content_seaice_melt(i, j))
         sst = hTr_T(i, j, nz)/max(h_layer(i, j, nz), h_min)
         massout = open_f*cp*sst*evap(i, j)
         heat_content_massin(i, j) = wet_mask(i, j)*massin
         heat_content_massout(i, j) = wet_mask(i, j)*massout
         Q_heat(i, j) = wet_mask(i, j)* &
                        (open_f*(Q_heat_const + q_sw(i, j) + q_lw(i, j) + q_lat(i, j) + &
                                 q_sens(i, j) + heat_added(i, j)) + &
                         heat_cavity(i, j) + massin + massout)
         Q_salt(i, j) = wet_mask(i, j)*(open_f*(Q_salt_const + salt_flux(i, j)) + &
                                        salt_cavity(i, j))
      end do
   end subroutine ocean_surfflux_assemble_cover_impl

end module rdb_ocean_surface_flux
