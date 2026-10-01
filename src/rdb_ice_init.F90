!! Sea-ice initial-condition path (PR 24 — analytic v1).
module rdb_ice_init
   !! Seeds `ocean_sea_ice_t`'s six category prognostics
   !! (`part_size`/`m_ice`/`m_snow`/`enth_ice`/`enth_snow`/`sal_ice`) from
   !! `&ocean_ice_ic_nml`, so a run can *start* with a live ice pack
   !! instead of growing one from frazil (`PLAN_PR24_ice_ic_path.md`).
   !!
   !! v1 ships two ANALYTIC seeding modes on top of the default no-op:
   !!   `"zero"`      (default) — early return, touches nothing.
   !!   `"uniform"`   — scalar concentration/thickness everywhere wet.
   !!   `"latitudes"` — SIS2 `initialize_concentration_from_latitudes`
   !!                   (`SIS_state_initialization.F90:489-531`): a 0/1
   !!                   step at `arctic_edge`/`antarctic_edge` off
   !!                   `metrics%geolatT`. Defaults (+-91) => no ice.
   !!
   !! **Host-side setup code, called ONCE, before `ocean_state_enter_data`
   !! (`rdb_driver.F90`, between `configure_ocean_land_mask` and
   !! `ocean_state_enter_data`).** No device kernel, no new state array,
   !! zero `do concurrent`.  A `do concurrent` here would either
   !! round-trip the (unmapped) host arrays through the device once per
   !! loop on `-stdpar=gpu` (measured 36.9 s regression elsewhere,
   !! `seed_h_layer_uniform_impl`) or read stale device memory — see
   !! `rdb_ocean_state`'s "Plain host loops, NOT do concurrent" precedent.
   !!
   !! **Mass, not thickness — and per-ICE-area, not per-cell, at ncat>1.**
   !! The prognostic is `m_ice = ICE_RHO_ICE*h_ice` [kg/m^2], and at
   !! `ncat>1` it is per unit ICE-COVERED area (SIS2 `mH_ice`) — NOT
   !! multiplied by `conc` (module docstring of `rdb_ice_state`,
   !! ":65-86"). At `ncat==1` (legacy lumped mode) `m_ice` is per unit
   !! CELL area instead, and `part_size` is left untouched
   !! (`part_size(:,:,0)=1` forever — the frozen PR 3b/3c contract).
   !!
   !! **Category allocation is an ITD fixed point, by construction.** The
   !! pack's mass is binned directly against the already-computed
   !! `ice%mh_lim(1:ncat+1)` (filled at slot `init` by
   !! `ice_itd_category_bounds`) via `ice_ic_target_category` — the SAME
   !! bounds `ice_adjust_categories` (`rdb_ice_itd`) uses to restore the
   !! ITD every thermo window. Consumers of the ITD MUST read
   !! `ice%mh_lim`, never `HLIM_DFLT_TABLE` (that table is private to
   !! `rdb_ice_state`) and never a re-derivation from `ncat` alone — PR-58
   !! may override `mh_lim`'s VALUES via `&ocean_ice_nml hlim`, and this
   !! module's IC must stay consistent with whatever `init` actually
   !! computed. `ice_adjust_categories` applied to an IC-seeded state is
   !! therefore a bit-exact no-op (test oracle, not a call site — see
   !! "What NOT to take" below).
   !!
   !! **Enthalpy, not temperature.** `enth_ice`/`enth_snow` are set from
   !! the namelist ice temperature `t_ice` through the EXACT
   !! `ice_enth_from_ts` inversion (`rdb_ice_enthalpy`), never assigned
   !! directly — the thermodynamic prognostic is specific enthalpy
   !! (J/kg), never temperature. `sal_ice = s_ice`; snow is fresh
   !! (`ice_enth_from_ts(t_ice, 0.0)`).
   !!
   !! **k-convention.** `enth_ice`/`sal_ice` are `(..., nk_ice)` BOTTOM-UP
   !! (k=1 = ice bottom/ocean side, k=nk_ice = ice top/atm side) — the
   !! opposite of `rdb_ice_column`'s internal top-down convention. A
   !! uniform-`t_ice` IC is k-symmetric so this cannot be gotten wrong
   !! here, but the first person to add a linear top-to-bottom
   !! temperature profile must respect it.
   !!
   !! **What NOT to take from SIS2 (D-list, deliberate divergences):**
   !!   D1  Do NOT call `ice_adjust_categories` (SIS2's "seed into cat 1,
   !!       then adjust_ice_categories" flow) — it is a `do concurrent`
   !!       device kernel and the IC runs before `enter_data`. Bin
   !!       directly on the host against `mh_lim`; use the kernel only as
   !!       a TEST ORACLE.
   !!   D2  Do NOT seed a placeholder mass (`mh_lim(c+1)`) into empty
   !!       categories the way SIS2's `ice_state_mass_init` does. Roundabout's
   !!       live convention is `m_ice = 0` for an empty category (what
   !!       `ice_frazil_uptake_multicat_impl` already produces, and what
   !!       the column entry gate keys on) — adopting SIS2's placeholder
   !!       would run a full Winton solve on every empty category of every
   !!       cell for no change in the answer (part-weighted to zero).
   !!   D3  Do NOT port `ICE_RELATIVE_TEMP_IC`/`ICE_RELATIVE_SALINITY`/
   !!       `spec_thermo_sal` — v1 ships one bulk salinity + one bulk
   !!       temperature. Deferred.
   !!   D4  Do NOT port `data_override` / file-backed ICs — that is PR-14
   !!       (v1.1), explicitly out of scope here (`"file"` is deliberately
   !!       NOT in the `conc_config` allowed list).
   !!
   !! `"zero"` (the default) is a stronger no-op than SIS2's own IC path
   !! (SIS2 sets the thermo fields unconditionally): it returns before
   !! touching a single array, including `sal_ice`/`enth_ice` (already
   !! sourced to `ICE_BULK_SALINITY`/`0.0` by `ocean_sea_ice_init`) — so
   !! the restart file bytes are unchanged too.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_column, only: ICE_RHO_ICE, ICE_RHO_SNOW, ICE_BULK_SALINITY
   use rdb_ice_enthalpy, only: ice_enth_from_ts
   implicit none
   private

   public :: ice_ic_params_t
   public :: ice_ic_params_from_config
   public :: ice_ic_parse_conc_config
   public :: ice_ic_target_category
   public :: ice_init_apply
   public :: ICE_IC_CONC_ZERO, ICE_IC_CONC_UNIFORM, ICE_IC_CONC_LATITUDES, ICE_IC_CONC_INVALID

   integer, parameter :: ICE_IC_CONC_ZERO = 0
      !! Default: no-op, byte-identical to today.
   integer, parameter :: ICE_IC_CONC_UNIFORM = 1
      !! Scalar concentration/thickness/salinity/temperature everywhere wet.
   integer, parameter :: ICE_IC_CONC_LATITUDES = 2
      !! SIS2 polar-cap analytic form off `metrics%geolatT`.
   integer, parameter :: ICE_IC_CONC_INVALID = -1
      !! Sentinel for an unrecognised `conc_config` string (fail-loud
      !! idiom — a `pure` function cannot `error stop`; `validate_config`
      !! turns this into a `logger%error`).

   type :: ice_ic_params_t
      !! Sea-ice IC scalar knob bundle (mirrors `ice_evp_params_t`).
      !! Built ONCE from `&ocean_ice_ic_nml` via `ice_ic_params_from_config`;
      !! every field a plain scalar (no derived-type deref inside the
      !! host loop, same discipline as the EVP params — though this
      !! module never runs on-device at all).
      integer  :: conc_config = ICE_IC_CONC_ZERO
         !! Parsed `conc_config` enum (`ICE_IC_CONC_*`).
      real(wp) :: conc = 0.0_wp
         !! Uniform-mode concentration [0,1] (nondim).
      real(wp) :: h_ice = 0.0_wp
         !! Ice thickness (m) where seeded — SIS2 `ICE_INIT_MASS` as a
         !! thickness.
      real(wp) :: h_snow = 0.0_wp
         !! Snow thickness (m) where seeded — SIS2 `SNOW_INIT_MASS` as a
         !! thickness.
      real(wp) :: t_ice = -4.0_wp
         !! Ice/snow temperature (degC) — SIS2 `ICE_TEMPERATURE_IC`
         !! default.
      real(wp) :: s_ice = ICE_BULK_SALINITY
         !! Ice bulk salinity (PSU) — SIS2 `ICE_SALINITY_IC` default,
         !! which matches `ocean_sea_ice_init`'s own `sal_ice` source
         !! value, so `"uniform"` at the default `s_ice` leaves `sal_ice`
         !! unchanged.
      real(wp) :: arctic_edge = 91.0_wp
         !! `"latitudes"` Arctic edge (degrees_north) — SIS2
         !! `ARCTIC_ICE_EDGE_IC` default (no cell qualifies).
      real(wp) :: antarctic_edge = -91.0_wp
         !! `"latitudes"` Antarctic edge (degrees_north) — SIS2
         !! `ANTARCTIC_ICE_EDGE_IC` default (no cell qualifies).
   end type ice_ic_params_t

contains

   ! =====================================================================
   ! Enum parse + params constructor
   ! =====================================================================

   pure function ice_ic_parse_conc_config(s) result(mode)
      !! `"zero"`/`"uniform"`/`"latitudes"` -> `ICE_IC_CONC_*`; any other
      !! string -> `ICE_IC_CONC_INVALID` (fail-loud idiom: a `pure`
      !! function cannot abort, so `validate_config` V8 turns the
      !! sentinel into a `logger%error`; `nml_enum allowed=` is the first
      !! line of defence).
      character(len=*), intent(in) :: s
      integer :: mode

      select case (trim(s))
      case ("zero")
         mode = ICE_IC_CONC_ZERO
      case ("uniform")
         mode = ICE_IC_CONC_UNIFORM
      case ("latitudes")
         mode = ICE_IC_CONC_LATITUDES
      case default
         mode = ICE_IC_CONC_INVALID
      end select
   end function ice_ic_parse_conc_config

   pure function ice_ic_params_from_config(conc_config, conc, h_ice, h_snow, &
                                           t_ice, s_ice, arctic_edge, &
                                           antarctic_edge) result(par)
      !! Small constructor — build once from `&ocean_ice_ic_nml`. Copies
      !! `ice_evp_params_from_config`'s shape; the one difference is
      !! `conc_config`, a STRING here (parsed internally), because the
      !! namelist knob is a `nml_enum` string, not an already-parsed
      !! integer.
      character(len=*), intent(in) :: conc_config
      real(wp), intent(in) :: conc, h_ice, h_snow, t_ice, s_ice
      real(wp), intent(in) :: arctic_edge, antarctic_edge
      type(ice_ic_params_t) :: par

      par%conc_config = ice_ic_parse_conc_config(conc_config)
      par%conc = conc
      par%h_ice = h_ice
      par%h_snow = h_snow
      par%t_ice = t_ice
      par%s_ice = s_ice
      par%arctic_edge = arctic_edge
      par%antarctic_edge = antarctic_edge
   end function ice_ic_params_from_config

   ! =====================================================================
   ! ITD category binning
   ! =====================================================================

   pure function ice_ic_target_category(m, mh_lim, ncat) result(c)
      !! The IC's target category for a pack of mass `m` [kg/m^2]
      !! (per-ICE-area, ncat>1 convention): the top bin `ncat` is
      !! UNBOUNDED above (`mh_lim(ncat+1)` is stored but never used as a
      !! cap — `ice_adjust_categories`'s upward pass stops at `c =
      !! ncat-1`), otherwise the unique `c` with
      !! `mh_lim(c) <= m < mh_lim(c+1)`. Declared `ncat` before the
      !! explicit-shape `mh_lim` that uses it (decl-order).
      integer, intent(in) :: ncat
         !! Number of ice thickness categories.
      real(wp), intent(in) :: m
         !! Pack mass (kg/m^2, per-ICE-area convention).
      real(wp), intent(in) :: mh_lim(ncat + 1)
         !! Category lower mass limits (kg/m^2), `ice%mh_lim` — NEVER a
         !! re-derivation from `ncat` alone (PR-58 may override its
         !! values via `&ocean_ice_nml hlim`).
      integer :: c

      integer :: k

      if (m >= mh_lim(ncat)) then
         c = ncat
         return
      end if
      c = 1
      do k = 1, ncat - 1
         if (m >= mh_lim(k) .and. m < mh_lim(k + 1)) then
            c = k
            exit
         end if
      end do
   end function ice_ic_target_category

   ! =====================================================================
   ! Outer shim
   ! =====================================================================

   pure subroutine ice_init_apply(grid, ms, metrics, ice, par)
      !! Outer shim (outer-shim + flat-impl pattern, mirrors
      !! `ice_adjust_categories`). No-op when the ice slot is not live
      !! (`is_init` gate) or when `par%conc_config == ICE_IC_CONC_ZERO`
      !! (the default-off bit-identity contract — nothing is written, not
      !! even `sal_ice`/`enth_ice`). `ice` is `intent(inout)` (not
      !! `out`): `h_lim`/`mh_lim` were filled at `init` and must survive.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
         !! READ-ONLY: only `wet_mask` is read.
      type(ocean_metrics_t), intent(in) :: metrics
         !! READ-ONLY: only `geolatT` is read (`"latitudes"` mode).
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ice_ic_params_t), intent(in) :: par

      if (.not. ice%is_init) return
      if (par%conc_config == ICE_IC_CONC_ZERO) return

      call ice_init_apply_impl(ms%wet_mask, metrics%geolatT, &
                               ice%part_size, ice%m_ice, ice%m_snow, &
                               ice%enth_ice, ice%enth_snow, ice%sal_ice, &
                               ice%mh_lim, &
                               par%conc_config, par%conc, par%h_ice, par%h_snow, &
                               par%t_ice, par%s_ice, par%arctic_edge, par%antarctic_edge, &
                               ice%ncat, ice%nk_ice, grid%nx_total, grid%ny_total)
   end subroutine ice_init_apply

   ! =====================================================================
   ! Flat impl — plain host loops, ZERO `do concurrent` (see module
   ! docstring: the IC runs before `ocean_state_enter_data`, so a device
   ! loop here would either round-trip unmapped host arrays or read
   ! stale device memory on `-stdpar=gpu, mem:separate`).
   ! =====================================================================

   pure subroutine ice_init_apply_impl(wet_mask, geolatT, part_size, m_ice, m_snow, &
                                       enth_ice, enth_snow, sal_ice, mh_lim, &
                                       conc_config, conc_uniform, h_ice, h_snow, &
                                       t_ice, s_ice, arctic_edge, antarctic_edge, &
                                       ncat, nk_ice, nx, ny)
      !! Per-cell seeder. Decl-order: integer dims before the
      !! explicit-shape arrays that use them.
      integer, intent(in) :: ncat, nk_ice, nx, ny
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: geolatT(nx, ny)
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk_ice)
      real(wp), intent(inout) :: enth_snow(nx, ny, ncat, 1)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk_ice)
      real(wp), intent(in) :: mh_lim(ncat + 1)
      integer, intent(in) :: conc_config
      real(wp), intent(in) :: conc_uniform, h_ice, h_snow, t_ice, s_ice
      real(wp), intent(in) :: arctic_edge, antarctic_edge

      integer :: i, j, c_star
      real(wp) :: conc, m_target, enth_i, enth_s

      ! Precompute once: t_ice/s_ice are uniform scalars, so the mapped
      ! enthalpy is the same at every seeded cell/category/layer
      ! (k-symmetric, see module docstring).
      m_target = ICE_RHO_ICE*h_ice
      enth_i = ice_enth_from_ts(t_ice, s_ice)
      enth_s = ice_enth_from_ts(t_ice, 0.0_wp)

      ! Seed the FULL array, ghosts included (CLAUDE.md "formula
      ! bathymetry setters must fill ghost rows" gotcha wearing a
      ! different hat): ice_transport_step and the EVP periodic wrap both
      ! read neighbour ghost cells, so an unseeded ghost band would bleed
      ! the pack out at the domain edge / be flatly wrong under EVP
      ! periodicity.
      do j = 1, ny
         do i = 1, nx
            if (wet_mask(i, j) <= 0.5_wp) cycle
            ! Land: leave the `init`-time defaults
            ! (part_size(:,:,0)=1, m_ice=m_snow=0) untouched.

            select case (conc_config)
            case (ICE_IC_CONC_UNIFORM)
               conc = conc_uniform
            case (ICE_IC_CONC_LATITUDES)
               if (geolatT(i, j) > arctic_edge .or. geolatT(i, j) < antarctic_edge) then
                  conc = 1.0_wp
               else
                  conc = 0.0_wp
               end if
            case default
               ! Unreachable: `ice_init_apply` already gated ZERO, and
               ! `validate_config` V8 rejects every other integer.
               conc = 0.0_wp
            end select

            if (conc <= 0.0_wp) cycle
            ! Open water: leave the `init`-time defaults untouched.

            if (ncat == 1) then
               ! Legacy lumped mode (module docstring, `rdb_ice_state`):
               ! `part_size` is NEVER maintained here; `m_ice`/`m_snow`
               ! are per-CELL area.
               m_ice(i, j, 1) = m_target
               m_snow(i, j, 1) = ICE_RHO_SNOW*h_snow
               sal_ice(i, j, 1, :) = s_ice
               enth_ice(i, j, 1, :) = enth_i
               enth_snow(i, j, 1, 1) = enth_s
            else
               ! SIS2 ITD mode: `part_size` is LIVE, `m_ice`/`m_snow` are
               ! per unit ICE-COVERED area (NOT x conc — see module
               ! docstring). Bin against the ALREADY-COMPUTED `mh_lim`
               ! (never a re-derivation) so the seeded state is an exact
               ! fixed point of `ice_adjust_categories`.
               c_star = ice_ic_target_category(m_target, mh_lim, ncat)
               part_size(i, j, 0) = 1.0_wp - conc
               part_size(i, j, c_star) = conc
               m_ice(i, j, c_star) = m_target
               m_snow(i, j, c_star) = ICE_RHO_SNOW*h_snow
               sal_ice(i, j, c_star, :) = s_ice
               enth_ice(i, j, c_star, :) = enth_i
               enth_snow(i, j, c_star, 1) = enth_s
            end if
         end do
      end do
   end subroutine ice_init_apply_impl

end module rdb_ice_init
