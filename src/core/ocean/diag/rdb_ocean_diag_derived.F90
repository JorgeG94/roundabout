!! Library of derived ocean diagnostics — static catalog binding a
!! diagnostic name to a fill procedure + CF-1.8 metadata, picked by
!! name from the namelist.  Every fill is a pure function of the
!! current `ocean_state_t`; outputs at cell centres (vgrid LAYER,
!! z-level remap via `remap_layer_to_z`).  Includes the opt-in sea-ice
!! drift diagnostics (`ice_speed`/`ice_u`/`ice_v`) — the canonical
!! `ice_conc`/`ice_thick` fills live in `rdb_ocean_diag_fills` (module-
!! cycle rule: this module USES that one).
module rdb_ocean_diag_derived
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp, GRAVITY
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_diag, only: ocean_diag_t, diag_fill_proc, diag_remap_proc, &
                             DIAG_OP_MEAN, DIAG_OP_INSTANT, DIAG_OP_UNSET, &
                             DIAG_VGRID_LAYER, DIAG_VGRID_SURFACE, DIAG_COORD_UNSET, &
                             DIAG_VGRID_DENSITY, diag_spec_t, parse_diag_spec
   use rdb_ocean_diag_fills, only: register_default_diags, is_canonical_diag_name, &
                                   coord_remap_proc, diag_mask_vanished_is_on, &
                                   canonical_diag_gate_hint
   use pic_logger, only: global_logger
   use rdb_error_ring, only: error_ring_push
   implicit none
   private

   public :: derived_entry_t
   public :: register_derived, apply_diag_selection
   public :: derived_catalog_size, derived_catalog_name, derived_catalog_requires
   public :: fill_h_layer, fill_rho_layer, fill_vorticity_z
   public :: fill_ke_total, fill_transport_x, fill_transport_y
   public :: fill_mld_density
   public :: fill_ice_speed, fill_ice_u, fill_ice_v
   public :: fill_melt, fill_melt_m_per_yr, fill_thermal_driving
   public :: fill_haline_driving, fill_tbdry, fill_sbdry, fill_tfreeze_ib
   public :: fill_exch_vel_t, fill_exch_vel_s, fill_ustar_shelf
   public :: fill_cavity_melt_status, fill_z_draft, fill_water_column
   public :: cavity_mask_impl, melt_m_per_yr_factor

   integer, parameter, public :: DERIVED_REQ_NONE = 0
      !! No prerequisite — the entry reads state that always exists.
   integer, parameter, public :: DERIVED_REQ_CAVITY_DYN = 1
      !! Needs `&ocean_cavity_dyn_nml enable` (the draft geometry).
   integer, parameter, public :: DERIVED_REQ_CAVITY_MELT = 2
      !! Needs `&ocean_cavity_melt_nml enable` (the melt slot), which in
      !! turn requires `&ocean_cavity_dyn_nml`.

   type :: derived_entry_t
      !! One entry in the static catalog.  Buffer layout is implicit
      !! in `is_layered`: layered → (nx, ny, nz_ml), 2D → (nx, ny, 1).
      character(len=64)  :: name = ""
      character(len=128) :: long_name = ""
      character(len=32)  :: units = ""
      character(len=64)  :: standard_name = ""
      procedure(diag_fill_proc), pointer, nopass :: fill => null()
      logical :: is_layered = .true.
      integer :: requires = DERIVED_REQ_NONE
         !! Prerequisite knob, if any.  Registering an entry whose
         !! prerequisite is off FAILS LOUD at configure — the same
         !! stance `register_derived` already takes on an unknown name,
         !! and for the same reason: a requested diagnostic that comes
         !! back as a plane of missing values, or worse as a plane of
         !! zeros, is indistinguishable from a physical answer.
   end type derived_entry_t

   integer, parameter :: N_CATALOG = 23
   type(derived_entry_t) :: CATALOG(N_CATALOG)
   logical :: catalog_initialised = .false.

   real(wp), parameter :: MLD_DENSITY_THRESHOLD = 0.03_wp
      !! De Boyer Montégut (2004) MLD criterion: Δσ_0 = 0.03 kg/m³ vs surface.

   real(wp), parameter :: RHO_FRESHWATER_ISOMIP = 1000.0_wp
      !! Freshwater density (kg/m³) the ISOMIP+ melt-rate CONVENTION
      !! divides by — Asay-Davis et al. (2016), *GMD* **9**, 2471–2497,
      !! §3.3: melt rates are reported in m/yr of ICE-EQUIVALENT
      !! freshwater.  It is a reporting convention, deliberately NOT the
      !! model's Boussinesq `rho_0` (1035) and NOT the melt law's own
      !! `const%rho_w`: quoting `melt_m_per_yr` against anything else
      !! puts the number 3.5 % off every published ISOMIP+ figure.
   real(wp), parameter :: SECONDS_PER_YEAR = 365.0_wp*86400.0_wp
      !! 365-day year, the ISOMIP+ protocol's own calendar.  A Julian
      !! year (365.25 d) would move every melt rate by 0.07 %, which is
      !! below the noise but not below the reader's attention.

contains

   subroutine ensure_catalog_initialised()
      !! Populate the catalog at runtime (procedure pointers can't be a
      !! parameter constructor).  Idempotent.
      if (catalog_initialised) return
      CATALOG(1) = derived_entry_t( &
                   name="h_layer", &
                   long_name="layer_thickness", &
                   units="m", &
                   standard_name="ocean_layer_thickness", &
                   fill=fill_h_layer, is_layered=.true.)
      CATALOG(2) = derived_entry_t( &
                   name="rho_layer", &
                   long_name="layer_in_situ_density", &
                   units="kg m-3", &
                   standard_name="sea_water_density", &
                   fill=fill_rho_layer, is_layered=.true.)
      CATALOG(3) = derived_entry_t( &
                   name="vorticity_z", &
                   long_name="vertical_relative_vorticity_at_cell_centre", &
                   units="s-1", &
                   standard_name="ocean_relative_vorticity", &
                   fill=fill_vorticity_z, is_layered=.true.)
      CATALOG(4) = derived_entry_t( &
                   name="ke_total", &
                   long_name="depth_integrated_kinetic_energy", &
                   units="m3 s-2", &
                   standard_name="ocean_kinetic_energy_content", &
                   fill=fill_ke_total, is_layered=.false.)
      CATALOG(5) = derived_entry_t( &
                   name="transport_x", &
                   long_name="depth_integrated_eastward_transport", &
                   units="m2 s-1", &
                   standard_name="ocean_volume_transport_x", &
                   fill=fill_transport_x, is_layered=.false.)
      CATALOG(6) = derived_entry_t( &
                   name="transport_y", &
                   long_name="depth_integrated_northward_transport", &
                   units="m2 s-1", &
                   standard_name="ocean_volume_transport_y", &
                   fill=fill_transport_y, is_layered=.false.)
      CATALOG(7) = derived_entry_t( &
                   name="mld_density", &
                   long_name="mixed_layer_depth_density_threshold", &
                   units="m", &
                   standard_name="ocean_mixed_layer_thickness_defined_by_sigma_t", &
                   fill=fill_mld_density, is_layered=.false.)
      CATALOG(8) = derived_entry_t( &
                   name="ice_speed", &
                   long_name="sea_ice_drift_speed_at_cell_centre", &
                   units="m s-1", &
                   standard_name="sea_ice_speed", &
                   fill=fill_ice_speed, is_layered=.false.)
      CATALOG(9) = derived_entry_t( &
                   name="ice_u", &
                   long_name="eastward_sea_ice_velocity_at_cell_centre", &
                   units="m s-1", &
                   standard_name="sea_ice_x_velocity", &
                   fill=fill_ice_u, is_layered=.false.)
      CATALOG(10) = derived_entry_t( &
                    name="ice_v", &
                    long_name="northward_sea_ice_velocity_at_cell_centre", &
                    units="m s-1", &
                    standard_name="sea_ice_y_velocity", &
                    fill=fill_ice_v, is_layered=.false.)
      ! ---- Ice-shelf cavity (P2c).  Eleven melt-interface fields, all
      ! gated on `&ocean_cavity_melt_nml enable`, plus two GEOMETRY
      ! fields that need only `&ocean_cavity_dyn_nml`.  Every one is
      ! NaN outside the cover, so a domain mean over the output is a
      ! mean over the CAVITY — see `cavity_mask_impl`.
      CATALOG(11) = derived_entry_t( &
                    name="melt", &
                    long_name="ice_shelf_basal_melt_mass_flux", &
                    units="kg m-2 s-1", &
                    standard_name="water_flux_into_sea_water_from_ice_shelf", &
                    fill=fill_melt, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(12) = derived_entry_t( &
                    name="melt_m_per_yr", &
                    long_name="ice_shelf_basal_melt_rate_isomip_convention", &
                    units="m yr-1", &
                    standard_name="ice_shelf_basal_melt_rate", &
                    fill=fill_melt_m_per_yr, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(13) = derived_entry_t( &
                    name="thermal_driving", &
                    long_name="far_field_thermal_driving_above_in_situ_freezing_point", &
                    units="degC", &
                    standard_name="ocean_thermal_driving", &
                    fill=fill_thermal_driving, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(14) = derived_entry_t( &
                    name="haline_driving", &
                    long_name="far_field_haline_driving_above_interface_salinity", &
                    units="g kg-1", &
                    standard_name="ocean_haline_driving", &
                    fill=fill_haline_driving, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(15) = derived_entry_t( &
                    name="tbdry", &
                    long_name="ice_ocean_interface_temperature", &
                    units="degC", &
                    standard_name="sea_water_temperature_at_ice_shelf_base", &
                    fill=fill_tbdry, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(16) = derived_entry_t( &
                    name="sbdry", &
                    long_name="ice_ocean_interface_salinity", &
                    units="g kg-1", &
                    standard_name="sea_water_salinity_at_ice_shelf_base", &
                    fill=fill_sbdry, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(17) = derived_entry_t( &
                    name="tfreeze_ib", &
                    long_name="in_situ_freezing_point_of_the_far_field_at_the_ice_base", &
                    units="degC", &
                    standard_name="sea_water_freezing_temperature", &
                    fill=fill_tfreeze_ib, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(18) = derived_entry_t( &
                    name="exch_vel_t", &
                    long_name="thermal_exchange_velocity_at_the_ice_base", &
                    units="m s-1", &
                    standard_name="ocean_thermal_exchange_velocity", &
                    fill=fill_exch_vel_t, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(19) = derived_entry_t( &
                    name="exch_vel_s", &
                    long_name="haline_exchange_velocity_at_the_ice_base", &
                    units="m s-1", &
                    standard_name="ocean_haline_exchange_velocity", &
                    fill=fill_exch_vel_s, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(20) = derived_entry_t( &
                    name="ustar_shelf", &
                    long_name="friction_velocity_of_the_ice_shelf_melt_law", &
                    units="m s-1", &
                    standard_name="ocean_friction_velocity_at_ice_shelf_base", &
                    fill=fill_ustar_shelf, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(21) = derived_entry_t( &
                    name="cavity_melt_status", &
                    long_name="per_column_basal_melt_solver_status_code", &
                    units="1", &
                    standard_name="ocean_basal_melt_solver_status", &
                    fill=fill_cavity_melt_status, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_MELT)
      CATALOG(22) = derived_entry_t( &
                    name="z_draft", &
                    long_name="ice_shelf_draft_below_the_geoid", &
                    units="m", &
                    standard_name="ice_shelf_draft", &
                    fill=fill_z_draft, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_DYN)
      CATALOG(23) = derived_entry_t( &
                    name="water_column", &
                    long_name="water_column_thickness_under_the_ice_shelf", &
                    units="m", &
                    standard_name="sea_water_column_thickness", &
                    fill=fill_water_column, is_layered=.false., &
                    requires=DERIVED_REQ_CAVITY_DYN)
      catalog_initialised = .true.
   end subroutine ensure_catalog_initialised

   function derived_catalog_size() result(n)
      !! Public only for the unit-test suite.
      integer :: n
      call ensure_catalog_initialised()
      n = N_CATALOG
   end function derived_catalog_size

   function derived_catalog_name(i) result(name)
      !! Public only for the unit-test suite.
      integer, intent(in) :: i
      character(len=64) :: name
      call ensure_catalog_initialised()
      name = CATALOG(i)%name
   end function derived_catalog_name

   function derived_catalog_requires(name) result(req)
      !! The `DERIVED_REQ_*` code a catalog entry is gated on, `-1` for
      !! an unknown name.  The DECISION `register_derived` fails loud on,
      !! exposed as a lookup so the suite can assert it without
      !! provoking `error stop` — the repo's standing pattern for
      !! testing a fail-loud rule.
      character(len=*), intent(in) :: name
      integer :: req
      integer :: i
      call ensure_catalog_initialised()
      req = -1
      do i = 1, N_CATALOG
         if (trim(CATALOG(i)%name) == trim(name)) then
            req = CATALOG(i)%requires
            return
         end if
      end do
   end function derived_catalog_requires

   subroutine register_derived(state, name, time_op, dt_out, coord)
      !! Register ONE derived diagnostic by catalog `name` (used by
      !! `apply_diag_selection`): look it up and forward to
      !! `state%diag%register(...)` with the catalog metadata + correct buffer
      !! shape.  Error-stops on an unknown name.  `coord` (a `DIAG_VGRID_*`)
      !! sets the output vgrid for a LAYERED diag + attaches the conservative
      !! remap (2D entries ignore it); default LAYER.  Remaps INTENSIVE.
      type(ocean_state_t), intent(inout), target :: state
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: time_op
      real(wp), intent(in), optional :: dt_out
      integer, intent(in), optional :: coord
      integer :: i, idx, nx, ny, nz, n3, ocoord
      type(derived_entry_t) :: entry
      procedure(diag_remap_proc), pointer :: remap

      call ensure_catalog_initialised()

      idx = 0
      do i = 1, N_CATALOG
         if (trim(CATALOG(i)%name) == trim(name)) then
            idx = i
            exit
         end if
      end do
      if (idx == 0) then
         call error_ring_push("unknown derived diagnostic '"//trim(name)// &
                              "'; valid names: "//catalog_name_list())
         call global_logger%error("unknown derived diagnostic '"//trim(name)// &
                                  "'; valid names: "//catalog_name_list())
         error stop "register_derived: unknown derived diagnostic name"
      end if
      entry = CATALOG(idx)

      ! Prerequisite gate.  FAIL LOUD, never a silent plane of zeros or
      ! of missing values: a user who asked for `melt` on a run with no
      ! melt slot has made a configuration error, and handing them a
      ! NetCDF variable full of `_FillValue` would let it reach a figure.
      select case (entry%requires)
      case (DERIVED_REQ_CAVITY_DYN)
         if (.not. state%metrics%use_cavity) then
            call derived_requires_fail(trim(entry%name), "&ocean_cavity_dyn_nml enable", &
                                       "the ice-shelf draft geometry")
         end if
      case (DERIVED_REQ_CAVITY_MELT)
         if (.not. state%cavity_flux%enable) then
            call derived_requires_fail(trim(entry%name), "&ocean_cavity_melt_nml enable", &
                                       "the basal-melt interface slot")
         end if
      case default
         continue
      end select

      nx = size(state%barotropic%h, 1)
      ny = size(state%barotropic%h, 2)
      nz = state%multilayer%nz_ml
      n3 = nz
      if (.not. entry%is_layered) n3 = 1

      ocoord = DIAG_VGRID_LAYER
      if (present(coord) .and. entry%is_layered) ocoord = coord
      call coord_remap_proc(ocoord, remap)

      if (associated(remap)) then
         call state%diag%register(name=trim(entry%name), units=trim(entry%units), &
                                  fill=entry%fill, n1=nx, n2=ny, n3=n3, &
                                  long_name=trim(entry%long_name), &
                                  standard_name=trim(entry%standard_name), &
                                  time_op=time_op, dt_out=dt_out, &
                                  output_vgrid=ocoord, remap=remap, is_extensive=.false., &
                                  has_missing=(diag_mask_vanished_is_on() .and. &
                                               ocoord /= DIAG_VGRID_DENSITY))
      else
         ! Every cavity entry writes the NaN sentinel outside the cover
         ! (`cavity_mask_impl`), so the NetCDF variable must advertise a
         ! `_FillValue` — unconditionally, not via
         ! `diag_mask_vanished_is_on()`, which gates the REMAP path's
         ! below-target fill and has nothing to do with a calving front.
         call state%diag%register(name=trim(entry%name), units=trim(entry%units), &
                                  fill=entry%fill, n1=nx, n2=ny, n3=n3, &
                                  long_name=trim(entry%long_name), &
                                  standard_name=trim(entry%standard_name), &
                                  time_op=time_op, dt_out=dt_out, &
                                  has_missing=(entry%requires /= DERIVED_REQ_NONE))
      end if
   end subroutine register_derived

   subroutine derived_requires_fail(name, knob, what)
      !! Fail loud on a derived diagnostic whose prerequisite knob is
      !! off.  Same three-step shape `register_derived` uses for an
      !! unknown name — error ring, logger, `error stop` — so the C ABI
      !! and the console both see it.
      character(len=*), intent(in) :: name
         !! Catalog name the user asked for.
      character(len=*), intent(in) :: knob
         !! The namelist key that has to be on.
      character(len=*), intent(in) :: what
         !! One phrase naming what that knob builds.
      character(len=:), allocatable :: msg
      msg = "&ocean_diag_nml diags requested '"//name//"', which reads "//what// &
            " and therefore requires "//knob//"=.true.  Refused rather than "// &
            "registered: the fill would be a plane of missing values, which is "// &
            "indistinguishable from a physical answer once it reaches a figure."
      call error_ring_push(msg)
      call global_logger%error(msg)
      error stop "register_derived: derived diagnostic requires a knob that is off"
   end subroutine derived_requires_fail

   subroutine apply_diag_selection(state, spec, dt_out, default_coord)
      !! Configure the ocean diagnostic set from the unified
      !! `&ocean_diag_nml diags` selection string.  Single production entry
      !! point: parses `spec` once, registers the canonical defaults
      !! (consulting the spec for per-diagnostic `:off` skips and
      !! `:cadence` / `:op` / `:coord` overrides), then registers any spec
      !! entries that name a derived-catalog diagnostic (with their overrides).
      !!
      !! `default_coord` (a `DIAG_VGRID_*`, default LAYER) is the global
      !! output vgrid (`&ocean_diag_nml vgrid`) applied to layered diagnostics
      !! that don't carry a per-diag `:coord`.
      !!
      !! An empty / blank `spec` + LAYER default reproduces the canonical set
      !! exactly (bit-identical).  A spec name that is neither canonical nor
      !! in the derived catalog fails loud via `register_derived`.  A `:off`
      !! on a non-canonical name that is not registered is a harmless no-op.
      !!
      !! Must be called AFTER the output-level setters (`set_output_z_levels`
      !! / `set_output_sigma_levels` / `set_output_zstar_levels` /
      !! `set_output_density_levels`) so non-layer diagnostics size their
      !! remap targets correctly.
      type(ocean_state_t), intent(inout), target :: state
      character(len=*), intent(in) :: spec
      real(wp), intent(in), optional :: dt_out
      integer, intent(in), optional :: default_coord
      type(diag_spec_t), allocatable :: specs(:)
      integer :: i, op, dcoord, coord
      real(wp) :: dtout, dto

      dtout = 3600.0_wp
      if (present(dt_out)) dtout = dt_out
      dcoord = DIAG_VGRID_LAYER
      if (present(default_coord)) dcoord = default_coord

      specs = parse_diag_spec(spec)

      ! Canonical set, with the spec's per-diagnostic skips / overrides.
      call register_default_diags(state, dt_out=dtout, specs=specs, default_coord=dcoord)

      ! Derived-catalog additions: any spec entry that is not a canonical
      ! diagnostic and is not turned off.
      do i = 1, size(specs)
         if (specs(i)%off) cycle
         if (is_canonical_diag_name(trim(specs(i)%name))) then
            if (.not. state%diag%is_registered(trim(specs(i)%name))) then
               call global_logger%warning( &
                  "&ocean_diag_nml diags requested '"//trim(specs(i)%name)// &
                  "', which is a canonical diagnostic but is NOT registered — it will "// &
                  "be absent from the output. "//gate_clause(trim(specs(i)%name)))
            end if
            cycle
         end if
         op = DIAG_OP_INSTANT
         if (specs(i)%time_op /= DIAG_OP_UNSET) op = specs(i)%time_op
         dto = dtout
         if (specs(i)%dt_out > 0.0_wp) dto = specs(i)%dt_out
         coord = dcoord
         if (specs(i)%coord /= DIAG_COORD_UNSET) coord = specs(i)%coord
         call register_derived(state, trim(specs(i)%name), time_op=op, dt_out=dto, coord=coord)
      end do
   end subroutine apply_diag_selection

   pure function gate_clause(name) result(clause)
      !! Render `canonical_diag_gate_hint(name)` into a human-readable
      !! sentence for the "requested but not registered" warning.  Phrased
      !! as a QUESTION, never a diagnosis: the hint is a hand-maintained
      !! mirror of `register_default_diags`'s gates (§11.2 of the PR-64
      !! plan) and can drift from the true cause.  An empty hint means the
      !! diagnostic is one of the four unconditional ones (SSH/u/v/KE) and
      !! SHOULD have registered — that is a bug, not a closed gate.
      character(len=*), intent(in) :: name
      character(len=96) :: clause
      character(len=64) :: hint

      hint = canonical_diag_gate_hint(name)
      if (len_trim(hint) > 0) then
         clause = "Is "//trim(hint)//" off?"
      else
         clause = "This diagnostic has no feature gate — please report this."
      end if
   end function gate_clause

   function catalog_name_list() result(list)
      !! Comma-separated list of catalog diagnostic names (for error text).
      character(len=:), allocatable :: list
      integer :: i
      call ensure_catalog_initialised()
      list = ""
      do i = 1, N_CATALOG
         if (i > 1) list = list//", "
         list = list//trim(CATALOG(i)%name)
      end do
   end function catalog_name_list

   ! ---------------------------------------------------------------------
   ! Fill procedures
   ! ---------------------------------------------------------------------

   ! Outer-shim + flat-impl: public fills are host-side `select type`
   ! shims that deref the multilayer slot then forward device-resident
   ! arrays to a `do concurrent` `_impl` — keeps polymorphic dispatch out
   ! of the device kernel.

   subroutine fill_h_layer(state_handle, buf)
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call copy3_impl(state%multilayer%h_layer, buf)
      end select
   end subroutine fill_h_layer

   subroutine fill_rho_layer(state_handle, buf)
      !! In-situ density per layer — direct read of the EOS slot (driver
      !! must have run `ocean_eos_compute` this step; not re-invoked here).
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call copy3_impl(state%multilayer%rho_layer, buf)
      end select
   end subroutine fill_rho_layer

   pure subroutine copy3_impl(src, buf)
      !! Shared device copy `buf = src` with shape clipping — every
      !! direct-copy derived field routes through here.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: src(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = min(size(buf, 1), size(src, 1))
      ny = min(size(buf, 2), size(src, 2))
      nz = min(size(buf, 3), size(src, 3))
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = src(i, j, k)
      end do
   end subroutine copy3_impl

   subroutine fill_vorticity_z(state_handle, buf)
      !! Relative vorticity ζ = ∂v/∂x − ∂u/∂y, averaged from the four
      !! surrounding C-grid corners onto the cell centre.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%metrics%is_init) then
            call zero3_impl(buf)
            return
         end if
         call fill_vorticity_z_impl(state%multilayer%u_face_x_layer, &
                                    state%multilayer%v_face_y_layer, &
                                    state%metrics%idxT, &
                                    state%metrics%idyT, &
                                    state%multilayer%nz_ml, buf)
      end select
   end subroutine fill_vorticity_z

   pure subroutine fill_vorticity_z_impl(u_face, v_face, idxT, idyT, nz_ml, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized u_face/v_face have nx+1/ny+1 dims; size() min-clips at call.
      !! Cell-centred relative vorticity using per-cell metric inverses
      !! `idxT`/`idyT` (= 1/dx, 1/dy on uniform Cartesian); the 0.25
      !! corner-average factor folds into the per-cell scale.
      real(wp), intent(in)    :: u_face(:, :, :), v_face(:, :, :)
      ! assumed-shape-ok: diag fill — cadence-bounded (once per output frame)
      real(wp), intent(in)    :: idxT(:, :), idyT(:, :)
      integer, intent(in)    :: nz_ml
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      real(wp) :: inv_dx, inv_dy, dvdx, dudy
      call zero3_impl(buf)
      nx = min(size(buf, 1), size(u_face, 1) - 1, size(idxT, 1))
      ny = min(size(buf, 2), size(v_face, 2) - 1, size(idyT, 2))
      nz = min(size(buf, 3), nz_ml)
      ! Interior cells only (i in [2, nx-1], j in [2, ny-1]) so the
      ! ±1 stencil stays in bounds; the rim stays at the zero seed.
      do concurrent(k=1:nz, j=2:ny - 1, i=2:nx - 1) &
         local(dvdx, dudy, inv_dx, inv_dy)
         inv_dx = 0.25_wp*idxT(i, j)
         inv_dy = 0.25_wp*idyT(i, j)
         dvdx = inv_dx*( &
                v_face(i + 1, j, k) - v_face(i - 1, j, k) + &
                v_face(i + 1, j + 1, k) - v_face(i - 1, j + 1, k))
         dudy = inv_dy*( &
                u_face(i, j + 1, k) - u_face(i, j - 1, k) + &
                u_face(i + 1, j + 1, k) - u_face(i + 1, j - 1, k))
         buf(i, j, k) = dvdx - dudy
      end do
   end subroutine fill_vorticity_z_impl

   pure subroutine zero3_impl(buf)
      !! Device-side zero of `buf` (seeds column sums / clean bail-out).
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: i, j, k, nx, ny, nz
      nx = size(buf, 1)
      ny = size(buf, 2)
      nz = size(buf, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = 0.0_wp
      end do
   end subroutine zero3_impl

   subroutine fill_ke_total(state_handle, buf)
      !! Depth-integrated kinetic energy at each cell:
      !!   KE_total(i, j) = Σ_k 0.5 · h_layer(k) · (u_c² + v_c²)
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_ke_total_impl(state%multilayer%h_layer, &
                                 state%multilayer%u_face_x_layer, &
                                 state%multilayer%v_face_y_layer, buf)
      end select
   end subroutine fill_ke_total

   pure subroutine fill_ke_total_impl(h_layer, u_face, v_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized arrays have nx+1/ny+1 dims; size() min-clips at call.
      real(wp), intent(in)    :: h_layer(:, :, :)
      real(wp), intent(in)    :: u_face(:, :, :), v_face(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded; face-sized dims
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      real(wp) :: uc, vc, col_ke
      nx = min(size(buf, 1), size(u_face, 1) - 1, size(v_face, 1))
      ny = min(size(buf, 2), size(u_face, 2), size(v_face, 2) - 1)
      nz = min(size(u_face, 3), size(v_face, 3), size(h_layer, 3))
      ! Per-cell column sum: outer DC over (j, i), inner serial over k.
      do concurrent(j=1:ny, i=1:nx) &
         local(uc, vc, col_ke, k)
         col_ke = 0.0_wp
         do k = 1, nz
            uc = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
            vc = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
            col_ke = col_ke + 0.5_wp*h_layer(i, j, k)*(uc*uc + vc*vc)
         end do
         buf(i, j, 1) = col_ke
      end do
   end subroutine fill_ke_total_impl

   subroutine fill_transport_x(state_handle, buf)
      !! Depth-integrated zonal transport (m²/s):
      !!   T_x(i, j) = Σ_k h_c · u_c   (cell-centre form)
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_transport_x_impl(state%multilayer%h_layer, &
                                    state%multilayer%u_face_x_layer, buf)
      end select
   end subroutine fill_transport_x

   pure subroutine fill_transport_x_impl(h_layer, u_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h_layer(:, :, :), u_face(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      real(wp) :: u_c, col
      nx = min(size(buf, 1), size(u_face, 1) - 1)
      ny = min(size(buf, 2), size(u_face, 2))
      nz = min(size(u_face, 3), size(h_layer, 3))
      do concurrent(j=1:ny, i=1:nx) &
         local(u_c, col, k)
         col = 0.0_wp
         do k = 1, nz
            u_c = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
            col = col + h_layer(i, j, k)*u_c
         end do
         buf(i, j, 1) = col
      end do
   end subroutine fill_transport_x_impl

   subroutine fill_transport_y(state_handle, buf)
      !! Depth-integrated meridional transport (m²/s) — mirror of x.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_transport_y_impl(state%multilayer%h_layer, &
                                    state%multilayer%v_face_y_layer, buf)
      end select
   end subroutine fill_transport_y

   pure subroutine fill_transport_y_impl(h_layer, v_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h_layer(:, :, :), v_face(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      real(wp) :: v_c, col
      nx = min(size(buf, 1), size(v_face, 1))
      ny = min(size(buf, 2), size(v_face, 2) - 1)
      nz = min(size(v_face, 3), size(h_layer, 3))
      do concurrent(j=1:ny, i=1:nx) &
         local(v_c, col, k)
         col = 0.0_wp
         do k = 1, nz
            v_c = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
            col = col + h_layer(i, j, k)*v_c
         end do
         buf(i, j, 1) = col
      end do
   end subroutine fill_transport_y_impl

   subroutine fill_mld_density(state_handle, buf)
      !! Mixed-layer depth via the de Boyer Montégut threshold.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: nz_ml
      select type (state => state_handle)
      class is (ocean_state_t)
         nz_ml = state%multilayer%nz_ml
         if (nz_ml < 1) then
            call zero3_impl(buf)
            return
         end if
         call fill_mld_density_impl(state%multilayer%h_layer, &
                                    state%multilayer%rho_layer, &
                                    nz_ml, buf)
      end select
      if (.false.) buf(1, 1, 1) = GRAVITY  ! keep GRAVITY import live
   end subroutine fill_mld_density

   pure subroutine fill_mld_density_impl(h_layer, rho_layer, nz_ml, buf)
      !! Per-column scan from surface (k=nz_ml) toward bed; first layer
      !! whose ρ exceeds (surface ρ + MLD_DENSITY_THRESHOLD) marks the
      !! MLD as the cumulative h-sum above it.  No crossing → MLD = full
      !! column depth.  Threshold met at the surface itself → MLD = 0.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h_layer(:, :, :), rho_layer(:, :, :)
      integer, intent(in)    :: nz_ml
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny
      real(wp) :: rho_surf, d_acc, mld
      logical :: crossed
      nx = min(size(buf, 1), size(rho_layer, 1))
      ny = min(size(buf, 2), size(rho_layer, 2))
      do concurrent(j=1:ny, i=1:nx) &
         local(rho_surf, d_acc, mld, crossed, k)
         rho_surf = rho_layer(i, j, nz_ml)
         d_acc = 0.0_wp
         mld = 0.0_wp
         crossed = .false.
         do k = nz_ml, 1, -1
            if (rho_layer(i, j, k) - rho_surf >= MLD_DENSITY_THRESHOLD) then
               mld = d_acc
               crossed = .true.
               exit
            end if
            d_acc = d_acc + h_layer(i, j, k)
         end do
         if (.not. crossed) mld = d_acc
         buf(i, j, 1) = mld
      end do
   end subroutine fill_mld_density_impl

   subroutine fill_ice_speed(state_handle, buf)
      !! Sea-ice drift speed |u_ice| at T-centres (m/s) — C-face pairs
      !! averaged to centre, 2-D twin of the ocean `fill_ke` stencil.
      !! Ice off / not init ⇒ zeros (never registered by default; opt-in
      !! via `&ocean_diag_nml diags`).
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%ice%is_init) then
            call zero3_impl(buf)
            return
         end if
         call fill_ice_speed_impl(state%ice%u_ice, state%ice%v_ice, buf)
      end select
   end subroutine fill_ice_speed

   pure subroutine fill_ice_speed_impl(u_ice, v_ice, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized u_ice/v_ice have nx+1/ny+1 dims; size() min-clips at call.
      real(wp), intent(in)    :: u_ice(:, :), v_ice(:, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: uc, vc
      nx = min(size(buf, 1), size(u_ice, 1) - 1, size(v_ice, 1))
      ny = min(size(buf, 2), size(u_ice, 2), size(v_ice, 2) - 1)
      do concurrent(j=1:ny, i=1:nx) local(uc, vc)
         uc = 0.5_wp*(u_ice(i, j) + u_ice(i + 1, j))
         vc = 0.5_wp*(v_ice(i, j) + v_ice(i, j + 1))
         buf(i, j, 1) = sqrt(uc*uc + vc*vc)
      end do
   end subroutine fill_ice_speed_impl

   subroutine fill_ice_u(state_handle, buf)
      !! Eastward sea-ice velocity at T-centres (m/s) — u-face pair
      !! averaged to centre, 2-D twin of `fill_u_centre`.  Ice off / not
      !! init ⇒ zeros.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%ice%is_init) then
            call zero3_impl(buf)
            return
         end if
         call fill_ice_u_impl(state%ice%u_ice, buf)
      end select
   end subroutine fill_ice_u

   pure subroutine fill_ice_u_impl(u_ice, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized u_ice(:,:) has nx+1 first dim, incompatible with cell-sized buf.
      real(wp), intent(in)    :: u_ice(:, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      nx = min(size(buf, 1), size(u_ice, 1) - 1)
      ny = min(size(buf, 2), size(u_ice, 2))
      do concurrent(j=1:ny, i=1:nx)
         buf(i, j, 1) = 0.5_wp*(u_ice(i, j) + u_ice(i + 1, j))
      end do
   end subroutine fill_ice_u_impl

   subroutine fill_ice_v(state_handle, buf)
      !! Northward sea-ice velocity at T-centres (m/s) — v-face pair
      !! averaged to centre, 2-D twin of `fill_v_centre`.  Ice off / not
      !! init ⇒ zeros.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%ice%is_init) then
            call zero3_impl(buf)
            return
         end if
         call fill_ice_v_impl(state%ice%v_ice, buf)
      end select
   end subroutine fill_ice_v

   ! ---------------------------------------------------------------------
   ! Ice-shelf cavity (P2c)
   ! ---------------------------------------------------------------------
   !
   ! Eleven melt-interface fields and two geometry fields.  Three rules
   ! run through all of them:
   !
   !   MISSING OUTSIDE THE CAVITY.  Every melt field is IEEE NaN where
   !   the column was not solved, the same sentinel `fill_tracer_impl`
   !   writes for a vanished layer, so a domain mean over the output is
   !   a mean over the CAVITY and not a cavity average diluted by a
   !   basin of zeros.  The console reduction already skips non-finite
   !   cells and reports how many it skipped, and the NetCDF writer tags
   !   the variable with `_FillValue` (`has_missing` above).  A zero
   !   would be much worse than absent here: zero melt is a LEGAL
   !   answer, and 0 m/yr averaged over open ocean is how a 30 m/yr
   !   shelf gets published as 2.
   !
   !   THE MASK IS `cavity_flux%active`, NOT `cover_frac`.  `active` is
   !   what the kernel was actually handed — covered AND wet AND "the
   !   far-field sample found mass" — so a covered column that was
   !   skipped reads missing rather than reading the zeros its outputs
   !   were left at.
   !
   !   NOTHING IS RE-DERIVED.  `exch_vel_t`/`exch_vel_s` come off the
   !   slot because under `hj99`/`yung25` they are implicit functions of
   !   the converged interface state; `tfreeze_ib` goes through the same
   !   `eos_freezing_point` handle the solve used.  A diagnostic that
   !   rebuilt either one would be a second, divergent implementation of
   !   the physics.

   pure subroutine cavity_mask_impl(src, active, buf)
      !! Copy a 2-D cavity field into the diag buffer, writing IEEE NaN
      !! wherever the column was not solved.  Every cavity fill that is
      !! a plain read of a slot array routes through here, so the
      !! missing-value convention cannot drift between them.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: src(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(src, 1), size(active, 1))
      ny = min(size(buf, 2), size(src, 2), size(active, 2))
      ! Sentinel computed once on the HOST before the device loop —
      ! `ieee_value` is a host intrinsic (the `fill_tracer_impl` pattern).
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = src(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cavity_mask_impl

   pure function melt_m_per_yr_factor() result(f)
      !! The ISOMIP+ melt-rate conversion, `kg m-2 s-1` → `m yr-1` of
      !! ice-equivalent freshwater: `f = SECONDS_PER_YEAR / rho_fw`.
      !!
      !! Asay-Davis et al. (2016), *GMD* **9**, 2471–2497 §3.3.  With
      !! `rho_fw = 1000 kg/m^3` and a 365-day year,
      !! `f = 31 536 000 / 1000 = 31 536 m yr-1 per kg m-2 s-1`, i.e. a
      !! melt flux of 1e-3 kg/m^2/s is 31.536 m/yr — hand-checkable, and
      !! the test checks it by hand.
      !!
      !! A function rather than a parameter so the test can assert the
      !! CONVERSION rather than re-typing the same two constants and
      !! asserting its own arithmetic.
      real(wp) :: f
      f = SECONDS_PER_YEAR/RHO_FRESHWATER_ISOMIP
   end function melt_m_per_yr_factor

   subroutine fill_melt(state_handle, buf)
      !! Basal melt mass flux (kg m-2 s-1), **positive = melting**.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%melt, state%cavity_flux%active, buf)
      end select
   end subroutine fill_melt

   subroutine fill_melt_m_per_yr(state_handle, buf)
      !! Basal melt rate in the ISOMIP+ reporting unit (m yr-1 of
      !! ice-equivalent freshwater) — see `melt_m_per_yr_factor`.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_scale_impl(state%cavity_flux%melt, state%cavity_flux%active, &
                                melt_m_per_yr_factor(), buf)
      end select
   end subroutine fill_melt_m_per_yr

   pure subroutine cavity_scale_impl(src, active, scale, buf)
      !! `cavity_mask_impl` with a constant multiplier folded in.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: src(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: scale
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(src, 1), size(active, 1))
      ny = min(size(buf, 2), size(src, 2), size(active, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = scale*src(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cavity_scale_impl

   subroutine fill_thermal_driving(state_handle, buf)
      !! `T* = T_far − T_f(S_far, p_top)` (degC) — the far-field
      !! temperature above the IN-SITU freezing point at the interface
      !! pressure.  Positive drives melting.  This is the single number
      !! the melt rate is roughly linear in, so it is the first thing to
      !! look at when a melt rate looks wrong.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call thermal_driving_impl(state%cavity_flux%t_far, state%cavity_flux%s_far, &
                                   state%multilayer%p_top, state%cavity_flux%active, &
                                   state%eos, buf)
      end select
   end subroutine fill_thermal_driving

   pure subroutine thermal_driving_impl(t_far, s_far, p_top, active, eos, buf)
      !! `T_far − eos_freezing_point(eos, S_far, p_top)`.  The liquidus
      !! comes off the SAME `eos_t` handle the solve used
      !! (`&ocean_eos_nml tfreeze_set`), never a local copy of the
      !! coefficients — the two sets differ by ~0.03 degC, which is
      !! enough to flip the sign of this field.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      use rdb_eos, only: eos_t, eos_freezing_point
      real(wp), intent(in)    :: t_far(:, :), s_far(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: p_top(:, :), active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(t_far, 1), size(p_top, 1), size(active, 1))
      ny = min(size(buf, 2), size(t_far, 2), size(p_top, 2), size(active, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = t_far(i, j) - eos_freezing_point(eos, s_far(i, j), p_top(i, j))
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine thermal_driving_impl

   subroutine fill_tfreeze_ib(state_handle, buf)
      !! `T_f(S_far, p_top)` (degC) — the in-situ freezing point of the
      !! FAR FIELD at the ice base.  Distinct from `tbdry`, which is the
      !! freezing point of the INTERFACE salinity `S_b`; their
      !! difference is the whole three-equation correction, so shipping
      !! both makes it visible.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call tfreeze_ib_impl(state%cavity_flux%s_far, state%multilayer%p_top, &
                              state%cavity_flux%active, state%eos, buf)
      end select
   end subroutine fill_tfreeze_ib

   pure subroutine tfreeze_ib_impl(s_far, p_top, active, eos, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      use rdb_eos, only: eos_t, eos_freezing_point
      real(wp), intent(in)    :: s_far(:, :), p_top(:, :), active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(s_far, 1), size(p_top, 1), size(active, 1))
      ny = min(size(buf, 2), size(s_far, 2), size(p_top, 2), size(active, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = eos_freezing_point(eos, s_far(i, j), p_top(i, j))
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine tfreeze_ib_impl

   subroutine fill_haline_driving(state_handle, buf)
      !! `S* = S_far − S_b` (g/kg) — the salinity contrast the haline
      !! exchange acts on.  Positive under melting (meltwater freshens
      !! the interface).
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_diff_impl(state%cavity_flux%s_far, state%cavity_flux%s_b, &
                               state%cavity_flux%active, buf)
      end select
   end subroutine fill_haline_driving

   pure subroutine cavity_diff_impl(a, b, active, buf)
      !! `a - b` with the cavity missing-value convention.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: a(:, :), b(:, :), active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(a, 1), size(b, 1), size(active, 1))
      ny = min(size(buf, 2), size(a, 2), size(b, 2), size(active, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = a(i, j) - b(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cavity_diff_impl

   subroutine fill_tbdry(state_handle, buf)
      !! Interface temperature `T_b` (degC) — on the liquidus at `S_b`.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%t_b, state%cavity_flux%active, buf)
      end select
   end subroutine fill_tbdry

   subroutine fill_sbdry(state_handle, buf)
      !! Interface salinity `S_b` (g/kg).
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%s_b, state%cavity_flux%active, buf)
      end select
   end subroutine fill_sbdry

   subroutine fill_exch_vel_t(state_handle, buf)
      !! Thermal exchange velocity `gamma_T` (m/s) the solve CONVERGED
      !! on — not `Gamma_T·u*` re-derived here.  Under `hj99`/`yung25`
      !! the two differ by the stratification suppression.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%gamma_t, state%cavity_flux%active, buf)
      end select
   end subroutine fill_exch_vel_t

   subroutine fill_exch_vel_s(state_handle, buf)
      !! Haline exchange velocity `gamma_S` (m/s) — see `fill_exch_vel_t`.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%gamma_s, state%cavity_flux%active, buf)
      end select
   end subroutine fill_exch_vel_s

   subroutine fill_ustar_shelf(state_handle, buf)
      !! Melt friction velocity `u*` (m/s), `max(sqrt(C_d(u²+v²+u_tide²)),
      !! u*_min)`.  This is the MELT law's `u*` and nothing else: it
      !! drives no momentum drag, and KPP/EPBL do not read it.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%cavity_flux%ustar, state%cavity_flux%active, buf)
      end select
   end subroutine fill_ustar_shelf

   subroutine fill_cavity_melt_status(state_handle, buf)
      !! Per-column `CAVITY_MELT_*` status code, as a real.  `0` is OK;
      !! anything else is a column that took the kernel's documented
      !! zero-melt safe state, and the console's warning line says how
      !! many there were.  Masked like the rest, so the plane cannot be
      !! read as "everything is fine" over open ocean.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_status_impl(state%cavity_flux%status, state%cavity_flux%active, buf)
      end select
   end subroutine fill_cavity_melt_status

   pure subroutine cavity_status_impl(status, active, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      integer, intent(in)     :: status(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: active(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(status, 1), size(active, 1))
      ny = min(size(buf, 2), size(status, 2), size(active, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            buf(i, j, 1) = real(status(i, j), wp)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine cavity_status_impl

   subroutine fill_z_draft(state_handle, buf)
      !! Ice-shelf draft (m, POSITIVE DOWN from the geoid) — geometry,
      !! so it needs only `&ocean_cavity_dyn_nml` and is masked by
      !! `cover_frac` rather than by the melt slot's `active`.  Beyond
      !! the calving front there is no draft, which is missing data, not
      !! a draft of zero.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call cavity_mask_impl(state%metrics%z_draft, state%metrics%cover_frac, buf)
      end select
   end subroutine fill_z_draft

   subroutine fill_water_column(state_handle, buf)
      !! Water-column thickness under the shelf, `D = bt_H_ref + bt_eta`
      !! (m) — the ONE expression every consumer of the column depth
      !! uses, which is exactly what the cavity datum
      !! (`bt_H_ref = b − z_draft`) buys.  Masked by `cover_frac`: the
      !! quantity is perfectly well defined in open water, but as a
      !! CAVITY diagnostic a domain mean should be the mean cavity
      !! thickness, not that diluted by the open basin.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call water_column_impl(state%dyn%bt_work%bt_H_ref, state%dyn%bt_work%bt_eta, &
                                state%metrics%cover_frac, buf)
      end select
   end subroutine fill_water_column

   pure subroutine water_column_impl(bt_H_ref, bt_eta, cover, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: bt_H_ref(:, :), bt_eta(:, :), cover(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      real(wp) :: qnan
      nx = min(size(buf, 1), size(bt_H_ref, 1), size(bt_eta, 1), size(cover, 1))
      ny = min(size(buf, 2), size(bt_H_ref, 2), size(bt_eta, 2), size(cover, 2))
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(j=1:ny, i=1:nx)
         if (cover(i, j) > 0.5_wp) then
            buf(i, j, 1) = bt_H_ref(i, j) + bt_eta(i, j)
         else
            buf(i, j, 1) = qnan
         end if
      end do
   end subroutine water_column_impl

   pure subroutine fill_ice_v_impl(v_ice, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized v_ice(:,:) has ny+1 second dim, incompatible with cell-sized buf.
      real(wp), intent(in)    :: v_ice(:, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      nx = min(size(buf, 1), size(v_ice, 1))
      ny = min(size(buf, 2), size(v_ice, 2) - 1)
      do concurrent(j=1:ny, i=1:nx)
         buf(i, j, 1) = 0.5_wp*(v_ice(i, j) + v_ice(i, j + 1))
      end do
   end subroutine fill_ice_v_impl

end module rdb_ocean_diag_derived
