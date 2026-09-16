!! File-backed surface forcing (PR-15) — the first production consumer
!! of the shared PR-14 NetCDF reader.
module rdb_ocean_data_forcing
   !! Binds time-varying NetCDF surface fields onto the ocean C-grid
   !! forcing slots: wind stress (`ocean_surface_stress_t%tau_x/tau_y`)
   !! and the surface heat / freshwater / salt fluxes
   !! (`ocean_surface_flux_t`).  Owns the `(file, variable) -> slot`
   !! mapping that `rdb_ocean_data_input` deliberately does not; owns no
   !! file handles, no time axis and no buffers of its own — those all
   !! live in the reader.
   !!
   !! **Two-call contract.**  `ocean_data_forcing_configure` runs once in
   !! the driver's configure phase, BEFORE `ocean_state_enter_data`
   !! (registration opens each file and allocates the reader's `f0`/`f1`
   !! bracket buffers, which `enter_data` then maps).
   !! `ocean_data_forcing_apply` runs once per outer step, immediately
   !! after the driver's `ocean_data_input_update_all(t_current)` — the
   !! reader's fail-loud freshness check enforces that ordering rather
   !! than letting a stale bracket blend silently.
   !!
   !! **This slot maps nothing.**  It holds registration `id`s and
   !! logicals that are read HOST-side only, so it adds no term to
   !! `ocean_state_enter_data` — stated explicitly because the standing
   !! rule is that a new ocean slot *does* wire in, and a silent omission
   !! there is a 150-1500x memcpy bug.  The arrays it writes into are
   !! mapped by their own owning slots.
   !!
   !! **Ghost cells are EXCHANGED, never extrapolated.**  The reader
   !! fills physical cells only and the consumer owns the halo.  For the
   !! stress pair that halo is filled by
   !! `ocean_seam_refresh_surface_stress` — MPI exchange, then periodic
   !! wrap on any axis the exchange did not own, then the tripolar fold —
   !! which also re-derives `stress_mag` (read by KPP/EPBL for `u_*`,
   !! and otherwise stale from the configure-time wind).
   !!
   !! An earlier revision zero-gradient-extended the blend into its
   !! ghosts.  That is wrong under decomposition and worth recording so
   !! it is not reintroduced: at an MPI seam the ghost belongs to the
   !! neighbour rank, so copying this rank's edge value there produces a
   !! decomposition-dependent answer that no single-rank test can see.
   !! Extrapolating forcing into a halo is a pattern MOM6 does not use
   !! anywhere; it exchanges the stress pair and relies on masking at
   !! true domain edges.
   !!
   !! The thermodynamic flux tags (`heat`, `salt`, `evap`, `lprec`) get
   !! NO ghost treatment at all, deliberately: they are applied strictly
   !! column-locally, so no kernel ever reads them in a ghost cell.  A
   !! future kernel that takes a horizontal gradient of one of them owns
   !! adding the exchange.
   !!
   !! **The namelist is the expert-level surface, not the intended one.**
   !! `&ocean_dataovr_nml` is a flat, fixed set of tags because that is
   !! what the strict schema can express; it is not the ergonomic way to
   !! describe a forcing dataset.  The binding itself is programmatic —
   !! `register_tag` is a thin wrapper over
   !! `ocean_data_input_register_2d` returning an opaque id — so a
   !! future Python/C driver should call the registration path DIRECTLY
   !! with its own field table rather than synthesising namelist text.
   !! Keep `ocean_data_forcing_configure` a pure translation of config
   !! to registrations, with no logic that a non-namelist caller would
   !! have to re-implement.
   !!
   !! **Heat/salt destination depends on `use_components`** — see
   !! `resolve_flux_targets`.  `ocean_surface_flux_assemble` fully
   !! overwrites `Q_heat`/`Q_salt` from the component set every thermo
   !! step, so under components the file must feed a component instead.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: ocean_dataovr_config_t, dataovr_entry_config_t
   use rdb_ocean_data_input, only: ocean_data_input_t, &
                                   ocean_data_input_register_2d, &
                                   ocean_data_input_update_2d, &
                                   DATA_OOR_ERROR, DATA_OOR_CLAMP
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_data_forcing_t
   public :: ocean_data_forcing_configure
   public :: ocean_data_forcing_apply

   type :: ocean_data_forcing_t
      !! Registration bookkeeping for the file-driven surface tags.  A
      !! zero `id_*` means that tag is not file-driven (blank `file` in
      !! `&ocean_dataovr_nml`) and its slot keeps whatever the
      !! configure-time scalar/formula path seeded.
      logical :: active = .false.
         !! `.true.` once at least one tag registered.  `.false.` makes
         !! `ocean_data_forcing_apply` an immediate return.
      integer :: id_tau_x = 0
      integer :: id_tau_y = 0
      integer :: id_heat = 0
      integer :: id_evap = 0
      integer :: id_lprec = 0
      integer :: id_salt = 0
      logical :: heat_to_component = .false.
         !! `.true.` => the heat tag blends into `heat_added` (components
         !! on); `.false.` => straight into `Q_heat`.  Resolved once at
         !! configure so the per-step path carries no mode test.
      logical :: salt_to_component = .false.
         !! Same for the salt tag: `salt_flux` vs `Q_salt`.
   end type ocean_data_forcing_t

contains

   ! ======================================================================
   ! Configure (once, before ocean_state_enter_data)
   ! ======================================================================

   subroutine ocean_data_forcing_configure(cfg, reader, grid, ss, sf, bc, this, ierr)
      !! Register every configured tag against the shared reader and
      !! resolve the heat/salt destinations.  A no-op when
      !! `enable = .false.` — nothing registers, so
      !! `ocean_data_input_update_all` stays a no-op and the run is
      !! bit-identical.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_bc_state_t), intent(in) :: bc
         !! Periodic / north-fold topology for the configure-time seam refresh.
      type(ocean_data_forcing_t), intent(out) :: this
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`/`OCEAN_STATUS_ERR_IO`, see
         !! `rdb_ocean_status`) on a malformed `&ocean_dataovr_nml` group
         !! or a registration failure (bad file/variable/dims) when
         !! present; absent behaves as today (`error stop`).

      integer :: oor, nx, ny, ng

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. cfg%enable) return

      call check_time_mode(cfg, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      oor = DATA_OOR_ERROR
      if (cfg%oor_clamp) oor = DATA_OOR_CLAMP

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost

      ! --- wind stress: face-shaped destinations ---
      ! `tau_x` is (nx_total+1, ny_total) on east faces, `tau_y` is
      ! (nx_total, ny_total+1) on north faces, landing at (nghost+1,
      ! nghost+1).
      !
      ! The file must supply nx_phys+1 x ny_phys values for `tau_x`
      ! (ny_phys+1 in y for `tau_y`) — one MORE than the cell count in
      ! the staggered direction.  An x-face array spans nx_phys+1 faces,
      ! and the trailing face is owned by THIS rank under the
      ! west/south-owns-the-seam rule, so no exchange can fill it in.
      ! Supplying only nx_phys would silently leave that face at zero,
      ! which halves the cell-centred `stress_mag` in the last column
      ! and hands a zero to the east neighbour's seam.
      call register_tag(cfg, cfg%tau_x, "tau_x", reader, grid, oor, &
                        nx + 1, ny, ng + 1, ng + 1, this%id_tau_x, nx_extra=1, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_tag(cfg, cfg%tau_y, "tau_y", reader, grid, oor, &
                        nx, ny + 1, ng + 1, ng + 1, this%id_tau_y, ny_extra=1, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      ! --- surface fluxes: cell-centred destinations ---
      call resolve_flux_targets(cfg, sf, this, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      call register_tag(cfg, cfg%heat, "heat", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_heat, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_tag(cfg, cfg%evap, "evap", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_evap, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_tag(cfg, cfg%lprec, "lprec", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_lprec, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_tag(cfg, cfg%salt, "salt", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_salt, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      this%active = (this%id_tau_x > 0 .or. this%id_tau_y > 0 .or. &
                     this%id_heat > 0 .or. this%id_evap > 0 .or. &
                     this%id_lprec > 0 .or. this%id_salt > 0)

      if (.not. this%active) then
         call logger%warning("ocean_dataovr: enable = .true. but no tag has a file — "// &
                             "every surface slot keeps its configure-time value")
         return
      end if

      ! Latches the consuming kernels gate on.  Set HERE, not per step:
      ! they are host-side gates on a device-mapped derived type, and a
      ! host write to a mapped component never reaches the device.
      ! Registration is for the whole run, so the latch is too.
      if (this%id_heat > 0) sf%has_heat = .true.
      if (this%id_salt > 0) sf%has_salt = .true.
      if (this%id_evap > 0 .or. this%id_lprec > 0) sf%has_mass_flux = .true.

      ! `stress_mag` must be consistent with whatever the first blend
      ! writes; the first `apply` re-derives it, but configure-time
      ! consumers that read it before the first step should not see a
      ! value derived from a wind field that is about to be replaced.
      if (this%id_tau_x > 0 .or. this%id_tau_y > 0) then
         call ocean_seam_refresh_surface_stress(ss, grid, bc, device_resident=.false.)
      end if
   end subroutine ocean_data_forcing_configure

   subroutine check_time_mode(cfg, ierr)
      !! Fail loud on a `cyclic` group with no period.  The reader makes
      !! the same check per field; catching it once here names the
      !! namelist group the user actually edited.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (trim(cfg%time_mode) == "cyclic" .and. cfg%cycle_period <= 0.0_wp) then
         call fail("&ocean_dataovr_nml: time_mode = 'cyclic' requires cycle_period > 0", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
      end if
   end subroutine check_time_mode

   subroutine resolve_flux_targets(cfg, sf, this, ierr)
      !! Pick the heat/salt destination for this run and reject the
      !! freshwater tags when the component set they need is absent.
      !!
      !! `ocean_surface_flux_assemble` rebuilds `Q_heat`/`Q_salt` from
      !! the components on every thermo step when `use_components` is
      !! on, so a file write straight into `Q_heat` would be overwritten
      !! before it was ever read.  With components off there is no
      !! component to write, so `Q_heat` is both correct and the only
      !! option.  `evap`/`lprec` exist only in the component set.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      type(ocean_surface_flux_t), intent(in) :: sf
      type(ocean_data_forcing_t), intent(inout) :: this
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      this%heat_to_component = sf%use_components
      this%salt_to_component = sf%use_components

      if (sf%use_components) return

      if (len_trim(cfg%evap%file) > 0 .or. len_trim(cfg%lprec%file) > 0) then
         call fail("&ocean_dataovr_nml: evap/lprec file forcing needs the surface-flux "// &
                   "component set — set &ocean_forcing_nml enable_components = .true. "// &
                   "(there is no non-component freshwater slot to write into)", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
      end if
   end subroutine resolve_flux_targets

   subroutine register_tag(cfg, e, tag, reader, grid, oor, dn1, dn2, di0, dj0, id, &
                           nx_extra, ny_extra, ierr)
      !! Register one tag if it names a file, else leave `id = 0`.
      !! Shared time-axis settings come from the group; `scale`/`add`
      !! are per tag.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      type(dataovr_entry_config_t), intent(in) :: e
      character(len=*), intent(in) :: tag
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: oor, dn1, dn2, di0, dj0
      integer, intent(out) :: id
      integer, intent(in), optional :: nx_extra, ny_extra
         !! Face-field slab widening; forwarded to the reader.
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      id = 0
      if (len_trim(e%file) == 0) return

      if (len_trim(e%var) == 0) then
         call fail("&ocean_dataovr_nml: "//tag//"_file is set ('"//trim(e%file)// &
                   "') but "//tag//"_var is blank — name the variable explicitly "// &
                   "(no guessed fallbacks for forcing)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      call ocean_data_input_register_2d(reader, trim(e%file), trim(e%var), grid, &
                                        dn1, dn2, di0, dj0, id, &
                                        time_mode=trim(cfg%time_mode), &
                                        cycle_period=cfg%cycle_period, &
                                        t_offset=cfg%t_offset, &
                                        scale=e%scale, add_offset=e%add_offset, &
                                        oor=oor, nx_extra=nx_extra, ny_extra=ny_extra, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      call logger%info("ocean_dataovr: "//tag//" <- '"//trim(e%var)//"' in "//trim(e%file)// &
                       " (mode "//trim(cfg%time_mode)//", id "//to_string(id)//")")
   end subroutine register_tag

   ! ======================================================================
   ! Per-step apply (immediately after ocean_data_input_update_all)
   ! ======================================================================

   subroutine ocean_data_forcing_apply(this, reader, grid, ss, sf, bc, t)
      !! Blend every active tag's current bracket into its slot.  `t` must
      !! be the same model time `ocean_data_input_update_all` was just
      !! called with — the reader enforces this and aborts otherwise.
      !!
      !! Only PHYSICAL cells are written here.  Stress ghosts are filled
      !! afterwards by `ocean_seam_refresh_surface_stress` (exchange /
      !! periodic wrap / fold), which also re-derives `stress_mag`.  The
      !! flux tags get no ghost treatment at all — see the module
      !! docstring for why that is correct rather than an omission.
      type(ocean_data_forcing_t), intent(in) :: this
      type(ocean_data_input_t), intent(in) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_bc_state_t), intent(in) :: bc
         !! Supplies the periodic / north-fold topology to the seam refresh.
      real(wp), intent(in) :: t

      integer :: nx, ny

      if (.not. this%active) return

      nx = grid%nx_total
      ny = grid%ny_total

      ! --- wind stress ---
      if (this%id_tau_x > 0) then
         call ocean_data_input_update_2d(reader, this%id_tau_x, t, nx + 1, ny, ss%tau_x)
      end if
      if (this%id_tau_y > 0) then
         call ocean_data_input_update_2d(reader, this%id_tau_y, t, nx, ny + 1, ss%tau_y)
      end if
      if (this%id_tau_x > 0 .or. this%id_tau_y > 0) then
         ! Mandatory, both halves: the ghosts because neighbour-reading
         ! kernels consume them, and `stress_mag` because KPP/EPBL would
         ! otherwise keep mixing on the configure-time wind.
         call ocean_seam_refresh_surface_stress(ss, grid, bc)
      end if

      ! --- surface fluxes (column-local; no ghost fill by design) ---
      if (this%id_heat > 0) then
         if (this%heat_to_component) then
            call ocean_data_input_update_2d(reader, this%id_heat, t, nx, ny, sf%heat_added)
         else
            call ocean_data_input_update_2d(reader, this%id_heat, t, nx, ny, sf%Q_heat)
         end if
      end if
      if (this%id_salt > 0) then
         if (this%salt_to_component) then
            call ocean_data_input_update_2d(reader, this%id_salt, t, nx, ny, sf%salt_flux)
         else
            call ocean_data_input_update_2d(reader, this%id_salt, t, nx, ny, sf%Q_salt)
         end if
      end if
      if (this%id_evap > 0) then
         call ocean_data_input_update_2d(reader, this%id_evap, t, nx, ny, sf%evap)
      end if
      if (this%id_lprec > 0) then
         call ocean_data_input_update_2d(reader, this%id_lprec, t, nx, ny, sf%lprec)
      end if
   end subroutine ocean_data_forcing_apply

end module rdb_ocean_data_forcing
