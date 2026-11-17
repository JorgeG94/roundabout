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
   !! The thermodynamic flux tags (`heat`, `salt`, `evap`, `lprec`) are
   !! applied column-locally, but the boundary-layer schemes read the
   !! ASSEMBLED `Q_heat`/`Q_salt` in the GHOST columns (their loops span
   !! the full array) and the momentum vertical viscosity then averages the
   !! ghost-column `kv` into the first owned face.  So the seam ghosts of the
   !! flux the schemes read are exchanged: under `use_components` the
   !! engine does it once per thermo step after the assembler, and with components off, where the
   !! tags write `Q_heat`/`Q_salt` themselves, this module exchanges them
   !! right after the write.
   !!
   !! **Runoff** (`lrunoff`, `frunoff`; kg m-2 s-1, into the ocean).  The files must already be
   !! on the model grid with every value on a WET cell (`tools/om_runoff_to_model.f90`
   !! builds them: conservative remap of a lon-lat product, then land values moved to the
   !! nearest wet coastal cell, global total preserved; Tsujino et al. 2018 for the
   !! JRA55-do runoff).  They land in `surface_flux%lrunoff/frunoff`, which the assembler
   !! folds into `fw_net` (real volume under `freshwater = "mass"`) or into the virtual
   !! salt flux (`"virtual"`); `"none"` refuses them.  `frunoff` also writes
   !! `q_fus_rof = -L_f*frunoff` (this routine is its only writer).  Runoff enthalpy follows
   !! the mass step's layer-temperature rule, so no `heat_content_*` is written.  A non-zero
   !! value on a land cell is FATAL (`check_runoff_on_land`, rank-uniform; run at the first
   !! read and whenever a new record has been read).
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
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_grid_reduce, only: grid_allreduce_sum_i8
   use rdb_ice_enthalpy, only: ICE_LAT_FUS
   use rdb_config, only: ocean_dataovr_config_t, dataovr_entry_config_t, &
                         ocean_restore_config_t
   use rdb_ocean_data_input, only: ocean_data_input_t, &
                                   ocean_data_input_register_2d, &
                                   ocean_data_input_update_2d, &
                                   DATA_OOR_ERROR, DATA_OOR_CLAMP
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress, &
                                   ocean_halo_exchange_centre_field
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_data_forcing_t
   public :: ocean_data_forcing_configure
   public :: ocean_data_forcing_register_sss
   public :: ocean_data_forcing_apply
   public :: runoff_on_land_is_fatal
   public :: runoff_land_count_impl
   public :: runoff_latent_heat_impl

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
      integer :: id_lrunoff = 0
         !! `lrunoff` (liquid runoff, kg m-2 s-1) -> `surface_flux%lrunoff`.
      integer :: id_frunoff = 0
         !! `frunoff` (frozen runoff / calving) -> `surface_flux%frunoff`, and
         !! `q_fus_rof = -L_f*frunoff` is derived from it each step.
      integer :: land_checked_l = -1
         !! `nreads` of the `lrunoff` field when `check_runoff_on_land` last ran (`-1` = never).
      integer :: land_checked_f = -1
         !! Same for `frunoff`.  The check runs only when a new record has been read.
      integer :: id_sss = 0
         !! File-based SSS restoring target (`&ocean_restore_nml sss_file`)
         !! -> `surface_flux%sss_target`; 0 = scalar target.  Registered by
         !! `ocean_data_forcing_register_sss`, independent of
         !! `&ocean_dataovr_nml enable`.
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

      call register_tag(cfg, cfg%lrunoff, "lrunoff", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_lrunoff, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call register_tag(cfg, cfg%frunoff, "frunoff", reader, grid, oor, &
                        nx, ny, ng + 1, ng + 1, this%id_frunoff, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if

      this%active = (this%id_tau_x > 0 .or. this%id_tau_y > 0 .or. &
                     this%id_heat > 0 .or. this%id_evap > 0 .or. &
                     this%id_lprec > 0 .or. this%id_salt > 0 .or. &
                     this%id_lrunoff > 0 .or. this%id_frunoff > 0)

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
      if (this%id_evap > 0 .or. this%id_lprec > 0 .or. &
          this%id_lrunoff > 0 .or. this%id_frunoff > 0) sf%has_mass_flux = .true.
      ! `stress_mag` must be consistent with whatever the first blend
      ! writes; the first `apply` re-derives it, but configure-time
      ! consumers that read it before the first step should not see a
      ! value derived from a wind field that is about to be replaced.
      if (this%id_tau_x > 0 .or. this%id_tau_y > 0) then
         call ocean_seam_refresh_surface_stress(ss, grid, bc, device_resident=.false.)
      end if
   end subroutine ocean_data_forcing_configure

   subroutine ocean_data_forcing_register_sss(rcfg, reader, grid, sf, this, ierr)
      !! Register the SSS-restoring target file (`&ocean_restore_nml
      !! sss_file`) against the shared reader.  Independent of
      !! `&ocean_dataovr_nml enable`; call AFTER `ocean_data_forcing_configure`
      !! (which resets `this`) and after `surface_flux%set_restore_file`
      !! (which allocates the target).  A blank `sss_file` is a no-op.
      !! The reader fails loud on a missing file, variable or grid mismatch.
      type(ocean_restore_config_t), intent(in) :: rcfg
      type(ocean_data_input_t), intent(inout) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in) :: sf
      type(ocean_data_forcing_t), intent(inout) :: this
      integer, intent(out), optional :: ierr

      character(len=8) :: tmode

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (len_trim(rcfg%sss_file) == 0) return
      if (.not. allocated(sf%sss_target)) then
         call fail("ocean_restore: sss_file is set but surface_flux%set_restore_file "// &
                   "has not run (no target field)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      tmode = "linear"
      if (rcfg%sss_time_cyclic) tmode = "cyclic"
      call ocean_data_input_register_2d(reader, trim(rcfg%sss_file), trim(rcfg%sss_var), &
                                        grid, grid%nx_total, grid%ny_total, &
                                        grid%nghost + 1, grid%nghost + 1, this%id_sss, &
                                        time_mode=trim(tmode), &
                                        cycle_period=rcfg%sss_cycle_period, &
                                        oor=DATA_OOR_ERROR, ierr=ierr)
      if (present(ierr)) then
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      this%active = .true.
      call logger%info("ocean_restore: SSS target <- '"//trim(rcfg%sss_var)//"' in "// &
                       trim(rcfg%sss_file)//" (mode "//trim(tmode)//", id "// &
                       to_string(this%id_sss)//")")
   end subroutine ocean_data_forcing_register_sss

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

      if (len_trim(cfg%evap%file) > 0 .or. len_trim(cfg%lprec%file) > 0 .or. &
          len_trim(cfg%lrunoff%file) > 0 .or. len_trim(cfg%frunoff%file) > 0) then
         call fail("&ocean_dataovr_nml: evap/lprec/lrunoff/frunoff file forcing needs the surface-flux "// &
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

   subroutine ocean_data_forcing_apply(this, reader, grid, ss, sf, bc, t, wet_mask)
      !! Blend every active tag's current bracket into its slot.  `t` must
      !! be the same model time `ocean_data_input_update_all` was just
      !! called with — the reader enforces this and aborts otherwise.
      !!
      !! Only PHYSICAL cells are written here.  Stress ghosts are filled
      !! afterwards by `ocean_seam_refresh_surface_stress` (exchange /
      !! periodic wrap / fold), which also re-derives `stress_mag`.  The
      !! flux tags' ghosts are refreshed too, because KPP/EPBL read the
      !! assembled flux in the ghost columns: with components off the
      !! `Q_heat`/`Q_salt` this routine writes are exchanged right after
      !! the write; with components on the tags fill components and the
      !! engine exchanges the ASSEMBLED `Q_heat`/`Q_salt` after the
      !! assembler (whenever the component set is on).
      type(ocean_data_forcing_t), intent(inout) :: this
         !! `inout` only for the runoff land-check read counters.
      type(ocean_data_input_t), intent(in) :: reader
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_bc_state_t), intent(in) :: bc
         !! Supplies the periodic / north-fold topology to the seam refresh.
      real(wp), intent(in) :: t
      real(wp), intent(in), optional :: wet_mask(grid%nx_total, grid%ny_total)
         !! `multilayer_state_t%wet_mask`.  Present, with a runoff tag registered, it
         !! arms the land-cell check (see below); the engine always passes it.

      integer :: nx, ny, n_rl, n_rf

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

      ! --- SSS-restoring target: the deposit kernel runs over the full
      ! array, so the seam ghosts are EXCHANGED (never extrapolated) ---
      if (this%id_sss > 0) then
         call ocean_data_input_update_2d(reader, this%id_sss, t, nx, ny, sf%sss_target)
         call ocean_halo_exchange_centre_field(sf%sss_target, grid, bc)
      end if

      ! --- surface fluxes (component destinations: ghosts refreshed by the
      ! engine after the assembler; direct Q_heat/Q_salt writes: exchanged here) ---
      if (this%id_heat > 0) then
         if (this%heat_to_component) then
            call ocean_data_input_update_2d(reader, this%id_heat, t, nx, ny, sf%heat_added)
         else
            call ocean_data_input_update_2d(reader, this%id_heat, t, nx, ny, sf%Q_heat)
            ! Written on physical cells only and read in the ghost columns
            ! by KPP/EPBL: exchange (never extrapolate).
            call ocean_halo_exchange_centre_field(sf%Q_heat, grid, bc)
         end if
      end if
      if (this%id_salt > 0) then
         if (this%salt_to_component) then
            call ocean_data_input_update_2d(reader, this%id_salt, t, nx, ny, sf%salt_flux)
         else
            call ocean_data_input_update_2d(reader, this%id_salt, t, nx, ny, sf%Q_salt)
            call ocean_halo_exchange_centre_field(sf%Q_salt, grid, bc)
         end if
      end if
      if (this%id_evap > 0) then
         call ocean_data_input_update_2d(reader, this%id_evap, t, nx, ny, sf%evap)
      end if
      if (this%id_lprec > 0) then
         call ocean_data_input_update_2d(reader, this%id_lprec, t, nx, ny, sf%lprec)
      end if

      ! --- runoff: model-grid mass fluxes (kg m-2 s-1, into the ocean).  Read on physical
      ! cells; the assembler snapshots them into `fw_net` (whose ghosts are exchanged), so no
      ! halo work belongs here.  `frunoff` also owes the latent heat of fusion `q_fus_rof`
      ! (its ONE owner is this routine).
      if (this%id_lrunoff > 0) then
         call ocean_data_input_update_2d(reader, this%id_lrunoff, t, nx, ny, sf%lrunoff)
      end if
      if (this%id_frunoff > 0) then
         call ocean_data_input_update_2d(reader, this%id_frunoff, t, nx, ny, sf%frunoff)
         call runoff_latent_heat_impl(nx, ny, ICE_LAT_FUS, sf%frunoff, sf%q_fus_rof)
      end if
      ! The land check reads the blended slot, whose land-ness is that of the two bracketing
      ! records, so it runs only when a record has been (re)read since the last check -- and
      ! once at the start (counters at -1).  The read counts are identical on every rank, so
      ! every rank takes the same branch and the collective inside stays matched.
      if ((this%id_lrunoff > 0 .or. this%id_frunoff > 0) .and. present(wet_mask)) then
         n_rl = -2
         n_rf = -2
         if (this%id_lrunoff > 0) n_rl = reader%fields(this%id_lrunoff)%nreads
         if (this%id_frunoff > 0) n_rf = reader%fields(this%id_frunoff)%nreads
         if (n_rl /= this%land_checked_l .or. n_rf /= this%land_checked_f) then
            call check_runoff_on_land(this, grid, sf, wet_mask)
            this%land_checked_l = n_rl
            this%land_checked_f = n_rf
         end if
      end if
   end subroutine ocean_data_forcing_apply

   subroutine check_runoff_on_land(this, grid, sf, wet_mask)
      !! Fatal, rank-uniform check, run at the first read and whenever a new record has been read: a runoff value on a LAND cell would vanish
      !! silently (`fw_net` is wet-masked, so the mass never enters the ocean and the global
      !! mass budget still "closes" against a source that was never spent).  The offline tool
      !! `tools/om_runoff_to_model.f90` moves such values to the nearest wet coastal cell;
      !! a file that skipped it is a setup error, not something to clip silently.  Counts
      !! over OWNED cells and sums once (int64), so every rank takes the same branch.
      type(ocean_data_forcing_t), intent(in) :: this
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(in) :: sf
      real(wp), intent(in) :: wet_mask(grid%nx_total, grid%ny_total)
      integer :: n_loc, n_f
      integer(int64) :: n_glob

      n_loc = 0
      if (this%id_lrunoff > 0) then
         call runoff_land_count_impl(grid%nx_total, grid%ny_total, grid%nghost, wet_mask, &
                                     sf%lrunoff, n_loc)
      end if
      if (this%id_frunoff > 0) then
         call runoff_land_count_impl(grid%nx_total, grid%ny_total, grid%nghost, wet_mask, &
                                     sf%frunoff, n_f)
         n_loc = n_loc + n_f
      end if
      call grid_allreduce_sum_i8(grid, int(n_loc, int64), n_glob)
      if (runoff_on_land_is_fatal(n_glob)) then
         call logger%error("ocean_dataovr: "//to_string(int(n_glob))// &
                           " land cell(s) carry a non-zero (or non-finite) lrunoff/frunoff")
         call fail("&ocean_dataovr_nml lrunoff/frunoff: "//to_string(int(n_glob))// &
                   " land cell(s) carry runoff.  That mass would vanish silently (the "// &
                   "assembled fresh-water flux is wet-masked).  Move it to the nearest wet "// &
                   "coastal cell offline with tools/om_runoff_to_model.f90.", &
                   code=OCEAN_STATUS_ERR_SETUP)
      end if
   end subroutine check_runoff_on_land

   pure function runoff_on_land_is_fatal(n_land) result(is_fatal)
      !! The verdict on the global count of land cells carrying runoff: any is fatal.
      integer(int64), intent(in) :: n_land
      logical :: is_fatal
      is_fatal = (n_land > 0_int64)
   end function runoff_on_land_is_fatal

   pure subroutine runoff_land_count_impl(nx, ny, nghost, wet_mask, runoff, n_land)
      !! Number of OWNED land cells (`wet_mask <= 0.5`) whose `runoff` is non-zero or
      !! non-finite.  One accumulator, no inner loop (reduction-kernel rule).
      integer, intent(in) :: nx, ny, nghost
      real(wp), intent(in) :: wet_mask(nx, ny), runoff(nx, ny)
      integer, intent(out) :: n_land
      integer :: i, j, c

      c = 0
      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) reduce(+:c)
         if (wet_mask(i, j) <= 0.5_wp) then
            if (.not. ieee_is_finite(runoff(i, j))) then
               c = c + 1
            else if (runoff(i, j) /= 0.0_wp) then
               c = c + 1
            end if
         end if
      end do
      n_land = c
   end subroutine runoff_land_count_impl

   pure subroutine runoff_latent_heat_impl(nx, ny, lfus, frunoff, q_fus_rof)
      !! `q_fus_rof = -L_f * frunoff` (W m-2, <= 0, positive down): the latent heat owed
      !! by frozen runoff / calving, over the whole array (the ghosts hold zero runoff).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: lfus
      real(wp), intent(in) :: frunoff(nx, ny)
      real(wp), intent(out) :: q_fus_rof(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         q_fus_rof(i, j) = -lfus*frunoff(i, j)
      end do
   end subroutine runoff_latent_heat_impl

end module rdb_ocean_data_forcing
