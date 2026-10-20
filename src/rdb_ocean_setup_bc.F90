submodule(rdb_ocean_setup) rdb_ocean_setup_bc
   !! R2 split: bc-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
   implicit none

contains

   module subroutine configure_ocean_bc(cfg, ocean_state, compute_rank, ierr)
      !! Populate `ocean_state%bc` edge tags + per-edge data values from
      !! `cfg%ocean%bc` (read from `&ocean_bc_nml`).  Run this after
      !! `configure_ocean_bt_split` and before `ocean_state_enter_data`
      !! so the BC state is set before the first dyn step.
      !!
      !! All defaults in `ocean_bc_config_t` resolve to OBC_WALL, so
      !! namelist files that omit `&ocean_bc_nml` are bit-identical to
      !! prior behaviour.
      !!
      !! After populating tags, calls `ocean_bc_validate_periodic` when any
      !! axis is periodic so the ghost-width and sponge incompatibility
      !! checks fire at setup time.
      type(config_t), intent(inout) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a boundary-condition configuration conflict when
         !! present; absent behaves as today (`error stop`).

      integer :: ntidal, n_tracers
      integer :: i, local_ierr

      associate (bc => ocean_state%bc, bc_cfg => cfg%ocean%bc)

         ! ---- Per-edge tag ----
         bc%west%bc_type = ocean_bc_type_from_string(bc_cfg%west)
         bc%east%bc_type = ocean_bc_type_from_string(bc_cfg%east)
         bc%south%bc_type = ocean_bc_type_from_string(bc_cfg%south)
         bc%north%bc_type = ocean_bc_type_from_string(bc_cfg%north)

         ! ---- Clamped Dirichlet values ----
         bc%west%clamped_eta = bc_cfg%west_clamped_eta
         bc%east%clamped_eta = bc_cfg%east_clamped_eta
         bc%south%clamped_eta = bc_cfg%south_clamped_eta
         bc%north%clamped_eta = bc_cfg%north_clamped_eta
         bc%west%clamped_u = bc_cfg%west_clamped_u
         bc%east%clamped_u = bc_cfg%east_clamped_u
         bc%south%clamped_v = bc_cfg%south_clamped_v
         bc%north%clamped_v = bc_cfg%north_clamped_v

         ! ---- Inflow tracer concentrations ----
         ! Populate clamped_tracer(idx_S) and clamped_tracer(idx_T) for
         ! each edge.  The bc_state_init has already allocated these arrays.
         n_tracers = bc%n_tracers
         if (n_tracers > 0 .and. allocated(bc%west%clamped_tracer)) then
            if (ocean_state%multilayer%idx_salinity > 0 .and. &
                ocean_state%multilayer%idx_salinity <= n_tracers) then
               i = ocean_state%multilayer%idx_salinity
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_S
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_S
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_S
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_S
            end if
            if (ocean_state%multilayer%idx_temperature > 0 .and. &
                ocean_state%multilayer%idx_temperature <= n_tracers) then
               i = ocean_state%multilayer%idx_temperature
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_T
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_T
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_T
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_T
            end if
            ! Pseudo-salt OBC mirror: same clamped inflow value as
            ! salinity (pseudo-salt receives exactly salinity's boundary
            ! fluxes, §3.2 of the plan).
            if (ocean_state%multilayer%idx_pseudo_salt > 0 .and. &
                ocean_state%multilayer%idx_pseudo_salt <= n_tracers) then
               i = ocean_state%multilayer%idx_pseudo_salt
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_S
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_S
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_S
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_S
            end if
         end if

         ! ---- Sponge ----
         bc%west%sponge_width = bc_cfg%sponge_width
         bc%east%sponge_width = bc_cfg%sponge_width
         bc%south%sponge_width = bc_cfg%sponge_width
         bc%north%sponge_width = bc_cfg%sponge_width
         bc%west%sponge_strength = bc_cfg%sponge_strength
         bc%east%sponge_strength = bc_cfg%sponge_strength
         bc%south%sponge_strength = bc_cfg%sponge_strength
         bc%north%sponge_strength = bc_cfg%sponge_strength
         bc%west%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%east%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%south%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%north%sponge_relax_tracers = bc_cfg%sponge_relax_tracers

         ! ---- Tidal constituents ----
         ntidal = min(bc_cfg%west_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%west%n_tidal_constituents = ntidal
         bc%west%tidal_amp(1:ntidal) = bc_cfg%west_tidal_amp(1:ntidal)
         bc%west%tidal_phase(1:ntidal) = bc_cfg%west_tidal_phase(1:ntidal)
         bc%west%tidal_omega(1:ntidal) = bc_cfg%west_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%east_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%east%n_tidal_constituents = ntidal
         bc%east%tidal_amp(1:ntidal) = bc_cfg%east_tidal_amp(1:ntidal)
         bc%east%tidal_phase(1:ntidal) = bc_cfg%east_tidal_phase(1:ntidal)
         bc%east%tidal_omega(1:ntidal) = bc_cfg%east_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%south_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%south%n_tidal_constituents = ntidal
         bc%south%tidal_amp(1:ntidal) = bc_cfg%south_tidal_amp(1:ntidal)
         bc%south%tidal_phase(1:ntidal) = bc_cfg%south_tidal_phase(1:ntidal)
         bc%south%tidal_omega(1:ntidal) = bc_cfg%south_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%north_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%north%n_tidal_constituents = ntidal
         bc%north%tidal_amp(1:ntidal) = bc_cfg%north_tidal_amp(1:ntidal)
         bc%north%tidal_phase(1:ntidal) = bc_cfg%north_tidal_phase(1:ntidal)
         bc%north%tidal_omega(1:ntidal) = bc_cfg%north_tidal_omega(1:ntidal)

         ! ---- OBC tidal nodal/astronomical correction (capability C3) ----
         ! When on, bake the 18.6-yr nodal factor f_c and the equilibrium+nodal
         ! phase (V_c + u_c) into each edge's tidal_fnodal/tidal_arg, using the
         ! SHARED &ocean_tides_nml reference epoch so the boundary tide stays
         ! phase-consistent with the interior body tide (C1).  Host setup —
         ! plain loops, no do concurrent.  Default off ⇒ f=1, arg=0 (the
         ! defaults) ⇒ legacy static-phase OBC sum bit-identical.
         bc%tidal_nodal = bc_cfg%obc_tidal_nodal
         if (bc%tidal_nodal) then
            block
               integer :: yr, mo, dy
               logical :: date_ok
               real(wp) :: dref, dnodal
               real(wp) :: v_all(TIDES_CATALOG_SIZE)
               real(wp) :: f_all(TIDES_CATALOG_SIZE), u_all(TIDES_CATALOG_SIZE)
               character(len=16) :: nodal_str

               ! Equilibrium argument V_c at the model-t=0 reference date.
               call parse_date_string(cfg%ocean%tides%ref_date, yr, mo, dy, date_ok)
               if (.not. date_ok) then
                  call fail("OBC tides: unparseable &ocean_tides_nml ref_date '"// &
                            trim(cfg%ocean%tides%ref_date)//"'", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
               dref = days_since_1900(yr, mo, dy)

               ! Nodal f/u at the nodal reference date ("" ⇒ ref_date).
               nodal_str = cfg%ocean%tides%nodal_ref_date
               if (len_trim(nodal_str) == 0) nodal_str = cfg%ocean%tides%ref_date
               call parse_date_string(nodal_str, yr, mo, dy, date_ok)
               if (.not. date_ok) then
                  call fail("OBC tides: unparseable &ocean_tides_nml nodal_ref_date '"// &
                            trim(nodal_str)//"'", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
               dnodal = days_since_1900(yr, mo, dy)

               call equilibrium_arguments(dref, v_all)
               call nodal_fu(dnodal, .true., f_all, u_all)

               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%west, f_all, u_all, v_all, "west", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%west, f_all, u_all, v_all, "west", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%east, f_all, u_all, v_all, "east", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%east, f_all, u_all, v_all, "east", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%south, f_all, u_all, v_all, "south", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%south, f_all, u_all, v_all, "south", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%north, f_all, u_all, v_all, "north", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%north, f_all, u_all, v_all, "north", compute_rank)
               end if
            end block
         end if

         ! ---- Open-edge tracer reservoirs (§1, v2) ----
         ! Cache length-scale knobs on bc so kernels can read them without
         ! going back to cfg.  Then allocate reservoir arrays when the feature
         ! is enabled and the edge is open-ish.  The RAW tag is compared on
         ! purpose: an OBC_SPONGE edge still gets a reservoir, which is inert
         ! there (the ghost fill is gated on is_open_ish, which excludes it,
         ! and the update sees the zero outer-face flux, so tres holds its
         ! seed).  Mapping through ocean_bc_outer_face_tag would drop the
         ! bc_tres_* fields from the restart registry.
         bc%res_lscale_out = bc_cfg%res_lscale_out
         bc%res_lscale_in = bc_cfg%res_lscale_in

         if ((bc%res_lscale_out > 0.0_wp .or. bc%res_lscale_in > 0.0_wp) .and. &
             bc%n_tracers > 0) then
            ! West
            if (bc%west%bc_type /= OBC_WALL .and. bc%west%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_west(bc%ny_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! East
            if (bc%east%bc_type /= OBC_WALL .and. bc%east%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_east(bc%ny_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! South
            if (bc%south%bc_type /= OBC_WALL .and. bc%south%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_south(bc%nx_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! North
            if (bc%north%bc_type /= OBC_WALL .and. bc%north%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_north(bc%nx_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if

            ! Seed reservoirs from the interior IC concentration.
            ! We seed from clamped_tracer (which was just populated above from
            ! west/east/south/north_inflow_S/T).  This is the best available
            ! "interior" estimate at setup time; the first-step update will
            ! refine it from the actual h_layer state.  Per-layer init is
            ! uniform in k (the IC is typically also uniform in k at setup time).
            if (allocated(bc%tres_west)) then
               do i = 1, bc%n_tracers
                  bc%tres_west(:, :, i) = bc%west%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_east)) then
               do i = 1, bc%n_tracers
                  bc%tres_east(:, :, i) = bc%east%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_south)) then
               do i = 1, bc%n_tracers
                  bc%tres_south(:, :, i) = bc%south%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_north)) then
               do i = 1, bc%n_tracers
                  bc%tres_north(:, :, i) = bc%north%clamped_tracer(i)
               end do
            end if

            if (compute_rank == 0) then
               call logger%info("OBC reservoirs: res_lscale_out=" &
                                //to_string(bc%res_lscale_out)//"m" &
                                //" res_lscale_in="//to_string(bc%res_lscale_in)//"m")
            end if
         end if

         ! ---- Per-layer Orlanski radiation (§2, v2) ----
         ! Parse the radiation_scheme string to an integer tag (0=anomaly, 1=orlanski).
         ! Allocate rx and u_prev arrays on radiating edges when the scheme is orlanski.
         ! Both arrays are initialised to 0; u_prev will be seeded from the IC velocity
         ! on the first ocean_obc_apply_baroclinic call (cold-start: dhdt = 0 ⇒ rx = 0).
         ! Not restart-registered — restarts cold (known gap).
         bc%orlanski_rx_max = bc_cfg%orlanski_rx_max
         bc%orlanski_gamma = bc_cfg%orlanski_gamma
         bc%nudge_tau_in = bc_cfg%nudge_tau_in
         bc%nudge_tau_out = bc_cfg%nudge_tau_out
         select case (trim(bc_cfg%radiation_scheme))
         case ("orlanski")
            bc%radiation_scheme = 1
         case default   ! "anomaly" or anything else: v1 path
            bc%radiation_scheme = 0
         end select

         if (bc%radiation_scheme == 1) then
            ! Raw tags: an OBC_SPONGE edge allocates these buffers but never
            ! uses them (the Orlanski kernels are gated on is_radiating).
            ! West
            if (bc%west%bc_type /= OBC_WALL .and. bc%west%bc_type /= OBC_PERIODIC .and. &
                bc%west%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_west(bc%ny_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_west(bc%ny_total, bc%nz_ml), source=0.0_wp)
            end if
            ! East
            if (bc%east%bc_type /= OBC_WALL .and. bc%east%bc_type /= OBC_PERIODIC .and. &
                bc%east%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_east(bc%ny_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_east(bc%ny_total, bc%nz_ml), source=0.0_wp)
            end if
            ! South
            if (bc%south%bc_type /= OBC_WALL .and. bc%south%bc_type /= OBC_PERIODIC .and. &
                bc%south%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_south(bc%nx_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_south(bc%nx_total, bc%nz_ml), source=0.0_wp)
            end if
            ! North
            if (bc%north%bc_type /= OBC_WALL .and. bc%north%bc_type /= OBC_PERIODIC .and. &
                bc%north%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_north(bc%nx_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_north(bc%nx_total, bc%nz_ml), source=0.0_wp)
            end if

            if (compute_rank == 0) then
               call logger%info("OBC Orlanski: rx_max=" &
                                //to_string(bc%orlanski_rx_max) &
                                //" gamma="//to_string(bc%orlanski_gamma))
            end if
         end if

         ! ---- Full Flather with exterior velocity (§4, v2) ----
         select case (trim(bc_cfg%flather_form))
         case ("full")
            bc%use_full_flather = .true.
         case default   ! "legacy" or anything else
            bc%use_full_flather = .false.
         end select
         bc%ext_u_west = bc_cfg%west_ext_u
         bc%ext_u_east = bc_cfg%east_ext_u
         bc%ext_v_south = bc_cfg%south_ext_v
         bc%ext_v_north = bc_cfg%north_ext_v

         ! ---- Periodic validation (also refreshes periodic_x/y flags) ----
         ! Only need to call if any axis could be periodic.
         if (bc%west%bc_type == OBC_PERIODIC .or. bc%east%bc_type == OBC_PERIODIC .or. &
             bc%south%bc_type == OBC_PERIODIC .or. bc%north%bc_type == OBC_PERIODIC) then
            if (present(ierr)) then
               call ocean_bc_validate_periodic(bc, ierr=local_ierr)
               if (local_ierr /= 0) then
                  ierr = local_ierr
                  return
               end if
            else
               call ocean_bc_validate_periodic(bc)
            end if
         else
            ! Derive periodic flags from final tags even for non-periodic runs.
            bc%periodic_x = .false.
            bc%periodic_y = .false.
         end if

         ! ---- Tripolar north-fold validation (also refreshes north_fold) ----
         ! Always call: it enforces "fold only on north" for every tag set,
         ! and "fold requires periodic w/e + nghost>=3" when the tag is set.
         if (present(ierr)) then
            call ocean_bc_validate_fold(bc, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call ocean_bc_validate_fold(bc)
         end if

         if (compute_rank == 0) then
            call logger%info("OBC:  west="//trim(bc_cfg%west)// &
                             " east="//trim(bc_cfg%east)// &
                             " south="//trim(bc_cfg%south)// &
                             " north="//trim(bc_cfg%north))
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_bc

   module subroutine configure_ocean_sponge(cfg, ocean_state, grid, compute_rank)
      !! Populate `ocean_state%sponge`'s per-cell `idamp_h`/`idamp_u`/
      !! `idamp_v` maps from `cfg%ocean%sponge` (`&ocean_sponge_nml`) + the
      !! already-configured `ocean_state%bc` edge tags. Run AFTER
      !! `configure_ocean_bc` (reads the edge tags + `has_*` flags) and
      !! BEFORE `ocean_state_enter_data`.
      !!
      !! `damp_source = "band"` (the only value implemented in v1): cosine
      !! ramp from every `OBC_SPONGE`-tagged edge, at the exact cell/u-face/
      !! v-face offsets the legacy `rdb_ocean_sponge::ocean_sponge_apply{,
      !! _tracers}` kernels use (§3.2 of the plan — reproduces today's band
      !! so `enable=.true., damp_source="band"` is *physically* the same
      !! sponge as today). Overlapping edges (corners) SUM their rates — the
      !! exact continuation of the legacy kernel's SEQUENTIAL per-edge decay
      !! (`exp(-a*dt)*exp(-b*dt) = exp(-(a+b)*dt)`).
      !!
      !! Per-edge `west_width`/`west_strength` etc. (sentinel `< 0` ⇒
      !! inherit the single global `&ocean_bc_nml sponge_width`/
      !! `sponge_strength` already fanned out onto `bc%<edge>%sponge_*` by
      !! `configure_ocean_bc`) close the audit gap that the per-edge tag
      !! fields were previously unreachable (no config path could ever make
      !! them differ).
      !!
      !! Does NOT snapshot the reference state — `ocean_sponge_snapshot_reference`
      !! does that, called MUCH earlier from the driver (right after the IC
      !! seed, before any restart read — see that subroutine's docstring for
      !! why the ordering is load-bearing).
      !!
      !! No-op when `.not. ocean_state%sponge%enable` (default) — the maps
      !! stay at their `init`-time zero fill (or unallocated, if `enable`
      !! was false at `init` time too).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      integer :: i0, i1, j0, j1, band, ramp
      real(wp) :: strength
      integer :: nonzero_h, nonzero_u, nonzero_v

      associate (sp => ocean_state%sponge, bc => ocean_state%bc, s_cfg => cfg%ocean%sponge)

         if (.not. sp%enable) return

         sp%relax_uv = s_cfg%relax_uv
         sp%relax_tracers = s_cfg%relax_tracers
         sp%relax_h = s_cfg%relax_h
         sp%damp_source = s_cfg%damp_source
         sp%target_source = s_cfg%target_source
         sp%lin_t_ref = s_cfg%lin_t_ref
         sp%lin_dt_dz = s_cfg%lin_dt_dz
         sp%lin_s_ref = s_cfg%lin_s_ref
         sp%lin_ds_dz = s_cfg%lin_ds_dz
         ramp = merge(SPONGE_RAMP_LINEAR, SPONGE_RAMP_COSINE, &
                      trim(s_cfg%ramp) == "linear")

         ! Latch the two registry indices the analytic `linear_z` refresh
         ! needs, so the per-step kernel never walks the tracer registry
         ! (the array-of-derived-types device-indirection rule).
         sp%idx_t = ocean_state%multilayer%idx_temperature
         sp%idx_s = ocean_state%multilayer%idx_salinity

         ! Latch the geopotential depth of the column TOP.  The draft is
         ! static by design, so this is a configure-time copy and the
         ! sponge slot stays self-contained at run time.  No cavity ⇒
         ! `z_draft` is the (1,1) placeholder and `z_top` keeps its
         ! init-time zero (the open-ocean free-surface datum).
         if (ocean_state%metrics%use_cavity .and. &
             size(ocean_state%metrics%z_draft, 1) == size(sp%z_top, 1) .and. &
             size(ocean_state%metrics%z_draft, 2) == size(sp%z_top, 2)) then
            sp%z_top = ocean_state%metrics%z_draft
         else
            sp%z_top = 0.0_wp
         end if

         i0 = grid%nghost + 1
         i1 = grid%nghost + grid%nx_phys
         j0 = grid%nghost + 1
         j1 = grid%nghost + grid%ny_phys

         if (trim(sp%damp_source) == "band") then
            ! ---- West edge ----
            if (bc%west%bc_type == OBC_SPONGE .and. bc%has_west) then
               band = merge(s_cfg%west_width, bc%west%sponge_width, s_cfg%west_width >= 0)
               strength = merge(s_cfg%west_strength, bc%west%sponge_strength, &
                                s_cfg%west_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_x(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + 1, band=band, &
                                         strength=strength, side=+1, j0=j0, j1=j1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- East edge ----
            if (bc%east%bc_type == OBC_SPONGE .and. bc%has_east) then
               band = merge(s_cfg%east_width, bc%east%sponge_width, s_cfg%east_width >= 0)
               strength = merge(s_cfg%east_strength, bc%east%sponge_strength, &
                                s_cfg%east_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_x(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + grid%nx_phys + 1, band=band, &
                                         strength=strength, side=-1, j0=j0, j1=j1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- South edge ----
            if (bc%south%bc_type == OBC_SPONGE .and. bc%has_south) then
               band = merge(s_cfg%south_width, bc%south%sponge_width, s_cfg%south_width >= 0)
               strength = merge(s_cfg%south_strength, bc%south%sponge_strength, &
                                s_cfg%south_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_y(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + 1, band=band, &
                                         strength=strength, side=+1, i0=i0, i1=i1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- North edge ----
            if (bc%north%bc_type == OBC_SPONGE .and. bc%has_north) then
               band = merge(s_cfg%north_width, bc%north%sponge_width, s_cfg%north_width >= 0)
               strength = merge(s_cfg%north_strength, bc%north%sponge_strength, &
                                s_cfg%north_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_y(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + grid%ny_phys + 1, band=band, &
                                         strength=strength, side=-1, i0=i0, i1=i1, &
                                         ramp=ramp)
               end if
            end if
         end if

         if (compute_rank == 0) then
            nonzero_h = count(sp%idamp_h > 0.0_wp)
            nonzero_u = count(sp%idamp_u > 0.0_wp)
            nonzero_v = count(sp%idamp_v > 0.0_wp)
            call logger%info("Sponge (map-driven): damp_source="//trim(sp%damp_source)// &
                             " target_source="//trim(sp%target_source)// &
                             " idamp_h max="//to_string(maxval(sp%idamp_h))//" 1/s"// &
                             " nonzero(h/u/v)="//to_string(nonzero_h)//"/"// &
                             to_string(nonzero_u)//"/"//to_string(nonzero_v))
         end if
      end associate
   end subroutine configure_ocean_sponge

   module subroutine configure_obc_edge_nodal(face, f_all, u_all, v_all, edge_name, compute_rank, ierr)
      !! Bake the C3 nodal/astronomical correction into one OBC edge, then
      !! fail loud (rank-0 error + `error stop`) if any of the edge's
      !! constituent frequencies matches no tide-catalog entry, else log the
      !! resolved constituent → (name, f_c, V_c+u_c) map on rank 0.  Host-side
      !! setup wrapper around the `pure` `obc_tide_nodal_fill` — the pure
      !! helper signals the failure via `fill_ierr`; the loud policy lives
      !! here.
      type(ocean_bc_face_tag_t), intent(inout) :: face
      real(wp), intent(in) :: f_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: u_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: v_all(TIDES_CATALOG_SIZE)
      character(len=*), intent(in) :: edge_name
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero when an OBC tidal constituent matches no tide-catalog
         !! entry, when present; absent behaves as today (`error stop`).

      integer :: fill_ierr, nc, ic

      call obc_tide_nodal_fill(face, f_all, u_all, v_all, fill_ierr)
      if (fill_ierr /= 0) then
         call fail("OBC tides: "//trim(edge_name)//" edge constituent " &
                   //to_string(fill_ierr)//" omega="//to_string(face%tidal_omega(fill_ierr)) &
                   //" rad/s matches no tide-catalog entry (rel tol 1e-4)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (compute_rank == 0) then
         do nc = 1, face%n_tidal_constituents
            ic = obc_match_constituent(face%tidal_omega(nc))
            call logger%info("OBC tide nodal ["//trim(edge_name)//"]: "//TIDE_NAME(ic) &
                             //" f="//to_string(face%tidal_fnodal(nc)) &
                             //" V+u="//to_string(face%tidal_arg(nc))//" rad")
         end do
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_obc_edge_nodal

end submodule rdb_ocean_setup_bc
