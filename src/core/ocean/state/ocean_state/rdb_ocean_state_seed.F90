submodule(rdb_ocean_state) rdb_ocean_state_seed
   !! R2b split: host-only seed/IC bodies, moved verbatim from
   !! rdb_ocean_state.F90 (analytical ICs, bathymetry/cavity seed
   !! dispatch, land-cell seeding; see
   !! src/core/ocean/state/rdb_ocean_state.F90 for the public interface).
   implicit none

contains

   module subroutine ocean_state_seed_from_cfg(state, grid, cfg, ierr, injected_b, injected_b_convention, &
                                               periodic_x, periodic_y)
      !! Populate the ocean prognostic state with an analytical IC
      !! derived from cfg scalars.  Bathymetry is set per `cfg%ocean%topo%topo_config`
      !! (`"flat"` → uniform `ocean_max_depth`; `"spoon"` → MOM6 spoon
      !! shape) — UNLESS `injected_b` is present, in which case it
      !! overrides `topo_config` entirely (P2.5 pre-create geometry
      !! injection).  Water column thickness `h = b` so the free surface
      !! starts at SSH = 0.  Layers split the local depth evenly
      !! (`h_layer(i,j,k) = b(i,j) / nz_ml`).  Tracers carry per-layer
      !! `h * Tr` at the configured initial T/S.  Velocities zeroed.
      !!
      !! Caller must have run `state%init(grid)` first so every slot's
      !! allocations are in place; this routine only writes into them
      !! on the host.  Run before `ocean_state_enter_data` so the GPU
      !! mapping captures the seeded values.
      !!
      !! Every position-aware formula fill below reads its global offsets
      !! and global extents from `grid` (`i_offset_global` /
      !! `j_offset_global`, `nx_global` / `ny_global`), so each rank seeds
      !! its own window on ONE global analytical IC.  There is deliberately
      !! no `decomp` argument and no optional offset quartet: an optional
      !! that defaults to the LOCAL extents is silently wrong under MPI, and
      !! forgetting to pass it shipped three separate bugs.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on an initial-condition configuration/data failure
         !! when present; absent behaves as today (`error stop`).
      real(wp), intent(in), optional :: injected_b(:, :)
         !! P2.5 pre-create geometry injection: interior-sized
         !! `(grid%nx_phys, grid%ny_phys)` bathymetry, RAW in the caller's
         !! own sign convention (see `injected_b_convention`) — REQUIRES
         !! `injected_b_convention`.  When present, this OVERRIDES
         !! `cfg%ocean%topo%topo_config` entirely: sign-normalised to
         !! positive-down depth + wet-fraction validated
         !! (`bathymetry_normalise_sign`), written into the interior of
         !! `state%barotropic%b`, and ghost-filled by constant
         !! extrapolation (`bathymetry_fill_ghosts_array`) — mirroring
         !! `topo_config='file'` rather than trusting the caller
         !! (CLAUDE.md formula-bathymetry ghost-fill gotcha).  Everything
         !! downstream in this routine (h=b, layer split, wet mask, T/S
         !! seed, IC overlays, the ALE z_ref table) then runs unchanged
         !! against the injected bathymetry.
      integer, intent(in), optional :: injected_b_convention
         !! REQUIRED alongside `injected_b` — no default (D6.2: sign is
         !! the single most dangerous argument in the geometry API). One
         !! of `BATHY_CONVENTION_DEPTH_POSITIVE_DOWN` /
         !! `_HEIGHT_POSITIVE_UP` (`rdb_ocean_bathymetry_inject`).
      logical, intent(in), optional :: periodic_x
         !! Grid topology the static geometry is wrapped with (see
         !! `seed_wrap_static_2d`).  Absent ⇒ derived from the
         !! `&ocean_bc_nml` west/east tags, exactly as `configure_ocean_bc`
         !! derives `bc%periodic_x`.  The engine passes it when a staged
         !! topology (`rdb_ocean_stage_topology`) overrides the tags.
      logical, intent(in), optional :: periodic_y
         !! As `periodic_x`, for south/north.

      integer :: nz_ml, idx_S, idx_T, i, j, k, nx, ny, local_ierr
      logical :: per_x, per_y, north_fold
      logical :: geographic
         !! Grid positions are in degrees (spherical/supergrid/tripolar):
         !! the input to `topo_length_to_grid_units` (#203).
      logical :: seeded_on_eta0_target
         !! The FOURTH or FIFTH `h_layer` seed branch (`zstar_full` under
         !! `zfixed_closed_faces`; MOM6 `zstar` under the knob or zinit) was
         !! taken ⇒ establish I1′ after the tracer IC.
      real(wp), allocatable :: water(:, :)
         !! Reference water-column thickness the IC seeds work on:
         !! `b − z_draft` under an ice shelf, a byte copy of `b`
         !! otherwise.  Host-only setup scratch, released on return.
      real(wp), allocatable :: eta_trim(:, :)
         !! Initial free-surface anomaly of the TRIMMED cavity IC
         !! (`&ocean_cavity_dyn_nml trim_ic_for_p_surf`); 0 when off.
      real(wp), allocatable :: h_col(:, :)
         !! Initial water column `water + eta_trim` the layer split seeds
         !! from; a byte copy of `water` when the trim is off.
      logical :: trim_ic
         !! `&ocean_cavity_dyn_nml trim_ic_for_p_surf` under an active cavity.

      nz_ml = state%multilayer%nz_ml
      idx_S = state%multilayer%idx_salinity
      idx_T = state%multilayer%idx_temperature
      nx = size(state%barotropic%b, 1)
      ny = size(state%barotropic%b, 2)

      ! Grid topology of the static geometry (the same rule
      ! `ocean_bc_state_init` applies to the tags, read here because the
      ! seed runs BEFORE `configure_ocean_bc`).
      per_x = ocean_bc_type_from_string(cfg%ocean%bc%west) == OBC_PERIODIC .and. &
              ocean_bc_type_from_string(cfg%ocean%bc%east) == OBC_PERIODIC
      per_y = ocean_bc_type_from_string(cfg%ocean%bc%south) == OBC_PERIODIC .and. &
              ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_PERIODIC
      if (present(periodic_x)) per_x = periodic_x
      if (present(periodic_y)) per_y = periodic_y
      north_fold = ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_TRIPOLAR_FOLD

      ! Bathymetry.  `slope_scale` (the spoon/seamount length scale) is a
      ! metres knob; `set_bathymetry_*` works in GRID coordinate units,
      ! which are metres on Cartesian grids but DEGREES on spherical/
      ! curvilinear ones.  Convert here so the length scale is consistent
      ! with the grid coordinate (without it, a metres scale against a
      ! degrees position collapses the seamount/spoon to a flat basin).
      if (present(injected_b)) then
         ! P2.5 pre-create geometry injection — bypasses topo_config
         ! entirely.  `injected_b_convention` has no default (D6.2).
         if (.not. present(injected_b_convention)) then
            call fail("ocean_state_seed_from_cfg: injected_b requires "// &
                      "injected_b_convention (no default — sign is the single most "// &
                      "dangerous argument; see D6.2 in the Python API design notes)", &
                      ierr, OCEAN_STATUS_ERR_BAD_SHAPE)
            return
         end if
         if (size(injected_b, 1) /= grid%nx_phys .or. size(injected_b, 2) /= grid%ny_phys) then
            call fail("ocean_state_seed_from_cfg: injected_b shape mismatch — got ("// &
                      to_string(size(injected_b, 1))//","//to_string(size(injected_b, 2))// &
                      "), expected the physical interior ("//to_string(grid%nx_phys)//","// &
                      to_string(grid%ny_phys)//")", ierr, OCEAN_STATUS_ERR_BAD_SHAPE)
            return
         end if
         block
            real(wp) :: b_local(grid%nx_phys, grid%ny_phys)
            integer :: ng
            b_local = injected_b
            call bathymetry_normalise_sign(b_local, injected_b_convention, ierr=local_ierr)
            if (local_ierr /= OCEAN_STATUS_OK) then
               if (present(ierr)) then
                  ierr = local_ierr
                  return
               end if
               error stop "ocean_state_seed_from_cfg: injected_b sign normalisation failed"
            end if
            ng = grid%nghost
            state%barotropic%b(ng + 1:ng + grid%nx_phys, ng + 1:ng + grid%ny_phys) = b_local
            call bathymetry_fill_ghosts_array(state%barotropic%b, grid)
         end block
      else
         ! Inline, NOT `grid_config_is_geographic(cfg%ocean%grid%grid_config)`:
         ! in this routine nvfortran 26.9 drops the hidden length on the first
         ! user call that passes this component, and the callee's `trim` then
         ! segfaults (#203). Intrinsic comparisons pass their lengths
         ! explicitly. Keep in step with `grid_config_is_geographic`.
         geographic = cfg%ocean%grid%grid_config == "spherical" .or. &
                      cfg%ocean%grid%grid_config == "supergrid" .or. &
                      cfg%ocean%grid%grid_config == "tripolar"
         select case (trim(cfg%ocean%topo%topo_config))
         case ("spoon")
            call set_bathymetry_spoon(state%barotropic%b, grid, &
                                      cfg%ocean%topo%max_depth, cfg%ocean%topo%edge_depth, &
                                      topo_length_to_grid_units(cfg%ocean%topo%slope_scale, &
                                                                geographic, &
                                                                cfg%ocean%grid%rad_earth))
         case ("seamount")
            call set_bathymetry_seamount(state%barotropic%b, grid, &
                                         cfg%ocean%topo%max_depth, cfg%ocean%topo%edge_depth, &
                                         topo_length_to_grid_units(cfg%ocean%topo%slope_scale, &
                                                                   geographic, &
                                                                   cfg%ocean%grid%rad_earth))
         case ("neverworld2")
            call set_bathymetry_neverworld2(state%barotropic%b, grid, &
                                            cfg%ocean%topo%max_depth, &
                                            cfg%ocean%topo%nl_continent_amp, &
                                            cfg%ocean%topo%nl_roughness_amp, &
                                            cfg%ocean%topo%nl_min_depth)
         case ("island")
            ! Flat basin at `max_depth` with a central LAND square (depth
            ! forced below LAND_DEPTH_THRESHOLD).  Produces the interior
            ! land block the static land-mask validation needs.  The square
            ! spans the central `island_frac` fraction of the physical
            ! domain in each direction (default via edge_depth/slope_scale
            ! re-use; see set_bathymetry_island).
            call set_bathymetry_island(state%barotropic%b, grid, &
                                       cfg%ocean%topo%max_depth, &
                                       cfg%ocean%topo%slope_scale)
         case ("double_drake")
            ! Ferreira et al. (2010) two-wall supercontinent with a
            ! reentrant southern channel — the land-mask showcase.
            call set_bathymetry_double_drake(state%barotropic%b, grid, &
                                             cfg%ocean%topo%max_depth, &
                                             cfg%ocean%topo%slope_scale)
         case ("isomip_plus")
            ! MISMIP+/ISOMIP+ analytic bedrock (Asay-Davis et al. 2016,
            ! Eqs. 1-4 + Table 1).  `max_depth` IS the deep clip
            ! (`-z_b,deep`, protocol 720 m) and `x_origin` places the
            ! model's west edge on the paper's absolute x axis (protocol
            ! 320 km).  The formula is written in METRES, so the setter
            ! is handed metres-per-grid-unit rather than a converted
            ! length: `topo_length_to_grid_units(1, ...)` is grid units
            ! per metre, and this is its reciprocal (exactly 1 on a
            ! Cartesian grid).
            call set_bathymetry_isomip_plus(state%barotropic%b, grid, &
                                            cfg%ocean%topo%max_depth, &
                                            cfg%ocean%topo%x_origin, &
                                            1.0_wp/topo_length_to_grid_units(1.0_wp, &
                                                                             geographic, &
                                                                             cfg%ocean%grid%rad_earth))
         case ("file")
            ! Real bathymetry from NetCDF.  File must be pre-projected onto
            ! the model's Cartesian grid (matching nx_phys × ny_phys); the
            ! loader fills ghost rows by constant extrapolation.  Sign
            ! convention: `b` is bottom depth positive-down (consistent with
            ! the spoon + flat branches).  GEBCO + ETOPO ship elevation
            ! positive-up — the preprocessing script should flip the sign
            ! before writing the NetCDF.
            if (len_trim(cfg%bathymetry_file) == 0) then
               call fail("ocean_state_seed_from_cfg: topo_config='file' requires "// &
                         "bathymetry_file", ierr, OCEAN_STATUS_ERR_IC_SEED)
               return
            end if
#ifndef RDB_NO_NETCDF
            if (present(ierr)) then
               call load_bathymetry_into_array(trim(cfg%bathymetry_file), &
                                               state%barotropic%b, grid, ierr=local_ierr)
               if (local_ierr /= OCEAN_STATUS_OK) then
                  ierr = local_ierr
                  return
               end if
            else
               call load_bathymetry_into_array(trim(cfg%bathymetry_file), &
                                               state%barotropic%b, grid)
            end if
#else
            ! File-loaded bathymetry needs the NetCDF reader in `rdb_bathymetry`,
            ! which isn't compiled in when RDB_ENABLE_NETCDF=OFF.  Fail fast
            ! with a descriptive error so the user knows to rebuild with NetCDF
            ! enabled rather than hitting a less-obvious runtime issue.
            call fail("ocean_state_seed_from_cfg: topo_config='file' requires "// &
                      "RDB_ENABLE_NETCDF=ON at build time (the NetCDF-backed "// &
                      "rdb_bathymetry reader is needed to load the file).", ierr, OCEAN_STATUS_ERR_IO)
            return
#endif
         case default  ! "flat"
            state%barotropic%b = cfg%ocean%topo%max_depth
         end select
      end if

      ! Make the bathymetry SEAM-CONSISTENT before anything reads it.  The
      ! file loader and the injected-array path fill the ghosts by constant
      ! extrapolation, and a formula setter evaluates its formula OUTSIDE
      ! the domain; on a periodic or folded edge neither is the value the
      ! seam needs.  Every field seeded below (the water column, h_layer,
      ! the wet mask, the tracers, the zstar_full z_ref table) and every
      ! setup-time consumer before the engine's first halo pass (the PGF's
      ! own bathymetry copy, bt_H_ref, the wet/dry and sponge setup) reads
      ! these ghosts — a stale seam ghost is a 1.9 m/s jet on the seam face
      ! of the 1-degree global grid within three hours.  The engine still
      ! re-wraps + halo-exchanges `b` later (the only fill a DECOMPOSED
      ! axis can get); on the local axes that is now a no-op.
      call seed_wrap_static_2d(state%barotropic%b, grid, per_x, per_y, north_fold)

      ! ---- Static ice-shelf cavity geometry (&ocean_cavity_dyn_nml, P5.1) ----
      ! ORDERING IS LOAD-BEARING and this is the only place it can go: the
      ! draft must exist before the water column `b − z_draft` seeds the
      ! layer split and the wet mask (grounding), both of which happen a
      ! few lines below, and long before any `configure_ocean_*` pass.
      ! The formula setters fill the FULL array including ghosts; the
      ! periodic/fold re-wrap + halo exchange that `barotropic%b` gets are
      ! applied to `z_draft` alongside it in `rdb_ocean_engine`.
      !
      ! `water` is the reference water-column thickness every seed below
      ! works on.  With the knob off it is a byte copy of `b`, so there is
      ! ONE code path and the default run is bit-identical.
      allocate (water(nx, ny))
      if (state%metrics%use_cavity) then
         call seed_cavity_draft(state, grid, cfg, ierr=local_ierr)
         if (local_ierr /= OCEAN_STATUS_OK) then
            if (present(ierr)) then
               ierr = local_ierr
               return
            end if
            error stop "ocean_state_seed_from_cfg: ice-shelf cavity geometry failed"
         end if
         ! Same seam rule as `b` just above, for the same reason: the
         ! draft and its cover fraction are bathymetry-class geometry, and
         ! the water column below reads their ghosts.
         call seed_wrap_static_2d(state%metrics%z_draft, grid, per_x, per_y, north_fold)
         call seed_wrap_static_2d(state%metrics%cover_frac, grid, per_x, per_y, north_fold)
         call cavity_water_column_impl(water, state%barotropic%b, &
                                       state%metrics%z_draft, nx, ny)
      else
         water = state%barotropic%b
      end if

      ! MOM6 TRIM_IC_FOR_P_SURF (`&ocean_cavity_dyn_nml
      ! trim_ic_for_p_surf`, default off).  The ice LOAD stays the
      ! Boussinesq-isostatic `rho_ref*g*z_draft`; each loaded column's
      ! initial TOP moves to the depth where the displaced water's own
      ! weight equals it, so the MOM6 barotropic split starts at rest (see
      ! `cavity_trim_eta_linear_impl`).  `water` itself — the datum and
      ! the grounding decision — is untouched: the trim is an initial
      ! `bt_eta`, not a geometry change.  The density is the linear EOS
      ! over the analytic zinit profile, both enforced by
      ! `validate_config`, and taken from `cfg` because they are exactly
      ! what `ocean_state_init` copied onto `state%eos`.
      allocate (eta_trim(nx, ny), source=0.0_wp)
      trim_ic = state%metrics%use_cavity .and. cfg%ocean%cavity_dyn%trim_ic_for_p_surf
      if (trim_ic) then
         block
            real(wp) :: rho_surf, drho_dz
            logical :: trim_ok
            rho_surf = cfg%ocean%ic%rho_0 + &
                       cfg%ocean%ic%beta_S*(cfg%ocean%zinit%lin_s_ref - cfg%ocean%ic%S_ref) - &
                       cfg%ocean%ic%alpha_T*(cfg%ocean%zinit%lin_t_ref - cfg%ocean%ic%T_ref)
            drho_dz = cfg%ocean%ic%beta_S*cfg%ocean%zinit%lin_ds_dz - &
                      cfg%ocean%ic%alpha_T*cfg%ocean%zinit%lin_dt_dz
            call cavity_trim_eta_linear_impl(eta_trim, trim_ok, state%metrics%z_draft, &
                                             water, cfg%ocean%cavity_dyn%h_min_cavity, &
                                             cfg%ocean%ic%rho_0, rho_surf, drho_dz, nx, ny)
            if (.not. trim_ok) then
               call fail("ocean_state_seed_from_cfg: &ocean_cavity_dyn_nml "// &
                         "trim_ic_for_p_surf found no admissible trim depth (the "// &
                         "initial density must be positive and stably stratified, "// &
                         "and the trimmed column must stay non-empty)", &
                         ierr, OCEAN_STATUS_ERR_IC_SEED)
               return
            end if
            call logger%info("ocean_cavity: trimmed the initial column under the ice "// &
                             "to the load (MOM6 TRIM_IC_FOR_P_SURF): eta in ["// &
                             to_string(minval(eta_trim))//", "// &
                             to_string(maxval(eta_trim))//"] m")
         end block
      end if
      allocate (h_col(nx, ny))
      if (trim_ic) then
         h_col = water + eta_trim
      else
         h_col = water
      end if

      ! Water column thickness h = b − z_draft → the free-surface anomaly
      ! `bt_eta = Σ h_layer − bt_H_ref` starts at zero (SSH = 0 without a
      ! cavity; the loaded equilibrium under one) — or at `eta_trim` under
      ! a trimmed cavity IC.
      state%barotropic%h = h_col
      state%barotropic%u_face_x = 0.0_wp
      state%barotropic%v_face_y = 0.0_wp
      state%barotropic%hu_face_x = 0.0_wp
      state%barotropic%hv_face_y = 0.0_wp

      ! Per-column layer thickness = depth / nz_ml.  Pass flat
      ! allocatables into a `_impl` helper so the `do concurrent` body
      ! reads/writes plain arrays — registry / array-of-DT indirection
      ! crashes NVHPC's device codegen with CUDA_ERROR_ILLEGAL_ADDRESS.
      ! With wet/dry ON, floor the seed to 2·H_VANISHED so an emerged
      ! intertidal rest-bed (b < 0 ⇒ now a LIVE wet_mask=1 column) never
      ! seeds a negative layer / negative barotropic depth D = Σ h_layer,
      ! and each seeded layer clears the strict `h_old > H_VANISHED`
      ! remap-drain gate (so the seeded S/T survives the first regrid).
      ! Tracers are seeded from THIS floored h_layer just below, so hTr
      ! stays consistent (const·2·H_VANISHED, finite T/S, no negative mass).
      ! `thickness_config = "uniform_z"` swaps the local-depth even split for
      ! MOM6's uniform-z-interface seed (flat resting isopycnals under a
      ! horizontally-uniform density stack).  Validated as mutually exclusive
      ! with wet/dry, so the two branches never both need the emerged-column
      ! floor.  Default "sigma" ⇒ byte-identical to the pre-knob path.
      !
      ! THIRD BRANCH — `VCOORD_Z_FIXED` under a cavity.  The running
      ! coordinate there is quasi-geopotential with inert fillers inside
      ! the ice, so a sigma-style seed is NOT on the coordinate: the very
      ! first ALE remap would relamp the whole column in one step, and a
      ! T/S profile that `&ocean_zinit_nml source="linear"` made exactly
      ! linear in geopotential z would come back through the PPM boundary
      ! closure NOT exactly linear — column by column, because the draft
      ! (and so the cut) differs column to column.  That difference IS a
      ! horizontal density gradient, i.e. exactly the spurious rest
      ! current this coordinate exists to remove.  So seed `h_layer`
      ! directly FROM the target (`η = 0`, which is the cavity datum's
      ! own resting state) and let the zinit overlay evaluate T/S at
      ! those layer centres: exact by construction, first remap an
      ! identity, no step-1 regrid shock.
      !
      ! The same holds WITHOUT a cavity whenever the T/S come from the
      ! GEOPOTENTIAL `&ocean_zinit_nml` overlay: a sigma-style seed puts
      ! layer `k` of a column of depth `H` at `(nz-k+1/2)*H/nz`, the zinit
      ! profile is evaluated THERE, and only the step-1 regrid moves it onto
      ! the `z_fixed` layers.  Anything snapshotted from the seed in between
      ! is then indexed on the wrong layers -- the `&ocean_sponge_nml
      ! target_source = "ic"` reference is (`ocean_sponge_snapshot_reference`
      ! runs right after this seed), so the sponge relaxed layer `k` of the
      ! 50-level tanh grid toward the zinit value hundreds of metres deeper,
      ! and toward a DIFFERENT depth in every column of a different `H`: a
      ! grid-scale, bathymetry-following horizontal density forcing across
      ! the whole band (measured: the band's surface layers pulled 5-6 degC
      ! cold, domain-mean T -0.55 degC in 5 days on the coastal-noise box).
      ! Seeding on the target makes the overlay exact and the snapshot
      ! consistent.  `z_fixed` WITHOUT zinit keeps the sigma-style seed: its
      ! `&tracer_nml` IC is defined PER LAYER INDEX, so moving the layers
      ! would change what that IC means.
      !
      ! FOURTH BRANCH — `VCOORD_ZSTAR_FULL` under `&vcoord_nml
      ! zfixed_closed_faces`.  The closed-face mask is built from the
      ! `ZSTAR_FULL` target at `η = 0`; a sigma-style seed is not on that
      ! coordinate, so step 1 would run the FULL sigma pressure gradient on
      ! the `b/nz` stack (on the 1-degree Southern Ocean, 0.67 m/s from
      ! rest in one step, at faces the mask does not even see) before the
      ! first regrid moved the layers, and a `target_source = "ic"` sponge
      ! would snapshot T/S on the wrong layers — the `z_fixed` seed bug
      ! above.  So seed `h_layer` from the SAME target the mask is built
      ! from (`ocean_vcoord_eta0_target`, after laying the `z_ref` table
      ! it walks), with or without zinit: a per-index `&tracer_nml` IC
      ! then means "per coordinate layer", which is the only reading under
      ! which the mask's live/filler pattern is the IC's.  Knob-gated, so
      ! every existing `zstar_full` namelist keeps its sigma-style seed
      ! byte for byte.
      !
      ! FIFTH BRANCH — `VCOORD_ZSTAR` (MOM6 z*) under `zfixed_closed_faces`
      ! OR `&ocean_zinit_nml`.  The same reasons as the fourth (the mask's
      ! pattern is the `η = 0` target's) and the third (a z-level IC is
      ! exact only at the running coordinate's own layer centres; a
      ! `target_source = "ic"` sponge snapshots what is seeded), on the
      ! `η = 0` z* target, which is the `z_fixed` one with `z_top = 0`.  A
      ! per-layer-index `&tracer_nml` IC without the knob keeps the
      ! sigma-style seed and is moved onto z* by the first regrid, as on
      ! `z_fixed`.  z* is refused under a cavity, so `eta_trim ≡ 0` here.
      seeded_on_eta0_target = .false.
      if (trim(cfg%thickness_config) == "uniform_z") then
         call seed_h_layer_uniform_z_impl(state%multilayer%h_layer, &
                                          water, nz_ml, &
                                          cfg%ocean%topo%max_depth, &
                                          cfg%ocean%isopycnal%angstrom_h)
      else if ((state%metrics%use_cavity .or. cfg%ocean%zinit%enable) .and. &
               parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_Z_FIXED .and. &
               cfg%ocean%topo%max_depth > 0.0_wp) then
         block
            real(wp) :: h_min_seed
            real(wp), allocatable :: z_top_seed(:, :)
            ! Column-top depth: the ice draft under a cavity, else `z = 0`
            ! (`metrics%z_draft` is only a `(1,1)` placeholder then).
            if (state%metrics%use_cavity) then
               z_top_seed = state%metrics%z_draft
            else
               allocate (z_top_seed(nx, ny), source=0.0_wp)
            end if
            ! `zstar_h_min` comes off the SLOT, not off `cfg`: there is
            ! one source of truth for the filler thickness and it is the
            ! one the running target builder will use.  `engine_setup`
            ! copies both `zstar_*` knobs onto the slot immediately BEFORE
            ! this seed (it has to — the seed's own tail calls
            ! `vcoord%build_zref_full`, which reads them), so the value is
            ! already the namelist's.  Fall back to `cfg` only for a
            ! caller that seeds a state whose vcoord slot was never
            ! initialised.
            h_min_seed = cfg%zstar_h_min
            if (state%vcoord%is_init) h_min_seed = state%vcoord%zstar_h_min
            ! A stretched nominal profile (`&vcoord_nml z_fixed_profile`)
            ! is installed on the slot by `engine_setup` before this seed,
            ! for the same one-source-of-truth reason as `zstar_h_min`.
            if (state%vcoord%is_init .and. state%vcoord%z_fixed_use_profile) then
               call ocean_vcoord_z_fixed_target(state%multilayer%h_layer, water, eta_trim, &
                                                z_top_seed, nx, ny, nz_ml, &
                                                cfg%ocean%topo%max_depth/real(nz_ml, wp), &
                                                .true., state%vcoord%z_fixed_zi, &
                                                state%vcoord%z_fixed_dz, h_min_seed)
            else
               call ocean_vcoord_z_fixed_target_uniform(state%multilayer%h_layer, water, &
                                                        eta_trim, z_top_seed, &
                                                        nx, ny, nz_ml, &
                                                        cfg%ocean%topo%max_depth/real(nz_ml, wp), &
                                                        h_min_seed)
            end if
         end block
      else if (cfg%zfixed_closed_faces .and. state%vcoord%is_init .and. &
               state%vcoord%coord_type == VCOORD_ZSTAR_FULL .and. &
               parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_ZSTAR_FULL) then
         call state%vcoord%build_zref_full(state%barotropic%b)
         call ocean_vcoord_eta0_target(state%vcoord, state%multilayer%h_layer, &
                                       water, nx, ny, nz_ml)
         seeded_on_eta0_target = .true.
      else if ((cfg%zfixed_closed_faces .or. cfg%ocean%zinit%enable) .and. &
               state%vcoord%is_init .and. &
               state%vcoord%coord_type == VCOORD_ZSTAR .and. &
               parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_ZSTAR .and. &
               state%vcoord%z_fixed_h_ref > 0.0_wp) then
         call ocean_vcoord_eta0_target(state%vcoord, state%multilayer%h_layer, &
                                       water, nx, ny, nz_ml)
         seeded_on_eta0_target = .true.
      else
         call seed_h_layer_uniform_impl(state%multilayer%h_layer, &
                                        h_col, nz_ml, &
                                        apply_wetdry_floor=cfg%ocean%wetdry%enable)
      end if
      ! Keep the barotropic water-column prognostic non-negative on the
      ! same emerged band (it is inert in the multilayer dyn-core — SSH is
      ! Σ h_layer − b — but is restart-registered as `bt_h`; floor it so a
      ! checkpoint never carries a negative rest depth).  Match the layer
      ! floor exactly: on an emerged column Σ h_layer = nz·2·H_VANISHED, so
      ! floor bt_h to the SAME value (not 0) — otherwise a restart that
      ! re-derives D from bt_h would disagree with Σ h_layer by nz·2·H_VANISHED
      ! on every emerged column.  Knob-off ⇒ h = b (byte-identical).
      if (cfg%ocean%wetdry%enable) then
         state%barotropic%h = max(water, &
                                  real(nz_ml, wp)*2.0_wp*H_VANISHED)
      end if
      state%multilayer%u_face_x_layer = 0.0_wp
      state%multilayer%v_face_y_layer = 0.0_wp
      state%multilayer%hu_face_x_layer = 0.0_wp
      state%multilayer%hv_face_y_layer = 0.0_wp

      ! Wet/dry mask: 1.0 where `b >= cutoff`, 0.0 elsewhere (land).
      ! Default cutoff = LAND_DEPTH_THRESHOLD (byte-identical for all
      ! non-wetdry configs).  When wetdry is enabled, pass `land_cutoff =
      ! -land_margin` so intertidal columns (bed above rest MSL but below the
      ! flood headroom) stay wet_mask=1 and the dynamic wd_wet_dyn gate
      ! handles their wetting/drying instead of the static land mask.
      if (cfg%ocean%wetdry%enable) then
         call seed_wet_mask_impl(state%multilayer%wet_mask, water, &
                                 land_cutoff=-cfg%ocean%wetdry%land_margin)
      else if (state%metrics%use_cavity) then
         ! GROUNDING.  A column with less than `h_min_cavity` of water
         ! under the ice is LAND — routed through the SAME wet-mask seed
         ! the bathymetry uses, so the static metric-zeroing land mask
         ! (`configure_ocean_land_mask`) and the finite land-state hold
         ! (`ocean_state_seed_land_cells`) follow for free.  Never a thin
         ! film of water under grounded ice.  Note the cutoff is applied
         ! to `b − z_draft`, so an ordinary land column (`b` below
         ! LAND_DEPTH_THRESHOLD, draft already zeroed there) is land for
         ! the same reason it always was.
         call seed_wet_mask_impl(state%multilayer%wet_mask, water, &
                                 land_cutoff=cfg%ocean%cavity_dyn%h_min_cavity)
      else
         call seed_wet_mask_impl(state%multilayer%wet_mask, water)
      end if

      ! Tracers carry per-layer h*Tr.  Multiply the (now spatially-
      ! varying) h_layer by the configured uniform scalar.
      if (idx_S > 0) then
         block
            real(wp) :: s_layer(nz_ml)
            real(wp) :: dS_dlayer
            logical :: stratify_s
            ! Linear S(z) when both surface + bottom are specified —
            ! EXACT mirror of the temperature branch below, including
            ! the gate (`both /= 0`), the vertical convention (k=1 is
            ! the bed = S_init_bottom, k=nz_ml the surface =
            ! S_init_surface), the single-layer fall-through, and the
            ! seed helper (so ghosts and land columns are filled the
            ! same way: hTr = S(k)·h_layer everywhere, land included,
            ! since h_layer is already 0/floored there).  The profile is
            ! linear in LAYER INDEX, which under the sigma-style
            ! `h_layer = b/nz_ml` seed is linear in layer-centre depth
            ! on every column — so a sloping bed gets the same endpoint
            ! values with a depth-proportional gradient.
            ! Note the stable polarity is the INVERSE of temperature:
            ! dense/salty water belongs at the bed, so a stable haline
            ! column has `S_init_bottom > S_init_surface`.
            ! Both-zero (the default) ⇒ uniform `initial_salinity`,
            ! byte-identical to the pre-knob path.
            stratify_s = (cfg%S_init_surface /= 0.0_wp) .and. &
                         (cfg%S_init_bottom /= 0.0_wp) .and. &
                         (nz_ml > 1)
            if (stratify_s) then
               dS_dlayer = (cfg%S_init_surface - cfg%S_init_bottom)/ &
                           real(nz_ml - 1, wp)
               do k = 1, nz_ml
                  s_layer(k) = cfg%S_init_bottom + dS_dlayer*real(k - 1, wp)
               end do
               call seed_tracer_stratified_impl( &
                  state%multilayer%tracers(idx_S)%hTr, &
                  state%multilayer%h_layer, s_layer, nz_ml)
            else
               call seed_tracer_uniform_impl( &
                  state%multilayer%tracers(idx_S)%hTr, &
                  state%multilayer%h_layer, &
                  cfg%initial_salinity, nz_ml)
            end if
         end block
      end if
      if (idx_T > 0) then
         block
            real(wp) :: t_layer(nz_ml)
            real(wp) :: dT_dlayer
            logical :: stratify
            ! Linear T(z) when both surface + bottom are specified.
            ! Convention: k=1 is the bed (T = T_init_bottom), k=nz_ml is
            ! the surface (T = T_init_surface).  Single-layer case
            ! falls through to T_init_surface (no profile to build).
            stratify = (cfg%T_init_surface /= 0.0_wp) .and. &
                       (cfg%T_init_bottom /= 0.0_wp) .and. &
                       (nz_ml > 1)
            if (stratify) then
               dT_dlayer = (cfg%T_init_surface - cfg%T_init_bottom)/ &
                           real(nz_ml - 1, wp)
               do k = 1, nz_ml
                  t_layer(k) = cfg%T_init_bottom + dT_dlayer*real(k - 1, wp)
               end do
            else
               t_layer = cfg%initial_temperature
            end if
            call seed_tracer_stratified_impl( &
               state%multilayer%tracers(idx_T)%hTr, &
               state%multilayer%h_layer, t_layer, nz_ml)
         end block
      end if

      ! Optional IC overlay applied after bathymetry + uniform h_layer
      ! + analytical T/S are in place.
      ! ierr threaded down ONLY when THIS routine's own ierr is present:
      ! otherwise each seed_*_ic helper must keep reaching its own
      ! `error stop` (specific text) rather than the generic wrapper
      ! message below (P0.1 review F2).
      if (present(ierr)) then
         select case (trim(cfg%ocean%ic%ic_config))
         case ("eady")
            call seed_eady_ic(state, grid, cfg, ierr=local_ierr)
         case ("geostrophic_adjustment")
            call seed_geostrophic_adjustment_ic(state, grid, cfg, ierr=local_ierr)
         case ("baroclinic_jet")
            call seed_baroclinic_jet_ic(state, grid, cfg, ierr=local_ierr)
         case default
            ! "" — keep the default analytical IC unchanged.
            local_ierr = 0
         end select
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         select case (trim(cfg%ocean%ic%ic_config))
         case ("eady")
            call seed_eady_ic(state, grid, cfg)
         case ("geostrophic_adjustment")
            call seed_geostrophic_adjustment_ic(state, grid, cfg)
         case ("baroclinic_jet")
            call seed_baroclinic_jet_ic(state, grid, cfg)
         case default
            ! "" — keep the default analytical IC unchanged.
         end select
      end if

      ! Z-level T/S IC overlay (capability A2).  Runs after bathymetry +
      ! uniform h_layer + wet_mask + analytical T/S are in place (it reads
      ! the seeded column depths + wet_mask) and intentionally OVERWRITES
      ! any analytical T/S.  Default off (`enable = .false.`) preserves
      ! bit-identity.  NetCDF-only: the reader lives in rdb_ocean_z_init,
      ! which only compiles with RDB_ENABLE_NETCDF=ON.
      !
      ! Under a cavity the overlay is handed `metrics%z_draft` so every
      ! layer centre's depth is measured from `z = 0` rather than from the
      ! ice base; see `seed_zinit_overlay`.
      if (cfg%ocean%zinit%enable) then
#ifndef RDB_NO_NETCDF
         ! A trimmed cavity IC moves the column top to `z_draft - eta_trim`.
         if (trim_ic) then
            if (present(ierr)) then
               call seed_zinit_overlay(state, grid, cfg, ierr=local_ierr, &
                                       z_top=state%metrics%z_draft - eta_trim)
               if (local_ierr /= 0) then
                  ierr = local_ierr
                  return
               end if
            else
               call seed_zinit_overlay(state, grid, cfg, &
                                       z_top=state%metrics%z_draft - eta_trim)
            end if
         else if (present(ierr)) then
            call seed_zinit_overlay(state, grid, cfg, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call seed_zinit_overlay(state, grid, cfg)
         end if
#else
         call fail("ocean_state_seed_from_cfg: ocean_zinit requires "// &
                   "RDB_ENABLE_NETCDF=ON at build time (the NetCDF-backed "// &
                   "rdb_ocean_z_init reader is needed to load the z-level T/S file).", ierr, OCEAN_STATUS_ERR_IO)
         return
#endif
      end if

      ! The on-target `zstar_full` / `zstar` seed (FOURTH / FIFTH BRANCH
      ! above) laid inert `zstar_h_min` fillers, and every IC writer
      ! since (`&tracer_nml` per layer index, the zinit overlay at the
      ! filler's own depth) gave them a concentration that is not their
      ! donor's.  Establish invariant I1′ NOW, with the one definition (host
      ! twin — before `enter_data`), so the sponge snapshot and the budget
      ! latch see the state every later step holds.  Left to the first
      ! in-step enforcement, the pooling moves that foreign content into the
      ! live partial cell above at step 1: a real horizontal density
      ! difference at a face the mask keeps open (measured 4.6e-6 m/s after
      ! two hours on `test_ocean_zstar_full_closed_faces`' resting
      ! staircase).  Column-conservative.  Before the pseudo-salt seed, so
      ! that copies the settled S.
      if (seeded_on_eta0_target) then
         call state%multilayer%enforce_vanished_content_host(nx, ny)
      end if

      ! Pseudo-salt seed — MUST run after every write to salinity's initial
      ! condition above (analytical IC, the eady/geostrophic/baroclinic-jet
      ! overlays, and critically the z-file overlay, which OVERWRITES any
      ! analytical S) and before ocean_state_enter_data (driver-ordered).
      ! Seeding before the z-file overlay would leave D(t=0) /= 0 in exactly
      ! the configuration (z-file IC) where the diagnostic matters most.
      ! Self-gates on idx_pseudo_salt <= 0 (knob off).
      call ocean_pseudo_salt_seed(state%multilayer)

      ! Direct per-layer density init — MOM6 `COORD_CONFIG="gprime"` analogue.
      ! When `ocean_layer_rho_init` is set (any non-sentinel value), write
      ! `ms%rho_layer(:,:,k) = ocean_layer_rho_init(k)` directly, bypassing
      ! the EOS path entirely.  Intended use case: `ocean_enable_thermodynamics
      ! = .false.` + this knob ⇒ static reduced-gravity stratification, mirroring
      ! MOM6's adiabatic gprime IC.  Without this, `rho_layer` stays at its
      ! alloc-time zero when thermo is off (the dyn step's EOS call is gated
      ! on `therm_active`), and FV-LITE's pressure-stack collapses to zero
      ! ⇒ no PGF, no Sverdrup balance, no WBC dynamics — only the wind +
      ! Coriolis side of the momentum equation produces anything.  Setting
      ! the per-layer ρ here is what turns the PGF back on under the
      ! adiabatic-stratified semantics MOM6 uses for `double_gyre`.
      !
      ! Convention: `k=1` bed → `k=nz_ml` surface (Roundabout bottom-up).  Pass
      ! the heavier value first.  Count of non-sentinel entries must match
      ! `nz_ml`; otherwise we abort with a descriptive error.
      !
      ! `rho_lightest >= 0` selects the linear density-range generator
      ! (MOM6 `COORD_CONFIG="linear"`): build `nz_ml` linearly-spaced
      ! densities and reuse the same write path, so the column scales by
      ! just bumping `nz_layers`.  Mutually exclusive with the explicit
      ! `layer_rho_init` list — set only one.
      if (cfg%ocean%ic%rho_lightest >= 0.0_wp) then
         if (count(cfg%ocean%ic%layer_rho_init >= 0.0_wp) > 0) then
            call fail("ocean_ic: rho_lightest (linear density-range) and "// &
                      "layer_rho_init (explicit list) are mutually exclusive — "// &
                      "set only one.", ierr, OCEAN_STATUS_ERR_IC_SEED)
            return
         end if
         if (present(ierr)) then
            call apply_layer_rho_init(state%multilayer, &
                                      ocean_linear_layer_density(cfg%ocean%ic%rho_lightest, &
                                                                 cfg%ocean%ic%rho_range, nz_ml), &
                                      nz_ml, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call apply_layer_rho_init(state%multilayer, &
                                      ocean_linear_layer_density(cfg%ocean%ic%rho_lightest, &
                                                                 cfg%ocean%ic%rho_range, nz_ml), &
                                      nz_ml)
         end if
      else
         if (present(ierr)) then
            call apply_layer_rho_init(state%multilayer, cfg%ocean%ic%layer_rho_init, nz_ml, &
                                      ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call apply_layer_rho_init(state%multilayer, cfg%ocean%ic%layer_rho_init, nz_ml)
         end if
      end if

      ! Populate the per-column z_ref table for VCOORD_ZSTAR_FULL.  Without
      ! this, `compute_target_h(ZSTAR_FULL)` walks an all-zero table and
      ! the first ALE remap collapses every layer to `zstar_h_min` — bug
      ! caught by the double_gyre adiabatic NK=2 setup (2026-05-26).
      ! Other vcoord types ignore z_ref so the call is safe regardless.
      call state%vcoord%build_zref_full(state%barotropic%b)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_state_seed_from_cfg

   pure subroutine seed_wrap_static_2d(fld, grid, per_x, per_y, north_fold)
      !! Fill the seam ghosts of a static, cell-centred 2-D geometry field
      !! (bathymetry, ice draft, cover fraction) from their periodic /
      !! north-fold images, on the host, at seed time.
      !!
      !! Only the axes this rank holds WHOLE are wrapped locally
      !! (`nx_phys == nx_global` / `ny_phys == ny_global`); on a decomposed
      !! axis the seam ghosts belong to a neighbour rank and are filled by
      !! the engine's init-time halo exchange, which runs later.  The fold
      !! is applied only on the rank that owns the north edge, and only when
      !! that rank holds the whole fold row: on an east-west split (px > 1)
      !! it needs the distributed exchange (`rdb_ocean_fold_exchange`), which
      !! is initialised later, so `engine_setup` folds these fields then
      !! (with the init folds).
      !! Every field seeded here FROM the unfolded north ghosts is repaired
      !! before anything reads it: h_layer and the tracers by the engine's
      !! cold-start fold, the stored wet mask by `configure_ocean_land_mask`
      !! (an elementwise map of the water column, so folding it reproduces
      !! the serial ghosts), the z_ref table, the PGF's bathymetry copy and
      !! bt_H_ref by their engine-side rebuild / re-fold after the fold; the
      !! barotropic CFL scan reads physical cells only.  Composition
      !! order matches the engine's init wrap: periodic first, so the fold
      !! reads the wrapped corner columns.
      real(wp), intent(inout) :: fld(:, :)
         !! Cell-centred field, shape `(nx_total, ny_total)`.
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: per_x
         !! West/east edges are periodic.
      logical, intent(in) :: per_y
         !! South/north edges are periodic.
      logical, intent(in) :: north_fold
         !! North edge is the tripolar fold.

      logical :: wrap_x, wrap_y, fold

      if (size(fld, 1) /= grid%nx_total .or. size(fld, 2) /= grid%ny_total) return
      wrap_x = per_x .and. grid%nx_phys == grid%nx_global
      wrap_y = per_y .and. grid%ny_phys == grid%ny_global
      fold = north_fold .and. grid%nx_phys == grid%nx_global .and. &
             grid%j_offset_global + grid%ny_phys == grid%ny_global
      if (wrap_x .or. wrap_y) then
         call ocean_periodic_wrap_centre_2d(fld, grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, grid%nghost, &
                                            wrap_x, wrap_y)
      end if
      if (fold) then
         call fold_north_centre(fld, grid%nx_total, grid%ny_total, &
                                grid%nx_phys, grid%ny_phys, grid%nghost)
      end if
   end subroutine seed_wrap_static_2d

#ifndef RDB_NO_NETCDF
   subroutine seed_zinit_overlay(state, grid, cfg, ierr, z_top)
      !! Dispatch the `&ocean_zinit_nml` T/S overlay across its two axes:
      !! the profile SOURCE (`"file"` — the pre-regridded NetCDF reader;
      !! `"linear"` — the analytic affine `lin_*` profile) and whether a
      !! cavity draft is present.
      !!
      !! The draft is the whole reason this is a separate routine.  Both
      !! seeders measure each layer centre's GEOPOTENTIAL depth from
      !! `z = 0`; under an ice shelf the column top is `z_draft` metres
      !! down, so `metrics%z_draft` has to reach them or a `T(z)` profile
      !! lands systematically too shallow (and, under a SLOPING lid,
      !! tilts the isopycnals with the ice base — not a state of rest).
      !! With the cavity off `metrics%z_draft` is the `(1, 1)`
      !! placeholder, so the argument is simply not passed and the
      !! arithmetic is bit-identical to the pre-cavity path.
      !!
      !! NOT `pure`: the seeders it dispatches to read files and log.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Threaded straight through to the seeder; an ABSENT `ierr`
         !! stays absent there, so each keeps its own `error stop` text.
      real(wp), intent(in), optional :: z_top(:, :)
         !! Column-top depth (m, positive down), FULL ghosted shape, in
         !! place of `metrics%z_draft`: the TRIMMED cavity IC
         !! (`&ocean_cavity_dyn_nml trim_ic_for_p_surf`) puts the top at
         !! `z_draft - eta_trim`.  Absent => `metrics%z_draft`.

      logical :: cav, conv

      cav = state%metrics%use_cavity
      ! File holding theta / SP into a CT / SA model: convert on read.
      ! (validate_config refuses the opposite pairing.)
      conv = trim(adjustl(cfg%ocean%eos%ts_convention)) == "cons_abs" .and. &
             trim(adjustl(cfg%ocean%zinit%input_convention)) == "pot_prac"
      if (present(z_top)) then
         if (trim(adjustl(cfg%ocean%zinit%source)) == "linear") then
            call seed_ts_linear_z(state%multilayer, cfg%ocean%zinit, ierr, z_top)
         else
            call seed_ts_from_zfile(state%multilayer, grid, cfg%ocean%zinit, ierr, z_top, &
                                    to_cons_abs=conv)
         end if
         return
      end if
      if (trim(adjustl(cfg%ocean%zinit%source)) == "linear") then
         if (cav) then
            call seed_ts_linear_z(state%multilayer, cfg%ocean%zinit, ierr, &
                                  state%metrics%z_draft)
         else
            call seed_ts_linear_z(state%multilayer, cfg%ocean%zinit, ierr)
         end if
      else
         if (cav) then
            call seed_ts_from_zfile(state%multilayer, grid, cfg%ocean%zinit, ierr, &
                                    state%metrics%z_draft, to_cons_abs=conv)
         else
            call seed_ts_from_zfile(state%multilayer, grid, cfg%ocean%zinit, ierr, &
                                    to_cons_abs=conv)
         end if
      end if
   end subroutine seed_zinit_overlay
#endif

   module subroutine apply_layer_rho_init(ms, layer_rho_init, nz_ml, ierr)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Helper: write per-layer density from `layer_rho_init(:)` when the
      !! caller has set at least one non-sentinel value.  Sentinel is
      !! `< 0`; we expect exactly `nz_ml` non-sentinel entries listed
      !! in `k=1..nz_ml` order (bed → surface).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: layer_rho_init(:)
      integer, intent(in) :: nz_ml
      integer, intent(out), optional :: ierr
         !! Non-zero when the non-sentinel count mismatches `nz_ml`, when
         !! present; absent behaves as today (`error stop`).
      integer :: i, j, k, n_set

      n_set = count(layer_rho_init >= 0.0_wp)
      if (n_set == 0) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return  ! sentinel everywhere — leave rho_layer alone
      end if

      if (n_set /= nz_ml) then
         call fail("ocean_layer_rho_init: count of non-sentinel "// &
                   "entries ("//to_string(n_set)//") does not match "// &
                   "nz_layers ("//to_string(nz_ml)//").", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      do k = 1, nz_ml
         do j = 1, size(ms%rho_layer, 2)
            do i = 1, size(ms%rho_layer, 1)
               ms%rho_layer(i, j, k) = layer_rho_init(k)
            end do
         end do
      end do

      ! Root-only: this is a once-per-run init note, not per-rank state — on a
      ! multi-rank run it would otherwise print nranks identical copies.
      if (comm_env_rank() == 0) then
         call logger%info("ocean_layer_rho_init: wrote per-layer ρ directly, "// &
                          "bypassing EOS init.  ρ(k=1, bed) = "// &
                          to_string(layer_rho_init(1))//" kg/m³, "// &
                          "ρ(k="//to_string(nz_ml)//", surf) = "// &
                          to_string(layer_rho_init(nz_ml))//" kg/m³.")
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine apply_layer_rho_init

   module subroutine seed_eady_ic(state, grid, cfg, ierr)
      !! Eady-front overlay IC.  Assumes flat bottom + uniform layer
      !! split already in place from `ocean_state_seed_from_cfg`.
      !!
      !! Sets:
      !!   T(i, j, k) = T_ref + dT_dz · z(k) + dT_dy · (y_phys(j) - y_mid)
      !!                                                       + ε(i, j, k)
      !!   u_face_x(i, j, k) = dU/dz · (z(k) - z_mid)
      !! with `dU/dz` from thermal-wind balance:
      !!     dU/dz = g · alpha_T / (rho_0 · f) · dT_dy
      !!
      !! Convention: k=1 is the bed, k=nz the surface.  z(k=1) =
      !! -H + dz/2, z(k=nz) = -dz/2.  y_phys uses cell centres.
      !!
      !! The perturbation ε is a uniform-random ±eady_pert_amp/2 noise
      !! seeded from `cfg%ocean%ic%eady_pert_seed`, applied to interior cells
      !! (j ∈ [ng+2, ng+ny_phys-1]) so wall-adjacent rows stay clean.
      !!
      !! Requires `coriolis_f /= 0`.  Errors if topo_config is not "flat".
      !!
      !! MPI: `y_mid` (the front centre) uses the GLOBAL meridional extent
      !! `grid%ny_global`, and the local row index is mapped to a global
      !! physical position with `grid%j_offset_global`.  Built from the
      !! LOCAL extents instead, every rank would put the front in the middle
      !! of its own tile.  On a single rank the offset is 0 and global ==
      !! local, so the seed is byte-identical.
      !!
      !! Not decomposition-invariant: the `eady_pert_amp` noise is drawn
      !! per-rank on the LOCAL array, so its realisation depends on the rank
      !! layout (the deterministic Eady profile underneath does not).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on an Eady-IC configuration conflict when present;
         !! absent behaves as today (`error stop`).

      integer :: i, j, k, nx_total, ny_total, nz_ml, ng
      integer :: ny_phys, j_phys, joff
      real(wp) :: H, dz, y_phys, y_mid, z_k, z_mid
      real(wp) :: dUdz, T_local, hk
      real(wp), allocatable :: noise(:, :, :)
      integer :: idx_T, n_seed, seed_val
      integer, allocatable :: seed_buf(:)

      if (trim(cfg%ocean%topo%topo_config) /= "flat") then
         call fail("seed_eady_ic: requires topo_config='flat'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (abs(cfg%coriolis_f) < tiny(1.0_wp)) then
         call fail("seed_eady_ic: requires coriolis_f /= 0", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = size(state%barotropic%b, 1)
      ny_total = size(state%barotropic%b, 2)
      nz_ml = state%multilayer%nz_ml
      ng = grid%nghost
      ! LOCAL ny_phys — used only to place the noise-free wall-adjacent
      ! rows, which is a per-rank array-bound question, not a coordinate.
      ny_phys = grid%ny_phys
      joff = grid%j_offset_global
      idx_T = state%multilayer%idx_temperature

      H = cfg%ocean%topo%max_depth
      dz = H/real(nz_ml, wp)
      y_mid = 0.5_wp*real(grid%ny_global, wp)*grid%dy
      z_mid = -0.5_wp*H

      ! Thermal-wind balance: ∂u_g/∂z = (g/(fρ₀))·∂ρ/∂y, with linear
      ! EOS ρ = ρ₀ − α(T−T_ref) giving ∂ρ/∂y = −α·∂T/∂y.  Therefore
      !     dU/dz = -g·α/(f·ρ₀) · dT/dy
      ! (cold-to-north dT/dy < 0 → eastward shear with z, surface jet
      !  westerly — the canonical mid-latitude convention).
      !
      ! Geostrophic balance for the y-independent part of -f·u would
      ! also demand a linearly tilted SSH (∂η/∂y = f·dU/dz·z_mid/g).
      ! Setting that η here AND refreshing h_layer accordingly produces
      ! ghost-row h_layer mismatch and tracer-mass corruption in the
      ! continuity-PPM path — the analytic balance doesn't exactly
      ! match the FV-lite PGF's discrete vertical integral.  Left for
      ! follow-up; for now the IC is hydrostatic-only and IG waves
      ! shed in the first ~24h as the system geostrophically adjusts.
      dUdz = -GRAVITY*cfg%ocean%ic%alpha_T/(cfg%ocean%ic%rho_0*cfg%coriolis_f) &
             *cfg%ocean%ic%eady_dT_dy

      ! Pre-generate noise on host (RNG isn't device-callable).  Fill
      ! the full array including ghosts; we zero ghost + bordering rows
      ! after to keep the perturbation cleanly inside the physical
      ! interior.
      allocate (noise(nx_total, ny_total, nz_ml))
      call random_seed(size=n_seed)
      allocate (seed_buf(n_seed))
      seed_val = cfg%ocean%ic%eady_pert_seed
      seed_buf = [(seed_val + 17*j, j=1, n_seed)]
      call random_seed(put=seed_buf)
      call random_number(noise)
      noise = cfg%ocean%ic%eady_pert_amp*(noise - 0.5_wp)
      ! Zero ghosts + first/last physical row to avoid wall-adjacent
      ! noise that would project onto the gravest mode.
      noise(:, 1:ng + 1, :) = 0.0_wp
      noise(:, ng + ny_phys:, :) = 0.0_wp

      ! Override T(i,j,k) with the Eady analytical profile + noise.
      if (idx_T > 0) then
         do k = 1, nz_ml
            z_k = -H + (real(k, wp) - 0.5_wp)*dz
            do j = 1, ny_total
               ! GLOBAL physical row index of this local row.
               j_phys = j - ng + joff
               y_phys = (real(j_phys, wp) - 0.5_wp)*grid%dy
               T_local = cfg%ocean%ic%eady_T_ref + cfg%ocean%ic%eady_dT_dz*z_k &
                         + cfg%ocean%ic%eady_dT_dy*(y_phys - y_mid)
               do i = 1, nx_total
                  hk = state%multilayer%h_layer(i, j, k)
                  state%multilayer%tracers(idx_T)%hTr(i, j, k) = &
                     (T_local + noise(i, j, k))*hk
               end do
            end do
         end do
      end if

      ! Thermal-wind-balanced u_face_x(z).  U depends on z only (no y, x).
      do k = 1, nz_ml
         z_k = -H + (real(k, wp) - 0.5_wp)*dz
         state%multilayer%u_face_x_layer(:, :, k) = dUdz*(z_k - z_mid)
         ! hu on x-faces: interpolate cell-centred h to the u-face (end faces
         ! take the adjacent cell), mirroring the production convention in
         ! face_depth_mean_u.  u_face/h_layer are staggered (nx+1 faces vs nx
         ! cells), so the old whole-array `u_face * h_layer` was non-conformant
         ! (silently OOB on non-checking compilers; LFortran flags it).
         associate (uf => state%multilayer%u_face_x_layer, &
                    hl => state%multilayer%h_layer, &
                    huf => state%multilayer%hu_face_x_layer)
            huf(1, :, k) = uf(1, :, k)*hl(1, :, k)
            huf(2:nx_total, :, k) = uf(2:nx_total, :, k) &
                                    *0.5_wp*(hl(1:nx_total - 1, :, k) + hl(2:nx_total, :, k))
            huf(nx_total + 1, :, k) = uf(nx_total + 1, :, k)*hl(nx_total, :, k)
         end associate
      end do

      deallocate (noise, seed_buf)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_eady_ic

   module subroutine seed_baroclinic_jet_ic(state, grid, cfg, ierr)
      !! Two-layer reduced-gravity baroclinic-instability IC — a
      !! geostrophically-balanced tanh jet in the upper (surface) layer of
      !! a spherical re-entrant channel, plus a front-localised sech²
      !! meander seed that the instability grows.  Reproduces the `bc_inst`
      !! spec (SIM_DETAILS.md §5), mapped to Roundabout's bottom-up layer
      !! convention (k=1 bed, k=nz surface).  The spec's layers FLIP:
      !!
      !!   * spec upper layer 1 (the jet) → Roundabout k=2 (surface)
      !!   * spec lower layer 2 (at rest) → Roundabout k=1 (bed)
      !!
      !! With `y = R·(φ − φ₀)` (m) about the jet centre and
      !! `xm = R·cosφ₀·(λ − λ_west)` (m):
      !!
      !!   ξ(x,y) = Δξ·tanh(y/L) + A_pert·sech²(y/L)·cos(k_x·xm)
      !!   η(y)   = −(g'/g_FS)·ξ                       (free-surface signature)
      !!   h_layer(:,:,1) = H_bed  + ξ                 (lower, k=1)
      !!   h_layer(:,:,2) = H_surf + (η − ξ)           (upper, k=2)
      !!   u_face_x_layer(:,:,2) = (g'/f)·(Δξ/L)·sech²(y/L)  (base tanh only)
      !!   u_face_x_layer(:,:,1) = 0 ;  v_face_y_layer = 0
      !!
      !! `g'` and `g_FS` are read from the SAME reduced-gravity parameters
      !! the gprime PGF consumes (`cfg%ocean%pgf%gprime_gint` /
      !! `gprime_gfs`), so the IC balance and the PGF that maintains it use
      !! one source of truth.  `f = 2Ω·sin(φ)` is the full variable
      !! Coriolis at the u-face (T-row) latitude — NOT the scalar
      !! `coriolis_f`.  The geostrophic jet uses the BASE tanh only; the
      !! perturbation is an unbalanced interface displacement the
      !! instability feeds on.
      !!
      !! Requires `nz_layers == 2`, `&ocean_pgf_nml form="gprime"`, and a
      !! spherical grid.  Fills the FULL arrays incl. ghost rows: φ/λ are
      !! continued analytically past the physical cells, which for a
      !! periodic-x channel wraps exactly because `k_x·L_x = 2π·pert_nx`
      !! (integer `pert_nx`), so the x-ghost columns match the physical
      !! ones.  Wall-y ghosts continue the smooth tanh/sech² profile.
      !!
      !! Host-compute (deterministic, RNG-free): plain host loops writing
      !! `state%multilayer` BEFORE `enter_data`, per the setup-code
      !! convention (a `do concurrent` here would round-trip unmapped
      !! arrays through the device).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a baroclinic-jet IC configuration conflict when
         !! present; absent behaves as today (`error stop`).
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`).  All four are load-bearing: the
      !! jet-centre latitude, the perturbation wavenumber AND both
      !! index→coordinate maps must describe the WHOLE domain, or every
      !! rank reproduces the entire jet inside its own tile (with
      !! `pert_nx = 3`, three wavelengths per TILE instead of three per
      !! DOMAIN).  On a single rank the offsets are 0 and global == local,
      !! so the seed is byte-identical.

      integer :: i, j, k, nx_total, ny_total, nz_ml, ng
      integer :: idx_S, idx_T
      integer :: ioff, joff, nxg, nyg
      real(wp) :: gprime, g_fs, H_total, H_surf, H_bed
      real(wp) :: jet_L, dxi, a_pert, kx, phi0, cos_phi0, Lx
      real(wp) :: rad_earth, omega, lon_west, lat_south, dlon, dlat
      real(wp) :: phi, y_arg, tanh_y, sech2, f_u, u_jet
      real(wp) :: lam_deg, xm, xi_full, eta_local

      nz_ml = state%multilayer%nz_ml
      if (nz_ml /= 2) then
         call fail("seed_baroclinic_jet_ic: requires nz_layers == 2 "// &
                   "(two-layer reduced gravity); got "//to_string(nz_ml), ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (trim(cfg%ocean%pgf%form) /= "gprime") then
         call fail("seed_baroclinic_jet_ic: requires &ocean_pgf_nml "// &
                   "form='gprime'; got '"//trim(cfg%ocean%pgf%form)//"'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (trim(cfg%ocean%grid%grid_config) /= "spherical") then
         call fail("seed_baroclinic_jet_ic: requires &ocean_grid_nml "// &
                   "grid_config='spherical'; got '"// &
                   trim(cfg%ocean%grid%grid_config)//"'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (cfg%ocean%grid%omega <= 0.0_wp) then
         call fail("seed_baroclinic_jet_ic: requires &ocean_grid_nml "// &
                   "omega > 0 (planetary Coriolis f = 2*Omega*sin(phi))", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = grid%nx_total
      ny_total = grid%ny_total
      ng = grid%nghost

      ! Reduced-gravity parameters — read from the gprime PGF's own config
      ! (single source of truth for g' and the free-surface gravity g_FS).
      gprime = cfg%ocean%pgf%gprime_gint
      g_fs = cfg%ocean%pgf%gprime_gfs

      ! Rest layer split: upper (surface, k=2) H_surf; lower (bed, k=1) H_bed.
      H_total = cfg%ocean%topo%max_depth
      H_surf = cfg%ocean%ic%upper_layer_rest
      H_bed = H_total - H_surf
      if (H_surf <= 0.0_wp .or. H_bed <= 0.0_wp) then
         call fail("seed_baroclinic_jet_ic: upper_layer_rest ("// &
                   to_string(H_surf)//" m) must lie in (0, max_depth="// &
                   to_string(H_total)//" m)", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      jet_L = cfg%ocean%ic%jet_half_width
      dxi = cfg%ocean%ic%interface_amp
      a_pert = cfg%ocean%ic%pert_amp_frac*dxi

      ! Spherical geometry.  dx/dy are dlon/dlat in DEGREES for spherical.
      rad_earth = cfg%ocean%grid%rad_earth
      omega = cfg%ocean%grid%omega
      lon_west = cfg%ocean%grid%lon_west
      lat_south = cfg%ocean%grid%lat_south
      dlon = grid%dx
      dlat = grid%dy
      ! Jet-centre latitude = domain centre; φ = (lat_south + ny/2·dlat).
      phi0 = (lat_south + 0.5_wp*real(nyg, wp)*dlat)*DEG2RAD
      cos_phi0 = cos(phi0)

      ! Zonal perturbation wavenumber.  L_x = R·cosφ₀·(lon span, rad);
      ! integer pert_nx ⇒ k_x·L_x = 2π·pert_nx makes the x-continuation
      ! periodic (ghost columns match the physical wrap).
      Lx = rad_earth*cos_phi0*(real(nxg, wp)*dlon*DEG2RAD)
      kx = TWO_PI*real(cfg%ocean%ic%pert_nx, wp)/Lx

      ! ---- h_layer (cell centres), full array incl. ghosts ----
      do j = 1, ny_total
         phi = (lat_south + (real(j - ng + joff, wp) - 0.5_wp)*dlat)*DEG2RAD
         y_arg = rad_earth*(phi - phi0)/jet_L
         tanh_y = tanh(y_arg)
         sech2 = 1.0_wp/cosh(y_arg)**2
         do i = 1, nx_total
            lam_deg = lon_west + (real(i - ng + ioff, wp) - 0.5_wp)*dlon
            xm = rad_earth*cos_phi0*((lam_deg - lon_west)*DEG2RAD)
            xi_full = dxi*tanh_y + a_pert*sech2*cos(kx*xm)
            eta_local = -(gprime/g_fs)*xi_full
            state%multilayer%h_layer(i, j, 1) = H_bed + xi_full            ! lower (bed)
            state%multilayer%h_layer(i, j, 2) = H_surf + (eta_local - xi_full)  ! upper (surface)
         end do
      end do

      ! ---- balanced zonal jet on the upper-layer u-faces (row-only) ----
      ! u depends on the T-row latitude only (base tanh; NOT the pert).  The
      ! u-face shares its T-row's latitude, so f = 2Ω·sin(φ_row).  Lower
      ! layer + all v-faces stay at rest.
      state%multilayer%u_face_x_layer = 0.0_wp
      state%multilayer%v_face_y_layer = 0.0_wp
      do j = 1, ny_total
         phi = (lat_south + (real(j - ng + joff, wp) - 0.5_wp)*dlat)*DEG2RAD
         y_arg = rad_earth*(phi - phi0)/jet_L
         sech2 = 1.0_wp/cosh(y_arg)**2
         f_u = 2.0_wp*omega*sin(phi)
         u_jet = (gprime/f_u)*(dxi/jet_L)*sech2
         do i = 1, nx_total + 1
            state%multilayer%u_face_x_layer(i, j, 2) = u_jet
         end do
      end do

      ! ---- hu on x-faces: u · face-interpolated h (mirror seed_eady_ic) ----
      ! u_face/h_layer are staggered (nx+1 faces vs nx cells); end faces take
      ! the adjacent cell.  Lower layer u=0 ⇒ hu=0.  hv stays zero (v=0).
      state%multilayer%hv_face_y_layer = 0.0_wp
      do k = 1, nz_ml
         associate (uf => state%multilayer%u_face_x_layer, &
                    hl => state%multilayer%h_layer, &
                    huf => state%multilayer%hu_face_x_layer)
            huf(1, :, k) = uf(1, :, k)*hl(1, :, k)
            huf(2:nx_total, :, k) = uf(2:nx_total, :, k) &
                                    *0.5_wp*(hl(1:nx_total - 1, :, k) + hl(2:nx_total, :, k))
            huf(nx_total + 1, :, k) = uf(nx_total + 1, :, k)*hl(nx_total, :, k)
         end associate
      end do

      ! ---- re-seed S/T uniform against the new h_layer ----
      ! The default seed set hTr = const·h_layer against the OLD uniform
      ! split; refresh so hTr/h_layer stays at the configured scalars.
      idx_S = state%multilayer%idx_salinity
      idx_T = state%multilayer%idx_temperature
      if (idx_S > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(idx_S)%hTr, &
            state%multilayer%h_layer, cfg%initial_salinity, nz_ml)
      end if
      if (idx_T > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(idx_T)%hTr, &
            state%multilayer%h_layer, cfg%initial_temperature, nz_ml)
      end if

      if (comm_env_rank() == 0) then
         call logger%info("seed_baroclinic_jet_ic: 2-layer reduced-gravity jet — "// &
                          "L="//to_string(jet_L)//" m, dxi="//to_string(dxi)// &
                          " m, g'="//to_string(gprime)//" m/s^2, H_surf="// &
                          to_string(H_surf)//" m, H_bed="//to_string(H_bed)// &
                          " m, n_x="//to_string(cfg%ocean%ic%pert_nx))
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_baroclinic_jet_ic

   module subroutine seed_geostrophic_adjustment_ic(state, grid, cfg, ierr)
      !! Rossby's classic geostrophic-adjustment problem.  Overlays a
      !! Gaussian SSH bump on a flat-bottom, single-layer (barotropic)
      !! state at rest:
      !!     η(x, y) = A · exp(-r² / L²)
      !!     h(i, j) = b(i, j) + η(i, j)
      !!     u = v = 0
      !! where r is the radial distance from the bump centre.
      !!
      !! Expected evolution (analytical):
      !!   * Length-scale ratio L / Rd with Rd = √(gH)/f governs the
      !!     split between radiated IG-wave energy and retained
      !!     balanced geostrophic ring.
      !!     L >> Rd → most energy retained, surface stays bumped with
      !!                 a balanced anticyclonic flow around it.
      !!     L << Rd → most energy radiates, the bump flattens out.
      !!     L ~ Rd  → partial, ~50/50.
      !!   * Time-series of η at the centre oscillates with period
      !!     2π/f, damping toward the residual value.
      !!   * IG-wave fronts propagate outward at c = √(gH).
      !!
      !! Diagnostic interpretation:
      !!   * Bump amplitude growing in time → numerical instability
      !!     (wrong PGF sign, bad continuity).
      !!   * Oscillation period != 2π/f → wrong Coriolis coupling.
      !!   * Energy not conserved (after damping out the IG transient)
      !!     → spurious sources / sinks in the bt step or Coriolis-adv.
      !!
      !! Requires `topo_config = "flat"` so the bump sits on a
      !! uniform reference depth.  Single-layer is recommended (nz=1)
      !! to isolate the barotropic adjustment; multi-layer works but
      !! the bump distributes uniformly across layers.
      !!
      !! MPI: the default bump centre uses the GLOBAL extents
      !! `grid%nx_global` / `grid%ny_global`, and local cell indices are
      !! mapped to global physical positions with `grid%i_offset_global` /
      !! `grid%j_offset_global`.  Built from the LOCAL extents instead,
      !! every rank would drop a full-amplitude bump in the middle of its
      !! own tile.  On a single rank the offsets are 0 and global == local,
      !! so the seed is byte-identical.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a geostrophic-adjustment IC configuration conflict
         !! when present; absent behaves as today (`error stop`).

      integer :: i, j, nx_total, ny_total, nz_ml, ng
      integer :: i_phys, j_phys, ioff, joff
      real(wp) :: x_phys, y_phys, x_c, y_c, r2, eta_local, inv_L2

      if (trim(cfg%ocean%topo%topo_config) /= "flat") then
         call fail("seed_geostrophic_adjustment_ic: requires topo_config='flat'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (cfg%ocean%ic%ga_length_scale <= 0.0_wp) then
         call fail("seed_geostrophic_adjustment_ic: ga_length_scale must be > 0", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = size(state%barotropic%b, 1)
      ny_total = size(state%barotropic%b, 2)
      nz_ml = state%multilayer%nz_ml
      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ! Default centre = GLOBAL basin midpoint (negative cfg values).
      x_c = cfg%ocean%ic%ga_x_center
      y_c = cfg%ocean%ic%ga_y_center
      if (x_c < 0.0_wp) x_c = 0.5_wp*real(grid%nx_global, wp)*grid%dx
      if (y_c < 0.0_wp) y_c = 0.5_wp*real(grid%ny_global, wp)*grid%dy

      inv_L2 = 1.0_wp/(cfg%ocean%ic%ga_length_scale*cfg%ocean%ic%ga_length_scale)

      ! Override h_total = b + η(x, y) and redistribute the per-layer
      ! thickness evenly.  Velocities + face fluxes left at zero from
      ! the default seed.  Velocities are unchanged at zero, so hu/hv
      ! stay at zero too.
      do j = 1, ny_total
         ! GLOBAL physical indices of this local cell.
         j_phys = j - ng + joff
         y_phys = (real(j_phys, wp) - 0.5_wp)*grid%dy
         do i = 1, nx_total
            i_phys = i - ng + ioff
            x_phys = (real(i_phys, wp) - 0.5_wp)*grid%dx
            r2 = (x_phys - x_c)**2 + (y_phys - y_c)**2
            eta_local = cfg%ocean%ic%ga_eta_amp*exp(-r2*inv_L2)
            state%barotropic%h(i, j) = state%barotropic%b(i, j) + eta_local
         end do
      end do

      call seed_h_layer_uniform_impl(state%multilayer%h_layer, &
                                     state%barotropic%h, nz_ml)
      ! Re-seed S, T so hTr/h_layer stays at the configured uniform
      ! values regardless of the η bump.
      if (state%multilayer%idx_salinity > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
            state%multilayer%h_layer, cfg%initial_salinity, nz_ml)
      end if
      if (state%multilayer%idx_temperature > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
            state%multilayer%h_layer, cfg%initial_temperature, nz_ml)
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_geostrophic_adjustment_ic

   pure subroutine seed_tracer_uniform_impl(hTr, h_layer, tracer_const, nz_ml)
      !! Seed `hTr = const * h_layer` for a uniform-IC tracer (salinity).
      ! assumed-shape-ok: init routine called once at startup; size(hTr,1/2)
      ! used to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in)    :: h_layer(:, :, :)  ! assumed-shape-ok: init routine; size() derives loop bounds
      real(wp), intent(in)    :: tracer_const
      integer, intent(in)    :: nz_ml
      integer :: i, j, k, nx, ny
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do k = 1, nz_ml
         do j = 1, ny
            do i = 1, nx
               hTr(i, j, k) = tracer_const*h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine seed_tracer_uniform_impl

   pure subroutine seed_tracer_stratified_impl(hTr, h_layer, t_layer, nz_ml)
      !! Seed `hTr(i,j,k) = t_layer(k) * h_layer(i,j,k)` for the
      !! per-layer stratified-IC tracer (temperature).  `t_layer(:)` is
      !! a small 1D array pre-computed on host.
      ! assumed-shape-ok: init routine called once at startup; size(hTr,1/2)
      ! used to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in)    :: h_layer(:, :, :)  ! assumed-shape-ok: init routine; size() derives loop bounds
      real(wp), intent(in)    :: t_layer(:)
      integer, intent(in)    :: nz_ml
      integer :: i, j, k, nx, ny
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do k = 1, nz_ml
         do j = 1, ny
            do i = 1, nx
               hTr(i, j, k) = t_layer(k)*h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine seed_tracer_stratified_impl

   pure subroutine seed_wet_mask_impl(wet_mask, b, land_cutoff)
      !! 1.0 where `b >= cutoff` (ocean), 0.0 elsewhere (land).
      !! Default (absent `land_cutoff`): `cutoff = LAND_DEPTH_THRESHOLD`
      !! (byte-identical to the pre-v2 path).  When `land_cutoff` is
      !! present, a column is land iff `b < land_cutoff`; the wetdry-aware
      !! seed passes `land_cutoff = -land_margin` so intertidal columns
      !! (bed above rest MSL but below the flood headroom) stay wet_mask=1
      !! and the dynamic wd_wet_dyn gate handles their wetting/drying.
      ! assumed-shape-ok: init routine called once at startup; size(b,1/2) used
      ! to derive loop bounds.
      real(wp), intent(inout) :: wet_mask(:, :)
      real(wp), intent(in)    :: b(:, :)  ! assumed-shape-ok: init routine; size(b,1/2) derives loop bounds
      real(wp), intent(in), optional :: land_cutoff
         !! When present, a column is land iff `b < land_cutoff` (used by
         !! the wetdry-aware seed: land_cutoff = -land_margin lets intertidal
         !! columns near/above rest MSL stay wet_mask=1). Absent ⇒ the default
         !! `b >= LAND_DEPTH_THRESHOLD` ocean test (byte-identical).
      integer :: i, j, nx, ny
      real(wp) :: cutoff
      nx = size(b, 1)
      ny = size(b, 2)
      cutoff = LAND_DEPTH_THRESHOLD
      if (present(land_cutoff)) cutoff = land_cutoff
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do j = 1, ny
         do i = 1, nx
            if (b(i, j) >= cutoff) then
               wet_mask(i, j) = 1.0_wp
            else
               wet_mask(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine seed_wet_mask_impl

   module subroutine ocean_state_seed_land_cells(state, grid)
      !! Hold land T-cells (`wet_mask==0`) at FINITE reference values so
      !! the masked dyn-core never evaluates `0*NaN` (a zeroed face metric
      !! times a NaN land contribution is still NaN).  Spec §13.3 /
      !! `land_mask_final_resolution.md` R-land-state:
      !!   * `h_layer` floored to `H_VANISHED` (never 0 — avoids 1/0);
      !!   * tracer content zeroed (`hTr = 0`, so `T = S = 0` — finite);
      !!   * layer + barotropic face velocities zeroed at land faces.
      !!
      !! The `(h, hTr) = (H_VANISHED, 0)` pair is the LAND-STATE CONTRACT:
      !! it is what every vanished-gated operator already holds a land
      !! column at, so a land column looks the same to the budget at the
      !! latch and at every later step.  `seed_land_tracer_hold_impl`
      !! carries the full statement and the history.
      !!
      !! Runs at SETUP, after `metrics_apply_land_mask` has derived
      !! `wet_u/wet_v`, BEFORE `ocean_state_enter_data`.  Plain host loops
      !! (pre-enter_data).  No-op for an all-wet domain (`wet_mask≡1`).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid

      integer :: nz_ml, t

      associate (ms => state%multilayer)
         nz_ml = size(ms%h_layer, 3)

         ! Tracer content on land is ZERO — the state every vanished-gated
         ! operator (ALE remap first) holds a land column at, and therefore
         ! the state the budget latch must see.  Order-independent of the
         ! thickness floor below: neither reads the other.  See
         ! `seed_land_tracer_hold_impl`'s land-state contract.
         if (allocated(ms%tracers)) then
            do t = 1, size(ms%tracers)
               call seed_land_tracer_hold_impl(ms%tracers(t)%hTr, &
                                               ms%wet_mask, nz_ml)
            end do
         end if

         ! h_layer floor on land.
         call seed_land_h_floor_impl(ms%h_layer, ms%wet_mask, nz_ml)

         ! On a WET column, the inert fillers a rigid-top coordinate leaves
         ! inside the ice (and below the bed) are seeded with the IC's own
         ! `c(z)*h_min` — the concentration at the filler's depth, not its
         ! donor's.  Establish invariant I1′ (`hTr = h*c_live`) HERE, with
         ! the one definition (`rdb_vl_merge_content`, host twin because
         ! this runs before `enter_data`), so the budget latch sees the
         ! state every later step holds.  Column-conservative, so the
         ! latch's content is unchanged either way; this is about the
         ! seeded state satisfying the invariant, not about the budget.
         ! (Under the previous rule, I1, the first regrid ZEROED these
         ! fillers un-budgeted — `-1.604E-06` relative at step 1 on
         ! `cavity_sloping_lid_rest_zfixed.nml` — which is why this fence
         ! was a zeroing hold.)  Fenced to `z_fixed` x cavity so every
         ! other namelist's step-0 state is bit-identical.
         if (state%metrics%use_cavity .and. state%vcoord%coord_type == VCOORD_Z_FIXED) then
            call ms%enforce_vanished_content_host(size(ms%h_layer, 1), size(ms%h_layer, 2))
         end if

         ! Zero layer face velocities at land faces.
         call seed_land_face_vel_impl(ms%u_face_x_layer, state%metrics%wet_u, nz_ml)
         call seed_land_face_vel_impl(ms%v_face_y_layer, state%metrics%wet_v, nz_ml)
      end associate

      ! Zero the layer mass-transports + barotropic face velocities at land
      ! faces (the transports ride masked metrics anyway, but resetting them
      ! keeps the held land state internally consistent).
      call seed_land_face_vel_impl(state%multilayer%hu_face_x_layer, &
                                   state%metrics%wet_u, nz_ml)
      call seed_land_face_vel_impl(state%multilayer%hv_face_y_layer, &
                                   state%metrics%wet_v, nz_ml)
      call seed_land_face_vel_2d_impl(state%barotropic%u_face_x, state%metrics%wet_u)
      call seed_land_face_vel_2d_impl(state%barotropic%v_face_y, state%metrics%wet_v)
      call seed_land_face_vel_2d_impl(state%barotropic%hu_face_x, state%metrics%wet_u)
      call seed_land_face_vel_2d_impl(state%barotropic%hv_face_y, state%metrics%wet_v)

      ! `grid` is part of the setup contract (shapes come from the state
      ! allocations); not re-read here.
      associate (unused => grid%nx_total)
      end associate
   end subroutine ocean_state_seed_land_cells

   pure subroutine seed_land_tracer_hold_impl(hTr, wet_mask, nz)
      !! Zero the extensive tracer content `hTr` on land T-cells
      !! (`wet_mask == 0`).  Wet cells untouched (bit-identical when
      !! `wet_mask ≡ 1`).
      !!
      !! ### The land-state contract: `h = H_VANISHED`, `hTr = 0`
      !!
      !! `seed_land_h_floor_impl` pins a land layer's thickness EXACTLY at
      !! `H_VANISHED`, which is the D4 *vanished* marker.  Every
      !! vanished-gated operator in the tree tests `h > H_VANISHED`
      !! (strictly), so a land layer is on the VANISHED side of every one
      !! of those gates — in particular the ALE remap's concentration
      !! step (`rdb_ocean_remap::ocean_remap_tracer_field`, `c = hTr/h` if
      !! `h > H_FLOOR` else `c = 0`), which therefore writes `hTr = 0` on
      !! every land column at the first regrid.
      !!
      !! So `hTr = 0` is not a choice made here — it is the land content
      !! the running solver holds.  The seed's job is to hand the budget
      !! latch the SAME land state that every later step will have, and
      !! any other seeded value is content the first regrid discards
      !! silently and un-budgeted, which shows up as a step change in the
      !! console `Error` residual between step 0 and step 1.
      !!
      !! This replaced a "recover `val = hTr/max(h_old, H_VANISHED)` and
      !! re-scale onto the floor" hold, which was wrong twice:
      !!
      !!   * the `max(...)` divisor makes the re-scale an exact IDENTITY
      !!     whenever `h_old <= H_VANISHED` — including `h_old < 0`, which
      !!     is what an ice-shelf column GROUNDED by `&ocean_cavity_dyn_nml
      !!     h_min_cavity` has (its water column `b - z_draft` is negative
      !!     by hundreds of metres).  Those columns kept a FULL-COLUMN,
      !!     negative `hTr` next to a floored `h`, i.e. an implied
      !!     concentration of order `-1e7` PSU, and contributed it to the
      !!     budget latch: 61 % of the initial salt content on
      !!     `validation_examples/ocean/isomip_plus/ocean0_idealised_draft.nml`;
      !!   * even where the re-scale DID work (ordinary land, `0 < h_old`),
      !!     the `val*H_VANISHED` it left is discarded by the first regrid,
      !!     a ~1e-8 relative step change in every land-bearing case.
      !!
      !! The `0*NaN` hazard the old hold existed to avoid is avoided the
      !! same way: `T = S = hTr/h = 0` is finite.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in) :: wet_mask(:, :)
      integer :: i, j, k, nx, ny
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (wet_mask(i, j) == 0.0_wp) hTr(i, j, k) = 0.0_wp
            end do
         end do
      end do
   end subroutine seed_land_tracer_hold_impl

   pure subroutine seed_land_h_floor_impl(h_layer, wet_mask, nz)
      !! Floor land-cell layer thickness to `H_VANISHED` (never 0 ⇒ no
      !! 1/0 in any per-layer divide).  Wet cells untouched.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h_layer(:, :, :)
      real(wp), intent(in) :: wet_mask(:, :)
      integer :: i, j, k, nx, ny
      nx = size(h_layer, 1)
      ny = size(h_layer, 2)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               ! vanished-ok: this IS the land-state contract's definition
               ! site (wet_mask==0 => h_layer=H_VANISHED), not a liveness read.
               if (wet_mask(i, j) == 0.0_wp) h_layer(i, j, k) = H_VANISHED
            end do
         end do
      end do
   end subroutine seed_land_h_floor_impl

   pure subroutine seed_land_face_vel_impl(fld, wet_face, nz)
      !! Zero a per-layer face field (`(nf1,nf2,nz)`) at land faces
      !! (`wet_face==0`).  `wet_face` is the matching `wet_u`/`wet_v`.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: fld(:, :, :)
      real(wp), intent(in) :: wet_face(:, :)
      integer :: i, j, k, nf1, nf2
      nf1 = size(fld, 1)
      nf2 = size(fld, 2)
      do k = 1, nz
         do j = 1, nf2
            do i = 1, nf1
               fld(i, j, k) = fld(i, j, k)*wet_face(i, j)
            end do
         end do
      end do
   end subroutine seed_land_face_vel_impl

   pure subroutine seed_land_face_vel_2d_impl(fld, wet_face)
      !! 2D (barotropic) variant of `seed_land_face_vel_impl`.
      real(wp), intent(inout) :: fld(:, :)
      real(wp), intent(in) :: wet_face(:, :)
      integer :: i, j, nf1, nf2
      nf1 = size(fld, 1)
      nf2 = size(fld, 2)
      do j = 1, nf2
         do i = 1, nf1
            fld(i, j) = fld(i, j)*wet_face(i, j)
         end do
      end do
   end subroutine seed_land_face_vel_2d_impl

   module subroutine seed_cavity_draft(state, grid, cfg, ierr)
      !! Fill `metrics%z_draft` (and its `cover_frac` companion) from
      !! `&ocean_cavity_dyn_nml`, then cross-validate the geometry against
      !! the seeded bathymetry.  Runs from `ocean_state_seed_from_cfg`
      !! IMMEDIATELY after the bathymetry and BEFORE the layer split and
      !! the wet-mask seed, which both read `b − z_draft`.
      !!
      !! Non-`pure` on purpose (the only cavity routine that is): it
      !! reports the grounded / over-land column counts through the
      !! logger and fails loud through the error ring.  The arithmetic it
      !! drives is in `rdb_ocean_cavity`, where every routine IS `pure`.
      !!
      !! UNITS.  The namelist carries metres; the setters work in GRID
      !! coordinate units (metres on Cartesian, DEGREES on
      !! spherical/curvilinear), so the box corners are converted with
      !! `topo_length_to_grid_units` and the dimensionless slope is
      !! converted the inverse way — the same trap `&ocean_topo_nml
      !! slope_scale` documents, where a metres length against a degrees
      !! position collapsed a seamount to a flat basin.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out) :: ierr

      integer :: draft_code, source_code, nx, ny, ng
      integer :: n_over_land, n_grounded, n_interior
      integer :: sign_code, local_ierr
      real(wp) :: per_metre, x0_g, x1_g, y0_g, y1_g, slope_g, amp
      logical :: geographic
         !! Grid positions are in degrees: the input to
         !! `topo_length_to_grid_units` (#203).
      real(wp) :: grounded_frac

      ierr = OCEAN_STATUS_OK
      nx = size(state%metrics%z_draft, 1)
      ny = size(state%metrics%z_draft, 2)
      ng = grid%nghost
      if (nx /= size(state%barotropic%b, 1) .or. ny /= size(state%barotropic%b, 2)) then
         call fail("&ocean_cavity_dyn_nml: z_draft is at its placeholder size — "// &
                   "metrics%use_cavity must be latched BEFORE metrics%init (it is "// &
                   "latched in ocean_state_init_from_config)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      draft_code = parse_cavity_draft_config(cfg%ocean%cavity_dyn%draft_config)
      source_code = parse_cavity_draft_source(cfg%ocean%cavity_dyn%draft_source)
      ! `validate_config` already refuses every spelling outside the v1
      ! envelope with a full explanation; this is the same gate one level
      ! down, for direct (test / API) callers that bypass it.
      if (draft_code /= CAVITY_DRAFT_NONE .and. draft_code /= CAVITY_DRAFT_FLAT &
          .and. draft_code /= CAVITY_DRAFT_LINEAR .and. draft_code /= CAVITY_DRAFT_FILE) then
         call fail("&ocean_cavity_dyn_nml draft_config='"// &
                   trim(cfg%ocean%cavity_dyn%draft_config)//"' is not available "// &
                   "(none|flat|linear|file)", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (source_code /= CAVITY_SOURCE_DRAFT .and. source_code /= CAVITY_SOURCE_THICKNESS) then
         call fail("&ocean_cavity_dyn_nml draft_source='"// &
                   trim(cfg%ocean%cavity_dyn%draft_source)//"' is not available "// &
                   "(draft|thickness; 'in_situ' isostasy is deferred)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ! GRID units per metre (1 on Cartesian; degrees-per-metre otherwise).
      ! The "no limit" sentinels pass through UNCONVERTED so they stay
      ! sentinels on every grid.
      ! Inline intrinsic comparisons, not a user call with `grid_config`:
      ! see the #203 note in `ocean_state_seed_from_cfg`.
      geographic = cfg%ocean%grid%grid_config == "spherical" .or. &
                   cfg%ocean%grid%grid_config == "supergrid" .or. &
                   cfg%ocean%grid%grid_config == "tripolar"
      per_metre = topo_length_to_grid_units(1.0_wp, geographic, cfg%ocean%grid%rad_earth)
      x0_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_x0, per_metre)
      x1_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_x1, per_metre)
      y0_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_y0, per_metre)
      y1_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_y1, per_metre)
      ! `draft_slope` is d(draft [m]) / d(x [m]); the setter wants
      ! d(draft [m]) / d(x [grid units]) = slope / (grid units per metre).
      slope_g = cfg%ocean%cavity_dyn%draft_slope/per_metre

      ! `draft_source = "thickness"`: the formula amplitude is an ice
      ! THICKNESS, converted by the Boussinesq-isostatic (flotation)
      ! relation `z_draft = rho_ice*h_ice/rho_0`.  Scaling the amplitude
      ! (and the slope with it) is exact because both setters are LINEAR
      ! in the amplitude — and it keeps one draft field, so everything
      ! downstream stays source-agnostic.
      amp = cfg%ocean%cavity_dyn%draft_depth
      if (source_code == CAVITY_SOURCE_THICKNESS) then
         amp = amp*cfg%ocean%cavity_dyn%rho_ice/state%eos%rho0
         slope_g = slope_g*cfg%ocean%cavity_dyn%rho_ice/state%eos%rho0
      end if

      select case (draft_code)
      case (CAVITY_DRAFT_FLAT)
         call set_draft_flat(state%metrics%z_draft, grid, amp, x0_g, x1_g, y0_g, y1_g)
      case (CAVITY_DRAFT_LINEAR)
         call set_draft_linear(state%metrics%z_draft, grid, amp, slope_g, &
                               x0_g, x1_g, y0_g, y1_g)
      case (CAVITY_DRAFT_FILE)
#ifndef RDB_NO_NETCDF
         ! Static 2-D NetCDF draft, through the PR-14 reader.  SINGLE
         ! RANK: the loader itself applies the global offset correctly,
         ! but the grounding statistics a few lines below are single-rank
         ! reductions and the whole cavity is fenced that way, so the
         ! restriction is asserted here rather than left implicit.
         if (grid%nx_phys /= grid%nx_global .or. grid%ny_phys /= grid%ny_global) then
            call fail("&ocean_cavity_dyn_nml draft_config='file' is single-rank "// &
                      "only (the cavity's grounding statistics are single-rank "// &
                      "reductions).  Run on one rank or use an analytic draft.", &
                      ierr, OCEAN_STATUS_ERR_IC_SEED)
            return
         end if
         sign_code = parse_cavity_draft_sign(cfg%ocean%cavity_dyn%draft_sign)
         if (sign_code == CAVITY_SIGN_INVALID) then
            call fail("&ocean_cavity_dyn_nml draft_sign='"// &
                      trim(adjustl(cfg%ocean%cavity_dyn%draft_sign))// &
                      "' is not recognised (depth|positive_down|elevation|"// &
                      "positive_up)", ierr, OCEAN_STATUS_ERR_IC_SEED)
            return
         end if
         ! Interior first (the reader writes the physical window only)...
         state%metrics%z_draft = 0.0_wp
         call ocean_data_input_load_static_2d( &
            trim(cfg%ocean%cavity_dyn%draft_file), &
            trim(cfg%ocean%cavity_dyn%draft_var), grid, nx, ny, &
            ng + 1, ng + 1, state%metrics%z_draft, ierr=local_ierr)
         if (local_ierr /= OCEAN_STATUS_OK) then
            ierr = local_ierr
            return
         end if
         ! ...sign-normalise onto DEPTH positive down...
         call cavity_draft_apply_sign(state%metrics%z_draft, nx, ny, sign_code)
         ! ...then fill the ghost band by constant extrapolation, the
         ! SAME routine and the same order the file bathymetry uses
         ! (`load_bathymetry_into_array` -> `fill_bathymetry_ghosts_array`).
         ! The periodic/fold re-wrap and the halo exchange that
         ! `metrics%z_draft` gets in `rdb_ocean_engine` run later and are
         ! shared with the formula path, so a file draft and a formula
         ! draft see an identical boundary treatment.
         call bathymetry_fill_ghosts_array(state%metrics%z_draft, grid)
#else
         call fail("&ocean_cavity_dyn_nml draft_config='file' requires "// &
                   "RDB_ENABLE_NETCDF=ON at build time (the static-2-D reader "// &
                   "lives in the NetCDF-gated rdb_ocean_data_input).", &
                   ierr, OCEAN_STATUS_ERR_IO)
         return
#endif
      case default  ! CAVITY_DRAFT_NONE
         state%metrics%z_draft = 0.0_wp
      end select

      ! Guard the field itself before anything derives geometry from it.
      ! Written as `.not. (z >= 0)` inside the helper so a NaN FAILS
      ! rather than sliding through a `z < 0` test that is false for NaN.
      if (.not. cavity_draft_is_finite_nonneg(state%metrics%z_draft, nx, ny)) then
         call fail("&ocean_cavity_dyn_nml: z_draft must be finite and >= 0 "// &
                   "everywhere (it is a DEPTH below z = 0, positive down)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ! NO ICE OVER LAND: zero the draft wherever the bathymetry already
      ! says land, so land columns keep the datum they always had
      ! (`bt_H_ref = b`) and the counted-once invariant stays exact on
      ! every column.
      call cavity_apply_land_exclusion(state%metrics%z_draft, state%barotropic%b, &
                                       nx, ny, n_over_land)

      call cavity_count_grounded(state%barotropic%b, state%metrics%z_draft, &
                                 cfg%ocean%cavity_dyn%h_min_cavity, ng, &
                                 grid%nx_phys, grid%ny_phys, nx, ny, &
                                 n_grounded, n_interior)
      grounded_frac = real(n_grounded, wp)/real(max(n_interior, 1), wp)
      if (grounded_frac > cfg%ocean%cavity_dyn%grounded_max_frac) then
         call fail("&ocean_cavity_dyn_nml: the prescribed draft grounds "// &
                   to_string(n_grounded)//" of "//to_string(n_interior)// &
                   " interior columns ("//to_string(grounded_frac)//"), above "// &
                   "grounded_max_frac = "// &
                   to_string(cfg%ocean%cavity_dyn%grounded_max_frac)// &
                   ".  Either the draft is too deep for this bathymetry or the "// &
                   "shelf box is in the wrong place (check the UNITS of "// &
                   "draft_x0/x1 — they are metres, converted to grid units).", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      call cavity_fill_cover_frac(state%metrics%cover_frac, state%metrics%z_draft, nx, ny)

      if (comm_env_rank() == 0) then
         call logger%info("Ice-shelf cavity: ON  draft_config='"// &
                          trim(cfg%ocean%cavity_dyn%draft_config)//"' source='"// &
                          trim(cfg%ocean%cavity_dyn%draft_source)//"' max draft = "// &
                          to_string(maxval(state%metrics%z_draft))//" m, "// &
                          "h_min_cavity = "// &
                          to_string(cfg%ocean%cavity_dyn%h_min_cavity)//" m")
         call logger%info("                  datum bt_H_ref = b - z_draft afloat, "// &
                          "0 where grounded; "// &
                          to_string(n_grounded)//" of "//to_string(n_interior)// &
                          " interior columns grounded (-> LAND via the wet mask)")
         if (n_over_land > 0) then
            call logger%warning("&ocean_cavity_dyn_nml: draft zeroed on "// &
                                to_string(n_over_land)//" column(s) whose bed is "// &
                                "already land (b < LAND_DEPTH_THRESHOLD) — no ice "// &
                                "over land.  Check the shelf box if that is a surprise.")
         end if
      end if
   end subroutine seed_cavity_draft

   pure function cavity_bound_to_grid(bound_m, per_metre) result(bound_grid)
      !! Convert ONE shelf-box bound from metres to grid coordinate units,
      !! leaving the `CAVITY_BOUND_INF` "no limit" sentinel alone.  Without
      !! the guard the sentinel would be scaled by the degrees-per-metre
      !! factor on a spherical grid and come out as a finite (if absurd)
      !! bound — harmless numerically, but it would stop meaning what it
      !! says, and the next reader would have to re-derive that.
      real(wp), intent(in) :: bound_m
      real(wp), intent(in) :: per_metre
         !! Grid units per metre (1 on Cartesian).
      real(wp) :: bound_grid
      if (abs(bound_m) >= CAVITY_BOUND_INF) then
         bound_grid = bound_m
      else
         bound_grid = bound_m*per_metre
      end if
   end function cavity_bound_to_grid

end submodule rdb_ocean_state_seed
