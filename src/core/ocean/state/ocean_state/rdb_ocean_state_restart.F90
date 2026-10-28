submodule(rdb_ocean_state) rdb_ocean_state_restart
   !! R2b split: host-only restart registry/metadata/read bodies,
   !! moved verbatim from rdb_ocean_state.F90 (restart WRITE stays in
   !! the parent -- it is the D->H pull side; see
   !! src/core/ocean/state/rdb_ocean_state.F90 for the public interface).
   implicit none

contains

   module subroutine ocean_state_build_restart_registry(state, grid, reg)
      !! Walk the ocean god state and register every field that must
      !! checkpoint for a bit-exact step-(N+1) resume (ROADMAP A1).
      !!
      !! Registered (always):
      !!   barotropic prognostics h, u_face_x, v_face_y;
      !!   multilayer prognostics h_layer, u/v_face_x/y_layer, every
      !!   tracers(:)%hTr (registry-of-registries — iterate the tracer
      !!   registry so passive tracers join automatically).
      !! Registered (conditionally — only when the slot's persistent
      !!   array is allocated):
      !!   vmix%bl_depth — KPP's lagged BL-depth seed; the V_t² wstar3
      !!     term reads the PREVIOUS step's value, so a restart between
      !!     steps reads it stale unless checkpointed.  Allocated
      !!     unconditionally (alloc tracks slot existence, not whether the
      !!     closure is active), so it always registers.
      !!   epbl%mld + epbl%kd_int, kshear%kd_int + kshear%tke_int.
      !!     These refresh only at thermo cadence but are *merged every
      !!     stage*, so a restart between refreshes reads them stale.
      !!     Also allocated unconditionally — register-always; the
      !!     enable-state itself is fixed at config and validated by the
      !!     schema metadata, so a flipped flag is a fresh run, not a
      !!     resume (flags cannot flip mid-run).
      !!   ml_rho_layer — diagnosed each stage from S/T, BUT when
      !!     `dt_therm_ratio > 1` it is NOT recomputed on the slow steps
      !!     that skip the EOS, so the value carried into the next step
      !!     depends on the last thermo step.  Cheapest correct fix:
      !!     checkpoint it (review #2).
      !!   bc%eta_old_chapman_{w,e,s,n} — Chapman radiation persistent
      !!     corner state (host scalars).  bc%tres_* tracer reservoirs —
      !!     registered when allocated (device-mapped, so they join the
      !!     update-self walk).
      !!
      !! Deliberately EXCLUDED:
      !!   sf%Q_heat / sf%Q_salt with the surface-flux component set OFF
      !!     (default) — configure-static: re-seeded from `&ocean_thermo_nml
      !!     q_heat / q_salt` on every resume, and the sea-ice couplers'
      !!     contributions are folded back in from the registered
      !!     `ice_salt_flux_diag` / `ice_heat_flux_diag` / `ice_sw_thru_diag`
      !!     (`engine_setup`), exactly the sums the couplers form.  With the
      !!     component set ON they ARE registered (`sf_Q_heat`/`sf_Q_salt`,
      !!     see below): `ocean_surface_flux_assemble` rebuilds them at the
      !!     END of each thermo step, and every step until the next one
      !!     reads that assembly -- carried state, with SST-dependent terms
      !!     no configure-time re-assembly could reproduce.
      !!   RK saves (h_layer0, hTr0, *_layer0) + bt accumulators
      !!     (bt_eta/ubt/...) — step-internal scratch, re-zeroed at the
      !!     top of every outer step (restart is step-aligned).
      !!   bt_H_ref — recomputed from `b` by configure_ocean_bt_split;
      !!     `b` is bathymetry (static), reconstructed at config time.
      !!   metrics%z_draft / cover_frac / p_ice_ref — the ice-shelf cavity
      !!     statics (`&ocean_cavity_dyn_nml`) fall under the SAME rule,
      !!     and for the same reason: the draft is a PRESCRIBED, static
      !!     geometry, rebuilt at configure by `seed_cavity_draft` from
      !!     the namelist before the restart read runs, exactly as `b` and
      !!     `bt_H_ref` are.  Checkpointing it would create the
      !!     file-vs-namelist ambiguity the derived-field exclusion exists
      !!     to avoid (does a saved draft beat an edited namelist?), and a
      !!     resume whose draft disagreed with its datum would be
      !!     silently wrong — `bt_H_ref = b − z_draft` ties the two
      !!     together, so they must be rebuilt together or not at all.
      !!     A TIME-VARYING draft (a coupled ice sheet) is a different
      !!     field with a different owner and would register itself, the
      !!     same way a time-varying surface-flux component must.
      !!   vcoord target_h / z_ref — recomputed per step (z_ref rebuilt
      !!     from `b` in the seed); never prognostic.
      !!   w_interface, mass_flux_* — diagnosed each stage from the
      !!     prognostics before first use.
      !!   diag accumulators — MEAN/MAX windows RESET on restart by
      !!     design (documented); the bit-exact gate covers the
      !!     prognostic trajectory, not mid-window diag aggregates.
      !!
      !! The dead `ocean_obc_t` scaffold (`state%obc`) is NOT registered:
      !! it is never enabled, never device-mapped, and has no enter_data —
      !! registering it (and issuing `update self` on its host-only
      !! buffers) was a latent GPU crash (review #3/#8).  When OBC lands,
      !! it registers its own live persistent state here.
      type(ocean_state_t), intent(in), target :: state
      type(hgrid_t), intent(in) :: grid
      type(restart_registry_t), intent(inout) :: reg

      integer :: it
      character(len=64) :: tag

      call reg%clear()

      ! Ghost-cell policy (v1): the structured prognostic + closure
      ! arrays are checkpointed as FULL local arrays (interior + ghosts),
      ! registered with ng=0 over the total extent.  Rationale — three
      ! kinds of ghost cell, distinguished for the MPI-readiness contract:
      !   * physical domain-edge (wall) ghosts are GENUINE owned boundary
      !     state — the slow-tendency operators (EOS / face-thickness /
      !     hvisc) read `h_layer` etc. at wall-adjacent ghost faces and
      !     there is NO stage-entry wall-mirror op to re-establish them,
      !     so owned-only loses information and breaks the bit-exact gate.
      !   * periodic / tripolar-fold seam ghosts are redundant — the
      !     stage-entry wrap (`ocean_periodic_wrap_state` / fold) rebuilds
      !     them from the interior every step; saving them is harmless.
      !   * inter-rank halo ghosts (E1, multi-rank) will be refilled by
      !     the halo exchange on resume; under the v1 SAME-decomp rule
      !     those columns round-trip identically anyway.
      ! So a full-array write is bit-exact today and forward-compatible
      ! with E1 (the halo columns simply get overwritten by exchange).
      ! The owned-slice machinery in restart_entry_t is retained for the
      ! ghostless edge buffers (OBC rings) and for a future owned-only
      ! mode once a wall-ghost-rebuild op exists.
      ! --- Barotropic prognostics ---
      call reg%register_2d("bt_h", state%barotropic%h, 0, &
                           size(state%barotropic%h, 1), size(state%barotropic%h, 2))
      call reg%register_2d("bt_u_face_x", state%barotropic%u_face_x, 0, &
                           size(state%barotropic%u_face_x, 1), size(state%barotropic%u_face_x, 2))
      call reg%register_2d("bt_v_face_y", state%barotropic%v_face_y, 0, &
                           size(state%barotropic%v_face_y, 1), size(state%barotropic%v_face_y, 2))

      if (state%use_multilayer) then
         ! --- Multilayer prognostics ---
         call register_full_3d(reg, "ml_h_layer", state%multilayer%h_layer)
         call register_full_3d(reg, "ml_u_face_x_layer", state%multilayer%u_face_x_layer)
         call register_full_3d(reg, "ml_v_face_y_layer", state%multilayer%v_face_y_layer)

         ! --- Tracers (registry-of-registries) ---
         if (allocated(state%multilayer%tracers)) then
            do it = 1, size(state%multilayer%tracers)
               write (tag, "(A,I0)") "tracer_hTr_", it
               call register_full_3d(reg, trim(tag), state%multilayer%tracers(it)%hTr)
            end do
         end if
         ! --- pred_corr step time-means (SPEC S1).  Under
         !     `&ocean_bt_nml split_scheme="pred_corr"` these are CARRIED
         !     STATE, not scratch: the predictor's Coriolis-advection and
         !     horizontal viscosity read u_av/v_av/h_av, and the dyn loop
         !     seeds them from the initial state only at
         !     `outer_step_count == 0`.  A resume restores outer_step_count
         !     from the file, so that seed does NOT re-fire -- without
         !     checkpointing them the first post-restart predictor reads
         !     the allocation zeros (a zero-thickness h_av) and the resume
         !     is not bit-exact.  Caught by `restart_bit_exact_*` the day
         !     pred_corr became the default.  `optional=.true.`: a pre-flip
         !     checkpoint has no such variables and must still resume --
         !     it will be an ssp_rk2 file, where they are unread.
         if (allocated(state%multilayer%u_av_layer)) then
            call reg%register_3d("ml_u_av_layer", state%multilayer%u_av_layer, 0, &
                                 size(state%multilayer%u_av_layer, 1), &
                                 size(state%multilayer%u_av_layer, 2), optional=.true.)
         end if
         if (allocated(state%multilayer%v_av_layer)) then
            call reg%register_3d("ml_v_av_layer", state%multilayer%v_av_layer, 0, &
                                 size(state%multilayer%v_av_layer, 1), &
                                 size(state%multilayer%v_av_layer, 2), optional=.true.)
         end if
         if (allocated(state%multilayer%h_av_layer)) then
            call reg%register_3d("ml_h_av_layer", state%multilayer%h_av_layer, 0, &
                                 size(state%multilayer%h_av_layer, 1), &
                                 size(state%multilayer%h_av_layer, 2), optional=.true.)
         end if
         ! --- pred_corr carried viscous tendency.  The predictor does NOT
         !     recompute the lateral viscosity: it reuses the previous
         !     step's corrector `du_visc`/`dv_visc` (MOM6 `diffu(u[n-1])`,
         !     `rdb_ocean_dyn`'s `is_pred` gate), so these buffers are
         !     CARRIED STATE under pred_corr, not scratch.  Unregistered,
         !     the first post-restart predictor read the device-zeroed
         !     payload -- one step with no lateral viscosity in any column
         !     -- and the resume was not bit-exact (found on the 1/4-degree
         !     Southern Ocean, 2026-10-01; the restart gate ran inviscid).
         !     Under ssp_rk2 they are recomputed before every read, so
         !     restoring them is inert.  `optional=.true.`: an older
         !     checkpoint still resumes, with the old one-step defect.
         if (allocated(state%hvisc%du_visc%data)) then
            call register_full_3d_opt(reg, "hvisc_du_visc", state%hvisc%du_visc%data)
         end if
         if (allocated(state%hvisc%dv_visc%data)) then
            call register_full_3d_opt(reg, "hvisc_dv_visc", state%hvisc%dv_visc%data)
         end if
      end if

      ! --- rho_layer (review #2): diagnosed each stage, but NOT
      !     recomputed on dt_therm-skipped slow steps, so the value
      !     carried forward is thermo-step-dependent.  Checkpoint it. ---
      if (state%use_multilayer .and. allocated(state%multilayer%rho_layer)) then
         call register_full_3d(reg, "ml_rho_layer", state%multilayer%rho_layer)
      end if

      ! --- KPP lagged BL-depth seed (review #1): the V_t² wstar3_lagged
      !     term reads the previous step's bl_depth.  Allocated
      !     unconditionally; required for a bit-exact KPP resume. ---
      if (allocated(state%vmix%bl_depth)) then
         call reg%register_2d("vmix_bl_depth", state%vmix%bl_depth, 0, &
                              size(state%vmix%bl_depth, 1), size(state%vmix%bl_depth, 2))
      end if

      ! --- PR-2 (bt-rem-from-av-rem): vmix%kv, closing the REAL root cause
      !     of compat row `restart_visc_rem`.  `visc_rem_precompute` runs
      !     at the START of every pred_corr stage (MOM6-order parity,
      !     PGF_BUG.md §9), BEFORE that stage's own `vmix_apply_in_stage`
      !     recomputes `vmix%kv` -- so it always reads the PREVIOUS
      !     stage's `kv`, one stage stale by design (mirrors MOM6:
      !     `vertvisc_coef`'s `visc%Kv_slow`/the shear-viscosity input is
      !     likewise carried across the predictor/corrector boundary, and
      !     MOM6 checkpoints exactly this class of field --
      !     `set_visc_register_restarts`
      !     registers `Kv_shear`/`Kd_shear`/`Kv_shear_Bu`/`MLD` for the
      !     same "read before recomputed" reason).  Unregistered, a warm
      !     restart re-seeds `kv` from the COLD background value
      !     (`vmix_seed_backgrounds`, `pp81_nu_bg`) rather than the
      !     spun-up profile, so the first resumed stage's visc_rem matrix
      !     (hence `F_bt`'s weighting under `bt_forcing_visc_rem`/
      !     `bt_renorm_visc_rem`, and `bt_rem_from_visc_rem`'s av_rem/
      !     bt_rem) differs from the continued run --
      !     every prognostic then drifts (measured 1e-11 relative by step
      !     24, `tests/regression/compat_expect.py::restart_visc_rem`).
      !     `kv` is allocated UNCONDITIONALLY by `ocean_vmix_init`
      !     (every vmix consumer shares it, not just visc_rem), so this
      !     registers on every run, not just visc_rem-chain ones --
      !     cheap and simple beats a consumer-flag gate here: one extra
      !     `(nx, ny, nz+1)` field, the same class of cost as
      !     `epbl_kd_int` two lines below.  `optional=.true.`: a
      !     pre-PR-2 checkpoint has no such field and must still resume
      !     (re-seeding `kv` cold, the pre-existing one-stage-stale
      !     defect, same compat posture as `bt_visc_rem_u/v` above). ---
      if (allocated(state%vmix%kv)) then
         call register_full_3d_opt(reg, "vmix_kv", state%vmix%kv)
      end if

      ! --- EPBL prev-MLD seed + kd_int (merged every stage) ---
      if (allocated(state%epbl%mld)) then
         call reg%register_2d("epbl_mld", state%epbl%mld, 0, &
                              size(state%epbl%mld, 1), size(state%epbl%mld, 2))
      end if
      if (allocated(state%epbl%kd_int)) then
         call register_full_3d(reg, "epbl_kd_int", state%epbl%kd_int)
      end if

      ! --- MEKE prognostic eddy-energy field (capability [5]): without
      !     this it would cold-reset to 0 each restart, losing the spun-up
      !     eddy field.  Allocated whenever the multilayer path is on. ---
      if (allocated(state%meke%meke)) then
         call reg%register_2d("meke", state%meke%meke, 0, &
                              size(state%meke%meke, 1), size(state%meke%meke, 2))
      end if

      ! --- Gent-McWilliams carried state.  The GM operator runs AFTER the
      !     dynamics of step n (`run_gm_step`) and leaves its PE release
      !     `gm_src`, which MEKE reads at the top of step n+1: carried, so
      !     a resume between the two must restore it (without it the first
      !     resumed MEKE step sourced from a cold GM work, ~3 % off at step
      !     24).  With `dt_therm_ratio > 1` the GM operator also runs on
      !     the steps BETWEEN thermo refreshes, reading the slopes and the
      !     VarMix(+MEKE) base KhTh the last thermo step left — carried as
      !     well (at ratio 1 they are recomputed before every read, so
      !     restoring them is inert).  Registered on EVERY GM run: this
      !     registry is built for the read BEFORE `configure_ocean_vmix`
      !     sets `dyn%dt_therm_ratio`, so the ratio cannot gate it (a
      !     ratio-gated registration wrote them and never read them back —
      !     `restart_engine_bit_exact_gm_meke_varmix_therm2`).  All
      !     optional: a checkpoint written before they were registered
      !     still resumes, with the old one-step defect. ---
      if (state%gm%enable .and. allocated(state%gm%gm_src)) then
         call reg%register_2d("gm_src", state%gm%gm_src, 0, &
                              size(state%gm%gm_src, 1), size(state%gm%gm_src, 2), optional=.true.)
         if (allocated(state%slopes%slope_x)) then
            call register_full_3d_opt(reg, "slopes_slope_x", state%slopes%slope_x)
            call register_full_3d_opt(reg, "slopes_slope_y", state%slopes%slope_y)
            call register_full_3d_opt(reg, "slopes_n2_u", state%slopes%n2_u)
            call register_full_3d_opt(reg, "slopes_n2_v", state%slopes%n2_v)
         end if
         if (state%varmix%enable .and. allocated(state%varmix%khth_u)) then
            call reg%register_2d("varmix_khth_u", state%varmix%khth_u, 0, &
                                 size(state%varmix%khth_u, 1), size(state%varmix%khth_u, 2), &
                                 optional=.true.)
            call reg%register_2d("varmix_khth_v", state%varmix%khth_v, 0, &
                                 size(state%varmix%khth_v, 1), size(state%varmix%khth_v, 2), &
                                 optional=.true.)
         end if
      end if

      ! --- Wet/dry hysteresis mask (docs/ocean_wetdry_plan.md): persistent
      !     front state (wet/dry/held-in-band).  optional — a restart
      !     written before the knob existed re-seeds from depth at
      !     configure instead of failing. ---
      if (allocated(state%dyn%bt_work%wd_wet_dyn)) then
         call reg%register_2d("wd_wet_dyn", state%dyn%bt_work%wd_wet_dyn, 0, &
                              size(state%dyn%bt_work%wd_wet_dyn, 1), &
                              size(state%dyn%bt_work%wd_wet_dyn, 2), optional=.true.)
      end if

      ! --- PR-1 viscous remnant γ (`bt_work%visc_rem_u/v`): REFRESHED,
      !     not re-derived from scratch, by the stage-end vdiff producer
      !     (`bt_visc_rem_producer`, D1 follow-up -- decoupled from the
      !     retired `bt_correction_visc_rem`) and the pre-substep
      !     `visc_rem_precompute` (`forcing_visc_rem`/`renorm_visc_rem`/
      !     `bt_rem_from_visc_rem`) -- a cold resume without this
      !     checkpoint would restart every post-restart stage from the
      !     `source=1.0` init value, which is one stage's worth of
      !     refresh behind a continued run (consumers read the PREVIOUS
      !     stage's γ by design -- see `bt_visc_rem_producer`'s
      !     docstring in `rdb_barotropic_workstate`).  Allocated
      !     unconditionally by
      !     `barotropic_workstate_t%init`.  `optional=.true.`: a pre-PR-1
      !     checkpoint has no such variables and must still resume -- it
      !     re-seeds at 1.0, same as a cold start (closes compat row
      !     `restart_visc_rem`). ---
      if (allocated(state%dyn%bt_work%visc_rem_u)) then
         call register_full_3d_opt(reg, "bt_visc_rem_u", state%dyn%bt_work%visc_rem_u)
      end if
      if (allocated(state%dyn%bt_work%visc_rem_v)) then
         call register_full_3d_opt(reg, "bt_visc_rem_v", state%dyn%bt_work%visc_rem_v)
      end if

      ! --- Sea-ice frazil bank (PR 1): un-spent supercooling heat the ice
      !     model (PR 3) will consume.  PERSISTENT (accumulates across
      !     steps), so a restart must carry it or banked energy would be
      !     silently discarded.  optional — only ice-on runs write it, and
      !     an ice-on resume from an older checkpoint re-seeds 0. ---
      if (allocated(state%ice%frazil_heat)) then
         call reg%register_2d("ice_frazil_heat", state%ice%frazil_heat, 0, &
                              size(state%ice%frazil_heat, 1), &
                              size(state%ice%frazil_heat, 2), optional=.true.)
      end if

      ! --- Sea-ice PR 3b: brine-rejection flux diag.  `salt_flux_diag` is
      !     REQUIRED for restart-exact Q_salt refill: `rdb_ice_ocean_coupler`
      !     rebuilds Q_salt from this field every thermo step, and the
      !     driver's configure-time resume fold re-applies it so the first
      !     post-resume window sees the same flux the uninterrupted run
      !     would have (Q_salt itself is configure-static and NOT
      !     registered).  `m_frozen_diag` is registered alongside for
      !     post-resume diagnostic continuity only.  optional — an
      !     ice-on resume from a pre-3b checkpoint re-seeds 0 (a one-window
      !     cold-start of the brine flux) rather than failing.
      !
      !     PR 3c extends the same contract to the melt side:
      !     `heat_flux_diag` is REQUIRED for restart-exact Q_heat refill
      !     (mirrors `salt_flux_diag` — `ice_ocean_heat_flux` rebuilds
      !     Q_heat from it every thermo step, and the driver resume fold
      !     re-applies it); `m_melt_diag` is diagnostic-continuity only,
      !     like `m_frozen_diag`.  PR 31 extends the same required-refill
      !     contract to `sw_thru_diag` — `ice_ocean_sw_flux` rebuilds the
      !     ocean shortwave from it every thermo step and the driver resume
      !     fold re-applies it, so it is registered alongside
      !     `heat_flux_diag`.  The coupleable atmospheric-forcing seam
      !     (`atm_sf0`/`atm_dsfdt`/`atm_sw_dn`) and the column-driver
      !     scratch (`fb`/`sst_seam`/`ssurf_seam`/`tfw_seam`/`tsurf_out`/
      !     `h2o_ocn_to_ice`/`h2o_ice_to_ocn`/`heat_to_ocn`/`sw_thru`) are
      !     ALL recomputed every thermo step — deliberately NOT
      !     registered (the per-category `sw_thru` stays unregistered; only
      !     its per-cell reduction `sw_thru_diag` is carried). ---
      if (allocated(state%ice%m_frozen_diag)) then
         call reg%register_2d("ice_m_frozen_diag", state%ice%m_frozen_diag, 0, &
                              size(state%ice%m_frozen_diag, 1), &
                              size(state%ice%m_frozen_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%salt_flux_diag)) then
         call reg%register_2d("ice_salt_flux_diag", state%ice%salt_flux_diag, 0, &
                              size(state%ice%salt_flux_diag, 1), &
                              size(state%ice%salt_flux_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_melt_diag)) then
         call reg%register_2d("ice_m_melt_diag", state%ice%m_melt_diag, 0, &
                              size(state%ice%m_melt_diag, 1), &
                              size(state%ice%m_melt_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_congel_win)) then
         call reg%register_2d("ice_m_congel_win", state%ice%m_congel_win, 0, &
                              size(state%ice%m_congel_win, 1), &
                              size(state%ice%m_congel_win, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_meltt_win)) then
         call reg%register_2d("ice_m_meltt_win", state%ice%m_meltt_win, 0, &
                              size(state%ice%m_meltt_win, 1), &
                              size(state%ice%m_meltt_win, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_meltb_win)) then
         call reg%register_2d("ice_m_meltb_win", state%ice%m_meltb_win, 0, &
                              size(state%ice%m_meltb_win, 1), &
                              size(state%ice%m_meltb_win, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_snoice_win)) then
         call reg%register_2d("ice_m_snoice_win", state%ice%m_snoice_win, 0, &
                              size(state%ice%m_snoice_win, 1), &
                              size(state%ice%m_snoice_win, 2), optional=.true.)
      end if
      if (allocated(state%ice%heat_flux_diag)) then
         call reg%register_2d("ice_heat_flux_diag", state%ice%heat_flux_diag, 0, &
                              size(state%ice%heat_flux_diag, 1), &
                              size(state%ice%heat_flux_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%sw_thru_diag)) then
         call reg%register_2d("ice_sw_thru_diag", state%ice%sw_thru_diag, 0, &
                              size(state%ice%sw_thru_diag, 1), &
                              size(state%ice%sw_thru_diag, 2), optional=.true.)
      end if

      ! --- Ice-shelf cavity basal melt (P2b): the TWO OWNED surface-flux
      !     components.  These ARE registered, and the exclusion rule
      !     above says exactly why: `Q_heat`/`Q_salt` are derived views
      !     and stay out, but "any filler that makes a COMPONENT
      !     time-varying MUST register it".  The melt rate is a function
      !     of the live state, so `heat_cavity`/`salt_cavity` are
      !     time-varying — and they are written at the END of outer step
      !     N and integrated on step N+1 (the ice coupler's documented
      !     one-step lag), so without them a warm restart would apply
      !     zero melt for its first step.  `optional=.true.`: an older
      !     checkpoint resumes with the zero seed rather than failing.
      !     The slot's own arrays (`cavity_flux%melt`, `t_b`, ...) are
      !     NOT registered — they are recomputed from state before first
      !     use, the same derived-field rule as `mass_flux_*`.
      if (allocated(state%surface_flux%heat_cavity)) then
         call reg%register_2d("sf_heat_cavity", state%surface_flux%heat_cavity, 0, &
                              size(state%surface_flux%heat_cavity, 1), &
                              size(state%surface_flux%heat_cavity, 2), optional=.true.)
      end if
      if (allocated(state%surface_flux%salt_cavity)) then
         call reg%register_2d("sf_salt_cavity", state%surface_flux%salt_cavity, 0, &
                              size(state%surface_flux%salt_cavity, 1), &
                              size(state%surface_flux%salt_cavity, 2), optional=.true.)
      end if

      ! --- Assembled net surface fluxes under the component set.  With
      !     `use_components` the assembler derives Q_heat/Q_salt at the END
      !     of each thermo step and the following steps read them, so they
      !     are carried state (an end-of-window checkpoint is read by the
      !     very next step); `sf_q_assembled` says the arrays hold an
      !     assembly (engine_setup then resumes them instead of the
      !     configure seed + ice fold, which sum the same terms in a
      !     different order and miss every non-ice component).  Optional:
      !     an older checkpoint resumes with the seed + fold as before.
      !     `set_components` runs BEFORE the restart read in engine_setup,
      !     so these (and sf_heat_cavity/sf_salt_cavity above) exist when
      !     the read walks the registry. ---
      if (state%surface_flux%use_components) then
         call reg%register_2d("sf_Q_heat", state%surface_flux%Q_heat, 0, &
                              size(state%surface_flux%Q_heat, 1), &
                              size(state%surface_flux%Q_heat, 2), optional=.true.)
         call reg%register_2d("sf_Q_salt", state%surface_flux%Q_salt, 0, &
                              size(state%surface_flux%Q_salt, 1), &
                              size(state%surface_flux%Q_salt, 2), optional=.true.)
         call reg%register_scalar("sf_q_assembled", state%surface_flux%q_assembled, &
                                  optional=.true.)
      end if

      ! --- Sea-ice PR 5: C-grid EVP dynamics prognostics.  u_ice/v_ice
      !     become PROGNOSTIC under dynamics=.true. (previously a v1
      !     interim sampler output — cheap to always register).
      !     str_d/str_t/str_s are the H&D stress tensor components;
      !     fxoc/fyoc are carried for DIAGNOSTIC CONTINUITY only (PR 63 —
      !     see below, they are no longer what makes tau_x/tau_y
      !     restart-exact).  All optional: a resume from a pre-5
      !     checkpoint re-seeds 0 (a one-window cold-start of the ice
      !     velocity/stress) rather than failing.  tau_a_x/tau_a_y are
      !     deliberately NOT registered (configure-time snapshot, rebuilt
      !     fresh from the wind-stress config on every run).
      !
      !     PR 63: tau_ocn_x/tau_ocn_y + tau_ocn_valid are what make the
      !     ice->ocean tau MEDIATION restart-exact.  ice_ocean_stress_flux
      !     mirrors its own blended output into tau_ocn_x/y every outer
      !     step; ice_ocean_stress_resume_apply COPIES them into
      !     surface_stress%tau_x/y at configure (after the wind re-seed),
      !     replacing the old deleted resume-fold RECONSTRUCT routine
      !     (which used the checkpoint's POST-thermo ci — wrong whenever a
      !     checkpoint step's thermo/transport changed ci after the blend,
      !     F4).  fxoc/fyoc above are no longer load-bearing for this: they
      !     are write-only-within-a-call accumulators
      !     (`ice_evp_dynamics` zeros them at entry), so the resume fold
      !     was their only cross-step-boundary reader — that reader is
      !     gone now.  tau_ocn_valid is register_scalar (device_mapped
      !     forced .false.) and gated on the SAME allocated(...) check as
      !     the arrays, so an ocean-only restart file never grows this
      !     variable. ---
      if (allocated(state%ice%u_ice)) then
         call reg%register_2d("ice_u_ice", state%ice%u_ice, 0, &
                              size(state%ice%u_ice, 1), &
                              size(state%ice%u_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%v_ice)) then
         call reg%register_2d("ice_v_ice", state%ice%v_ice, 0, &
                              size(state%ice%v_ice, 1), &
                              size(state%ice%v_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_d)) then
         call reg%register_2d("ice_str_d", state%ice%str_d, 0, &
                              size(state%ice%str_d, 1), &
                              size(state%ice%str_d, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_t)) then
         call reg%register_2d("ice_str_t", state%ice%str_t, 0, &
                              size(state%ice%str_t, 1), &
                              size(state%ice%str_t, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_s)) then
         call reg%register_2d("ice_str_s", state%ice%str_s, 0, &
                              size(state%ice%str_s, 1), &
                              size(state%ice%str_s, 2), optional=.true.)
      end if
      if (allocated(state%ice%fxoc)) then
         call reg%register_2d("ice_fxoc", state%ice%fxoc, 0, &
                              size(state%ice%fxoc, 1), &
                              size(state%ice%fxoc, 2), optional=.true.)
      end if
      if (allocated(state%ice%fyoc)) then
         call reg%register_2d("ice_fyoc", state%ice%fyoc, 0, &
                              size(state%ice%fyoc, 1), &
                              size(state%ice%fyoc, 2), optional=.true.)
      end if
      ! PR 63: the ice->ocean tau mediation seam (see the comment block
      ! above).  The allocated(...) gate is mandatory for all THREE
      ! entries, including the scalar — tau_ocn_valid lives on the type
      ! whether or not the ice slot is live; gating it on the array's
      ! allocation (rather than registering it unconditionally) keeps an
      ! ocean-only restart file free of a stray ice_tau_ocn_valid
      ! variable.
      if (allocated(state%ice%tau_ocn_x)) then
         call reg%register_2d("ice_tau_ocn_x", state%ice%tau_ocn_x, 0, &
                              size(state%ice%tau_ocn_x, 1), &
                              size(state%ice%tau_ocn_x, 2), optional=.true.)
         call reg%register_2d("ice_tau_ocn_y", state%ice%tau_ocn_y, 0, &
                              size(state%ice%tau_ocn_y, 1), &
                              size(state%ice%tau_ocn_y, 2), optional=.true.)
         call reg%register_scalar("ice_tau_ocn_valid", state%ice%tau_ocn_valid, &
                                  optional=.true.)
      end if

      ! --- Sea-ice Winton column prognostics (PR 3a): part_size, m_ice,
      !     m_snow are rank-3 (register directly); enth_ice/sal_ice
      !     (rank-4, per-k bottom-up) and enth_snow (rank-4, nk=1) are
      !     registered as per-k rank-3 slices — the registry is
      !     rank-2/3 only.  All optional: a resume from a pre-3a
      !     checkpoint warns-and-seeds the init values (part_size all
      !     open water, m_ice/m_snow/enth_* = 0, sal_ice =
      !     ICE_BULK_SALINITY) rather than failing. ---
      if (allocated(state%ice%part_size)) then
         call reg%register_3d("ice_part_size", state%ice%part_size, 0, &
                              size(state%ice%part_size, 1), &
                              size(state%ice%part_size, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_ice)) then
         call reg%register_3d("ice_m_ice", state%ice%m_ice, 0, &
                              size(state%ice%m_ice, 1), &
                              size(state%ice%m_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_snow)) then
         call reg%register_3d("ice_m_snow", state%ice%m_snow, 0, &
                              size(state%ice%m_snow, 1), &
                              size(state%ice%m_snow, 2), optional=.true.)
      end if
      if (allocated(state%ice%enth_ice)) then
         do it = 1, size(state%ice%enth_ice, 4)
            call reg%register_3d("ice_enth_ice_k"//to_string(it), &
                                 state%ice%enth_ice(:, :, :, it), 0, &
                                 size(state%ice%enth_ice, 1), &
                                 size(state%ice%enth_ice, 2), optional=.true.)
         end do
      end if
      if (allocated(state%ice%sal_ice)) then
         do it = 1, size(state%ice%sal_ice, 4)
            call reg%register_3d("ice_sal_ice_k"//to_string(it), &
                                 state%ice%sal_ice(:, :, :, it), 0, &
                                 size(state%ice%sal_ice, 1), &
                                 size(state%ice%sal_ice, 2), optional=.true.)
         end do
      end if
      if (allocated(state%ice%enth_snow)) then
         do it = 1, size(state%ice%enth_snow, 4)
            call reg%register_3d("ice_enth_snow_k"//to_string(it), &
                                 state%ice%enth_snow(:, :, :, it), 0, &
                                 size(state%ice%enth_snow, 1), &
                                 size(state%ice%enth_snow, 2), optional=.true.)
         end do
      end if

      ! --- Tidal-mixing + law-of-the-wall BBL kd_int (merged every stage).
      !     Their column solves refresh once per thermo interval, so a
      !     resume between refreshes (always the case under
      !     `diabatic_first = .false.`, where the refresh is the interval's
      !     last step) must carry the interval's diffusivity: a cold zero
      !     would drop the closure for the rest of the interval.
      !     `optional`: a pre-fix checkpoint has no such field. ---
      if (allocated(state%vmix_tidal%kd_int)) then
         call register_full_3d_opt(reg, "tidal_kd_int", state%vmix_tidal%kd_int)
      end if
      if (allocated(state%vdiff%bbl_mix%kd_int)) then
         call register_full_3d_opt(reg, "bbl_kd_int", state%vdiff%bbl_mix%kd_int)
      end if

      ! --- Kappa-shear kd_int / tke_int (merged every stage) ---
      if (allocated(state%kshear%kd_int)) then
         call register_full_3d(reg, "kshear_kd_int", state%kshear%kd_int)
      end if
      if (allocated(state%kshear%tke_int)) then
         call register_full_3d(reg, "kshear_tke_int", state%kshear%tke_int)
      end if

      ! --- Kappa-shear VERTEX corner diffusivity (c069 restart fix):
      !     `kd_corner` is the Pass-B -> Pass-C carrier `kappa_shear_compute`
      !     writes and `vdiff_apply_momentum`'s `kv_corner_source` reads --
      !     but it is only refreshed at STAGE 1 of a thermo step
      !     (`vmix_apply_in_stage`, "stage-invariant by design"), while
      !     `visc_rem_precompute` -- called at the START of every stage,
      !     BEFORE that refresh -- reads it on EVERY stage, including
      !     stage 1 itself (MOM6-order parity, PGF_BUG.md Sec.9:
      !     `vertvisc_coef` runs before `btstep`).  On a continued run
      !     stage 1's `visc_rem_precompute` sees the PREVIOUS step's
      !     converged `kd_corner`, still live in process memory.  On a
      !     warm restart `kd_corner` is a plain workspace -- allocated
      !     zero by `init_vertex`, never written before the first resumed
      !     stage's `visc_rem_precompute` call -- so that FIRST call
      !     silently drops the vertex kappa-shear's corner viscosity from
      !     the remnant matrix, perturbing `bt_visc_rem_u/v` at the
      !     `kappa_trunc` scale (~1e-9) and from there every consumer of
      !     the BT-correction weighting (compat row c069).  `optional
      !     =.true.`: a pre-fix checkpoint has no such field and must
      !     still resume (re-seeding `kd_corner` at 0 for the first
      !     stage only, same as the pre-existing defect). ---
      if (allocated(state%kshear%kd_corner)) then
         call register_full_3d_opt(reg, "kshear_kd_corner", state%kshear%kd_corner)
      end if

      ! --- Live persistent BC state (review #3) ---
      ! Chapman radiation corner scalars: host-side, always present.
      call reg%register_scalar("bc_eta_old_chapman_w", state%bc%eta_old_chapman_w)
      call reg%register_scalar("bc_eta_old_chapman_e", state%bc%eta_old_chapman_e)
      call reg%register_scalar("bc_eta_old_chapman_s", state%bc%eta_old_chapman_s)
      call reg%register_scalar("bc_eta_old_chapman_n", state%bc%eta_old_chapman_n)
      ! Tracer reservoirs: rank-3 (edge, layer, tracer); device-mapped
      ! when allocated, so they join the update-self walk by default.
      if (allocated(state%bc%tres_west)) then
         call register_full_3d(reg, "bc_tres_west", state%bc%tres_west)
      end if
      if (allocated(state%bc%tres_east)) then
         call register_full_3d(reg, "bc_tres_east", state%bc%tres_east)
      end if
      if (allocated(state%bc%tres_south)) then
         call register_full_3d(reg, "bc_tres_south", state%bc%tres_south)
      end if
      if (allocated(state%bc%tres_north)) then
         call register_full_3d(reg, "bc_tres_north", state%bc%tres_north)
      end if
   end subroutine ocean_state_build_restart_registry

   subroutine register_full_3d(reg, tag, arr)
      !! Register a rank-3 field as a FULL local array (interior +
      !! ghosts) — ng=0 over the total extent.  See the ghost-cell
      !! policy note in `ocean_state_build_restart_registry`.
      type(restart_registry_t), intent(inout) :: reg
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: arr(:, :, :)
      call reg%register_3d(tag, arr, 0, size(arr, 1), size(arr, 2))
   end subroutine register_full_3d

   subroutine register_full_3d_opt(reg, tag, arr)
      !! `register_full_3d`, but OPTIONAL on read: a checkpoint written
      !! before the field was registered still resumes.
      type(restart_registry_t), intent(inout) :: reg
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: arr(:, :, :)
      call reg%register_3d(tag, arr, 0, size(arr, 1), size(arr, 2), optional=.true.)
   end subroutine register_full_3d_opt

#ifndef RDB_NO_NETCDF
   module subroutine ocean_state_fill_restart_metadata(state, grid, meta)
      !! Build the grid/vcoord/tracer fingerprint validated on resume
      !! (review #7).
      type(ocean_state_t), intent(in) :: state
      type(hgrid_t), intent(in) :: grid
      type(ocean_restart_metadata_t), intent(out) :: meta
      meta%nz_ml = state%multilayer%nz_ml
      meta%nghost = grid%nghost
      meta%vcoord_type = state%vcoord%coord_type
      meta%vcoord_name = vcoord_type_name(state%vcoord%coord_type)
      meta%idx_salinity = state%multilayer%idx_salinity
      meta%idx_temperature = state%multilayer%idx_temperature
      if (allocated(state%multilayer%tracers)) then
         meta%n_tracers = size(state%multilayer%tracers)
      else
         meta%n_tracers = 0
      end if
   end subroutine ocean_state_fill_restart_metadata
#endif

   pure function vcoord_type_name(code) result(name)
      !! Human-readable tag for a VCOORD_* enum (mismatch messages only;
      !! the integer enum is what is validated).
      use rdb_constants, only: VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                               VCOORD_ZSTAR_SIGMA, VCOORD_ZSTAR_FULL, &
                               VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN
      integer, intent(in) :: code
      character(len=32) :: name
      select case (code)
      case (VCOORD_SIGMA)
         name = "sigma"
      case (VCOORD_ZSIGMA)
         name = "zsigma"
      case (VCOORD_ZSTAR)
         name = "zstar"
      case (VCOORD_ZSTAR_SIGMA)
         name = "zstar_sigma"
      case (VCOORD_ZSTAR_FULL)
         name = "zstar_full"
      case (VCOORD_EULERIAN_Z)
         name = "eulerian_z"
      case (VCOORD_LAGRANGIAN)
         name = "lagrangian"
      case default
         name = "code_"//to_string(code)
      end select
   end function vcoord_type_name

   module subroutine ocean_state_restart_read(state, grid, decomp, filename, t, step, ierr)
      !! Read a per-rank ocean restart into the (host) prognostic arrays.
      !! MUST run BEFORE `ocean_state_enter_data` — the subsequent H->D
      !! copy carries the restored values to the device.  Restores
      !! `dyn%outer_step_count` (dt_therm alignment).
      !!
      !! Error policy (review #4): when the caller PASSES `ierr`, a
      !! decomp/schema/grid mismatch returns it non-zero (the test path).
      !! When the caller OMITS `ierr` (the production driver path), a
      !! mismatch is FATAL — `error stop` naming the kind + the file —
      !! so a mismatched resume dies loudly instead of silently
      !! cold-starting and then clobbering the good checkpoint at the next
      !! cadence.
      type(ocean_state_t), intent(inout), target :: state
      type(hgrid_t), intent(in) :: grid
      type(decomp_t), intent(in) :: decomp
      character(len=*), intent(in) :: filename
      real(wp), intent(out) :: t
      integer, intent(out) :: step
      integer, intent(out), optional :: ierr
#ifndef RDB_NO_NETCDF
      type(restart_registry_t) :: reg
      type(ocean_restart_metadata_t) :: meta
      integer :: osc, local_ierr

      call ocean_state_build_restart_registry(state, grid, reg)
      call ocean_state_fill_restart_metadata(state, grid, meta)
      call ocean_restart_read_local(filename, reg, decomp, meta, t, step, osc, &
                                    ierr=local_ierr)
      if (present(ierr)) then
         ! Translate the INTERNAL 1/2/3 restart_mismatch_kind convention
         ! (private to `ocean_restart_check_decomp`/`ocean_restart_read_local`)
         ! into a named `OCEAN_STATUS_ERR_RESTART_*` code before it reaches
         ! the caller — the raw integers collide with
         ! `OCEAN_STATUS_ERR_CONFIG_PARSE`/`_CONFIG_VALIDATE`/`_SETUP`
         ! (P0.1 review F4).
         if (local_ierr /= 0) then
            ierr = restart_mismatch_status_code(local_ierr)
         else
            ierr = OCEAN_STATUS_OK
         end if
      else if (local_ierr /= 0) then
         call fail("Ocean restart resume aborted: "// &
                   restart_mismatch_kind(local_ierr)//" mismatch reading "// &
                   trim(filename)//".  A mismatched resume would cold-start "// &
                   "and then overwrite the good checkpoint; refusing.", &
                   code=restart_mismatch_status_code(local_ierr))
      end if
      if (local_ierr == 0) then
         state%dyn%outer_step_count = osc
         ! PR-2 (bt-rem-from-av-rem review): record whether THIS read
         ! actually found vmix_kv (vs. an older checkpoint without it, or
         ! this call never running at all on a mismatch) -- read by
         ! `configure_ocean_lateral`'s config-copy to decide whether to
         ! skip reseeding `kv` (see `ocean_vmix_t%kv_from_restart`).
         state%vmix%kv_from_restart = reg%entry_found("vmix_kv")
      end if
#else
      t = 0.0_wp
      step = 0
      call fail("ocean restart requires RDB_ENABLE_NETCDF=ON", ierr, OCEAN_STATUS_ERR_IO)
      return
#endif
   end subroutine ocean_state_restart_read

   pure function restart_mismatch_kind(code) result(kind)
      !! Map the check-decomp ierr code (INTERNAL 1/2/3 convention —
      !! see `ocean_restart_check_decomp`) to a human-readable kind for
      !! the legacy `error stop` message.
      integer, intent(in) :: code
      character(len=:), allocatable :: kind
      select case (code)
      case (1)
         kind = "schema-version"
      case (2)
         kind = "decomposition"
      case (3)
         kind = "grid/vcoord/tracer"
      case default
         kind = "unknown"
      end select
   end function restart_mismatch_kind

   pure function restart_mismatch_status_code(code) result(status_code)
      !! Map the check-decomp ierr code (INTERNAL 1/2/3 convention) onto
      !! the PUBLIC, collision-free `OCEAN_STATUS_ERR_RESTART_*` codes
      !! returned through `ocean_state_restart_read`'s `ierr` (P0.1
      !! review F4).
      integer, intent(in) :: code
      integer :: status_code
      select case (code)
      case (1)
         status_code = OCEAN_STATUS_ERR_RESTART_SCHEMA
      case (2)
         status_code = OCEAN_STATUS_ERR_RESTART_DECOMP
      case (3)
         status_code = OCEAN_STATUS_ERR_RESTART_GRID
      case default
         status_code = OCEAN_STATUS_ERR_RESTART_GRID
      end select
   end function restart_mismatch_status_code

end submodule rdb_ocean_state_restart
