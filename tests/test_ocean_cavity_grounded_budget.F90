!! THE GROUNDED-CAVITY BUDGET GATE: a land column holds nothing, and it
!! holds nothing from step 0.
module test_ocean_cavity_grounded_budget
   !! ### What a budget is, and what this suite refuses to accept
   !!
   !! The console prints, per conserved quantity, a RELATIVE residual
   !!
   !! ```
   !!   Error = [ (total(t) - total(0)) + out - src ] / total(0)
   !! ```
   !!
   !! and the statement this suite gates is the maintainer's: **that
   !! residual must sit at round-off, at EVERY step, starting at step 0.
   !! When the domain total changes, the change must be a TRACKED SOURCE
   !! — never an unexplained residual.**  A jump between step 0 and step
   !! 1 is not "the initial transient"; it is content that entered or
   !! left the total without passing through a term the budget names.
   !!
   !! So the bound here is `1e-12` RELATIVE, asserted at every one of the
   !! 20 steps INCLUDING step 0, and it is never relaxed to accommodate a
   !! result.  Re-latching the baseline after the first step, subtracting
   !! the jump, or dropping cells from the integral after the fact would
   !! each make this suite pass while the defect it was written for was
   !! still there.
   !!
   !! ### The defect this suite was written for
   !!
   !! An ice-shelf column that GROUNDS (`b - z_draft < &ocean_cavity_dyn_nml
   !! h_min_cavity`) becomes LAND through `seed_wet_mask_impl`, exactly as
   !! an ordinary shallow-bathymetry column does.  Its water column
   !! `b - z_draft` is NEGATIVE — by hundreds of metres for a real draft
   !! over a real bed — so its seeded `h_layer` is negative too.
   !!
   !! `ocean_state_seed_land_cells` then floored `h_layer` to
   !! `H_VANISHED` (correct: the mass total closed) but tried to hold the
   !! tracer at its seeded concentration with
   !! `val = hTr/max(h_old, H_VANISHED)`, then `hTr = val*H_VANISHED`.
   !! For `h_old <= H_VANISHED` that pair is an exact ALGEBRAIC IDENTITY:
   !! `hTr` came out of it UNCHANGED, i.e. a full-column `S*(b - z_draft)`
   !! sitting next to a floored `h` — an implied concentration of order
   !! `-1e7` PSU.  The budget latch integrated it; the first ALE regrid
   !! then discarded it (a land layer is pinned AT the `H_VANISHED`
   !! vanish marker, so the remap's `h > H_FLOOR` concentration gate is
   !! false and it writes `hTr = 0`), and the whole thing showed up as a
   !! step change in `Error` between step 0 and step 1.
   !!
   !! On `validation_examples/ocean/isomip_plus/ocean0_idealised_draft.nml`
   !! (39.4 % of interior columns grounded) that step change was **0.609
   !! of the initial salt content** and `-0.080` of the heat, while mass
   !! closed at `-5.0e-14`.  It scaled with the grounded columns' summed
   !! DEPTH DEFICIT `sum(z_draft - b)`, not with their AREA, which is why
   !! two runs at 39.4 % and 14 % grounded gave 6.09e-01 and 1.87e-02
   !! rather than anything proportional to the fraction.
   !!
   !! ### The land-state contract this suite pins
   !!
   !! ```
   !!   a LAND T-cell holds   h_layer = H_VANISHED   and   hTr = 0
   !! ```
   !!
   !! — for every land cell, however it became land, at `t = 0` and at
   !! every step after.  `h` is pinned AT the D4 vanish marker, so a land
   !! layer is on the vanished side of every `h > H_VANISHED` gate in the
   !! tree; `hTr = 0` is the content those gates already give it.  The
   !! seed's job is to hand the budget latch that same state, so
   !! `grounded_holds_the_ordinary_land_state` asserts it BITWISE against
   !! an ordinary-land column in the same run, before and after 20 steps.
   !!
   !! ### `mem:separate` discipline
   !!
   !! The engine maps the whole state with `engine_enter_data`, and every
   !! budget term here is formed by the PRODUCTION reducers
   !! (`compute_total_h` / `compute_total_tracer`, which carry their own
   !! `!$acc parallel loop ... present(...)`) rather than by a host loop,
   !! so the totals are correct on the GPU build without a single extra
   !! transfer.  The three places that genuinely need HOST data — the
   !! land/wet cell-class comparisons — pull it with an explicit
   !! `!$acc update self` of COMPONENT arrays, never of an aggregate
   !! derived type, and the tracer registry's array-of-derived-type
   !! indirection is dereferenced through an `associate` first (the
   !! `window_totals` pattern in `test_ocean_conservation_salt_heat`).
   !! All of it is inert on the host build.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_console_stats, only: compute_total_h, compute_total_tracer, &
                                      ocean_budget_src, ocean_budget_out, &
                                      ocean_budget_stage_weight, &
                                      ocean_salt_src_sum, ocean_heat_src_sum
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_cavity_grounded_budget_tests

   integer, parameter :: N_STEPS = 20
      !! Steps of the FULL split solver every budget gate runs for.  The
      !! residual is asserted after each one AND at step 0.
   real(wp), parameter :: DT = 300.0_wp
   real(wp), parameter :: BUDGET_TOL = 1.0e-12_wp
      !! The gate.  A relative residual, so it is dimensionless and the
      !! same number serves mass, salt and heat.  Not a tuned tolerance:
      !! a closed budget in double precision drifts at `~n_steps * eps *
      !! cancellation`, which for these domains is 1e-15..1e-13.  If this
      !! ever needs raising, something stopped being tracked.
   real(wp), parameter :: BED = 400.0_wp
      !! Flat bed (m).  Shallow on purpose: it makes the grounded
      !! columns' depth deficit a LARGE fraction of the wet content, so
      !! the pre-fix defect is a 33 % residual rather than a subtle one.
   real(wp), parameter :: DRAFT = 600.0_wp
      !! Uniform ice draft (m) inside the shelf box.  `BED - DRAFT` is
      !! `-200 m`: the shelf box grounds, and every column in it goes
      !! land with a NEGATIVE seeded thickness — the case the old hold
      !! turned into an identity.
   real(wp), parameter :: MELT_BED = 800.0_wp
      !! Flat bed (m) for the MELT variant only — see `melt_topo`.
   real(wp), parameter :: DX = 2000.0_wp

   ! Shelf box, chosen to coincide CELL-FOR-CELL with the land square
   ! `set_bathymetry_island` carves at `slope_scale = 0.25`: that setter
   ! lands physical indices `nx/2-4+1 .. nx/2+4` (5..12 for NX = 16), and
   ! `set_draft_flat` compares CELL CENTRES `(i_phys - 0.5)*dx`, so
   ! `[4*dx, 12*dx]` covers centres `4.5*dx .. 11.5*dx` — indices 5..12,
   ! and nothing else.  The two runs in `grounded_matches_plain_land`
   ! therefore make the SAME 64 cells land by two different routes.
   real(wp), parameter :: BOX_LO = 4.0_wp*DX
   real(wp), parameter :: BOX_HI = 12.0_wp*DX
   real(wp), parameter :: ISLAND_HALF_FRAC = 0.25_wp

contains

   subroutine collect_ocean_cavity_grounded_budget_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("grounded_cavity_budget_closes_from_step_zero", &
                               test_budget_closes_melt_off), &
                  new_unittest("grounded_cavity_budget_closes_with_melt", &
                               test_budget_closes_melt_on), &
                  new_unittest("grounded_holds_the_ordinary_land_state", &
                               test_grounded_is_ordinary_land), &
                  new_unittest("grounded_matches_plain_land_in_every_wet_column", &
                               test_grounded_matches_plain_land) &
                  ]
   end subroutine collect_ocean_cavity_grounded_budget_tests

   ! ------------------------------------------------------------------
   ! Namelists
   ! ------------------------------------------------------------------

   function base_nml(topo, tfreeze) result(nml)
      !! Everything the configurations share.  No wind, no sponge, walls
      !! all round: the only things that may move the salt and heat
      !! totals are the ones the budget names.
      character(len=*), intent(in) :: topo
         !! The whole `&ocean_topo_nml` group — `flat` for the cavity
         !! runs, `island` for the cavity-free twin.  Passed as a group
         !! rather than appended, because a namelist may not carry the
         !! same group twice.
      character(len=*), intent(in) :: tfreeze
         !! `&ocean_eos_nml tfreeze_set` — `isomip` is REQUIRED by the
         !! melt package (the two shipped liquidi differ by enough to
         !! flip the sign of the melt rate), `seaice` elsewhere.
      character(len=:), allocatable :: nml
      character(len=32) :: dx_s

      write (dx_s, '(F12.2)') DX
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 16, ny = 16, nghost = 2, dx = "// &
            trim(adjustl(dx_s))//", dy = "//trim(adjustl(dx_s))//" /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0 /"//new_line("a")// &
            "&physics_nml coriolis_f = -1.409e-4, wind_stress_x = 0.0, "// &
            "wind_stress_y = 0.0 /"//new_line("a")// &
            topo// &
            ! Stratified, so the wet columns carry a live baroclinic
            ! pressure gradient and the run is not a vacuous rest state.
            !
            ! T is kept STRICTLY POSITIVE (1..3 degC) on purpose.  The
            ! heat total is `sum h*T*area*rho`, so a profile straddling
            ! 0 degC would make the reference `total(0)` a near-cancelled
            ! sum and turn a RELATIVE residual into a measure of that
            ! cancellation rather than of conservation.  A budget gate
            ! must not be reading its own denominator's round-off.
            ! (It is also comfortably above the in-situ freezing point at
            ! these depths, so the melt variant melts rather than freezes.)
            "&tracer_nml initial_salinity = 34.5, T_init_bottom = 1.0, "// &
            "T_init_surface = 3.0 /"//new_line("a")// &
            ! fv_mom6 + p_top_in_bc: REQUIRED for a varying draft, and
            ! carried by the cavity-free twin too so the two runs differ
            ! in the cavity group alone.
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//new_line("a")// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = '"//tfreeze//"' /"// &
            new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', "// &
            "north = 'wall' /"//new_line("a")// &
            ! n_inner pinned: `auto_n_inner` derives it from max(b), which
            ! is the same 400 m in all these runs, but pinning it keeps
            ! the comparison independent of that derivation.
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .false., "// &
            "n_inner = 12 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function base_nml

   function flat_topo() result(grp)
      character(len=:), allocatable :: grp
      character(len=32) :: bed_s
      write (bed_s, '(F12.2)') BED
      grp = "&ocean_topo_nml topo_config = 'flat', max_depth = "// &
            trim(adjustl(bed_s))//", wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//new_line("a")
   end function flat_topo

   function island_topo() result(grp)
      !! The cavity-free twin's bed: the SAME 64 cells made land by
      !! ordinary shallow bathymetry (`b = 0`) instead of by grounding.
      character(len=:), allocatable :: grp
      character(len=32) :: bed_s, hf_s
      write (bed_s, '(F12.2)') BED
      write (hf_s, '(F8.4)') ISLAND_HALF_FRAC
      grp = "&ocean_topo_nml topo_config = 'island', max_depth = "// &
            trim(adjustl(bed_s))//", slope_scale = "//trim(adjustl(hf_s))// &
            ", wind_config = 'constant', taux_magnitude = 0.0 /"//new_line("a")
   end function island_topo

   function cavity_grounded_box() result(grp)
      !! The fully grounded shelf: a uniform `DRAFT` over the central
      !! box, zero outside it.  `BED - DRAFT = -200 m` inside, so every
      !! column of the box grounds with a NEGATIVE seeded thickness —
      !! the case the old tracer hold turned into an identity.
      character(len=:), allocatable :: grp
      character(len=32) :: d_s, lo_s, hi_s

      write (d_s, '(F12.2)') DRAFT
      write (lo_s, '(F16.2)') BOX_LO
      write (hi_s, '(F16.2)') BOX_HI
      grp = "&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', "// &
            "draft_depth = "//trim(adjustl(d_s))// &
            ", draft_x0 = "//trim(adjustl(lo_s))// &
            ", draft_x1 = "//trim(adjustl(hi_s))// &
            ", draft_y0 = "//trim(adjustl(lo_s))// &
            ", draft_y1 = "//trim(adjustl(hi_s))// &
            ", h_min_cavity = 20.0 /"//new_line("a")
   end function cavity_grounded_box

   function melt_topo() result(grp)
      !! The melt case's own bed: deeper (`MELT_BED`) than the grounded-box
      !! runs', because a shelf that both GROUNDS and floats over a real
      !! water column needs vertical room for both, and squeezing them
      !! into a 400 m basin only buys a draft slope steep enough to blow
      !! the barotropic loop up — a numerically violent geometry proves
      !! nothing about a budget.
      character(len=:), allocatable :: grp
      character(len=32) :: bed_s
      write (bed_s, '(F12.2)') MELT_BED
      grp = "&ocean_topo_nml topo_config = 'flat', max_depth = "// &
            trim(adjustl(bed_s))//", wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//new_line("a")
   end function melt_topo

   function cavity_melt_shelf() result(grp)
      !! A shelf with BOTH kinds of column under it, over the `MELT_BED`
      !! basin: a linear draft of 800 m at `x = 0` thinning eastward at
      !! `8e-3 m/m`.  With `h_min_cavity = 100 m` a column GROUNDS while
      !! the draft exceeds `MELT_BED - 100 = 700 m`, so the western six
      !! columns are land and the rest of the basin floats under ice over
      !! 100-250 m of water.
      !!
      !! There is deliberately NO CALVING FRONT inside the domain: the
      !! shelf box is left open on every side, so the lid runs to the
      !! east wall.  A front would put a several-hundred-metre STEP in
      !! the ice load across one face, which at this resolution is a
      !! barotropic shock, not a cavity — the run blows up in two steps
      !! and the budget it would have measured never happens.  A front
      !! is a real configuration and it is tested elsewhere
      !! (`test_ocean_cavity_draft`, and the ISOMIP+ case, where the
      !! draft has thinned to 90 m by the time it reaches one); this
      !! gate is about conservation and buys nothing by carrying it.
      !!
      !! Both kinds must be present for the melt gate to mean anything:
      !! the grounded columns are the defect's habitat, and the covered
      !! WET columns are the only place the basal-melt kernel delivers a
      !! source at all.  A fully grounded box would give a melt run with
      !! zero melt and a gate that passes on nothing.
      !!
      !! `h_min_cavity` is 100 m rather than the 20 m of the other cases
      !! for a stated reason: it puts the grounding line where the water
      !! column is still thick enough to carry four honest layers, rather
      !! than leaving a ring of 20 m columns under 700 m of ice at the
      !! grounding line — which is a thin-layer stress test, not a
      !! conservation one.
      character(len=:), allocatable :: grp

      grp = "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 800.0, draft_slope = -8.0e-3, draft_x0 = 0.0"// &
            ", h_min_cavity = 100.0 /"//new_line("a")// &
            ! The basal-melt package: a heat flux and a VIRTUAL salt
            ! flux, both delivered through the ordinary
            ! `heat_budget_surface` / `salt_budget_surface` contributors,
            ! so the `src` term the console prints is the only place they
            ! can appear.
            "&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0 /"// &
            new_line("a")// &
            "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"// &
            new_line("a")// &
            ! A real basin has lateral friction; a sloping-lid cavity run
            ! with none is a grid-noise experiment.  Stated rather than
            ! tuned: 100 m2/s at dx = 2 km is a viscous CFL of
            ! nu*dt/dx^2 = 7.5e-3, four decades inside the limit.
            "&ocean_hvisc_nml nu_h = 100.0 /"//new_line("a")
   end function cavity_melt_shelf

   ! ------------------------------------------------------------------
   ! Budget bookkeeping — the console's own terms, through its own helpers
   ! ------------------------------------------------------------------

   subroutine totals(engine, total_mass, total_salt, total_heat, &
                     src_salt, src_heat, out_salt, out_heat, mass_out)
      !! Form exactly what `ocean_console_stats_report` forms: the
      !! area-weighted domain totals and the `src` / `out` budget terms,
      !! through the PRODUCTION reducers and the production
      !! `ocean_budget_*` helpers.  Asserting on a re-derivation would
      !! gate a second implementation; this gates the one that prints.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: total_mass, total_salt, total_heat
      real(wp), intent(out) :: src_salt, src_heat, out_salt, out_heat
      real(wp), intent(out) :: mass_out

      real(wp) :: total_h, bud_w

      associate (ms => engine%state%multilayer, mt => engine%state%metrics, &
                 ng => engine%grid%nghost)
         bud_w = ocean_budget_stage_weight(engine%state%dyn%split_scheme == &
                                           SPLIT_SCHEME_PRED_CORR)

         total_h = compute_total_h(ms%h_layer, mt%areaT, ng)
         total_mass = total_h*RHO_WATER
         total_salt = compute_total_tracer(ms%tracers(ms%idx_salinity)%hTr, &
                                           mt%areaT, ng)*RHO_WATER
         total_heat = compute_total_tracer(ms%tracers(ms%idx_temperature)%hTr, &
                                           mt%areaT, ng)*RHO_WATER

         src_salt = ocean_budget_src( &
                    ocean_salt_src_sum( &
                    compute_total_tracer(ms%salt_budget_surface, mt%areaT, ng), &
                    compute_total_tracer(ms%salt_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_salt = ocean_budget_out( &
                    compute_total_tracer(ms%salt_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_tracer(ms%salt_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)

         src_heat = ocean_budget_src( &
                    ocean_heat_src_sum( &
                    compute_total_tracer(ms%heat_budget_surface, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_geothermal, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_heat = ocean_budget_out( &
                    compute_total_tracer(ms%heat_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)

         mass_out = ms%mass_out
      end associate
   end subroutine totals

   pure function residual(total, ref, out_term, src_term) result(rel)
      !! `rdb_console_stats::emit_drift_line`'s residual, relative: what
      !! the `Error` column prints.
      real(wp), intent(in) :: total, ref, out_term, src_term
      real(wp) :: rel
      if (abs(ref) > 1.0e-30_wp) then
         rel = ((total - ref) + out_term - src_term)/ref
      else
         rel = 0.0_wp
      end if
   end function residual

   ! ------------------------------------------------------------------
   ! Engine driving
   ! ------------------------------------------------------------------

   subroutine start_engine(nml, engine, cfg, ok)
      !! Config -> validate -> setup -> map.  One place, so the three
      !! budget gates cannot drift apart in how they start a run.
      character(len=*), intent(in) :: nml
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      integer :: ierr

      ok = .false.
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= 0) return
      call validate_config(cfg, ierr)
      if (ierr /= 0) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine start_engine

   subroutine advance(engine, cfg, t, ierr)
      !! One outer step in `driver_run_ocean`'s MANDATED ORDER:
      !! dyn-core, then the sea-ice block, then the surface-flux
      !! assembly + diag step.  The melt package writes its two owned
      !! flux components inside `engine_step_finalize`, so a gate that
      !! skipped it would never see the melt source at all.
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t
      integer, intent(out) :: ierr

      call engine_step(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_finalize(engine, DT, t, ierr=ierr)
   end subroutine advance

   subroutine stop_engine(engine)
      type(ocean_engine_t), intent(inout) :: engine
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine stop_engine

   ! ------------------------------------------------------------------
   ! Cell-class snapshots
   ! ------------------------------------------------------------------

   subroutine pull_state(engine, h, hs, ht, wet)
      !! Host copies of the four fields the cell-class comparisons read.
      !! COMPONENT arrays only, and the registry indirection is resolved
      !! through `associate` before the `update self` names it.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), allocatable, intent(out) :: h(:, :, :), hs(:, :, :), ht(:, :, :)
      real(wp), allocatable, intent(out) :: wet(:, :)

      associate (ms => engine%state%multilayer)
         associate (hl => ms%h_layer, wm => ms%wet_mask, &
                    hsal => ms%tracers(ms%idx_salinity)%hTr, &
                    htmp => ms%tracers(ms%idx_temperature)%hTr)
            !$acc update self(hl, wm, hsal, htmp) if_present
            allocate (h, source=hl)
            allocate (hs, source=hsal)
            allocate (ht, source=htmp)
            allocate (wet, source=wm)
         end associate
      end associate
   end subroutine pull_state

   ! ------------------------------------------------------------------
   ! The gates
   ! ------------------------------------------------------------------

   subroutine budget_run(error, nml, label, live_source)
      !! Step a partly grounded cavity for `N_STEPS` and assert the three
      !! relative residuals at EVERY step, step 0 included.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml
      character(len=*), intent(in) :: label
      logical, intent(in) :: live_source
         !! `.true.` when the configuration is supposed to DELIVER a
         !! tracked source (the melt variant).  A closed budget with an
         !! identically-zero source is the same assertion as "nothing
         !! happened", so the melt gate also demands that the source and
         !! the change in the total are both non-zero — otherwise it
         !! would pass on a melt package that silently did nothing.

      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr
      real(wp) :: t
      real(wp) :: m0, s0, h0
      real(wp) :: mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out
      real(wp) :: r_m, r_s, r_h, worst_m, worst_s, worst_h

      call start_engine(nml, engine, cfg, ok)
      call check(error, ok, label//": engine starts")
      if (allocated(error)) return

      gate: block
         ! --- step 0: the latch.  The reference the residual is measured
         ! against is THIS state, so a defect that only corrupts the
         ! latch reports as a jump at step 1 and never afterwards.
         call totals(engine, m0, s0, h0, src_s, src_h, out_s, out_h, m_out)
         call check(error, ieee_is_finite(s0) .and. ieee_is_finite(h0), &
                    label//": step-0 totals are finite")
         if (allocated(error)) exit gate
         call check(error, abs(s0) > 0.0_wp, label//": step-0 salt total is non-zero")
         if (allocated(error)) exit gate

         worst_m = 0.0_wp
         worst_s = 0.0_wp
         worst_h = 0.0_wp
         t = 0.0_wp
         do n = 1, N_STEPS
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": step succeeds")
            if (allocated(error)) exit gate
            t = t + DT

            call totals(engine, mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out)
            r_m = residual(mt_, m0, m_out*RHO_WATER, 0.0_wp)
            r_s = residual(st_, s0, out_s, src_s)
            r_h = residual(ht_, h0, out_h, src_h)
            worst_m = max(worst_m, abs(r_m))
            worst_s = max(worst_s, abs(r_s))
            worst_h = max(worst_h, abs(r_h))
         end do

         ! The gate.  `worst_*` is the maximum over all 20 steps, and the
         ! step-0 residual is identically zero by construction (the latch
         ! is the reference), so asserting the maximum asserts every step.
         call check(error, worst_m <= BUDGET_TOL, &
                    label//": mass residual stays at round-off over "// &
                    "20 steps of the full split solver")
         if (allocated(error)) exit gate
         call check(error, worst_s <= BUDGET_TOL, &
                    label//": salt residual stays at round-off — the change "// &
                    "in the salt total equals the tracked source, every step")
         if (allocated(error)) exit gate
         call check(error, worst_h <= BUDGET_TOL, &
                    label//": heat residual stays at round-off — the change "// &
                    "in the heat total equals the tracked source, every step")
         if (allocated(error)) exit gate

         ! Non-vacuity.  With melt on, `src` must actually be a number
         ! and the totals must actually have moved by it; a budget that
         ! closes because nothing was delivered proves nothing about a
         ! budget that has to account for something.
         if (live_source) then
            call check(error, abs(src_s) > 0.0_wp .and. abs(src_h) > 0.0_wp, &
                       label//": the melt package delivers a non-zero tracked "// &
                       "salt AND heat source (else the closure is vacuous)")
            if (allocated(error)) exit gate
            call check(error, abs(st_ - s0) > 0.0_wp .and. abs(ht_ - h0) > 0.0_wp, &
                       label//": and the domain totals moved — the change in "// &
                       "each total IS that accumulated source, to 1e-12")
         end if
      end block gate

      call stop_engine(engine)
   end subroutine budget_run

   subroutine test_budget_closes_melt_off(error)
      !! Melt OFF: no tracked source at all, so the salt and heat totals
      !! must be CONSTANT to round-off across 20 steps.  This is the
      !! variant the reported defect showed up on — pre-fix the salt
      !! residual is 3.3e-01 from step 1 and flat thereafter.
      type(error_type), allocatable, intent(out) :: error
      call budget_run(error, base_nml(flat_topo(), "seaice")//cavity_grounded_box(), &
                      label="grounded/melt-off", live_source=.false.)
   end subroutine test_budget_closes_melt_off

   subroutine test_budget_closes_melt_on(error)
      !! Melt ON: the basal-melt package delivers a heat flux and a
      !! VIRTUAL salt flux through `surface_flux%heat_cavity` /
      !! `salt_cavity`, which ride the ordinary
      !! `heat_budget_surface` / `salt_budget_surface` contributors.  The
      !! assertion is the same one — `(total(t) - total(0)) + out - src`
      !! at round-off — which for a live source says the change in the
      !! total EQUALS the accumulated tracked source, not merely that the
      !! total is unchanged.
      type(error_type), allocatable, intent(out) :: error
      call budget_run(error, base_nml(melt_topo(), "isomip")//cavity_melt_shelf(), &
                      label="grounded/melt-on", live_source=.true.)
   end subroutine test_budget_closes_melt_on

   subroutine test_grounded_is_ordinary_land(error)
      !! The land-state contract, BITWISE, at `t = 0` and after 20 steps.
      !!
      !! This run has both kinds of land in one domain: the 64 GROUNDED
      !! columns of the shelf box (seeded from a NEGATIVE `b - z_draft`)
      !! and — nothing else, since the bed is flat and wet everywhere
      !! else.  So the reference it is compared against is the CONTRACT
      !! itself: `h == H_VANISHED` exactly and `hTr == 0` exactly on
      !! every land cell.  `grounded_matches_plain_land` then closes the
      !! loop by making the same cells land the ordinary way and finding
      !! the same numbers.
      type(error_type), allocatable, intent(out) :: error

      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr, i, j, k, n_land
      real(wp) :: t
      real(wp), allocatable :: h(:, :, :), hs(:, :, :), ht(:, :, :), wet(:, :)
      logical :: h_ok, tr_ok

      call start_engine(base_nml(flat_topo(), "seaice")//cavity_grounded_box(), &
                        engine, cfg, ok)
      call check(error, ok, "land-state: engine starts")
      if (allocated(error)) return

      gate: block
         call pull_state(engine, h, hs, ht, wet)
         n_land = 0
         h_ok = .true.
         tr_ok = .true.
         do j = 1, size(wet, 2)
            do i = 1, size(wet, 1)
               if (wet(i, j) /= 0.0_wp) cycle
               n_land = n_land + 1
               do k = 1, size(h, 3)
                  if (h(i, j, k) /= H_VANISHED) h_ok = .false.
                  if (hs(i, j, k) /= 0.0_wp) tr_ok = .false.
                  if (ht(i, j, k) /= 0.0_wp) tr_ok = .false.
               end do
            end do
         end do
         call check(error, n_land > 0, "land-state: the shelf box actually grounds")
         if (allocated(error)) exit gate
         call check(error, h_ok, "land-state t=0: every land layer is exactly H_VANISHED")
         if (allocated(error)) exit gate
         call check(error, tr_ok, "land-state t=0: every land layer holds exactly "// &
                    "zero salt and heat content (the seed hands the budget latch "// &
                    "the state every vanished-gated operator holds land at)")
         if (allocated(error)) exit gate

         t = 0.0_wp
         do n = 1, N_STEPS
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "land-state: step succeeds")
            if (allocated(error)) exit gate
            t = t + DT
         end do

         deallocate (h, hs, ht, wet)
         call pull_state(engine, h, hs, ht, wet)
         h_ok = .true.
         tr_ok = .true.
         do j = 1, size(wet, 2)
            do i = 1, size(wet, 1)
               if (wet(i, j) /= 0.0_wp) cycle
               do k = 1, size(h, 3)
                  if (h(i, j, k) /= H_VANISHED) h_ok = .false.
                  if (hs(i, j, k) /= 0.0_wp) tr_ok = .false.
                  if (ht(i, j, k) /= 0.0_wp) tr_ok = .false.
               end do
            end do
         end do
         call check(error, h_ok, "land-state after 20 steps: land thickness unmoved")
         if (allocated(error)) exit gate
         call check(error, tr_ok, "land-state after 20 steps: land content still "// &
                    "exactly zero — the seed state IS the running state")
      end block gate

      call stop_engine(engine)
   end subroutine test_grounded_is_ordinary_land

   subroutine test_grounded_matches_plain_land(error)
      !! The same 64 cells made land two ways — GROUNDED under a 600 m
      !! draft over a 400 m bed, and ordinary `b = 0` island bathymetry —
      !! must give the same ocean.
      !!
      !! The two runs differ on the land block in `bt_H_ref` (`-200 m`
      !! against `0 m`) and in `p_top` (`rho_ref*g*600` against `0`), and
      !! they agree everywhere else.  If land is genuinely inert — every
      !! face metric zeroed, every land layer on the vanished side of
      !! every gate — those differences reach no wet cell at all and the
      !! agreement is BITWISE.  That is the assertion, and it is the
      !! strongest available: a tolerance here would hide exactly the
      !! leak the test exists to detect.
      !!
      !! Non-vacuity is asserted separately: the stratified column under
      !! the flat lid is not at rest, so `max|hS|` over the wet cells
      !! moves by many decades more than the bound.
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: hA(:, :, :), hsA(:, :, :), htA(:, :, :), wetA(:, :)
      real(wp), allocatable :: hB(:, :, :), hsB(:, :, :), htB(:, :, :), wetB(:, :)
      real(wp) :: d_h, d_s, d_t, sig
      logical :: ok, mask_same

      call run_and_snapshot(base_nml(flat_topo(), "seaice")//cavity_grounded_box(), &
                            hA, hsA, htA, wetA, ok)
      call check(error, ok, "twin A (grounded cavity) runs")
      if (allocated(error)) return
      call run_and_snapshot(base_nml(island_topo(), "seaice"), hB, hsB, htB, wetB, ok)
      call check(error, ok, "twin B (ordinary island land) runs")
      if (allocated(error)) return

      mask_same = all(wetA == wetB)
      call check(error, mask_same, "the two routes make the SAME cells land "// &
                 "(the comparison is otherwise meaningless)")
      if (allocated(error)) return
      call check(error, any(wetA == 0.0_wp) .and. any(wetA /= 0.0_wp), &
                 "the domain really has both land and ocean in it")
      if (allocated(error)) return

      call wet_diff(hA, hB, wetA, d_h)
      call wet_diff(hsA, hsB, wetA, d_s)
      call wet_diff(htA, htB, wetA, d_t)
      sig = wet_scale(hsA, wetA)

      call check(error, d_h == 0.0_wp, "wet-column thickness is BITWISE equal "// &
                 "between grounded-via-cavity and the same cells as plain land")
      if (allocated(error)) return
      call check(error, d_s == 0.0_wp, "wet-column salt content is BITWISE equal — "// &
                 "no land state, however extreme, reaches a wet cell")
      if (allocated(error)) return
      call check(error, d_t == 0.0_wp, "wet-column heat content is BITWISE equal")
      if (allocated(error)) return
      call check(error, sig > 1.0_wp, "non-vacuity: the wet columns carry real content")
   end subroutine test_grounded_matches_plain_land

   subroutine run_and_snapshot(nml, h, hs, ht, wet, ok)
      !! Start, step `N_STEPS`, snapshot, tear down.  Sequential by
      !! construction: one engine at a time, each fully destroyed before
      !! the next is built.
      character(len=*), intent(in) :: nml
      real(wp), allocatable, intent(out) :: h(:, :, :), hs(:, :, :), ht(:, :, :)
      real(wp), allocatable, intent(out) :: wet(:, :)
      logical, intent(out) :: ok

      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      integer :: n, ierr
      real(wp) :: t

      call start_engine(nml, engine, cfg, ok)
      if (.not. ok) return
      t = 0.0_wp
      do n = 1, N_STEPS
         call advance(engine, cfg, t, ierr)
         if (ierr /= OCEAN_STATUS_OK) then
            ok = .false.
            call stop_engine(engine)
            return
         end if
         t = t + DT
      end do
      call pull_state(engine, h, hs, ht, wet)
      call stop_engine(engine)
   end subroutine run_and_snapshot

   pure subroutine wet_diff(a, b, wet, d)
      !! `max|a - b|` over WET cells only.  Land is excluded on purpose:
      !! the two runs are SUPPOSED to disagree there (different
      !! `bt_H_ref`), and the claim under test is that the disagreement
      !! stays there.
      real(wp), intent(in) :: a(:, :, :), b(:, :, :), wet(:, :)
      real(wp), intent(out) :: d
      integer :: i, j, k
      d = 0.0_wp
      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               if (wet(i, j) == 0.0_wp) cycle
               d = max(d, abs(a(i, j, k) - b(i, j, k)))
            end do
         end do
      end do
   end subroutine wet_diff

   pure function wet_scale(a, wet) result(s)
      !! `max|a|` over wet cells — the signal the bitwise claim is
      !! measured against.
      real(wp), intent(in) :: a(:, :, :), wet(:, :)
      real(wp) :: s
      integer :: i, j, k
      s = 0.0_wp
      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               if (wet(i, j) == 0.0_wp) cycle
               s = max(s, abs(a(i, j, k)))
            end do
         end do
      end do
   end function wet_scale

end module test_ocean_cavity_grounded_budget
