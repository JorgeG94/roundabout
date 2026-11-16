module test_ocean_cavity_engine
   !! Engine-driven ice-shelf cavity cases, built from namelists and stepped
   !! through the production `engine_*` sequence on ONE rank, with the
   !! interesting geometry placed ON the seams of the MPI bit-identity gate
   !! (`tests/mpi/test_ocean_decomp_bitid_mpi.F90`: 26 x 18, nghost = 3, so the
   !! 2 x 2 split lines are i = 13 | 14 and j = 9 | 10).
   !!
   !! The geometry (both cases): a Gaussian seamount centred on the domain
   !! centre (the corner of the 2 x 2 split) pierces the ice base, so a
   !! GROUNDED disc straddles both split lines; the shelf is the box
   !! `x <= draft_x1, y <= draft_y1` with a draft that slopes in x, so the
   !! calving front is an L (an x-normal leg crossing the j split line and a
   !! y-normal leg crossing the i split line) and is not symmetric; open
   !! water lies north and east of it.  Virtual-freshwater melt and
   !! quadratic top drag are ON.
   !!
   !! What the suite asserts, per case and split scheme:
   !!
   !!  * every prognostic stays finite over `N_STEPS` (the gate's 48) and
   !!    over twice that;
   !!  * the mass / salt / heat budget residuals the console prints stay at
   !!    round-off (`BUDGET_TOL`, the same 1e-12 `test_ocean_cavity_grounded_budget`
   !!    holds), step by step;
   !!  * NON-VACUITY: on BOTH sides of BOTH split lines (two columns / rows
   !!    each side) there are ice-covered wet cells, open-water wet cells,
   !!    grounded cells, cells with a non-zero melt rate, covered cells with a
   !!    non-zero top-drag stress and, for `z_fixed`, columns whose `k_top < nz`
   !!    (layers outcropping into the ice).  Any zero count fails the test, so
   !!    a geometry that drifts off the seams can not pass quietly.
   !!
   !! `cavity_case_nml` is public so the bit-identity gate can paste the same
   !! text: it opens with the gate's own `common` groups (nz_layers is the
   !! only difference: 8 for z_fixed) and, like the gate's other cases, leaves
   !! `&output_nml output_to_file = .false.` to the caller.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, RHO_WATER
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_console_stats, only: compute_total_h, compute_total_tracer, compute_total_content, &
                                      ocean_budget_src, ocean_budget_out, &
                                      ocean_budget_stage_weight, &
                                      ocean_salt_src_sum, ocean_heat_src_sum
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_cavity_engine_tests, cavity_case_nml

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: N_STEPS = 48
      !! The bit-identity gate's step count.
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: BUDGET_TOL = 1.0e-12_wp
      !! `test_ocean_cavity_grounded_budget`'s tolerance, unchanged.

   logical, save :: comm_inited = .false.
      !! One-shot guard for the comm env bring-up (shared per-test main finalises it).

   type :: band_counts_t
      !! Non-vacuity counters of one band (two columns or rows beside a split line).
      integer :: covered = 0
      integer :: open = 0
      integer :: grounded = 0
      integer :: melting = 0
      integer :: dragged = 0
      integer :: outcrop = 0
   end type band_counts_t

contains

   subroutine collect_ocean_cavity_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_zfixed_pred_corr", test_zfixed_pred_corr), &
                  new_unittest("cavity_zfixed_ssp_rk2", test_zfixed_ssp_rk2), &
                  new_unittest("cavity_sigma_pred_corr", test_sigma_pred_corr), &
                  new_unittest("cavity_sigma_ssp_rk2", test_sigma_ssp_rk2) &
                  ]
   end subroutine collect_ocean_cavity_engine_tests

   function cavity_case_nml(label, scheme, px, py) result(nml)
      !! Namelist of one case ("cavity_zfixed" | "cavity_sigma") on a px x py grid.
      character(len=*), intent(in) :: label, scheme
      integer, intent(in) :: px, py
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: spx, spy
      character(len=:), allocatable :: vc, nzl, cav, vm

      write (spx, '(i0)') px
      write (spy, '(i0)') py
      if (label == "cavity_zfixed") then
         ! nz = 8 over 1000 m: h_nominal = 125 m, so h_min_cavity = 250 m
         ! (validate_config: >= 2 h_nominal).  The draft (~300-500 m) cuts
         ! the top two to four nominal layers into fillers (k_top < nz).
         nzl = "8"
         vc = "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true. /"//NL
         cav = ", h_min_cavity = 250.0 /"//NL
         vm = "&ocean_vmix_nml use_kpp = .false., pp81_nu0 = 0.0, pp81_nu_bg = 1.0e-3, "// &
              "pp81_kappa_bg = 5.0e-5 /"//NL
      else
         nzl = "4"
         vc = "&vcoord_nml vcoord_type = 'sigma' /"//NL
         cav = ", h_min_cavity = 100.0 /"//NL
         vm = ""
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = "//nzl//" /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_surface = 3.0, T_init_bottom = 1.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 26, ny = 18, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = -1.409e-4 /"//NL// &
            vc//vm// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 400.0, "// &
            "slope_scale = 80000.0, wind_config = 'constant', taux_magnitude = 0.05 /"//NL// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 500.0, draft_slope = -4.0e-4, draft_x0 = 0.0, "// &
            "draft_x1 = 380000.0, draft_y1 = 260000.0"//cav// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0 /"//NL// &
            "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"//NL// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//NL// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
   end function cavity_case_nml

   ! ------------------------------------------------------------------
   ! Engine driving
   ! ------------------------------------------------------------------

   subroutine start_engine(nml, engine, cfg, ok)
      character(len=*), intent(in) :: nml
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      integer :: ierr

      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
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
      !! One outer step in the driver's order: dyn-core, ice block, finalize.
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

   subroutine totals(engine, total_mass, total_salt, total_heat, &
                     src_salt, src_heat, out_salt, out_heat, mass_out)
      !! The console's own budget terms, through its own reducers.
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
         total_salt = compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, &
                                           mt%areaT, ng)*RHO_WATER
         total_heat = compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, &
                                           mt%areaT, ng)*RHO_WATER

         src_salt = ocean_budget_src( &
                    ocean_salt_src_sum( &
                    compute_total_content(ms%salt_budget_surface, mt%areaT, ng), &
                    compute_total_content(ms%salt_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_salt = ocean_budget_out( &
                    compute_total_content(ms%salt_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_content(ms%salt_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)

         src_heat = ocean_budget_src( &
                    ocean_heat_src_sum( &
                    compute_total_content(ms%heat_budget_surface, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_geothermal, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_heat = ocean_budget_out( &
                    compute_total_content(ms%heat_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_content(ms%heat_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)

         mass_out = ms%mass_out
      end associate
   end subroutine totals

   pure function residual(total, ref, out_term, src_term) result(rel)
      !! The console's `Error` column: relative drift not explained by out/src.
      real(wp), intent(in) :: total, ref, out_term, src_term
      real(wp) :: rel
      if (abs(ref) > 1.0e-30_wp) then
         rel = ((total - ref) + out_term - src_term)/ref
      else
         rel = 0.0_wp
      end if
   end function residual

   function state_is_finite(engine) result(fin)
      !! All prognostics finite on the host copy (COMPONENT arrays only).
      type(ocean_engine_t), intent(inout) :: engine
      logical :: fin
      integer :: it

      fin = .true.
      !$acc update self(engine%state%multilayer%h_layer)
      !$acc update self(engine%state%multilayer%u_face_x_layer)
      !$acc update self(engine%state%multilayer%v_face_y_layer)
      fin = fin .and. all(ieee_is_finite(engine%state%multilayer%h_layer)) &
            .and. all(ieee_is_finite(engine%state%multilayer%u_face_x_layer)) &
            .and. all(ieee_is_finite(engine%state%multilayer%v_face_y_layer))
      do it = 1, size(engine%state%multilayer%tracers)
         !$acc update self(engine%state%multilayer%tracers(it)%c)
         fin = fin .and. all(ieee_is_finite(engine%state%multilayer%tracers(it)%c))
      end do
   end function state_is_finite

   subroutine count_bands(engine, zfixed, counts)
      !! Fill the four bands: (1) i-split west side, (2) i-split east side,
      !! (3) j-split south side, (4) j-split north side, each two columns /
      !! rows wide, running the whole other direction.
      type(ocean_engine_t), intent(inout) :: engine
      logical, intent(in) :: zfixed
      type(band_counts_t), intent(out) :: counts(4)

      integer :: ib, i, j, ng, nz, i0, i1, j0, j1
      logical :: wet, cov

      ng = engine%grid%nghost
      nz = engine%state%multilayer%nz_ml
      ! The test always runs melt and top drag, so all three arrays are mapped.
      !$acc update self(engine%state%multilayer%wet_mask)
      !$acc update self(engine%state%cavity_flux%melt)
      !$acc update self(engine%state%tdrag%stress_top)
      associate (ms => engine%state%multilayer, mt => engine%state%metrics, &
                 cf => engine%state%cavity_flux, td => engine%state%tdrag)
         associate (wm => ms%wet_mask, melt => cf%melt, st => td%stress_top)
            do ib = 1, 4
               select case (ib)
               case (1)
                  i0 = ng + NX_G/2 - 1; i1 = ng + NX_G/2; j0 = ng + 1; j1 = ng + NY_G
               case (2)
                  i0 = ng + NX_G/2 + 1; i1 = ng + NX_G/2 + 2; j0 = ng + 1; j1 = ng + NY_G
               case (3)
                  i0 = ng + 1; i1 = ng + NX_G; j0 = ng + NY_G/2 - 1; j1 = ng + NY_G/2
               case default
                  i0 = ng + 1; i1 = ng + NX_G; j0 = ng + NY_G/2 + 1; j1 = ng + NY_G/2 + 2
               end select
               do j = j0, j1
                  do i = i0, i1
                     wet = wm(i, j) > 0.5_wp
                     cov = mt%cover_frac(i, j) > 0.0_wp
                     if (cov .and. .not. wet .and. mt%z_draft(i, j) > 0.0_wp) &
                        counts(ib)%grounded = counts(ib)%grounded + 1
                     if (.not. wet) cycle
                     if (cov) then
                        counts(ib)%covered = counts(ib)%covered + 1
                        if (abs(melt(i, j)) > 0.0_wp) counts(ib)%melting = counts(ib)%melting + 1
                        if (st(i, j) > 0.0_wp) counts(ib)%dragged = counts(ib)%dragged + 1
                        if (zfixed .and. ms%k_top(i, j) < nz) counts(ib)%outcrop = counts(ib)%outcrop + 1
                     else
                        counts(ib)%open = counts(ib)%open + 1
                     end if
                  end do
               end do
            end do
         end associate
      end associate
   end subroutine count_bands

   subroutine check_counts(error, label, counts, zfixed)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      type(band_counts_t), intent(in) :: counts(4)
      logical, intent(in) :: zfixed
      character(len=*), parameter :: BAND(4) = [character(len=14) :: "i-split west", &
                                                                                   "i-split east", "j-split south", "j-split north"]
      integer :: ib
      character(len=96) :: msg

      do ib = 1, 4
         write (msg, '(5a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') trim(label), " ", trim(BAND(ib)), &
            ": cov/open/gnd/melt/drag/outcrop = ", "", counts(ib)%covered, "/", counts(ib)%open, "/", &
            counts(ib)%grounded, "/", counts(ib)%melting, "/", counts(ib)%dragged, "/", counts(ib)%outcrop
         write (*, '(a)') trim(msg)
         call check(error, counts(ib)%covered > 0, trim(label)//" "//trim(BAND(ib))//": ice-covered wet cells")
         if (allocated(error)) return
         call check(error, counts(ib)%open > 0, trim(label)//" "//trim(BAND(ib))//": open-water wet cells")
         if (allocated(error)) return
         call check(error, counts(ib)%grounded > 0, trim(label)//" "//trim(BAND(ib))//": grounded cells")
         if (allocated(error)) return
         call check(error, counts(ib)%melting > 0, trim(label)//" "//trim(BAND(ib))//": non-zero melt")
         if (allocated(error)) return
         call check(error, counts(ib)%dragged > 0, trim(label)//" "//trim(BAND(ib))//": non-zero top drag")
         if (allocated(error)) return
         if (zfixed) then
            call check(error, counts(ib)%outcrop > 0, &
                       trim(label)//" "//trim(BAND(ib))//": k_top < nz columns")
            if (allocated(error)) return
         end if
      end do
   end subroutine check_counts

   subroutine run_case(error, label, scheme)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, scheme

      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      type(band_counts_t) :: counts(4)
      logical :: ok
      integer :: n, ierr
      real(wp) :: t, m0, s0, h0, mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out
      real(wp) :: worst_m, worst_s, worst_h, r

      call start_engine(cavity_case_nml(label, scheme, 1, 1)// &
                        "&output_nml output_to_file = .false. /", engine, cfg, ok)
      call check(error, ok, label//"/"//scheme//": engine starts (validate_config + engine_setup)")
      if (allocated(error)) return

      gate: block
         call totals(engine, m0, s0, h0, src_s, src_h, out_s, out_h, m_out)
         worst_m = 0.0_wp; worst_s = 0.0_wp; worst_h = 0.0_wp
         t = 0.0_wp
         do n = 1, 2*N_STEPS
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//"/"//scheme//": step succeeds")
            if (allocated(error)) exit gate
            t = t + DT
            call totals(engine, mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out)
            r = residual(mt_, m0, m_out*RHO_WATER, 0.0_wp)
            worst_m = max(worst_m, abs(r))
            r = residual(st_, s0, out_s, src_s)
            worst_s = max(worst_s, abs(r))
            r = residual(ht_, h0, out_h, src_h)
            worst_h = max(worst_h, abs(r))
            if (n == N_STEPS) then
               call check(error, state_is_finite(engine), label//"/"//scheme//": finite at N_STEPS")
               if (allocated(error)) exit gate
               call count_bands(engine, label == "cavity_zfixed", counts)
               call check_counts(error, label//"/"//scheme, counts, label == "cavity_zfixed")
               if (allocated(error)) exit gate
            end if
         end do
         call check(error, state_is_finite(engine), label//"/"//scheme//": finite at 2 x N_STEPS")
         if (allocated(error)) exit gate
         write (*, '(4a,3(a,es10.3))') trim(label), "/", trim(scheme), " worst budget residuals", &
            " mass ", worst_m, " salt ", worst_s, " heat ", worst_h
         call check(error, worst_m <= BUDGET_TOL, label//"/"//scheme//": mass budget at round-off")
         if (allocated(error)) exit gate
         call check(error, worst_s <= BUDGET_TOL, label//"/"//scheme//": salt budget at round-off")
         if (allocated(error)) exit gate
         call check(error, worst_h <= BUDGET_TOL, label//"/"//scheme//": heat budget at round-off")
         if (allocated(error)) exit gate
         call check(error, abs(src_s) > 0.0_wp .and. abs(src_h) > 0.0_wp, &
                    label//"/"//scheme//": melt delivers a non-zero tracked source")
      end block gate

      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_case

   subroutine test_zfixed_pred_corr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "cavity_zfixed", "pred_corr")
   end subroutine test_zfixed_pred_corr

   subroutine test_zfixed_ssp_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "cavity_zfixed", "ssp_rk2")
   end subroutine test_zfixed_ssp_rk2

   subroutine test_sigma_pred_corr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "cavity_sigma", "pred_corr")
   end subroutine test_sigma_pred_corr

   subroutine test_sigma_ssp_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "cavity_sigma", "ssp_rk2")
   end subroutine test_sigma_ssp_rk2

end module test_ocean_cavity_engine
