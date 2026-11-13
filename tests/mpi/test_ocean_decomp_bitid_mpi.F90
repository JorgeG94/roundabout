!! Dycore bit-identity across domain decompositions, through the
!! production engine.
!!
!! "MPI works" means a decomposed run IS the serial run: every owned value
!! of every prognostic field bit-for-bit, not an integral that agrees to
!! round-off.  A round-off-level difference is where every seam bug
!! starts (a one-sided stencil at a seam face, an allreduce feeding back
!! into the state, a seam treated as a wall, a global-index formula
!! evaluated on local indices); integrals hide it for hundreds of steps.
!!
!! On EACH launch (ctest runs np = 1, 2, 4) the binary, for every case
!! below:
!!   1. runs the SERIAL reference (px = py = 1, compute_size = 1) on every
!!      rank — deterministic, no collective that could differ — through
!!      `engine_setup` -> `engine_enter_data` -> `engine_step` x N_STEPS;
!!   2. runs every px x py factorisation of the launched rank count
!!      (np = 2: 2x1, 1x2; np = 4: 4x1, 1x4, 2x2) on the same namelist;
!!   3. compares, BITWISE, each rank's OWNED window (ghosts excluded; a
!!      staggered face array includes both of its tile's edge faces, so a
!!      seam face is checked on both ranks that compute it) of every field
!!      the restart registry checkpoints — h, u, v, every tracer hTr, the
!!      barotropic prognostics, the KPP/EPBL/kappa-shear persistent state,
!!      rho, the OBC reservoirs — plus the BT end-of-step eta, against the
!!      matching window of the reference.
!!
!! Dev filter: `RDB_BITID_CASES=<substring>` runs only the cases whose label
!! contains it (e.g. `RDB_BITID_CASES=drain`).
!!
!! Cases (each under `pred_corr` AND `ssp_rk2`, 48 steps):
!!   * island_basin — closed Cartesian basin with an interior land block
!!     (static land mask), beta plane, `2gyre` wind, sigma;
!!   * weno7_pv — island_basin on the Sadourny enstrophy path with the weno7
!!     PV face interpolation at nghost = 5, the stencil radius + 1 its
!!     configure gate asks for (`pv_adv_required_nghost`);
!!   * periodic_channel_zstar — re-entrant channel over a seamount on z*
!!     (ALE remap every step), with porous barriers;
!!   * visc_rem_zstar — closed, cooled Cartesian spoon basin on z* (bed
!!     fillers on the shallow rim) with the visc_rem chain (`visc_rem_chain`
!!     + the implicit bed-drag fold): av_rem/bt_rem and the wt_u/renorm
!!     weights are per-face column sums built on the full face extent
!!     including ghosts, so a decomposition-sensitive seam in them shows
!!     here (D1 follow-up, 2026-10: this case used to drive the
!!     visc_rem-WEIGHTED BT-correction fold directly via the now-retired
!!     `correction_visc_rem`, whose ghost writes needed an explicit
!!     refresh before the vertical mixing read them; that fold has no
!!     MOM6 counterpart and is gone, so this specific ghost-staleness
!!     class has no decomposition coverage any more -- same posture as
!!     `accel_visc_rem`'s kernel, single-rank-tested only);
!!   * periodic_sponge — re-entrant channel, periodic west/east, with a
!!     relaxing sponge band on the closed north edge; besides the usual
!!     bitwise field comparison, its closed-budget mass/salt/heat totals
!!     are also compared bit-for-bit against the single-rank SERIAL
!!     reference (`check_periodic_sponge_serial_out`) — the regression
!!     gate for the periodic-seam / sponge-seam-ghost bug (see its
!!     docstring);
!!   * open_obc — tidal west edge (eta target) + Flather east edge,
!!     zstar_sigma, over a seamount;
!!   * obc_epbl — Flather east edge, sigma, EPBL, cooled surface, 24 x 16 at
!!     nghost = 4 (the compat-matrix cell): at -np 4 the 4x1 east tile is 14
!!     columns wide and the open-edge ghost columns it evolves sit in the
!!     vector remainder of the continuity loops; with FMA contraction on
!!     their last-bit rounding differed from the serial run and reached the
!!     owned boundary face through the horizontal viscosity (first owned
!!     difference at step 12, 6.7e-16, 3.6e-13 by step 24).  CPU builds now
!!     compile with contraction off (`cmake/compiler_flags.cmake`).  The
!!     case depends on the tile width (14) relative to the vector width
!!     (AVX2): on an AVX-512 target it may not catch a regression;
!!   * spherical — lon-lat sector, planetary Coriolis, spoon basin, Wright
!!     EOS, `2gyre` wind;
!!   * obc_radiation_sponge — Orlanski-radiating open south edge with tracer
!!     reservoirs, clamped north inflow, legacy relaxing sponge band west;
!!   * obc_sponge_east — the same with the sponge band on the EAST edge (west
!!     a wall), the combination once documented as drifting at round-off;
!!   * closures — the spherical case with EPBL (instead of KPP), Fox-Kemper
!!     MLE, GM + MEKE, Redi, kappa-shear, tidal mixing, convective
!!     adjustment, geothermal heating and tracer hdiff;
!!   * file_readers — the per-rank windowed readers (a periodic 360-degree
!!     MOM6 supergrid, a C-order bathymetry file with land, a z-level T/S
!!     IC), all written by rank 0 first, with the global-1-degree physics
!!     set (z_fixed + closed faces, fv_mom6, energy Coriolis, Wright),
!!     plus Redi at a NON-zero `khtr` — its open-window pairing reads the
!!     neighbour column across every seam face, so a stale ghost or a
!!     one-sided window shows here (the `closures` case enables Redi at
!!     the default `khtr = 0`, which builds the coefficients but applies
!!     no flux); with two coordinate twins on the same files:
!!     `file_readers_zstar_full` and `file_readers_zstar` (MOM6 z* + closed
!!     faces — the z_fixed staircase dilated per column every regrid);
!!   * sss_restore - the island basin with FILE sea-surface-salinity restoring
!!     (per-rank windowed two-record target, halo-exchanged, piston salt flux,
!!     net-zero shift through the EFP allreduce); needs NetCDF;
!!   * manizza - the island basin with the chlorophyll-dependent N-band
!!     shortwave (`sw_scheme = 'manizza05'`): per-rank windowed two-record
!!     chlorophyll, its seam-ghost exchange, the per-band opacities; needs
!!     NetCDF;
!!   * balanced_vortex / balanced_vortex_t0 - the isolated balanced cold-core
!!     vortex IC (`ic_config = 'balanced_vortex'`) on a closed beta-plane, centred
!!     on the 2x1/1x2 seams and the 2x2 corner; `_t0` takes zero steps so the
!!     seeded state itself (incl. ghosts) is compared;
!!   * drain_ratio / drain_ratio_h3 - windowed tracer drain (ratio 4, ppm / ppm_h3
!!     + diabatic_first = .false.) on a periodic channel, >1 sub-cycle pass;
!!   * drain_ratio_weno5 / _weno7 / _weno9 - the same drain with the WENO swept
!!     reconstruction at nghost 3 / 4 / 5 (rung ladder keyed on global edges);
!!   * dtbt_reset / dtbt_reset_period - `&ocean_bt_nml dtbt` re-derived every step / hourly;
!!     `dtbt` sits at an n_inner threshold so n_inner flips 14 <-> 15 mid-run (27x19 trench)
!!   * drain_ratio_fold / drain_ratio_fold_h3 - the same drain across the tripolar fold
!!     (distributed fold exchange for px > 1);
!!   * tracer_ppm_h3 - re-entrant channel (periodic x) with
!!     `tracer_advection_scheme = 'ppm_h3'` on the per-step path;
!!   * the OM3-parity feature cases, each the spherical spoon basin plus ONE
!!     feature: `gm_ebt_struct` (EBT vertical structure of KhTh), `gm_fgnv`
!!     (FGNV streamfunction), `meke_gm_work` (MEKE source from the GM work),
!!     `hvisc_thin_layer` (stress-form viscosity with the thin-layer ratio,
!!     harmonic corner thickness and T/corner-point Laplacian floors),
!!     `hvisc_biharm_stress` (the biharmonic folded into the stress tensor
!!     with the joint Laplacian + biharmonic bound; a second strain stencil
!!     on top of the Laplacian's),
!!     `hvisc_om3_floors` (latitude / velocity-scale / Coriolis-bounded
!!     viscosity), `vmix_no_interior` (`interior_closure = 'none'`),
!!     `redi_cfl_fold` (Redi diffusive-CFL subcycling, `n_itt > 1`, on the
!!     tripolar fold), `tidal_files` (tidal `tideamp`/`h2` read from a file),
!!     `tidal_polzin` (the Polzin vertical profile) and `bbl_lotw` (the
!!     law-of-the-wall bottom-boundary-layer tracer mixing) and
!!     `hvisc_kh_points` (stress-tensor viscosity with the static Laplacian
!!     floors built at the T centres and Bu corners), `diabatic_end`
!!     (`diabatic_first = .false.`, `dt_therm_ratio = 2`: the thermo tracer
!!     applies + ALE remap close the interval, pred_corr only; the windowed
!!     drain is not in play there); the two tidal cases
!!     need NetCDF;
!!   * sea_ice — sea ice on the ocean's decomposition: Winton thermo +
!!     ITD, EVP dynamics (CFL clip, `project_ci`) and the ice->ocean
!!     stress blend in a cooled periodic channel, `dt_therm_ratio = 2`;
!!     `sea_ice_transport` is the same with the category transport on.
!!     Every ice registry field is compared like the ocean's.  `run_one` calls
!!     `engine_step_ice` between the step and the finalize, as the driver
!!     does (a no-op for every other case).  Both ice cases also carry a
!!     WRITE/RESUME leg on the 4x1 and 2x2 layouts (2x1 on the 2-rank
!!     launch; `check_resume_leg`):
!!     the same decomposed run, checkpointed to per-rank restart files at
!!     step N_CKPT (odd, so mid thermo window) and resumed through
!!     `engine_setup`'s warm-restart path, must equal the straight
!!     decomposed run bitwise on every registry field -- FULL local arrays,
!!     ghosts included (a resume lands on the same decomposition).
!!   * tripolar — the analytic tripolar cap closed by the north fold, under a
!!     constant wind, with the double-Drake land reaching the fold line; its
!!     x splits (2x1, 4x1 = 8/8/7/7, 2x2) fold through the distributed fold
!!     exchange, its 1xN splits through the local kernels.
!!   * sea_ice_fold / sea_ice_fold_nodyn (A7, #143) -- sea ice (thermo + ITD +
!!     EVP + CAS transport + the stress blend) on the tripolar cap, so the
!!     category state, the flux diagnostics, the EVP ice velocity and
!!     stresses and the blended stress all fold across the north seam.  The
!!     serial reference must carry nonzero ice to the seam with an exactly
!!     antisymmetric v_ice fold row, every rank layout (2x1, 4x1 = 8/8/7/7,
!!     2x2 through the distributed fold; 1xN through the local kernels) must
!!     equal it bitwise with zero fold-buffer regrows, and the write/resume
!!     leg runs under fold + MPI.  `_nodyn` has `dynamics = .false.` so a
!!     thermo / transport seam bug is separable from an EVP one.
!!   * visc_rem_chain — the island_basin topology with the visc_rem chain
!!     live (`&ocean_vdiff_nml hvel_mom6 + bbl_glue`, `&ocean_bt_nml
!!     visc_rem_chain`) AND the frhat port's `frhat_scheme = "hybrid"`
!!     explicit (MOM6 `btcalc` HVEL_SCHEME=HYBRID; both now the PR-4
!!     defaults, spelled out here anyway so the case stays a decomp-bitid
!!     gate on the combination regardless of future default changes) --
!!     av_rem/bt_rem and the wt_u/renorm weights are per-face column sums
!!     built on the FULL face extent including ghosts, now via the HYBRID
!!     face-thickness closure instead of the plain arithmetic mean, so a
!!     decomposition-sensitive seam in EITHER the chain or the frhat sweep
!!     shows here exactly like every other compared field; the
!!     BT-correction fold itself stays UNIFORM (D1 follow-up — MOM6 never
!!     weights it).
!!   * cavity_zfixed / cavity_sigma / cavity_sigma_ratio2 / cavity_periodic /
!!     cavity_melt_only / cavity_combo (EPBL + hj99 + implicit_top_drag + tracer ratio 2) / cavity_gall — the ice-shelf cavity (dyn + virtual
!!     melt + quadratic top drag, `cavity_melt_only` without the drag) with a
!!     grounded disc and an L-shaped calving front on the 2x2 seams
!!     (`cavity_nml`).  On top of the registry fields the gate compares the
!!     cavity's own outputs (melt, u*, t_b, stress_shelf, stress_top) on owned
!!     cells, checks the SEAM-GHOST images of the four face twins
!!     `cover_u/v`, `k_top_u/v` against the serial reference, and fails when a
!!     non-vacuity count of the serial reference is 0.  `cavity_gall` has a
!!     tile that is entirely grounded on the 4x1 split; its global grounded
!!     fraction is under `grounded_max_frac`, so every layout must ACCEPT it,
!!     and the same draft against a lower limit must be REFUSED on every rank
!!     without a hang (`check_cavity_grounded_refusal`).
!! All are stratified with a boundary-layer scheme on, so the tiles exchange real
!! flow and real tracer structure.  26 x 18 cells (tripolar: 30 x 24),
!! nghost = 3 (weno7_pv: 26 x 24, nghost = 5): every factorisation above is
!! uneven somewhere.
!!
!! The barotropic march-in (`&ocean_bt_nml bt_halo > 0`) is NOT covered: it
!! is opt-in precisely because it is not bit-identical to the serial run over
!! variable bathymetry or with open boundaries (see `resolve_bt_halo`).
!!
!! Each bug this test found, and the fix that made it pass, is named in the
!! commit that fixed it; a regression anywhere prints the first mismatching
!! cell of every mismatching field on every rank.
#ifdef RDB_ENABLE_MPI
program test_ocean_decomp_bitid_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_error_ring, only: error_ring_clear, error_ring_count, error_ring_get
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, &
                               engine_teardown
   use rdb_ocean_state, only: ocean_state_build_restart_registry, ocean_state_restart_write
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm, &
                           comm_env_push_compute_comm, comm_env_pop_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use rdb_console_stats, only: console_stats_t, conservation_budget_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, ocean_budget_stage_weight
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_fold_exchange, only: ocean_fold_exchange_grow_count
#ifndef RDB_NO_NETCDF
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_def_var_3d, nc_enddef, nc_put_var_2d, rdb_def_var_1d, &
                            rdb_put_var_1d, output_rank_filename, ensure_directory_exists, &
                            nc_put_att
   use netcdf, only: nf90_put_var
#endif
   implicit none

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: N_STEPS = 48
   character(len=*), parameter :: DTBT_FRAC = "-0.84678989"
      !! `dtbt` of the `dtbt_reset*` cases (fraction of the largest stable substep).
   integer, parameter :: N_CKPT = 25
      !! Checkpoint step of the write/resume leg: odd, so with
      !! `dt_therm_ratio = 2` it lands mid thermo window.
   character(len=*), parameter :: RST_DIR = "bitid_resume_rst"
   real(wp), parameter :: DT = 900.0_wp
   integer, parameter :: MAXF = 128
   character(len=*), parameter :: SG_FILE = "bitid_supergrid.nc"
   character(len=*), parameter :: BATHY_FILE = "bitid_bathy.nc"
   character(len=*), parameter :: ZINIT_FILE = "bitid_zinit.nc"
   character(len=*), parameter :: SSS_FILE = "bitid_sss.nc"
   character(len=*), parameter :: TIDE_FILE = "bitid_tide.nc"
   character(len=*), parameter :: CHL_FILE = "bitid_chl.nc"
   character(len=*), parameter :: ATM_WARM = "bitid_atm_warm.nc"
   character(len=*), parameter :: ATM_COLD = "bitid_atm_cold.nc"
   character(len=*), parameter :: ATM_CRSE = "bitid_atm_coarse.nc"
   character(len=*), parameter :: ATM_WARM_LIST = "bitid_atm_warm_files.txt"
   character(len=*), parameter :: ATM_CRSE_LIST = "bitid_atm_coarse_files.txt"
   character(len=*), parameter :: ATM_COLD_LIST = "bitid_atm_cold_files.txt"

   type :: field_t
      !! One compared field: a host copy of the whole local array (2-D
      !! fields are stored with a unit third extent) + its tag.
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snap_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
      integer :: io = 0, jo = 0, nxl = 0, nyl = 0, ng = 0
      integer :: ngrow = 0
         !! Mid-run regrows of the fold-exchange buffers (must stay 0).
      integer :: drain_iter = 0
         !! Largest windowed-drain pass count used (0 when the drain never ran).
      integer :: ncav(6) = 0
         !! Cavity non-vacuity counts over the OWNED cells: covered wet, open wet,
         !! grounded, melting, top-drag stress > 0, k_top < nz.
      type(conservation_budget_t) :: bud
         !! Closed-budget mass/salt/heat out+src totals at the final step —
         !! EFP-reduced (`reproducing_sums`, default on), so these must be
         !! bit-identical between the serial reference and every
         !! decomposition on a periodic configuration (see `compare_budget`).
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail, ic
   integer :: n_steps_cur = N_STEPS
      !! Steps `run_one` takes; 0 for the `*_t0` labels (compare the seeded IC itself).
   character(len=256) :: case_filter = ""
      !! `RDB_BITID_CASES`: run only the cases whose name contains this substring.
   type(comm_t) :: comm
   character(len=16), parameter :: SCHEMES(2) = [character(len=16) :: "pred_corr", "ssp_rk2"]
   character(len=24), parameter :: CASES(61) = [character(len=24) :: &
                                                "island_basin", "weno7_pv", &
                                                "periodic_channel_zstar", &
                                                "visc_rem_zstar", &
                                                "periodic_sponge", &
                                                "open_obc", "spherical", "obc_radiation_sponge", "obc_sponge_east", &
                                                "closures", "file_readers", &
                                                "file_readers_zstar_full", &
                                                "file_readers_zstar", "sea_ice", &
                                                "sea_ice_transport", "tripolar", &
                                                "visc_rem_chain", "sea_ice_fold", &
                                                "sea_ice_fold_nodyn", "sss_restore", &
                                                "gm_ebt_struct", "gm_fgnv", "meke_gm_work", &
                                                "hvisc_om3_floors", "vmix_no_interior", &
                                                "redi_cfl_fold", "tidal_files", &
                                                "tidal_polzin", "bbl_lotw", "tracer_ppm_h3", &
                                                "hvisc_kh_points", "diabatic_end", &
                                                "hvisc_thin_layer", "hvisc_biharm_stress", &
                                                "balanced_vortex_t0", "balanced_vortex", &
                                                "drain_ratio", "drain_ratio_h3", &
                                                "drain_ratio_weno5", "drain_ratio_weno7", &
                                                "drain_ratio_weno9", "drain_ratio_fold", &
                                                "drain_ratio_fold_h3", "eulz", "eulz_periodic", &
                                                "eulz_tripolar", "eulz_epbl_mle", "bulk_flux", &
                                                "bulk_flux_ice", "bulk_flux_fold", "manizza", "dtbt_reset", &
                                                "dtbt_reset_period", &
                                                "cavity_zfixed", "cavity_sigma", &
                                                "cavity_sigma_ratio2", "cavity_gall", "cavity_melt_only", "cavity_periodic", &
                                                "cavity_combo", "obc_epbl"]

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0
#ifndef RDB_NO_NETCDF
   if (rank == 0) call write_input_files()
   call comm%barrier()
#endif

   ! Dev filter: RDB_BITID_CASES=<substring> runs only matching case labels.
   call get_environment_variable("RDB_BITID_CASES", case_filter)
   do ic = 1, size(CASES)
      if (len_trim(case_filter) > 0) then
         if (index(trim(CASES(ic)), trim(case_filter)) == 0) cycle
      end if
#ifdef RDB_NO_NETCDF
      if (trim(CASES(ic)) == "file_readers" .or. &
          trim(CASES(ic)) == "file_readers_zstar_full" .or. &
          trim(CASES(ic)) == "sss_restore" .or. &
          trim(CASES(ic)) == "manizza" .or. &
          trim(CASES(ic)) == "tidal_files" .or. &
          trim(CASES(ic)) == "tidal_polzin" .or. &
          index(trim(CASES(ic)), "bulk_flux") == 1 .or. &
          trim(CASES(ic)) == "file_readers_zstar") cycle
#endif
      ! eulerian_z is refused under pred_corr (validate_config): ssp_rk2 only.
      if (CASES(ic) (1:4) /= "eulz") call run_case(trim(CASES(ic)), trim(SCHEMES(1)))
      ! diabatic_first = .false. is refused under ssp_rk2 (validate_config).
      ! (drain_ratio_h3 pins it too.)
      if (trim(CASES(ic)) == "diabatic_end" .or. trim(CASES(ic)) == "drain_ratio_h3" .or. &
          trim(CASES(ic)) == "drain_ratio_fold_h3") cycle
      call run_case(trim(CASES(ic)), trim(SCHEMES(2)))
   end do

   if (nprocs >= 2) call check_single_rank_fences()
   if (len_trim(case_filter) == 0 .or. index("cavity_grounded_refused", trim(case_filter)) > 0) &
      call check_cavity_grounded_refusal()

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_decomp_bitid_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_decomp_bitid_mpi: all checks PASSED (nprocs=", &
            nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   function case_nml(label, scheme, px, py) result(nml)
      !! The namelist of one case on a px x py process grid.
      character(len=*), intent(in) :: label, scheme
      integer, intent(in) :: px, py
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: spx, spy, snx, sny, sng, snz
      character(len=:), allocatable :: common, bt_extra, bt_head, ice_dyn, hv_extra, drain_extra
      character(len=16) :: rad_earth_s, wind_s, wind_y_s
      character(len=8) :: sponge_w, sponge_e
      character(len=64) :: wl

      ! The visc_rem chain rides on the common &ocean_bt_nml group (one
      ! group per namelist).  D1 follow-up (2026-10): the weighted
      ! BT-correction fold (`correction_visc_rem`) is RETIRED -- MOM6
      ! never weights it, and it is what NaNs the 1-degree Southern Ocean
      ! z* open-step case under bbl_glue -- so both cases below now
      ! exercise `visc_rem_chain` (producer + bt_rem_from_av_rem + wt_u
      ! forcing + renorm_visc_rem, UNIFORM BT-correction fold) instead.
      ice_dyn = "dynamics = .true., evp_sub_steps = 30, cfl_trunc = 0.01, "// &
                "project_ci = .true., "
      if (label == "sea_ice_fold_nodyn") ice_dyn = "dynamics = .false., "
      bt_extra = ""
      bt_head = "auto_n_inner = .true., "
      ! Per-step (period 0) / hourly (period 3600 s) re-derivation of n_inner from
      ! the live thicknesses; the fraction sits at an n_inner threshold (see case).
      if (label == "dtbt_reset") bt_head = "dtbt = "//DTBT_FRAC//", dtbt_reset_period = 0.0, "
      if (label == "dtbt_reset_period") bt_head = "dtbt = "//DTBT_FRAC//", dtbt_reset_period = 3600.0, "
      if (label == "visc_rem_zstar") bt_extra = ", visc_rem_chain = .true."
      if (label == "visc_rem_chain") then
         bt_extra = ", visc_rem_chain = .true., frhat_scheme = 'hybrid'"
      end if
      hv_extra = ""
      ! dx = 250 m puts the fixed nu_h = 200 far past the viscous CFL: clamp it.
      if (index(label, "drain") == 1) hv_extra = ", bound_kh = .true."
      drain_extra = ""
      if (index(label, "drain_ratio_weno") == 1) then
         drain_extra = ", tracer_recon = '"//label(len_trim(label) - 4:len_trim(label))//"'"
      end if
      rad_earth_s = "6.378e6"
      wind_s = "0.05"
      wind_y_s = "0.0"
      if (index(label, "drain_ratio_fold") == 1) then
         rad_earth_s = "2.0e4"
         wind_s = "0.3"
         wind_y_s = "0.3"   ! meridional stress drives flow INTO the fold row
      end if
      if (label == "drain_ratio_h3" .or. label == "drain_ratio_fold_h3") then
         drain_extra = ", tracer_advection_scheme = 'ppm_h3', diabatic_first = .false."
      end if
      if (label == "hvisc_om3_floors") then
         hv_extra = ", kh_sin_lat = 2000.0, kh_pwr_of_sine = 4.0, ah_vel_scale = 0.01, "// &
                    "bound_coriolis_biharm = .true., bound_coriolis_vel = 6.0"
      end if
      if (label == "hvisc_kh_points") then
         ! Stress-tensor path with the static Laplacian floors built AT the T
         ! centres and Bu corners (local latitude, local harmonic length).
         hv_extra = ", stress_tensor = .true., kh_bg_at_points = .true., "// &
                    "kh_vel_scale = 0.01, kh_sin_lat = 2000.0, kh_pwr_of_sine = 4.0"
      end if
      if (label == "hvisc_thin_layer") then
         hv_extra = ", stress_tensor = .true., bound_thin_layer = .true., "// &
                    "kh_bg_at_points = .true., kh_vel_scale = 0.01, "// &
                    "kh_sin_lat = 2000.0, kh_pwr_of_sine = 4.0"
      end if
      if (label == "hvisc_biharm_stress") then
         hv_extra = ", stress_tensor = .true., biharm_stress = .true., nu_4 = 1.0e9"
      end if
      ! biharm_stress chains four stencil levels: validate_config needs nghost >= 4.
      sng = "3"
      if (label == "hvisc_biharm_stress") sng = "4"
      ! WENO drain rungs need nghost >= 3 / 4 / 5 (weno5 / 7 / 9).
      if (label == "drain_ratio_weno7") sng = "4"
      if (label == "drain_ratio_weno9") sng = "5"
      write (spx, '(i0)') px
      write (spy, '(i0)') py
      write (snx, '(i0)') NX_G
      write (sny, '(i0)') NY_G
      snz = "4"
      if (label == "manizza") snz = "8"
      ! weno9 needs nghost = 5: widen the grid so the np = 4 layouts (1x4 ->
      ! ny/4 rows, 4x1 -> nx/4 columns) stay at least nghost wide.
      if (label == "drain_ratio_weno9") then
         snx = "32"
         sny = "24"
      end if
      ! Odd extents put the trench centre in the middle of a cell, so the CFL-controlling
      ! (deepest) cell is unique and sits beside the 2-way seam on both axes.
      if (index(label, "dtbt_reset") == 1) then
         snx = "27"
         sny = "19"
      end if
      common = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = "//trim(snz)//" /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml "//bt_head//"split_scheme = '"//scheme//"'"// &
               bt_extra//" /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true."//hv_extra//" /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               ""

      if (index(label, "cavity_") == 1) then
         nml = cavity_nml(label, scheme, spx, spy)//"&output_nml output_to_file = .false. /"//NL
         return
      end if

      select case (label)
      case ("island_basin")
         ! Closed Cartesian basin, interior land block (static land mask),
         ! beta plane, eastward wind onto the island.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("sss_restore")
         ! The island basin with FILE sea-surface-salinity restoring: the
         ! per-rank windowed read of a two-record (t = 0 and 86400 s) target,
         ! its halo exchange, the piston salt flux, and the net-zero
         ! (global-mean shift) adjustment through the EFP allreduce.  The
         ! target varies in both directions so a wrong window cannot match.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_restore_nml enable_restore_salt = .true., piston_s = 0.5, "// &
               "sss_file = '"//SSS_FILE//"', sss_var = 'salt', sss_as_sflux = .true., "// &
               "sss_adjust_net_zero = .true. /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("manizza")
         ! The island basin with the Manizza (2005) chlorophyll shortwave: the
         ! per-rank windowed read of a two-record spatially non-uniform,
         ! non-symmetric chlorophyll climatology, its seam-ghost exchange, and
         ! the per-band opacities feeding the deposition kernel and the KPP
         ! B_0 coupling.  A wrong window or a stale ghost changes the heating.
         ! Shallow (100 m) and 8 layers deep so the blue AND red bands cross
         ! several interfaces: the opacities, not just the surface layer, drive
         ! the heating.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 100.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = 400.0, "// &
               "sw_pen_frac = 1.0, sw_scheme = 'manizza05', "// &
               "sw_chl_file = '"//CHL_FILE//"', sw_chl_var = 'CHL' /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("visc_rem_chain")
         ! The island basin with the visc_rem chain live end to end: MOM6
         ! BOTTOMDRAGLAW glue (hvel_mom6 + bbl_glue) feeding the producer,
         ! the wt_u BT-forcing weight, the continuity renormaliser and
         ! bt_rem = av_rem**(1/n_inner) substep damping (all via
         ! `visc_rem_chain`, bt_extra above) -- the BT-correction fold
         ! itself stays UNIFORM (D1 follow-up).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_vdiff_nml hvel_mom6 = .true., bbl_glue = .true. /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("balanced_vortex", "balanced_vortex_t0")
         ! Isolated balanced cold-core vortex (ic_config = 'balanced_vortex') on a
         ! closed beta-plane, centred on the 2x1/1x2 seams and the 2x2 four-way
         ! corner (x0 = 13 dx, y0 = 9 dx), R = 3 dx so the core also straddles the
         ! 4x1 (x = 7, 14, 20 dx) and 1x4 (y = 5, 10, 14 dx) seams.  Strong beta
         ! so f(y) differs rank to rank.  The `_t0` twin takes zero steps: the
         ! seeded IC (h, T, face velocities incl. ghosts) alone must match.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'flat', max_depth = 1000.0, "// &
               "coriolis_beta = 1.0e-10, coriolis_y_ref = 180000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
               "&ocean_eos_nml eos = 'linear' /"//NL// &
               "&ocean_ic_nml ic_config = 'balanced_vortex', alpha_T = 0.17, "// &
               "rho_0 = 1035.0, vortex_x0 = 260000.0, vortex_y0 = 180000.0, "// &
               "vortex_radius = 60000.0, vortex_hv = 300.0, vortex_dT = -2.0, "// &
               "vortex_T_surf = 20.0, vortex_dT_dz = 1.5e-2 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("weno7_pv")
         ! The island basin on the Sadourny enstrophy path with the weno7 PV
         ! face interpolation, at its gate's halo (stencil radius 4 + 1 =
         ! nghost 5).  At nghost = 4 (radius) the stencil's edge corners are
         ! built from the outermost ghost ring of u/v and the decomposed run
         ! drifts at the last bit (compat-matrix row `decomp_weno_pv`).  Its
         ! own ny = 24, so every 1xN tile (1x4: 6 rows) is wider than the
         ! 5-deep ghost band it sends.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = 24, nghost = 5, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.15, wind_config = '2gyre', taux_magnitude = 0.1, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_coriolis_nml form = 'sadourny', pv_adv_scheme = 'weno7' /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("eulz", "eulz_periodic")
         ! Plain eulerian_z (per-step vertical advection pins h to the
         ! reference; no EPBL, no MLE): the vertical advection runs after the
         ! post-continuity exchange, so a ghost column must be refreshed from
         ! its owner afterwards.  "eulz" = closed seamount basin
         ! (walls: seam ghosts only); "eulz_periodic" = the periodic channel
         ! of periodic_channel_zstar (wrap ghosts, also on one rank).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = -1.0e-4, wind_stress_x = 0.1 /"//NL// &
               "&vcoord_nml vcoord_type = 'eulerian_z' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL
         if (label == "eulz_periodic") then
            nml = nml//"&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
                  "north = 'wall' /"//NL
         else
            nml = nml//"&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
         end if
      case ("periodic_channel_zstar")
         ! Re-entrant channel over a seamount, z* (ALE remap every step).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = -1.0e-4, wind_stress_x = 0.1 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_porous_nml enable = .true. /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("visc_rem_zstar")
         ! Closed, cooled Cartesian spoon basin, MOM6 z* (the 300 m rim
         ! carries bed fillers), with the visc_rem chain (producer +
         ! forcing/renorm/bt_rem_from -- D1 follow-up: no longer the
         ! weighted BT-correction fold, which is retired).  visc_rem
         ! needs the implicit drag fold, which refuses an HBBL-distributed
         ! drag: bed-only drag.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
               "&ocean_topo_nml topo_config = 'spoon', max_depth = 3000.0, "// &
               "edge_depth = 300.0, slope_scale = 300000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//NL// &
               "&ocean_vdiff_nml implicit_drag = .true. /"//NL// &
               "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3, hbbl = 0.0 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      case ("periodic_sponge")
         ! Re-entrant channel, periodic west/east, with a relaxing sponge
         ! band on the closed north edge -- the shape of the real Southern
         ! Ocean 1-degree configuration the periodic-seam / sponge-seam-ghost
         ! bug (this test's `bud` capture) was found on: a periodic axis
         ! carrying a physical-span-only relaxation on the orthogonal one.
         ! Legacy per-edge band sponge (not map-driven): simplest namelist
         ! that reaches `ocean_sponge_apply`/`ocean_sponge_apply_tracers`
         ! and the seam-ghost refresh gated on `sponge_seam` in
         ! `rdb_ocean_dyn.F90`.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.2 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'sponge', sponge_width = 3, sponge_strength = 1.0e-2, "// &
               "sponge_relax_tracers = .true. /"//NL
      case ("drain_ratio", "drain_ratio_h3", "drain_ratio_weno5", "drain_ratio_weno7", &
            "drain_ratio_weno9")
         ! Windowed tracer-advection drain (`dt_tracer_advect_ratio = 4`, 12
         ! windows in N_STEPS) on a re-entrant seamount channel (periodic
         ! x: the drain reads tile and wrap ghosts), a small dx so the
         ! accumulated window Courant number exceeds 0.5 and the sub-cycle
         ! runs more than one pass.  `drain_ratio_h3` swaps in the
         ! `ppm_h3` drain reconstruction and the OM3 interval order
         ! (`diabatic_first = .false.`, pred_corr only).
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = "//trim(sng)//", "// &
               "dx = 250.0, dy = 250.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 1.0 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 1000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
               "&ocean_vmix_nml dt_therm_ratio = 4, dt_tracer_advect_ratio = 4"// &
               ""//trim(drain_extra)//" /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("dtbt_reset", "dtbt_reset_period")
         ! Re-entrant channel with a Gaussian TRENCH (peak_depth > max_depth) at the
         ! domain centre, which on the odd 27 x 19 grid is the middle of cell (14, 10):
         ! the unique deepest, hence CFL-controlling, cell.  It is the first/last
         ! column and row of a 2x1 / 1x2 tile, so its limit reads h across the seam.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.2 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 3000.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("tracer_ppm_h3")
         ! Re-entrant channel (periodic west/east) over a seamount, sigma, with
         ! the ppm_h3 horizontal tracer advection (Huynh edges + limiter, a
         ! 4-point stencil: radius 2 < nghost 3).  Periodic x makes the
         ! per-step path read tracer ghosts across the seam on every tile
         ! boundary and across the wrap.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.2 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_vmix_nml tracer_advection_scheme = 'ppm_h3' /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("open_obc")
         ! Tidal west edge (eta target) + Flather open east edge.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.05, wind_stress_y = 0.02 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar_sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 200.0, slope_scale = 40000.0 /"//NL// &
               "&ocean_bc_nml west = 'tidal', east = 'open', south = 'wall', north = 'wall', "// &
               "west_n_tidal = 1, west_tidal_amp = 0.5, west_tidal_phase = 0.0, "// &
               "west_tidal_omega = 1.4051890e-4 /"//NL
      case ("obc_epbl")
         ! Flather east edge, sigma, EPBL, cooled surface: the open-edge ghost
         ! columns are evolved by the continuity chain and the owned boundary
         ! face reads them through the horizontal viscosity, so with FMA
         ! contraction on, a vector body vs a differently-contracted remainder
         ! made a ghost cell's rounding depend on its position in the TILE
         ! (nvfortran CPU, 4x1, AVX2; the mismatch list starts with
         ! bt_visc_rem_u/v and epbl_kd_int only because it is alphabetical --
         ! 56 fields differ).  Caught by the compat matrix DECOMP leg.  Written
         ! out in full (not on `common`): the rounding event needs this
         ! stratification, forcing and tile width (24 x 16, nghost 4, 10
         ! layers).  "14 columns" is STORAGE (6 owned + 2 x 4 ghosts).  CI pins
         ! `-march=x86-64-v3` / `-tp=haswell` (AVX2), where that width puts the
         ! ghosts in the remainder; on an AVX-512 build (`-march=native`) the
         ! tile may sit differently in the vector width and this case may then
         ! not catch a regression.  That caveat does not apply to CI.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 21600.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 10 /"//NL// &
               "&tracer_nml initial_temperature = 9.0, initial_salinity = 34.8, "// &
               "T_init_surface = 16.0, T_init_bottom = 3.0, "// &
               "S_init_surface = 34.6, S_init_bottom = 35.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, nu_4 = 1.0e11 /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&ocean_eos_nml eos = 'linear' /"//NL// &
               "&ocean_ic_nml T_ref = 10.0, S_ref = 35.0, alpha_T = 0.2, beta_S = 0.77, "// &
               "rho_0 = 1027.0 /"//NL// &
               "&ocean_bdrag_nml form = 'quadratic', cd = 0.003, bg_vel = 0.05, hbbl = 0.0 /"//NL// &
               "&ocean_vdiff_nml implicit_drag = .true. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 4, dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.0 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'spoon', max_depth = 2000.0, "// &
               "edge_depth = 300.0, slope_scale = 150000.0, wind_config = 'constant', "// &
               "coriolis_beta = 2.0e-11, coriolis_y_ref = 160000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//NL// &
               "&ocean_vmix_nml use_kpp = .false. /"//NL// &
               "&ocean_epbl_nml enable = .true., mld_use_prev_guess = .true., "// &
               "mstar_scheme = 'om4' /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'open', south = 'wall', north = 'wall' /"//NL
      case ("tripolar", "redi_cfl_fold", "drain_ratio_fold", "drain_ratio_fold_h3", "eulz_tripolar")
         ! Analytic tripolar cap (59N-83N, phi_join 74N) closed by the north
         ! fold, periodic in x, under a constant wind, with the double-Drake
         ! land walls reaching the fold line at asymmetric columns (1 and
         ! nx/4 mirror to nx and 3nx/4+1).  Its own x extent, nx = 30 at
         ! 12 deg (360/26 is inexact): every split of the launched rank
         ! count runs, so 2x1 and 4x1 (8/8/7/7 — uneven) and 2x2 fold
         ! through the distributed fold exchange, 1xN through the local
         ! kernels.  `compare` is size-agnostic.
         nml = common// &
               "&grid_nml nx = 30, ny = 24, nghost = 3, dx = 12.0, dy = 1.0 /"//NL// &
               "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, lat_south = 59.0, "// &
               "phi_join = 74.0, lon_pole = 0.0, rad_earth = "//trim(rad_earth_s)//", "// &
               "coriolis_scheme = 'planetary' /"//NL// &
               "&physics_nml wind_stress_x = "//trim(wind_s)//", wind_stress_y = "//trim(wind_y_s)//" /"//NL// &
               "&vcoord_nml vcoord_type = '"//trim(merge("eulerian_z", "sigma     ", label == "eulz_tripolar"))//"' /"//NL// &
               "&ocean_topo_nml topo_config = 'double_drake', max_depth = 1000.0, "// &
               "slope_scale = 0.2 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'tripolar_fold' /"//NL
         if (label == "drain_ratio_fold" .or. label == "drain_ratio_fold_h3") then
            ! Windowed tracer drain (ratio 4) across the tripolar fold.  A
            ! 1/100-size planet (rad_earth) brings the accumulated window
            ! Courant number above 0.5 so the sub-cycle runs >1 pass, and
            ! the cap rows (fold seam) carry the flow; thermodynamics on.
            nml = nml// &
                  "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
                  "&ocean_vmix_nml dt_therm_ratio = 4, dt_tracer_advect_ratio = 4 /"//NL
         end if
         if (label == "redi_cfl_fold") then
            ! Redi at a diffusivity far above the explicit-diffusion limit
            ! with `check_diffusive_cfl`: the cell CFL (a cross-rank max) sets
            ! n_itt > 1 passes, each followed by a tracer-halo refresh that
            ! must cross the fold and the periodic seam.
            nml = nml// &
                  "&ocean_eos_nml eos = 'wright' /"//NL// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_redi_nml enable = .true., khtr = 2.0e6, "// &
                  "check_diffusive_cfl = .true. /"//NL
         end if
      case ("spherical", "closures", "gm_ebt_struct", "gm_fgnv", "meke_gm_work", &
            "hvisc_om3_floors", "vmix_no_interior", "tidal_files", "tidal_polzin", &
            "bbl_lotw", "hvisc_kh_points", "diabatic_end", "hvisc_thin_layer", &
            "hvisc_biharm_stress", "eulz_epbl_mle")
         ! Lon-lat sector, spoon basin, Wright EOS.  "closures" adds the
         ! lateral and vertical parameterisation set on top (EPBL instead of
         ! KPP, Fox-Kemper MLE, GM + MEKE, Redi, kappa-shear, tidal mixing,
         ! convective adjustment, geothermal heating, tracer hdiff).  The
         ! OM3-parity labels (the `select case` after the block below) each
         ! add ONE feature to this basin.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = "//trim(sng)//", "// &
               "dx = 1.0, dy = 1.0 /"//NL// &
               "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, lat_south = 20.0, "// &
               "rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
               "&vcoord_nml vcoord_type = '"// &
               trim(merge("zstar     ", merge("eulerian_z", "sigma     ", label == "eulz_epbl_mle"), &
                          label == "diabatic_end"))// &
               "' /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_topo_nml topo_config = 'spoon', max_depth = 3000.0, "// &
               "edge_depth = 300.0, slope_scale = 300000.0, wind_config = '2gyre', "// &
               "taux_magnitude = 0.1 /"//NL// &
               "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
         if (label == "closures") then
            nml = nml// &
                  "&ocean_vmix_nml use_kpp = .false. /"//NL// &
                  "&ocean_epbl_nml enable = .true. /"//NL// &
                  "&ocean_foxkemper_nml enable = .true. /"//NL// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_gm_nml enable = .true. /"//NL// &
                  "&ocean_meke_nml enable = .true. /"//NL// &
                  "&ocean_redi_nml enable = .true., khtr = 100.0 /"//NL// &
                  "&ocean_kappa_shear_nml enable = .true. /"//NL// &
                  "&ocean_tidal_mixing_nml enable = .true., e_uniform = 1.0e-3 /"//NL// &
                  "&ocean_conv_nml enable = .true. /"//NL// &
                  "&ocean_geothermal_nml enable = .true. /"//NL// &
                  "&ocean_hdiff_nml kappa_h = 100.0 /"//NL
         end if
         select case (label)
         case ("eulz_epbl_mle")
            ! eulerian_z (legacy per-stage vertical advection + h rescale)
            ! with EPBL and Fox-Kemper MLE: the pair that once drifted at
            ! round-off under ssp_rk2.
            nml = nml// &
                  "&ocean_vmix_nml use_kpp = .false. /"//NL// &
                  "&ocean_epbl_nml enable = .true., mstar_scheme = 'om4', "// &
                  "mld_use_prev_guess = .true. /"//NL// &
                  "&ocean_foxkemper_nml enable = .true., use_mom_mixrate = .true., "// &
                  "mld_decay_time = 86400.0 /"//NL
         case ("gm_ebt_struct")
            ! GM with the EBT vertical structure of KhTh (the structure is the
            ! wave-speed slot's `ebt_struct`, a per-cell field read at both
            ! cells of every face).
            nml = nml// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_wavespeed_nml enable = .true., use_ebt = .true. /"//NL// &
                  "&ocean_gm_nml enable = .true., khth = 500.0, use_ebt_struct = .true. /"//NL
         case ("gm_fgnv")
            ! GM with the FGNV streamfunction (a tridiagonal BVP per face
            ! column on the cg1 of the two bounding cells).
            nml = nml// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_wavespeed_nml enable = .true. /"//NL// &
                  "&ocean_gm_nml enable = .true., khth = 500.0, use_fgnv = .true., "// &
                  "fgnv_filter_scale = 0.1 /"//NL
         case ("meke_gm_work")
            ! MEKE sourced from the work of the final bolus transports, fed
            ! back into the VarMix KhTh (`khth_fac`) so the source reaches
            ! the state: with the feedback off the MEKE field alone changes.
            nml = nml// &
                  "&ocean_slopes_nml enable = .true. /"//NL// &
                  "&ocean_wavespeed_nml enable = .true. /"//NL// &
                  "&ocean_varmix_nml enable = .true., khth = 500.0 /"//NL// &
                  "&ocean_gm_nml enable = .true., khth = 500.0, gm_src_form = 'work' /"//NL// &
                  "&ocean_meke_nml enable = .true., gmcoeff = 0.15, damping = 1.0e-6, "// &
                  "khth_fac = 1.0, alpha_grid = 1.0 /"//NL
         case ("vmix_no_interior")
            ! No interior (shear / PP81) closure: KPP boundary layer only.
            nml = nml// &
                  "&ocean_vmix_nml interior_closure = 'none' /"//NL
         case ("tidal_files", "tidal_polzin")
            ! Tidal mixing with `tideamp` and `h2` read from a whole-grid
            ! file (each rank reads its own band); `tidal_polzin` adds the
            ! Polzin stretched-height vertical profile.
            nml = nml// &
                  "&ocean_tidal_mixing_nml enable = .true., e_compute = .true., "// &
                  "tideamp_source = 'file', tideamp_file = '"//TIDE_FILE//"', "// &
                  "h2_source = 'file', h2_file = '"//TIDE_FILE//"', "// &
                  "kd_max = 0.1, energy_limiter = 'reinject', bed_share = 'lump'"
            if (label == "tidal_polzin") nml = nml//", profile = 'polzin'"
            nml = nml//" /"//NL
         case ("bbl_lotw")
            ! Law-of-the-wall BBL tracer mixing: needs the per-face viscous
            ! BBL (bbl_glue + a bottom drag with hbbl > 0) under it.
            nml = nml// &
                  "&ocean_vdiff_nml hvel_mom6 = .true., bbl_glue = .true. /"//NL// &
                  "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3, hbbl = 10.0, "// &
                  "bg_vel = 0.1 /"//NL// &
                  "&ocean_bbl_mixing_nml enable = .true. /"//NL
         case ("diabatic_end")
            ! OM3 ordering of the thermo interval: surface heat flux +
            ! KPP vertical tracer diffusion applied at the interval end,
            ! then the ALE remap (pred_corr only; the loop skips ssp_rk2).
            nml = nml// &
                  "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
                  "&ocean_vmix_nml dt_therm_ratio = 2, diabatic_first = .false. /"//NL
         end select
      case ("obc_radiation_sponge", "obc_sponge_east")
         ! Orlanski-radiating open south edge with tracer reservoirs, a
         ! clamped (inflow) north edge, and a relaxing sponge band west.
         sponge_w = "sponge"
         sponge_e = "wall"
         if (label == "obc_sponge_east") then
            sponge_w = "wall"
            sponge_e = "sponge"
         end if
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 500.0, slope_scale = 40000.0, wind_config = '2gyre', "// &
               "taux_magnitude = 0.05 /"//NL// &
               "&ocean_bc_nml west = '"//trim(sponge_w)//"', east = '"//trim(sponge_e)// &
               "', south = 'open', "// &
               "north = 'clamped', north_clamped_v = -0.01, sponge_width = 3, "// &
               "sponge_strength = 1.0e-4, sponge_relax_tracers = .true., "// &
               "radiation_scheme = 'orlanski', res_lscale_out = 20000.0, "// &
               "res_lscale_in = 20000.0 /"//NL
      case ("file_readers", "file_readers_zstar_full", "file_readers_zstar")
         ! The three per-rank windowed readers (supergrid, bathymetry,
         ! z-level T/S IC) on files the test writes, with the global-1-degree
         ! physics set: z_fixed + closed partial-step faces, fv_mom6 PGF,
         ! energy Coriolis, Wright EOS, periodic in x — plus Redi with a
         ! real diffusivity on the closed-face open-window path, and
         ! Fox-Kemper MLE (EPBL supplies its MLD, so EPBL replaces KPP
         ! here) on the closed-face open-column path.  The z_fixed levels
         ! are tanh-stretched (100 m at the surface) so the mixed layer
         ! spans several of the four layers: with uniform 750 m layers the
         ! ML is the top layer alone, mu(0) - mu(-1) = 0, and MLE moves
         ! nothing.
         nml = common// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3 /"//NL// &
               "&ocean_grid_nml grid_config = 'supergrid', supergrid_file = '"//SG_FILE// &
               "', coriolis_scheme = 'planetary', rad_earth = 6.371e6 /"//NL// &
               "&physics_nml wind_stress_x = 0.08, wind_stress_y = 0.0 /"//NL// &
               "&ocean_topo_nml topo_config = 'file', max_depth = 3000.0 /"//NL// &
               "&output_nml bathymetry_file = '"//BATHY_FILE//"', output_to_file = .false. /"//NL// &
               "&ocean_zinit_nml enable = .true., source = 'file', file = '"//ZINIT_FILE//"' /"//NL// &
               "&ocean_pgf_nml form = 'fv_mom6' /"//NL// &
               "&ocean_coriolis_nml form = 'sadourny_energy' /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3, hbbl = 10.0, bg_vel = 0.1 /"//NL// &
               "&ocean_redi_nml enable = .true., khtr = 600.0 /"//NL// &
               "&ocean_vmix_nml use_kpp = .false. /"//NL// &
               "&ocean_epbl_nml enable = .true. /"//NL// &
               "&ocean_foxkemper_nml enable = .true., ce = 0.08 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("sea_ice", "sea_ice_transport")
         ! Sea ice on the ocean's decomposition (thermo + EVP dynamics +
         ! the ice->ocean stress blend; category transport off): a cold,
         ! cooled re-entrant channel over a seamount under an oblique wind,
         ! partial ice cover, `dt_therm_ratio = 2` so the EVP (every step)
         ! and thermo (every other step) cadences both run, and the CFL
         ! clip + `project_ci` on (`sea_ice_transport`: the same with the
         ! category transport on, 2 advective substeps, across the periodic
         ! seam and every rank seam).  `cfl_trunc = 0.01` is deliberately
         ! tight (bound ~0.1 m/s against a ~0.17 m/s drift) so the final
         ! clip really fires: the post-clip exchange and the rank-summed
         ! truncation count are exercised on live values, not zeros.  The periodic edge
         ! puts the ice across a wrap seam on every factorisation, the
         ! interior seams across a rank seam.  Its own tracer/thermo
         ! groups: a column a few hundredths of a degree above the
         ! liquidus (T_f(34) = -1.836 degC).  The basal flux hands the ice
         ! the surface layer's whole above-freezing heat content each
         ! window, so a warmer column (the common block's 4-20 degC, or
         ! even -1.2 degC at the surface) melts every physical cell out at
         ! the first thermo window and leaves the EVP, the transport and
         ! the blend nothing to move.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
               "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
               "dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.1, wind_stress_y = 0.04 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, "// &
               "edge_depth = 800.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
               "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
               "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
               "restore_lambda = 20.0, sw_down = 0.0, dynamics = .true., "// &
               "evp_sub_steps = 30, cfl_trunc = 0.01, project_ci = .true., "// &
               "transport = "//merge(".true. ", ".false.", label == "sea_ice_transport")// &
               ", adv_substeps = 2 /"//NL// &
               "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL
      case ("sea_ice_fold", "sea_ice_fold_nodyn")
         ! Sea ice on the tripolar cap (A7, #143): the `tripolar` grid / land /
         ! fold + the `sea_ice_transport` ice block (cold column just above
         ! the liquidus, partial uniform cover, CFL clip live, `project_ci`,
         ! 2 advective substeps), under an oblique wind so the ice drifts
         ! into the fold row.  `_nodyn` switches the EVP off.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
               "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 30, ny = 24, nghost = 3, dx = 12.0, dy = 1.0 /"//NL// &
               "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, lat_south = 59.0, "// &
               "phi_join = 74.0, lon_pole = 0.0, rad_earth = 6.378e6, "// &
               "coriolis_scheme = 'planetary' /"//NL// &
               "&physics_nml wind_stress_x = 0.1, wind_stress_y = 0.04 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'double_drake', max_depth = 1000.0, "// &
               "slope_scale = 0.2 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
               "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
               "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
               "restore_lambda = 20.0, sw_down = 0.0, "// &
               ice_dyn// &
               "transport = .true., adv_substeps = 2 /"//NL// &
               "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'tripolar_fold' /"//NL
      case ("bulk_flux", "bulk_flux_ice", "bulk_flux_fold")
         ! Atmosphere-driven bulk forcing (`&ocean_atm_forcing_nml` +
         ! `&ocean_bulk_flux_nml`): nine prescribed atmospheric fields from
         ! files the test writes (spatially NON-uniform, non-symmetric in
         ! both directions), regridded per rank, turned into stress + heat +
         ! freshwater fluxes at the cell centres and averaged to faces
         ! across every rank seam.  `bulk_flux`: lon-lat spoon basin, walls.
         ! `bulk_flux_ice`: a cold re-entrant lon-lat channel with sea ice
         ! (EVP + category transport) whose open-water fraction scales the
         ! bulk fluxes and whose surface stress is blended.  `bulk_flux_fold`:
         ! the tripolar cap, so the stress crosses the north fold.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&mpi_nml px = "//trim(spx)//", py = "//trim(spy)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0, start_date = '2000-06-21' /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
               "&ocean_forcing_nml enable_components = .true. /"//NL// &
               "&ocean_bulk_flux_nml enable = .true., n_iter = 5, z_wind = 10.0, "// &
               "z_ta = 10.0 /"//NL
         if (label == "bulk_flux_ice") then
            nml = nml// &
                  "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
                  "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
                  "&ocean_atm_forcing_nml enable = .true., time_mode = 'interannual', "// &
                  "filelist_tas = '"//ATM_COLD_LIST//"', filelist_huss = '"//ATM_COLD_LIST// &
                  "', filelist_uas = '"//ATM_COLD_LIST//"', filelist_vas = '"//ATM_COLD_LIST// &
                  "', filelist_psl = '"//ATM_COLD_LIST//"', filelist_rlds = '"//ATM_COLD_LIST// &
                  "', filelist_rsds = '"//ATM_COLD_LIST//"', filelist_prra = '"//ATM_COLD_LIST// &
                  "', filelist_prsn = '"//ATM_COLD_LIST//"' /"//NL// &
                  "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
                  "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, "// &
                  "restore_lambda = 20.0, dynamics = .true., evp_sub_steps = 30, "// &
                  "cfl_trunc = 0.01, project_ci = .true., transport = .true., "// &
                  "adv_substeps = 2 /"//NL// &
                  "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL
         else
            wl = ATM_WARM_LIST
            if (label == "bulk_flux_fold") wl = ATM_CRSE_LIST
            nml = nml// &
                  "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
                  "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
                  "&ocean_atm_forcing_nml enable = .true., time_mode = 'interannual', "// &
                  "filelist_tas = '"//trim(wl)//"', filelist_huss = '"//trim(wl)// &
                  "', filelist_uas = '"//trim(wl)//"', filelist_vas = '"//trim(wl)// &
                  "', filelist_psl = '"//trim(wl)//"', filelist_rlds = '"//trim(wl)// &
                  "', filelist_rsds = '"//trim(wl)//"', filelist_prra = '"//trim(wl)// &
                  "', filelist_prsn = '"//trim(wl)//"' /"//NL
         end if
         if (label == "bulk_flux_fold") then
            ! Coarse (30 x 10 degree) atmosphere source.  The default
            ! conservative flux regrid (`fluxes_regrid` is not set) overflows
            ! its per-target-cell candidate cap (`REGRID_NMAX_CAP`) on the
            ! tripolar cap with a fine (5-degree) source, even on one rank:
            ! a serial limitation independent of the decomposition.  The
            ! 10-degree latitude spacing keeps a source centre (85N) poleward
            ! of `phi_join`, so the bilinear state regrid has a meridional
            ! gradient over the whole cap, not only in longitude.
            nml = nml// &
                  "&grid_nml nx = 30, ny = 24, nghost = 3, dx = 12.0, dy = 1.0 /"//NL// &
                  "&ocean_grid_nml grid_config = 'tripolar', lon_west = 0.0, "// &
                  "lat_south = 59.0, phi_join = 74.0, lon_pole = 0.0, "// &
                  "rad_earth = 6.378e6, coriolis_scheme = 'planetary' /"//NL// &
                  "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
                  "&ocean_topo_nml topo_config = 'double_drake', max_depth = 1000.0, "// &
                  "slope_scale = 0.2 /"//NL// &
                  "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
                  "north = 'tripolar_fold' /"//NL
         else
            nml = nml// &
                  "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
                  "dx = 1.0, dy = 1.0 /"//NL// &
                  "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, "// &
                  "lat_south = 20.0, rad_earth = 6.371e6, "// &
                  "coriolis_scheme = 'planetary' /"//NL// &
                  "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
                  "&ocean_topo_nml topo_config = 'spoon', max_depth = 3000.0, "// &
                  "edge_depth = 300.0, slope_scale = 300000.0 /"//NL
            if (label == "bulk_flux_ice") then
               nml = nml//"&ocean_bc_nml west = 'periodic', east = 'periodic', "// &
                     "south = 'wall', north = 'wall' /"//NL
            else
               nml = nml//"&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', "// &
                     "north = 'wall' /"//NL
            end if
         end if
      case default
         error stop "test_ocean_decomp_bitid_mpi: unknown case"
      end select
      ! The closed partial-step faces on both geometric coordinates.  The
      ! `zstar_full` twin's 1500 m fine zone (nz/3 = 1 layer) leaves every
      ! column shallower than that — the northern shelf and the ridge crest
      ! — as one partial cell over three fillers, so its mask (built per
      ! rank from the exchanged bathymetry) and its on-target seed both
      ! straddle the rank seams.
      if (label == "file_readers") then
         nml = nml//"&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true., "// &
               "z_fixed_profile = 'tanh', z_fixed_dz_top = 100.0 /"//NL
      else if (label == "file_readers_zstar_full") then
         nml = nml//"&vcoord_nml vcoord_type = 'zstar_full', zstar_h_surf_target = 1500.0, "// &
               "zfixed_closed_faces = .true. /"//NL
      else if (label == "file_readers_zstar") then
         ! MOM6 z*: the z_fixed staircase (the same fillers, the same mask)
         ! with every live layer dilated by the column's (H + eta)/H each
         ! regrid, and the IC seeded on the eta = 0 target.
         nml = nml//"&vcoord_nml vcoord_type = 'zstar', zfixed_closed_faces = .true. /"//NL
      end if
      if (label /= "file_readers" .and. label /= "file_readers_zstar_full" .and. &
          label /= "file_readers_zstar") then
         nml = nml//"&output_nml output_to_file = .false. /"//NL
      end if
   end function case_nml

   function cavity_nml(label, scheme, spx, spy, gfrac) result(nml)
      !! Ice-shelf cavity cases (dyn + virtual melt + quadratic top drag, all
      !! walls).  The text is `cavity_case_nml` of `test_ocean_cavity_engine`
      !! (which asserts the geometry sits ON the 2x2 seams: a grounded disc
      !! from a centred seamount, an L-shaped sloping calving front).
      !!   cavity_zfixed       z_fixed + closed faces, nz 8, PP81;
      !!   cavity_sigma        sigma, nz 4, KPP on (live under-ice u*);
      !!   cavity_sigma_ratio2 cavity_sigma with `dt_therm_ratio = 2`;
      !!   cavity_melt_only cavity_sigma with the top drag OFF (melt-only under-ice u*);
      !!   cavity_periodic cavity_sigma re-entrant in x (periodic west/east): the shelf
      !!                       box sits against the wrap and the face twins wrap;
      !!   cavity_combo       cavity_sigma with the cheap composition knobs that have no case of
      !!                       their own: EPBL (KPP off), the `hj99` exchange law, the vdiff
      !!                       `implicit_top_drag` fold and `dt_tracer_advect_ratio = 2`;
      !!   cavity_gall         a flat bed under a deep sloping draft with no box
      !!                       limit: the west 14 of 26 columns are grounded
      !!                       (global 14/26 = 0.54), so a 4x1 west tile is 100 %
      !!                       grounded; `grounded_max_frac` is `gfrac`.
      character(len=*), intent(in) :: label, scheme, spx, spy
      real(wp), intent(in), optional :: gfrac
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: vc, nzl, cav, vm, topo, draft, tdrag
      character(len=24) :: sfrac
      character(len=8) :: xedge
      character(len=48) :: vdf, xlaw

      vm = ""
      vdf = ""
      xlaw = "const_gamma"
      xedge = "wall"
      if (label == "cavity_periodic") xedge = "periodic"
      tdrag = "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"//NL
      if (label == "cavity_melt_only") tdrag = ""
      if (label == "cavity_zfixed") then
         nzl = "8"
         vc = "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true. /"//NL
         cav = ", h_min_cavity = 250.0 /"//NL
         vm = "&ocean_vmix_nml use_kpp = .false., pp81_nu0 = 0.0, pp81_nu_bg = 1.0e-3, "// &
              "pp81_kappa_bg = 5.0e-5 /"//NL
      else
         nzl = "4"
         vc = "&vcoord_nml vcoord_type = 'sigma' /"//NL
         cav = ", h_min_cavity = 100.0 /"//NL
         if (label == "cavity_sigma_ratio2") vm = "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL
         if (label == "cavity_combo") then
            vm = "&ocean_vmix_nml use_kpp = .false., dt_therm_ratio = 2, dt_tracer_advect_ratio = 2 /"//NL// &
                 "&ocean_epbl_nml enable = .true. /"//NL
            vdf = "&ocean_vdiff_nml implicit_top_drag = .true. /"//NL
            xlaw = "hj99"
         end if
      end if
      topo = "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 400.0, "// &
             "slope_scale = 80000.0, wind_config = 'constant', taux_magnitude = 0.05 /"//NL
      draft = "draft_depth = 500.0, draft_slope = -4.0e-4, draft_x0 = 0.0, "// &
              "draft_x1 = 380000.0, draft_y1 = 260000.0"
      if (index(label, "cavity_gall") == 1) then
         topo = "&ocean_topo_nml topo_config = 'flat', max_depth = 1000.0, "// &
                "wind_config = 'constant', taux_magnitude = 0.05 /"//NL
         draft = "draft_depth = 1400.0, draft_slope = -2.2e-3, draft_x0 = 0.0"
         sfrac = "0.6"
         if (present(gfrac)) write (sfrac, '(f6.3)') gfrac
         cav = ", h_min_cavity = 100.0, grounded_max_frac = "//trim(adjustl(sfrac))//" /"//NL
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
            vc//vm//topo// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            draft//cav// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = '"//trim(xlaw)//"', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0 /"//NL// &
            tdrag//vdf// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//NL// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//NL// &
            "&ocean_bc_nml west = '"//xedge//"', east = '"//xedge//"', south = 'wall', "// &
            "north = 'wall' /"//NL
   end function cavity_nml

   subroutine run_one(nml, csize, crank, snap, ok)
      !! Configure, step N_STEPS, snapshot every registry field to host.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: csize, crank
      type(snap_t), intent(out) :: snap
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(console_stats_t) :: cstats
      type(comm_t) :: real_comm, self_comm
      logical :: pushed_self_comm
      integer :: ierr, n
      real(wp) :: t

      ok = .false.
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr, compute_rank=crank, compute_size=csize)
      if (ierr /= OCEAN_STATUS_OK) return

      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      ierr = OCEAN_STATUS_OK
      do n = 1, n_steps_cur
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         ! Driver order: the sea-ice block (a no-op unless `&ocean_ice_nml
         ! enable`) between the dyn-core advance and the finalize.
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         ! Closed-budget mass/salt/heat out+src at the final step, EFP-reduced
         ! (reproducing_sums default on) so it is order-invariant across rank
         ! counts -- the periodic-seam / sponge-seam-ghost regression this test
         ! case exists for (see `compare_budget` and
         ! `check_periodic_sponge_serial_out`).  A fresh `cstats` per call is
         ! fine: `budget_out` is computed independently of the t=0 reference
         ! `cstats` latches internally.
         !
         ! HAZARD (found in review, fixed by the comm_env push/pop seam):
         ! `ocean_console_stats_report` takes no `compute_size` -- it is
         ! COLLECTIVE over whatever `comm_env_compute_comm()` returns, the
         ! REAL job communicator, regardless of the logical `csize`/`crank`
         ! the caller is emulating.  The serial REFERENCE call below always
         ! passes `csize=1, crank=0` to `engine_setup` so every physical rank
         ! independently computes a full, non-decomposed copy of the SAME
         ! problem -- calling the collective report unmodified from inside
         ! that call would allreduce across every physical rank's
         ! independent copy (an `nprocs`-fold bogus sum for the
         ! mass/salt/heat totals feeding `bud`), not the single logical rank
         ! the reference is standing in for.
         !
         ! Fix: when `csize == 1` (the serial-reference call, on ANY real
         ! `nprocs`), push a per-rank SELF-communicator
         ! (`comm_env_push_compute_comm`, `src/comm/rdb_comm_env.F90`) around
         ! the report so its collective becomes an identity op over each
         ! rank's own independent copy instead of allreducing `nprocs`
         ! copies together, then pop it back off.  When `csize == nprocs`
         ! (a real decomposed `px x py` run), the real communicator already
         ! IS the right scope and no override is needed.  Either way
         ! `snap%bud` now comes out correctly populated (`*_active=.true.`)
         ! on every launch size, which is what lets
         ! `check_periodic_sponge_serial_out` compare the reference against
         ! every decomposed factorization bit-for-bit at nprocs=1, 2 AND 4
         ! (previously it could only run the old, thresholded check at
         ! nprocs==1, where the reference's own trivial `csize=1` collective
         ! happened to already be correct).
         pushed_self_comm = .false.
         if (csize == 1) then
            real_comm = comm_env_compute_comm()
            self_comm = real_comm%split_by(real_comm%rank())
            call comm_env_push_compute_comm(self_comm)
            pushed_self_comm = .true.
         else if (csize /= nprocs) then
            error stop "run_one: csize must be 1 (serial reference) or nprocs (decomposed run)"
         end if
         call ocean_console_stats_report(cstats, engine%grid, engine%state%metrics, &
                                         engine%state%multilayer, t, DT, max(n_steps_cur, 1), &
                                         compute_rank=crank, reproducing_sums=.true., &
                                         budget_stage_weight=ocean_budget_stage_weight( &
                                         engine%state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR), &
                                         budget_out=snap%bud)
         if (pushed_self_comm) then
            call comm_env_pop_compute_comm()
            call self_comm%finalize()
         end if
         call take_snapshot(engine, snap)
         snap%drain_iter = engine%state%continuity%drain_iter_max
         ok = .true.
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_one

   subroutine push_field(s, tag, a)
      !! Append one 2-D field (stored with a unit third extent).
      type(snap_t), intent(inout) :: s
      character(len=*), intent(in) :: tag
      real(wp), intent(in) :: a(:, :)
      s%n = s%n + 1
      if (s%n > MAXF) error stop "test_ocean_decomp_bitid_mpi: raise MAXF"
      s%f(s%n)%tag = tag
      allocate (s%f(s%n)%a(size(a, 1), size(a, 2), 1))
      s%f(s%n)%a(:, :, 1) = a
   end subroutine push_field

   subroutine add_cavity_fields(engine, s)
      !! The ice-shelf cavity's own outputs and static face twins, which the
      !! restart registry does not carry: the slot arrays (melt rate, friction
      !! velocity, interface temperature, under-ice stress, drag stress),
      !! compared on owned cells like every other field, and the four face
      !! twins `cover_u/v`, `k_top_u/v` (tag `img:`), additionally compared over
      !! the SEAM GHOSTS by `check_ghost_images`.  Also fills the non-vacuity
      !! counts.  Component arrays only are pulled from the device.
      type(ocean_engine_t), intent(inout), target :: engine
      type(snap_t), intent(inout) :: s
      integer :: i, j, ng, nz

      !$acc update self(engine%state%cavity_flux%melt, engine%state%cavity_flux%ustar, engine%state%cavity_flux%t_b)
      !$acc update self(engine%state%surface_stress%stress_shelf, engine%state%tdrag%stress_top)
      !$acc update self(engine%state%tdrag%cover_u, engine%state%tdrag%cover_v)
      !$acc update self(engine%state%multilayer%wet_mask)
      call push_field(s, "cav:melt", engine%state%cavity_flux%melt)
      call push_field(s, "cav:ustar", engine%state%cavity_flux%ustar)
      call push_field(s, "cav:t_b", engine%state%cavity_flux%t_b)
      call push_field(s, "cav:stress_shelf", engine%state%surface_stress%stress_shelf)
      call push_field(s, "cav:stress_top", engine%state%tdrag%stress_top)
      call push_field(s, "img:cover_u", engine%state%tdrag%cover_u)
      call push_field(s, "img:cover_v", engine%state%tdrag%cover_v)
      call push_field(s, "img:k_top_u", real(engine%state%multilayer%k_top_u, wp))
      call push_field(s, "img:k_top_v", real(engine%state%multilayer%k_top_v, wp))

      ng = engine%grid%nghost
      nz = engine%state%multilayer%nz_ml
      s%ncav = 0
      do j = ng + 1, ng + engine%grid%ny_phys
         do i = ng + 1, ng + engine%grid%nx_phys
            associate (wet => engine%state%multilayer%wet_mask(i, j) > 0.5_wp, &
                       cov => engine%state%metrics%cover_frac(i, j) > 0.0_wp)
               if (wet .and. cov) s%ncav(1) = s%ncav(1) + 1
               if (wet .and. .not. cov) s%ncav(2) = s%ncav(2) + 1
               if (.not. wet .and. engine%state%metrics%z_draft(i, j) > 0.0_wp) &
                  s%ncav(3) = s%ncav(3) + 1
               if (wet .and. cov .and. engine%state%cavity_flux%melt(i, j) /= 0.0_wp) &
                  s%ncav(4) = s%ncav(4) + 1
               if (wet .and. cov .and. engine%state%tdrag%stress_top(i, j) > 0.0_wp) &
                  s%ncav(5) = s%ncav(5) + 1
               if (engine%state%multilayer%k_top(i, j) < nz) s%ncav(6) = s%ncav(6) + 1
            end associate
         end do
      end do
   end subroutine add_cavity_fields

   subroutine take_snapshot(engine, s)
      !! Host copy of every restart-registry field (device -> host first)
      !! plus the BT end-of-step eta.
      type(ocean_engine_t), intent(inout), target :: engine
      type(snap_t), intent(inout) :: s
      type(restart_registry_t) :: reg
      integer :: e

      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      s%n = 0
      do e = 1, reg%n
         associate (en => reg%entries(e))
            if (en%rank == 0) cycle     ! rank-local host scalars (Chapman corners)
            if (en%device_mapped) then
               if (en%rank == 2) then
                  !$acc update self(en%p2)
               else
                  !$acc update self(en%p3)
               end if
            end if
            s%n = s%n + 1
            if (s%n > MAXF) error stop "test_ocean_decomp_bitid_mpi: raise MAXF"
            s%f(s%n)%tag = en%tag
            if (en%rank == 2) then
               allocate (s%f(s%n)%a(size(en%p2, 1), size(en%p2, 2), 1))
               s%f(s%n)%a(:, :, 1) = en%p2
            else
               s%f(s%n)%a = en%p3
            end if
         end associate
      end do
      !$acc update self(engine%state%dyn%bt_work%bt_eta_end)
      s%n = s%n + 1
      s%f(s%n)%tag = "bt_eta_end"
      allocate (s%f(s%n)%a(size(engine%state%dyn%bt_work%bt_eta_end, 1), &
                           size(engine%state%dyn%bt_work%bt_eta_end, 2), 1))
      s%f(s%n)%a(:, :, 1) = engine%state%dyn%bt_work%bt_eta_end
      if (engine%state%metrics%use_cavity) call add_cavity_fields(engine, s)
      s%io = engine%grid%i_offset_global
      s%jo = engine%grid%j_offset_global
      s%nxl = engine%grid%nx_phys
      s%nyl = engine%grid%ny_phys
      s%ng = engine%grid%nghost
      s%ngrow = ocean_fold_exchange_grow_count()
   end subroutine take_snapshot

#ifndef RDB_NO_NETCDF
   subroutine check_resume_leg(nml, straight, tag)
      !! The write/resume leg: the same decomposed run, checkpointed at
      !! N_CKPT to per-rank files (the driver's `<dir>/restart_rank_NNNNNN.nc`
      !! convention, through `ocean_state_restart_write`), torn down, and
      !! resumed through `engine_setup(restart_file=<dir>)` -- the warm
      !! restart the driver runs -- for the remaining N_STEPS - N_CKPT
      !! steps.  Every registry field (and the BT end-of-step eta) must
      !! equal the STRAIGHT decomposed run's bitwise over the FULL local
      !! array, ghosts included.  Counts into `n_fail` (rank-summed).
      character(len=*), intent(in) :: nml, tag
      type(snap_t), intent(in) :: straight
      type(snap_t) :: res
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      character(len=512) :: fname
      integer :: ierr, step0, glob(2), nb
      real(wp) :: t
      logical :: ok

      ok = .false.
      nb = 0
      if (rank == 0) call ensure_directory_exists(RST_DIR)
      call comm%barrier()
      fname = output_rank_filename(RST_DIR, "restart", rank)

      ! ---- leg 1: N_CKPT steps, checkpoint ----
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr == OCEAN_STATUS_OK) call validate_config(cfg, ierr)
      if (ierr == OCEAN_STATUS_OK) &
         call engine_setup(engine, cfg, ierr, compute_rank=rank, compute_size=nprocs)
      if (ierr == OCEAN_STATUS_OK) then
         call engine_enter_data(engine, cfg)
         t = 0.0_wp
         call advance(engine, cfg, t, N_CKPT, ierr)
         if (ierr == OCEAN_STATUS_OK) &
            call ocean_state_restart_write(engine%state, engine%grid, engine%decomp, &
                                           trim(fname), t, N_CKPT, ierr=ierr)
         call engine_exit_data(engine)
         call engine_teardown(engine)
      end if
      call comm%barrier()

      ! ---- leg 2: warm restart from the per-rank files, finish ----
      if (ierr == OCEAN_STATUS_OK) then
         call read_config_from_string(nml, cfg, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call validate_config(cfg, ierr)
         if (ierr == OCEAN_STATUS_OK) &
            call engine_setup(engine, cfg, ierr, compute_rank=rank, compute_size=nprocs, &
                              restart_file=RST_DIR, t_restart=t, step_restart=step0)
         if (ierr == OCEAN_STATUS_OK) then
            if (step0 == N_CKPT .and. engine%warm_restart) then
               call engine_enter_data(engine, cfg)
               call advance(engine, cfg, t, N_STEPS - N_CKPT, ierr)
               if (ierr == OCEAN_STATUS_OK) then
                  call take_snapshot(engine, res)
                  ok = .true.
               end if
               call engine_exit_data(engine)
            end if
            call engine_teardown(engine)
         end if
      end if
      call delete_file(trim(fname))

      if (ok) call compare_full(res, straight, nb, tag)
      glob = [nb, merge(0, 1, ok)]
      call allreduce(comm, glob, op=MPI_SUM)
      if (rank == 0) then
         if (sum(glob) == 0) then
            write (*, '(3a,i0,a)') "case ", tag, " write/resume: IDENTICAL (", &
               straight%n, " fields, ghosts included)"
         else
            write (*, '(3a,i0,a,i0)') "FAIL ", tag, " write/resume: mismatches=", glob(1), &
               " leg_failures=", glob(2)
         end if
      end if
      if (sum(glob) > 0) n_fail = n_fail + 1
   end subroutine check_resume_leg

   subroutine advance(engine, cfg, t, nsteps, ierr)
      !! The driver's step sequence, `nsteps` times from `t`.
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      integer, intent(out) :: ierr
      integer :: n
      ierr = OCEAN_STATUS_OK
      do n = 1, nsteps
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + DT
      end do
   end subroutine advance

   subroutine compare_full(a, b, nbad, report)
      !! Bitwise, every element of every field (ghosts included) -- the
      !! two runs share a decomposition.
      type(snap_t), intent(in) :: a, b
      integer, intent(out) :: nbad
      character(len=*), intent(in) :: report
      integer :: e, nb
      nbad = 0
      if (a%n /= b%n) then
         write (*, '(a,i0,a,i0,a,i0)') "  rank ", rank, ": field count differs ", a%n, " vs ", b%n
         nbad = 1
         return
      end if
      do e = 1, a%n
         if (a%f(e)%tag /= b%f(e)%tag .or. any(shape(a%f(e)%a) /= shape(b%f(e)%a))) then
            nbad = nbad + 1
            cycle
         end if
         nb = count(transfer(a%f(e)%a, 0_int64, size(a%f(e)%a)) /= &
                    transfer(b%f(e)%a, 0_int64, size(b%f(e)%a)))
         if (nb > 0) write (*, '(a,i0,5a,i0,a,es10.3)') "  rank ", rank, " ", report, &
            " write/resume ", trim(a%f(e)%tag), ": mismatches=", nb, &
            "  max|diff|=", maxval(abs(a%f(e)%a - b%f(e)%a))
         nbad = nbad + nb
      end do
   end subroutine compare_full

   subroutine delete_file(fname)
      character(len=*), intent(in) :: fname
      integer :: u, ios
      logical :: exists
      inquire (file=fname, exist=exists)
      if (.not. exists) return
      open (newunit=u, file=fname, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine delete_file
#endif

   subroutine compare(dec, ref, nbad, nfin, report)
      !! Bitwise comparison of each field's owned window.  A face array is
      !! one wider than the cell array along its stagger; its owned window
      !! then includes both edge faces of the tile.
      type(snap_t), intent(in) :: dec, ref
      integer, intent(out) :: nbad, nfin
      character(len=*), intent(in) :: report
      integer :: e, i, j, k, ex, ey, nb
      character(len=160) :: first

      nbad = 0
      nfin = 0
      if (dec%n /= ref%n) then
         write (*, '(a,i0,a,i0,a,i0)') "  rank ", rank, ": field count differs ", dec%n, " vs ", ref%n
         nbad = 1
         return
      end if
      do e = 1, dec%n
         if (dec%f(e)%tag /= ref%f(e)%tag) then
            write (*, '(a,i0,4a)') "  rank ", rank, ": field order differs ", &
               trim(dec%f(e)%tag), " vs ", trim(ref%f(e)%tag)
            nbad = nbad + 1
            cycle
         end if
         if (carried_tendency(dec%f(e)%tag)) cycle
         associate (a => dec%f(e)%a, b => ref%f(e)%a)
            ex = size(a, 1) - (dec%nxl + 2*dec%ng)
            ey = size(a, 2) - (dec%nyl + 2*dec%ng)
            nb = 0
            first = ""
            do k = 1, size(a, 3)
               do j = dec%ng + 1, dec%ng + dec%nyl + ey
                  do i = dec%ng + 1, dec%ng + dec%nxl + ex
                     if (.not. ieee_is_finite(b(i + dec%io, j + dec%jo, k))) nfin = nfin + 1
                     if (transfer(a(i, j, k), 0_int64) /= &
                         transfer(b(i + dec%io, j + dec%jo, k), 0_int64)) then
                        if (nb == 0) write (first, '(a,3(i0,1x),a,es24.16,a,es24.16)') &
                           "first at local (i,j,k)=", i, j, k, ": ", a(i, j, k), " vs ", &
                           b(i + dec%io, j + dec%jo, k)
                        nb = nb + 1
                     end if
                  end do
               end do
            end do
         end associate
         if (nb > 0) then
            write (*, '(a,i0,5a,i0,2a)') "  rank ", rank, " ", report, " ", &
               trim(dec%f(e)%tag), ": mismatches=", nb, "  ", trim(first)
         end if
         nbad = nbad + nb
      end do
   end subroutine compare

   pure logical function is_ice_fold_case(label)
      !! The sea-ice-on-the-tripolar-cap cases (A7).
      character(len=*), intent(in) :: label
      is_ice_fold_case = label == "sea_ice_fold" .or. label == "sea_ice_fold_nodyn"
   end function is_ice_fold_case

   subroutine check_ice_seam(ref, tag)
      !! On the SERIAL reference of an ice-on-fold case: (1) ice must have
      !! reached the seam (nonzero `ice_m_ice` in the last physical row, and
      !! under EVP nonzero `ice_v_ice` in the top rows), else bit-identity of
      !! the folds would be vacuous; (2) the `ice_v_ice` fold row must be
      !! exactly antisymmetric under the 180-degree cap rotation
      !! (`v(i) + v(2*ng + nx + 1 - i) == 0`, Murray 1996 bipolar fold).
      !! Bitwise equality with the reference then carries both properties to
      !! every rank layout.  Failures count into `n_fail`.
      type(snap_t), intent(in) :: ref
      character(len=*), intent(in) :: tag
      integer :: e, i, im, jf, nbad_anti
      real(wp) :: seam_mice, seam_v
      logical :: have_v

      seam_mice = 0.0_wp
      seam_v = 0.0_wp
      nbad_anti = 0
      have_v = .false.
      jf = ref%ng + ref%nyl + 1
      do e = 1, ref%n
         if (trim(ref%f(e)%tag) == "ice_m_ice") then
            seam_mice = maxval(abs(ref%f(e)%a(ref%ng + 1:ref%ng + ref%nxl, &
                                              ref%ng + ref%nyl, :)))
         else if (trim(ref%f(e)%tag) == "ice_v_ice") then
            have_v = .true.
            seam_v = maxval(abs(ref%f(e)%a(ref%ng + 1:ref%ng + ref%nxl, jf - 3:jf, 1)))
            do i = ref%ng + 1, ref%ng + ref%nxl
               im = 2*ref%ng + ref%nxl + 1 - i
               if (ref%f(e)%a(i, jf, 1) + ref%f(e)%a(im, jf, 1) /= 0.0_wp) nbad_anti = nbad_anti + 1
            end do
         end if
      end do
      if (.not. (seam_mice > 0.0_wp)) then
         write (*, '(3a,i0)') "FAIL ", tag, ": no ice reaches the fold seam, rank ", rank
         n_fail = n_fail + 1
      end if
      if (have_v .and. .not. (seam_v > 0.0_wp)) then
         write (*, '(3a,i0)') "FAIL ", tag, ": v_ice is zero at the fold seam, rank ", rank
         n_fail = n_fail + 1
      end if
      if (nbad_anti > 0) then
         write (*, '(3a,i0,a,i0)') "FAIL ", tag, ": v_ice fold row not antisymmetric, cells=", &
            nbad_anti, " rank ", rank
         n_fail = n_fail + 1
      end if
   end subroutine check_ice_seam

   pure logical function carried_tendency(tag)
      !! The `pred_corr` predictor's carried viscous tendency
      !! (`hvisc_du_visc` / `hvisc_dv_visc`) is restart state -- the
      !! predictor reuses it -- but it is not compared here.  It is not
      !! decomposition-invariant to the last bit on a CPU build: nvfortran
      !! vectorises the hvisc kernels, whether a face lands in the vector
      !! body or the scalar remainder (with a different FMA contraction)
      !! depends on its position in the TILE, and the tendency is a
      !! difference of fluxes, so the reordering survives cancellation.
      !! Measured on nvfortran 26.5 CPU + HPC-X, 4x1 and 2x2, with FMA
      !! contraction ON: up to a relative 3e-15 at the faces at and beside
      !! a tile seam; the GPU build's device kernels have no remainder
      !! loop and are bitwise.  It CAN reach the prognostics: under
      !! z* + visc_rem (`visc_rem_zstar`) the FMA difference drifted eta
      !! 1e-8 in 24 steps on 4x1.  CPU builds therefore compile with
      !! contraction off, globally (`-Mnofma` / `-ffp-contract=off`,
      !! cmake/compiler_flags.cmake), under which every prognostic field
      !! compared here is BITWISE; that is what catches a real
      !! decomposition bug in the viscosity.  This carried array itself
      !! stays excluded from the compare (whether its seam-face residual
      !! vanishes with contraction off was not re-measured).  Restarts resume on the same decomposition and
      !! restore these arrays verbatim.
      character(len=*), intent(in) :: tag
      carried_tendency = trim(tag) == "hvisc_du_visc" .or. trim(tag) == "hvisc_dv_visc"
   end function carried_tendency

   subroutine run_case(label, scheme)
      character(len=*), intent(in) :: label, scheme
      type(snap_t) :: ref, dec
      logical :: ok_ref, ok_dec
      integer :: px, nbad, nfin, glob(3), nbud, nb_img
      character(len=96) :: tag
      type(conservation_budget_t) :: buds(8)

      n_steps_cur = N_STEPS
      if (len(label) > 3) then
         if (label(len(label) - 2:) == "_t0") n_steps_cur = 0
      end if
      call run_one(case_nml(label, scheme, 1, 1), 1, 0, ref, ok_ref)
      if (.not. ok_ref) then
         write (*, '(5a,i0)') "FAIL ", label, "/", scheme, ": serial reference failed on rank ", rank
         n_fail = n_fail + 1
         return
      end if

      if (is_ice_fold_case(label)) call check_ice_seam(ref, trim(label)//"/"//scheme)
      if (index(label, "cavity_") == 1) call check_cavity_counters(ref, label, trim(label)//"/"//scheme)

      nbud = 0
      do px = nprocs, 1, -1
         if (mod(nprocs, px) /= 0) cycle
         write (tag, '(4a,i0,a,i0)') label, "/", scheme, " ", px, "x", nprocs/px
         call run_one(case_nml(label, scheme, px, nprocs/px), nprocs, rank, dec, ok_dec)
         if (.not. ok_dec) then
            write (*, '(3a,i0)') "FAIL ", trim(tag), ": decomposed run failed on rank ", rank
            glob = [1, 0, 1]
         else
            call compare(dec, ref, nbad, nfin, trim(tag))
            if (index(label, "cavity_") == 1) then
               call check_ghost_images(dec, ref, nb_img, trim(tag))
               nbad = nbad + nb_img
            end if
            ! Windowed-drain cases: the sub-cycle must run >1 pass (else the
            ! pass-count reduction is never exercised) and every layout must
            ! have used the same count as the serial reference.
            ! (`cavity_combo` waives ONLY the ">1 pass" requirement: its quiescent
            ! flow needs a single pass per window.  It must still have run the
            ! drain at all (ref >= 1) with the same count on every layout; the
            ! windowed accumulation is compared bitwise.)
            if (index(case_nml(label, scheme, 1, 1), "dt_tracer_advect_ratio") > 0) then
               if (ref%drain_iter < merge(1, 2, label == "cavity_combo") .or. &
                   dec%drain_iter /= ref%drain_iter) then
                  write (*, '(3a,i0,a,i0)') "FAIL ", trim(tag), ": drain passes ref=", &
                     ref%drain_iter, " dec=", dec%drain_iter
                  nbad = nbad + 1
               else if (rank == 0) then
                  write (*, '(3a,i0,a)') "case ", trim(tag), ": drain max_iter=", ref%drain_iter, " (ok)"
               end if
            end if
            ! A7: the fold-exchange reserve must cover every ice group.
            if (is_ice_fold_case(label) .or. index(label, "drain_ratio_fold") == 1) &
               nbad = nbad + dec%ngrow
            glob = [nbad, nfin, 0]
#ifndef RDB_NO_NETCDF
            ! Cavity: only `cavity_sigma_ratio2` (melt + top drag, checkpoint
            ! mid thermo window).  The ratio-1 cavity cases (`cavity_zfixed`,
            ! `cavity_sigma`, `cavity_gall`, `cavity_melt_only`) are NOT in the
            ! list: measured on ONE rank (np 1, same failure decomposed), their
            ! resume differs from the straight run (tracers ~1e-3, u/v ~3e-8,
            ! `sf_heat_cavity` ~9e-2) while the same namelist with the melt off
            ! (dyn + draft only) resumes bit-identically and the ratio-2 twin
            ! does too: a serial warm-restart defect of the melt path at
            ! `dt_therm_ratio = 1`, not a decomposition defect.
            if ((label == "sea_ice" .or. label == "sea_ice_transport" .or. &
                 label == "bulk_flux_ice" .or. label == "dtbt_reset" .or. &
                 label == "cavity_sigma_ratio2" .or. is_ice_fold_case(label)) .and. &
                (nprocs/px == 1 .or. px == nprocs/px) .and. &
                (nprocs > 1 .or. label == "cavity_sigma_ratio2")) then
               ! 4x1 and 2x2 at np = 4 (and 2x1 at np = 2).
               call check_resume_leg(case_nml(label, scheme, px, nprocs/px), dec, &
                                     trim(tag))
            end if
#endif
            if (nbud < size(buds)) then
               nbud = nbud + 1
               buds(nbud) = dec%bud
            end if
         end if
         call allreduce(comm, glob, op=MPI_SUM)
         if (rank == 0) then
            if (glob(1) == 0 .and. glob(3) == 0 .and. glob(2) == 0) then
               write (*, '(3a,i0,a)') "case ", trim(tag), ": IDENTICAL (", ref%n, " fields)"
            else
               write (*, '(3a,i0,a,i0,a,i0)') "FAIL ", trim(tag), ": mismatches=", glob(1), &
                  " nonfinite_ref=", glob(2), " run_failures=", glob(3)
            end if
         end if
         if (sum(glob) > 0) n_fail = n_fail + 1
      end do
      ! The real regression gate: the single-rank reference's budget
      ! against every decomposed factorization, bit-for-bit, at whatever
      ! `nprocs` this launch is (1, 2, or 4 under ctest) -- see the
      ! subroutine's own docstring.
      call check_periodic_sponge_serial_out(label, scheme, ref%bud, buds(1:nbud))
      ! Every collected `dec%bud` above has `csize == nprocs` (always true
      ! for the decomposed loop), so its collective was safe regardless of
      ! `nprocs`; comparing them to each other is meaningful once nprocs > 1
      ! (>= 2 distinct factorizations).
      if (nbud >= 2) call compare_budget(buds(1:nbud), label, scheme)
   end subroutine run_case

   subroutine check_ghost_images(dec, ref, nbad, tag)
      !! Seam-ghost image check for the `img:` face twins: EVERY element of the
      !! tile's local array (ghost ring and the outermost face included) whose
      !! global index lies inside the serial reference must equal the reference
      !! there, bitwise.  Indices outside the reference (beyond a physical edge)
      !! have no image and are skipped.
      type(snap_t), intent(in) :: dec, ref
      integer, intent(out) :: nbad
      character(len=*), intent(in) :: tag
      integer :: e, i, j, gi, gj, nb, nchecked
      character(len=160) :: first

      nbad = 0
      nchecked = 0
      do e = 1, dec%n
         if (dec%f(e)%tag(1:4) /= "img:") cycle
         nb = 0
         first = ""
         do j = 1, size(dec%f(e)%a, 2)
            do i = 1, size(dec%f(e)%a, 1)
               gi = i + dec%io
               gj = j + dec%jo
               if (gi < 1 .or. gj < 1 .or. gi > size(ref%f(e)%a, 1) .or. &
                   gj > size(ref%f(e)%a, 2)) cycle
               nchecked = nchecked + 1
               if (transfer(dec%f(e)%a(i, j, 1), 0_int64) /= &
                   transfer(ref%f(e)%a(gi, gj, 1), 0_int64)) then
                  if (nb == 0) write (first, '(a,2(i0,1x),a,es12.4,a,es12.4)') &
                     "first at local (i,j)=", i, j, ": ", dec%f(e)%a(i, j, 1), " vs ", &
                     ref%f(e)%a(gi, gj, 1)
                  nb = nb + 1
               end if
            end do
         end do
         if (nb > 0) write (*, '(a,i0,5a,i0,2a)') "  rank ", rank, " ", tag, " ghost image ", &
            trim(dec%f(e)%tag), ": mismatches=", nb, "  ", trim(first)
         nbad = nbad + nb
      end do
      if (nchecked == 0) then
         write (*, '(3a,i0)') "FAIL ", tag, ": ghost-image check compared nothing on rank ", rank
         nbad = nbad + 1
      end if
   end subroutine check_ghost_images

   subroutine check_cavity_counters(ref, label, tag)
      !! Non-vacuity of a cavity case on the serial reference: every count the
      !! case relies on must be non-zero (a geometry that drifts off the cavity
      !! would otherwise pass trivially).  `cavity_gall` is covered end to end,
      !! so it has no open water.
      type(snap_t), intent(in) :: ref
      character(len=*), intent(in) :: label, tag
      integer :: c
      logical :: need(6)

      need = .true.
      if (index(label, "cavity_gall") == 1) need(2) = .false.
      if (label /= "cavity_zfixed") need(6) = .false.
      if (label == "cavity_melt_only") need(5) = .false.
      do c = 1, 6
         if (need(c) .and. ref%ncav(c) == 0) then
            write (*, '(3a,i0,a)') "FAIL ", tag, ": cavity non-vacuity count #", c, " is 0"
            n_fail = n_fail + 1
         end if
      end do
      if (rank == 0) write (*, '(3a,6(1x,i0))') "case ", tag, ": cavity counts (covered, open, grounded, "// &
         "melting, dragged, k_top<nz) =", ref%ncav
   end subroutine check_cavity_counters

   subroutine check_cavity_grounded_refusal()
      !! A draft that grounds more of the domain than `grounded_max_frac` must
      !! be REFUSED, identically, on every rank of every layout (and on the
      !! serial reference): the verdict is taken on the GLOBAL integer count.
      !! Before it was, each tile compared its own fraction, so one rank
      !! refused while the others carried on and the job hung.  The probe has
      !! global fraction 14/26 = 0.54 against a limit of 0.4 (the accepted
      !! twin, `cavity_gall`, uses 0.6).  The run is guarded by the
      !! explicit ctest TIMEOUT on this test: a hang is a failure.
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: px, ierr, nbad, glob(1), im
      logical :: why_found
      character(len=96) :: tag

      do px = nprocs, 1, -1
         if (mod(nprocs, px) /= 0) cycle
         nbad = 0
         write (tag, '(a,i0,a,i0)') "cavity_grounded_refused ", px, "x", nprocs/px
         nml = cavity_nml("cavity_gall", "pred_corr", int_str(px), int_str(nprocs/px), 0.4_wp)// &
               "&output_nml output_to_file = .false. /"//new_line("a")
         call read_config_from_string(nml, cfg, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call validate_config(cfg, ierr)
         if (ierr /= OCEAN_STATUS_OK) then
            write (*, '(3a,i0)') "FAIL ", trim(tag), ": validate_config refused the probe, rank ", rank
            nbad = 1
         else
            call error_ring_clear()
            call engine_setup(engine, cfg, ierr, compute_rank=rank, compute_size=nprocs)
            call engine_teardown(engine)
            if (ierr == OCEAN_STATUS_OK) then
               write (*, '(3a,i0)') "FAIL ", trim(tag), ": ACCEPTED, rank ", rank
               nbad = 1
            else
               ! Refused is not enough: it must be the grounded-fraction
               ! refusal ("the prescribed draft grounds N of M ..."), on this rank.
               why_found = .false.
               do im = 0, error_ring_count() - 1
                  if (index(error_ring_get(im), "draft grounds") > 0) why_found = .true.
               end do
               if (.not. why_found) then
                  write (*, '(3a,i0)') "FAIL ", trim(tag), &
                     ": refused, but NOT for the grounded-fraction reason, rank ", rank
                  nbad = 1
               end if
            end if
         end if
         glob(1) = nbad
         call allreduce(comm, glob, op=MPI_SUM)
         if (rank == 0 .and. glob(1) == 0) &
            write (*, '(3a)') "case ", trim(tag), ": refused on every rank (ok)"
         if (glob(1) > 0) n_fail = n_fail + 1
      end do
   end subroutine check_cavity_grounded_refusal

   function int_str(n) result(str)
      integer, intent(in) :: n
      character(len=16) :: str
      write (str, '(i0)') n
   end function int_str

   subroutine compare_budget(buds, label, scheme)
      !! Bit-identity check across every decomposed-factorization budget of
      !! one case/scheme collected by `run_case` (2x1 vs 1x2 for nprocs=2;
      !! 4x1 vs 2x2 vs 1x4 for nprocs=4) -- every entry has `csize ==
      !! nprocs` so its `ocean_console_stats_report` collective was taken
      !! over the REAL communicator correctly (see the hazard note in
      !! `run_one`); no unsafe reference-vs-decomposed comparison here.
      !! The EFP reproducing-sums path (default on) makes the closed
      !! mass/salt/heat out+src totals order-invariant across rank counts,
      !! so any two factorizations of the SAME real launch must land on the
      !! identical bit pattern.
      !!
      !! NOTE what this does and does not catch: every valid `px x py`
      !! factorization of nprocs >= 2 has `px > 1 .or. py > 1`, so the OLD
      !! buggy gate (`ocean_halo_is_decomposed_x() .or. ..._y()`) already
      !! fired the seam-ghost refresh for EVERY entry compared here, both
      !! before and after the fix -- a mismatch here is a regression in
      !! the general refresh mechanism, not a re-detection of the original
      !! single-rank bug.  Only a genuinely undecomposed run (the serial
      !! reference, compared against these decomposed factorizations by
      !! `check_periodic_sponge_serial_out` below) exercises the code path
      !! the fix changed.
      type(conservation_budget_t), intent(in) :: buds(:)
      character(len=*), intent(in) :: label, scheme
      integer :: k, nbad

      nbad = 0
      do k = 2, size(buds)
         if (transfer(buds(k)%mass_out, 0_int64) /= transfer(buds(1)%mass_out, 0_int64) .or. &
             transfer(buds(k)%salt_out, 0_int64) /= transfer(buds(1)%salt_out, 0_int64) .or. &
             transfer(buds(k)%heat_out, 0_int64) /= transfer(buds(1)%heat_out, 0_int64) .or. &
             transfer(buds(k)%mass_src, 0_int64) /= transfer(buds(1)%mass_src, 0_int64) .or. &
             transfer(buds(k)%salt_src, 0_int64) /= transfer(buds(1)%salt_src, 0_int64) .or. &
             transfer(buds(k)%heat_src, 0_int64) /= transfer(buds(1)%heat_src, 0_int64)) then
            nbad = nbad + 1
         end if
      end do
      if (nbad > 0) then
         n_fail = n_fail + nbad
         if (rank == 0) write (*, '(5a,i0,a)') "FAIL budget ", trim(label), "/", trim(scheme), &
            ": ", nbad, " decomposed factorization(s) disagree on the closed mass/salt/heat budget"
      else if (rank == 0) then
         write (*, '(4a)') "case ", trim(label), "/", trim(scheme), &
            " budget: IDENTICAL across factorizations"
      end if
   end subroutine compare_budget

   subroutine check_periodic_sponge_serial_out(label, scheme, ref_bud, dec_buds)
      !! The real regression gate for the periodic-seam / sponge bug, on
      !! "periodic_sponge" (periodic west/east, a relaxing sponge band on
      !! the closed north edge): the closed-budget mass/salt/heat totals
      !! (`out` AND `src`, AND their `*_active` gating flags) of the
      !! genuinely undecomposed SERIAL REFERENCE must equal, BIT FOR BIT,
      !! the same totals from every decomposed `px x py` factorization of
      !! the SAME real launch. The EFP
      !! reproducing-sums path (`reproducing_sums=.true.`, the default) is
      !! what makes that an EXACT identity on any compiler/toolchain/rank
      !! count rather than a round-off-level agreement -- no threshold,
      !! calibrated or otherwise, is needed or wanted here.
      !!
      !! Before the fix (the post-sponge seam-ghost refresh in
      !! `rdb_ocean_dyn.F90`, gated on
      !! `ocean_halo_is_decomposed_x() .or. ..._y()`), a real single-rank
      !! launch (`nprocs == 1`) skipped that refresh entirely -- both
      !! `ocean_halo_is_decomposed_{x,y}()` are false with no MPI seam at
      !! all -- so the reference's periodic-seam ghost columns kept the
      !! UN-relaxed value after the sponge touched only the tile's
      !! physical cells, and the corrector's advection read the stale
      !! ghosts before the next exchange: the reference's closed-budget
      !! totals came out WRONG BY ORDERS OF MAGNITUDE relative to any
      !! decomposed run of the same problem (the real Southern Ocean
      !! 1-degree case: Salt out 1.6e8 on 1 rank vs -81 on 2 ranks). That
      !! is exactly what this bitwise comparison catches, with no need to
      !! calibrate a margin against noise: reference and decomposed
      !! budgets of the SAME physical problem are either identical (fixed)
      !! or wildly different (broken), never "close but outside a
      !! threshold".
      !!
      !! `sponge_width` (set in `case_nml`) is deliberately narrow (3
      !! cells, `<=` the smallest per-rank tile at any tested `px x py`):
      !! a wider band that crosses an uneven y-decomposition tile boundary
      !! (found at `sponge_width = 6`, `1x4`: the legacy band sponge does
      !! not handle a band split across ranks) is a SEPARATE, pre-existing
      !! defect this case must not also exercise.
      !!
      !! `run_one`'s comm-override seam (`comm_env_push_compute_comm`,
      !! `src/comm/rdb_comm_env.F90`) is what makes `ref_bud` safe to trust
      !! here at EVERY `nprocs` the suite runs at (1, 2, 4), not just
      !! `nprocs == 1` as before: the reference engine always runs with
      !! `csize == 1` on every physical rank, and the console's collective
      !! is now scoped to a per-rank self-communicator for that call, so
      !! it is never contaminated by the real job size.
      !!
      !! At `nprocs == 1` there is no decomposed factorization to compare
      !! against (the `px x py` loop only ever produces 1x1, the SAME run
      !! as the reference) -- there this can only assert the reference
      !! budget came out finite and marked active, which is still a
      !! meaningful smoke check (a NaN/inactive budget here is itself a
      !! regression).
      character(len=*), intent(in) :: label, scheme
      type(conservation_budget_t), intent(in) :: ref_bud
      type(conservation_budget_t), intent(in) :: dec_buds(:)
      integer :: k, nbad

      if (trim(label) /= "periodic_sponge") return

      if (.not. ref_bud%mass_active) then
         n_fail = n_fail + 1
         write (*, '(4a,i0)') "FAIL serial-out ", trim(label), "/", trim(scheme), &
            ": reference budget was never computed (mass_active=.false.) on rank ", rank
         return
      end if
      if (.not. (ieee_is_finite(ref_bud%mass_out) .and. ieee_is_finite(ref_bud%salt_out) .and. &
                 ieee_is_finite(ref_bud%heat_out) .and. ieee_is_finite(ref_bud%mass_src) .and. &
                 ieee_is_finite(ref_bud%salt_src) .and. ieee_is_finite(ref_bud%heat_src))) then
         n_fail = n_fail + 1
         write (*, '(4a,i0)') "FAIL serial-out ", trim(label), "/", trim(scheme), &
            ": reference budget has a non-finite mass/salt/heat term on rank ", rank
         return
      end if

      if (size(dec_buds) == 0) then
         if (rank == 0) write (*, '(4a)') "case ", trim(label), "/", trim(scheme), &
            " serial-out: reference budget finite (no decomposed layout to compare at nprocs=1)"
         return
      end if

      nbad = 0
      do k = 1, size(dec_buds)
         if (dec_buds(k)%mass_active .neqv. ref_bud%mass_active .or. &
             dec_buds(k)%salt_active .neqv. ref_bud%salt_active .or. &
             dec_buds(k)%heat_active .neqv. ref_bud%heat_active) then
            nbad = nbad + 1
            cycle
         end if
         if (.not. dec_buds(k)%mass_active) then
            nbad = nbad + 1
            cycle
         end if
         if (transfer(dec_buds(k)%mass_out, 0_int64) /= transfer(ref_bud%mass_out, 0_int64) .or. &
             transfer(dec_buds(k)%salt_out, 0_int64) /= transfer(ref_bud%salt_out, 0_int64) .or. &
             transfer(dec_buds(k)%heat_out, 0_int64) /= transfer(ref_bud%heat_out, 0_int64) .or. &
             transfer(dec_buds(k)%mass_src, 0_int64) /= transfer(ref_bud%mass_src, 0_int64) .or. &
             transfer(dec_buds(k)%salt_src, 0_int64) /= transfer(ref_bud%salt_src, 0_int64) .or. &
             transfer(dec_buds(k)%heat_src, 0_int64) /= transfer(ref_bud%heat_src, 0_int64)) then
            nbad = nbad + 1
         end if
      end do

      if (nbad > 0) then
         n_fail = n_fail + nbad
         if (rank == 0) write (*, '(5a,i0,a,i0,a)') "FAIL serial-out ", trim(label), "/", trim(scheme), &
            ": ", nbad, " of ", size(dec_buds), " decomposed factorization(s) disagree "// &
            "bit-for-bit with the single-rank reference's closed mass/salt/heat budget"
      else if (rank == 0) then
         write (*, '(4a)') "case ", trim(label), "/", trim(scheme), &
            " serial-out: reference budget bit-identical to every decomposed factorization"
      end if
   end subroutine check_periodic_sponge_serial_out

   subroutine check_single_rank_fences()
      !! The single-rank features must be REFUSED at configure on more than
      !! one rank -- with the process grid left unset (px = py = 1, which
      !! the engine auto-factors and which `validate_config`'s px*py fences
      !! cannot see), so the engine-side gate is the one exercised.
      character(len=*), parameter :: NL = new_line("a")
      ! The cavity entries need a namelist `validate_config` ACCEPTS (FV_MOM6 PGF
      ! with `p_top_in_bc`, a sigma column, an analytic draft, virtual-freshwater
      ! melt with a matching top-drag coefficient): a bare `enable = .true.` is
      ! refused by `validate_config` before `engine_setup` runs, which made
      ! that fence vacuous.  The validate-accepts check below keeps every
      ! entry honest.  The cavity itself RUNS decomposed (the `cavity_*` cases);
      ! each entry below is one narrow combination that has no case yet.
      character(len=*), parameter :: CAV = &
                                     "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
                                     "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//new_line("a")// &
                                     "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//new_line("a")// &
                                     "&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
                                     "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"//new_line("a")
      character(len=*), parameter :: CAV_DYN = &
                                     "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
                                     "draft_depth = 300.0, draft_slope = 0.0, draft_x0 = 0.0, h_min_cavity = 10.0"
      character(len=*), parameter :: CAV_MELT = &
                                     "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
                                     "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
                                     "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0"
      character(len=1400), parameter :: KNOBS(8) = [character(len=1400) :: &
                                                    "&ocean_wetdry_nml enable = .true. /", &
                                                    "&ocean_bc_nml east = 'chapman' /", &
                                                    ! cavity + freshwater = 'mass'
                                                    CAV//CAV_DYN//" /"//new_line("a")// &
                                                    CAV_MELT//", freshwater = 'mass' /", &
                                                    ! cavity + a file draft
                                                    CAV//"&ocean_cavity_dyn_nml enable = .true., "// &
                                                    "draft_config = 'file', draft_file = 'draft.nc', "// &
                                                    "draft_var = 'draft', h_min_cavity = 10.0 /", &
                                                    ! cavity + bulk forcing
                                                    CAV//CAV_DYN//" /"//new_line("a")// &
                                                  "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, lat_south = 0.0 /"// &
                                                    new_line("a")// &
                                                    "&ocean_atm_forcing_nml enable = .true., filelist_tas = 'a.txt', "// &
                                                    "filelist_huss = 'a.txt', filelist_uas = 'a.txt', "// &
                                                    "filelist_vas = 'a.txt', filelist_psl = 'a.txt', "// &
                                                    "filelist_rlds = 'a.txt', filelist_rsds = 'a.txt', "// &
                                                    "filelist_prra = 'a.txt', filelist_prsn = 'a.txt' /"// &
                                                    new_line("a")// &
                                                    "&ocean_bulk_flux_nml enable = .true. /", &
                                                    ! cavity + an open edge
                                                    CAV//CAV_DYN//" /"//new_line("a")// &
                                                    "&ocean_bc_nml east = 'open' /", &
                                                    ! cavity + a sponge edge
                                                    CAV//CAV_DYN//" /"//new_line("a")// &
                                                    "&ocean_bc_nml north = 'sponge' /", &
                                                    ! cavity + the tripolar fold
                                                    CAV//CAV_DYN//" /"//new_line("a")// &
                                                    "&ocean_grid_nml grid_config = 'tripolar', "// &
                                                    "lat_south = 59.0, phi_join = 74.0, "// &
                                                    "coriolis_scheme = 'planetary' /"//new_line("a")// &
                                                    "&ocean_bc_nml west = 'periodic', "// &
                                                    "east = 'periodic', north = 'tripolar_fold' /"]
      ! One distinctive substring of the message each entry must produce, so a
      ! fence that is deleted while a LATER, more generic refusal still fires
      ! (the bulk-forcing one, the validate envelope) cannot pass unnoticed.
      character(len=64), parameter :: WHY(8) = [character(len=64) :: &
                                                "&ocean_wetdry_nml enable = .true. is single-rank", &
                                                "'chapman' edges are single-rank", &
                                                "freshwater = 'mass' with a cavity is unproven", &
                                                "draft_config = 'file' is unproven", &
                                                "enable with &ocean_bulk_flux_nml enable is unproven", &
                                                "is not 'wall' or 'periodic'", &
                                                "is not 'wall' or 'periodic'", &
                                                "'tripolar_fold' is unproven"]
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ik, ierr, im
      logical :: why_found, engine_ran

      do ik = 1, size(KNOBS)
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&grid_nml nx = 26, ny = 18, nghost = 3, "// &
               merge("dx = 1.0, dy = 1.0        ", "dx = 20000.0, dy = 20000.0", &
                     index(KNOBS(ik), "'spherical'") > 0)//" /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 900.0, start_date = '1958-01-01' /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'ssp_rk2' /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&output_nml output_to_file = .false. /"//NL// &
               trim(KNOBS(ik))//NL
         call read_config_from_string(nml, cfg, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call validate_config(cfg, ierr)
         if (ierr /= OCEAN_STATUS_OK .and. index(KNOBS(ik), "cavity_dyn") > 0) then
            write (*, '(a,i0)') "  fence entry ", ik
            write (*, '(a,i0)') "FAIL fence: the cavity namelist is refused by validate_config, "// &
               "so the engine fence is never reached (vacuous), rank ", rank
            n_fail = n_fail + 1
            cycle
         end if
         call error_ring_clear()
         engine_ran = ierr == OCEAN_STATUS_OK
         if (ierr == OCEAN_STATUS_OK) then
            call engine_setup(engine, cfg, ierr, compute_rank=rank, compute_size=nprocs)
            call engine_teardown(engine)
         end if
         why_found = .false.
         do im = 0, error_ring_count() - 1
            if (index(error_ring_get(im), trim(WHY(ik))) > 0) why_found = .true.
         end do
         ! Entries `validate_config` itself refuses (wet/dry) never reach the engine
         ! fence: status-only, no reason to read.
         if (engine_ran .and. ierr /= OCEAN_STATUS_OK .and. .not. why_found) then
            if (error_ring_count() > 0) write (*, '(2a)') "  first ring message: ", &
               trim(error_ring_get(0))
            write (*, '(3a,i0)') "FAIL fence: '", trim(KNOBS(ik)), &
               "' was refused, but NOT for the expected reason (no error-ring message "// &
               "contains '"//trim(WHY(ik))//"'), rank ", rank
            n_fail = n_fail + 1
            cycle
         end if
         if (ierr == OCEAN_STATUS_OK) then
            write (*, '(3a,i0)') "FAIL fence: '", trim(KNOBS(ik)), &
               "' was ACCEPTED on multi-rank, rank ", rank
            n_fail = n_fail + 1
         else if (rank == 0) then
            write (*, '(3a)') "case fence: '", trim(KNOBS(ik)), &
               "' refused on multi-rank for the expected reason (ok)"
         end if
      end do
   end subroutine check_single_rank_fences

#ifndef RDB_NO_NETCDF
   subroutine write_input_files()
      !! Rank 0 writes the three whole-grid input files the "file_readers"
      !! case reads through the per-rank windowed readers: a periodic lon-lat
      !! MOM6 supergrid (360 deg in x, 20-56 N), a bathymetry with a land
      !! block stored C-ORDER (Fortran dims (y, x) — the reader's transpose
      !! path) and a z-level T/S initial condition (x, y, z).  Every field
      !! varies in both directions, so a tile that read the wrong window
      !! cannot match the serial run.
      integer, parameter :: NZS = 8
      real(wp), parameter :: DEG = 3.14159265358979323846_wp/180.0_wp
      real(wp), parameter :: RE = 6.371e6_wp, LAT0 = 20.0_wp, DLAT = 2.0_wp
      real(wp) :: dlon
      integer :: ncid, dnxp, dnyp, dnx, dny, dx_, dy_, dz_, v1, v2, v3, v4, v5, ierr
      integer :: m, n, i, j, k
      real(wp), allocatable :: sx(:, :), sy(:, :), sdx(:, :), sdy(:, :), sar(:, :)
      real(wp), allocatable :: byx(:, :), tt(:, :, :), ss(:, :, :)
      real(wp) :: zs(NZS), lat, x, y

      dlon = 360.0_wp/real(NX_G, wp)
      allocate (sx(2*NX_G + 1, 2*NY_G + 1), sy(2*NX_G + 1, 2*NY_G + 1))
      allocate (sdx(2*NX_G, 2*NY_G + 1), sdy(2*NX_G + 1, 2*NY_G), sar(2*NX_G, 2*NY_G))
      do n = 1, 2*NY_G + 1
         lat = LAT0 + real(n - 1, wp)*0.5_wp*DLAT
         do m = 1, 2*NX_G + 1
            sx(m, n) = real(m - 1, wp)*0.5_wp*dlon
            sy(m, n) = lat
         end do
         do m = 1, 2*NX_G
            sdx(m, n) = RE*cos(lat*DEG)*0.5_wp*dlon*DEG
         end do
      end do
      sdy = RE*0.5_wp*DLAT*DEG
      do n = 1, 2*NY_G
         do m = 1, 2*NX_G
            sar(m, n) = sdx(m, n)*sdy(m, n)
         end do
      end do
      call nc_create_file(SG_FILE, ncid)
      call nc_def_dim(ncid, "nxp", 2*NX_G + 1, dnxp)
      call nc_def_dim(ncid, "nyp", 2*NY_G + 1, dnyp)
      call nc_def_dim(ncid, "nx", 2*NX_G, dnx)
      call nc_def_dim(ncid, "ny", 2*NY_G, dny)
      call nc_def_var_2d(ncid, "x", [dnxp, dnyp], v1)
      call nc_def_var_2d(ncid, "y", [dnxp, dnyp], v2)
      call nc_def_var_2d(ncid, "dx", [dnx, dnyp], v3)
      call nc_def_var_2d(ncid, "dy", [dnxp, dny], v4)
      call nc_def_var_2d(ncid, "area", [dnx, dny], v5)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, v1, sx)
      call nc_put_var_2d(ncid, v2, sy)
      call nc_put_var_2d(ncid, v3, sdx)
      call nc_put_var_2d(ncid, v4, sdy)
      call nc_put_var_2d(ncid, v5, sar)
      call nc_close(ncid)

      ! Bathymetry, positive-down, C-order: byx(j, i).  A shelf rising to
      ! the north, a Gaussian ridge, and a small land block.
      allocate (byx(NY_G, NX_G))
      do i = 1, NX_G
         do j = 1, NY_G
            x = real(i, wp)/real(NX_G, wp)
            y = real(j, wp)/real(NY_G, wp)
            byx(j, i) = 3000.0_wp - 1800.0_wp*y**2 - &
                        900.0_wp*exp(-((x - 0.3_wp)**2 + (y - 0.5_wp)**2)/0.02_wp)
            if (i >= 17 .and. i <= 19 .and. j >= 7 .and. j <= 10) byx(j, i) = 0.0_wp
         end do
      end do
      call nc_create_file(BATHY_FILE, ncid)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_var_2d(ncid, "depth", [dy_, dx_], v1)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, v1, byx)
      call nc_close(ncid)

      ! z-level T/S (x, y, z), positive-down source depths.
      allocate (tt(NX_G, NY_G, NZS), ss(NX_G, NY_G, NZS))
      do k = 1, NZS
         zs(k) = real(k - 1, wp)*450.0_wp
      end do
      do k = 1, NZS
         do j = 1, NY_G
            do i = 1, NX_G
               x = real(i, wp)/real(NX_G, wp)
               y = real(j, wp)/real(NY_G, wp)
               tt(i, j, k) = 22.0_wp - 18.0_wp*y - zs(k)/250.0_wp + 1.5_wp*sin(6.2831853_wp*x)
               tt(i, j, k) = max(tt(i, j, k), -1.5_wp)
               ss(i, j, k) = 34.5_wp + 0.6_wp*y - 0.2_wp*cos(6.2831853_wp*x) + zs(k)*1.0e-4_wp
            end do
         end do
      end do
      call nc_create_file(ZINIT_FILE, ncid)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "z", NZS, dz_)
      call nc_def_var_3d(ncid, "temp", [dx_, dy_, dz_], v1)
      call nc_def_var_3d(ncid, "salt", [dx_, dy_, dz_], v2)
      call rdb_def_var_1d(ncid, "z_src", dz_, v3)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, v1, tt)
      ierr = nf90_put_var(ncid, v2, ss)
      call rdb_put_var_1d(ncid, v3, zs)
      call nc_close(ncid)

      ! SSS restoring target (x, y, time): two records spanning the run,
      ! varying in both directions.
      deallocate (tt)
      allocate (tt(NX_G, NY_G, 2))
      do j = 1, NY_G
         do i = 1, NX_G
            x = real(i, wp)/real(NX_G, wp)
            y = real(j, wp)/real(NY_G, wp)
            tt(i, j, 1) = 34.0_wp + 1.5_wp*y + 0.4_wp*sin(6.2831853_wp*x)
            tt(i, j, 2) = tt(i, j, 1) + 0.3_wp*x - 0.2_wp*y
         end do
      end do
      call nc_create_file(SSS_FILE, ncid)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "time", 2, dz_)
      call nc_def_var_3d(ncid, "salt", [dx_, dy_, dz_], v1)
      call rdb_def_var_1d(ncid, "time", dz_, v2)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, v1, tt)
      call rdb_put_var_1d(ncid, v2, [0.0_wp, 86400.0_wp])
      call nc_close(ncid)

      ! Internal-tide fields (x, y): a tidal speed and a roughness variance
      ! that differ from each other and vary in both directions, so a tile
      ! that read the wrong window cannot match the serial run.
      deallocate (tt)
      allocate (tt(NX_G, NY_G, 2))
      do j = 1, NY_G
         do i = 1, NX_G
            x = real(i, wp)/real(NX_G, wp)
            y = real(j, wp)/real(NY_G, wp)
            tt(i, j, 1) = 0.02_wp + 0.04_wp*x*(1.0_wp - y)
            tt(i, j, 2) = 2000.0_wp + 6000.0_wp*y*(1.0_wp - 0.5_wp*x)
         end do
      end do
      call nc_create_file(TIDE_FILE, ncid)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_var_2d(ncid, "tideamp", [dx_, dy_], v1)
      call nc_def_var_2d(ncid, "h2", [dx_, dy_], v2)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, v1, tt(:, :, 1))
      call nc_put_var_2d(ncid, v2, tt(:, :, 2))
      call nc_close(ncid)

      ! Chlorophyll climatology (x, y, time): two records, strongly
      ! non-uniform and non-symmetric in both directions (a low-chl corner
      ! against a bloom), so a wrong window or a stale seam ghost changes the
      ! band opacities.
      tt(:, :, 2) = 0.0_wp
      do j = 1, NY_G
         do i = 1, NX_G
            x = real(i, wp)/real(NX_G, wp)
            y = real(j, wp)/real(NY_G, wp)
            tt(i, j, 1) = 0.01_wp + 0.25_wp*x**2*(1.0_wp - 0.7_wp*y) + &
                          0.12_wp*abs(sin(9.0_wp*x + 4.0_wp*y))
            tt(i, j, 2) = 1.5_wp*tt(i, j, 1) + 0.03_wp*y
         end do
      end do
      call nc_create_file(CHL_FILE, ncid)
      call nc_def_dim(ncid, "x", NX_G, dx_)
      call nc_def_dim(ncid, "y", NY_G, dy_)
      call nc_def_dim(ncid, "time", 2, dz_)
      call nc_def_var_3d(ncid, "CHL", [dx_, dy_, dz_], v1)
      call rdb_def_var_1d(ncid, "time", dz_, v2)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, v1, tt)
      call rdb_put_var_1d(ncid, v2, [0.0_wp, 172800.0_wp])
      call nc_close(ncid)
      call write_atm_file(ATM_WARM, ATM_WARM_LIST, 288.0_wp, 4.0e-3_wp, 72, 36)
      call write_atm_file(ATM_COLD, ATM_COLD_LIST, 255.0_wp, 6.0e-4_wp, 72, 36)
      ! 30 x 10 degree source for the tripolar cap (see the bulk_flux_fold case).
      call write_atm_file(ATM_CRSE, ATM_CRSE_LIST, 288.0_wp, 4.0e-3_wp, 12, 18)
   end subroutine write_input_files

   subroutine write_atm_file(path, listname, tas0, huss0, nlon, nlat)
      !! One global lon-lat file (5 degrees, two time records, all nine CMOR
      !! names) of spatially varying, east-west and north-south ASYMMETRIC
      !! prescribed atmospheric state (affine + sinusoid in lon/lat), plus
      !! the one-line filelist the namelist points at.
      character(len=*), intent(in) :: path, listname
      real(wp), intent(in) :: tas0, huss0
      integer, intent(in) :: nlon, nlat
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      character(len=4), parameter :: VN(9) = [character(len=4) :: "tas", "huss", "uas", &
                                                                                       "vas", "psl", "rlds", "rsds", "prra", "prsn"]
      real(wp) :: lon(nlon), lat(nlat), lonb(2, nlon), latb(2, nlat), tax(2), f(nlon, nlat, 2)
      real(wp) :: sx, sy, dlo, dla
      integer :: ncid, dlon_, dlat_, dbnd, dtim, vlon, vlat, vlonb, vlatb, vtim, vf(9)
      integer :: v, i, j, k, ierr, u

      dlo = 360.0_wp/real(nlon, wp)
      dla = 180.0_wp/real(nlat, wp)
      do i = 1, nlon
         lonb(1, i) = dlo*real(i - 1, wp)
         lonb(2, i) = dlo*real(i, wp)
         lon(i) = 0.5_wp*(lonb(1, i) + lonb(2, i))
      end do
      do j = 1, nlat
         latb(1, j) = -90.0_wp + dla*real(j - 1, wp)
         latb(2, j) = -90.0_wp + dla*real(j, wp)
         lat(j) = 0.5_wp*(latb(1, j) + latb(2, j))
      end do
      tax = [0.0_wp, 10.0_wp]
      call nc_create_file(path, ncid)
      call nc_def_dim(ncid, "lon", nlon, dlon_)
      call nc_def_dim(ncid, "lat", nlat, dlat_)
      call nc_def_dim(ncid, "bnds", 2, dbnd)
      call nc_def_dim(ncid, "time", 2, dtim)
      call rdb_def_var_1d(ncid, "lon", dlon_, vlon)
      call rdb_def_var_1d(ncid, "lat", dlat_, vlat)
      call nc_def_var_2d(ncid, "lon_bnds", [dbnd, dlon_], vlonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dbnd, dlat_], vlatb)
      call rdb_def_var_1d(ncid, "time", dtim, vtim)
      call nc_put_att(ncid, vtim, "units", "days since 2000-06-20 00:00:00")
      do v = 1, 9
         call nc_def_var_3d(ncid, trim(VN(v)), [dlon_, dlat_, dtim], vf(v))
      end do
      call nc_enddef(ncid)
      call rdb_put_var_1d(ncid, vlon, lon)
      call rdb_put_var_1d(ncid, vlat, lat)
      ierr = nf90_put_var(ncid, vlonb, lonb)
      ierr = nf90_put_var(ncid, vlatb, latb)
      call rdb_put_var_1d(ncid, vtim, tax)
      do v = 1, 9
         do k = 1, 2
            do j = 1, nlat
               do i = 1, nlon
                  sx = sin(2.0_wp*PI*lon(i)/37.0_wp + 0.3_wp)
                  sy = cos(2.0_wp*PI*lat(j)/29.0_wp - 0.5_wp)
                  select case (v)
                  case (1)
                     f(i, j, k) = tas0 + 4.0_wp*sx + 0.15_wp*(lat(j) - 20.0_wp) + &
                                  0.05_wp*tax(k)
                  case (2)
                     f(i, j, k) = huss0*(1.0_wp + 0.2_wp*sy + 0.1_wp*sx)
                  case (3)
                     f(i, j, k) = 6.0_wp + 4.0_wp*sx + 0.05_wp*(lat(j) - 20.0_wp)
                  case (4)
                     f(i, j, k) = -3.0_wp + 3.0_wp*sy + 0.04_wp*lon(i) + 0.02_wp*tax(k)
                  case (5)
                     f(i, j, k) = 1.01e5_wp - 300.0_wp*sx + 5.0_wp*lat(j)
                  case (6)
                     f(i, j, k) = 330.0_wp + 15.0_wp*sy + 20.0_wp*sx
                  case (7)
                     f(i, j, k) = 180.0_wp + 40.0_wp*sx*sy + 0.5_wp*lat(j)
                  case (8)
                     f(i, j, k) = 2.0e-5_wp*(1.0_wp + 0.5_wp*sx)
                  case default
                     f(i, j, k) = 1.0e-6_wp*(1.0_wp + 0.5_wp*sy)
                  end select
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vf(v), f)
      end do
      call nc_close(ncid)
      open (newunit=u, file=listname, status='replace', action='write')
      write (u, '(a)') trim(path)
      close (u)
   end subroutine write_atm_file
#endif

end program test_ocean_decomp_bitid_mpi
#else
program test_ocean_decomp_bitid_mpi
   implicit none
   write (*, '(a)') "test_ocean_decomp_bitid_mpi: skipped (RDB_ENABLE_MPI=OFF)"
end program test_ocean_decomp_bitid_mpi
#endif
