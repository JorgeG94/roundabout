!! Tests for the ocean bulk-flux consumer (`rdb_ocean_bulk_forcing`,
!! OM3 forcing wave-3, C2c-3a) -- PLAN_C2c.md §5 tests 4-5 (ocean half).
!!
!! Every synthetic `atm_state_t` here is built BY HAND (the nine field
!! arrays allocated and filled directly, `enable = .true.` set
!! directly) -- no NetCDF reader, no filelist, matching the brief's
!! "synthetic in-memory atm state is fine, no files needed". Likewise
!! `bulk_forcing_t` is hand-built rather than run through
!! `bulk_forcing_configure` (which needs a full parsed `config_t`):
!! every field it reads at apply time (`cfg_atm`, `cfg_bulk`,
!! `taux_c`/`tauy_c`/`open_water` scratch) is set directly.
!!
!! Tests:
!!   bulk_forcing_matches_column_and_rotates -- §5 test 4: uniform atm
!!     + uniform SST + zero current, ice disabled (open_water == 1) ==>
!!     every physical cell's q_sens/q_lat/evap/taux_c/tauy_c equals a
!!     DIRECT `bulk_flux_column` call BITWISE; the interior stress FACE
!!     equals the (uniform) cell value; q_sw == (1-albedo)*rsds; q_lw
!!     matches `bulk_flux_longwave_up` exactly. A second sub-case sets
!!     `angle_dx = 30 deg` uniformly and checks the output stress is
!!     the EXACT 30-degree rotation (to 1e-13) of the angle=0 case's
!!     stress -- a property of this bulk formula (stress is parallel
!!     to the rotated relative-wind vector, whose magnitude rotation
!!     does not change), not merely a code-structure echo.
!!   bulk_forcing_open_water_halves -- §5 test 5 (ocean half): 50% ice
!!     cover (`part_size(:,:,0) = 0.5`) halves every ocean component
!!     (q_sens, q_lat, evap, taux_c, tauy_c, q_sw, q_lw, lprec, fprec)
!!     relative to the ice-disabled run, to round-off; `p_surf_atm` is
!!     NOT scaled (it is a pressure-loading term, not a flux).
module test_ocean_bulk_forcing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, DEG2RAD
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_atm_state, only: atm_state_t
   use rdb_ocean_bulk_flux, only: bulk_flux_column, bulk_flux_longwave_up, bulk_flux_qsat_ice, &
                                  bulk_flux_config_t, BULK_T0_KELVIN, BULK_STEFAN_BOLTZMANN, &
                                  BULK_GAS_CONST_AIR, BULK_TVQ_AIR, BULK_CP_AIR, &
                                  BULK_LATENT_SUBLIMATION, BULK_QSAT_ICE_COEF_B, BULK_WIND_MIN
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_bulk_forcing, only: bulk_forcing_t, bulk_forcing_apply
   use rdb_ice_atm_forcing, only: ice_atm_forcing_from_atm
   implicit none
   private

   public :: collect_ocean_bulk_forcing_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 3, NYP = 3
   real(wp), parameter :: TOL = 1.0e-12_wp
   real(wp), parameter :: TOL_ROT = 1.0e-13_wp

   ! Uniform synthetic atmosphere + ocean state.
   real(wp), parameter :: UAS = 5.0_wp, VAS = 3.0_wp
   real(wp), parameter :: TAS = 290.0_wp, HUSS = 0.0080_wp, PSL = 101000.0_wp
   real(wp), parameter :: RLDS = 300.0_wp, RSDS = 200.0_wp
   real(wp), parameter :: PRRA = 2.0e-5_wp, PRSN = 1.0e-6_wp
   real(wp), parameter :: SST_DEGC = 15.0_wp, H0 = 10.0_wp
   real(wp), parameter :: Z_WIND = 10.0_wp, Z_TA = 2.0_wp
   integer, parameter :: N_ITER = 5
   real(wp), parameter :: EMISSIVITY = 0.97_wp, OCEAN_ALBEDO = 0.066_wp

contains

   subroutine collect_ocean_bulk_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bulk_forcing_matches_column", test_matches_column), &
                  new_unittest("bulk_forcing_rotation_exact", test_rotation_exact), &
                  new_unittest("bulk_forcing_open_water_halves", test_open_water_halves), &
                  new_unittest("ice_atm_forcing_seb_linearization", test_ice_seb_linearization) &
                  ]
   end subroutine collect_ocean_bulk_forcing_tests

   ! =================================================================
   ! Test support
   ! =================================================================

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine init_halo_single_rank(grid)
      !! `ocean_halo_centre`/`ocean_halo_face_x/y` module state -- shared
      !! singleton, so each test owns one init/destroy pair. Non-periodic,
      !! single rank: ghosts are a true no-op, which is fine here since
      !! every assertion below reads an INTERIOR face/cell.
      type(hgrid_t), intent(in) :: grid
      type(decomp_t) :: decomp
      logical, save :: comm_inited = .false.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
      call decomp_init(decomp, grid%nx_phys, grid%ny_phys, 1, 1, 0)
      call ocean_halo_init(decomp, grid%nghost, .false., .false.)
   end subroutine init_halo_single_rank

   subroutine make_atm(grid, this)
      !! Hand-built nine-field synthetic atmospheric state -- spatially
      !! uniform, every field including ghosts (no consumer here reads
      !! a ghost cell, but filling them removes any doubt).
      type(hgrid_t), intent(in) :: grid
      type(atm_state_t), intent(inout) :: this
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (this%tas(nx, ny), source=TAS)
      allocate (this%huss(nx, ny), source=HUSS)
      allocate (this%uas(nx, ny), source=UAS)
      allocate (this%vas(nx, ny), source=VAS)
      allocate (this%psl(nx, ny), source=PSL)
      allocate (this%rlds(nx, ny), source=RLDS)
      allocate (this%prra(nx, ny), source=PRRA)
      allocate (this%prsn(nx, ny), source=PRSN)
      allocate (this%rsds(nx, ny), source=RSDS)
      allocate (this%rsds_rec(nx, ny), source=RSDS)
      allocate (this%cosz(nx, ny), source=0.0_wp)
      allocate (this%cosz_mean(nx, ny), source=0.0_wp)
      ! Cached geometry (`atm_state_configure`'s own job in production):
      ! `ice_atm_forcing_from_atm` reads `this%nx_phys/ny_phys/nghost` --
      ! unlike `bulk_forcing_apply`, which takes `grid` separately and
      ! never looks at these -- so a hand-built atm must set them too.
      this%nx_total = nx
      this%ny_total = ny
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%nghost = grid%nghost
      this%enable = .true.
   end subroutine make_atm

   subroutine make_ms(grid, ms)
      !! Single-layer (`nz_ml = 1`), uniform SST, zero current.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      ms%nz_ml = 1
      call ms%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = SST_DEGC*H0
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine make_ms

   subroutine make_bulk_forcing(grid, wind_relative, this)
      !! Hand-built workspace -- bypasses `bulk_forcing_configure` (no
      !! `config_t` needed): sets exactly the fields `bulk_forcing_apply`
      !! reads.
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: wind_relative
      type(bulk_forcing_t), intent(inout) :: this
      this%enable = .true.
      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      this%cfg_atm%wind_relative = wind_relative
      this%cfg_atm%ocean_albedo = OCEAN_ALBEDO
      this%cfg_bulk%n_iter = N_ITER
      this%cfg_bulk%z_wind = Z_WIND
      this%cfg_bulk%z_ta = Z_TA
      this%cfg_bulk%emissivity = EMISSIVITY
      allocate (this%taux_c(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (this%tauy_c(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (this%open_water(grid%nx_total, grid%ny_total), source=1.0_wp)
   end subroutine make_bulk_forcing

   ! =================================================================
   ! Test 4 — matches a direct bulk_flux_column call + rotates exactly
   ! =================================================================

   subroutine test_matches_column(error)
      !! angle_dx = 0 (default): every physical cell's output equals a
      !! DIRECT `bulk_flux_column` call bitwise (no trig recombination
      !! is involved at angle=0 -- `uas`/`vas` pass straight through --
      !! so this holds to the exact last bit, thr=0).
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf

      real(wp) :: taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar
      real(wp) :: sst_k, lw_up, q_sw_ref, q_lw_ref
      integer :: ip, jp, face_x, face_y

      checks: block
         call make_grid(grid)
         call make_cartesian_metrics(metrics, grid)
         call init_halo_single_rank(grid)

         ip = grid%nghost + 2 ! middle physical column -- both faces interior
         jp = grid%nghost + 2 ! middle physical row
         face_x = ip ! west face of cell ip (both neighbours physical)
         face_y = jp ! south face of cell jp

         call make_ms(grid, ms)
         call make_atm(grid, atm)
         call make_bulk_forcing(grid, .true., bf)
         call ss%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)

         call run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

         ! Direct reference: bulk_flux_column with zero current, the
         ! same uniform atm + SST.
         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, SST_DEGC, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar)

         call check(error, sf%q_sens(ip, jp), sens_ref, "q_sens == direct bulk_flux_column", thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, sf%q_lat(ip, jp), lat_ref, "q_lat == direct bulk_flux_column", thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, sf%evap(ip, jp), -evap_ref, "evap == -direct bulk_flux_column (MOM6 sign)", &
                    thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, bf%taux_c(ip, jp), taux_ref, "taux_c == direct bulk_flux_column", thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, bf%tauy_c(ip, jp), tauy_ref, "tauy_c == direct bulk_flux_column", thr=0.0_wp)
         if (allocated(error)) exit checks

         ! Interior stress FACE == the (uniform) cell value.
         call check(error, ss%tau_x(face_x, jp), taux_ref, "tau_x face == uniform cell value", thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, ss%tau_y(ip, face_y), tauy_ref, "tau_y face == uniform cell value", thr=0.0_wp)
         if (allocated(error)) exit checks

         ! q_sw / q_lw.
         q_sw_ref = (1.0_wp - OCEAN_ALBEDO)*RSDS
         call check(error, sf%q_sw(ip, jp), q_sw_ref, "q_sw == (1-albedo)*rsds", thr=0.0_wp)
         if (allocated(error)) exit checks
         sst_k = SST_DEGC + BULK_T0_KELVIN
         lw_up = bulk_flux_longwave_up(sst_k, EMISSIVITY)
         q_lw_ref = RLDS - lw_up
         call check(error, sf%q_lw(ip, jp), q_lw_ref, "q_lw == rlds - bulk_flux_longwave_up", thr=0.0_wp)
         if (allocated(error)) exit checks

         call ocean_halo_destroy()
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_matches_column

   subroutine test_rotation_exact(error)
      !! angle_dx = 30 deg uniformly: the output stress is the EXACT
      !! 30-degree rotation (to round-off, `TOL_ROT = 1e-13`) of the
      !! angle=0 stress -- a property of this bulk formula (stress is
      !! parallel to the rotated relative wind, whose SPEED -- hence
      !! cd/ch/ce -- is rotation-invariant), not merely a code-structure
      !! echo. Every comparison here uses `TOL_ROT`, not an exact
      !! `thr=0`: `cos(30deg)*uas(i,j) + sin(30deg)*vas(i,j)` and this
      !! test's own identically-written scalar recombination are two
      !! SEPARATE compilation units, and `-O3 -march=native` is free to
      !! contract either into a differently-rounded FMA sequence -- a
      !! round-off artifact, not a correctness bug (the angle=0 test
      !! above, which involves no trig recombination at all, holds to
      !! the exact bit).
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf

      real(wp) :: taux_ref, tauy_ref, d3, d4, d5, cd, ch, ce, ustar
      real(wp) :: taux_rot_direct, tauy_rot_direct, d1, d2
      real(wp) :: taux_rot_expected, tauy_rot_expected
      real(wp) :: cosA, sinA, angle
      integer :: ip, jp

      checks: block
         call make_grid(grid)
         call make_cartesian_metrics(metrics, grid)
         call init_halo_single_rank(grid)

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         angle = 30.0_wp*DEG2RAD
         metrics%angle_dx = angle
         !$omp target update to(metrics%angle_dx)

         call make_ms(grid, ms)
         call make_atm(grid, atm)
         call make_bulk_forcing(grid, .true., bf)
         call ss%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)

         call run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

         ! Reference #1: bulk_flux_column with the UNROTATED wind.
         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, SST_DEGC, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, d3, d4, d5, cd, ch, ce, ustar)

         ! Reference #2: bulk_flux_column with the ROTATED wind.
         cosA = cos(angle)
         sinA = sin(angle)
         d1 = cosA*UAS + sinA*VAS
         d2 = -sinA*UAS + cosA*VAS
         call bulk_flux_column(d1, d2, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, SST_DEGC, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_rot_direct, tauy_rot_direct, d3, d4, d5, cd, ch, ce, ustar)

         call check(error, bf%taux_c(ip, jp), taux_rot_direct, &
                    "rotated taux_c == direct bulk_flux_column on rotated wind", thr=TOL_ROT)
         if (allocated(error)) exit checks
         call check(error, bf%tauy_c(ip, jp), tauy_rot_direct, &
                    "rotated tauy_c == direct bulk_flux_column on rotated wind", thr=TOL_ROT)
         if (allocated(error)) exit checks

         ! Property check: rotating the ANGLE=0 stress by 30 deg
         ! reproduces the rotated-wind column answer to 1e-13.
         taux_rot_expected = cosA*taux_ref + sinA*tauy_ref
         tauy_rot_expected = -sinA*taux_ref + cosA*tauy_ref
         call check(error, taux_rot_expected, taux_rot_direct, &
                    "stress rotates by exactly 30 deg (taux)", thr=TOL_ROT)
         if (allocated(error)) exit checks
         call check(error, tauy_rot_expected, tauy_rot_direct, &
                    "stress rotates by exactly 30 deg (tauy)", thr=TOL_ROT)
         if (allocated(error)) exit checks

         call ocean_halo_destroy()
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_rotation_exact

   ! =================================================================
   ! Test 5 (ocean half) — open-water weighting halves every component
   ! =================================================================

   subroutine test_open_water_halves(error)
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice_off, ice_on
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf_off, bf_on
      type(ocean_surface_stress_t) :: ss_off, ss_on
      type(ocean_surface_flux_t) :: sf_off, sf_on
      integer :: ip, jp

      checks: block
         call make_grid(grid)
         call make_cartesian_metrics(metrics, grid)
         call init_halo_single_rank(grid)

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         call make_ms(grid, ms)
         call make_atm(grid, atm)

         ! ---- ice disabled: the reference run ----
         call make_bulk_forcing(grid, .true., bf_off)
         call ss_off%init(grid)
         call sf_off%init(grid)
         call sf_off%set_components(grid, .true.)
         call run_apply(grid, metrics, bc, atm, ms, ice_off, bf_off, ss_off, sf_off)

         ! ---- ice enabled, 50% open water ----
         ice_on%ncat = 1
         call ice_on%init(grid)
         ice_on%enable = .true.
         ice_on%part_size(:, :, 0) = 0.5_wp
         ice_on%part_size(:, :, 1) = 0.5_wp

         call make_bulk_forcing(grid, .true., bf_on)
         call ss_on%init(grid)
         call sf_on%init(grid)
         call sf_on%set_components(grid, .true.)
         call run_apply(grid, metrics, bc, atm, ms, ice_on, bf_on, ss_on, sf_on)

         call check(error, sf_on%q_sens(ip, jp), 0.5_wp*sf_off%q_sens(ip, jp), &
                    "50% ice halves q_sens", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%q_lat(ip, jp), 0.5_wp*sf_off%q_lat(ip, jp), &
                    "50% ice halves q_lat", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%evap(ip, jp), 0.5_wp*sf_off%evap(ip, jp), &
                    "50% ice halves evap", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, bf_on%taux_c(ip, jp), 0.5_wp*bf_off%taux_c(ip, jp), &
                    "50% ice halves taux_c", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, bf_on%tauy_c(ip, jp), 0.5_wp*bf_off%tauy_c(ip, jp), &
                    "50% ice halves tauy_c", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%q_sw(ip, jp), 0.5_wp*sf_off%q_sw(ip, jp), &
                    "50% ice halves q_sw", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%q_lw(ip, jp), 0.5_wp*sf_off%q_lw(ip, jp), &
                    "50% ice halves q_lw", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%lprec(ip, jp), 0.5_wp*sf_off%lprec(ip, jp), &
                    "50% ice halves lprec", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%fprec(ip, jp), 0.5_wp*sf_off%fprec(ip, jp), &
                    "50% ice halves fprec", thr=TOL)
         if (allocated(error)) exit checks

         ! p_surf_atm is a pressure-loading term, not a flux -- NOT scaled.
         call check(error, sf_on%p_surf_atm(ip, jp), sf_off%p_surf_atm(ip, jp), &
                    "p_surf_atm is not open-water weighted", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, sf_on%p_surf_atm(ip, jp), PSL, &
                    "p_surf_atm == psl", thr=TOL)
         if (allocated(error)) exit checks

         call ocean_halo_destroy()
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_open_water_halves

   subroutine run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)
      !! `bulk_forcing_apply` with the production device residency: the
      !! engine maps every object before stepping, and the GPU build makes
      !! no implicit host<->device copies (mem:separate), so the stress
      !! seam refresh in particular needs device-resident state. Map,
      !! apply, bring back what the tests read (component arrays only,
      !! never an aggregate), unmap. Directives are inert on CPU builds.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_bc_state_t), intent(in) :: bc
      type(atm_state_t), intent(inout) :: atm
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(bulk_forcing_t), intent(inout) :: bf
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf

      !$omp target enter data map(to: atm, ms, bf, ss, sf)
      call atm%enter_data()
      call ms%enter_data()
      call bf%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      if (ice%enable) then
         !$omp target enter data map(to: ice)
         call ice%enter_data()
      end if

      call bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

      !$omp target update from(sf%q_sens, sf%q_lat, sf%evap, sf%q_sw, sf%q_lw)
      !$omp target update from(sf%lprec, sf%fprec, sf%p_surf_atm)
      !$omp target update from(bf%taux_c, bf%tauy_c, ss%tau_x, ss%tau_y)

      if (ice%enable) then
         call ice%exit_data()
         !$omp target exit data map(delete: ice)
      end if
      call sf%exit_data()
      call ss%exit_data()
      call bf%exit_data()
      call ms%exit_data()
      call atm%exit_data()
      !$omp target exit data map(delete: atm, ms, bf, ss, sf)
   end subroutine run_apply

   ! =================================================================
   ! Test 5 (ice half) — SEB linearisation matches a hand computation
   ! =================================================================

   subroutine test_ice_seb_linearization(error)
      !! PLAN_C2c.md S5 test 5 (ice half): `ice_atm_forcing_from_atm`'s
      !! `atm_sf0`/`atm_dsfdt` match an independently-written evaluation
      !! of the LY09 over-ice SEB linearisation at `T_s = 0 degC`, built
      !! from the SAME public primitives (`bulk_flux_qsat_ice`,
      !! `BULK_STEFAN_BOLTZMANN`, ...) the production kernel uses but
      !! written out as its own expression here (not a call to a shared
      !! "reference" routine, since none exists -- the formula IS the
      !! implementation). `atm_sw_dn`/`atm_fprec` are plain pass-through
      !! checks (`rsds`/`prsn`).
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(atm_state_t) :: atm
      type(ocean_sea_ice_t) :: ice
      type(bulk_flux_config_t) :: cfg_bulk

      real(wp), parameter :: ICE_BULK_CH = 1.63e-3_wp, ICE_BULK_CE = 1.63e-3_wp
      real(wp) :: rho_air, wv, qsat0, dqsat0_dt, lw_t3, sens_coef, lat_coef
      real(wp) :: sf0_ref, dsfdt_ref
      integer :: ip, jp

      checks: block
         call make_grid(grid)
         call make_atm(grid, atm)

         ice%ncat = 1
         call ice%init(grid)

         cfg_bulk%emissivity = EMISSIVITY

         ! Production device residency (mem:separate: no implicit copies):
         ! map both objects, run the seam, bring the four outputs back.
         !$omp target enter data map(to: atm, ice)
         call atm%enter_data()
         call ice%enter_data()
         call ice_atm_forcing_from_atm(ice, atm, cfg_bulk)
         !$omp target update from(ice%atm_sf0, ice%atm_dsfdt, ice%atm_sw_dn, ice%atm_fprec)
         call ice%exit_data()
         call atm%exit_data()
         !$omp target exit data map(delete: atm, ice)

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         ! Independent evaluation of the plan's S4 formula.
         rho_air = PSL/(BULK_GAS_CONST_AIR*TAS*(1.0_wp + BULK_TVQ_AIR*HUSS))
         wv = max(sqrt(UAS**2 + VAS**2), BULK_WIND_MIN)
         qsat0 = bulk_flux_qsat_ice(BULK_T0_KELVIN, rho_air)
         dqsat0_dt = qsat0*BULK_QSAT_ICE_COEF_B/(BULK_T0_KELVIN*BULK_T0_KELVIN)
         lw_t3 = EMISSIVITY*BULK_STEFAN_BOLTZMANN*BULK_T0_KELVIN**3
         sens_coef = rho_air*BULK_CP_AIR*ICE_BULK_CH*wv
         lat_coef = rho_air*BULK_LATENT_SUBLIMATION*ICE_BULK_CE*wv

         sf0_ref = lw_t3*BULK_T0_KELVIN - RLDS + sens_coef*(BULK_T0_KELVIN - TAS) + &
                   lat_coef*(qsat0 - HUSS)
         dsfdt_ref = 4.0_wp*lw_t3 + sens_coef + lat_coef*dqsat0_dt

         call check(error, ice%atm_sf0(ip, jp), sf0_ref, "atm_sf0 == hand-computed SEB intercept", &
                    thr=1.0e-12_wp, rel=.true.)
         if (allocated(error)) exit checks
         call check(error, ice%atm_dsfdt(ip, jp), dsfdt_ref, "atm_dsfdt == hand-computed SEB slope", &
                    thr=1.0e-12_wp, rel=.true.)
         if (allocated(error)) exit checks
         call check(error, ice%atm_sw_dn(ip, jp), RSDS, "atm_sw_dn == rsds", thr=TOL)
         if (allocated(error)) exit checks
         call check(error, ice%atm_fprec(ip, jp), PRSN, "atm_fprec == prsn", thr=TOL)
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_ice_seb_linearization

end module test_ocean_bulk_forcing
