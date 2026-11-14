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
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_assemble
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux, ice_ocean_heat_flux, ice_ocean_sw_flux
   use rdb_ocean_bulk_forcing, only: bulk_forcing_t, bulk_forcing_apply
   use rdb_teos10_convert, only: pt_from_ct
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
   real(wp), parameter :: READER_HEAT = 11.0_wp, READER_SALT = 5.0e-6_wp

contains

   subroutine collect_ocean_bulk_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bulk_forcing_matches_column", test_matches_column), &
                  new_unittest("bulk_forcing_cons_abs_sst_is_theta", test_cons_abs_sst), &
                  new_unittest("bulk_forcing_sst_first_live", test_sst_first_live), &
                  new_unittest("bulk_forcing_rotation_exact", test_rotation_exact), &
                  new_unittest("bulk_forcing_open_water_halves", test_open_water_halves), &
                  new_unittest("bulk_forcing_stress_fold_line", test_stress_fold_line), &
                  new_unittest("bulk_forcing_stress_no_fold_unchanged", test_stress_no_fold_line), &
                  new_unittest("ice_atm_forcing_seb_linearization", test_ice_seb_linearization), &
                  new_unittest("ice_enabled_icefree_keeps_shortwave", test_icefree_keeps_sw), &
                  new_unittest("ice_partial_cover_open_plus_transmitted_sw", test_partial_cover_sw), &
                  new_unittest("ice_and_reader_flux_both_arrive", test_reader_and_ice_both_arrive), &
                  new_unittest("ice_coupler_no_accumulation_across_window", test_no_accumulation) &
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
      ms%tracers(ms%idx_temperature)%c = SST_DEGC
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

   subroutine test_sst_first_live(error)
      !! The bulk SST is read from the first LIVE layer from the top.  Two
      !! layers; the top one has collapsed (`h = 0`) and holds a stale
      !! T = -40, the live one below holds `SST_DEGC`.  The sensible heat
      !! and longwave must equal a direct `bulk_flux_column` /
      !! `bulk_flux_longwave_up` at `SST_DEGC`, bit for bit.
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
      integer :: ip, jp

      checks: block
         call make_grid(grid)
         call make_cartesian_metrics(metrics, grid)
         call init_halo_single_rank(grid)
         ip = grid%nghost + 2
         jp = grid%nghost + 2

         ms%nz_ml = 2
         call ms%init(grid)
         ms%h_layer = H0
         ms%h_layer(:, :, 2) = 0.0_wp
         ms%tracers(ms%idx_temperature)%c = SST_DEGC
         ms%tracers(ms%idx_temperature)%c(:, :, 2) = -40.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         call make_atm(grid, atm)
         call make_bulk_forcing(grid, .true., bf)
         call ss%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)

         call run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, SST_DEGC, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar)
         call check(error, sf%q_sens(ip, jp), sens_ref, "q_sens must read SST from the first live layer", &
                    thr=0.0_wp)
         if (allocated(error)) exit checks
         call check(error, sf%q_lw(ip, jp), &
                    RLDS - bulk_flux_longwave_up(SST_DEGC + BULK_T0_KELVIN, EMISSIVITY), &
                    "q_lw must read SST from the first live layer", thr=0.0_wp)
         call ocean_halo_destroy()
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_sst_first_live

   subroutine test_cons_abs_sst(error)
      !! Under `cons_abs` the model top-layer T is a Conservative
      !! Temperature: the SST reaching `bulk_flux_column` and the longwave
      !! emission is `theta = pt_from_ct(SA, CT)` (IOC/SCOR/IAPSO 2010;
      !! McDougall & Barker 2011).  Compared with a DIRECT call at theta,
      !! and shown to differ from the unconverted CT call.  The default
      !! (`cons_abs = .false.`) stays the exact-bit result of
      !! `bulk_forcing_matches_column`.
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

      real(wp), parameter :: SA_TOP = 35.2_wp
      real(wp) :: theta, taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar
      real(wp) :: sens_ct, lw_ref, lw_ct
      integer :: ip, jp

      checks: block
         call make_grid(grid)
         call make_cartesian_metrics(metrics, grid)
         call init_halo_single_rank(grid)
         ip = grid%nghost + 2
         jp = grid%nghost + 2

         call make_ms(grid, ms)
         ms%tracers(ms%idx_salinity)%c = SA_TOP
         call make_atm(grid, atm)
         call make_bulk_forcing(grid, .true., bf)
         bf%cons_abs = .true.
         call ss%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)

         call run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

         theta = pt_from_ct(SA_TOP, SST_DEGC)
         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, theta, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar)
         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, SST_DEGC, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, sens_ct, lat_ref, evap_ref, cd, ch, ce, ustar)
         call bulk_flux_column(UAS, VAS, 0.0_wp, 0.0_wp, TAS, HUSS, PSL, theta, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_ref, tauy_ref, sens_ref, lat_ref, evap_ref, cd, ch, ce, ustar)
         lw_ref = RLDS - bulk_flux_longwave_up(theta + BULK_T0_KELVIN, EMISSIVITY)
         lw_ct = RLDS - bulk_flux_longwave_up(SST_DEGC + BULK_T0_KELVIN, EMISSIVITY)

         call check(error, abs(theta - SST_DEGC) > 1.0e-2_wp, "theta must differ from CT here")
         if (allocated(error)) exit checks
         call check(error, sf%q_sens(ip, jp), sens_ref, "q_sens uses theta", thr=1.0e-10_wp, rel=.true.)
         if (allocated(error)) exit checks
         call check(error, sf%q_lat(ip, jp), lat_ref, "q_lat uses theta", thr=1.0e-10_wp, rel=.true.)
         if (allocated(error)) exit checks
         call check(error, sf%q_lw(ip, jp), lw_ref, "q_lw uses theta", thr=1.0e-10_wp, rel=.true.)
         if (allocated(error)) exit checks
         call check(error, abs(sf%q_lw(ip, jp) - lw_ct) > 1.0e-3_wp, "q_lw must NOT use the raw CT")
         if (allocated(error)) exit checks
         call check(error, abs(sf%q_sens(ip, jp) - sens_ct) > 1.0e-6_wp, "q_sens must NOT use the raw CT")
         if (allocated(error)) exit checks

         call ocean_halo_destroy()
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_cons_abs_sst

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
         !$acc update device(metrics%angle_dx)

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

   ! =================================================================
   ! Tripolar fold line -- the meridional stress on the fold-line face
   ! =================================================================

   subroutine fold_stress_case(fold, tau_y_fold, tauy_c_last, ni, jf)
      !! Run `bulk_forcing_apply` on a (6 x 3) single-rank tile with a
      !! wind field that is a TRUE VECTOR field continuous across the
      !! tripolar fold: the 180-degree rotation of the local frame negates
      !! both components at the mirror cell `i' = ni+1-i`, so
      !! `uas(i') = -uas(i)` and `vas(i') = -vas(i)` on the last row (the
      !! stress is odd in the wind).  Returns the fold-line `tau_y` row
      !! (storage row `jf = nghost+nj+1`) and the last physical row of the
      !! cell-centred `tauy_c`.
      logical, intent(in) :: fold
      integer, intent(in) :: ni, jf
      real(wp), intent(out) :: tau_y_fold(ni), tauy_c_last(ni)

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer :: i, j
      real(wp) :: d

      call grid%init(ni, NYP, NGHOST, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call init_halo_single_rank(grid)
      bc%north_fold = fold
      call make_ms(grid, ms)
      call make_atm(grid, atm)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            d = real(i - NGHOST, wp) - 0.5_wp*real(ni + 1, wp)
            atm%uas(i, j) = 2.0_wp*d + 0.1_wp*d**3
            atm%vas(i, j) = 3.0_wp*d - 0.2_wp*d**3 + 1.0_wp*d*real(j, wp)/real(NYP, wp)
         end do
      end do
      call make_bulk_forcing(grid, .true., bf)
      call ss%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      ! uas/vas are host-set: run_apply's atm%enter_data() copies them in (no update before the map).
      call run_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)
      do i = 1, ni
         tau_y_fold(i) = ss%tau_y(NGHOST + i, jf)
         tauy_c_last(i) = bf%tauy_c(NGHOST + i, NGHOST + NYP)
      end do
      call ocean_halo_destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine fold_stress_case

   subroutine test_stress_fold_line(error)
      !! The fold-line `tau_y` face (the NORTH face of the last physical
      !! T-row, storage row `nghost+nj+1`) must be the two-sided average
      !! across the fold, the mirror cell's stress entering with the vector
      !! sign: `0.5*(tauy_c(i) - tauy_c(i'))`.  Antisymmetric by itself.
      !! Before the fix the north ghost row of `tauy_c` was never folded
      !! (zero), giving `0.5*tauy_c(i)` on the east half.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NI = 6
      real(wp) :: tau_y(NI), tc(NI), want, odd_err
      integer :: i, ip

      call fold_stress_case(.true., tau_y, tc, NI, NGHOST + NYP + 1)
      checks: block
         odd_err = 0.0_wp
         do i = 1, NI
            odd_err = max(odd_err, abs(tc(i) + tc(NI + 1 - i)))
         end do
         call check(error, odd_err < 1.0e-14_wp, .true., "precondition: tauy_c odd across the fold")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(tc)) > 1.0e-3_wp, .true., "precondition: nonzero stress")
         if (allocated(error)) exit checks
         do i = 1, NI
            ip = NI + 1 - i
            want = 0.5_wp*(tc(i) - tc(ip))
            call check(error, abs(tau_y(i) - want) < 1.0e-14_wp, .true., &
                       "fold-line tau_y == two-sided average 0.5*(tauy_c(i) - tauy_c(i'))")
            if (allocated(error)) exit checks
            call check(error, tau_y(i) + tau_y(ip), 0.0_wp, "fold-line tau_y antisymmetric", thr=0.0_wp)
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_stress_fold_line

   subroutine test_stress_no_fold_line(error)
      !! Without a fold the last-row face is the one-sided edge average against the
      !! zero ghost, `0.5*tauy_c(i, last)` -- unchanged by the fold fix.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NI = 6
      real(wp) :: tau_y(NI), tc(NI)
      integer :: i

      call fold_stress_case(.false., tau_y, tc, NI, NGHOST + NYP + 1)
      do i = 1, NI
         call check(error, tau_y(i), 0.5_wp*tc(i), "no fold: tau_y(jf) = 0.5*(tauy_c(last) + 0 ghost)", &
                    thr=0.0_wp)
         if (allocated(error)) return
      end do
   end subroutine test_stress_no_fold_line

   ! =================================================================
   ! Ice coupler x bulk consumer x data-forcing reader: write order
   ! =================================================================
   ! The real per-step sequence of `engine_step` / `engine_step_ice` /
   ! `engine_step_finalize`, at unit level: reader write (every step) ->
   ! `bulk_forcing_apply` (every step) -> on a thermo step (`mod(step,
   ! ratio) == 0`) the three ice couplers -> `ocean_surface_flux_assemble`
   ! (same gate).  Before the fix the couplers full-overwrote `q_sw` /
   ! `heat_added` / `salt_flux`, so the assembled `Q_heat` carried NO open-
   ! water shortwave and none of the reader's heat/salt.

   subroutine run_sequence(ratio, nsteps, ice_on, cover, swt, hice, sice, reader, &
                           qh, qs, comp)
      !! `comp` = [q_sw, sw_ice, heat_added, heat_ice, salt_flux, salt_ice]
      !! at the probe cell after the last step; `qh`/`qs` = `Q_heat`/`Q_salt`
      !! there after every step (carried between thermo steps, as in the
      !! engine).  `reader` writes `heat_added = READER_HEAT`,
      !! `salt_flux = READER_SALT` every step.
      integer, intent(in) :: ratio, nsteps
      logical, intent(in) :: ice_on, reader
      real(wp), intent(in) :: cover, swt, hice, sice
      real(wp), intent(out) :: qh(nsteps), qs(nsteps), comp(6)

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(atm_state_t) :: atm
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer :: ip, jp, step
      logical :: thermo

      call make_grid(grid)
      call make_cartesian_metrics(metrics, grid)
      call init_halo_single_rank(grid)
      ip = grid%nghost + 2
      jp = grid%nghost + 2
      call make_ms(grid, ms)
      ms%wet_mask = 1.0_wp
      call make_atm(grid, atm)
      call make_bulk_forcing(grid, .true., bf)
      call ss%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      sf%has_heat = .true.
      sf%has_salt = .true.
      if (ice_on) then
         ice%ncat = 1
         call ice%init(grid)
         ice%enable = .true.
         ice%part_size(:, :, 0) = 1.0_wp - cover
         ice%part_size(:, :, 1) = cover
         ice%sw_thru_diag = swt
         ice%heat_flux_diag = hice
         ice%salt_flux_diag = sice
      end if

      !$acc enter data copyin(atm, ms, bf, ss, sf)
      call atm%enter_data()
      call ms%enter_data()
      call bf%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      if (ice_on) then
         !$acc enter data copyin(ice)
         call ice%enter_data()
      end if

      do step = 1, nsteps
         thermo = (mod(step, ratio) == 0)
         if (reader) then
            sf%heat_added = READER_HEAT
            sf%salt_flux = READER_SALT
            !$acc update device(sf%heat_added, sf%salt_flux)
         end if
         call bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)
         if (thermo) then
            if (ice_on) then
               call ice_ocean_brine_flux(sf, ice)
               call ice_ocean_heat_flux(sf, ice)
               call ice_ocean_sw_flux(sf, ice)
            end if
            call ocean_surface_flux_assemble(grid, sf, ms, active=.true.)
         end if
         !$acc update self(sf%Q_heat, sf%Q_salt)
         qh(step) = sf%Q_heat(ip, jp)
         qs(step) = sf%Q_salt(ip, jp)
      end do

      !$acc update self(sf%q_sw, sf%sw_ice, sf%heat_added, sf%heat_ice, sf%salt_flux, sf%salt_ice)
      comp = [sf%q_sw(ip, jp), sf%sw_ice(ip, jp), sf%heat_added(ip, jp), &
              sf%heat_ice(ip, jp), sf%salt_flux(ip, jp), sf%salt_ice(ip, jp)]

      if (ice_on) then
         call ice%exit_data()
         !$acc exit data delete(ice)
      end if
      call sf%exit_data()
      call ss%exit_data()
      call bf%exit_data()
      call ms%exit_data()
      call atm%exit_data()
      !$acc exit data delete(atm, ms, bf, ss, sf)
      call ocean_halo_destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_sequence

   subroutine test_icefree_keeps_sw(error)
      !! Ice ENABLED but no ice present (all diags 0, open fraction 1): the
      !! ocean must receive exactly what it does with ice DISABLED -- the
      !! open-water shortwave included -- at `dt_therm_ratio` 1 and 3.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS = 6
      real(wp) :: qh_on(NS), qs_on(NS), qh_off(NS), qs_off(NS), c_on(6), c_off(6)
      integer :: ratio, k, ratios(2)

      ratios = [1, 3]
      do k = 1, 2
         ratio = ratios(k)
         call run_sequence(ratio, NS, .true., 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, .false., &
                           qh_on, qs_on, c_on)
         call run_sequence(ratio, NS, .false., 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, .false., &
                           qh_off, qs_off, c_off)
         call check(error, c_off(1) > 1.0_wp, .true., "precondition: bulk open-water q_sw > 0")
         if (allocated(error)) return
         call check(error, c_on(1), c_off(1), "ice enabled: q_sw survives the ice coupler", thr=0.0_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(qh_on - qh_off)), 0.0_wp, &
                    "ice-free ice-enabled Q_heat == ice-disabled Q_heat (bitwise)", thr=0.0_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(qs_on - qs_off)), 0.0_wp, &
                    "ice-free ice-enabled Q_salt == ice-disabled Q_salt (bitwise)", thr=0.0_wp)
         if (allocated(error)) return
      end do
   end subroutine test_icefree_keeps_sw

   subroutine test_partial_cover_sw(error)
      !! 30 % ice cover: the cell gets the bulk open-water shortwave
      !! (open fraction x (1 - albedo) x rsds, written by the bulk consumer)
      !! PLUS the ice-transmitted shortwave (per cell area) -- each once.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS = 3
      real(wp), parameter :: COVER = 0.3_wp, SWT = 12.5_wp
      real(wp) :: qh(NS), qs(NS), c(6), qh0(NS), qs0(NS), c0(6), sw_open

      call run_sequence(3, NS, .true., COVER, SWT, 0.0_wp, 0.0_wp, .false., qh, qs, c)
      call run_sequence(3, NS, .false., 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, .false., qh0, qs0, c0)
      sw_open = (1.0_wp - COVER)*(1.0_wp - OCEAN_ALBEDO)*RSDS
      call check(error, c(1), sw_open, "q_sw == open fraction x (1-albedo) x rsds", thr=TOL)
      if (allocated(error)) return
      call check(error, c(2), SWT, "sw_ice == transmitted shortwave", thr=0.0_wp)
      if (allocated(error)) return
      ! Every bulk band scales with the open fraction, so the whole-cell
      ! heat is (1-a)*Q_heat(ice-free) + the transmitted shortwave, once.
      call check(error, abs(qh(NS) - ((1.0_wp - COVER)*qh0(NS) + SWT)) &
                 < 1.0e-12_wp*max(abs(qh(NS)), 1.0_wp), .true., &
                 "Q_heat == (1-a)*Q_heat_icefree + transmitted shortwave")
   end subroutine test_partial_cover_sw

   subroutine test_reader_and_ice_both_arrive(error)
      !! A data-forcing reader writing `heat_added`/`salt_flux` every step
      !! AND the ice contributing `heat_flux_diag`/`salt_flux_diag`: the
      !! assembled fluxes carry BOTH, once each; neither slot is clobbered.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS = 6
      real(wp), parameter :: HICE = -7.0_wp, SICE = 3.0e-6_wp
      real(wp) :: qh(NS), qs(NS), c(6), qh_r(NS), qs_r(NS), c_r(6)
      real(wp) :: qh_i(NS), qs_i(NS), c_i(6), qh_n(NS), qs_n(NS), c_n(6)

      call run_sequence(3, NS, .true., 0.0_wp, 0.0_wp, HICE, SICE, .true., qh, qs, c)
      call run_sequence(3, NS, .true., 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, .true., qh_r, qs_r, c_r)
      call run_sequence(3, NS, .true., 0.0_wp, 0.0_wp, HICE, SICE, .false., qh_i, qs_i, c_i)
      call run_sequence(3, NS, .true., 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, .false., qh_n, qs_n, c_n)
      call check(error, c(3), READER_HEAT, "reader heat_added survives the ice coupler", thr=0.0_wp)
      if (allocated(error)) return
      call check(error, c(5), READER_SALT, "reader salt_flux survives the ice coupler", thr=0.0_wp)
      if (allocated(error)) return
      call check(error, c(4), HICE, "heat_ice == heat_flux_diag", thr=0.0_wp)
      if (allocated(error)) return
      call check(error, c(6), SICE, "salt_ice == salt_flux_diag", thr=0.0_wp)
      if (allocated(error)) return
      ! reader + ice == (reader alone) + (ice alone) - (neither), to round-off
      call check(error, abs(qh(NS) - (qh_r(NS) + qh_i(NS) - qh_n(NS))) &
                 < 1.0e-12_wp*max(abs(qh(NS)), 1.0_wp), .true., "Q_heat carries reader AND ice")
      if (allocated(error)) return
      call check(error, abs(qh(NS) - qh_n(NS) - (READER_HEAT + HICE)) &
                 < 1.0e-12_wp*max(abs(qh(NS)), 1.0_wp), .true., "Q_heat step = reader + ice heat")
      if (allocated(error)) return
      call check(error, abs(qs(NS) - qs_n(NS) - (READER_SALT + SICE)) &
                 < 1.0e-12_wp*max(abs(qs(NS)), 1.0e-6_wp), .true., "Q_salt step = reader + ice salt")
   end subroutine test_reader_and_ice_both_arrive

   subroutine test_no_accumulation(error)
      !! Constant forcing, `dt_therm_ratio` 3 over 4 windows: the assembled
      !! flux is the SAME at every thermo step (the ice slots are re-assigned
      !! before each assembly, never added to), is carried unchanged through
      !! the non-thermo steps in between, and ice + reader + bulk are all in.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NS = 12
      real(wp) :: qh(NS), qs(NS), c(6)
      integer :: w

      call run_sequence(3, NS, .true., 0.2_wp, 9.0_wp, -4.0_wp, 2.0e-6_wp, .true., qh, qs, c)
      do w = 2, NS/3
         call check(error, qh(3*w), qh(3), "Q_heat identical at every thermo step", thr=0.0_wp)
         if (allocated(error)) return
         call check(error, qs(3*w), qs(3), "Q_salt identical at every thermo step", thr=0.0_wp)
         if (allocated(error)) return
         call check(error, qh(3*w - 1), qh(3*w - 3), "Q_heat carried through non-thermo steps", &
                    thr=0.0_wp)
         if (allocated(error)) return
      end do
      call check(error, c(2), 9.0_wp, "sw_ice is the single-window value, not a running sum", thr=0.0_wp)
   end subroutine test_no_accumulation

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

      !$acc enter data copyin(atm, ms, bf, ss, sf)
      call atm%enter_data()
      call ms%enter_data()
      call bf%enter_data()
      call ss%enter_data()
      call sf%enter_data()
      if (ice%enable) then
         !$acc enter data copyin(ice)
         call ice%enter_data()
      end if

      call bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)

      !$acc update self(sf%q_sens, sf%q_lat, sf%evap, sf%q_sw, sf%q_lw)
      !$acc update self(sf%lprec, sf%fprec, sf%p_surf_atm)
      !$acc update self(bf%taux_c, bf%tauy_c, ss%tau_x, ss%tau_y)

      if (ice%enable) then
         call ice%exit_data()
         !$acc exit data delete(ice)
      end if
      call sf%exit_data()
      call ss%exit_data()
      call bf%exit_data()
      call ms%exit_data()
      call atm%exit_data()
      !$acc exit data delete(atm, ms, bf, ss, sf)
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
         !$acc enter data copyin(atm, ice)
         call atm%enter_data()
         call ice%enter_data()
         call ice_atm_forcing_from_atm(ice, atm, cfg_bulk)
         !$acc update self(ice%atm_sf0, ice%atm_dsfdt, ice%atm_sw_dn, ice%atm_fprec)
         call ice%exit_data()
         call atm%exit_data()
         !$acc exit data delete(atm, ice)

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
