!! Conservation + quiescence tests on genuinely NON-UNIFORM (spherical)
!! metrics for the M2 slice-2a curvilinear kernel conversions
!! (continuity-PPM, Coriolis-adv, barotropic substep + coupling).
!!
!! These are the tests uniform-Cartesian grids cannot provide: when every
!! metric is equal, a kernel that still secretly uses `1/dx` or `dx*dy`
!! looks correct.  On a lon-lat sector `dxT` varies with latitude, so a
!! surviving scalar-metric form breaks mass conservation here.
!!
!! Cases:
!!   * T3 — spherical mass conservation.  Closed-basin spherical sector,
!!     non-trivial h + a velocity perturbation, ~50 continuity steps:
!!     Σ h·areaT conserved to round-off.  FAILS on the old `1/dx` form.
!!   * T6 — quiescent spherical.  Resting stratified column on the
!!     spherical sector through a short ocean_dyn_step sequence:
!!     velocities stay 0 to round-off (catches metric/PGF inconsistency
!!     introduced by the curvilinear conversion).
!!   * bitwise-on-uniform sanity — the NEW continuity form on uniform
!!     Cartesian metrics preserves a non-uniform lake at rest bit-for-bit
!!     (a guard that the conservative rewrite did not loosen the
!!     constancy-preservation property the existing tests rely on).
module test_ocean_metrics_conservation
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_fill_spherical, metrics_fill_coriolis, &
                                metrics_finalize, CORIOLIS_SCHEME_PLANETARY
   use ocean_test_metrics, only: make_cartesian_metrics, make_spherical_metrics, &
                                 destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, continuity_compute_fluxes, &
                             continuity_apply_fluxes
   use rdb_coriolis_adv, only: coriolis_adv_t, coriolis_adv_compute_tendencies
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t, tracer_hdiff
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_eos, only: eos_t
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pressure_force, only: ocean_pressure_force_compute
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   implicit none
   private

   public :: collect_ocean_metrics_conservation_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   ! Spherical sector parameters: 1-degree cells starting at 30N, so the
   ! cos(lat) factor swings ~14% across the basin (dxT non-uniform).
   real(wp), parameter :: LON_W = 0.0_wp
   real(wp), parameter :: LAT_S = 30.0_wp
   real(wp), parameter :: DLON = 1.0_wp
   real(wp), parameter :: DLAT = 1.0_wp
   real(wp), parameter :: REARTH = 6.378e6_wp

contains

   subroutine collect_ocean_metrics_conservation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("T3_spherical_mass_conservation", test_spherical_mass), &
                  new_unittest("T3b_spherical_hdiff_tracer_conservation", &
                               test_spherical_hdiff_tracer), &
                  new_unittest("T6_quiescent_spherical", test_quiescent_spherical), &
                  new_unittest("M3_sector_vs_beta_plane_convergence", &
                               test_sector_vs_beta_plane), &
                  new_unittest("bitwise_uniform_lake_at_rest", test_bitwise_uniform) &
                  ]
   end subroutine collect_ocean_metrics_conservation_tests

   ! =================================================================
   ! T3 — spherical mass conservation
   ! =================================================================

   subroutine test_spherical_mass(error)
      !! Closed-basin spherical sector (16×12 cells), 2-layer h with a
      !! latitude-varying perturbation + a divergent velocity field.
      !! Step the multilayer continuity kernel 50× and require the total
      !! mass Σ h·areaT to be conserved to round-off.  Because areaT
      !! varies with latitude (and the conservative form divides flux
      !! differences by iareaT), this exercises the genuinely
      !! non-uniform metric path — the old `1/dx` divergence would leak.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      integer, parameter :: NXP = 16, NYP = 12, NZC = 2
      real(wp), parameter :: DT = 30.0_wp
      integer, parameter :: N_STEPS = 50
      integer :: i, j, k, nx, ny, step
      real(wp) :: mass0, mass1, rel
      checks: block

         call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
         call make_spherical_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, REARTH)
         ms%nz_ml = NZC
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZC)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Non-trivial h_layer per layer (still strictly positive).
         do k = 1, NZC
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = 500.0_wp + 50.0_wp*real(k, wp) + &
                                        10.0_wp*sin(0.3_wp*real(i, wp)) + &
                                        7.0_wp*cos(0.2_wp*real(j, wp))
               end do
            end do
         end do
         ! Interior divergent velocity perturbation; walls stay zero so
         ! the basin is closed (no flux across the physical boundary).
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZC
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 2, NGHOST + NXP
                  ms%u_face_x_layer(i, j, k) = 0.05_wp*sin(0.4_wp*real(i, wp))
               end do
            end do
            do j = NGHOST + 2, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP
                  ms%v_face_y_layer(i, j, k) = 0.03_wp*cos(0.5_wp*real(j, wp))
               end do
            end do
         end do

         mass0 = total_mass(ms, metrics, grid, NZC)

         call enter(ms, ct)
         do step = 1, N_STEPS
            call continuity_compute_fluxes(grid, metrics, ct, ms)
            call continuity_apply_fluxes(ms, DT)
         end do
         call leave(ms, ct)

         mass1 = total_mass(ms, metrics, grid, NZC)
         rel = abs(mass1 - mass0)/abs(mass0)
         call check(error, rel < 1.0e-13_wp, &
                    "T3: spherical Σ h·areaT not conserved to round-off")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call ct%destroy()
      call ms%destroy()
   end subroutine test_spherical_mass

   ! =================================================================
   ! T3b — spherical tracer-mass conservation under horizontal diffusion
   ! =================================================================

   subroutine test_spherical_hdiff_tracer(error)
      !! Closed-basin spherical sector (16×12 cells, 2 layers).  Seed a
      !! non-uniform salinity field and run the conservative
      !! horizontal-diffusion kernel 50×.  The total tracer mass
      !! Σ hTr·areaT must be conserved to ~1e-13: the width-weighted
      !! face flux + iareaT divergence telescopes exactly over the
      !! closed walls only if the metric form is correct (the old
      !! `1/dx`/`1/dy` form leaks on the latitude-varying areaT).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      integer, parameter :: NXP = 16, NYP = 12, NZC = 2
      real(wp), parameter :: DT = 30.0_wp, KAPPA = 50.0_wp
      integer, parameter :: N_STEPS = 50
      integer :: i, j, k, nx, ny, step, it_S
      real(wp) :: tr0, tr1, rel
      checks: block

         call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
         call make_spherical_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, REARTH)
         ms%nz_ml = NZC
         call ms%init(grid)
         call hd%init(grid, nz_ml=NZC)
         hd%kappa_h = KAPPA
         nx = grid%nx_total
         ny = grid%ny_total
         it_S = ms%idx_salinity

         ! Strictly-positive h and a non-uniform salinity field.
         do k = 1, NZC
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = 500.0_wp + 50.0_wp*real(k, wp)
                  ms%tracers(it_S)%hTr(i, j, k) = &
                     (35.0_wp + 2.0_wp*sin(0.4_wp*real(i, wp)) + &
                      1.5_wp*cos(0.3_wp*real(j, wp)))*ms%h_layer(i, j, k)
               end do
            end do
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         ! The kernel's closed walls are the ARRAY edges (F = 0 at
         ! i=1,nx+1 / j=1,ny+1), not the physical-cell boundary — the
         ! diffusion flux is non-zero across the interior/ghost edge,
         ! so the conserved domain is the full array.
         tr0 = total_tracer_full(ms, metrics, it_S, NZC)

         call enter_hd(ms, hd)
         do step = 1, N_STEPS
            call tracer_hdiff(grid, metrics, hd, ms, DT)
         end do
         call leave_hd(ms, hd)

         tr1 = total_tracer_full(ms, metrics, it_S, NZC)
         rel = abs(tr1 - tr0)/abs(tr0)
         call check(error, rel < 1.0e-13_wp, &
                    "T3b: spherical Σ hTr·areaT not conserved under hdiff")

      end block checks
      call destroy_cartesian_metrics(metrics)
      call hd%destroy()
      call ms%destroy()
   end subroutine test_spherical_hdiff_tracer

   ! =================================================================
   ! T6 — quiescent spherical
   ! =================================================================

   subroutine test_quiescent_spherical(error)
      !! Resting stratified state on the spherical sector.  Uniform per-
      !! layer thickness + uniform per-layer S/T (so the layer densities
      !! are horizontally uniform and the PGF gradient is identically
      !! zero), zero velocity, run several full ocean_dyn_step_split
      !! steps.  Velocities must stay zero to round-off: any
      !! metric/PGF/Coriolis inconsistency introduced by the curvilinear
      !! conversion injects spurious motion that a uniform grid hides.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      integer, parameter :: NXP = 12, NYP = 10
      real(wp), parameter :: H0 = 200.0_wp
      real(wp), parameter :: F_C = 7.0e-5_wp
      real(wp), parameter :: DT = 10.0_wp
      integer, parameter :: N_INNER = 5, N_OUTER = 3
      real(wp) :: S_k(NZ), T_k(NZ)
      real(wp) :: max_du, max_dv
      integer :: k, step
      checks: block

         call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         call pgf%init(grid, nz_ml=NZ)
         call hv%init(grid, nz_ml=NZ)
         call bd%init(grid, nz_ml=NZ)
         call ss%init(grid, nz_ml=NZ)
         call va%init(grid, nz_ml=NZ)
         call hd%init(grid, nz_ml=NZ)
         call vd%init(grid, nz_ml=NZ)
         call vmix%init(grid, nz_ml=NZ)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZ)

         ! Stable stratification (bottom k=1 saltiest/coldest = densest),
         ! anchored on the EOS reference so the rest state is a true PGF
         ! balance (mirrors the Cartesian stratified-rest dyn-split test).
         S_k = [eos%S_ref + 1.0_wp, eos%S_ref, eos%S_ref - 1.0_wp]
         T_k = [eos%T_ref - 2.0_wp, eos%T_ref, eos%T_ref + 2.0_wp]
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*H0
         end do
         ! Split driver needs the total column reference depth.
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         call make_spherical_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, REARTH)
         !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ms%enter_data(); call ct%enter_data(); call cor%enter_data()
         call pgf%enter_data(); call hv%enter_data(); call bd%enter_data()
         call ss%enter_data(); call va%enter_data(); call hd%enter_data()
         call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()

         do step = 1, N_OUTER
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do

         call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
         call hd%exit_data(); call va%exit_data(); call ss%exit_data()
         call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
         call cor%exit_data(); call ct%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)

         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         call check(error, max_du < 1.0e-10_wp, &
                    "T6: quiescent spherical injected u")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp, &
                    "T6: quiescent spherical injected v")

      end block checks
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine test_quiescent_spherical

   ! =================================================================
   ! M3 — sector vs equivalent beta-plane convergence
   ! =================================================================

   subroutine test_sector_vs_beta_plane(error)
      !! Deliverable-2 / design §5 M3: the SAME physical state evaluated
      !! on (a) a small spherical sector and (b) the equivalent Cartesian
      !! beta-plane must agree to a tolerance set by the metric/f
      !! variation across the sector — NOT roundoff.  The point is to
      !! catch gross metric (cos-lat) or Coriolis (f sign/scale) errors;
      !! a wrong cos(lat) or a missing beta term breaks this at the tens-
      !! of-percent level, far above the derived ~3% physical band.
      !!
      !! We compare the PGF + Coriolis tendencies on identical ICs (rather
      !! than a full dyn-step) — the geometry-sensitive operators in
      !! isolation, free of the substep/remap entanglement, which keeps
      !! the comparison interpretable and the derived tolerance honest.
      !!
      !! Setup: 2°×2° sector centred on 45N (16×16 cells, dlon=dlat=0.125°).
      !!   spherical:  metrics_fill_spherical + planetary f = 2Ω sin(lat).
      !!   beta-plane: dx = R·cos(45°)·dλ, dy = R·dφ (constants),
      !!     f0 = 2Ω sin(45°), β = 2Ω cos(45°)/R, y from the sector
      !!     centre latitude (the f-fill y_ref placed at lat_ref).
      !!
      !! Tolerance derivation:
      !!   * f field: the beta-plane is the 1st-order Taylor expansion of
      !!     2Ω sin(lat) about 45°.  The truncation error is 2nd order:
      !!     |Δf|/|f| ≈ ½·tan(45°)·(Δφ)² with |Δφ|_max = 1° = 0.01745 rad
      !!     ⇒ ~1.5e-4.  We compare the f_corner FIELDS directly and
      !!     assert agreement < 1e-3 (with margin) — this isolates the
      !!     Coriolis-parameter generalisation from the metric variation.
      !!   * metrics: spherical dxCu ∝ cos(lat); beta-plane uses the
      !!     constant cos(45°).  cos(44°)/cos(45°)-1 ≈ +2.47e-2,
      !!     cos(46°)/cos(45°)-1 ≈ -2.45e-2 ⇒ idxCu varies by up to ~2.5%
      !!     near the meridional edges.  BOTH the PGF u-tendency
      !!     (= ΔP·idxCu) AND the Coriolis u-tendency (ζ circulation +
      !!     KE gradient also carry idxCu/iareaBu) inherit this ~2.5%
      !!     band, so we assert each agrees to 3.0e-2 — well above
      !!     roundoff, far below the breakage scale of a gross cos-lat or
      !!     f error (tens of percent).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: m_sph, m_bp
      type(multilayer_state_t) :: ms_s, ms_b
      type(coriolis_adv_t) :: cor_s, cor_b
      type(ocean_pressure_force_t) :: pgf_s, pgf_b
      type(eos_t) :: eos
      integer, parameter :: NXP = 16, NYP = 16, NZC = 3
      real(wp), parameter :: LAT_REF = 45.0_wp
      real(wp), parameter :: DLL = 0.125_wp           ! dlon = dlat (deg)
      real(wp), parameter :: OMEGA = 7.292115e-5_wp
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: D2R = PI/180.0_wp
      real(wp), parameter :: F_TOL = 1.0e-3_wp        ! f field: ~beta trunc + margin
      real(wp), parameter :: TEND_TOL = 3.0e-2_wp     ! tendencies: ~cos-lat span over 2°
      real(wp) :: lat_s, f0, beta, y_ref, dx_bp, dy_bp
      real(wp) :: rel_cor, rel_pgf, denom
      real(wp), allocatable :: f_centre_dump(:, :)
      integer :: nx, ny
      checks: block

         lat_s = LAT_REF - real(NYP, wp)*0.5_wp*DLL   ! south edge latitude
         call grid%init(NXP, NYP, NGHOST, DLL, DLL)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (f_centre_dump(nx, ny), source=0.0_wp)

         ! ---- (a) spherical sector + planetary f ----
         call m_sph%init(grid)
         call metrics_fill_spherical(m_sph, grid, LON_W, lat_s, DLL, DLL, REARTH)
         call metrics_finalize(m_sph)
         call cor_s%init(grid, nz_ml=NZC)
         call metrics_fill_coriolis(m_sph, CORIOLIS_SCHEME_PLANETARY, &
                                    0.0_wp, 0.0_wp, 0.0_wp, OMEGA, grid, &
                                    cor_s%f_corner, f_centre_dump)

         ! ---- (b) equivalent Cartesian beta-plane ----
         ! dy uses R·dφ; dx uses R·cos(lat_ref)·dλ.  grid%dy is set to
         ! dy_bp so the beta-plane f-fill (which measures y = (j-…)·grid%dy)
         ! lines up with the spherical latitude spacing.
         dy_bp = REARTH*DLL*D2R
         dx_bp = REARTH*cos(LAT_REF*D2R)*DLL*D2R
         grid%dx = dx_bp
         grid%dy = dy_bp
         f0 = 2.0_wp*OMEGA*sin(LAT_REF*D2R)
         beta = 2.0_wp*OMEGA*cos(LAT_REF*D2R)/REARTH
         ! y_ref places f0 at lat_ref: the corner fill uses
         ! y = (j-1-ng)·dy, and lat_ref corner row is j = ng+1+NYP/2.
         y_ref = real(NYP, wp)*0.5_wp*dy_bp
         call m_bp%init(grid)
         call metrics_fill_cartesian(m_bp, grid, dx_bp, dy_bp)
         call metrics_finalize(m_bp)
         call cor_b%init(grid, nz_ml=NZC)
         call cor_b%set_beta_plane(grid, f0, beta, y_ref)

         ! ---- identical baroclinic ICs on both states ----
         call eos%init(grid)
         call pgf_s%init(grid, nz_ml=NZC)
         call pgf_b%init(grid, nz_ml=NZC)
         ms_s%nz_ml = NZC; ms_b%nz_ml = NZC
         call ms_s%init(grid); call ms_b%init(grid)
         call seed_baroclinic_ic(ms_s, eos, grid, NZC)
         call seed_baroclinic_ic(ms_b, eos, grid, NZC)

         ! ---- run EOS + PGF + Coriolis on the spherical state ----
         call eval_pgf_cor(grid, m_sph, eos, pgf_s, cor_s, ms_s)
         call eval_pgf_cor(grid, m_bp, eos, pgf_b, cor_b, ms_b)

         ! ---- (1) f field: isolated from metrics, ~beta-truncation band ----
         denom = maxval(abs(cor_b%f_corner))
         if (denom <= 0.0_wp) denom = 1.0_wp
         rel_cor = maxval(abs(cor_s%f_corner - cor_b%f_corner))/denom
         call check(error, rel_cor < F_TOL, &
                    "M3: planetary f vs beta-plane f rel-diff exceeds the "// &
                    "beta-truncation band (~1.5e-4 expected)")
         if (allocated(error)) exit checks
         ! and the f fields must NOT be identical — beta-plane is only an
         ! approximation, so a small (but non-zero) difference is expected.
         call check(error, rel_cor > 1.0e-7_wp, &
                    "M3: planetary f bit-identical to beta-plane f — "// &
                    "sphericity of f not exercised")
         if (allocated(error)) exit checks

         ! ---- (2) PGF u-tendency: ~cos-lat metric band ----
         denom = maxval(abs(pgf_s%dpdx_face%data))
         if (denom <= 0.0_wp) denom = 1.0_wp
         rel_pgf = maxval(abs(pgf_s%dpdx_face%data - pgf_b%dpdx_face%data))/denom
         call check(error, rel_pgf < TEND_TOL, &
                    "M3: PGF u-tendency sector vs beta-plane rel-diff "// &
                    "exceeds the cos-lat metric-variation band")
         if (allocated(error)) exit checks
         ! Must not be byte-identical (metric distinction really present).
         call check(error, rel_pgf > 1.0e-6_wp, &
                    "M3: sector and beta-plane PGF identical — metric "// &
                    "variation not exercised")
         if (allocated(error)) exit checks

         ! ---- (3) Coriolis u-tendency: same ~cos-lat metric band ----
         denom = maxval(abs(cor_s%pv_flux_x%data))
         if (denom <= 0.0_wp) denom = 1.0_wp
         rel_cor = maxval(abs(cor_s%pv_flux_x%data - cor_b%pv_flux_x%data))/denom
         call check(error, rel_cor < TEND_TOL, &
                    "M3: Coriolis u-tendency sector vs beta-plane rel-diff "// &
                    "exceeds the metric-variation band")

      end block checks
      call pgf_s%destroy(); call pgf_b%destroy()
      call cor_s%destroy(); call cor_b%destroy()
      call eos%destroy()
      call ms_s%destroy(); call ms_b%destroy()
      call m_sph%destroy(); call m_bp%destroy()
      if (allocated(f_centre_dump)) deallocate (f_centre_dump)
   end subroutine test_sector_vs_beta_plane

   subroutine seed_baroclinic_ic(ms, eos, grid, nz)
      !! Latitude-varying stratified S/T + a smooth velocity field, so
      !! both the PGF (density gradients) and the Coriolis (advective +
      !! planetary) tendencies are non-trivial and identical between the
      !! two grids (same array indices ⇒ same physical column).
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      integer :: i, j, k, nx, ny
      real(wp) :: s_ijk, t_ijk
      nx = grid%nx_total
      ny = grid%ny_total
      ms%h_layer = 200.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               ! Stable stratification (densest at the bed, k=1) + BOTH a
               ! zonal and a meridional S/T gradient — the zonal gradient
               ! drives a non-trivial east-face PGF (∂/∂x), which is the
               ! component the cos-lat metric difference acts on.
               s_ijk = eos%S_ref + real(nz - k, wp)*0.5_wp + &
                       0.02_wp*real(i, wp) + 0.01_wp*real(j, wp)
               t_ijk = eos%T_ref - real(nz - k, wp)*1.0_wp - &
                       0.03_wp*real(i, wp) - 0.02_wp*real(j, wp)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = s_ijk*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_ijk*ms%h_layer(i, j, k)
            end do
         end do
      end do
      ! Smooth interior velocity for the Coriolis advection term.
      do k = 1, nz
         do j = NGHOST + 1, ny - NGHOST
            do i = NGHOST + 2, nx - NGHOST
               ms%u_face_x_layer(i, j, k) = 0.10_wp*sin(0.3_wp*real(j, wp))
            end do
            do i = NGHOST + 1, nx - NGHOST
               ms%v_face_y_layer(i, j, k) = 0.08_wp*cos(0.25_wp*real(i, wp))
            end do
         end do
      end do
   end subroutine seed_baroclinic_ic

   subroutine eval_pgf_cor(grid, metrics, eos, pgf, cor, ms)
      !! Device round-trip: EOS → PGF → Coriolis on one (metrics, state)
      !! pair, leaving the tendency buffers updated on the host.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(eos_t), intent(inout) :: eos
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(coriolis_adv_t), intent(inout) :: cor
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
      !$acc enter data copyin(ms, pgf, cor)
      call ms%enter_data(); call pgf%enter_data(); call cor%enter_data()
      call ocean_eos_compute(eos, ms)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
      !$acc update self(pgf%dpdx_face%data, cor%pv_flux_x%data)
      call cor%exit_data(); call pgf%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, pgf, cor)
      call metrics%exit_data()
      !$acc exit data delete(metrics)
   end subroutine eval_pgf_cor

   ! =================================================================
   ! Bitwise-on-uniform sanity
   ! =================================================================

   subroutine test_bitwise_uniform(error)
      !! The conservative-curvilinear continuity rewrite must still
      !! preserve a non-uniform lake at rest bit-for-bit on uniform
      !! Cartesian metrics (zero velocity ⇒ zero flux ⇒ h unchanged).
      !! This guards the constancy-preservation property at the same
      !! 1e-14 tolerance the existing continuity tests use — the
      !! rewrite must not have loosened it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), allocatable :: h_ic(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: max_diff
      checks: block

         call grid%init(12, 10, NGHOST, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (h_ic(nx, ny, NZ))
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  h_ic(i, j, k) = 100.0_wp + 5.0_wp*real(k, wp) + &
                                  0.5_wp*real(i, wp) + 0.3_wp*real(j, wp)
               end do
            end do
         end do
         ms%h_layer = h_ic
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call enter(ms, ct)
         call continuity_compute_fluxes(grid, metrics, ct, ms)
         call continuity_apply_fluxes(ms, 1.0_wp)
         call leave(ms, ct)

         max_diff = maxval(abs(ms%h_layer - h_ic))
         call check(error, max_diff < 1.0e-14_wp, &
                    "bitwise uniform: lake at rest drifted under new form")

      end block checks
      call destroy_cartesian_metrics(metrics)
      if (allocated(h_ic)) deallocate (h_ic)
      call ct%destroy()
      call ms%destroy()
   end subroutine test_bitwise_uniform

   ! =================================================================
   ! Helpers
   ! =================================================================

   function total_mass(ms, metrics, grid, nz) result(m)
      !! Σ over interior cells + layers of h_layer·areaT (host-side).
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp) :: m
      integer :: i, j, k, ng, nxp, nyp
      ng = grid%nghost
      nxp = grid%nx_phys
      nyp = grid%ny_phys
      m = 0.0_wp
      do k = 1, nz
         do j = ng + 1, ng + nyp
            do i = ng + 1, ng + nxp
               m = m + ms%h_layer(i, j, k)*metrics%areaT(i, j)
            end do
         end do
      end do
   end function total_mass

   function total_tracer_full(ms, metrics, it, nz) result(m)
      !! Σ over the FULL array (incl. ghosts) of hTr·areaT.  The hdiff
      !! kernel's closed walls are the array edges, so the conserved
      !! domain is the whole array, not the physical interior.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: it, nz
      real(wp) :: m
      integer :: i, j, k
      m = 0.0_wp
      do k = 1, nz
         do j = 1, size(metrics%areaT, 2)
            do i = 1, size(metrics%areaT, 1)
               m = m + ms%tracers(it)%hTr(i, j, k)*metrics%areaT(i, j)
            end do
         end do
      end do
   end function total_tracer_full

   subroutine enter_hd(ms, hd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hd)
      call hd%enter_data()
   end subroutine enter_hd

   subroutine leave_hd(ms, hd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      call hd%exit_data()
      !$acc exit data delete(hd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine leave_hd

   subroutine enter(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms, ct)
      call ms%enter_data()
      call ct%enter_data()
      ! metrics already entered by make_*_metrics.
   end subroutine enter

   subroutine leave(ms, ct)
      !! Unmap ms + ct only.  The metrics slot is left mapped + allocated
      !! so the caller can still read `metrics%areaT` on the host to
      !! compute the post-run mass; the caller destroys it afterwards.
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct)
   end subroutine leave

end module test_ocean_metrics_conservation
