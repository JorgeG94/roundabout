!! Validation gates for C-grid EVP sea-ice dynamics (PR 5): `rdb_ice_evp`.
!!
!! Harness: `grid%init` + `ocean_test_metrics::make_cartesian_metrics`
!! (all-wet Cartesian) + raw test-allocated arrays for the core seam
!! (`ice_evp_dynamics`), plus a live `ocean_state_t` for the restart-
!! roundtrip and disabled-bitident gates.  GPU mem:separate discipline
!! throughout (canonical template: `test_open_boundary_out_closes` in
!! `test_ocean_conservation_salt_heat.F90`): unconditional
!! `!$acc enter data copyin(...)` before every kernel call, `!$acc update
!! device` after every host re-fill, `!$acc update self` before every host
!! assertion, `exit data` at teardown.
!!
!! All goldens verified against `tmp_local_artifacts/evp_pr5_goldens.log`
!! + `proto_evp_validate.py`; SIS2 defaults unless stated (EC=2,
!! p0=2.75e4, c0=20, cdw=3.24e-3, rho_ocean=1030, del_sh_min_scale=2,
!! tdamp=-0.2, evp_sub_steps=432).  Walls are made by the core's own
!! non-periodic edge policy (ghost mask 0) — tests do NOT allocate land
!! cells unless testing interior land.
module test_ocean_ice_evp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_evp, only: ice_evp_dynamics, ice_evp_params_t, ice_evp_mi_ratio_point, &
                          evp_truncate_final_impl
   use rdb_ice_state, only: ocean_sea_ice_t, evp_workspace_t
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux, ice_ocean_stress_resume_apply, &
                                    ice_ocean_stress_cleanup
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data, &
                              ocean_state_restart_write, ocean_state_restart_read
   use rdb_decomp, only: decomp_t, decomp_init
   implicit none
   private

   public :: collect_ocean_ice_evp_tests

   integer, parameter :: NGHOST = 3
   real(wp), parameter :: PI = 3.14159265358979323846_wp

   ! SIS2 defaults (module-wide unless a gate overrides).
   real(wp), parameter :: P0_DEFAULT = 2.75e4_wp
   real(wp), parameter :: C0_DEFAULT = 20.0_wp
   real(wp), parameter :: EC_DEFAULT = 2.0_wp
   real(wp), parameter :: CDW_DEFAULT = 3.24e-3_wp
   real(wp), parameter :: RHO_OCEAN_DEFAULT = 1030.0_wp
   real(wp), parameter :: DEL_SH_MIN_SCALE_DEFAULT = 2.0_wp
   real(wp), parameter :: TDAMP_DEFAULT = -0.2_wp
   integer, parameter :: EVP_SUB_STEPS_DEFAULT = 432
   real(wp), parameter :: CDW_TINY = 1.0e-30_wp
      !! PR 62: ballistic-limit trick (drag_u -> 0) so the wind term alone
      !! is analytically isolable — see `test_a_face_ballistic_scaling`.

contains

   subroutine collect_ocean_ice_evp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nansen_free_drift", test_nansen_free_drift), &
                  new_unittest("wall_arrest", test_wall_arrest), &
                  new_unittest("shear_2d", test_shear_2d), &
                  new_unittest("limit_cadence", test_limit_cadence), &
                  new_unittest("ice_margin", test_ice_margin), &
                  new_unittest("cavitating", test_cavitating), &
                  new_unittest("inertial_coriolis", test_inertial_coriolis), &
                  new_unittest("mi_ratio_harmonic", test_mi_ratio_harmonic), &
                  new_unittest("gather_ncat_equivalence", test_gather_ncat_equivalence), &
                  new_unittest("tau_coupling", test_tau_coupling), &
#ifndef RDB_NO_NETCDF
                  new_unittest("evp_restart_roundtrip", test_evp_restart_roundtrip), &
#endif
                  new_unittest("disabled_bitident", test_disabled_bitident), &
                  new_unittest("a_face_full_cover_bitident", test_a_face_full_cover_bitident), &
                  new_unittest("a_face_ballistic_scaling", test_a_face_ballistic_scaling), &
                  new_unittest("a_face_momentum_closes", test_a_face_momentum_closes), &
                  new_unittest("a_face_wind_only_would_leak", test_a_face_wind_only_would_leak), &
                  new_unittest("a_face_no_ghost_drift", test_a_face_no_ghost_drift), &
                  new_unittest("a_face_disabled_bitident", test_a_face_disabled_bitident), &
                  new_unittest("trunc_disabled_bitident", test_trunc_disabled_bitident), &
                  new_unittest("trunc_bound_donor_asymmetry", test_trunc_bound_donor_asymmetry), &
                  new_unittest("trunc_final_clips_and_counts", test_trunc_final_clips_and_counts), &
                  new_unittest("trunc_below_bound_bitident", test_trunc_below_bound_bitident), &
                  new_unittest("trunc_uses_transport_dt", test_trunc_uses_transport_dt), &
                  new_unittest("project_ci_zero_divergence_bitident", &
                               test_project_ci_zero_divergence_bitident), &
                  new_unittest("project_ci_stiffens_convergence", &
                               test_project_ci_stiffens_convergence), &
                  new_unittest("project_ci_extreme_divergence_finite", &
                               test_project_ci_extreme_divergence_finite) &
                  ]
   end subroutine collect_ocean_ice_evp_tests

   ! =====================================================================
   ! Shared params builder
   ! =====================================================================

   pure function default_params(evp_sub_steps, ec, cdw, tdamp, a_face_stress, &
                                cfl_trunc, cfl_trunc_dyn_its, project_ci) result(par)
      integer, intent(in), optional :: evp_sub_steps
      real(wp), intent(in), optional :: ec, cdw, tdamp
      logical, intent(in), optional :: a_face_stress
      real(wp), intent(in), optional :: cfl_trunc
      logical, intent(in), optional :: cfl_trunc_dyn_its, project_ci
      type(ice_evp_params_t) :: par

      par%p0 = P0_DEFAULT
      par%c0 = C0_DEFAULT
      par%ec = EC_DEFAULT
      par%cdw = CDW_DEFAULT
      par%rho_ocean = RHO_OCEAN_DEFAULT
      par%del_sh_min_scale = DEL_SH_MIN_SCALE_DEFAULT
      par%tdamp = TDAMP_DEFAULT
      par%evp_sub_steps = EVP_SUB_STEPS_DEFAULT
      par%a_face_stress = .false.
      par%cfl_trunc = 0.0_wp
      par%cfl_trunc_dyn_its = .false.
      par%project_ci = .false.
      if (present(evp_sub_steps)) par%evp_sub_steps = evp_sub_steps
      if (present(ec)) par%ec = ec
      if (present(cdw)) par%cdw = cdw
      if (present(tdamp)) par%tdamp = tdamp
      if (present(a_face_stress)) par%a_face_stress = a_face_stress
      if (present(cfl_trunc)) par%cfl_trunc = cfl_trunc
      if (present(cfl_trunc_dyn_its)) par%cfl_trunc_dyn_its = cfl_trunc_dyn_its
      if (present(project_ci)) par%project_ci = project_ci
   end function default_params

   ! =====================================================================
   ! Gate 1: Nansen free drift
   ! =====================================================================

   subroutine test_nansen_free_drift(error)
      !! Double-periodic 4x4, uniform mi/ci, f=0, uniform tau_ax sweep.
      !! Golden |u| = sqrt(tau/(rho_ocean*cdw)); face-uniformity; str_d
      !! spatial uniformity; v==0; fxoc~=tau_ax at steady state.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 48
      integer, parameter :: N_TAU = 4
      real(wp), parameter :: TAU_VALUES(N_TAU) = [0.05_wp, 0.1_wp, 0.2_wp, 0.5_wp]
      real(wp), parameter :: U_ANALYTIC(N_TAU) = [0.122403513677564_wp, 0.173104709124932_wp, &
                                                  0.244807027355129_wp, 0.387073896828676_wp]

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, itau, n, i, j
      real(wp) :: u_mean, u_min, u_max, str_min, str_max, worst_tau_err

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%init(nx, ny)
      call ws%enter_data()

      checks: block
         do itau = 1, N_TAU
            mis = 0.0_wp
            mice = 0.0_wp
            ci = 0.0_wp
            uo = 0.0_wp
            vo = 0.0_wp
            tau_ax = 0.0_wp
            tau_ay = 0.0_wp
            ui = 0.0_wp
            vi = 0.0_wp
            str_d = 0.0_wp
            str_t = 0.0_wp
            str_s = 0.0_wp
            f_corner = 0.0_wp

            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP
                  mis(i, j) = 3.0_wp*ICE_RHO_ICE
                  mice(i, j) = 3.0_wp*ICE_RHO_ICE
                  ci(i, j) = 1.0_wp
               end do
            end do
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP + 1
                  tau_ax(i, j) = TAU_VALUES(itau)
               end do
            end do

            !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
            !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

            do n = 1, N_OUTER
               call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                     tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                     fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
            end do

            !$acc update self(ui, vi, str_d, fxoc)
            !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
            !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

            u_mean = sum(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)) &
                     /real((NXP + 1)*NYP, wp)
            u_min = minval(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP))
            u_max = maxval(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP))
            str_min = minval(str_d(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP))
            str_max = maxval(str_d(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP))

            call check(error, abs(u_mean - U_ANALYTIC(itau)) <= 1.0e-10_wp*U_ANALYTIC(itau), &
                       "nansen_free_drift: |u| mismatch")
            if (allocated(error)) exit checks
            call check(error, (u_max - u_min) <= 1.0e-9_wp, "nansen_free_drift: face spread")
            if (allocated(error)) exit checks
            call check(error, (str_max - str_min) <= 1.0e-6_wp*par%p0, &
                       "nansen_free_drift: str_d not spatially uniform")
            if (allocated(error)) exit checks
            call check(error, all(abs(vi(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP + 1)) &
                                  <= 1.0e-14_wp), "nansen_free_drift: v /= 0")
            if (allocated(error)) exit checks

            worst_tau_err = 0.0_wp
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP + 1
                  worst_tau_err = max(worst_tau_err, abs(fxoc(i, j) - tau_ax(i, j)))
               end do
            end do
            call check(error, worst_tau_err <= 1.0e-8_wp*TAU_VALUES(itau), &
                       "nansen_free_drift: fxoc != tau_ax at steady state")
            if (allocated(error)) exit checks
         end do
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nansen_free_drift

   ! =====================================================================
   ! Gate 2: wall arrest
   ! =====================================================================

   subroutine test_wall_arrest(error)
      !! x-WALLS + y-PERIODIC channel, nx=20 phys, ny=4 phys, dx=dy=25000,
      !! uniform mi/ci, tau_ax=+0.1, f=0, 100 outer steps of 3600s.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 20, NYP = 4
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 100
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: pres_mice, l_domain, x_i, sigma_ana, worst_ramp_err, max_u, u_free
      real(wp) :: str_d_row_spread, str_d_j0, str_d_jk

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         ui = 0.0_wp
         vi = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .false., .true., ws)
         end do

         !$acc update self(ui, vi, str_d, str_t, str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         pres_mice = P0_DEFAULT/ICE_RHO_ICE*exp(0.0_wp)*MI_CONST
         l_domain = real(NXP, wp)*DX

         worst_ramp_err = 0.0_wp
         do i = NGHOST + 1, NGHOST + NXP
            x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
            sigma_ana = -0.5_wp*pres_mice - TAU_A*(x_i - 0.5_wp*l_domain)
            do j = NGHOST + 1, NGHOST + NYP
               worst_ramp_err = max(worst_ramp_err, &
                                    abs((str_d(i, j) + str_t(i, j)) - sigma_ana))
            end do
         end do
         call check(error, worst_ramp_err/(TAU_A*l_domain) <= 2.0e-5_wp, &
                    "wall_arrest: stress ramp mismatch")
         if (allocated(error)) exit checks

         max_u = maxval(abs(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)))
         u_free = sqrt(TAU_A/(RHO_OCEAN_DEFAULT*CDW_DEFAULT))
         call check(error, max_u <= 1.0e-2_wp*u_free, "wall_arrest: max|u| too large")
         if (allocated(error)) exit checks

         str_d_row_spread = 0.0_wp
         do i = NGHOST + 1, NGHOST + NXP
            str_d_j0 = str_d(i, NGHOST + 1)
            do j = NGHOST + 1, NGHOST + NYP
               str_d_jk = str_d(i, j)
               str_d_row_spread = max(str_d_row_spread, abs(str_d_jk - str_d_j0))
            end do
         end do
         call check(error, str_d_row_spread <= 1.0e-12_wp, &
                    "wall_arrest: str_d not y-uniform (periodic-y wrap)")
         if (allocated(error)) exit checks

         ! "Exactly 0" up to float noise: 43200 subcycles of semi-implicit
         ! relaxation on quantities that are analytically 0 (no y-variation,
         ! no v-forcing anywhere) accumulate ~1e-13-scale round-off, not
         ! literal bitwise zero.
         call check(error, maxval(abs(str_s)) <= 1.0e-9_wp*par%p0, "wall_arrest: str_s /= 0")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vi)) <= 1.0e-12_wp, "wall_arrest: v /= 0")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_wall_arrest

   ! =====================================================================
   ! Gate 3: 2-D pure shear (corner machinery + str_s clamp)
   ! =====================================================================

   subroutine test_shear_2d(error)
      !! Fully-walled box nx=ny=20, prescribed uniform simple shear,
      !! re-imposed every subcycle (n_sub=1, dt=1800/432).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXY = 20
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: MI_CONST = 2.0_wp*ICE_RHO_ICE
      real(wp), parameter :: GAMMA_DOT = 1.0e-5_wp
      real(wp), parameter :: DT_SLOW_NOMINAL = 1800.0_wp
      integer, parameter :: N_SUB_NOMINAL = 432
      integer, parameter :: N_REIMPOSE = N_SUB_NOMINAL*10
      integer, parameter :: MARGIN = 3

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j, i_mid, j_mid
      real(wp) :: pres_mice_val, pressure, str_d_exp, str_t_exp, str_s_exp
      real(wp) :: worst_str_d_ratio, worst_str_t_ratio, worst_str_s_ratio
      real(wp) :: pres_avg
      real(wp) :: dt_eff

      call grid%init(NXY, NXY, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(evp_sub_steps=1)
      dt_eff = DT_SLOW_NOMINAL/real(N_SUB_NOMINAL, wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         ui = 0.0_wp
         vi = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do

         call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_REIMPOSE
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, dt_eff, par, .false., .false., ws)
            ui = 0.0_wp
            call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)
            !$acc update device(ui, vi)
         end do

         !$acc update self(str_d, str_t, str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         pres_mice_val = P0_DEFAULT/ICE_RHO_ICE
         pressure = pres_mice_val*MI_CONST
         str_d_exp = -0.5_wp*pressure
         str_t_exp = 0.0_wp
         str_s_exp = 0.5_wp*pressure/EC_DEFAULT

         i_mid = NGHOST + NXY/2
         j_mid = NGHOST + NXY/2

         call check(error, abs(str_d(i_mid, j_mid) - str_d_exp) <= 1.0e-13_wp*pressure, &
                    "shear_2d: str_d closed-form mismatch")
         if (allocated(error)) exit checks
         call check(error, abs(str_t(i_mid, j_mid) - str_t_exp) <= 1.0e-13_wp*pressure, &
                    "shear_2d: str_t closed-form mismatch")
         if (allocated(error)) exit checks
         call check(error, abs(str_s(i_mid, j_mid) - str_s_exp) <= 1.0e-13_wp*pressure, &
                    "shear_2d: str_s closed-form mismatch")
         if (allocated(error)) exit checks

         ! Yield bound everywhere.
         worst_str_d_ratio = 0.0_wp
         worst_str_t_ratio = 0.0_wp
         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               worst_str_d_ratio = max(worst_str_d_ratio, -str_d(i, j)/pressure)
               worst_str_t_ratio = max(worst_str_t_ratio, &
                                       abs(EC_DEFAULT*str_t(i, j))/(0.5_wp*pressure))
            end do
         end do
         call check(error, worst_str_d_ratio <= 1.0_wp + 1.0e-12_wp, &
                    "shear_2d: str_d yield-bound violation")
         if (allocated(error)) exit checks
         call check(error, worst_str_t_ratio <= 1.0_wp + 1.0e-12_wp, &
                    "shear_2d: str_t yield-bound violation")
         if (allocated(error)) exit checks

         worst_str_s_ratio = 0.0_wp
         do j = NGHOST + 1 + MARGIN, NGHOST + NXY + 1 - MARGIN
            do i = NGHOST + 1 + MARGIN, NGHOST + NXY + 1 - MARGIN
               pres_avg = 0.25_wp*((pres_mice_val*mice(i - 1, j - 1) + &
                                    pres_mice_val*mice(i, j)) + &
                                   (pres_mice_val*mice(i, j - 1) + &
                                    pres_mice_val*mice(i - 1, j)))
               worst_str_s_ratio = max(worst_str_s_ratio, &
                                       abs(EC_DEFAULT*str_s(i, j))/(0.5_wp*pres_avg))
            end do
         end do
         call check(error, worst_str_s_ratio <= 1.0_wp + 1.0e-12_wp, &
                    "shear_2d: str_s yield-bound violation")
         if (allocated(error)) exit checks

         ! Requirement-(4) discriminator (prototype row bookkeeping: prototype
         ! I=10/i=10 (0-based) -> rdb index NGHOST+11 (0-based -> 1-based
         ! +1, then +NGHOST for the ghost prefix) for BOTH the cell index
         ! (str_t) and the corner index (str_s) -- see SPEC SS1. J=0 <-> south
         ! physical boundary corner row (rdb NGHOST+1), J=1 <-> first
         ! interior corner row (rdb NGHOST+2); j=0 <-> first physical cell
         ! row (rdb NGHOST+1).
         call check(error, abs(str_s(NGHOST + 11, NGHOST + 2) - 7199.4493087499_wp) &
                    <= 1.0e-6_wp*7199.4493087499_wp, &
                    "shear_2d: requirement-4 str_s(first interior row) mismatch")
         if (allocated(error)) exit checks
         call check(error, abs(str_s(NGHOST + 11, NGHOST + 1)) <= 1.0e-9_wp, &
                    "shear_2d: requirement-4 str_s(wall corner row) /= 0")
         if (allocated(error)) exit checks
         call check(error, abs(str_t(NGHOST + 11, NGHOST + 1) - (-6147.6502176361_wp)) &
                    <= 1.0e-6_wp*6147.6502176361_wp, &
                    "shear_2d: requirement-4 str_t(first physical row) mismatch")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_shear_2d

   subroutine impose_shear(vi, gamma_dot, dx, nxy, nghost, nx, ny)
      !! `ui=0` everywhere (caller pre-zeroed); `vi(i,J) = gamma_dot*x_i`
      !! for interior rows, 0 at the physical N/S boundary v-faces.
      integer, intent(in) :: nxy, nghost, nx, ny
      real(wp), intent(in) :: gamma_dot, dx
      real(wp), intent(inout) :: vi(nx, ny + 1)
      integer :: i, j
      real(wp) :: x_i

      do j = 1, ny + 1
         do i = 1, nx
            vi(i, j) = 0.0_wp
         end do
      end do
      do i = nghost + 1, nghost + nxy
         x_i = real(i - nghost - 1, wp)*dx
         do j = nghost + 2, nghost + nxy
            vi(i, j) = gamma_dot*x_i
         end do
      end do
   end subroutine impose_shear

   subroutine impose_divergence(ui, alpha, dx, nxy, nghost, nx, ny)
      !! PR 36 (cases 8/9): a uniform 1-D divergence field --
      !! `ui(i,j) = alpha*x_i` at every physical u-face (`vi==0`
      !! everywhere, caller's responsibility). On a uniform Cartesian
      !! grid this gives `sh_Dd = du/dx + dv/dy = alpha` exactly at any
      !! T-cell whose two bracketing u-faces are both in the filled
      !! range (true for any cell away from the domain edge -- callers
      !! check an interior cell, e.g. `i_mid = nghost + nxy/2`).
      !! `alpha < 0` => CONVERGENT (u decreases with x => flow points
      !! inward => `sh_Dd < 0`); `alpha > 0` => DIVERGENT.
      integer, intent(in) :: nxy, nghost, nx, ny
      real(wp), intent(in) :: alpha, dx
      real(wp), intent(inout) :: ui(nx + 1, ny)
      integer :: i, j
      real(wp) :: x_i

      do j = 1, ny
         do i = 1, nx + 1
            ui(i, j) = 0.0_wp
         end do
      end do
      do i = nghost + 1, nghost + nxy + 1
         x_i = real(i - nghost - 1, wp)*dx
         do j = nghost + 1, nghost + nxy
            ui(i, j) = alpha*x_i
         end do
      end do
   end subroutine impose_divergence

   function run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui_in, vi_in, dt_slow, par, &
                                  i_mid, j_mid) result(str_d_val)
      !! PR 36 (cases 8/9) shared helper: run ONE fresh `ice_evp_dynamics`
      !! call (`str_d`/`str_t`/`str_s` start at 0) on the caller's
      !! prescribed `ui_in`/`vi_in` and return `str_d(i_mid,j_mid)`
      !! host-side. Self-contained GPU mem:separate span (owns its own
      !! `evp_workspace_t`, maps/unmaps every array it touches).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
      real(wp), intent(in) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), intent(in) :: uo(:, :), vo(:, :), tau_ax(:, :), tau_ay(:, :)
      real(wp), intent(in) :: ui_in(:, :), vi_in(:, :)
      real(wp), intent(in) :: dt_slow
      type(ice_evp_params_t), intent(in) :: par
      integer, intent(in) :: i_mid, j_mid
      real(wp) :: str_d_val

      type(evp_workspace_t) :: ws
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :)
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      allocate (ui, source=ui_in)
      allocate (vi, source=vi_in)
      allocate (str_d(nx, ny), source=0.0_wp)
      allocate (str_t(nx, ny), source=0.0_wp)
      allocate (str_s(nx + 1, ny + 1), source=0.0_wp)
      allocate (fxoc(nx + 1, ny), source=0.0_wp)
      allocate (fyoc(nx, ny + 1), source=0.0_wp)

      call ws%init(nx, ny)
      call ws%enter_data()
      !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

      call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                            tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                            fxoc, fyoc, dt_slow, par, .false., .false., ws)

      !$acc update self(str_d)
      !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
      !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%exit_data()
      call ws%destroy()

      str_d_val = str_d(i_mid, j_mid)
      deallocate (ui, vi, str_d, str_t, str_s, fxoc, fyoc)
   end function run_one_substep_str_d

   ! =====================================================================
   ! Gate 4: limit_stresses cadence discriminator
   ! =====================================================================

   subroutine test_limit_cadence(error)
      !! Requirements (2)+(3) discriminator: opposing wind bands drive
      !! str_s past the yield surface between calls (per-call clamp, not
      !! per-subcycle); a subsequent zero-forcing call must decay back
      !! under the bound.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXY = 20
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: TAU_A = 0.2_wp
      real(wp), parameter :: MI_CONST = 2.0_wp*ICE_RHO_ICE
      real(wp), parameter :: DT_SLOW = 1800.0_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j
      real(wp) :: worst_ratio_call1, worst_ratio_final

      call grid%init(NXY, NXY, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         ui = 0.0_wp
         vi = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         ! Opposing wind bands: + for the upper half of physical rows.
         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY + 1
               if (j > NGHOST + NXY/2) then
                  tau_ax(i, j) = TAU_A
               else
                  tau_ax(i, j) = -TAU_A
               end if
            end do
         end do

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .false., .false., ws)
         !$acc update self(str_s, mice)
         worst_ratio_call1 = worst_str_s_ec_ratio(str_s, mice, nx, ny, NGHOST, NXY)
         call check(error, worst_ratio_call1 > 1.1_wp, &
                    "limit_cadence: call 1 must overshoot the yield bound (per-call clamp)")
         if (allocated(error)) exit checks

         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .false., .false., ws)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .false., .false., ws)

         ! Zero forcing/velocity, then ONE call with n_sub=1.
         ui = 0.0_wp
         vi = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         !$acc update device(ui, vi, tau_ax, tau_ay)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW/432.0_wp, &
                               default_params(evp_sub_steps=1), .false., .false., ws)

         !$acc update self(str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         worst_ratio_final = worst_str_s_ec_ratio(str_s, mice, nx, ny, NGHOST, NXY)
         call check(error, worst_ratio_final <= 1.0_wp + 1.0e-12_wp, &
                    "limit_cadence: post-clamp ratio exceeds the yield bound")
         if (allocated(error)) exit checks
         call check(error, worst_ratio_final < 1.0_wp, &
                    "limit_cadence: post-clamp ratio did not decay below 1")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_limit_cadence

   function worst_str_s_ec_ratio(str_s, mice, nx, ny, nghost, nxy) result(worst)
      !! Matches `ice_limit_stresses`'s masked-area-weighted pres_avg
      !! EXACTLY (uniform areaT on this Cartesian grid, so area weights
      !! cancel to a plain mask-weighted average): a wall-adjacent corner
      !! sees only its <=2 physically-wet T-cell neighbours, never the
      !! ghost-row zero-mice cells diluting the average (that dilution is
      !! a distinct bug from what the corner clamp actually computes).
      integer, intent(in) :: nx, ny, nghost, nxy
      real(wp), intent(in) :: str_s(nx + 1, ny + 1)
      real(wp), intent(in) :: mice(nx, ny)
      real(wp) :: worst
      integer :: i, j, i_worst, j_worst
      real(wp) :: pres_avg, pres_mice_val, ratio, sum_mask
      real(wp) :: mask_sw, mask_se, mask_nw, mask_ne
      integer :: i_lo, i_hi, j_lo, j_hi

      i_lo = nghost + 1
      i_hi = nghost + nxy
      j_lo = nghost + 1
      j_hi = nghost + nxy

      pres_mice_val = P0_DEFAULT/ICE_RHO_ICE
      worst = 0.0_wp
      i_worst = -1
      j_worst = -1
      do j = nghost + 1, nghost + nxy + 1
         do i = nghost + 1, nghost + nxy + 1
            mask_sw = wet_mask_at(i - 1, j - 1, i_lo, i_hi, j_lo, j_hi)
            mask_se = wet_mask_at(i, j - 1, i_lo, i_hi, j_lo, j_hi)
            mask_nw = wet_mask_at(i - 1, j, i_lo, i_hi, j_lo, j_hi)
            mask_ne = wet_mask_at(i, j, i_lo, i_hi, j_lo, j_hi)
            sum_mask = (mask_sw + mask_ne) + (mask_nw + mask_se)
            if (sum_mask > 0.0_wp) then
               pres_avg = ((mask_sw*pres_mice_val*mice(i - 1, j - 1) + &
                            mask_ne*pres_mice_val*mice(i, j)) + &
                           (mask_nw*pres_mice_val*mice(i - 1, j) + &
                            mask_se*pres_mice_val*mice(i, j - 1)))/sum_mask
               if (pres_avg > 0.0_wp) then
                  ratio = abs(EC_DEFAULT*str_s(i, j))/(0.5_wp*pres_avg)
                  if (ratio > worst) then
                     worst = ratio
                     i_worst = i
                     j_worst = j
                  end if
               end if
            end if
         end do
      end do
   end function worst_str_s_ec_ratio

   pure function wet_mask_at(i, j, i_lo, i_hi, j_lo, j_hi) result(m)
      integer, intent(in) :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: m
      if (i >= i_lo .and. i <= i_hi .and. j >= j_lo .and. j <= j_hi) then
         m = 1.0_wp
      else
         m = 0.0_wp
      end if
   end function wet_mask_at

   ! =====================================================================
   ! Gate 5: ice margin (no ice-edge special-casing)
   ! =====================================================================

   subroutine test_ice_margin(error)
      !! x-walled + y-periodic channel, 20x4: west 10 cells iced, east 10
      !! ice-free. Ice-free cells relax str_d to (near) 0 emergently.
      !! Knob-off (legacy, `a_face_stress=.false.`) path. See
      !! `test_a_face_no_ghost_drift` for the `&ocean_ice_nml
      !! a_face_stress=.true.` behaviour at the same ice edge.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 20, NYP = 4
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 20
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      integer, parameter :: I_SEED = NGHOST + 15
         !! An ice-free cell (east half starts at physical col 11) seeded
         !! with a spurious str_d.

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: u_free, worst_far_east_err, edge_speed

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         ui = 0.0_wp
         vi = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = 0.0_wp

         ! West 10 physical cells iced; east 10 ice-free.
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + 10
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do
         ! Seed a spurious str_d in an ice-free cell.
         str_d(I_SEED, NGHOST + 2) = -1000.0_wp

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .false., .true., ws)
         end do

         !$acc update self(ui, vi, str_d)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         ! Ice-free cells that started at 0 stay exactly 0 (the seeded
         ! cell I_SEED is checked separately below via decay).
         block
            logical :: ok_zero
            ok_zero = .true.
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 11, NGHOST + NXP
                  if (i == I_SEED .and. j == NGHOST + 2) cycle
                  if (str_d(i, j) /= 0.0_wp) ok_zero = .false.
               end do
            end do
            call check(error, ok_zero, &
                       "ice_margin: an ice-free cell that started at 0 drifted off 0")
         end block
         if (allocated(error)) exit checks

         call check(error, abs(str_d(I_SEED, NGHOST + 2)) < 1.0e-6_wp, &
                    "ice_margin: seeded str_d did not decay")
         if (allocated(error)) exit checks

         call check(error,.not. any(ieee_is_nan_grid(str_d, nx, ny)), &
                    "ice_margin: NaN in str_d")
         if (allocated(error)) exit checks
         call check(error,.not. any(ieee_is_nan_grid(ui, nx + 1, ny)), &
                    "ice_margin: NaN in ui")
         if (allocated(error)) exit checks

         ! "Far east": the last INTERIOR face before the east wall (face
         ! NGHOST+NXP+1 is the wall itself, masked to 0 by mask_u -- not
         ! a meaningful probe of the ice-free interior's drift speed).
         ! PR 62: this is the UNWEIGHTED (a_face_stress=.false.) form's
         ! massless-slab GHOST DRIFT (see the PR-62 plan / rdb_ice_evp's D7
         ! docstring) -- NOT physical ice-free-water drift.  A massless
         ! slab pushed by the full wind and braked by the full quadratic
         ! drag has a finite terminal velocity == the Nansen free-drift
         ! speed, regenerated every substep, even though there is no ice
         ! mass at this face.  `test_a_face_no_ghost_drift` is the knob-on
         ! twin where the same face is exactly `uo` (== 0 here).
         u_free = sqrt(TAU_A/(RHO_OCEAN_DEFAULT*CDW_DEFAULT))
         worst_far_east_err = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            worst_far_east_err = max(worst_far_east_err, &
                                     abs(ui(NGHOST + NXP, j) - u_free))
         end do
         call check(error, worst_far_east_err <= 1.0e-3_wp, &
                    "ice_margin (a_face_stress=.false., legacy): ice-free face "// &
                    "not at the massless-slab ghost drift")
         if (allocated(error)) exit checks

         edge_speed = ui(NGHOST + 11, NGHOST + 2)
         call check(error, edge_speed > 0.5_wp*u_free, &
                    "ice_margin: iced pack does not drift east across the open margin")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_ice_margin

   pure elemental function ieee_is_nan_grid(x, nx, ny) result(isnan_val)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: x
      logical :: isnan_val
      isnan_val = (x /= x)
   end function ieee_is_nan_grid

   ! =====================================================================
   ! Gate 6: cavitating fluid (EC=0)
   ! =====================================================================

   subroutine test_cavitating(error)
      !! EC=0 wall-arrest config: str_t/str_s stay exactly 0; str_d stays
      !! within the yield bound; the ramp still holds (looser tolerance);
      !! max|u| stays small.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 20, NYP = 4
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 100
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: pres_mice, l_domain, x_i, sigma_ana, worst_ramp_err, max_u, u_free
      real(wp) :: worst_str_d_ratio

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(ec=0.0_wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         ui = 0.0_wp
         vi = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .false., .true., ws)
         end do

         !$acc update self(ui, str_d, str_t, str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, all(str_t == 0.0_wp), "cavitating: str_t /= 0")
         if (allocated(error)) exit checks
         call check(error, all(str_s == 0.0_wp), "cavitating: str_s /= 0")
         if (allocated(error)) exit checks

         pres_mice = P0_DEFAULT/ICE_RHO_ICE
         worst_str_d_ratio = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               worst_str_d_ratio = max(worst_str_d_ratio, &
                                       -str_d(i, j)/(pres_mice*MI_CONST))
            end do
         end do
         call check(error, worst_str_d_ratio <= 1.0_wp + 1.0e-12_wp, &
                    "cavitating: str_d yield-bound violation")
         if (allocated(error)) exit checks

         l_domain = real(NXP, wp)*DX
         worst_ramp_err = 0.0_wp
         do i = NGHOST + 1, NGHOST + NXP
            x_i = (real(i - NGHOST, wp) - 0.5_wp)*DX
            sigma_ana = -0.5_wp*pres_mice*MI_CONST - TAU_A*(x_i - 0.5_wp*l_domain)
            do j = NGHOST + 1, NGHOST + NYP
               worst_ramp_err = max(worst_ramp_err, abs((str_d(i, j) + str_t(i, j)) - sigma_ana))
            end do
         end do
         call check(error, worst_ramp_err/(TAU_A*l_domain) <= 5.0e-5_wp, &
                    "cavitating: stress ramp mismatch")
         if (allocated(error)) exit checks

         max_u = maxval(abs(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)))
         u_free = sqrt(TAU_A/(RHO_OCEAN_DEFAULT*CDW_DEFAULT))
         call check(error, max_u <= 2.0e-2_wp*u_free, "cavitating: max|u| too large")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_cavitating

   ! =====================================================================
   ! Gate 7: inertial Coriolis
   ! =====================================================================

   subroutine test_inertial_coriolis(error)
      !! Double-periodic 4x4, uniform mi/ci, f_corner=1.4e-4, cdw~0, no
      !! tau. Collapses exactly to a scalar quasi-implicit rotation
      !! recurrence; assert every face matches at every outer step.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 12
      integer, parameter :: N_SUB = 432
      real(wp), parameter :: F_CORIOLIS = 1.4e-4_wp
      real(wp), parameter :: CDW_TINY = 1.0e-30_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      real(wp), parameter :: U_IC = 0.3_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j, sub
      real(wp) :: u_scalar, v_scalar, u_new, v_new, dt_sub, dtf, k_factor, m_neglect
      real(wp) :: worst_err

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(cdw=CDW_TINY)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp
         mice = 0.0_wp
         ci = 0.0_wp
         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.0_wp
         tau_ay = 0.0_wp
         str_d = 0.0_wp
         str_t = 0.0_wp
         str_s = 0.0_wp
         f_corner = F_CORIOLIS

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         ui = U_IC
         vi = 0.0_wp

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         u_scalar = U_IC
         v_scalar = 0.0_wp
         dt_sub = DT_SLOW/real(N_SUB, wp)
         dtf = dt_sub*F_CORIOLIS
         m_neglect = ICE_RHO_ICE*1.0e-30_wp
         k_factor = MI_CONST/(MI_CONST + m_neglect)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .true., ws)

            do sub = 1, N_SUB
               u_new = (u_scalar + dtf*v_scalar)/(1.0_wp + dtf*dtf)*k_factor
               v_new = (v_scalar - dtf*u_scalar)/(1.0_wp + dtf*dtf)*k_factor
               u_scalar = u_new
               v_scalar = v_new
            end do

            !$acc update self(ui, vi)
            worst_err = 0.0_wp
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP + 1
                  worst_err = max(worst_err, abs(ui(i, j) - u_scalar))
               end do
            end do
            do j = NGHOST + 1, NGHOST + NYP + 1
               do i = NGHOST + 1, NGHOST + NXP
                  worst_err = max(worst_err, abs(vi(i, j) - v_scalar))
               end do
            end do
            call check(error, worst_err <= 1.0e-12_wp*max(abs(u_scalar), abs(v_scalar), 1.0_wp), &
                       "inertial_coriolis: face mismatch vs scalar recurrence")
            if (allocated(error)) exit checks
         end do

         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_inertial_coriolis

   ! =====================================================================
   ! Gate 8: mi_ratio_A_q harmonic-mean unit test
   ! =====================================================================

   subroutine test_mi_ratio_harmonic(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: AREA = 25000.0_wp*25000.0_wp
      real(wp), parameter :: M_NEGLECT = 905.0_wp*1.0e-30_wp
      real(wp) :: m2, m4, ratio, sum_area

      m2 = M_NEGLECT*M_NEGLECT
      m4 = m2*m2

      checks: block
         ! Uniform mis=1810 -> ratio*sum_area == 1.
         sum_area = 4.0_wp*AREA
         ratio = ice_evp_mi_ratio_point(1810.0_wp, 1810.0_wp, 1810.0_wp, 1810.0_wp, &
                                        1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        AREA, AREA, AREA, AREA, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        m2, m4)
         call check(error, abs(ratio*sum_area - 1.0_wp) <= 1.0e-14_wp, &
                    "mi_ratio_harmonic: uniform case")
         if (allocated(error)) exit checks

         ! Quad [905, 1810, 2715, 3620] SW/SE/NW/NE -> ratio*sum_area = 0.896.
         ratio = ice_evp_mi_ratio_point(905.0_wp, 1810.0_wp, 2715.0_wp, 3620.0_wp, &
                                        1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        AREA, AREA, AREA, AREA, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        m2, m4)
         call check(error, abs(ratio*sum_area - 0.896000000000_wp) <= 1.0e-12_wp, &
                    "mi_ratio_harmonic: quad case")
         if (allocated(error)) exit checks

         ! MIZ quad [9.05, 2715, 2715, 2715] -> ratio*sum_area = 0.889873257116.
         ratio = ice_evp_mi_ratio_point(9.05_wp, 2715.0_wp, 2715.0_wp, 2715.0_wp, &
                                        1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        AREA, AREA, AREA, AREA, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
                                        m2, m4)
         call check(error, abs(ratio*sum_area - 0.889873257116_wp) <= 1.0e-9_wp, &
                    "mi_ratio_harmonic: MIZ quad case")
         if (allocated(error)) exit checks

         ! Straight coast: only the SW/SE u-face-below wet, 2 wet T-cells.
         ! mask_q=0, (mask_u_below+mask_u_above)+(mask_v_left+mask_v_right)
         ! <= 1.5 -> ratio = 1/sum_area over the masked area (2 wet cells).
         sum_area = 2.0_wp*AREA
         ratio = ice_evp_mi_ratio_point(1810.0_wp, 1810.0_wp, 0.0_wp, 0.0_wp, &
                                        1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                        AREA, AREA, AREA, AREA, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                                        m2, m4)
         call check(error, abs(ratio - 1.0_wp/sum_area) <= 1.0e-14_wp*abs(1.0_wp/sum_area), &
                    "mi_ratio_harmonic: straight-coast branch")
         if (allocated(error)) exit checks

         ! All-land -> 0.
         ratio = ice_evp_mi_ratio_point(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                        0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                        AREA, AREA, AREA, AREA, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                        m2, m4)
         call check(error, ratio == 0.0_wp, "mi_ratio_harmonic: all-land branch")
      end block checks
   end subroutine test_mi_ratio_harmonic

   ! =====================================================================
   ! Gate 9: ncat gather equivalence
   ! =====================================================================

   subroutine test_gather_ncat_equivalence(error)
      !! Same physical ice as ncat=1 lumped vs ncat=3 fully-covered:
      !! gathered mis/mice/ci agree bitwise, so one EVP call gives
      !! bit-identical u_ice.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      real(wp), parameter :: M_TOTAL = 3.0_wp*ICE_RHO_ICE
      real(wp), parameter :: S_TOTAL = 0.5_wp*ICE_RHO_ICE

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      type(ocean_sea_ice_t) :: ice1, ice3
      real(wp), allocatable :: mis1(:, :), mice1(:, :), ci1(:, :)
      real(wp), allocatable :: mis3(:, :), mice3(:, :), ci3(:, :)
      real(wp), allocatable :: uo(:, :), vo(:, :), tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui1(:, :), vi1(:, :), str_d1(:, :), str_t1(:, :), str_s1(:, :)
      real(wp), allocatable :: ui3(:, :), vi3(:, :), str_d3(:, :), str_t3(:, :), str_s3(:, :)
      real(wp), allocatable :: fxoc1(:, :), fyoc1(:, :), fxoc3(:, :), fyoc3(:, :)
      real(wp), allocatable :: f_corner(:, :)
      integer :: nx, ny, i, j

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      ice1%enable = .true.
      ice1%ncat = 1
      ice1%nk_ice = 1
      call ice1%init(grid)
      ice3%enable = .true.
      ice3%ncat = 3
      ice3%nk_ice = 1
      call ice3%init(grid)

      allocate (mis1(nx, ny), mice1(nx, ny), ci1(nx, ny))
      allocate (mis3(nx, ny), mice3(nx, ny), ci3(nx, ny))
      allocate (uo(nx + 1, ny), vo(nx, ny + 1), tau_ax(nx + 1, ny), tau_ay(nx, ny + 1))
      allocate (ui1(nx + 1, ny), vi1(nx, ny + 1))
      allocate (str_d1(nx, ny), str_t1(nx, ny), str_s1(nx + 1, ny + 1))
      allocate (ui3(nx + 1, ny), vi3(nx, ny + 1))
      allocate (str_d3(nx, ny), str_t3(nx, ny), str_s3(nx + 1, ny + 1))
      allocate (fxoc1(nx + 1, ny), fyoc1(nx, ny + 1), fxoc3(nx + 1, ny), fyoc3(nx, ny + 1))
      allocate (f_corner(nx + 1, ny + 1))

      checks: block
         ice1%m_ice = 0.0_wp
         ice1%m_snow = 0.0_wp
         ice3%m_ice = 0.0_wp
         ice3%m_snow = 0.0_wp
         ice3%part_size = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               ice1%m_ice(i, j, 1) = M_TOTAL
               ice1%m_snow(i, j, 1) = S_TOTAL
               ! ncat=3 fully covered: Σ part=1, Σ part*m = M, uniform
               ! per-cat intensive masses.
               ice3%part_size(i, j, 0) = 0.0_wp
               ice3%part_size(i, j, 1) = 1.0_wp/3.0_wp
               ice3%part_size(i, j, 2) = 1.0_wp/3.0_wp
               ice3%part_size(i, j, 3) = 1.0_wp/3.0_wp
               ice3%m_ice(i, j, :) = M_TOTAL
               ice3%m_snow(i, j, :) = S_TOTAL
            end do
         end do

         uo = 0.0_wp
         vo = 0.0_wp
         tau_ax = 0.1_wp
         tau_ay = 0.0_wp
         ui1 = 0.0_wp
         vi1 = 0.0_wp
         str_d1 = 0.0_wp
         str_t1 = 0.0_wp
         str_s1 = 0.0_wp
         ui3 = 0.0_wp
         vi3 = 0.0_wp
         str_d3 = 0.0_wp
         str_t3 = 0.0_wp
         str_s3 = 0.0_wp
         f_corner = 0.0_wp

         !$acc enter data copyin(ice1)
         call ice1%enter_data()
         !$acc enter data copyin(ice3)
         call ice3%enter_data()
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis1, mice1, ci1, mis3, mice3, ci3)
         !$acc enter data copyin(uo, vo, tau_ax, tau_ay, f_corner)
         !$acc enter data copyin(ui1, vi1, str_d1, str_t1, str_s1, fxoc1, fyoc1)
         !$acc enter data copyin(ui3, vi3, str_d3, str_t3, str_s3, fxoc3, fyoc3)

         call gather_and_run(grid, metrics, f_corner, ice1, uo, vo, tau_ax, tau_ay, &
                             ui1, vi1, str_d1, str_t1, str_s1, fxoc1, fyoc1, &
                             mis1, mice1, ci1, DT_SLOW, par, ws)
         call gather_and_run(grid, metrics, f_corner, ice3, uo, vo, tau_ax, tau_ay, &
                             ui3, vi3, str_d3, str_t3, str_s3, fxoc3, fyoc3, &
                             mis3, mice3, ci3, DT_SLOW, par, ws)

         !$acc update self(mis1, mice1, ci1, mis3, mice3, ci3, ui1, ui3)
         !$acc exit data delete(mis1, mice1, ci1, mis3, mice3, ci3)
         !$acc exit data delete(uo, vo, tau_ax, tau_ay, f_corner)
         !$acc exit data delete(ui1, vi1, str_d1, str_t1, str_s1, fxoc1, fyoc1)
         !$acc exit data delete(ui3, vi3, str_d3, str_t3, str_s3, fxoc3, fyoc3)
         call ice1%exit_data()
         !$acc exit data delete(ice1)
         call ice3%exit_data()
         !$acc exit data delete(ice3)

         call check(error, all(mis1 == mis3), "gather_ncat_equivalence: mis mismatch")
         if (allocated(error)) exit checks
         call check(error, all(mice1 == mice3), "gather_ncat_equivalence: mice mismatch")
         if (allocated(error)) exit checks
         call check(error, all(ci1 == ci3), "gather_ncat_equivalence: ci mismatch")
         if (allocated(error)) exit checks
         call check(error, all(ui1 == ui3), "gather_ncat_equivalence: u_ice not bit-identical")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call ice1%destroy()
      call ice3%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_gather_ncat_equivalence

   subroutine gather_and_run(grid, metrics, f_corner, ice, uo, vo, tau_ax, tau_ay, &
                             ui, vi, str_d, str_t, str_s, fxoc, fyoc, mis, mice, ci, &
                             dt_slow, par, ws)
      use rdb_ice_state, only: ice_cell_concentration_impl
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in) :: f_corner(:, :)
      type(ocean_sea_ice_t), intent(in) :: ice
      real(wp), intent(in) :: uo(:, :), vo(:, :), tau_ax(:, :), tau_ay(:, :)
      real(wp), intent(inout) :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), intent(inout) :: fxoc(:, :), fyoc(:, :)
      real(wp), intent(inout) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), intent(in) :: dt_slow
      type(ice_evp_params_t), intent(in) :: par
      type(evp_workspace_t), intent(inout) :: ws

      call ice_cell_concentration_impl(metrics%wet_T, ice%part_size, ice%m_ice, ice%m_snow, &
                                       mis, mice, ci, ice%ncat, ice%nx_total, ice%ny_total)
      call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, tau_ax, tau_ay, &
                            ui, vi, str_d, str_t, str_s, fxoc, fyoc, dt_slow, par, &
                            .true., .true., ws)
   end subroutine gather_and_run

   ! =====================================================================
   ! Gate 10: tau coupling blend
   ! =====================================================================

   subroutine test_tau_coupling(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_stress_t) :: stress
      integer :: nx, ny, i, j

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ice%enable = .true.
      ice%ncat = 1
      ice%nk_ice = 1
      call ice%init(grid)
      call stress%init(grid, nz_ml=1)

      checks: block
         ! Full cover (ci=1 everywhere): tau_x == fxoc bitwise.
         ice%m_ice = 0.0_wp
         ice%m_snow = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               ice%m_ice(i, j, 1) = 3.0_wp*ICE_RHO_ICE
            end do
         end do
         ice%tau_a_x = 0.3_wp
         ice%tau_a_y = 0.0_wp
         ice%fxoc = 0.7_wp
         ice%fyoc = -0.2_wp
         stress%tau_x = 0.0_wp
         stress%tau_y = 0.0_wp

         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(stress)
         call stress%enter_data()
         ! metrics%wet_T is read on-device by ice_cell_concentration_impl inside
         ! ice_ocean_stress_flux -> it MUST be device-present (mem:separate).
         ! `make_cartesian_metrics` above ALREADY mapped it (parent + arrays),
         ! so do NOT re-map here: a second copyin(metrics)+enter_data() would
         ! leave the metrics arrays at presentcount 2, and if any check below
         ! trips `exit checks` the paired in-block unmap is skipped while the
         ! post-block `destroy_cartesian_metrics` removes only ONE reference ->
         ! idxCu/idyCu leak on-device with their host block freed, then collide
         ! with a later subtest's fresh allocation (present-table FATAL).

         call ice_ocean_stress_flux(metrics, stress, ice)

         ! Interior faces only (NGHOST+2 .. NGHOST+NXP): both neighbouring
         ! T-cells are physical + fully covered there.  The two boundary
         ! faces (NGHOST+1, NGHOST+NXP+1) straddle a ghost T-cell that was
         ! never given m_ice ⇒ ci=0 there ⇒ a genuine (and correct) 0.5
         ! blend, not a full-cover a_u=1 — those are exercised by the
         ! half-cover case below instead.
         !$acc update self(stress%tau_x, stress%tau_y)
         call check(error, all(abs(stress%tau_x(NGHOST + 2:NGHOST + NXP, &
                                                NGHOST + 1:NGHOST + NYP) - 0.7_wp) &
                               <= epsilon(1.0_wp)*10.0_wp), &
                    "tau_coupling: full-cover tau_x /= fxoc")
         if (allocated(error)) exit checks
         call check(error, all(abs(stress%tau_y(NGHOST + 1:NGHOST + NXP, &
                                                NGHOST + 2:NGHOST + NYP) - (-0.2_wp)) &
                               <= epsilon(1.0_wp)*10.0_wp), &
                    "tau_coupling: full-cover tau_y /= fyoc")
         if (allocated(error)) exit checks

         ! Zero cover: tau_x == tau_a_x bitwise (every face — a=0 exactly
         ! for EVERY face when ci=0 everywhere, including boundaries).
         ice%m_ice = 0.0_wp
         !$acc update device(ice%m_ice)
         call ice_ocean_stress_flux(metrics, stress, ice)
         !$acc update self(stress%tau_x, stress%tau_y)
         call check(error, all(stress%tau_x(NGHOST + 1:NGHOST + NXP + 1, &
                                            NGHOST + 1:NGHOST + NYP) == 0.3_wp), &
                    "tau_coupling: zero-cover tau_x /= tau_a_x")
         if (allocated(error)) exit checks

         ! Half cover: exact 0.5/0.5 blend on an interior u-face where BOTH
         ! neighbours are iced on one side, ice-free the other, giving a_u=0.5.
         ice%m_ice = 0.0_wp
         ice%m_ice(NGHOST + 1, NGHOST + 2, 1) = 3.0_wp*ICE_RHO_ICE
         ! (NGHOST+2, NGHOST+2) stays ice-free -> a_u at face (NGHOST+2,NGHOST+2) = 0.5.
         !$acc update device(ice%m_ice)
         call ice_ocean_stress_flux(metrics, stress, ice)
         !$acc update self(stress%tau_x, stress%tau_y)
         call check(error, abs(stress%tau_x(NGHOST + 2, NGHOST + 2) - &
                               (0.5_wp*0.3_wp + 0.5_wp*0.7_wp)) <= 1.0e-13_wp, &
                    "tau_coupling: half-cover blend mismatch")
         if (allocated(error)) exit checks

         ! Resume apply reproduces the last-written blend bitwise (PR 63 —
         ! the carry replaces the old bitwise-reconstruct-on-host gate:
         ! ice_ocean_stress_flux already mirrored stress%tau_x/y into
         ! ice%tau_ocn_x/y above, so a plain COPY back must match exactly).
         block
            real(wp), allocatable :: tau_x_device(:, :), tau_y_device(:, :)
            allocate (tau_x_device, source=stress%tau_x)
            allocate (tau_y_device, source=stress%tau_y)
            stress%tau_x = -999.0_wp
            stress%tau_y = -999.0_wp
            !$acc update device(stress%tau_x, stress%tau_y)
            ! Pull the device-mirrored ice_tau_mirror_impl output (the
            ! LAST value ice_ocean_stress_flux wrote) to the host — the
            ! resume apply is host-only and reads ice%tau_ocn_x/y as plain
            ! host arrays (component arrays only, never `update self(ice)`
            ! — commit 72152870).
            !$acc update self(ice%tau_ocn_x, ice%tau_ocn_y)
            ! `ice_ocean_stress_resume_apply` is a PLAIN-HOST routine (it
            ! writes stress%tau_x/tau_y on the host, never on the device —
            ! that is the whole point of the resume path). Do NOT
            ! `!$acc update self` after it: that would pull the stale
            ! device values (-999 above) back over the host apply result
            ! and the compare below would fail.
            call ice_ocean_stress_resume_apply(stress, ice)
            call check(error, all(stress%tau_x == tau_x_device), &
                       "tau_coupling: resume apply != device blend (tau_x)")
            if (allocated(error)) exit checks
            call check(error, all(stress%tau_y == tau_y_device), &
                       "tau_coupling: resume apply != device blend (tau_y)")
         end block

         !$acc exit data delete(stress)
         call stress%exit_data()
         !$acc exit data delete(ice)
         call ice%exit_data()
         ! metrics was mapped once by `make_cartesian_metrics`; it is unmapped
         ! once by `destroy_cartesian_metrics` below (which runs on EVERY exit
         ! path, in or out of `checks`). Do not unmap it here.
      end block checks

      call ice_ocean_stress_cleanup()
      call ice%destroy()
      call stress%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_tau_coupling

   ! =====================================================================
   ! Gate 11: restart round-trip
   ! =====================================================================

   subroutine test_evp_restart_roundtrip(error)
      !! Ice-enabled ocean_state with dynamics=.true.: fill the seven PR-5
      !! fields with non-trivial values, registry write -> clear -> read,
      !! bit-exact (model on test_ocean_restart.F90).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 6, NYP = 5
      real(wp), parameter :: DX = 1000.0_wp
      character(len=*), parameter :: FN = "test_ocean_ice_evp_rt.nc"

      type(hgrid_t) :: grid
      type(ocean_state_t) :: a_state, b_state
      type(decomp_t) :: decomp
      integer :: nx, ny, i, j
      real(wp) :: t_read
      integer :: step_read
      real(wp), allocatable :: u_ice_a(:, :), v_ice_a(:, :)
      real(wp), allocatable :: str_d_a(:, :), str_t_a(:, :), str_s_a(:, :)
      real(wp), allocatable :: fxoc_a(:, :), fyoc_a(:, :)

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx = grid%nx_total
      ny = grid%ny_total
      call decomp_init(decomp, NXP, NYP, 1, 1, 0)
      call cleanup_test_file(FN)

      a_state%multilayer%nz_ml = 1
      a_state%ice%enable = .true.
      a_state%ice%ncat = 1
      a_state%ice%nk_ice = 1
      a_state%ice%dynamics = .true.
      call a_state%init(grid)
      call metrics_seed_cartesian(a_state, grid)

      checks: block
         do j = 1, ny
            do i = 1, nx
               a_state%ice%u_ice(min(i, nx), j) = 0.01_wp*real(i, wp) + 0.001_wp*real(j, wp)
               a_state%ice%str_d(i, j) = -100.0_wp - real(i + j, wp)
               a_state%ice%str_t(i, j) = 5.0_wp + real(i, wp)
               a_state%ice%fxoc(min(i, nx), j) = 0.001_wp*real(i - j, wp)
               a_state%ice%fyoc(i, min(j, ny)) = -0.002_wp*real(i, wp)
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               a_state%ice%v_ice(i, j) = 0.02_wp*real(i, wp) - 0.003_wp*real(j, wp)
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx + 1
               a_state%ice%str_s(i, j) = 200.0_wp + real(i, wp) - real(j, wp)
            end do
         end do
         allocate (u_ice_a, source=a_state%ice%u_ice)
         allocate (v_ice_a, source=a_state%ice%v_ice)
         allocate (str_d_a, source=a_state%ice%str_d)
         allocate (str_t_a, source=a_state%ice%str_t)
         allocate (str_s_a, source=a_state%ice%str_s)
         allocate (fxoc_a, source=a_state%ice%fxoc)
         allocate (fyoc_a, source=a_state%ice%fyoc)

         call ocean_state_enter_data(a_state)
         !$acc update device(a_state%ice%u_ice, a_state%ice%v_ice)
         !$acc update device(a_state%ice%str_d, a_state%ice%str_t, a_state%ice%str_s)
         !$acc update device(a_state%ice%fxoc, a_state%ice%fyoc)

         call ocean_state_restart_write(a_state, grid, decomp, FN, 0.0_wp, 0)
         call ocean_state_exit_data(a_state)
         call a_state%destroy()

         b_state%multilayer%nz_ml = 1
         b_state%ice%enable = .true.
         b_state%ice%ncat = 1
         b_state%ice%nk_ice = 1
         b_state%ice%dynamics = .true.
         call b_state%init(grid)
         call metrics_seed_cartesian(b_state, grid)
         call ocean_state_restart_read(b_state, grid, decomp, FN, t_read, step_read)
         call ocean_state_enter_data(b_state)
         !$acc update self(b_state%ice%u_ice, b_state%ice%v_ice)
         !$acc update self(b_state%ice%str_d, b_state%ice%str_t, b_state%ice%str_s)
         !$acc update self(b_state%ice%fxoc, b_state%ice%fyoc)

         call check(error, all(b_state%ice%u_ice == u_ice_a), "restart: u_ice mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%v_ice == v_ice_a), "restart: v_ice mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%str_d == str_d_a), "restart: str_d mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%str_t == str_t_a), "restart: str_t mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%str_s == str_s_a), "restart: str_s mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%fxoc == fxoc_a), "restart: fxoc mismatch")
         if (allocated(error)) exit checks
         call check(error, all(b_state%ice%fyoc == fyoc_a), "restart: fyoc mismatch")
         if (allocated(error)) exit checks

         call ocean_state_exit_data(b_state)
      end block checks

      call b_state%destroy()
      call cleanup_test_file(FN)
   end subroutine test_evp_restart_roundtrip

   subroutine metrics_seed_cartesian(state, grid)
      !! Fill the ocean_state's own metrics slot with a Cartesian fill
      !! (mirrors `ocean_test_metrics::make_cartesian_metrics` without a
      !! separate metrics object — `state%metrics` is already allocated
      !! by `state%init`).
      use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine metrics_seed_cartesian

   subroutine cleanup_test_file(fn)
      character(len=*), intent(in) :: fn
      integer :: u, ios
      character(len=256) :: msg
      logical :: exists
      inquire (file=fn, exist=exists)
      if (exists) then
         open (newunit=u, file=fn, status="old", action="readwrite", iostat=ios, iomsg=msg)
         if (ios == 0) close (u, status="delete")
      end if
   end subroutine cleanup_test_file

   ! =====================================================================
   ! Gate 12: disabled bit-identity
   ! =====================================================================

   subroutine test_disabled_bitident(error)
      !! dynamics=.false.: ice_evp_step is a no-op on an uninitialised or
      !! dynamics-off slot (kernel-level defence-in-depth guard).
      use rdb_ice_evp, only: ice_evp_step
      use rdb_multilayer_state, only: multilayer_state_t
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 1000.0_wp
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_sea_ice_t) :: ice_uninit, ice_off
      type(multilayer_state_t) :: ms
      type(ice_evp_params_t) :: par
      real(wp), allocatable :: f_corner(:, :)
      integer :: nx, ny

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params()

      ! Case A: uninitialised slot (is_init == .false.).
      allocate (f_corner(nx + 1, ny + 1), source=0.0_wp)
      !$acc enter data copyin(f_corner)
      call ice_evp_step(grid, metrics, f_corner, ice_uninit, ms, 3600.0_wp, par, &
                        .false., .false.)
      call check(error,.not. ice_uninit%is_init, &
                 "disabled_bitident: uninit slot must stay uninit after a no-op call")

      ! Case B: is_init but dynamics=.false. (transport-only ice slot).
      ice_off%enable = .true.
      ice_off%ncat = 1
      ice_off%nk_ice = 1
      ice_off%dynamics = .false.
      call ice_off%init(grid)
      ms%nz_ml = 1
      call ms%init(grid)
      !$acc enter data copyin(ice_off)
      call ice_off%enter_data()
      !$acc enter data copyin(ms)
      call ms%enter_data()

      block
         real(wp), allocatable :: u_ice0(:, :), v_ice0(:, :)
         allocate (u_ice0, source=ice_off%u_ice)
         allocate (v_ice0, source=ice_off%v_ice)
         call ice_evp_step(grid, metrics, f_corner, ice_off, ms, 3600.0_wp, par, &
                           .false., .false.)
         !$acc update self(ice_off%u_ice, ice_off%v_ice)
         if (.not. allocated(error)) then
            call check(error, all(ice_off%u_ice == u_ice0), &
                       "disabled_bitident: u_ice mutated by a dynamics=.false. call")
         end if
         if (.not. allocated(error)) then
            call check(error, all(ice_off%v_ice == v_ice0), &
                       "disabled_bitident: v_ice mutated by a dynamics=.false. call")
         end if
      end block

      !$acc exit data delete(f_corner)
      !$acc exit data delete(ms)
      call ms%exit_data()
      !$acc exit data delete(ice_off)
      call ice_off%exit_data()

      call ice_off%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_disabled_bitident

   ! =====================================================================
   ! PR 62: a_face_stress -- EVP stress weighting (momentum conservation)
   ! =====================================================================

   subroutine test_a_face_full_cover_bitident(error)
      !! Case 1: a fully-covered pack (`ci==1` everywhere, `test_nansen_
      !! free_drift`'s exact fully-periodic uniform-cover setup) must be
      !! the IDENTITY under `a_face_stress`: `a=1` reproduces the knob-off
      !! answer to a TIGHT relative tolerance (1e-10, matching this file's
      !! `test_nansen_free_drift` convention), on/off.  NOT exact `==`:
      !! empirically, on this GPU build, routing the SAME face's wind/drag
      !! through the extra `a_fac`/`tau_eff`/`drag_eff` locals (even at
      !! `a_fac==1.0_wp` exactly, itself IEEE-exact) versus the bare
      !! literals accumulates a ~1e-14 RELATIVE drift over the ~5000
      !! chained substeps in this test (432 subcycles/outer * N_OUTER) --
      !! the `7d283fe8` GPU FMA/codegen-non-associativity gotcha
      !! (CLAUDE.md Gotchas) applies even with the off arm fully separated
      !! and textually untouched (verified: the off arm alone reproduces
      !! itself bit-for-bit run-to-run; only crossing through the on arm's
      !! differently-shaped code introduces the drift).  The REAL default-
      !! off bit-identity gate is structural (the off arm's source is
      !! byte-for-byte the pre-PR expression) plus the full `ctest -R
      !! rdb` regression, which encodes tight pre-PR golden values this
      !! PR does not move.  Also catches an `a_u` built from the UNWRAPPED
      !! `ci` instead of `ci_w` -- the periodic seam would then read
      !! a_u=0.5, not 1, and this case would fail by an O(1) margin, far
      !! outside any FMA-drift tolerance.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 12
      real(wp), parameter :: TAU_A = 0.2_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: ui_off(:, :), vi_off(:, :), str_d_off(:, :), fxoc_off(:, :)
      integer :: nx, ny, n, i, j
      logical :: a_face

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         do n = 1, 2
            a_face = (n == 2)
            mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
            uo = 0.0_wp; vo = 0.0_wp
            tau_ax = 0.0_wp; tau_ay = 0.0_wp
            ui = 0.0_wp; vi = 0.0_wp
            str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
            f_corner = 0.0_wp
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP
                  mis(i, j) = 3.0_wp*ICE_RHO_ICE
                  mice(i, j) = 3.0_wp*ICE_RHO_ICE
                  ci(i, j) = 1.0_wp
               end do
            end do
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 1, NGHOST + NXP + 1
                  tau_ax(i, j) = TAU_A
               end do
            end do

            par = default_params(a_face_stress=a_face)
            call ws%init(nx, ny)
            call ws%enter_data()
            !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
            !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
            block
               integer :: n_inner
               do n_inner = 1, N_OUTER
                  call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                        tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                        fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
               end do
            end block
            !$acc update self(ui, vi, str_d, fxoc)
            !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
            !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
            call ws%exit_data()
            call ws%destroy()

            if (n == 1) then
               allocate (ui_off, source=ui)
               allocate (vi_off, source=vi)
               allocate (str_d_off, source=str_d)
               allocate (fxoc_off, source=fxoc)
            end if
         end do

         block
            real(wp), parameter :: TOL_REL = 1.0e-10_wp
            real(wp) :: scale_u, scale_str, scale_fxoc

            scale_u = max(maxval(abs(ui_off)), tiny(1.0_wp))
            scale_str = max(maxval(abs(str_d_off)), tiny(1.0_wp))
            scale_fxoc = max(maxval(abs(fxoc_off)), tiny(1.0_wp))

            call check(error, maxval(abs(ui - ui_off)) <= TOL_REL*scale_u, &
                       "a_face_full_cover_bitident: ui differs beyond FMA-drift tolerance")
            if (allocated(error)) exit checks
            call check(error, maxval(abs(vi - vi_off)) <= 1.0e-13_wp, &
                       "a_face_full_cover_bitident: vi differs")
            if (allocated(error)) exit checks
            call check(error, maxval(abs(str_d - str_d_off)) <= TOL_REL*scale_str, &
                       "a_face_full_cover_bitident: str_d differs beyond FMA-drift tolerance")
            if (allocated(error)) exit checks
            call check(error, maxval(abs(fxoc - fxoc_off)) <= TOL_REL*scale_fxoc, &
                       "a_face_full_cover_bitident: fxoc differs beyond FMA-drift tolerance")
         end block
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_full_cover_bitident

   subroutine test_a_face_ballistic_scaling(error)
      !! Case 2 (analytical): periodic uniform box, `f_corner=0` (cor=0),
      !! `uo=vo=0`, uniform `tau_ax`, uniform `mis=mice`, `ci==0.5`,
      !! `cdw=CDW_TINY` (drag_u -> 0).  One outer step from rest.  Pins
      !! the `a_u` FACTOR, PLACEMENT (numerator wind term), and that it
      !! reaches the kernel at all: a factor-of-2 discriminator between
      !! the on/off arms that a dead knob would fail.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      real(wp), parameter :: CI_CONST = 0.5_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j
      real(wp) :: ui_on, ui_off, expect_on, expect_off, tol

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         ! ---- a_face_stress = .true. ----
         mis = MI_CONST; mice = MI_CONST; ci = CI_CONST
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = TAU_A; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp
         par = default_params(cdw=CDW_TINY, a_face_stress=.true.)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         !$acc update self(ui)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()
         ui_on = ui(NGHOST + 2, NGHOST + 2)

         ! ---- a_face_stress = .false. ----
         mis = MI_CONST; mice = MI_CONST; ci = CI_CONST
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = TAU_A; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp
         par = default_params(cdw=CDW_TINY, a_face_stress=.false.)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         !$acc update self(ui)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()
         ui_off = ui(NGHOST + 2, NGHOST + 2)

         expect_on = DT_SLOW*CI_CONST*TAU_A/MI_CONST
         expect_off = DT_SLOW*TAU_A/MI_CONST
         tol = 1.0e-12_wp

         call check(error, abs(ui_on - expect_on) <= tol*expect_on, &
                    "a_face_ballistic_scaling: on-arm ui != dt*a*tau_a/mi")
         if (allocated(error)) exit checks
         call check(error, abs(ui_off - expect_off) <= tol*expect_off, &
                    "a_face_ballistic_scaling: off-arm ui != dt*tau_a/mi")
         if (allocated(error)) exit checks
         call check(error, abs(ui_off - 2.0_wp*ui_on) <= tol*ui_off, &
                    "a_face_ballistic_scaling: off/on is not the exact factor-of-2")
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_ballistic_scaling

   subroutine test_a_face_momentum_closes(error)
      !! Case 3 (analytical, the headline): periodic uniform box,
      !! `f_corner=0`, `ci==0.5`, DEFAULT cdw (drag live), from rest, ONE
      !! outer step (deliberately a TRANSIENT -- today's unweighted form
      !! is ALSO exactly conservative at steady free drift, so a converged
      !! run would not discriminate).  Drives the REAL coupler kernel
      !! (`ice_ocean_stress_flux_impl`) and checks the momentum the
      !! atmosphere supplied == the momentum the ice took + the momentum
      !! the ocean got, exactly, at fractional cover, in a transient --
      !! the whole PR in one assertion.  The knob-off control's residual
      !! is pinned QUANTITATIVELY to the predicted leak
      !! `dt_slow*(1-a)*(tau_a-fxoc)` (PR-62 plan Sec 3.2/3.3).
      use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux_impl
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      real(wp), parameter :: CI_CONST = 0.5_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: tau_x(:, :), tau_y(:, :)
      integer :: nx, ny, i, j
      real(wp) :: residual, predicted_leak, tol, worst_residual, worst_leak_mismatch

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      allocate (tau_x(nx + 1, ny), tau_y(nx, ny + 1))
      tol = 1.0e-12_wp*DT_SLOW*TAU_A

      checks: block
         ! ---- ON arm ----
         mis = MI_CONST; mice = MI_CONST; ci = CI_CONST
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = TAU_A; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp
         par = default_params(a_face_stress=.true.)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         !$acc update self(ui, fxoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         tau_x = 0.0_wp; tau_y = 0.0_wp
         call ice_ocean_stress_flux_impl(tau_x, tau_y, tau_ax, tau_ay, fxoc, fyoc, ci, nx, ny)

         worst_residual = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               residual = MI_CONST*ui(i, j) + DT_SLOW*tau_x(i, j) - DT_SLOW*tau_ax(i, j)
               worst_residual = max(worst_residual, abs(residual))
            end do
         end do
         call check(error, worst_residual <= tol, &
                    "a_face_momentum_closes: budget does not close with a_face_stress=.true.")
         if (allocated(error)) exit checks

         ! ---- OFF arm (control) ----
         mis = MI_CONST; mice = MI_CONST; ci = CI_CONST
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = TAU_A; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp
         par = default_params(a_face_stress=.false.)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         !$acc update self(ui, fxoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         tau_x = 0.0_wp; tau_y = 0.0_wp
         call ice_ocean_stress_flux_impl(tau_x, tau_y, tau_ax, tau_ay, fxoc, fyoc, ci, nx, ny)

         worst_leak_mismatch = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               residual = MI_CONST*ui(i, j) + DT_SLOW*tau_x(i, j) - DT_SLOW*tau_ax(i, j)
               predicted_leak = DT_SLOW*(1.0_wp - CI_CONST)*(tau_ax(i, j) - fxoc(i, j))
               worst_leak_mismatch = max(worst_leak_mismatch, abs(residual - predicted_leak))
            end do
         end do
         call check(error, worst_leak_mismatch <= tol, &
                    "a_face_momentum_closes: off-arm leak does not match (1-a)*(tau_a-fxoc)")
         if (allocated(error)) exit checks

         ! The off-arm residual must be a REAL leak (not accidentally ~0):
         ! guards against a degenerate transient that would not discriminate.
         call check(error, abs(residual) > 1.0e-6_wp*DT_SLOW*TAU_A, &
                    "a_face_momentum_closes: off-arm control is not a real leak "// &
                    "(transient too close to steady state -- strengthen it, do not delete)")
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_momentum_closes

   subroutine test_a_face_wind_only_would_leak(error)
      !! Case 4: guards the DESIGN, not the code.  Purely arithmetic (no
      !! new kernel call -- reruns case 3's off-arm setup to get a
      !! realistic transient `tau_a`/`fxoc`/`a`), evaluates the three
      !! Sec 3.2 rows of the PR-62 plan, and asserts weighting the wind
      !! ALONE leaks MORE at steady free drift than weighting neither --
      !! i.e. it is WORSE than doing nothing.  The next person who reads
      !! the task title ("a_face WIND weighting") and "simplifies" by
      !! deleting `drag_eff` should have this test fail on them.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      real(wp), parameter :: CI_CONST = 0.5_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j
      real(wp) :: leak_today, leak_wind_only, leak_both, tau_a_here, fxoc_here

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = MI_CONST; mice = MI_CONST; ci = CI_CONST
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = TAU_A; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp
         par = default_params(a_face_stress=.false.)

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                               tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                               fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         !$acc update self(fxoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         tau_a_here = TAU_A
         fxoc_here = fxoc(NGHOST + 2, NGHOST + 2)

         ! Sec 3.2 table, per unit cell area, per step:
         !   today (neither weighted):  leak = (1-a)*(tau_a - fxoc)
         !   wind only:                 leak = -(1-a)*a*tau_a
         !   both (this PR):            leak = 0
         leak_today = (1.0_wp - CI_CONST)*(tau_a_here - fxoc_here)
         leak_wind_only = -(1.0_wp - CI_CONST)*CI_CONST*tau_a_here
         leak_both = 0.0_wp

         call check(error, leak_both == 0.0_wp, &
                    "a_face_wind_only_would_leak: 'both' row is not exactly zero")
         if (allocated(error)) exit checks

         ! At STEADY free drift (fxoc == tau_a) today's row is EXACTLY zero
         ! while the wind-only row is a PERMANENT -(1-a)*a*tau_a: the
         ! wind-only fix converts a transient-only leak into one that never
         ! closes.  Verify the steady-state comparison directly (not just
         ! the recorded transient, whose today-leak also happens to be
         ! nonzero and would not by itself show wind-only is a regression
         ! at EVERY point in time, only at this one).
         block
            real(wp) :: leak_today_steady, leak_wind_only_steady
            leak_today_steady = (1.0_wp - CI_CONST)*(tau_a_here - tau_a_here)
            leak_wind_only_steady = -(1.0_wp - CI_CONST)*CI_CONST*tau_a_here
            call check(error, leak_today_steady == 0.0_wp, &
                       "a_face_wind_only_would_leak: today's steady-state leak != 0")
            if (allocated(error)) exit checks
            call check(error, abs(leak_wind_only_steady) > 0.0_wp, &
                       "a_face_wind_only_would_leak: wind-only steady-state leak == 0 "// &
                       "(it must be PERMANENT, not zero)")
            if (allocated(error)) exit checks
            call check(error, abs(leak_wind_only_steady) > abs(leak_today_steady), &
                       "a_face_wind_only_would_leak: wind-only is not worse than today "// &
                       "at steady free drift")
         end block
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_wind_only_would_leak

   subroutine test_a_face_no_ghost_drift(error)
      !! Case 5: `test_ice_margin`'s exact setup (20x4, west 10 iced
      !! ci=1, east 10 ice-free ci=0), knob ON.  Kills the ghost: every
      !! genuinely ice-free face (i in [NGHOST+12, NGHOST+NXP] -- BOTH
      !! neighbours ice-free, a_u==0 exactly; face NGHOST+11 is the mixed
      !! margin, a_u=0.5) has ui==uo==0 EXACTLY and fxoc==0 EXACTLY -- only
      !! reachable via the `a_fac > 0` BRANCH (a `max(a_u,eps)` floor or a
      !! naive `a_u*drag_u` gives 1e30/NaN here instead, PR-62 plan Sec
      !! 3.4).  The margin face's steady free-drift speed is UNCHANGED
      !! (a cancels out of the free-drift equilibrium `tau_a=drag(u)*u`).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 20, NYP = 4
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 20
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: u_free, edge_speed

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(a_face_stress=.true.)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + 10
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .false., .true., ws)
         end do

         !$acc update self(ui, str_d, fxoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error,.not. any(ieee_is_nan_grid(ui, nx + 1, ny)), &
                    "a_face_no_ghost_drift: NaN in ui")
         if (allocated(error)) exit checks
         call check(error,.not. any(ieee_is_nan_grid(fxoc, nx + 1, ny)), &
                    "a_face_no_ghost_drift: NaN in fxoc")
         if (allocated(error)) exit checks

         block
            logical :: ok_zero_u, ok_zero_fxoc
            ok_zero_u = .true.
            ok_zero_fxoc = .true.
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 12, NGHOST + NXP
                  if (ui(i, j) /= uo(i, j)) ok_zero_u = .false.
                  if (fxoc(i, j) /= 0.0_wp) ok_zero_fxoc = .false.
               end do
            end do
            call check(error, ok_zero_u, &
                       "a_face_no_ghost_drift: an ice-free face is not exactly at uo")
            if (allocated(error)) exit checks
            call check(error, ok_zero_fxoc, &
                       "a_face_no_ghost_drift: an ice-free face has nonzero fxoc")
            if (allocated(error)) exit checks
            call check(error, all(abs(ui(NGHOST + 12:NGHOST + NXP, &
                                         NGHOST + 1:NGHOST + NYP)) < huge(1.0_wp)), &
                       "a_face_no_ghost_drift: 1e30-scale blowup at an ice-free face")
            if (allocated(error)) exit checks
         end block
         if (allocated(error)) exit checks

         u_free = sqrt(TAU_A/(RHO_OCEAN_DEFAULT*CDW_DEFAULT))
         edge_speed = ui(NGHOST + 11, NGHOST + 2)
         call check(error, edge_speed > 0.5_wp*u_free, &
                    "a_face_no_ghost_drift: margin face steady drift changed with the knob on")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_no_ghost_drift

   subroutine test_a_face_disabled_bitident(error)
      !! Case 6: the house default-off gate, localised to the kernel (the
      !! full `ctest -R rdb` byte-identity is the real proof -- every
      !! existing EVP test already runs `default_params()`, whose type
      !! default is `a_face_stress=.false.`).  Explicitly builds
      !! `a_face_stress=.false.` and reruns `test_ice_margin`'s decay/
      !! ghost-drift golden and `test_nansen_free_drift`'s analytic |u|
      !! golden -- both must hold to the SAME tolerances as those tests.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 20, NYP = 4
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 20
      real(wp), parameter :: TAU_A = 0.1_wp
      real(wp), parameter :: MI_CONST = 3.0_wp*ICE_RHO_ICE
      integer, parameter :: I_SEED = NGHOST + 15

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: u_free, worst_far_east_err, edge_speed

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(a_face_stress=.false.)
      call check(error,.not. par%a_face_stress, &
                 "a_face_disabled_bitident: default_params() did not default a_face_stress off")
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + 10
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do
         str_d(I_SEED, NGHOST + 2) = -1000.0_wp

         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .false., .true., ws)
         end do

         !$acc update self(ui, str_d)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         ! Same golden as test_ice_margin: the massless-slab ghost drift.
         u_free = sqrt(TAU_A/(RHO_OCEAN_DEFAULT*CDW_DEFAULT))
         worst_far_east_err = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            worst_far_east_err = max(worst_far_east_err, &
                                     abs(ui(NGHOST + NXP, j) - u_free))
         end do
         call check(error, worst_far_east_err <= 1.0e-3_wp, &
                    "a_face_disabled_bitident: ghost drift golden moved with a_face_stress "// &
                    "explicitly off")
         if (allocated(error)) exit checks

         edge_speed = ui(NGHOST + 11, NGHOST + 2)
         call check(error, edge_speed > 0.5_wp*u_free, &
                    "a_face_disabled_bitident: iced-pack drift golden moved")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_a_face_disabled_bitident

   ! =====================================================================
   ! PR 36: CFL velocity truncation + PROJECT_ICE_CONCENTRATION
   ! =====================================================================

   subroutine test_trunc_disabled_bitident(error)
      !! Case 1 (PR-36): the house default-off gate, localised to the
      !! kernel (the real proof is the full byte-identical `ctest -R
      !! rdb`). `cfl_trunc=0` however spelled -- the type-level
      !! default (unset) vs an explicit `cfl_trunc=0.0` -- must give
      !! BIT-IDENTICAL `ui`/`vi`/`str_d`/`str_t`/`str_s`/`fxoc`/`fyoc`
      !! and `n_trunc==0`, driven with a tau strong enough that the
      !! analytic free-drift speed is well ABOVE what a `cfl_trunc=0.5`
      !! bound would allow -- proving the clip genuinely did not
      !! silently engage.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 20
      real(wp), parameter :: TAU_A = 5.0_wp
      real(wp), parameter :: WOULD_BE_BOUND = 0.95_wp*0.5_wp*DX/DT_SLOW
         !! `0.95*cfl_trunc*areaT/(dt_slow*dy_cu)` at `cfl_trunc=0.5` on
         !! this uniform grid (`areaT=DX*DX`, `dy_cu=DX`) -- documents
         !! that the drive is past the bound the clip WOULD apply if it
         !! were on.

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par1, par2
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: ui1(:, :), vi1(:, :), str_d1(:, :), str_t1(:, :), str_s1(:, :)
      real(wp), allocatable :: fxoc1(:, :), fyoc1(:, :)
      integer :: nx, ny, n, i, j, n_trunc1, n_trunc2
      real(wp) :: u_mean

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par1 = default_params()
      par2 = default_params(cfl_trunc=0.0_wp, cfl_trunc_dyn_its=.false., project_ci=.false.)
      call check(error, par1%cfl_trunc == 0.0_wp .and. (.not. par1%project_ci), &
                 "trunc_disabled_bitident: default_params() did not default "// &
                 "cfl_trunc/project_ci off")
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         f_corner = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = 3.0_wp*ICE_RHO_ICE
               mice(i, j) = 3.0_wp*ICE_RHO_ICE
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         ! ---- run 1: par1 = default_params() (type-level default off) ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par1, .true., .true., ws, &
                                  n_trunc=n_trunc1)
         end do
         !$acc update self(ui, vi, str_d, str_t, str_s, fxoc, fyoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         allocate (ui1, source=ui); allocate (vi1, source=vi)
         allocate (str_d1, source=str_d); allocate (str_t1, source=str_t)
         allocate (str_s1, source=str_s)
         allocate (fxoc1, source=fxoc); allocate (fyoc1, source=fyoc)

         ! ---- run 2: par2 = explicit cfl_trunc=0.0 ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         fxoc = 0.0_wp; fyoc = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par2, .true., .true., ws, &
                                  n_trunc=n_trunc2)
         end do
         !$acc update self(ui, vi, str_d, str_t, str_s, fxoc, fyoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, n_trunc1 == 0 .and. n_trunc2 == 0, &
                    "trunc_disabled_bitident: n_trunc must be 0 with cfl_trunc off "// &
                    "(however spelled)")
         if (allocated(error)) exit checks
         call check(error, all(ui == ui1), &
                    "trunc_disabled_bitident: ui differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(vi == vi1), &
                    "trunc_disabled_bitident: vi differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(str_d == str_d1), &
                    "trunc_disabled_bitident: str_d differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(str_t == str_t1), &
                    "trunc_disabled_bitident: str_t differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(str_s == str_s1), &
                    "trunc_disabled_bitident: str_s differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(fxoc == fxoc1), &
                    "trunc_disabled_bitident: fxoc differs by spelling of off")
         if (allocated(error)) exit checks
         call check(error, all(fyoc == fyoc1), &
                    "trunc_disabled_bitident: fyoc differs by spelling of off")
         if (allocated(error)) exit checks

         u_mean = sum(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)) &
                  /real((NXP + 1)*NYP, wp)
         call check(error, u_mean > WOULD_BE_BOUND, &
                    "trunc_disabled_bitident: drive must exceed the would-be "// &
                    "cfl_trunc=0.5 bound to prove the clip really did not fire")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      if (allocated(ui1)) deallocate (ui1, vi1, str_d1, str_t1, str_s1, fxoc1, fyoc1)
      call destroy_cartesian_metrics(metrics)
   end subroutine test_trunc_disabled_bitident

   subroutine test_trunc_bound_donor_asymmetry(error)
      !! Case 2 (PR-36, THE DISCRIMINATOR): direct call to the published
      !! `evp_truncate_final_impl` on a hand-built NON-UNIFORM `areaT`.
      !! Every other EVP test uses `make_cartesian_metrics` (uniform
      !! `areaT`), under which a symmetrised (WRONG) donor bound is
      !! indistinguishable from the correct asymmetric one -- this is
      !! the only gate that can tell them apart. Also pins the `0.95`
      !! back-off, the +/- sign, and the `mi > m_neglect` massless-face
      !! count gate (SIS2 :1466,1469).
      !!
      !! Minimal grid: nghost=1, nx_phys=ny_phys=1 -> nx=ny=3, ONE
      !! physical T-cell (2,2) flanked by ghost cells. u-face(3,2) (east
      !! face of the physical cell) has WEST donor areaT(2,2)=1e8, EAST
      !! donor areaT(3,2)=4e8; v-face(2,3) (north face) mirrors with
      !! areaT(2,2)=1e8 south donor, areaT(2,3)=4e8 north donor.
      !! u-face(2,2) (west face) has BOTH donors at the default 1e8 --
      !! the massless-face discriminator (E) reuses that symmetry.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3, NGH = 1, NXP = 1, NYP = 1
      real(wp), parameter :: CFL_TRUNC = 0.5_wp, DT_TR = 600.0_wp
      real(wp), parameter :: AREA_LO = 1.0e8_wp, AREA_HI = 4.0e8_wp, FACE_LEN = 1.0e4_wp
      real(wp), parameter :: M_NEGLECT_TEST = ICE_RHO_ICE*1.0e-30_wp
      real(wp), parameter :: U_HI = 0.95_wp*CFL_TRUNC*AREA_LO/(DT_TR*FACE_LEN)
      real(wp), parameter :: U_LO = -0.95_wp*CFL_TRUNC*AREA_HI/(DT_TR*FACE_LEN)
      real(wp), allocatable :: areaT(:, :), dy_cu(:, :), dx_cv(:, :)
      real(wp), allocatable :: mi_u(:, :), mi_v(:, :), ui(:, :), vi(:, :)
      integer :: n_trunc

      allocate (areaT(NX, NY), source=AREA_LO)
      areaT(3, 2) = AREA_HI
      areaT(2, 3) = AREA_HI
      allocate (dy_cu(NX + 1, NY), source=FACE_LEN)
      allocate (dx_cv(NX, NY + 1), source=FACE_LEN)
      allocate (mi_u(NX + 1, NY), source=1.0_wp)
      allocate (mi_v(NX, NY + 1), source=1.0_wp)
      allocate (ui(NX + 1, NY), source=0.0_wp)
      allocate (vi(NX, NY + 1), source=0.0_wp)

      !$acc enter data copyin(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi)

      checks: block
         ! (A) u positive clip: face (3,2), west donor 1e8.
         ui = 0.0_wp; vi = 0.0_wp
         ui(3, 2) = 100.0_wp
         !$acc update device(ui, vi)
         call evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                      CFL_TRUNC, DT_TR, M_NEGLECT_TEST, NGH, NXP, NYP, &
                                      NX, NY, n_trunc)
         !$acc update self(ui, vi)
         call check(error, abs(ui(3, 2) - U_HI) <= 1.0e-13_wp*abs(U_HI), &
                    "trunc_bound_donor_asymmetry: +u clip against WEST donor mismatch")
         if (allocated(error)) exit checks
         call check(error, n_trunc == 1, &
                    "trunc_bound_donor_asymmetry: +u case must count exactly 1 face")
         if (allocated(error)) exit checks

         ! (B) u negative clip: face (3,2), EAST donor 4e8.
         ui = 0.0_wp; vi = 0.0_wp
         ui(3, 2) = -100.0_wp
         !$acc update device(ui, vi)
         call evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                      CFL_TRUNC, DT_TR, M_NEGLECT_TEST, NGH, NXP, NYP, &
                                      NX, NY, n_trunc)
         !$acc update self(ui, vi)
         call check(error, abs(ui(3, 2) - U_LO) <= 1.0e-13_wp*abs(U_LO), &
                    "trunc_bound_donor_asymmetry: -u clip against EAST donor mismatch")
         if (allocated(error)) exit checks
         call check(error, n_trunc == 1, &
                    "trunc_bound_donor_asymmetry: -u case must count exactly 1 face")
         if (allocated(error)) exit checks

         ! (C) v positive clip: face (2,3), SOUTH donor areaT(2,2)=1e8.
         ui = 0.0_wp; vi = 0.0_wp
         vi(2, 3) = 100.0_wp
         !$acc update device(ui, vi)
         call evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                      CFL_TRUNC, DT_TR, M_NEGLECT_TEST, NGH, NXP, NYP, &
                                      NX, NY, n_trunc)
         !$acc update self(ui, vi)
         call check(error, abs(vi(2, 3) - U_HI) <= 1.0e-13_wp*abs(U_HI), &
                    "trunc_bound_donor_asymmetry: +v clip against SOUTH donor mismatch")
         if (allocated(error)) exit checks
         call check(error, n_trunc == 1, &
                    "trunc_bound_donor_asymmetry: +v case must count exactly 1 face")
         if (allocated(error)) exit checks

         ! (D) v negative clip: face (2,3), NORTH donor areaT(2,3)=4e8.
         ui = 0.0_wp; vi = 0.0_wp
         vi(2, 3) = -100.0_wp
         !$acc update device(ui, vi)
         call evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                      CFL_TRUNC, DT_TR, M_NEGLECT_TEST, NGH, NXP, NYP, &
                                      NX, NY, n_trunc)
         !$acc update self(ui, vi)
         call check(error, abs(vi(2, 3) - U_LO) <= 1.0e-13_wp*abs(U_LO), &
                    "trunc_bound_donor_asymmetry: -v clip against NORTH donor mismatch")
         if (allocated(error)) exit checks
         call check(error, n_trunc == 1, &
                    "trunc_bound_donor_asymmetry: -v case must count exactly 1 face")
         if (allocated(error)) exit checks

         ! (E) massless-face discriminator: face (2,2), symmetric default
         ! donors (both 1e8) -- clips to the SAME magnitude as (A), but
         ! mi_u(2,2)=0 means it must NOT be counted (SIS2 :1466,1469).
         ui = 0.0_wp; vi = 0.0_wp
         ui(2, 2) = 100.0_wp
         mi_u(2, 2) = 0.0_wp
         !$acc update device(ui, vi, mi_u)
         call evp_truncate_final_impl(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi, &
                                      CFL_TRUNC, DT_TR, M_NEGLECT_TEST, NGH, NXP, NYP, &
                                      NX, NY, n_trunc)
         !$acc update self(ui, vi)
         call check(error, abs(ui(2, 2) - U_HI) <= 1.0e-13_wp*abs(U_HI), &
                    "trunc_bound_donor_asymmetry: massless face must still be clipped")
         if (allocated(error)) exit checks
         call check(error, n_trunc == 0, &
                    "trunc_bound_donor_asymmetry: massless face must NOT be counted")
      end block checks

      !$acc exit data delete(areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi)
      deallocate (areaT, dy_cu, dx_cv, mi_u, mi_v, ui, vi)
   end subroutine test_trunc_bound_donor_asymmetry

   subroutine test_trunc_final_clips_and_counts(error)
      !! Case 3 (PR-36): through `ice_evp_dynamics`, `cfl_trunc=0.5`, a
      !! tau large enough to drive free drift past the bound -- on
      !! return every physical face satisfies the CLIPPED range and
      !! `n_trunc` is > 0. The integration gate: the clip is actually
      !! reached from the production seam and the returned velocity is
      !! transport-safe. (The `mi_u=0`-not-counted discriminator is
      !! pinned directly on the kernel by case 2(E) above -- reaching
      !! that exact state through the full `ice_evp_dynamics` seam would
      !! require a face already zeroed by `evp_zero_massless_velocity_impl`,
      !! which is never subsequently re-driven past the bound.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 5
      real(wp), parameter :: TAU_A = 5.0_wp
      real(wp), parameter :: BOUND = 0.95_wp*0.5_wp*DX/DT_SLOW

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, n, i, j, n_trunc

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(cfl_trunc=0.5_wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)
      call ws%init(nx, ny)
      call ws%enter_data()

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = 3.0_wp*ICE_RHO_ICE
               mice(i, j) = 3.0_wp*ICE_RHO_ICE
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .true., ws, &
                                  n_trunc=n_trunc)
         end do

         !$acc update self(ui, vi)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, n_trunc > 0, &
                    "trunc_final_clips_and_counts: a forced-runaway drive must clip >0 faces")
         if (allocated(error)) exit checks

         call check(error, all(abs(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)) &
                               <= BOUND*(1.0_wp + 1.0e-9_wp)), &
                    "trunc_final_clips_and_counts: clipped ui exceeds the 0.95*bound ceiling")
         if (allocated(error)) exit checks
         call check(error, all(abs(vi(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP + 1)) &
                               <= BOUND*(1.0_wp + 1.0e-9_wp)), &
                    "trunc_final_clips_and_counts: clipped vi exceeds the 0.95*bound ceiling")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_trunc_final_clips_and_counts

   subroutine test_trunc_below_bound_bitident(error)
      !! Case 4 (PR-36): the no-false-trigger gate. `cfl_trunc=0.5` ON,
      !! but a drift far below the bound: `ui`/`vi` (physical faces) and
      !! `str_d`/`str_t`/`str_s`/`fxoc`/`fyoc` (full array -- neither is
      !! touched by the clip machinery) are bit-identical to
      !! `cfl_trunc=0`, and `n_trunc==0`. A sign error, a missing
      !! `dt_tr` in the denominator, or a lost `areaT` would clamp every
      !! ordinary drift to a wrong value and this is what catches it.
      !! (Physical-face-only comparison for `ui`/`vi`: the FINAL clip's
      !! mandatory periodic re-wrap runs whenever `cfl_trunc>0`,
      !! refreshing the ghost rows even when no face actually clips --
      !! an intentional PR-36 divergence in GHOST bookkeeping only, not
      !! a claim about the physical interior.)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 48
      real(wp), parameter :: TAU_A = 0.05_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par_off, par_on
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: ui_off(:, :), vi_off(:, :), str_d_off(:, :)
      real(wp), allocatable :: str_t_off(:, :), str_s_off(:, :), fxoc_off(:, :), fyoc_off(:, :)
      integer :: nx, ny, n, i, j, n_trunc

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par_off = default_params()
      par_on = default_params(cfl_trunc=0.5_wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         f_corner = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = 3.0_wp*ICE_RHO_ICE
               mice(i, j) = 3.0_wp*ICE_RHO_ICE
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         ! ---- run OFF ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par_off, .true., .true., ws)
         end do
         !$acc update self(ui, vi, str_d, str_t, str_s, fxoc, fyoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         allocate (ui_off, source=ui); allocate (vi_off, source=vi)
         allocate (str_d_off, source=str_d); allocate (str_t_off, source=str_t)
         allocate (str_s_off, source=str_s)
         allocate (fxoc_off, source=fxoc); allocate (fyoc_off, source=fyoc)

         ! ---- run ON: cfl_trunc=0.5, but drift stays below bound ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         fxoc = 0.0_wp; fyoc = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par_on, .true., .true., ws, &
                                  n_trunc=n_trunc)
         end do
         !$acc update self(ui, vi, str_d, str_t, str_s, fxoc, fyoc)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, n_trunc == 0, &
                    "trunc_below_bound_bitident: an ordinary drift must not clip (n_trunc=0)")
         if (allocated(error)) exit checks
         call check(error, all(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP) == &
                               ui_off(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)), &
                    "trunc_below_bound_bitident: ui (physical) moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(vi(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP + 1) == &
                               vi_off(NGHOST + 1:NGHOST + NXP, NGHOST + 1:NGHOST + NYP + 1)), &
                    "trunc_below_bound_bitident: vi (physical) moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(str_d == str_d_off), &
                    "trunc_below_bound_bitident: str_d moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(str_t == str_t_off), &
                    "trunc_below_bound_bitident: str_t moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(str_s == str_s_off), &
                    "trunc_below_bound_bitident: str_s moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(fxoc == fxoc_off), &
                    "trunc_below_bound_bitident: fxoc moved with cfl_trunc on")
         if (allocated(error)) exit checks
         call check(error, all(fyoc == fyoc_off), &
                    "trunc_below_bound_bitident: fyoc moved with cfl_trunc on")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      if (allocated(ui_off)) deallocate (ui_off, vi_off, str_d_off, str_t_off, str_s_off, &
                                         fxoc_off, fyoc_off)
      call destroy_cartesian_metrics(metrics)
   end subroutine test_trunc_below_bound_bitident

   subroutine test_trunc_uses_transport_dt(error)
      !! Case 5 (PR-36): the sharpest structural finding (Sec 2.5),
      !! asserted. Identical state, `cfl_trunc=0.5`, driven hard enough
      !! to saturate the clip every outer call: `dt_transport` ABSENT
      !! must equal `dt_transport=dt_slow` bit-for-bit, and
      !! `dt_transport=4*dt_slow` must clip to EXACTLY 1/4 the bound
      !! (`evp_truncate_final_impl` overwrites a clipped face
      !! unconditionally with a value that depends only on
      !! `areaT`/`dy_cu`/`cfl_trunc`/`dt_tr` -- never on the momentum
      !! solve's trajectory -- so once saturated the ratio is exact).
      !! Without this, `dt_transport` is a dead argument the driver
      !! could silently stop passing, and at `dt_therm_ratio > 1` the
      !! bound would be too loose by exactly that factor.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 4, NYP = 4
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DT_SLOW = 3600.0_wp
      integer, parameter :: N_OUTER = 5
      real(wp), parameter :: TAU_A = 5.0_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: ui_a(:, :), vi_a(:, :)
      integer :: nx, ny, n, i, j

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par = default_params(cfl_trunc=0.5_wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         f_corner = 0.0_wp
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP
               mis(i, j) = 3.0_wp*ICE_RHO_ICE
               mice(i, j) = 3.0_wp*ICE_RHO_ICE
               ci(i, j) = 1.0_wp
            end do
         end do
         do j = NGHOST + 1, NGHOST + NYP
            do i = NGHOST + 1, NGHOST + NXP + 1
               tau_ax(i, j) = TAU_A
            end do
         end do

         ! ---- run A: dt_transport ABSENT (=> dt_slow) ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .true., ws)
         end do
         !$acc update self(ui, vi)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()
         allocate (ui_a, source=ui); allocate (vi_a, source=vi)

         ! ---- run B: dt_transport = 4*dt_slow ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         fxoc = 0.0_wp; fyoc = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .true., ws, &
                                  dt_transport=4.0_wp*DT_SLOW)
         end do
         !$acc update self(ui, vi)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()

         call check(error, &
                    all(abs(ui(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP) - &
                            0.25_wp*ui_a(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP)) &
                        <= 1.0e-13_wp*abs(ui_a(NGHOST + 1:NGHOST + NXP + 1, &
                                               NGHOST + 1:NGHOST + NYP)) + 1.0e-13_wp), &
                    "trunc_uses_transport_dt: dt_transport=4*dt_slow must clip to 1/4 the bound")
         if (allocated(error)) exit checks

         ! ---- run C: dt_transport = dt_slow (explicit) -> matches A bit-for-bit ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         fxoc = 0.0_wp; fyoc = 0.0_wp
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_OUTER
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, DT_SLOW, par, .true., .true., ws, &
                                  dt_transport=DT_SLOW)
         end do
         !$acc update self(ui, vi)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, all(ui == ui_a), &
                    "trunc_uses_transport_dt: explicit dt_transport=dt_slow must match "// &
                    "ABSENT bit-for-bit")
         if (allocated(error)) exit checks
         call check(error, all(vi == vi_a), &
                    "trunc_uses_transport_dt: explicit dt_transport=dt_slow (v) must "// &
                    "match ABSENT")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      if (allocated(ui_a)) deallocate (ui_a, vi_a)
      call destroy_cartesian_metrics(metrics)
   end subroutine test_trunc_uses_transport_dt

   subroutine test_project_ci_zero_divergence_bitident(error)
      !! Case 7 (PR-36): `PROJECT_ICE_CONCENTRATION` under PURE SHEAR
      !! (`sh_Dd == 0` by construction, `impose_shear`/`test_shear_2d`'s
      !! setup) is an EXACT no-op: `exp(-t_cum*0) = 1 => ci_proj == ci`,
      !! so `pres_mice` is untouched regardless of `project_ci`. Run at
      !! `ci=0.9` (unsaturated -- `max(1-ci,0) /= 0`) so the assertion
      !! isn't a vacuous "already-saturated" case. A dropped minus, a
      !! `sh_Dt` used for `sh_Dd`, or a `t_cum` that is `dt` instead of
      !! `n*dt` all survive a naive "does it stiffen" test but die here.
      !!
      !! `impose_shear` zeros `v` at the physical N/S WALL faces
      !! (`j=nghost+1`/`j=nghost+nxy+1`), so `sh_Dd` is exactly 0 only in
      !! the DEEP interior (`j` in `[nghost+2, nghost+nxy-1]`) -- the two
      !! rows immediately touching a wall see `v_wall=0` against
      !! `v_interior=gamma_dot*x_i`, a genuine nonzero divergence there
      !! (confirmed by `test_shear_2d`'s own near-wall requirement-4
      !! checks). The bit-identity assertion below is restricted to that
      !! interior J range accordingly -- it is not a weaker test, it is
      !! the CORRECT domain of the "sh_Dd==0" claim.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXY = 20
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: CI_VAL = 0.9_wp
      real(wp), parameter :: MI_CONST = 2.0_wp*ICE_RHO_ICE
      real(wp), parameter :: GAMMA_DOT = 1.0e-5_wp
      real(wp), parameter :: DT_SLOW_NOMINAL = 1800.0_wp
      integer, parameter :: N_SUB_NOMINAL = 432
      integer, parameter :: N_REIMPOSE = 200
      integer, parameter :: J_LO = NGHOST + 2, J_HI = NGHOST + NXY - 1
         !! Deep-interior row range where sh_Dd==0 exactly (excludes the
         !! two wall-adjacent rows).

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par_off, par_on
      type(evp_workspace_t) :: ws
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      real(wp), allocatable :: str_d_off(:, :), str_t_off(:, :), str_s_off(:, :)
      integer :: nx, ny, n, i, j
      real(wp) :: dt_eff

      call grid%init(NXY, NXY, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      par_off = default_params(evp_sub_steps=1, project_ci=.false.)
      par_on = default_params(evp_sub_steps=1, project_ci=.true.)
      dt_eff = DT_SLOW_NOMINAL/real(N_SUB_NOMINAL, wp)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp; ci = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         ui = 0.0_wp; vi = 0.0_wp
         str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         f_corner = 0.0_wp

         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
               ci(i, j) = CI_VAL
            end do
         end do

         call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)

         ! ---- run OFF ----
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_REIMPOSE
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, dt_eff, par_off, .false., .false., ws)
            ui = 0.0_wp
            call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)
            !$acc update device(ui, vi)
         end do
         !$acc update self(str_d, str_t, str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         call ws%exit_data()
         call ws%destroy()
         allocate (str_d_off, source=str_d); allocate (str_t_off, source=str_t)
         allocate (str_s_off, source=str_s)

         ! ---- run ON ----
         ui = 0.0_wp; vi = 0.0_wp; str_d = 0.0_wp; str_t = 0.0_wp; str_s = 0.0_wp
         call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)
         call ws%init(nx, ny)
         call ws%enter_data()
         !$acc enter data copyin(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc enter data copyin(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)
         do n = 1, N_REIMPOSE
            call ice_evp_dynamics(grid, metrics, f_corner, mis, mice, ci, uo, vo, &
                                  tau_ax, tau_ay, ui, vi, str_d, str_t, str_s, &
                                  fxoc, fyoc, dt_eff, par_on, .false., .false., ws)
            ui = 0.0_wp
            call impose_shear(vi, GAMMA_DOT, DX, NXY, NGHOST, nx, ny)
            !$acc update device(ui, vi)
         end do
         !$acc update self(str_d, str_t, str_s)
         !$acc exit data delete(mis, mice, ci, uo, vo, tau_ax, tau_ay)
         !$acc exit data delete(ui, vi, str_d, str_t, str_s, fxoc, fyoc, f_corner)

         call check(error, all(str_d(NGHOST + 1:NGHOST + NXY, J_LO:J_HI) == &
                               str_d_off(NGHOST + 1:NGHOST + NXY, J_LO:J_HI)), &
                    "project_ci_zero_divergence_bitident: str_d moved under pure shear")
         if (allocated(error)) exit checks
         call check(error, all(str_t(NGHOST + 1:NGHOST + NXY, J_LO:J_HI) == &
                               str_t_off(NGHOST + 1:NGHOST + NXY, J_LO:J_HI)), &
                    "project_ci_zero_divergence_bitident: str_t moved under pure shear")
         if (allocated(error)) exit checks
         call check(error, all(str_s(NGHOST + 1:NGHOST + NXY + 1, J_LO + 1:J_HI) == &
                               str_s_off(NGHOST + 1:NGHOST + NXY + 1, J_LO + 1:J_HI)), &
                    "project_ci_zero_divergence_bitident: str_s moved under pure shear")
      end block checks

      call ws%exit_data()
      call ws%destroy()
      if (allocated(str_d_off)) deallocate (str_d_off, str_t_off, str_s_off)
      call destroy_cartesian_metrics(metrics)
   end subroutine test_project_ci_zero_divergence_bitident

   subroutine test_project_ci_stiffens_convergence(error)
      !! Case 8 (PR-36, analytical, signed) -- the whole feature in one
      !! case, including its limit. `evp_sub_steps=1` makes `t_cum =
      !! dt_slow` exactly and the algebra closed-form: with `str_d`
      !! starting at 0 (a fresh call) and `zeta = 0.5*pres_mice*mice/
      !! del_sh` (`del_sh` purely kinematic -- independent of
      !! `pres_mice`/`ci`), `str_d` after the single substep is LINEAR
      !! in `pres_mice`, so `str_d(project_ci=on)/str_d(off) ==
      !! pres_mice_proj/pres_mice0` EXACTLY (in this regime, away from
      !! the `del_sh_min_pr` floor branch). (a)+(b) pin BOTH the
      !! magnitude (to the analytic ratio, 1e-12 rel) AND both signs --
      !! a symmetric "it changed" assertion would pass with the sign
      !! flipped, which would make converging ice WEAKER and actively
      !! destabilise the pack. (c) is Sec 2.4's conclusion written into
      !! the suite: at `ci=1` (a saturated jam), `max(1-ci_proj,0)=0`
      !! either way, so `project_ci` is PROVABLY inert at the wall.
      !!
      !! `GAMMA`/`DT_SLOW` are chosen so `del_sh` (purely kinematic,
      !! `~GAMMA*sqrt(1+1/EC^2)`) dominates `del_sh_min_pr*pres_mice`
      !! (the floor branch, where zeta -- and hence the ratio -- is NOT
      !! linear in `pres_mice`) by ~2 orders of magnitude: at the
      !! SIS2-default `del_sh_min_scale=2`/`tdamp<0` rule, `tdamp_eff =
      !! 3*dt_slow` (since `evp_sub_steps=1` => `dt=dt_slow`), so
      !! `del_sh_min_pr = 2*del_sh_min_scale*dt_slow/(3*dxharm^2)` grows
      !! LINEARLY with `dt_slow` -- a smaller `dt_slow` (100s, not the
      !! 1800s SIS2-nominal value) keeps the floor comfortably
      !! subdominant while `GAMMA` stays small enough for a genuinely
      !! UNSATURATED `ci_proj` (a graded ratio, not a degenerate 0/1 clip).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXY = 20
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: MI_CONST = 2.0_wp*ICE_RHO_ICE
      real(wp), parameter :: GAMMA = 2.0e-4_wp
      real(wp), parameter :: DT_SLOW = 100.0_wp
      real(wp), parameter :: CI_VAL = 0.9_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par_off, par_on
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j, i_mid, j_mid
      real(wp) :: str_d_off_conv, str_d_on_conv, str_d_off_div, str_d_on_div
      real(wp) :: str_d_off_full, str_d_on_full
      real(wp) :: ci_proj_conv, ci_proj_div, ratio_conv, ratio_div

      call grid%init(NXY, NXY, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      i_mid = NGHOST + NXY/2
      j_mid = NGHOST + NXY/2

      par_off = default_params(evp_sub_steps=1, project_ci=.false.)
      par_on = default_params(evp_sub_steps=1, project_ci=.true.)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         vi = 0.0_wp
         f_corner = 0.0_wp
         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
            end do
         end do

         ! ---- convergent field: sh_Dd = -GAMMA ----
         ci = CI_VAL
         call impose_divergence(ui, -GAMMA, DX, NXY, NGHOST, nx, ny)
         str_d_off_conv = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                                uo, vo, tau_ax, tau_ay, ui, vi, &
                                                DT_SLOW, par_off, i_mid, j_mid)
         str_d_on_conv = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                               uo, vo, tau_ax, tau_ay, ui, vi, &
                                               DT_SLOW, par_on, i_mid, j_mid)

         ci_proj_conv = CI_VAL*exp(GAMMA*DT_SLOW)
            !! -t_cum*sh_Dd = -DT_SLOW*(-GAMMA) = +GAMMA*DT_SLOW.
         ratio_conv = exp(-par_off%c0*(max(1.0_wp - ci_proj_conv, 0.0_wp) - &
                                       max(1.0_wp - CI_VAL, 0.0_wp)))

         call check(error, ratio_conv > 1.0_wp, &
                    "project_ci_stiffens_convergence: sanity -- convergence ratio must be > 1")
         if (allocated(error)) exit checks
         call check(error, abs(str_d_on_conv) > abs(str_d_off_conv), &
                    "project_ci_stiffens_convergence: convergence must STIFFEN "// &
                    "(|str_d| larger)")
         if (allocated(error)) exit checks
         call check(error, abs(str_d_on_conv - ratio_conv*str_d_off_conv) &
                    <= 1.0e-12_wp*abs(ratio_conv*str_d_off_conv), &
                    "project_ci_stiffens_convergence: convergence ratio mismatch")
         if (allocated(error)) exit checks

         ! ---- divergent field: sh_Dd = +GAMMA ----
         call impose_divergence(ui, GAMMA, DX, NXY, NGHOST, nx, ny)
         str_d_off_div = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                               uo, vo, tau_ax, tau_ay, ui, vi, &
                                               DT_SLOW, par_off, i_mid, j_mid)
         str_d_on_div = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                              uo, vo, tau_ax, tau_ay, ui, vi, &
                                              DT_SLOW, par_on, i_mid, j_mid)

         ci_proj_div = CI_VAL*exp(-GAMMA*DT_SLOW)
         ratio_div = exp(-par_off%c0*(max(1.0_wp - ci_proj_div, 0.0_wp) - &
                                      max(1.0_wp - CI_VAL, 0.0_wp)))

         call check(error, ratio_div < 1.0_wp, &
                    "project_ci_stiffens_convergence: sanity -- divergence ratio must be < 1")
         if (allocated(error)) exit checks
         call check(error, abs(str_d_on_div) < abs(str_d_off_div), &
                    "project_ci_stiffens_convergence: divergence must WEAKEN "// &
                    "(|str_d| smaller)")
         if (allocated(error)) exit checks
         call check(error, abs(str_d_on_div - ratio_div*str_d_off_div) &
                    <= 1.0e-12_wp*abs(ratio_div*str_d_off_div), &
                    "project_ci_stiffens_convergence: divergence ratio mismatch")
         if (allocated(error)) exit checks

         ! ---- (c) ci=1 (saturated jam) + convergent field: bit-identical ----
         ci = 1.0_wp
         call impose_divergence(ui, -GAMMA, DX, NXY, NGHOST, nx, ny)
         str_d_off_full = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                                uo, vo, tau_ax, tau_ay, ui, vi, &
                                                DT_SLOW, par_off, i_mid, j_mid)
         str_d_on_full = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                               uo, vo, tau_ax, tau_ay, ui, vi, &
                                               DT_SLOW, par_on, i_mid, j_mid)
         call check(error, str_d_on_full == str_d_off_full, &
                    "project_ci_stiffens_convergence: ci=1 (saturated) must be "// &
                    "bit-identical")
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_project_ci_stiffens_convergence

   subroutine test_project_ci_extreme_divergence_finite(error)
      !! Case 9 (PR-36): `t_cum*|sh_Dd| > 709` overflows `exp` to
      !! `+Inf` => `ci_proj = Inf` => `max(1-Inf,0) = max(-Inf,0) = 0`
      !! => `exp(-c0*0) = 1` => `pres_mice = p0_rho`, the CORRECT
      !! saturated value -- IEEE launders the overflow benignly, no
      !! guard needed (SIS2 has none). Assert `str_d` stays FINITE (not
      !! NaN) and equals the `ci=1` (fully saturated, `project_ci=off`)
      !! reference to 1e-13 rel -- both use `pres_mice=p0_rho` in their
      !! single substep, so an exact match is the correct outcome, not
      !! a coincidence.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXY = 20
      real(wp), parameter :: DX = 25000.0_wp
      real(wp), parameter :: MI_CONST = 2.0_wp*ICE_RHO_ICE
      real(wp), parameter :: GAMMA_EXTREME = 1.0_wp
         !! `DT_SLOW*GAMMA_EXTREME = 1800 >> 709` -- guarantees `exp` overflow.
      real(wp), parameter :: DT_SLOW = 1800.0_wp
      real(wp), parameter :: CI_VAL = 0.9_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ice_evp_params_t) :: par_on, par_ref
      real(wp), allocatable :: mis(:, :), mice(:, :), ci(:, :), uo(:, :), vo(:, :)
      real(wp), allocatable :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable :: ui(:, :), vi(:, :), str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable :: fxoc(:, :), fyoc(:, :), f_corner(:, :)
      integer :: nx, ny, i, j, i_mid, j_mid
      real(wp) :: str_d_extreme, str_d_ref

      call grid%init(NXY, NXY, NGHOST, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      i_mid = NGHOST + NXY/2
      j_mid = NGHOST + NXY/2

      par_on = default_params(evp_sub_steps=1, project_ci=.true.)
      par_ref = default_params(evp_sub_steps=1, project_ci=.false.)

      call alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                            str_d, str_t, str_s, fxoc, fyoc, f_corner)

      checks: block
         mis = 0.0_wp; mice = 0.0_wp
         uo = 0.0_wp; vo = 0.0_wp
         tau_ax = 0.0_wp; tau_ay = 0.0_wp
         f_corner = 0.0_wp
         do j = NGHOST + 1, NGHOST + NXY
            do i = NGHOST + 1, NGHOST + NXY
               mis(i, j) = MI_CONST
               mice(i, j) = MI_CONST
            end do
         end do
         call impose_divergence(ui, -GAMMA_EXTREME, DX, NXY, NGHOST, nx, ny)
         vi = 0.0_wp

         ci = CI_VAL
         str_d_extreme = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                               uo, vo, tau_ax, tau_ay, ui, vi, &
                                               DT_SLOW, par_on, i_mid, j_mid)

         call check(error, str_d_extreme == str_d_extreme, &
                    "project_ci_extreme_divergence_finite: str_d must not be NaN")
            !! `(x /= x)` is the portable NaN test; `check(x==x)` fails on NaN.
         if (allocated(error)) exit checks
         call check(error, abs(str_d_extreme) < huge(1.0_wp), &
                    "project_ci_extreme_divergence_finite: str_d must not be +/-Inf")
         if (allocated(error)) exit checks

         ci = 1.0_wp
         str_d_ref = run_one_substep_str_d(grid, metrics, f_corner, mis, mice, ci, &
                                           uo, vo, tau_ax, tau_ay, ui, vi, &
                                           DT_SLOW, par_ref, i_mid, j_mid)

         call check(error, abs(str_d_extreme - str_d_ref) <= 1.0e-13_wp*abs(str_d_ref), &
                    "project_ci_extreme_divergence_finite: must equal the "// &
                    "p0_rho-saturated reference (ci=1, project_ci=off)")
      end block checks

      call destroy_cartesian_metrics(metrics)
   end subroutine test_project_ci_extreme_divergence_finite

   ! =====================================================================
   ! Shared allocation helper
   ! =====================================================================

   subroutine alloc_evp_arrays(nx, ny, mis, mice, ci, uo, vo, tau_ax, tau_ay, ui, vi, &
                               str_d, str_t, str_s, fxoc, fyoc, f_corner)
      integer, intent(in) :: nx, ny
      real(wp), allocatable, intent(out) :: mis(:, :), mice(:, :), ci(:, :)
      real(wp), allocatable, intent(out) :: uo(:, :), vo(:, :)
      real(wp), allocatable, intent(out) :: tau_ax(:, :), tau_ay(:, :)
      real(wp), allocatable, intent(out) :: ui(:, :), vi(:, :)
      real(wp), allocatable, intent(out) :: str_d(:, :), str_t(:, :), str_s(:, :)
      real(wp), allocatable, intent(out) :: fxoc(:, :), fyoc(:, :)
      real(wp), allocatable, intent(out) :: f_corner(:, :)

      allocate (mis(nx, ny), mice(nx, ny), ci(nx, ny))
      allocate (uo(nx + 1, ny), vo(nx, ny + 1))
      allocate (tau_ax(nx + 1, ny), tau_ay(nx, ny + 1))
      allocate (ui(nx + 1, ny), vi(nx, ny + 1))
      allocate (str_d(nx, ny), str_t(nx, ny), str_s(nx + 1, ny + 1))
      allocate (fxoc(nx + 1, ny), fyoc(nx, ny + 1))
      allocate (f_corner(nx + 1, ny + 1))
   end subroutine alloc_evp_arrays

end module test_ocean_ice_evp
