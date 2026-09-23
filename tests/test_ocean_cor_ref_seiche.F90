!! ROTATING CLOSED-BASIN NON-GROWTH — the barotropic Coriolis-reference
!! gate for the `pred_corr` outer split scheme.
!!
!! ## What is guarded
!!
!! The split-explicit driver hands the barotropic substep a frozen
!! forcing `F_bt_u/v` that already contains the depth mean of the SLOW
!! layer Coriolis-advection tendency, and the substep then integrates
!! its own LIVE `(ζ+f)·v − ∇KE` on top.  `subtract_fast_cor_ref`
!! removes the double count by subtracting a REFERENCE `(ζ+f)·v − ∇KE`
!! evaluated at a barotropic velocity `bt_work%cor_ref_u/v` (MOM6
!! `Cor_ref_u`, built there from `ubt_Cor`).
!!
!! That reference velocity is not free.  It must be the depth mean of
!! the SAME layer velocity the slow tendency was evaluated on:
!!
!!   * `ssp_rk2`  — slow tendencies on the prognostic `u^n`  ⇒ the
!!     stage-entry `bt_ubt/bt_vbt`;
!!   * `pred_corr`  — slow tendencies on the time-mean `u_av/v_av` ⇒ the
!!     depth mean of `u_av/v_av`.
!!
!! Use `u^n` under `pred_corr` and the two no longer cancel: the leftover
!! `f × (v̄_av − v̄^n)` is injected into EVERY barotropic substep as a
!! near-constant forcing.  In a closed rotating basin it projects onto
!! the gravest Poincaré seiche and pumps it exponentially — five orders
!! of magnitude and then a NaN on a realistic configuration, with a
!! growth rate ∝ dt and RISING with substep count (the fingerprint of a
!! fixed per-substep forcing rather than an inner-loop instability).
!!
!! ## Why this shape of test
!!
!! The property is energetic, not pointwise: the basin is closed,
!! unforced, undamped, inviscid and uniform-density, so `KE + PE` is a
!! conserved quantity of the continuous system.  Any sustained GAIN is
!! the integrator manufacturing energy.  A 51-case A/B over shipped
!! namelists missed this defect entirely because no shipped case is a
!! long rotating closed basin — hence a dedicated one here.
!!
!! `ssp_rk2` runs first on the identical configuration as a control.
!! It has never had this defect (both sides sit on `u^n`, so they
!! cancel exactly), and requiring it to clear the same bar proves the
!! bar is achievable on this case rather than merely tight.
!!
!! ## The energy / HK transport forms on masked walls
!!
!! `energy_form_masked_walls_no_growth` runs the same basin with the
!! PRODUCTION wall convention (`mask_wall_velocity`, `wet_u = wet_v = 0`
!! on the four wall faces) under `&ocean_coriolis_nml form =
!! "sadourny_energy"` and `"sadourny_hk"`, from a jet that is NOT
!! masked at the wall faces.  `u_av` used to be seeded from that
!! unmasked state and was never rewritten at a masked face (the
!! renormaliser's `u_cor` is its only writer, and it skips walls and
!! land), so the wall value lived for the whole run.  The fast-loop
!! reference (`set_cor_ref_velocity` → `subtract_fast_cor_ref`) read it;
!! the transport forms did not (their `vh` rides `dx_cv = 0`); the
!! difference `a = −(f/4)·(v̄_av(i−1) + v̄_av(i))` forced every substep
!! of the two wall-adjacent u-rows.  Measured before the fix (gfortran,
!! 1200 steps): KE+PE ×126 (energy), ×129 (HK), with the time-integrated
!! work of `a` accounting for the gain to 0.1 %; the enstrophy form
!! read the stale value on BOTH sides, cancelled, and stayed at 0.89.
!! After (`mask_time_mean_velocities` on the `pred_corr` step-0 seed):
!! 0.892 for all three forms.
module test_ocean_cor_ref_seiche
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t, PV_VARIANT_SADOURNY_ENERGY, &
                               PV_VARIANT_SADOURNY_HK
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_LITE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, &
                            SPLIT_SCHEME_PRED_CORR, SPLIT_SCHEME_SSP_RK2
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_cor_ref_seiche_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 32
   integer, parameter :: NYP = 8
   integer, parameter :: NZ = 2
   integer, parameter :: N_INNER = 24
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: H0 = 200.0_wp        !! flat total depth (m)
   real(wp), parameter :: GRAV = 9.80665_wp
   real(wp), parameter :: JET_CELLS = 8.0_wp   !! seed wavelength (grid cells)
   real(wp), parameter :: PI_L = 3.14159265358979324_wp

   ! Deliberately SMALL seed amplitude: the reference residual is a
   ! LINEAR instability, so its energy GROWTH FACTOR does not care about
   ! the seed (measured: the same factor at V0 = 5e-3 and at 1e-3),
   ! while the nonlinear grid-scale cascade of an undamped zero-
   ! viscosity jet very much does.  Seeding small buys a long clean
   ! window in which the only thing that can grow is the defect.
   ! Rotation is stronger than Earth's and the outer step is long
   ! because the residual scales like `(f·dt)²` per outer step — that is
   ! what makes the signal decisive in ~1200 cheap steps.  The substep
   ! count still resolves the gravity CFL: c = sqrt(g·H0) = 44 m/s,
   ! dt_inner = 300/24 = 12.5 s, c·dt_inner·sqrt(2)/dx = 0.78.
   real(wp), parameter :: F0 = 2.0e-4_wp       !! f-plane Coriolis (1/s)
   real(wp), parameter :: DT = 300.0_wp        !! outer step (s)
   real(wp), parameter :: V0 = 1.0e-3_wp       !! seed jet amplitude (m/s)
   integer, parameter :: N_STEPS = 1200

   real(wp), parameter :: ENERGY_GROWTH_BAR = 2.0_wp
      !! Bar for `(KE+PE)_end / (KE+PE)_0`.  Measured on this case
      !! (gfortran host build, 1200 steps):
      !!   ssp_rk2 control                                    0.614
      !!   pred_corr, reference from `u_av`   (correct)         0.875
      !!   pred_corr, reference from `u^n`    (the defect)     126.3
      !! The bar sits BETWEEN the regimes — 2.3x of headroom below it,
      !! 63x above.  It is not a tuned threshold, and it must never be
      !! widened to make this pass: a ratio creeping up means the
      !! residual is back.

   real(wp), parameter :: OMEGA_SEICHE = PI_L*44.2865_wp/(real(NXP, wp)*DX)
      !! Gravest seiche frequency `π·c/L` (1/s), `c = √(g·H0) = 44.29 m/s`
      !! written as a literal (an intrinsic `sqrt` in a constant
      !! expression is F2008 but not every toolchain here folds it).
   real(wp), parameter :: ENERGY_NONINCREASE_BAR = &
                          1.0_wp + 0.5_wp*OMEGA_SEICHE*DT/real(N_INNER, wp)
      !! DERIVED bar `1 + ω·dt_inner/2 = 1.027` for the masked-wall
      !! transport-form gate — the derivation of
      !! `test_ocean_zfixed_cor_ref`'s bar on this basin: closed,
      !! unforced, inviscid, uniform density, so the only discrete energy
      !! sources are the time schemes.  `pred_corr` at `pc_be = 0.6`
      !! gives `|G|² = 1 − 0.2·θ² + 0.36·θ⁴ < 1` for every slow linear
      !! mode (`θ = f·dt = 0.06`); the forward-backward substep is
      !! neutral (`c·dt_inner·√2/dx = 0.78 < 2`) but conserves a MODIFIED
      !! energy, from which the sampled `KE+PE` of a wave at `ω` departs
      !! by at most `ω·dt_inner/2`, the gravest seiche being the largest
      !! `ω` that holds a finite share of this seed.  The seed is linear
      !! (`V0/(f·L) = 1.6e-4`), so no nonlinear allowance is owed.

contains

   subroutine collect_ocean_cor_ref_seiche_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cor_ref_rotating_basin_no_growth", test_rotating_basin_no_growth), &
                  new_unittest("energy_form_masked_walls_no_growth", test_energy_form_masked_walls) &
                  ]
   end subroutine collect_ocean_cor_ref_seiche_tests

   subroutine run_basin(split_scheme, finite, energy_ratio, pv_variant, solid_walls, stale_av)
      !! Integrate `N_STEPS` outer steps of the closed f-plane basin
      !! under `split_scheme` and report `(KE+PE)_end / (KE+PE)_0`.
      !!
      !! Everything that could hide (or fake) the signal is switched
      !! off: no wind, no drag, no lateral viscosity, no tracer
      !! diffusion, no vertical mixing, no thermodynamics, uniform
      !! density (so the PGF is identically zero) and a LAGRANGIAN
      !! vertical coordinate (so there is no remap).  `ocean_dyn_step_split`
      !! is called without a `bc`, which defaults every edge to
      !! `OBC_WALL` — a genuinely closed basin, which is what supports
      !! the seiche.
      integer, intent(in) :: split_scheme
      logical, intent(out) :: finite
      real(wp), intent(out) :: energy_ratio
      integer, intent(in), optional :: pv_variant
         !! Coriolis form (`PV_VARIANT_*`).  Absent ⇒ the `coriolis_adv_t`
         !! default (enstrophy), as the original case runs.
      logical, intent(in), optional :: solid_walls
         !! `.true.` ⇒ the production wall convention (`wet_u = wet_v = 0`
         !! on the wall faces).  Absent ⇒ the legacy unmasked walls, as
         !! the original case runs.
      real(wp), intent(out), optional :: stale_av
         !! max `|u_av|`, `|v_av|` over the faces with `wet = 0` after the
         !! run — the land contract says exactly zero.

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
      type(ocean_vcoord_t) :: vc

      integer :: i, j, k, ig, i0, i1, j0, j1, step
      real(wp) :: x_f, vjet, energy0, energy1
      logical :: walls

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = F0
      call cor%init(grid, nz_ml=NZ)
      if (present(pv_variant)) cor%pv_variant = pv_variant
      walls = .false.
      if (present(solid_walls)) walls = solid_walls
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_FV_LITE
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_LAGRANGIAN
      dyn%split_scheme = split_scheme

      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP

      ! Flat two-layer column, uniform density, depth-uniform
      ! (barotropic) sinusoidal v-seed in x.
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            ms%h_layer(i, j, 1) = 0.5_wp*H0
            ms%h_layer(i, j, 2) = 0.5_wp*H0
            ms%rho_layer(i, j, :) = eos%rho0
            do k = 1, NZ
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            x_f = real(i - ig, wp)
            vjet = V0*sin(2.0_wp*PI_L*x_f/JET_CELLS)
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = vjet
            dyn%bt_work%bt_H_ref(i, j) = H0
         end do
      end do
      ms%u_face_x_layer(grid%nx_total + 1, :, :) = 0.0_wp
      ms%v_face_y_layer(:, grid%ny_total + 1, :) = 0.0_wp

      energy0 = basin_energy(ms, i0, i1, j0, j1)

      ! The seed above is deliberately NOT masked: with `solid_walls` the
      ! wall faces carry the jet at step 0, which is what the transport-
      ! form gate needs (see the module docstring).
      call make_cartesian_metrics(metrics, grid, solid_walls=walls)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
      end do

      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      finite = .true.
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%v_face_y_layer(i, j, k))) finite = .false.
            end do
         end do
      end do

      if (present(stale_av)) then
         !$acc update self(ms%u_av_layer, ms%v_av_layer)
         stale_av = 0.0_wp
         do k = 1, NZ
            stale_av = max(stale_av, maxval(abs(ms%u_av_layer(:, :, k))*(1.0_wp - metrics%wet_u)))
            stale_av = max(stale_av, maxval(abs(ms%v_av_layer(:, :, k))*(1.0_wp - metrics%wet_v)))
         end do
      end if

      energy1 = basin_energy(ms, i0, i1, j0, j1)
      if (finite .and. energy0 > 0.0_wp .and. ieee_is_finite(energy1)) then
         energy_ratio = energy1/energy0
      else
         energy_ratio = huge(1.0_wp)
      end if

      call vc%exit_data()
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data()
      call ss%exit_data(); call bd%exit_data(); call hv%exit_data()
      call pgf%exit_data(); call cor%exit_data(); call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine run_basin

   pure function basin_energy(ms, i0, i1, j0, j1) result(energy)
      !! Depth-integrated `KE + PE` over the interior, per unit area and
      !! per unit density: `Σ_k h·(u² + v²)/2 + g·η²/2`, with
      !! `η = Σ_k h − H0` (flat bottom, `bt_H_ref = H0`).  Face
      !! velocities are squared in place — a monotone energy DIAGNOSTIC,
      !! not a discretely conserved energy, which is all a non-growth
      !! assertion needs.
      type(multilayer_state_t), intent(in) :: ms
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: energy
      integer :: i, j, k
      real(wp) :: eta, h_col

      energy = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            h_col = 0.0_wp
            do k = 1, NZ
               h_col = h_col + ms%h_layer(i, j, k)
               energy = energy + 0.5_wp*ms%h_layer(i, j, k)* &
                        (ms%u_face_x_layer(i, j, k)**2 + ms%v_face_y_layer(i, j, k)**2)
            end do
            eta = h_col - H0
            energy = energy + 0.5_wp*GRAV*eta*eta
         end do
      end do
   end function basin_energy

   subroutine test_rotating_basin_no_growth(error)
      !! `ssp_rk2` control first (must itself clear the bar), then
      !! `pred_corr` on the identical configuration.
      type(error_type), allocatable, intent(out) :: error
      logical :: fin_pc, fin_rk2
      real(wp) :: ratio_pc, ratio_rk2
      character(len=320) :: msg

      call run_basin(SPLIT_SCHEME_SSP_RK2, fin_rk2, ratio_rk2)
      call check(error, fin_rk2, "ssp_rk2 rotating basin: the control went non-finite")
      if (allocated(error)) return
      write (msg, '("ssp_rk2 rotating basin control: KE+PE ratio = ", es11.3, &
            &" (the control must itself not grow, else the case is unusable)")') ratio_rk2
      call check(error, ratio_rk2 < ENERGY_GROWTH_BAR, trim(msg))
      if (allocated(error)) return

      call run_basin(SPLIT_SCHEME_PRED_CORR, fin_pc, ratio_pc)
      call check(error, fin_pc, &
                 "pred_corr rotating basin: went non-finite (the seiche blew up)")
      if (allocated(error)) return

      write (msg, '("pred_corr rotating basin: KE+PE ratio = ", es11.3, " vs ssp_rk2 ", &
            &es11.3, ". Unforced + undamped, so growth is manufactured energy;", &
            &" suspect the fast-loop Coriolis reference (set_cor_ref_velocity).")') &
         ratio_pc, ratio_rk2
      call check(error, ratio_pc < ENERGY_GROWTH_BAR, trim(msg))
   end subroutine test_rotating_basin_no_growth

   subroutine test_energy_form_masked_walls(error)
      !! `pred_corr` + the TRANSPORT Coriolis forms (`sadourny_energy`,
      !! `sadourny_hk`) on the production masked walls, from a seed that
      !! is not masked at the wall faces: `KE+PE` must stay under the
      !! derived `ENERGY_NONINCREASE_BAR`, and `u_av` must be exactly
      !! zero on every masked face (a product with `wet = 0`, not a
      !! cancellation, so exact zero is the right assertion).
      !!
      !! Fails before / passes after, measured (gfortran 15.1 Release,
      !! 1200 steps): energy 1.262E+02 → 0.892, HK 1.287E+02 → 0.892;
      !! wall `v_av` 1.0E-03 (the seed amplitude, for the whole run) → 0.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: forms(2) = [PV_VARIANT_SADOURNY_ENERGY, PV_VARIANT_SADOURNY_HK]
      character(len=15), parameter :: names(2) = ["sadourny_energy", "sadourny_hk    "]
      logical :: fin
      real(wp) :: ratio, stale
      integer :: n
      character(len=640) :: msg

      do n = 1, size(forms)
         call run_basin(SPLIT_SCHEME_PRED_CORR, fin, ratio, pv_variant=forms(n), &
                        solid_walls=.true., stale_av=stale)
         call check(error, fin, "pred_corr + "//trim(names(n))//" masked-wall basin went non-finite")
         if (allocated(error)) return
         write (msg, '("pred_corr + ",a," on masked walls: KE+PE ratio = ",es11.3, &
               &" over ",I0," steps; derived bar ",f6.4,". Unforced + inviscid: growth is ", &
               &"manufactured. Suspect a reference/slow-Coriolis mismatch at the walls ", &
               &"(u_av vs the metric-masked transports). Do NOT widen this bar.")') &
            trim(names(n)), ratio, N_STEPS, ENERGY_NONINCREASE_BAR
         call check(error, ratio <= ENERGY_NONINCREASE_BAR, trim(msg))
         if (allocated(error)) return
         write (msg, '("pred_corr + ",a,": u_av on a masked face = ",es11.3, &
               &" after the run; the land contract says exactly 0. A non-zero value is ", &
               &"read by the fast-loop Coriolis reference but not by the transport forms.")') &
            trim(names(n)), stale
         call check(error, stale == 0.0_wp, trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_energy_form_masked_walls

end module test_ocean_cor_ref_seiche
