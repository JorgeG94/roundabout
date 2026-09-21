!! BAROTROPIC CONSISTENCY of the partial-step z-level face closure
!! (`&vcoord_nml zfixed_closed_faces`) — a closed-basin seiche over a
!! STAIRCASE bed under `VCOORD_Z_FIXED`.
!!
!! ### What this guards, and why a seiche
!!
!! Closing a face removes transport CAPACITY from the column, not water.
!! The barotropic mode therefore has to be told about it in three places
!! at once, and they have to agree:
!!
!!   1. `derive_bt_from_layers` must build `ubt` as the OPEN-column depth
!!      mean, because the fast loop transports on
!!      `ubt·FA·dy_cu_bt` with `dy_cu_bt` narrowed by the open fraction —
!!      so `ubt` there already MEANS "open-column mean";
!!   2. `face_depth_mean_*` must depth-average the slow tendencies under
!!      the SAME weights, or the `dt·F_bt` the fold subtracts back out is
!!      not the quantity the fast loop integrated;
!!   3. `apply_bt_correction` must distribute `Δu` over the OPEN layers,
!!      `wt = h_face·open·vr / h_bar_o`, so the open-column mean is
!!      shifted by exactly `Δu` and a CLOSED layer receives nothing.
!!
!! The throwaway spike did (3) as fold-uniformly-then-mask, which leaves
!! the layer depth mean short of `ubt_end` by `Δu·(Σ_closed h)/(Σ_k h)`
!! — a cancellation, not a construction, and invisible in a case that
!! never leaves rest.  A SEICHE is the cheapest state that is not at rest
!! and has no energy source: an unforced, undamped closed basin sloshing
!! at `√(gH)`.  If any of the three disagrees, the disagreement is a
!! per-face velocity bias that the next step re-derives `ubt` from, and
!! it compounds — which shows up as energy growth even though nothing is
!! doing work.
!!
!! ### The geometry
!!
!! A closed basin, flat-bottomed in the western half and a STEP shallower
!! in the eastern half, run on `z_fixed` with `h_nominal` chosen so the
!! step costs the shallow columns exactly one nominal layer.  The mask
!! then closes the bed layer at every face inside and bordering the
!! shallow half — the test asserts that a partially closed face exists
!! before it asserts anything about it, so a future change that silently
!! stops closing faces cannot make this pass vacuously.
!!
!! Unforced and undamped: no wind, no drag, no lateral viscosity, no
!! tracer diffusion, no vertical mixing, no thermodynamics, uniform
!! density.
!!
!! ### The assertions
!!
!! **`fold_preserves_open_depth_mean`** is the transport identity, and it
!! is asserted where it is EXACT: on the fold itself.  A hand-built face
!! column with a known mask is folded and masked, `ubt` is re-derived from
!! the result by the very routine the next step will use, and
!!
!! ```
!! Σ_k h_face·open·u_k = (Σ_k h_face·open)·ubt_end
!! ```
!! is checked to a DERIVED round-off bound — never bit-zero: both sides
!! are cancellations of `nz` products, the compiler is free to contract
!! `h*u` into an FMA on one side and not the other, and a test that
!! asserts `== 0` across that is a test of the optimiser.  The SPIKE's
!! combination — full-column `derive_bt_from_layers` against a
!! uniform-then-masked fold — is then run on the SAME case, and the test
!! asserts it MISSES by exactly `Δu·(Σ_closed h)/(Σ_k h)`: a gate that
!! cannot discriminate is not a gate.  One further exact statement pins
!! the constructive difference: the open-layer fold leaves a closed
!! layer at zero BEFORE the mask runs (`wt = 0` there), where the uniform
!! fold leaves `Δu` in it and needs the mask to take it back out.
!!
!! The integration cases then assert, over a real run:
!!
!!   * **a closed face carries exactly zero normal velocity at the end of
!!     EVERY step** — with `==`, because that one IS exact (a multiply by
!!     a literal 0).  It is an end-of-STEP statement, not end-of-stage:
!!     the ALE remap runs after the last stage's mask and is itself a
!!     velocity writer, which is why the mask is re-asserted after it.
!!   * **no energy growth over >= 20 seiche periods**, and
!!   * **the measured period within a few % of the analytic `2L/√(gH)`** —
!!     the barotropic mode has to see the OPEN depth of the staircase, and
!!     a mode that does not is a mode with the wrong phase speed.
!!
!! The per-step END-OF-STEP comparison of `mean_open(u)` against
!! `bt_ubt_end` is deliberately NOT asserted: the ALE remap runs between
!! the fold and the step boundary and redistributes momentum along the
!! face column, so the two are not equal there to round-off and pretending
!! otherwise would need a tolerance nobody could derive.
module test_ocean_zfixed_bt_seiche
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_closed_faces_alloc
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
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
                            SPLIT_SCHEME_PRED_CORR, SPLIT_SCHEME_SSP_RK2, &
                            mask_layer_velocities
   use rdb_barotropic_coupling, only: apply_bt_correction, derive_bt_from_layers
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_Z_FIXED, &
                               ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_closed_face_masks
   use testdrive, only: error_type, check, new_unittest, unittest_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_zfixed_bt_seiche_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 40
   integer, parameter :: NYP = 4
   integer, parameter :: NZ = 4
   integer, parameter :: N_INNER = 20
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: GRAV = 9.80665_wp

   real(wp), parameter :: H_DEEP = 400.0_wp
      !! Western (deep) column depth (m).  With `h_nominal = 100 m` the
      !! deep column is 4 full nominal layers.
   real(wp), parameter :: H_SHELF = 300.0_wp
      !! Eastern (shelf) column depth (m) — exactly one nominal layer
      !! shallower, so the `z_fixed` target vanishes the BED layer there
      !! and every face touching the shelf closes at `k = 1`.
   real(wp), parameter :: H_NOM = 100.0_wp
   real(wp), parameter :: H_MIN = 1.0e-4_wp

   real(wp), parameter :: ETA0 = 0.05_wp
      !! Seiche amplitude (m).  Small enough that the mode stays linear
      !! (`ETA0/H_DEEP ~ 1e-4`), large enough to sit far above round-off.
   real(wp), parameter :: PI_L = 3.14159265358979324_wp
   real(wp), parameter :: DT = 60.0_wp

   integer, parameter :: N_PERIODS = 22
      !! At least the 20 the gate asks for, with headroom.

   real(wp), parameter :: ENERGY_GROWTH_BAR = 1.30_wp
      !! Bar for `(KE+PE)_end / (KE+PE)_0` over `N_PERIODS` periods.
      !! MEASURED on this case (gfortran 15.1 host build, 22 periods):
      !!   `pred_corr`  0.990   — the mode propagates for the whole run
      !!                          and is very slightly damped
      !!   `ssp_rk2`    1.8E-07 — the SSP average damps THIS mode away
      !!                          almost completely.  That is a property
      !!                          of the outer scheme, not of the mask,
      !!                          and it is why the bar is one-sided: the
      !!                          `ssp_rk2` leg is here to prove the mask
      !!                          does not depend on the split, not to
      !!                          carry the energy statement.
      !! The bar is 31 % of headroom over the `pred_corr` leg, which is
      !! the one that actually carries the mode.  It must never be widened
      !! to make this pass: a ratio creeping toward it means the fold,
      !! `derive_bt_from_layers` and `face_depth_mean_*` have stopped
      !! agreeing about what `ubt` is, and that disagreement compounds
      !! every step.

   real(wp), parameter :: PERIOD_TOL = 0.05_wp
      !! Fractional tolerance on the measured seiche period against the
      !! analytic `2L/√(g·H_eff)`.  MEASURED 1.27 % apart under both outer
      !! schemes, so 5 % is ~4x of headroom.  A tolerance at all (rather
      !! than a tight one) because `H_eff` over a two-step
      !! bed is not a single number — the analytic value below uses the
      !! width-weighted mean depth, which is the standard first-order
      !! estimate for a stepped basin, not an exact eigenvalue.

contains

   subroutine collect_ocean_zfixed_bt_seiche_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fold_preserves_open_depth_mean", test_fold_identity), &
                  new_unittest("staircase_seiche_pred_corr", test_pred_corr), &
                  new_unittest("staircase_seiche_ssp_rk2", test_ssp_rk2) &
                  ]
   end subroutine collect_ocean_zfixed_bt_seiche_tests

   subroutine run_seiche(split_scheme, finite, energy_ratio, period_meas, &
                         n_closed, closed_u_max)
      !! Integrate the staircase basin and report everything the
      !! assertions need.
      integer, intent(in) :: split_scheme
      logical, intent(out) :: finite
      real(wp), intent(out) :: energy_ratio
         !! `(KE+PE)_end / (KE+PE)_0`.
      real(wp), intent(out) :: period_meas
         !! Seiche period (s) measured from the sign changes of the
         !! basin-ends SSH difference.
      integer, intent(out) :: n_closed
         !! Number of CLOSED interior u-face entries in the mask — the
         !! anti-vacuity guard.
      real(wp), intent(out) :: closed_u_max
         !! Worst `|u|` ever seen on a CLOSED u-face over the whole run.

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

      integer :: i, j, k, ig, i0, i1, j0, j1, step, n_steps, n_cross
      integer :: nx_t, ny_t
      real(wp) :: energy0, energy1, h_bed, x_rel, eta_seed
      real(wp) :: h_eff, period_analytic, t_prev_cross, t_first_cross, t_last_cross
      real(wp) :: d_prev, d_now
      real(wp), allocatable :: tgt(:, :, :), tot_h(:, :), eta0f(:, :), z_top(:, :)

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 0.0_wp
      call cor%init(grid, nz_ml=NZ)
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
      vc%coord_type = VCOORD_Z_FIXED
      vc%z_fixed_h_ref = H_NOM*real(NZ, wp)
      vc%zstar_h_min = H_MIN
      vc%zfixed_closed_faces = .true.
      dyn%split_scheme = split_scheme

      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      vd%zlevel_faces = .true.
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)

      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP

      call make_cartesian_metrics(metrics, grid)

      ! ---- the staircase bed, and the z_fixed target it implies ----
      allocate (tot_h(nx_t, ny_t), source=H_DEEP)
      allocate (eta0f(nx_t, ny_t), source=0.0_wp)
      allocate (z_top(nx_t, ny_t), source=0.0_wp)
      allocate (tgt(nx_t, ny_t, NZ), source=0.0_wp)
      do j = 1, ny_t
         do i = 1, nx_t
            if (i > ig + NXP/2) then
               tot_h(i, j) = H_SHELF
            else
               tot_h(i, j) = H_DEEP
            end if
         end do
      end do
      call ocean_vcoord_z_fixed_target(tgt, tot_h, eta0f, z_top, &
                                       nx_t, ny_t, NZ, H_NOM, H_MIN)

      call metrics_closed_faces_alloc(metrics, grid, NZ)
      call ocean_vcoord_closed_face_masks(metrics%open_u, metrics%open_v, &
                                          tgt, nx_t, ny_t, NZ, H_VANISHED)
      metrics%use_closed_faces = .true.

      n_closed = 0
      do k = 1, NZ
         do j = j0, j1
            do i = i0, i1 + 1
               if (metrics%open_u(i, j, k) == 0.0_wp) n_closed = n_closed + 1
            end do
         end do
      end do

      ! ---- seed: the gravest seiche mode, h from the target ----
      ! eta = ETA0*cos(pi*x/L) at rest; the basin then sloshes.  The
      ! anomaly rides the FIRST LIVE layer (the partial bed cell absorbs
      ! eta under z_fixed anyway, but the surface layer is the honest
      ! place to put a free-surface anomaly).
      do j = 1, ny_t
         do i = 1, nx_t
            x_rel = (real(i - ig, wp) - 0.5_wp)/real(NXP, wp)
            eta_seed = ETA0*cos(PI_L*x_rel)
            do k = 1, NZ
               ms%h_layer(i, j, k) = tgt(i, j, k)
            end do
            ms%h_layer(i, j, NZ) = ms%h_layer(i, j, NZ) + eta_seed
            ms%rho_layer(i, j, :) = eos%rho0
            do k = 1, NZ
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            ms%u_face_x_layer(i, j, :) = 0.0_wp
            ms%v_face_y_layer(i, j, :) = 0.0_wp
            dyn%bt_work%bt_H_ref(i, j) = tot_h(i, j)
         end do
      end do
      ms%u_face_x_layer(nx_t + 1, :, :) = 0.0_wp
      ms%v_face_y_layer(:, ny_t + 1, :) = 0.0_wp

      ! Width-weighted mean depth — the standard first-order estimate for
      ! a stepped basin's gravest mode.
      h_eff = 0.5_wp*(H_DEEP + H_SHELF)
      period_analytic = 2.0_wp*real(NXP, wp)*DX/sqrt(GRAV*h_eff)
      n_steps = nint(real(N_PERIODS, wp)*period_analytic/DT)

      energy0 = basin_energy(ms, dyn, i0, i1, j0, j1)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      closed_u_max = 0.0_wp
      n_cross = 0
      t_first_cross = 0.0_wp
      t_last_cross = 0.0_wp
      t_prev_cross = 0.0_wp
      d_prev = ssh_tilt(ms, dyn, i0, i1, j0, j1)

      do step = 1, n_steps
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
         ! Host read-back: the assertions below are HOST scans, so pull the
         ! prognostics and the barotropic velocity down every step.  The
         ! COMPONENT arrays only -- never the aggregate derived type (that
         ! overwrites the host descriptors with device addresses).
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         call scan_closed_faces(ms, metrics, i0, i1, j0, j1, closed_u_max)
         d_now = ssh_tilt(ms, dyn, i0, i1, j0, j1)
         if (d_prev*d_now < 0.0_wp) then
            n_cross = n_cross + 1
            t_prev_cross = real(step, wp)*DT
            if (n_cross == 1) t_first_cross = t_prev_cross
            t_last_cross = t_prev_cross
         end if
         d_prev = d_now
      end do

      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      finite = .true.
      do j = j0, j1
         do i = i0, i1
            do k = 1, NZ
               if (.not. ieee_is_finite(ms%h_layer(i, j, k))) finite = .false.
               if (.not. ieee_is_finite(ms%u_face_x_layer(i, j, k))) finite = .false.
            end do
         end do
      end do

      energy1 = basin_energy(ms, dyn, i0, i1, j0, j1)
      if (finite .and. energy0 > 0.0_wp .and. ieee_is_finite(energy1)) then
         energy_ratio = energy1/energy0
      else
         energy_ratio = huge(1.0_wp)
      end if

      ! Two sign changes per period.
      if (n_cross >= 3) then
         period_meas = 2.0_wp*(t_last_cross - t_first_cross)/real(n_cross - 1, wp)
      else
         period_meas = 0.0_wp
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
      deallocate (tgt, tot_h, eta0f, z_top)
   end subroutine run_seiche

   pure subroutine scan_closed_faces(ms, metrics, i0, i1, j0, j1, closed_u_max)
      !! Worst `|u|` on a CLOSED u-face.  A closed face is a z-level WALL
      !! and the mask multiplies by a literal 0, so the only acceptable
      !! answer is EXACTLY zero; anything else means a velocity writer ran
      !! after the last mask.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: i0, i1, j0, j1
      real(wp), intent(inout) :: closed_u_max
      integer :: i, j, k

      do j = j0, j1
         do i = i0 + 1, i1
            do k = 1, NZ
               if (metrics%open_u(i, j, k) == 0.0_wp) then
                  closed_u_max = max(closed_u_max, abs(ms%u_face_x_layer(i, j, k)))
               end if
            end do
         end do
      end do
   end subroutine scan_closed_faces

   pure function ssh_tilt(ms, dyn, i0, i1, j0, j1) result(tilt)
      !! `eta(west end) − eta(east end)`, the gravest-mode amplitude.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_dyn_t), intent(in) :: dyn
      integer, intent(in) :: i0, i1, j0, j1
      real(wp) :: tilt
      integer :: j, k
      real(wp) :: w, e

      w = 0.0_wp
      e = 0.0_wp
      do j = j0, j1
         do k = 1, NZ
            w = w + ms%h_layer(i0, j, k)
            e = e + ms%h_layer(i1, j, k)
         end do
         w = w - dyn%bt_work%bt_H_ref(i0, j)
         e = e - dyn%bt_work%bt_H_ref(i1, j)
      end do
      tilt = w - e
   end function ssh_tilt

   pure function basin_energy(ms, dyn, i0, i1, j0, j1) result(energy)
      !! Depth-integrated `KE + PE` per unit area and per unit density:
      !! `Σ_k h·(u² + v²)/2 + g·η²/2`, with `η = Σ_k h − bt_H_ref`.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_dyn_t), intent(in) :: dyn
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
            eta = h_col - dyn%bt_work%bt_H_ref(i, j)
            energy = energy + 0.5_wp*GRAV*eta*eta
         end do
      end do
   end function basin_energy

   subroutine one_scheme(error, split_scheme, label)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: split_scheme
      character(len=*), intent(in) :: label
      logical :: fin
      real(wp) :: ratio, period_meas, closed_u_max
      real(wp) :: h_eff, period_analytic, rel
      integer :: n_closed
      character(len=400) :: msg

      call run_seiche(split_scheme, fin, ratio, period_meas, n_closed, closed_u_max)

      ! ---- anti-vacuity: the mask must actually close something ----
      write (msg, '(A,": the staircase closed ",I0," interior u-face entries. ", &
            &"Zero means the geometry stopped producing a staircase and every ", &
            &"assertion below is vacuous -- fix the GEOMETRY, not the bar.")') &
         label, n_closed
      call check(error, n_closed > 0, trim(msg))
      if (allocated(error)) return

      write (msg, '(A,": went non-finite over the seiche run")') label
      call check(error, fin, trim(msg))
      if (allocated(error)) return

      ! ---- a closed face is a WALL: exactly zero, asserted with == ----
      write (msg, '(A,": max |u| on a CLOSED face = ",es11.3,". A closed face is a ", &
            &"z-level WALL and mask_layer_velocities multiplies by a literal 0, so ", &
            &"this is EXACT -- anything non-zero means a fold, a remap or a ", &
            &"tendency is writing velocity after the mask ran.")') label, closed_u_max
      call check(error, closed_u_max == 0.0_wp, trim(msg))
      if (allocated(error)) return

      ! ---- no energy growth ----
      write (msg, '(A,": KE+PE ratio = ",es11.3," over ",I0," seiche periods. ", &
            &"Unforced and undamped, so any growth is MANUFACTURED energy -- ", &
            &"the fold, derive_bt_from_layers and face_depth_mean_* must all ", &
            &"weight by h_face*open or their disagreement compounds per step.")') &
         label, ratio, N_PERIODS
      call check(error, ratio < ENERGY_GROWTH_BAR, trim(msg))
      if (allocated(error)) return

      ! ---- the period is the physical one ----
      h_eff = 0.5_wp*(H_DEEP + H_SHELF)
      period_analytic = 2.0_wp*real(NXP, wp)*DX/sqrt(GRAV*h_eff)
      rel = abs(period_meas - period_analytic)/period_analytic
      write (msg, '(A,": seiche period ",es11.3," s vs analytic 2L/sqrt(g*H_eff) = ", &
            &es11.3," s (",f6.2,"% apart). A period this far off means the ", &
            &"barotropic mode is not seeing the open depth of the staircase.")') &
         label, period_meas, period_analytic, 100.0_wp*rel
      call check(error, period_meas > 0.0_wp .and. rel < PERIOD_TOL, trim(msg))
   end subroutine one_scheme

   subroutine test_fold_identity(error)
      !! THE transport identity, asserted where it is exact: on the
      !! barotropic fold itself.
      !!
      !! One u-face column, `nz = 4`, layer 1 CLOSED but carrying 40 m of
      !! real water (the mask, not the thickness, is what closes a face —
      !! a test whose closed layer was a `1e-4` filler could not see the
      !! defect at all).  A known pre-fold profile and a known `Δu` go in,
      !! and the identity
      !!
      !! ```
      !! Σ_k h_face·open·u_k = (Σ_k h_face·open)·ubt_end
      !! ```
      !!
      !! is checked by re-deriving `ubt` from the folded layers with
      !! `derive_bt_from_layers(..., metrics)` — the very routine the next
      !! step will use — and comparing it against the `ubt_end` the fast
      !! loop handed over.  Round-off bound, DERIVED, never bit-zero: both
      !! sides are cancellations of `nz` products and the compiler is free
      !! to contract `h*u` into an FMA on one side and not the other.
      !!
      !! Then the SPIKE's combination is reproduced — full-column
      !! `derive_bt_from_layers` against a uniform-then-masked fold — and
      !! the test asserts it MISSES by exactly
      !! `Δu·(Σ_closed h)/(Σ_k h)`.  A gate that cannot discriminate is not
      !! a gate.
      !!
      !! One more exact statement, and it is the constructive difference
      !! between the two folds: the OPEN-layer fold leaves a closed
      !! layer's velocity at zero BEFORE the mask runs, because `wt = 0`
      !! there.  The uniform fold leaves `Δu` in it and relies on the mask
      !! to take it back out — preservation by cancellation, which is what
      !! the spike write-up itself said a production version must not do.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_dyn_t) :: dyn
      integer :: i, j, k, ig, iface, jface, nx_t, ny_t
      real(wp) :: du_want, resid, bound, umax, mean_open0, mean_full0
      real(wp) :: ubt_open, ubt_full, shortfall_want, shortfall_got
      real(wp) :: u_closed_open_fold, u_closed_uniform_fold
      character(len=420) :: msg
      real(wp), parameter :: SAFETY = 64.0_wp
      real(wp), parameter :: H_CLOSED = 40.0_wp
      real(wp), parameter :: H_TOTAL = 400.0_wp

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      call make_cartesian_metrics(metrics, grid)
      call metrics_closed_faces_alloc(metrics, grid, NZ)
      metrics%open_u = 1.0_wp
      metrics%open_v = 1.0_wp
      metrics%open_u(:, :, 1) = 0.0_wp
      metrics%open_v(:, :, 1) = 0.0_wp
      metrics%use_closed_faces = .true.

      iface = ig + 5
      jface = ig + 2

      do j = 1, ny_t
         do i = 1, nx_t
            ms%h_layer(i, j, 1) = H_CLOSED
            ms%h_layer(i, j, 2) = 100.0_wp
            ms%h_layer(i, j, 3) = 120.0_wp
            ms%h_layer(i, j, 4) = 140.0_wp
            dyn%bt_work%bt_H_ref(i, j) = H_TOTAL
         end do
      end do
      call seed_profile(ms, nx_t, ny_t)

      call column_means(ms, metrics, iface, jface, mean_open0, mean_full0, umax)

      ! The fold's Δu is `ubt_end − ubt_at_n − dt·F_bt`.  Anchor `ubt_at_n`
      ! on the PRE-FOLD open mean and zero the forcing, so the identity
      ! under test — post-fold open mean == ubt_end — is a statement about
      ! the FOLD and not about whatever the slow tendencies did.
      du_want = 0.037_wp
      dyn%bt_work%ubt_at_n = mean_open0
      dyn%bt_work%F_bt_u = 0.0_wp
      dyn%bt_work%vbt_at_n = 0.0_wp
      dyn%bt_work%F_bt_v = 0.0_wp
      dyn%bt_work%bt_ubt_end = mean_open0 + du_want
      dyn%bt_work%bt_vbt_end = 0.0_wp
      dyn%bt_work%bt_eta_end = 0.0_wp

      ! ---- the OPEN-layer fold ----
      call apply_bt_correction(dyn%bt_work, ms, 1.0_wp, skip_h_rescale=.true., &
                               metrics=metrics)
      u_closed_open_fold = ms%u_face_x_layer(iface, jface, 1)
      call mask_layer_velocities(grid, metrics, ms)
      call derive_bt_from_layers(grid, dyn%bt_work, ms, metrics)
      ubt_open = dyn%bt_work%bt_ubt(iface, jface)

      resid = abs(ubt_open - (mean_open0 + du_want))
      bound = SAFETY*real(NZ, wp)*epsilon(1.0_wp)*max(abs(mean_open0 + du_want), umax)
      write (msg, '("open-layer fold: sum_k h*open*u / sum_k h*open = ",es20.13, &
            &" but ubt_end = ",es20.13," (residual ",es11.3," vs derived bound ", &
            &es11.3,"). The folded layer velocities and the barotropic velocity ", &
            &"disagree about the OPEN column -- suspect apply_bt_correction''s ", &
            &"open branch or derive_bt_from_layers'' open weighting.")') &
         ubt_open, mean_open0 + du_want, resid, bound
      call check(error, resid <= bound, trim(msg))
      if (allocated(error)) then
         call teardown(ms, metrics)
         return
      end if

      ! ---- a CLOSED layer receives nothing, BEFORE the mask ----
      write (msg, '("the open-layer fold left ",es11.3," m/s in a CLOSED layer ", &
            &"before mask_layer_velocities ran. It must be EXACTLY zero: wt = 0 ", &
            &"there by construction, and relying on the mask to take it back out ", &
            &"is preservation by cancellation -- the thing the spike write-up ", &
            &"said a production version must not do.")') u_closed_open_fold
      call check(error, u_closed_open_fold == 0.0_wp, trim(msg))
      if (allocated(error)) then
         call teardown(ms, metrics)
         return
      end if

      ! ---- the SPIKE's combination, for contrast ----
      call seed_profile(ms, nx_t, ny_t)
      ! Uniform Δu into every layer (no `metrics` => no open branch)...
      call apply_bt_correction(dyn%bt_work, ms, 1.0_wp, skip_h_rescale=.true.)
      u_closed_uniform_fold = ms%u_face_x_layer(iface, jface, 1)
      ! ...then the mask takes it back out of the closed ones...
      call mask_layer_velocities(grid, metrics, ms)
      ! ...and `ubt` is re-derived over the FULL column, as the spike did.
      call derive_bt_from_layers(grid, dyn%bt_work, ms)
      ubt_full = dyn%bt_work%bt_ubt(iface, jface)

      write (msg, '("the uniform fold should leave Du = ",es11.3," m/s in the ", &
            &"CLOSED layer before masking; measured ",es11.3,". If it is already ", &
            &"zero the two folds have become the same and this contrast is ", &
            &"vacuous.")') du_want, u_closed_uniform_fold
      call check(error, abs(u_closed_uniform_fold - du_want) <= 1.0e-14_wp, trim(msg))
      if (allocated(error)) then
         call teardown(ms, metrics)
         return
      end if

      ! The spike's `ubt` is the FULL-column mean of a field whose closed
      ! layers were just zeroed, so it falls short of `ubt_end` by exactly
      ! `Du*(sum_closed h)/(sum_k h)` -- plus the pre-fold profile's own
      ! open-vs-full mean difference, which is why the comparison is made
      ! against `mean_full0 + Du*(1 - f)`.
      shortfall_want = du_want*(H_CLOSED/H_TOTAL)
      shortfall_got = (mean_full0 + du_want) - ubt_full
      write (msg, '("the spike combination (full-column derive_bt + uniform ", &
            &"fold + mask) should MISS ubt_end by Du*(sum_closed h)/(sum_k h) = ", &
            &es11.3,"; measured ",es11.3,". If the miss has vanished, the ", &
            &"open weighting is no longer the thing being tested.")') &
         shortfall_want, shortfall_got
      call check(error, abs(shortfall_got - shortfall_want) <= 1.0e-12_wp, trim(msg))

      call teardown(ms, metrics)
   end subroutine test_fold_identity

   pure subroutine seed_profile(ms, nx_t, ny_t)
      !! A depth-varying pre-fold velocity, with the CLOSED layer at rest
      !! (which is where `mask_layer_velocities` would have left it at the
      !! end of the previous stage).
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx_t, ny_t
      integer :: i, j, k

      do k = 1, NZ
         do j = 1, ny_t
            do i = 1, nx_t + 1
               ms%u_face_x_layer(i, j, k) = 0.01_wp*real(k, wp) - 0.02_wp
            end do
         end do
      end do
      ms%v_face_y_layer = 0.0_wp
      ms%u_face_x_layer(:, :, 1) = 0.0_wp
   end subroutine seed_profile

   pure subroutine column_means(ms, metrics, iface, jface, mean_open, mean_full, umax)
      !! The OPEN-weighted and FULL-column depth means at one u-face, plus
      !! the velocity magnitude the round-off bound is derived from.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: iface, jface
      real(wp), intent(out) :: mean_open, mean_full, umax
      integer :: k
      real(wp) :: h_face, w, num_o, den_o, num_f, den_f

      num_o = 0.0_wp
      den_o = 0.0_wp
      num_f = 0.0_wp
      den_f = 0.0_wp
      umax = 0.0_wp
      do k = 1, NZ
         h_face = 0.5_wp*(ms%h_layer(iface - 1, jface, k) + ms%h_layer(iface, jface, k))
         w = h_face*metrics%open_u(iface, jface, k)
         num_o = num_o + w*ms%u_face_x_layer(iface, jface, k)
         den_o = den_o + w
         num_f = num_f + h_face*ms%u_face_x_layer(iface, jface, k)
         den_f = den_f + h_face
         umax = max(umax, abs(ms%u_face_x_layer(iface, jface, k)))
      end do
      mean_open = num_o/den_o
      mean_full = num_f/den_f
   end subroutine column_means

   subroutine teardown(ms, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine teardown

   subroutine test_pred_corr(error)
      type(error_type), allocatable, intent(out) :: error
      call one_scheme(error, SPLIT_SCHEME_PRED_CORR, "pred_corr staircase seiche")
   end subroutine test_pred_corr

   subroutine test_ssp_rk2(error)
      !! The mask must not depend on the outer time split.
      type(error_type), allocatable, intent(out) :: error
      call one_scheme(error, SPLIT_SCHEME_SSP_RK2, "ssp_rk2 staircase seiche")
   end subroutine test_ssp_rk2

end module test_ocean_zfixed_bt_seiche
