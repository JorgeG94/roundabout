!! THE CLOSED-FACE Coriolis-reference gate: `&vcoord_nml
!! zfixed_closed_faces` under `VCOORD_Z_FIXED`, with ROTATION, on the
!! `pred_corr` outer split.
!!
!! ### Why this test exists, given the two tests either side of it
!!
!! `tests/test_ocean_cor_ref_seiche.F90` already guards the `pred_corr`
!! barotropic Coriolis reference — and it CANNOT see the defect this
!! file exists for, because it is ALL-OPEN.  With `φ ≡ 1` a full-column
!! depth mean and an open-column depth mean are the same number, so
!! `set_cor_ref_velocity` omitting the closed-face weights is invisible
!! there.  `tests/test_ocean_zfixed_bt_seiche.F90` has the staircase and
!! the mask but runs at `f = 0`, and the residual this guards carries a
!! factor `f`: it is dead in that case too.  The defect needed BOTH —
!! closed faces AND rotation — and so needs a third case.
!!
!! ### The defect, stated
!!
!! A CLOSED layer carries exactly zero velocity (`mask_layer_velocities`)
!! but a NON-ZERO `h_face`.  Writing
!! `φ_u = Σ_k h_face·open_u / Σ_k h_face` for the open fraction of a
!! u-face, a FULL-column depth mean of `u_av` returns
!!
!! ```
!! Σ_k u_av·h_face / Σ_k h_face  =  φ_u · [ Σ_k u_av·h_face·open / Σ_k h_face·open ]
!!                               =  φ_u · ū_open
!! ```
!!
!! while the fast loop integrates its live `(ζ+f)·v̄ − ∇KE` on
!! `bt_ubt = ū_open` and `F_bt` is the OPEN-weighted mean of the slow
!! `cor%pv_flux_*`.  `subtract_fast_cor_ref` then removes `f·φ_v·v̄`
!! where it must remove `f·v̄`, leaving
!!
!! ```
!! Δa_u = + f·(1 − φ_v)·v̄ ,        Δa_v = − f·(1 − φ_u)·ū
!! ```
!!
!! forcing EVERY barotropic substep **proportionally to the barotropic
!! velocity itself** — an amplifier, not a seed, which is why the answer
!! is exponential and why `f` appears steeply.  On a partial-step face
!! `φ` is O(0.5), not O(1 − 1e-4).
!!
!! ### The three assertions, what each one gates, and which two carry it
!!
!! **`cor_ref_is_open_column_mean`** — the DIRECT statement, on a
!! hand-built column with a known mask: `cor_ref_u` from
!! `set_cor_ref_velocity(..., from_u_av = .true.)` must equal
!! `Σ_k h_face·open·u_av / Σ_k h_face·open` to a DERIVED round-off bound
!! — and it asserts in the same breath that the open and full means
!! DIFFER on this column, so a future all-open column cannot make it
!! pass vacuously.  Measured: `1.1111111111111E-02` (correct) against
!! `1.0000000000000E-02` (the full-column answer) on the same column.
!!
!! **`bt_coriolis_residual_is_roundoff`** — the MECHANISM, in its own
!! units.  On a state that is EXACTLY barotropic on the open layers the
!! two branches of `set_cor_ref_velocity` must produce the same
!! reference, so `subtract_fast_cor_ref` must remove the same forcing
!! either way; the difference IS the spurious acceleration the substep
!! then integrates.  Measured: `2.500E-07 m/s²` before the fix —
!! matching the derived `f·(1−φ)·V = 2e-4 · 0.25 · 5e-3` to four figures
!! — against a `5.7E-20` round-off bound after it.
!!
!! Those two are the gate.  The third is not, and says so:
!!
!! **`rotating_staircase_no_growth`** — the integration case (closed
!! rotating basin, no wind, drag, viscosity, tracer diffusion, vertical
!! mixing or thermodynamics, uniform density) asserting that `KE + PE`
!! does not grow over `N_PERIODS` seiche periods.  It is a CATASTROPHE
!! guard, NOT the discriminator, and the distinction is measured rather
!! than assumed: at this horizon it reads ~1.0 both with and without the
!! fix, and pushed out to 2000 outer steps it reads 25x fixed against
!! 55x unfixed — still growing on BOTH.  The reason is a SECOND,
!! independent amplifier on the `pred_corr` × closed-faces × rotation
!! path that this fix does not touch: it is linear (seed-independent
!! from `V0 = 1e-3` down to `1e-6`), absent at `f = 0` (ratio 0.97),
!! absent with the mask all-open (1.11 at 1700 steps) and absent under
!! `ssp_rk2` (0.69).  It is a separate finding and a separate slice; it
!! is recorded here so that nobody reads this assertion as evidence the
!! closed-face barotropic path is clean.  Do NOT tighten this bar
!! expecting it to catch the reference defect — the two above do that,
!! exactly and cheaply.
!!
!! `ssp_rk2` runs first as a control.  It takes the OTHER branch of
!! `set_cor_ref_velocity` — a plain COPY of `bt_ubt`, which IS the
!! open-column mean — so it cannot have the defect by construction, and
!! requiring it to clear the same bar proves the bar is achievable on
!! this case rather than merely tight.
module test_ocean_zfixed_cor_ref
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_closed_faces_alloc
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: set_cor_ref_velocity, subtract_fast_cor_ref, &
                                      derive_bt_from_layers
   use rdb_ocean_boundary_types, only: OBC_WALL
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, ocean_porous_refresh, &
                            SPLIT_SCHEME_PRED_CORR, SPLIT_SCHEME_SSP_RK2
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_Z_FIXED, &
                               ocean_vcoord_z_fixed_target
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: collect_ocean_zfixed_cor_ref_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 32
   integer, parameter :: NYP = 16
   integer, parameter :: NZ = 4
   integer, parameter :: N_INNER = 20
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: GRAV = 9.80665_wp

   real(wp), parameter :: H_NOM = 100.0_wp
   real(wp), parameter :: H_MIN = 1.0e-4_wp
   real(wp), parameter :: H0 = H_NOM*real(NZ, wp)
      !! FLAT total depth (m) — `NZ` full nominal layers, so the `z_fixed`
      !! target is exactly the live column, no fillers anywhere, and the
      !! ALE regrid has nothing to move.  The bed-layer wall is a
      !! DIAGONAL LEDGE imposed on `open_u`/`open_v` directly (see
      !! `run_basin`) rather than implied by a stepped bathymetry.  Both
      !! halves of that sentence are deliberate:
      !!
      !!   * the CLOSED layer must carry REAL water.  A layer closed
      !!     because it is an inert `H_MIN = 1e-4` filler leaves
      !!     `φ = Σ h·open / Σ h = 1 − 3e-7`, and the residual is `∝ (1−φ)`
      !!     — it is unreachable there.  Here the closed layer is 100 m of
      !!     water and `φ = 0.75`.  (On the real ISOMIP+ geometry the same
      !!     thing happens through PARTIAL steps: a face with 100 m live on
      !!     one side and a filler on the other contributes `h_face = 50 m`
      !!     to the FULL sum and nothing to the open one.)
      !!   * `φ_u` and `φ_v` must DIFFER, because the residual is
      !!     `Δa_u = +f·(1−φ_v)·v̄` and `Δa_v = −f·(1−φ_u)·ū`; a spatially
      !!     uniform `φ` retunes the inertial frequency instead of
      !!     growing.  A DIAGONAL ledge gives the u-mask and the v-mask
      !!     different column pairs all along it.
      !!
      !! A real stepped bed was tried first and REJECTED, measured: a
      !! terraced `z_fixed` bathymetry carries an f-dependent,
      !! `pred_corr`-only exponential mode of its own (e-folding ~430
      !! outer steps at this `f` and `dt`), present in IDENTICAL strength
      !! with and without the fix under test — 115x vs 123x over 2300
      !! steps — so it swamps this signal and gates nothing.  Keeping the
      !! bed FLAT removes it (ratio 0.97 at `f = 0`, and see the measured
      !! table below) and leaves the mask as the only difference from
      !! `test_ocean_cor_ref_seiche`, which is exactly what this file is
      !! for.  That mode is a separate finding, not this one.

   real(wp), parameter :: F0 = 2.0e-4_wp
      !! f-plane Coriolis (1/s).  Deliberately stronger than Earth's, for
      !! the same reason `test_ocean_cor_ref_seiche` does it: the residual
      !! is `∝ f` per substep and compounds, so a large `f` makes the
      !! signal decisive in a few thousand cheap steps instead of a
      !! simulated month.
   real(wp), parameter :: DT = 300.0_wp
      !! Outer step (s).  Gravity CFL: `c = √(g·H0) = 62.6 m/s`,
      !! `dt_inner = 300/20 = 15 s`, `c·dt_inner·√2/dx = 0.66`.
   real(wp), parameter :: V0 = 1.0e-3_wp
      !! Seed jet amplitude (m/s).  Small on purpose: the residual is a
      !! LINEAR amplifier, so its GROWTH FACTOR does not care about the
      !! seed (verified: the growth curve is unchanged from `1e-3` down
      !! to `1e-6`), while the nonlinear grid-scale cascade of an
      !! undamped zero-viscosity jet very much does.
   real(wp), parameter :: JET_CELLS = 8.0_wp
   real(wp), parameter :: PI_L = 3.14159265358979324_wp

   integer, parameter :: N_PERIODS = 22
      !! Seiche periods to integrate, `T = 2L/√(g·H_eff)` on the mean
      !! depth — at least the 20 the gate asks for, with headroom, and
      !! the same clock the sibling staircase-seiche test uses.  The
      !! barotropic mode is what the residual feeds on (it is `∝ v̄`), so
      !! the seiche period is the honest clock.  The inertial period
      !! (`2π/f = 31416 s`) is 13x longer, and running that far lets the
      !! undamped zero-viscosity jet cascade to the grid scale, which
      !! manufactures energy of its own and swamps the signal — measured:
      !! at 20 inertial periods BOTH legs read a ratio of ~31.

   real(wp), parameter :: ENERGY_GROWTH_BAR = 2.0_wp
      !! Bar for `(KE+PE)_end / (KE+PE)_0` over `N_PERIODS` seiche
      !! periods, for the CATASTROPHE guard only.  MEASURED on this case
      !! (gfortran 15.1 Release host build, 150 outer steps): `pred_corr`
      !! 1.01, `ssp_rk2` 0.80 — and 1.01 again with the reference defect
      !! reintroduced, which is exactly why this assertion is documented
      !! as a guard and not as the gate (see the module docstring).  It
      !! must never be widened to make a run pass.

contains

   subroutine collect_ocean_zfixed_cor_ref_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cor_ref_is_open_column_mean", test_cor_ref_open_mean), &
                  new_unittest("bt_coriolis_residual_is_roundoff", test_bt_cor_residual), &
                  new_unittest("rotating_staircase_no_growth", test_rotating_staircase) &
                  ]
   end subroutine collect_ocean_zfixed_cor_ref_tests

   ! ---------------------------------------------------------------
   !  1. The direct statement, on a hand-built column.
   ! ---------------------------------------------------------------

   subroutine test_cor_ref_open_mean(error)
      !! `set_cor_ref_velocity(..., from_u_av = .true.)` — the
      !! `pred_corr` branch — must return the OPEN-column depth mean of
      !! `u_av_layer`, not the full-column one.
      !!
      !! One u-face column, `nz = 4`, layer 1 CLOSED but carrying 40 m of
      !! real water: the MASK closes a face, not the thickness, so a test
      !! whose closed layer was a `1e-4` filler could not see the defect
      !! at all (`φ` would be `1 − 1e-6`).  The closed layer's `u_av` is
      !! left at zero, which is where `mask_layer_velocities` leaves the
      !! prognostic velocity it is the time-mean of.
      !!
      !! Bound is DERIVED, never bit-zero: both sides are cancellations of
      !! `nz` products and the compiler is free to contract `h*u` into an
      !! FMA on one side and not the other.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      integer :: i, j, k, ig, iface, jface, nx_t, ny_t
      real(wp) :: h_face, w, num_o, den_o, num_f, den_f
      real(wp) :: mean_open, mean_full, umax, got, resid, bound
      character(len=420) :: msg
      real(wp), parameter :: SAFETY = 64.0_wp
      real(wp), parameter :: H_CLOSED = 40.0_wp

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=NZ)
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
         end do
      end do
      ! `u_av` on the OPEN layers, zero on the closed one.
      do k = 1, NZ
         do j = 1, ny_t
            do i = 1, nx_t + 1
               ms%u_av_layer(i, j, k) = 0.01_wp*real(k, wp) - 0.02_wp
            end do
         end do
      end do
      ms%u_av_layer(:, :, 1) = 0.0_wp
      ms%v_av_layer = 0.0_wp

      num_o = 0.0_wp; den_o = 0.0_wp
      num_f = 0.0_wp; den_f = 0.0_wp
      umax = 0.0_wp
      do k = 1, NZ
         h_face = 0.5_wp*(ms%h_layer(iface - 1, jface, k) + ms%h_layer(iface, jface, k))
         w = h_face*metrics%open_u(iface, jface, k)
         num_o = num_o + w*ms%u_av_layer(iface, jface, k)
         den_o = den_o + w
         num_f = num_f + h_face*ms%u_av_layer(iface, jface, k)
         den_f = den_f + h_face
         umax = max(umax, abs(ms%u_av_layer(iface, jface, k)))
      end do
      mean_open = num_o/den_o
      mean_full = num_f/den_f

      ! ---- anti-vacuity: the two means must actually differ ----
      write (msg, '("the hand-built column has open mean ",es20.13," and full-column ", &
            &"mean ",es20.13,". If these agree the column is effectively all-open ", &
            &"and every assertion below is vacuous -- fix the COLUMN, not the bar.")') &
         mean_open, mean_full
      call check(error, abs(mean_open - mean_full) > 1.0e-6_wp, trim(msg))
      if (allocated(error)) then
         call teardown(ms, metrics)
         return
      end if

      call set_cor_ref_velocity(grid, bt_work, ms, .true., metrics)
      got = bt_work%cor_ref_u(iface, jface)

      resid = abs(got - mean_open)
      bound = SAFETY*real(NZ, wp)*epsilon(1.0_wp)*max(abs(mean_open), umax)
      write (msg, '("cor_ref_u = ",es20.13," but the OPEN-column mean of u_av is ", &
            &es20.13," (residual ",es11.3," vs derived bound ",es11.3,"). The ", &
            &"full-column mean here is ",es20.13," -- if cor_ref_u matched THAT, ", &
            &"set_cor_ref_velocity is taking the full-column depth mean and the ", &
            &"uncancelled f*(1-phi)*vbar forces every barotropic substep.")') &
         got, mean_open, resid, bound, mean_full
      call check(error, resid <= bound, trim(msg))

      call teardown(ms, metrics)
   end subroutine test_cor_ref_open_mean

   subroutine teardown(ms, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
   end subroutine teardown

   ! ---------------------------------------------------------------
   !  2. The MECHANISM, measured in its own units (m/s^2).
   ! ---------------------------------------------------------------

   subroutine test_bt_cor_residual(error)
      !! The spurious barotropic acceleration itself, on a state where
      !! the correct answer is EXACTLY zero.
      !!
      !! Take a column whose layer velocities are a pure barotropic
      !! `(U0, V0)` on the OPEN layers and zero on the closed ones — i.e.
      !! exactly what `mask_layer_velocities` leaves, and exactly a state
      !! whose open-column depth mean is `(U0, V0)`.  Put the same field
      !! in `u_av/v_av` (the `pred_corr` reference source) and in
      !! `u_face_x/v_face_y` (from which `derive_bt_from_layers` builds
      !! `bt_ubt/bt_vbt`).  Then the two branches of
      !! `set_cor_ref_velocity` MUST agree: `pred_corr`'s depth mean of
      !! `u_av` and `ssp_rk2`'s plain copy of `bt_ubt` are the same
      !! number, `(U0, V0)`, and so `subtract_fast_cor_ref` must remove
      !! the SAME forcing either way.
      !!
      !! Running it both ways and differencing gives
      !! `max|ΔF_bt_u|` = the spurious per-substep acceleration the
      !! defect injects, in m/s² — which is `f·(1−φ_v)·V0` by the
      !! derivation in the module docstring.  This is the amplifier
      !! measured directly, and unlike an energy ratio it cannot be
      !! swamped by anything else in the chain.
      !!
      !! Bound is DERIVED from the magnitude of the terms being
      !! differenced, never bit-zero: `subtract_fast_cor_ref` evaluates
      !! `(ζ+f)·v̄ − ∇KE` twice over the same stencil and the compiler may
      !! contract differently in the two passes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      integer :: i, j, k, ig, nx_t, ny_t, i0, i1, j0, j1
      integer :: n_closed
      real(wp) :: dmax, bound, expect
      real(wp), allocatable :: f_corner(:, :), fu_pc(:, :), fv_pc(:, :)
      logical, allocatable :: shelf(:, :)
      character(len=460) :: msg
      real(wp), parameter :: SAFETY = 256.0_wp
      real(wp), parameter :: U0 = 3.0e-3_wp
      real(wp), parameter :: VV0 = 5.0e-3_wp
      real(wp), parameter :: PHI = 0.75_wp
         !! `Σ h·open / Σ h` on a ledge face: 3 open layers of 100 m out
         !! of 4.  Quoted in the failure message so the measured residual
         !! can be read against `f·(1−φ)·V0` on sight.

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      i0 = ig + 1; i1 = ig + NXP
      j0 = ig + 1; j1 = ig + NYP
      ms%nz_ml = NZ
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=NZ)
      call make_cartesian_metrics(metrics, grid)
      call metrics_closed_faces_alloc(metrics, grid, NZ)

      allocate (shelf(nx_t, ny_t))
      allocate (f_corner(nx_t + 1, ny_t + 1), source=F0)
      allocate (fu_pc(size(bt_work%F_bt_u_fast, 1), size(bt_work%F_bt_u_fast, 2)))
      allocate (fv_pc(size(bt_work%F_bt_v_fast, 1), size(bt_work%F_bt_v_fast, 2)))

      call build_ledge_mask(metrics, shelf, nx_t, ny_t, ig, n_closed)

      do k = 1, NZ
         do j = 1, ny_t
            do i = 1, nx_t
               ms%h_layer(i, j, k) = H_NOM
            end do
         end do
      end do
      ! A PURE barotropic state on the open layers.
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_t
            do i = 1, nx_t + 1
               ms%u_face_x_layer(i, j, k) = U0*metrics%open_u(i, j, k)
            end do
         end do
         do j = 1, ny_t + 1
            do i = 1, nx_t
               ms%v_face_y_layer(i, j, k) = VV0*metrics%open_v(i, j, k)
            end do
         end do
      end do
      ms%u_av_layer = ms%u_face_x_layer
      ms%v_av_layer = ms%v_face_y_layer

      call derive_bt_from_layers(grid, bt_work, ms, metrics)

      ! ---- pred_corr: reference from the depth mean of u_av ----
      bt_work%F_bt_u_fast = 0.0_wp
      bt_work%F_bt_v_fast = 0.0_wp
      call set_cor_ref_velocity(grid, bt_work, ms, .true., metrics)
      call subtract_fast_cor_ref(grid, metrics, bt_work, f_corner, &
                                 OBC_WALL, OBC_WALL, OBC_WALL, OBC_WALL, &
                                 .true., .true., .true., .true.)
      fu_pc = bt_work%F_bt_u_fast
      fv_pc = bt_work%F_bt_v_fast

      ! ---- ssp_rk2: reference is a plain copy of bt_ubt ----
      bt_work%F_bt_u_fast = 0.0_wp
      bt_work%F_bt_v_fast = 0.0_wp
      call set_cor_ref_velocity(grid, bt_work, ms, .false., metrics)
      call subtract_fast_cor_ref(grid, metrics, bt_work, f_corner, &
                                 OBC_WALL, OBC_WALL, OBC_WALL, OBC_WALL, &
                                 .true., .true., .true., .true.)

      dmax = 0.0_wp
      do j = j0, j1
         do i = i0, i1 + 1
            dmax = max(dmax, abs(fu_pc(i, j) - bt_work%F_bt_u_fast(i, j)))
         end do
      end do
      do j = j0, j1 + 1
         do i = i0, i1
            dmax = max(dmax, abs(fv_pc(i, j) - bt_work%F_bt_v_fast(i, j)))
         end do
      end do

      write (msg, '("the ledge closed ",I0," interior u-face entries. Zero means ", &
            &"phi == 1 everywhere and this assertion is vacuous -- fix the MASK, ", &
            &"not the bound.")') n_closed
      call check(error, n_closed > 0, trim(msg))
      if (allocated(error)) then
         call teardown(ms, metrics)
         deallocate (shelf, f_corner, fu_pc, fv_pc)
         return
      end if

      expect = F0*(1.0_wp - PHI)*max(U0, VV0)
      bound = SAFETY*epsilon(1.0_wp)*F0*max(U0, VV0)
      write (msg, '("the pred_corr and ssp_rk2 Coriolis references differ by ", &
            &es11.3," m/s^2 on an exactly barotropic open-column state, where ", &
            &"both must be (U0,V0) and the difference must be round-off (derived ", &
            &"bound ",es11.3,"). That difference IS the spurious acceleration the ", &
            &"barotropic substep integrates every inner step; f*(1-phi)*V for this ", &
            &"ledge is ",es11.3," m/s^2. Suspect set_cor_ref_velocity taking the ", &
            &"FULL-column depth mean of u_av instead of the OPEN-column one.")') &
         dmax, bound, expect
      call check(error, dmax <= bound, trim(msg))

      call teardown(ms, metrics)
      deallocate (shelf, f_corner, fu_pc, fv_pc)
   end subroutine test_bt_cor_residual

   pure subroutine build_ledge_mask(metrics, shelf, nx_t, ny_t, ig, n_closed)
      !! The DIAGONAL bed-layer ledge, shared by both integration-shaped
      !! cases.  Built with the same "closed if the layer is dead on
      !! EITHER side" rule `ocean_vcoord_closed_face_masks` uses, so the
      !! u-mask and the v-mask are taken from DIFFERENT column pairs and
      !! `φ_u /= φ_v` all along the diagonal — which is the asymmetry the
      !! residual needs (a uniform `φ` would only retune the inertial
      !! frequency).
      integer, intent(in) :: nx_t, ny_t, ig
      type(ocean_metrics_t), intent(inout) :: metrics
      logical, intent(out) :: shelf(nx_t, ny_t)
      integer, intent(out) :: n_closed
      integer :: i, j, k

      do j = 1, ny_t
         do i = 1, nx_t
            shelf(i, j) = (real(min(max(i - ig, 1), NXP), wp)/real(NXP, wp) &
                           + real(min(max(j - ig, 1), NYP), wp)/real(NYP, wp)) > 1.0_wp
         end do
      end do
      metrics%open_u = 1.0_wp
      metrics%open_v = 1.0_wp
      do j = 1, ny_t
         do i = 2, nx_t
            if (shelf(i - 1, j) .or. shelf(i, j)) metrics%open_u(i, j, 1) = 0.0_wp
         end do
      end do
      do j = 2, ny_t
         do i = 1, nx_t
            if (shelf(i, j - 1) .or. shelf(i, j)) metrics%open_v(i, j, 1) = 0.0_wp
         end do
      end do
      metrics%use_closed_faces = .true.

      n_closed = 0
      do k = 1, NZ
         do j = ig + 1, ig + NYP
            do i = ig + 1, ig + NXP + 1
               if (metrics%open_u(i, j, k) == 0.0_wp) n_closed = n_closed + 1
            end do
         end do
      end do
   end subroutine build_ledge_mask

   ! ---------------------------------------------------------------
   !  2. The energetic statement, over a real run.
   ! ---------------------------------------------------------------

   subroutine run_basin(split_scheme, finite, energy_ratio, n_closed)
      !! Integrate the rotating staircase basin and report
      !! `(KE+PE)_end / (KE+PE)_0` plus the anti-vacuity mask census.
      integer, intent(in) :: split_scheme
      logical, intent(out) :: finite
      real(wp), intent(out) :: energy_ratio
      integer, intent(out) :: n_closed
         !! Number of CLOSED interior u-face entries in the mask.

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

      integer :: i, j, k, ig, i0, i1, j0, j1, step, n_steps, nx_t, ny_t
      real(wp) :: energy0, energy1, x_f, vjet, period_seiche, h_eff
      real(wp), allocatable :: tgt(:, :, :), tot_h(:, :), eta0f(:, :), z_top(:, :)
      logical, allocatable :: shelf(:, :)

      ig = NGHOST
      call grid%init(NXP, NYP, NGHOST, DX, DX)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = F0
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

      ! ---- flat bed, z_fixed target, and a DIAGONAL bed-layer ledge ----
      allocate (tot_h(nx_t, ny_t), source=H0)
      allocate (eta0f(nx_t, ny_t), source=0.0_wp)
      allocate (z_top(nx_t, ny_t), source=0.0_wp)
      allocate (tgt(nx_t, ny_t, NZ), source=0.0_wp)
      allocate (shelf(nx_t, ny_t), source=.false.)
      h_eff = H0
      call ocean_vcoord_z_fixed_target(tgt, tot_h, eta0f, z_top, &
                                       nx_t, ny_t, NZ, H_NOM, H_MIN)

      call metrics_closed_faces_alloc(metrics, grid, NZ)
      call build_ledge_mask(metrics, shelf, nx_t, ny_t, ig, n_closed)

      ! ---- seed: a depth-uniform sinusoidal v-jet on the OPEN layers ----
      ! Masked at seed time so the initial state is the one
      ! `mask_layer_velocities` would leave: a closed face is a z-level
      ! WALL and carries exactly zero normal velocity.
      do j = 1, ny_t
         do i = 1, nx_t
            do k = 1, NZ
               ms%h_layer(i, j, k) = tgt(i, j, k)
            end do
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
            do k = 1, NZ
               ms%v_face_y_layer(i, j, k) = vjet*metrics%open_v(i, j, k)
            end do
            dyn%bt_work%bt_H_ref(i, j) = tot_h(i, j)
         end do
      end do
      ms%u_face_x_layer(nx_t + 1, :, :) = 0.0_wp
      ms%v_face_y_layer(:, ny_t + 1, :) = 0.0_wp

      period_seiche = 2.0_wp*real(NXP, wp)*DX/sqrt(GRAV*h_eff)
      n_steps = nint(real(N_PERIODS, wp)*period_seiche/DT)

      energy0 = basin_energy(ms, dyn, i0, i1, j0, j1)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()

      do step = 1, n_steps
         ! The engine refreshes the barotropic face widths from the LIVE
         ! `h` once per outer step (`rdb_ocean_engine`), and a test that
         ! drives `ocean_dyn_step_split` directly has to do it too: without
         ! it `dy_cu_bt`/`dx_cv_bt` keep their un-narrowed `dy_cu`/`dx_cv`
         ! seed, the fast loop transports on the FULL column while the
         ! layers transport on the OPEN one, and THAT inconsistency is an
         ! amplifier of its own -- measured at x25 over 2000 steps here,
         ! which is enough to swamp the residual under test.
         call ocean_porous_refresh(grid, metrics, ms)
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc)
      end do

      ! Host read-back of the COMPONENT arrays only -- never the aggregate
      ! derived type (that overwrites the host descriptors with device
      ! addresses and the next host read segfaults).
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

      energy1 = basin_energy(ms, dyn, i0, i1, j0, j1)
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
      deallocate (tgt, tot_h, eta0f, z_top, shelf)
   end subroutine run_basin

   pure function basin_energy(ms, dyn, i0, i1, j0, j1) result(energy)
      !! Depth-integrated `KE + PE` per unit area and per unit density:
      !! `Σ_k h·(u² + v²)/2 + g·η²/2`, with `η = Σ_k h − bt_H_ref`.  A
      !! monotone energy DIAGNOSTIC, not a discretely conserved energy,
      !! which is all a non-growth assertion needs.
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

   subroutine test_rotating_staircase(error)
      !! The CATASTROPHE guard, not the discriminator -- read the module
      !! docstring before touching this bar.  `ssp_rk2` control first (it
      !! takes the copy-of-`bt_ubt` branch and so cannot have the defect),
      !! then `pred_corr` on the identical configuration.
      !!
      !! What it does gate: that a closed-face rotating basin still
      !! completes, stays finite, and does not blow up over
      !! `N_PERIODS` seiche periods.  What it does NOT gate is the
      !! Coriolis reference -- measured, it reads the same ratio with and
      !! without the reference defect, because a SECOND `pred_corr` x
      !! closed-faces x rotation amplifier of comparable strength lives on
      !! this path and neither test nor fix addresses it yet.
      type(error_type), allocatable, intent(out) :: error
      logical :: fin_pc, fin_rk2
      real(wp) :: ratio_pc, ratio_rk2
      integer :: n_closed_pc, n_closed_rk2
      character(len=440) :: msg

      call run_basin(SPLIT_SCHEME_SSP_RK2, fin_rk2, ratio_rk2, n_closed_rk2)

      write (msg, '("the staircase closed ",I0," interior u-face entries. Zero ", &
            &"means the geometry stopped producing a staircase, phi == 1 and ", &
            &"every assertion below is vacuous -- fix the GEOMETRY, not the bar.")') &
         n_closed_rk2
      call check(error, n_closed_rk2 > 0, trim(msg))
      if (allocated(error)) return

      call check(error, fin_rk2, "ssp_rk2 rotating staircase: the control went non-finite")
      if (allocated(error)) return
      write (msg, '("ssp_rk2 rotating staircase control: KE+PE ratio = ",es11.3, &
            &" (the control must itself not grow, else the case is unusable)")') ratio_rk2
      call check(error, ratio_rk2 < ENERGY_GROWTH_BAR, trim(msg))
      if (allocated(error)) return

      call run_basin(SPLIT_SCHEME_PRED_CORR, fin_pc, ratio_pc, n_closed_pc)
      call check(error, fin_pc, &
                 "pred_corr rotating staircase: went non-finite (the mode blew up)")
      if (allocated(error)) return

      write (msg, '("pred_corr rotating staircase: KE+PE ratio = ",es11.3," vs ssp_rk2 ", &
            &es11.3," over ",I0," seiche periods. Closed basin, unforced and ", &
            &"undamped, so growth is MANUFACTURED energy. This is the CATASTROPHE ", &
            &"guard, not the Coriolis-reference gate -- that is ", &
            &"cor_ref_is_open_column_mean and bt_coriolis_residual_is_roundoff. Do ", &
            &"NOT widen this bar; find what started manufacturing energy.")') &
         ratio_pc, ratio_rk2, N_PERIODS
      call check(error, ratio_pc < ENERGY_GROWTH_BAR, trim(msg))
   end subroutine test_rotating_staircase

end module test_ocean_zfixed_cor_ref
