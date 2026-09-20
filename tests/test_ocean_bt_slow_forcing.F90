!! THE `F_slow` MEMBERSHIP CONTRACT — a layer velocity tendency applied
!! inside the split solver's corrected set reaches the DEPTH MEAN exactly
!! once, whether or not its depth mean was summed into `F_slow`.
!!
!! ## Why this test exists
!!
!! `sum_slow_tendencies_into_F_slow` builds the frozen barotropic forcing
!! from a HARD-CODED five-term list (PGF, Coriolis-advection, horizontal
!! viscosity, bottom drag, surface stress).  The side-wall CHANNEL drag
!! (`ocean_channel_drag_apply_tendencies`) is applied to the layers one
!! line after the bottom drag and is NOT in that list, which invites the
!! reading that its depth-mean effect is "lost" or enters the barotropic
!! mode one stage late.
!!
!! It is neither, and the reason is structural: `apply_bt_correction`
!! adds an INCREMENT
!!
!!     Δu = u_bt^end − u_bt^n − dt·F_bt
!!
!! to every layer — it does not REPLACE the layer depth mean with
!! `u_bt^end`.  Write `T_k` for the tendencies that are in `F_slow` and
!! `D_k` for one that is applied but omitted.  The stage does
!!
!!     u_k^{n+1} = u_k^n + dt·(T_k + D_k) + (u_bt^end − u_bt^n − dt·F_bt)
!!
!! with `F_bt = ⟨T⟩_h` and `⟨u^n⟩_h = u_bt^n` (`derive_bt_from_layers`).
!! Taking the thickness-weighted depth mean:
!!
!!     ⟨u^{n+1}⟩ = u_bt^n + dt·⟨T⟩ + dt·⟨D⟩ + u_bt^end − u_bt^n − dt·F_bt
!!               = u_bt^end + dt·⟨D⟩.
!!
!! So `⟨D⟩` is applied EXACTLY ONCE, on top of the barotropic solution —
!! not lost, not double counted.  The mirror algebra for an in-`F_slow`
!! term is just as tight: the fast loop integrates `⟨T⟩` over the
!! substeps and the `−dt·F_bt` guard takes it straight back out, so a
!! term's `F_slow` membership is depth-mean NEUTRAL to leading order.
!! What membership actually buys is second order: the substep's live
!! η / ζ / KE / `bt_rem` trajectory, and the time-mean transports
!! `bt_uhbt` handed to continuity, see the forcing.  Omitting a term is a
!! first-order-in-dt operator SPLIT, not a missing term.
!!
!! ## The two arms
!!
!! A doubly-periodic, flat-bottom, uniform-density, non-rotating box
!! carrying a uniform barotropic current, with every closure off but one
!! drag.  Uniformity makes `∇·(hu) ≡ 0`, `η ≡ 0`, `ζ ≡ 0`, `∇KE ≡ 0`, so
!! the PGF, the Coriolis-advection and the whole barotropic response are
!! identically zero and the ONE active drag owns the entire answer.
!!
!! * `channel_drag_depth_mean_decay` — channel drag, which is NOT in
!!   `F_slow`.  `wet_q ≡ 0` makes the blocked perimeter fraction 1 at
!!   every face, so `λ = cdrag_side·|u|/dyCu` is uniform and the implicit
!!   apply `u ← u/(1 + dt·λ)` has a closed form through the outer scheme.
!! * `bottom_drag_depth_mean_decay` — bed-only linear bottom drag, which
!!   IS in `F_slow`, on the identical box.  Its depth mean is `−r·u/nz`
!!   and only `k=1` moves.
!!
!! Both arms are checked against the outer scheme's exact scalar
!! recursion under BOTH split schemes.  If an omitted tendency's depth
!! mean were dropped, arm 1 would not decay at all (every layer is
!! dragged identically, so the whole signal is barotropic); if it were
!! double counted it would decay at twice the rate.  The measured value
!! is neither — it is the analytic one to round-off.
!!
!! **Contract for a new tendency**: read the `sum_slow_tendencies_into_F_slow`
!! docstring before adding one.  Summing it into `F_slow` is the default
!! (it keeps the fast loop's η and transports consistent with the
!! forcing); omitting it is legal and costs one operator-split lag, but
!! never a lost or duplicated term.
module test_ocean_bt_slow_forcing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_LITE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, BDRAG_LINEAR
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, &
                                       ocean_bc_state_enter_data, ocean_bc_state_exit_data, &
                                       ocean_bc_state_set_topology
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, &
                            SPLIT_SCHEME_PRED_CORR, SPLIT_SCHEME_SSP_RK2
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_LAGRANGIAN
   implicit none
   private

   public :: collect_ocean_bt_slow_forcing_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NXP = 8
   integer, parameter :: NYP = 8
   integer, parameter :: NZ = 4
   integer, parameter :: N_INNER = 24
   integer, parameter :: N_STEPS = 40
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: H0 = 400.0_wp        !! flat total depth (m)
   real(wp), parameter :: U0 = 1.0_wp          !! uniform barotropic current (m/s)
   real(wp), parameter :: DT = 200.0_wp        !! outer step (s)

   real(wp), parameter :: CDRAG_SIDE = 0.4_wp
      !! Deliberately far above a physical side-drag coefficient: the
      !! per-stage decay `dt·cdrag_side·U0/dy = 0.08` has to be large
      !! enough that "decayed at the analytic rate", "never decayed" and
      !! "decayed twice" are decades apart after `N_STEPS`.
   real(wp), parameter :: R_LINEAR = 4.0e-4_wp
      !! Linear bed-drag rate (1/s); `dt·r = 0.08`, matching arm 1's
      !! per-stage strength so the two arms are directly comparable.

   real(wp), parameter :: TOL_REL = 1.0e-11_wp
      !! Relative agreement demanded against the closed form.  The
      !! configuration makes every other operator identically zero, so
      !! the only slack is floating-point accumulation over `N_STEPS`
      !! outer steps and `N_INNER` barotropic substeps each.

contains

   subroutine collect_ocean_bt_slow_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("channel_drag_depth_mean_decay", test_channel_drag_decay), &
                  new_unittest("bottom_drag_depth_mean_decay", test_bottom_drag_decay) &
                  ]
   end subroutine collect_ocean_bt_slow_forcing_tests

   subroutine run_box(split_scheme, side_drag, u_layer)
      !! Integrate `N_STEPS` outer steps of the doubly-periodic uniform
      !! box and report the per-layer velocity at an interior u-face.
      !! `side_drag = .true.` selects the channel (side-wall) drag — the
      !! tendency that is NOT in `F_slow`; `.false.` selects the bed-only
      !! linear bottom drag, which is.  Exactly one is ever active.
      integer, intent(in) :: split_scheme
      logical, intent(in) :: side_drag
      real(wp), intent(out) :: u_layer(NZ)

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
      type(ocean_bc_state_t) :: bc

      integer :: i, j, k, step, i_probe, j_probe, ierr

      call grid%init(NXP, NYP, NGHOST, DX, DX)
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
      vc%coord_type = VCOORD_LAGRANGIAN
      dyn%split_scheme = split_scheme

      ! Everything off but the one drag under test.
      dyn%enable_thermodynamics = .false.
      vmix%use_closure = .false.
      vmix%use_kpp = .false.
      vd%K_v_tracer = 0.0_wp
      vd%K_v_momentum = 0.0_wp
      hv%nu_h = 0.0_wp
      hd%kappa_h = 0.0_wp
      bd%c_drag = 0.0_wp
      bd%r_linear = 0.0_wp
      bd%hbbl = 0.0_wp
      bd%variant = BDRAG_LINEAR
      bd%channel_drag = .false.
      bd%cdrag_side = 0.0_wp
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
      if (side_drag) then
         bd%channel_drag = .true.
         bd%cdrag_side = CDRAG_SIDE
      else
         bd%r_linear = R_LINEAR
      end if

      ! Uniform flat box: every cell (ghosts included) identical, so the
      ! divergence, the density gradient, the vorticity and the KE
      ! gradient are all identically zero.
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            do k = 1, NZ
               ms%h_layer(i, j, k) = H0/real(NZ, wp)
               if (allocated(ms%tracers)) then
                  if (ms%idx_salinity > 0) &
                     ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
                  if (ms%idx_temperature > 0) &
                     ms%tracers(ms%idx_temperature)%hTr(i, j, k) = eos%T_ref*ms%h_layer(i, j, k)
               end if
            end do
            ms%rho_layer(i, j, :) = eos%rho0
            dyn%bt_work%bt_H_ref(i, j) = H0
         end do
      end do
      ms%u_face_x_layer(:, :, :) = U0
      ms%v_face_y_layer(:, :, :) = 0.0_wp

      call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=size(ms%tracers))
      call ocean_bc_state_set_topology(bc, .true., .true., ierr)

      call make_cartesian_metrics(metrics, grid)
      ! Blocked side perimeter everywhere: `wet_q = 0` drives the channel
      ! drag's `f_blocked` to 1 at every face, which makes λ uniform and
      ! the decay analytic.  `wet_u` / `wet_v` stay 1, so the land-face
      ! masks (`mask_bt_rem`, `mask_layer_velocities`) remain no-ops, and
      ! with `f_0 = 0` and a uniform flow the only other `wet_q` consumers
      ! (Coriolis-advection's corner vorticity, the Leith/Smagorinsky
      ! strain) are multiplying zero.
      metrics%wet_q(:, :) = 0.0_wp
      !$acc update device(metrics%wet_q)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ct%enter_data(); call cor%enter_data(); call pgf%enter_data()
      call hv%enter_data(); call bd%enter_data(); call ss%enter_data()
      call va%enter_data(); call hd%enter_data()
      call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()
      call vc%enter_data()
      call ocean_bc_state_enter_data(bc)

      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, vcoord=vc, bc=bc)
      end do

      !$acc update self(ms%u_face_x_layer)
      i_probe = NGHOST + NXP/2
      j_probe = NGHOST + NYP/2
      do k = 1, NZ
         u_layer(k) = ms%u_face_x_layer(i_probe, j_probe, k)
      end do

      call ocean_bc_state_exit_data(bc)
      call vc%exit_data()
      call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
      call hd%exit_data(); call va%exit_data()
      call ss%exit_data(); call bd%exit_data(); call hv%exit_data()
      call pgf%exit_data(); call cor%exit_data(); call ct%exit_data()
      !$acc exit data delete(ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
      call ms%destroy()
   end subroutine run_box

   pure function channel_drag_reference(split_scheme) result(u_end)
      !! Closed form for arm 1, carried through the outer scheme.
      !! Per stage the side drag is the frozen-rate backward-Euler
      !! `u ← u/(1 + dt_vel·λ)` with `λ = cdrag_side·|u_stage_entry|/dy`
      !! and `f_blocked = 1`.  Every other operator contributes zero, and
      !! `apply_bt_correction`'s Δ is `u_bt^end − u_bt^n − dt·F_bt = 0`
      !! because the substep's only forcing IS `F_bt` and η never moves.
      !!
      !!   `ssp_rk2`   — two identical stages at `dt`, then the SSP
      !!                 average `u^{n+1} = ½(u^n + u^{(2)})`.
      !!   `pred_corr` — the predictor's provisional velocity is
      !!                 discarded by `restore_state`, and the corrector
      !!                 re-derives λ from the restored `u^n`, so the
      !!                 whole step is ONE `dt` application.
      integer, intent(in) :: split_scheme
      real(wp) :: u_end
      real(wp) :: u, u1, u2, c_side
      integer :: step

      c_side = CDRAG_SIDE/DX
      u = U0
      do step = 1, N_STEPS
         if (split_scheme == SPLIT_SCHEME_PRED_CORR) then
            u = u/(1.0_wp + DT*c_side*abs(u))
         else
            u1 = u/(1.0_wp + DT*c_side*abs(u))
            u2 = u1/(1.0_wp + DT*c_side*abs(u1))
            u = 0.5_wp*(u + u2)
         end if
      end do
      u_end = u
   end function channel_drag_reference

   pure subroutine bottom_drag_reference(split_scheme, u_bed, u_int)
      !! Closed form for arm 2.  The bed-only linear tendency is explicit,
      !! `du_drag(:,:,1) = −r·u_1`, and is zero for `k ≥ 2`.  Its depth
      !! mean `F_bt = −r·u_1/nz` is what the fast loop integrates and what
      !! `apply_bt_correction` subtracts back out, so Δ is again zero and
      !! only `k = 1` moves.
      integer, intent(in) :: split_scheme
      real(wp), intent(out) :: u_bed, u_int
      real(wp) :: b, b1, b2
      integer :: step

      b = U0
      do step = 1, N_STEPS
         if (split_scheme == SPLIT_SCHEME_PRED_CORR) then
            b = b*(1.0_wp - DT*R_LINEAR)
         else
            b1 = b*(1.0_wp - DT*R_LINEAR)
            b2 = b1*(1.0_wp - DT*R_LINEAR)
            b = 0.5_wp*(b + b2)
         end if
      end do
      u_bed = b
      ! Interior layers see no tendency at all, and the SSP average of an
      ! unchanged field is the field, under either scheme.
      u_int = U0
   end subroutine bottom_drag_reference

   subroutine test_channel_drag_decay(error)
      !! Arm 1 — the tendency that is NOT in `F_slow`.  Every layer is
      !! dragged identically, so the ENTIRE signal is barotropic: a
      !! dropped depth mean shows up as no decay at all, a double count as
      !! decay at twice the rate.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_pc(NZ), u_rk2(NZ), ref_pc, ref_rk2
      integer :: k
      character(len=360) :: msg

      ref_pc = channel_drag_reference(SPLIT_SCHEME_PRED_CORR)
      ref_rk2 = channel_drag_reference(SPLIT_SCHEME_SSP_RK2)

      call run_box(SPLIT_SCHEME_PRED_CORR, .true., u_pc)
      call run_box(SPLIT_SCHEME_SSP_RK2, .true., u_rk2)

      write (msg, '("channel drag (NOT in F_slow), pred_corr: u = ", es22.15, &
            &" vs analytic ", es22.15, " (undecayed ", es12.5, ").  A wrong value here ", &
            &"means the omitted tendency no longer reaches the depth mean exactly once ", &
            &"— see sum_slow_tendencies_into_F_slow.")') u_pc(1), ref_pc, U0
      call check(error, abs(u_pc(1) - ref_pc) <= TOL_REL*abs(ref_pc), trim(msg))
      if (allocated(error)) return

      write (msg, '("channel drag (NOT in F_slow), ssp_rk2: u = ", es22.15, &
            &" vs analytic ", es22.15, " (undecayed ", es12.5, ").")') &
         u_rk2(1), ref_rk2, U0
      call check(error, abs(u_rk2(1) - ref_rk2) <= TOL_REL*abs(ref_rk2), trim(msg))
      if (allocated(error)) return

      ! Uniform drag ⇒ no baroclinic structure may appear.  The bt
      ! correction distributes Δ uniformly, so any spread between layers
      ! would be the barotropic bookkeeping leaking into the layers.
      do k = 2, NZ
         write (msg, '("channel drag: layer ", i0, " velocity ", es22.15, &
               &" differs from layer 1 ", es22.15, " — a uniformly applied drag must ", &
               &"stay depth-uniform.")') k, u_pc(k), u_pc(1)
         call check(error, abs(u_pc(k) - u_pc(1)) <= TOL_REL*abs(u_pc(1)), trim(msg))
         if (allocated(error)) return
         write (msg, '("channel drag (ssp_rk2): layer ", i0, " velocity ", es22.15, &
               &" differs from layer 1 ", es22.15, ".")') k, u_rk2(k), u_rk2(1)
         call check(error, abs(u_rk2(k) - u_rk2(1)) <= TOL_REL*abs(u_rk2(1)), trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_channel_drag_decay

   subroutine test_bottom_drag_decay(error)
      !! Arm 2 — the mirror, with a tendency that IS in `F_slow`.  It too
      !! reaches the depth mean exactly once: the fast loop integrates
      !! `⟨D⟩` and `−dt·F_bt` takes it straight back out, so the bed layer
      !! carries the whole decay and the interior layers do not move.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_pc(NZ), u_rk2(NZ), bed_pc, int_pc, bed_rk2, int_rk2
      integer :: k
      character(len=360) :: msg

      call bottom_drag_reference(SPLIT_SCHEME_PRED_CORR, bed_pc, int_pc)
      call bottom_drag_reference(SPLIT_SCHEME_SSP_RK2, bed_rk2, int_rk2)

      call run_box(SPLIT_SCHEME_PRED_CORR, .false., u_pc)
      call run_box(SPLIT_SCHEME_SSP_RK2, .false., u_rk2)

      write (msg, '("bottom drag (IN F_slow), pred_corr: bed u = ", es22.15, &
            &" vs analytic ", es22.15, ".  Too much decay means the depth mean is ", &
            &"counted twice (once by the layer apply, once left in via F_bt).")') &
         u_pc(1), bed_pc
      call check(error, abs(u_pc(1) - bed_pc) <= TOL_REL*abs(bed_pc), trim(msg))
      if (allocated(error)) return

      write (msg, '("bottom drag (IN F_slow), ssp_rk2: bed u = ", es22.15, &
            &" vs analytic ", es22.15, ".")') u_rk2(1), bed_rk2
      call check(error, abs(u_rk2(1) - bed_rk2) <= TOL_REL*abs(bed_rk2), trim(msg))
      if (allocated(error)) return

      do k = 2, NZ
         write (msg, '("bottom drag: interior layer ", i0, " moved to ", es22.15, &
               &" from ", es22.15, " — a bed-only tendency must leave the interior ", &
               &"layers untouched once its depth mean is removed by apply_bt_correction.")') &
            k, u_pc(k), int_pc
         call check(error, abs(u_pc(k) - int_pc) <= TOL_REL*abs(int_pc), trim(msg))
         if (allocated(error)) return
         write (msg, '("bottom drag (ssp_rk2): interior layer ", i0, " moved to ", es22.15, &
               &" from ", es22.15, ".")') k, u_rk2(k), int_rk2
         call check(error, abs(u_rk2(k) - int_rk2) <= TOL_REL*abs(int_rk2), trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_bottom_drag_decay

end module test_ocean_bt_slow_forcing
