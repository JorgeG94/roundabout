!! O2 acceptance gate: 2-rank vs 1-rank trajectory agreement for the ocean
!! split-RK2 dyn-core.
!!
!! Proves the O2 multi-rank wiring by running the same small double-gyre-like
!! problem in two ways within ONE binary, dispatched by MPI rank count.
!!
!! The binary is registered by ctest at np=1 AND np=2 (two legs).
!! On EACH launch it:
!!   1. Builds the SERIAL REFERENCE on a full NX_G x NY_G px=1 grid.
!!      (All ranks compute the same reference — deterministic.)
!!   2. Destroys the reference halo topology.
!!   3. Builds the DECOMPOSED run on the local subdomain.
!!      On np=1: subdomain == global, trivially identical.
!!      On np=2: x-split (px=2, py=1); the seam bisects the domain.
!!   4. Asserts:
!!      (a) CONSERVATION:  global mass after N steps == global mass at step 0
!!          to relative < 1e-6.  Checked on BOTH serial reference and decomposed.
!!      (b) FINITENESS: global mass and global KE are finite.
!!      (c) DECOMPOSITION-INVARIANCE: decomposed global mass / global KE
!!          agree with the serial reference to relative < 1e-6 (mass) and
!!          < 1e-4 (KE).  On np=1 this is trivially bit-exact; on np=2
!!          it proves the seam exchange is wired correctly.
!!
!! Grid: NX_G=16, NY_G=8, NZ=2, nghost=3.
!! Physics: f-plane Coriolis, uniform flat bathymetry H0=200 m,
!!          2-gyre wind stress tau_x(j) = TAUX*(1 - cos(2π*(j-0.5)/NY_G))
!!          seeded per-rank via global j index.  No explicit viscosity,
!!          no tracer diffusion (defaults off).  The wind drives real flow
!!          that crosses the x-seam within N_STEPS=100 steps, testing O2 wiring.
!! N_STEPS=100 gives a broken seam exchange time to grow far beyond the KE
!! gate (a fully broken exchange measured 7.5e-3 at 100 steps vs 9.5e-6 at 20;
!! the gate is now 1e-13 against measured round-off ~9e-16 — see AGREE_KE_TOL).
!!
!! nghost = 3 on ALL legs (O4): the continuity/tracer PPM reconstruction is
!! 5-point, so a rank-seam face needs THREE ghost columns for full-order
!! reconstruction of the ghost-side donor cell.  At nghost = 2 the PPM
!! local-array-edge first-order fallback fires AT the seam face and the
!! wall/island legs drift from the serial reference at O(1e-6) relative KE
!! per 100 steps (measured growth 3.4e-9 @ N=1 → 5.7e-8 @ N=10 → 2.0e-7
!! @ N=50 → 8.9e-7 @ N=100 pre-fix) while the periodic leg — which already
!! ran nghost = 3 — was bit-exact.  Post-fix (nghost = 3 + the has_*-gated
!! seam-face flux renormalisation in rdb_continuity) all three legs agree
!! with serial at round-off: wall 5.2e-16, island 8.5e-16, periodic 6.4e-16.
!!
!! Architecture note: ocean_halo_init is a module-level singleton.
!!   We init it for the serial reference (px=1), run+measure, destroy; then
!!   re-init for the decomposed run (actual rank count), run+measure, destroy.
#ifdef RDB_ENABLE_MPI
program test_ocean_dyn_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, GRAVITY, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_finalize, metrics_apply_land_mask
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, ocean_dyn_enable_bt_wide
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, ocean_bc_state_set_edges, &
                                       OBC_PERIODIC, ocean_bc_validate_periodic
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy, ocean_halo_reserve, &
                             ocean_halo_centre, ocean_halo_is_decomposed
   use rdb_ocean_halo_counters, only: oh_counters_reset, oh_counters_get, oh_counters_msgs
   use rdb_halo, only: halo_allreduce_sum
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   ! --- Problem parameters (global) ---
   integer, parameter :: NX_G = 16
      !! Global physical cells in x
   integer, parameter :: NY_G = 8
      !! Global physical cells in y
   integer, parameter :: NZ = 2
      !! Number of vertical layers
   integer, parameter :: NG = 3
      !! Ghost cell width.  3 (not 2) — the continuity/tracer PPM 5-point
      !! reconstruction needs 3 ghost columns for full-order reconstruction
      !! at a rank-seam face (see header).  Serial runs are insensitive
      !! (the array-edge fallback coincides with zero-flux walls there).
   real(wp), parameter :: DX = 50000.0_wp
      !! Cell spacing in x (m)
   real(wp), parameter :: DY = 50000.0_wp
      !! Cell spacing in y (m)
   real(wp), parameter :: H0 = 200.0_wp
      !! Total column depth per layer split evenly
   real(wp), parameter :: F_0 = 1.0e-4_wp
      !! f-plane Coriolis (s^-1)
   real(wp), parameter :: TAUX = 0.1_wp
      !! Wind stress amplitude (N/m^2)
   integer, parameter :: N_STEPS = 100
      !! Outer steps to integrate
   integer, parameter :: N_INNER = 4
      !! Fixed inner substeps (deterministic — no auto_n_inner)
   real(wp), parameter :: DT = 600.0_wp
      !! Outer time step (s)
   real(wp), parameter :: MASS_TOL = 1.0e-12_wp
      !! Conservation tolerance (relative).  Tightened from 1e-6 to 1e-12 in
      !! O3 after the intra-continuity Lie-split seam exchange landed: the
      !! measured decomposed conservation drift is ~3.7e-16 (round-off), so a
      !! 1e-12 gate leaves ~4 orders of margin while being a genuinely tight
      !! band.  Robust: holds bit-for-bit whether or not the Lie-split seam
      !! exchange is present (this flat-bottom all-wall config's column-mass
      !! metric is dominated by the BT-transport renormalisation, which pins
      !! sum_k(h) = H_ref + eta_end to machine precision either way; the seam
      !! exchange corrects PER-LAYER seam transport, verified separately by
      !! test_halo_ocean_mpi's bit-exact primitive round-trip).
   real(wp), parameter :: AGREE_MASS_TOL = 1.0e-12_wp
      !! Decomposition-invariance tolerance for mass (relative).  Tightened
      !! from 1e-6 to 1e-12 in O3; measured decomp-vs-serial mass reldiff is
      !! ~8.5e-16 (round-off).  See MASS_TOL for the robustness rationale.
   real(wp), parameter :: AGREE_KE_TOL = 1.0e-13_wp
      !! Decomposition-invariance tolerance for KE (relative).  Tightened
      !! from 1e-4 to 1e-13 in O4 after the seam-face PPM (nghost = 3) +
      !! has_*-gated flux-renorm fixes landed: measured post-fix KE reldiff
      !! is 5.2e-16 (wall) / 8.5e-16 (island) / 6.4e-16 (periodic) at
      !! N_STEPS = 100 — a 1e-13 gate leaves ~2.5 orders of margin while
      !! catching any seam regression at round-off scale.
   real(wp), parameter :: TRACER_CLOSURE_TOL = 1.0e-12_wp
      !! Tracer CONSERVATION tolerance (relative drift of the global salt /
      !! heat sum vs its own step-0 value).  This is the gate that the
      !! seamount heat-leak postmortem exposed as MISSING: the old test
      !! checked mass + KE but never salt/heat closure.  At the shipped
      !! nghost = 3 default the measured drift on this config is ~1e-15
      !! (round-off) for BOTH the uniform-S and structured-T tracers, so a
      !! 1e-12 gate leaves ~3 orders of margin.  TEETH: reverting the default
      !! to nghost = 2 makes the seam-face PPM fall to first order and the
      !! STRUCTURED-T heat sum drifts + grows (~1e-11 by N_STEPS on this
      !! config; ~5e-11/day on seamount_bench_full) — heat FAILS this gate
      !! while the UNIFORM-S salt stays clean, exactly the salt-vs-heat
      !! discriminator (verified by temporarily forcing NG = 2).
   real(wp), parameter :: TRACER_AGREE_TOL = 1.0e-12_wp
      !! Decomposition-invariance tolerance for the global salt / heat sums
      !! (decomp vs serial, relative).  Same round-off margin as
      !! TRACER_CLOSURE_TOL; catches any future seam regression that moves
      !! tracer transport off the serial trajectory.

   integer :: rank, nprocs, n_fail, total_fail
   type(comm_t) :: comm

   real(wp) :: ref_mass0, ref_massN, ref_ke
   real(wp) :: dec_mass0, dec_massN, dec_ke
   real(wp) :: ref_salt0, ref_saltN, ref_heat0, ref_heatN
   real(wp) :: dec_salt0, dec_saltN, dec_heat0, dec_heatN

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   n_fail = 0

   ! --- Step 1: SERIAL REFERENCE (all ranks, same global grid, px=1) ---
   call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                             ref_salt0, ref_saltN, ref_heat0, ref_heatN)

   ! --- Step 2: DECOMPOSED RUN (rank's local subdomain) ---
   call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                       dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- Assertions ---

   ! (b) Finiteness
   call assert_finite(ref_massN, "serial reference final mass", rank)
   call assert_finite(ref_ke, "serial reference final KE", rank)
   call assert_finite(dec_massN, "decomposed final mass", rank)
   call assert_finite(dec_ke, "decomposed final KE", rank)

   ! (a) Conservation — serial reference
   block
      real(wp) :: drift
      drift = abs(ref_massN - ref_mass0)/max(abs(ref_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL conservation: rank ", rank, &
            " serial ref mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (a) Conservation — decomposed run
   block
      real(wp) :: drift
      drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL conservation: rank ", rank, &
            " decomposed mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance: mass
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4)') "MEASURED mass agreement reldiff=", rel_diff
      if (rel_diff >= AGREE_MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL mass agreement: rank ", rank, &
            " decomp=", dec_massN, " ref=", ref_massN, &
            " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance: KE
   block
      real(wp) :: rel_diff
      ! Avoid divide-by-zero for near-zero KE (trivial IC); use absolute floor.
      rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
         "MEASURED KE agreement reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
      if (rel_diff >= AGREE_KE_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL KE agreement: rank ", rank, &
            " decomp=", dec_ke, " ref=", ref_ke, &
            " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
      end if
   end block

   ! (a)+(c) Salt/heat closure + agreement — wall leg (the gate gap the
   ! seamount heat-leak postmortem exposed).
   call check_tracer_closure_and_agreement("wall", rank, &
                                           ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                           dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- Land (island straddling the seam) leg ---
   ! The KEY assertion is decomposition-invariance with the island:
   ! the seam x-faces between two land cells must be masked consistently
   ! on both the serial (px=1) and decomposed (px=nprocs) runs.  This is
   ! only true if the wet_mask seam exchange in metrics_apply_land_mask
   ! (Part A) is wired correctly.  Labelled "island" for diagnosability.
   call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                             ref_salt0, ref_saltN, ref_heat0, ref_heatN, use_island=.true.)
   call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                       dec_salt0, dec_saltN, dec_heat0, dec_heatN, use_island=.true.)

   ! (b) Finiteness — island leg
   call assert_finite(ref_massN, "island serial reference final mass", rank)
   call assert_finite(ref_ke, "island serial reference final KE", rank)
   call assert_finite(dec_massN, "island decomposed final mass", rank)
   call assert_finite(dec_ke, "island decomposed final KE", rank)

   ! (a) Conservation — island serial reference
   block
      real(wp) :: drift
      drift = abs(ref_massN - ref_mass0)/max(abs(ref_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL island conservation: rank ", rank, &
            " serial ref mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (a) Conservation — island decomposed run
   block
      real(wp) :: drift
      drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL island conservation: rank ", rank, &
            " decomposed mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — island mass
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4)') "MEASURED mass agreement reldiff=", rel_diff
      if (rel_diff >= AGREE_MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL island mass agreement: rank ", rank, &
            " decomp=", dec_massN, " ref=", ref_massN, &
            " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — island KE
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
         "MEASURED KE agreement reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
      if (rel_diff >= AGREE_KE_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL island KE agreement: rank ", rank, &
            " decomp=", dec_ke, " ref=", ref_ke, &
            " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
      end if
   end block

   ! (a)+(c) Salt/heat closure + agreement — island leg
   call check_tracer_closure_and_agreement("island", rank, &
                                           ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                           dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- Periodic-x (Eady channel) leg ---
   ! Exercises the multi-rank periodic wraparound: rank 0's west ghost is
   ! filled from the last rank's east interior THROUGH the halo module (a
   ! rank-to-rank message under decomposition, a local wrap at px=1).  The
   ! wall/island legs never cross the domain-wrap seam.  Decomposition-
   ! invariance here proves the periodic rank link is wired correctly.
   call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                             ref_salt0, ref_saltN, ref_heat0, ref_heatN, use_periodic=.true.)
   call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                       dec_salt0, dec_saltN, dec_heat0, dec_heatN, use_periodic=.true.)

   ! (b) Finiteness — periodic leg
   call assert_finite(ref_massN, "periodic serial reference final mass", rank)
   call assert_finite(ref_ke, "periodic serial reference final KE", rank)
   call assert_finite(dec_massN, "periodic decomposed final mass", rank)
   call assert_finite(dec_ke, "periodic decomposed final KE", rank)

   ! (a) Conservation — periodic serial reference
   block
      real(wp) :: drift
      drift = abs(ref_massN - ref_mass0)/max(abs(ref_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL periodic conservation: rank ", rank, &
            " serial ref mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (a) Conservation — periodic decomposed run
   block
      real(wp) :: drift
      drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL periodic conservation: rank ", rank, &
            " decomposed mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — periodic mass
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4)') "MEASURED mass agreement reldiff=", rel_diff
      if (rel_diff >= AGREE_MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL periodic mass agreement: rank ", rank, &
            " decomp=", dec_massN, " ref=", ref_massN, &
            " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — periodic KE
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
         "MEASURED KE agreement reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
      if (rel_diff >= AGREE_KE_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL periodic KE agreement: rank ", rank, &
            " decomp=", dec_ke, " ref=", ref_ke, &
            " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
      end if
   end block

   ! (a)+(c) Salt/heat closure + agreement — periodic leg
   call check_tracer_closure_and_agreement("periodic", rank, &
                                           ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                           dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- POISONED WALL leg ---
   ! Purpose: prove that the ghost-poison mechanism is transparent when
   ! every poisoned cell is overwritten by the subsequent exchange.  The
   ! poisoned decomposed run must agree with the UNPOISONED serial reference
   ! to the same tolerances as the unpoisoned decomposed run.
   !
   ! On a wall leg with no periodic axes, the band-selection rule gives:
   !   on np=1  : has_*=.true., periodic_*=.false. => no edges qualify => no-op
   !   on np=2  : MPI seam edges (.not. has_*=.true.) qualify; those ghost
   !              bands are poisoned then overwritten by ocean_halo_exchange_ml_state
   !              at the start of each stage.
   ! Serial reference: runs with poison=OFF (it is the reference).
   call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                             ref_salt0, ref_saltN, ref_heat0, ref_heatN)
   call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                       dec_salt0, dec_saltN, dec_heat0, dec_heatN, &
                       use_poison=.true.)

   ! (b) Finiteness — poisoned wall leg
   call assert_finite(dec_massN, "poison-wall decomposed final mass", rank)
   call assert_finite(dec_ke, "poison-wall decomposed final KE", rank)

   ! (a) Conservation — poisoned wall decomposed run
   block
      real(wp) :: drift
      drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL poison-wall conservation: rank ", rank, &
            " decomposed mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — poisoned wall mass + KE
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4)') &
         "MEASURED poison-wall mass agreement reldiff=", rel_diff
      if (rel_diff >= AGREE_MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL poison-wall mass agreement: rank ", rank, &
            " decomp=", dec_massN, " ref=", ref_massN, &
            " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
      end if
   end block
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
         "MEASURED poison-wall KE reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
      if (rel_diff >= AGREE_KE_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL poison-wall KE agreement: rank ", rank, &
            " decomp=", dec_ke, " ref=", ref_ke, &
            " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
      end if
   end block
   call check_tracer_closure_and_agreement("poison-wall", rank, &
                                           ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                           dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- POISONED PERIODIC leg ---
   ! On a periodic-x run bc%periodic_x = .true., so ALL ranks' west and east
   ! ghost bands are poisoned (single-rank: local wrap refills them; multi-rank:
   ! halo exchange + periodic wrap refill them).  Proves the wrap/exchange chain
   ! covers the bands on BOTH single-rank and multi-rank.
   call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                             ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                             use_periodic=.true.)
   call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                       dec_salt0, dec_saltN, dec_heat0, dec_heatN, &
                       use_periodic=.true., use_poison=.true.)

   ! (b) Finiteness — poisoned periodic leg
   call assert_finite(dec_massN, "poison-periodic decomposed final mass", rank)
   call assert_finite(dec_ke, "poison-periodic decomposed final KE", rank)

   ! (a) Conservation — poisoned periodic decomposed run
   block
      real(wp) :: drift
      drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
      if (drift >= MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4)') &
            "FAIL poison-periodic conservation: rank ", rank, &
            " decomposed mass drift=", drift, " >= tol=", MASS_TOL
      end if
   end block

   ! (c) Decomposition-invariance — poisoned periodic mass + KE
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4)') &
         "MEASURED poison-periodic mass agreement reldiff=", rel_diff
      if (rel_diff >= AGREE_MASS_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL poison-periodic mass agreement: rank ", rank, &
            " decomp=", dec_massN, " ref=", ref_massN, &
            " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
      end if
   end block
   block
      real(wp) :: rel_diff
      rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
      if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
         "MEASURED poison-periodic KE reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
      if (rel_diff >= AGREE_KE_TOL) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            "FAIL poison-periodic KE agreement: rank ", rank, &
            " decomp=", dec_ke, " ref=", ref_ke, &
            " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
      end if
   end block
   call check_tracer_closure_and_agreement("poison-periodic", rank, &
                                           ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                           dec_salt0, dec_saltN, dec_heat0, dec_heatN)

   ! --- BT_HALO=4 wall leg ---
   ! Proves the wide-halo march-in path produces the same trajectory as the
   ! serial reference (bt_halo=0).  This is exact — BIT-IDENTICAL, not
   ! merely within tolerance: at seams the march keeps the physical
   ! interior clean by construction (valid wide-band data within the
   ! 2-cell/substep creep budget), and at physical edges the substep's
   ! per-side effective-edge insets re-anchor v1's array-edge closures at
   ! the normal-array edge, so the band evolves exactly as v1's (see the
   ! physical-edge-insets bullet in the BT wide-halo march-in design notes; the
   ! naive un-inset widening measured KE reldiff 8.4 after 100 steps).
   ! At np=1 every reldiff below measures exactly 0.0; at np=2 the
   ! reldiffs equal the v1 wall-leg values digit-for-digit.
   ! Skip when ng_wide = NG + 4 = 7 would
   ! exceed min(nx_phys, ny_phys): on np=2 nxl=8 >= 7 so it runs; on np>2
   ! nxl=4 < 7 so we skip.  np=1 serial ref always has nx_phys=16 >= 7.
   !
   ! Counter expectations per outer step (2 RK2 stages; num_cycles=2):
   !   entry_exchange (per stage): +1 bt_group, +1 centre_2d, +3 face_x_2d, +3 face_y_2d
   !   in-loop wide exchange n=2  (per stage): +1 bt_group  (mod(2,2)==0 .and. 2<4)
   !   exit seam freshen     (per stage): +1 bt_group
   !   => per step: bt_group=6, centre_2d=2, face_x_2d=6, face_y_2d=6; bt_u_mid=0
   !   Over N_STEPS=100: bt_group=600, centre_2d=200, face_x_2d=600, face_y_2d=600
   !   Msgs (np=2 wall, x-only): 68 Isends/step × 100 = 6800 — see the
   !   per-primitive derivation in check_exchange_counts.
   if (nprocs <= 2) then
      call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                                ref_salt0, ref_saltN, ref_heat0, ref_heatN)
      call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                          dec_salt0, dec_saltN, dec_heat0, dec_heatN, &
                          use_bt_halo=4)

      ! (b) Finiteness — bt_halo wall leg
      call assert_finite(dec_massN, "bt_halo wall decomposed final mass", rank)
      call assert_finite(dec_ke, "bt_halo wall decomposed final KE", rank)

      ! (a) Conservation — bt_halo wall decomposed
      block
         real(wp) :: drift
         drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
         if (drift >= MASS_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4)') &
               "FAIL bt_halo conservation: rank ", rank, &
               " decomposed mass drift=", drift, " >= tol=", MASS_TOL
         end if
      end block

      ! (c) Decomposition-invariance — bt_halo wall mass + KE
      block
         real(wp) :: rel_diff
         rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
         if (rank == 0) write (*, '(a,es12.4)') &
            "MEASURED bt_halo wall mass agreement reldiff=", rel_diff
         if (rel_diff >= AGREE_MASS_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
               "FAIL bt_halo wall mass agreement: rank ", rank, &
               " decomp=", dec_massN, " ref=", ref_massN, &
               " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
         end if
      end block
      block
         real(wp) :: rel_diff
         rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
         if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
            "MEASURED bt_halo wall KE reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
         if (rel_diff >= AGREE_KE_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
               "FAIL bt_halo wall KE agreement: rank ", rank, &
               " decomp=", dec_ke, " ref=", ref_ke, &
               " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
         end if
      end block
      call check_tracer_closure_and_agreement("bt_halo-wall", rank, &
                                              ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                              dec_salt0, dec_saltN, dec_heat0, dec_heatN)

      ! --- POISONED BT_HALO=4 wall leg ---
      ! Ghost-band sentinels beyond ng_wide must be overwritten by entry_exchange;
      ! the poisoned run must agree with the unpoisoned serial reference.
      call run_serial_reference(ref_mass0, ref_massN, ref_ke, &
                                ref_salt0, ref_saltN, ref_heat0, ref_heatN)
      call run_decomposed(dec_mass0, dec_massN, dec_ke, &
                          dec_salt0, dec_saltN, dec_heat0, dec_heatN, &
                          use_bt_halo=4, use_poison=.true.)

      ! (b) Finiteness — poisoned bt_halo wall leg
      call assert_finite(dec_massN, "poison-bt_halo decomposed final mass", rank)
      call assert_finite(dec_ke, "poison-bt_halo decomposed final KE", rank)

      ! (a) Conservation — poisoned bt_halo decomposed
      block
         real(wp) :: drift
         drift = abs(dec_massN - dec_mass0)/max(abs(dec_mass0), 1.0_wp)
         if (drift >= MASS_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4)') &
               "FAIL poison-bt_halo conservation: rank ", rank, &
               " decomposed mass drift=", drift, " >= tol=", MASS_TOL
         end if
      end block

      ! (c) Decomposition-invariance — poisoned bt_halo mass + KE
      block
         real(wp) :: rel_diff
         rel_diff = abs(dec_massN - ref_massN)/max(abs(ref_massN), 1.0_wp)
         if (rank == 0) write (*, '(a,es12.4)') &
            "MEASURED poison-bt_halo mass agreement reldiff=", rel_diff
         if (rel_diff >= AGREE_MASS_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
               "FAIL poison-bt_halo mass agreement: rank ", rank, &
               " decomp=", dec_massN, " ref=", ref_massN, &
               " reldiff=", rel_diff, " >= tol=", AGREE_MASS_TOL
         end if
      end block
      block
         real(wp) :: rel_diff
         rel_diff = abs(dec_ke - ref_ke)/max(abs(ref_ke), 1.0_wp)
         if (rank == 0) write (*, '(a,es12.4,a,es12.4,a,es12.4)') &
            "MEASURED poison-bt_halo KE reldiff=", rel_diff, " ref_ke=", ref_ke, " dec_ke=", dec_ke
         if (rel_diff >= AGREE_KE_TOL) then
            n_fail = n_fail + 1
            write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
               "FAIL poison-bt_halo KE agreement: rank ", rank, &
               " decomp=", dec_ke, " ref=", ref_ke, &
               " reldiff=", rel_diff, " >= tol=", AGREE_KE_TOL
         end if
      end block
      call check_tracer_closure_and_agreement("poison-bt_halo", rank, &
                                              ref_salt0, ref_saltN, ref_heat0, ref_heatN, &
                                              dec_salt0, dec_saltN, dec_heat0, dec_heatN)
   end if

   ! --- Report ---
   comm = comm_env_compute_comm()
   call comm%barrier()
   if (n_fail > 0) then
      write (*, '(a,i0,a,i0,a,i0,a)') &
         "Rank ", rank, ": ", n_fail, " check(s) FAILED (nprocs=", nprocs, ")"
   else
      write (*, '(a,i0,a,i0,a)') &
         "Rank ", rank, ": all checks PASSED (nprocs=", nprocs, ")"
   end if
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   ! =======================================================================
   ! Island land-mask helper: init + fill cartesian metrics, seed a small
   ! rectangular island that straddles the np=2 x-seam (global cells
   ! i_global in [NX_G/2-1, NX_G/2+1], j_global in [NY_G/2-1, NY_G/2+1]),
   ! apply the land mask (which triggers the seam exchange in Part A),
   ! then map to device.
   !
   ! Call AFTER ocean_halo_init and ms%init; BEFORE ms%enter_data and the
   ! !$acc enter data copyin(ms, ...) call.
   ! =======================================================================
   subroutine make_cartesian_metrics_island(metrics, grid, ms, decomp, &
                                            nx_g_in, ny_g_in, &
                                            per_x, per_y, north_fold)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nx_g_in, ny_g_in
      logical, intent(in) :: per_x, per_y, north_fold

      integer :: i, j, k, i_global, j_global
      integer :: nxt, nyt, ng_loc, nz_loc

      nxt = grid%nx_total
      nyt = grid%ny_total
      ng_loc = grid%nghost
      nz_loc = ms%nz_ml

      ! Seed wet_mask so the test has TEETH for the Part A seam exchange.
      !
      ! Geometry: a vertical land strip on rank 0's LAST physical column
      ! (global i = NX_G/2 = 8 for NX_G=16), a few j-rows tall.  On the
      ! np=2 x-split rank 0 owns i_global 1..8, rank 1 owns 9..16; the seam
      ! u-face is between global cells 8 and 9.  The land sits ONLY on
      ! rank 0's side of the seam — rank 1 is all-wet near the seam.
      !
      ! Why this bites: the seam u-face at i_global=9 has
      ! `wet_u = wet_T(8)*wet_T(9)`.  On RANK 1 the cell i_global=8 is a WEST
      ! GHOST.  With the Part A seam exchange that ghost carries rank 0's
      ! real wet_T(8)=0 (land) => wet_u=0*1=0, the face is correctly CLOSED.
      ! WITHOUT the exchange the ghost is constant-extrapolated from rank 1's
      ! own edge (wet_T(9)=1) => wet_u=1*1=1, so rank 1 wrongly OPENS a face
      ! that abuts rank 0's land — flow leaks into the land shadow and the
      ! decomposed run diverges from the serial reference.  The serial run
      ! (px=1) always sees the true wet_T(8)=0, so serial-vs-decomposed
      ! agreement is the discriminator.
      !
      ! CRITICAL for teeth: seed only the PHYSICAL region, then
      ! constant-extrapolate the ghosts (mimicking the production bathymetry
      ! fill).  If we seeded ghosts by global index the seam ghost would be
      ! correct WITHOUT any exchange and the test would have no teeth.
      ms%wet_mask = 1.0_wp
      do j = ng_loc + 1, ng_loc + grid%ny_phys
         j_global = (j - ng_loc) + decomp%j_start - 1
         do i = ng_loc + 1, ng_loc + grid%nx_phys
            i_global = (i - ng_loc) + decomp%i_start - 1
            if (i_global == nx_g_in/2 .and. &
                j_global >= ny_g_in/2 - 1 .and. j_global <= ny_g_in/2 + 1) then
               ms%wet_mask(i, j) = 0.0_wp
            end if
         end do
      end do
      ! Constant-extrapolate wet_mask into the ghost band (production
      ! bathymetry-fill pattern): west/east ghost columns copy the nearest
      ! physical column, south/north ghost rows copy the nearest physical row.
      ! This deliberately puts the WRONG value in the seam ghost so only the
      ! Part A exchange can correct it.
      do j = 1, nyt
         do i = 1, ng_loc
            ms%wet_mask(i, j) = ms%wet_mask(ng_loc + 1, j)
            ms%wet_mask(nxt - ng_loc + i, j) = ms%wet_mask(nxt - ng_loc, j)
         end do
      end do
      do i = 1, nxt
         do j = 1, ng_loc
            ms%wet_mask(i, j) = ms%wet_mask(i, ng_loc + 1)
            ms%wet_mask(i, nyt - ng_loc + j) = ms%wet_mask(i, nyt - ng_loc)
         end do
      end do

      ! Floor h_layer to H_VANISHED at the land cells (same single-column
      ! strip as the wet_mask seeding above) so no layer has zero thickness;
      ! zero the layer velocities there.  The zeroed face metrics prevent
      ! flux across the land faces regardless.  Physical region only (ghosts
      ! are handled by the initial halo exchange / const-extrap).
      do k = 1, nz_loc
         do j = ng_loc + 1, ng_loc + grid%ny_phys
            j_global = (j - ng_loc) + decomp%j_start - 1
            do i = ng_loc + 1, ng_loc + grid%nx_phys
               i_global = (i - ng_loc) + decomp%i_start - 1
               if (i_global == nx_g_in/2 .and. &
                   j_global >= ny_g_in/2 - 1 .and. j_global <= ny_g_in/2 + 1) then
                  ms%h_layer(i, j, k) = H_VANISHED
                  ms%u_face_x_layer(i, j, k) = 0.0_wp
                  ms%v_face_y_layer(i, j, k) = 0.0_wp
               end if
            end do
         end do
      end do

      ! Build uniform Cartesian metrics, derive inverses.
      call metrics%init(grid)
      call metrics_fill_cartesian(metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(metrics)
      ! Multi-rank seam ghost fill of wet_mask BEFORE the land mask (O3): the
      ! seam ghost columns must carry the neighbour rank's real wet_T so the
      ! derived wet_u/wet_v/wet_q close the seam faces.  Mirrors what
      ! configure_ocean_land_mask does in production, incl. the is_decomposed
      ! guard (px=1 serial reference ⇒ no-op skip).  ocean_halo_init has
      ! already run in the calling run_* subroutine.
      ! Init-time host path: wet_mask is not yet device-mapped (this runs
      ! before ms%enter_data / the metrics copyin below).
      if (ocean_halo_is_decomposed()) then
         call ocean_halo_centre(ms%wet_mask, device_resident=.false.)
      end if
      call metrics_apply_land_mask(metrics, ms%wet_mask, grid, per_x, per_y, north_fold)
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
   end subroutine make_cartesian_metrics_island

   ! =======================================================================
   ! Serial reference: full NX_G x NY_G grid, px=py=1, rank=0 everywhere.
   ! Every rank builds the SAME identical serial run.
   ! =======================================================================
   subroutine run_serial_reference(mass0, massN, ke_out, salt0, saltN, heat0, heatN, &
                                   use_island, use_periodic)
      real(wp), intent(out) :: mass0, massN, ke_out
      real(wp), intent(out) :: salt0, saltN, heat0, heatN
         !! Global interior salt / heat sums (Σ hS·area / Σ hT·area) at step
         !! 0 and step N.  Salt is a UNIFORM tracer, heat is STRUCTURED — the
         !! salt-vs-heat discriminator for the decomposed-PPM seam leak.
      logical, intent(in), optional :: use_island
         !! When .true., seed a small island straddling the x-seam and apply
         !! the land mask before stepping.  Default .false. = all-wet.
      logical, intent(in), optional :: use_periodic
         !! When .true., apply periodic-x boundary conditions (west <-> east wrap).
         !! Default .false. = all-wall.

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
      type(ocean_bc_state_t) :: bc
      type(decomp_t) :: decomp

      integer :: step, ng_use
      real(wp) :: local_mass0, local_massN, local_ke
      logical :: island_run, periodic_run

      island_run = .false.
      if (present(use_island)) island_run = use_island
      periodic_run = .false.
      if (present(use_periodic)) periodic_run = use_periodic

      ! nghost = 3 on every leg: PPM 5-point seam-face stencil (see header).
      ng_use = NG

      ! Full global grid on px=py=1 (rank 0 everywhere — serial reference)
      call decomp_init(decomp, NX_G, NY_G, 1, 1, 0)
      call ocean_halo_init(decomp, ng_use, periodic_run, .false.)
      call ocean_halo_reserve(NZ, 0)

      call grid%init(NX_G, NY_G, ng_use, DX, DY)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = F_0
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

      ! Set boundary conditions.  Periodic-x: west+east both OBC_PERIODIC;
      ! validate derives periodic_x = .true. and checks nghost >= 3.
      ! Wall legs: all-WALL, single-rank so all edges physical.
      call ocean_bc_state_init(bc, grid, nz_ml=NZ)
      if (periodic_run) then
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
         ! Serial reference (px=1): west and east are the periodic wrap edges.
         call ocean_bc_state_set_edges(bc, .true., .true., .true., .true.)
      else
         call ocean_bc_state_set_edges(bc, .true., .true., .true., .true.)
      end if

      ! Seed IC and wind stress (full global domain)
      call seed_ic_and_wind(ms, ss, eos, decomp, NX_G, NY_G, grid%nx_total, grid%ny_total, ng_use, NZ)

      dyn%bt_work%bt_H_ref = real(NZ, wp)*(H0/real(NZ, wp))  ! = H0

      ! Initial halo fill
      call ocean_halo_exchange_ml_state(ms, device_resident=.false.)

      if (island_run) then
         call make_cartesian_metrics_island(metrics, grid, ms, decomp, &
                                            NX_G, NY_G, .false., .false., .false.)
      else
         call make_cartesian_metrics(metrics, grid)
      end if
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      call dyn%enter_data()

      ! Measure initial mass + tracer sums
      !$acc update self(ms%h_layer)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      local_mass0 = interior_mass_sum(ms%h_layer, NX_G, NY_G, ng_use, NZ, DX, DY)
      mass0 = local_mass0  ! serial — no allreduce needed; same on all ranks
      salt0 = interior_tracer_sum(ms%tracers(ms%idx_salinity)%hTr, NX_G, NY_G, ng_use, NZ, DX, DY)
      heat0 = interior_tracer_sum(ms%tracers(ms%idx_temperature)%hTr, NX_G, NY_G, ng_use, NZ, DX, DY)

      ! Step N_STEPS
      do step = 1, N_STEPS
         call ocean_dyn_step_split( &
            grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
            va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
      end do

      ! Pull state to host for measurement
      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)

      local_massN = interior_mass_sum(ms%h_layer, NX_G, NY_G, ng_use, NZ, DX, DY)
      local_ke = interior_ke_sum(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                                 NX_G, NY_G, ng_use, NZ, DX, DY)
      massN = local_massN
      ke_out = local_ke
      saltN = interior_tracer_sum(ms%tracers(ms%idx_salinity)%hTr, NX_G, NY_G, ng_use, NZ, DX, DY)
      heatN = interior_tracer_sum(ms%tracers(ms%idx_temperature)%hTr, NX_G, NY_G, ng_use, NZ, DX, DY)

      ! Cleanup
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
      call dyn%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call eos%destroy(); call ms%destroy()

      call ocean_halo_destroy()

   end subroutine run_serial_reference

   ! =======================================================================
   ! Decomposed run: each rank operates on its local subdomain.
   ! On np=1 this is identical to the serial reference.
   ! On np=2 the x-seam splits NX_G=16 into two halves.
   ! =======================================================================
   subroutine run_decomposed(mass0, massN, ke_out, salt0, saltN, heat0, heatN, &
                             use_island, use_periodic, use_poison, use_bt_halo)
      real(wp), intent(out) :: mass0, massN, ke_out
      real(wp), intent(out) :: salt0, saltN, heat0, heatN
         !! Global (allreduced) interior salt / heat sums at step 0 and N.
      logical, intent(in), optional :: use_island
         !! When .true., seed a small island straddling the x-seam and apply
         !! the land mask before stepping.  Default .false. = all-wet.
      logical, intent(in), optional :: use_periodic
         !! When .true., apply periodic-x boundary conditions (west <-> east wrap).
         !! Default .false. = all-wall.
      logical, intent(in), optional :: use_poison
         !! When .true., enable ghost-band sentinel-NaN (`dyn%poison_ghosts`).
         !! Every exchange-covered ghost is overwritten each step; the run must
         !! produce IDENTICAL results to the unpoisoned reference.  Default .false.
      integer, intent(in), optional :: use_bt_halo
         !! When > 0, enable the wide-halo BT march-in with this width.
         !! Default 0 = v1 per-substep exchange (bit-identical).

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
      type(ocean_bc_state_t) :: bc
      type(decomp_t) :: decomp

      integer :: nxl, nyl, step, ng_use
      integer :: bt_halo_val
      real(wp) :: local_mass0, local_massN, local_ke
      real(wp) :: global_mass0, global_massN, global_ke
      real(wp) :: local_salt0, local_saltN, local_heat0, local_heatN
      real(wp) :: global_salt0, global_saltN, global_heat0, global_heatN
      logical :: island_run, periodic_run, poison_run

      island_run = .false.
      if (present(use_island)) island_run = use_island
      periodic_run = .false.
      if (present(use_periodic)) periodic_run = use_periodic
      poison_run = .false.
      if (present(use_poison)) poison_run = use_poison
      bt_halo_val = 0
      if (present(use_bt_halo)) bt_halo_val = use_bt_halo

      ! nghost = 3 on every leg: PPM 5-point seam-face stencil (see header).
      ng_use = NG

      ! Decompose over actual rank count: x-split (px=nprocs, py=1)
      call decomp_init(decomp, NX_G, NY_G, nprocs, 1, rank)
      call ocean_halo_init(decomp, ng_use, periodic_run, .false.)
      ! Reserve worst-case: NZ layers + wide bt_halo if used.
      call ocean_halo_reserve(NZ, merge(NG + bt_halo_val, 0, bt_halo_val > 0))

      nxl = decomp%nx_local
      nyl = decomp%ny_local

      call grid%init(nxl, nyl, ng_use, DX, DY)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = F_0
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

      ! Set boundary conditions.  Periodic-x: west+east both OBC_PERIODIC;
      ! validate derives periodic_x = .true. and checks nghost >= 3.
      ! Wall legs: gate seam edges off (MPI seam is never a physical wall).
      call ocean_bc_state_init(bc, grid, nz_ml=NZ)
      if (periodic_run) then
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
         ! Under periodic-x decomposition every rank's west and east are
         ! periodic-neighbour links (rank 0's west wraps to the last rank,
         ! etc.).  has_west/has_east from decomp_init are .true. only at
         ! rank-column 0 (west) and px-1 (east); the periodic halo covers
         ! those wrap links; the bc wall-closure is gated off at seams by
         ! has_*.  South/north remain physical walls.
         call ocean_bc_state_set_edges(bc, decomp%has_west, decomp%has_east, &
                                       .true., .true.)
      else
         call ocean_bc_state_set_edges(bc, decomp%has_west, decomp%has_east, &
                                       decomp%has_south, decomp%has_north)
      end if

      ! Seed IC per rank using GLOBAL index (i_global = i_start + i_phys - 1)
      call seed_ic_and_wind(ms, ss, eos, decomp, NX_G, NY_G, grid%nx_total, grid%ny_total, ng_use, NZ)

      dyn%bt_work%bt_H_ref = H0  ! = NZ * (H0/NZ) = H0
      ! Poison-ghost knob: when enabled, sentinel-NaN the exchange-covered
      ! ghost bands at each outer-step start.  Set BEFORE enter_data so the
      ! flag is visible to the step loop on both CPU and GPU paths.
      dyn%poison_ghosts = poison_run

      ! Initial halo fill (fills seam ghosts from neighbours; periodic wrap for
      ! the domain-wrap edges under periodic_run)
      call ocean_halo_exchange_ml_state(ms, device_resident=.false.)

      if (island_run) then
         call make_cartesian_metrics_island(metrics, grid, ms, decomp, &
                                            NX_G, NY_G, .false., .false., .false.)
      else
         call make_cartesian_metrics(metrics, grid)
      end if
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      ! Set bt_halo BEFORE enter_data so the GPU-mapped struct has bt_halo=0
      ! (we want 0 there — the CPU dispatch checks the HOST copy after
      ! enable_bt_wide re-sets it; bt_wide arrays are attached separately by
      ! enable_bt_wide).  But we DO need the HOST copy to carry bt_halo_val
      ! before the step loop runs, so set it here on the host side.
      dyn%bt_halo = bt_halo_val
      call dyn%enter_data()
      ! Wide-halo march-in: allocate + GPU-attach bt_wide AFTER enter_data.
      ! No-op when bt_halo_val = 0.
      if (bt_halo_val > 0) then
         call ocean_dyn_enable_bt_wide(dyn, grid, DX, DY, &
                                       0.0_wp, 0.0_wp, 6.378e6_wp, &
                                       "cartesian", F_0, 0.0_wp, 0.0_wp, &
                                       "beta_plane")
      end if

      ! Measure initial mass + tracer sums (allreduce over ranks)
      !$acc update self(ms%h_layer)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
      local_mass0 = interior_mass_sum(ms%h_layer, nxl, nyl, ng_use, NZ, DX, DY)
      call halo_allreduce_sum(local_mass0, global_mass0)
      mass0 = global_mass0
      local_salt0 = interior_tracer_sum(ms%tracers(ms%idx_salinity)%hTr, nxl, nyl, ng_use, NZ, DX, DY)
      local_heat0 = interior_tracer_sum(ms%tracers(ms%idx_temperature)%hTr, nxl, nyl, ng_use, NZ, DX, DY)
      call halo_allreduce_sum(local_salt0, global_salt0)
      call halo_allreduce_sum(local_heat0, global_heat0)
      salt0 = global_salt0
      heat0 = global_heat0

      ! Step N_STEPS.  Counter reset AFTER the init-time ml_state fill so
      ! the expectations in check_exchange_counts cover exactly the step loop.
      call oh_counters_reset()
      do step = 1, N_STEPS
         call ocean_dyn_step_split( &
            grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
            va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
      end do

      ! EXCHANGE-COUNT GATE (deterministic, load-independent): the semantic
      ! halo-exchange counts over the step loop must match the hand-derived
      ! exchange topology of the split-RK2 step (the RK2 step comm-topology map).
      ! Any change to the exchange set (added / removed / re-cadenced
      ! exchange) fails here and must update the expectations consciously.
      call check_exchange_counts(nprocs > 1, periodic_run, bt_halo_val)

      ! Pull state to host for measurement
      !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
      !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
      !$acc update self(ms%tracers(ms%idx_temperature)%hTr)

      local_massN = interior_mass_sum(ms%h_layer, nxl, nyl, ng_use, NZ, DX, DY)
      local_ke = interior_ke_sum(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                                 nxl, nyl, ng_use, NZ, DX, DY)
      call halo_allreduce_sum(local_massN, global_massN)
      call halo_allreduce_sum(local_ke, global_ke)
      massN = global_massN
      ke_out = global_ke
      local_saltN = interior_tracer_sum(ms%tracers(ms%idx_salinity)%hTr, nxl, nyl, ng_use, NZ, DX, DY)
      local_heatN = interior_tracer_sum(ms%tracers(ms%idx_temperature)%hTr, nxl, nyl, ng_use, NZ, DX, DY)
      call halo_allreduce_sum(local_saltN, global_saltN)
      call halo_allreduce_sum(local_heatN, global_heatN)
      saltN = global_saltN
      heatN = global_heatN

      ! Cleanup
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      call destroy_cartesian_metrics(metrics)
      call ocean_bc_state_destroy(bc)
      call dyn%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call eos%destroy(); call ms%destroy()

      call ocean_halo_destroy()

   end subroutine run_decomposed

   ! =======================================================================
   ! Exchange-count gate: the semantic halo-exchange counters accumulated
   ! over the N_STEPS step loop must equal the hand-derived exchange
   ! topology of the split-RK2 step (the RK2 step comm-topology map):
   !
   !   ml_state  = 3 * 2 * N_STEPS     (stage entry + post-continuity +
   !                                    stage tail, per RK2 stage)
   !   centre_3d = 3 * 2 * N_STEPS     (continuity mid-Lie-split: h + the
   !                                    2 registered tracers, per stage)
   !   bt_group  = N_INNER * 2 * N_STEPS  (one grouped BT exchange per
   !                                       fast-loop substep)
   !   bt_u_mid  = N_INNER * 2 * N_STEPS on a REAL decomposition (px > 1),
   !               0 single-rank (the mid-substep u exchange is gated on
   !               ocean_halo_is_decomposed())
   !   every other category = 0
   !
   ! Deterministic and load-independent: counts are identical on all ranks
   ! by construction (all call sites unconditional or uniformly gated), so
   ! the gate runs on every rank.  Fires on each leg (wall/island/periodic —
   ! the exchange set is leg-invariant).
   ! =======================================================================
   subroutine check_exchange_counts(decomposed, periodic, bt_halo_in)
      logical, intent(in) :: decomposed
         !! .true. when the run has a real multi-rank seam (px > 1)
      logical, intent(in) :: periodic
         !! .true. when the run uses periodic-x (both wrap directions active)
      integer, intent(in), optional :: bt_halo_in
         !! bt_halo width; when > 0 the march-in counter schedule applies.
         !! Default 0 = v1 per-substep exchange.

      integer(int64) :: c_bt_group, c_bt_u_mid, c_ml_state
      integer(int64) :: c_centre_2d, c_centre_3d
      integer(int64) :: c_face_x_2d, c_face_x_3d, c_face_y_2d, c_face_y_3d
      integer(int64) :: expect_ml, expect_c3d, expect_btg, expect_btu
      integer(int64) :: expect_c2d, expect_fx2d, expect_fy2d
      integer(int64) :: expect_msgs
      integer :: bh

      bh = 0
      if (present(bt_halo_in)) bh = bt_halo_in

      expect_ml = int(3*2*N_STEPS, int64)
      expect_c3d = int(3*2*N_STEPS, int64)

      if (bh > 0) then
         ! Wide-halo march-in schedule (num_cycles = bh/2 = 2 for bh=4,
         ! N_INNER=4: in-loop wide exchange fires once per stage at n=2).
         ! Per outer step (2 stages):
         !   bt_group: entry_exchange(1) + in-loop(1) + exit-freshen(1) = 3 per stage = 6
         !   bt_u_mid: 0 (suppressed)
         !   centre_2d: entry_exchange(1) per stage = 2
         !   face_x_2d: entry_exchange(3) per stage = 6
         !   face_y_2d: entry_exchange(3) per stage = 6
         expect_btg = int(3*2*N_STEPS, int64)
         expect_btu = 0_int64
         expect_c2d = int(1*2*N_STEPS, int64)
         expect_fx2d = int(3*2*N_STEPS, int64)
         expect_fy2d = int(3*2*N_STEPS, int64)
         ! Msgs (np=2 wall rank 0, py=1, 1 active x-direction).  Each
         ! exchange posts ONE Isend per PRIMITIVE per direction — a grouped
         ! bt exchange carries 3 prims (eta, ubt, vbt) = 3 Isends, and a
         ! face_y array still exchanges its x-seam columns (py=1 only kills
         ! the y-direction sends, not the x-seam of v-located arrays).
         !   ml_state + continuity: 36 Isends/step (unchanged from v1)
         !   entry_exchange (per stage): bt_group_wide(3) + centre_wide(1)
         !     + 3×face_x_wide(3) + 3×face_y_wide(3) = 10 Isends; ×2 stages = 20
         !   in-loop wide (per stage): bt_group_wide(3) Isends; ×2 stages = 6
         !   exit freshen (per stage): bt_group(3) Isends; ×2 stages = 6
         !   total: 36 + 20 + 6 + 6 = 68 Isends/step × N_STEPS = 6800
         block
            integer(int64) :: isends_per_step
            isends_per_step = 36_int64 + 20_int64 + 6_int64 + 6_int64
            expect_msgs = isends_per_step*int(N_STEPS, int64)
         end block
      else
         ! v1 path.
         expect_btg = int(N_INNER*2*N_STEPS, int64)
         expect_btu = 0_int64
         if (decomposed) expect_btu = int(N_INNER*2*N_STEPS, int64)
         expect_c2d = 0_int64
         expect_fx2d = 0_int64
         expect_fy2d = 0_int64
         ! MPI isend count formula (nprocs==2 rank 0, py=1):
         !   Rank 0 = west edge rank.  X-directions active: 1 (wall) or 2 (periodic).
         !   Per outer step per direction:
         !     ml_state:    3 calls/stage × 2 stages × 5 prims (h, ux, vy, S, T) = 30
         !     continuity:  1 call/stage  × 2 stages × 3 prims (h, S, T)         = 6
         !     bt_group:    N_INNER/stage  × 2 stages × 3 prims (η, ubt, vbt)    = N_INNER*6
         !     bt_u_mid:    N_INNER/stage  × 2 stages × 1 prim  (ubt face_x)     = N_INNER*2
         !   = (36 + N_INNER*8) × N_STEPS × n_x_dirs
         block
            integer(int64) :: isends_per_step_per_dir, n_x_dirs
            n_x_dirs = 1_int64
            if (periodic .and. decomposed) n_x_dirs = 2_int64
            isends_per_step_per_dir = 36_int64 + int(N_INNER, int64)*8_int64
            expect_msgs = isends_per_step_per_dir*int(N_STEPS, int64)*n_x_dirs
         end block
      end if

      call oh_counters_get(c_bt_group, c_bt_u_mid, c_ml_state, &
                           c_centre_2d, c_centre_3d, &
                           c_face_x_2d, c_face_x_3d, c_face_y_2d, c_face_y_3d)

      call assert_count("ml_state", c_ml_state, expect_ml)
      call assert_count("centre_3d", c_centre_3d, expect_c3d)
      call assert_count("bt_group", c_bt_group, expect_btg)
      call assert_count("bt_u_mid", c_bt_u_mid, expect_btu)
      call assert_count("centre_2d", c_centre_2d, expect_c2d)
      call assert_count("face_x_2d", c_face_x_2d, expect_fx2d)
      call assert_count("face_x_3d", c_face_x_3d, 0_int64)
      call assert_count("face_y_2d", c_face_y_2d, expect_fy2d)
      call assert_count("face_y_3d", c_face_y_3d, 0_int64)

      ! MPI message count assertion (rank 0, nprocs==2 only — per-rank count
      ! varies by position for nprocs>2, so we only gate the 2-rank topology
      ! where rank 0 has exactly 1 active X-direction (wall) or 2 (periodic)).
      ! For bt_halo > 0 only pin the wall leg (np=2, not periodic) for now.
      if (nprocs == 2 .and. (.not. periodic .or. bh == 0)) then
         call assert_count("msgs", oh_counters_msgs(), expect_msgs)
      end if

   end subroutine check_exchange_counts

   subroutine assert_count(label, got, expected)
      !! Fail (bump host-associated n_fail) unless got == expected.
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: got, expected

      if (got /= expected) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,a,a,i0,a,i0)') &
            "FAIL exchange count: rank ", rank, " ", label, &
            " got ", got, " expected ", expected
      end if
   end subroutine assert_count

   ! =======================================================================
   ! Seed IC and wind stress.
   !
   ! IC: uniform h_layer = H0/NZ per layer; zero velocity;
   !     salinity / temperature at EOS reference values.
   ! Wind: 2-gyre tau_x(j_global) = TAUX*(1 - cos(2*pi*(j_global-0.5)/NY_G))
   !       which drives real flow.  The IC is zero-velocity so the first steps
   !       generate u crossing the x-seam.  tau_y = 0.
   !
   ! The IC is trivially decomposition-consistent (uniform h, zero vel).
   ! The wind stress uses the global j index so each rank seeds its own band.
   ! =======================================================================
   subroutine seed_ic_and_wind(ms, ss, eos, decomp, nx_g, ny_g, nxt, nyt, ng, nz_in)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(eos_t), intent(in) :: eos
      type(decomp_t), intent(in) :: decomp
      integer, intent(in) :: nx_g, ny_g, nxt, nyt, ng, nz_in

      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: i, j, k, i_global, j_global
      real(wp) :: h_per_layer, tau, amp, t_struct

      h_per_layer = H0/real(nz_in, wp)

      ! Non-uniform h with a sinusoidal x-y perturbation across the GLOBAL
      ! domain.  The perturbation is centered on x so it has non-zero gradient
      ! at the x-seam between the two ranks (i_global = NX_G/2 = 8).
      ! Without correct ghost-zone halos, the PPM continuity kernel reads
      ! stale h at the seam and the mass diverges from the serial reference
      ! within the first step.
      amp = 0.05_wp*h_per_layer   ! 5% amplitude — small enough to stay positive
      do k = 1, nz_in
         do j = 1, nyt
            do i = 1, nxt
               ! i_global: physical cells start at ng+1; i=ng+1 → i_global=i_start
               i_global = (i - ng) + decomp%i_start - 1
               j_global = (j - ng) + decomp%j_start - 1
               ms%h_layer(i, j, k) = h_per_layer + amp* &
                                     sin(2.0_wp*PI*real(i_global, wp)/real(nx_g, wp))* &
                                     cos(PI*real(j_global, wp)/real(ny_g, wp))
               ! Salinity: UNIFORM (S_ref everywhere).  Temperature: STRUCTURED
               ! in x — a sinusoid crossing the x-seam (i_global = NX_G/2).
               ! This is the salt-vs-heat discriminator that exposed the ng=2
               ! decomposed-PPM seam leak: horizontal tracer advection through
               ! a seam-inconsistent flux conserves a UNIFORM field exactly
               ! (Σ flux·C cancels) but LEAKS a STRUCTURED one (rank-A-out /=
               ! rank-B-in with T differing across the seam).  So a healthy
               ! seam keeps BOTH salt and heat at round-off; a broken seam
               ! (e.g. nghost=2) leaves salt clean while heat drifts + grows.
               t_struct = eos%T_ref*(1.0_wp + 0.20_wp* &
                                     sin(2.0_wp*PI*real(i_global, wp)/real(nx_g, wp)))
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = eos%S_ref*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_struct*ms%h_layer(i, j, k)
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      ! Wind stress: 2-gyre tau_x using global j index
      ss%tau_x = 0.0_wp
      ss%tau_y = 0.0_wp
      do j = ng + 1, ng + ny_g
         j_global = (j - ng) + decomp%j_start - 1
         tau = TAUX*(1.0_wp - cos(2.0_wp*PI*(real(j_global, wp) - 0.5_wp)/real(ny_g, wp)))
         do i = 1, nxt
            ss%tau_x(i, j) = tau
         end do
      end do

   end subroutine seed_ic_and_wind

   ! =======================================================================
   ! Compute interior mass sum on host.
   ! Interior = i in [ng+1, ng+nx_phys], j in [ng+1, ng+ny_phys], k in [1,nz].
   ! =======================================================================
   pure function interior_mass_sum(h_layer, nx_phys, ny_phys, ng, nz_in, dx, dy) result(total)
      integer, intent(in) :: nx_phys, ny_phys, ng, nz_in
      real(wp), intent(in) :: h_layer(ng + nx_phys + ng, ng + ny_phys + ng, nz_in)
      real(wp), intent(in) :: dx, dy
      real(wp) :: total
      integer :: i, j, k
      real(wp) :: area
      area = dx*dy
      total = 0.0_wp
      do k = 1, nz_in
         do j = ng + 1, ng + ny_phys
            do i = ng + 1, ng + nx_phys
               total = total + h_layer(i, j, k)*area
            end do
         end do
      end do
   end function interior_mass_sum

   ! =======================================================================
   ! Compute interior tracer sum on host:  sum_interior( hTr ) * area.
   ! Mirrors interior_mass_sum for a single tracer's hTr slice (PSU·m /
   ! °C·m per cell before the area weight).  Used for the salt/heat
   ! closure + agreement checks.
   ! =======================================================================
   pure function interior_tracer_sum(hTr, nx_phys, ny_phys, ng, nz_in, dx, dy) result(total)
      integer, intent(in) :: nx_phys, ny_phys, ng, nz_in
      real(wp), intent(in) :: hTr(ng + nx_phys + ng, ng + ny_phys + ng, nz_in)
      real(wp), intent(in) :: dx, dy
      real(wp) :: total
      integer :: i, j, k
      real(wp) :: area
      area = dx*dy
      total = 0.0_wp
      do k = 1, nz_in
         do j = ng + 1, ng + ny_phys
            do i = ng + 1, ng + nx_phys
               total = total + hTr(i, j, k)*area
            end do
         end do
      end do
   end function interior_tracer_sum

   ! =======================================================================
   ! Compute interior KE proxy sum on host:
   !   KE = sum_interior( 0.5 * h_layer * u_face^2 ) * area
   ! using the u-face at i (x-face) and v-face at j (y-face), averaged
   ! to cell centres as a simple proxy: u_c ~ 0.5*(u(i,j)+u(i+1,j)).
   ! This is deterministic across decompositions as long as the seam
   ! halos are correct (which is what we're testing).
   ! =======================================================================
   pure function interior_ke_sum(h_layer, u_face_x, v_face_y, &
                                 nx_phys, ny_phys, ng, nz_in, dx, dy) result(total)
      integer, intent(in) :: nx_phys, ny_phys, ng, nz_in
      real(wp), intent(in) :: h_layer(ng + nx_phys + ng, ng + ny_phys + ng, nz_in)
      real(wp), intent(in) :: u_face_x(ng + nx_phys + ng + 1, ng + ny_phys + ng, nz_in)
      real(wp), intent(in) :: v_face_y(ng + nx_phys + ng, ng + ny_phys + ng + 1, nz_in)
      real(wp), intent(in) :: dx, dy
      real(wp) :: total
      integer :: i, j, k
      real(wp) :: uc, vc, area
      area = dx*dy
      total = 0.0_wp
      do k = 1, nz_in
         do j = ng + 1, ng + ny_phys
            do i = ng + 1, ng + nx_phys
               uc = 0.5_wp*(u_face_x(i, j, k) + u_face_x(i + 1, j, k))
               vc = 0.5_wp*(v_face_y(i, j, k) + v_face_y(i, j + 1, k))
               total = total + 0.5_wp*h_layer(i, j, k)*(uc*uc + vc*vc)*area
            end do
         end do
      end do
   end function interior_ke_sum

   ! =======================================================================
   ! Simple finiteness check; increments n_fail.
   ! =======================================================================
   subroutine assert_finite(val, label, rnk)
      real(wp), intent(in) :: val
      character(len=*), intent(in) :: label
      integer, intent(in) :: rnk
      if (.not. (val == val) .or. abs(val) > huge(val)*0.5_wp) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,a,a,es12.4)') &
            "FAIL finiteness: rank ", rnk, " ", label, " = ", val
      end if
   end subroutine assert_finite

   ! =======================================================================
   ! Salt + heat closure (drift vs step 0) AND decomposition-invariance
   ! (decomp vs serial) checks for one leg.  Prints the MEASURED reldiffs
   ! (rank 0) so a regression is diagnosable from the log.  Salt is a
   ! UNIFORM tracer (should stay clean under any seam bug); heat is
   ! STRUCTURED (leaks through a seam-inconsistent flux) — this is the
   ! salt-vs-heat discriminator that exposed the ng=2 decomposed-PPM leak.
   ! =======================================================================
   subroutine check_tracer_closure_and_agreement(leg, rnk, &
                                                 ref_s0, ref_sN, ref_h0, ref_hN, &
                                                 dec_s0, dec_sN, dec_h0, dec_hN)
      character(len=*), intent(in) :: leg
      integer, intent(in) :: rnk
      real(wp), intent(in) :: ref_s0, ref_sN, ref_h0, ref_hN
      real(wp), intent(in) :: dec_s0, dec_sN, dec_h0, dec_hN
      real(wp) :: d

      ! (a) CONSERVATION — serial reference (drift vs its own step 0)
      d = abs(ref_sN - ref_s0)/max(abs(ref_s0), 1.0_wp)
      call one_check(leg//" ref salt closure", d, TRACER_CLOSURE_TOL, rnk)
      d = abs(ref_hN - ref_h0)/max(abs(ref_h0), 1.0_wp)
      call one_check(leg//" ref heat closure", d, TRACER_CLOSURE_TOL, rnk)

      ! (a) CONSERVATION — decomposed run (the leg that catches the seam leak)
      d = abs(dec_sN - dec_s0)/max(abs(dec_s0), 1.0_wp)
      if (rnk == 0) write (*, '(a,a,a,es12.4)') &
         "MEASURED ", leg, " decomp salt closure drift=", d
      call one_check(leg//" decomp salt closure", d, TRACER_CLOSURE_TOL, rnk)
      d = abs(dec_hN - dec_h0)/max(abs(dec_h0), 1.0_wp)
      if (rnk == 0) write (*, '(a,a,a,es12.4)') &
         "MEASURED ", leg, " decomp heat closure drift=", d
      call one_check(leg//" decomp heat closure", d, TRACER_CLOSURE_TOL, rnk)

      ! (c) DECOMPOSITION-INVARIANCE — decomp final sum vs serial final sum
      d = abs(dec_sN - ref_sN)/max(abs(ref_sN), 1.0_wp)
      if (rnk == 0) write (*, '(a,a,a,es12.4)') &
         "MEASURED ", leg, " salt agreement reldiff=", d
      call one_check(leg//" salt agreement", d, TRACER_AGREE_TOL, rnk)
      d = abs(dec_hN - ref_hN)/max(abs(ref_hN), 1.0_wp)
      if (rnk == 0) write (*, '(a,a,a,es12.4)') &
         "MEASURED ", leg, " heat agreement reldiff=", d
      call one_check(leg//" heat agreement", d, TRACER_AGREE_TOL, rnk)
   end subroutine check_tracer_closure_and_agreement

   subroutine one_check(label, measured, tol, rnk)
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: measured, tol
      integer, intent(in) :: rnk
      if (measured >= tol) then
         n_fail = n_fail + 1
         write (*, '(a,i0,a,a,a,es12.4,a,es12.4)') &
            "FAIL ", rnk, " ", label, ": measured=", measured, " >= tol=", tol
      end if
   end subroutine one_check

end program test_ocean_dyn_mpi
#endif
