!! Unit tests for the Adcroft (2013) porous-barrier fit (`rdb_ocean_porous`)
!! and its seam into the ocean continuity-PPM transport.
!!
!! The three curve tests pin the fit at its two endpoints and in between:
!! a face whose subgrid sill is FULLY SUBMERGED must give the full width
!! back EXACTLY, a face that is FULLY BLOCKED must give exactly zero, and
!! a partially blocked face must reproduce a hand-computed intermediate in
!! each of the three branches of the fit (`m < 1/2`, `m = 1/2`, `m > 1/2`).
!!
!! The kernel tests do the same on the layer-averaged open-area fractions
!! the transport actually consumes, including the flat-bottom degeneracy
!! (every fraction exactly 1 — the property the default-off bit-identity
!! argument rests on) and the masking-depth gate.
!!
!! The conservation test is the one that catches a wrong seam INTO
!! continuity: narrowing a face width may redistribute mass but must never
!! create or destroy it, because the flux divergence telescopes.  It is
!! paired with a non-vacuity check — the porous run must actually differ
!! from the un-narrowed one, or the conservation assertion proves nothing.
!!
!! The device test poisons the HOST copy of the mapped inputs so the two
!! memories disagree; only a kernel that read the DEVICE copy can pass.
!! See `tests/test_ocean_diag_reduce.F90` for the full rationale — for
!! explicit-shape `do concurrent` dummies nvfortran emits an implicit
!! `copyin`, so an unmapped buffer stages STALE HOST data rather than
!! faulting and a plain map-and-compare would pass vacuously.
module test_ocean_porous
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, continuity_compute_fluxes, &
                             continuity_apply_fluxes, continuity_zonal_flux, &
                             continuity_meridional_flux
   use rdb_coriolis_adv, only: coriolis_adv_t, &
                               coriolis_adv_compute_tendencies_hk, &
                               coriolis_adv_compute_tendencies_sadourny_energy
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear_interior
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_porous
   use rdb_config, only: config_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_finalize, metrics_porous_alloc
   use rdb_ocean_porous, only: porous_open_width, porous_cum_area, &
                               porous_eta_face, porous_update_face_areas, &
                               porous_fill_stats_resolved, &
                               porous_stats_are_ordered, &
                               parse_porous_source, parse_porous_eta_interp, &
                               POROUS_SOURCE_RESOLVED, POROUS_SOURCE_FILE, &
                               POROUS_ETA_MAX, POROUS_ETA_MIN, &
                               POROUS_ETA_ARITH
   implicit none
   private

   public :: collect_ocean_porous_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

   real(wp), allocatable :: last_bt_u(:, :), last_bt_v(:, :)
      !! Column-integrated open fraction (unit widths) from the most
      !! recent `run_update`, so a case can assert on it without
      !! threading another argument through every caller.

   real(wp), parameter :: TOL_ULP = 1.0e-14_wp
      !! Relative bound for the per-face narrowing identities.  They hold
      !! EXACTLY in the code — the transport is literally multiplied by the
      !! open fraction — but the comparison cannot be written as an
      !! equality: under `-fast` nvfortran contracts `got - plain*por` into
      !! an FMA, so the difference measured is the (nonzero, half-ulp)
      !! rounding error of `plain*por` rather than a real discrepancy.  A
      !! few-ulp bound is still a hard mutation gate: dropping a porous
      !! factor changes the transport by 25-100%, not by 1e-14.

   real(wp), parameter :: POISON = -7.77e5_wp
      !! Host-side poison for the device-residency test.  Chosen far
      !! outside any physical bottom elevation / thickness used here, so a
      !! kernel that reads it produces an unmistakably wrong answer.

contains

   subroutine collect_ocean_porous_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("curve_submerged_gives_full_width", &
                               test_curve_submerged), &
                  new_unittest("curve_blocked_gives_zero_width", &
                               test_curve_blocked), &
                  new_unittest("curve_partial_hand_computed", &
                               test_curve_partial), &
                  new_unittest("cum_area_reproduces_mean_height", &
                               test_cum_area_mean), &
                  new_unittest("cum_area_integrates_open_width", &
                               test_cum_area_quadrature), &
                  new_unittest("eta_face_interp_options", &
                               test_eta_face_interp), &
                  new_unittest("enum_parsing", test_enum_parsing), &
                  new_unittest("kernel_flat_bottom_is_inert", &
                               test_kernel_flat_bottom), &
                  new_unittest("kernel_sill_partial_and_full", &
                               test_kernel_sill_partial), &
                  new_unittest("kernel_sill_fully_blocks_bottom", &
                               test_kernel_sill_blocked), &
                  new_unittest("kernel_masking_depth_leaves_open", &
                               test_kernel_masking), &
                  new_unittest("stats_resolved_flat_is_degenerate", &
                               test_stats_resolved_flat), &
                  new_unittest("continuity_conserves_mass_with_porous", &
                               test_continuity_conservation), &
                  new_unittest("bt_width_carries_column_fraction", &
                               test_bt_width_column), &
                  new_unittest("degenerate_sill_needs_no_division", &
                               test_degenerate_sill), &
                  new_unittest("renormalised_split_flux_matches_uhbt", &
                               test_renorm_split), &
                  new_unittest("coriolis_transport_is_narrowed", &
                               test_coriolis_narrowed), &
                  new_unittest("coriolis_energy_transport_is_narrowed", &
                               test_coriolis_energy_narrowed), &
                  new_unittest("face_areas_are_device_resident", &
                               test_device_resident), &
                  new_unittest("stats_resolved_simpson_weighting", &
                               test_stats_simpson), &
                  new_unittest("stats_resolved_ignores_land_corners", &
                               test_stats_land_gating), &
                  new_unittest("stats_ordering_validator", &
                               test_stats_ordering), &
                  new_unittest("along_face_uniform_ridge_is_not_walled", &
                               test_ridge_not_walled), &
                  new_unittest("vanished_column_blocks_bt_width", &
                               test_vanished_column), &
                  new_unittest("single_layer_column_is_handled", &
                               test_nz1_column), &
                  new_unittest("eta_interp_max_blocks_least", &
                               test_eta_interp_monotone), &
                  new_unittest("split_path_narrows_both_directions", &
                               test_split_narrows_both), &
                  new_unittest("fused_path_narrows_both_directions", &
                               test_fused_narrows_both), &
                  new_unittest("renormalisation_keeps_blocked_layer_at_zero", &
                               test_renorm_blocked_layer), &
                  new_unittest("bt_substep_consumes_bt_widths", &
                               test_bt_substep_consumes), &
                  new_unittest("configure_porous_sign_conventions", &
                               test_configure_signs) &
                  ]
   end subroutine collect_ocean_porous_tests

   ! =================================================================
   ! Curve — the three points that pin the fit
   ! =================================================================

   subroutine test_curve_submerged(error)
      !! A subgrid sill entirely BELOW the interface leaves the face
      !! fully open: the width fraction is EXACTLY one (not merely close),
      !! because the fit short-circuits above `d_max` rather than
      !! evaluating the monomial.  Checked in all three `m` branches.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp
      real(wp) :: davg(3)
      integer :: n

      davg = [-90.0_wp, -75.0_wp, -60.0_wp]   ! m = 0.2, 0.5, 0.8
      checks: block
         do n = 1, 3
            call check(error, porous_open_width(DMIN, DMAX, davg(n), -49.0_wp) == 1.0_wp, &
                       "submerged sill must give EXACTLY the full width back")
            if (allocated(error)) exit checks
            call check(error, porous_open_width(DMIN, DMAX, davg(n), 0.0_wp) == 1.0_wp, &
                       "sea-surface interface must give EXACTLY the full width")
            if (allocated(error)) exit checks
            ! Just above d_max the cumulative area must be continuous with
            ! the interior branch value (d_max - d_min)*(1 - m).
            call check(error, abs(porous_cum_area(DMIN, DMAX, davg(n), DMAX) - &
                                  (DMAX - DMIN)*(1.0_wp - (davg(n) - DMIN)/(DMAX - DMIN))) &
                       < 1.0e-12_wp, "cum_area must be continuous at d_max")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_curve_submerged

   subroutine test_curve_blocked(error)
      !! An interface at or below the DEEPEST along-face point sees a
      !! fully blocked face: EXACTLY zero open width and zero cumulative
      !! area.  Again exact, in all three `m` branches.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp
      real(wp) :: davg(3)
      integer :: n

      davg = [-90.0_wp, -75.0_wp, -60.0_wp]
      checks: block
         do n = 1, 3
            call check(error, porous_open_width(DMIN, DMAX, davg(n), DMIN) == 0.0_wp, &
                       "interface AT the deepest point must give EXACTLY zero width")
            if (allocated(error)) exit checks
            call check(error, porous_open_width(DMIN, DMAX, davg(n), -120.0_wp) == 0.0_wp, &
                       "interface below the sill must give EXACTLY zero width")
            if (allocated(error)) exit checks
            call check(error, porous_cum_area(DMIN, DMAX, davg(n), -120.0_wp) == 0.0_wp, &
                       "cumulative area below the sill must be EXACTLY zero")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_curve_blocked

   subroutine test_curve_partial(error)
      !! Hand-computed intermediates, one per branch of the fit, at the
      !! midpoint of the sill range (`zeta = 1/2`):
      !!
      !!   m = 1/2 (d_avg = -75): w = zeta            = 0.5 exactly
      !!   m = 1/5 (d_avg = -90): a = 4,   w = zeta**(1/4)
      !!   m = 4/5 (d_avg = -60): a = 1/4, w = 1 - (1-zeta)**(1/4)
      !!
      !! with `0.5**0.25 = 0.8408964152537145`.  The cumulative areas at
      !! the same point follow from the closed forms and are pinned too.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp, ETA = -75.0_wp
      real(wp), parameter :: Q = 0.8408964152537145_wp      ! 0.5**0.25
      real(wp), parameter :: P = 0.4204482076268558_wp      ! 0.5**1.25
      real(wp), parameter :: TOL = 1.0e-12_wp

      checks: block
         ! --- m = 1/2: the linear branch, exactly representable ---
         call check(error, porous_open_width(DMIN, DMAX, -75.0_wp, ETA) == 0.5_wp, &
                    "m=1/2 open width at the sill midpoint must be EXACTLY 0.5")
         if (allocated(error)) exit checks
         ! A = (d_max-d_min)*0.5*zeta**2 = 50*0.5*0.25 = 6.25
         call check(error, porous_cum_area(DMIN, DMAX, -75.0_wp, ETA) == 6.25_wp, &
                    "m=1/2 cumulative area must be EXACTLY 6.25")
         if (allocated(error)) exit checks

         ! --- m = 1/5: w = zeta**(1/4) ---
         call check(error, abs(porous_open_width(DMIN, DMAX, -90.0_wp, ETA) - Q) < TOL, &
                    "m=1/5 open width must match the hand-computed 0.5**0.25")
         if (allocated(error)) exit checks
         ! A = 50*(1-m)*zeta**(1/(1-m)) = 50*0.8*0.5**1.25 = 40*P
         call check(error, abs(porous_cum_area(DMIN, DMAX, -90.0_wp, ETA) - 40.0_wp*P) < TOL, &
                    "m=1/5 cumulative area must match 40*0.5**1.25")
         if (allocated(error)) exit checks

         ! --- m = 4/5: w = 1 - (1-zeta)**(1/4) ---
         call check(error, abs(porous_open_width(DMIN, DMAX, -60.0_wp, ETA) - &
                               (1.0_wp - Q)) < TOL, &
                    "m=4/5 open width must match the hand-computed 1 - 0.5**0.25")
         if (allocated(error)) exit checks
         ! A = 50*(zeta - m + m*(1-zeta)**(1/m)) = 50*(0.5 - 0.8 + 0.8*0.5**1.25)
         call check(error, abs(porous_cum_area(DMIN, DMAX, -60.0_wp, ETA) - &
                               50.0_wp*(0.5_wp - 0.8_wp + 0.8_wp*P)) < TOL, &
                    "m=4/5 cumulative area must match the closed form")
         if (allocated(error)) exit checks

         ! Monotonicity in d_avg: a shallower MEAN sill blocks more.
         call check(error, porous_open_width(DMIN, DMAX, -90.0_wp, ETA) > &
                    porous_open_width(DMIN, DMAX, -75.0_wp, ETA), &
                    "a deeper mean sill must leave MORE of the face open")
         if (allocated(error)) exit checks
         call check(error, porous_open_width(DMIN, DMAX, -75.0_wp, ETA) > &
                    porous_open_width(DMIN, DMAX, -60.0_wp, ETA), &
                    "a shallower mean sill must leave LESS of the face open")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_curve_partial

   subroutine test_cum_area_mean(error)
      !! The defining property of the three-parameter fit: averaged over
      !! the whole sill range the open fraction equals the fraction of the
      !! range that lies above the MEAN height,
      !!   (A(d_max) - A(d_min)) / (d_max - d_min) = (d_max - d_avg)/(d_max - d_min).
      !! If the fit did not reproduce the prescribed `d_avg` this is the
      !! assertion that would catch it.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp
      real(wp) :: davg(5), got, want
      integer :: n

      davg = [-95.0_wp, -90.0_wp, -75.0_wp, -60.0_wp, -55.0_wp]
      checks: block
         do n = 1, 5
            got = (porous_cum_area(DMIN, DMAX, davg(n), DMAX) - &
                   porous_cum_area(DMIN, DMAX, davg(n), DMIN))/(DMAX - DMIN)
            want = (DMAX - davg(n))/(DMAX - DMIN)
            call check(error, abs(got - want) < 1.0e-12_wp, &
                       "the fit must reproduce the prescribed mean sill height")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_cum_area_mean

   subroutine test_cum_area_quadrature(error)
      !! `porous_cum_area` must be the exact vertical integral of
      !! `porous_open_width` — that is what makes the layer-averaged
      !! fraction a difference of two cumulative areas rather than a
      !! quadrature.  Verified against composite Simpson over the interior
      !! of the sill range (the endpoints carry an infinite slope for
      !! `m /= 1/2`, so the comparison is made over `[d_min + s, eta]`
      !! with the analytic offset subtracted).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp
      integer, parameter :: NQ = 20000
      real(wp) :: davg(3), lo, hi, hstep, quad, got, x
      integer :: n, q

      davg = [-90.0_wp, -75.0_wp, -60.0_wp]
      ! DELIBERATELY ASYMMETRIC about the sill midpoint: on a symmetric
      ! interval the m=1/2 quadratic and linear forms give the same
      ! difference, so a `0.5*zeta*zeta` -> `0.5*zeta` slip would slide
      ! straight through.
      lo = -92.0_wp
      hi = -57.0_wp
      hstep = (hi - lo)/real(NQ, wp)
      checks: block
         do n = 1, 3
            quad = porous_open_width(DMIN, DMAX, davg(n), lo) + &
                   porous_open_width(DMIN, DMAX, davg(n), hi)
            do q = 1, NQ - 1
               x = lo + real(q, wp)*hstep
               if (mod(q, 2) == 1) then
                  quad = quad + 4.0_wp*porous_open_width(DMIN, DMAX, davg(n), x)
               else
                  quad = quad + 2.0_wp*porous_open_width(DMIN, DMAX, davg(n), x)
               end if
            end do
            quad = quad*hstep/3.0_wp
            got = porous_cum_area(DMIN, DMAX, davg(n), hi) - &
                  porous_cum_area(DMIN, DMAX, davg(n), lo)
            call check(error, abs(got - quad) < 1.0e-6_wp*max(1.0_wp, abs(quad)), &
                       "cum_area must be the vertical integral of open_width")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_cum_area_quadrature

   subroutine test_eta_face_interp(error)
      !! The four interface-at-velocity-point rules.  MAX (the default)
      !! picks the SHALLOWER interface — the most blocking choice.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: ZA = -80.0_wp, ZB = -40.0_wp

      checks: block
         call check(error, porous_eta_face(ZA, ZB, POROUS_ETA_MAX) == ZB, &
                    "MAX must pick the shallower interface")
         if (allocated(error)) exit checks
         call check(error, porous_eta_face(ZA, ZB, POROUS_ETA_MIN) == ZA, &
                    "MIN must pick the deeper interface")
         if (allocated(error)) exit checks
         call check(error, porous_eta_face(ZA, ZB, POROUS_ETA_ARITH) == -60.0_wp, &
                    "ARITHMETIC must be the mean of the two interfaces")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_eta_face_interp

   subroutine test_enum_parsing(error)
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_porous_source("resolved") == POROUS_SOURCE_RESOLVED, &
                    "source='resolved' must parse")
         if (allocated(error)) exit checks
         call check(error, parse_porous_source("file") == POROUS_SOURCE_FILE, &
                    "source='file' must parse (configure then fails loud)")
         if (allocated(error)) exit checks
         call check(error, parse_porous_source("nonsense") < 0, &
                    "an unknown source must be rejected, not silently defaulted")
         if (allocated(error)) exit checks
         call check(error, parse_porous_eta_interp("harmonic") >= 0, &
                    "eta_interp='harmonic' must parse")
         if (allocated(error)) exit checks
         call check(error, parse_porous_eta_interp("nonsense") < 0, &
                    "an unknown eta_interp must be rejected")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_enum_parsing

   ! =================================================================
   ! Kernel — layer-averaged open-area fractions
   ! =================================================================

   subroutine run_update(nx, ny, interp, mask_depth, bed, h_layer, &
                         dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, &
                         por_u, por_v)
      !! Map every argument, run the device kernel, pull the answers back.
      !! `create` on the outputs is deliberate: they are pure kernel
      !! products, so there is nothing to stage in.
      integer, intent(in) :: nx, ny, interp
      real(wp), intent(in) :: mask_depth
      real(wp), intent(in) :: bed(nx, ny), h_layer(nx, ny, NZ)
      real(wp), intent(in) :: dmin_u(nx + 1, ny), dmax_u(nx + 1, ny), davg_u(nx + 1, ny)
      real(wp), intent(in) :: dmin_v(nx, ny + 1), dmax_v(nx, ny + 1), davg_v(nx, ny + 1)
      real(wp), intent(out) :: por_u(nx + 1, ny, NZ), por_v(nx, ny + 1, NZ)
      real(wp) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp) :: dy_cu_bt(nx + 1, ny), dx_cv_bt(nx, ny + 1)

      ! Unit widths make `dy_cu_bt` read directly as the column-integrated
      ! open fraction.
      dy_cu = 1.0_wp
      dx_cv = 1.0_wp
      !$acc enter data copyin(bed, h_layer, dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v)
      !$acc enter data copyin(dy_cu, dx_cv)
      !$acc enter data create(por_u, por_v, dy_cu_bt, dx_cv_bt)
      call porous_update_face_areas(nx, ny, NZ, interp, mask_depth, bed, h_layer, &
                                    dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, &
                                    dy_cu, dx_cv, por_u, por_v, dy_cu_bt, dx_cv_bt)
      !$acc update self(por_u, por_v, dy_cu_bt, dx_cv_bt)
      !$acc exit data delete(bed, h_layer, dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v)
      !$acc exit data delete(dy_cu, dx_cv)
      !$acc exit data delete(por_u, por_v, dy_cu_bt, dx_cv_bt)
      last_bt_u = dy_cu_bt
      last_bt_v = dx_cv_bt
   end subroutine run_update

   subroutine test_kernel_flat_bottom(error)
      !! A flat along-face seafloor makes `d_min = d_max = d_avg`, the fit
      !! degenerates to the binary open/closed step, and every layer above
      !! the bed comes back EXACTLY fully open.  This is the property the
      !! default-off bit-identity argument leans on: with the resolved-
      !! bathymetry source a flat basin is a literal no-op.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)

      bed = -100.0_wp
      h_layer = 25.0_wp
      call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                      dmin_u, dmax_u, davg_u, &
                                      dmin_v, dmax_v, davg_v)
      call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      checks: block
         call check(error, all(por_u == 1.0_wp), &
                    "flat bottom must leave EVERY u-face fraction exactly 1")
         if (allocated(error)) exit checks
         call check(error, all(por_v == 1.0_wp), &
                    "flat bottom must leave EVERY v-face fraction exactly 1")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_kernel_flat_bottom

   subroutine test_kernel_sill_partial(error)
      !! Hand-computed layer fractions across a linear (`m = 1/2`) sill.
      !!
      !! Bed at -100 m, four 25 m layers, so the face interfaces sit at
      !! -100, -75, -50, -25, 0.  The sill spans `d_min = -100` to
      !! `d_max = -50` with `d_avg = -75`, giving `A(eta) = 50*0.5*zeta**2`
      !! inside the range and `A = eta + 75` above it:
      !!
      !!   layer 1 (-100 -> -75): (6.25 - 0)/25    = 0.25
      !!   layer 2 ( -75 -> -50): (25 - 6.25)/25   = 0.75
      !!   layer 3 ( -50 -> -25): (50 - 25)/25     = 1     (sill submerged)
      !!   layer 4 ( -25 ->   0): (75 - 50)/25     = 1
      !!
      !! Layers 3 and 4 must come back EXACTLY 1 — the "full width given
      !! back" endpoint — and layers 1 and 2 exactly 0.25 / 0.75, all three
      !! values being exactly representable in binary.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)
      real(wp) :: want(NZ)
      integer :: k

      bed = -100.0_wp
      h_layer = 25.0_wp
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp
      call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      want = [0.25_wp, 0.75_wp, 1.0_wp, 1.0_wp]
      checks: block
         do k = 1, NZ
            ! Interior faces only: the array-edge faces (i=1, i=nx+1) have
            ! no adjacent cell pair and are held fully open by design.
            call check(error, all(por_u(2:NXT, :, k) == want(k)), &
                       "u-face layer fraction must match the hand-computed value")
            if (allocated(error)) exit checks
            call check(error, all(por_v(:, 2:NYT, k) == want(k)), &
                       "v-face layer fraction must match the hand-computed value")
            if (allocated(error)) exit checks
         end do
      end block checks
   end subroutine test_kernel_sill_partial

   subroutine test_kernel_sill_blocked(error)
      !! A sill that sits entirely ABOVE the two lowest layers blocks them
      !! completely: exactly zero open area, no leak.
      !!
      !! Bed -100 m, four 25 m layers (interfaces -100, -75, -50, -25, 0),
      !! sill `d_min = -40`, `d_max = -30`, `d_avg = -35`:
      !!   layer 1 (-100 -> -75): both interfaces below d_min => 0 exactly
      !!   layer 2 ( -75 -> -50): both below d_min             => 0 exactly
      !!   layer 3 ( -50 -> -25): (A(-25) - 0)/25 = 10/25      = 0.4
      !!   layer 4 ( -25 ->   0): (35 - 10)/25                 = 1
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)

      bed = -100.0_wp
      h_layer = 25.0_wp
      dmin_u = -40.0_wp; dmax_u = -30.0_wp; davg_u = -35.0_wp
      dmin_v = -40.0_wp; dmax_v = -30.0_wp; davg_v = -35.0_wp
      call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      checks: block
         call check(error, all(por_u(2:NXT, :, 1) == 0.0_wp), &
                    "a fully blocked bed layer must give EXACTLY zero open area")
         if (allocated(error)) exit checks
         call check(error, all(por_u(2:NXT, :, 2) == 0.0_wp), &
                    "a fully blocked second layer must give EXACTLY zero open area")
         if (allocated(error)) exit checks
         call check(error, all(abs(por_u(2:NXT, :, 3) - 0.4_wp) < 1.0e-13_wp), &
                    "the straddling layer must match the hand-computed 0.4")
         if (allocated(error)) exit checks
         call check(error, all(por_u(2:NXT, :, 4) == 1.0_wp), &
                    "the layer above the sill must give EXACTLY the full width")
         if (allocated(error)) exit checks
         call check(error, all(por_v(:, 2:NYT, 1) == 0.0_wp), &
                    "v-faces must block identically")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_kernel_sill_blocked

   subroutine test_kernel_masking(error)
      !! `masking_depth` gate: a face whose MEAN along-face height is at or
      !! above the gate is left fully open regardless of its sill, so the
      !! parameterization never narrows shelf faces the grid resolves.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)

      bed = -100.0_wp
      h_layer = 25.0_wp
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp

      ! Gate height -80 m: d_avg = -75 is ABOVE it, so nothing is narrowed.
      call run_update(NXT, NYT, POROUS_ETA_MAX, -80.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)
      checks: block
         call check(error, all(por_u == 1.0_wp) .and. all(por_v == 1.0_wp), &
                    "a face shallower than masking_depth must stay fully open")
         if (allocated(error)) exit checks

         ! Gate height -70 m: d_avg = -75 is BELOW it, so the sill applies.
         call run_update(NXT, NYT, POROUS_ETA_MAX, -70.0_wp, bed, h_layer, &
                         dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)
         call check(error, all(por_u(2:NXT, :, 1) == 0.25_wp), &
                    "a face deeper than masking_depth must be narrowed")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_kernel_masking

   subroutine test_stats_resolved_flat(error)
      !! The resolved-bathymetry source on a sloping bed must bracket the
      !! along-face samples (`d_min <= d_avg <= d_max`) and collapse to a
      !! single value on a flat bed.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      integer :: i, j

      checks: block
         bed = -100.0_wp
         call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call check(error, all(dmin_u(2:NXT, :) == -100.0_wp) .and. &
                    all(dmax_u(2:NXT, :) == -100.0_wp) .and. &
                    all(davg_u(2:NXT, :) == -100.0_wp), &
                    "flat bed must give a degenerate (single-valued) face statistic")
         if (allocated(error)) exit checks

         ! A cross-face ridge: bed varies with j, so the along-face samples
         ! of a u-face genuinely differ.
         do j = 1, NYT
            do i = 1, NXT
               bed(i, j) = -100.0_wp + 10.0_wp*real(j, wp)
            end do
         end do
         call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call check(error, all(dmin_u(2:NXT, 2:NYT - 1) <= davg_u(2:NXT, 2:NYT - 1)) .and. &
                    all(davg_u(2:NXT, 2:NYT - 1) <= dmax_u(2:NXT, 2:NYT - 1)), &
                    "resolved statistics must satisfy d_min <= d_avg <= d_max")
         if (allocated(error)) exit checks
         call check(error, any(dmax_u(2:NXT, 2:NYT - 1) > dmin_u(2:NXT, 2:NYT - 1)), &
                    "a cross-face ridge must produce a non-degenerate u-face sill")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_resolved_flat

   ! =================================================================
   ! Conservation — the seam into continuity
   ! =================================================================

   subroutine test_continuity_conservation(error)
      !! Narrowing face widths redistributes mass but must not create or
      !! destroy any: with closed walls the flux divergence telescopes, so
      !! `sum(h_layer)` (uniform Cartesian => uniform areaT) is invariant.
      !! A seam that narrowed the OUTGOING face of a cell but not the
      !! matching INCOMING face of its neighbour would break exactly here.
      !!
      !! Paired with a non-vacuity check: the porous run must actually
      !! DIFFER from the un-narrowed run, otherwise conservation proves
      !! nothing.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8, N_STEPS = 40
      real(wp), parameter :: DT = 0.05_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: h_plain(:, :, :)
      real(wp) :: mass0, mass1, drift, spread
      integer :: i, j, k, nx, ny, step

      checks: block
         ! ---- Reference run: porous OFF ----
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         do step = 1, N_STEPS
            call continuity_compute_fluxes(grid, metrics, ct, ms)
            call continuity_apply_fluxes(ms, DT)
         end do
         call unmap_state(ms, ct)
         allocate (h_plain(nx, ny, NZ))
         h_plain = ms%h_layer
         call ct%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         ! ---- Porous ON, same IC ----
         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         mass0 = sum(ms%h_layer)
         call map_state(ms, ct)
         do step = 1, N_STEPS
            call continuity_compute_fluxes(grid, metrics, ct, ms)
            call continuity_apply_fluxes(ms, DT)
         end do
         call unmap_state(ms, ct)
         mass1 = sum(ms%h_layer)

         drift = abs(mass1 - mass0)/mass0
         call check(error, drift < 1.0e-12_wp, &
                    "porous barriers must not create or destroy mass")
         if (allocated(error)) exit checks

         spread = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  spread = max(spread, abs(ms%h_layer(i, j, k) - h_plain(i, j, k)))
               end do
            end do
         end do
         call check(error, spread > 1.0e-6_wp, &
                    "the porous run must differ from the un-narrowed run "// &
                    "(otherwise the conservation check is vacuous)")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(h_plain)) deallocate (h_plain)
      call ct%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_continuity_conservation

   subroutine seed_state(ms, nx, ny)
      !! Non-trivial thickness + a smooth wall-vanishing velocity field, so
      !! the transport is genuinely active in the interior and exactly zero
      !! at the closed walls.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny
      integer :: i, j, k
      real(wp) :: sx, sy

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 25.0_wp + 0.5_wp*real(k, wp) + &
                                     2.0_wp*sin(real(i, wp)*0.4_wp)*cos(real(j, wp)*0.3_wp)
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         do j = 2, ny - 1
            do i = 2, nx
               sx = sin(3.14159265358979_wp*real(i - 1, wp)/real(nx - 1, wp))
               ms%u_face_x_layer(i, j, k) = 0.3_wp*sx*sin(real(j, wp)*0.5_wp)
            end do
         end do
         do j = 2, ny
            do i = 2, nx - 1
               sy = sin(3.14159265358979_wp*real(j - 1, wp)/real(ny - 1, wp))
               ms%v_face_y_layer(i, j, k) = 0.2_wp*sy*cos(real(i, wp)*0.4_wp)
            end do
         end do
      end do
   end subroutine seed_state

   pure function all_wet(nx, ny) result(w)
      !! An all-wet `wet_T` mask for the fill-statistics helper.  Most
      !! cases have no land; `test_stats_resolved_ignores_land_corners`
      !! builds its own.
      integer, intent(in) :: nx, ny
      real(wp) :: w(nx, ny)
      w = 1.0_wp
   end function all_wet

   subroutine build_metrics(metrics, grid, porous, blocked)
      !! Uniform-Cartesian metrics, optionally with porous barriers armed
      !! on a mid-depth sill.  The porous arrays are grown BEFORE
      !! `enter_data`, which is the ordering the production configure path
      !! also obeys (a realloc after the map would leave the device
      !! pointing at freed host memory).
      !!
      !! `blocked` swaps the mid-depth sill (fractions 0.25/0.75/1/1) for
      !! a high one that CLOSES the two bed layers outright
      !! (0/0/0.4/1) — a `por = 0` layer is the only configuration in
      !! which a dropped porous factor in the renormaliser's per-layer
      !! increment shows up as a hard zero rather than a redistribution.
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: porous
      logical, intent(in), optional :: blocked

      logical :: use_blocked

      use_blocked = .false.
      if (present(blocked)) use_blocked = blocked

      call metrics%init(grid)
      call metrics_fill_cartesian(metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(metrics)
      if (porous) then
         call metrics_porous_alloc(metrics, grid, NZ)
         metrics%porous_eta_interp = POROUS_ETA_MAX
         metrics%porous_mask_depth = 0.0_wp
         metrics%por_bed = -100.0_wp
         if (use_blocked) then
            metrics%por_dmin_u = -40.0_wp
            metrics%por_dmax_u = -30.0_wp
            metrics%por_davg_u = -35.0_wp
            metrics%por_dmin_v = -40.0_wp
            metrics%por_dmax_v = -30.0_wp
            metrics%por_davg_v = -35.0_wp
         else
            metrics%por_dmin_u = -100.0_wp
            metrics%por_dmax_u = -50.0_wp
            metrics%por_davg_u = -75.0_wp
            metrics%por_dmin_v = -100.0_wp
            metrics%por_dmax_v = -50.0_wp
            metrics%por_davg_v = -75.0_wp
         end if
         metrics%use_porous = .true.
      end if
      !$acc enter data copyin(metrics)
      call metrics%enter_data()
      if (porous) then
         ! Fill the fractions once on the device from the seeded state's
         ! layer thicknesses; the production path does this per outer step
         ! via `ocean_porous_refresh`.
         call metrics_refresh(metrics, grid)
      end if
   end subroutine build_metrics

   subroutine metrics_refresh(metrics, grid)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), allocatable :: h_uniform(:, :, :)
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (h_uniform(nx, ny, NZ), source=25.0_wp)
      !$acc enter data copyin(h_uniform)
      call porous_update_face_areas(nx, ny, NZ, metrics%porous_eta_interp, &
                                    metrics%porous_mask_depth, metrics%por_bed, &
                                    h_uniform, &
                                    metrics%por_dmin_u, metrics%por_dmax_u, &
                                    metrics%por_davg_u, &
                                    metrics%por_dmin_v, metrics%por_dmax_v, &
                                    metrics%por_davg_v, &
                                    metrics%dy_cu, metrics%dx_cv, &
                                    metrics%por_face_area_u, metrics%por_face_area_v, &
                                    metrics%dy_cu_bt, metrics%dx_cv_bt)
      !$acc exit data delete(h_uniform)
      deallocate (h_uniform)
   end subroutine metrics_refresh

   subroutine teardown_metrics(metrics)
      type(ocean_metrics_t), intent(inout) :: metrics
      call metrics%exit_data()
      !$acc exit data delete(metrics)
      call metrics%destroy()
   end subroutine teardown_metrics

   subroutine map_state(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
      ! `enter_data` maps the mass fluxes `create` (the production step
      ! recomputes them on-device), so the device copies start as GARBAGE
      ! and any face the kernel under test does not write keeps it.  A
      ! face-by-face comparison of two runs would then differ on those
      ! faces for reasons that have nothing to do with the code under
      ! test, so push the host zeros over explicitly.
      ms%mass_flux_x_layer = 0.0_wp
      ms%mass_flux_y_layer = 0.0_wp
      !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
   end subroutine map_state

   subroutine unmap_state(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc update self(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine unmap_state

   ! =================================================================
   ! Device residency
   ! =================================================================

   subroutine test_bt_width_column(error)
      !! The barotropic transport width must carry the COLUMN-INTEGRATED
      !! open fraction, i.e. the thickness-weighted mean of the per-layer
      !! fractions.  Without it the barotropic solve is porous-blind and
      !! the per-layer renormalisation to `uhbt` hands the blocked
      !! transport straight back — the barrier would only redistribute
      !! transport in the vertical, never reduce it.
      !!
      !! Same configuration as `kernel_sill_partial_and_full`
      !! (fractions 0.25, 0.75, 1, 1 over four equal 25 m layers), so the
      !! thickness-weighted mean is exactly (0.25+0.75+1+1)/4 = 0.75.
      !! `run_update` uses unit widths, so `dy_cu_bt` IS that fraction.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)
      real(wp) :: want

      bed = -100.0_wp
      h_layer = 25.0_wp
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp
      call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      want = 0.75_wp
      checks: block
         call check(error, all(abs(last_bt_u(2:NXT, :) - want) < 1.0e-13_wp), &
                    "the BT u-width must carry the column-integrated open fraction")
         if (allocated(error)) exit checks
         call check(error, all(abs(last_bt_v(:, 2:NYT) - want) < 1.0e-13_wp), &
                    "the BT v-width must carry the column-integrated open fraction")
         if (allocated(error)) exit checks
         ! And it must equal the thickness-weighted mean of the per-layer
         ! fractions, computed here from the kernel's own output.
         call check(error, abs(last_bt_u(3, 3) - sum(por_u(3, 3, :))/real(NZ, wp)) &
                    < 1.0e-13_wp, &
                    "BT width must be the thickness-weighted mean of the layer fractions")
         if (allocated(error)) exit checks

         ! Flat bottom => degenerate stats => BT width untouched (exactly 1
         ! at unit widths).  This is the bit-identity property.
         call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                         dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)
         call check(error, all(last_bt_u == 1.0_wp) .and. all(last_bt_v == 1.0_wp), &
                    "a flat bottom must leave the BT widths EXACTLY un-narrowed")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_bt_width_column

   subroutine test_degenerate_sill(error)
      !! `a = (1-m)/m` divides by `m`, so a face whose mean height ties the
      !! deepest sample (`m = 0`) or the shallowest (`m = 1`) must be
      !! short-circuited, not evaluated.  Both collapse to a clean step.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -50.0_wp

      checks: block
         ! m = 0: step at d_min, everything above is open.
         call check(error, porous_open_width(DMIN, DMAX, DMIN, -75.0_wp) == 1.0_wp, &
                    "m=0 must be a step at d_min (fully open above), not a division")
         if (allocated(error)) exit checks
         call check(error, abs(porous_cum_area(DMIN, DMAX, DMIN, -75.0_wp) - 25.0_wp) &
                    < 1.0e-13_wp, "m=0 cumulative area must be linear from d_min")
         if (allocated(error)) exit checks
         ! Continuity with the above-d_max branch.
         call check(error, abs(porous_cum_area(DMIN, DMAX, DMIN, DMAX) - &
                               (DMAX - DMIN)) < 1.0e-13_wp, &
                    "m=0 must join the above-d_max branch continuously")
         if (allocated(error)) exit checks

         ! m = 1: step at d_max, nothing below it is open.
         call check(error, porous_open_width(DMIN, DMAX, DMAX, -75.0_wp) == 0.0_wp, &
                    "m=1 must be a step at d_max (closed below), not a division")
         if (allocated(error)) exit checks
         call check(error, porous_cum_area(DMIN, DMAX, DMAX, -75.0_wp) == 0.0_wp, &
                    "m=1 cumulative area must be exactly zero below d_max")
         if (allocated(error)) exit checks

         ! Fully degenerate face (d_max == d_min): the binary step.
         call check(error, porous_open_width(DMIN, DMIN, DMIN, -99.0_wp) == 1.0_wp, &
                    "a flat face must be fully open above the bed")
         if (allocated(error)) exit checks
         call check(error, porous_open_width(DMIN, DMIN, DMIN, -101.0_wp) == 0.0_wp, &
                    "a flat face must be fully closed below the bed")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_degenerate_sill

   subroutine test_renorm_split(error)
      !! Covers the DIRECTION-SPLIT continuity path and its barotropic
      !! renormalisation — the `wk = w*por(i,j,k)` conversion inside
      !! `renormalise_zonal_flux_to_uhbt`, which the fused-path
      !! conservation test never reaches.
      !!
      !! The renormaliser drives `sum_k mass_flux_x_layer = uhbt` at every
      !! face.  That constraint must hold WITH the narrowed face areas —
      !! if the porous factor were dropped from `sum_h` or from the
      !! per-layer increment, the corrected fluxes would miss the target.
      !! The porous run must also redistribute the transport differently
      !! from the un-narrowed one, or the check is vacuous.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8
      real(wp), parameter :: DT = 0.05_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: uhbt(:, :), flux_plain(:, :, :)
      real(wp) :: worst, worst_plain, spread, col
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (uhbt(nx + 1, ny), source=0.0_wp)
         do j = 2, ny - 1
            do i = 3, nx - 1
               uhbt(i, j) = 4.0_wp*sin(real(i, wp)*0.3_wp)*cos(real(j, wp)*0.2_wp)
            end do
         end do

         ! ---- un-narrowed reference ----
         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_zonal_flux(grid, metrics, ct, ms, DT, uhbt=uhbt)
         call unmap_state(ms, ct)
         allocate (flux_plain(nx + 1, ny, NZ))
         flux_plain = ms%mass_flux_x_layer
         worst_plain = 0.0_wp
         do j = 2, ny - 1
            do i = 3, nx - 1
               col = 0.0_wp
               do k = 1, NZ
                  col = col + flux_plain(i, j, k)
               end do
               worst_plain = max(worst_plain, abs(col - uhbt(i, j)))
            end do
         end do
         call ct%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         ! ---- porous ----
         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_zonal_flux(grid, metrics, ct, ms, DT, uhbt=uhbt)
         call unmap_state(ms, ct)

         worst = 0.0_wp
         spread = 0.0_wp
         do j = 2, ny - 1
            do i = 3, nx - 1
               col = 0.0_wp
               do k = 1, NZ
                  col = col + ms%mass_flux_x_layer(i, j, k)
                  spread = max(spread, abs(ms%mass_flux_x_layer(i, j, k) - &
                                           flux_plain(i, j, k)))
               end do
               worst = max(worst, abs(col - uhbt(i, j)))
            end do
         end do

         ! The oracle is the UN-NARROWED run's OWN residual, not an
         ! absolute tolerance: the renormaliser stops on a relative
         ! `RENORM_TOL` and is CFL-bracketed, so it never promises an
         ! exact hit (the plain residual here is already above 1e-8).
         ! What porous barriers must not do is make the constraint any
         ! worse — and dropping the factor from `sum_h` or from the
         ! per-layer increment would do exactly that.  Measured: the
         ! porous residual comes in at or below the plain one, so the
         ! bound below is tight, not decorative.
         call check(error, worst <= worst_plain + 1.0e-12_wp, &
                    "narrowed per-layer fluxes must still sum to uhbt as "// &
                    "tightly as the un-narrowed ones (the porous factor has "// &
                    "to be inside the renormalisation)")
         if (allocated(error)) exit checks
         call check(error, spread > 1.0e-6_wp, &
                    "the porous split-path run must differ from the un-narrowed one")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(uhbt)) deallocate (uhbt)
      if (allocated(flux_plain)) deallocate (flux_plain)
      call ct%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_renorm_split

   subroutine test_coriolis_narrowed(error)
      !! The Coriolis/advection TRANSPORT form must carry the same
      !! narrowed face width continuity does — two mass-flux definitions
      !! that disagree would inject spurious PV.  Checked exactly: the
      !! porous transport is the un-narrowed one times `por_face_area_u`.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: flux_plain(:, :, :), por(:, :, :)
      real(wp) :: worst, spread
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call cor%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor)
         call cor%enter_data()
         call coriolis_adv_compute_tendencies_hk(grid, metrics, cor, ms, &
                                                 ms%u_face_x_layer, &
                                                 ms%v_face_y_layer, ms%h_layer)
         !$acc update self(cor%mass_flux_u%data)
         call cor%exit_data()
         !$acc exit data delete(cor)
         call ms%exit_data()
         !$acc exit data delete(ms)
         allocate (flux_plain(nx + 1, ny, NZ))
         flux_plain = cor%mass_flux_u%data
         call cor%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call cor%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor)
         call cor%enter_data()
         call coriolis_adv_compute_tendencies_hk(grid, metrics, cor, ms, &
                                                 ms%u_face_x_layer, &
                                                 ms%v_face_y_layer, ms%h_layer)
         !$acc update self(cor%mass_flux_u%data)
         !$acc update self(metrics%por_face_area_u)
         call cor%exit_data()
         !$acc exit data delete(cor)
         call ms%exit_data()
         !$acc exit data delete(ms)

         allocate (por(nx + 1, ny, NZ))
         por = metrics%por_face_area_u
         worst = 0.0_wp
         spread = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  worst = max(worst, abs(cor%mass_flux_u%data(i, j, k) - &
                                         flux_plain(i, j, k)*por(i, j, k)))
                  spread = max(spread, abs(cor%mass_flux_u%data(i, j, k) - &
                                           flux_plain(i, j, k)))
               end do
            end do
         end do

         call check(error, worst < 1.0e-12_wp, &
                    "the Coriolis transport must be the un-narrowed one "// &
                    "times the open-area fraction")
         if (allocated(error)) exit checks
         call check(error, spread > 1.0e-6_wp, &
                    "the Coriolis transport must actually change (non-vacuity)")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(flux_plain)) deallocate (flux_plain)
      if (allocated(por)) deallocate (por)
      call cor%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_coriolis_narrowed

   subroutine test_device_resident(error)
      !! The recompute must read the DEVICE copies of the bed, the layer
      !! thicknesses and the sill statistics.
      !!
      !! For explicit-shape `do concurrent` dummies nvfortran emits an
      !! implicit `copyin(...) [if not already present]`, so an UNMAPPED
      !! buffer is silently staged from the host rather than faulting — a
      !! plain map-and-compare would pass whether or not the map existed.
      !! Here the HOST copies are overwritten with POISON after the map, so
      !! the two memories disagree and only a kernel bound to the resident
      !! device copies can reproduce the hand-computed fractions.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)
      real(wp) :: dy_cu(NXT + 1, NYT), dx_cv(NXT, NYT + 1)
      real(wp) :: dy_cu_bt(NXT + 1, NYT), dx_cv_bt(NXT, NYT + 1)
      real(wp) :: want(NZ)
      integer :: k

      ! Same configuration as `kernel_sill_partial`: fractions 0.25, 0.75, 1, 1.
      bed = -100.0_wp
      h_layer = 25.0_wp
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp
      want = [0.25_wp, 0.75_wp, 1.0_wp, 1.0_wp]

      checks: block
         ! The poison only discriminates if it is nowhere near the real data.
         call check(error, POISON < minval(bed) .and. POISON < minval(dmin_u), &
                    "test setup: POISON must sit below every real input")
         if (allocated(error)) exit checks

         dy_cu = 1.0_wp
         dx_cv = 1.0_wp
         !$acc enter data copyin(bed, h_layer, dmin_u, dmax_u, davg_u)
         !$acc enter data copyin(dmin_v, dmax_v, davg_v, dy_cu, dx_cv)
         !$acc enter data create(por_u, por_v, dy_cu_bt, dx_cv_bt)
#ifdef RDB_GPU_OFFLOAD
         call poison_2d(bed)
         call poison_3d(h_layer)
         call poison_2d(dmin_u)
         call poison_2d(dmax_u)
         call poison_2d(davg_u)
         call poison_2d(dy_cu)
#endif
         call porous_update_face_areas(NXT, NYT, NZ, POROUS_ETA_MAX, 0.0_wp, &
                                       bed, h_layer, dmin_u, dmax_u, davg_u, &
                                       dmin_v, dmax_v, davg_v, dy_cu, dx_cv, &
                                       por_u, por_v, dy_cu_bt, dx_cv_bt)
         !$acc update self(por_u, por_v, dy_cu_bt)
         !$acc exit data delete(bed, h_layer, dmin_u, dmax_u, davg_u)
         !$acc exit data delete(dmin_v, dmax_v, davg_v, dy_cu, dx_cv)
         !$acc exit data delete(por_u, por_v, dy_cu_bt, dx_cv_bt)

         do k = 1, NZ
            call check(error, all(por_u(2:NXT, :, k) == want(k)), &
                       "the open fractions must come from the mapped DEVICE "// &
                       "inputs, not the poisoned host copies")
            if (allocated(error)) exit checks
         end do
         ! `dy_cu` was poisoned too, so the BT width is a second, independent
         ! witness that the kernel bound to device memory.
         call check(error, all(abs(dy_cu_bt(2:NXT, :) - 0.75_wp) < 1.0e-13_wp), &
                    "the BT width must come from the mapped DEVICE dy_cu, "// &
                    "not the poisoned host copy")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_device_resident

   ! =================================================================
   ! Along-face statistics — weighting + gating
   ! =================================================================

   subroutine test_stats_simpson(error)
      !! The along-face mean must be the SIMPSON `(1,4,1)/6` combination
      !! of the three samples, not a plain `(1,1,1)/3` average.
      !!
      !! A LINEAR along-face profile cannot see the difference — there
      !! `s_lo + s_hi = 2*s_mid` and both weightings collapse to `s_mid`.
      !! The bed here is QUADRATIC in `j`, where the two forms differ by
      !! exactly `C/6` for `b(j) = -100 - C*j**2` (Simpson gives
      !! `b(j) - C/6`, the arithmetic mean `b(j) - C/3`), so the test both
      !! pins the right value and rejects the wrong one by a wide margin.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp), parameter :: C = 2.0_wp
      real(wp) :: bed(NXT, NYT)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: want_simpson, want_arith, bj
      integer :: i, j

      do j = 1, NYT
         do i = 1, NXT
            bed(i, j) = -100.0_wp - C*real(j, wp)**2
         end do
      end do
      call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                      dmin_u, dmax_u, davg_u, &
                                      dmin_v, dmax_v, davg_v)

      ! u-face (4,4): samples run in j, so the quadratic curvature is
      ! along the face and the weighting is visible.
      bj = bed(4, 4)
      want_simpson = bj - C/6.0_wp
      want_arith = bj - C/3.0_wp
      checks: block
         call check(error, abs(davg_u(4, 4) - want_simpson) < 1.0e-12_wp, &
                    "d_avg must be the Simpson (1,4,1)/6 along-face mean")
         if (allocated(error)) exit checks
         call check(error, abs(davg_u(4, 4) - want_arith) > 1.0e-3_wp, &
                    "test setup: Simpson and the plain average must differ here, "// &
                    "or the weighting is untested")
         if (allocated(error)) exit checks
         ! v-faces sample in i, where this bed is uniform, so they are the
         ! degenerate case and pin the midpoint fallback instead.
         call check(error, davg_v(4, 4) == 0.5_wp*(bed(4, 3) + bed(4, 4)), &
                    "a v-face over an along-face-uniform bed must be the midpoint")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_simpson

   subroutine test_stats_land_gating(error)
      !! A LAND cell in a corner stencil must not enter the along-face
      !! average.  Without the gate its elevation pulls `d_max` up and the
      !! fit blocks the deep layers of a face the grid resolves as fully
      !! open — measured as ~37% spurious blockage for one 4000 m land
      !! diagonal beside a 4000 m column.
      !!
      !! Non-vacuous by construction: the SAME bed with an all-wet mask
      !! must produce the non-degenerate (blocking) statistic, so the test
      !! fails both if the gate disappears and if it never fired.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), wet(NXT, NYT)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)

      bed = -4000.0_wp
      bed(5, 5) = 0.0_wp          ! land: elevation at the datum
      wet = 1.0_wp
      wet(5, 5) = 0.0_wp

      checks: block
         ! ---- all-wet mask: the land elevation DOES contaminate ----
         call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call check(error, dmax_u(5, 4) > dmin_u(5, 4), &
                    "test setup: un-gated, the land corner must make the face "// &
                    "statistic non-degenerate (otherwise nothing is being fixed)")
         if (allocated(error)) exit checks

         ! Every gated corner in turn — the land cell (5,5) sits in the
         ! NORTH corner stencil of u-face (5,4), the SOUTH corner stencil
         ! of u-face (5,6), the EAST corner stencil of v-face (4,5) and the
         ! WEST corner stencil of v-face (6,5).  All four gates must fire,
         ! so no single one can be removed unnoticed.
         call check(error, dmax_u(5, 4) > dmin_u(5, 4) .and. &
                    dmax_u(5, 6) > dmin_u(5, 6) .and. &
                    dmax_v(4, 5) > dmin_v(4, 5) .and. &
                    dmax_v(6, 5) > dmin_v(6, 5), &
                    "test setup: un-gated, the land corner must contaminate all "// &
                    "four neighbouring faces")
         if (allocated(error)) exit checks

         ! ---- wet-gated: the contaminated corner is dropped ----
         call porous_fill_stats_resolved(NXT, NYT, bed, wet, &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call check(error, dmax_u(5, 4) == dmin_u(5, 4) .and. &
                    davg_u(5, 4) == dmin_u(5, 4), &
                    "a wet-wet u-face with a land diagonal to the NORTH must fall "// &
                    "back to the two-cell midpoint, not average the land in")
         if (allocated(error)) exit checks
         call check(error, dmax_u(5, 6) == dmin_u(5, 6) .and. &
                    davg_u(5, 6) == dmin_u(5, 6), &
                    "same for a land diagonal to the SOUTH of a u-face")
         if (allocated(error)) exit checks
         call check(error, dmax_v(4, 5) == dmin_v(4, 5) .and. &
                    davg_v(4, 5) == dmin_v(4, 5), &
                    "same for a land diagonal to the EAST of a v-face")
         if (allocated(error)) exit checks
         call check(error, dmax_v(6, 5) == dmin_v(6, 5) .and. &
                    davg_v(6, 5) == dmin_v(6, 5), &
                    "same for a land diagonal to the WEST of a v-face")
         if (allocated(error)) exit checks
         call check(error, dmax_u(5, 4) == -4000.0_wp, &
                    "the fallback must be the two-cell midpoint depth")
         if (allocated(error)) exit checks
         ! Faces well away from the land are untouched by the gate.
         call check(error, dmax_u(3, 4) == -4000.0_wp .and. dmin_u(3, 4) == -4000.0_wp, &
                    "a face far from land must be unaffected")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_land_gating

   subroutine test_stats_ordering(error)
      !! The reader-boundary invariant `d_min <= d_avg <= d_max`.  The
      !! resolved filler cannot break it; a file-backed source could, and
      !! `configure_ocean_porous` refuses to run on statistics that do.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N1 = 5, N2 = 4
      real(wp) :: dmin(N1, N2), dmax(N1, N2), davg(N1, N2)

      dmin = -100.0_wp
      dmax = -50.0_wp
      davg = -75.0_wp
      checks: block
         call check(error, porous_stats_are_ordered(N1, N2, dmin, dmax, davg), &
                    "well-ordered statistics must be accepted")
         if (allocated(error)) exit checks
         davg(3, 2) = -40.0_wp        ! above d_max
         call check(error,.not. porous_stats_are_ordered(N1, N2, dmin, dmax, davg), &
                    "d_avg above d_max must be rejected")
         if (allocated(error)) exit checks
         davg(3, 2) = -120.0_wp       ! below d_min
         call check(error,.not. porous_stats_are_ordered(N1, N2, dmin, dmax, davg), &
                    "d_avg below d_min must be rejected")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_stats_ordering

   ! =================================================================
   ! Degenerate / vanishing geometry
   ! =================================================================

   subroutine test_ridge_not_walled(error)
      !! A one-cell ridge with deep basins either side, running PARALLEL
      !! to the u-faces, gives every along-face sample the same value.
      !! The bare fit's limit there is a STEP at the two-cell mean height
      !! — a hard wall on the deep layers of a face the grid resolves as
      !! open, manufactured by the proxy rather than measured.  The kernel
      !! must recognise the degenerate statistic and block NOTHING.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)
      real(wp) :: mean_h
      integer :: i, j, k

      do j = 1, NYT
         do i = 1, NXT
            if (i == 5) then
               bed(i, j) = -500.0_wp        ! the ridge crest
            else
               bed(i, j) = -3000.0_wp       ! deep basin either side
            end if
         end do
      end do
      do k = 1, NZ
         do j = 1, NYT
            do i = 1, NXT
               h_layer(i, j, k) = -bed(i, j)/real(NZ, wp)
            end do
         end do
      end do

      call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                      dmin_u, dmax_u, davg_u, &
                                      dmin_v, dmax_v, davg_v)
      ! POROUS_ETA_MIN, not the MAX default: MAX takes the SHALLOWER of the
      ! two adjacent interfaces, which on this ridge sits above the sampled
      ! sill top and leaves the face open whatever the fit says — so a MAX
      ! run cannot tell a degenerate-face fallback from a no-op.  MIN is
      ! the blocking rule and puts the manufactured wall on display.
      call run_update(NXT, NYT, POROUS_ETA_MIN, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      mean_h = 0.5_wp*(-3000.0_wp - 500.0_wp)
      checks: block
         ! The statistic really is degenerate at the ridge faces.
         call check(error, dmax_u(5, 4) == dmin_u(5, 4) .and. &
                    dmax_u(5, 4) == mean_h, &
                    "test setup: the ridge face statistic must collapse to the "// &
                    "two-cell mean height")
         if (allocated(error)) exit checks
         ! The bare fit at that degenerate statistic IS a wall — this is
         ! what the kernel must refuse to apply.
         call check(error, porous_open_width(mean_h, mean_h, mean_h, -2000.0_wp) &
                    == 0.0_wp, &
                    "test setup: the degenerate fit walls below the mean height, "// &
                    "so the kernel fallback is doing real work")
         if (allocated(error)) exit checks

         call check(error, all(por_u == 1.0_wp), &
                    "a u-face whose along-face statistic is degenerate must be "// &
                    "left FULLY OPEN, not walled at the two-cell mean depth")
         if (allocated(error)) exit checks
         call check(error, all(last_bt_u == 1.0_wp), &
                    "the BT u-width must be un-narrowed on a degenerate face")
         if (allocated(error)) exit checks
         ! The v-faces CROSS the ridge, so their along-face samples (which
         ! run in i) genuinely differ and the scheme narrows them.  That is
         ! the intended behaviour, and asserting it keeps the u-face check
         ! above from being a "porous barriers do nothing here" tautology.
         call check(error, any(por_v < 1.0_wp), &
                    "v-faces, whose samples run ACROSS the ridge, must still be "// &
                    "narrowed — the degenerate fallback must not disable the "// &
                    "scheme wholesale")
         if (allocated(error)) exit checks

         ! ---- Same ridge rotated 90 degrees ----
         ! Now the ridge runs in i, so it is the V-faces whose samples are
         ! along-face uniform.  Without this half the v-face fallback would
         ! be dead code no test reaches.
         do j = 1, NYT
            do i = 1, NXT
               if (j == 5) then
                  bed(i, j) = -500.0_wp
               else
                  bed(i, j) = -3000.0_wp
               end if
            end do
         end do
         do k = 1, NZ
            do j = 1, NYT
               do i = 1, NXT
                  h_layer(i, j, k) = -bed(i, j)/real(NZ, wp)
               end do
            end do
         end do
         call porous_fill_stats_resolved(NXT, NYT, bed, all_wet(NXT, NYT), &
                                         dmin_u, dmax_u, davg_u, &
                                         dmin_v, dmax_v, davg_v)
         call run_update(NXT, NYT, POROUS_ETA_MIN, 0.0_wp, bed, h_layer, &
                         dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)
         call check(error, dmax_v(4, 5) == dmin_v(4, 5) .and. &
                    dmax_v(4, 5) == mean_h, &
                    "test setup: the rotated ridge must make the V-face statistic "// &
                    "degenerate at the two-cell mean height")
         if (allocated(error)) exit checks
         call check(error, all(por_v == 1.0_wp), &
                    "a v-face whose along-face statistic is degenerate must be "// &
                    "left FULLY OPEN, not walled at the two-cell mean depth")
         if (allocated(error)) exit checks
         call check(error, all(last_bt_v == 1.0_wp), &
                    "the BT v-width must be un-narrowed on a degenerate face")
         if (allocated(error)) exit checks
         call check(error, any(por_u < 1.0_wp), &
                    "u-faces crossing the rotated ridge must still be narrowed")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_ridge_not_walled

   subroutine test_vanished_column(error)
      !! A wholly VANISHED column blocks every layer (`por = 0`); the
      !! barotropic width must go to zero with them.  Leaving it
      !! un-narrowed would let the BT mode transport at FULL width across
      !! a face every layer has closed — the two modes disagreeing about
      !! the same face.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 8, NYT = 8
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, NZ)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, NZ), por_v(NXT, NYT + 1, NZ)

      bed = -100.0_wp
      h_layer = 0.0_wp             ! every layer vanished
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp
      call run_update(NXT, NYT, POROUS_ETA_MAX, 0.0_wp, bed, h_layer, &
                      dmin_u, dmax_u, davg_u, dmin_v, dmax_v, davg_v, por_u, por_v)

      checks: block
         call check(error, all(por_u(2:NXT, :, :) == 0.0_wp), &
                    "every layer of a vanished column must be fully blocked")
         if (allocated(error)) exit checks
         call check(error, all(last_bt_u(2:NXT, :) == 0.0_wp), &
                    "the BT u-width must be zero where every layer is blocked")
         if (allocated(error)) exit checks
         call check(error, all(last_bt_v(:, 2:NYT) == 0.0_wp), &
                    "the BT v-width must be zero where every layer is blocked")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_vanished_column

   subroutine test_nz1_column(error)
      !! `nz = 1`: the column recurrence has a single layer, so the layer
      !! fraction and the column fraction are the SAME number.  Exercised
      !! separately because every other kernel case here is built on the
      !! module's `NZ = 4`.
      !!
      !! Bed -100 m, one 100 m layer (interfaces -100 -> 0), sill
      !! -100/-50/-75: `A(0) - A(-100) = 75`, so both fractions are 0.75.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXT = 6, NYT = 6
      real(wp) :: bed(NXT, NYT), h_layer(NXT, NYT, 1)
      real(wp) :: dmin_u(NXT + 1, NYT), dmax_u(NXT + 1, NYT), davg_u(NXT + 1, NYT)
      real(wp) :: dmin_v(NXT, NYT + 1), dmax_v(NXT, NYT + 1), davg_v(NXT, NYT + 1)
      real(wp) :: por_u(NXT + 1, NYT, 1), por_v(NXT, NYT + 1, 1)
      real(wp) :: dy_cu(NXT + 1, NYT), dx_cv(NXT, NYT + 1)
      real(wp) :: dy_cu_bt(NXT + 1, NYT), dx_cv_bt(NXT, NYT + 1)

      bed = -100.0_wp
      h_layer = 100.0_wp
      dmin_u = -100.0_wp; dmax_u = -50.0_wp; davg_u = -75.0_wp
      dmin_v = -100.0_wp; dmax_v = -50.0_wp; davg_v = -75.0_wp
      dy_cu = 1.0_wp
      dx_cv = 1.0_wp

      !$acc enter data copyin(bed, h_layer, dmin_u, dmax_u, davg_u)
      !$acc enter data copyin(dmin_v, dmax_v, davg_v, dy_cu, dx_cv)
      !$acc enter data create(por_u, por_v, dy_cu_bt, dx_cv_bt)
      call porous_update_face_areas(NXT, NYT, 1, POROUS_ETA_MAX, 0.0_wp, &
                                    bed, h_layer, dmin_u, dmax_u, davg_u, &
                                    dmin_v, dmax_v, davg_v, dy_cu, dx_cv, &
                                    por_u, por_v, dy_cu_bt, dx_cv_bt)
      !$acc update self(por_u, por_v, dy_cu_bt, dx_cv_bt)
      !$acc exit data delete(bed, h_layer, dmin_u, dmax_u, davg_u)
      !$acc exit data delete(dmin_v, dmax_v, davg_v, dy_cu, dx_cv)
      !$acc exit data delete(por_u, por_v, dy_cu_bt, dx_cv_bt)

      checks: block
         call check(error, all(abs(por_u(2:NXT, :, 1) - 0.75_wp) < 1.0e-13_wp), &
                    "nz=1: the single layer fraction must be the hand-computed 0.75")
         if (allocated(error)) exit checks
         call check(error, all(abs(dy_cu_bt(2:NXT, :) - 0.75_wp) < 1.0e-13_wp), &
                    "nz=1: the BT width must equal the single layer fraction")
         if (allocated(error)) exit checks
         call check(error, all(abs(por_v(:, 2:NYT, 1) - 0.75_wp) < 1.0e-13_wp), &
                    "nz=1: the v-face fraction must match too")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_nz1_column

   subroutine test_eta_interp_monotone(error)
      !! `w` is monotone increasing in the interface height, so the rule
      !! that returns the HIGHER interface (MAX, the default) is the LEAST
      !! blocking and MIN the most.  Pinned because the docstrings used to
      !! claim the opposite.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DMIN = -100.0_wp, DMAX = -30.0_wp, DAVG = -65.0_wp
      real(wp), parameter :: ZA = -80.0_wp, ZB = -40.0_wp
      real(wp) :: w_max, w_min, w_arith

      w_max = porous_open_width(DMIN, DMAX, DAVG, porous_eta_face(ZA, ZB, POROUS_ETA_MAX))
      w_min = porous_open_width(DMIN, DMAX, DAVG, porous_eta_face(ZA, ZB, POROUS_ETA_MIN))
      w_arith = porous_open_width(DMIN, DMAX, DAVG, porous_eta_face(ZA, ZB, POROUS_ETA_ARITH))

      checks: block
         call check(error, w_max > w_arith .and. w_arith > w_min, &
                    "MAX must leave MORE of the face open than ARITHMETIC, which "// &
                    "must leave more open than MIN")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_eta_interp_monotone

   ! =================================================================
   ! Per-face narrowing — the seams a conservation check cannot see
   ! =================================================================

   subroutine test_split_narrows_both(error)
      !! The PRODUCTION direction-split path must narrow BOTH directions.
      !!
      !! A telescoping `sum(h_layer)` check cannot catch a missing
      !! narrowing: mass is conserved whatever per-face multiplier the
      !! transport carries, because the same flux array is added to one
      !! cell and subtracted from its neighbour.  The assertion has to be
      !! per FACE — the narrowed transport must be EXACTLY the un-narrowed
      !! one times the open fraction, face by face, layer by layer, in x
      !! and in y.
      !!
      !! Driven WITHOUT `uhbt`/`vhbt` so no renormalisation runs: this
      !! isolates the narrowing seam from the barotropic constraint.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8
      real(wp), parameter :: DT = 0.05_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fx_plain(:, :, :), fy_plain(:, :, :)
      real(wp), allocatable :: por_u(:, :, :), por_v(:, :, :)
      real(wp) :: worst_x, worst_y, spread_x, spread_y, mag_x, mag_y, expect
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         ! ---- un-narrowed reference ----
         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_zonal_flux(grid, metrics, ct, ms, DT)
         call continuity_meridional_flux(grid, metrics, ct, ms, DT)
         call unmap_state(ms, ct)
         allocate (fx_plain(nx + 1, ny, NZ))
         allocate (fy_plain(nx, ny + 1, NZ))
         fx_plain = ms%mass_flux_x_layer
         fy_plain = ms%mass_flux_y_layer
         call ct%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         ! ---- porous ----
         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_zonal_flux(grid, metrics, ct, ms, DT)
         call continuity_meridional_flux(grid, metrics, ct, ms, DT)
         !$acc update self(metrics%por_face_area_u, metrics%por_face_area_v)
         call unmap_state(ms, ct)
         allocate (por_u(nx + 1, ny, NZ))
         allocate (por_v(nx, ny + 1, NZ))
         por_u = metrics%por_face_area_u
         por_v = metrics%por_face_area_v

         worst_x = 0.0_wp; spread_x = 0.0_wp; mag_x = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ! `expect` is a rounded LOCAL on purpose: written inline,
                  ! nvfortran's -fast contracts `a - b*c` into an FMA and the
                  ! comparison then measures the (nonzero) rounding error of
                  ! `b*c` rather than the difference under test.
                  expect = fx_plain(i, j, k)*por_u(i, j, k)
                  worst_x = max(worst_x, abs(ms%mass_flux_x_layer(i, j, k) - expect))
                  mag_x = max(mag_x, abs(fx_plain(i, j, k)))
                  spread_x = max(spread_x, abs(ms%mass_flux_x_layer(i, j, k) - &
                                               fx_plain(i, j, k)))
               end do
            end do
         end do
         worst_y = 0.0_wp; spread_y = 0.0_wp; mag_y = 0.0_wp
         do k = 1, NZ
            do j = 1, ny + 1
               do i = 1, nx
                  expect = fy_plain(i, j, k)*por_v(i, j, k)
                  worst_y = max(worst_y, abs(ms%mass_flux_y_layer(i, j, k) - expect))
                  mag_y = max(mag_y, abs(fy_plain(i, j, k)))
                  spread_y = max(spread_y, abs(ms%mass_flux_y_layer(i, j, k) - &
                                               fy_plain(i, j, k)))
               end do
            end do
         end do

         call check(error, worst_x <= TOL_ULP*mag_x, &
                    "split-path ZONAL transport must be EXACTLY the un-narrowed "// &
                    "one times the u-face open fraction")
         if (allocated(error)) exit checks
         call check(error, worst_y <= TOL_ULP*mag_y, &
                    "split-path MERIDIONAL transport must be EXACTLY the "// &
                    "un-narrowed one times the v-face open fraction")
         if (allocated(error)) exit checks
         call check(error, spread_x > 1.0e-6_wp, &
                    "the zonal comparison must be non-vacuous")
         if (allocated(error)) exit checks
         call check(error, spread_y > 1.0e-6_wp, &
                    "the meridional comparison must be non-vacuous")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(fx_plain)) deallocate (fx_plain)
      if (allocated(fy_plain)) deallocate (fy_plain)
      if (allocated(por_u)) deallocate (por_u)
      if (allocated(por_v)) deallocate (por_v)
      call ct%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_split_narrows_both

   subroutine test_fused_narrows_both(error)
      !! The same per-face assertion for the FUSED `continuity_compute_fluxes`
      !! path — the unsplit oracle the conservation test drives.  Its
      !! `sum(h_layer)` telescopes off one shared flux array, so only a
      !! face-by-face comparison can see a direction whose narrowing was
      !! dropped.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fx_plain(:, :, :), fy_plain(:, :, :)
      real(wp), allocatable :: por_u(:, :, :), por_v(:, :, :)
      real(wp) :: worst_x, worst_y, spread_x, spread_y, mag_x, mag_y, expect
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_compute_fluxes(grid, metrics, ct, ms)
         call unmap_state(ms, ct)
         allocate (fx_plain(nx + 1, ny, NZ))
         allocate (fy_plain(nx, ny + 1, NZ))
         fx_plain = ms%mass_flux_x_layer
         fy_plain = ms%mass_flux_y_layer
         call ct%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_compute_fluxes(grid, metrics, ct, ms)
         !$acc update self(metrics%por_face_area_u, metrics%por_face_area_v)
         call unmap_state(ms, ct)
         allocate (por_u(nx + 1, ny, NZ))
         allocate (por_v(nx, ny + 1, NZ))
         por_u = metrics%por_face_area_u
         por_v = metrics%por_face_area_v

         worst_x = 0.0_wp; spread_x = 0.0_wp; mag_x = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ! `expect` is a rounded LOCAL on purpose: written inline,
                  ! nvfortran's -fast contracts `a - b*c` into an FMA and the
                  ! comparison then measures the (nonzero) rounding error of
                  ! `b*c` rather than the difference under test.
                  expect = fx_plain(i, j, k)*por_u(i, j, k)
                  worst_x = max(worst_x, abs(ms%mass_flux_x_layer(i, j, k) - expect))
                  mag_x = max(mag_x, abs(fx_plain(i, j, k)))
                  spread_x = max(spread_x, abs(ms%mass_flux_x_layer(i, j, k) - &
                                               fx_plain(i, j, k)))
               end do
            end do
         end do
         worst_y = 0.0_wp; spread_y = 0.0_wp; mag_y = 0.0_wp
         do k = 1, NZ
            do j = 1, ny + 1
               do i = 1, nx
                  expect = fy_plain(i, j, k)*por_v(i, j, k)
                  worst_y = max(worst_y, abs(ms%mass_flux_y_layer(i, j, k) - expect))
                  mag_y = max(mag_y, abs(fy_plain(i, j, k)))
                  spread_y = max(spread_y, abs(ms%mass_flux_y_layer(i, j, k) - &
                                               fy_plain(i, j, k)))
               end do
            end do
         end do

         call check(error, worst_x <= TOL_ULP*mag_x, &
                    "fused-path ZONAL transport must be EXACTLY the un-narrowed "// &
                    "one times the u-face open fraction")
         if (allocated(error)) exit checks
         call check(error, worst_y <= TOL_ULP*mag_y, &
                    "fused-path MERIDIONAL transport must be EXACTLY the "// &
                    "un-narrowed one times the v-face open fraction")
         if (allocated(error)) exit checks
         call check(error, spread_x > 1.0e-6_wp .and. spread_y > 1.0e-6_wp, &
                    "both fused-path comparisons must be non-vacuous")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(fx_plain)) deallocate (fx_plain)
      if (allocated(fy_plain)) deallocate (fy_plain)
      if (allocated(por_u)) deallocate (por_u)
      if (allocated(por_v)) deallocate (por_v)
      call ct%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_fused_narrows_both

   subroutine test_renorm_blocked_layer(error)
      !! The barotropic renormalisation must carry the open fraction in the
      !! PER-LAYER INCREMENT, not only in the `sum_h` denominator.
      !!
      !! `Σ_k flux_k = uhbt` cannot catch a dropped factor: the Newton loop
      !! is self-consistent, so whichever weight it uses it converges on the
      !! same total and only the DISTRIBUTION changes.  A layer the barrier
      !! has closed outright is the discriminator — its increment must be
      !! identically zero, so the flux stays at the zero the narrowing pass
      !! left.  With the factor dropped it acquires `du*h_face*dy_cu`.
      !!
      !! Run on BOTH renormaliser branches (Newton and the legacy
      !! single-step form), because each carries its own copy of the
      !! weight.
      type(error_type), allocatable, intent(out) :: error

      call renorm_blocked_case(error, legacy=.false.)
      if (allocated(error)) return
      call renorm_blocked_case(error, legacy=.true.)
   end subroutine test_renorm_blocked_layer

   subroutine renorm_blocked_case(error, legacy)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: legacy
      integer, parameter :: NXP = 12, NYP = 8
      real(wp), parameter :: DT = 0.05_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: uhbt(:, :)
      real(wp) :: worst_blocked, open_mag
      integer :: i, j, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (uhbt(nx + 1, ny), source=0.0_wp)
         do j = 2, ny - 1
            do i = 3, nx - 1
               uhbt(i, j) = 4.0_wp*sin(real(i, wp)*0.3_wp)*cos(real(j, wp)*0.2_wp)
            end do
         end do

         ! `blocked=.true.` puts the sill ABOVE the two bed layers, so
         ! por(:,:,1) = por(:,:,2) = 0 exactly.
         call build_metrics(metrics, grid, porous=.true., blocked=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         ct%renorm_legacy_single_step = legacy
         call seed_state(ms, nx, ny)
         call map_state(ms, ct)
         call continuity_zonal_flux(grid, metrics, ct, ms, DT, uhbt=uhbt)
         !$acc update self(metrics%por_face_area_u)
         call unmap_state(ms, ct)

         call check(error, all(metrics%por_face_area_u(3:nx - 1, 2:ny - 1, 1) == 0.0_wp), &
                    "test setup: the bed layer must be fully blocked")
         if (allocated(error)) exit checks

         worst_blocked = 0.0_wp
         open_mag = 0.0_wp
         do j = 2, ny - 1
            do i = 3, nx - 1
               worst_blocked = max(worst_blocked, abs(ms%mass_flux_x_layer(i, j, 1)))
               worst_blocked = max(worst_blocked, abs(ms%mass_flux_x_layer(i, j, 2)))
               open_mag = max(open_mag, abs(ms%mass_flux_x_layer(i, j, 4)))
            end do
         end do

         call check(error, worst_blocked == 0.0_wp, &
                    "a layer the barrier closes must stay at EXACTLY zero "// &
                    "transport through the renormalisation (the open fraction "// &
                    "has to be in the per-layer increment, not just in sum_h)")
         if (allocated(error)) exit checks
         call check(error, open_mag > 1.0e-6_wp, &
                    "the renormalisation must actually be moving transport in "// &
                    "the open layers (otherwise the zero above is vacuous)")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(uhbt)) deallocate (uhbt)
      call ct%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine renorm_blocked_case

   subroutine test_coriolis_energy_narrowed(error)
      !! `sadourny_energy` is the OTHER transport-form Coriolis scheme.  It
      !! builds its own `mass_flux_u/v` and must narrow them exactly as the
      !! HK form does — a scheme that skipped it would feed the PV flux a
      !! mass transport continuity does not agree with.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 12, NYP = 8
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fu_plain(:, :, :), fv_plain(:, :, :)
      real(wp), allocatable :: por_u(:, :, :), por_v(:, :, :)
      real(wp) :: worst_u, worst_v, spread, mag, expect
      integer :: i, j, k, nx, ny

      checks: block
         call grid%init(NXP, NYP, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         call build_metrics(metrics, grid, porous=.false.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call cor%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor)
         call cor%enter_data()
         ! Same `create`-mapped-garbage trap as `map_state`: the transport
         ! buffers are kernel products, so zero the device copies before a
         ! face-by-face comparison of two runs.
         cor%mass_flux_u%data = 0.0_wp
         cor%mass_flux_v%data = 0.0_wp
         !$acc update device(cor%mass_flux_u%data, cor%mass_flux_v%data)
         call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                              ms%u_face_x_layer, &
                                                              ms%v_face_y_layer, ms%h_layer)
         !$acc update self(cor%mass_flux_u%data, cor%mass_flux_v%data)
         call cor%exit_data()
         !$acc exit data delete(cor)
         call ms%exit_data()
         !$acc exit data delete(ms)
         allocate (fu_plain(nx + 1, ny, NZ))
         allocate (fv_plain(nx, ny + 1, NZ))
         fu_plain = cor%mass_flux_u%data
         fv_plain = cor%mass_flux_v%data
         call cor%destroy()
         call ms%destroy()
         call teardown_metrics(metrics)

         call build_metrics(metrics, grid, porous=.true.)
         ms%nz_ml = NZ
         call ms%init(grid)
         call cor%init(grid, nz_ml=NZ)
         call seed_state(ms, nx, ny)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor)
         call cor%enter_data()
         ! Same `create`-mapped-garbage trap as `map_state`: the transport
         ! buffers are kernel products, so zero the device copies before a
         ! face-by-face comparison of two runs.
         cor%mass_flux_u%data = 0.0_wp
         cor%mass_flux_v%data = 0.0_wp
         !$acc update device(cor%mass_flux_u%data, cor%mass_flux_v%data)
         call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                              ms%u_face_x_layer, &
                                                              ms%v_face_y_layer, ms%h_layer)
         !$acc update self(cor%mass_flux_u%data, cor%mass_flux_v%data)
         !$acc update self(metrics%por_face_area_u, metrics%por_face_area_v)
         call cor%exit_data()
         !$acc exit data delete(cor)
         call ms%exit_data()
         !$acc exit data delete(ms)

         allocate (por_u(nx + 1, ny, NZ))
         allocate (por_v(nx, ny + 1, NZ))
         por_u = metrics%por_face_area_u
         por_v = metrics%por_face_area_v

         worst_u = 0.0_wp; worst_v = 0.0_wp; spread = 0.0_wp; mag = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  expect = fu_plain(i, j, k)*por_u(i, j, k)
                  worst_u = max(worst_u, abs(cor%mass_flux_u%data(i, j, k) - expect))
                  mag = max(mag, abs(fu_plain(i, j, k)))
                  spread = max(spread, abs(cor%mass_flux_u%data(i, j, k) - &
                                           fu_plain(i, j, k)))
               end do
            end do
         end do
         do k = 1, NZ
            do j = 1, ny + 1
               do i = 1, nx
                  expect = fv_plain(i, j, k)*por_v(i, j, k)
                  worst_v = max(worst_v, abs(cor%mass_flux_v%data(i, j, k) - expect))
                  mag = max(mag, abs(fv_plain(i, j, k)))
                  spread = max(spread, abs(cor%mass_flux_v%data(i, j, k) - &
                                           fv_plain(i, j, k)))
               end do
            end do
         end do

         call check(error, worst_u <= TOL_ULP*mag .and. worst_v <= TOL_ULP*mag, &
                    "sadourny_energy transports must be EXACTLY the un-narrowed "// &
                    "ones times the open-area fractions, in both directions")
         if (allocated(error)) exit checks
         call check(error, spread > 1.0e-6_wp, &
                    "the sadourny_energy comparison must be non-vacuous")
         if (allocated(error)) exit checks
      end block checks

      if (allocated(fu_plain)) deallocate (fu_plain)
      if (allocated(fv_plain)) deallocate (fv_plain)
      if (allocated(por_u)) deallocate (por_u)
      if (allocated(por_v)) deallocate (por_v)
      call cor%destroy()
      call ms%destroy()
      call teardown_metrics(metrics)
   end subroutine test_coriolis_energy_narrowed

   ! =================================================================
   ! The barotropic solve must CONSUME the narrowed widths
   ! =================================================================

   subroutine test_bt_substep_consumes(error)
      !! The kernel half of the porous-aware barotropic fix is guarded by
      !! `bt_width_carries_column_fraction`; this is the CONSUMING half.
      !! Reverting the call site to `metrics%dy_cu` / `dx_cv` leaves every
      !! other test green, so the assertion has to be made on the solve
      !! itself:
      !!
      !!   run A: dy_cu = W,   dy_cu_bt = W/2
      !!   run B: dy_cu = W,   dy_cu_bt = W
      !!   run C: dy_cu = W/2, dy_cu_bt = W/2
      !!
      !! A must DIFFER from B (the BT widths are read at all) and be
      !! BITWISE EQUAL to C (the slow-path widths are not).  A call site
      !! wired to `dy_cu` flips both.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: eta_a, eta_b, eta_c

      call bt_width_probe(error, bt_scale=0.5_wp, slow_scale=1.0_wp, eta_probe=eta_a)
      if (allocated(error)) return
      call bt_width_probe(error, bt_scale=1.0_wp, slow_scale=1.0_wp, eta_probe=eta_b)
      if (allocated(error)) return
      call bt_width_probe(error, bt_scale=0.5_wp, slow_scale=0.5_wp, eta_probe=eta_c)
      if (allocated(error)) return

      checks: block
         call check(error, eta_a /= eta_b, &
                    "the barotropic solve must CONSUME dy_cu_bt/dx_cv_bt "// &
                    "(halving them changed nothing)")
         if (allocated(error)) exit checks
         call check(error, eta_a == eta_c, &
                    "the barotropic solve must NOT consume the slow-path "// &
                    "dy_cu/dx_cv (changing them changed the answer)")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_bt_substep_consumes

   subroutine bt_width_probe(error, bt_scale, slow_scale, eta_probe)
      !! One short nonlinear barotropic substep run off a Gaussian SSH
      !! bump, with the BT transport widths and the slow-path widths
      !! scaled independently.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: bt_scale, slow_scale
      real(wp), intent(out) :: eta_probe
      integer, parameter :: NXP = 16, NYP = 12
      real(wp), parameter :: DXY = 1000.0_wp, H_REF = 500.0_wp
      real(wp), parameter :: DT_INNER = 2.0_wp
      integer, parameter :: N_STEPS = 4
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      integer :: i, j, nx, ny, ip, jp

      call grid%init(NXP, NYP, NGHOST, DXY, DXY)
      nx = grid%nx_total
      ny = grid%ny_total
      call build_metrics(metrics, grid, porous=.false.)
      ! Scale the two width families independently, then push both to the
      ! device (they were mapped `copyin` by `metrics%enter_data`).
      metrics%dy_cu_bt = bt_scale*metrics%dy_cu
      metrics%dx_cv_bt = bt_scale*metrics%dx_cv
      metrics%dy_cu = slow_scale*metrics%dy_cu
      metrics%dx_cv = slow_scale*metrics%dx_cv
      !$acc update device(metrics%dy_cu_bt, metrics%dx_cv_bt)
      !$acc update device(metrics%dy_cu, metrics%dx_cv)

      call dyn%init(grid, nz_ml=NZ)
      call cor%init(grid, nz_ml=NZ)
      cor%f_corner = 0.0_wp
      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp
      dyn%bt_work%bt_rem_u = 1.0_wp
      dyn%bt_work%bt_rem_v = 1.0_wp
      dyn%bt_work%F_bt_u_fast = 0.0_wp
      dyn%bt_work%F_bt_v_fast = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            dyn%bt_work%bt_eta(i, j) = 0.5_wp* &
                                       exp(-((real(i, wp) - 0.5_wp*real(nx, wp))**2 + &
                                             (real(j, wp) - 0.5_wp*real(ny, wp))**2)/9.0_wp)
         end do
      end do

      !$acc enter data copyin(dyn, cor)
      call dyn%enter_data()
      call cor%enter_data()
      !$acc update device(dyn%bt_work%bt_eta, dyn%bt_work%bt_H_ref)
      !$acc update device(dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt)
      !$acc update device(dyn%bt_work%bt_rem_u, dyn%bt_work%bt_rem_v)
      !$acc update device(dyn%bt_work%F_bt_u_fast, dyn%bt_work%F_bt_v_fast)
      call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                 cor%f_corner, N_STEPS, DT_INNER)
      !$acc update self(dyn%bt_work%bt_eta_end)
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor)

      ip = grid%nghost + NXP/2 + 2
      jp = grid%nghost + NYP/2
      eta_probe = dyn%bt_work%bt_eta_end(ip, jp)
      call check(error, eta_probe == eta_probe, "bt_eta_end probe is NaN")

      call dyn%destroy()
      call cor%destroy()
      call teardown_metrics(metrics)
   end subroutine bt_width_probe

   ! =================================================================
   ! configure_ocean_porous — the sign conventions
   ! =================================================================

   subroutine test_configure_signs(error)
      !! `configure_ocean_porous` performs the two sign flips the whole
      !! scheme depends on, and both can be inverted without any kernel
      !! test noticing:
      !!   * `porous_mask_depth = -masking_depth` (a namelist DEPTH,
      !!     positive down, becomes a kernel gate HEIGHT, positive up),
      !!   * `por_bed = -barotropic%b` (`b` is a reference column DEPTH,
      !!     positive down; the fit works in heights).
      !! Flipping either turns a deep-sill gate into a shallow one, or the
      !! seafloor into a mountain range 2H above the surface.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 8, NYP = 6
      real(wp), parameter :: DXY = 1000.0_wp, MASK_DEPTH = 500.0_wp
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      integer :: i, j, nx, ny

      checks: block
         cfg%sim_type = "ocean"
         cfg%nx = NXP; cfg%ny = NYP; cfg%dx = DXY; cfg%dy = DXY
         cfg%nz_layers = NZ; cfg%nghost = NGHOST
         cfg%ocean%porous%enable = .true.
         cfg%ocean%porous%source = "resolved"
         cfg%ocean%porous%eta_interp = "max"
         ! Deliberately non-zero and non-symmetric: a sign flip on zero is
         ! invisible.
         cfg%ocean%porous%masking_depth = MASK_DEPTH

         call grid%init(NXP, NYP, NGHOST, DXY, DXY)
         nx = grid%nx_total
         ny = grid%ny_total
         state%multilayer%nz_ml = NZ
         call state%init(grid)
         ! Reference column DEPTH, positive DOWN, with a ramp so the sign
         ! test is not degenerate.
         do j = 1, ny
            do i = 1, nx
               state%barotropic%b(i, j) = 3000.0_wp + 10.0_wp*real(i, wp)
            end do
         end do

         call configure_ocean_porous(cfg, state, grid, compute_rank=0)

         call check(error, state%metrics%use_porous, &
                    "configure_ocean_porous did not arm the master switch")
         if (allocated(error)) exit checks
         call check(error, state%metrics%porous_mask_depth == -MASK_DEPTH, &
                    "masking_depth must be negated into a gate HEIGHT "// &
                    "(positive-down depth -> positive-up height)")
         if (allocated(error)) exit checks
         call check(error, all(state%metrics%por_bed == -state%barotropic%b), &
                    "por_bed must be the NEGATED reference depth")
         if (allocated(error)) exit checks
         call check(error, maxval(state%metrics%por_bed) < 0.0_wp, &
                    "every bed HEIGHT must be below the sea surface (negative)")
         if (allocated(error)) exit checks
         ! Full-size arrays, not the (1,1)/(1,1,1) placeholders.
         call check(error, size(state%metrics%por_face_area_u, 1) == nx + 1 .and. &
                    size(state%metrics%por_face_area_u, 3) == NZ, &
                    "the open-fraction arrays must be grown to full face size")
         if (allocated(error)) exit checks
         ! And the statistics must inherit the height convention.
         call check(error, maxval(state%metrics%por_davg_u(2:nx, 2:ny - 1)) < 0.0_wp, &
                    "the along-face statistics must be heights, not depths")
         if (allocated(error)) exit checks
      end block checks
      call state%destroy()
   end subroutine test_configure_signs

   subroutine poison_2d(buf)
      !! Overwrite the HOST copy of an already-mapped array.  Plain
      !! sequential loops — never `do concurrent`, never array syntax — so
      !! the store cannot itself be offloaded onto the copy under test.
      ! assumed-shape-ok: host-only helper, no `do concurrent` here.
      real(wp), intent(inout) :: buf(:, :)
      integer :: i, j
      do j = 1, size(buf, 2)
         do i = 1, size(buf, 1)
            buf(i, j) = POISON
         end do
      end do
   end subroutine poison_2d

   subroutine poison_3d(buf)
      ! assumed-shape-ok: host-only helper, no `do concurrent` here.
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: i, j, k
      do k = 1, size(buf, 3)
         do j = 1, size(buf, 2)
            do i = 1, size(buf, 1)
               buf(i, j, k) = POISON
            end do
         end do
      end do
   end subroutine poison_3d

end module test_ocean_porous
