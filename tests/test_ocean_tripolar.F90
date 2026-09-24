!! Integration tests for a RUNNING tripolar ocean configuration (M4c).
!! The pure fold seam operators are unit-tested in test_ocean_fold.F90;
!! here we drive the FULL split-RK2 dyn loop on a small Murray bipolar-cap
!! grid with the north-fold tag active and assert the seam stays
!! physically consistent under dynamics.
!!
!! Tests:
!!   T1 quiescent gate: resting uniform-(T,S) tripolar basin, a few full
!!      dyn steps -> max|u|,|v| < 1e-10 (tier-1.5 trap incl. the seam).
!!   T2 v-seam antisymmetry: seed a smooth flow, step, assert v on the
!!      fold row satisfies v(i) = -v(i') to roundoff after each step.
!!   T3 cross-seam coherence: a tracer anomaly adjacent to the seam stays
!!      single-valued + globally conserved (no doubling/ghosting).
!!   T4 bit-identity guard: a non-tripolar (cartesian) bc has north_fold
!!      = .false. so it takes none of the new branches.
!!   T5/T6 near-pole metrics + quiescent gate.
!!   T7/T8 cross-fold conservation (pred_corr / ssp_rk2): a bump + salt
!!      blob seeded ON the fold; mass, salt and the boundary outflux are
!!      conserved to round-off, the fold-line v and mass flux pair off
!!      exactly, and the north ghosts are exact fold images (finding B1).
!!   T9 rest over topography + a land island straddling the fold, with a
!!      flat density interface: stays at rest.
module test_ocean_tripolar
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, PI
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_tripolar, metrics_finalize, &
                                metrics_apply_land_mask
   use rdb_ocean_fold, only: fold_north_centre
   use pic_strings, only: to_string
   use ocean_test_metrics, only: make_tripolar_metrics, make_cartesian_metrics, &
                                 destroy_cartesian_metrics
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split, SPLIT_SCHEME_PRED_CORR, &
                            SPLIT_SCHEME_SSP_RK2
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_PERIODIC, &
                                       OBC_WALL, OBC_TRIPOLAR_FOLD, &
                                       ocean_bc_validate_periodic, &
                                       ocean_bc_validate_fold
   implicit none
   private

   public :: collect_ocean_tripolar_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NZ = 2
   integer, parameter :: NXP = 32, NYP = 20
   real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 47.0_wp
   real(wp), parameter :: DLON = 360.0_wp/real(NXP, wp)
   real(wp), parameter :: DLAT = 1.5_wp
   real(wp), parameter :: REARTH = 6.378e6_wp
   real(wp), parameter :: PHI_JOIN = 65.0_wp, LON_POLE = 0.0_wp
   real(wp), parameter :: H_TOTAL = 4000.0_wp

contains

   subroutine collect_ocean_tripolar_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("tripolar_quiescent_gate", test_quiescent), &
                  new_unittest("tripolar_v_seam_antisymmetry", test_v_antisym), &
                  new_unittest("tripolar_cross_seam_coherence", test_cross_seam), &
                  new_unittest("tripolar_bit_identity_guard", test_bit_identity), &
                  new_unittest("tripolar_near_pole_metrics_finite", test_near_pole_metrics), &
                  new_unittest("tripolar_near_pole_quiescent_gate", test_near_pole_quiescent), &
                  new_unittest("tripolar_cross_fold_conservation_pred_corr", &
                               test_cross_fold_conservation_pc), &
                  new_unittest("tripolar_cross_fold_conservation_ssp_rk2", &
                               test_cross_fold_conservation_ssp), &
                  new_unittest("tripolar_rest_topography", test_rest_topography) &
                  ]
   end subroutine collect_ocean_tripolar_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      ! dx/dy reinterpreted as dlon/dlat (degrees) for the tripolar fill.
      call grid%init(NXP, NYP, NGHOST, DLON, DLAT)
   end subroutine make_grid

   subroutine make_bc_fold(bc, grid)
      type(ocean_bc_state_t), intent(out) :: bc
      type(hgrid_t), intent(in) :: grid
      call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
      bc%west%bc_type = OBC_PERIODIC
      bc%east%bc_type = OBC_PERIODIC
      bc%north%bc_type = OBC_TRIPOLAR_FOLD
      call ocean_bc_validate_periodic(bc)
      call ocean_bc_validate_fold(bc)
   end subroutine make_bc_fold

   subroutine init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
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
   end subroutine init_all

   subroutine map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
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
   end subroutine map_in

   subroutine map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
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
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine destroy_all

   ! Seed a resting uniform-(T,S) basin: flat layers, zero velocity,
   ! uniform tracers, bt_H_ref = column total.  f from corner planetary
   ! is irrelevant at rest.
   subroutine seed_rest(grid, ms, dyn)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_dyn_t), intent(inout) :: dyn
      integer :: nx, ny
      nx = grid%nx_total
      ny = grid%ny_total
      ms%h_layer = H_TOTAL/real(NZ, wp)
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*(H_TOTAL/real(NZ, wp))
      ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*(H_TOTAL/real(NZ, wp))
      dyn%bt_work%bt_H_ref = H_TOTAL
      dyn%bt_work%g_bt = 9.81_wp
   end subroutine seed_rest

   ! -----------------------------------------------------------------
   ! T1: quiescent gate
   ! -----------------------------------------------------------------
   subroutine test_quiescent(error)
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
      type(ocean_bc_state_t) :: bc
      integer :: step
      real(wp) :: maxu, maxv

      checks: block
         call make_grid(grid)
         call make_bc_fold(bc, grid)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                    REARTH, PHI_JOIN, LON_POLE)
         call seed_rest(grid, ms, dyn)
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         do step = 1, 8
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 600.0_wp, 20, bc=bc)
         end do

         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         maxu = maxval(abs(ms%u_face_x_layer))
         maxv = maxval(abs(ms%v_face_y_layer))

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         call check(error, maxu < 1.0e-10_wp, "T1: max|u| should be < 1e-10 at rest")
         if (allocated(error)) exit checks
         call check(error, maxv < 1.0e-10_wp, "T1: max|v| should be < 1e-10 at rest")
      end block checks
   end subroutine test_quiescent

   ! -----------------------------------------------------------------
   ! T2: v-seam antisymmetry preserved under dynamics
   ! -----------------------------------------------------------------
   subroutine test_v_antisym(error)
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
      type(ocean_bc_state_t) :: bc
      integer :: step, i, k, ng, ni, j_fold, isum, ip
      real(wp) :: asym, val

      checks: block
         call make_grid(grid)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                    REARTH, PHI_JOIN, LON_POLE)
         call seed_rest(grid, ms, dyn)
         ! Seed a smooth zonal flow so the dynamics are non-trivial.
         ng = grid%nghost
         ni = grid%nx_phys
         block
            integer :: ii, jj, kk
            real(wp) :: yy
            do kk = 1, NZ
               do jj = 1, grid%ny_total
                  yy = real(jj - ng, wp)/real(grid%ny_phys, wp)
                  do ii = 1, grid%nx_total + 1
                     ms%u_face_x_layer(ii, jj, kk) = 0.05_wp*sin(3.14159265_wp*yy)
                  end do
               end do
            end do
         end block
         call make_bc_fold(bc, grid)
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         ! Fold line = north face of the last T-row = v storage row
         ! ng+ny_phys+1 (v is the SOUTH face of T(i,j)).  This test used to
         ! check row ng+ny_phys — the last row's south face, an ordinary
         ! interior face that the old (MOM6 north-face) fold antisymmetrised
         ! by mistake — so it passed while the real seam face leaked mass.
         j_fold = ng + grid%ny_phys + 1
         isum = 2*ng + ni + 1

         do step = 1, 4
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 600.0_wp, 20, bc=bc)
         end do

         !$acc update self(ms%v_face_y_layer)
         ! After the stage-entry fold of the NEXT step would re-project; but
         ! the post-continuity fold inside the LAST step already left the
         ! fold row antisymmetric.  Check v(i,j_fold) = -v(i',j_fold).
         asym = 0.0_wp
         do k = 1, NZ
            do i = ng + 1, ng + ni
               ip = isum - i
               if (ip == i) then
                  val = abs(ms%v_face_y_layer(i, j_fold, k))
               else
                  val = abs(ms%v_face_y_layer(i, j_fold, k) + &
                            ms%v_face_y_layer(ip, j_fold, k))
               end if
               asym = max(asym, val)
            end do
         end do

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         call check(error, asym < 1.0e-12_wp, &
                    "T2: v on the fold row must be antisymmetric to roundoff")
      end block checks
   end subroutine test_v_antisym

   ! -----------------------------------------------------------------
   ! T3: cross-seam tracer coherence + global conservation
   ! -----------------------------------------------------------------
   subroutine test_cross_seam(error)
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
      type(ocean_bc_state_t) :: bc
      integer :: step, i, j, k, ng, idxS
      real(wp) :: salt0, salt1, smin, smax, rel

      checks: block
         call make_grid(grid)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                    REARTH, PHI_JOIN, LON_POLE)
         call seed_rest(grid, ms, dyn)
         idxS = ms%idx_salinity
         ng = grid%nghost
         ! Place a salinity anomaly in the interior rows adjacent to the seam,
         ! and seed a meridional flow that crosses the fold.
         block
            integer :: ii, jj, kk
            do kk = 1, NZ
               do jj = ng + grid%ny_phys - 2, ng + grid%ny_phys
                  do ii = ng + 1, ng + grid%nx_phys
                     ms%tracers(idxS)%hTr(ii, jj, kk) = 36.5_wp*(H_TOTAL/real(NZ, wp))
                  end do
               end do
            end do
            ! Seed a fold-CONSISTENT meridional flow: v antisymmetric about
            ! the seam (v(i') = -v(i)) so the cross-seam transport balances
            ! column-for-column.  A flow that crosses the fold but respects
            ! its antisymmetry is the physically-admissible cross-seam field;
            ! a uniform (non-antisymmetric) v would inject a spurious seam
            ! convergence that is a property of the IC, not the scheme.
            block
               integer :: ii, jj, kk, iip
               real(wp) :: vmag
               do kk = 1, NZ
                  do jj = 1, grid%ny_total + 1
                     do ii = ng + 1, ng + grid%nx_phys
                        iip = 2*ng + grid%nx_phys + 1 - ii
                        ! West half flows north, east-conjugate flows south.
                        if (ii <= iip) then
                           vmag = 0.03_wp
                        else
                           vmag = -0.03_wp
                        end if
                        ms%v_face_y_layer(ii, jj, kk) = vmag
                     end do
                  end do
               end do
            end block
         end block
         call make_bc_fold(bc, grid)
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         salt0 = interior_sum(ms%tracers(idxS)%hTr, grid)

         do step = 1, 6
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 300.0_wp, 20, bc=bc)
         end do

         !$acc update self(ms%tracers(idxS)%hTr)
         salt1 = interior_sum(ms%tracers(idxS)%hTr, grid)
         ! No anomalous extrema growth (no seam doubling): interior salinity
         ! concentration must stay within the seeded [35, 36.5] band.
         smin = huge(1.0_wp)
         smax = -huge(1.0_wp)
         do k = 1, NZ
            do j = ng + 1, ng + grid%ny_phys
               do i = ng + 1, ng + grid%nx_phys
                  smin = min(smin, ms%tracers(idxS)%hTr(i, j, k))
                  smax = max(smax, ms%tracers(idxS)%hTr(i, j, k))
               end do
            end do
         end do

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         rel = abs(salt1 - salt0)/abs(salt0)
         ! Global cross-seam conservation: the fold halo-fill + antisymmetric
         ! v projection make the cross-seam transport balance column-for-
         ! column, so the interior integral is conserved to the documented
         ! hTr/h roundoff-amplification floor (~1e-8 here) — orders of
         ! magnitude below any doubling/ghosting artefact (which would be
         ! O(1)).
         call check(error, rel < 1.0e-6_wp, &
                    "T3: total interior salt must be conserved across the seam")
         if (allocated(error)) exit checks
         ! No anomalous extrema growth (the doubling/ghosting failure mode
         ! produces a ~2x concentration spike at the seam).  The PPM scheme
         ! is not strictly monotone, so allow a small relative band slack.
         call check(error, smax <= 36.5_wp*(H_TOTAL/real(NZ, wp))*(1.0_wp + 1.0e-4_wp), &
                    "T3: no anomalous salinity maximum growth at the seam")
         if (allocated(error)) exit checks
         call check(error, smin >= 34.5_wp*(H_TOTAL/real(NZ, wp)), &
                    "T3: no anomalous salinity minimum at the seam")
      end block checks
   end subroutine test_cross_seam

   pure function interior_sum(fld, grid) result(s)
      real(wp), intent(in) :: fld(:, :, :)
      type(hgrid_t), intent(in) :: grid
      real(wp) :: s
      integer :: i, j, k, ng
      ng = grid%nghost
      s = 0.0_wp
      do k = 1, size(fld, 3)
         do j = ng + 1, ng + grid%ny_phys
            do i = ng + 1, ng + grid%nx_phys
               s = s + fld(i, j, k)
            end do
         end do
      end do
   end function interior_sum

   ! -----------------------------------------------------------------
   ! T4: bit-identity guard — non-tripolar bc never folds
   ! -----------------------------------------------------------------
   subroutine test_bit_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc_cart, bc_periodic

      checks: block
         call grid%init(16, 12, NGHOST, 1.0_wp, 1.0_wp)
         ! A plain closed-wall cartesian bc: north_fold must be .false.
         call ocean_bc_state_init(bc_cart, grid, nz_ml=NZ, n_tracers=2)
         call ocean_bc_validate_fold(bc_cart)
         call check(error,.not. bc_cart%north_fold, &
                    "T4: cartesian/wall bc must have north_fold = .false.")
         call ocean_bc_state_destroy(bc_cart)
         if (allocated(error)) exit checks

         ! A periodic-x bc with a closed north wall: still no fold.
         call ocean_bc_state_init(bc_periodic, grid, nz_ml=NZ, n_tracers=2)
         bc_periodic%west%bc_type = OBC_PERIODIC
         bc_periodic%east%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc_periodic)
         call ocean_bc_validate_fold(bc_periodic)
         call check(error,.not. bc_periodic%north_fold, &
                    "T4: periodic-x + wall-north bc must have north_fold = .false.")
         call ocean_bc_state_destroy(bc_periodic)
      end block checks
   end subroutine test_bit_identity

   ! -----------------------------------------------------------------
   ! T5: near-pole metrics — finite, positive, actually shrinking
   ! -----------------------------------------------------------------
   ! The M4c suite above (T1-T4) uses LAT_S=47, NYP=20, DLAT=1.5 => tops
   ! out at 77N where cells are still ~37 km — it tests the fold
   ! MECHANICS but never probes the small-cell regime where a real global
   ! run fails (FINDINGS.md section 4: "This is the whole argument for
   ! the end-to-end global test bed: the fold mechanics are unit-tested,
   ! the polar metrics are not.").  T5/T6 close that gap: a grid whose
   ! physical rows actually reach into the bipolar cap's near-pole
   ! neighbourhood.
   subroutine test_near_pole_metrics(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP2 = 48, NYP2 = 40, NGHOST2 = 3
      real(wp), parameter :: LAT_S2 = 47.0_wp, DLAT2 = 1.5_wp
      real(wp), parameter :: DLON2 = 360.0_wp/real(NXP2, wp)
      real(wp), parameter :: PHI_JOIN2 = 65.0_wp, LON_POLE2 = 0.0_wp
      real(wp), parameter :: REARTH2 = 6.378e6_wp
      integer :: i, j, ng, i0, i1, j0, j1
      logical :: all_finite, all_positive
      real(wp) :: dx_min, dx_at_south

      checks: block
         call grid%init(NXP2, NYP2, NGHOST2, DLON2, DLAT2)
         call make_tripolar_metrics(metrics, grid, 0.0_wp, LAT_S2, DLON2, DLAT2, &
                                    REARTH2, PHI_JOIN2, LON_POLE2)
         ! Physical rows span LAT_S2=47N through PHI_JOIN2=65N (ordinary
         ! lon-lat) and on into the Murray bipolar cap — the top row's
         ! bipolar y-coordinate reaches the immediate neighbourhood of
         ! the two cap poles (lon_pole, lon_pole+180).  NYP2=40 rows
         ! (vs. the M4c suite's 20) is sized to actually land in the
         ! pathological small-cell regime a global run hits immediately.
         ng = grid%nghost
         i0 = ng + 1; i1 = ng + grid%nx_phys
         j0 = ng + 1; j1 = ng + grid%ny_phys

         all_finite = .true.
         all_positive = .true.
         dx_min = huge(1.0_wp)
         do j = j0, j1
            do i = i0, i1
               if (.not. (ieee_is_finite(metrics%dxT(i, j)) .and. &
                          ieee_is_finite(metrics%dyT(i, j)) .and. &
                          ieee_is_finite(metrics%areaT(i, j)))) all_finite = .false.
               if (.not. (metrics%dxT(i, j) > 0.0_wp .and. metrics%dyT(i, j) > 0.0_wp &
                          .and. metrics%areaT(i, j) > 0.0_wp)) all_positive = .false.
               dx_min = min(dx_min, metrics%dxT(i, j), metrics%dyT(i, j))
            end do
         end do

         call destroy_cartesian_metrics(metrics)

         call check(error, all_finite, &
                    "T5: near-pole T-cell metrics (dxT/dyT/areaT) must be finite")
         if (allocated(error)) exit checks
         call check(error, all_positive, &
                    "T5: near-pole T-cell metrics (dxT/dyT/areaT) must be strictly positive")
         if (allocated(error)) exit checks
         ! Regression net for FINDINGS.md section 3/4: the cap must
         ! actually shrink cells toward the pole, not hold them at the
         ! coarse equatorward lon-lat spacing — dx_min must be well
         ! below the nominal cell size at the domain's southern
         ! (ordinary lon-lat) edge.
         dx_at_south = DLON2*(PI/180.0_wp)*REARTH2*cos(LAT_S2*PI/180.0_wp)
         call check(error, dx_min < dx_at_south, &
                    "T5: cells must shrink toward the cap, not stay at the "// &
                    "equatorward lon-lat spacing (dx_min reported via the "// &
                    "check message on failure)")
      end block checks
   end subroutine test_near_pole_metrics

   ! -----------------------------------------------------------------
   ! T6: near-pole quiescent gate — a resting basin including the cap
   ! poles themselves must stay finite and at rest.  Weak in the same
   ! sense T1 is weak (uniform density + zero velocity keeps every
   ! tendency exactly zero regardless of the fold/cap being right) —
   ! but it is NOT vacuous for a geometry singularity: a degenerate
   ! metric (zero area / corner-length) at the exact cap pole would
   ! still produce 0/0 in a reciprocal-metric kernel even with zero
   ! velocity, independent of any forcing.  The forced/dynamical
   ! near-pole failure mode FINDINGS.md section 7 leaves open (does an
   ! aquaplanet's missing land at the cap poles — vs. the fold itself —
   ! explain the day 4-8 failure) is NOT settled by this test; it needs
   ! real forcing + long integration, out of scope for a unit test.
   subroutine test_near_pole_quiescent(error)
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
      type(ocean_bc_state_t) :: bc
      integer, parameter :: NXP2 = 48, NYP2 = 40, NGHOST2 = 3
      real(wp), parameter :: LAT_S2 = 47.0_wp, DLAT2 = 1.5_wp
      real(wp), parameter :: DLON2 = 360.0_wp/real(NXP2, wp)
      real(wp), parameter :: PHI_JOIN2 = 65.0_wp, LON_POLE2 = 0.0_wp
      real(wp), parameter :: REARTH2 = 6.378e6_wp
      integer :: step
      real(wp) :: maxu, maxv

      checks: block
         call grid%init(NXP2, NYP2, NGHOST2, DLON2, DLAT2)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         bc%north%bc_type = OBC_TRIPOLAR_FOLD
         call ocean_bc_validate_periodic(bc)
         call ocean_bc_validate_fold(bc)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call make_tripolar_metrics(metrics, grid, 0.0_wp, LAT_S2, DLON2, DLAT2, &
                                    REARTH2, PHI_JOIN2, LON_POLE2)
         call seed_rest(grid, ms, dyn)
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         do step = 1, 8
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 600.0_wp, 20, bc=bc)
         end do

         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         maxu = maxval(abs(ms%u_face_x_layer))
         maxv = maxval(abs(ms%v_face_y_layer))

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         call check(error, ieee_is_finite(maxu) .and. ieee_is_finite(maxv), &
                    "T6: near-pole quiescent basin must stay FINITE (no NaN at the cap)")
         if (allocated(error)) exit checks
         call check(error, maxu < 1.0e-8_wp, &
                    "T6: near-pole quiescent basin max|u| should stay ~0")
         if (allocated(error)) exit checks
         call check(error, maxv < 1.0e-8_wp, &
                    "T6: near-pole quiescent basin max|v| should stay ~0")
      end block checks
   end subroutine test_near_pole_quiescent

   ! -----------------------------------------------------------------
   ! T7/T8: CROSS-FOLD CONSERVATION (global-tripolar finding B1 gate).
   ! -----------------------------------------------------------------
   ! A surface-height bump and a salt blob are seeded ON the fold (the last
   ! physical T-row, x = NXP/4), so the gravity wave and the advected blob
   ! cross the north seam every step.  Over NSTEP steps:
   !   * total mass  Σ h·areaT      and salt Σ hTr·areaT  are conserved to
   !     round-off (a closed periodic+fold basin has no boundary);
   !   * the accumulated boundary outflux `ms%mass_out` stays at round-off
   !     (the fold line telescopes).  Measured on this case, 40 steps:
   !     mass / salt / outflux drift 2e-15 / 3e-15 / 1e-21 relative; with
   !     the pre-fix (MOM6 north-face) indexing all three are 6e-6 — the
   !     B1 leak through the mis-placed seam face;
   !   * after every step the fold-line v row is EXACTLY antisymmetric and
   !     the last stage's meridional mass flux on it pairs off exactly
   !     (without the flux projection in `continuity_tracer_step_split` the
   !     two slots differ by ~1e-12 relative), while carrying a real
   !     cross-fold transport (see the f_max check);
   !   * the north ghost rows of h / S are the exact fold images (the IC is
   !     seeded seam-consistent, as the engine's init wrap leaves it).
   ! Both split schemes: the pred_corr time-means (u_av/v_av/h_av) are
   ! folded by their own seam fill, which ssp_rk2 never reads.
   subroutine test_cross_fold_conservation_pc(error)
      type(error_type), allocatable, intent(out) :: error
      call run_cross_fold_conservation(error, SPLIT_SCHEME_PRED_CORR)
   end subroutine test_cross_fold_conservation_pc

   subroutine test_cross_fold_conservation_ssp(error)
      type(error_type), allocatable, intent(out) :: error
      call run_cross_fold_conservation(error, SPLIT_SCHEME_SSP_RK2)
   end subroutine test_cross_fold_conservation_ssp

   subroutine run_cross_fold_conservation(error, scheme)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
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
      integer, parameter :: NSTEP = 40
      integer :: step, i, j, k, ng, ni, nj, jf, isum, jsum, idxS, p, pm
      real(wp) :: mass0, mass1, salt0, salt1, rel_mass, rel_salt, rel_out
      real(wp) :: v_asym, f_asym, f_max, ghost_err, r2, x0
      real(wp) :: hl

      checks: block
         call make_grid(grid)
         call make_bc_fold(bc, grid)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         dyn%split_scheme = scheme
         call make_tripolar_metrics(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                    REARTH, PHI_JOIN, LON_POLE)
         call seed_rest(grid, ms, dyn)
         ng = grid%nghost
         ni = grid%nx_phys
         nj = grid%ny_phys
         jf = ng + nj + 1                  ! v / corner fold-line row
         isum = 2*ng + ni + 1
         jsum = 2*ng + 2*nj + 1
         idxS = ms%idx_salinity
         hl = H_TOTAL/real(NZ, wp)
         ! Bump + blob centred on the top physical row at x = ni/4 (the
         ! fold image of that point is x = 3ni/4 on the same row).  Seeded
         ! on every row up to the fold (periodic ghost columns from the
         ! wrapped coordinate), then the north ghost rows are folded -- a
         ! seam-consistent IC, exactly what the engine's init wrap leaves.
         x0 = real(ni, wp)/4.0_wp
         do j = 1, ng + nj
            do i = 1, grid%nx_total
               r2 = ((modulo(real(i - ng, wp) - 0.5_wp, real(ni, wp)) - x0)**2 + &
                     (real(j - ng, wp) - 0.5_wp - real(nj, wp))**2)/9.0_wp
               ms%h_layer(i, j, NZ) = hl + 5.0_wp*exp(-r2)
               do k = 1, NZ
                  ms%tracers(idxS)%hTr(i, j, k) = (35.0_wp + exp(-r2))*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 10.0_wp*ms%h_layer(i, j, k)
               end do
            end do
         end do
         call fold_north_centre(ms%h_layer, grid%nx_total, grid%ny_total, NZ, ni, nj, ng)
         call fold_north_centre(ms%tracers(idxS)%hTr, grid%nx_total, grid%ny_total, NZ, ni, nj, ng)
         call fold_north_centre(ms%tracers(ms%idx_temperature)%hTr, grid%nx_total, grid%ny_total, &
                                NZ, ni, nj, ng)
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         mass0 = weighted_interior_sum(ms%h_layer, metrics%areaT, grid)
         salt0 = weighted_interior_sum(ms%tracers(idxS)%hTr, metrics%areaT, grid)
         ms%mass_out = 0.0_wp

         v_asym = 0.0_wp
         f_asym = 0.0_wp
         f_max = 0.0_wp
         ghost_err = 0.0_wp
         do step = 1, NSTEP
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 600.0_wp, 20, bc=bc)
            !$acc update self(ms%v_face_y_layer, ms%mass_flux_y_layer, ms%h_layer)
            !$acc update self(ms%tracers(idxS)%hTr)
            do k = 1, NZ
               do i = ng + 1, ng + ni
                  p = i - ng
                  pm = ni + 1 - p
                  if (p >= pm) cycle
                  v_asym = max(v_asym, abs(ms%v_face_y_layer(i, jf, k) + &
                                           ms%v_face_y_layer(ng + pm, jf, k)))
                  f_asym = max(f_asym, abs(ms%mass_flux_y_layer(i, jf, k) + &
                                           ms%mass_flux_y_layer(ng + pm, jf, k)))
                  f_max = max(f_max, abs(ms%mass_flux_y_layer(i, jf, k)))
               end do
               do j = ng + nj + 1, grid%ny_total
                  do i = 1, grid%nx_total
                     ghost_err = max(ghost_err, &
                                     abs(ms%h_layer(i, j, k) - ms%h_layer(isum - i, jsum - j, k)), &
                                     abs(ms%tracers(idxS)%hTr(i, j, k) - &
                                         ms%tracers(idxS)%hTr(isum - i, jsum - j, k)))
                  end do
               end do
            end do
         end do

         mass1 = weighted_interior_sum(ms%h_layer, metrics%areaT, grid)
         salt1 = weighted_interior_sum(ms%tracers(idxS)%hTr, metrics%areaT, grid)
         rel_mass = abs(mass1 - mass0)/mass0
         rel_salt = abs(salt1 - salt0)/salt0
         ! mass_out is in kg (ρ_water·m³); normalise by the same total.
         rel_out = abs(ms%mass_out)/(mass0*1000.0_wp)

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         ! The bump drives ~8e7 m^3/s through the fold line.  A fold the
         ! barotropic fast loop treats as a north WALL (the BT dispatch used
         ! to zero the fold-line vbt) lets only ~3e5 m^3/s of baroclinic
         ! residual through, conserving mass perfectly while blocking the
         ! flow -- so conservation alone cannot certify the seam.
         call check(error, f_max > 1.0e7_wp, &
                    "T7: the bump must drive a real cross-fold transport, max|F_fold| = "// &
                    to_string(f_max))
         if (allocated(error)) exit checks
         call check(error, v_asym == 0.0_wp, &
                    "T7: fold-line v must be exactly antisymmetric, max|v(i)+v(i')| = "// &
                    to_string(v_asym))
         if (allocated(error)) exit checks
         call check(error, f_asym == 0.0_wp, &
                    "T7: fold-line mass flux must pair off exactly, max|F(i)+F(i')| = "// &
                    to_string(f_asym))
         if (allocated(error)) exit checks
         call check(error, ghost_err == 0.0_wp, &
                    "T7: north ghost rows must be the exact fold image, err = "// &
                    to_string(ghost_err))
         if (allocated(error)) exit checks
         call check(error, rel_mass < 1.0e-13_wp, &
                    "T7: total mass must be conserved across the fold, rel = "//to_string(rel_mass))
         if (allocated(error)) exit checks
         call check(error, rel_out < 1.0e-13_wp, &
                    "T7: boundary outflux of a closed fold basin must be round-off, rel = "// &
                    to_string(rel_out))
         if (allocated(error)) exit checks
         call check(error, rel_salt < 1.0e-13_wp, &
                    "T7: total salt must be conserved across the fold, rel = "//to_string(rel_salt))
      end block checks
   end subroutine run_cross_fold_conservation

   ! -----------------------------------------------------------------
   ! T9: rest over topography on the tripolar grid.
   ! -----------------------------------------------------------------
   ! Two layers, lighter water on top (T = 15 over T = 5), a FLAT interface
   ! at 500 m, and a bed that varies everywhere: a meridional ridge plus a
   ! seamount straddling the fold line at x = 3ni/4, and a land island
   ! straddling the fold at x = ni/4 (so its fold image — the same island —
   ! closes the fold-line faces on both sides).  With a flat interface and
   ! flat free surface every pressure gradient vanishes, so a correct fold
   ! (metric ghosts, masks, the seam treatment of the BT fast loop) keeps the
   ! basin at rest.  An off-by-one fold pairs the fold-line faces with the
   ! wrong columns' metrics / masks and drives flow through the island.
   subroutine test_rest_topography(error)
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
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: wet(:, :), depth(:, :), h_top0(:, :)
      integer :: step, i, j, ng, ni, nj, nxt, nyt
      real(wp) :: x, y, maxu, maxv, dh_top
      real(wp), parameter :: H_TOP = 500.0_wp

      checks: block
         call make_grid(grid)
         call make_bc_fold(bc, grid)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ng = grid%nghost
         ni = grid%nx_phys
         nj = grid%ny_phys
         nxt = grid%nx_total
         nyt = grid%ny_total
         allocate (wet(nxt, nyt), depth(nxt, nyt), h_top0(nxt, nyt))
         ! Bed + island on EVERY storage cell from its (wrapped) coordinate;
         ! the ghosts are then made exact fold images below.
         do j = 1, nyt
            do i = 1, nxt
               x = modulo(real(i - ng, wp) - 0.5_wp, real(ni, wp))
               y = real(j - ng, wp) - 0.5_wp
               depth(i, j) = 3000.0_wp + 800.0_wp*cos(2.0_wp*PI*x/real(ni, wp)) &
                             - 600.0_wp*real(j - ng, wp)/real(nj, wp) &
                             - 1200.0_wp*exp(-((x - 0.75_wp*ni)**2 + (y - nj)**2)/4.0_wp)
               wet(i, j) = 1.0_wp
               if (abs(x - 0.25_wp*ni) < 1.6_wp .and. y > nj - 2.0_wp) wet(i, j) = 0.0_wp
            end do
         end do
         call fold_north_centre(depth, nxt, nyt, ni, nj, ng)
         call fold_north_centre(wet, nxt, nyt, ni, nj, ng)
         ! Metrics with the island masked (host edit BEFORE the device map).
         call metrics%init(grid)
         call metrics_fill_tripolar(metrics, grid, LON_W, LAT_S, DLON, DLAT, &
                                    REARTH, PHI_JOIN, LON_POLE)
         call metrics_finalize(metrics)
         call metrics_apply_land_mask(metrics, wet, grid, periodic_x=.true., &
                                      periodic_y=.false., north_fold=.true.)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()

         call seed_rest(grid, ms, dyn)
         do j = 1, nyt
            do i = 1, nxt
               ms%h_layer(i, j, NZ) = H_TOP
               ms%h_layer(i, j, 1) = depth(i, j) - H_TOP
               ms%tracers(ms%idx_temperature)%hTr(i, j, NZ) = 15.0_wp*H_TOP
               ms%tracers(ms%idx_temperature)%hTr(i, j, 1) = 5.0_wp*ms%h_layer(i, j, 1)
               ms%tracers(ms%idx_salinity)%hTr(i, j, :) = 35.0_wp*ms%h_layer(i, j, :)
               dyn%bt_work%bt_H_ref(i, j) = depth(i, j)
               h_top0(i, j) = H_TOP
            end do
         end do
         call map_in(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         do step = 1, 16
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, &
                                      bd, ss, va, hd, vd, vmix, ms, 600.0_wp, 20, bc=bc)
         end do

         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
         maxu = maxval(abs(ms%u_face_x_layer(ng + 1:ng + ni + 1, ng + 1:ng + nj, :)))
         maxv = maxval(abs(ms%v_face_y_layer(ng + 1:ng + ni, ng + 1:ng + nj + 1, :)))
         dh_top = maxval(abs(ms%h_layer(ng + 1:ng + ni, ng + 1:ng + nj, NZ) - &
                             h_top0(ng + 1:ng + ni, ng + 1:ng + nj)))

         call map_out(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call destroy_cartesian_metrics(metrics)
         call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call ocean_bc_state_destroy(bc)

         call check(error, ieee_is_finite(maxu) .and. ieee_is_finite(maxv), &
                    "T9: rest over topography must stay finite")
         if (allocated(error)) exit checks
         call check(error, maxu < 1.0e-10_wp .and. maxv < 1.0e-10_wp, &
                    "T9: rest over topography must stay at rest, max|u|,|v| = "// &
                    to_string(maxu)//", "//to_string(maxv))
         if (allocated(error)) exit checks
         call check(error, dh_top < 1.0e-8_wp, &
                    "T9: the flat interface must stay flat, max|dh_top| = "//to_string(dh_top))
      end block checks
   end subroutine test_rest_topography

   function weighted_interior_sum(fld, area, grid) result(s)
      !! Σ_interior Σ_k fld·area (host copy).
      real(wp), intent(in) :: fld(:, :, :), area(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp) :: s
      integer :: i, j, k, ng
      ng = grid%nghost
      s = 0.0_wp
      do k = 1, size(fld, 3)
         do j = ng + 1, ng + grid%ny_phys
            do i = ng + 1, ng + grid%nx_phys
               s = s + fld(i, j, k)*area(i, j)
            end do
         end do
      end do
   end function weighted_interior_sum

end module test_ocean_tripolar
