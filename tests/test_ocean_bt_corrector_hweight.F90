!! Unit tests for the h-weighted BT corrector
!! (rdb_barotropic_coupling::apply_bt_correction with
!! `use_h_weighted = .true.`).
!!
!! h-weighted distributes the BT-mode Δu by `h_face(k) / ⟨h⟩_h`
!! per-layer instead of uniformly.  Thin layers get a small share;
!! thick layers a correspondingly larger one; depth-mean preserved.
!!
!! Tests:
!!   1. uniform h → bit-identical to uniform Δu
!!   2. stratified h, Δu_bar > 0 → bed gets less, surf more
!!   3. depth-mean preservation (the mass-flux constraint)
module test_ocean_bt_corrector_hweight
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: apply_bt_correction
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_is_finite
   implicit none
   private

   public :: collect_bt_corrector_hweight_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_bt_corrector_hweight_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("uniform_h_matches_uniform_delta", test_uniform_h_match), &
                  new_unittest("stratified_h_thin_gets_less", test_stratified_h), &
                  new_unittest("depth_mean_preserved", test_depth_mean), &
                  new_unittest("visc_rem_unity_matches_h_only", test_vr_unity_matches), &
                  new_unittest("visc_rem_bed_damped_biases_against_bed", test_vr_bed_damped), &
                  new_unittest("skip_nonfinite_fold_input", test_skip_nonfinite) &
                  ]
   end subroutine collect_bt_corrector_hweight_tests

   subroutine build_state(grid, ms, bt_work, nz_ml)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      integer, intent(in) :: nz_ml
      call grid%init(4, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_ml
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=nz_ml)
   end subroutine build_state

   subroutine cleanup(ms, bt_work)
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(inout) :: bt_work
      call bt_work%destroy()
      call ms%destroy()
   end subroutine cleanup

   subroutine test_uniform_h_match(error)
      !! Uniform h_layer → h-weighted = uniform Δu, bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt_a, bt_b
      real(wp) :: max_diff
      integer :: k
      checks: block
         call build_state(grid, ms_a, bt_a, 3)
         call build_state(grid, ms_b, bt_b, 3)

         ms_a%h_layer = 100.0_wp
         ms_b%h_layer = 100.0_wp
         ms_a%u_face_x_layer = 0.5_wp
         ms_b%u_face_x_layer = 0.5_wp
         ms_a%v_face_y_layer = 0.3_wp
         ms_b%v_face_y_layer = 0.3_wp

         bt_a%bt_ubt_end = 1.0_wp; bt_a%bt_vbt_end = 0.4_wp
         bt_a%ubt_at_n = 0.5_wp; bt_a%vbt_at_n = 0.3_wp
         bt_a%F_bt_u = 0.0_wp; bt_a%F_bt_v = 0.0_wp
         bt_a%bt_H_ref = 300.0_wp; bt_a%bt_eta_end = 0.0_wp

         bt_b%bt_ubt_end = 1.0_wp; bt_b%bt_vbt_end = 0.4_wp
         bt_b%ubt_at_n = 0.5_wp; bt_b%vbt_at_n = 0.3_wp
         bt_b%F_bt_u = 0.0_wp; bt_b%F_bt_v = 0.0_wp
         bt_b%bt_H_ref = 300.0_wp; bt_b%bt_eta_end = 0.0_wp

         call apply_bt_correction(bt_a, ms_a, 100.0_wp, skip_h_rescale=.true.)
         call apply_bt_correction(bt_b, ms_b, 100.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true.)

         max_diff = maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer))
         call check(error, max_diff < 1.0e-12_wp, &
                    "uniform h: h-weighted should equal uniform-Δu")
      end block checks
      call cleanup(ms_a, bt_a); call cleanup(ms_b, bt_b)
   end subroutine test_uniform_h_match

   subroutine test_stratified_h(error)
      !! Stratified h: h_bed thin (k=1), h_surf thick (k=3).
      !! With Δu_bar > 0, bed (thin) gets a SMALLER share than surf.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      real(wp) :: du_bed, du_surf, ratio
      integer :: i, j
      checks: block
         call build_state(grid, ms, bt, 3)
         ! Layer thicknesses: k=1 thin (10m), k=2 medium (100m), k=3 thick (200m)
         ms%h_layer(:, :, 1) = 10.0_wp
         ms%h_layer(:, :, 2) = 100.0_wp
         ms%h_layer(:, :, 3) = 200.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         bt%bt_ubt_end = 1.0_wp; bt%bt_vbt_end = 0.0_wp
         bt%ubt_at_n = 0.0_wp; bt%vbt_at_n = 0.0_wp
         bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 310.0_wp; bt%bt_eta_end = 0.0_wp

         call apply_bt_correction(bt, ms, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true.)

         ! Pick an interior u-face.
         i = NGHOST + 2; j = NGHOST + 2
         du_bed = ms%u_face_x_layer(i, j, 1)
         du_surf = ms%u_face_x_layer(i, j, 3)

         ! Sanity: both non-zero, same sign as Δu_bar (positive)
         call check(error, du_bed > 0.0_wp .and. du_surf > 0.0_wp, &
                    "stratified h: both layers should get positive Δu")
         if (allocated(error)) exit checks
         ! Bed (thin) MUST get less than surf (thick)
         call check(error, du_bed < du_surf, &
                    "stratified h: bed should get smaller Δu than surf")
         if (allocated(error)) exit checks
         ! Ratio of corrections matches ratio of thicknesses
         ratio = du_bed/du_surf
         ! Expected: h_bed/h_surf = 10/200 = 0.05
         call check(error, abs(ratio - 10.0_wp/200.0_wp) < 1.0e-12_wp, &
                    "stratified h: Δu_bed/Δu_surf should equal h_bed/h_surf")
      end block checks
      call cleanup(ms, bt)
   end subroutine test_stratified_h

   subroutine test_depth_mean(error)
      !! After h-weighted correction, the depth-mean of u_layer must
      !! match bt_ubt_end to roundoff.  This is the mass-flux constraint.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      real(wp) :: u_mean, h_total, err_max
      integer :: i, j, k, nz
      checks: block
         call build_state(grid, ms, bt, 4)
         nz = ms%nz_ml
         ms%h_layer(:, :, 1) = 5.0_wp     ! thin bed
         ms%h_layer(:, :, 2) = 30.0_wp
         ms%h_layer(:, :, 3) = 150.0_wp
         ms%h_layer(:, :, 4) = 800.0_wp   ! thick surf
         ms%u_face_x_layer = 0.1_wp        ! some background flow
         ms%v_face_y_layer = -0.05_wp

         bt%bt_ubt_end = 0.7_wp; bt%bt_vbt_end = -0.2_wp
         bt%ubt_at_n = 0.1_wp; bt%vbt_at_n = -0.05_wp
         bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 985.0_wp; bt%bt_eta_end = 0.0_wp

         call apply_bt_correction(bt, ms, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true.)

         ! Compute depth-mean at one interior u-face.
         err_max = 0.0_wp
         do j = NGHOST + 1, NGHOST + 4
            do i = NGHOST + 2, NGHOST + 4
               u_mean = 0.0_wp
               h_total = 0.0_wp
               do k = 1, nz
                  u_mean = u_mean + ms%u_face_x_layer(i, j, k)* &
                           0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
                  h_total = h_total + 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
               end do
               u_mean = u_mean/h_total
               err_max = max(err_max, abs(u_mean - bt%bt_ubt_end(i, j)))
            end do
         end do

         call check(error, err_max < 1.0e-12_wp, &
                    "h-weighted: depth-mean of u_layer should equal bt_ubt_end")
      end block checks
      call cleanup(ms, bt)
   end subroutine test_depth_mean

   ! -----------------------------------------------------------------
   ! visc_rem joint weight (MOM6 frhatu · visc_rem) — extends the
   ! h-weighted corrector with a per-layer viscous-damping bias.
   ! -----------------------------------------------------------------

   subroutine test_vr_unity_matches(error)
      !! With `visc_rem ≡ 1.0` (the init value here — this test hand-fills
      !! rather than running vdiff), `use_visc_rem=.true.` must give a
      !! bit-identical result to `use_visc_rem=.false.`.  This pins the
      !! §3.3 row-sum invariant from the PRODUCER side (PR-19,
      !! `rdb_ocean_vdiff.F90`): a no-flux-top / no-flux-bottom viscous
      !! operator cannot remove a uniform acceleration, so γ ≡ 1 exactly
      !! for ANY interior viscosity/thickness/dt as long as the vdiff
      !! operator carries no drag — see `test_ocean_visc_rem.F90`'s
      !! `no_drag_remnant_is_exactly_one` for the producer-side proof.
      !! Here we only need the CONSUMER half: vr≡1 must reduce the
      !! h·visc_rem joint weight to the plain h-weighted path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt_a, bt_b
      real(wp) :: max_diff
      integer :: k
      checks: block
         call build_state(grid, ms_a, bt_a, 3)
         call build_state(grid, ms_b, bt_b, 3)

         ! Stratified h so the h-weighted path actually does something
         do k = 1, 3
            ms_a%h_layer(:, :, k) = real(k*50, wp)
            ms_b%h_layer(:, :, k) = real(k*50, wp)
         end do
         ms_a%u_face_x_layer = 0.5_wp; ms_b%u_face_x_layer = 0.5_wp
         ms_a%v_face_y_layer = 0.3_wp; ms_b%v_face_y_layer = 0.3_wp

         bt_a%bt_ubt_end = 1.0_wp; bt_a%bt_vbt_end = 0.4_wp
         bt_a%ubt_at_n = 0.5_wp; bt_a%vbt_at_n = 0.3_wp
         bt_a%F_bt_u = 0.0_wp; bt_a%F_bt_v = 0.0_wp
         bt_a%bt_H_ref = 300.0_wp; bt_a%bt_eta_end = 0.0_wp

         bt_b%bt_ubt_end = 1.0_wp; bt_b%bt_vbt_end = 0.4_wp
         bt_b%ubt_at_n = 0.5_wp; bt_b%vbt_at_n = 0.3_wp
         bt_b%F_bt_u = 0.0_wp; bt_b%F_bt_v = 0.0_wp
         bt_b%bt_H_ref = 300.0_wp; bt_b%bt_eta_end = 0.0_wp

         ! Both runs: h-weighted ON.  bt_b also sets use_visc_rem ON
         ! but visc_rem stays at its init value 1.0 ⇒ should match.
         call apply_bt_correction(bt_a, ms_a, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true.)
         call apply_bt_correction(bt_b, ms_b, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true., use_visc_rem=.true.)

         max_diff = max(maxval(abs(ms_a%u_face_x_layer - ms_b%u_face_x_layer)), &
                        maxval(abs(ms_a%v_face_y_layer - ms_b%v_face_y_layer)))
         call check(error, max_diff < 1.0e-12_wp, &
                    "visc_rem ≡ 1: use_visc_rem=.true. must be bit-identical to .false.")
      end block checks
      call cleanup(ms_a, bt_a); call cleanup(ms_b, bt_b)
   end subroutine test_vr_unity_matches

   subroutine test_vr_bed_damped(error)
      !! Halving visc_rem at the bed (k=1) biases the per-layer
      !! Δu weight against the bed: bed gets LESS share, the other
      !! layers compensate.  Uniform h so the only asymmetry comes
      !! from visc_rem.  Depth-mean still preserved (norm absorbs
      !! the vr factor).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt_a, bt_b
      real(wp) :: ubt_mean_a, ubt_mean_b, du_bed_a, du_bed_b
      integer :: k, i_probe, j_probe
      checks: block
         call build_state(grid, ms_a, bt_a, 3)
         call build_state(grid, ms_b, bt_b, 3)
         i_probe = grid%nghost + 2; j_probe = grid%nghost + 2

         ! Uniform h, uniform initial u
         ms_a%h_layer = 100.0_wp; ms_b%h_layer = 100.0_wp
         ms_a%u_face_x_layer = 0.5_wp; ms_b%u_face_x_layer = 0.5_wp
         ms_a%v_face_y_layer = 0.0_wp; ms_b%v_face_y_layer = 0.0_wp

         bt_a%bt_ubt_end = 1.0_wp; bt_a%bt_vbt_end = 0.0_wp
         bt_a%ubt_at_n = 0.5_wp; bt_a%vbt_at_n = 0.0_wp
         bt_a%F_bt_u = 0.0_wp; bt_a%F_bt_v = 0.0_wp
         bt_a%bt_H_ref = 300.0_wp; bt_a%bt_eta_end = 0.0_wp

         bt_b%bt_ubt_end = 1.0_wp; bt_b%bt_vbt_end = 0.0_wp
         bt_b%ubt_at_n = 0.5_wp; bt_b%vbt_at_n = 0.0_wp
         bt_b%F_bt_u = 0.0_wp; bt_b%F_bt_v = 0.0_wp
         bt_b%bt_H_ref = 300.0_wp; bt_b%bt_eta_end = 0.0_wp

         ! bt_b: damp bed via visc_rem; bt_a stays at vr=1
         bt_b%visc_rem_u(:, :, 1) = 0.5_wp
         ! Layers 2 + 3 stay at the init 1.0

         call apply_bt_correction(bt_a, ms_a, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true.)
         call apply_bt_correction(bt_b, ms_b, 1.0_wp, skip_h_rescale=.true., &
                                  use_h_weighted=.true., use_visc_rem=.true.)

         ! Bed-layer Δu — bt_b should have SMALLER delta than bt_a
         ! since the visc_rem factor reduces the bed weight.
         du_bed_a = ms_a%u_face_x_layer(i_probe, j_probe, 1) - 0.5_wp
         du_bed_b = ms_b%u_face_x_layer(i_probe, j_probe, 1) - 0.5_wp
         call check(error, du_bed_b < du_bed_a - 1.0e-6_wp, &
                    "visc_rem(bed)=0.5: bed gets smaller Δu than h-only path")
         if (allocated(error)) exit checks

         ! Depth-mean preservation: u_avg should still equal bt_ubt_end = 1.0
         ! for both runs.  With uniform h, the depth-mean = (u1+u2+u3)/3.
         ubt_mean_a = (ms_a%u_face_x_layer(i_probe, j_probe, 1) &
                       + ms_a%u_face_x_layer(i_probe, j_probe, 2) &
                       + ms_a%u_face_x_layer(i_probe, j_probe, 3))/3.0_wp
         ubt_mean_b = (ms_b%u_face_x_layer(i_probe, j_probe, 1) &
                       + ms_b%u_face_x_layer(i_probe, j_probe, 2) &
                       + ms_b%u_face_x_layer(i_probe, j_probe, 3))/3.0_wp
         call check(error, abs(ubt_mean_a - 1.0_wp) < 1.0e-12_wp, &
                    "h-only run: depth-mean must equal bt_ubt_end")
         if (allocated(error)) exit checks
         call check(error, abs(ubt_mean_b - 1.0_wp) < 1.0e-12_wp, &
                    "visc_rem run: depth-mean must STILL equal bt_ubt_end")
      end block checks
      call cleanup(ms_a, bt_a); call cleanup(ms_b, bt_b)
   end subroutine test_vr_bed_damped

   subroutine test_skip_nonfinite(error)
      !! Hot-state armour: a NON-FINITE bt-correction input (the BT substep
      !! loop reaches Inf on at-floor columns in a supercritical state; the
      !! fold's `finite − Inf` would mint NaN into the layer velocity) at ONE
      !! face must (a) SKIP the fold write there (velocity left as-is for the
      !! truncation NaN-catch backstop — not silently zeroed), (b) leave the
      !! neighbour face's correction untouched, (c) be counted on `n_nonfin`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      real(wp) :: pinf, u_bad_before, u_good_before, u_bad_after, u_good_after
      integer :: nnf, ib, jb, ig, jg, k
      character(len=160) :: msg
      checks: block
         call build_state(grid, ms, bt, 3)
         ms%h_layer = 100.0_wp
         ms%u_face_x_layer = 0.5_wp
         ms%v_face_y_layer = 0.0_wp
         bt%bt_ubt_end = 1.0_wp; bt%bt_vbt_end = 0.0_wp
         bt%ubt_at_n = 0.5_wp; bt%vbt_at_n = 0.0_wp
         bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
         bt%bt_H_ref = 300.0_wp; bt%bt_eta_end = 0.0_wp

         pinf = ieee_value(0.0_wp, ieee_positive_inf)
         ib = NGHOST + 2; jb = NGHOST + 2   ! the blown-up face
         ig = NGHOST + 3; jg = NGHOST + 2   ! a finite neighbour
         bt%bt_ubt_end(ib, jb) = pinf       ! ⇒ Δu = Inf − 0.5 − 0 = Inf

         u_bad_before = ms%u_face_x_layer(ib, jb, 1)   ! 0.5
         u_good_before = ms%u_face_x_layer(ig, jg, 1)  ! 0.5

         call apply_bt_correction(bt, ms, 1.0_wp, skip_h_rescale=.true., n_nonfin=nnf)

         u_bad_after = ms%u_face_x_layer(ib, jb, 1)
         u_good_after = ms%u_face_x_layer(ig, jg, 1)

         ! (a) The blown-up face is SKIPPED — unchanged and still finite (NOT
         ! Inf/NaN written), at EVERY layer.
         write (msg, '("skip: bad face u=", es12.4, " (want unchanged ", es12.4, ")")') &
            u_bad_after, u_bad_before
         call check(error, u_bad_after == u_bad_before .and. ieee_is_finite(u_bad_after), trim(msg))
         if (allocated(error)) exit checks
         do k = 1, 3
            call check(error, ieee_is_finite(ms%u_face_x_layer(ib, jb, k)), &
                       "skip: blown-up face left non-finite at some layer")
            if (allocated(error)) exit checks
         end do
         ! (b) The neighbour is corrected normally: Δu = 1.0 − 0.5 − 0 = 0.5.
         call check(error, abs(u_good_after - (u_good_before + 0.5_wp)) < 1.0e-12_wp, &
                    "skip: finite neighbour face not corrected")
         if (allocated(error)) exit checks
         ! (c) Counted loudly (one non-finite u-face; all v faces finite).
         write (msg, '("skip: n_nonfin=", i0, " (want 1)")') nnf
         call check(error, nnf == 1, trim(msg))
      end block checks
      call cleanup(ms, bt)
   end subroutine test_skip_nonfinite

end module test_ocean_bt_corrector_hweight
