!! The barotropic-correction FOLD WEIGHT
!! (`rdb_barotropic_coupling::apply_bt_correction`).
!!
!! The fold adds `Δ·wt_k` to every layer velocity, `Δ` the barotropic
!! increment, with
!!
!!     wt_k = open_k·vr_k / ⟨vr⟩_h,     ⟨vr⟩_h = Σ h_o·vr / Σ h_o,
!!
!! `vr = visc_rem` when the kernel's `use_visc_rem` dummy is `.true.`, else 1,
!! and `open ≡ 1` unless `&vcoord_nml zfixed_closed_faces`.  Without visc_rem
!! that is the uniform fold, MOM6's `accel_layer_u(I,j,k) = u_accel_bt(I,j)`.
!! `&ocean_bt_nml correction_visc_rem` — the namelist path that used to drive
!! `use_visc_rem` — is RETIRED (D1 follow-up, 2026-10): MOM6 never weights
!! this fold, and roundabout's weighted version is what NaNs the 1-degree
!! Southern Ocean z* open-step case under `bbl_glue`.  `visc_rem_chain` uses
!! the UNIFORM branch below. The tests here call `apply_bt_correction`
!! directly with a raw `use_visc_rem` argument, so the weighted-fold KERNEL
!! path stays covered even though no live namelist reaches it any more.
!!
!! ### Why this file exists
!!
!! The retired `&ocean_bt_nml correction_h_weighted` used `wt = h/⟨h⟩_h`.
!! Both weights preserve the depth mean, so a depth-mean test cannot tell
!! them apart.  ENERGY can: for any depth-mean-preserving weight the
!! column KE change is
!!
!!     ΔKE = H·(ū·Δ + Δ²/2)  +  ½Δ²·Σ h(wt−1)²  +  Δ·Σ h·u′·(wt−1),
!!
!! the first term the barotropic KE change the fast loop accounts for,
!! the second a positive-definite source (`½Δ²H(κ−1)` for the h weight,
!! `κ = Σh³Σh/(Σh²)²`), the third a shear feedback.  On the stretched
!! 50-level tanh `z_fixed` stack the h weight has `κ − 1 ≈ 0.2`, and that
!! feedback grew the 1-degree Southern Ocean and the coastal-noise box to a
!! non-finite state.  The uniform fold has only the first term.
!!
!! Tests:
!!   * `uniform_fold_energy_stretched_column` — on a 500 m column cut from
!!     the shipped 50-level tanh stack (κ−1 asserted ≥ 0.15): depth mean
!!     moved by exactly Δ, every layer moved by exactly Δ, and the column
!!     KE change equals the barotropic `H(ūΔ + Δ²/2)` to round-off, while
!!     the retired h weight's extra energy on the SAME column is asserted
!!     to be decades above that round-off (the column discriminates).
!!   * `closed_faces_open_column` — the same with the bed layers CLOSED:
!!     closed layers untouched, OPEN-column mean moved by Δ, OPEN-column
!!     KE change barotropic; with and without visc_rem.
!!   * `visc_rem_unity_is_uniform` — `visc_rem ≡ 1` reproduces the uniform
!!     fold bit for bit.
!!   * `visc_rem_biases_against_bed` — damped bed gets less, depth mean
!!     preserved exactly, on the stretched column.
!!   * `skip_nonfinite_fold_input` — a non-finite Δ is skipped and counted,
!!     on the uniform and the visc_rem branch.
!!   * `h_weighted_refused` / `visc_rem_chain_standalone` — the retired
!!     knob is refused with its reason on the error ring; the visc_rem
!!     chain configures without it.
module test_ocean_bt_correction_weight
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: apply_bt_correction
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_vcoord, only: z_fixed_nominal_dz, ZFIXED_PROFILE_TANH, ZFIXED_DZ_OK
   use rdb_config, only: config_t, read_config_from_string, validate_config, &
                         ocean_bt_forcing_visc_rem_on, ocean_bt_renorm_visc_rem_on, &
                         ocean_bt_rem_from_visc_rem_on, ocean_bt_visc_rem_producer_on
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_error_ring, only: error_ring_clear, error_ring_get, error_ring_count, &
                             ERROR_RING_MSG_LEN
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
   implicit none
   private

   public :: collect_ocean_bt_correction_weight_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ_STACK = 50
      !! The shipped OM_1deg / southern_ocean_1deg / coastal_noise_box stack.
   real(wp), parameter :: H_STACK = 6500.0_wp, DZ_TOP = 2.0_wp
   real(wp), parameter :: COL_DEPTH = 500.0_wp
      !! A shelf-break column: 20 live layers, 2 m … ~90 m, κ − 1 ≈ 0.21.
   real(wp), parameter :: H_FILL = 1.0e-4_wp
      !! Inert filler thickness below the bed.
   real(wp), parameter :: DELTA = 0.05_wp
      !! Barotropic increment (m/s) the fast loop hands the fold.

contains

   subroutine collect_ocean_bt_correction_weight_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("uniform_fold_energy_stretched_column", test_uniform_energy), &
                  new_unittest("closed_faces_open_column", test_closed_faces), &
                  new_unittest("visc_rem_unity_is_uniform", test_vr_unity), &
                  new_unittest("visc_rem_biases_against_bed", test_vr_bed), &
                  new_unittest("skip_nonfinite_fold_input", test_skip_nonfinite), &
                  new_unittest("h_weighted_refused", test_h_weighted_refused), &
                  new_unittest("visc_rem_chain_standalone", test_visc_rem_chain), &
                  new_unittest("visc_rem_chain_on_helpers", test_chain_on_helpers), &
                  new_unittest("visc_rem_chain_switch_configures", test_chain_switch), &
                  new_unittest("correction_visc_rem_retired", test_correction_visc_rem_retired), &
                  new_unittest("visc_rem_chain_substep_drag_refused", test_chain_substep_drag), &
                  new_unittest("visc_rem_chain_strong_drag_accepted", test_chain_strong_drag), &
                  new_unittest("accel_visc_rem_retired", test_accel_visc_rem_retired) &
                  ]
   end subroutine collect_ocean_bt_correction_weight_tests

   ! ------------------------------------------------------------------
   ! Fixtures
   ! ------------------------------------------------------------------

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

   subroutine stretched_column(h, n_live)
      !! The 50-level tanh `z_fixed` stack cut at `COL_DEPTH`, BOTTOM-UP
      !! (`h(1)` the bed): full nominal layers down to the bed, a partial
      !! bed cell, `H_FILL` fillers below — `ocean_vcoord_z_fixed_target`'s
      !! shape at `η = 0`.  `n_live` is the number of live layers.
      real(wp), intent(out) :: h(NZ_STACK)
      integer, intent(out) :: n_live
      real(wp) :: dz(NZ_STACK), dz_list(1), z
      integer :: n, ierr
      dz_list = 0.0_wp
      call z_fixed_nominal_dz(ZFIXED_PROFILE_TANH, NZ_STACK, H_STACK, dz_list, DZ_TOP, &
                              0.5_wp, 0.25_wp, dz, ierr)
      if (ierr /= ZFIXED_DZ_OK) error stop "stretched_column: z_fixed_nominal_dz failed"
      h = H_FILL
      z = 0.0_wp
      n_live = 0
      do n = 1, NZ_STACK                     ! n = 1 is the SURFACE layer
         if (z >= COL_DEPTH) exit
         n_live = n_live + 1
         h(NZ_STACK + 1 - n) = min(dz(n), COL_DEPTH - z)
         z = z + dz(n)
      end do
   end subroutine stretched_column

   pure function kappa(h) result(kap)
      !! `κ = Σh³Σh/(Σh²)²` — 1 iff the layers are equally thick.
      real(wp), intent(in) :: h(:)
      real(wp) :: kap
      kap = sum(h**3)*sum(h)/sum(h**2)**2
   end function kappa

   pure function sheared_profile(nz) result(u)
      !! A non-trivial baroclinic background (surface-intensified, sign
      !! change at depth) so the shear-feedback term is NOT zero.
      integer, intent(in) :: nz
      real(wp) :: u(nz)
      integer :: k
      do k = 1, nz
         u(k) = 0.02_wp + 0.25_wp*(real(k, wp)/real(nz, wp))**3 - 0.05_wp*sin(0.7_wp*real(k, wp))
      end do
   end function sheared_profile

   subroutine seed(ms, bt, h, u, ubt_n, dlt)
      !! Every column/face the same: `h`, layer velocity `u` (u- and
      !! v-faces), barotropic state with `ubt_end − ubt_at_n = dlt`,
      !! `F_bt = 0`.
      type(multilayer_state_t), intent(inout) :: ms
      type(barotropic_workstate_t), intent(inout) :: bt
      real(wp), intent(in) :: h(:), u(:), ubt_n, dlt
      integer :: k
      do k = 1, size(h)
         ms%h_layer(:, :, k) = h(k)
         ms%u_face_x_layer(:, :, k) = u(k)
         ms%v_face_y_layer(:, :, k) = -u(k)
      end do
      bt%ubt_at_n = ubt_n; bt%bt_ubt_end = ubt_n + dlt
      bt%vbt_at_n = -ubt_n; bt%bt_vbt_end = -ubt_n - dlt
      bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
      bt%bt_H_ref = sum(h); bt%bt_eta_end = 0.0_wp
   end subroutine seed

   ! ------------------------------------------------------------------
   ! Tests
   ! ------------------------------------------------------------------

   subroutine test_uniform_energy(error)
      !! The ENERGY statement on a stretched column.  All columns carry the
      !! same `h`, so `h_face = h` exactly and the face column is the cell
      !! column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      real(wp) :: h(NZ_STACK), u0(NZ_STACK), u1(NZ_STACK), w(NZ_STACK)
      real(wp) :: hh, ubar, ke0, ke1, dke_bt, tol, kap, e_hw, max_layer_err
      integer :: n_live, i, j
      character(len=200) :: msg
      call stretched_column(h, n_live)
      kap = kappa(h)
      u0 = sheared_profile(NZ_STACK)
      hh = sum(h)
      ubar = sum(h*u0)/hh
      checks: block
         write (msg, '("fixture: kappa-1 = ", es10.3, " on ", i0, " live layers (want >= 0.15)")') &
            kap - 1.0_wp, n_live
         call check(error, kap - 1.0_wp >= 0.15_wp .and. n_live == 20, trim(msg))
         if (allocated(error)) exit checks

         call build_state(grid, ms, bt, NZ_STACK)
         call make_cartesian_metrics(metrics, grid)
         call seed(ms, bt, h, u0, ubar, DELTA)
         call apply_bt_correction(bt, ms, 1.0_wp, metrics, skip_h_rescale=.true.)

         i = NGHOST + 2; j = NGHOST + 2
         u1 = ms%u_face_x_layer(i, j, :)
         ! Every layer, the 2 m surface layer and the fillers included, moves
         ! by exactly Δ (the uniform fold, MOM6's u_accel_bt).
         max_layer_err = maxval(abs((u1 - u0) - DELTA))
         write (msg, '("uniform fold: max |du_k - Delta| = ", es10.3)') max_layer_err
         call check(error, max_layer_err < 1.0e-15_wp, trim(msg))
         if (allocated(error)) exit checks
         ! Depth mean moved by exactly Δ.
         write (msg, '("depth mean: ", es22.15, " want ", es22.15)') sum(h*u1)/hh, ubar + DELTA
         call check(error, abs(sum(h*u1)/hh - (ubar + DELTA)) < 1.0e-14_wp, trim(msg))
         if (allocated(error)) exit checks

         ! Column KE change == the barotropic KE change.
         ke0 = 0.5_wp*sum(h*u0**2)
         ke1 = 0.5_wp*sum(h*u1**2)
         dke_bt = hh*(ubar*DELTA + 0.5_wp*DELTA**2)
         tol = 1.0e-12_wp*ke1
         write (msg, '("uniform fold: dKE - dKE_bt = ", es10.3, " (tol ", es10.3, ")")') &
            (ke1 - ke0) - dke_bt, tol
         call check(error, abs((ke1 - ke0) - dke_bt) < tol, trim(msg))
         if (allocated(error)) exit checks

         ! The column DISCRIMINATES: the retired h weight's positive-definite
         ! source on it, `½Δ²Σh(w−1)² = ½Δ²H(κ−1)`, is decades above `tol`.
         w = h/(sum(h*h)/hh)
         e_hw = 0.5_wp*DELTA**2*sum(h*(w - 1.0_wp)**2)
         write (msg, '("retired h weight source ", es10.3, " (want > 1e4 * tol = ", es10.3, ")")') &
            e_hw, 1.0e4_wp*tol
         call check(error, e_hw > 1.0e4_wp*tol, trim(msg))
      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt)
   end subroutine test_uniform_energy

   subroutine test_closed_faces(error)
      !! `zfixed_closed_faces`: the filler layers below the bed are CLOSED
      !! at every face.  The closed layers must not move; the OPEN column
      !! mean moves by exactly Δ and its KE change is barotropic (uniform
      !! fold); with visc_rem the open mean is still exact.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      real(wp) :: h(NZ_STACK), u0(NZ_STACK), u1(NZ_STACK), op(NZ_STACK)
      real(wp) :: ho, ubar_o, dke, dke_bt
      integer :: n_live, n_closed, i, j, pass
      character(len=200) :: msg
      call stretched_column(h, n_live)
      n_closed = NZ_STACK - n_live
      op = 1.0_wp
      op(1:n_closed) = 0.0_wp
      u0 = sheared_profile(NZ_STACK)
      u0(1:n_closed) = 0.0_wp                ! masked closed layers
      ho = sum(h*op)
      ubar_o = sum(h*op*u0)/ho
      do pass = 1, 2                         ! 1: uniform, 2: visc_rem
         checks: block
            call build_state(grid, ms, bt, NZ_STACK)
            call make_cartesian_metrics(metrics, grid, nz_closed=NZ_STACK)
            do j = 1, NZ_STACK
               metrics%open_u(:, :, j) = op(j)
               metrics%open_v(:, :, j) = op(j)
            end do
            !$acc update device(metrics%open_u, metrics%open_v)
            metrics%use_closed_faces = .true.
            call seed(ms, bt, h, u0, ubar_o, DELTA)
            if (pass == 2) then
               do j = 1, NZ_STACK                ! damped near the bed
                  bt%visc_rem_u(:, :, j) = max(0.05_wp, min(1.0_wp, 0.3_wp + 0.05_wp*real(j - n_closed, wp)))
                  bt%visc_rem_v(:, :, j) = bt%visc_rem_u(:, :, j)
               end do
            end if
            call apply_bt_correction(bt, ms, 1.0_wp, metrics, skip_h_rescale=.true., &
                                     use_visc_rem=(pass == 2))
            i = NGHOST + 2; j = NGHOST + 2
            u1 = ms%u_face_x_layer(i, j, :)
            write (msg, '("pass ", i0, ": closed layers moved by ", es10.3)') pass, &
               maxval(abs(u1(1:n_closed)))
            call check(error, all(u1(1:n_closed) == 0.0_wp), trim(msg))
            if (allocated(error)) exit checks
            write (msg, '("pass ", i0, ": open mean err ", es10.3)') pass, &
               sum(h*op*u1)/ho - (ubar_o + DELTA)
            call check(error, abs(sum(h*op*u1)/ho - (ubar_o + DELTA)) < 1.0e-14_wp, trim(msg))
            if (allocated(error)) exit checks
            if (pass == 1) then
               dke = 0.5_wp*sum(h*op*(u1**2 - u0**2))
               dke_bt = ho*(ubar_o*DELTA + 0.5_wp*DELTA**2)
               write (msg, '("open-column dKE - dKE_bt = ", es10.3)') dke - dke_bt
               call check(error, abs(dke - dke_bt) < 1.0e-12_wp*abs(dke_bt), trim(msg))
               if (allocated(error)) exit checks
            end if
         end block checks
         call destroy_cartesian_metrics(metrics)
         call cleanup(ms, bt)
         if (allocated(error)) return
      end do
   end subroutine test_closed_faces

   subroutine test_vr_unity(error)
      !! `visc_rem ≡ 1` (the producer's no-drag answer, exact by the
      !! operator's row sums — `test_ocean_visc_rem`) must reproduce the
      !! uniform fold BIT FOR BIT: `⟨vr⟩_h = Σh/Σh = 1` exactly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(barotropic_workstate_t) :: bt_a, bt_b
      type(ocean_metrics_t) :: metrics
      real(wp) :: h(NZ_STACK), u0(NZ_STACK)
      integer :: n_live
      call stretched_column(h, n_live)
      u0 = sheared_profile(NZ_STACK)
      call build_state(grid, ms_a, bt_a, NZ_STACK)
      call build_state(grid, ms_b, bt_b, NZ_STACK)
      call make_cartesian_metrics(metrics, grid)
      call seed(ms_a, bt_a, h, u0, 0.1_wp, DELTA)
      call seed(ms_b, bt_b, h, u0, 0.1_wp, DELTA)
      bt_b%visc_rem_u = 1.0_wp; bt_b%visc_rem_v = 1.0_wp
      call apply_bt_correction(bt_a, ms_a, 1.0_wp, metrics, skip_h_rescale=.true.)
      call apply_bt_correction(bt_b, ms_b, 1.0_wp, metrics, skip_h_rescale=.true., &
                               use_visc_rem=.true.)
      call check(error, all(ms_a%u_face_x_layer == ms_b%u_face_x_layer) .and. &
                 all(ms_a%v_face_y_layer == ms_b%v_face_y_layer), &
                 "visc_rem = 1: the weighted fold must equal the uniform fold bit for bit")
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms_a, bt_a); call cleanup(ms_b, bt_b)
   end subroutine test_vr_unity

   subroutine test_vr_bed(error)
      !! A damped bed (`visc_rem(bed) = 0.4`) gets LESS than Δ, the
      !! undamped layers more, and the depth mean still moves by exactly Δ
      !! — on the stretched column, where the old `h·vr` weight would also
      !! have starved the 2 m surface layer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      real(wp) :: h(NZ_STACK), u0(NZ_STACK), du(NZ_STACK), hh, ubar, vr_bar
      integer :: n_live, kbed, i, j
      character(len=200) :: msg
      call stretched_column(h, n_live)
      kbed = NZ_STACK - n_live + 1             ! the partial bed cell
      u0 = sheared_profile(NZ_STACK)
      hh = sum(h)
      ubar = sum(h*u0)/hh
      checks: block
         call build_state(grid, ms, bt, NZ_STACK)
         call make_cartesian_metrics(metrics, grid)
         call seed(ms, bt, h, u0, ubar, DELTA)
         bt%visc_rem_u = 1.0_wp; bt%visc_rem_v = 1.0_wp
         bt%visc_rem_u(:, :, 1:kbed) = 0.4_wp
         bt%visc_rem_v(:, :, 1:kbed) = 0.4_wp
         call apply_bt_correction(bt, ms, 1.0_wp, metrics, skip_h_rescale=.true., &
                                  use_visc_rem=.true.)
         i = NGHOST + 2; j = NGHOST + 2
         du = ms%u_face_x_layer(i, j, :) - u0
         vr_bar = (0.4_wp*sum(h(1:kbed)) + sum(h(kbed + 1:)))/hh
         write (msg, '("bed du = ", es12.5, " want ", es12.5)') du(kbed), DELTA*0.4_wp/vr_bar
         call check(error, abs(du(kbed) - DELTA*0.4_wp/vr_bar) < 1.0e-14_wp, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("surface du = ", es12.5, " want ", es12.5)') du(NZ_STACK), DELTA/vr_bar
         call check(error, abs(du(NZ_STACK) - DELTA/vr_bar) < 1.0e-14_wp .and. &
                    du(NZ_STACK) > DELTA, trim(msg))
         if (allocated(error)) exit checks
         write (msg, '("visc_rem fold: depth-mean err ", es10.3)') sum(h*du)/hh - DELTA
         call check(error, abs(sum(h*du)/hh - DELTA) < 1.0e-14_wp, trim(msg))
      end block checks
      call destroy_cartesian_metrics(metrics)
      call cleanup(ms, bt)
   end subroutine test_vr_bed

   subroutine test_skip_nonfinite(error)
      !! Hot-state armour, on BOTH fold branches: a non-finite Δ at one
      !! face (a blown-up BT loop) is SKIPPED there — velocity left as-is
      !! for the truncation NaN-catch, never NaN-ed — the neighbour is
      !! corrected normally, and the face is counted on `n_nonfin`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt
      type(ocean_metrics_t) :: metrics
      integer :: nnf, ib, jb, ig, jg, pass
      character(len=160) :: msg
      do pass = 1, 2
         checks: block
            call build_state(grid, ms, bt, 3)
            call make_cartesian_metrics(metrics, grid)
            ms%h_layer(:, :, 1) = 10.0_wp
            ms%h_layer(:, :, 2) = 100.0_wp
            ms%h_layer(:, :, 3) = 200.0_wp
            ms%u_face_x_layer = 0.5_wp
            ms%v_face_y_layer = 0.0_wp
            bt%bt_ubt_end = 1.0_wp; bt%bt_vbt_end = 0.0_wp
            bt%ubt_at_n = 0.5_wp; bt%vbt_at_n = 0.0_wp
            bt%F_bt_u = 0.0_wp; bt%F_bt_v = 0.0_wp
            bt%bt_H_ref = 310.0_wp; bt%bt_eta_end = 0.0_wp
            bt%visc_rem_u = 1.0_wp; bt%visc_rem_v = 1.0_wp
            ib = NGHOST + 2; jb = NGHOST + 2
            ig = NGHOST + 3; jg = NGHOST + 2
            bt%bt_ubt_end(ib, jb) = ieee_value(0.0_wp, ieee_positive_inf)
            call apply_bt_correction(bt, ms, 1.0_wp, metrics, skip_h_rescale=.true., &
                                     use_visc_rem=(pass == 2), n_nonfin=nnf)
            write (msg, '("pass ", i0, ": blown-up face not skipped")') pass
            call check(error, all(ms%u_face_x_layer(ib, jb, :) == 0.5_wp), trim(msg))
            if (allocated(error)) exit checks
            write (msg, '("pass ", i0, ": finite neighbour not corrected")') pass
            call check(error, all(abs(ms%u_face_x_layer(ig, jg, :) - 1.0_wp) < 1.0e-14_wp), trim(msg))
            if (allocated(error)) exit checks
            write (msg, '("pass ", i0, ": n_nonfin=", i0, " (want 1)")') pass, nnf
            call check(error, nnf == 1, trim(msg))
         end block checks
         call destroy_cartesian_metrics(metrics)
         call cleanup(ms, bt)
         if (allocated(error)) return
      end do
   end subroutine test_skip_nonfinite

   subroutine test_h_weighted_refused(error)
      !! `correction_h_weighted = .true.` is REFUSED at configure, and the
      !! specific reason — energy non-conservation, no MOM6 counterpart,
      !! `visc_rem_chain` as the replacement — is on the error ring, not
      !! only the generic rollup.  The same namelist without the key
      !! configures (so the refusal is the knob, not the fixture).
      type(error_type), allocatable, intent(out) :: error
      character(len=ERROR_RING_MSG_LEN) :: msg
      logical :: found
      integer :: n
      call expect_status(error, base_nml(""), .true., "the fixture namelist")
      if (allocated(error)) return
      call error_ring_clear()
      call expect_status(error, base_nml("correction_h_weighted = .true."), .false., &
                         "correction_h_weighted = .true.")
      if (allocated(error)) return
      found = .false.
      do n = 0, error_ring_count() - 1
         msg = error_ring_get(n)
         if (index(msg, "correction_h_weighted is RETIRED") > 0) then
            found = .true.
            exit
         end if
      end do
      call check(error, found, "the retirement reason must be on the error ring")
      if (allocated(error)) return
      call check(error, index(msg, "energy-non-conserving") > 0 .and. &
                 index(msg, "MOM6 has no h-weighted fold") > 0 .and. &
                 index(msg, "visc_rem_chain") > 0, &
                 "the refusal must name the defect, MOM6 and the replacement: "//trim(msg))
   end subroutine test_h_weighted_refused

   subroutine test_visc_rem_chain(error)
      !! D1 (revised once MOM6 settled the BT-correction fold question):
      !! `correction_visc_rem` (the weighted fold) alone is now REFUSED.
      !! The three real consumers are each SELF-SUFFICIENT — no longer
      !! "require" a separate producer knob — and compose freely.
      type(error_type), allocatable, intent(out) :: error
      call expect_status(error, base_nml("correction_visc_rem = .true."), .false., &
                         "correction_visc_rem alone (retired)")
      if (allocated(error)) return
      call expect_status(error, base_nml("forcing_visc_rem = .true."), .true., &
                         "forcing_visc_rem alone (self-sufficient)")
      if (allocated(error)) return
      call expect_status(error, base_nml("renorm_visc_rem = .true."), .true., &
                         "renorm_visc_rem alone (self-sufficient)")
      if (allocated(error)) return
      call expect_status(error, base_nml("bt_rem_from_visc_rem = .true."), .true., &
                         "bt_rem_from_visc_rem alone (self-sufficient)")
      if (allocated(error)) return
      call expect_status(error, base_nml("forcing_visc_rem = .true., renorm_visc_rem = .true., "// &
                                         "bt_rem_from_visc_rem = .true."), .true., &
                         "forcing_visc_rem + renorm_visc_rem + bt_rem_from_visc_rem")
   end subroutine test_visc_rem_chain

   subroutine test_chain_on_helpers(error)
      !! PR-3 (D1): the `ocean_bt_*_on` helpers in `rdb_config` are
      !! `visc_rem_chain .OR. <the direct knob>` — never a superset, never
      !! a subset — and `ocean_bt_visc_rem_producer_on` is the OR of the
      !! three real consumers plus the (retired, direct-cfg-only) raw
      !! `correction_visc_rem` field.  Truth-table, no namelist I/O.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      ! PR-4 (the flip): `visc_rem_chain` now defaults ON at the type level,
      ! so this truth-table test (exercising every combination by hand,
      ! starting from "neither set") pins it OFF explicitly first.
      cfg%ocean%bt%visc_rem_chain = .false.

      ! Neither set: every helper false.
      call check(error,.not. ocean_bt_visc_rem_producer_on(cfg), "neither: producer off")
      if (allocated(error)) return
      call check(error,.not. ocean_bt_forcing_visc_rem_on(cfg), "neither: forcing off")
      if (allocated(error)) return
      call check(error,.not. ocean_bt_renorm_visc_rem_on(cfg), "neither: renorm off")
      if (allocated(error)) return
      call check(error,.not. ocean_bt_rem_from_visc_rem_on(cfg), "neither: bt_rem off")
      if (allocated(error)) return

      ! The retired raw knob alone (direct cfg construction, bypassing the
      ! nml refusal): the PRODUCER still runs (so a test exercising it
      ! directly gets live visc_rem), but none of the three real
      ! consumers turn on.
      cfg%ocean%bt%correction_visc_rem = .true.
      call check(error, ocean_bt_visc_rem_producer_on(cfg), "direct correction_visc_rem: producer on")
      if (allocated(error)) return
      call check(error,.not. ocean_bt_forcing_visc_rem_on(cfg), "direct correction_visc_rem: forcing stays off")
      if (allocated(error)) return
      cfg%ocean%bt%correction_visc_rem = .false.

      ! Each real consumer alone also turns the producer on.
      cfg%ocean%bt%forcing_visc_rem = .true.
      call check(error, ocean_bt_visc_rem_producer_on(cfg), "forcing_visc_rem alone: producer on")
      if (allocated(error)) return
      cfg%ocean%bt%forcing_visc_rem = .false.

      ! Chain only: producer + the three real consumers, NEVER the
      ! retired weighted fold.
      cfg%ocean%bt%visc_rem_chain = .true.
      call check(error, ocean_bt_visc_rem_producer_on(cfg), "chain: producer on")
      if (allocated(error)) return
      call check(error, ocean_bt_forcing_visc_rem_on(cfg), "chain: forcing on")
      if (allocated(error)) return
      call check(error, ocean_bt_renorm_visc_rem_on(cfg), "chain: renorm on")
      if (allocated(error)) return
      call check(error, ocean_bt_rem_from_visc_rem_on(cfg), "chain: bt_rem on")
      if (allocated(error)) return
      call check(error,.not. cfg%ocean%bt%correction_visc_rem, &
                 "chain: the retired raw weighted-fold field itself stays untouched")
   end subroutine test_chain_on_helpers

   subroutine test_chain_switch(error)
      !! `visc_rem_chain = .true.` alone configures exactly like setting
      !! the three real consumer knobs by hand (the `test_visc_rem_chain`
      !! case two lines up), with no other knob touched — in particular
      !! it does NOT turn on the retired weighted BT-correction fold.
      type(error_type), allocatable, intent(out) :: error
      call expect_status(error, base_nml("visc_rem_chain = .true."), .true., &
                         "visc_rem_chain alone")
   end subroutine test_chain_switch

   subroutine test_correction_visc_rem_retired(error)
      !! `&ocean_bt_nml correction_visc_rem = .true.` is REFUSED at
      !! configure, naming MOM6's uniform `accel_layer_u`, the NaN
      !! isolation and `visc_rem_chain` as the replacement on the error
      !! ring, not only the generic rollup.
      type(error_type), allocatable, intent(out) :: error
      character(len=ERROR_RING_MSG_LEN) :: msg
      logical :: found
      integer :: n
      call error_ring_clear()
      call expect_status(error, base_nml("correction_visc_rem = .true."), .false., &
                         "correction_visc_rem = .true.")
      if (allocated(error)) return
      found = .false.
      do n = 0, error_ring_count() - 1
         msg = error_ring_get(n)
         if (index(msg, "correction_visc_rem is RETIRED") > 0) then
            found = .true.
            exit
         end if
      end do
      call check(error, found, "the retirement reason must be on the error ring")
      if (allocated(error)) return
      call check(error, index(msg, "accel_layer_u") > 0 .and. &
                 index(msg, "UNIFORMLY") > 0 .and. &
                 index(msg, "visc_rem_chain") > 0, &
                 "the refusal must name MOM6's uniform fold and the replacement: "//trim(msg))
   end subroutine test_correction_visc_rem_retired

   subroutine test_chain_substep_drag(error)
      !! D2: `visc_rem_chain` composes with `substep_drag` exactly like
      !! `bt_rem_from_visc_rem` already does — mutually exclusive,
      !! fail-loud (double-counted bed drag).
      type(error_type), allocatable, intent(out) :: error
      call expect_status(error, base_nml("visc_rem_chain = .true., substep_drag = .true."), &
                         .false., "visc_rem_chain + substep_drag")
   end subroutine test_chain_substep_drag

   subroutine test_chain_strong_drag(error)
      !! D3: `strong_drag` (MOM6 BT_STRONG_DRAG, opt-in, default off) is
      !! reachable through `visc_rem_chain` exactly as it already is
      !! through `bt_rem_from_visc_rem`.
      type(error_type), allocatable, intent(out) :: error
      call expect_status(error, base_nml("visc_rem_chain = .true., strong_drag = .true."), &
                         .true., "visc_rem_chain + strong_drag")
   end subroutine test_chain_strong_drag

   subroutine test_accel_visc_rem_retired(error)
      !! PR-3: `&ocean_vdiff_nml accel_visc_rem = .true.` is REFUSED at
      !! configure (no MOM6 state-update equivalent — see its docstring
      !! in `rdb_config.F90`), naming the real MOM6 mechanisms on the
      !! error ring.  The underlying kernels + their own direct unit
      !! tests (`tests/test_ocean_accel_visc_rem.F90`) are untouched —
      !! only the configure-time path from a namelist is refused here.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml
      character(len=ERROR_RING_MSG_LEN) :: msg
      logical :: found
      integer :: n

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 400.0 /"//new_line("a")// &
            "&ocean_vdiff_nml accel_visc_rem = .true. /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")

      call error_ring_clear()
      call expect_status(error, nml, .false., "accel_visc_rem = .true.")
      if (allocated(error)) return
      found = .false.
      do n = 0, error_ring_count() - 1
         msg = error_ring_get(n)
         if (index(msg, "accel_visc_rem is RETIRED") > 0) then
            found = .true.
            exit
         end if
      end do
      call check(error, found, "the retirement reason must be on the error ring")
      if (allocated(error)) return
      call check(error, index(msg, "renorm_visc_rem") > 0 .and. &
                 index(msg, "rescale_strong_drag") > 0, &
                 "the refusal must name the real MOM6 mechanisms: "//trim(msg))
   end subroutine test_accel_visc_rem_retired

   pure function base_nml(bt_body) result(nml)
      !! A minimal in-envelope `z_fixed` namelist; `bt_body` is the extra
      !! `&ocean_bt_nml` entries under test.
      character(len=*), intent(in) :: bt_body
      character(len=:), allocatable :: nml, bt_l
      bt_l = "auto_n_inner = .false., n_inner = 8"
      if (len_trim(bt_body) > 0) bt_l = bt_l//", "//bt_body
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 400.0 /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true. /"//new_line("a")// &
            "&ocean_bt_nml "//bt_l//" /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function base_nml

   subroutine expect_status(error, nml, valid, what)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, what
      logical, intent(in) :: valid
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "namelist must PARSE: "//what)
      if (allocated(error)) return
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      if (valid) then
         call check(error, ierr == OCEAN_STATUS_OK, "must CONFIGURE: "//what)
      else
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, "must be REFUSED: "//what)
      end if
   end subroutine expect_status

end module test_ocean_bt_correction_weight
