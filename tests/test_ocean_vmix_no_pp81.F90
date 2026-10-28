module test_ocean_vmix_no_pp81
   !! `&ocean_vmix_nml interior_closure = "none"` (no PP81 interior, the
   !! ACCESS-OM3 stance) and `bkgnd_prandtl` (momentum background tied to
   !! the tracer background in use).
   !!
   !! Under "none" the interior `kv`/`kt` are reset to zero each stage and the
   !! `vmix_assemble` floor supplies the base, so on a column with no other
   !! contributor the assembled interior field IS the background: the scalar
   !! `kv_bg`/`kt_bg`/`ks_bg`, or the Henyey, Wright & Flatte (1986) /
   !! Harrison & Hallberg (2008) `max(bkgnd_kd_min, kt_bg*L(phi))`.  Every
   !! full-dispatch case drives `ocean_dyn_step` (so the real
   !! `vmix_apply_in_stage` branch is exercised) on a strongly SHEARED, weakly
   !! stratified column where Pacanowski & Philander (1981) gives a kv far above
   !! the background -- the PP81 run on the same column is the non-vacuity
   !! witness for every "== background" assertion.
   !!
   !! Cases:
   !!   * none_scalar_background: kt/kv/ks == kt_bg/kv_bg/ks_bg bit-exactly;
   !!     PP81 on the same column gives max kv > 20 kv_bg.
   !!   * none_henyey_three_regimes: kt follows max(kd_min, kt_bg L(phi)) on
   !!     columns in the floor-bound, factor-bound and poleward-clamped
   !!     regimes (all three must occur); kv stays the unscaled kv_bg.
   !!   * prandtl_ties_kv_to_kt: bkgnd_prandtl=1, kv_bg=0 gives kv == kt (scalar and
   !!     Henyey, bit-exactly), and prandtl=2, kv_bg=1e-5 gives 2 kt + 1e-5.
   !!   * none_kappa_shear: with kappa-shear on and ONE unstable interface
   !!     (Ri << 0.25), kv/kt exceed the background there and only there.
   !!   * parse / defaults: the selector parser, the reserved spellings, the
   !!     defaults (pp81, prandtl off).  The configure refusals terminate
   !!     the process, so the predicate they share is what is asserted.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_assemble, vmix_split_kd_heat_salt, &
                             vmix_seed_interior_none, henyey_lat_factor_impl, &
                             parse_vmix_interior_closure, &
                             vmix_interior_closure_is_implemented, &
                             VMIX_INTERIOR_PP81, VMIX_INTERIOR_NONE, &
                             VMIX_INTERIOR_LARGE94, VMIX_INTERIOR_CVMIX
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t, kappa_shear_compute, &
                                    kappa_shear_merge_into_kv_kt
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step
   implicit none
   private

   public :: collect_ocean_vmix_no_pp81_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 100.0_wp
   real(wp), parameter :: KT_BG = 3.0e-5_wp
   real(wp), parameter :: KS_BG = 5.0e-5_wp
   real(wp), parameter :: KV_BG = 2.0e-4_wp
   integer, parameter :: NZK = 12, KJ = 6
   real(wp), parameter :: DZ = 5.0_wp, N2 = 3.0e-5_wp

contains

   subroutine collect_ocean_vmix_no_pp81_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vmix_none_scalar_background", test_none_scalar), &
                  new_unittest("vmix_none_henyey_three_regimes", test_none_henyey), &
                  new_unittest("vmix_prandtl_ties_kv_to_kt", test_prandtl), &
                  new_unittest("vmix_none_kappa_shear", test_none_kappa_shear), &
                  new_unittest("vmix_none_seed_henyey", test_seed_henyey), &
                  new_unittest("vmix_interior_closure_parse", test_parse_and_defaults) &
                  ]
   end subroutine collect_ocean_vmix_no_pp81_tests

   subroutine run_step(grid, metrics, ms, vmix)
      !! One `ocean_dyn_step` with every slot device-mapped; pulls kv/kt/ks
      !! back to the host.  The caller owns `ms`, `vmix` and `metrics`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn

      call ct%init(grid, nz_ml=NZ)
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
      !$acc enter data copyin(cor)
      call cor%enter_data()
      !$acc enter data copyin(pgf)
      call pgf%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
      !$acc enter data copyin(bd)
      call bd%enter_data()
      !$acc enter data copyin(ss)
      call ss%enter_data()
      !$acc enter data copyin(va)
      call va%enter_data()
      !$acc enter data copyin(hd)
      call hd%enter_data()
      !$acc enter data copyin(vd)
      call vd%enter_data()
      !$acc enter data copyin(vmix)
      call vmix%enter_data()

      call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, 1.0_wp)

      !$acc update self(vmix%kv, vmix%kt, vmix%ks)
      call vmix%exit_data()
      !$acc exit data delete(vmix)
      call vd%exit_data()
      !$acc exit data delete(vd)
      call hd%exit_data()
      !$acc exit data delete(hd)
      call va%exit_data()
      !$acc exit data delete(va)
      call ss%exit_data()
      !$acc exit data delete(ss)
      call bd%exit_data()
      !$acc exit data delete(bd)
      call hv%exit_data()
      !$acc exit data delete(hv)
      call pgf%exit_data()
      !$acc exit data delete(pgf)
      call cor%exit_data()
      !$acc exit data delete(cor)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)

      call dyn%destroy(); call eos%destroy(); call vd%destroy(); call hd%destroy()
      call va%destroy(); call ss%destroy(); call bd%destroy(); call hv%destroy()
      call pgf%destroy(); call cor%destroy(); call ct%destroy()
   end subroutine run_step

   subroutine make_sheared_column(grid, ms, vmix, closure, lat_south, dlat)
      !! Fresh grid/state/vmix/metrics-ready column: u increases 0.3 m/s per
      !! layer (S^2 ~ 9e-6 1/s^2), T warms 0.05 K per layer toward the surface
      !! (N^2 ~ 1e-6 1/s^2), so Ri ~ 0.1 and PP81 mixes hard.  `vmix` carries the
      !! scalar test backgrounds and the chosen interior closure.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: closure
      real(wp), intent(in), optional :: lat_south, dlat
      type(eos_t) :: eos
      integer :: k

      call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      ms%h_layer = H0
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%u_face_x_layer(:, :, k) = 0.3_wp*real(k, wp)
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (eos%T_ref + 0.05_wp*real(k - 2, wp))*H0
      end do
      call eos%destroy()
      vmix%use_closure = .true.
      vmix%use_kpp = .false.
      vmix%interior_closure = closure
      vmix%pp81_nu_bg = KV_BG
      vmix%pp81_kappa_bg = KT_BG
      call vmix%seed_backgrounds()
      vmix%ks_bg = KS_BG
   end subroutine make_sheared_column

   subroutine test_none_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(multilayer_state_t) :: ms_p
      type(ocean_vmix_t) :: vmix_p
      integer :: k
      real(wp) :: kv_pp81_max
      checks: block
         ! ---- PP81 witness on the identical column ----
         call make_sheared_column(grid, ms_p, vmix_p, VMIX_INTERIOR_PP81)
         call make_spherical_metrics(metrics, grid, 0.0_wp, 10.0_wp, 2.0_wp, 2.0_wp, 6.371e6_wp)
         call run_step(grid, metrics, ms_p, vmix_p)
         kv_pp81_max = maxval(vmix_p%kv(NGHOST + 1:grid%nx_total - NGHOST, &
                                        NGHOST + 1:grid%ny_total - NGHOST, 2:NZ))
         call check(error, kv_pp81_max > 20.0_wp*KV_BG, &
                    "test setup vacuous: PP81 kv on the sheared column is not far above the background")
         if (allocated(error)) exit checks

         ! ---- interior_closure = none on the same column ----
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         call run_step(grid, metrics, ms, vmix)
         do k = 2, NZ
            call check(error, all(vmix%kt(:, :, k) == KT_BG), &
                       "none: interior kt is not exactly kt_bg")
            if (allocated(error)) exit checks
            call check(error, all(vmix%kv(:, :, k) == KV_BG), &
                       "none: interior kv is not exactly kv_bg (PP81 ran or state accumulated)")
            if (allocated(error)) exit checks
            ! ks is derived from kt, then floored at ks_bg (> kt_bg here).
            call check(error, all(vmix%ks(:, :, k) == KS_BG), &
                       "none: interior ks is not exactly ks_bg")
            if (allocated(error)) exit checks
         end do
         call check(error, all(vmix%kv(:, :, 1) == 0.0_wp) .and. all(vmix%kv(:, :, NZ + 1) == 0.0_wp), &
                    "none: closed-BC boundary interfaces must stay zero")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call vmix%destroy(); call ms%destroy(); call vmix_p%destroy(); call ms_p%destroy()
   end subroutine test_none_scalar

   subroutine test_none_henyey(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: KD_MIN = 5.0e-6_wp, N0_2OM = 20.0_wp, MAX_LAT = 40.0_wp
      integer :: i, j, k, n_floor, n_factor, n_clamp
      real(wp) :: lat, lfac, want
      checks: block
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         ! Rows near -2, 10, 22, 34, 46, 58 degN: floor-bound, factor-bound,
         ! and (poleward of MAX_LAT) clamped.
         call make_spherical_metrics(metrics, grid, 0.0_wp, -8.0_wp, 2.0_wp, 12.0_wp, 6.371e6_wp)
         vmix%bkgnd_henyey = .true.
         vmix%bkgnd_kd_min = KD_MIN
         vmix%bkgnd_henyey_n0_2omega = N0_2OM
         vmix%bkgnd_henyey_max_lat = MAX_LAT
         call run_step(grid, metrics, ms, vmix)

         n_floor = 0; n_factor = 0; n_clamp = 0
         do j = NGHOST + 1, grid%ny_total - NGHOST
            i = NGHOST + 2
            lat = metrics%geolatT(i, j)
            lfac = henyey_lat_factor_impl(lat, N0_2OM, MAX_LAT)
            want = max(KD_MIN, KT_BG*lfac)
            if (abs(lat) > MAX_LAT) then
               n_clamp = n_clamp + 1
               call check(error, vmix%kt(i, j, 2) == KD_MIN, &
                          "none+Henyey: poleward of max_lat kt must sit on the kd_min floor")
            else if (KT_BG*lfac < KD_MIN) then
               n_floor = n_floor + 1
               call check(error, vmix%kt(i, j, 2) == KD_MIN, &
                          "none+Henyey: equatorial column must sit on the kd_min floor")
            else
               n_factor = n_factor + 1
               call check(error, abs(vmix%kt(i, j, 2) - want) <= 1.0e-12_wp*want, &
                          "none+Henyey: kt must equal kt_bg*L(phi) in the factor-bound regime")
            end if
            if (allocated(error)) exit checks
            do k = 2, NZ
               call check(error, vmix%kt(i, j, k) == vmix%kt(i, j, 2), &
                          "none+Henyey: kt must be uniform down the column")
               if (allocated(error)) exit checks
               call check(error, vmix%kv(i, j, k) == KV_BG, &
                          "none+Henyey: kv floor is the unscaled kv_bg with prandtl off")
               if (allocated(error)) exit checks
            end do
         end do
         call check(error, n_floor >= 1 .and. n_factor >= 2 .and. n_clamp >= 1, &
                    "test setup: all three Henyey regimes (floor, factor, clamp) must occur")
         if (allocated(error)) exit checks
         ! The factor-bound value must differ from the unscaled kt_bg, so a
         ! missing Henyey branch cannot pass.
         j = NGHOST + 2
         call check(error, abs(vmix%kt(NGHOST + 2, j, 2) - KT_BG) > 1.0e-6_wp, &
                    "test setup: Henyey factor at the probe row is ~1, assertion would not discriminate")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call vmix%destroy(); call ms%destroy()
   end subroutine test_none_henyey

   subroutine test_prandtl(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer :: i, j, k
      checks: block
         ! ---- scalar background, prandtl=1, kv_bg=0 (the OM3 setting) ----
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         call make_spherical_metrics(metrics, grid, 0.0_wp, -8.0_wp, 2.0_wp, 12.0_wp, 6.371e6_wp)
         vmix%kv_bg = 0.0_wp
         vmix%bkgnd_prandtl = 1.0_wp
         call run_step(grid, metrics, ms, vmix)
         do k = 2, NZ
            call check(error, all(vmix%kv(:, :, k) == vmix%kt(:, :, k)) .and. &
                       all(vmix%kv(:, :, k) == KT_BG), &
                       "prandtl=1, kv_bg=0, scalar: kv must equal kt_bg exactly (no shear term)")
            if (allocated(error)) exit checks
         end do
         call vmix%destroy(); call ms%destroy()

         ! ---- Henyey: kv follows the latitude-dependent tracer background ----
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         vmix%kv_bg = 0.0_wp
         vmix%bkgnd_prandtl = 1.0_wp
         vmix%bkgnd_henyey = .true.
         vmix%bkgnd_kd_min = 5.0e-6_wp
         vmix%bkgnd_henyey_max_lat = 40.0_wp
         call run_step(grid, metrics, ms, vmix)
         do j = NGHOST + 1, grid%ny_total - NGHOST
            i = NGHOST + 2
            call check(error, vmix%kv(i, j, 2) == vmix%kt(i, j, 2), &
                       "prandtl=1 + Henyey: kv must equal the latitude-dependent kt exactly")
            if (allocated(error)) exit checks
         end do
         call check(error, abs(vmix%kv(NGHOST + 2, NGHOST + 2, 2) - KT_BG) > 1.0e-6_wp, &
                    "prandtl=1 + Henyey: kv must differ from the unscaled kt_bg (Henyey carried into momentum)")
         if (allocated(error)) exit checks
         call vmix%destroy(); call ms%destroy()

         ! ---- general: kv = prandtl*kt + kv_bg ----
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         ! `vmix` is a re-used object: switch the Henyey scheme of the previous
         ! case off explicitly.
         vmix%bkgnd_henyey = .false.
         vmix%kv_bg = 1.0e-5_wp
         vmix%bkgnd_prandtl = 2.0_wp
         call run_step(grid, metrics, ms, vmix)
         ! Tolerance, not equality: the kernel's prandtl*kt + kv_bg may be FMA-contracted.
         call check(error, all(abs(vmix%kv(:, :, 2:NZ) - (2.0_wp*KT_BG + 1.0e-5_wp)) < 1.0e-18_wp), &
                    "prandtl=2, kv_bg=1e-5: kv must be 2*kt_bg + kv_bg")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call vmix%destroy(); call ms%destroy()
   end subroutine test_prandtl

   subroutine ks_pipeline(du, kv_col, kt_col, kv_bg_in, kt_bg_in)
      !! Production stage order under interior_closure="none" -- reset,
      !! kappa-shear merge, heat/salt split, assemble -- on a uniformly
      !! stratified column with ONE velocity jump `du` between layers
      !! KJ and KJ+1 (`du = 0` is the resting control).  Returns the
      !! assembled column at a probe cell.
      real(wp), intent(in) :: du, kv_bg_in, kt_bg_in
      real(wp), intent(out) :: kv_col(NZK + 1), kt_col(NZK + 1)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_kappa_shear_t) :: ks
      real(wp), parameter :: RHO0 = 1025.0_wp, ALPHA = RHO0*2.0e-4_wp, DT = 1800.0_wp
      real(wp) :: dtdz
      integer :: k

      call grid%init(4, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZK
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZK)
      vmix%use_closure = .true.
      vmix%interior_closure = VMIX_INTERIOR_NONE
      vmix%pp81_nu_bg = kv_bg_in
      vmix%pp81_kappa_bg = kt_bg_in
      call vmix%seed_backgrounds()
      call ks%init(grid, nz_ml=NZK)
      ks%enable = .true.
      ks%eos%variant = EOS_VARIANT_LINEAR
      ks%eos%alpha_T = ALPHA
      ks%eos%beta_S = RHO0*7.6e-4_wp
      ks%eos%rho0 = RHO0
      ks%rho0 = RHO0
      call ks%set_f_centre(grid, 1.0e-4_wp, 0.0_wp, 0.0_wp)

      dtdz = N2*RHO0/(GRAVITY*ALPHA)
      ms%h_layer = DZ
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZK
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (20.0_wp + dtdz*DZ*real(k - NZK, wp))*DZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*DZ
         if (k > KJ) ms%u_face_x_layer(:, :, k) = du
      end do

      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vmix)
      call vmix%enter_data()
      !$acc enter data copyin(ks)
      call ks%enter_data()
      call vmix_seed_interior_none(grid, vmix, ms)
      call kappa_shear_compute(grid, ks, ms, DT)
      call kappa_shear_merge_into_kv_kt(ks, grid%nx_total, grid%ny_total, NZK + 1, vmix%kv, vmix%kt)
      call vmix_split_kd_heat_salt(grid, vmix, ms)
      call vmix_assemble(grid, vmix, ms)
      !$acc update self(vmix%kv, vmix%kt)
      call ks%exit_data()
      !$acc exit data delete(ks)
      call vmix%exit_data()
      !$acc exit data delete(vmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      kv_col = vmix%kv(grid%nx_total/2, grid%ny_total/2, :)
      kt_col = vmix%kt(grid%nx_total/2, grid%ny_total/2, :)
      call ks%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine ks_pipeline

   subroutine test_none_kappa_shear(error)
      !! Kappa-shear (Jackson, Hallberg & Legg 2008) on top of the "none"
      !! interior.  Resting control: Ri = infinity everywhere, so kv/kt are
      !! EXACTLY the background.  One 0.2 m/s jump (Ri = N^2 dz^2/du^2 =
      !! 0.0019 << 0.25) at interface KJ+1: the mixing peaks THERE and falls
      !! off monotonically on both sides.  (The prognostic TKE diffuses into
      !! neighbouring interfaces, so the added mixing is not confined to the
      !! single unstable interface -- that is the scheme, not the interior
      !! closure; the control column proves nothing else leaks in.)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kv0(NZK + 1), kt0(NZK + 1), kv1(NZK + 1), kt1(NZK + 1)
      real(wp) :: kvs(NZK + 1), kts(NZK + 1)
      integer :: k
      checks: block
         call check(error, N2*DZ*DZ/(0.2_wp**2) < 0.25_wp, "test setup: jump interface must have Ri < 0.25")
         if (allocated(error)) exit checks
         call ks_pipeline(0.0_wp, kv0, kt0, KV_BG, KT_BG)
         do k = 2, NZK
            call check(error, kv0(k) == KV_BG .and. kt0(k) == KT_BG, &
                       "none+kappa-shear control: a resting column must be exactly the background")
            if (allocated(error)) exit checks
         end do
         call ks_pipeline(0.2_wp, kv1, kt1, KV_BG, KT_BG)
         call ks_pipeline(0.2_wp, kvs, kts, 0.0_wp, 0.0_wp)
         call check(error, kv1(KJ + 1) > 10.0_wp*KV_BG .and. kt1(KJ + 1) > 10.0_wp*KT_BG, &
                    "none+kappa-shear: kv/kt at the Ri<0.25 interface must far exceed the background")
         if (allocated(error)) exit checks
         ! ADDITIVE, not max(): background + shear.  kts/kvs is the same column
         ! run with a zero background, i.e. the pure kappa-shear diffusivity.
         call check(error, kts(KJ + 1) > 10.0_wp*KT_BG, "test setup: shear kt must dominate the background")
         if (allocated(error)) exit checks
         do k = 2, NZK
            call check(error, abs(kt1(k) - (KT_BG + kts(k))) <= 1.0e-15_wp*kt1(k), &
                       "none+kappa-shear: kt must be kt_bg + kt_shear (additive), not max()")
            if (allocated(error)) exit checks
            call check(error, abs(kv1(k) - (KV_BG + kvs(k))) <= 1.0e-15_wp*kv1(k), &
                       "none+kappa-shear: kv must be kv_bg + kv_shear (additive), not max()")
            if (allocated(error)) exit checks
         end do
         call check(error, maxloc(kv1(2:NZK), 1) + 1 == KJ + 1, &
                    "none+kappa-shear: the mixing peak must sit at the unstable interface")
         if (allocated(error)) exit checks
         do k = 2, NZK
            if (k == KJ + 1) cycle
            call check(error, kv1(k) < kv1(KJ + 1) .and. kt1(k) < kt1(KJ + 1), &
                       "none+kappa-shear: every other interface must mix less than the peak")
            if (allocated(error)) exit checks
            call check(error, kv1(k) >= KV_BG .and. kt1(k) >= KT_BG, &
                       "none+kappa-shear: mixing must never fall below the background")
            if (allocated(error)) exit checks
         end do
         call check(error, kv1(1) == 0.0_wp .and. kv1(NZK + 1) == 0.0_wp, &
                    "none+kappa-shear: closed-BC boundary interfaces must stay zero")
      end block checks
   end subroutine test_none_kappa_shear

   subroutine test_seed_henyey(error)
      !! The seed ITSELF (before any floor) is the Henyey background at
      !! three latitudes: equatorial (kd_min floor), mid-latitude (factor)
      !! and poleward of max_lat (kd_min floor), kv = prandtl kt + kv_bg,
      !! boundaries zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: KD_MIN = 5.0e-6_wp, N0 = 20.0_wp, MAXL = 40.0_wp
      real(wp), parameter :: LATS(3) = [1.0_wp, 45.0_wp, 60.0_wp]
      real(wp), allocatable :: geolat(:, :)
      real(wp) :: want
      integer :: m
      checks: block
         call make_sheared_column(grid, ms, vmix, VMIX_INTERIOR_NONE)
         vmix%bkgnd_henyey = .true.
         vmix%bkgnd_kd_min = KD_MIN
         vmix%bkgnd_henyey_n0_2omega = N0
         vmix%bkgnd_henyey_max_lat = MAXL
         vmix%bkgnd_prandtl = 2.0_wp
         vmix%kv_bg = 1.0e-5_wp
         vmix%kt = -1.0_wp; vmix%kv = -1.0_wp
         allocate (geolat(grid%nx_total, grid%ny_total))
         geolat = 1.0_wp
         geolat(:, 4:5) = 45.0_wp
         geolat(:, 6:) = 60.0_wp
         !$acc enter data copyin(ms, vmix, geolat)
         call ms%enter_data()
         call vmix%enter_data()
         call vmix_seed_interior_none(grid, vmix, ms, geolat=geolat)
         !$acc update self(vmix%kv, vmix%kt)
         call vmix%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vmix, geolat)
         do m = 1, 3
            want = max(KD_MIN, KT_BG*henyey_lat_factor_impl(LATS(m), N0, MAXL))
            if (m == 3) want = KD_MIN
            call check(error, abs(vmix%kt(3, merge(1, merge(4, 6, m == 2), m == 1), 2) - want) <= 1.0e-14_wp*want, &
                       "seed kt must be the Henyey background (floor, factor, clamp)")
            if (allocated(error)) exit checks
            call check(error, abs(vmix%kv(3, merge(1, merge(4, 6, m == 2), m == 1), 2) &
                                  - (2.0_wp*want + 1.0e-5_wp)) <= 1.0e-14_wp*want, &
                       "seed kv must be prandtl*kt_background + kv_bg")
            if (allocated(error)) exit checks
         end do
         call check(error, all(vmix%kt(:, :, 1) == 0.0_wp) .and. all(vmix%kt(:, :, NZ + 1) == 0.0_wp) &
                    .and. all(vmix%kv(:, :, 1) == 0.0_wp) .and. all(vmix%kv(:, :, NZ + 1) == 0.0_wp), &
                    "seed must leave the boundary interfaces at zero")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_seed_henyey

   subroutine test_parse_and_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vmix_t) :: vmix
      checks: block
         call check(error, parse_vmix_interior_closure("pp81") == VMIX_INTERIOR_PP81, "pp81 must parse")
         if (allocated(error)) exit checks
         call check(error, parse_vmix_interior_closure("none") == VMIX_INTERIOR_NONE, "none must parse")
         if (allocated(error)) exit checks
         call check(error, parse_vmix_interior_closure("large94") < 0, &
                    "reserved large94 must be refused")
         if (allocated(error)) exit checks
         call check(error, parse_vmix_interior_closure("cvmix") < 0, "reserved cvmix must be refused")
         if (allocated(error)) exit checks
         call check(error, parse_vmix_interior_closure("PP81") < 0 .and. &
                    parse_vmix_interior_closure("") < 0 .and. &
                    parse_vmix_interior_closure("bogus") < 0, &
                    "unknown spellings must be refused")
         if (allocated(error)) exit checks
         call check(error, vmix_interior_closure_is_implemented(VMIX_INTERIOR_NONE) .and. &
                    vmix_interior_closure_is_implemented(VMIX_INTERIOR_PP81) .and. &
                    .not. vmix_interior_closure_is_implemented(VMIX_INTERIOR_LARGE94) .and. &
                    .not. vmix_interior_closure_is_implemented(VMIX_INTERIOR_CVMIX), &
                    "only pp81 and none are implemented")
         if (allocated(error)) exit checks
         call check(error, vmix%interior_closure == VMIX_INTERIOR_PP81, "default interior closure must be pp81")
         if (allocated(error)) exit checks
         call check(error, vmix%bkgnd_prandtl <= 0.0_wp, "bkgnd_prandtl must default to off (<= 0)")
      end block checks
   end subroutine test_parse_and_defaults

end module test_ocean_vmix_no_pp81
