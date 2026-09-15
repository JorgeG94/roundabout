!! Unit tests for the C7 depth-varying background mixing profiles
!! (Bryan & Lewis 1979 + the Henyey, Wright & Flatte 1986 latitude
!! factor, in the constant-N0 simplification of Harrison & Hallberg
!! 2008) in `rdb_ocean_vmix`.
!!
!! `bkgnd_profile` (Bryan-Lewis) and `bkgnd_henyey` are MUTUALLY
!! EXCLUSIVE background schemes, matching the reference formulation
!! (which FATALs when a second background scheme is selected):
!!
!!   * `bkgnd_profile` — `vmix_assemble` replaces the SCALAR kt_bg/ks_bg
!!     additive floor with a per-interface Bryan-Lewis depth profile
!!     (kv floor = bkgnd_prandtl * kd_bg).
!!   * `bkgnd_henyey`  — the SCALAR kt_bg/ks_bg floors are scaled by a
!!     latitude-only factor `L(phi)` (`henyey_lat_factor_impl`) and
!!     floored at `bkgnd_kd_min`:  max(kd_min, kt_bg*L(phi)).  The
!!     momentum floor `kv_bg` is deliberately NOT scaled.
!!
!! Both default off ⇒ the scalar floor path is used verbatim
!! (bit-identical).
!!
!! MUTATION NOTE.  Every Henyey case below deliberately avoids phi = 30
!! deg as the *discriminating* latitude: `L(30 deg) == 1` by
!! construction, so a test that only probes 30 deg cannot tell "factor
!! applied" from "factor never applied" and lets the whole branch be
!! deleted while staying green.  The load-bearing latitudes here are
!! 45 deg, 90 deg and their southern mirrors, and the knob values are
!! deliberately NON-default so a kernel that hardcoded them fails.
!!
!! Cases:
!!   * Bryan-Lewis SHAPE: a column's kd_bg follows the atan between the
!!     surface (Kd_sfc) and deep (Kd_deep) asymptotes, hits the midpoint
!!     (Kd_sfc+Kd_deep)/2 at z=z0, and is monotonic with depth.
!!   * DEFAULT-OFF bit-identity: with bkgnd_profile off the assembly
!!     background == the scalar kv_bg/kt_bg/ks_bg path, exactly.
!!   * Profile-ON floor: an interface below the local Bryan-Lewis floor is
!!     raised to it; kv floor = bkgnd_prandtl * kd_bg.
!!   * Henyey knob defaults: N0_2Omega=20, max_lat=95, master switch off.
!!   * Henyey ANALYTICAL: the closed-form factor at phi = 0/30/45/90 deg
!!     against a Python `math.acosh` oracle (not merely "small"/"less
!!     than 1" -- the pole value is > 1 under the default N0_2Omega),
!!     then the same latitudes through the clip kernel where the
!!     equatorial column lands on the `kd_min` FLOOR rather than zero.
!!   * Henyey HEMISPHERIC SYMMETRY: L(-phi) == L(+phi) bitwise at three
!!     latitudes, and the poleward clamp fires in the SOUTH too -- the
!!     `abs()` on both `sin(lat)` and the clamp comparison is the only
!!     thing standing between the Southern Ocean and a NEGATIVE floor
!!     that `max(kt, floor)` would swallow without a trace.
!!   * Henyey TRACER FLOORS: with kt_bg /= ks_bg the SAME factor scales
!!     both tracer floors -- checked with the two backgrounds an order of
!!     magnitude apart, so scaling only one cannot reproduce both -- and
!!     the momentum floor kv_bg is left unscaled (the documented scope).
!!   * Henyey KD_MIN floor: `max(kd_min, kt_bg*L)` with the MOM6-shaped
!!     `0.01*kt_bg` default resolved from the negative sentinel, and the
!!     explicit-knob path.  The equatorial column is the discriminator.
!!   * Henyey max_lat clamp: poleward of the cutoff the factor collapses
!!     to a tiny POSITIVE floor (~1.2e-9), NOT to the exact zero the
!!     equator gives.  The oracle is DERIVED from the public
!!     `HENYEY_MIN_SINLAT` rather than hardcoded.
!!   * Henyey DEVICE RESIDENCY + slot knobs: `geolat` mapped at 45 deg
!!     with a NON-default `bkgnd_henyey_n0_2omega` on the slot, host copy
!!     poisoned on GPU builds -- an unmapped/stale-host read would give
!!     the poisoned answer, and a kernel ignoring the slot knob would
!!     give the default-N0 answer (pattern: test_ocean_diag_reduce.F90).
!!   * Henyey SLOT max_lat, southern latitude, through `vmix_assemble`.
!!   * Henyey METRICS THREADING: a full `ocean_dyn_step` on a SPHERICAL
!!     grid -- the only case that proves `vmix_apply_in_stage` forwards
!!     `metrics%geolatT` and not `geolonT` / `geolatBu`, which on this
!!     grid give demonstrably different factors.
module test_ocean_bkgnd_mixing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, PI
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_assemble, vmix_bkgnd_fill_impl, &
                             vmix_assemble_clip_henyey_impl, henyey_lat_factor_impl, &
                             vmix_resolve_kd_min, bkgnd_henyey_conflicts_profile, &
                             HENYEY_MIN_SINLAT, HENYEY_KD_MIN_FRAC
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step
   implicit none
   private

   public :: collect_ocean_bkgnd_mixing_tests

   integer, parameter :: NGHOST = 2

   ! Host-poison sentinel for the device-residency test.  Any latitude
   ! value clearly distinct from the mapped reference works; 0 deg is the
   ! natural choice (its Henyey factor is exactly 0, maximally
   ! distinguishable from any non-equatorial factor).
   real(wp), parameter :: POISON_LAT = 0.0_wp

   ! Python `math` oracles for L(phi) at the DEFAULT n0_2omega = 20:
   !   L(45)  = |sin45|*acosh(20/|sin45|) / (0.5*acosh(40))
   real(wp), parameter :: L45_N20 = 1.3023092477420493_wp
   real(wp), parameter :: L90_N20 = 1.6834153338247386_wp
   ! ... and at the NON-default n0_2omega = 5 used to prove the slot knob
   ! actually reaches the kernel.
   real(wp), parameter :: L45_N5 = 1.2492726512763956_wp

contains

   subroutine collect_ocean_bkgnd_mixing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bkgnd_bryan_lewis_shape", test_bryan_lewis_shape), &
                  new_unittest("bkgnd_default_off_bit_identity", test_default_off_bit_identity), &
                  new_unittest("bkgnd_profile_floor_applies", test_profile_floor_applies), &
                  new_unittest("bkgnd_henyey_defaults", test_henyey_defaults), &
                  new_unittest("bkgnd_henyey_factor_analytical", test_henyey_factor_analytical), &
                  new_unittest("bkgnd_henyey_hemispheric_symmetry", test_henyey_hemispheric_symmetry), &
                  new_unittest("bkgnd_henyey_scales_tracer_floors", test_henyey_scales_tracer_floors), &
                  new_unittest("bkgnd_henyey_kd_min_floor", test_henyey_kd_min_floor), &
                  new_unittest("bkgnd_schemes_mutually_exclusive", test_schemes_mutually_exclusive), &
                  new_unittest("bkgnd_henyey_max_lat_clamp", test_henyey_max_lat_clamp), &
                  new_unittest("bkgnd_henyey_device_resident", test_henyey_device_resident), &
                  new_unittest("bkgnd_henyey_slot_clamp_south", test_henyey_slot_clamp_south), &
                  new_unittest("bkgnd_henyey_metrics_threading", test_henyey_metrics_threading) &
                  ]
   end subroutine collect_ocean_bkgnd_mixing_tests

   subroutine make_state(grid, ms, vmix, nz, layer_h)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: nz
      real(wp), intent(in) :: layer_h
      call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      call vmix%init(grid, nz_ml=nz)
      ms%h_layer = layer_h
   end subroutine make_state

   pure function rel_err(got, want) result(e)
      !! Relative error, safe for the ~1e-14 magnitudes the poleward-clamp
      !! cases produce (an absolute tolerance there would be vacuous).
      real(wp), intent(in) :: got, want
      real(wp) :: e
      e = abs(got - want)/max(abs(want), tiny(1.0_wp))
   end function rel_err

   ! -----------------------------------------------------------------

   subroutine test_bryan_lewis_shape(error)
      !! Drive `vmix_bkgnd_fill_impl` directly on a uniform-thickness
      !! column and check the atan shape: midpoint at z=z0, bounded by the
      !! asymptotes, monotonically increasing with depth.  This routine has
      !! NO Henyey arguments at all any more (the latitude-scaled variant
      !! is a separate `_impl`), so it is also the bit-identity proof for
      !! the plain path: assertions/tolerances unchanged from before Henyey.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 40
      ! Layer thickness chosen so an interface lands exactly at z0.
      real(wp), parameter :: KD_SFC = 1.0e-5_wp, KD_DEEP = 1.3e-4_wp
      real(wp), parameter :: Z0 = 1000.0_wp, DELTA = 222.0_wp, HZ = 100.0_wp
      real(wp) :: mid, depth, expected
      integer :: i, j, k, k_at_z0
      checks: block
         call make_state(grid, ms, vmix, NZ, HZ)
         vmix%bkgnd_kd_sfc = KD_SFC
         vmix%bkgnd_kd_deep = KD_DEEP
         vmix%bkgnd_z0 = Z0
         vmix%bkgnd_delta = DELTA

         i = 5; j = 4
         ! Host-side fill (the kernel runs on host data here; the device
         ! round-trip is exercised in the floor test below).
         call vmix_bkgnd_fill_impl(grid%nx_total, grid%ny_total, NZ + 1, &
                                   vmix%kd_bg, ms%h_layer, KD_SFC, KD_DEEP, Z0, DELTA)

         ! Interface k has depth = sum of layers k..NZ = (NZ-k+1)*HZ.
         ! z0=1000, HZ=100 ⇒ interface where depth=1000 is k = NZ-10+1 = 31.
         k_at_z0 = NZ - nint(Z0/HZ) + 1
         depth = real(NZ - k_at_z0 + 1, wp)*HZ
         call check(error, abs(depth - Z0) < 1.0e-9_wp, "test setup: no interface at z0")
         if (allocated(error)) exit checks

         ! (1) midpoint at z=z0.
         mid = 0.5_wp*(KD_SFC + KD_DEEP)
         call check(error, abs(vmix%kd_bg(i, j, k_at_z0) - mid) < 1.0e-12_wp, &
                    "Bryan-Lewis not at midpoint at z=z0")
         if (allocated(error)) exit checks

         ! (2) bounded by the asymptotes at every interior interface.
         do k = 2, NZ
            call check(error, vmix%kd_bg(i, j, k) >= KD_SFC - 1.0e-12_wp .and. &
                       vmix%kd_bg(i, j, k) <= KD_DEEP + 1.0e-12_wp, &
                       "Bryan-Lewis value outside [Kd_sfc, Kd_deep]")
            if (allocated(error)) exit checks
         end do

         ! (3) monotonic increase with DEPTH (k decreasing ⇒ deeper ⇒ larger).
         do k = NZ, 3, -1
            call check(error, vmix%kd_bg(i, j, k - 1) >= vmix%kd_bg(i, j, k) - 1.0e-15_wp, &
                       "Bryan-Lewis not monotonic with depth")
            if (allocated(error)) exit checks
         end do

         ! (4) exact formula at one interior interface (k = NZ-5, depth=600).
         depth = real(NZ - (NZ - 5) + 1, wp)*HZ
         expected = KD_SFC + (KD_DEEP - KD_SFC)* &
                    (0.5_wp + atan((depth - Z0)/DELTA)/PI)
         call check(error, abs(vmix%kd_bg(i, j, NZ - 5) - expected) < 1.0e-15_wp, &
                    "Bryan-Lewis formula mismatch at interior interface")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_bryan_lewis_shape

   subroutine run_assemble(grid, ms, vmix, geolat, poison_geolat)
      !! Device round-trip: enter_data, assemble, pull kv/kt/ks/kd_bg back.
      !! `geolat`, when present, is mapped alongside ms/vmix; when
      !! `poison_geolat` is also `.true.`, its HOST copy is overwritten
      !! AFTER the device copyin on GPU builds only (the residency
      !! discriminator -- see test_ocean_diag_reduce.F90).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      real(wp), intent(inout), optional :: geolat(grid%nx_total, grid%ny_total)
         !! Explicit-shape (not assumed-shape) so the actual is passed by
         !! plain base address with no chance of a compiler-made contiguous
         !! temporary — a temporary would be absent from the device present
         !! table and turn the residency check below into a false pass.
      logical, intent(in), optional :: poison_geolat
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      !$acc update device(vmix%kv, vmix%kt, vmix%ks)
      if (present(geolat)) then
         !$acc enter data copyin(geolat)
#ifdef RDB_GPU_OFFLOAD
         if (present(poison_geolat)) then
            if (poison_geolat) call poison_host_2d(geolat)
         end if
#endif
         call vmix_assemble(grid, vmix, ms, geolat=geolat)
         !$acc exit data delete(geolat)
      else
         call vmix_assemble(grid, vmix, ms)
      end if
      !$acc update self(vmix%kv, vmix%kt, vmix%ks, vmix%kd_bg)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_assemble

   subroutine clip_henyey(grid, vmix, nz, k_bg, kd_min, n0_2omega, max_lat, geolat)
      !! Zero kv/kt/ks, then drive `vmix_assemble_clip_henyey_impl` with
      !! `kt_bg == ks_bg == k_bg` and NO momentum floor (`kv_bg = 0`), so
      !! the assembled kt/ks ARE the Henyey tracer floor
      !! `max(kd_min, k_bg*L(phi))` with nothing else mixed in.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: vmix
      integer, intent(in) :: nz
      real(wp), intent(in) :: k_bg, kd_min, n0_2omega, max_lat
      real(wp), intent(in) :: geolat(grid%nx_total, grid%ny_total)
      vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
      call vmix_assemble_clip_henyey_impl(grid%nx_total, grid%ny_total, nz + 1, &
                                          vmix%kv, vmix%kt, vmix%ks, &
                                          0.0_wp, k_bg, k_bg, huge(1.0_wp), huge(1.0_wp), &
                                          n0_2omega, max_lat, kd_min, geolat)
   end subroutine clip_henyey

   subroutine poison_host_2d(arr)
      !! Overwrite the HOST copy of an already-mapped 2D buffer.  A plain
      !! sequential loop, never `do concurrent` and never array syntax, so
      !! there is no chance of the store being offloaded to the device
      !! (which would poison the very copy under test).
      ! assumed-shape-ok: host-only helper, no `do concurrent` here.
      real(wp), intent(inout) :: arr(:, :)
      integer :: i, j
      do j = 1, size(arr, 2)
         do i = 1, size(arr, 1)
            arr(i, j) = POISON_LAT
         end do
      end do
   end subroutine poison_host_2d

   subroutine test_default_off_bit_identity(error)
      !! With bkgnd_profile off (default), the assembly background floor
      !! is the scalar kv_bg/kt_bg/ks_bg path — bit-for-bit unchanged from
      !! the pre-C7 behaviour.  A field above the scalar floor survives
      !! exactly; a sub-floor value is raised to the SCALAR floor.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 5
      real(wp), allocatable :: kv0(:, :, :), kt0(:, :, :), ks0(:, :, :)
      integer :: k
      checks: block
         call make_state(grid, ms, vmix, NZ, 5.0_wp)
         ! bkgnd_profile defaults .false. — do not touch it.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0e-3_wp + 1.0e-4_wp*real(k, wp)
            vmix%kt(:, :, k) = 2.0e-4_wp + 1.0e-5_wp*real(k, wp)
            vmix%ks(:, :, k) = 3.0e-4_wp + 1.0e-5_wp*real(k, wp)
         end do
         kv0 = vmix%kv; kt0 = vmix%kt; ks0 = vmix%ks

         call run_assemble(grid, ms, vmix)

         call check(error, all(vmix%kv == kv0), "kv changed with bkgnd off")
         if (allocated(error)) exit checks
         call check(error, all(vmix%kt == kt0), "kt changed with bkgnd off")
         if (allocated(error)) exit checks
         call check(error, all(vmix%ks == ks0), "ks changed with bkgnd off")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_default_off_bit_identity

   subroutine test_profile_floor_applies(error)
      !! With bkgnd_profile on (bkgnd_henyey left OFF) and a sub-floor
      !! closure output, kt/ks are raised to the local Bryan-Lewis kd_bg
      !! and kv to prandtl*kd_bg.  Exercises the full device assembly
      !! path; also the default-off bit-identity proof for `bkgnd_henyey`
      !! at the `vmix_assemble` level (unchanged from before Henyey).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 30
      real(wp), parameter :: KD_SFC = 1.0e-5_wp, KD_DEEP = 1.3e-4_wp
      real(wp), parameter :: Z0 = 1000.0_wp, DELTA = 222.0_wp
      real(wp), parameter :: PR = 2.0_wp, HZ = 100.0_wp
      real(wp) :: depth, expected_kd
      integer :: i, j, k, kk
      checks: block
         call make_state(grid, ms, vmix, NZ, HZ)
         vmix%bkgnd_profile = .true.
         vmix%bkgnd_kd_sfc = KD_SFC
         vmix%bkgnd_kd_deep = KD_DEEP
         vmix%bkgnd_z0 = Z0
         vmix%bkgnd_delta = DELTA
         vmix%bkgnd_prandtl = PR
         ! Closure output deliberately well below any background floor.
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         do k = 2, NZ
            vmix%kv(:, :, k) = 1.0e-12_wp
            vmix%kt(:, :, k) = 1.0e-12_wp
            vmix%ks(:, :, k) = 1.0e-12_wp
         end do

         call run_assemble(grid, ms, vmix)

         i = 5; j = 4; kk = NZ - 5
         depth = real(NZ - kk + 1, wp)*HZ
         expected_kd = KD_SFC + (KD_DEEP - KD_SFC)* &
                       (0.5_wp + atan((depth - Z0)/DELTA)/PI)
         call check(error, abs(vmix%kt(i, j, kk) - expected_kd) < 1.0e-15_wp, &
                    "kt not floored to Bryan-Lewis kd_bg")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%ks(i, j, kk) - expected_kd) < 1.0e-15_wp, &
                    "ks not floored to Bryan-Lewis kd_bg")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kv(i, j, kk) - PR*expected_kd) < 1.0e-15_wp, &
                    "kv not floored to prandtl*kd_bg")
         if (allocated(error)) exit checks
         ! Boundary interfaces untouched (closed BC).
         call check(error, abs(vmix%kv(i, j, 1)) < 1.0e-15_wp, "bed interface altered")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%kv(i, j, NZ + 1)) < 1.0e-15_wp, "surface interface altered")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_profile_floor_applies

   subroutine test_henyey_defaults(error)
      !! The knob + its two constants exist on the slot with the shipped
      !! defaults, master switch off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      checks: block
         call make_state(grid, ms, vmix, 4, 5.0_wp)
         call check(error,.not. vmix%bkgnd_henyey, "bkgnd_henyey should default off")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%bkgnd_henyey_n0_2omega - 20.0_wp) < 1.0e-12_wp, &
                    "bkgnd_henyey_n0_2omega default should be 20.0")
         if (allocated(error)) exit checks
         call check(error, abs(vmix%bkgnd_henyey_max_lat - 95.0_wp) < 1.0e-12_wp, &
                    "bkgnd_henyey_max_lat default should be 95.0 (inert clamp)")
         if (allocated(error)) exit checks
         ! `bkgnd_kd_min` ships as a NEGATIVE "unset" sentinel that
         ! `vmix_resolve_kd_min` turns into the reference 0.01*Kd default.
         call check(error, vmix%bkgnd_kd_min < 0.0_wp, &
                    "bkgnd_kd_min default should be the negative unset sentinel")
         if (allocated(error)) exit checks
         call check(error, abs(vmix_resolve_kd_min(vmix%bkgnd_kd_min, vmix%kt_bg) - &
                               HENYEY_KD_MIN_FRAC*vmix%kt_bg) < 1.0e-20_wp, &
                    "the unset sentinel did not resolve to HENYEY_KD_MIN_FRAC*kt_bg")
         if (allocated(error)) exit checks
         ! A non-negative value is taken literally -- kd_min = 0 must stay 0
         ! (a legal way to ask for no floor at all), not fall back to 0.01*Kd.
         call check(error, vmix_resolve_kd_min(0.0_wp, vmix%kt_bg) == 0.0_wp, &
                    "kd_min = 0 was overridden by the sentinel default")
      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_henyey_defaults

   subroutine test_henyey_factor_analytical(error)
      !! Hand-verified closed-form latitude factor (Python `math.acosh`
      !! oracle) at four latitudes, asserted BOTH on the scalar kernel
      !! `henyey_lat_factor_impl` and end-to-end through
      !! `vmix_assemble_clip_henyey_impl` with `kt_bg == 1`, so the
      !! assembled tracer floor is `max(kd_min, L(phi))` directly.
      !!
      !! 45 deg and 90 deg are the load-bearing values: `L(30 deg) == 1`
      !! by construction, so 30 deg alone cannot distinguish "factor
      !! applied" from "factor never applied".
      !!
      !! RE-BASELINED EQUATOR ORACLE.  This case used to assert that the
      !! assembled background at the equator was EXACTLY 0, which encoded a
      !! deliberate deviation: rdb multiplied `L(phi)` straight in with
      !! no minimum-diffusivity floor, so the equatorial and poleward limits
      !! collapsed toward zero instead of the documented "returned to the
      !! MINIMUM diffusivity".  The floor now ships (`max(Kd_min, Kd*L)`,
      !! MOM6 `KD_MIN`), so the assembled equatorial value is `KD_MIN` --
      !! here the explicitly-set 1.0e-7 -- and the old `== 0` oracle would
      !! now be asserting the absence of the floor.  The RAW factor is still
      !! exactly 0 at the equator; that assertion stays, one level down, on
      !! `henyey_lat_factor_impl` itself.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      real(wp), parameter :: K = 1.0_wp
      real(wp), parameter :: KD_MIN = 1.0e-7_wp
      real(wp), allocatable :: geolat(:, :)
      integer :: i, j
      checks: block
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         i = 5; j = 4
         allocate (geolat(grid%nx_total, grid%ny_total), source=0.0_wp)

         ! --- scalar kernel, direct ---
         ! The RAW factor carries no floor -- that lives in the clip kernel.
         call check(error, henyey_lat_factor_impl(0.0_wp, 20.0_wp, 95.0_wp) == 0.0_wp, &
                    "L(0 deg) must be EXACTLY 0 (the raw factor carries no floor)")
         if (allocated(error)) exit checks
         call check(error, abs(henyey_lat_factor_impl(30.0_wp, 20.0_wp, 95.0_wp) - 1.0_wp) &
                    < 1.0e-14_wp, "L(30 deg) must be 1 (the normalisation latitude)")
         if (allocated(error)) exit checks
         call check(error, rel_err(henyey_lat_factor_impl(45.0_wp, 20.0_wp, 95.0_wp), &
                                   L45_N20) < 1.0e-13_wp, &
                    "L(45 deg) mismatched the Python oracle")
         if (allocated(error)) exit checks
         call check(error, rel_err(henyey_lat_factor_impl(90.0_wp, 20.0_wp, 95.0_wp), &
                                   L90_N20) < 1.0e-13_wp, &
                    "L(90 deg) mismatched the Python oracle (it is ~1.68, NOT < 1)")
         if (allocated(error)) exit checks
         ! The knob must actually change the answer -- a kernel that
         ! hardcoded n0_2omega = 20 would return L45_N20 here.
         call check(error, rel_err(henyey_lat_factor_impl(45.0_wp, 5.0_wp, 95.0_wp), &
                                   L45_N5) < 1.0e-13_wp, &
                    "L(45 deg) ignored the n0_2omega argument")
         if (allocated(error)) exit checks

         ! --- through the clip kernel (kt_bg == 1 ⇒ floor == max(kd_min, L)) ---
         ! phi = 0 deg: sin(0) = 0 EXACTLY -> factor EXACTLY 0, not NaN/Inf,
         ! so the assembled floor is KD_MIN itself.
         geolat(i, j) = 0.0_wp
         call clip_henyey(grid, vmix, NZ, K, KD_MIN, 20.0_wp, 95.0_wp, geolat)
         call check(error, vmix%kt(i, j, 2) == KD_MIN, &
                    "at the equator L=0, so the tracer floor must be exactly KD_MIN "// &
                    "(re-baselined from the old no-floor oracle of exactly 0)")
         if (allocated(error)) exit checks
         call check(error, vmix%ks(i, j, 2) == KD_MIN, &
                    "the salinity floor missed KD_MIN at the equator")
         if (allocated(error)) exit checks

         geolat(i, j) = 45.0_wp
         call clip_henyey(grid, vmix, NZ, K, KD_MIN, 20.0_wp, 95.0_wp, geolat)
         call check(error, rel_err(vmix%kt(i, j, 2), L45_N20) < 1.0e-13_wp, &
                    "clip kernel did not apply L(45 deg)")
         if (allocated(error)) exit checks

         geolat(i, j) = 90.0_wp
         call clip_henyey(grid, vmix, NZ, K, KD_MIN, 20.0_wp, 95.0_wp, geolat)
         call check(error, rel_err(vmix%kt(i, j, 2), L90_N20) < 1.0e-13_wp, &
                    "clip kernel did not apply L(90 deg)")
         if (allocated(error)) exit checks
         ! The floor must not be masking the factor at these latitudes:
         ! L(45)*K and L(90)*K are both ~1e7 x KD_MIN.
         call check(error, vmix%kt(i, j, 2) > 1.0e3_wp*KD_MIN, &
                    "test setup: KD_MIN is large enough to mask the latitude factor")
      end block checks
      call vmix%destroy(); call ms%destroy()
      if (allocated(geolat)) deallocate (geolat)
   end subroutine test_henyey_factor_analytical

   subroutine test_henyey_hemispheric_symmetry(error)
      !! The factor is symmetric about the equator: L(-phi) == L(+phi),
      !! bitwise.  Without the `abs()` on `sin(lat)` the Southern
      !! Hemisphere would produce a NEGATIVE kd_bg, which the downstream
      !! `max(kt, kd_bg)` floor swallows silently -- no NaN, no crash,
      !! just no background mixing south of the equator.
      !!
      !! Second half: the poleward clamp compares |lat| against max_lat.
      !! Without that `abs()` the Southern Ocean is never clamped, so
      !! L(-45) with max_lat = 20 would come back as the unclamped 1.302
      !! instead of the ~1.2e-9 floor.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: floor_expected
      checks: block
         call check(error, henyey_lat_factor_impl(-30.0_wp, 20.0_wp, 95.0_wp) == &
                    henyey_lat_factor_impl(30.0_wp, 20.0_wp, 95.0_wp), &
                    "L(-30) /= L(+30): the abs() on sin(lat) is missing")
         if (allocated(error)) exit checks
         call check(error, henyey_lat_factor_impl(-45.0_wp, 20.0_wp, 95.0_wp) == &
                    henyey_lat_factor_impl(45.0_wp, 20.0_wp, 95.0_wp), &
                    "L(-45) /= L(+45): the abs() on sin(lat) is missing")
         if (allocated(error)) exit checks
         call check(error, henyey_lat_factor_impl(-90.0_wp, 20.0_wp, 95.0_wp) == &
                    henyey_lat_factor_impl(90.0_wp, 20.0_wp, 95.0_wp), &
                    "L(-90) /= L(+90): the abs() on sin(lat) is missing")
         if (allocated(error)) exit checks
         ! Positivity everywhere (the sign bug's actual symptom).
         call check(error, henyey_lat_factor_impl(-45.0_wp, 20.0_wp, 95.0_wp) > 0.0_wp, &
                    "L(-45) is not positive: a negative kd_bg reaches the floor")
         if (allocated(error)) exit checks

         ! Southern poleward clamp.
         floor_expected = HENYEY_MIN_SINLAT*acosh(20.0_wp/HENYEY_MIN_SINLAT) &
                          /(0.5_wp*acosh(40.0_wp))
         call check(error, rel_err(henyey_lat_factor_impl(-45.0_wp, 20.0_wp, 20.0_wp), &
                                   floor_expected) < 1.0e-13_wp, &
                    "the poleward clamp did not fire at -45 deg (missing abs() on the cutoff)")
      end block checks
   end subroutine test_henyey_hemispheric_symmetry

   subroutine test_henyey_scales_tracer_floors(error)
      !! The SAME latitude factor must scale BOTH tracer background floors.
      !! `kt_bg` and `ks_bg` are set an order of magnitude apart, so a
      !! kernel that scaled only one of them (or reused one floor for both)
      !! cannot reproduce `kt_bg*L` and `ks_bg*L` simultaneously.
      !!
      !! Second half pins the documented SCOPE: the momentum floor `kv_bg`
      !! is deliberately NOT scaled.  Roundabout's scalar background carries
      !! `kv_bg` as an independent momentum floor rather than
      !! `prandtl*kt_bg`, so there is no single background `Kd` for the
      !! reference `Kv_bkgnd = PRANDTL_BKGND*Kd` tie to reproduce.  Asserted
      !! rather than merely commented, so changing it is a deliberate act
      !! that re-baselines a test.
      !!
      !! 45 deg (not 30 deg, where L == 1 by construction) so "factor never
      !! applied" is distinguishable; `kd_min` is 0 so the floor cannot mask
      !! the factor at either background.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 6
      real(wp), parameter :: KT_BG = 1.0e-5_wp, KS_BG = 1.3e-4_wp
      real(wp), parameter :: KV_BG = 1.0e-4_wp
      real(wp), parameter :: LAT = 45.0_wp
      real(wp), allocatable :: geolat(:, :)
      integer :: i, j, k
      checks: block
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         i = 5; j = 4; k = 3
         allocate (geolat(grid%nx_total, grid%ny_total), source=LAT)

         ! Sanity: the two tracer backgrounds really are far apart,
         ! otherwise "scaled kt, reused it for ks" would stay green.
         call check(error, KS_BG > 5.0_wp*KT_BG, &
                    "test setup: kt_bg and ks_bg are not far enough apart")
         if (allocated(error)) exit checks

         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         call vmix_assemble_clip_henyey_impl(grid%nx_total, grid%ny_total, NZ + 1, &
                                             vmix%kv, vmix%kt, vmix%ks, &
                                             KV_BG, KT_BG, KS_BG, &
                                             huge(1.0_wp), huge(1.0_wp), &
                                             20.0_wp, 95.0_wp, 0.0_wp, geolat)

         call check(error, rel_err(vmix%kt(i, j, k), KT_BG*L45_N20) < 1.0e-13_wp, &
                    "the temperature floor is not kt_bg*L(45 deg)")
         if (allocated(error)) exit checks
         call check(error, rel_err(vmix%ks(i, j, k), KS_BG*L45_N20) < 1.0e-13_wp, &
                    "the salinity floor is not ks_bg*L(45 deg)")
         if (allocated(error)) exit checks
         ! Momentum floor left at the plain scalar kv_bg -- documented scope.
         call check(error, vmix%kv(i, j, k) == KV_BG, &
                    "kv was scaled by the latitude factor: the Henyey scaling is "// &
                    "documented to touch the TRACER floors only")
      end block checks
      call vmix%destroy(); call ms%destroy()
      if (allocated(geolat)) deallocate (geolat)
   end subroutine test_henyey_scales_tracer_floors

   subroutine test_henyey_kd_min_floor(error)
      !! The `max(Kd_min, Kd*L(phi))` minimum-diffusivity floor (MOM6
      !! `KD_MIN`), through the full `vmix_assemble` device path.
      !!
      !! Equatorial latitude: `L(0 deg) == 0` exactly, so the floor is the
      !! ONLY thing keeping the background off zero -- which makes this the
      !! discriminating case for the floor's existence.  Two configurations:
      !! an explicit `bkgnd_kd_min`, and the negative sentinel resolving to
      !! the reference `0.01*kt_bg` default.
      !! The two configurations use SEPARATE slot variables on purpose: a
      !! `destroy` + `init` round trip does NOT restore a derived type's
      !! default component initialisation, so reusing one slot would carry
      !! the explicit `bkgnd_kd_min` into the sentinel case and test nothing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid, grid_d
      type(multilayer_state_t) :: ms, ms_d
      type(ocean_vmix_t) :: vmix     !! explicit-kd_min case
      type(ocean_vmix_t) :: vmix_d   !! shipped-default (sentinel) case
      integer, parameter :: NZ = 4
      real(wp), parameter :: KT_BG = 2.0e-5_wp
      real(wp), parameter :: KD_MIN = 3.0e-7_wp
      real(wp), allocatable :: geolat(:, :), geolat_d(:, :)
      integer :: i, j
      checks: block
         i = 5; j = 4

         ! --- explicit kd_min ---
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         vmix%bkgnd_henyey = .true.
         vmix%kt_bg = KT_BG; vmix%ks_bg = KT_BG
         vmix%bkgnd_kd_min = KD_MIN
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         allocate (geolat(grid%nx_total, grid%ny_total), source=0.0_wp)

         call run_assemble(grid, ms, vmix, geolat=geolat)

         call check(error, vmix%kt(i, j, 2) == KD_MIN, &
                    "the equatorial background is not the explicit bkgnd_kd_min: "// &
                    "either the floor is missing (it would be 0) or kd_min was ignored")
         if (allocated(error)) exit checks
         call check(error, vmix%ks(i, j, 2) == KD_MIN, &
                    "the equatorial salinity background is not bkgnd_kd_min")
         if (allocated(error)) exit checks
         ! It really is a floor and not "kt_bg passed through unscaled".
         call check(error, vmix%kt(i, j, 2) < KT_BG, &
                    "the equatorial background equals kt_bg: the latitude factor "// &
                    "was never applied")
         if (allocated(error)) exit checks

         ! --- negative sentinel ⇒ HENYEY_KD_MIN_FRAC*kt_bg ---
         call make_state(grid_d, ms_d, vmix_d, NZ, 10.0_wp)
         vmix_d%bkgnd_henyey = .true.
         vmix_d%kt_bg = KT_BG; vmix_d%ks_bg = KT_BG
         ! bkgnd_kd_min left at its shipped negative sentinel on purpose.
         vmix_d%kv = 0.0_wp; vmix_d%kt = 0.0_wp; vmix_d%ks = 0.0_wp
         allocate (geolat_d(grid_d%nx_total, grid_d%ny_total), source=0.0_wp)

         call run_assemble(grid_d, ms_d, vmix_d, geolat=geolat_d)

         call check(error, rel_err(vmix_d%kt(i, j, 2), HENYEY_KD_MIN_FRAC*KT_BG) &
                    < 1.0e-13_wp, &
                    "the unset bkgnd_kd_min sentinel did not resolve to "// &
                    "HENYEY_KD_MIN_FRAC*kt_bg at the equator")
      end block checks
      call vmix%destroy(); call ms%destroy()
      call vmix_d%destroy(); call ms_d%destroy()
      if (allocated(geolat)) deallocate (geolat)
      if (allocated(geolat_d)) deallocate (geolat_d)
   end subroutine test_henyey_kd_min_floor

   subroutine test_schemes_mutually_exclusive(error)
      !! Bryan-Lewis and Henyey are MUTUALLY EXCLUSIVE background schemes,
      !! matching the reference code (`check_bkgnd_scheme` FATALs when a
      !! second scheme is selected).
      !!
      !! The rule is asserted on the `bkgnd_henyey_conflicts_profile`
      !! PREDICATE rather than end-to-end, because both enforcement sites
      !! (`validate_config`'s logger + the `configure_ocean_vmix`
      !! `error stop` backstop) terminate the process and are unreachable
      !! from a unit test without a subprocess harness.  Extracting the
      !! rule into a `pure` predicate is what makes it testable at all --
      !! with the condition inlined at the guard, deleting it was a
      !! mutation the whole suite survived.  Both call sites go through
      !! this function, so a wrong rule here is a wrong rule everywhere.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, bkgnd_henyey_conflicts_profile(.true., .true.), &
                    "Henyey + Bryan-Lewis must CONFLICT (they are mutually exclusive "// &
                    "background schemes, not a composition)")
         if (allocated(error)) exit checks
         ! Every other combination is legal -- in particular Henyey ALONE,
         ! which the pre-parity code refused outright.
         call check(error,.not. bkgnd_henyey_conflicts_profile(.true., .false.), &
                    "Henyey alone must be legal (it scales the SCALAR background; "// &
                    "it no longer requires bkgnd_profile)")
         if (allocated(error)) exit checks
         call check(error,.not. bkgnd_henyey_conflicts_profile(.false., .true.), &
                    "Bryan-Lewis alone must be legal")
         if (allocated(error)) exit checks
         call check(error,.not. bkgnd_henyey_conflicts_profile(.false., .false.), &
                    "both off must be legal (the default scalar background)")
      end block checks
   end subroutine test_schemes_mutually_exclusive

   subroutine test_henyey_max_lat_clamp(error)
      !! Setting `henyey_max_lat` below the test latitude collapses the
      !! factor to the tiny POSITIVE floor the `HENYEY_MIN_SINLAT` guard
      !! produces -- which is NOT the exact 0 the equator gives (the
      !! equator keeps its true |sin phi| = 0 in the outer multiplication;
      !! the clamp substitutes the floor there too).  The oracle is
      !! DERIVED from the public `HENYEY_MIN_SINLAT` so it cannot rot if
      !! that constant is ever retuned; its current value is 1.21933e-9.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      real(wp), parameter :: K = 1.0_wp
      real(wp), allocatable :: geolat(:, :)
      real(wp) :: expected
      integer :: i, j
      checks: block
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         i = 5; j = 4
         allocate (geolat(grid%nx_total, grid%ny_total), source=45.0_wp)
         expected = HENYEY_MIN_SINLAT*acosh(20.0_wp/HENYEY_MIN_SINLAT) &
                    /(0.5_wp*acosh(40.0_wp))

         ! kd_min = 0 so the clamped value reaches the assembled floor
         ! unmasked -- this case is about the CLAMP, not about the floor
         ! (which `test_henyey_kd_min_floor` owns).
         call clip_henyey(grid, vmix, NZ, K, 0.0_wp, 20.0_wp, 20.0_wp, geolat)

         call check(error, rel_err(vmix%kt(i, j, 2), expected) < 1.0e-13_wp, &
                    "max_lat clamp did not reset the factor to the min_sinlat floor value")
         if (allocated(error)) exit checks
         ! It is a floor, not a zero -- distinct from the equatorial value.
         call check(error, vmix%kt(i, j, 2) > 0.0_wp, &
                    "the poleward clamp produced 0, not the min_sinlat floor")
      end block checks
      call vmix%destroy(); call ms%destroy()
      if (allocated(geolat)) deallocate (geolat)
   end subroutine test_henyey_max_lat_clamp

   subroutine test_henyey_device_resident(error)
      !! GPU residency discriminator (pattern: test_ocean_diag_reduce.F90)
      !! AND slot-knob wiring, in one case.
      !!
      !! `geolat` is mapped to device at 45 deg, then the HOST copy is
      !! poisoned to 0 deg (factor == 0) on GPU builds only: a kernel that
      !! bound to an implicit per-launch copyin instead of the mapped
      !! device buffer would come back with the equatorial answer.
      !!
      !! `bkgnd_henyey_n0_2omega` is set to a NON-default 5.0 on the slot,
      !! so the expected answer (L45_N5) differs from the default-knob
      !! answer (L45_N20) by ~4% -- a kernel that hardcoded the shipped
      !! default instead of reading the slot fails here.  45 deg (not the
      !! normalisation latitude 30 deg) so "factor never applied" is also
      !! distinguishable.
      !!
      !! `bkgnd_kd_min = 0` so the floor cannot mask any of that; the floor
      !! itself is `test_henyey_kd_min_floor`'s business.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      real(wp), parameter :: KD_CONST = 5.0e-5_wp
      real(wp), allocatable :: geolat(:, :)
      integer :: i, j
      checks: block
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         vmix%bkgnd_henyey = .true.
         vmix%bkgnd_henyey_n0_2omega = 5.0_wp   ! NON-default on purpose
         vmix%kt_bg = KD_CONST; vmix%ks_bg = KD_CONST
         vmix%bkgnd_kd_min = 0.0_wp
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         allocate (geolat(grid%nx_total, grid%ny_total), source=45.0_wp)

         call run_assemble(grid, ms, vmix, geolat=geolat, poison_geolat=.true.)

         i = 5; j = 4
         call check(error, rel_err(vmix%kt(i, j, 2), KD_CONST*L45_N5) < 1.0e-13_wp, &
                    "Henyey kt floor is not KD_CONST*L(45 deg, N0_2Omega=5): either the "// &
                    "poisoned host latitude was read, the slot knob was ignored, or "// &
                    "the factor was never applied")
         if (allocated(error)) exit checks
         ! Explicitly rule out the two silent failure modes.
         call check(error, vmix%kt(i, j, 2) /= 0.0_wp, &
                    "Henyey kt floor came back at the equatorial (poison) value")
         if (allocated(error)) exit checks
         call check(error, rel_err(vmix%kt(i, j, 2), KD_CONST*L45_N20) > 1.0e-3_wp, &
                    "Henyey kt floor matched the DEFAULT n0_2omega: the slot knob is ignored")
      end block checks
      call vmix%destroy(); call ms%destroy()
      if (allocated(geolat)) deallocate (geolat)
   end subroutine test_henyey_device_resident

   subroutine test_henyey_slot_clamp_south(error)
      !! `bkgnd_henyey_max_lat` reaches the kernel from the slot, and the
      !! clamp fires at a SOUTHERN latitude — through the full
      !! `vmix_assemble` device path rather than the bare `_impl`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      integer, parameter :: NZ = 4
      real(wp), parameter :: KD_CONST = 5.0e-5_wp
      real(wp), allocatable :: geolat(:, :)
      real(wp) :: expected
      integer :: i, j
      checks: block
         call make_state(grid, ms, vmix, NZ, 10.0_wp)
         vmix%bkgnd_henyey = .true.
         vmix%bkgnd_henyey_max_lat = 20.0_wp    ! NON-default on purpose
         vmix%kt_bg = KD_CONST; vmix%ks_bg = KD_CONST
         vmix%bkgnd_kd_min = 0.0_wp   ! floor off: this case is about the clamp
         vmix%kv = 0.0_wp; vmix%kt = 0.0_wp; vmix%ks = 0.0_wp
         allocate (geolat(grid%nx_total, grid%ny_total), source=-45.0_wp)
         expected = KD_CONST*HENYEY_MIN_SINLAT*acosh(20.0_wp/HENYEY_MIN_SINLAT) &
                    /(0.5_wp*acosh(40.0_wp))

         call run_assemble(grid, ms, vmix, geolat=geolat)

         i = 5; j = 4
         call check(error, rel_err(vmix%kt(i, j, 2), expected) < 1.0e-12_wp, &
                    "slot bkgnd_henyey_max_lat did not clamp at -45 deg")
      end block checks
      call vmix%destroy(); call ms%destroy()
      if (allocated(geolat)) deallocate (geolat)
   end subroutine test_henyey_slot_clamp_south

   subroutine test_henyey_metrics_threading(error)
      !! END-TO-END: one full `ocean_dyn_step` on a SPHERICAL grid with
      !! Henyey on and the closure output driven below any floor, so the
      !! assembled `kt(i,j,k)` must equal
      !! `KD_CONST * L(metrics%geolatT(i,j))` at every column.
      !!
      !! This is the only case that covers the `metrics%geolatT` argument
      !! threaded through `vmix_apply_in_stage` into `vmix_assemble`.
      !! Everything else in this file calls `vmix_assemble` directly with
      !! a latitude array the test itself supplies, so a production
      !! dispatch that forwarded `geolonT` (same shape) or `geolatBu`
      !! (corner-shaped, sequence-associated) would stay green.  The final
      !! assertion proves the discrimination: on this grid the longitude
      !! field gives a materially different factor at the probe column.
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
      integer, parameter :: NZ = 3
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: KD_CONST = 5.0e-5_wp
      real(wp), parameter :: N0_2OM = 20.0_wp, MAX_LAT = 95.0_wp
      real(wp) :: expected, lon_factor
      integer :: i, j, k
      checks: block
         call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
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
         call dyn%init(grid)

         ! Lon-lat sector: latitudes 10..~30 degN, longitudes 0..~24 degE
         ! -- geolatT, geolonT and geolatBu all hold DIFFERENT numbers at
         ! any given (i, j), which is what makes the check discriminating.
         call make_spherical_metrics(metrics, grid, 0.0_wp, 10.0_wp, &
                                     2.0_wp, 2.0_wp, 6.371e6_wp)

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do

         ! use_closure ON so vmix_apply_in_stage actually reaches
         ! vmix_assemble; KPP left off (no surface-flux slot here).
         vmix%use_closure = .true.
         vmix%use_kpp = .false.
         vmix%bkgnd_henyey = .true.
         vmix%kt_bg = KD_CONST; vmix%ks_bg = KD_CONST
         vmix%bkgnd_kd_min = 0.0_wp   ! floor off: kt == KD_CONST * L(lat)
         vmix%bkgnd_henyey_n0_2omega = N0_2OM
         vmix%bkgnd_henyey_max_lat = MAX_LAT
         ! PP81 interior background must not out-floor the Henyey result:
         ! `vmix_compute_pp81` writes pp81_kappa_bg + kappa0*factor into kt,
         ! and the assembly floor is a max(), so a PP81 background above
         ! KD_CONST*L would mask the very thing under test.
         vmix%pp81_kappa_bg = 0.0_wp
         vmix%pp81_nu_bg = 0.0_wp
         vmix%pp81_nu0 = 0.0_wp

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

         call ocean_dyn_step(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, va, hd, vd, vmix, ms, DT)

         !$acc update self(vmix%kt, vmix%ks)
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

         ! Every physical column must carry its own latitude's factor.
         do j = NGHOST + 1, grid%ny_total - NGHOST
            do i = NGHOST + 1, grid%nx_total - NGHOST
               expected = KD_CONST*henyey_lat_factor_impl(metrics%geolatT(i, j), &
                                                          N0_2OM, MAX_LAT)
               call check(error, rel_err(vmix%kt(i, j, 2), expected) < 1.0e-12_wp, &
                          "kt does not match KD_CONST*L(metrics%geolatT): the "// &
                          "production dispatch is not threading geolatT")
               if (allocated(error)) exit checks
            end do
         end do

         ! Discrimination proof: at the probe column the LONGITUDE field
         ! would give a visibly different factor, so the assertion above
         ! genuinely pins geolatT specifically.
         i = NGHOST + 3; j = NGHOST + 2
         lon_factor = KD_CONST*henyey_lat_factor_impl(metrics%geolonT(i, j), N0_2OM, MAX_LAT)
         call check(error, rel_err(lon_factor, KD_CONST* &
                                   henyey_lat_factor_impl(metrics%geolatT(i, j), &
                                                          N0_2OM, MAX_LAT)) > 1.0e-2_wp, &
                    "test setup: geolonT and geolatT give the same factor here, so the "// &
                    "threading assertion above would not discriminate")
      end block checks
      call destroy_cartesian_metrics(metrics)
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy()
   end subroutine test_henyey_metrics_threading

end module test_ocean_bkgnd_mixing
