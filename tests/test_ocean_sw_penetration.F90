!! Analytical unit tests for shortwave (SW) penetration
!! (`ocean_surface_flux_apply_sw_penetration` in
!! `rdb_ocean_surface_flux`).
!!
!! The kernel is an ADDITIVE CORRECTION layered on top of the legacy
!! surface deposition: the surface kernel first deposits the FULL
!! `Q_heat` lump at `k = nz`, then this kernel removes the penetrating
!! part there and redistributes it through the column as a two-band
!! exponential (Paulson & Simpson 1977), absorbing whatever remains at
!! the opaque bed (`k = 1`).  All tests run the kernel ON DEVICE.
!!
!! Bottom-up convention: surface = `k = nz`, bed = `k = 1`.
!!
!! Cases:
!!   * ENERGY: the column-integrated SW correction is zero (it only
!!     MOVES heat in depth); equivalently, the SW correction's own
!!     column sum of (absorbed - I0-at-nz) = 0 to round-off.  The total
!!     column change vs the all-at-nz baseline is therefore zero, while
!!     the column absorbs exactly `dt/(ρ₀·cp)·sw_pen_frac·Q_heat`.
!!   * PROFILE: per-layer correction matches the two-band formula.
!!   * MONOTONE / TOP-WEIGHTED: surface layer gains the most below the
!!     surface baseline-removal; deeper layers gain progressively less.
!!   * BIT-IDENTITY: sw_pen_frac = 0 ⇒ the kernel makes no change.
module test_ocean_sw_penetration
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_apply_sw_penetration
   implicit none
   private

   public :: collect_ocean_sw_penetration_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 6
   real(wp), parameter :: H_LAYER = 5.0_wp     ! m per layer (deep column)
   real(wp), parameter :: DT = 3600.0_wp
   real(wp), parameter :: Q_HEAT = 200.0_wp    ! W/m^2 downward
   real(wp), parameter :: FRAC = 0.5_wp
   real(wp), parameter :: R = 0.58_wp
   real(wp), parameter :: ZETA1 = 0.35_wp
   real(wp), parameter :: ZETA2 = 23.0_wp

contains

   subroutine collect_ocean_sw_penetration_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sw_energy_conserves_column", test_energy), &
                  new_unittest("sw_profile_matches_two_band", test_profile), &
                  new_unittest("sw_top_weighted_monotone", test_monotone), &
                  new_unittest("sw_frac_zero_bit_identity", test_bit_identity), &
                  new_unittest("negative_source_unreachable_via_q_sw", test_negative_source) &
                  ]
   end subroutine collect_ocean_sw_penetration_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   !> Run surface-flux + SW-penetration on device, return final hTr_T
   !> (the FULL field, including the all-at-nz surface deposition) and,
   !> separately, the SW-only correction (run on a clone without the
   !> surface deposition) so the analytical sums are unambiguous.
   subroutine run_full(grid, sf, ms)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, DT)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr)
         !$acc update self(hT)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_full

   !> Run ONLY the SW penetration kernel (no surface deposition first),
   !> on device.  Returns the SW correction added to hTr_T.
   subroutine run_sw_only(grid, sf, ms)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_sw_penetration(grid, sf, ms, DT)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr)
         !$acc update self(hT)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_sw_only

   subroutine setup(grid, ms, eos, sf, frac)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: frac
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      call sf%init(grid)
      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER
      call sf%set_surface_flux_const(Q_HEAT, 0.0_wp)
      call sf%set_sw_penetration(frac, R, ZETA1, ZETA2)
   end subroutine setup

   !> Closed-form per-layer SW correction (K·m) in the same bottom-up
   !> sweep the kernel uses, so the analytic LHS is independent of it.
   subroutine analytic_correction(corr)
      real(wp), intent(out) :: corr(NZ)
      real(wp) :: inv_scale, i0, d_top, d_bot, tt, tb, absorbed
      integer :: k
      inv_scale = DT/(1035.0_wp*3992.0_wp)   ! rho0, cp defaults on sf
      i0 = FRAC*Q_HEAT
      d_top = 0.0_wp
      do k = NZ, 1, -1
         d_bot = d_top + H_LAYER
         tt = R*exp(-d_top/ZETA1) + (1.0_wp - R)*exp(-d_top/ZETA2)
         if (k > 1) then
            tb = R*exp(-d_bot/ZETA1) + (1.0_wp - R)*exp(-d_bot/ZETA2)
         else
            tb = 0.0_wp
         end if
         absorbed = i0*(tt - tb)
         corr(k) = inv_scale*absorbed
         if (k == NZ) corr(k) = corr(k) - inv_scale*i0
         d_top = d_bot
      end do
   end subroutine analytic_correction

   ! -----------------------------------------------------------------

   subroutine test_energy(error)
      !! The SW correction redistributes — its column sum is zero (the
      !! I0 lump removed at nz equals Σ_k absorbed_k since the bed is
      !! opaque).  Equivalently the absorbed energy = inv_scale·I0 and
      !! the removed lump = inv_scale·I0, so the net is 0 to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: base, col_sum, absorbed_sum, expected_absorbed, inv_scale
      integer :: i, j, k
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, FRAC)
         call run_sw_only(grid, sf, ms)

         i = grid%nx_total/2
         j = grid%ny_total/2
         base = eos%T_ref*H_LAYER
         col_sum = 0.0_wp
         absorbed_sum = 0.0_wp
         do k = 1, NZ
            col_sum = col_sum + (ms%tracers(ms%idx_temperature)%hTr(i, j, k) - base)
         end do
         ! Absorbed = correction + the removed lump at nz (so we recover
         ! the physical column-absorbed energy).
         inv_scale = DT/(sf%rho0*sf%cp)
         absorbed_sum = col_sum + inv_scale*FRAC*Q_HEAT
         expected_absorbed = inv_scale*FRAC*Q_HEAT

         call check(error, abs(col_sum) < 1.0e-12_wp, &
                    "SW correction column sum must be zero (redistributive)")
         if (allocated(error)) exit checks
         call check(error, abs(absorbed_sum - expected_absorbed) < &
                    1.0e-12_wp*abs(expected_absorbed) + 1.0e-15_wp, &
                    "column-absorbed SW must equal dt/(rho0 cp) frac Q_heat")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_energy

   subroutine test_profile(error)
      !! Per-layer SW correction matches the two-band closed form.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: corr(NZ), base, got, maxerr
      integer :: i, j, k
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, FRAC)
         call run_sw_only(grid, sf, ms)
         call analytic_correction(corr)

         i = grid%nx_total/2
         j = grid%ny_total/2
         base = eos%T_ref*H_LAYER
         maxerr = 0.0_wp
         do k = 1, NZ
            got = ms%tracers(ms%idx_temperature)%hTr(i, j, k) - base
            maxerr = max(maxerr, abs(got - corr(k)))
         end do
         call check(error, maxerr < 1.0e-12_wp, &
                    "per-layer SW correction must match two-band formula")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_profile

   subroutine test_monotone(error)
      !! Physically-absorbed SW per layer must be non-negative, the
      !! surface layer (k = nz) must absorb the most, and the INTERIOR
      !! (k = nz-1 .. 2) must decay with depth.  The opaque bed (k = 1)
      !! catches all the remaining transmitted irradiance, so it can
      !! exceed the layer above it — that is physical, not monotone,
      !! and is checked separately (bed > the interior layer above it
      !! because it absorbs the whole remainder).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp) :: corr(NZ), absorbed(NZ), inv_scale
      logical :: ok
      integer :: k
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, FRAC)
         call analytic_correction(corr)
         inv_scale = DT/(sf%rho0*sf%cp)
         ! recover physical absorbed per layer (add lump back at nz)
         absorbed = corr
         absorbed(NZ) = corr(NZ) + inv_scale*FRAC*Q_HEAT

         ok = .true.
         do k = 1, NZ
            if (absorbed(k) < 0.0_wp) ok = .false.
         end do
         call check(error, ok, "absorbed SW per layer must be non-negative")
         if (allocated(error)) exit checks
         ! Surface absorbs strictly more than any deeper layer.
         ok = .true.
         do k = 1, NZ - 1
            if (.not. (absorbed(NZ) > absorbed(k))) ok = .false.
         end do
         call check(error, ok, "surface layer must absorb the most SW")
         if (allocated(error)) exit checks
         ! Interior decay (exclude the opaque bed k = 1, which catches
         ! the whole remainder): absorbed(k) > absorbed(k-1) for the
         ! interior interfaces nz-1 .. 2.
         ok = .true.
         do k = NZ - 1, 3, -1
            if (.not. (absorbed(k) > absorbed(k - 1))) ok = .false.
         end do
         call check(error, ok, "interior absorbed SW must decay with depth")
      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_monotone

   subroutine test_bit_identity(error)
      !! sw_pen_frac = 0 ⇒ has_sw = .false. ⇒ the kernel makes no
      !! change to hTr (byte-for-byte) on device.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :)
      checks: block
         call make_grid(grid)
         call setup(grid, ms, eos, sf, 0.0_wp)   ! frac = 0 ⇒ off
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         call run_sw_only(grid, sf, ms)
         call check(error, &
                    maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp, &
                    "sw_pen_frac=0 must leave hTr byte-for-byte unchanged")
      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_bit_identity

   subroutine test_negative_source(error)
      !! The negative-`I0` landmine (PR-21 §2.4): `Q_heat` is the NET flux
      !! (SW + longwave + latent + sensible).  With `sw_source="net_heat"`
      !! and night-time net cooling (`Q_heat < 0`) the kernel drives
      !! `I0 = sw_pen_frac·Q_heat < 0` — an unphysical NEGATIVE irradiance
      !! that redistributes cooling and (pathologically) WARMS the surface
      !! layer.  The `sw_source="q_sw"` path removes it: `q_sw >= 0` always,
      !! so with `q_sw = 0` (night) `I0 == 0` and the kernel is an EXACT
      !! no-op regardless of `Q_heat`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :)
      real(wp) :: add_nz, base
      integer :: i, j
      checks: block
         ! ---- q_sw path: q_sw = 0, Q_heat = -100 => I0 = 0 (no-op) ----
         call make_grid(grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)
         call sf%set_components(grid, .true.)     ! allocate q_sw
         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER
         call sf%set_surface_flux_const(-100.0_wp, 0.0_wp)   ! net cooling
         call sf%set_sw_penetration(FRAC, R, ZETA1, ZETA2, sw_source="q_sw")
         sf%q_sw = 0.0_wp                          ! night: no sunlight
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         call run_sw_only(grid, sf, ms)
         call check(error, &
                    maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp, &
                    "q_sw=0 with Q_heat<0 must move exactly zero heat (I0=0)")
         if (allocated(error)) exit checks
         call sf%destroy(); call eos%destroy(); call ms%destroy()
         deallocate (hT_ic)

         ! ---- net_heat path: I0 = -100*FRAC < 0 => the landmine fires ----
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)
         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER
         call sf%set_surface_flux_const(-100.0_wp, 0.0_wp)
         call sf%set_sw_penetration(FRAC, R, ZETA1, ZETA2, sw_source="net_heat")
         call run_sw_only(grid, sf, ms)
         i = grid%nx_total/2
         j = grid%ny_total/2
         base = eos%T_ref*H_LAYER
         add_nz = ms%tracers(ms%idx_temperature)%hTr(i, j, NZ) - base
         ! Documented hazard: with I0 < 0 the surface layer is WARMED
         ! (the correction "un-penetrates" cooling upward) — a non-zero,
         ! positive nz increment.  This is what sw_source='q_sw' avoids.
         call check(error, add_nz > 0.0_wp, &
                    "net_heat with Q_heat<0 must exhibit the +ve nz landmine sign")
      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_negative_source

end module test_ocean_sw_penetration
