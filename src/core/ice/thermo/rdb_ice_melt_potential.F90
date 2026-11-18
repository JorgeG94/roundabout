!! Ocean-side melt-potential provider for the sea ice (`frzmlt`).
module rdb_ice_melt_potential
   !! One signed heat-flux field per cell for the sea ice, `frzmlt` [W m-2]:
   !!
   !!   frzmlt > 0   freezing potential: the ocean already made frazil, the
   !!                standing frazil BANK `frazil_heat` [J m-2] divided by the
   !!                step;
   !!   frzmlt < 0   melt potential: heat in the upper ocean above the
   !!                freezing point, available to melt ice.
   !!
   !! The two signs are EXCLUSIVE per cell and the bank wins: a positive bank
   !! gives `frzmlt = bank/dt` however warm the water below is (that heat is
   !! not available to melt ice that the same cell is being asked to grow).
   !!
   !! Melt potential (the ocean-model formulation of Adcroft et al. 2019,
   !! "The GFDL Global Ocean and Sea Ice Model OM4.0", GMD 12, 3167-3204, and
   !! of the OM3/CESM coupling, where the ocean hands the sea ice the heat
   !! needed to bring its mixed layer back to the freezing point):
   !!
   !!   MP = rho * Cp * sum_k dh_k * (T_k - T_f(S_k, p = 0))     [J m-2]
   !!
   !! summed downward from the first live layer at or below the surface
   !! (`k = nz`; vanished layers skipped) over exactly
   !! `min(hfreeze, h_bl)` metres: a layer that straddles that depth is
   !! counted in part (`dh_k = min(h_k, remaining)`), and the sum stops at
   !! the bed.  `h_bl` is the surface boundary-layer depth of the active
   !! scheme: KPP `vmix%bl_depth`, EPBL `epbl%mld`; a zero depth gives zero
   !! melt potential (a negative one is clamped to zero).  ONLY with no
   !! boundary-layer scheme (PP81 only; `h_bl` absent) is the fixed `hfreeze`
   !! used.  `MP` is stored raw (it can be negative: a cold column);
   !! `frzmlt = -MP/dt` only where `MP > 0` and the bank is not positive.
   !!
   !! Freezing point and constants.  `T_f` is `eos_freezing_point(eos, S, 0)`
   !! of the OCEAN EOS (the liquidus the frazil kernel clamps to, the
   !! configured `tfreeze_set`, `cons_abs` conventions included), evaluated
   !! at pressure ZERO exactly as the frazil kernel's surface clamp does --
   !! no second liquidus is written.  `rho * Cp` is `RHO_WATER * eos%c_p`,
   !! the product that filled the bank, so the bank and the melt potential
   !! are in the same joules.
   !!
   !! Non-finite input (T, S or h of a layer the sum touches, `h_bl`, the
   !! bank) makes BOTH outputs of that cell NaN, never a silent zero; the
   !! consumer's finiteness check (the ice error word of the import step)
   !! turns that into its fatal bit.  NaN is produced from a scalar the host
   !! shim builds with `ieee_value`: device code never calls it.
   !!
   !! CONTRACT FOR THE COUPLING PR (so heat is counted ONCE).
   !!  1. This module is a PROVIDER.  It changes no ocean state and debits
   !!     nothing.
   !!  2. The ocean is debited only through the ice's `fhocn` heat flux (what
   !!     the ice actually used to melt or grow ice), never by `-frzmlt`: the
   !!     melt potential is an upper bound on what the ice MAY take, not a
   !!     flux that happened.
   !!  3. Frazil heat was already paid by the ocean when the bank was filled
   !!     (the surface layer was warmed to `T_f` then).  New ice formed from
   !!     the bank must never put that heat into `fhocn`; the bank is zeroed
   !!     by the energy spent.
   !!  4. Two freezing points exist and are NOT interchangeable: the ocean
   !!     EOS liquidus (`eos_freezing_point`, here and in the frazil clamp)
   !!     and the ice's linear `-depressT * S` (`Tf` on the ice import).  The
   !!     ocean side (this kernel, the bank) uses the first; the ice's basal
   !!     and lateral heat balance uses the second.  A cell whose SST sits
   !!     between them has a small melt potential here and no thermodynamic
   !!     drive there; the ice side must clamp, not reconcile.
   !!  5. Provider cadence: `frzmlt = bank/dt` is a rate only if `dt` is the
   !!     window over which the bank accumulated.  The bank is filled every
   !!     outer step and spent (zeroed) by the frazil uptake at thermo
   !!     cadence, so the engine calls this once per thermo step, with the
   !!     thermo-window step, BEFORE the uptake.  A caller with another
   !!     bank-spending cadence must pass that window.
   !!  6. Ghost cells of `frzmlt` / `melt_potential` are NOT filled (owned
   !!     cells only).  A consumer that averages `frzmlt` to faces or reads a
   !!     neighbour must halo-exchange it first.
   !!  7. The two output arrays exist on the ice slot only when
   !!     `&ocean_ice_nml hfreeze > 0`; it is refused under an ice-shelf
   !!     cavity (`&ocean_cavity_dyn_nml enable`), where melt potential must
   !!     not be offered to sea ice.
   !!
   !! Off by default (`&ocean_ice_nml hfreeze = -1`): the kernel is then never
   !! called and nothing changes.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED, H_DIV_EPS, GRAVITY, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   implicit none
   private

   public :: ice_melt_potential

contains

   pure subroutine ice_melt_potential(grid, eos, ms, frazil_heat, hfreeze, dt, &
                                      melt_potential, frzmlt, h_bl)
      !! Outer shim: dereference the tracer registry on the host, build the
      !! NaN scalar, forward bare explicit-shape arrays to the device kernel.
      !! No-op when T or S is unregistered.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: frazil_heat(:, :)
         !! Frazil bank [J m-2] (`ocean_sea_ice_t%frazil_heat`).
      real(wp), intent(in) :: hfreeze
         !! Depth cap [m], > 0.
      real(wp), intent(in) :: dt
         !! Step [s], > 0.
      real(wp), intent(inout) :: melt_potential(:, :)
         !! Out: raw melt potential [J m-2].
      real(wp), intent(inout) :: frzmlt(:, :)
         !! Out: signed `frzmlt` [W m-2].
      real(wp), intent(in), optional :: h_bl(:, :)
         !! Boundary-layer depth [m]; absent => fixed `hfreeze`.

      integer :: idx_T, idx_S
      real(wp) :: nan

      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (idx_T <= 0 .or. idx_S <= 0) return
      nan = ieee_value(1.0_wp, ieee_quiet_nan)

      if (present(h_bl)) then
         call ice_melt_potential_impl(ms%tracers(idx_T)%c, ms%tracers(idx_S)%c, &
                                      ms%h_layer, ms%wet_mask, frazil_heat, h_bl, &
                                      .true., hfreeze, dt, nan, eos, &
                                      melt_potential, frzmlt, grid%nghost, &
                                      ms%nz_ml, grid%nx_total, grid%ny_total)
      else
         call ice_melt_potential_impl(ms%tracers(idx_T)%c, ms%tracers(idx_S)%c, &
                                      ms%h_layer, ms%wet_mask, frazil_heat, frazil_heat, &
                                      .false., hfreeze, dt, nan, eos, &
                                      melt_potential, frzmlt, grid%nghost, &
                                      ms%nz_ml, grid%nx_total, grid%ny_total)
      end if
   end subroutine ice_melt_potential

   pure subroutine ice_melt_potential_impl(c_T, c_S, h_layer, wet_mask, bank, h_bl, &
                                           use_hbl, hfreeze, dt, nan, eos, &
                                           melt_potential, frzmlt, nghost, nz, nx, ny)
      !! Device kernel over PHYSICAL cells (ghosts untouched).  One thread
      !! per column; the layer loop is sequential inside it.  Dry cells and
      !! land write zeros.
      integer, intent(in) :: nghost, nz, nx, ny
      real(wp), intent(in) :: c_T(nx, ny, nz)
      real(wp), intent(in) :: c_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: bank(nx, ny)
      real(wp), intent(in) :: h_bl(nx, ny)
      logical, intent(in) :: use_hbl
      real(wp), intent(in) :: hfreeze, dt, nan
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: melt_potential(nx, ny)
      real(wp), intent(inout) :: frzmlt(nx, ny)

      integer :: i, j, k, i_lo, i_hi, j_lo, j_hi
      real(wp) :: depth, rem, h, dh, t_k, s_k, mp, b, rcp
      logical :: bad

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      rcp = RHO_WATER*eos%c_p

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(k, depth, rem, h, dh, t_k, s_k, mp, b, bad)
         mp = 0.0_wp
         b = 0.0_wp
         melt_potential(i, j) = 0.0_wp
         frzmlt(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            bad = .false.
            depth = hfreeze
            if (use_hbl) then
               if (.not. ieee_is_finite(h_bl(i, j))) then
                  bad = .true.
               else
                  depth = min(hfreeze, max(h_bl(i, j), 0.0_wp))
               end if
            end if
            b = bank(i, j)
            if (.not. ieee_is_finite(b)) bad = .true.
            rem = depth
            ! Scan from the surface layer down; a vanished layer is skipped
            ! (the first-live rule), but a NON-FINITE thickness is an error,
            ! not a skipped layer, so the scan checks it itself.
            do k = nz, 1, -1
               if (rem <= 0.0_wp) exit
               h = h_layer(i, j, k)
               if (.not. ieee_is_finite(h)) then
                  bad = .true.
                  exit
               end if
               if (.not. rdb_vl_is_live(h)) cycle
               t_k = c_T(i, j, k)
               s_k = c_S(i, j, k)
               if (.not. (ieee_is_finite(t_k) .and. ieee_is_finite(s_k))) then
                  bad = .true.
                  exit
               end if
               dh = min(h, rem)
               mp = mp + dh*(t_k - eos_freezing_point(eos, s_k, 0.0_wp))
               rem = rem - dh
            end do
            mp = rcp*mp
            if (bad) then
               melt_potential(i, j) = nan
               frzmlt(i, j) = nan
            else
               melt_potential(i, j) = mp
               if (b > 0.0_wp) then
                  frzmlt(i, j) = b/dt
               else if (mp > 0.0_wp) then
                  frzmlt(i, j) = -mp/dt
               end if
            end if
         end if
      end do
   end subroutine ice_melt_potential_impl

#include "rdb_vanished_layer.inc"

end module rdb_ice_melt_potential
