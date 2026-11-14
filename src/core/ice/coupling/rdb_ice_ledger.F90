!! Host mirror of the sea-ice conservation ledger (roundabout-original).
module rdb_ice_ledger
   !! `ice_ledger_mirror_t` pulls, at a diagnostic cadence, the global sums of the
   !! conserved ice quantities and of the accumulated ocean/atmosphere fluxes, and
   !! reports the per-step closure residual of each line of MASTER_PLAN C.5.  It
   !! holds scalars only; the sums are small `reduce(+)` kernels over the OWNED
   !! cells (`rdb_ice_ledger_sums`), so the only host<-device traffic is a
   !! handful of scalars, and only when `open`/`close` are called.  It never
   !! changes the state or the fluxes.
   !!
   !! **Quantities** (area-weighted totals over the owned cells, cell areas `A`
   !! supplied by the caller):
   !! `M = sum A*(rhoi*vice + rhos*vsno)` [kg], `Ssalt = sum A*rhoi*vice*Sref*1e-3`
   !! [kg], `E = sum A*sum_n( sum_k qicen*vicen/nilyr + sum_k qsnon*vsnon/nslyr )`
   !! [J].
   !!
   !! **Lines** (Icepack signs, see `rdb_ice_flux_state`; the fluxes are the cell
   !! accumulators, so the start-of-step `aicen_init` weight is already inside
   !! them and is not applied again).  With `evap > 0` = condensation (mass
   !! gain), `fsens`, `flat` positive downward, `flwout` negative and the
   !! `flwout` merge already net of the reflected incoming longwave:
   !! - water: `M1 - M0 = W_ext + dt*sum A*(evap - fresh)` (`fresh` already
   !!   contains `evap`: pure sublimation, `evap < 0`, has `fresh = 0` and loses
   !!   mass);
   !! - salt:  `S1 - S0 = -dt*sum A*fsalt`;
   !! - heat:  `E1 - E0 = H_ext + dt*sum A*(fswabs - fswthru + flwout + fsens - fhocn)`.
   !! `flat` is NOT in the heat line: the latent heat leaves with the vapour mass
   !! (upstream `einp = (fsurfn - flatn + fswint - fhocnn - fsnow*Lfresh)*dt`,
   !! where `fsurfn` contains `flatn`, so it cancels; `fswabs = fswsfc + fswint +
   !! fswthru`).  `W_ext` and `H_ext` [kg, J] are the external sources the caller
   !! integrates over the step:
   !! `W_ext = dt*sum A*(fsnow*aice_init + frain*aice)` (snowfall on ice plus the
   !! rain that passes through `fresh`);
   !! `H_ext = dt*sum A*(aice_init*flw - Lfresh*fsnow*aice_init) + E_frz`, where
   !! `aice_init*flw` is the incoming longwave weighted by the start-of-step area
   !! (the merge stores `flwout - (1-emis)*flw` per area, so `sum_n aicen_init*
   !! (emis*flw + flwoutn) = flwout_cell + aice_init*flw`), `-Lfresh*fsnow*aice_init`
   !! is the falling-snow enthalpy and `E_frz` the bank energy spent by new ice
   !! (`frzmlt_pos`; the ocean already gave it up).  Where the plan's C.5 text
   !! differs (a second `aice_init` weight on the atmosphere terms, `flat` in
   !! the heat line) upstream wins.
   !! The residual of a line is `LHS - RHS` (signed); its scale is the sum of the
   !! absolute values of every term, and `rel_residual` is `|res|/scale`.  There
   !! are no dropped-flux terms: this port conserves, so a non-zero residual is a
   !! bug, never an accounted loss.
   !!
   !! **Cumulative interface sums.**  `cum_fresh/fsalt/fhocn/fswthru`
   !! accumulate `dt*sum A*flux` over every closed step: the left-hand side of the
   !! interface closure (what the ocean must have received).
   !!
   !! **Reproducibility / multi-rank.**  The sums are `reduce(+)` (round-off
   !! accurate, not order-reproducible) over the local tile.  `open` takes the rank
   !! count; with more than one rank the ledger is INVALID (`close` leaves
   !! `valid = .false.`, NaN residuals, `closed` false) because the tile sums would
   !! silently under-count.  The multi-rank PR replaces `ice_ledger_reduce` (the
   !! single collective hook) with an all-reduce of reproducing (EFP) sums and
   !! drops the refusal.  There is no MPI here.
   !!
   !! **Non-finite input** propagates through the sums; `finite` is then
   !! `.false.`, every residual is NaN and `closed` is `.false.`.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t, ice_slots_agree
   use rdb_ice_ledger_sums, only: ice_ledger_mass_salt_sums_impl, ice_ledger_energy_sum_impl, &
                                  ice_ledger_ocean_flux_sums_impl, ice_ledger_atm_flux_sums_impl
   implicit none
   private

   public :: ice_ledger_mirror_t
   public :: ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT

   integer, parameter :: ICE_LEDGER_WATER = 1
      !! Line selector: water.
   integer, parameter :: ICE_LEDGER_SALT = 2
      !! Line selector: salt.
   integer, parameter :: ICE_LEDGER_HEAT = 3
      !! Line selector: heat.
   integer, parameter :: NSUM = 6
      !! Flux integrals pulled at `close`.

   type :: ice_ledger_mirror_t
      !! Scalar ledger state (see the module docstring).
      logical :: is_open = .false.
         !! `open` has been called and `close` has not.
      logical :: valid = .false.
         !! The last `close` produced residuals (it followed an `open`).
      logical :: finite = .true.
         !! Every total of the last `close` was finite.
      real(wp) :: dt = 0.0_wp
         !! Step length [s] given to `open`.
      real(wp) :: mass0 = 0.0_wp
         !! `M` at `open` [kg].
      real(wp) :: salt0 = 0.0_wp
         !! `Ssalt` at `open` [kg].
      real(wp) :: energy0 = 0.0_wp
         !! `E` at `open` [J].
      real(wp) :: d_mass = 0.0_wp
         !! `M1 - M0` of the last closed step [kg].
      real(wp) :: d_salt = 0.0_wp
         !! `Ssalt1 - Ssalt0` [kg].
      real(wp) :: d_energy = 0.0_wp
         !! `E1 - E0` [J].
      real(wp) :: i_fresh = 0.0_wp
         !! `dt*sum A*fresh` of the last step [kg].
      real(wp) :: i_evap = 0.0_wp
         !! `dt*sum A*evap` [kg].
      real(wp) :: i_fsalt = 0.0_wp
         !! `dt*sum A*fsalt` [kg].
      real(wp) :: i_fhocn = 0.0_wp
         !! `dt*sum A*fhocn` [J].
      real(wp) :: i_fswthru = 0.0_wp
         !! `dt*sum A*fswthru` [J].
      real(wp) :: i_atm = 0.0_wp
         !! `dt*sum A*(fswabs + flwout + fsens)` [J].
      real(wp) :: res(3) = 0.0_wp
         !! Signed residual `LHS - RHS` of each line (`ICE_LEDGER_*`).
      real(wp) :: scale(3) = 0.0_wp
         !! Sum of the absolute values of every term of each line.
      real(wp) :: cum_fresh = 0.0_wp
         !! Cumulative `dt*sum A*fresh` over closed steps [kg].
      real(wp) :: cum_fsalt = 0.0_wp
         !! Cumulative `dt*sum A*fsalt` [kg].
      real(wp) :: cum_fhocn = 0.0_wp
         !! Cumulative `dt*sum A*fhocn` [J].
      real(wp) :: cum_fswthru = 0.0_wp
         !! Cumulative `dt*sum A*fswthru` [J].
      logical :: multirank = .false.
         !! `open` was given more than one rank: the ledger is refused.
      integer :: nsteps = 0
         !! Number of closed steps.
   contains
      procedure, non_overridable :: open => ledger_open
      procedure, non_overridable :: close => ledger_close
      procedure, non_overridable :: rel_residual => ledger_rel_residual
      procedure, non_overridable :: closed => ledger_closed
      procedure, non_overridable :: reset => ledger_reset
   end type ice_ledger_mirror_t

contains

   subroutine ice_ledger_reduce(v)
      !! Collective hook: sum `v` over the ranks.  Single-rank: identity.
      real(wp), intent(inout) :: v(:)
      if (size(v) < 0) v = 0.0_wp   ! single-rank identity; a multi-rank PR all-reduces here
   end subroutine ice_ledger_reduce

   subroutine state_totals(cice, area, mass, salt, energy)
      !! Local-then-global totals of `M`, `Ssalt`, `E` of the slot.
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(in) :: area(cice%nx, cice%ny)
      real(wp), intent(out) :: mass, salt, energy
      integer :: is, ie, js, je
      real(wp) :: v(3)

      is = cice%nghost + 1
      ie = cice%nghost + cice%nx_phys
      js = cice%nghost + 1
      je = cice%nghost + cice%ny_phys
      call ice_ledger_mass_salt_sums_impl(cice%nx, cice%ny, cice%ncat, is, ie, js, je, area, &
                                          cice%vicen, cice%vsnon, cice%p%c%rhoi, cice%p%c%rhos, &
                                          cice%p%thermo%ice_ref_salinity, v(1), v(2))
      call ice_ledger_energy_sum_impl(cice%nx, cice%ny, cice%ncat, cice%nilyr, cice%nslyr, &
                                      is, ie, js, je, area, cice%vicen, cice%vsnon, &
                                      cice%qicen, cice%qsnon, v(3))
      call ice_ledger_reduce(v)
      mass = v(1)
      salt = v(2)
      energy = v(3)
   end subroutine state_totals

   subroutine ledger_open(this, cice, area, dt, nranks)
      !! Snapshot `M`, `Ssalt`, `E` at the start of an ice step of length `dt` [s]
      !! run on `nranks` ranks.  `area(nx,ny)` are the cell areas [m2]; the state
      !! (and, on an offload build, `area`) must be device-mapped.  `nranks > 1`
      !! marks the ledger invalid (see the module docstring).
      class(ice_ledger_mirror_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      real(wp), intent(in) :: area(cice%nx, cice%ny)
      real(wp), intent(in) :: dt
      integer, intent(in) :: nranks

      this%dt = dt
      this%multirank = nranks > 1
      call state_totals(cice, area, this%mass0, this%salt0, this%energy0)
      this%is_open = .true.
   end subroutine ledger_open

   subroutine ledger_close(this, cice, flux, area, w_ext, h_ext)
      !! Close the step: pull the end-of-step totals and the flux integrals and
      !! evaluate the three residual lines.  `w_ext` [kg] and `h_ext` [J] are the
      !! external sources of the module docstring (default 0).  A `close` without
      !! an `open` leaves `valid = .false.`.
      class(ice_ledger_mirror_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_flux_t), intent(in) :: flux
      real(wp), intent(in) :: area(cice%nx, cice%ny)
      real(wp), intent(in), optional :: w_ext, h_ext

      integer :: is, ie, js, je
      real(wp) :: m1, s1, e1, we, he, dt
      real(wp) :: v(NSUM)

      this%valid = .false.
      this%finite = .false.
      this%res = ieee_value(1.0_wp, ieee_quiet_nan)
      if (.not. this%is_open) return
      this%is_open = .false.
      if (this%multirank .or. .not. ice_slots_agree(flux, cice)) return
      this%finite = .true.
      we = 0.0_wp
      he = 0.0_wp
      if (present(w_ext)) we = w_ext
      if (present(h_ext)) he = h_ext
      dt = this%dt
      is = cice%nghost + 1
      ie = cice%nghost + cice%nx_phys
      js = cice%nghost + 1
      je = cice%nghost + cice%ny_phys

      call state_totals(cice, area, m1, s1, e1)
      call ice_ledger_ocean_flux_sums_impl(cice%nx, cice%ny, is, ie, js, je, area, &
                                           flux%fresh, flux%fsalt, flux%fhocn, v(1), v(2), v(3))
      call ice_ledger_atm_flux_sums_impl(cice%nx, cice%ny, is, ie, js, je, area, &
                                         flux%fswthru, flux%evap, flux%fswabs, flux%flwout, &
                                         flux%fsens, v(4), v(5), v(6))
      call ice_ledger_reduce(v)

      this%d_mass = m1 - this%mass0
      this%d_salt = s1 - this%salt0
      this%d_energy = e1 - this%energy0
      this%i_fresh = dt*v(1)
      this%i_fsalt = dt*v(2)
      this%i_fhocn = dt*v(3)
      this%i_fswthru = dt*v(4)
      this%i_evap = dt*v(5)
      this%i_atm = dt*v(6)

      this%res(ICE_LEDGER_WATER) = this%d_mass - (we + (this%i_evap - this%i_fresh))
      this%scale(ICE_LEDGER_WATER) = abs(this%d_mass) + abs(we) + abs(this%i_evap) + abs(this%i_fresh)
      this%res(ICE_LEDGER_SALT) = this%d_salt + this%i_fsalt
      this%scale(ICE_LEDGER_SALT) = abs(this%d_salt) + abs(this%i_fsalt)
      this%res(ICE_LEDGER_HEAT) = this%d_energy - (he + ((this%i_atm - this%i_fswthru) - this%i_fhocn))
      this%scale(ICE_LEDGER_HEAT) = abs(this%d_energy) + abs(he) + abs(this%i_atm) &
                                    + abs(this%i_fhocn) + abs(this%i_fswthru)

      this%finite = all(ieee_is_finite([this%res, this%scale]))
      if (.not. this%finite) this%res = ieee_value(1.0_wp, ieee_quiet_nan)
      this%cum_fresh = this%cum_fresh + this%i_fresh
      this%cum_fsalt = this%cum_fsalt + this%i_fsalt
      this%cum_fhocn = this%cum_fhocn + this%i_fhocn
      this%cum_fswthru = this%cum_fswthru + this%i_fswthru
      this%nsteps = this%nsteps + 1
      this%valid = .true.
   end subroutine ledger_close

   pure function ledger_rel_residual(this, line) result(r)
      !! `|res|/max(scale, tiny)` of `line` (`ICE_LEDGER_*`) for the last closed
      !! step; NaN when the last close was non-finite or never happened.
      class(ice_ledger_mirror_t), intent(in) :: this
      integer, intent(in) :: line
      real(wp) :: r

      r = ieee_value(1.0_wp, ieee_quiet_nan)
      if (.not. this%valid .or. .not. this%finite) return
      if (line < 1 .or. line > 3) return
      r = abs(this%res(line))/max(this%scale(line), tiny(1.0_wp))
   end function ledger_rel_residual

   pure function ledger_closed(this, tol) result(ok)
      !! True when the last step closed and every line's relative residual is a
      !! finite number `<= tol`.
      class(ice_ledger_mirror_t), intent(in) :: this
      real(wp), intent(in) :: tol
      logical :: ok
      integer :: l
      real(wp) :: r

      ok = this%valid .and. this%finite
      if (.not. ok) return
      do l = 1, 3
         r = this%rel_residual(l)
         if (.not. ieee_is_finite(r)) then
            ok = .false.
         else if (r > tol) then
            ok = .false.
         end if
      end do
   end function ledger_closed

   subroutine ledger_reset(this)
      !! Clear the cumulative sums and the step counter (restart of a run).
      class(ice_ledger_mirror_t), intent(inout) :: this
      select type (this)
      type is (ice_ledger_mirror_t)
         this = ice_ledger_mirror_t()
      end select
   end subroutine ledger_reset

end module rdb_ice_ledger
