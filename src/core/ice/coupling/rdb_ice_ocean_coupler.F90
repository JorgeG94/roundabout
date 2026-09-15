!! Ice -> ocean brine-rejection coupling: refresh the surface salt-flux
!! field from the sea-ice slot's last uptake (SIS2 port, PR 3b).
module rdb_ice_ocean_coupler
   !! Bridges `ocean_sea_ice_t%salt_flux_diag` (filled by
   !! `rdb_ice_frazil_uptake%ice_frazil_uptake`) into
   !! `ocean_surface_flux_t%Q_salt` — the field the REAL production
   !! kernel `ocean_surface_flux_apply_tracers`
   !! (`rdb_ocean_surface_flux.F90`) reads every thermo step.
   !!
   !! **TRAP #2 — sign convention** (resolved against the CODE, not the
   !! (now-fixed) stale module-header docstring of
   !! `rdb_ocean_surface_flux`): `apply_surface_src_2d_impl` does
   !! `hTr_S(i,j,nz) += (dt/rho0) * Q_salt(i,j) * wet_mask(i,j)`, i.e.
   !! `dS/dt = +Q_salt/(rho0*h_top)` — Q_salt is POSITIVE-SALINIFIES.
   !! Brine rejection (freezing seawater leaves excess salt behind) must
   !! RAISE the ocean surface salinity, so `salt_flux_diag` (already
   !! computed with that sign in `ice_frazil_uptake_impl`) is added
   !! directly, unnegated.
   !!
   !! OWNERSHIP: with ice on and `&ocean_forcing_nml enable_components =
   !! .false.` (the default), this coupler OWNS the `Q_salt` field-fill —
   !! any future field-fill path (Area-A3 data-override) composing with
   !! sea ice must add into `salt_flux_diag` or fold in here, not
   !! bypass this call. **With `enable_components = .true.` (PR-12), the
   !! coupler instead OWNS the `salt_flux` / `heat_added` COMPONENTS** —
   !! `ocean_surface_flux_assemble` derives the net `Q_salt`/`Q_heat`
   !! from `Q_salt_const + salt_flux` / `Q_heat_const + heat_added` at
   !! the same point in the step, so the two modes produce identical
   !! numbers (the components-on ice bit-identity gate,
   !! `tests/test_ocean_surface_forcing_type.F90:ice_components_bitident`).
   !!
   !! Q_heat was NOT written by the PR-3b coupler: the frazil latent heat
   !! was already credited to the ocean by the PR-1 surface clamp (see
   !! `rdb_ice_frazil_uptake` module docstring); writing it there would
   !! have double-counted. PR 3c adds `ice_ocean_heat_flux` below, which
   !! refreshes Q_heat from the melt-side `heat_flux_diag` (filled by
   !! `rdb_ice_thermo_driver` — the column's grow/melt exchange, NOT the
   !! frazil bank) — a disjoint energy pathway from the frazil clamp, so
   !! no double-count.  Under `enable_components`, `heat_flux_diag`
   !! lands in `heat_added` instead — MOM6's slot for a net, already-
   !! summed heat term that isn't further decomposable (PR-12 plan §5.4).
   !!
   !! PR 31 adds `ice_ocean_sw_flux`, which OWNS the shortwave the ice
   !! transmits to the ocean (`ocean_sea_ice_t%sw_thru_diag`, W/m^2,
   !! +down).  Because `q_sw` is the shortwave SHARE OF the net heat
   !! (Q_heat CONTAINS q_sw, it is not Q_heat + q_sw), the delivery is
   !! component-mode-dependent and the coupler owns exactly ONE path in
   !! each mode so the energy reaches `Q_heat` once, never twice:
   !!   * components ON  — writes the `q_sw` COMPONENT; the assembler
   !!     sums it into `Q_heat`.  `heat_flux_diag` (-> `heat_added`)
   !!     carries only the NON-shortwave heat, so there is no overlap.
   !!   * components OFF — `q_sw` is unallocated and there is no
   !!     assembler, so the shortwave is ADDED directly into the
   !!     components-off `Q_heat` (which `ice_ocean_heat_flux`
   !!     full-overwrote to `Q_heat_const + heat_flux_diag` immediately
   !!     before, per the driver's mandated order), giving
   !!     `Q_heat_const + heat_flux_diag + sw_thru_diag`.
   !! Either way the shortwave is NEVER folded into `heat_flux_diag`
   !! (`rdb_ice_thermo_driver`), which is what keeps it out of
   !! `heat_added` and prevents a double count.
   use rdb_constants, only: wp
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ice_state, only: ocean_sea_ice_t, ice_cell_concentration_impl
   implicit none
   private

   public :: ice_ocean_brine_flux
   public :: ice_ocean_heat_flux
   public :: ice_ocean_sw_flux
   public :: ice_ocean_stress_flux
   public :: ice_ocean_stress_resume_apply
   public :: ice_ocean_stress_flux_impl
   public :: ice_ocean_stress_cleanup

   ! ---- F6: persistent device scratch for the per-step tau blend, so
   ! `ice_ocean_stress_flux` never per-step allocate + enter/exit-data
   ! churns (the repo forbids per-step device alloc/map).  Lazily sized
   ! by `stress_scratch_ensure`, released by `ice_ocean_stress_cleanup`
   ! (driver teardown, next to `ice_evp_cleanup`). ----
   logical, save :: stress_scratch_ready = .false.
   integer, save :: ss_nx = 0, ss_ny = 0
   real(wp), allocatable, save :: ci_scratch(:, :)
   real(wp), allocatable, save :: mis_scratch(:, :), mice_scratch(:, :)
      !! `mis_scratch`/`mice_scratch` are throwaway outputs the shared
      !! `ice_cell_concentration_impl` requires; only `ci_scratch` feeds
      !! the blend.

contains

   pure subroutine ice_ocean_brine_flux(sf, ice)
      !! Refresh the ocean surface salt-flux field from the ice slot.
      !! **Components off** (default):
      !!   Q_salt(i,j) = Q_salt_const + salt_flux_diag(i,j)
      !! POSITIVE SALINIFIES (TRAP #2 above — matches
      !! `apply_surface_src_2d_impl`'s `hTr_S += dt/rho0*Q_salt`). This
      !! is a FULL OVERWRITE from the configure-time constant plus the
      !! last uptake's rate: no accumulation drift, and a thermo window
      !! with no freezing resets the field back to the background
      !! (`salt_flux_diag == 0` when the bank was empty or the cell was
      !! dry/land/vanished — see `ice_frazil_uptake_impl`'s unconditional
      !! diag zeroing).
      !! **Components on** (`sf%use_components`, PR-12): write the ice's
      !! own COMPONENT instead — `sf%salt_flux(i,j) = salt_flux_diag(i,j)`
      !! (no `Q_salt_const` term: `ocean_surface_flux_assemble` adds it).
      !! Sets `has_salt` either way (§13 item 6 of the PR-12 plan: any
      !! filler that can produce non-zero net salt must set the latch).
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_sea_ice_t), intent(in) :: ice

      sf%has_salt = .true.
      if (sf%use_components) then
         call ice_ocean_brine_flux_components_impl(sf%salt_flux, ice%salt_flux_diag, &
                                                   ice%nx_total, ice%ny_total)
      else
         call ice_ocean_brine_flux_impl(sf%Q_salt, ice%salt_flux_diag, &
                                        sf%Q_salt_const, ice%nx_total, ice%ny_total)
      end if
   end subroutine ice_ocean_brine_flux

   pure subroutine ice_ocean_brine_flux_components_impl(salt_flux, salt_flux_diag, nx, ny)
      !! Components-on branch: full-array overwrite of the `salt_flux`
      !! COMPONENT (not `Q_salt`) — same sign convention, no `_const`
      !! term (the assembler adds it).  Explicit-shape + decl-order.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: salt_flux(nx, ny)
      real(wp), intent(in) :: salt_flux_diag(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         salt_flux(i, j) = salt_flux_diag(i, j)
      end do
   end subroutine ice_ocean_brine_flux_components_impl

   pure subroutine ice_ocean_brine_flux_impl(Q_salt, salt_flux_diag, Q_salt_const, nx, ny)
      !! Full-array overwrite (including ghosts — they get
      !! `Q_salt_const + 0`, the same value the configure-time seed
      !! already gave them, so this is a no-op there). Explicit-shape
      !! dummies + decl-order (integer dims before the arrays that use
      !! them) so NVHPC stdpar compiles a device kernel against static
      !! bounds. Runs on the device-resident `Q_salt` (mapped by
      !! `sf%enter_data`) and `salt_flux_diag` (mapped by `ice%enter_data`).
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: Q_salt(nx, ny)
      real(wp), intent(in) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: Q_salt_const
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         Q_salt(i, j) = Q_salt_const + salt_flux_diag(i, j)
      end do
   end subroutine ice_ocean_brine_flux_impl

   pure subroutine ice_ocean_heat_flux(sf, ice)
      !! Refresh the ocean surface heat-flux field from the ice slot.
      !! **Components off** (default):
      !!   Q_heat(i,j) = Q_heat_const + heat_flux_diag(i,j)
      !! POSITIVE DOWN into the ocean (the apply-tracers convention:
      !! d(hT) = Q_heat*dt/(rho0*cp)). Full overwrite from the const + last
      !! window's rate — a window with no ice exchange resets to background.
      !! **Components on** (`sf%use_components`, PR-12): write
      !! `sf%heat_added(i,j) = heat_flux_diag(i,j)` instead — MOM6's slot
      !! for a net, already-summed heat term (no `_const` term: the
      !! assembler adds it).  Sets `has_heat` either way.
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_sea_ice_t), intent(in) :: ice

      sf%has_heat = .true.
      if (sf%use_components) then
         call ice_ocean_heat_flux_components_impl(sf%heat_added, ice%heat_flux_diag, &
                                                  ice%nx_total, ice%ny_total)
      else
         call ice_ocean_heat_flux_impl(sf%Q_heat, ice%heat_flux_diag, &
                                       sf%Q_heat_const, ice%nx_total, ice%ny_total)
      end if
   end subroutine ice_ocean_heat_flux

   pure subroutine ice_ocean_heat_flux_components_impl(heat_added, heat_flux_diag, nx, ny)
      !! Components-on branch: full-array overwrite of the `heat_added`
      !! COMPONENT (not `Q_heat`) — no `_const` term (the assembler adds
      !! it).  Explicit-shape + decl-order.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: heat_added(nx, ny)
      real(wp), intent(in) :: heat_flux_diag(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         heat_added(i, j) = heat_flux_diag(i, j)
      end do
   end subroutine ice_ocean_heat_flux_components_impl

   pure subroutine ice_ocean_heat_flux_impl(Q_heat, heat_flux_diag, Q_heat_const, nx, ny)
      !! Full-array overwrite (including ghosts — same no-op-there
      !! reasoning as `ice_ocean_brine_flux_impl`). Explicit-shape dummies
      !! + decl-order. Runs on the device-resident `Q_heat` (mapped by
      !! `sf%enter_data`) and `heat_flux_diag` (mapped by `ice%enter_data`).
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: Q_heat(nx, ny)
      real(wp), intent(in) :: heat_flux_diag(nx, ny)
      real(wp), intent(in) :: Q_heat_const
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         Q_heat(i, j) = Q_heat_const + heat_flux_diag(i, j)
      end do
   end subroutine ice_ocean_heat_flux_impl

   pure subroutine ice_ocean_sw_flux(sf, ice)
      !! Refresh the ocean-surface shortwave from the ice slot's
      !! `sw_thru_diag` (PR 31) — the shortwave that penetrated the ice to
      !! the water below (W/m^2, >= 0, POSITIVE DOWN into the ocean).
      !!
      !! `q_sw` is the shortwave SHARE OF the net surface heat flux
      !! (`Q_heat` CONTAINS `q_sw`, it is not `Q_heat + q_sw`), so the
      !! delivery is BRANCH-DEPENDENT and this routine owns exactly ONE
      !! path in each mode — the energy reaches `Q_heat` once, never twice:
      !!
      !! **Components on** (`sf%use_components`, PR-12): full-array
      !! overwrite of the `q_sw` COMPONENT (no `_const` term — the
      !! assembler adds none for shortwave; `q_sw` IS the whole shortwave
      !! summand it folds into `Q_heat`).  `heat_flux_diag` (which the
      !! heat coupler routed to `heat_added`) carries only the non-SW
      !! heat, so there is no overlap.  Sets `has_q_sw`.
      !!
      !! **Components off** (default): `q_sw` is unallocated and no
      !! assembler runs, so the shortwave is ADDED into `Q_heat` here.
      !! MUST run AFTER `ice_ocean_heat_flux` (driver mandated order),
      !! which full-overwrote `Q_heat = Q_heat_const + heat_flux_diag`
      !! this window — so the net is `Q_heat_const + heat_flux_diag +
      !! sw_thru_diag`, recomputed fresh each thermo window (no
      !! accumulation drift).  `has_heat` is already latched by the heat
      !! coupler, so the apply-tracers kernel fires; `has_q_sw` is left
      !! untouched (its backing array is unallocated in this mode, and
      !! nothing consumes it — PR-21's penetration reader is components-on).
      !!
      !! No new `validate_config` guard is needed: the `q_sw` write is
      !! reached ONLY under `use_components`, so an ice-on / components-off
      !! run never touches the unallocated `q_sw` array.
      type(ocean_surface_flux_t), intent(inout) :: sf
      type(ocean_sea_ice_t), intent(in) :: ice

      if (sf%use_components) then
         sf%has_q_sw = .true.
         call ice_ocean_sw_flux_components_impl(sf%q_sw, ice%sw_thru_diag, &
                                                ice%nx_total, ice%ny_total)
      else
         call ice_ocean_sw_flux_add_impl(sf%Q_heat, ice%sw_thru_diag, &
                                         ice%nx_total, ice%ny_total)
      end if
   end subroutine ice_ocean_sw_flux

   pure subroutine ice_ocean_sw_flux_components_impl(q_sw, sw_thru_diag, nx, ny)
      !! Components-on branch: full-array overwrite of the `q_sw` COMPONENT
      !! (including ghosts — they get `0`).  No `_const` term.
      !! Explicit-shape + decl-order.  Runs on the device-resident `q_sw`
      !! (mapped by `sf%enter_data`) and `sw_thru_diag` (mapped by
      !! `ice%enter_data`).
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: q_sw(nx, ny)
      real(wp), intent(in) :: sw_thru_diag(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         q_sw(i, j) = sw_thru_diag(i, j)
      end do
   end subroutine ice_ocean_sw_flux_components_impl

   pure subroutine ice_ocean_sw_flux_add_impl(Q_heat, sw_thru_diag, nx, ny)
      !! Components-off branch: ADD the ice-transmitted shortwave into the
      !! net `Q_heat` (including ghosts — they get `+0`, since the ice
      !! column gates `sw_thru = 0` off-ice and `sw_thru_diag` is
      !! correspondingly 0 there).  This is an ADD, not an overwrite, and
      !! is well-defined precisely because `ice_ocean_heat_flux` ran first
      !! this window and full-overwrote `Q_heat` — so `Q_heat` holds
      !! `Q_heat_const + heat_flux_diag` (SW-free) when this kernel lands
      !! on it, and the net is recomputed fresh every window.  Explicit-
      !! shape + decl-order.  Both arrays device-resident.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: Q_heat(nx, ny)
      real(wp), intent(in) :: sw_thru_diag(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         Q_heat(i, j) = Q_heat(i, j) + sw_thru_diag(i, j)
      end do
   end subroutine ice_ocean_sw_flux_add_impl

   subroutine stress_scratch_ensure(nx, ny)
      !! Lazily (re)allocate + device-map the per-step tau-blend scratch
      !! (F6). Size change tears down and rebuilds. Never called inside a
      !! per-substep loop (this whole coupler runs once per outer step).
      integer, intent(in) :: nx, ny
      if (stress_scratch_ready .and. ss_nx == nx .and. ss_ny == ny) return
      if (stress_scratch_ready) call ice_ocean_stress_cleanup()
      allocate (ci_scratch(nx, ny), source=0.0_wp)
      allocate (mis_scratch(nx, ny), source=0.0_wp)
      allocate (mice_scratch(nx, ny), source=0.0_wp)
      !$acc enter data create(ci_scratch, mis_scratch, mice_scratch)
      ss_nx = nx
      ss_ny = ny
      stress_scratch_ready = .true.
   end subroutine stress_scratch_ensure

   subroutine ice_ocean_stress_cleanup()
      !! Release the persistent tau-blend scratch. Idempotent (safe on an
      !! already-clean workspace — the driver calls it unconditionally at
      !! ocean teardown, next to `ice_evp_cleanup`).
      if (.not. stress_scratch_ready) return
      !$acc exit data delete(ci_scratch, mis_scratch, mice_scratch)
      deallocate (ci_scratch, mis_scratch, mice_scratch)
      ss_nx = 0
      ss_ny = 0
      stress_scratch_ready = .false.
   end subroutine ice_ocean_stress_cleanup

   subroutine ice_ocean_stress_flux(metrics, stress, ice)
      !! Ice->ocean momentum-mediation blend (PR 5), the momentum mirror
      !! of `ice_ocean_brine_flux`: FULL overwrite each outer step from
      !! the pristine wind snapshot + the lagged EVP drag, weighted by
      !! ice concentration at each face.
      !!   a_u(i,j) = 0.5*(ci(i-1,j) + ci(i,j))
      !!   tau_x(i,j) = (1-a_u)*tau_a_x(i,j) + a_u*fxoc(i,j)
      !! ditto y (`a_v(i,j) = 0.5*(ci(i,j-1)+ci(i,j))`).  `ci` is
      !! re-gathered via `ice_cell_concentration_impl` (shared with
      !! `rdb_ice_evp`, so this module never depends on the EVP kernel).
      !!
      !! **Momentum-budget caveat (F5, D7).**  This blend gives the ocean
      !! `(1-a)*tau_a + a*fxoc` per face.  With `&ocean_ice_nml
      !! a_face_stress=.false.` (default, legacy) the ice absorbs the FULL
      !! wind stress `tau_a` and sheds the FULL drag `fxoc` — so at
      !! fractional cover `a in (0,1)` the coupled system sees a spurious
      !! net input `(1-a)*(tau_a - fxoc)` per face.  It is exact only at
      !! `a in {0,1}` and at steady free drift (`fxoc == tau_a`); a
      !! generic transient at fractional `ci` (ITD) leaks momentum.
      !! `&ocean_ice_nml a_face_stress=.true.` (PR 62,
      !! `evp_u_momentum_impl`/`evp_v_momentum_impl` in `rdb_ice_evp`)
      !! CLOSES this identically, for every `a` and in every transient, by
      !! weighting BOTH the wind AND the ice-ocean drag by the same `a_u`
      !! *inside the ice momentum balance* (this blend needs no code
      !! change — it already matches SIS2's `set_ocean_top_stress_Cgrid`
      !! exactly).  Weighting the wind alone — the naive reading of "weight
      !! the wind the ICE feels by `a_face`" — is WORSE than doing nothing:
      !! it converts this leak (zero at steady free drift) into a
      !! PERMANENT one, `-(1-a)*a*tau_a`, nonzero at steady state forever.
      !! Do not implement that half-measure.
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_surface_stress_t), intent(inout) :: stress
      type(ocean_sea_ice_t), intent(inout) :: ice

      call stress_scratch_ensure(ice%nx_total, ice%ny_total)
      call ice_cell_concentration_impl(metrics%wet_T, ice%part_size, ice%m_ice, ice%m_snow, &
                                       mis_scratch, mice_scratch, ci_scratch, ice%ncat, &
                                       ice%nx_total, ice%ny_total)
      call ice_ocean_stress_flux_impl(stress%tau_x, stress%tau_y, ice%tau_a_x, ice%tau_a_y, &
                                      ice%fxoc, ice%fyoc, ci_scratch, ice%nx_total, ice%ny_total)
      ! PR 63: mirror the exact blended value into the restart-carried
      ! tau_ocn_x/y (device copy — ice_tau_mirror_impl is NOT part of the
      ! blend math above, so ice_ocean_stress_flux_impl stays byte-identical
      ! by inspection and the tau_coupling gate is unaffected).  The host
      ! scalar write needs no signature change: this routine is already
      ! non-pure (stress_scratch_ensure) with ice intent(inout).
      call ice_tau_mirror_impl(ice%tau_ocn_x, ice%tau_ocn_y, stress%tau_x, stress%tau_y, &
                               ice%nx_total, ice%ny_total)
      ice%tau_ocn_valid = 1.0_wp
   end subroutine ice_ocean_stress_flux

   pure subroutine ice_ocean_stress_flux_impl(tau_x, tau_y, tau_a_x, tau_a_y, fxoc, fyoc, ci, &
                                              nx, ny)
      !! Face-blend kernel. `a_u`/`a_v` interpolate `ci` onto the u/v
      !! faces (Adcroft-style simple average — no mask needed since `ci`
      !! is already 0 on land). Explicit-shape + decl-order.
      !!
      !! PUBLIC TEST SEAM (PR 62): exported so a test can drive the
      !! real coupler blend on plain arrays (matching a raw `ci_evp_dynamics`
      !! call bit-for-bit) without building a full `ocean_sea_ice_t` and
      !! reverse-engineering a fractional `ci` through the category ITD
      !! gather — mirrors `ice_evp_dynamics`'s own public-seam rationale.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: tau_x(nx + 1, ny)
      real(wp), intent(inout) :: tau_y(nx, ny + 1)
      real(wp), intent(in) :: tau_a_x(nx + 1, ny)
      real(wp), intent(in) :: tau_a_y(nx, ny + 1)
      real(wp), intent(in) :: fxoc(nx + 1, ny)
      real(wp), intent(in) :: fyoc(nx, ny + 1)
      real(wp), intent(in) :: ci(nx, ny)
      integer :: i, j
      real(wp) :: a_u, a_v

      do concurrent(j=1:ny, i=1:nx + 1) local(a_u)
         if (i == 1) then
            a_u = 0.5_wp*ci(1, j)
         else if (i == nx + 1) then
            a_u = 0.5_wp*ci(nx, j)
         else
            a_u = 0.5_wp*(ci(i - 1, j) + ci(i, j))
         end if
         tau_x(i, j) = (1.0_wp - a_u)*tau_a_x(i, j) + a_u*fxoc(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx) local(a_v)
         if (j == 1) then
            a_v = 0.5_wp*ci(i, 1)
         else if (j == ny + 1) then
            a_v = 0.5_wp*ci(i, ny)
         else
            a_v = 0.5_wp*(ci(i, j - 1) + ci(i, j))
         end if
         tau_y(i, j) = (1.0_wp - a_v)*tau_a_y(i, j) + a_v*fyoc(i, j)
      end do
   end subroutine ice_ocean_stress_flux_impl

   pure subroutine ice_tau_mirror_impl(tau_ocn_x, tau_ocn_y, tau_x, tau_y, nx, ny)
      !! PR 63. Device copy of the blend's OUTPUT into the restart-carried
      !! mirror fields — deliberately NOT fused into
      !! `ice_ocean_stress_flux_impl` (that kernel stays untouched so the
      !! `tau_coupling` bitwise gate is unaffected by this PR).
      !! Explicit-shape + decl-order (integer dims before the arrays that
      !! use them).
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: tau_ocn_x(nx + 1, ny)
      real(wp), intent(out) :: tau_ocn_y(nx, ny + 1)
      real(wp), intent(in) :: tau_x(nx + 1, ny)
      real(wp), intent(in) :: tau_y(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         tau_ocn_x(i, j) = tau_x(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         tau_ocn_y(i, j) = tau_y(i, j)
      end do
   end subroutine ice_tau_mirror_impl

   pure subroutine ice_ocean_stress_resume_apply(stress, ice)
      !! Configure-time resume apply (PR 63).  `tau_x`/`tau_y` were just
      !! re-seeded from the wind-stress config by `configure_ocean_forcing`;
      !! if the checkpoint carried a blend (`tau_ocn_valid > 0.5`),
      !! overwrite them with it — the EXACT stress the uninterrupted run
      !! would have handed the ocean at this step boundary.  Otherwise
      !! leave the pristine wind (D8: a fresh run's first outer step
      !! drives the ocean with pure wind — and, for the first time, this
      !! is now literally true: no `(1-a)*tau_a` residual from an
      !! uninitialised `fxoc`).
      !!
      !! REPLACES the old configure-time resume-fold routine (deleted).
      !! That routine RECONSTRUCTED the blend from the checkpointed `ci`
      !! — which is the POST-thermo/post-transport concentration, not the
      !! PRE-thermo `ci`
      !! the uninterrupted run actually blended against (F4) — so it was
      !! wrong whenever a checkpoint step's thermo/transport changed `ci`
      !! after the blend.  This routine COPIES the blend's own recorded
      !! output instead, which is exact by construction and stays exact
      !! under any future change to the blend formula (SIS2 `Ice%flux_u`/
      !! `flux_v`, `ice_type.F90:240-243` — carried, not reconstructed).
      !!
      !! Host-only whole-array assignment; runs BEFORE `enter_data`
      !! (mirrors the deleted routine's host-only contract, but for a
      !! trivial reason now — a plain copy needs no device kernel at all,
      !! so the GPU-hazard this used to carry — do-concurrent kernels over
      !! unmapped host arrays at configure time, see rdb_ice_evp's D8 note
      !! and CLAUDE.md — is gone by construction, not by discipline).
      type(ocean_surface_stress_t), intent(inout) :: stress
      type(ocean_sea_ice_t), intent(in) :: ice
      if (ice%tau_ocn_valid <= 0.5_wp) return
      stress%tau_x = ice%tau_ocn_x
      stress%tau_y = ice%tau_ocn_y
   end subroutine ice_ocean_stress_resume_apply

end module rdb_ice_ocean_coupler
