!! Configure-time numerical-stability audit for the ocean dyn-core.
module rdb_ocean_stability_audit
   !! Motivating failure (`tmp_local_artifacts/global_run/FINDINGS.md`,
   !! 2026-09-11): a global tripolar aquaplanet NaN'd at outer step 7. The
   !! ONLY diagnostic on offer was a bare non-finite-face count — "producer
   !! 0/0 upstream, investigate". A human needed several runs + a
   !! bisection to find the actual cause: `nu_h = 2.0e4` with `dt = 900 s`
   !! at the ~3 km polar cells gives a viscous-diffusion number of 1.65
   !! against an explicit-Laplacian bound of 0.125 — 13x over — and
   !! `bound_kh` (the per-cell runtime clamp that would have protected
   !! against exactly this) was never enabled.
   !!
   !! Nothing at configure time said any of that. This module is the fix:
   !! a small set of checks run once, AFTER the real per-cell metric
   !! arrays exist (`ocean_metrics_t`, filled by `configure_ocean_metrics`
   !! + land-masked by `configure_ocean_land_mask`) — never off the
   !! nominal `&grid_nml dx`/`dy`, which are DEGREES on spherical/tripolar
   !! grids and are in any case the NOMINAL spacing, not the smallest
   !! actual cell (a spherical/tripolar grid's smallest cell can be an
   !! order of magnitude below nominal near a pole).
   !!
   !! Each check reports the computed number, the limit, the knob(s)
   !! responsible, and a concrete fix — never a bare "X exceeded".
   !!
   !! Severity: a VIOLATED HARD STABILITY BOUND (viscous CFL, the tracer
   !! diffusive number) is an ERROR — returned via the P0 `ierr` status
   !! (`OCEAN_STATUS_ERR_SETUP`), never `error stop` (this module always
   !! has `ierr` to report through; `configure_ocean_metrics` et al. do
   !! the same). A MARGINAL/QUALITY issue (the Munk-layer resolution
   !! criterion, `ah_max` silently clamping `nu_h`) is a WARNING — logged,
   !! run proceeds. The viscous-CFL check is itself downgraded from ERROR
   !! to an informational WARNING when the run already carries automatic
   !! runtime protection (`bound_kh`, or the `stress_tensor` operator's
   !! own always-on per-cell CFL limiter) — the raw `nu_h*dt/dx^2` number
   !! is no longer what the kernel actually uses in that case, so a hard
   !! configure-time failure would be a false positive (see
   !! `docs/CLOSURE_MATRIX.md` / `rdb_ocean_horizontal_viscosity.F90`
   !! module header for `bound_kh` / `stress_tensor` semantics).
   use rdb_constants, only: wp, PI
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, parse_grid_config, &
                                parse_coriolis_scheme, GRID_CONFIG_CARTESIAN, &
                                CORIOLIS_SCHEME_PLANETARY
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_error_ring, only: fail
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   implicit none
   private

   public :: ocean_stability_audit
   public :: ocean_viscous_cfl_number
   public :: ocean_viscous_cfl_max_nu_h
   public :: ocean_viscous_cfl_max_dt
   public :: ocean_diffusive_number
   public :: ocean_munk_delta_m
   public :: ocean_munk_required_nu_h
   public :: ocean_viscous_cfl_limit
   public :: ocean_diffusive_number_limit
   public :: ocean_munk_min_cells

   real(wp), parameter :: VISCOUS_CFL_LIMIT = 0.125_wp
      !! Single-axis viscous-diffusion stability bound `nu_h*dt/dx_min^2`.
      !! Taken directly from the constant this codebase ALREADY uses at
      !! runtime for exactly this quantity: the `bound_kh` per-face clamp
      !! (`rdb_ocean_horizontal_viscosity.F90`) limits the harmonic
      !! viscosity to `bound_coef*0.125/(dt*(idx^2+idy^2))`, documented
      !! there as "~1/4 of the forward-Euler stability limit" of 0.5 (the
      !! same two-axis sum-form bound `rdb_ocean_hdiff_tracer.F90`'s
      !! `kappa_h` check already uses, see `DIFFUSIVE_NUMBER_LIMIT`
      !! below). This audit collapses that two-axis form to the single
      !! worst axis (`nu_h*dt/dx_min^2` rather than
      !! `nu_h*dt*(1/dx_min^2+1/dy_min^2)`) so the reported number matches
      !! the plain "viscous CFL" a user computes by hand, at `bound_coef=1`
      !! parity with the runtime clamp's own margin — i.e. this check trips
      !! at exactly the `nu_h` that would need `bound_kh`'s protection.
   real(wp), parameter :: DIFFUSIVE_NUMBER_LIMIT = 0.5_wp
      !! Two-axis explicit forward-Euler Laplacian stability bound
      !! `kappa_h*dt_therm*(1/dx^2+1/dy^2) <= 0.5` — unchanged from the
      !! existing `rdb_config.F90` check this module absorbs (only the
      !! length scale changes: the real per-cell metric minimum, not
      !! nominal `dx`/`dy`).
   real(wp), parameter :: MUNK_MIN_CELLS = 2.0_wp
      !! Minimum number of grid cells the Munk sidewall boundary layer
      !! `delta_M = (nu_h/beta)^(1/3)` must span; below this the wall
      !! carries grid-scale (2-delta) noise instead of a resolved
      !! boundary-layer profile (see
      !! `validation_examples/ocean/acc_channel/acc_channel.nml`, where
      !! this exact criterion is documented and was hand-derived).

contains

   pure function ocean_viscous_cfl_limit() result(lim)
      !! Accessor for `VISCOUS_CFL_LIMIT` — tests reference this instead
      !! of duplicating the literal.
      real(wp) :: lim
      lim = VISCOUS_CFL_LIMIT
   end function ocean_viscous_cfl_limit

   pure function ocean_diffusive_number_limit() result(lim)
      !! Accessor for `DIFFUSIVE_NUMBER_LIMIT`.
      real(wp) :: lim
      lim = DIFFUSIVE_NUMBER_LIMIT
   end function ocean_diffusive_number_limit

   pure function ocean_munk_min_cells() result(n)
      !! Accessor for `MUNK_MIN_CELLS`.
      real(wp) :: n
      n = MUNK_MIN_CELLS
   end function ocean_munk_min_cells

   pure function ocean_viscous_cfl_number(nu_h, dt, dx_min) result(cfl)
      !! `nu_h*dt/dx_min^2` — the single-axis viscous-diffusion stability
      !! number checked against `VISCOUS_CFL_LIMIT`. `dx_min` MUST be the
      !! smallest actual cell edge in the domain (e.g.
      !! `metrics_dx_min`/`ocean_stability_min_cell`), never a nominal
      !! `&grid_nml dx`/`dy` (degrees on non-Cartesian grids). `dx_min<=0`
      !! (degenerate/unset grid) returns 0 (no constraint expressible).
      real(wp), intent(in) :: nu_h, dt, dx_min
      real(wp) :: cfl
      if (dx_min > 0.0_wp) then
         cfl = nu_h*dt/dx_min**2
      else
         cfl = 0.0_wp
      end if
   end function ocean_viscous_cfl_number

   pure function ocean_viscous_cfl_max_nu_h(dt, dx_min, limit) result(nu_h_max)
      !! Largest `nu_h` (m^2/s) that keeps `ocean_viscous_cfl_number` at
      !! or below `limit`, at fixed `dt`/`dx_min` — the "reduce nu_h
      !! below ..." half of the audit's suggested fix.
      real(wp), intent(in) :: dt, dx_min, limit
      real(wp) :: nu_h_max
      if (dt > 0.0_wp) then
         nu_h_max = limit*dx_min**2/dt
      else
         nu_h_max = 0.0_wp
      end if
   end function ocean_viscous_cfl_max_nu_h

   pure function ocean_viscous_cfl_max_dt(nu_h, dx_min, limit) result(dt_max)
      !! Largest `dt` (s) that keeps `ocean_viscous_cfl_number` at or
      !! below `limit`, at fixed `nu_h`/`dx_min` — the "reduce dt below
      !! ..." half of the audit's suggested fix.
      real(wp), intent(in) :: nu_h, dx_min, limit
      real(wp) :: dt_max
      if (nu_h > 0.0_wp) then
         dt_max = limit*dx_min**2/nu_h
      else
         dt_max = huge(1.0_wp)
      end if
   end function ocean_viscous_cfl_max_dt

   pure function ocean_diffusive_number(kappa_h, dt_therm, dx_min) result(dnum)
      !! Two-axis explicit forward-Euler diffusive number
      !! `kappa_h*dt_therm*(1/dx_min^2+1/dy_min^2)`, conservatively
      !! evaluated at the SAME worst-case `dx_min` on both axes (matches
      !! the "use the minimum cell" instruction; exact on an isotropic
      !! worst cell, strictly more conservative than using the true
      !! per-axis pair). `dx_min<=0` returns 0.
      real(wp), intent(in) :: kappa_h, dt_therm, dx_min
      real(wp) :: dnum
      if (dx_min > 0.0_wp) then
         dnum = kappa_h*dt_therm*(2.0_wp/dx_min**2)
      else
         dnum = 0.0_wp
      end if
   end function ocean_diffusive_number

   pure function ocean_munk_delta_m(nu_h, beta) result(delta_m)
      !! Munk boundary-layer width `delta_M = (nu_h/beta)^(1/3)`.
      !! `beta<=0` (f-plane — no meridional PV gradient, no Munk
      !! boundary layer) or `nu_h<=0` returns `huge(1.0_wp)` (no
      !! constraint — never trips the >= 2-cell criterion).
      real(wp), intent(in) :: nu_h, beta
      real(wp) :: delta_m
      if (beta > 0.0_wp .and. nu_h > 0.0_wp) then
         delta_m = (nu_h/beta)**(1.0_wp/3.0_wp)
      else
         delta_m = huge(1.0_wp)
      end if
   end function ocean_munk_delta_m

   pure function ocean_munk_required_nu_h(beta, dx, n_cells) result(nu_h_req)
      !! `nu_h` (m^2/s) needed for `delta_M` to span exactly `n_cells` of
      !! width `dx` — the audit's suggested fix for a Munk-layer warning.
      real(wp), intent(in) :: beta, dx, n_cells
      real(wp) :: nu_h_req
      nu_h_req = beta*(n_cells*dx)**3
   end function ocean_munk_required_nu_h

   pure subroutine ocean_stability_min_cell(metrics, grid, dx_min, i_at, j_at, is_x)
      !! Smallest actual cell edge over the PHYSICAL domain (excludes
      !! ghosts), taken over BOTH `dxT` and `dyT`, with its (i,j)
      !! location and which axis (`is_x`) it came from — for actionable
      !! messages ("near j=110"). Host-side, configure time; the grid
      !! sizes here are at most a few 10^5 cells (a global tripolar
      !! config), trivial to scan once.
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: dx_min
      integer, intent(out) :: i_at, j_at
      logical, intent(out) :: is_x
      integer :: i, j, ng, i0, i1, j0, j1

      ng = grid%nghost
      i0 = ng + 1
      i1 = ng + grid%nx_phys
      j0 = ng + 1
      j1 = ng + grid%ny_phys
      dx_min = huge(1.0_wp)
      i_at = i0
      j_at = j0
      is_x = .true.
      do j = j0, j1
         do i = i0, i1
            if (metrics%dxT(i, j) < dx_min) then
               dx_min = metrics%dxT(i, j)
               i_at = i
               j_at = j
               is_x = .true.
            end if
            if (metrics%dyT(i, j) < dx_min) then
               dx_min = metrics%dyT(i, j)
               i_at = i
               j_at = j
               is_x = .false.
            end if
         end do
      end do
   end subroutine ocean_stability_min_cell

   pure subroutine ocean_munk_worst_case(cfg, metrics, grid, nu_h, ratio_min, beta_at, dx_at, j_at)
      !! Worst-case (smallest) `delta_M/dx` ratio over the physical
      !! domain, honouring a latitude-varying `beta` under
      !! `coriolis_scheme='planetary'` on a non-Cartesian grid (`beta =
      !! 2*omega*cos(lat)/R`, maximal — hence `delta_M` MINIMAL, the
      !! worst case — at the most equatorward row) and a constant `beta`
      !! (`&ocean_topo_nml coriolis_beta`) everywhere else. `ratio_min` is
      !! `huge(1.0_wp)` (no constraint) when `beta<=0` everywhere (an
      !! f-plane run has no Munk boundary layer to resolve).
      type(config_t), intent(in) :: cfg
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: nu_h
      real(wp), intent(out) :: ratio_min, beta_at, dx_at
      integer, intent(out) :: j_at

      real(wp), parameter :: DEG2RAD = PI/180.0_wp
      logical :: planetary
      integer :: i, j, ng, i0, i1, j0, j1
      real(wp) :: beta_ij, dx_ij, delta_ij, ratio_ij

      planetary = (parse_grid_config(cfg%ocean%grid%grid_config) /= GRID_CONFIG_CARTESIAN &
                   .and. parse_coriolis_scheme(cfg%ocean%grid%coriolis_scheme) == CORIOLIS_SCHEME_PLANETARY)

      ratio_min = huge(1.0_wp)
      beta_at = 0.0_wp
      dx_at = 0.0_wp
      j_at = 0

      if (nu_h <= 0.0_wp) return

      ng = grid%nghost
      i0 = ng + 1
      i1 = ng + grid%nx_phys
      j0 = ng + 1
      j1 = ng + grid%ny_phys

      if (.not. planetary) then
         ! Uniform beta (beta_plane, or the f-plane beta<=0 no-op).
         beta_ij = cfg%ocean%topo%coriolis_beta
         if (beta_ij <= 0.0_wp) return
         do j = j0, j1
            do i = i0, i1
               dx_ij = min(metrics%dxT(i, j), metrics%dyT(i, j))
               delta_ij = ocean_munk_delta_m(nu_h, beta_ij)
               ratio_ij = delta_ij/dx_ij
               if (ratio_ij < ratio_min) then
                  ratio_min = ratio_ij
                  beta_at = beta_ij
                  dx_at = dx_ij
                  j_at = j
               end if
            end do
         end do
      else
         do j = j0, j1
            do i = i0, i1
               beta_ij = 2.0_wp*cfg%ocean%grid%omega*cos(metrics%geolatT(i, j)*DEG2RAD)/ &
                         cfg%ocean%grid%rad_earth
               if (beta_ij <= 0.0_wp) cycle
               dx_ij = min(metrics%dxT(i, j), metrics%dyT(i, j))
               delta_ij = ocean_munk_delta_m(nu_h, beta_ij)
               ratio_ij = delta_ij/dx_ij
               if (ratio_ij < ratio_min) then
                  ratio_min = ratio_ij
                  beta_at = beta_ij
                  dx_at = dx_ij
                  j_at = j
               end if
            end do
         end do
      end if
   end subroutine ocean_munk_worst_case

   subroutine ocean_stability_audit(cfg, metrics, grid, rank, ierr)
      !! Run all configure-time stability checks. Must run AFTER
      !! `configure_ocean_metrics` + `configure_ocean_land_mask` (needs
      !! the real filled `ocean_metrics_t`), before `ocean_state_enter_data`.
      !! `ierr` present -> `OCEAN_STATUS_ERR_SETUP` on any hard-bound
      !! violation (never `error stop`); warnings always just log,
      !! whatever `ierr` does. Rank-0-only logging (mirrors every other
      !! `configure_ocean_*` info/warning line); the ERROR path itself
      !! always fires (every rank must agree the config is broken).
      type(config_t), intent(in) :: cfg
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: rank
      integer, intent(out), optional :: ierr

      real(wp) :: dx_min, dt, dt_therm, nu_h, kappa_h, ah_max
      integer :: i_at, j_at
      logical :: is_x
      logical :: auto_protected

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      dt = cfg%dt_fixed
      nu_h = cfg%ocean%hvisc%nu_h
      kappa_h = cfg%ocean%hdiff%kappa_h
      ah_max = cfg%ocean%hvisc%ah_max

      ! dt_fixed <= 0 => adaptive-CFL dt, not known at configure time
      ! (matches the pre-existing kappa_h check's own gate).
      if (dt <= 0.0_wp) return

      call ocean_stability_min_cell(metrics, grid, dx_min, i_at, j_at, is_x)
      if (dx_min <= 0.0_wp .or. dx_min == huge(1.0_wp)) return

      ! ---- Check 1: viscous CFL (nu_h*dt/dx_min^2) ----
      auto_protected = cfg%ocean%hvisc%bound_kh .or. cfg%ocean%hvisc%stress_tensor
      if (nu_h > 0.0_wp) then
         block
            real(wp) :: cfl, nu_h_max, dt_max
            cfl = ocean_viscous_cfl_number(nu_h, dt, dx_min)
            if (cfl > VISCOUS_CFL_LIMIT) then
               nu_h_max = ocean_viscous_cfl_max_nu_h(dt, dx_min, VISCOUS_CFL_LIMIT)
               dt_max = ocean_viscous_cfl_max_dt(nu_h, dx_min, VISCOUS_CFL_LIMIT)
               if (auto_protected) then
                  if (rank == 0) then
                     call logger%warning("&ocean_hvisc_nml nu_h = "//to_string(nu_h)// &
                                         " with dt = "//to_string(dt)//"s gives viscous CFL "// &
                                         to_string(cfl)//" at the smallest cell (dx = "// &
                                         to_string(dx_min/1000.0_wp)//" km, near i="//to_string(i_at)// &
                                         " j="//to_string(j_at)//"); limit is "// &
                                         to_string(VISCOUS_CFL_LIMIT)//" — but bound_kh/stress_tensor "// &
                                         "is enabled, so the runtime per-cell clamp will keep the "// &
                                         "EFFECTIVE viscosity within bound. Informational only.")
                  end if
               else
                  call fail("&ocean_hvisc_nml nu_h = "//to_string(nu_h)// &
                            " with dt = "//to_string(dt)//"s gives viscous CFL "// &
                            to_string(cfl)//" at the smallest cell (dx = "// &
                            to_string(dx_min/1000.0_wp)//" km, near i="//to_string(i_at)// &
                            " j="//to_string(j_at)//"); limit is "//to_string(VISCOUS_CFL_LIMIT)// &
                            ". Set &ocean_hvisc_nml bound_kh = .true. for a per-cell clamp, "// &
                            "or reduce nu_h below "//to_string(nu_h_max)// &
                            ", or reduce dt below "//to_string(dt_max)//"s.", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
            end if
         end block
      end if

      ! ---- Check 2: tracer diffusive number (kappa_h) ----
      ! Absorbed + fixed from the old rdb_config.F90 check (Cartesian-only,
      ! nominal dx/dy): now uses the REAL minimum cell, on ANY grid type.
      if (kappa_h > 0.0_wp) then
         block
            real(wp) :: dnum
            dt_therm = dt*real(cfg%ocean%vmix%dt_therm_ratio, wp)
            dnum = ocean_diffusive_number(kappa_h, dt_therm, dx_min)
            if (dnum > DIFFUSIVE_NUMBER_LIMIT) then
               call fail("&ocean_hdiff_nml kappa_h = "//to_string(kappa_h)// &
                         " violates the explicit forward-Euler diffusive stability bound "// &
                         "kappa_h*dt_therm*(1/dx_min^2+1/dy_min^2) <= "// &
                         to_string(DIFFUSIVE_NUMBER_LIMIT)//" (computed "//to_string(dnum)// &
                         " using dt_therm = "//to_string(dt_therm)//"s, dx_min = "// &
                         to_string(dx_min)//"m near i="//to_string(i_at)//" j="//to_string(j_at)// &
                         " — the ACTUAL smallest cell, not nominal dx/dy). Reduce kappa_h "// &
                         "below "//to_string(DIFFUSIVE_NUMBER_LIMIT*dx_min**2/(2.0_wp*dt_therm))// &
                         ", or reduce dt/dt_therm_ratio.", ierr, OCEAN_STATUS_ERR_SETUP)
               return
            end if
         end block
      end if

      ! ---- Check 3: Munk-layer resolution (WARNING) ----
      if (nu_h > 0.0_wp) then
         block
            real(wp) :: ratio_min, beta_at, dx_at, nu_h_req
            integer :: j_munk
            call ocean_munk_worst_case(cfg, metrics, grid, nu_h, ratio_min, beta_at, dx_at, j_munk)
            if (ratio_min < MUNK_MIN_CELLS .and. ratio_min < huge(1.0_wp)) then
               nu_h_req = ocean_munk_required_nu_h(beta_at, dx_at, MUNK_MIN_CELLS)
               if (rank == 0) then
                  call logger%warning("Munk sidewall boundary layer under-resolved: "// &
                                      "delta_M = (nu_h/beta)^(1/3) = "//to_string(ocean_munk_delta_m(nu_h, beta_at))// &
                                      "m spans only "//to_string(ratio_min)//" cells (dx = "// &
                                      to_string(dx_at)//"m, beta = "//to_string(beta_at)//" 1/(m*s) near j="// &
                                      to_string(j_munk)//"); need >= "//to_string(MUNK_MIN_CELLS)// &
                                      " cells or the wall carries grid-scale (2-delta) noise. Raise "// &
                                      "&ocean_hvisc_nml nu_h to at least "//to_string(nu_h_req)// &
                                      " m2/s (and raise ah_max to match — see the next check).")
               end if
            end if
         end block
      end if

      ! ---- Check 4: ah_max clamping nu_h inert (WARNING) ----
      if (nu_h > 0.0_wp .and. ah_max > 0.0_wp .and. ah_max < nu_h) then
         if (rank == 0) then
            call logger%warning("&ocean_hvisc_nml ah_max = "//to_string(ah_max)// &
                                " is BELOW nu_h = "//to_string(nu_h)//" — ah_max is a hard CAP on the "// &
                                "per-face viscosity, so the configured nu_h is silently clamped down to "// &
                                "ah_max and never actually applied. Set ah_max >= nu_h (e.g. ah_max = "// &
                                to_string(nu_h)//" or higher).")
         end if
      end if

      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_stability_audit

end module rdb_ocean_stability_audit
