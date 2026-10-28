!! File-based sea-surface-salinity restoring as a salt flux, with a
!! per-cell clamp and a global net-zero adjustment.
module rdb_ocean_sss_restore
   !! Surface salinity restoring toward a time-varying 2-D target field
   !! (`surface_flux%sss_target`, blended from `&ocean_restore_nml sss_file`
   !! by `ocean_data_forcing_apply`), in the salt-flux form used by
   !! climatological ocean models (Griffies et al. 2009, "Coordinated
   !! Ocean-ice Reference Experiments (COREs)", Ocean Modelling 26, sec.
   !! 4.1; the net-zero adjustment is the CORE-II salinity-restoring
   !! practice).  Per wet cell, with `S` the surface-layer salinity:
   !!
   !!   1. `dS = clamp(S_tgt - S, +-max_delta)`
   !!   2. `F  = rho0 * w * dS * open`        (`w` = piston velocity, m/s)
   !!   3. global net removal (optional), over PHYSICAL cells, weighted by
   !!      the cell area `A`:
   !!        scaling:  `P = sum A max(F,0)`, `N = sum A min(F,0)`;
   !!                  if `P > 0 > N` the larger-magnitude side is scaled
   !!                  by `min(P,-N)/max(P,-N)` so `sum A F = 0` and no
   !!                  cell changes sign;
   !!        mean:     `F <- F - (sum A F / sum A open) * open`.
   !!   4. deposit `d(hS_top) = (dt/rho0) F` into the first live layer
   !!      (`k_top`) and mirror it into `salt_budget_surface`, exactly as
   !!      `apply_surface_src_2d_impl` does for `Q_salt` (so the closed
   !!      salt budget sees it as a source).
   !!
   !! **Units.**  `F` is in the model's `Q_salt` convention, the one
   !! `dS/dt = Q_salt/(rho0 h)` is written in: salinity in g/kg times
   !! kg/m^2/s.  In SI salt mass that is `1e-3 * F` [kg salt m^-2 s^-1];
   !! the 1e3 is absorbed because the tracer is carried in g/kg.  The
   !! resulting tendency is the piston law `dS/dt = w dS / h`, identical to
   !! the scalar-target path.
   !!
   !! **Why `F` is never stored.**  `F` is a pure function of the live
   !! surface salinity, the target, and scalars, so the kernels recompute
   !! it in each of the (at most) three passes - the two global sums and
   !! the deposit.  That avoids a persistent workspace (whose lifecycle and
   !! `intent(in)` access through the dynamical core's `sf` dummy would
   !! both be awkward) at the price of one extra trivially cheap
   !! evaluation per thermo step.
   !!
   !! **Decomposition independence.**  `P`, `N` and the open area are
   !! accumulated in EFP fixed point (Hallberg & Adcroft 2014) and
   !! combined across ranks in ONE `halo_allreduce_efp_list`, so the
   !! scaling factor is bit-identical for any rank count.  The target's
   !! seam ghosts are exchanged by the caller, so the full-array deposit
   !! matches the neighbour tile.
   !!
   !! A cell whose target is non-finite (a masked file value) or whose
   !! open fraction is zero takes no part: `F = 0` there and it is
   !! excluded from every sum.  `mem:separate`: every array is owned and
   !! mapped by `ms` / `sf` / `metrics`.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_console_stats, only: efp_decompose_impl
   use rdb_efp, only: efp_t, EFP_DIGITS, efp_carry, efp_to_real
   use rdb_halo, only: halo_allreduce_efp_list, halo_allreduce_max
   use rdb_teos10_convert, only: sr_from_sp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_sss_restore_apply
   public :: sss_cell_flux
   public :: sss_net_factors
   public :: sss_flux_sums
   public :: sss_flux_maxabs
   public :: sss_deposit_impl

   integer, parameter, public :: SSS_NET_NONE = 0
   integer, parameter, public :: SSS_NET_MEAN = 1
   integer, parameter, public :: SSS_NET_SCALE = 2

   integer, parameter :: SUM_POS = 1
   integer, parameter :: SUM_NEG = 2
   integer, parameter :: SUM_AREA = 3

   integer, save :: n_calls = 0
      !! Calls so far, for the sparse progress log line.

contains

   pure subroutine sss_cell_flux(hTr, h, tgt_in, w_open, coef, max_delta, h_min, &
                                 convert_sp, f, pw)
      !! Restoring salt flux of one cell (steps 1-2 above).  `coef` is
      !! `rho0 * w`.  `pw` is the participation weight: the open fraction
      !! for a participating cell, `0` for a cell that is closed or has a
      !! non-finite target (those return `f = 0`).
      !$acc routine seq
      real(wp), intent(in) :: hTr
         !! Surface-layer `S * h` [g/kg m].
      real(wp), intent(in) :: h
         !! Surface-layer thickness [m].
      real(wp), intent(in) :: tgt_in
         !! Target salinity as read (SP when `convert_sp`).
      real(wp), intent(in) :: w_open
         !! Open fraction: `wet_mask * (1 - cover_frac)`.
      real(wp), intent(in) :: coef
         !! `rho0 * w`.
      real(wp), intent(in) :: max_delta
         !! Clamp on `|S_tgt - S|` [g/kg].
      real(wp), intent(in) :: h_min
         !! Floor on `h` for the surface concentration [m].
      logical, intent(in) :: convert_sp
      real(wp), intent(out) :: f
      real(wp), intent(out) :: pw
      real(wp) :: tgt, ds

      f = 0.0_wp
      pw = 0.0_wp
      if (.not. (w_open > 0.0_wp)) return
      if (.not. ieee_is_finite(tgt_in)) return
      tgt = tgt_in
      if (convert_sp) tgt = sr_from_sp(tgt_in)
      ! vanished-ok: surface concentration against a caller-supplied
      ! restoring floor, as the scalar restoring path does.
      ds = tgt - hTr/max(h, h_min)
      ds = sign(min(abs(ds), max_delta), ds)
      f = coef*ds*w_open
      pw = w_open
   end subroutine sss_cell_flux

   pure subroutine sss_net_factors(p, n, fpos, fneg)
      !! Scaling factors that make `P + N = 0` without moving the zero
      !! contour: the larger-magnitude side is scaled by
      !! `min(P,-N)/max(P,-N)`.  Unchanged (1, 1) unless `P > 0 > N`.
      real(wp), intent(in) :: p
         !! `sum A max(F,0)` (>= 0).
      real(wp), intent(in) :: n
         !! `sum A min(F,0)` (<= 0).
      real(wp), intent(out) :: fpos
         !! Multiplier for positive fluxes.
      real(wp), intent(out) :: fneg
         !! Multiplier for negative fluxes.
      fpos = 1.0_wp
      fneg = 1.0_wp
      if (p > 0.0_wp .and. n < 0.0_wp) then
         if (p > -n) then
            fpos = -n/p
         else if (p < -n) then
            fneg = -p/n
         end if
      end if
   end subroutine sss_net_factors

   function sss_flux_sums(mode, nx, ny, nz, ng, hTr, h_layer, tgt, wet, cover, use_cover, &
                          areaT, k_top, coef, max_delta, h_min, convert_sp) result(tot)
      !! Rank-LOCAL EFP area integral over physical cells of one
      !! integrand: `mode = SUM_POS` -> `A max(F,0)`, `SUM_NEG` ->
      !! `A min(F,0)`, `SUM_AREA` -> `A * participation`.  One device
      !! reduction; combine across ranks with `halo_allreduce_efp_list`.
      integer, intent(in) :: mode, nx, ny, nz, ng
      real(wp), intent(in) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: tgt(nx, ny)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: cover(nx, ny)
      logical, intent(in) :: use_cover
      real(wp), intent(in) :: areaT(nx, ny)
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(in) :: coef, max_delta, h_min
      logical, intent(in) :: convert_sp
      type(efp_t) :: tot
      integer :: i, j, kt
      integer(int64) :: e1, e2, e3, e4, e5, e6, epoison
      integer(int64) :: d1, d2, d3, d4, d5, d6, dpoison
      integer(int64) :: bins(EFP_DIGITS)
      real(wp) :: f, pw, ow, val

      e1 = 0_int64
      e2 = 0_int64
      e3 = 0_int64
      e4 = 0_int64
      e5 = 0_int64
      e6 = 0_int64
      epoison = 0_int64
      !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6,epoison) &
      !$acc&         private(f, pw, ow, val, kt, d1, d2, d3, d4, d5, d6, dpoison) &
      !$acc&         present(hTr, h_layer, tgt, wet, cover, areaT, k_top)
      do j = ng + 1, ny - ng
         do i = ng + 1, nx - ng
            kt = k_top(i, j)
            ow = wet(i, j)
            if (use_cover) ow = ow*(1.0_wp - cover(i, j))
            call sss_cell_flux(hTr(i, j, kt), h_layer(i, j, kt), tgt(i, j), ow, &
                               coef, max_delta, h_min, convert_sp, f, pw)
            select case (mode)
            case (SUM_POS)
               val = areaT(i, j)*max(f, 0.0_wp)
            case (SUM_NEG)
               val = areaT(i, j)*min(f, 0.0_wp)
            case default
               val = areaT(i, j)*pw
            end select
            call efp_decompose_impl(real(val, real64), d1, d2, d3, d4, d5, d6, dpoison)
            e1 = e1 + d1
            e2 = e2 + d2
            e3 = e3 + d3
            e4 = e4 + d4
            e5 = e5 + d5
            e6 = e6 + d6
            epoison = epoison + dpoison
         end do
      end do
      bins = [e1, e2, e3, e4, e5, e6]
      call efp_carry(bins)
      tot%v = bins
      tot%poison = epoison
   end function sss_flux_sums

   function sss_flux_maxabs(nx, ny, nz, ng, hTr, h_layer, tgt, wet, cover, use_cover, &
                            k_top, coef, max_delta, h_min, convert_sp) result(fmax)
      !! Rank-local `max |F|` over physical cells (before the net
      !! adjustment); a diagnostic for the sparse progress line.
      integer, intent(in) :: nx, ny, nz, ng
      real(wp), intent(in) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: tgt(nx, ny)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: cover(nx, ny)
      logical, intent(in) :: use_cover
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(in) :: coef, max_delta, h_min
      logical, intent(in) :: convert_sp
      real(wp) :: fmax
      integer :: i, j, kt
      real(wp) :: f, pw, ow
      fmax = 0.0_wp
      !$acc parallel loop collapse(2) reduction(max:fmax) private(f, pw, ow, kt) &
      !$acc&         present(hTr, h_layer, tgt, wet, cover, k_top)
      do j = ng + 1, ny - ng
         do i = ng + 1, nx - ng
            kt = k_top(i, j)
            ow = wet(i, j)
            if (use_cover) ow = ow*(1.0_wp - cover(i, j))
            call sss_cell_flux(hTr(i, j, kt), h_layer(i, j, kt), tgt(i, j), ow, &
                               coef, max_delta, h_min, convert_sp, f, pw)
            fmax = max(fmax, abs(f))
         end do
      end do
   end function sss_flux_maxabs

   pure subroutine sss_deposit_impl(nx, ny, nz, hTr, budget, h_layer, tgt, wet, cover, &
                                    use_cover, k_top, coef, max_delta, h_min, convert_sp, &
                                    net_mode, fpos, fneg, shift, inv_scale)
      !! Recompute `F`, apply the net adjustment, and deposit
      !! `inv_scale * F` (`inv_scale = dt/rho0`) into the first live layer
      !! and the salt budget contributor.  Over the full array (ghosts
      !! included - the target's ghosts were exchanged).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: tgt(nx, ny)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: cover(nx, ny)
      logical, intent(in) :: use_cover
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(in) :: coef, max_delta, h_min
      logical, intent(in) :: convert_sp
      integer, intent(in) :: net_mode
      real(wp), intent(in) :: fpos, fneg
         !! Scaling factors (`SSS_NET_SCALE`).
      real(wp), intent(in) :: shift
         !! Area-mean flux to subtract (`SSS_NET_MEAN`).
      real(wp), intent(in) :: inv_scale
      integer :: i, j, kt
      real(wp) :: f, pw, ow, inc
      do concurrent(j=1:ny, i=1:nx) local(f, pw, ow, inc, kt)
         kt = k_top(i, j)
         ow = wet(i, j)
         if (use_cover) ow = ow*(1.0_wp - cover(i, j))
         call sss_cell_flux(hTr(i, j, kt), h_layer(i, j, kt), tgt(i, j), ow, &
                            coef, max_delta, h_min, convert_sp, f, pw)
         if (net_mode == SSS_NET_SCALE) then
            if (f > 0.0_wp) then
               f = f*fpos
            else if (f < 0.0_wp) then
               f = f*fneg
            end if
         else if (net_mode == SSS_NET_MEAN) then
            f = f - shift*pw
         end if
         inc = inv_scale*f
         hTr(i, j, kt) = hTr(i, j, kt) + inc
         budget(i, j, kt) = budget(i, j, kt) + inc
      end do
   end subroutine sss_deposit_impl

   subroutine ocean_sss_restore_apply(grid, metrics, sf, ms, dt, active)
      !! One thermo-step application of the file-based SSS restoring.
      !! No-op unless `sf%has_restore_S_file` (and the thermo cadence gate,
      !! when given).  Collective across compute ranks when the net
      !! adjustment is on (one packed EFP allreduce).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
         !! `areaT`, and `cover_frac` when `use_cavity`.
      type(ocean_surface_flux_t), intent(in) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active

      integer :: nx, ny, nz, ng, idx_S
      logical :: use_cover
      real(wp) :: coef, p, n, a_open, fpos, fneg, shift, fmax_loc, fmax
      type(efp_t) :: loc(3), glob(3)

      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. sf%has_restore_S_file) return
      if (.not. allocated(ms%tracers)) return
      idx_S = ms%idx_salinity
      if (idx_S <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      ng = grid%nghost
      use_cover = metrics%use_cavity
      coef = sf%rho0*sf%sss_piston
      fpos = 1.0_wp
      fneg = 1.0_wp
      shift = 0.0_wp
      p = 0.0_wp
      n = 0.0_wp

      if (sf%sss_net_mode /= SSS_NET_NONE) then
         call sums_all(sf, ms, metrics, idx_S, use_cover, nx, ny, nz, ng, coef, loc)
         call halo_allreduce_efp_list(loc, glob, 3)
         p = efp_to_real(glob(SUM_POS))
         n = efp_to_real(glob(SUM_NEG))
         a_open = efp_to_real(glob(SUM_AREA))
         if (sf%sss_net_mode == SSS_NET_SCALE) then
            call sss_net_factors(p, n, fpos, fneg)
         else if (a_open > 0.0_wp) then
            shift = (p + n)/a_open
         end if
      end if

      call deposit_dispatch(sf, ms, metrics, idx_S, use_cover, nx, ny, nz, coef, &
                            fpos, fneg, shift, dt/sf%rho0)

      n_calls = n_calls + 1
      if (n_calls == 1 .or. mod(n_calls, 500) == 0) then
         ! Sparse progress line; the max is global, the line is rank-0 only.
         if (use_cover) then
            fmax_loc = sss_flux_maxabs(nx, ny, nz, ng, ms%tracers(idx_S)%hTr, ms%h_layer, &
                                       sf%sss_target, ms%wet_mask, metrics%cover_frac, &
                                       .true., ms%k_top, coef, sf%sss_max_delta, sf%h_min, &
                                       sf%sss_convert_sp)
         else
            fmax_loc = sss_flux_maxabs(nx, ny, nz, ng, ms%tracers(idx_S)%hTr, ms%h_layer, &
                                       sf%sss_target, ms%wet_mask, ms%wet_mask, &
                                       .false., ms%k_top, coef, sf%sss_max_delta, sf%h_min, &
                                       sf%sss_convert_sp)
         end if
         call halo_allreduce_max(fmax_loc, fmax)
         if (grid%i_offset_global == 0 .and. grid%j_offset_global == 0) then
            call logger%info("sss_restore: call "//to_string(n_calls)//"  max|F| = "// &
                             to_string(fmax)//"  P = "//to_string(p)//"  N = "// &
                             to_string(n)//"  fpos = "//to_string(fpos)//"  fneg = "// &
                             to_string(fneg))
         end if
      end if
   end subroutine ocean_sss_restore_apply

   subroutine sums_all(sf, ms, metrics, idx_S, use_cover, nx, ny, nz, ng, coef, loc)
      !! The three local EFP integrals, with the tracer-registry deref on
      !! the host (outer-shim rule).
      type(ocean_surface_flux_t), intent(in) :: sf
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: idx_S, nx, ny, nz, ng
      logical, intent(in) :: use_cover
      real(wp), intent(in) :: coef
      type(efp_t), intent(out) :: loc(3)
      integer :: m
      do m = 1, 3
         if (use_cover) then
            loc(m) = sss_flux_sums(m, nx, ny, nz, ng, ms%tracers(idx_S)%hTr, ms%h_layer, &
                                   sf%sss_target, ms%wet_mask, metrics%cover_frac, .true., &
                                   metrics%areaT, ms%k_top, coef, sf%sss_max_delta, &
                                   sf%h_min, sf%sss_convert_sp)
         else
            loc(m) = sss_flux_sums(m, nx, ny, nz, ng, ms%tracers(idx_S)%hTr, ms%h_layer, &
                                   sf%sss_target, ms%wet_mask, ms%wet_mask, .false., &
                                   metrics%areaT, ms%k_top, coef, sf%sss_max_delta, &
                                   sf%h_min, sf%sss_convert_sp)
         end if
      end do
   end subroutine sums_all

   subroutine deposit_dispatch(sf, ms, metrics, idx_S, use_cover, nx, ny, nz, coef, &
                               fpos, fneg, shift, inv_scale)
      !! Host-side shim for `sss_deposit_impl` (registry deref on host).
      type(ocean_surface_flux_t), intent(in) :: sf
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: idx_S, nx, ny, nz
      logical, intent(in) :: use_cover
      real(wp), intent(in) :: coef, fpos, fneg, shift, inv_scale
      if (use_cover) then
         call sss_deposit_impl(nx, ny, nz, ms%tracers(idx_S)%hTr, ms%salt_budget_surface, &
                               ms%h_layer, sf%sss_target, ms%wet_mask, metrics%cover_frac, &
                               .true., ms%k_top, coef, sf%sss_max_delta, sf%h_min, &
                               sf%sss_convert_sp, sf%sss_net_mode, fpos, fneg, shift, &
                               inv_scale)
      else
         call sss_deposit_impl(nx, ny, nz, ms%tracers(idx_S)%hTr, ms%salt_budget_surface, &
                               ms%h_layer, sf%sss_target, ms%wet_mask, ms%wet_mask, &
                               .false., ms%k_top, coef, sf%sss_max_delta, sf%h_min, &
                               sf%sss_convert_sp, sf%sss_net_mode, fpos, fneg, shift, &
                               inv_scale)
      end if
   end subroutine deposit_dispatch

end module rdb_ocean_sss_restore
