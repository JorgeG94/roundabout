!! Per-segment layer-KE attribution probe ("kediag" rebuilt).
module rdb_ocean_ke_probe
   !! Rebuild of the transient per-term KE attribution meter used to
   !! localise the split-corrector anti-damping term
   !! (LAGRANGIAN_PGF_BUG.md §R).  Samples the h-weighted layer
   !! kinetic energy between the split stage's tendency applies and
   !! prints `KE_ATTR` rows: the dKE between consecutive samples IS the
   !! energy injected/removed by the segment just executed, so one
   !! instrumented run names the anti-damping term directly.
   !!
   !! Three regions per sample, one small device pass each:
   !!   full — all physical cells;
   !!   rim  — cells within `rim_band` of any GLOBAL domain wall (the slow
   !!          exponential is rim-trapped);
   !!   east — cells within `east_band` of the GLOBAL east wall (the mode-1
   !!          internal Kelvin wave is wall-trapped there).
   !!
   !! Multi-rank: every region is defined in GLOBAL cell indices
   !! (`i_offset_global`/`nx_global`).  The KE is summed PER OWNED FACE
   !! (see `ke_region_one`): each u/v face is read by exactly one rank, so
   !! NO u/v ghost is ever read -- on a PERIODIC axis the east (north) wrap
   !! slot is a copy of face 1 that the full-range applies advance on their
   !! own, so it is not read either; face 1 then also belongs to cell
   !! `nx_g` (`ny_g`) through the wrap.  The only ghosts the probe relies on
   !! are those of `h_layer` (the face thickness straddles a seam, and at a
   !! periodic face 1 it needs cell `nx_g`), current at both sample points
   !! (the entry sample follows the stage-entry exchange; at stage close
   !! `h_layer` was last set by the continuity chain, which ends in a halo
   !! exchange and periodic wrap).  Caveat: under `eulerian_z` + `ssp_rk2`
   !! the vertical advection updates ghost thickness locally AFTER that
   !! exchange, so the seam-face thickness at the stage-close sample is
   !! layout-dependent there (known defect, fixed separately).  Tripolar
   !! fold: the v row at `jg = ny_g+1` (the fold line, owned by the north
   !! rank row) is read AS STORED, like the one-rank run; whether it is
   !! refreshed before the stage-close sample is not verified.  Sums are reproducing EFP sums,
   !! combined with `halo_allreduce_efp_list` (collective: all ranks call
   !! whenever the step-window gate is open) and printed once by rank 0,
   !! so the digits are identical on any rank count.
   !!
   !! Gated by `&ocean_debug_nml ke_attr` with an optional
   !! `[ke_attr_start_step, ke_attr_end_step]` outer-step window.
   !! Default off ⇒ bit-identical (untaken branches only); model state is
   !! never written.  `ke_probe_coradv_split` (KE_ATTR_PV rows) has NO
   !! caller in the stepping loop today (it is only imported there); it is
   !! decomposition-safe as a routine but additionally needs the Coriolis
   !! scratch (`ke_centre`) ghosts to be current when called.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use rdb_constants, only: wp
   use rdb_efp, only: efp_t, efp_carry, efp_to_real, EFP_DIGITS
   use rdb_halo, only: halo_allreduce_efp_list
   use rdb_comm_env, only: comm_env_rank
   use rdb_ocean_console_stats, only: efp_decompose_impl, efp_summands_guard
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   implicit none
   private

   public :: ke_probe_t, ke_probe_sample, ke_probe_active
   public :: ke_probe_coradv_split
   public :: ke_regions_efp_local, coradv_split_efp_local, N_KE, N_SPLIT

   integer, parameter :: N_KE = 3
      !! Regions per KE sample: full, rim, east.
   integer, parameter :: N_SPLIT = 6
      !! PV/KE-work split values: (PV, KE) x (full, rim, east).

   type :: ke_probe_t
      !! Config + running state for the KE attribution meter.  Lives on
      !! `ocean_dyn_t`; host-only (never device-mapped).
      logical :: enable = .false.
         !! Master gate — `&ocean_debug_nml ke_attr`.
      integer :: start_step = 0
         !! First outer step to sample (0 = from the start).
      integer :: end_step = 0
         !! Last outer step to sample (0 = no upper bound).
      integer :: rim_band = 8
         !! Rim-region width in cells (distance from any physical wall).
      integer :: east_band = 16
         !! East-wall band width in cells (the Kelvin-wave region; 16
         !! covers the ~6-cell internal-radius wall trapping).
      logical :: have_prev = .false.
         !! .true. once a sample has landed; first row prints dKE = KE.
      logical :: header_done = .false.
         !! Column-header row emitted.
      real(wp) :: prev_full = 0.0_wp
      real(wp) :: prev_rim = 0.0_wp
      real(wp) :: prev_east = 0.0_wp
         !! KE at the previous sample, per region.
   end type ke_probe_t

contains

   pure function ke_probe_active(probe, step) result(active)
      !! Gate: enabled AND inside the step window.
      type(ke_probe_t), intent(in) :: probe
      integer, intent(in) :: step
      logical :: active
      active = probe%enable
      if (active .and. probe%start_step > 0) active = step >= probe%start_step
      if (active .and. probe%end_step > 0) active = step <= probe%end_step
   end function ke_probe_active

   subroutine ke_probe_sample(grid, ms, probe, label, stage, step, periodic_x, periodic_y)
      !! Sample the three-region layer KE and print the attribution row.
      !! Call AFTER the segment being attributed; the printed dKE is
      !! "this sample minus the previous one" — i.e. the segment's
      !! energy contribution.  Drains all device queues first so the
      !! async velocity-apply chain (queue 1) has landed.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(ke_probe_t), intent(inout) :: probe
      character(len=*), intent(in) :: label
         !! Segment tag, e.g. "coriolis_adv", "bt_correction".
      integer, intent(in) :: stage, step
      logical, intent(in) :: periodic_x, periodic_y
         !! GLOBAL periodic-axis flags (bc state tags, rank-uniform; not `has_*`).

      real(wp) :: ke_full, ke_rim, ke_east
      real(wp) :: d_full, d_rim, d_east
      type(efp_t) :: loc(N_KE), glob(N_KE)

      if (.not. ke_probe_active(probe, step)) return

      !$acc wait
      call ke_regions_efp_local(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                                grid%nx_total, grid%ny_total, ms%nz_ml, &
                                grid%nghost, grid%nx_phys, grid%ny_phys, &
                                grid%i_offset_global, grid%j_offset_global, &
                                grid%nx_global, grid%ny_global, &
                                probe%rim_band, probe%east_band, periodic_x, periodic_y, loc)
      call halo_allreduce_efp_list(loc, glob, N_KE)
      ke_full = real(efp_to_real(glob(1)), wp)
      ke_rim = real(efp_to_real(glob(2)), wp)
      ke_east = real(efp_to_real(glob(3)), wp)

      if (probe%have_prev) then
         d_full = ke_full - probe%prev_full
         d_rim = ke_rim - probe%prev_rim
         d_east = ke_east - probe%prev_east
      else
         d_full = ke_full
         d_rim = ke_rim
         d_east = ke_east
      end if

      if (comm_env_rank() == 0) then
         if (.not. probe%header_done) then
            write (*, "(a)") "# KE_ATTR columns: step stage segment "// &
               "dKE_full KE_full dKE_rim KE_rim dKE_east KE_east   [m3/s2]"
            write (*, "(a)") "# dKE = energy injected (+) / removed (-) by the segment just run"
            probe%header_done = .true.
         end if
         write (*, "(a,1x,i6,1x,a,i1,1x,a16,6(1x,es13.5))") &
            "KE_ATTR", step, "s", stage, adjustl(label), &
            d_full, ke_full, d_rim, ke_rim, d_east, ke_east
      end if

      probe%prev_full = ke_full
      probe%prev_rim = ke_rim
      probe%prev_east = ke_east
      probe%have_prev = .true.
   end subroutine ke_probe_sample

   subroutine ke_probe_coradv_split(grid, metrics, ms, cor, probe, stage, step, dt, periodic_x, periodic_y)
      !! Sub-attribute the coriolis_adv segment: split the just-applied
      !! tendency into its PV-flux part and its −∇KE part and print each
      !! part's work integral per region (`KE_ATTR_PV` rows).  Call
      !! immediately AFTER `coriolis_adv_apply_tendencies` (and after
      !! the "coriolis_adv" KE_ATTR sample).
      !!
      !! Method: the kernel left `pv_flux_*` = F_tot = F_pv − ∇KE and
      !! `ke_centre` in its scratch buffers.  Per face,
      !! F_ke = −∇KE (recomputed from ke_centre), F_pv = F_tot − F_ke,
      !! and the work of each part over the apply is
      !! `W_part = h_face · u_mid · dt · F_part` with
      !! u_mid = u_new − ½·dt·F_tot (u_new is the post-apply face
      !! velocity), so W_pv + W_ke equals the segment's dKE up to the
      !! face/cell weighting convention (a few % — use the split, not
      !! the sum, as the signal).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
      type(coriolis_adv_t), intent(in) :: cor
      type(ke_probe_t), intent(inout) :: probe
      integer, intent(in) :: stage, step
      real(wp), intent(in) :: dt
      logical, intent(in) :: periodic_x, periodic_y
         !! GLOBAL periodic-axis flags (bc state tags).

      type(efp_t) :: loc(N_SPLIT), glob(N_SPLIT)

      if (.not. ke_probe_active(probe, step)) return

      !$acc wait
      call coradv_split_efp_local(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                                  cor%pv_flux_x%data, cor%pv_flux_y%data, &
                                  cor%ke_centre%data, metrics%idxCu, metrics%idyCv, &
                                  grid%nx_total, grid%ny_total, ms%nz_ml, &
                                  grid%nghost, grid%nx_phys, grid%ny_phys, &
                                  grid%i_offset_global, grid%j_offset_global, &
                                  grid%nx_global, grid%ny_global, &
                                  probe%rim_band, probe%east_band, dt, periodic_x, periodic_y, loc)
      call halo_allreduce_efp_list(loc, glob, N_SPLIT)

      if (comm_env_rank() == 0) then
         write (*, "(a,1x,i6,1x,a,i1,6(1x,es13.5))") &
            "KE_ATTR_PV", step, "s", stage, &
            real(efp_to_real(glob(1)), wp), real(efp_to_real(glob(4)), wp), &
            real(efp_to_real(glob(2)), wp), real(efp_to_real(glob(5)), wp), &
            real(efp_to_real(glob(3)), wp), real(efp_to_real(glob(6)), wp)
      end if
   end subroutine ke_probe_coradv_split

   subroutine ke_regions_efp_local(h, u, v, nx, ny, nz, nghost, nx_phys, ny_phys, &
                                   i_off, j_off, nx_g, ny_g, rim_w, east_w, &
                                   periodic_x, periodic_y, tot)
      !! RANK-LOCAL EFP partial sums of the h-weighted layer KE over the
      !! three regions (`tot(1:3)` = full, rim, east); the caller combines
      !! them across ranks.  Per cell: KE = ¼(h_W·u_W² + h_E·u_E²) +
      !! ¼(h_S·v_S² + h_N·v_N²) with h_face the centred face thickness
      !! (area factor dropped; attribution only needs a consistent
      !! measure).  Regions are tested on GLOBAL cell indices
      !! (`i_off`/`j_off` = the tile's global offsets, `nx_g`/`ny_g` the
      !! global extents), so they are decomposition-independent.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: i_off, j_off, nx_g, ny_g, rim_w, east_w
      logical, intent(in) :: periodic_x, periodic_y
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      type(efp_t), intent(out) :: tot(N_KE)

      integer :: r
      do r = 1, N_KE
         call ke_region_one(h, u, v, nx, ny, nz, nghost, nx_phys, ny_phys, &
                            i_off, j_off, nx_g, ny_g, rim_w, east_w, periodic_x, periodic_y, r, tot(r))
      end do
   end subroutine ke_regions_efp_local

   pure function in_region(region, ic, jc, nx_g, ny_g, rim_w, east_w) result(r)
      !! Is GLOBAL cell (ic, jc) in the region (1 full, 2 rim, 3 east)?
      !! Cells outside the domain are in no region.
      !$acc routine seq
      integer, intent(in) :: region, ic, jc, nx_g, ny_g, rim_w, east_w
      logical :: r
      r = (ic >= 1 .and. ic <= nx_g .and. jc >= 1 .and. jc <= ny_g)
      if (r .and. region == 2) then
         r = (ic - 1 < rim_w .or. nx_g - ic < rim_w .or. jc - 1 < rim_w .or. ny_g - jc < rim_w)
      else if (r .and. region == 3) then
         r = (nx_g - ic < east_w)
      end if
   end function in_region

   subroutine ke_region_one(h, u, v, nx, ny, nz, nghost, nx_phys, ny_phys, &
                            i_off, j_off, nx_g, ny_g, rim_w, east_w, &
                            periodic_x, periodic_y, region, total)
      !! One region (1 full, 2 rim, 3 east) of `ke_regions_efp_local`,
      !! summed PER OWNED FACE.  The per-cell KE
      !! `¼(h_W u_W² + h_E u_E²) + ¼(h_S v_S² + h_N v_N²)` is regrouped as
      !! `Σ_faces ¼·h_face·u_face²·n`, `n` = number of the face's two
      !! adjacent cells that lie in the region (0..2), so the global sum is
      !! algebraically the per-cell one.  A face is read only by the rank
      !! that owns it: west faces (global index 1..nx_g) of its cells, plus
      !! the global east boundary face (`ig = nx_g+1`) on the east-most
      !! rank; v likewise with south faces / the global north face.  No u/v
      !! GHOST is ever read (they are stale at the stage-close sample).
      !! PERIODIC axis (`periodic_x`/`periodic_y`, the GLOBAL edge tags): the
      !! east (north) boundary face slot is a wrap copy of face 1, advanced
      !! on its own between exchanges, so it is NOT read; face 1 is the
      !! east face of cell `nx_g` (`ny_g`), so the region test wraps the
      !! non-existent cell 0 to `nx_g` and face 1 gets its weight from
      !! cells `nx_g` and 1.  Only owned prognostic faces are read.
      !! `h` ghosts are read (the face thickness straddles a seam; at a
      !! periodic face 1 it needs the h of cell `nx_g`, the west-most
      !! rank's west ghost, a wrap/halo copy that is current for the same
      !! reason as every other h ghost): `h` is current at both sample
      !! points -- the stage-entry sample follows the stage-entry
      !! exchange, and at stage close `h_layer` was last updated by the
      !! continuity chain, which ends in a halo exchange + periodic wrap
      !! (`ssp_rk2`), or not at all before the sample (`pred_corr` runs
      !! continuity after it).  Not so under `eulerian_z` + `ssp_rk2`: the
      !! vertical advection then updates ghost thickness locally after
      !! that exchange, so the seam-face thickness is layout-dependent
      !! there (known defect, fixed separately).  Tripolar fold: the v row
      !! `jg = ny_g+1` is the fold line, read as stored (refresh before the
      !! stage-close sample not verified).
      !! One k-slab at a time, seven EFP accumulators, carried between slabs.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: i_off, j_off, nx_g, ny_g, rim_w, east_w, region
      logical, intent(in) :: periodic_x, periodic_y
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      type(efp_t), intent(out) :: total

      integer :: i, j, k, ig, jg, nw, icw, jcs
      integer(int64) :: e1, e2, e3, e4, e5, e6, epoison
      integer(int64) :: d1, d2, d3, d4, d5, d6, dpoison
      integer(int64) :: slab_e(EFP_DIGITS)
      real(real64) :: val
      real(wp) :: hfc

      call efp_summands_guard(2*(nx_phys + 1), ny_phys + 1, "ke_region_one")
      total%v = 0_int64
      total%poison = 0_int64
      do k = 1, nz
         e1 = 0_int64
         e2 = 0_int64
         e3 = 0_int64
         e4 = 0_int64
         e5 = 0_int64
         e6 = 0_int64
         epoison = 0_int64
         !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6,epoison) &
         !$acc&         private(val, hfc, ig, jg, nw, icw, jcs, d1, d2, d3, d4, d5, d6, dpoison) &
         !$acc&         present(h, u, v)
         do j = nghost + 1, nghost + ny_phys + 1
            do i = nghost + 1, nghost + nx_phys + 1
               ig = i - nghost + i_off
               jg = j - nghost + j_off
               ! u face (west face of cell ig): owned rows only.
               val = 0.0_real64
               if (j <= nghost + ny_phys .and. &
                   (i <= nghost + nx_phys .or. (ig == nx_g + 1 .and. .not. periodic_x))) then
                  nw = 0
                  icw = ig - 1
                  if (periodic_x .and. icw < 1) icw = nx_g
                  if (in_region(region, icw, jg, nx_g, ny_g, rim_w, east_w)) nw = nw + 1
                  if (in_region(region, ig, jg, nx_g, ny_g, rim_w, east_w)) nw = nw + 1
                  if (nw > 0) then
                     hfc = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
                     val = real(0.25_wp*nw*hfc*u(i, j, k)*u(i, j, k), real64)
                  end if
               end if
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6, dpoison)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
               epoison = epoison + dpoison
               ! v face (south face of cell jg): owned columns only.
               val = 0.0_real64
               if (i <= nghost + nx_phys .and. &
                   (j <= nghost + ny_phys .or. (jg == ny_g + 1 .and. .not. periodic_y))) then
                  nw = 0
                  jcs = jg - 1
                  if (periodic_y .and. jcs < 1) jcs = ny_g
                  if (in_region(region, ig, jcs, nx_g, ny_g, rim_w, east_w)) nw = nw + 1
                  if (in_region(region, ig, jg, nx_g, ny_g, rim_w, east_w)) nw = nw + 1
                  if (nw > 0) then
                     hfc = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
                     val = real(0.25_wp*nw*hfc*v(i, j, k)*v(i, j, k), real64)
                  end if
               end if
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6, dpoison)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
               epoison = epoison + dpoison
            end do
         end do
         slab_e = [e1, e2, e3, e4, e5, e6]
         call efp_carry(slab_e)
         total%v = total%v + slab_e
         total%poison = total%poison + epoison
         call efp_carry(total%v)
      end do
   end subroutine ke_region_one

   subroutine coradv_split_efp_local(h, u, v, fx, fy, kec, idxCu, idyCv, &
                                     nx, ny, nz, nghost, nx_phys, ny_phys, &
                                     i_off, j_off, nx_g, ny_g, rim_w, east_w, dt, periodic_x, periodic_y, tot)
      !! RANK-LOCAL EFP partials of the per-face PV / -grad(KE) work
      !! split: `tot(1:3)` = PV work (full, rim, east), `tot(4:6)` = KE
      !! work (full, rim, east).  Faces are attributed to their west/south
      !! cell's region (global indices) and owned by the cell on their
      !! east/north side; the global west/south wall faces are excluded
      !! (on a periodic axis they are face 1 of the wrap: included, with
      !! the region taken from cell `nx_g` / `ny_g`), so a face is counted
      !! once on any decomposition.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: i_off, j_off, nx_g, ny_g, rim_w, east_w
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      real(wp), intent(in) :: fx(nx + 1, ny, nz)
      real(wp), intent(in) :: fy(nx, ny + 1, nz)
      real(wp), intent(in) :: kec(nx, ny, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: dt
      logical, intent(in) :: periodic_x, periodic_y
      type(efp_t), intent(out) :: tot(N_SPLIT)

      integer :: part, region
      do part = 1, 2
         do region = 1, N_KE
            call coradv_one(h, u, v, fx, fy, kec, idxCu, idyCv, nx, ny, nz, &
                            nghost, nx_phys, ny_phys, i_off, j_off, nx_g, ny_g, &
                            rim_w, east_w, dt, periodic_x, periodic_y, part, region, tot(N_KE*(part - 1) + region))
         end do
      end do
   end subroutine coradv_split_efp_local

   subroutine coradv_one(h, u, v, fx, fy, kec, idxCu, idyCv, nx, ny, nz, &
                         nghost, nx_phys, ny_phys, i_off, j_off, nx_g, ny_g, &
                         rim_w, east_w, dt, periodic_x, periodic_y, part, region, total)
      !! One (part: 1 PV / 2 KE, region: 1 full / 2 rim / 3 east) EFP sum
      !! of `coradv_split_efp_local`; seven accumulators per k-slab.
      !! Reads only OWNED faces (u west faces ig=2..nx_g of its cells, v
      !! south faces jg=2..ny_g; ig=1 / jg=1 too on a periodic axis, with the
      !! region from the wrapped cell nx_g / ny_g) and, at a seam, the `h` and `ke_centre`
      !! ghost of the neighbouring cell -- `kec` is the Coriolis kernel's
      !! scratch, current only if the caller samples right after it
      !! (this routine has no caller in the stepping loop today).
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: i_off, j_off, nx_g, ny_g, rim_w, east_w, part, region
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      real(wp), intent(in) :: fx(nx + 1, ny, nz)
      real(wp), intent(in) :: fy(nx, ny + 1, nz)
      real(wp), intent(in) :: kec(nx, ny, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: dt
      logical, intent(in) :: periodic_x, periodic_y
      type(efp_t), intent(out) :: total

      integer :: i, j, k, ig, jg, icw, jcs
      integer(int64) :: e1, e2, e3, e4, e5, e6, epoison
      integer(int64) :: d1, d2, d3, d4, d5, d6, dpoison
      integer(int64) :: slab_e(EFP_DIGITS)
      real(real64) :: val
      real(wp) :: h_face, f_tot, f_ke, f_pv, u_mid, w
      logical :: in_reg

      call efp_summands_guard(2*nx_phys, ny_phys, "coradv_one")
      total%v = 0_int64
      total%poison = 0_int64
      do k = 1, nz
         e1 = 0_int64
         e2 = 0_int64
         e3 = 0_int64
         e4 = 0_int64
         e5 = 0_int64
         e6 = 0_int64
         epoison = 0_int64
         !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6,epoison) &
         !$acc&         private(val, h_face, f_tot, f_ke, f_pv, u_mid, w, ig, jg, icw, jcs, in_reg, d1, d2, d3, d4, d5, d6, dpoison) &
         !$acc&         present(h, u, v, fx, fy, kec, idxCu, idyCv)
         do j = nghost + 1, nghost + ny_phys
            do i = nghost + 1, nghost + nx_phys
               ig = i - nghost + i_off
               jg = j - nghost + j_off
               ! u-face (west face of this cell); region from the west cell ig-1.
               val = 0.0_real64
               if (ig > 1 .or. periodic_x) then
                  icw = ig - 1
                  if (icw < 1) icw = nx_g
                  in_reg = .true.
                  if (region == 2) then
                     in_reg = (icw - 1 < rim_w .or. nx_g - icw < rim_w .or. &
                               jg - 1 < rim_w .or. ny_g - jg < rim_w)
                  else if (region == 3) then
                     in_reg = (nx_g - icw < east_w)
                  end if
                  if (in_reg) then
                     h_face = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
                     f_tot = fx(i, j, k)
                     f_ke = -(kec(i, j, k) - kec(i - 1, j, k))*idxCu(i, j)
                     f_pv = f_tot - f_ke
                     u_mid = u(i, j, k) - 0.5_wp*dt*f_tot
                     w = h_face*u_mid*dt*f_pv
                     if (part == 2) w = h_face*u_mid*dt*f_ke
                     val = real(w, real64)
                  end if
               end if
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6, dpoison)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
               epoison = epoison + dpoison
               ! v-face (south face of this cell); region from the south cell jg-1.
               val = 0.0_real64
               if (jg > 1 .or. periodic_y) then
                  jcs = jg - 1
                  if (jcs < 1) jcs = ny_g
                  in_reg = .true.
                  if (region == 2) then
                     in_reg = (ig - 1 < rim_w .or. nx_g - ig < rim_w .or. &
                               jcs - 1 < rim_w .or. ny_g - jcs < rim_w)
                  else if (region == 3) then
                     in_reg = (nx_g - ig < east_w)
                  end if
                  if (in_reg) then
                     h_face = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
                     f_tot = fy(i, j, k)
                     f_ke = -(kec(i, j, k) - kec(i, j - 1, k))*idyCv(i, j)
                     f_pv = f_tot - f_ke
                     u_mid = v(i, j, k) - 0.5_wp*dt*f_tot
                     w = h_face*u_mid*dt*f_pv
                     if (part == 2) w = h_face*u_mid*dt*f_ke
                     val = real(w, real64)
                  end if
               end if
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6, dpoison)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
               epoison = epoison + dpoison
            end do
         end do
         slab_e = [e1, e2, e3, e4, e5, e6]
         call efp_carry(slab_e)
         total%v = total%v + slab_e
         total%poison = total%poison + epoison
         call efp_carry(total%v)
      end do
   end subroutine coradv_one

end module rdb_ocean_ke_probe
