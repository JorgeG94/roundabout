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
   !! Three regions per sample, one device pass (three `reduction(+:)`
   !! accumulators):
   !!   full — all physical cells;
   !!   rim  — cells within `rim_band` of any physical wall (the slow
   !!          exponential is rim-trapped);
   !!   east — cells within `east_band` of the east wall (the mode-1
   !!          internal Kelvin wave is wall-trapped there).
   !!
   !! Gated by `&ocean_bt_nml debug_ke_attr` with an optional
   !! `[ke_attr_start_step, ke_attr_end_step]` outer-step window.
   !! Default off ⇒ bit-identical (untaken branches only).  Single-rank
   !! debug tool: no halo reduction — multi-rank sums are per-rank.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   implicit none
   private

   public :: ke_probe_t, ke_probe_sample, ke_probe_active
   public :: ke_probe_coradv_split

   type :: ke_probe_t
      !! Config + running state for the KE attribution meter.  Lives on
      !! `ocean_dyn_t`; host-only (never device-mapped).
      logical :: enable = .false.
         !! Master gate — `&ocean_bt_nml debug_ke_attr`.
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

   subroutine ke_probe_sample(grid, ms, probe, label, stage, step)
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

      real(wp) :: ke_full, ke_rim, ke_east
      real(wp) :: d_full, d_rim, d_east

      if (.not. ke_probe_active(probe, step)) return

      !$acc wait
      call ke_regions_impl(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                           grid%nx_total, grid%ny_total, ms%nz_ml, &
                           grid%nghost, grid%nx_phys, grid%ny_phys, &
                           probe%rim_band, probe%east_band, &
                           ke_full, ke_rim, ke_east)

      if (probe%have_prev) then
         d_full = ke_full - probe%prev_full
         d_rim = ke_rim - probe%prev_rim
         d_east = ke_east - probe%prev_east
      else
         d_full = ke_full
         d_rim = ke_rim
         d_east = ke_east
      end if

      if (.not. probe%header_done) then
         write (*, "(a)") "# KE_ATTR columns: step stage segment "// &
            "dKE_full KE_full dKE_rim KE_rim dKE_east KE_east   [m3/s2]"
         write (*, "(a)") "# dKE = energy injected (+) / removed (-) by the segment just run"
         probe%header_done = .true.
      end if
      write (*, "(a,1x,i6,1x,a,i1,1x,a16,6(1x,es13.5))") &
         "KE_ATTR", step, "s", stage, adjustl(label), &
         d_full, ke_full, d_rim, ke_rim, d_east, ke_east

      probe%prev_full = ke_full
      probe%prev_rim = ke_rim
      probe%prev_east = ke_east
      probe%have_prev = .true.
   end subroutine ke_probe_sample

   subroutine ke_probe_coradv_split(grid, metrics, ms, cor, probe, stage, step, dt)
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

      real(wp) :: wpv_full, wpv_rim, wpv_east
      real(wp) :: wke_full, wke_rim, wke_east

      if (.not. ke_probe_active(probe, step)) return

      !$acc wait
      call coradv_split_impl(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                             cor%pv_flux_x%data, cor%pv_flux_y%data, &
                             cor%ke_centre%data, metrics%idxCu, metrics%idyCv, &
                             grid%nx_total, grid%ny_total, ms%nz_ml, &
                             grid%nghost, grid%nx_phys, grid%ny_phys, &
                             probe%rim_band, probe%east_band, dt, &
                             wpv_full, wpv_rim, wpv_east, &
                             wke_full, wke_rim, wke_east)

      write (*, "(a,1x,i6,1x,a,i1,6(1x,es13.5))") &
         "KE_ATTR_PV", step, "s", stage, &
         wpv_full, wke_full, wpv_rim, wke_rim, wpv_east, wke_east
   end subroutine ke_probe_coradv_split

   subroutine coradv_split_impl(h, u, v, fx, fy, kec, idxCu, idyCv, &
                                nx, ny, nz, nghost, nx_phys, ny_phys, &
                                rim_w, east_w, dt, &
                                wpv_full, wpv_rim, wpv_east, &
                                wke_full, wke_rim, wke_east)
      !! Device pass: per-face PV/∇KE work split, three regions.
      !! Faces attributed to their west/south cell's region mask.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: rim_w, east_w
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      real(wp), intent(in) :: fx(nx + 1, ny, nz)
      real(wp), intent(in) :: fy(nx, ny + 1, nz)
      real(wp), intent(in) :: kec(nx, ny, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: dt
      real(wp), intent(out) :: wpv_full, wpv_rim, wpv_east
      real(wp), intent(out) :: wke_full, wke_rim, wke_east

      integer :: i, j, k, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h_face, f_tot, f_ke, f_pv, u_mid, w_pv, w_ke
      logical :: in_rim, in_east

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      wpv_full = 0.0_wp
      wpv_rim = 0.0_wp
      wpv_east = 0.0_wp
      wke_full = 0.0_wp
      wke_rim = 0.0_wp
      wke_east = 0.0_wp

      ! u-faces: interior faces i in [i_lo+1, i_hi]; region from the
      ! west cell (i-1 ≥ i_lo).
      do concurrent(k=1:nz, j=j_lo:j_hi, i=i_lo + 1:i_hi) &
         reduce(+:wpv_full, wpv_rim, wpv_east, wke_full, wke_rim, wke_east)
         h_face = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
         f_tot = fx(i, j, k)
         f_ke = -(kec(i, j, k) - kec(i - 1, j, k))*idxCu(i, j)
         f_pv = f_tot - f_ke
         u_mid = u(i, j, k) - 0.5_wp*dt*f_tot
         w_pv = h_face*u_mid*dt*f_pv
         w_ke = h_face*u_mid*dt*f_ke
         wpv_full = wpv_full + w_pv
         wke_full = wke_full + w_ke
         in_rim = (i - 1 - i_lo < rim_w .or. i_hi - (i - 1) < rim_w .or. &
                   j - j_lo < rim_w .or. j_hi - j < rim_w)
         in_east = (i_hi - (i - 1) < east_w)
         if (in_rim) then
            wpv_rim = wpv_rim + w_pv
            wke_rim = wke_rim + w_ke
         end if
         if (in_east) then
            wpv_east = wpv_east + w_pv
            wke_east = wke_east + w_ke
         end if
      end do

      ! v-faces: interior faces j in [j_lo+1, j_hi]; region from the
      ! south cell.
      do concurrent(k=1:nz, j=j_lo + 1:j_hi, i=i_lo:i_hi) &
         reduce(+:wpv_full, wpv_rim, wpv_east, wke_full, wke_rim, wke_east)
         h_face = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
         f_tot = fy(i, j, k)
         f_ke = -(kec(i, j, k) - kec(i, j - 1, k))*idyCv(i, j)
         f_pv = f_tot - f_ke
         u_mid = v(i, j, k) - 0.5_wp*dt*f_tot
         w_pv = h_face*u_mid*dt*f_pv
         w_ke = h_face*u_mid*dt*f_ke
         wpv_full = wpv_full + w_pv
         wke_full = wke_full + w_ke
         in_rim = (i - i_lo < rim_w .or. i_hi - i < rim_w .or. &
                   j - 1 - j_lo < rim_w .or. j_hi - (j - 1) < rim_w)
         in_east = (i_hi - i < east_w)
         if (in_rim) then
            wpv_rim = wpv_rim + w_pv
            wke_rim = wke_rim + w_ke
         end if
         if (in_east) then
            wpv_east = wpv_east + w_pv
            wke_east = wke_east + w_ke
         end if
      end do
   end subroutine coradv_split_impl

   subroutine ke_regions_impl(h, u, v, nx, ny, nz, nghost, nx_phys, ny_phys, &
                              rim_w, east_w, ke_full, ke_rim, ke_east)
      !! h-weighted layer KE over the three regions, one device pass.
      !! Per cell: KE = Σ_k [ ¼(h_W·u_W² + h_E·u_E²) + ¼(h_S·v_S² + h_N·v_N²) ]
      !! with h_face the centred face thickness (m³/s² per unit area;
      !! uniform-Δx Cartesian, so the area factor is a constant and is
      !! dropped — attribution only needs a consistent measure).
      !! Explicit OpenACC reduction — a `sum()` intrinsic on a
      !! present-mapped array silently runs host-side under NVHPC
      !! non-managed mode and returns the stale host shadow.
      integer, intent(in) :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer, intent(in) :: rim_w, east_w
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: u(nx + 1, ny, nz)
      real(wp), intent(in) :: v(nx, ny + 1, nz)
      real(wp), intent(out) :: ke_full, ke_rim, ke_east

      integer :: i, j, k, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h_w, h_e, h_s, h_n, ke_cell

      i_lo = nghost + 1
      i_hi = nghost + nx_phys
      j_lo = nghost + 1
      j_hi = nghost + ny_phys

      ke_full = 0.0_wp
      ke_rim = 0.0_wp
      ke_east = 0.0_wp
      do concurrent(k=1:nz, j=j_lo:j_hi, i=i_lo:i_hi) reduce(+:ke_full, ke_rim, ke_east)
         h_w = 0.5_wp*(h(i - 1, j, k) + h(i, j, k))
         h_e = 0.5_wp*(h(i, j, k) + h(i + 1, j, k))
         h_s = 0.5_wp*(h(i, j - 1, k) + h(i, j, k))
         h_n = 0.5_wp*(h(i, j, k) + h(i, j + 1, k))
         ke_cell = 0.25_wp*(h_w*u(i, j, k)*u(i, j, k) + &
                            h_e*u(i + 1, j, k)*u(i + 1, j, k)) + &
                   0.25_wp*(h_s*v(i, j, k)*v(i, j, k) + &
                            h_n*v(i, j + 1, k)*v(i, j + 1, k))
         ke_full = ke_full + ke_cell
         if (i - i_lo < rim_w .or. i_hi - i < rim_w .or. &
             j - j_lo < rim_w .or. j_hi - j < rim_w) then
            ke_rim = ke_rim + ke_cell
         end if
         if (i_hi - i < east_w) ke_east = ke_east + ke_cell
      end do
   end subroutine ke_regions_impl

end module rdb_ocean_ke_probe
