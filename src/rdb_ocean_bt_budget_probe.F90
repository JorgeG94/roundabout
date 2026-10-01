!! Per-region BT-mode energy + per-term power diagnostic probe.
module rdb_ocean_bt_budget_probe
   !! `print_bt_budget` walks the per-kernel slow-tendency arrays before
   !! they sum into `F_slow_u/v`, decomposes the BT-mode work rate by
   !! source term (PGF, Coriolis-adv, hvisc, bottom drag, surface
   !! stress), splits the basin into a N/S region triple (subpolar / jet
   !! / subtropical for double_gyre), and prints a tabular per-region
   !! snapshot. Diagnoses which slow-term over/under-energizes a gyre.
   !!
   !! Gated by `bt_work%debug_bt_budget` (namelist `ocean_debug_bt_budget`).
   !! Default off ⇒ existing runs unchanged. Best paired with serial
   !! multicore builds so prints land in causal order.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   implicit none
   private

   public :: print_bt_budget

contains

   subroutine print_bt_budget(grid, ms, bt_work, pgf, cor, hv, bd, ss, &
                              t_value, t_unit, stage_tag, header)
      !! Compute + print the per-region BT-mode budget snapshot.
      !!
      !! Power terms `P_<term>` are computed as the cell-centred
      !! u_bt·F_u + v_bt·F_v dot product, then averaged over the
      !! cells in each region.  Units: m²/s³ (acceleration × velocity).
      !! Sign convention: positive = the term ADDS energy to the BT
      !! mode (the wind P_ss should be positive in the gyres; hvisc /
      !! drag should be negative).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(barotropic_workstate_t), intent(in) :: bt_work
      type(ocean_pressure_force_t), intent(in) :: pgf
      type(coriolis_adv_t), intent(in) :: cor
      type(ocean_horizontal_viscosity_t), intent(in) :: hv
      type(ocean_bottom_drag_t), intent(in) :: bd
      type(ocean_surface_stress_t), intent(in) :: ss
      real(wp), intent(in) :: t_value
      character(len=*), intent(in) :: t_unit, stage_tag
      logical, intent(in) :: header

      integer :: nx, ny, nz, nghost
      integer :: jp_lo, jp_mid_lo, jp_mid_hi, jp_hi
      integer :: nz_top
      real(wp) :: ke_sub, ke_jet, ke_pol
      real(wp) :: eta_sub_min, eta_sub_max, eta_jet_min, eta_jet_max
      real(wp) :: eta_pol_min, eta_pol_max
      real(wp) :: usurf_sub, usurf_jet, usurf_pol
      real(wp) :: ubed_sub, ubed_jet, ubed_pol
      real(wp) :: vsurf_sub, vsurf_jet, vsurf_pol
      real(wp) :: vbed_sub, vbed_jet, vbed_pol
      real(wp) :: p_pgf_sub, p_pgf_jet, p_pgf_pol
      real(wp) :: p_cor_sub, p_cor_jet, p_cor_pol
      real(wp) :: p_hv_sub, p_hv_jet, p_hv_pol
      real(wp) :: p_bd_sub, p_bd_jet, p_bd_pol
      real(wp) :: p_ss_sub, p_ss_jet, p_ss_pol

      ! ---- Region selection: thirds along j (NS) within the
      ! physical domain.  The double_gyre setup puts the spoon's
      ! east-west margin in the top third (subpolar gyre).
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nghost = grid%nghost
      nz_top = nz                            ! Surface layer index
      jp_lo = nghost + 1                     ! First interior j
      jp_hi = nghost + grid%ny_phys          ! Last interior j
      jp_mid_lo = nghost + grid%ny_phys/3 + 1
      jp_mid_hi = nghost + 2*grid%ny_phys/3

      ! ---- η ranges + |u|/|v| extremes per region ----
      call region_eta_uv(grid, ms, bt_work, &
                         jp_lo, jp_mid_lo - 1, &
                         eta_sub_min, eta_sub_max, &
                         usurf_sub, ubed_sub, vsurf_sub, vbed_sub, ke_sub)
      call region_eta_uv(grid, ms, bt_work, &
                         jp_mid_lo, jp_mid_hi, &
                         eta_jet_min, eta_jet_max, &
                         usurf_jet, ubed_jet, vsurf_jet, vbed_jet, ke_jet)
      call region_eta_uv(grid, ms, bt_work, &
                         jp_mid_hi + 1, jp_hi, &
                         eta_pol_min, eta_pol_max, &
                         usurf_pol, ubed_pol, vsurf_pol, vbed_pol, ke_pol)

      ! ---- Per-term BT power per region ----
      call region_power(grid, ms, bt_work, &
                        pgf%dpdx_face%data, pgf%dpdy_face%data, &
                        jp_lo, jp_mid_lo - 1, p_pgf_sub)
      call region_power(grid, ms, bt_work, &
                        pgf%dpdx_face%data, pgf%dpdy_face%data, &
                        jp_mid_lo, jp_mid_hi, p_pgf_jet)
      call region_power(grid, ms, bt_work, &
                        pgf%dpdx_face%data, pgf%dpdy_face%data, &
                        jp_mid_hi + 1, jp_hi, p_pgf_pol)

      call region_power(grid, ms, bt_work, &
                        cor%pv_flux_x%data, cor%pv_flux_y%data, &
                        jp_lo, jp_mid_lo - 1, p_cor_sub)
      call region_power(grid, ms, bt_work, &
                        cor%pv_flux_x%data, cor%pv_flux_y%data, &
                        jp_mid_lo, jp_mid_hi, p_cor_jet)
      call region_power(grid, ms, bt_work, &
                        cor%pv_flux_x%data, cor%pv_flux_y%data, &
                        jp_mid_hi + 1, jp_hi, p_cor_pol)

      call region_power(grid, ms, bt_work, &
                        hv%du_visc%data, hv%dv_visc%data, &
                        jp_lo, jp_mid_lo - 1, p_hv_sub)
      call region_power(grid, ms, bt_work, &
                        hv%du_visc%data, hv%dv_visc%data, &
                        jp_mid_lo, jp_mid_hi, p_hv_jet)
      call region_power(grid, ms, bt_work, &
                        hv%du_visc%data, hv%dv_visc%data, &
                        jp_mid_hi + 1, jp_hi, p_hv_pol)

      ! Bottom drag: explicit-tendency path writes `-r·u·(h_in_bbl/h)`
      ! to `du_drag/dv_drag` directly. (Implicit-rate variant handled by
      ! region_power_drag_implicit when used.)
      call region_power(grid, ms, bt_work, &
                        bd%du_drag%data, bd%dv_drag%data, &
                        jp_lo, jp_mid_lo - 1, p_bd_sub)
      call region_power(grid, ms, bt_work, &
                        bd%du_drag%data, bd%dv_drag%data, &
                        jp_mid_lo, jp_mid_hi, p_bd_jet)
      call region_power(grid, ms, bt_work, &
                        bd%du_drag%data, bd%dv_drag%data, &
                        jp_mid_hi + 1, jp_hi, p_bd_pol)

      call region_power(grid, ms, bt_work, &
                        ss%du_stress%data, ss%dv_stress%data, &
                        jp_lo, jp_mid_lo - 1, p_ss_sub)
      call region_power(grid, ms, bt_work, &
                        ss%du_stress%data, ss%dv_stress%data, &
                        jp_mid_lo, jp_mid_hi, p_ss_jet)
      call region_power(grid, ms, bt_work, &
                        ss%du_stress%data, ss%dv_stress%data, &
                        jp_mid_hi + 1, jp_hi, p_ss_pol)

      ! ---- Print ----
      if (header) then
         write (*, "(a)") "# BT-BUDGET PROBE"
         write (*, "(a)") "# columns: t stage region <KE>[m2/s2] eta_min eta_max "// &
            "|u_surf|max |u_bed|max |v_surf|max |v_bed|max "// &
            "P_pgf P_cor P_hv P_bd P_ss   [m2/s3]"
         write (*, "(a)") "# regions: SUB=subtropical (low j), JET=middle, POL=subpolar (high j)"
      end if

      call emit_row(t_value, t_unit, stage_tag, "SUB", ke_sub, &
                    eta_sub_min, eta_sub_max, &
                    usurf_sub, ubed_sub, vsurf_sub, vbed_sub, &
                    p_pgf_sub, p_cor_sub, p_hv_sub, p_bd_sub, p_ss_sub)
      call emit_row(t_value, t_unit, stage_tag, "JET", ke_jet, &
                    eta_jet_min, eta_jet_max, &
                    usurf_jet, ubed_jet, vsurf_jet, vbed_jet, &
                    p_pgf_jet, p_cor_jet, p_hv_jet, p_bd_jet, p_ss_jet)
      call emit_row(t_value, t_unit, stage_tag, "POL", ke_pol, &
                    eta_pol_min, eta_pol_max, &
                    usurf_pol, ubed_pol, vsurf_pol, vbed_pol, &
                    p_pgf_pol, p_cor_pol, p_hv_pol, p_bd_pol, p_ss_pol)
   end subroutine print_bt_budget

   subroutine region_eta_uv(grid, ms, bt_work, j_lo, j_hi, &
                            eta_min, eta_max, &
                            usurf_max, ubed_max, vsurf_max, vbed_max, ke_mean)
      !! Walk cells (interior x-range, j in [j_lo, j_hi]); collect
      !! η extremes, |u|/|v| extremes at surface + bed, and the
      !! cell-area-averaged BT KE.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(barotropic_workstate_t), intent(in) :: bt_work
      integer, intent(in) :: j_lo, j_hi
      real(wp), intent(out) :: eta_min, eta_max
      real(wp), intent(out) :: usurf_max, ubed_max, vsurf_max, vbed_max
      real(wp), intent(out) :: ke_mean

      integer :: i, j, ip_lo, ip_hi, nz_top, n
      real(wp) :: u_cell, v_cell, ke_sum

      ip_lo = grid%nghost + 1
      ip_hi = grid%nghost + grid%nx_phys
      nz_top = ms%nz_ml

      eta_min = huge(1.0_wp)
      eta_max = -huge(1.0_wp)
      usurf_max = 0.0_wp
      ubed_max = 0.0_wp
      vsurf_max = 0.0_wp
      vbed_max = 0.0_wp
      ke_sum = 0.0_wp
      n = 0
      do j = j_lo, j_hi
         do i = ip_lo, ip_hi
            eta_min = min(eta_min, bt_work%bt_eta(i, j))
            eta_max = max(eta_max, bt_work%bt_eta(i, j))
            ubed_max = max(ubed_max, abs(ms%u_face_x_layer(i, j, 1)), &
                           abs(ms%u_face_x_layer(i + 1, j, 1)))
            usurf_max = max(usurf_max, abs(ms%u_face_x_layer(i, j, nz_top)), &
                            abs(ms%u_face_x_layer(i + 1, j, nz_top)))
            vbed_max = max(vbed_max, abs(ms%v_face_y_layer(i, j, 1)), &
                           abs(ms%v_face_y_layer(i, j + 1, 1)))
            vsurf_max = max(vsurf_max, abs(ms%v_face_y_layer(i, j, nz_top)), &
                            abs(ms%v_face_y_layer(i, j + 1, nz_top)))
            u_cell = 0.5_wp*(bt_work%bt_ubt(i, j) + bt_work%bt_ubt(i + 1, j))
            v_cell = 0.5_wp*(bt_work%bt_vbt(i, j) + bt_work%bt_vbt(i, j + 1))
            ke_sum = ke_sum + 0.5_wp*(u_cell*u_cell + v_cell*v_cell)
            n = n + 1
         end do
      end do
      if (n > 0) then
         ke_mean = ke_sum/real(n, wp)
      else
         ke_mean = 0.0_wp
         eta_min = 0.0_wp
         eta_max = 0.0_wp
      end if
   end subroutine region_eta_uv

   subroutine region_power(grid, ms, bt_work, F_u_3d, F_v_3d, j_lo, j_hi, p_mean)
      !! Cell-centred BT power per region: P = ⟨u_bt·F_u + v_bt·F_v⟩.
      !! `F_u_3d` is per-layer at u-faces, shape `(nx+1, ny, nz)`.
      !! Depth-averages with the centred face thickness as weight to
      !! get the BT-mode contribution, then dots with `bt_ubt/bt_vbt`
      !! at the same face, then averages east+west (north+south) into
      !! the cell.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(barotropic_workstate_t), intent(in) :: bt_work
      real(wp), intent(in) :: F_u_3d(:, :, :), F_v_3d(:, :, :)
      integer, intent(in) :: j_lo, j_hi
      real(wp), intent(out) :: p_mean

      integer :: i, j, k, nz, ip_lo, ip_hi, n
      real(wp) :: h_face, fu_W, fu_E, fv_S, fv_N
      real(wp) :: w_W, w_E, w_S, w_N
      real(wp) :: power_sum, p_cell

      ip_lo = grid%nghost + 1
      ip_hi = grid%nghost + grid%nx_phys
      nz = ms%nz_ml

      power_sum = 0.0_wp
      n = 0
      do j = j_lo, j_hi
         do i = ip_lo, ip_hi
            ! Depth-mean each tendency at the four cell faces.  Weight
            ! by the centred face thickness so the BT mode is what we
            ! actually project onto.
            fu_W = 0.0_wp
            w_W = 0.0_wp
            fu_E = 0.0_wp
            w_E = 0.0_wp
            fv_S = 0.0_wp
            w_S = 0.0_wp
            fv_N = 0.0_wp
            w_N = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
               fu_W = fu_W + F_u_3d(i, j, k)*h_face
               w_W = w_W + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i + 1, j, k))
               fu_E = fu_E + F_u_3d(i + 1, j, k)*h_face
               w_E = w_E + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
               fv_S = fv_S + F_v_3d(i, j, k)*h_face
               w_S = w_S + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i, j + 1, k))
               fv_N = fv_N + F_v_3d(i, j + 1, k)*h_face
               w_N = w_N + h_face
            end do
            if (w_W > 0.0_wp) fu_W = fu_W/w_W
            if (w_E > 0.0_wp) fu_E = fu_E/w_E
            if (w_S > 0.0_wp) fv_S = fv_S/w_S
            if (w_N > 0.0_wp) fv_N = fv_N/w_N
            ! Power per face: u_bt·F.  Average into the cell.
            p_cell = 0.5_wp*(bt_work%bt_ubt(i, j)*fu_W + bt_work%bt_ubt(i + 1, j)*fu_E) + &
                     0.5_wp*(bt_work%bt_vbt(i, j)*fv_S + bt_work%bt_vbt(i, j + 1)*fv_N)
            power_sum = power_sum + p_cell
            n = n + 1
         end do
      end do
      if (n > 0) then
         p_mean = power_sum/real(n, wp)
      else
         p_mean = 0.0_wp
      end if
   end subroutine region_power

   subroutine region_power_drag_implicit(grid, ms, bt_work, rate_u, rate_v, &
                                         j_lo, j_hi, p_mean)
      !! Bottom drag in implicit mode: the slow tendency for u at a
      !! face is `−rate_u(i,j,k)·u_face_layer(i,j,k)` (1/s × m/s
      !! → m/s²).  We construct that on the fly, depth-mean by
      !! face thickness, dot with `bt_ubt` (and v counterpart),
      !! sum the cell-centred result over the region.  Sign comes
      !! out negative: drag removes BT-mode energy.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(barotropic_workstate_t), intent(in) :: bt_work
      real(wp), intent(in) :: rate_u(:, :, :), rate_v(:, :, :)
      integer, intent(in) :: j_lo, j_hi
      real(wp), intent(out) :: p_mean

      integer :: i, j, k, nz, ip_lo, ip_hi, n
      real(wp) :: h_face
      real(wp) :: fu_W, fu_E, fv_S, fv_N
      real(wp) :: w_W, w_E, w_S, w_N
      real(wp) :: power_sum, p_cell

      ip_lo = grid%nghost + 1
      ip_hi = grid%nghost + grid%nx_phys
      nz = ms%nz_ml

      power_sum = 0.0_wp
      n = 0
      do j = j_lo, j_hi
         do i = ip_lo, ip_hi
            fu_W = 0.0_wp
            w_W = 0.0_wp
            fu_E = 0.0_wp
            w_E = 0.0_wp
            fv_S = 0.0_wp
            w_S = 0.0_wp
            fv_N = 0.0_wp
            w_N = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
               fu_W = fu_W + (-rate_u(i, j, k)*ms%u_face_x_layer(i, j, k))*h_face
               w_W = w_W + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i + 1, j, k))
               fu_E = fu_E + (-rate_u(i + 1, j, k)*ms%u_face_x_layer(i + 1, j, k))*h_face
               w_E = w_E + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
               fv_S = fv_S + (-rate_v(i, j, k)*ms%v_face_y_layer(i, j, k))*h_face
               w_S = w_S + h_face
               h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i, j + 1, k))
               fv_N = fv_N + (-rate_v(i, j + 1, k)*ms%v_face_y_layer(i, j + 1, k))*h_face
               w_N = w_N + h_face
            end do
            if (w_W > 0.0_wp) fu_W = fu_W/w_W
            if (w_E > 0.0_wp) fu_E = fu_E/w_E
            if (w_S > 0.0_wp) fv_S = fv_S/w_S
            if (w_N > 0.0_wp) fv_N = fv_N/w_N
            p_cell = 0.5_wp*(bt_work%bt_ubt(i, j)*fu_W + bt_work%bt_ubt(i + 1, j)*fu_E) + &
                     0.5_wp*(bt_work%bt_vbt(i, j)*fv_S + bt_work%bt_vbt(i, j + 1)*fv_N)
            power_sum = power_sum + p_cell
            n = n + 1
         end do
      end do
      if (n > 0) then
         p_mean = power_sum/real(n, wp)
      else
         p_mean = 0.0_wp
      end if
   end subroutine region_power_drag_implicit

   subroutine emit_row(t_value, t_unit, stage_tag, region_tag, ke, &
                       eta_min, eta_max, usurf, ubed, vsurf, vbed, &
                       p_pgf, p_cor, p_hv, p_bd, p_ss)
      real(wp), intent(in) :: t_value, ke, eta_min, eta_max
      real(wp), intent(in) :: usurf, ubed, vsurf, vbed
      real(wp), intent(in) :: p_pgf, p_cor, p_hv, p_bd, p_ss
      character(len=*), intent(in) :: t_unit, stage_tag, region_tag
      write (*, "(a, f9.4, 1x, a, 1x, a, 1x, a, 13(1x, es12.4))") &
         "BT_BUDGET ", t_value, trim(t_unit), trim(stage_tag), trim(region_tag), &
         ke, eta_min, eta_max, usurf, ubed, vsurf, vbed, &
         p_pgf, p_cor, p_hv, p_bd, p_ss
   end subroutine emit_row

end module rdb_ocean_bt_budget_probe
