!! Largest stable barotropic substep (OM3-style `DTBT` definition).
module rdb_barotropic_dtbt
   !! The forward-backward barotropic substep (continuity on the BEBT-
   !! projected face velocity, then the momentum update with the new
   !! elevation and an explicit Coriolis term) is stable for
   !!
   !!     dt_bt^2 <= dt_max(i,j)^2  at every wet cell, with
   !!
   !!     1/dt_max^2 = 1/2 (1 + 2 bebt) [ (1/A_T) ( g D_E dy_E/dx_E + g D_W dy_W/dx_W
   !!                                             + g D_N dx_N/dy_N + g D_S dx_S/dy_S )
   !!                                     + f_NE^2 + f_NW^2 + f_SE^2 + f_SW^2 ].
   !!
   !! The bracketed gravity-wave sum is the squared highest wave frequency the
   !! cell's four faces can carry; the `1/2 (1 + 2 bebt)` factor is the
   !! stability limit of the projected scheme (Shchepetkin & McWilliams 2005;
   !! Hallberg & Adcroft 2009), and the corner `f^2` terms are the explicit
   !! Coriolis rotation.
   !!
   !! Mapping onto the BT arrays: `D_f dy_f` is the BT continuity face AREA the
   !! substep itself uses -- the face thickness `0.5 (D_a + D_b)` of the two
   !! adjacent column depths `D = sum_k h_layer` times `dy_cu_bt` / `dx_cv_bt`
   !! -- divided by the centre-to-centre distance `1/idxCu` / `1/idyCv`.
   !! Porous narrowing enters only through the `_bt` widths the caller passes:
   !! the PER-STEP re-derivation (`dtbt_reset_period`) sees the narrowing the
   !! porous refresh applied, but the SETUP-time value is formed before that
   !! refresh, from the un-narrowed widths -- the smaller, conservative
   !! `dt_max`.  A face is open only when BOTH adjacent
   !! cells are wet.  `g` is the substep gravity `bt_work%g_bt` (full gravity
   !! unless a reduced-gravity / SAL correction was configured).
   !!
   !! The step actually used is `|dtbt| * min_cells dt_max` for `dtbt < 0`
   !! (`dtbt = 0` means `-0.98`), or `dtbt` seconds for `dtbt > 0`;
   !! `n_inner = ceiling(dt / step)` and the substep is `dt / n_inner`.
   use rdb_constants, only: wp
   implicit none
   private

   public :: bt_dtbt_cell_inv2
   public :: bt_dtbt_inv2_terms
   public :: bt_dtbt_dt_max
   public :: bt_dtbt_safe_step
   public :: bt_dtbt_n_inner

contains

   pure function bt_dtbt_inv2_terms(cE, cW, cN, cS, iarea, f_sw, f_se, f_nw, f_ne, &
                                    g, bebt) result(inv2)
      !! `1/dt_max^2` of one cell from its four face conductances
      !! `c = D_face dy_face / dx_face` (m), `1/A_T` (1/m^2) and the four
      !! corner Coriolis parameters.
      !$acc routine seq
      real(wp), intent(in) :: cE, cW, cN, cS
         !! Face conductances `D dy/dx` (E, W) and `D dx/dy` (N, S) (m).
      real(wp), intent(in) :: iarea
         !! Inverse T-cell area (1/m^2).
      real(wp), intent(in) :: f_sw, f_se, f_nw, f_ne
         !! Corner Coriolis parameters (1/s).
      real(wp), intent(in) :: g
         !! Substep gravity (m/s^2).
      real(wp), intent(in) :: bebt
         !! `&ocean_bt_nml bebt`.
      real(wp) :: inv2
      inv2 = 0.5_wp*(1.0_wp + 2.0_wp*bebt)* &
             (iarea*g*(((cE + cW) + cN) + cS) + &
              (((f_ne*f_ne + f_nw*f_nw) + f_se*f_se) + f_sw*f_sw))
   end function bt_dtbt_inv2_terms

   pure function bt_dtbt_cell_inv2(i, j, nx, ny, nz, h, wet, dy_cu, dx_cv, idx_cu, idy_cv, &
                                   iarea_t, f_corner, g, bebt) result(inv2)
      !! `1/dt_max(i,j)^2` for T cell `(i,j)` from the layer thicknesses `h`
      !! (column depth `D = sum_k h`).  0 for a land cell.  Needs one ghost
      !! ring of `h` and `wet`.
      !$acc routine seq
      integer, intent(in) :: nx, ny, nz
      integer, intent(in) :: i, j
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idx_cu(nx + 1, ny)
      real(wp), intent(in) :: idy_cv(nx, ny + 1)
      real(wp), intent(in) :: iarea_t(nx, ny)
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1)
      real(wp), intent(in) :: g, bebt
      real(wp) :: inv2
      integer :: k
      real(wp) :: dc, de, dw, dn, ds, cE, cW, cN, cS

      inv2 = 0.0_wp
      if (wet(i, j) <= 0.5_wp) return
      dc = 0.0_wp
      de = 0.0_wp
      dw = 0.0_wp
      dn = 0.0_wp
      ds = 0.0_wp
      do k = 1, nz
         dc = dc + h(i, j, k)
         de = de + h(i + 1, j, k)
         dw = dw + h(i - 1, j, k)
         dn = dn + h(i, j + 1, k)
         ds = ds + h(i, j - 1, k)
      end do
      cE = 0.0_wp
      cW = 0.0_wp
      cN = 0.0_wp
      cS = 0.0_wp
      if (wet(i + 1, j) > 0.5_wp) cE = 0.5_wp*(dc + de)*dy_cu(i + 1, j)*idx_cu(i + 1, j)
      if (wet(i - 1, j) > 0.5_wp) cW = 0.5_wp*(dw + dc)*dy_cu(i, j)*idx_cu(i, j)
      if (wet(i, j + 1) > 0.5_wp) cN = 0.5_wp*(dc + dn)*dx_cv(i, j + 1)*idy_cv(i, j + 1)
      if (wet(i, j - 1) > 0.5_wp) cS = 0.5_wp*(ds + dc)*dx_cv(i, j)*idy_cv(i, j)
      inv2 = bt_dtbt_inv2_terms(max(cE, 0.0_wp), max(cW, 0.0_wp), max(cN, 0.0_wp), &
                                max(cS, 0.0_wp), iarea_t(i, j), &
                                f_corner(i, j), f_corner(i + 1, j), &
                                f_corner(i, j + 1), f_corner(i + 1, j + 1), g, bebt)
   end function bt_dtbt_cell_inv2

   subroutine bt_dtbt_dt_max(nx, ny, nz, i0, i1, j0, j1, h, wet, dy_cu, dx_cv, idx_cu, &
                             idy_cv, iarea_t, f_corner, g, bebt, on_device, dt_max)
      !! Domain-local `dt_max = min over wet cells in [i0:i1, j0:j1]` of the
      !! cell limit; `huge` if no wet cell has a nonzero limit (the identity of
      !! the cross-rank `min`).  `on_device = .true.` runs the `do concurrent`
      !! reduction (every array device-present); `.false.` is a plain host loop
      !! for configure time, before the state is mapped.  The device branch
      !! spells the cell body out inline: the same arithmetic reached through
      !! the nested `routine seq` call `bt_dtbt_cell_inv2` returned a wrong
      !! (too large) `dt_max` on nvfortran `-stdpar=gpu` (155 s against the
      !! host 74.6 s on the 1-degree grid) while every input array matched the
      !! host copy; `test_ocean_bt_cfl` pins the device branch against the host
      !! one.
      integer, intent(in) :: nx, ny, nz, i0, i1, j0, j1
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: dy_cu(nx + 1, ny)
      real(wp), intent(in) :: dx_cv(nx, ny + 1)
      real(wp), intent(in) :: idx_cu(nx + 1, ny)
      real(wp), intent(in) :: idy_cv(nx, ny + 1)
      real(wp), intent(in) :: iarea_t(nx, ny)
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1)
      real(wp), intent(in) :: g, bebt
      logical, intent(in) :: on_device
      real(wp), intent(out) :: dt_max
      integer :: i, j, k
      real(wp) :: inv2_max, v, dc, de, dw, dn, ds, cE, cW, cN, cS

      inv2_max = 0.0_wp
      if (on_device) then
         do concurrent(j=j0:j1, i=i0:i1) local(k, v, dc, de, dw, dn, ds, cE, cW, cN, cS) &
            reduce(max:inv2_max)
            ! Cell body inlined on purpose: calling the `!$acc routine seq`
            ! helper `bt_dtbt_cell_inv2` from inside this device reduction
            ! gave a wrong (too large) `dt_max` on nvfortran -stdpar=gpu.
            ! The host branch below still uses the helper.
            v = 0.0_wp
            if (wet(i, j) > 0.5_wp) then
               dc = 0.0_wp
               de = 0.0_wp
               dw = 0.0_wp
               dn = 0.0_wp
               ds = 0.0_wp
               do k = 1, nz
                  dc = dc + h(i, j, k)
                  de = de + h(i + 1, j, k)
                  dw = dw + h(i - 1, j, k)
                  dn = dn + h(i, j + 1, k)
                  ds = ds + h(i, j - 1, k)
               end do
               cE = merge(0.5_wp*(dc + de)*dy_cu(i + 1, j)*idx_cu(i + 1, j), 0.0_wp, wet(i + 1, j) > 0.5_wp)
               cW = merge(0.5_wp*(dw + dc)*dy_cu(i, j)*idx_cu(i, j), 0.0_wp, wet(i - 1, j) > 0.5_wp)
               cN = merge(0.5_wp*(dc + dn)*dx_cv(i, j + 1)*idy_cv(i, j + 1), 0.0_wp, wet(i, j + 1) > 0.5_wp)
               cS = merge(0.5_wp*(ds + dc)*dx_cv(i, j)*idy_cv(i, j), 0.0_wp, wet(i, j - 1) > 0.5_wp)
               v = 0.5_wp*(1.0_wp + 2.0_wp*bebt)* &
                   (iarea_t(i, j)*g*(((max(cE, 0.0_wp) + max(cW, 0.0_wp)) + max(cN, 0.0_wp)) + max(cS, 0.0_wp)) + &
                    (((f_corner(i + 1, j + 1)**2 + f_corner(i, j + 1)**2) + f_corner(i + 1, j)**2) + f_corner(i, j)**2))
            end if
            inv2_max = max(inv2_max, v)
         end do
      else
         do j = j0, j1
            do i = i0, i1
               v = bt_dtbt_cell_inv2(i, j, nx, ny, nz, h, wet, dy_cu, dx_cv, idx_cu, idy_cv, &
                                     iarea_t, f_corner, g, bebt)
               inv2_max = max(inv2_max, v)
            end do
         end do
      end if
      if (inv2_max > 0.0_wp) then
         dt_max = 1.0_wp/sqrt(inv2_max)
      else
         dt_max = huge(1.0_wp)
      end if
   end subroutine bt_dtbt_dt_max

   pure function bt_dtbt_safe_step(dtbt, dt_max) result(step)
      !! Largest substep to use: `|dtbt| * dt_max` for `dtbt < 0`
      !! (`dtbt = 0` means `-0.98`), `dtbt` itself for `dtbt > 0`.
      real(wp), intent(in) :: dtbt
         !! `&ocean_bt_nml dtbt`.
      real(wp), intent(in) :: dt_max
         !! Global minimum stable substep (s); `huge` if no wet cell.
      real(wp) :: step
      if (dtbt > 0.0_wp) then
         step = dtbt
      else if (dtbt == 0.0_wp) then
         step = 0.98_wp*dt_max
      else
         step = (-dtbt)*dt_max
      end if
   end function bt_dtbt_safe_step

   pure function bt_dtbt_n_inner(dt, dtbt, dt_max) result(n_inner)
      !! `n_inner = max(1, ceiling(dt / step))`.
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: dtbt
      real(wp), intent(in) :: dt_max
      integer :: n_inner
      real(wp) :: step
      step = bt_dtbt_safe_step(dtbt, min(dt_max, 1.0e300_wp))
      if (dtbt <= 0.0_wp .and. dt_max >= huge(1.0_wp)) then
         n_inner = 1      ! no wet cell anywhere
      else
         n_inner = max(1, ceiling(dt/step))
      end if
   end function bt_dtbt_n_inner

end module rdb_barotropic_dtbt
