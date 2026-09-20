!! Barotropic ↔ multilayer coupling kernels for the split-explicit ocean
!! driver — the bt↔ml bridges run around each barotropic-substep call:
!! `derive_bt_from_layers` (layers → bt state, depth-averaged),
!! `sum_slow_tendencies_into_F_slow` (per-kernel scratch → composite F_slow),
!! `face_depth_mean_u/_v` (3D slow tendency → 2D face forcing),
!! `apply_bt_correction` (bt time-mean → per-layer correction + h rescale).
module rdb_barotropic_coupling
   use rdb_constants, only: wp, GRAVITY
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_MOM6
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: OBC_PERIODIC
   implicit none
   private

   public :: derive_bt_from_layers
   public :: compute_h_face_upstream
   public :: sum_slow_tendencies_into_F_slow
   public :: subtract_fast_cor_ref
   public :: set_cor_ref_velocity
   public :: face_depth_mean_u
   public :: face_depth_mean_v
   public :: face_depth_mean_rem_u
   public :: face_depth_mean_rem_v
   public :: apply_bt_correction
   public :: snapshot_eta_PF
   public :: compute_pbce
   public :: compute_gtot_faces
   public :: compute_e_anom
   public :: compute_bt_rem
   public :: reset_bt_rem
   public :: compute_bt_rem_wave_drag
   public :: mask_bt_rem
   ! BT_cont_type producer — fills BTCL_u/v on bt_work
   public :: set_local_BT_cont_types

   ! Floor below which a layer counts as "vanished" for FA accumulation.
   real(wp), parameter :: BTC_H_NEGLECT = 1.0e-10_wp
   ! Volume-CFL safety factor (vol_CFL = 0.5). NOTE: too loose for the
   ! upstream-h-sum producer — the closure stays in the cubic-near-zero
   ! branch (≈ naive u·h_face); real saturation needs a PPM-perturbation FA.
   real(wp), parameter :: BTC_VOL_CFL = 0.5_wp

contains

   pure subroutine derive_bt_from_layers(grid, bt_work, ms)
      !! Populate `bt_eta`, `bt_ubt`, `bt_vbt` from the current multilayer
      !! state. `bt_H_ref` must already be set.
      !!   bt_eta = Σ_k h_layer − H_ref;  bt_ubt = Σ_k(u·h_face)/Σ_k h_face.
      !! Face thickness averages the two abutting columns (wall faces use the
      !! single cell). With `use_upstream_h_face`, the interior face-h is the
      !! first-order upwind pick — consistent with `compute_h_face_upstream`.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz, nx_face, ny_face
      logical :: use_upstream
      real(wp) :: total_h, hu_sum, h_face_sum, h_face, hv_sum

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      use_upstream = bt_work%use_upstream_h_face

      do concurrent(j=1:ny, i=1:nx) local(k, total_h)
         total_h = 0.0_wp
         do k = 1, nz
            total_h = total_h + ms%h_layer(i, j, k)
         end do
         bt_work%bt_eta(i, j) = total_h - bt_work%bt_H_ref(i, j)
      end do

      do concurrent(j=1:ny, i=1:nx_face) &
         local(k, hu_sum, h_face_sum, h_face)
         hu_sum = 0.0_wp
         h_face_sum = 0.0_wp
         do k = 1, nz
            if (i == 1) then
               h_face = ms%h_layer(1, j, k)
            else if (i == nx_face) then
               h_face = ms%h_layer(nx, j, k)
            else if (use_upstream) then
               if (ms%u_face_x_layer(i, j, k) >= 0.0_wp) then
                  h_face = ms%h_layer(i - 1, j, k)
               else
                  h_face = ms%h_layer(i, j, k)
               end if
            else
               h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
            end if
            hu_sum = hu_sum + ms%u_face_x_layer(i, j, k)*h_face
            h_face_sum = h_face_sum + h_face
         end do
         if (h_face_sum > 0.0_wp) then
            bt_work%bt_ubt(i, j) = hu_sum/h_face_sum
         else
            bt_work%bt_ubt(i, j) = 0.0_wp
         end if
      end do

      do concurrent(j=1:ny_face, i=1:nx) &
         local(k, hv_sum, h_face_sum, h_face)
         hv_sum = 0.0_wp
         h_face_sum = 0.0_wp
         do k = 1, nz
            if (j == 1) then
               h_face = ms%h_layer(i, 1, k)
            else if (j == ny_face) then
               h_face = ms%h_layer(i, ny, k)
            else if (use_upstream) then
               if (ms%v_face_y_layer(i, j, k) >= 0.0_wp) then
                  h_face = ms%h_layer(i, j - 1, k)
               else
                  h_face = ms%h_layer(i, j, k)
               end if
            else
               h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
            end if
            hv_sum = hv_sum + ms%v_face_y_layer(i, j, k)*h_face
            h_face_sum = h_face_sum + h_face
         end do
         if (h_face_sum > 0.0_wp) then
            bt_work%bt_vbt(i, j) = hv_sum/h_face_sum
         else
            bt_work%bt_vbt(i, j) = 0.0_wp
         end if
      end do
   end subroutine derive_bt_from_layers

   pure subroutine compute_h_face_upstream(grid, bt_work, ms)
      !! Per-face upstream column-sum thickness `h_face_up_x/y(I,j) =
      !! Σ_k h_layer(I_upstream,j,k)` used by the BT chain when
      !! `use_upstream_h_face = .true.`. First-order upwind pick by face-velocity
      !! sign (sampled at the top of the outer step). Wall faces use the single
      !! available cell. No-op when the knob is off.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz, nx_face, ny_face
      real(wp) :: h_sum, h_k

      if (.not. bt_work%use_upstream_h_face) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_face = size(bt_work%h_face_up_x, 1)
      ny_face = size(bt_work%h_face_up_y, 2)

      ! East-face: upstream pick from u_face_x_layer sign.
      do concurrent(j=1:ny, i=1:nx_face) local(k, h_sum, h_k)
         h_sum = 0.0_wp
         do k = 1, nz
            if (i == 1) then
               h_k = ms%h_layer(1, j, k)
            else if (i == nx_face) then
               h_k = ms%h_layer(nx, j, k)
            else if (ms%u_face_x_layer(i, j, k) >= 0.0_wp) then
               h_k = ms%h_layer(i - 1, j, k)
            else
               h_k = ms%h_layer(i, j, k)
            end if
            h_sum = h_sum + h_k
         end do
         bt_work%h_face_up_x(i, j) = h_sum
      end do

      ! North-face: upstream pick from v_face_y_layer sign.
      do concurrent(j=1:ny_face, i=1:nx) local(k, h_sum, h_k)
         h_sum = 0.0_wp
         do k = 1, nz
            if (j == 1) then
               h_k = ms%h_layer(i, 1, k)
            else if (j == ny_face) then
               h_k = ms%h_layer(i, ny, k)
            else if (ms%v_face_y_layer(i, j, k) >= 0.0_wp) then
               h_k = ms%h_layer(i, j - 1, k)
            else
               h_k = ms%h_layer(i, j, k)
            end if
            h_sum = h_sum + h_k
         end do
         bt_work%h_face_up_y(i, j) = h_sum
      end do
   end subroutine compute_h_face_upstream

   pure subroutine sum_slow_tendencies_into_F_slow(bt_work, pgf, cor, hv, bd, ss, ms)
      !! Sum the per-kernel slow-tendency scratch buffers into a
      !! single (`bt_work%F_slow_u`, `bt_work%F_slow_v`) field per face per
      !! layer.  Reads `pgf%dpdx_face`, `cor%pv_flux_x`, `hv%du_visc`,
      !! `bd%du_drag`, `ss%du_stress` and their v counterparts —
      !! all already at the matching u-face / v-face shape.  Each
      !! kernel must have run its compute step before this is called.
      !!
      !! ## What this list IS (and what it is not)
      !!
      !! The five terms are the ADDITIVE layer-tendency buffers that
      !! `run_stage_split` applies between the forcing assembly and
      !! `apply_bt_correction`.  The depth mean of this sum becomes
      !! `F_bt_u/v`, the frozen forcing the barotropic substep integrates,
      !! and the SAME `F_bt` is subtracted again inside the correction.
      !!
      !! It is tempting to read the list as a completeness requirement —
      !! "a tendency missing from here never reaches the barotropic mode".
      !! It is not, because `apply_bt_correction` adds an INCREMENT rather
      !! than REPLACING the layer depth mean.  Writing `T_k` for the
      !! summed tendencies and `D_k` for one that is applied but omitted,
      !! a stage does
      !!
      !!     u_k^{n+1} = u_k^n + dt·(T_k + D_k)
      !!                       + (u_bt^end − u_bt^n − dt·F_bt),
      !!
      !! and with `F_bt = ⟨T⟩_h` plus `⟨u^n⟩_h = u_bt^n`
      !! (`derive_bt_from_layers`) its thickness-weighted depth mean is
      !!
      !!     ⟨u^{n+1}⟩ = u_bt^end + dt·⟨D⟩.
      !!
      !! So `⟨D⟩` is applied EXACTLY ONCE, on top of the barotropic
      !! solution — never lost, never double counted.  The mirror holds
      !! for a summed term: the fast loop integrates `⟨T⟩` across the
      !! substeps and the `−dt·F_bt` guard takes it straight back out, so
      !! membership is depth-mean NEUTRAL to leading order.
      !!
      !! What membership actually buys is SECOND order, and it is real:
      !! a summed term shapes the substep's live η / ζ / KE / `bt_rem`
      !! trajectory and therefore the time-mean transports `bt_uhbt` that
      !! the slow continuity renormalises to.  Omitting a term is a
      !! first-order-in-dt operator SPLIT — the depth mean is applied
      !! after the fast loop instead of inside it.
      !!
      !! ## The contract for a new tendency
      !!
      !! A new layer velocity tendency applied inside the corrected set
      !! MUST come with a decision about this sum, recorded here:
      !!
      !! * **Summed (the default).**  An additive `du/dt` buffer that the
      !!   barotropic mode should feel DURING the substeps — anything that
      !!   changes the depth-mean momentum on the fast timescale, or whose
      !!   transport must be consistent with `bt_uhbt`.  Add it to both
      !!   `do concurrent` bodies below; keep them branch-free and
      !!   explicit-shape.
      !! * **Deliberately omitted.**  A term with no additive tendency
      !!   buffer to sum (a multiplicative / implicit operator), or one
      !!   whose depth-mean lag is physically irrelevant.  The side-wall
      !!   CHANNEL drag (`ocean_channel_drag_apply_tendencies`) is the one
      !!   such term today: it is a frozen-rate BACKWARD-EULER factor
      !!   `u ← u/(1 + dt·λ_side)`, not a `du/dt` field, so there is no
      !!   buffer to add — and by the algebra above its depth mean still
      !!   lands exactly once.
      !!
      !! Tendencies applied by separate operator-split steps AFTER
      !! `apply_bt_correction` (the baroclinic OBC, the sponges, the
      !! implicit vertical friction `vdiff_apply_momentum` and the
      !! `&ocean_vdiff_nml implicit_stress` / `implicit_drag` folds) are
      !! outside the corrected set and must NOT be summed here.  Note the
      !! deliberate asymmetry that follows: `bd%du_drag` / `ss%du_stress`
      !! stay in this sum even when their explicit applies are skipped by
      !! those folds — the fast loop adds their depth mean and the
      !! correction removes it, so the fold's later application is still
      !! the only one.
      !!
      !! The gate is `tests/test_ocean_bt_slow_forcing.F90`
      !! (`channel_drag_depth_mean_decay` for an omitted term,
      !! `bottom_drag_depth_mean_decay` for a summed one): on a
      !! doubly-periodic uniform box each decays at the outer scheme's
      !! exact analytic rate under both `pred_corr` and `ssp_rk2`.
      !!
      !! Terms that reach the barotropic mode by a DIFFERENT seam have no
      !! business here either: `&ocean_bt_nml substep_drag` and the linear
      !! wave drag multiply `bt_rem_u/v` inside the substep, the porous
      !! barriers narrow the substep transport widths, and the tide / SAL
      !! / surface-pressure loads arrive as `eta_forcing`.
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(ocean_pressure_force_t), intent(in) :: pgf
      type(coriolis_adv_t), intent(in) :: cor
      type(ocean_horizontal_viscosity_t), intent(in) :: hv
      type(ocean_bottom_drag_t), intent(in) :: bd
      type(ocean_surface_stress_t), intent(in) :: ss
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nu, nv, nx, ny, nz

      nu = size(bt_work%F_slow_u, 1)
      nv = size(bt_work%F_slow_v, 2)
      nx = size(bt_work%F_slow_v, 1)
      ny = size(bt_work%F_slow_u, 2)
      nz = ms%nz_ml

      do concurrent(k=1:nz, j=1:ny, i=1:nu)
         bt_work%F_slow_u(i, j, k) = pgf%dpdx_face%data(i, j, k) + &
                                     cor%pv_flux_x%data(i, j, k) + &
                                     hv%du_visc%data(i, j, k) + &
                                     bd%du_drag%data(i, j, k) + &
                                     ss%du_stress%data(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:nv, i=1:nx)
         bt_work%F_slow_v(i, j, k) = pgf%dpdy_face%data(i, j, k) + &
                                     cor%pv_flux_y%data(i, j, k) + &
                                     hv%dv_visc%data(i, j, k) + &
                                     bd%dv_drag%data(i, j, k) + &
                                     ss%dv_stress%data(i, j, k)
      end do
   end subroutine sum_slow_tendencies_into_F_slow

   pure subroutine subtract_fast_cor_ref(grid, metrics, bt_work, f_corner, &
                                         bc_w, bc_e, bc_s, bc_n, &
                                         has_w, has_e, has_s, has_n)
      !! Subtract the fast-loop Coriolis + vector-invariant advection,
      !! evaluated at the reference barotropic velocity
      !! `bt_work%cor_ref_u`/`cor_ref_v` (filled by
      !! `set_cor_ref_velocity`), from the substep forcing
      !! `F_bt_u_fast`/`F_bt_v_fast`.
      !!
      !! The reference velocity is NOT free: it must be the depth mean
      !! of the same layer velocity the slow `cor%pv_flux_*` inside
      !! `F_bt_u/v` was evaluated on, or the difference survives as a
      !! near-constant per-substep forcing.  Under `ssp_rk2` that is the
      !! stage-entry `bt_ubt/bt_vbt` (what this routine used to read
      !! directly); under `pred_corr` it is the depth mean of `u_av/v_av`.
      !! See `set_cor_ref_velocity` and the `cor_ref_u` docstring in
      !! `rdb_barotropic_workstate`.
      !!
      !! Why: `F_bt_u` is the depth mean of ALL slow layer tendencies —
      !! including the layer Coriolis-advection (`cor%pv_flux_*`), whose
      !! depth mean is ≈ the fast solver's own `(ζ+f)·v − ∇KE` at the
      !! stage-entry state.  The substep then integrates its own LIVE
      !! `(ζ+f)·v − ∇KE` on top, so without this subtraction the
      !! barotropic Coriolis/advection is integrated TWICE — the exact
      !! analogue of the PGF double-count the `F_bt_u_fast = F_bt_u −
      !! depth_mean(PGF)` subtraction already guards against ("√(gH)
      !! inflates to √(2gH)").  The Coriolis double-count is what pumps
      !! the exponential wall/corner barotropic mode on shelf rims under
      !! `VCOORD_LAGRANGIAN` (600² double-gyre h-guard trap).  MOM6
      !! removes it with a reference Coriolis/advection (`Cor_ref_u/v`)
      !! subtracted inside its barotropic substep loop; this is that
      !! subtraction, folded into the forcing so the substep kernel is
      !! untouched.
      !!
      !! After this, at τ=0 the substep's net Coriolis/advection
      !! contribution is zero and only the ANOMALY that develops over
      !! the substeps is integrated — so `Δu = ubt_end − ubt_at_n −
      !! dt·F_bt_u` hands the layers the genuine fast anomaly instead
      !! of an extra `dt·f·v̄` rotation per stage.
      !!
      !! Uses `bt_work%bt_zeta_corner` and `bt_work%bt_ke_centre` as
      !! scratch — both are per-substep scratch the substep recomputes
      !! from scratch before reading (Pass 1 precedes Pass 2b/2c).
      !! ζ/KE formulas and wall closures mirror the substep exactly
      !! (`rdb_barotropic_substep` Pass 1 + the corner-ζ closure), so
      !! the τ=0 cancellation holds at walls too — where the unstable
      !! mode lives.  bt_halo = 0 index convention (the wide-halo
      !! march-in path receives the same interior-computed forcing via
      !! `copy_in`; its seam rows differ at O(ghost) — acceptable, the
      !! halo exchange owns them).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      real(wp), intent(in) :: f_corner(grid%nx_total + 1, grid%ny_total + 1)
      integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
         !! Per-edge OBC tags (OBC_WALL when no bc present).
      logical, intent(in) :: has_w, has_e, has_s, has_n
         !! Physical-edge flags (.false. at an MPI seam).

      integer :: i, j, nx, ny
      real(wp) :: zeta_at_u, f_at_u, v_at_u, ke_grad_x
      real(wp) :: zeta_at_v, f_at_v, u_at_v, ke_grad_y
      real(wp) :: w_nl
         !! Mirror of the substep's live-nonlinear weight
         !! (`bt_work%substep_zeta_ke`): the reference MUST subtract
         !! exactly what the substep re-integrates live — full
         !! `(ζ+f)·v − ∇KE` when live (1), planetary `f·v̄` only when
         !! the substep is planetary-only (0).

      nx = grid%nx_total
      ny = grid%ny_total
      w_nl = merge(1.0_wp, 0.0_wp, bt_work%substep_zeta_ke)

      ! ---- ζ at corners from the stage-entry bt state (substep Pass-1
      ! formula + closure).  Corner ring + physical-wall lines → 0.
      do concurrent(j=1:ny + 1, i=1:nx + 1)
         if (i >= 2 .and. i <= nx .and. j >= 2 .and. j <= ny) then
            bt_work%bt_zeta_corner(i, j) = &
               ((bt_work%cor_ref_v(i, j)*metrics%dyCv(i, j) &
                 - bt_work%cor_ref_v(i - 1, j)*metrics%dyCv(i - 1, j)) - &
                (bt_work%cor_ref_u(i, j)*metrics%dxCu(i, j) &
                 - bt_work%cor_ref_u(i, j - 1)*metrics%dxCu(i, j - 1)))* &
               metrics%iareaBu(i, j)
         else
            bt_work%bt_zeta_corner(i, j) = 0.0_wp
         end if
         if (bc_w /= OBC_PERIODIC .and. has_w .and. i == grid%nghost + 1) then
            bt_work%bt_zeta_corner(i, j) = 0.0_wp
         end if
         if (bc_e /= OBC_PERIODIC .and. has_e .and. i == grid%nghost + grid%nx_phys + 1) then
            bt_work%bt_zeta_corner(i, j) = 0.0_wp
         end if
         if (bc_s /= OBC_PERIODIC .and. has_s .and. j == grid%nghost + 1) then
            bt_work%bt_zeta_corner(i, j) = 0.0_wp
         end if
         if (bc_n /= OBC_PERIODIC .and. has_n .and. j == grid%nghost + grid%ny_phys + 1) then
            bt_work%bt_zeta_corner(i, j) = 0.0_wp
         end if
      end do

      ! ---- KE at centres (substep Pass-1 formula). ----
      do concurrent(j=1:ny, i=1:nx)
         bt_work%bt_ke_centre(i, j) = 0.25_wp*metrics%iareaT(i, j)*( &
                                      metrics%areaCu(i, j)*bt_work%cor_ref_u(i, j)**2 + &
                                      metrics%areaCu(i + 1, j)*bt_work%cor_ref_u(i + 1, j)**2 + &
                                      metrics%areaCv(i, j)*bt_work%cor_ref_v(i, j)**2 + &
                                      metrics%areaCv(i, j + 1)*bt_work%cor_ref_v(i, j + 1)**2)
      end do

      ! ---- Subtract the u-face reference (substep Pass-2b operand). ----
      do concurrent(j=1:ny, i=2:nx) local(zeta_at_u, f_at_u, v_at_u, ke_grad_x)
         zeta_at_u = w_nl*0.5_wp*(bt_work%bt_zeta_corner(i, j) + bt_work%bt_zeta_corner(i, j + 1))
         f_at_u = 0.5_wp*(f_corner(i, j) + f_corner(i, j + 1))
         if (j > 1 .and. j < ny) then
            v_at_u = 0.25_wp*(bt_work%cor_ref_v(i - 1, j) + bt_work%cor_ref_v(i - 1, j + 1) + &
                              bt_work%cor_ref_v(i, j) + bt_work%cor_ref_v(i, j + 1))
         else if (j == 1) then
            v_at_u = 0.5_wp*(bt_work%cor_ref_v(i - 1, j + 1) + bt_work%cor_ref_v(i, j + 1))
         else
            v_at_u = 0.5_wp*(bt_work%cor_ref_v(i - 1, j) + bt_work%cor_ref_v(i, j))
         end if
         ke_grad_x = w_nl*(bt_work%bt_ke_centre(i, j) - bt_work%bt_ke_centre(i - 1, j))*metrics%idxCu(i, j)
         bt_work%F_bt_u_fast(i, j) = bt_work%F_bt_u_fast(i, j) - &
                                     ((zeta_at_u + f_at_u)*v_at_u - ke_grad_x)
      end do

      ! ---- Subtract the v-face reference (substep Pass-2c operand). ----
      do concurrent(j=2:ny, i=1:nx) local(zeta_at_v, f_at_v, u_at_v, ke_grad_y)
         zeta_at_v = w_nl*0.5_wp*(bt_work%bt_zeta_corner(i, j) + bt_work%bt_zeta_corner(i + 1, j))
         f_at_v = 0.5_wp*(f_corner(i, j) + f_corner(i + 1, j))
         if (i > 1 .and. i < nx) then
            u_at_v = 0.25_wp*(bt_work%cor_ref_u(i, j - 1) + bt_work%cor_ref_u(i + 1, j - 1) + &
                              bt_work%cor_ref_u(i, j) + bt_work%cor_ref_u(i + 1, j))
         else if (i == 1) then
            u_at_v = 0.5_wp*(bt_work%cor_ref_u(i + 1, j - 1) + bt_work%cor_ref_u(i + 1, j))
         else
            u_at_v = 0.5_wp*(bt_work%cor_ref_u(i, j - 1) + bt_work%cor_ref_u(i, j))
         end if
         ke_grad_y = w_nl*(bt_work%bt_ke_centre(i, j) - bt_work%bt_ke_centre(i, j - 1))*metrics%idyCv(i, j)
         bt_work%F_bt_v_fast(i, j) = bt_work%F_bt_v_fast(i, j) - &
                                     (-(zeta_at_v + f_at_v)*u_at_v - ke_grad_y)
      end do
   end subroutine subtract_fast_cor_ref

   pure subroutine set_cor_ref_velocity(grid, bt_work, ms, from_u_av)
      !! Fill `bt_work%cor_ref_u/v` — the barotropic velocity at which
      !! `subtract_fast_cor_ref` evaluates the Coriolis/advection
      !! reference it removes from the substep forcing (MOM6
      !! `ubt_Cor`/`vbt_Cor`).
      !!
      !! The reference MUST be the depth mean of the same layer
      !! velocity whose Coriolis-advection tendency (`cor%pv_flux_*`)
      !! was depth-averaged into `F_bt_u/v`, under the same weights.
      !! Otherwise the two do not cancel at τ=0 and the residual
      !! `f × (v̄_ref − v̄_slow)` enters EVERY barotropic substep as a
      !! near-constant forcing.  In a closed rotating basin that
      !! residual projects onto the gravest Poincaré seiche and pumps
      !! it exponentially (e-folding ~0.6 d on a 240 km f-plane square
      !! at dt = 300 s; growth rate ∝ dt and rising with `n_inner` —
      !! the fingerprint of a fixed per-substep forcing, not an inner
      !! loop instability).
      !!
      !! * `from_u_av = .false.` (`ssp_rk2`) — the slow tendencies were
      !!   evaluated on the prognostic `u^n`, whose depth mean is the
      !!   stage-entry `bt_ubt/bt_vbt` from `derive_bt_from_layers`.
      !!   A plain copy, so the arithmetic downstream is bit-identical
      !!   to reading `bt_ubt/bt_vbt` directly.
      !! * `from_u_av = .true.` (`pred_corr`) — the slow tendencies were
      !!   evaluated on the time-mean `u_av/v_av` (`run_stage_split`
      !!   step 2), so take ITS depth mean, weighted exactly as the
      !!   forcing depth-mean was (h, or h·visc_rem when `&ocean_bt_nml
      !!   forcing_visc_rem`).  MOM6 does the same by construction:
      !!   `ubt_Cor = Σ_k wt_u·U_Cor` with `U_Cor = u_av`, the velocity
      !!   its `CorAdCalc` used.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      logical, intent(in) :: from_u_av
         !! `.true.` under `split_scheme = "pred_corr"`.

      integer :: i, j, nu, nv, nx, ny
      logical :: use_av

      use_av = from_u_av
      if (use_av) use_av = allocated(ms%u_av_layer) .and. allocated(ms%v_av_layer)

      if (use_av) then
         if (bt_work%bt_forcing_visc_rem) then
            call face_depth_mean_rem_u(grid, ms%u_av_layer, ms%h_layer, &
                                       bt_work%visc_rem_u, bt_work%cor_ref_u, ms%nz_ml)
            call face_depth_mean_rem_v(grid, ms%v_av_layer, ms%h_layer, &
                                       bt_work%visc_rem_v, bt_work%cor_ref_v, ms%nz_ml)
         else
            call face_depth_mean_u(grid, ms%u_av_layer, ms%h_layer, bt_work%cor_ref_u, ms%nz_ml)
            call face_depth_mean_v(grid, ms%v_av_layer, ms%h_layer, bt_work%cor_ref_v, ms%nz_ml)
         end if
      else
         nu = size(bt_work%bt_ubt, 1)
         ny = size(bt_work%bt_ubt, 2)
         nx = size(bt_work%bt_vbt, 1)
         nv = size(bt_work%bt_vbt, 2)
         do concurrent(j=1:ny, i=1:nu)
            bt_work%cor_ref_u(i, j) = bt_work%bt_ubt(i, j)
         end do
         do concurrent(j=1:nv, i=1:nx)
            bt_work%cor_ref_v(i, j) = bt_work%bt_vbt(i, j)
         end do
      end if
   end subroutine set_cor_ref_velocity

   pure subroutine face_depth_mean_u(grid, F_3d, h_layer, F_mean_2d, nz)
      !! Depth-average a u-face 3D field, weighted by the face
      !! thickness (= mean of the two abutting cell columns'
      !! `h_layer` values).  Writes to a 2D field at the same u-face
      !! shape.  Wall faces (i=1, nx+1) fall back to the single
      !! available cell.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx+1,ny,nz) / (nx,ny+1,nz);
      ! a single (nx,ny,nz) explicit-shape triplet would mis-bound the face axis.
      ! size() is used to derive loop bounds from the actual face dimension.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      integer :: i, j, k, nu, ny, nx_cells
      real(wp) :: h_face, num, denom

      nu = size(F_3d, 1)
      ny = size(F_3d, 2)
      nx_cells = grid%nx_total

      do concurrent(j=1:ny, i=1:nu) local(k, h_face, num, denom)
         num = 0.0_wp
         denom = 0.0_wp
         do k = 1, nz
            if (i == 1) then
               h_face = h_layer(1, j, k)
            else if (i == nu) then
               h_face = h_layer(nx_cells, j, k)
            else
               h_face = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            end if
            num = num + F_3d(i, j, k)*h_face
            denom = denom + h_face
         end do
         if (denom > 0.0_wp) then
            F_mean_2d(i, j) = num/denom
         else
            F_mean_2d(i, j) = 0.0_wp
         end if
      end do
   end subroutine face_depth_mean_u

   pure subroutine face_depth_mean_rem_u(grid, F_3d, h_layer, rem, F_mean_2d, nz)
      !! `face_depth_mean_u` with MOM6 `wt_u` weighting (`&ocean_bt_nml
      !! forcing_visc_rem`): the weight is
      !! `h_face·visc_rem(k)` instead of `h_face`, so layers the implicit
      !! vertical-friction solve will immediately damp (grounded sliver
      !! stacks under the BBL glue, visc_rem → 0) contribute nothing to
      !! the barotropic forcing.  Without this the spurious grounded-layer
      !! PGF's depth-mean drives the fast loop ballistically even after
      !! the layer velocities themselves are glued (PGF_BUG.md §9).
      !! Denominator falls back to zero-output on an all-remnant-zero
      !! column (the substep should not force an immobilized column).
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx+1,ny,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(in) :: rem(:, :, :)      ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      integer :: i, j, k, nu, ny, nx_cells
      real(wp) :: h_face, wt, num, denom

      nu = size(F_3d, 1)
      ny = size(F_3d, 2)
      nx_cells = grid%nx_total

      do concurrent(j=1:ny, i=1:nu) local(k, h_face, wt, num, denom)
         num = 0.0_wp
         denom = 0.0_wp
         do k = 1, nz
            if (i == 1) then
               h_face = h_layer(1, j, k)
            else if (i == nu) then
               h_face = h_layer(nx_cells, j, k)
            else
               h_face = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
            end if
            wt = h_face*min(max(rem(i, j, k), 0.0_wp), 1.0_wp)
            num = num + F_3d(i, j, k)*wt
            denom = denom + wt
         end do
         if (denom > 0.0_wp) then
            F_mean_2d(i, j) = num/denom
         else
            F_mean_2d(i, j) = 0.0_wp
         end if
      end do
   end subroutine face_depth_mean_rem_u

   pure subroutine face_depth_mean_rem_v(grid, F_3d, h_layer, rem, F_mean_2d, nz)
      !! Symmetric v-face counterpart of `face_depth_mean_rem_u`.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx,ny+1,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(in) :: rem(:, :, :)      ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      integer :: i, j, k, nx, nv, ny_cells
      real(wp) :: h_face, wt, num, denom

      nx = size(F_3d, 1)
      nv = size(F_3d, 2)
      ny_cells = grid%ny_total

      do concurrent(j=1:nv, i=1:nx) local(k, h_face, wt, num, denom)
         num = 0.0_wp
         denom = 0.0_wp
         do k = 1, nz
            if (j == 1) then
               h_face = h_layer(i, 1, k)
            else if (j == nv) then
               h_face = h_layer(i, ny_cells, k)
            else
               h_face = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            end if
            wt = h_face*min(max(rem(i, j, k), 0.0_wp), 1.0_wp)
            num = num + F_3d(i, j, k)*wt
            denom = denom + wt
         end do
         if (denom > 0.0_wp) then
            F_mean_2d(i, j) = num/denom
         else
            F_mean_2d(i, j) = 0.0_wp
         end if
      end do
   end subroutine face_depth_mean_rem_v

   pure subroutine face_depth_mean_v(grid, F_3d, h_layer, F_mean_2d, nz)
      !! Symmetric v-face counterpart of `face_depth_mean_u`.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx,ny+1,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      integer :: i, j, k, nx, nv, ny_cells
      real(wp) :: h_face, num, denom

      nx = size(F_3d, 1)
      nv = size(F_3d, 2)
      ny_cells = grid%ny_total

      do concurrent(j=1:nv, i=1:nx) local(k, h_face, num, denom)
         num = 0.0_wp
         denom = 0.0_wp
         do k = 1, nz
            if (j == 1) then
               h_face = h_layer(i, 1, k)
            else if (j == nv) then
               h_face = h_layer(i, ny_cells, k)
            else
               h_face = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
            end if
            num = num + F_3d(i, j, k)*h_face
            denom = denom + h_face
         end do
         if (denom > 0.0_wp) then
            F_mean_2d(i, j) = num/denom
         else
            F_mean_2d(i, j) = 0.0_wp
         end if
      end do
   end subroutine face_depth_mean_v

   pure subroutine apply_bt_correction(bt_work, ms, dt, skip_h_rescale, use_h_weighted, &
                                       grid, use_bc_pgf, use_visc_rem, metrics, scale, n_nonfin)
      !! Replace the bt mode in the per-layer face velocities with the
      !! barotropic-substep end-step value, adding Δu = u_bt_end − u_bt_at_n −
      !! dt·F_bt_u to every layer (same for v). Split-explicit convention
      !! (Hallberg 2009): momentum uses the END-of-step barotropic velocity;
      !! layer continuity earlier used the time-mean transports. Both legs of the
      !! corrector use the same end-step anchor (mismatched anchors overshoot the
      !! gravity-wave phase speed). Also rescales `h_layer` uniformly so the
      !! column total matches `H_ref + η_end`. `hTr` is deliberately NOT rescaled
      !! (would break exact tracer mass conservation; T = hTr/h drifts by
      !! O((η_end−η*_slow)/H) per step).
      !!
      !! `skip_h_rescale` — disable the h-rescale (Lagrangian vcoord, where slow
      !!   continuity's Σh_layer is authoritative and the ALE remap relayers).
      !!   Default `.false.`.
      !! `use_h_weighted` — distribute Δu ∝ h_face(k)/⟨h⟩_h (⟨h⟩_h = Σh²/Σh)
      !!   instead of uniformly; preserves depth-mean by construction. No-op for
      !!   uniform-h columns. Default `.false.` (bit-identical).
      !! `use_bc_pgf` — add the per-layer baroclinic-PGF retro-correction
      !!   Δu_bc = -dt·((pbce(R,k)-gtot_W(R))·e_anom(R) -
      !!   (pbce(L,k)-gtot_E(L))·e_anom(L))/dx. Depth-mean zero by construction,
      !!   so the mass-flux invariant survives. Requires `grid`+`metrics` present
      !!   and pbce/gtot_*/e_anom populated by the caller. No-op when omitted.
      type(barotropic_workstate_t), intent(in) :: bt_work
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: skip_h_rescale
      logical, intent(in), optional :: use_h_weighted
      type(hgrid_t), intent(in), optional :: grid
      logical, intent(in), optional :: use_bc_pgf
      logical, intent(in), optional :: use_visc_rem
      type(ocean_metrics_t), intent(in), optional :: metrics
         !! Curvilinear horizontal metrics — required only when
         !! `use_bc_pgf = .true.` (the per-face bc-PGF retro-correction
         !! divides the e_anom gradient by `idxCu`/`idyCv`).
      real(wp), intent(in), optional :: scale
         !! Multiplier on the Δu correction (default 1, bit-identical).
         !! The pred_corr PREDICTOR passes `BE` so the provisional velocity
         !! is `up = u + dt_pred·(u_bc_accel + u_accel_bt)` with
         !! `dt_pred = BE·dt` (SPEC §2 P8) — the tendency applies are
         !! scaled by BE at their call sites, and this scales the
         !! barotropic-increment leg to match.
      integer, intent(out), optional :: n_nonfin
         !! Count of faces whose barotropic-correction Δ (`bt_*_end − *_at_n −
         !! dt·F_bt`) is NON-FINITE.  In a supercritical hot state the BT
         !! substep loop can reach Inf on at-floor columns, and the fold's
         !! `finite − Inf` mints NaN into the layer velocity; the guard SKIPS
         !! the fold write for such a face (leaving its velocity as-is for the
         !! truncation's NaN-catch backstop) and this counts it loudly.  0 on
         !! a healthy run.

      integer :: i, j, k, nu, nv, nx, ny, nz, nfin
      real(wp) :: delta_u, delta_v, total_h_old, total_h_new, ratio
      real(wp) :: du_scale
      real(wp) :: h_face, sum_h, sum_h2, h_bar_h, wt, vr_k
      real(wp) :: du_bc, dv_bc
      logical :: do_rescale, do_h_weighted, do_bc_pgf, do_visc_rem

      do_rescale = .true.
      if (present(skip_h_rescale)) do_rescale = .not. skip_h_rescale
      do_h_weighted = .false.
      if (present(use_h_weighted)) do_h_weighted = use_h_weighted
      do_bc_pgf = .false.
      if (present(use_bc_pgf)) do_bc_pgf = use_bc_pgf
      do_visc_rem = .false.
      if (present(use_visc_rem)) do_visc_rem = use_visc_rem
      du_scale = 1.0_wp
      if (present(scale)) du_scale = scale
      if (do_bc_pgf .and. .not. present(grid)) then
         error stop "apply_bt_correction: use_bc_pgf=.true. requires grid"
      end if
      if (do_bc_pgf .and. .not. present(metrics)) then
         error stop "apply_bt_correction: use_bc_pgf=.true. requires metrics"
      end if

      nu = size(ms%u_face_x_layer, 1)
      ny = size(ms%u_face_x_layer, 2)
      nx = size(ms%v_face_y_layer, 1)
      nv = size(ms%v_face_y_layer, 2)
      nz = ms%nz_ml

      ! Loud count of non-finite fold Δ (2D, once per face) — the BT loop can
      ! reach Inf on at-floor columns in a supercritical state.  Read-only
      ! reduction (write out of the reduction loop, per the truncation fix).
      ! (No `present()` clause: this kernel's fold uses stdpar `do concurrent`
      ! managed memory, so callers do not acc-map bt_work; present_or_copyin
      ! reads the device copy in production and copies-in host data in the
      ! unmapped unit tests.)
      if (present(n_nonfin)) then
         nfin = 0
         do concurrent(j=1:ny, i=1:nu) reduce(+:nfin)
            if (.not. ieee_is_finite(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) &
                                     - dt*bt_work%F_bt_u(i, j))) nfin = nfin + 1
         end do
         do concurrent(j=1:nv, i=1:nx) reduce(+:nfin)
            if (.not. ieee_is_finite(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) &
                                     - dt*bt_work%F_bt_v(i, j))) nfin = nfin + 1
         end do
         n_nonfin = nfin
      end if

      if (.not. do_h_weighted) then
         ! Uniform Δu distribution — every layer gets the same Δu.  The finite
         ! guard skips a face whose Δ is non-finite (Inf/NaN from a blown-up BT
         ! loop) so the fold never mints NaN into the layer velocity.
         do concurrent(k=1:nz, j=1:ny, i=1:nu) local(delta_u)
            delta_u = du_scale*(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) - dt*bt_work%F_bt_u(i, j))
            if (ieee_is_finite(delta_u)) then
               ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u
            end if
         end do
         do concurrent(k=1:nz, j=1:nv, i=1:nx) local(delta_v)
            delta_v = du_scale*(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) - dt*bt_work%F_bt_v(i, j))
            if (ieee_is_finite(delta_v)) then
               ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v
            end if
         end do
      else
         ! H-weighted Δu distribution. Per-layer weight = h_face(k)·vr_k / h_bar_h
         ! with h_bar_h = (Σ h²·vr)/Σ h, vr_k = visc_rem(k) when do_visc_rem else
         ! 1.0 (MOM6 wt_u = frhatu·visc_rem). Preserves the depth-mean; vr≡1
         ! reduces bit-identically to the h-only branch.
         do concurrent(j=1:ny, i=1:nu) &
            local(k, delta_u, sum_h, sum_h2, h_face, h_bar_h, wt, vr_k)
            delta_u = du_scale*(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) - dt*bt_work%F_bt_u(i, j))
            sum_h = 0.0_wp
            sum_h2 = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(ms%h_layer(max(1, i - 1), j, k) + ms%h_layer(min(nu - 1, i), j, k))
               if (do_visc_rem) then
                  vr_k = bt_work%visc_rem_u(i, j, k)
               else
                  vr_k = 1.0_wp
               end if
               sum_h = sum_h + h_face
               sum_h2 = sum_h2 + h_face*h_face*vr_k
            end do
            if (sum_h2 > 0.0_wp .and. ieee_is_finite(delta_u)) then
               h_bar_h = sum_h2/sum_h
               do k = 1, nz
                  h_face = 0.5_wp*(ms%h_layer(max(1, i - 1), j, k) + ms%h_layer(min(nu - 1, i), j, k))
                  if (do_visc_rem) then
                     vr_k = bt_work%visc_rem_u(i, j, k)
                  else
                     vr_k = 1.0_wp
                  end if
                  wt = h_face*vr_k/h_bar_h
                  ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u*wt
               end do
            end if
         end do
         do concurrent(j=1:nv, i=1:nx) &
            local(k, delta_v, sum_h, sum_h2, h_face, h_bar_h, wt, vr_k)
            delta_v = du_scale*(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) - dt*bt_work%F_bt_v(i, j))
            sum_h = 0.0_wp
            sum_h2 = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(ms%h_layer(i, max(1, j - 1), k) + ms%h_layer(i, min(nv - 1, j), k))
               if (do_visc_rem) then
                  vr_k = bt_work%visc_rem_v(i, j, k)
               else
                  vr_k = 1.0_wp
               end if
               sum_h = sum_h + h_face
               sum_h2 = sum_h2 + h_face*h_face*vr_k
            end do
            if (sum_h2 > 0.0_wp .and. ieee_is_finite(delta_v)) then
               h_bar_h = sum_h2/sum_h
               do k = 1, nz
                  h_face = 0.5_wp*(ms%h_layer(i, max(1, j - 1), k) + ms%h_layer(i, min(nv - 1, j), k))
                  if (do_visc_rem) then
                     vr_k = bt_work%visc_rem_v(i, j, k)
                  else
                     vr_k = 1.0_wp
                  end if
                  wt = h_face*vr_k/h_bar_h
                  ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v*wt
               end do
            end if
         end do
      end if

      ! bc-PGF additive per-layer correction. L = west/south cell, R =
      ! east/north cell. Each (pbce(k) - gtot_face) term is depth-mean-zero per
      ! column, so the column-mean velocity is unchanged (mass-flux invariant
      ! preserved). Interior faces only (wall faces already BT-zeroed).
      if (do_bc_pgf) then
         do concurrent(k=1:nz, j=1:ny, i=2:nu - 1) local(du_bc)
            du_bc = -dt*((bt_work%pbce(i, j, k) - bt_work%gtot_W(i, j)) &
                         *bt_work%e_anom(i, j) &
                         - (bt_work%pbce(i - 1, j, k) - bt_work%gtot_E(i - 1, j)) &
                         *bt_work%e_anom(i - 1, j))*metrics%idxCu(i, j)
            if (ieee_is_finite(du_bc)) ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + du_bc
         end do
         do concurrent(k=1:nz, j=2:nv - 1, i=1:nx) local(dv_bc)
            dv_bc = -dt*((bt_work%pbce(i, j, k) - bt_work%gtot_S(i, j)) &
                         *bt_work%e_anom(i, j) &
                         - (bt_work%pbce(i, j - 1, k) - bt_work%gtot_N(i, j - 1)) &
                         *bt_work%e_anom(i, j - 1))*metrics%idyCv(i, j)
            if (ieee_is_finite(dv_bc)) ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + dv_bc
         end do
      end if

      if (do_rescale) then
         do concurrent(j=1:size(ms%h_layer, 2), i=1:size(ms%h_layer, 1)) &
            local(k, total_h_old, total_h_new, ratio)
            total_h_old = 0.0_wp
            do k = 1, nz
               total_h_old = total_h_old + ms%h_layer(i, j, k)
            end do
            total_h_new = bt_work%bt_H_ref(i, j) + bt_work%bt_eta_end(i, j)
            if (total_h_old > 0.0_wp) then
               ratio = total_h_new/total_h_old
               do k = 1, nz
                  ms%h_layer(i, j, k) = ms%h_layer(i, j, k)*ratio
               end do
            end if
         end do
      end if
   end subroutine apply_bt_correction

   pure subroutine snapshot_eta_PF(bt_work)
      !! Snapshot `bt_eta` into `eta_PF` — the free-surface height the slow PGF
      !! sees this stage. Later differenced by `compute_e_anom`.
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer :: i, j, nx, ny
      nx = size(bt_work%bt_eta, 1)
      ny = size(bt_work%bt_eta, 2)
      do concurrent(j=1:ny, i=1:nx)
         bt_work%eta_PF(i, j) = bt_work%bt_eta(i, j)
      end do
   end subroutine snapshot_eta_PF

   pure subroutine compute_pbce(grid, bt_work, pgf, ms)
      !! Per-layer pressure-anomaly gravity coefficient (m/s²): the response of
      !! layer k's pressure to a unit change in η. Montgomery form, bottom-up
      !! convention (k=1 bed, k=nz surface):
      !!     pbce(:,:,nz) = g·ρ_ref/ρ_0
      !!     do k = nz-1, 1, -1
      !!        g_prime_K = g·(rho_layer(k+1) − rho_layer(k))/ρ_0
      !!        pbce(:,:,k) = pbce(:,:,k+1) + g_prime_K·(e_top_of_k − e_bed)/H
      !! Uniform-density column ⇒ pbce−gtot ≡ 0 ⇒ bc-PGF correction a no-op.
      !! Reads `pgf%e_face`; requires `pgf%variant == OPGF_VARIANT_FV_MOM6`
      !! (error_stop otherwise — other variants don't fill e_face).
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(ocean_pressure_force_t), intent(in) :: pgf
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: g_surf, inv_rho0, h_col, g_prime_K, e_above, e_bed

      if (pgf%variant /= OPGF_VARIANT_FV_MOM6) then
         error stop "compute_pbce: requires ocean_pgf_form = 'fv_mom6' "// &
            "(pgf%e_face is not populated by other variants)."
      end if

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      inv_rho0 = 1.0_wp/pgf%rho0
      g_surf = GRAVITY*pgf%rho_ref*inv_rho0

      do concurrent(j=1:ny, i=1:nx) local(k, h_col, g_prime_K, e_above, e_bed)
         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + ms%h_layer(i, j, k)
         end do
         if (h_col <= 0.0_wp) h_col = 1.0_wp     ! dry cell — pbce won't be used
         e_bed = pgf%e_face%data(i, j, 1)
         bt_work%pbce(i, j, nz) = g_surf
         do k = nz - 1, 1, -1
            g_prime_K = GRAVITY*(ms%rho_layer(i, j, k + 1) - ms%rho_layer(i, j, k))*inv_rho0
            e_above = pgf%e_face%data(i, j, k + 1)
            bt_work%pbce(i, j, k) = bt_work%pbce(i, j, k + 1) &
                                    + g_prime_K*(e_above - e_bed)/h_col
         end do
      end do
   end subroutine compute_pbce

   pure subroutine compute_gtot_faces(grid, bt_work, ms)
      !! Face-centred depth-weighted column averages of `pbce` (gtot_E/W/N/S).
      !! Wall cells fall back to `pbce(:,:,nz)`. By construction
      !! Σ_k h_face(k)·(pbce(k) − gtot_face) = 0 per column, making the bc-PGF
      !! Δu correction depth-mean zero.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: h_face, h_sum, p_sum

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      do concurrent(j=1:ny, i=1:nx - 1) local(k, h_face, h_sum, p_sum)
         h_sum = 0.0_wp
         p_sum = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i + 1, j, k))
            h_sum = h_sum + h_face
            p_sum = p_sum + h_face*bt_work%pbce(i, j, k)
         end do
         if (h_sum > 0.0_wp) then
            bt_work%gtot_E(i, j) = p_sum/h_sum
         else
            bt_work%gtot_E(i, j) = bt_work%pbce(i, j, nz)
         end if
      end do
      do concurrent(j=1:ny)
         bt_work%gtot_E(nx, j) = bt_work%pbce(nx, j, nz)
      end do

      do concurrent(j=1:ny, i=2:nx) local(k, h_face, h_sum, p_sum)
         h_sum = 0.0_wp
         p_sum = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
            h_sum = h_sum + h_face
            p_sum = p_sum + h_face*bt_work%pbce(i, j, k)
         end do
         if (h_sum > 0.0_wp) then
            bt_work%gtot_W(i, j) = p_sum/h_sum
         else
            bt_work%gtot_W(i, j) = bt_work%pbce(i, j, nz)
         end if
      end do
      do concurrent(j=1:ny)
         bt_work%gtot_W(1, j) = bt_work%pbce(1, j, nz)
      end do

      do concurrent(j=1:ny - 1, i=1:nx) local(k, h_face, h_sum, p_sum)
         h_sum = 0.0_wp
         p_sum = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i, j, k) + ms%h_layer(i, j + 1, k))
            h_sum = h_sum + h_face
            p_sum = p_sum + h_face*bt_work%pbce(i, j, k)
         end do
         if (h_sum > 0.0_wp) then
            bt_work%gtot_N(i, j) = p_sum/h_sum
         else
            bt_work%gtot_N(i, j) = bt_work%pbce(i, j, nz)
         end if
      end do
      do concurrent(i=1:nx)
         bt_work%gtot_N(i, ny) = bt_work%pbce(i, ny, nz)
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, h_face, h_sum, p_sum)
         h_sum = 0.0_wp
         p_sum = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
            h_sum = h_sum + h_face
            p_sum = p_sum + h_face*bt_work%pbce(i, j, k)
         end do
         if (h_sum > 0.0_wp) then
            bt_work%gtot_S(i, j) = p_sum/h_sum
         else
            bt_work%gtot_S(i, j) = bt_work%pbce(i, j, nz)
         end if
      end do
      do concurrent(i=1:nx)
         bt_work%gtot_S(i, 1) = bt_work%pbce(i, 1, nz)
      end do
   end subroutine compute_gtot_faces

   pure subroutine compute_e_anom(bt_work)
      !! SSH anomaly = 0.5·(bt_eta_end + bt_eta) − eta_PF: the part of η the BT
      !! substep produced beyond what the slow PGF saw. Zero at steady state.
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer :: i, j, nx, ny
      nx = size(bt_work%e_anom, 1)
      ny = size(bt_work%e_anom, 2)
      do concurrent(j=1:ny, i=1:nx)
         bt_work%e_anom(i, j) = 0.5_wp*(bt_work%bt_eta_end(i, j) + bt_work%bt_eta(i, j)) &
                                - bt_work%eta_PF(i, j)
      end do
   end subroutine compute_e_anom

   pure subroutine compute_bt_rem(grid, bt_work, ms, r_linear, hbbl, dt_inner)
      !! Per-face multiplicative damping factor for the BT-substep velocity
      !! update (linear-drag branch):
      !!     bt_rem_face = Htot_face / (Htot_face + r·hbbl·dt_inner)
      !! applied as ubt_new = bt_rem_u·(ubt_old + dt_inner·forces) each inner
      !! step. When the `bt_substep_drag` knob is off this must NOT be called and
      !! the workspace stays at 1 (no-op, bit-identical).
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: r_linear   !! Linear drag rate at the bed (1/s)
      real(wp), intent(in) :: hbbl       !! BBL thickness over which drag acts (m)
      real(wp), intent(in) :: dt_inner   !! BT-substep dt (s)

      integer :: i, j, k, nx, ny, nz
      real(wp) :: htot_face, drag_dt

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      drag_dt = r_linear*hbbl*dt_inner   ! product is in metres

      do concurrent(j=1:ny, i=2:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
         end do
         if (htot_face > 0.0_wp) then
            bt_work%bt_rem_u(i, j) = htot_face/(htot_face + drag_dt)
         else
            bt_work%bt_rem_u(i, j) = 1.0_wp
         end if
      end do
      do concurrent(j=1:ny)
         bt_work%bt_rem_u(1, j) = 1.0_wp
         bt_work%bt_rem_u(nx + 1, j) = 1.0_wp
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
         end do
         if (htot_face > 0.0_wp) then
            bt_work%bt_rem_v(i, j) = htot_face/(htot_face + drag_dt)
         else
            bt_work%bt_rem_v(i, j) = 1.0_wp
         end if
      end do
      do concurrent(i=1:nx)
         bt_work%bt_rem_v(i, 1) = 1.0_wp
         bt_work%bt_rem_v(i, ny + 1) = 1.0_wp
      end do
   end subroutine compute_bt_rem

   pure subroutine reset_bt_rem(grid, bt_work)
      !! bt_rem_u/v ≡ 1 (the init value). `bt_rem_u/v` is otherwise reset
      !! only by `compute_bt_rem`, which only runs when `bt_substep_drag`
      !! is on. `compute_bt_rem_wave_drag` MULTIPLIES into `bt_rem_u/v`,
      !! so when wave drag is on and `bt_substep_drag` is off, something
      !! must still reset it to 1 each stage — otherwise it compounds
      !! geometrically across outer steps (bt_rem = R^n after n stages),
      !! silently annihilating the barotropic mode. See
      !! `src/core/ocean/README.md` for the multiplicative-accumulator
      !! contract this establishes.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_work%bt_rem_u(i, j) = 1.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_work%bt_rem_v(i, j) = 1.0_wp
      end do
   end subroutine reset_bt_rem

   pure subroutine compute_bt_rem_wave_drag(grid, bt_work, ms, dt_inner)
      !! MULTIPLIES the Egbert & Ray (2001) / Jayne & St Laurent (2001)
      !! linear (Rayleigh) barotropic wave drag into `bt_rem_u/v`:
      !!     bt_rem_u *= Htot_face / (Htot_face + lwd_drag_u·dt_inner)
      !! `lwd_drag_u/v` is a static, face-resident piston velocity [m/s]
      !! built once at configure by `configure_ocean_wave_drag`. Uses the
      !! IDENTICAL `Htot_face` expression as `compute_bt_rem` (reuse, not
      !! a second `H_tot`). Composes with `substep_drag` exactly as MOM6
      !! composes `lin_drag_u` with the viscous remnant. `Htot_face <= 0`
      !! ⇒ leave `bt_rem` unmodified (MOM6's guard).
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt_inner   !! BT-substep dt (s)

      integer :: i, j, k, nx, ny, nz
      real(wp) :: htot_face

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      do concurrent(j=1:ny, i=2:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
         end do
         if (htot_face > 0.0_wp) then
            bt_work%bt_rem_u(i, j) = bt_work%bt_rem_u(i, j)* &
                                     (htot_face/(htot_face + bt_work%lwd_drag_u(i, j)*dt_inner))
         end if
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
         end do
         if (htot_face > 0.0_wp) then
            bt_work%bt_rem_v(i, j) = bt_work%bt_rem_v(i, j)* &
                                     (htot_face/(htot_face + bt_work%lwd_drag_v(i, j)*dt_inner))
         end if
      end do
   end subroutine compute_bt_rem_wave_drag

   pure subroutine mask_bt_rem(grid, metrics, bt_work)
      !! Fold the static land face masks into the BT-substep damping factor
      !! (bt_rem_u(land)=0 ⇒ no velocity across a land face). Runs every outer
      !! step AFTER `compute_bt_rem` (which resets bt_rem each step, so the mask
      !! must be re-applied). All-wet ⇒ wet_u/v≡1 ⇒ no-op (bit-identical).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_work%bt_rem_u(i, j) = metrics%wet_u(i, j)*bt_work%bt_rem_u(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_work%bt_rem_v(i, j) = metrics%wet_v(i, j)*bt_work%bt_rem_v(i, j)
      end do
   end subroutine mask_bt_rem

   pure subroutine set_local_BT_cont_types(grid, metrics, bt_work, ms, dt_outer)
      !! Populate `bt_work%BTCL_u/v` — the per-face flux-closure coefficients
      !! consumed by `find_uhbt` — from the current `h_layer`. Upstream-h-sum
      !! approach: FA_u_W0=FA_u_E0=Σ_k h_face (centred); FA_u_WW=Σ_k h_layer(west)
      !! and FA_u_EE=Σ_k h_layer(east) (saturated-regime upstream draw); uBT_WW/EE
      !! = ±VOL_CFL·dx/dt_outer pin the saturation velocity; uh_crv*/uh_** are the
      !! C¹-matching coefficients (Hallberg & Adcroft 2009). No-op when
      !! `use_bt_cont_type = .false.` (BTCL_u/v unallocated).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt_outer

      real(wp), parameter :: C1_3 = 1.0_wp/3.0_wp
      integer :: i, j, k, nx, ny, nz
      real(wp) :: fa_centre, fa_up_W, fa_up_E, fa_up_N, fa_up_S
      real(wp) :: h_face, u_cfl_x, u_cfl_y, inv_ucfl_x2, inv_ucfl_y2
      real(wp) :: inv_dt
      ! u_cfl_x/y + inv_ucfl_x2/y2 are written per-iteration as DC locals
      ! (declared here so the `local()` clause can name them).

      ! Boundary u/v-faces are not written — they keep their type-default zero;
      ! find_uhbt(0, zero_BTC) = 0, matching the wall-zero ubt they carry.

      if (.not. bt_work%use_bt_cont_type) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      inv_dt = 1.0_wp/dt_outer

      ! ---- u-faces (i in 2..nx) ----
      ! Per-face CFL velocity BTC_VOL_CFL·dxCu/dt_outer.
      do concurrent(j=1:ny, i=2:nx) &
         local(k, fa_centre, fa_up_W, fa_up_E, h_face, u_cfl_x, inv_ucfl_x2)
         u_cfl_x = BTC_VOL_CFL*metrics%dxCu(i, j)*inv_dt
         inv_ucfl_x2 = 0.0_wp
         if (u_cfl_x > 0.0_wp) inv_ucfl_x2 = 1.0_wp/(u_cfl_x*u_cfl_x)
         fa_centre = 0.0_wp
         fa_up_W = 0.0_wp
         fa_up_E = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
            if (h_face > BTC_H_NEGLECT) fa_centre = fa_centre + h_face
            if (ms%h_layer(i - 1, j, k) > BTC_H_NEGLECT) then
               fa_up_W = fa_up_W + ms%h_layer(i - 1, j, k)
            end if
            if (ms%h_layer(i, j, k) > BTC_H_NEGLECT) then
               fa_up_E = fa_up_E + ms%h_layer(i, j, k)
            end if
         end do

         bt_work%BTCL_u(i, j)%FA_u_W0 = fa_centre
         bt_work%BTCL_u(i, j)%FA_u_E0 = fa_centre
         bt_work%BTCL_u(i, j)%FA_u_WW = fa_up_W
         bt_work%BTCL_u(i, j)%FA_u_EE = fa_up_E
         bt_work%BTCL_u(i, j)%uBT_WW = u_cfl_x
         bt_work%BTCL_u(i, j)%uBT_EE = -u_cfl_x
         bt_work%BTCL_u(i, j)%uh_crvW = C1_3*(fa_up_W - fa_centre)*inv_ucfl_x2
         bt_work%BTCL_u(i, j)%uh_crvE = C1_3*(fa_up_E - fa_centre)*inv_ucfl_x2
         bt_work%BTCL_u(i, j)%uh_WW = u_cfl_x*C1_3*(2.0_wp*fa_centre + fa_up_W)
         bt_work%BTCL_u(i, j)%uh_EE = -u_cfl_x*C1_3*(2.0_wp*fa_centre + fa_up_E)
      end do

      ! ---- v-faces (j in 2..ny) ----
      do concurrent(j=2:ny, i=1:nx) &
         local(k, fa_centre, fa_up_N, fa_up_S, h_face, u_cfl_y, inv_ucfl_y2)
         u_cfl_y = BTC_VOL_CFL*metrics%dyCv(i, j)*inv_dt
         inv_ucfl_y2 = 0.0_wp
         if (u_cfl_y > 0.0_wp) inv_ucfl_y2 = 1.0_wp/(u_cfl_y*u_cfl_y)
         fa_centre = 0.0_wp
         fa_up_N = 0.0_wp
         fa_up_S = 0.0_wp
         do k = 1, nz
            h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
            if (h_face > BTC_H_NEGLECT) fa_centre = fa_centre + h_face
            if (ms%h_layer(i, j - 1, k) > BTC_H_NEGLECT) then
               fa_up_S = fa_up_S + ms%h_layer(i, j - 1, k)
            end if
            if (ms%h_layer(i, j, k) > BTC_H_NEGLECT) then
               fa_up_N = fa_up_N + ms%h_layer(i, j, k)
            end if
         end do

         ! Convention: vBT_SS > 0 → southward-draw saturation (v > 0 flow,
         ! upstream is cell j-1, the "S" column).  vBT_NN < 0 → northward-
         ! draw saturation.  Mirrors MOM6's local_BT_cont_v_type.
         bt_work%BTCL_v(i, j)%FA_v_S0 = fa_centre
         bt_work%BTCL_v(i, j)%FA_v_N0 = fa_centre
         bt_work%BTCL_v(i, j)%FA_v_SS = fa_up_S
         bt_work%BTCL_v(i, j)%FA_v_NN = fa_up_N
         bt_work%BTCL_v(i, j)%vBT_SS = u_cfl_y
         bt_work%BTCL_v(i, j)%vBT_NN = -u_cfl_y
         bt_work%BTCL_v(i, j)%vh_crvS = C1_3*(fa_up_S - fa_centre)*inv_ucfl_y2
         bt_work%BTCL_v(i, j)%vh_crvN = C1_3*(fa_up_N - fa_centre)*inv_ucfl_y2
         bt_work%BTCL_v(i, j)%vh_SS = u_cfl_y*C1_3*(2.0_wp*fa_centre + fa_up_S)
         bt_work%BTCL_v(i, j)%vh_NN = -u_cfl_y*C1_3*(2.0_wp*fa_centre + fa_up_N)
      end do

      ! Boundary u-faces (i=1, i=nx_face) and v-faces (j=1, j=ny_face)
      ! stay at their type-default zero — wall faces in our setup carry
      ! ubt=0 and find_uhbt(0, anything)=0, so no flux through them.
   end subroutine set_local_BT_cont_types

end module rdb_barotropic_coupling
