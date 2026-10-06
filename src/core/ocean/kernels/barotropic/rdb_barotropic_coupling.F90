!! Barotropic ↔ multilayer coupling kernels for the split-explicit ocean
!! driver — the bt↔ml bridges run around each barotropic-substep call:
!! `derive_bt_from_layers` (layers → bt state, depth-averaged),
!! `sum_slow_tendencies_into_F_slow` (per-kernel scratch → composite F_slow),
!! `face_depth_mean_u/_v` (3D slow tendency → 2D face forcing),
!! `apply_bt_correction` (bt time-mean → per-layer correction + h rescale).
module rdb_barotropic_coupling
   use rdb_constants, only: wp, GRAVITY, H_DIV_EPS, FRHAT_HYBRID
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_MOM6, &
                                       OPGF_VARIANT_GPRIME
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_top_drag, only: ocean_top_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_boundary_types, only: OBC_PERIODIC, OBC_TRIPOLAR_FOLD
   implicit none
   private

   public :: derive_bt_from_layers
   public :: compute_h_face_upstream
   public :: sum_slow_tendencies_into_F_slow
   public :: add_top_drag_into_F_slow
   public :: subtract_fast_cor_ref
   public :: set_cor_ref_velocity
   public :: face_depth_mean_u
   public :: face_depth_mean_v
   public :: face_depth_mean_rem_u
   public :: face_depth_mean_rem_v
   public :: apply_bt_correction
   public :: snapshot_eta_PF
   public :: set_fast_forcing_eta_pf
   public :: pgf_free_surface_gravity
   public :: compute_pbce
   public :: compute_gtot_faces
   public :: compute_e_anom
   public :: compute_bt_rem
   public :: reset_bt_rem
   public :: compute_bt_rem_from_visc_rem
   public :: compute_bt_rem_wave_drag
   public :: mask_bt_rem
   ! BT_cont_type producer — fills BTCL_u/v on bt_work
   public :: set_local_BT_cont_types
   ! frhat port (src/shared_module_utilities/rdb_frhat_face.inc) -- exported
   ! for its own direct unit tests; every production call site reaches it
   ! through the routines above instead.
   public :: frhat_h_face_step

   ! Floor below which a layer counts as "vanished" for FA accumulation.
   real(wp), parameter :: BTC_H_NEGLECT = 1.0e-10_wp
   ! Volume-CFL safety factor (vol_CFL = 0.5). NOTE: too loose for the
   ! upstream-h-sum producer — the closure stays in the cubic-near-zero
   ! branch (≈ naive u·h_face); real saturation needs a PPM-perturbation FA.
   real(wp), parameter :: BTC_VOL_CFL = 0.5_wp
   ! MOM6 `wt_u`'s round-off guard (MOM_barotropic.F90 module parameter
   ! `subroundoff`, :467) — only ever used inside the `wt_u` floor below.
   real(wp), parameter :: VISC_REM_SUBROUNDOFF = 1.0e-30_wp

contains

   pure subroutine derive_bt_from_layers(grid, bt_work, ms, metrics)
      !! Populate `bt_eta`, `bt_ubt`, `bt_vbt` from the current multilayer
      !! state. `bt_H_ref` must already be set.
      !!   bt_eta = Σ_k h_layer − H_ref;  bt_ubt = Σ_k(u·h_face)/Σ_k h_face.
      !! Face thickness averages the two abutting columns (wall faces use the
      !! single cell). With `use_upstream_h_face`, the interior face-h is the
      !! first-order upwind pick — consistent with `compute_h_face_upstream`.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED — and required on purpose.  Under `&vcoord_nml
         !! zfixed_closed_faces` a CLOSED layer carries no transport, so
         !! it must not dilute the face mean either: the weight becomes
         !! `h_face·open` and `ubt` is the OPEN-column depth mean.  That
         !! is not a refinement, it is a consistency requirement: the
         !! barotropic substep transports on `ubt·FA·dy_cu_bt` with
         !! `dy_cu_bt` narrowed by the open fraction, so the `ubt` the
         !! fast loop integrates ALREADY means "the open-column mean".
         !! Deriving `ubt_at_n` from the full column would make the
         !! fold's `Δu = ubt_end − ubt_at_n − dt·F_bt` a difference
         !! between two different quantities — diluted by `0.5·h_live`
         !! per one-sided-filler layer, which at a partial face is not
         !! small.
         !!
         !! This dummy was OPTIONAL for one release and that is exactly
         !! how the `pred_corr` Coriolis-reference defect shipped: one
         !! call site in `set_cor_ref_velocity` omitted it, silently took
         !! the full-column branch, and the un-cancelled `f·(1−φ)·v̄`
         !! forced every barotropic substep.  An optional argument that
         !! silently changes the physics is a defect CLASS, not a defect;
         !! passing it is now mandatory so a new call site cannot quietly
         !! take the wrong branch.
         !!
         !! `metrics%use_closed_faces = .false.` (the default) ⇒ the
         !! ORIGINAL loops run, textually unchanged (byte-identical).  The
         !! open branch is written out in full rather than folded into the
         !! original with a runtime `if` so that the default-path kernel
         !! never NAMES `open_u`/`open_v` at all — with the knob off those
         !! are `(1,1,1)` placeholders, and a placeholder indexed inside a
         !! `do concurrent` is exactly the kind of thing that works on the
         !! host and faults under `mem:separate`.

      integer :: i, j, k, nx, ny, nz, nx_face, ny_face, il, ir, jl, jr
      logical :: use_upstream, use_open
      real(wp) :: total_h, hu_sum, h_face_sum, h_face, hv_sum
      real(wp) :: href_l, href_r, e_prev
      integer :: scheme

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_face = size(ms%u_face_x_layer, 1)
      ny_face = size(ms%v_face_y_layer, 2)
      use_upstream = bt_work%use_upstream_h_face
      use_open = metrics%use_closed_faces
      scheme = bt_work%frhat_scheme

      do concurrent(j=1:ny, i=1:nx) local(k, total_h)
         total_h = 0.0_wp
         do k = 1, nz
            total_h = total_h + ms%h_layer(i, j, k)
         end do
         bt_work%bt_eta(i, j) = total_h - bt_work%bt_H_ref(i, j)
      end do

      ! `il`/`ir` are the two columns a u-face abuts (equal at a wall, where
      ! `frhat_h_face_step` degenerates to the single available column —
      ! see that routine's docstring).  `use_upstream` (no MOM6 frhat
      ! counterpart — a roundabout-only alternative orthogonal to
      ! centred/hybrid) bypasses the frhat sweep entirely and is untouched.
      if (use_open) then
         do concurrent(j=1:ny, i=1:nx_face) &
            local(k, il, ir, href_l, href_r, e_prev, hu_sum, h_face_sum, h_face)
            il = max(1, i - 1)
            ir = min(nx, i)
            href_l = bt_work%bt_H_ref(il, j)
            href_r = bt_work%bt_H_ref(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            hu_sum = 0.0_wp
            h_face_sum = 0.0_wp
            do k = 1, nz
               if (use_upstream .and. i > 1 .and. i < nx_face) then
                  if (ms%u_face_x_layer(i, j, k) >= 0.0_wp) then
                     h_face = ms%h_layer(i - 1, j, k)
                  else
                     h_face = ms%h_layer(i, j, k)
                  end if
               else
                  call frhat_h_face_step(ms%h_layer(il, j, k), ms%h_layer(ir, j, k), &
                                         href_l, href_r, scheme, e_prev, h_face)
               end if
               h_face = h_face*metrics%open_u(i, j, k)
               hu_sum = hu_sum + ms%u_face_x_layer(i, j, k)*h_face
               h_face_sum = h_face_sum + h_face
            end do
            if (h_face_sum > 0.0_wp) then
               bt_work%bt_ubt(i, j) = hu_sum/h_face_sum
            else
               bt_work%bt_ubt(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:ny, i=1:nx_face) &
            local(k, il, ir, href_l, href_r, e_prev, hu_sum, h_face_sum, h_face)
            il = max(1, i - 1)
            ir = min(nx, i)
            href_l = bt_work%bt_H_ref(il, j)
            href_r = bt_work%bt_H_ref(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            hu_sum = 0.0_wp
            h_face_sum = 0.0_wp
            do k = 1, nz
               if (use_upstream .and. i > 1 .and. i < nx_face) then
                  if (ms%u_face_x_layer(i, j, k) >= 0.0_wp) then
                     h_face = ms%h_layer(i - 1, j, k)
                  else
                     h_face = ms%h_layer(i, j, k)
                  end if
               else
                  call frhat_h_face_step(ms%h_layer(il, j, k), ms%h_layer(ir, j, k), &
                                         href_l, href_r, scheme, e_prev, h_face)
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
      end if

      if (use_open) then
         do concurrent(j=1:ny_face, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, hv_sum, h_face_sum, h_face)
            jl = max(1, j - 1)
            jr = min(ny, j)
            href_l = bt_work%bt_H_ref(i, jl)
            href_r = bt_work%bt_H_ref(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            hv_sum = 0.0_wp
            h_face_sum = 0.0_wp
            do k = 1, nz
               if (use_upstream .and. j > 1 .and. j < ny_face) then
                  if (ms%v_face_y_layer(i, j, k) >= 0.0_wp) then
                     h_face = ms%h_layer(i, j - 1, k)
                  else
                     h_face = ms%h_layer(i, j, k)
                  end if
               else
                  call frhat_h_face_step(ms%h_layer(i, jl, k), ms%h_layer(i, jr, k), &
                                         href_l, href_r, scheme, e_prev, h_face)
               end if
               h_face = h_face*metrics%open_v(i, j, k)
               hv_sum = hv_sum + ms%v_face_y_layer(i, j, k)*h_face
               h_face_sum = h_face_sum + h_face
            end do
            if (h_face_sum > 0.0_wp) then
               bt_work%bt_vbt(i, j) = hv_sum/h_face_sum
            else
               bt_work%bt_vbt(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:ny_face, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, hv_sum, h_face_sum, h_face)
            jl = max(1, j - 1)
            jr = min(ny, j)
            href_l = bt_work%bt_H_ref(i, jl)
            href_r = bt_work%bt_H_ref(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            hv_sum = 0.0_wp
            h_face_sum = 0.0_wp
            do k = 1, nz
               if (use_upstream .and. j > 1 .and. j < ny_face) then
                  if (ms%v_face_y_layer(i, j, k) >= 0.0_wp) then
                     h_face = ms%h_layer(i, j - 1, k)
                  else
                     h_face = ms%h_layer(i, j, k)
                  end if
               else
                  call frhat_h_face_step(ms%h_layer(i, jl, k), ms%h_layer(i, jr, k), &
                                         href_l, href_r, scheme, e_prev, h_face)
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
      end if
   end subroutine derive_bt_from_layers

   pure subroutine compute_h_face_upstream(grid, bt_work, ms, metrics)
      !! Per-face upstream column-sum thickness `h_face_up_x/y(I,j) =
      !! Σ_k h_layer(I_upstream,j,k)` used by the BT chain when
      !! `use_upstream_h_face = .true.`. First-order upwind pick by face-velocity
      !! sign (sampled at the top of the outer step). Wall faces use the single
      !! available cell. No-op when the knob is off.
      !!
      !! Under `&vcoord_nml zfixed_closed_faces` the sum is over the OPEN
      !! column instead — see `metrics` and `h_face_upstream_open_impl`.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED, for the reason `derive_bt_from_layers` gives: an
         !! optional dummy that silently selects the full-column branch is
         !! how the `pred_corr` Coriolis-reference defect shipped.
         !!
         !! `metrics%use_closed_faces = .false.` (the default) ⇒ the
         !! ORIGINAL full-column loops below run, textually unchanged
         !! (byte-identical), and neither `open_*` nor `dy_cu_bt` is
         !! named.  `.true.` ⇒ the open-column builder
         !! `h_face_upstream_open_impl`, whose docstring derives why the
         !! full-column sum is not merely imprecise there but makes the
         !! barotropic transport disagree with the renormalised layer
         !! transports by a factor up to ~2 at a staircase face.

      integer :: i, j, k, nx, ny, nz, nx_face, ny_face
      real(wp) :: h_sum, h_k

      if (.not. bt_work%use_upstream_h_face) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_face = size(bt_work%h_face_up_x, 1)
      ny_face = size(bt_work%h_face_up_y, 2)

      if (metrics%use_closed_faces) then
         ! Two calls, one per porous state: with porous OFF the
         ! `por_face_area_*` arrays are the `(1,1,1)` placeholder and must
         ! NOT reach the callee's explicit-shape dummy (nvfortran builds
         ! the `do concurrent` data clause from the loop bounds, so a
         ! placeholder aborts under `mem:separate` even when the branch
         ! that indexes it is never taken).  The mask itself is the inert
         ! stand-in — right shape, already mapped, `intent(in)` at both
         ! dummies — the same device the `ocean_porous_refresh` call uses.
         if (metrics%use_porous) then
            call h_face_upstream_open_impl(nx, ny, nz, .true., ms%h_layer, &
                                           ms%u_face_x_layer, ms%v_face_y_layer, &
                                           metrics%open_u, metrics%open_v, &
                                           metrics%por_face_area_u, &
                                           metrics%por_face_area_v, &
                                           metrics%dy_cu, metrics%dx_cv, &
                                           metrics%dy_cu_bt, metrics%dx_cv_bt, &
                                           bt_work%h_face_up_x, bt_work%h_face_up_y)
         else
            call h_face_upstream_open_impl(nx, ny, nz, .false., ms%h_layer, &
                                           ms%u_face_x_layer, ms%v_face_y_layer, &
                                           metrics%open_u, metrics%open_v, &
                                           metrics%open_u, metrics%open_v, &
                                           metrics%dy_cu, metrics%dx_cv, &
                                           metrics%dy_cu_bt, metrics%dx_cv_bt, &
                                           bt_work%h_face_up_x, bt_work%h_face_up_y)
         end if
         return
      end if

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

   pure subroutine h_face_upstream_open_impl(nx, ny, nz, use_por, h_layer, &
                                             u_face, v_face, open_u, open_v, &
                                             por_u, por_v, dy_cu, dx_cv, &
                                             dy_cu_bt, dx_cv_bt, h_up_x, h_up_y)
      !! OPEN-column upstream face thickness, for `upstream_h_face` under
      !! `&vcoord_nml zfixed_closed_faces`.
      !!
      !! ### What the barotropic transport has to equal
      !!
      !! Everything else on the barotropic path under closed faces already
      !! means "the OPEN column": `derive_bt_from_layers` builds `ubt` as
      !! `Σ_k u_k·h_up,k·open_k / Σ_k h_up,k·open_k` (upstream pick under
      !! this knob), and the renormaliser hands the fast loop's `uhbt` to
      !! the layers with weight `dy_cu·por·open` on the upstream PPM
      !! thickness.  For a depth-uniform open-layer velocity `u` the layers
      !! therefore carry `dy_cu·Σ_k h_up,k·por_k·open_k·u`, and the fast
      !! loop must transport EXACTLY that, or the renormaliser's `du` is
      !! not zero: `u_av = r·u` with `r` = (barotropic face depth)/(open
      !! upstream depth), and the slow tendencies of the next stage are
      !! evaluated on a velocity the barotropic mode never had.
      !!
      !! ### Why the full-column sum is wrong by O(1), not by round-off
      !!
      !! The fast loop transports `h_face_up·ubt·dy_cu_bt`, and
      !! `dy_cu_bt = dy_cu·φ_c` carries the CENTRED open fraction
      !! `φ_c = Σ h_c·por·open / Σ h_c` (`closed_faces_update_bt_widths`).
      !! The full-column upstream sum `H_up` times `φ_c` is the open depth
      !! only when `H_up = Σ h_c` — a flat face.  At a staircase face
      !! between a deep column `H_D` and a shallow one `H_S` (`h_c` of a
      !! layer live on one side only is `≈ h/2`):
      !!
      !! ```
      !! φ_c        = H_S / ((H_D + H_S)/2)
      !! from deep:    r = H_D·φ_c / H_S = 2·H_D/(H_D + H_S)    -> 2
      !! from shallow: r = H_S·φ_c / H_S = 2·H_S/(H_D + H_S)    -> 0
      !! ```
      !!
      !! (prototype: `python_prototypes/bt_upstream_zfixed/`).  Measured on
      !! the 1-degree Southern Ocean (`zfixed_audit/bt_upstream_h_face`):
      !! barotropic velocity 2.7 m/s by step 13 against 0.84 m/s knob-off,
      !! then the `maxvel` clamp, then NaN at step 309; with the open
      !! column it runs the 10 days at `En 5.709E-04` against the knob-off
      !! `5.493E-04`.
      !!
      !! ### What is stored
      !!
      !! The substep multiplies `h_face_up` by the NARROWED width
      !! `dy_cu_bt`, which already carries `φ_c`.  Storing the open sum
      !! `s = Σ_k h_up,k·por_k·open_k` itself would count the open fraction
      !! twice, so the stored value is the full-column EQUIVALENT
      !!
      !! ```
      !! h_face_up = s · dy_cu / dy_cu_bt     ⇒   h_face_up·dy_cu_bt = s·dy_cu
      !! ```
      !!
      !! which makes the fast-loop transport the open-column upstream
      !! transport to round-off, on BOTH substep kernels, without touching
      !! either (the hot nonlinear kernel receives `dy_cu_bt` as its only
      !! width).  Read against the WIDTH the substep will actually use, so
      !! it is exact whichever step `dy_cu_bt` was refreshed at.  A face
      !! with `dy_cu_bt = 0` (every layer closed, or land) transports
      !! nothing whatever is stored, and stores `s` (`= 0` when every layer
      !! is closed).
      !!
      !! Porous barriers (`use_por`) enter `s` the way they enter the
      !! renormaliser's weight; `por_u`/`por_v` are never indexed when
      !! `.false.` (the caller hands over the mask as an inert,
      !! device-present stand-in — see `compute_h_face_upstream`).
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: use_por
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: u_face(nx + 1, ny, nz), v_face(nx, ny + 1, nz)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz), open_v(nx, ny + 1, nz)
      real(wp), intent(in) :: por_u(nx + 1, ny, nz), por_v(nx, ny + 1, nz)
      real(wp), intent(in) :: dy_cu(nx + 1, ny), dx_cv(nx, ny + 1)
      real(wp), intent(in) :: dy_cu_bt(nx + 1, ny), dx_cv_bt(nx, ny + 1)
      real(wp), intent(inout) :: h_up_x(nx + 1, ny), h_up_y(nx, ny + 1)

      integer :: i, j, k
      real(wp) :: s, h_k, w

      do concurrent(j=1:ny, i=1:nx + 1) local(k, s, h_k, w)
         s = 0.0_wp
         do k = 1, nz
            if (i == 1) then
               h_k = h_layer(1, j, k)
            else if (i == nx + 1) then
               h_k = h_layer(nx, j, k)
            else if (u_face(i, j, k) >= 0.0_wp) then
               h_k = h_layer(i - 1, j, k)
            else
               h_k = h_layer(i, j, k)
            end if
            w = open_u(i, j, k)
            if (use_por) w = w*por_u(i, j, k)
            s = s + h_k*w
         end do
         if (dy_cu_bt(i, j) > 0.0_wp) then
            h_up_x(i, j) = s*(dy_cu(i, j)/dy_cu_bt(i, j))
         else
            h_up_x(i, j) = s
         end if
      end do

      do concurrent(j=1:ny + 1, i=1:nx) local(k, s, h_k, w)
         s = 0.0_wp
         do k = 1, nz
            if (j == 1) then
               h_k = h_layer(i, 1, k)
            else if (j == ny + 1) then
               h_k = h_layer(i, ny, k)
            else if (v_face(i, j, k) >= 0.0_wp) then
               h_k = h_layer(i, j - 1, k)
            else
               h_k = h_layer(i, j, k)
            end if
            w = open_v(i, j, k)
            if (use_por) w = w*por_v(i, j, k)
            s = s + h_k*w
         end do
         if (dx_cv_bt(i, j) > 0.0_wp) then
            h_up_y(i, j) = s*(dx_cv(i, j)/dx_cv_bt(i, j))
         else
            h_up_y(i, j) = s
         end if
      end do
   end subroutine h_face_upstream_open_impl

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

   pure subroutine add_top_drag_into_F_slow(bt_work, td, ms)
      !! Add the ice-shelf top-drag tendency into the already-summed
      !! slow forcing.  Separate from `sum_slow_tendencies_into_F_slow`
      !! (rather than a sixth term in it) for one reason: the top-drag
      !! slot is OPTIONAL all the way down the driver chain, and the sum
      !! above must stay a single unconditional kernel with no `present`
      !! branch inside its `do concurrent`.
      !!
      !! **Why it has to be here at all.**  `F_slow` is depth-meaned into
      !! `F_bt`, the barotropic substep integrates `F_bt`, and
      !! `apply_bt_correction` subtracts `dt*F_bt` back out of the layer
      !! update.  A layer tendency that is applied to the layers but NOT
      !! in `F_slow` is therefore (a) invisible to the fast mode — a
      !! barotropic cavity flow would feel no top friction at all inside
      !! the substep loop — and (b) not subtracted by the correction, so
      !! its damping re-enters the barotropic state one stage late as an
      !! uncorrected residue.  Bottom drag is in the sum for exactly this
      !! reason; the side-wall (channel) drag is NOT, and is the standing
      !! counter-example of the bug this avoids.
      !!
      !! The buffer is added whether or not the explicit apply runs: when
      !! the drag is folded into the vdiff `k = nz` diagonal the layers
      !! get it implicitly AFTER the barotropic correction, while the fast
      !! mode still needs the explicit estimate — the same split the
      !! bottom drag's `implicit_drag` path already takes.
      !!
      !! No-op (and no kernel launch) when the slot is disabled: its
      !! buffers are `(1,1,1)` placeholders then.
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(ocean_top_drag_t), intent(in) :: td
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nu, nv, nx, ny, nz

      if (.not. td%enable) return

      nu = size(bt_work%F_slow_u, 1)
      nv = size(bt_work%F_slow_v, 2)
      nx = size(bt_work%F_slow_v, 1)
      ny = size(bt_work%F_slow_u, 2)
      nz = ms%nz_ml

      do concurrent(k=1:nz, j=1:ny, i=1:nu)
         bt_work%F_slow_u(i, j, k) = bt_work%F_slow_u(i, j, k) + &
                                     td%du_drag%data(i, j, k)
      end do
      do concurrent(k=1:nz, j=1:nv, i=1:nx)
         bt_work%F_slow_v(i, j, k) = bt_work%F_slow_v(i, j, k) + &
                                     td%dv_drag%data(i, j, k)
      end do
   end subroutine add_top_drag_into_F_slow

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
      !! analogue of the PGF double-count `set_fast_forcing_eta_pf`
      !! already guards against by shedding the free-surface term the
      !! slow PGF carries ("√(gH) inflates to √(2gH)").  The Coriolis double-count is what pumps
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
         ! The tripolar fold line is a seam (interior corners), not a wall.
         if (bc_n /= OBC_PERIODIC .and. bc_n /= OBC_TRIPOLAR_FOLD .and. has_n .and. &
             j == grid%nghost + grid%ny_phys + 1) then
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

   pure subroutine set_cor_ref_velocity(grid, bt_work, ms, from_u_av, metrics, n_inner)
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
      !!
      !! **Under `&vcoord_nml zfixed_closed_faces` "the same weights" is
      !! load-bearing and was, for one release, wrong here.**  A closed
      !! layer carries exactly zero velocity (`mask_layer_velocities`)
      !! but a non-zero `h_face`, so with
      !! `φ = Σ_k h_face·open / Σ_k h_face` a FULL-column mean of `u_av`
      !! returns `φ·ū_open`, not `ū_open`.  The fast loop meanwhile
      !! integrates its live `(ζ+f)·v̄ − ∇KE` on `bt_ubt = ū_open`, so
      !! `subtract_fast_cor_ref` would remove `f·φ·v̄` where it must
      !! remove `f·v̄`, leaving `Δa_u = +f·(1−φ)·v̄` forcing EVERY
      !! barotropic substep *proportionally to the barotropic velocity* —
      !! an amplifier, not a seed.  Measured on ISOMIP+ Ocean0 (melt off,
      !! 30 d): `En` `1.295E-06 → 5.711E-08` and barotropic KE ×195
      !! smaller once `metrics` is passed, matching the `ssp_rk2` twin
      !! (whose `.false.` branch below is a plain copy of `bt_ubt`, and so
      !! could never have the defect) to 2.7 %.  `φ` is O(0.5) on a
      !! partial-step face, not O(1 − 1e-4).
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      logical, intent(in) :: from_u_av
         !! `.true.` under `split_scheme = "pred_corr"`.
      type(ocean_metrics_t), intent(in) :: metrics
         !! The closed-face mask carrier, REQUIRED — see the paragraph
         !! above, and `face_depth_mean_u`'s own `metrics` docstring.
         !! Knob off ⇒ the depth means take their original branch and
         !! this is byte-identical.
      integer, intent(in), optional :: n_inner
         !! Barotropic substep count, forwarded to `face_depth_mean_rem_u/v`'s
         !! MOM6 `wt_u` floor when `bt_forcing_visc_rem` is on.  Optional
         !! (defaults to 1) ONLY so call sites that never set
         !! `forcing_visc_rem` (that branch is then never taken) need not
         !! be touched — a real `forcing_visc_rem` run must pass the true
         !! value or the floor's `Instep` is wrong.

      integer :: i, j, nu, nv, nx, ny, n_inner_use
      logical :: use_av

      n_inner_use = 1
      if (present(n_inner)) n_inner_use = n_inner
      use_av = from_u_av
      if (use_av) use_av = allocated(ms%u_av_layer) .and. allocated(ms%v_av_layer)

      if (use_av) then
         if (bt_work%bt_forcing_visc_rem) then
            call face_depth_mean_rem_u(grid, ms%u_av_layer, ms%h_layer, &
                                       bt_work%visc_rem_u, bt_work%cor_ref_u, ms%nz_ml, metrics, &
                                       n_inner_use, bt_work%bt_H_ref, bt_work%frhat_scheme)
            call face_depth_mean_rem_v(grid, ms%v_av_layer, ms%h_layer, &
                                       bt_work%visc_rem_v, bt_work%cor_ref_v, ms%nz_ml, metrics, &
                                       n_inner_use, bt_work%bt_H_ref, bt_work%frhat_scheme)
         else
            call face_depth_mean_u(grid, ms%u_av_layer, ms%h_layer, bt_work%cor_ref_u, &
                                   ms%nz_ml, metrics, bt_work%bt_H_ref, bt_work%frhat_scheme)
            call face_depth_mean_v(grid, ms%v_av_layer, ms%h_layer, bt_work%cor_ref_v, &
                                   ms%nz_ml, metrics, bt_work%bt_H_ref, bt_work%frhat_scheme)
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

   pure subroutine face_depth_mean_u(grid, F_3d, h_layer, F_mean_2d, nz, metrics, href, scheme)
      !! Depth-average a u-face 3D field, weighted by the per-layer face
      !! thickness `frhat_h_face_step` returns (`&ocean_bt_nml
      !! frhat_scheme`; `FRHAT_ARITHMETIC` = the plain mean of the two
      !! abutting cell columns' `h_layer` values, today's default; the
      !! argument is otherwise unchanged).  Writes to a 2D field at the
      !! same u-face shape.  Wall faces (i=1, nx+1) fall back to the
      !! single available cell (`frhat_h_face_step` degenerates there for
      !! either scheme — see its docstring).
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx+1,ny,nz) / (nx,ny+1,nz);
      ! a single (nx,ny,nz) explicit-shape triplet would mis-bound the face axis.
      ! size() is used to derive loop bounds from the actual face dimension.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED.  Under `&vcoord_nml zfixed_closed_faces` every
         !! depth mean in the split chain MUST use the same weights
         !! `derive_bt_from_layers` and `apply_bt_correction` use —
         !! `h_face·open` — or the `dt·F_bt` the fold subtracts back out
         !! is not the quantity the fast loop integrated, and the
         !! difference survives as a permanent per-face bias.
         !!
         !! Not optional, deliberately: this dummy WAS optional, and the
         !! one call site that omitted it (`set_cor_ref_velocity`) turned
         !! the `pred_corr` Coriolis reference into `φ·ū_open` and pumped
         !! the ISOMIP+ Ocean0 barotropic mode by ×22.7 in energy over 30
         !! days.  A caller that has no mask still has an
         !! `ocean_metrics_t` to hand; `use_closed_faces = .false.` ⇒ the
         !! ORIGINAL loop runs, textually unchanged (byte-identical).
      real(wp), intent(in) :: href(:, :)  ! assumed-shape-ok: cell-centred (nx,ny); size() unused, indexed by column
         !! Reference column depth (m), `barotropic_workstate_t%bt_H_ref` —
         !! REQUIRED for the same reason `metrics` is: every call site has
         !! a `bt_work` to hand, and an omitted/mismatched `href` would
         !! silently change `FRHAT_HYBRID`'s answer. Unused under
         !! `FRHAT_ARITHMETIC` but still dereferenced (no branch to skip
         !! it), so it must be a real, fully-sized array on every call,
         !! never a placeholder.
      integer, intent(in) :: scheme
         !! A `FRHAT_*` constant (`rdb_constants`; `&ocean_bt_nml
         !! frhat_scheme`, via `bt_work%frhat_scheme`).
      integer :: i, j, k, nu, ny, nx_cells, il, ir
      real(wp) :: h_face, num, denom, href_l, href_r, e_prev
      logical :: use_open

      nu = size(F_3d, 1)
      ny = size(F_3d, 2)
      nx_cells = grid%nx_total
      use_open = metrics%use_closed_faces

      if (use_open) then
         do concurrent(j=1:ny, i=1:nu) local(k, il, ir, href_l, href_r, e_prev, h_face, num, denom)
            il = max(1, i - 1)
            ir = min(nx_cells, i)
            href_l = href(il, j)
            href_r = href(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(il, j, k), h_layer(ir, j, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               h_face = h_face*metrics%open_u(i, j, k)
               num = num + F_3d(i, j, k)*h_face
               denom = denom + h_face
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:ny, i=1:nu) local(k, il, ir, href_l, href_r, e_prev, h_face, num, denom)
            il = max(1, i - 1)
            ir = min(nx_cells, i)
            href_l = href(il, j)
            href_r = href(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(il, j, k), h_layer(ir, j, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               num = num + F_3d(i, j, k)*h_face
               denom = denom + h_face
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      end if
   end subroutine face_depth_mean_u

   pure subroutine face_depth_mean_rem_u(grid, F_3d, h_layer, rem, F_mean_2d, nz, metrics, n_inner, &
                                         href, scheme)
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
      !!
      !! `rem` is run through MOM6's exact `wt_u` floor before it weights
      !! anything (`MOM_barotropic.F90:1082-1101`): `vr = min(rem, 1)`,
      !! `vr = max(vr, 1 - 0.5·Instep/(vr + subroundoff))`,
      !! `vr = max(vr, 0)`, `Instep = 1/n_inner` — NOT roundabout's old
      !! plain `[0,1]` clamp, which let a near-zero `visc_rem` on a
      !! many-substep column weight the forcing far closer to zero than
      !! MOM6 ever lets it (the floor's whole job is to keep the
      !! `Instep`-th root finite on exactly this kind of thin cell; see
      !! the plan's "Thin-cell floors" risk note).  `ieee_is_finite`-
      !! guarded per the CLAUDE.md NaN-clamp gotcha: a non-finite `rem`
      !! passes through unmasked rather than being laundered to 0 or 1.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx+1,ny,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(in) :: rem(:, :, :)      ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED (see `face_depth_mean_u`): the weight becomes
         !! `h_face·visc_rem·open`.  It must match
         !! `derive_bt_from_layers` and `apply_bt_correction` or the
         !! `dt·F_bt` the fold subtracts is not what the fast loop
         !! integrated.  Note a CLOSED layer's `visc_rem` is ~1, not 0 —
         !! the closed-face vdiff decoupling leaves it uncoupled, so
         !! `visc_rem` alone does NOT stand in for the mask here.
         !! Knob off ⇒ the ORIGINAL loop, byte-identical.
      integer, intent(in) :: n_inner
         !! Barotropic substep count (MOM6 `nstep`); `Instep = 1/n_inner`
         !! in the `wt_u` floor.  `max(n_inner, 1)` guards the unsplit
         !! (`n_inner = 0`) configuration.
      real(wp), intent(in) :: href(:, :)  ! assumed-shape-ok: cell-centred (nx,ny); see face_depth_mean_u
      integer, intent(in) :: scheme
         !! A `FRHAT_*` constant.  See `face_depth_mean_u`.
      integer :: i, j, k, nu, ny, nx_cells, il, ir
      real(wp) :: h_face, wt, num, denom, vr, instep, href_l, href_r, e_prev
      logical :: use_open

      nu = size(F_3d, 1)
      ny = size(F_3d, 2)
      nx_cells = grid%nx_total

      use_open = metrics%use_closed_faces
      instep = 1.0_wp/real(max(n_inner, 1), wp)

      if (use_open) then
         do concurrent(j=1:ny, i=1:nu) &
            local(k, il, ir, href_l, href_r, e_prev, h_face, wt, num, denom, vr)
            il = max(1, i - 1)
            ir = min(nx_cells, i)
            href_l = href(il, j)
            href_r = href(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(il, j, k), h_layer(ir, j, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               if (ieee_is_finite(rem(i, j, k))) then
                  ! Clamp to [0, 1] BEFORE MOM6's floor: MOM6's remnant is
                  ! non-negative by construction, but a round-off -1e-17 here
                  ! would make 1 - 0.5*Instep/vr ~ +5e15 and survive max(., 0).
                  vr = max(min(rem(i, j, k), 1.0_wp), 0.0_wp)
                  vr = max(vr, 1.0_wp - 0.5_wp*instep/(vr + VISC_REM_SUBROUNDOFF))
                  vr = max(vr, 0.0_wp)
               else
                  vr = rem(i, j, k)
               end if
               wt = metrics%open_u(i, j, k)*h_face*vr
               num = num + F_3d(i, j, k)*wt
               denom = denom + wt
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:ny, i=1:nu) &
            local(k, il, ir, href_l, href_r, e_prev, h_face, wt, num, denom, vr)
            il = max(1, i - 1)
            ir = min(nx_cells, i)
            href_l = href(il, j)
            href_r = href(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(il, j, k), h_layer(ir, j, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               if (ieee_is_finite(rem(i, j, k))) then
                  ! Clamp to [0, 1] BEFORE MOM6's floor: MOM6's remnant is
                  ! non-negative by construction, but a round-off -1e-17 here
                  ! would make 1 - 0.5*Instep/vr ~ +5e15 and survive max(., 0).
                  vr = max(min(rem(i, j, k), 1.0_wp), 0.0_wp)
                  vr = max(vr, 1.0_wp - 0.5_wp*instep/(vr + VISC_REM_SUBROUNDOFF))
                  vr = max(vr, 0.0_wp)
               else
                  vr = rem(i, j, k)
               end if
               wt = h_face*vr
               num = num + F_3d(i, j, k)*wt
               denom = denom + wt
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      end if
   end subroutine face_depth_mean_rem_u

   pure subroutine face_depth_mean_rem_v(grid, F_3d, h_layer, rem, F_mean_2d, nz, metrics, n_inner, &
                                         href, scheme)
      !! Symmetric v-face counterpart of `face_depth_mean_rem_u` — same
      !! MOM6 `wt_u` floor, `ieee_is_finite`-guarded the same way.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx,ny+1,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(in) :: rem(:, :, :)      ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED (see `face_depth_mean_u`): the weight becomes
         !! `h_face·visc_rem·open`.  It must match
         !! `derive_bt_from_layers` and `apply_bt_correction` or the
         !! `dt·F_bt` the fold subtracts is not what the fast loop
         !! integrated.  Note a CLOSED layer's `visc_rem` is ~1, not 0 —
         !! the closed-face vdiff decoupling leaves it uncoupled, so
         !! `visc_rem` alone does NOT stand in for the mask here.
         !! Knob off ⇒ the ORIGINAL loop, byte-identical.
      integer, intent(in) :: n_inner
         !! Barotropic substep count (MOM6 `nstep`); see `face_depth_mean_rem_u`.
      real(wp), intent(in) :: href(:, :)  ! assumed-shape-ok: cell-centred (nx,ny); see face_depth_mean_u
      integer, intent(in) :: scheme
         !! A `FRHAT_*` constant.  See `face_depth_mean_u`.
      integer :: i, j, k, nx, nv, ny_cells, jl, jr
      real(wp) :: h_face, wt, num, denom, vr, instep, href_l, href_r, e_prev
      logical :: use_open

      nx = size(F_3d, 1)
      nv = size(F_3d, 2)
      ny_cells = grid%ny_total

      use_open = metrics%use_closed_faces
      instep = 1.0_wp/real(max(n_inner, 1), wp)

      if (use_open) then
         do concurrent(j=1:nv, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, h_face, wt, num, denom, vr)
            jl = max(1, j - 1)
            jr = min(ny_cells, j)
            href_l = href(i, jl)
            href_r = href(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(i, jl, k), h_layer(i, jr, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               if (ieee_is_finite(rem(i, j, k))) then
                  ! Clamp to [0, 1] BEFORE MOM6's floor: MOM6's remnant is
                  ! non-negative by construction, but a round-off -1e-17 here
                  ! would make 1 - 0.5*Instep/vr ~ +5e15 and survive max(., 0).
                  vr = max(min(rem(i, j, k), 1.0_wp), 0.0_wp)
                  vr = max(vr, 1.0_wp - 0.5_wp*instep/(vr + VISC_REM_SUBROUNDOFF))
                  vr = max(vr, 0.0_wp)
               else
                  vr = rem(i, j, k)
               end if
               wt = metrics%open_v(i, j, k)*h_face*vr
               num = num + F_3d(i, j, k)*wt
               denom = denom + wt
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:nv, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, h_face, wt, num, denom, vr)
            jl = max(1, j - 1)
            jr = min(ny_cells, j)
            href_l = href(i, jl)
            href_r = href(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(i, jl, k), h_layer(i, jr, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               if (ieee_is_finite(rem(i, j, k))) then
                  ! Clamp to [0, 1] BEFORE MOM6's floor: MOM6's remnant is
                  ! non-negative by construction, but a round-off -1e-17 here
                  ! would make 1 - 0.5*Instep/vr ~ +5e15 and survive max(., 0).
                  vr = max(min(rem(i, j, k), 1.0_wp), 0.0_wp)
                  vr = max(vr, 1.0_wp - 0.5_wp*instep/(vr + VISC_REM_SUBROUNDOFF))
                  vr = max(vr, 0.0_wp)
               else
                  vr = rem(i, j, k)
               end if
               wt = h_face*vr
               num = num + F_3d(i, j, k)*wt
               denom = denom + wt
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      end if
   end subroutine face_depth_mean_rem_v

   pure subroutine face_depth_mean_v(grid, F_3d, h_layer, F_mean_2d, nz, metrics, href, scheme)
      !! Symmetric v-face counterpart of `face_depth_mean_u`.
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: face arrays have shape (nx,ny+1,nz); a single (nx,ny,nz)
      ! explicit-shape triplet would mis-bound the face axis.
      real(wp), intent(in) :: F_3d(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      real(wp), intent(out) :: F_mean_2d(:, :)  ! assumed-shape-ok: face-sized array; size() derives loop bounds
      integer, intent(in) :: nz
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED.  See `face_depth_mean_u` for the argument and for
         !! why it is not optional; knob off ⇒ the ORIGINAL loop,
         !! byte-identical.
      real(wp), intent(in) :: href(:, :)  ! assumed-shape-ok: cell-centred (nx,ny); see face_depth_mean_u
      integer, intent(in) :: scheme
         !! A `FRHAT_*` constant.  See `face_depth_mean_u`.
      integer :: i, j, k, nx, nv, ny_cells, jl, jr
      real(wp) :: h_face, num, denom, href_l, href_r, e_prev
      logical :: use_open

      nx = size(F_3d, 1)
      nv = size(F_3d, 2)
      ny_cells = grid%ny_total
      use_open = metrics%use_closed_faces

      if (use_open) then
         do concurrent(j=1:nv, i=1:nx) local(k, jl, jr, href_l, href_r, e_prev, h_face, num, denom)
            jl = max(1, j - 1)
            jr = min(ny_cells, j)
            href_l = href(i, jl)
            href_r = href(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(i, jl, k), h_layer(i, jr, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               h_face = h_face*metrics%open_v(i, j, k)
               num = num + F_3d(i, j, k)*h_face
               denom = denom + h_face
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:nv, i=1:nx) local(k, jl, jr, href_l, href_r, e_prev, h_face, num, denom)
            jl = max(1, j - 1)
            jr = min(ny_cells, j)
            href_l = href(i, jl)
            href_r = href(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            num = 0.0_wp
            denom = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(h_layer(i, jl, k), h_layer(i, jr, k), &
                                      href_l, href_r, scheme, e_prev, h_face)
               num = num + F_3d(i, j, k)*h_face
               denom = denom + h_face
            end do
            if (denom > 0.0_wp) then
               F_mean_2d(i, j) = num/denom
            else
               F_mean_2d(i, j) = 0.0_wp
            end if
         end do
      end if
   end subroutine face_depth_mean_v

   pure subroutine apply_bt_correction(bt_work, ms, dt, metrics, skip_h_rescale, &
                                       grid, use_bc_pgf, use_visc_rem, scale, n_nonfin, n_inner)
      !! Replace the bt mode in the per-layer face velocities with the
      !! barotropic-substep end-step value, adding `Δu·wt_k` to every layer,
      !! `Δu = u_bt_end − u_bt_at_n − dt·F_bt_u` (same for v). Split-explicit
      !! convention (Hallberg 2009): momentum uses the END-of-step barotropic
      !! velocity; layer continuity earlier used the time-mean transports.
      !! Both legs of the corrector use the same end-step anchor (mismatched
      !! anchors overshoot the gravity-wave phase speed). Also rescales
      !! `h_layer` uniformly so the column total matches `H_ref + η_end`.
      !! `hTr` is deliberately NOT rescaled (would break exact tracer mass
      !! conservation; T = hTr/h drifts by O((η_end−η*_slow)/H) per step).
      !!
      !! ### The fold weight
      !!
      !! `wt_k = open_k·vr_k / ⟨vr⟩_h`,  `⟨vr⟩_h = Σ_k h_o·vr_k / Σ_k h_o`,
      !! `h_o = h_face·open`, so `Σ_k h_o·wt_k = Σ_k h_o` and the OPEN-column
      !! depth mean moves by exactly `Δu`.  `vr_k = visc_rem(k)` with
      !! `use_visc_rem` (MOM6 `visc_rem_u`; the drag-aware weighting),
      !! else 1, and `open ≡ 1` unless `metrics%use_closed_faces`.  With
      !! neither it is the uniform fold, `wt ≡ 1` — MOM6's own barotropic
      !! acceleration, `accel_layer_u(I,j,k) = u_accel_bt(I,j)`
      !! (`MOM_barotropic.F90`), the same increment in every layer.
      !!
      !! The weight does NOT carry `h`.  An earlier opt-in h-weighted form
      !! (`wt = h_face/⟨h⟩_h`, `&ocean_bt_nml correction_h_weighted`, now
      !! retired and refused by `validate_config`) adds, beyond the
      !! barotropic `ΔKE`, a positive-definite source `½Δ²·H·(κ−1)`,
      !! `κ = Σh³Σh/(Σh²)² ≥ 1`, plus a shear feedback `Δ·Σ h·u′·(wt−1)`
      !! that the fold keeps feeding: on a stretched `z_fixed` stack
      !! (`κ−1 ≈ 0.2`) the feedback ran 15-28x the source and grew the
      !! 1-degree Southern Ocean and the coastal-noise box ~5-8x until
      !! they went non-finite.  `frhatu·visc_rem` is MOM6's AVERAGING
      !! weight (BT_force, ubt), never its distribution.
      !!
      !! `skip_h_rescale` — disable the h-rescale (Lagrangian vcoord, where slow
      !!   continuity's Σh_layer is authoritative and the ALE remap relayers).
      !!   Default `.false.`.
      !! `use_bc_pgf` — add the per-layer baroclinic-PGF retro-correction
      !!   Δu_bc = -dt·((pbce(R,k)-gtot_W(R))·e_anom(R) -
      !!   (pbce(L,k)-gtot_E(L))·e_anom(L))/dx. Depth-mean zero by construction,
      !!   so the mass-flux invariant survives. Requires `grid` present
      !!   and pbce/gtot_*/e_anom populated by the caller. No-op when omitted.
      type(barotropic_workstate_t), intent(in) :: bt_work
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      type(ocean_metrics_t), intent(in) :: metrics
         !! Curvilinear horizontal metrics — read by the bc-PGF
         !! retro-correction (`use_bc_pgf = .true.` divides the e_anom
         !! gradient by `idxCu`/`idyCv`), and the carrier of the z-level
         !! closed-face mask.
         !!
         !! **REQUIRED**, like the same argument of
         !! `derive_bt_from_layers`, `face_depth_mean_u/v` and
         !! `set_cor_ref_velocity`: an OPTIONAL `metrics` that silently
         !! selects the full-column fold when omitted is the defect class
         !! that shipped the `set_cor_ref_velocity` bug.  A caller that
         !! wants the full-column fold passes a metrics object whose
         !! `use_closed_faces` latch is `.false.` (the default).
         !!
         !! When `metrics%use_closed_faces` is set a CLOSED layer gets
         !! `wt = 0` and receives nothing, BY CONSTRUCTION — the weight
         !! carries `open` and `⟨vr⟩_h` is the OPEN-column mean, so the
         !! OPEN-column depth mean (the column `derive_bt_from_layers` and
         !! `face_depth_mean_*` weight by once the knob is on) shifts by
         !! exactly `Δu`.  This replaces the spike's fold-then-mask: there
         !! `Δu` was added uniformly to every layer and
         !! `mask_layer_velocities` removed it again from the closed ones,
         !! so the layer depth mean fell short of `ubt_end` by
         !! `Δu·(Σ_closed h)/(Σ_k h)` — preserved by CANCELLATION rather
         !! than by construction.
      logical, intent(in), optional :: skip_h_rescale
      type(hgrid_t), intent(in), optional :: grid
      logical, intent(in), optional :: use_bc_pgf
      logical, intent(in), optional :: use_visc_rem
         !! Weight the fold by `visc_rem/⟨visc_rem⟩_h`.  RETIRED as a live
         !! namelist path (`&ocean_bt_nml correction_visc_rem` is
         !! fail-loud at configure, D1 follow-up — MOM6's `accel_layer_u`
         !! never weights this fold) but the kernel dispatch stays, for
         !! its own direct unit tests.  Default `.false.` ⇒ the uniform
         !! fold, which `visc_rem_chain` also uses.
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
      integer, intent(in), optional :: n_inner
         !! Barotropic substeps per outer step.  REQUIRED when
         !! `bt_work%bt_rescale_strong_drag` is on (PR-2, MOM6
         !! `RESCALE_STRONG_DRAG`, `MOM_barotropic.F90:1989-1997`):
         !! `bt_strong_drag`'s rational-approximation `bt_rem` does not
         !! satisfy `bt_rem**n_in == av_rem` exactly (unlike the plain
         !! power form, which does by construction), so the Δu/Δv
         !! correction is rescaled by `min(bt_rem**n_in/av_rem, 1.0)`
         !! before being distributed into the layers — keeping the
         !! correction consistent with the TRUE depth-mean remnant.
         !! Ignored when `bt_rescale_strong_drag` is off.

      integer :: i, j, k, nu, nv, nx, ny, nz, nfin
      integer :: n_in
         !! Local copy of the optional `n_inner`.  Never read an absent
         !! optional inside `do concurrent`: ifx lowers the loop to an
         !! OpenMP region and loads every captured scalar at region entry,
         !! so a null reference segfaults even on a branch never taken.
      real(wp) :: delta_u, delta_v, total_h_old, total_h_new, ratio
      real(wp) :: du_scale
      real(wp) :: h_face, sum_h, sum_hvr, vr_bar, wt, vr_k
      real(wp) :: du_bc, dv_bc
      real(wp) :: href_l, href_r, e_prev
      integer :: il, ir, jl, jr, frhat_scheme
      logical :: do_rescale, do_bc_pgf, do_visc_rem, do_open, do_bt_rescale

      frhat_scheme = bt_work%frhat_scheme
      do_rescale = .true.
      if (present(skip_h_rescale)) do_rescale = .not. skip_h_rescale
      do_bc_pgf = .false.
      if (present(use_bc_pgf)) do_bc_pgf = use_bc_pgf
      do_visc_rem = .false.
      if (present(use_visc_rem)) do_visc_rem = use_visc_rem
      du_scale = 1.0_wp
      if (present(scale)) du_scale = scale
      do_open = metrics%use_closed_faces
      n_in = 1
      if (present(n_inner)) n_in = n_inner
      do_bt_rescale = bt_work%bt_rescale_strong_drag
      if (do_bt_rescale .and. .not. present(n_inner)) then
         error stop "apply_bt_correction: rescale_strong_drag requires n_inner"
      end if
      if (do_bc_pgf .and. .not. present(grid)) then
         error stop "apply_bt_correction: use_bc_pgf=.true. requires grid"
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

      if (do_open) then
         ! OPEN-LAYER fold (`&vcoord_nml zfixed_closed_faces`):
         ! `wt = open·vr/⟨vr⟩_h` over `h_o = h_face·open`.  A CLOSED layer
         ! receives nothing.  Without visc_rem `wt = open` — written as
         ! exactly that (not `open·1/1`) so the no-visc_rem closed-face
         ! answer is the one this branch always gave.
         do concurrent(j=1:ny, i=1:nu) &
            local(k, il, ir, href_l, href_r, e_prev, delta_u, sum_h, sum_hvr, h_face, vr_bar, wt, vr_k)
            delta_u = du_scale*(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) - dt*bt_work%F_bt_u(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_u(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_u(i, j))) then
                  delta_u = delta_u*min(bt_work%bt_rem_u(i, j)**n_in/bt_work%av_rem_u(i, j), 1.0_wp)
               end if
            end if
            il = max(1, i - 1)
            ir = min(nu - 1, i)
            href_l = bt_work%bt_H_ref(il, j)
            href_r = bt_work%bt_H_ref(ir, j)
            e_prev = -0.5_wp*(href_l + href_r)
            sum_h = 0.0_wp
            sum_hvr = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(ms%h_layer(il, j, k), ms%h_layer(ir, j, k), &
                                      href_l, href_r, frhat_scheme, e_prev, h_face)
               h_face = h_face*metrics%open_u(i, j, k)
               if (do_visc_rem) then
                  vr_k = bt_work%visc_rem_u(i, j, k)
               else
                  vr_k = 1.0_wp
               end if
               sum_h = sum_h + h_face
               sum_hvr = sum_hvr + h_face*vr_k
            end do
            if (sum_h > 0.0_wp .and. ieee_is_finite(delta_u)) then
               vr_bar = 1.0_wp
               if (do_visc_rem .and. sum_hvr > 0.0_wp) vr_bar = sum_hvr/sum_h
               do k = 1, nz
                  if (do_visc_rem) then
                     wt = metrics%open_u(i, j, k)*bt_work%visc_rem_u(i, j, k)/vr_bar
                  else
                     wt = metrics%open_u(i, j, k)
                  end if
                  ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u*wt
               end do
            end if
         end do
         do concurrent(j=1:nv, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, delta_v, sum_h, sum_hvr, h_face, vr_bar, wt, vr_k)
            delta_v = du_scale*(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) - dt*bt_work%F_bt_v(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_v(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_v(i, j))) then
                  delta_v = delta_v*min(bt_work%bt_rem_v(i, j)**n_in/bt_work%av_rem_v(i, j), 1.0_wp)
               end if
            end if
            jl = max(1, j - 1)
            jr = min(nv - 1, j)
            href_l = bt_work%bt_H_ref(i, jl)
            href_r = bt_work%bt_H_ref(i, jr)
            e_prev = -0.5_wp*(href_l + href_r)
            sum_h = 0.0_wp
            sum_hvr = 0.0_wp
            do k = 1, nz
               call frhat_h_face_step(ms%h_layer(i, jl, k), ms%h_layer(i, jr, k), &
                                      href_l, href_r, frhat_scheme, e_prev, h_face)
               h_face = h_face*metrics%open_v(i, j, k)
               if (do_visc_rem) then
                  vr_k = bt_work%visc_rem_v(i, j, k)
               else
                  vr_k = 1.0_wp
               end if
               sum_h = sum_h + h_face
               sum_hvr = sum_hvr + h_face*vr_k
            end do
            if (sum_h > 0.0_wp .and. ieee_is_finite(delta_v)) then
               vr_bar = 1.0_wp
               if (do_visc_rem .and. sum_hvr > 0.0_wp) vr_bar = sum_hvr/sum_h
               do k = 1, nz
                  if (do_visc_rem) then
                     wt = metrics%open_v(i, j, k)*bt_work%visc_rem_v(i, j, k)/vr_bar
                  else
                     wt = metrics%open_v(i, j, k)
                  end if
                  ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v*wt
               end do
            end if
         end do
      else if (.not. do_visc_rem) then
         ! Uniform Δu distribution — every layer gets the same Δu.  The finite
         ! guard skips a face whose Δ is non-finite (Inf/NaN from a blown-up BT
         ! loop) so the fold never mints NaN into the layer velocity.
         do concurrent(k=1:nz, j=1:ny, i=1:nu) local(delta_u)
            delta_u = du_scale*(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) - dt*bt_work%F_bt_u(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_u(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_u(i, j))) then
                  delta_u = delta_u*min(bt_work%bt_rem_u(i, j)**n_in/bt_work%av_rem_u(i, j), 1.0_wp)
               end if
            end if
            if (ieee_is_finite(delta_u)) then
               ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u
            end if
         end do
         do concurrent(k=1:nz, j=1:nv, i=1:nx) local(delta_v)
            delta_v = du_scale*(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) - dt*bt_work%F_bt_v(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_v(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_v(i, j))) then
                  delta_v = delta_v*min(bt_work%bt_rem_v(i, j)**n_in/bt_work%av_rem_v(i, j), 1.0_wp)
               end if
            end if
            if (ieee_is_finite(delta_v)) then
               ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v
            end if
         end do
      else
         ! visc_rem-weighted fold on the full column: `wt = vr/⟨vr⟩_h`.
         ! A column with no thickness, or with `Σ h·vr = 0` (every layer
         ! fully damped), takes the uniform increment: `wt ≡ 1` is the
         ! only weight with the right depth mean when `⟨vr⟩_h` is
         ! undefined.
         do concurrent(j=1:ny, i=1:nu) &
            local(k, il, ir, href_l, href_r, e_prev, delta_u, sum_h, sum_hvr, h_face, vr_bar, wt)
            delta_u = du_scale*(bt_work%bt_ubt_end(i, j) - bt_work%ubt_at_n(i, j) - dt*bt_work%F_bt_u(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_u(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_u(i, j))) then
                  delta_u = delta_u*min(bt_work%bt_rem_u(i, j)**n_in/bt_work%av_rem_u(i, j), 1.0_wp)
               end if
            end if
            if (ieee_is_finite(delta_u)) then
               il = max(1, i - 1)
               ir = min(nu - 1, i)
               href_l = bt_work%bt_H_ref(il, j)
               href_r = bt_work%bt_H_ref(ir, j)
               e_prev = -0.5_wp*(href_l + href_r)
               sum_h = 0.0_wp
               sum_hvr = 0.0_wp
               do k = 1, nz
                  call frhat_h_face_step(ms%h_layer(il, j, k), ms%h_layer(ir, j, k), &
                                         href_l, href_r, frhat_scheme, e_prev, h_face)
                  sum_h = sum_h + h_face
                  sum_hvr = sum_hvr + h_face*bt_work%visc_rem_u(i, j, k)
               end do
               if (sum_h > 0.0_wp .and. sum_hvr > 0.0_wp) then
                  vr_bar = sum_hvr/sum_h
                  do k = 1, nz
                     wt = bt_work%visc_rem_u(i, j, k)/vr_bar
                     ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u*wt
                  end do
               else
                  do k = 1, nz
                     ms%u_face_x_layer(i, j, k) = ms%u_face_x_layer(i, j, k) + delta_u
                  end do
               end if
            end if
         end do
         do concurrent(j=1:nv, i=1:nx) &
            local(k, jl, jr, href_l, href_r, e_prev, delta_v, sum_h, sum_hvr, h_face, vr_bar, wt)
            delta_v = du_scale*(bt_work%bt_vbt_end(i, j) - bt_work%vbt_at_n(i, j) - dt*bt_work%F_bt_v(i, j))
            if (do_bt_rescale) then
               if (bt_work%av_rem_v(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_v(i, j))) then
                  delta_v = delta_v*min(bt_work%bt_rem_v(i, j)**n_in/bt_work%av_rem_v(i, j), 1.0_wp)
               end if
            end if
            if (ieee_is_finite(delta_v)) then
               jl = max(1, j - 1)
               jr = min(nv - 1, j)
               href_l = bt_work%bt_H_ref(i, jl)
               href_r = bt_work%bt_H_ref(i, jr)
               e_prev = -0.5_wp*(href_l + href_r)
               sum_h = 0.0_wp
               sum_hvr = 0.0_wp
               do k = 1, nz
                  call frhat_h_face_step(ms%h_layer(i, jl, k), ms%h_layer(i, jr, k), &
                                         href_l, href_r, frhat_scheme, e_prev, h_face)
                  sum_h = sum_h + h_face
                  sum_hvr = sum_hvr + h_face*bt_work%visc_rem_v(i, j, k)
               end do
               if (sum_h > 0.0_wp .and. sum_hvr > 0.0_wp) then
                  vr_bar = sum_hvr/sum_h
                  do k = 1, nz
                     wt = bt_work%visc_rem_v(i, j, k)/vr_bar
                     ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v*wt
                  end do
               else
                  do k = 1, nz
                     ms%v_face_y_layer(i, j, k) = ms%v_face_y_layer(i, j, k) + delta_v
                  end do
               end if
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

   pure function pgf_free_surface_gravity(pgf) result(g_pf)
      !! The gravity of the free-surface term the slow layer PGF CARRIES,
      !! i.e. `−∂⟨PGF⟩/∂(∇η)` for a uniform-density column (m/s²):
      !!
      !! * MONT, FV_LITE, FV_WRIGHT are built from the FREE SURFACE down
      !!   (`M(nz) ≡ 0`, `p_edge(nz+1) = 0`, surface-relative `z`), so they
      !!   carry NO `−g·∇η` at all — 0.  Their layer PGF is purely
      !!   baroclinic.
      !! * FV_MOM6 closes its anomaly stack with `pa(nz+1) = ρ_ref·g·η_geo`
      !!   (MOM6 `PressureForce_FV`) and the MOM6 `GFS_scale` correction
      !!   removes `(1 − gfs_scale)·g·ρ_surf/ρ₀·∇η` — so, with the
      !!   surface density at `ρ_ref`, `gfs_scale·g·ρ_ref/ρ₀`.
      !! * GPRIME's top layer is `−g_FS·∇η` by construction — `g_FS`.
      !!
      !! This is what the barotropic substep's own `−g_bt·∇η` duplicates,
      !! and therefore the only part of `⟨PGF⟩` the fast forcing may shed
      !! (`set_fast_forcing_eta_pf`).
      type(ocean_pressure_force_t), intent(in) :: pgf
      real(wp) :: g_pf
      select case (pgf%variant)
      case (OPGF_VARIANT_FV_MOM6)
         g_pf = pgf%gfs_scale*GRAVITY*pgf%rho_ref/pgf%rho0
      case (OPGF_VARIANT_GPRIME)
         g_pf = pgf%gprime_gfs
      case default
         g_pf = 0.0_wp
      end select
   end function pgf_free_surface_gravity

   pure subroutine set_fast_forcing_eta_pf(grid, metrics, bt_work, nx, ny, g_pf, eta_seam, use_seam)
      !! The barotropic substep's frozen forcing under the MOM6 split
      !! (`&ocean_bt_nml bc_pgf_forcing`, default on):
      !!
      !!     F_bt_u_fast = F_bt_u + g_pf·(η_PF(i) − η_PF(i−1))·idxCu
      !!
      !! `F_bt` holds the depth mean of the FULL slow layer PGF.  The
      !! substep integrates `−g_bt·∇η` live, so the one thing the forcing
      !! must shed is the free-surface term the slow PGF itself carries,
      !! `−g_pf·∇η_PF` (`pgf_free_surface_gravity`: 0 for the
      !! surface-relative MONT/FV_LITE/FV_WRIGHT forms, ≈ `g_bt` for
      !! FV_MOM6/GPRIME), evaluated at `η_PF = bt_eta` at stage entry —
      !! the free surface the slow PGF of this stage was built on
      !! (`derive_bt_from_layers` fills it from the same `h_layer`).  The
      !! substep then sees `⟨PGF_bc⟩ − g_bt·∇η`: the depth-mean BAROCLINIC
      !! pressure gradient plus its own free surface.  MOM6: `BT_force`
      !! (Σ wt·bc_accel, PFu included) with `btloop_find_PF`'s
      !! `−gtot·∇(η − eta_PF)`.
      !!
      !! The legacy split (`F_bt − ⟨PGF⟩`) shed the WHOLE depth-mean PGF —
      !! its baroclinic part too, which is the bottom-pressure gradient of
      !! a sloping density field (JEBAR).  Because `apply_bt_correction`
      !! subtracts `dt·F_bt` from every layer, the layers lost it as well:
      !! the depth mean ended every stage at `u_bt^end`, which never felt
      !! it.
      !!
      !! `eta_seam`/`use_seam`: when the surface-pressure load ALSO enters
      !! the slow PGF (`&ocean_pgf_nml p_top_in_bc` with the psurf seam on),
      !! `⟨PGF⟩` carries `−∇p_surf/ρ₀ = +g·∇η_ib` and the substep carries
      !! it again through `eta_forcing`; shedding `g_pf·∇η_ib` here counts
      !! it once.  Otherwise pass any full-size array and `.false.`.
      !!
      !! Array-edge faces (`i = 1`, `nx + 1`; `j = 1`, `ny + 1`) have no
      !! η on one side; the substep never reads their forcing (they are
      !! overwritten by the boundary dispatch / halo), so they take `F_bt`.
      integer, intent(in) :: nx, ny
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      real(wp), intent(in) :: g_pf
      real(wp), intent(in) :: eta_seam(nx, ny)
      logical, intent(in) :: use_seam
      integer :: i, j, nu, nv
      real(wp) :: d_eta

      if (.false.) nu = grid%nx_total
      nu = nx + 1
      nv = ny + 1
      do concurrent(j=1:ny, i=1:nu) local(d_eta)
         if (i == 1 .or. i == nu) then
            bt_work%F_bt_u_fast(i, j) = bt_work%F_bt_u(i, j)
         else
            d_eta = bt_work%bt_eta(i, j) - bt_work%bt_eta(i - 1, j)
            if (use_seam) d_eta = d_eta - (eta_seam(i, j) - eta_seam(i - 1, j))
            bt_work%F_bt_u_fast(i, j) = bt_work%F_bt_u(i, j) + g_pf*d_eta*metrics%idxCu(i, j)
         end if
      end do
      do concurrent(j=1:nv, i=1:nx) local(d_eta)
         if (j == 1 .or. j == nv) then
            bt_work%F_bt_v_fast(i, j) = bt_work%F_bt_v(i, j)
         else
            d_eta = bt_work%bt_eta(i, j) - bt_work%bt_eta(i, j - 1)
            if (use_seam) d_eta = d_eta - (eta_seam(i, j) - eta_seam(i, j - 1))
            bt_work%F_bt_v_fast(i, j) = bt_work%F_bt_v(i, j) + g_pf*d_eta*metrics%idyCv(i, j)
         end if
      end do
   end subroutine set_fast_forcing_eta_pf

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
      !! (other variants don't fill e_face).  `validate_config` and
      !! `configure_ocean_pgf` refuse `correction_bc_pgf` with any other
      !! form, so the `error stop` below is a backstop for direct callers.
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

   pure subroutine compute_bt_rem(grid, bt_work, ms, metrics, r_linear, hbbl, dt_inner)
      !! Per-face multiplicative damping factor for the BT-substep velocity
      !! update (linear-drag branch):
      !!     bt_rem_face = Htot_face / (Htot_face + r·hbbl·dt_inner)
      !! applied as ubt_new = bt_rem_u·(ubt_old + dt_inner·forces) each inner
      !! step. When the `bt_substep_drag` knob is off this must NOT be called and
      !! the workspace stays at 1 (no-op, bit-identical).
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED (see `derive_bt_from_layers`).  Under
         !! `&vcoord_nml zfixed_closed_faces` `Htot_face` is the OPEN-column
         !! face depth `Σ_k h_face·open` — the column `ubt` is the mean of
         !! and the one the fast loop transports on — so the bed stress
         !! `r·hbbl·ubt` is spread over the depth that actually carries the
         !! barotropic momentum.  The full-column depth would under-damp a
         !! partially closed face by its closed fraction.  `.false.` (the
         !! default) ⇒ the original loops, textually unchanged; `open_*`
         !! is never named.
      real(wp), intent(in) :: r_linear   !! Linear drag rate at the bed (1/s)
      real(wp), intent(in) :: hbbl       !! BBL thickness over which drag acts (m)
      real(wp), intent(in) :: dt_inner   !! BT-substep dt (s)

      integer :: i, j, k, nx, ny, nz
      real(wp) :: htot_face, drag_dt

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      drag_dt = r_linear*hbbl*dt_inner   ! product is in metres

      if (metrics%use_closed_faces) then
         call bt_rem_open_impl(nx, ny, nz, ms%h_layer, metrics%open_u, &
                               metrics%open_v, drag_dt, &
                               bt_work%bt_rem_u, bt_work%bt_rem_v)
         return
      end if

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

   pure subroutine compute_bt_rem_from_visc_rem(grid, bt_work, ms, metrics, n_inner)
      !! PR-2 (bt-rem-from-av-rem): build `bt_rem_u/v` from the SAME
      !! viscous remnant the layered momentum solve uses, MOM6
      !! `MOM_barotropic.F90:1553-1580`.  Dispatched the same way as
      !! `compute_bt_rem` — a RESETTER, mutually exclusive at configure
      !! with `bt_substep_drag` (D2, double-counted bed drag) and with
      !! `bt_halo > 0` (`validate_config`) — so this and `compute_bt_rem`/
      !! `reset_bt_rem` never both run for the same stage; `src/core/
      !! ocean/README.md`'s "exactly one resets, everything else
      !! MULTIPLIES" contract gets this as its third resetter.
      !!
      !! Two steps:
      !!
      !! 1. `av_rem_u/v := Σ_k frhat_k·visc_rem_k`, `frhat_k` the PLAIN
      !!    face-thickness fraction (`h_face_k / Σ_k h_face_k`) —
      !!    **not** `face_depth_mean_rem_u`'s `wt_u = h_face·visc_rem`
      !!    weighting (that one is MOM6's FORCING weight, `forcing_
      !!    visc_rem`/PR-3 scope; this is the plain depth mean MOM6 calls
      !!    `frhatu`). `frhat_k/Σ_k h_face_k` is EXACTLY
      !!    `face_depth_mean_u`'s own weight (num = Σ F·h_face, denom =
      !!    Σ h_face), so `av_rem_u = face_depth_mean_u(visc_rem_u,
      !!    h_layer)` — no separate kernel needed; this reuses the SAME
      !!    `h_face`/`metrics%open_u` branches `derive_bt_from_layers`
      !!    builds `ubt` with (the `metrics` REQUIRED-argument contract:
      !!    see that routine's docstring), so `av_rem` is the depth mean
      !!    over the SAME column the fast loop actually transports on.
      !!    `face_depth_mean_u` already returns `0` on a dry/fully-closed
      !!    face (`denom <= 0`), which is exactly the MOM6 "av_rem = 0 on
      !!    a massless column" edge case.
      !!
      !!    NOTE this `frhat` is roundabout's own: the two-abutting-cell
      !!    arithmetic-mean `h_face` `face_depth_mean_u`/`derive_bt_from_
      !!    layers`/`apply_bt_correction` already share, which is what
      !!    SELF-CONSISTENCY across the BT chain requires here — not
      !!    necessarily MOM6's own `frhatu`, which comes from `BT_cont`'s
      !!    face thicknesses (a different, flux-bounded construction).
      !!    Auditing that parity (or documenting the deliberate
      !!    divergence) is PR-3 scope, not this one.
      !!
      !! 2. `bt_rem = av_rem**(1/n_inner)` where `av_rem > 0` (MOM6
      !!    `Instep = 1/nstep`), else `0` — no `max(..., eps)` floor
      !!    substitute (CLAUDE.md: the thin-cell floor is the `av_rem >
      !!    0` MASK itself, ported exactly).  `bt_strong_drag` (MOM6
      !!    `BT_STRONG_DRAG`) swaps in the rational approximation
      !!    `n_inner·av_rem/(1+(n_inner-1)·av_rem)` instead.  Land/dry
      !!    faces are left to the existing `mask_bt_rem` call that always
      !!    runs last in the dispatch (same posture as `compute_bt_rem`,
      !!    which also does not self-mask) — MOM6's own `mask2dCu`
      !!    multiply is therefore redundant with, not additional to, that
      !!    final mask pass.
      !!
      !! Ghosts: both steps run over the FULL face extent (`size(...,1)`)
      !! including ghost columns/rows, matching `face_depth_mean_u`'s own
      !! convention — `visc_rem_u/v`'s ghosts are halo-valid after PR-1's
      !! `visc_rem_halo_refresh`, so `av_rem`/`bt_rem` are correct on
      !! every face the substep loop reads, not just the owned interior.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(inout) :: bt_work
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED — see `derive_bt_from_layers`/`face_depth_mean_u`.
      integer, intent(in) :: n_inner
         !! Barotropic substeps per outer step (MOM6 `nstep`; must be
         !! >= 1 — `auto_n_inner`/the namelist floor already enforce
         !! that).  `Instep = 1/n_inner`.

      integer :: i, j, nu, nv, ny_u, nx_v
      real(wp) :: instep
      real(wp) :: rn

      call face_depth_mean_u(grid, bt_work%visc_rem_u, ms%h_layer, bt_work%av_rem_u, &
                             ms%nz_ml, metrics, bt_work%bt_H_ref, bt_work%frhat_scheme)
      call face_depth_mean_v(grid, bt_work%visc_rem_v, ms%h_layer, bt_work%av_rem_v, &
                             ms%nz_ml, metrics, bt_work%bt_H_ref, bt_work%frhat_scheme)

      nu = size(bt_work%av_rem_u, 1)
      ny_u = size(bt_work%av_rem_u, 2)
      nx_v = size(bt_work%av_rem_v, 1)
      nv = size(bt_work%av_rem_v, 2)
      instep = 1.0_wp/real(n_inner, wp)
      rn = real(n_inner, wp)

      if (bt_work%bt_strong_drag) then
         do concurrent(j=1:ny_u, i=1:nu)
            bt_work%bt_rem_u(i, j) = 0.0_wp
            if (bt_work%av_rem_u(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_u(i, j))) then
               bt_work%bt_rem_u(i, j) = (rn*bt_work%av_rem_u(i, j))/ &
                                        (1.0_wp + (rn - 1.0_wp)*bt_work%av_rem_u(i, j))
            end if
         end do
         do concurrent(j=1:nv, i=1:nx_v)
            bt_work%bt_rem_v(i, j) = 0.0_wp
            if (bt_work%av_rem_v(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_v(i, j))) then
               bt_work%bt_rem_v(i, j) = (rn*bt_work%av_rem_v(i, j))/ &
                                        (1.0_wp + (rn - 1.0_wp)*bt_work%av_rem_v(i, j))
            end if
         end do
      else
         do concurrent(j=1:ny_u, i=1:nu)
            bt_work%bt_rem_u(i, j) = 0.0_wp
            if (bt_work%av_rem_u(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_u(i, j))) then
               bt_work%bt_rem_u(i, j) = bt_work%av_rem_u(i, j)**instep
            end if
         end do
         do concurrent(j=1:nv, i=1:nx_v)
            bt_work%bt_rem_v(i, j) = 0.0_wp
            if (bt_work%av_rem_v(i, j) > 0.0_wp .and. ieee_is_finite(bt_work%av_rem_v(i, j))) then
               bt_work%bt_rem_v(i, j) = bt_work%av_rem_v(i, j)**instep
            end if
         end do
      end if
   end subroutine compute_bt_rem_from_visc_rem

   pure subroutine bt_rem_open_impl(nx, ny, nz, h_layer, open_u, open_v, drag_dt, &
                                    bt_rem_u, bt_rem_v)
      !! `compute_bt_rem` under `&vcoord_nml zfixed_closed_faces`: the
      !! same `H/(H + r·hbbl·dt_inner)` with `H = Σ_k h_face·open` (the
      !! OPEN-column centred face depth, the weight `face_depth_mean_*`
      !! uses).  A face whose every layer is closed has `H = 0` and keeps
      !! `bt_rem = 1`, exactly like a dry face on the original path (it
      !! carries no barotropic transport: `dy_cu_bt = 0` there).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz), open_v(nx, ny + 1, nz)
      real(wp), intent(in) :: drag_dt
         !! `r_linear·hbbl·dt_inner` (m).
      real(wp), intent(inout) :: bt_rem_u(nx + 1, ny), bt_rem_v(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: htot_face

      do concurrent(j=1:ny, i=2:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + &
                        0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))*open_u(i, j, k)
         end do
         if (htot_face > 0.0_wp) then
            bt_rem_u(i, j) = htot_face/(htot_face + drag_dt)
         else
            bt_rem_u(i, j) = 1.0_wp
         end if
      end do
      do concurrent(j=1:ny)
         bt_rem_u(1, j) = 1.0_wp
         bt_rem_u(nx + 1, j) = 1.0_wp
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + &
                        0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))*open_v(i, j, k)
         end do
         if (htot_face > 0.0_wp) then
            bt_rem_v(i, j) = htot_face/(htot_face + drag_dt)
         else
            bt_rem_v(i, j) = 1.0_wp
         end if
      end do
      do concurrent(i=1:nx)
         bt_rem_v(i, 1) = 1.0_wp
         bt_rem_v(i, ny + 1) = 1.0_wp
      end do
   end subroutine bt_rem_open_impl

   pure subroutine compute_bt_rem_wave_drag(grid, bt_work, ms, metrics, dt_inner)
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
      type(ocean_metrics_t), intent(in) :: metrics
         !! REQUIRED.  Under `&vcoord_nml zfixed_closed_faces`
         !! `Htot_face = Σ_k h_face·open` — the same OPEN-column depth
         !! `compute_bt_rem` uses there, so the "identical `Htot_face`"
         !! contract above holds on both paths: the piston velocity
         !! `r_H` damps the column that carries `ubt`, at the rate
         !! `r_H/H_open`.  `.false.` ⇒ the original loops, unchanged.
      real(wp), intent(in) :: dt_inner   !! BT-substep dt (s)

      integer :: i, j, k, nx, ny, nz
      real(wp) :: htot_face

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (metrics%use_closed_faces) then
         call bt_rem_wave_drag_open_impl(nx, ny, nz, ms%h_layer, metrics%open_u, &
                                         metrics%open_v, bt_work%lwd_drag_u, &
                                         bt_work%lwd_drag_v, dt_inner, &
                                         bt_work%bt_rem_u, bt_work%bt_rem_v)
         return
      end if

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

   pure subroutine bt_rem_wave_drag_open_impl(nx, ny, nz, h_layer, open_u, open_v, &
                                              lwd_u, lwd_v, dt_inner, bt_rem_u, bt_rem_v)
      !! `compute_bt_rem_wave_drag` under `&vcoord_nml zfixed_closed_faces`:
      !! MULTIPLIES `H/(H + r_H·dt_inner)` into `bt_rem` with the OPEN-column
      !! face depth `H = Σ_k h_face·open` (the `bt_rem_open_impl` depth).
      !! `H <= 0` (every layer closed) ⇒ unmodified, MOM6's guard.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: open_u(nx + 1, ny, nz), open_v(nx, ny + 1, nz)
      real(wp), intent(in) :: lwd_u(nx + 1, ny), lwd_v(nx, ny + 1)
      real(wp), intent(in) :: dt_inner
      real(wp), intent(inout) :: bt_rem_u(nx + 1, ny), bt_rem_v(nx, ny + 1)
      integer :: i, j, k
      real(wp) :: htot_face

      do concurrent(j=1:ny, i=2:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + &
                        0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))*open_u(i, j, k)
         end do
         if (htot_face > 0.0_wp) then
            bt_rem_u(i, j) = bt_rem_u(i, j)*(htot_face/(htot_face + lwd_u(i, j)*dt_inner))
         end if
      end do

      do concurrent(j=2:ny, i=1:nx) local(k, htot_face)
         htot_face = 0.0_wp
         do k = 1, nz
            htot_face = htot_face + &
                        0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))*open_v(i, j, k)
         end do
         if (htot_face > 0.0_wp) then
            bt_rem_v(i, j) = bt_rem_v(i, j)*(htot_face/(htot_face + lwd_v(i, j)*dt_inner))
         end if
      end do
   end subroutine bt_rem_wave_drag_open_impl

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

#include "rdb_frhat_face.inc"

end module rdb_barotropic_coupling
