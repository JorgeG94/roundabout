!! Ideal-age tracer kernels for the ocean path.
module rdb_ocean_ideal_age
   !! Passive ideal-age tracer: interior tendency `dA/dt = 1` (s/s) at
   !! thermo cadence, surface-layer value held at a Dirichlet condition
   !! `A = A_young(t)` (0 by default). With the standard advection +
   !! diffusion this gives "time since this parcel last touched the
   !! surface" — a diagnostic of circulation pathways and spurious
   !! diapycnal mixing.
   !!
   !! Storage: `hTr_age` is `h_layer · age` (m · s), like S/T; the
   !! diagnostic reports `age = hTr_age / h_layer`. `k = 1` bed,
   !! `k = nz` surface.
   !!
   !! **Two-entry-point call contract** (see `src/core/ocean/dynamics/
   !! split_rk2/rdb_ocean_dyn.F90`, `run_stage[_split]` /
   !! `ocean_dyn_step[_split]`):
   !!
   !!   - `ocean_ideal_age_apply` — the interior-aging SOURCE term.
   !!     Called once per RK2 STAGE, thermo-cadence gated
   !!     (`dt = therm_dt`, `active = dyn%is_thermo_step()`). The RK2
   !!     average turns the per-stage `+therm_dt` into exactly one
   !!     `+therm_dt` per outer step — by design, not by accident.
   !!   - `ocean_ideal_age_reset_surface` — the surface Dirichlet BC.
   !!     Called exactly ONCE per outer step, AFTER `rk2_average`,
   !!     AFTER `continuity_tracer_drain`, and AFTER the ALE remap —
   !!     i.e. it is the last operator to touch `hTr_age(:, :, nz)`
   !!     each step. A reset applied inside an RK2 stage is halved by
   !!     the subsequent average (`0.5*(A0 + 0) = 0.5*A0`, not `0`) —
   !!     that was the historical bug this split fixes. A reset applied
   !!     before the ALE remap is overwritten by the remap's vertical
   !!     redistribution of subsurface age into the new top cell.
   !!
   !! Both entry points self-gate on `ms%idx_age <= 0` (no-op when the
   !! ideal-age tracer isn't registered) and keep the outer-shim +
   !! flat-impl split: the `ms%tracers(idx)%hTr` registry dereference
   !! happens on the host, only flat arrays cross into `do concurrent`
   !! (two-level derived-type indirection segfaults on the GPU).
   !!
   !! Two deliberate divergences from the MOM6 reference
   !! (`ideal_age_example.F90`):
   !!   - **Units are SECONDS, not years** (`hTr_age` units `s`,
   !!     `rdb_multilayer_state.F90`); MOM6 works in years
   !!     (`Isecs_per_year`). Do not convert — the diagnostic registry
   !!     and the console status line both assume seconds.
   !!   - **Reset-after-remap.** MOM6's `ideal_age_tracer_column_physics`
   !!     runs in the diabatic driver, itself followed by `ALE_main`.
   !!     Roundabout places the reset explicitly after its own ALE remap
   !!     call, which is strictly stronger (idempotent, exact) and
   !!     matches the tracer's definition as a hard Dirichlet condition
   !!     rather than a relaxation.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   implicit none
   private

   public :: ocean_ideal_age_apply
   public :: ocean_ideal_age_age_step
   public :: ocean_ideal_age_reset_surface
   public :: ocean_ideal_age_reset_step
   public :: ocean_ideal_age_young_val

contains

   subroutine ocean_ideal_age_apply(grid, ms, dt, active)
      !! Driver-facing entry point for the interior-aging SOURCE term.
      !! Self-gates on `ms%idx_age <= 0` (no-op when the ideal-age
      !! tracer isn't registered), else delegates to the flat-impl
      !! `ocean_ideal_age_age_step` kernel. Call once per RK2 stage
      !! with `dt = therm_dt` and `active = dyn%is_thermo_step()` — see
      !! the module header for the full contract.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      logical, intent(in), optional :: active
         !! Optional gate (thermo cadence). Absent => kernel runs;
         !! present-and-false => early return.

      if (present(active)) then
         if (.not. active) return
      end if
      if (ms%idx_age <= 0) return
      call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                    dt, grid%nx_total, grid%ny_total, ms%nz_ml)
   end subroutine ocean_ideal_age_apply

   pure subroutine ocean_ideal_age_age_step(hTr_age, h_layer, dt, nx, ny, nz)
      !! Interior aging, `k = 1 .. nz-1` (subsurface layers only —
      !! `k = nz` is the surface and is owned exclusively by
      !! `ocean_ideal_age_reset_step`):
      !!   `hTr_age(i,j,k) += dt * h_layer(i,j,k)`
      !! `dt` here is the caller's `therm_dt`, so this is a per-RK2-
      !! stage source; two stages + `rk2_average_field_3d` net exactly
      !! one `+dt` per outer step.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: hTr_age(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      do concurrent(k=1:nz - 1, j=1:ny, i=1:nx)
         hTr_age(i, j, k) = hTr_age(i, j, k) + dt*h_layer(i, j, k)
      end do
   end subroutine ocean_ideal_age_age_step

   subroutine ocean_ideal_age_reset_surface(grid, ms, young_eff, active)
      !! Driver-facing entry point for the surface Dirichlet BC.
      !! Self-gates on `ms%idx_age <= 0`, else delegates to the
      !! flat-impl `ocean_ideal_age_reset_step` kernel. Call exactly
      !! ONCE per outer step, after `rk2_average` + the ALE remap —
      !! see the module header for the full ordering contract.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: young_eff
         !! Host-computed surface value (s) — see `ocean_ideal_age_young_val`.
         !! Computed ONCE on the host per outer step and passed by value;
         !! never call `exp()` inside a `do concurrent` body.
      logical, intent(in), optional :: active
         !! Optional gate (thermo cadence). Absent => kernel runs;
         !! present-and-false => early return.

      if (present(active)) then
         if (.not. active) return
      end if
      if (ms%idx_age <= 0) return
      call ocean_ideal_age_reset_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                      young_eff, grid%nx_total, grid%ny_total, ms%nz_ml)
   end subroutine ocean_ideal_age_reset_surface

   pure subroutine ocean_ideal_age_reset_step(hTr_age, h_layer, young_eff, nx, ny, nz)
      !! Surface Dirichlet BC, `k = nz` only:
      !!   `hTr_age(i,j,nz) = young_eff * h_layer(i,j,nz)`
      !! `young_eff` is `A_young(t)` (s), a host scalar from
      !! `ocean_ideal_age_young_val` — this writes a CONCENTRATION
      !! times thickness, not a raw concentration (today's hard-coded
      !! `young_eff = 0` reset never exercised the `h_layer` factor).
      !! Writes the full `1:nx, 1:ny` extent including ghosts/land, same
      !! as the historical kernel; land `h_layer -> H_VANISHED` makes
      !! the product negligible and no kernel reads land age.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: hTr_age(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: young_eff

      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         hTr_age(i, j, nz) = young_eff*h_layer(i, j, nz)
      end do
   end subroutine ocean_ideal_age_reset_step

   pure function ocean_ideal_age_young_val(young_val, sfc_growth_rate, t) result(young)
      !! Host-scalar evaluation of the surface-band age value `A_young(t)`
      !! (s), MOM6 `ideal_age_example.F90`'s `young_val` computation
      !! (`:380-385`), taken exactly in spirit including the
      !! `growth_rate == 0` short-circuit (no `exp` call, no `t`
      !! dependence) — this is what keeps the default bit-identical:
      !!   `young = young_val`                     if sfc_growth_rate == 0
      !!   `young = young_val * exp(sfc_growth_rate*t)`  otherwise
      !! `t` is model time (s) since run start. `young_val` (s) and
      !! `sfc_growth_rate` (1/s) are host scalars from `ocean_dyn_t`;
      !! call this ONCE per outer step on the host, never inside a
      !! `do concurrent` body.
      real(wp), intent(in) :: young_val
      real(wp), intent(in) :: sfc_growth_rate
      real(wp), intent(in) :: t
      real(wp) :: young

      if (sfc_growth_rate == 0.0_wp) then
         young = young_val
      else
         young = young_val*exp(sfc_growth_rate*t)
      end if
   end function ocean_ideal_age_young_val

end module rdb_ocean_ideal_age
