!! Ice -> ocean adapters of the CICE-port sea ice: surface-flux fields and the face stress blend (roundabout-original).
module rdb_ice_cice_adapters
   !! Turn the ice's window-mean cell fluxes (`ice_flux_avg_t`, or the fluxes of one
   !! step handed over directly) into the ocean's surface forcing fields
   !! (MASTER_PLAN C.4), in the ocean's units and signs.  Each adapter is ONE small
   !! device kernel over the OWNED cells (ghosts keep their seed); the ice-side
   !! signs are Icepack's (`rdb_ice_flux_state`).
   !!
   !! | ice flux (cell, per s) | ocean field | formula | sign / unit on the ocean side |
   !! |---|---|---|---|
   !! | `fhocn` [W m-2, + into ocean] | `ice_heat` (comp. on) / `Q_heat` (off) | `fhocn`; `Q_heat_const + fhocn + fswthru` when off | + warms, W m-2 |
   !! | `fswthru` [W m-2, + down] | `ice_sw` (comp. on) / `Q_heat` (off) | `fswthru` | + warms, W m-2; `Q_heat` CONTAINS the shortwave |
   !! | `fsalt`, `fresh` [kg m-2 s-1] | `ice_salt` (comp. on) / `Q_salt` (off) | `1000*fsalt - S_surf*fresh` (+ `Q_salt_const` when off) | + salinifies, psu kg m-2 s-1 |
   !! | `fresh` | optional output of `ice_cice_apply_window` | `fresh` | + into ocean, kg m-2 s-1 |
   !! | `fswthru_vis` | not adapted here (the visible-band consumer is a later PR); stays in the averager's `acc_fswthru_vis` | | |
   !!
   !! **Write order and ownership of the ocean slots (components on).**  Per outer step the
   !! ocean's shared components are written by `engine_step` (data forcing, atm state, bulk
   !! forcing: `q_sw`, `q_lw`, `q_lat`, `q_sens`, `evap`, `heat_added` by a reader), every
   !! step, as FULL overwrites; the ice hand-off runs at the thermo step (after them);
   !! `ocean_surface_flux_assemble` (`engine_step_finalize`, thermo step) sums them.  An
   !! adapter writing `q_sw`/`heat_added`/`salt_flux` would therefore erase the open-water
   !! shortwave and a reader's heat at every hand-off.  So the CICE ice owns its OWN
   !! components `ice_heat`, `ice_sw`, `ice_salt` (`ocean_surface_flux_t%set_ice_components`,
   !! allocated only under this scheme), overwritten at each hand-off and ADDED by the
   !! assembler after the shared sum; nobody else writes them, nothing needs zeroing, and
   !! under `scheme = "sis2"` they do not exist.  The new scheme requires
   !! `enable_components`; the components-off branch below (`Q_heat`/`Q_salt` overwrite, the
   !! SIS2 couplers' form) exists for the equivalence test and for runs without the bulk
   !! forcing.  The under-ice shortwave is added to `Q_heat` at the surface layer; it is not
   !! penetrated (`sw_source = "q_sw"` reads `q_sw` only) until the visible-band PR.
   !!
   !! **Seam ghosts.**  The adapters write OWNED cells; the ocean's surface-flux application
   !! reads a seam ghost, so `ice_cice_apply_window` ends with the centre-field halo
   !! exchange (`ocean_halo_exchange_centre_field`: two-pass, periodic wrap, tripolar fold)
   !! of every field it wrote (needs `ocean_halo_init`; a no-op on one non-periodic rank).
   !!
   !! **Virtual-salt form (interim, until the real-mass PR).**  The ocean is
   !! volume conserving, so a fresh-water flux enters as a salt flux that has the
   !! same effect on the surface salinity: `dS/dt = Q_salt/(rho0*h)` with
   !! `Q_salt = 1000*fsalt - S_surf*fresh`.  `1000*fsalt` is the salt the ice
   !! moved (kg salt -> psu kg), `S_surf*fresh` the dilution (`fresh > 0`) or
   !! concentration (`fresh < 0`) of the ocean's own surface salinity by the water
   !! the ice moved.  `S_surf` is the FIRST LIVE layer's salinity in the OCEAN
   !! convention (`ice_ocean_in_t%sss_ocn`) at the hand-off, one ice step old: the
   !! averager carries fluxes, not salinities, so with the window mean the
   !! identity `applied*window = 1000*sum(fsalt*dt) - S_surf*sum(fresh*dt)`
   !! holds for the `S_surf` of the hand-off.  It also equals the SIS2 virtual-salt form
   !! `(m*S_surf - m*S_ice)/dt` ONLY when the ocean carries Practical Salinity: under the
   !! absolute-salinity option (`cons_abs`) the dilution uses the ocean's own salinity
   !! variable (right for the ocean's salt CONTENT), which differs from the Practical value by
   !! the Reference-Salinity factor, so that identity is off there.  Signs by construction: ice growth
   !! (`dhi > 0`): `fresh = -rhoi*dhi/dt`, `fsalt = -Sref*1e-3*rhoi*dhi/dt`, so
   !! `Q_salt = rhoi*dhi*(S_surf - Sref)/dt > 0` (brine rejection salinifies);
   !! melt (`dhi < 0`): `Q_salt = -rhoi*|dhi|*(S_surf - Sref)/dt < 0` (freshens).
   !! The water itself is NOT given to the ocean (no sea-level change; `fresh`
   !! enters only through the dilution term); it is returned as `fresh_mass_out`
   !! so the real-mass PR can take it, at which point `Q_salt` becomes
   !! `1000*fsalt` alone.  The ocean's `seaice_melt`/`lprec` mass components are
   !! deliberately NOT written (it would count the water twice).
   !!
   !! **`update_ocn_f`.**  The frazil mass and salt reach `fresh`/`fsalt` only when
   !! the ice runs with `update_ocn_f = .true.` (the OM3 value); the adapters are
   !! agnostic and pass whatever the ice emitted.  With `.false.` the ice gains
   !! mass and salt the ocean does not lose: the coupling does NOT compensate (the
   !! plan keeps the knob as legal and documented), and the ledger reports the
   !! unpaid `rhoi*frazil_vol` as the water/salt residual.  The HEAT side is
   !! independent of the knob: frazil heat never enters `fhocn` (the ocean paid it
   !! when the bank was filled).
   !!
   !! **Open-water weighting of the ocean's own fluxes** is not an adapter: the
   !! ocean's bulk forcing already weights by `ocean_sea_ice_t%aice0`
   !! (D1.6), and `ice_aggregates_fill` copies the CICE aggregates into it when
   !! `ocean_sea_ice_t%cice` is initialised (the open-fraction hook, inert under
   !! `scheme = "sis2"`).
   !!
   !! **Stress.**  `ice_cice_stress_blend` is the retarget point of the face stress
   !! blend: the same `(1 - a)*tau_atm + a*tau_ice` as the SIS2 coupler with
   !! `a` the face mean of the CICE aggregate `cice%aice` and `tau_ice` the
   !! per-ice-area stress ON the ocean (`-strocn/ai`, MASTER_PLAN C.4) from the
   !! dynamics.  **`tau_a` must be the UNWEIGHTED wind stress.**  The bulk forcing already
   !! sends `aice0 x wind stress` to the faces now that the open fraction is live; blending
   !! that again by `(1 - a)` weights the wind by open water twice.  The wiring PR snapshots
   !! the stress BEFORE the open-water weighting and passes the snapshot.  The blend is
   !! mirrored into `ice_cice_stress_mirror_t` (restart tags `ice_cice_tau_x/y`, the valid
   !! scalar `ice_cice_tau_valid`, registered only under this scheme) and
   !! `ice_cice_stress_resume_apply` restores it at configure, so the first ocean step after a
   !! restart sees the blended stress, not the wind.
   !!
   !! **Hand-off order for the wiring PRs, every outer step:** ice step;
   !! `ice_flux_avg_accumulate`; the stress blend; and at `is_thermo_step()` only:
   !! `ice_cice_apply_window` (adapters, then `ice_flux_avg_reset`).
   !! `ice_cice_apply_window` refuses an EMPTY window (`n_acc = 0`) instead of
   !! zeroing the ocean's forcing.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, ocean_surface_stress_refresh_mag
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress, ocean_halo_exchange_centre_field
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_avg, only: ice_flux_avg_t, ice_flux_avg_reset
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux_impl
   implicit none
   private

   public :: ice_cice_to_surface_flux, ice_cice_apply_window, ice_cice_stress_blend
   public :: ice_cice_stress_mirror_t, ice_cice_stress_resume_apply
   public :: ICE_ADAPT_ERR_NO_COMPONENT, ICE_ADAPT_ERR_WINDOW
   public :: ice_adapt_heat_impl, ice_adapt_set_impl, ice_adapt_salt_impl
   public :: ICE_ADAPT_OK, ICE_ADAPT_ERR_SHAPE, ICE_ADAPT_ERR_EMPTY
   public :: PSU_PER_MASS_FRACTION

   real(wp), parameter :: PSU_PER_MASS_FRACTION = 1000.0_wp
      !! [psu per (kg salt / kg water)]: the ocean carries salinity in ppt-like units.

   integer, parameter :: ICE_ADAPT_OK = 0
      !! Call ran.
   integer, parameter :: ICE_ADAPT_ERR_SHAPE = 1
      !! A slot is uninitialised or the shapes disagree; nothing ran.
   integer, parameter :: ICE_ADAPT_ERR_NO_COMPONENT = 3
      !! Components are on but the ice components were not allocated; nothing ran.
   integer, parameter :: ICE_ADAPT_ERR_WINDOW = 4
      !! The accumulated time differs from the expected thermo-window length; nothing ran.
   integer, parameter :: ICE_ADAPT_ERR_EMPTY = 2
      !! The averaging window holds no ice step; nothing ran, the ocean fields are untouched.

   type :: ice_cice_stress_mirror_t
      !! Restart-carried copy of the blended face stress (the new scheme's `tau_ocn_*`).
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
         !! Total extent in x including ghosts.
      integer :: ny = 0
         !! Total extent in y including ghosts.
      real(wp) :: valid = 0.0_wp
         !! Host scalar: 1 once a blend was recorded or restored (never device-mapped).
      real(wp), allocatable :: tau_x(:, :)
         !! [N m-2] blended stress on the u faces `(nx+1, ny)`.
      real(wp), allocatable :: tau_y(:, :)
         !! [N m-2] blended stress on the v faces `(nx, ny+1)`.
   contains
      procedure, non_overridable :: init => mirror_init
      procedure, non_overridable :: destroy => mirror_destroy
      procedure, non_overridable :: enter_data => mirror_enter_data
      procedure, non_overridable :: exit_data => mirror_exit_data
      procedure, non_overridable :: bytes => mirror_bytes
      procedure, non_overridable :: register_restart => mirror_register_restart
   end type ice_cice_stress_mirror_t

contains

   subroutine ice_cice_to_surface_flux(grid, sf, fresh, fsalt, fhocn, fswthru, sss_ocn, ierr)
      !! Write the ocean surface forcing from the (window-mean or single-step) ice
      !! fluxes, the three adapters of the module table, in the two modes of the
      !! SIS2 couplers (`sf%use_components`): components OFF overwrites `Q_heat`
      !! (`Q_heat_const + fhocn`, then `+ fswthru`) and `Q_salt`; components ON
      !! overwrites the ice's OWN `ice_heat` (`fhocn`), `ice_sw` (`fswthru`) and
      !! `ice_salt`, summed by the assembler (`ICE_ADAPT_ERR_NO_COMPONENT` if they were
      !! not allocated).  Latches `has_heat`, `has_salt`.  No halo exchange here (see
      !! `ice_cice_apply_window`).  OWNED cells only.  All arrays are `(nx_total, ny_total)`; mapped on an
      !! offload build.
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: fresh(:, :), fsalt(:, :), fhocn(:, :), fswthru(:, :), sss_ocn(:, :)
      integer, intent(out) :: ierr

      integer :: nx, ny, is, ie, js, je

      ierr = ICE_ADAPT_OK
      nx = grid%nx_total
      ny = grid%ny_total
      if (.not. sf%is_init .or. any([size(fresh, 1), size(fsalt, 1), size(fhocn, 1), size(fswthru, 1), &
                                     size(sss_ocn, 1)] /= nx) .or. &
          any([size(fresh, 2), size(fsalt, 2), size(fhocn, 2), size(fswthru, 2), size(sss_ocn, 2)] /= ny)) then
         ierr = ICE_ADAPT_ERR_SHAPE
         return
      end if
      is = grid%nghost + 1
      ie = grid%nghost + grid%nx_phys
      js = grid%nghost + 1
      je = grid%nghost + grid%ny_phys
      if (sf%use_components) then
         if (.not. sf%use_ice_components) then
            ierr = ICE_ADAPT_ERR_NO_COMPONENT
            return
         end if
         call ice_adapt_set_impl(nx, ny, is, ie, js, je, fhocn, sf%ice_heat)
         call ice_adapt_set_impl(nx, ny, is, ie, js, je, fswthru, sf%ice_sw)
         call ice_adapt_salt_impl(nx, ny, is, ie, js, je, fsalt, fresh, sss_ocn, 0.0_wp, sf%ice_salt)
      else
         call ice_adapt_heat_impl(nx, ny, is, ie, js, je, fhocn, fswthru, sf%Q_heat_const, sf%Q_heat)
         call ice_adapt_salt_impl(nx, ny, is, ie, js, je, fsalt, fresh, sss_ocn, sf%Q_salt_const, sf%Q_salt)
      end if
      sf%has_heat = .true.
      sf%has_salt = .true.
   end subroutine ice_cice_to_surface_flux

   subroutine ice_cice_apply_window(grid, bc, avg, sf, sss_ocn, ierr, fresh_mass_out, window_expected, t_window)
      !! The thermo-step hand-off: adapt the window mean held by `avg`
      !! (`ice_cice_to_surface_flux` on `avg%acc_*`), optionally copy the mean
      !! `fresh` [kg m-2 s-1, + into ocean] into `fresh_mass_out(nx_total, ny_total)`,
      !! exchange the seam ghosts of every field written (centre exchange + periodic
      !! wrap + tripolar fold), then `ice_flux_avg_reset`.  An empty window
      !! (`avg%n_acc < 1`) returns `ICE_ADAPT_ERR_EMPTY` and touches nothing.
      !! `t_window` returns the seconds the window accumulated; with `window_expected`
      !! (the thermo-window length, `dt_therm_ratio*dt`) a mismatch beyond 1e-9 relative
      !! returns `ICE_ADAPT_ERR_WINDOW` and touches nothing (a skipped or doubled ice step
      !! would otherwise be applied for the wrong duration).
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
         !! Supplies `north_fold` for the exchange; `ocean_halo_init` must have run.
      type(ice_flux_avg_t), intent(inout) :: avg
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: sss_ocn(:, :)
      integer, intent(out) :: ierr
      real(wp), intent(inout), optional :: fresh_mass_out(:, :)
      real(wp), intent(in), optional :: window_expected
      real(wp), intent(out), optional :: t_window

      ierr = ICE_ADAPT_OK
      if (present(t_window)) t_window = 0.0_wp
      if (.not. avg%is_init .or. avg%nx /= grid%nx_total .or. avg%ny /= grid%ny_total) then
         ierr = ICE_ADAPT_ERR_SHAPE
         return
      end if
      if (avg%n_acc < 1.0_wp) then
         ierr = ICE_ADAPT_ERR_EMPTY
         return
      end if
      if (present(t_window)) t_window = avg%t_acc
      if (present(window_expected)) then
         if (.not. (abs(avg%t_acc - window_expected) <= 1.0e-9_wp*abs(window_expected))) then
            ierr = ICE_ADAPT_ERR_WINDOW
            return
         end if
      end if
      if (present(fresh_mass_out)) then
         if (size(fresh_mass_out, 1) /= grid%nx_total .or. size(fresh_mass_out, 2) /= grid%ny_total) then
            ierr = ICE_ADAPT_ERR_SHAPE
            return
         end if
      end if
      call ice_cice_to_surface_flux(grid, sf, avg%acc_fresh, avg%acc_fsalt, avg%acc_fhocn, avg%acc_fswthru, &
                                    sss_ocn, ierr)
      if (ierr /= ICE_ADAPT_OK) return
      if (present(fresh_mass_out)) then
         call ice_adapt_set_impl(grid%nx_total, grid%ny_total, grid%nghost + 1, grid%nghost + grid%nx_phys, &
                                 grid%nghost + 1, grid%nghost + grid%ny_phys, avg%acc_fresh, fresh_mass_out)
         call ocean_halo_exchange_centre_field(fresh_mass_out, grid, bc)
      end if
      if (sf%use_components) then
         call ocean_halo_exchange_centre_field(sf%ice_heat, grid, bc)
         call ocean_halo_exchange_centre_field(sf%ice_sw, grid, bc)
         call ocean_halo_exchange_centre_field(sf%ice_salt, grid, bc)
      else
         call ocean_halo_exchange_centre_field(sf%Q_heat, grid, bc)
         call ocean_halo_exchange_centre_field(sf%Q_salt, grid, bc)
      end if
      call ice_flux_avg_reset(avg)
   end subroutine ice_cice_apply_window

   pure subroutine ice_adapt_heat_impl(nx, ny, is, ie, js, je, fhocn, fswthru, Q_const, Q_heat)
      !! `Q_heat = (Q_const + fhocn) + fswthru` [W m-2, + warms the ocean]: the
      !! operand order of the SIS2 heat then shortwave couplers.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: fhocn(nx, ny), fswthru(nx, ny)
      real(wp), intent(in) :: Q_const
      real(wp), intent(inout) :: Q_heat(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         Q_heat(i, j) = (Q_const + fhocn(i, j)) + fswthru(i, j)
      end do
   end subroutine ice_adapt_heat_impl

   pure subroutine ice_adapt_set_impl(nx, ny, is, ie, js, je, src, dst)
      !! `dst = src` over the owned cells: `heat_added <- fhocn`, `q_sw <- fswthru`,
      !! the mass flux `<- fresh`.  Signs unchanged (ice sign = ocean sign).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: src(nx, ny)
      real(wp), intent(inout) :: dst(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         dst(i, j) = src(i, j)
      end do
   end subroutine ice_adapt_set_impl

   pure subroutine ice_adapt_salt_impl(nx, ny, is, ie, js, je, fsalt, fresh, sss_ocn, Q_const, Q_salt)
      !! `Q_salt = Q_const + (1000*fsalt - sss_ocn*fresh)` [psu kg m-2 s-1,
      !! + salinifies]: the virtual-salt form of the module docstring.  A
      !! non-finite input gives a non-finite `Q_salt` (no clamp).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: fsalt(nx, ny), fresh(nx, ny), sss_ocn(nx, ny)
      real(wp), intent(in) :: Q_const
      real(wp), intent(inout) :: Q_salt(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         Q_salt(i, j) = Q_const + (PSU_PER_MASS_FRACTION*fsalt(i, j) - sss_ocn(i, j)*fresh(i, j))
      end do
   end subroutine ice_adapt_salt_impl

   subroutine ice_cice_stress_blend(stress, tau_a_x, tau_a_y, tau_ice_x, tau_ice_y, cice, mirror, grid, bc)
      !! Face stress blend of the CICE ice (retarget of `ice_ocean_stress_flux`):
      !! `tau = (1 - a)*tau_atm + a*tau_ice`, `a` the face mean of `cice%aice` (half
      !! weight at the array edge faces, as the SIS2 blend), a FULL overwrite of
      !! `stress%tau_x/tau_y`, then `stress_mag` (and, with `grid` and `bc`, the seam
      !! ghosts) is re-derived and the blend is copied into `mirror`.
      !! **`tau_a_*` MUST BE THE UNWEIGHTED WIND STRESS** (`(nx+1, ny)` / `(nx, ny+1)`
      !! snapshot taken BEFORE the bulk forcing's open-water weighting), else the wind is
      !! weighted by open water twice; `tau_ice_*` is the per-ice-area stress on the ocean.
      !! Nothing calls this under `scheme = "sis2"`.
      type(ocean_surface_stress_t), intent(inout) :: stress
      real(wp), intent(in) :: tau_a_x(:, :), tau_a_y(:, :), tau_ice_x(:, :), tau_ice_y(:, :)
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_cice_stress_mirror_t), intent(inout) :: mirror
      type(hgrid_t), intent(in), optional :: grid
      type(ocean_bc_state_t), intent(in), optional :: bc

      call ice_ocean_stress_flux_impl(stress%tau_x, stress%tau_y, tau_a_x, tau_a_y, tau_ice_x, tau_ice_y, &
                                      cice%aice, cice%nx, cice%ny)
      if (present(grid) .and. present(bc)) then
         call ocean_seam_refresh_surface_stress(stress, grid, bc)
      else
         call ocean_surface_stress_refresh_mag(stress)
      end if
      ! Mirror AFTER the seam refresh, so a checkpoint carries the exchanged ghosts.
      call ice_cice_tau_mirror_impl(mirror%tau_x, mirror%tau_y, stress%tau_x, stress%tau_y, cice%nx, cice%ny)
      mirror%valid = 1.0_wp
   end subroutine ice_cice_stress_blend

   pure subroutine ice_cice_tau_mirror_impl(m_x, m_y, tau_x, tau_y, nx, ny)
      !! Device copy of the blend into the restart-carried mirror.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: m_x(nx + 1, ny), m_y(nx, ny + 1)
      real(wp), intent(in) :: tau_x(nx + 1, ny), tau_y(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         m_x(i, j) = tau_x(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         m_y(i, j) = tau_y(i, j)
      end do
   end subroutine ice_cice_tau_mirror_impl

   subroutine ice_cice_stress_resume_apply(stress, mirror)
      !! Configure-time resume (host, BEFORE `enter_data`): if the checkpoint carried a
      !! blend (`mirror%valid > 0.5`) overwrite the freshly re-seeded `stress%tau_x/y` with
      !! it and re-derive `stress_mag`; otherwise leave the wind (a fresh run).
      type(ocean_surface_stress_t), intent(inout) :: stress
      type(ice_cice_stress_mirror_t), intent(in) :: mirror

      if (mirror%valid <= 0.5_wp) return
      stress%tau_x = mirror%tau_x
      stress%tau_y = mirror%tau_y
      call ocean_surface_stress_refresh_mag(stress)
   end subroutine ice_cice_stress_resume_apply

   subroutine mirror_init(this, grid)
      !! Allocate the two face arrays (zero), not valid.
      class(ice_cice_stress_mirror_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (this%is_init) call this%destroy()
      this%nx = grid%nx_total
      this%ny = grid%ny_total
      allocate (this%tau_x(this%nx + 1, this%ny), source=0.0_wp)
      allocate (this%tau_y(this%nx, this%ny + 1), source=0.0_wp)
      this%valid = 0.0_wp
      this%is_init = .true.
      this%mapped = .false.
   end subroutine mirror_init

   subroutine mirror_destroy(this)
      class(ice_cice_stress_mirror_t), intent(inout) :: this
      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%tau_x)) deallocate (this%tau_x)
      if (allocated(this%tau_y)) deallocate (this%tau_y)
      this%valid = 0.0_wp
   end subroutine mirror_destroy

   subroutine mirror_enter_data(this)
      class(ice_cice_stress_mirror_t), intent(inout) :: this
      select type (this)
      type is (ice_cice_stress_mirror_t)
         if (.not. this%is_init .or. this%mapped) return
         !$acc enter data copyin(this%tau_x, this%tau_y)
         this%mapped = .true.
      end select
   end subroutine mirror_enter_data

   subroutine mirror_exit_data(this)
      class(ice_cice_stress_mirror_t), intent(inout) :: this
      select type (this)
      type is (ice_cice_stress_mirror_t)
         if (.not. this%is_init .or. .not. this%mapped) return
         !$acc exit data delete(this%tau_x, this%tau_y)
         this%mapped = .false.
      end select
   end subroutine mirror_exit_data

   pure function mirror_bytes(this) result(nbytes)
      class(ice_cice_stress_mirror_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = arr_bytes(this%tau_x) + arr_bytes(this%tau_y)
   end function mirror_bytes

   subroutine mirror_register_restart(this, reg)
      !! Optional restart entries `ice_cice_tau_x/y` and the host scalar `ice_cice_tau_valid`.
      !! Call ONLY under the new scheme.  The caller pulls the device arrays first.
      class(ice_cice_stress_mirror_t), intent(inout) :: this
      type(restart_registry_t), intent(inout) :: reg
      if (.not. this%is_init) return
      call reg%register_2d("ice_cice_tau_x", this%tau_x, 0, size(this%tau_x, 1), size(this%tau_x, 2), optional=.true.)
      call reg%register_2d("ice_cice_tau_y", this%tau_y, 0, size(this%tau_y, 1), size(this%tau_y, 2), optional=.true.)
      call reg%register_scalar("ice_cice_tau_valid", this%valid, optional=.true.)
   end subroutine mirror_register_restart

end module rdb_ice_cice_adapters
