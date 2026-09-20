!! Ice-shelf basal melt, COUPLED: the slot, the far-field sampler and the
!! once-per-thermo-step driver that turns `rdb_ocean_cavity_melt`'s scalar
!! three-equation kernel into two surface-flux components the ocean
!! integrates.
!!
!! `&ocean_cavity_melt_nml enable` (default `.false.`) is the gate; the
!! geometry it stands on (`z_draft`, `cover_frac`, `p_ice_ref`) is
!! `&ocean_cavity_dyn_nml`'s and is required.  Knob off ⇒ the slot's
!! arrays stay at their `(1,1)` placeholders, no kernel is launched, and
!! every path is byte-identical.
!!
!! ## What this module owns
!!
!! 1. `ocean_cavity_flux_t` — the 2-D device-resident state of the
!!    interface: the sampled far field, `u*`, the interface `(T_b, S_b)`,
!!    the canonical melt mass flux, the ocean → interface heat flux, and
!!    the per-column solver status.
!! 2. `cavity_far_field_impl` — the far-field sample, in METRES below the
!!    ice base, never "layer `nz`".
!! 3. `ocean_cavity_flux_step` — the driver: sample, solve (through
!!    `cavity_melt_columns_2d`, which must own the `do concurrent`; see
!!    its docstring for the nvlink reason), deliver, count, fail loud.
!!
!! ## Far-field sampling (why metres)
!!
!! The melt rate a three-equation law returns is roughly linear in the
!! thermal driving it is handed, and the thermal driving depends on HOW
!! FAR FROM THE ICE the model sampled it.  That is the dominant
!! resolution artefact in the subject, not a refinement: Gwyther et al.
!! (2020), Burchard et al. (2022) Table 2 p. 15 (the all-bulk melt-rate
!! error GROWS under refinement — -8 % at 6.66 m, -30 % at 0.015 m) and
!! Yung et al. (2026) p. 2074 all identify it.  Sampling "the top layer"
!! would therefore make the melt rate a function of the vertical
!! coordinate's cell thickness — which is precisely the quantity a
!! coordinate study must hold fixed.  So `far_field_depth` is in metres
!! and the sampler is a thickness-weighted mean over the layers that
!! span it, with a PARTIAL last layer.
!!
!! Vanished layers (`h <= H_VANISHED`, the D4 skip/merge marker) carry no
!! mass and are skipped rather than clamped.  A column whose sample finds
!! no mass at all is dropped from the cover mask — zero melt, `OK`
!! status — not solved on a fabricated state.
!!
!! ## The three pressures, and which one this is
!!
!! The liquidus is evaluated at `multilayer_state_t%p_top` — THE
!! interface pressure, the same field the FV_MOM6 pressure-stack
!! boundary condition and the in-situ EOS read (`src/core/ocean/README.md`,
!! the `p_top` seam contract).  It is NOT `eos%p_ref` (a scalar potential
!! -density reference, deliberately horizontally uniform) and NOT the
!! `eta_forcing` seam (a gradient).  One pressure, three consumers, and
!! this module is a CONSUMER — it writes nothing to `p_top`.  The sole
!! producer is `configure_ocean_cavity`, which assembles
!! `p_top = metrics%p_ice_ref + sf%p_surf` (re-assembled per outer step
!! in `ocean_dyn_step_split` when the psurf seam makes `sf%p_surf`
!! live); `configure_ocean_cavity_melt` asserts that it ran rather than
!! seeding a second copy, because a second writer running after it would
!! silently drop `sf%p_surf`.  A consequence worth stating: an
!! atmospheric load under the shelf depresses the freezing point exactly
!! as the ice load does, with no extra wiring.
!!
!! ## Sign and units — the table that decides whether the answer is right
!!
!! | quantity | symbol | unit | sign |
!! |---|---|---|---|
!! | melt | `melt` (`m_mass`) | kg/m^2/s | **> 0 melting** |
!! | kernel heat flux | `q_ocean` | W/m^2 | **> 0 warms the INTERFACE** (the ocean cools) |
!! | delivered heat | `sf%heat_cavity` | W/m^2 | positive DOWN into the ocean |
!! | delivered salt | `sf%salt_cavity` | same as `sf%salt_flux` | positive SALINIFIES |
!!
!! **Heat.**  `q_ocean = rho_w*c_w*gamma_t*(T_w - T_b) > 0` is the
!! turbulent flux the interface takes FROM the ocean, so the component
!! the ocean must be given is its negative:
!!
!!     heat_cavity = -q_ocean
!!
!! The assembler adds it to `Q_heat` with a plain `+`, and
!! `apply_surface_src_2d_impl` deposits `dt/(rho0*cp) * Q_heat` at
!! `k = nz` — positive warms.  Warm water under a shelf therefore COOLS
!! the top of the column, which is the whole point.
!!
!! **Salt — and why it is not the same construction.**  Melting adds
!! freshwater MASS the Boussinesq column does not yet carry (that is
!! Phase 3).  Emulate the dilution with a salt flux and the fixed-mass
!! equivalent is exact, not approximate.  For a column of mass `M` per
!! unit area receiving mass `m` at salinity `S_i`:
!!
!!     d(M*S)/dt = m*S_i ,   dM/dt = m
!!  => M*dS/dt = m*S_i - S*m = -m*(S - S_i)
!!
!! so the virtual flux that reproduces `dS/dt` at FIXED `M` is
!!
!!     salt_cavity = -m_mass*(S_far - s_ice)
!!
!! referenced to the SAME far-field salinity the solve was handed.
!! Melting (`m_mass > 0`, `S_far > s_ice`) freshens.
!!
!! The heat twin of that `-S*dM/dt` dilution term is `-m*c_w*(T_w - T_b)`
!! and it is DELIBERATELY DROPPED.  The asymmetry is physical, not
!! sloppiness: `S_far - s_ice ~ 34 g/kg` is O(1) of the salinity itself,
!! while `T_w - T_b` is a few hundredths of a degree, so the dropped heat
!! term is `m/(rho_w*gamma_t) ~ 1e-3` of `q_ocean` — a 0.1 % correction
!! that belongs with the real-mass PR, whereas the salt term is the
!! entire meltwater buoyancy signal and cannot wait for it.
!!
!! NOTE the prototype's `column_demo.py` removes `f_turb = rho_w*gamma_s*
!! (S_w - S_b) = m_mass*(S_b - s_ice)` instead.  That is the salt flux
!! ACROSS the interface, not the dilution equivalent — it omits the
!! meltwater that in reality returns to the column carrying `S_b`, and so
!! under-freshens by `(S_b - s_ice)/(S_far - s_ice)`.  That demo is
!! explicitly not part of the oracle (its README: "for reading"); the
!! identity above is what the ocean's fixed-mass salinity equation
!! requires, and it is what the budget test asserts.
!!
!! ## Budgets
!!
!! Both components ride `Q_heat`/`Q_salt`, so they are integrated by
!! `ocean_surface_flux_apply_tracers` and land in the EXISTING
!! `ms%heat_budget_surface` / `ms%salt_budget_surface` contributors —
!! which the console already folds in with the correct
!! `ocean_budget_stage_weight`.  No separate frazil-style accumulator is
!! needed, and no full-weight/half-weight decision is taken here: the
!! source is applied by the same per-stage kernel as every other surface
!! flux, so it carries the same weight by construction.
!!
!! ## `mem:separate`
!!
!! Every array this module touches is device-resident and mapped by the
!! slot's `enter_data`.  Nothing is copied host↔device per step except
!! the four status COUNTS, which come back as `do concurrent ... reduce`
!! results (host scalars), exactly like `continuity_t%n_limited_step`.
!!
!! ## Citations (papers, never another model's source)
!!
!!   * Holland, D. M. and Jenkins, A. (1999): J. Phys. Oceanogr. 29,
!!     1787-1800.
!!   * Jenkins, A., Nicholls, K. W. and Corr, H. F. J. (2010): J. Phys.
!!     Oceanogr. 40, 2298-2312.
!!   * Asay-Davis, X. S. et al. (2016): Geosci. Model Dev. 9, 2471-2497
!!     (ISOMIP+).
!!   * Burchard, H. et al. (2022): Ocean Modelling 179, 102119.
!!   * Yung, C. K. et al. (2025): The Cryosphere 19, 5827-5861.
!!   * Yung, C. K. et al. (2026): The Cryosphere 20, 2053-2088.
module rdb_ocean_cavity_flux
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t
   use rdb_mem_report, only: arr_bytes
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_cavity_melt, only: ocean_cavity_const_t, ocean_cavity_exchange_t, &
                                    ocean_cavity_ice_t, cavity_melt_columns_2d, &
                                    CAVITY_MELT_OK, CAVITY_MELT_NONFINITE_INPUT, &
                                    CAVITY_MELT_NONFINITE_STATE, &
                                    CAVITY_MELT_NOT_CONVERGED, &
                                    CAVITY_MELT_NO_PHYSICAL_ROOT
   implicit none
   private

   public :: ocean_cavity_flux_t
   public :: ocean_cavity_flux_step
   public :: cavity_far_field_impl
   public :: cavity_flux_fill_impl
   public :: cavity_status_counts_impl
   public :: cavity_melt_status_is_fatal

   type :: ocean_cavity_flux_t
      !! Ice-shelf basal-melt slot: the 2-D interface state plus the
      !! parameter bundles the kernel is called with.
      !!
      !! Every array is `(nx_total, ny_total)` when `enable`, `(1,1)`
      !! otherwise — the `z_draft` gating convention, latched in
      !! `ocean_state_init_from_config` BEFORE `init` so the allocation
      !! gate can read it.
      logical :: is_init = .false.
         !! Set by `init` after every allocation succeeds; cleared by
         !! `destroy` first.  Always test this, never `allocated(...)`.
      logical :: enable = .false.
         !! `&ocean_cavity_melt_nml enable`, latched before `init`.

      ! ---- Knobs the kernel is called with (host scalars) ----
      real(wp) :: far_field_depth = 10.0_wp
         !! Thickness (m) below the ice base the far field is averaged
         !! over.  See the module docstring: metres, not layers.
      real(wp) :: cdrag_top = 2.5e-3_wp
         !! Top drag coefficient for the MELT friction velocity only.
      real(wp) :: u_tide = 1.0e-2_wp
         !! RMS tidal velocity (m/s) in the melt `u*`.
      real(wp) :: ustar_min = 1.0e-4_wp
         !! Friction-velocity floor (m/s).
      real(wp) :: s_ice = 0.0_wp
         !! Ice salinity (g/kg).
      type(ocean_cavity_exchange_t) :: par
         !! Exchange-law selector + coefficients.  Flat POD; its
         !! `f_cor` member is unused here (the driver substitutes the
         !! per-column `f_cor` array).
      type(ocean_cavity_ice_t) :: ice
         !! Ice-conduction selector + `T_ice`.  Flat POD.
      type(ocean_cavity_const_t) :: const
         !! Thermodynamic + turbulence constants.  Flat POD, left at the
         !! ISOMIP+ protocol values (Asay-Davis et al. (2016) Table 4):
         !! `rho_w`, `c_w`, `alpha_T`, `beta_S` here are the MELT LAW's
         !! own calibrated constants, NOT copies of the model's
         !! Boussinesq `rho_0` — the same "out of scope on purpose"
         !! category as `RHO_WATER` and `&ocean_ice_nml rho_ocean` in
         !! `src/core/ocean/README.md`'s reference-density table.
         !! Overriding them would break ISOMIP+ comparability, which is
         !! the reason this path exists.

      ! ---- Static, filled once at configure ----
      real(wp), allocatable :: f_cor(:, :)
         !! Cell-centred Coriolis parameter (1/s), from the same
         !! `fill_coriolis_centre` every other slot uses.  Read by
         !! `hj99` only, as `|f|`; the law does not exist at `f = 0`.

      ! ---- Per-step interface state (all device-resident) ----
      real(wp), allocatable :: active(:, :)
         !! Composed solve mask (0/1): `cover_frac` AND wet AND "the
         !! far-field sample found mass".  This is what the kernel's
         !! `cover` argument is given.
      real(wp), allocatable :: t_far(:, :)
         !! Far-field temperature (degC), thickness-weighted over
         !! `far_field_depth`.
      real(wp), allocatable :: s_far(:, :)
         !! Far-field salinity (g/kg), same average.
      real(wp), allocatable :: u_far(:, :)
         !! Far-field x velocity at the cell CENTRE (m/s), same average.
      real(wp), allocatable :: v_far(:, :)
         !! Far-field y velocity at the cell CENTRE (m/s), same average.
      real(wp), allocatable :: ustar(:, :)
         !! Melt friction velocity (m/s).
      real(wp), allocatable :: t_b(:, :)
         !! Interface temperature (degC), on the liquidus.
      real(wp), allocatable :: s_b(:, :)
         !! Interface salinity (g/kg).
      real(wp), allocatable :: melt(:, :)
         !! **Canonical** melt mass flux (kg/m^2/s), > 0 melting.
      real(wp), allocatable :: q_ocean(:, :)
         !! Turbulent heat flux ocean → interface (W/m^2), > 0 cools the
         !! ocean.
      real(wp), allocatable :: gamma_t(:, :)
         !! Thermal exchange velocity (m/s) the column's solve converged
         !! on — NOT `par%gamma_t_coeff`, which is the dimensionless
         !! coefficient.  Under `hj99` / `yung25` it is an implicit
         !! function of the converged interface state, so it is stored
         !! rather than re-derived: a diagnostic that rebuilt it from
         !! `u*` alone would report the NEUTRAL value instead of the
         !! stratification-suppressed one.  Zero where inactive.
      real(wp), allocatable :: gamma_s(:, :)
         !! Haline exchange velocity (m/s), same convention.
      integer, allocatable :: status(:, :)
         !! `CAVITY_MELT_*` per column; `CAVITY_MELT_OK` where inactive.

      ! ---- Status counters (host scalars, reduced off the device) ----
      integer :: n_nonfinite_step = 0
         !! Active columns returning `NONFINITE_INPUT`/`NONFINITE_STATE`
         !! this step.  **FATAL** — the driver fails loud.
      integer :: n_not_converged_step = 0
         !! Active columns returning `NOT_CONVERGED` this step (zero melt
         !! applied there — the kernel's documented safe state).
      integer :: n_no_root_step = 0
         !! Active columns returning `NO_PHYSICAL_ROOT` this step.
      integer :: n_other_step = 0
         !! Active columns returning any other non-OK status
         !! (`BAD_INPUT`, `NO_CORIOLIS`, `LAW_DOMAIN`, ...).
      integer(int64) :: n_not_converged_total = 0_int64
         !! Running total, drained to the console like
         !! `continuity_t%n_limited_total`.  int64 because a
         !! long run with one stubborn column can exceed 2e9.
      integer(int64) :: n_no_root_total = 0_int64
         !! Running total of `NO_PHYSICAL_ROOT`.
      integer(int64) :: n_other_total = 0_int64
         !! Running total of every other non-OK, non-fatal status.
   contains
      procedure, non_overridable :: init => ocean_cavity_flux_init
      procedure, non_overridable :: destroy => ocean_cavity_flux_destroy
      procedure, non_overridable :: enter_data => ocean_cavity_flux_enter_data
      procedure, non_overridable :: exit_data => ocean_cavity_flux_exit_data
      procedure, non_overridable :: bytes => ocean_cavity_flux_bytes
   end type ocean_cavity_flux_t

contains

   ! ======================================================================
   ! Lifecycle
   ! ======================================================================

   subroutine ocean_cavity_flux_init(this, grid)
      !! Allocate the slot.  Gated on `enable` (latched by
      !! `ocean_state_init_from_config` before this runs), so a run
      !! without a cavity pays fourteen `(1,1)` placeholders.
      class(ocean_cavity_flux_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny

      if (this%enable) then
         nx = grid%nx_total
         ny = grid%ny_total
      else
         nx = 1
         ny = 1
      end if

      allocate (this%f_cor(nx, ny), source=0.0_wp)
      allocate (this%active(nx, ny), source=0.0_wp)
      allocate (this%t_far(nx, ny), source=0.0_wp)
      allocate (this%s_far(nx, ny), source=0.0_wp)
      allocate (this%u_far(nx, ny), source=0.0_wp)
      allocate (this%v_far(nx, ny), source=0.0_wp)
      allocate (this%ustar(nx, ny), source=0.0_wp)
      allocate (this%t_b(nx, ny), source=0.0_wp)
      allocate (this%s_b(nx, ny), source=0.0_wp)
      allocate (this%melt(nx, ny), source=0.0_wp)
      allocate (this%q_ocean(nx, ny), source=0.0_wp)
      allocate (this%gamma_t(nx, ny), source=0.0_wp)
      allocate (this%gamma_s(nx, ny), source=0.0_wp)
      allocate (this%status(nx, ny), source=CAVITY_MELT_OK)
      this%is_init = .true.
   end subroutine ocean_cavity_flux_init

   subroutine ocean_cavity_flux_destroy(this)
      !! Release the slot.  `is_init` is cleared FIRST.
      class(ocean_cavity_flux_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%f_cor)) deallocate (this%f_cor)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%t_far)) deallocate (this%t_far)
      if (allocated(this%s_far)) deallocate (this%s_far)
      if (allocated(this%u_far)) deallocate (this%u_far)
      if (allocated(this%v_far)) deallocate (this%v_far)
      if (allocated(this%ustar)) deallocate (this%ustar)
      if (allocated(this%t_b)) deallocate (this%t_b)
      if (allocated(this%s_b)) deallocate (this%s_b)
      if (allocated(this%melt)) deallocate (this%melt)
      if (allocated(this%q_ocean)) deallocate (this%q_ocean)
      if (allocated(this%gamma_t)) deallocate (this%gamma_t)
      if (allocated(this%gamma_s)) deallocate (this%gamma_s)
      if (allocated(this%status)) deallocate (this%status)
   end subroutine ocean_cavity_flux_destroy

   subroutine ocean_cavity_flux_enter_data(this)
      !! Type-bound wrapper — delegates to the non-polymorphic impl so
      !! the device-attach map base is the heap object, not a
      !! polymorphic stack box.
      class(ocean_cavity_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_cavity_flux_t)
         call ocean_cavity_flux_enter_data_impl(this)
      end select
   end subroutine ocean_cavity_flux_enter_data

   subroutine ocean_cavity_flux_enter_data_impl(this)
      !! `copyin` (not `create`) throughout: `f_cor` carries a
      !! configure-time host fill that MUST reach the device, and the
      !! rest carry the zero `init` promised on both toolchains.
      type(ocean_cavity_flux_t), intent(inout) :: this
      !$acc enter data copyin(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&                  this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&                  this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&                  this%gamma_s, this%status)
      !$acc update device(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&              this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&              this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&              this%gamma_s, this%status)
   end subroutine ocean_cavity_flux_enter_data_impl

   subroutine ocean_cavity_flux_exit_data(this)
      class(ocean_cavity_flux_t), intent(inout) :: this
      select type (this)
      type is (ocean_cavity_flux_t)
         call ocean_cavity_flux_exit_data_impl(this)
      end select
   end subroutine ocean_cavity_flux_exit_data

   subroutine ocean_cavity_flux_exit_data_impl(this)
      type(ocean_cavity_flux_t), intent(inout) :: this
      !$acc exit data delete(this%f_cor, this%active, this%t_far, this%s_far, &
      !$acc&                 this%u_far, this%v_far, this%ustar, this%t_b, &
      !$acc&                 this%s_b, this%melt, this%q_ocean, this%gamma_t, &
      !$acc&                 this%gamma_s, this%status)
   end subroutine ocean_cavity_flux_exit_data_impl

   pure function ocean_cavity_flux_bytes(this) result(nbytes)
      !! Counted allocatable footprint (0 when unallocated).  One
      !! `arr_bytes` term per array — add one here when an array joins
      !! the type.
      class(ocean_cavity_flux_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%f_cor) &
               + arr_bytes(this%active) &
               + arr_bytes(this%t_far) &
               + arr_bytes(this%s_far) &
               + arr_bytes(this%u_far) &
               + arr_bytes(this%v_far) &
               + arr_bytes(this%ustar) &
               + arr_bytes(this%t_b) &
               + arr_bytes(this%s_b) &
               + arr_bytes(this%melt) &
               + arr_bytes(this%q_ocean) &
               + arr_bytes(this%gamma_t) &
               + arr_bytes(this%gamma_s) &
               + arr_bytes(this%status)
   end function ocean_cavity_flux_bytes

   ! ======================================================================
   ! Far-field sampling
   ! ======================================================================

   pure subroutine cavity_far_field_impl(nx, ny, nz, far_depth, cover, wet_mask, &
                                         h_layer, hTr_T, hTr_S, u_face_x, v_face_y, &
                                         active, t_far, s_far, u_far, v_far)
      !! Thickness-weighted mean of `(T, S, u, v)` over `far_depth`
      !! METRES below the ice base, with a PARTIAL last layer.
      !!
      !! Bottom-up stack (`k = nz` is the surface layer, i.e. the one
      !! against the ice base), so the walk runs `k = nz, 1, -1` and
      !! stops when the budget is spent.  `far_depth >= column
      !! thickness` ⇒ the whole column, which is the correct limit and
      !! not an error.
      !!
      !! Vanished layers (`h <= H_VANISHED`) are SKIPPED, not clamped —
      !! the D4 taxonomy's skip/merge marker.  They carry no mass, so
      !! including them would be dividing a zero tracer load by a
      !! near-zero thickness.
      !!
      !! Velocities are centred from the C-grid faces per layer BEFORE
      !! the vertical average (`u_c = 1/2 (u_{i} + u_{i+1})`), not after:
      !! the two commute only for a uniform column, and the order that
      !! matches "the speed the ice base feels" is the one that averages
      !! the centred, per-layer velocity.
      !!
      !! `active` is the composed solve mask handed to the melt kernel:
      !! covered AND wet AND the sample found mass.  A covered column
      !! whose entire sample is vanished gets `active = 0` and zeroed
      !! outputs — no melt, no status noise.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      integer, intent(in) :: nz
         !! Layer count; `k = nz` is the surface / ice-base layer.
      real(wp), intent(in) :: far_depth
         !! Sampling thickness (m), > 0.
      real(wp), intent(in) :: cover(nx, ny)
         !! Ice-cover fraction (v1 binary).
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Static wet (1) / land (0) mask.
      real(wp), intent(in) :: h_layer(nx, ny, nz)
         !! Layer thickness (m).
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
         !! `h*T` (degC m).
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
         !! `h*S` ((g/kg) m).
      real(wp), intent(in) :: u_face_x(nx + 1, ny, nz)
         !! Zonal face velocity (m/s).
      real(wp), intent(in) :: v_face_y(nx, ny + 1, nz)
         !! Meridional face velocity (m/s).
      real(wp), intent(out) :: active(nx, ny)
         !! Composed solve mask, 0 or 1.
      real(wp), intent(out) :: t_far(nx, ny)
         !! Sampled temperature (degC); 0 where inactive.
      real(wp), intent(out) :: s_far(nx, ny)
         !! Sampled salinity (g/kg); 0 where inactive.
      real(wp), intent(out) :: u_far(nx, ny)
         !! Sampled centred x velocity (m/s); 0 where inactive.
      real(wp), intent(out) :: v_far(nx, ny)
         !! Sampled centred y velocity (m/s); 0 where inactive.
      integer :: i, j, k
      real(wp) :: acc_h, acc_t, acc_s, acc_u, acc_v, remain, h, w, inv_h

      do concurrent(j=1:ny, i=1:nx) &
         local(k, acc_h, acc_t, acc_s, acc_u, acc_v, remain, h, w, inv_h)
         acc_h = 0.0_wp
         acc_t = 0.0_wp
         acc_s = 0.0_wp
         acc_u = 0.0_wp
         acc_v = 0.0_wp
         remain = far_depth
         if (cover(i, j) > 0.5_wp .and. wet_mask(i, j) > 0.5_wp) then
            do k = nz, 1, -1
               if (remain <= 0.0_wp) exit
               h = h_layer(i, j, k)
               if (h <= H_VANISHED) cycle
               w = min(h, remain)
               inv_h = 1.0_wp/h
               acc_h = acc_h + w
               acc_t = acc_t + w*hTr_T(i, j, k)*inv_h
               acc_s = acc_s + w*hTr_S(i, j, k)*inv_h
               acc_u = acc_u + w*0.5_wp*(u_face_x(i, j, k) + u_face_x(i + 1, j, k))
               acc_v = acc_v + w*0.5_wp*(v_face_y(i, j, k) + v_face_y(i, j + 1, k))
               remain = remain - w
            end do
         end if
         if (acc_h > 0.0_wp) then
            active(i, j) = 1.0_wp
            t_far(i, j) = acc_t/acc_h
            s_far(i, j) = acc_s/acc_h
            u_far(i, j) = acc_u/acc_h
            v_far(i, j) = acc_v/acc_h
         else
            active(i, j) = 0.0_wp
            t_far(i, j) = 0.0_wp
            s_far(i, j) = 0.0_wp
            u_far(i, j) = 0.0_wp
            v_far(i, j) = 0.0_wp
         end if
      end do
   end subroutine cavity_far_field_impl

   ! ======================================================================
   ! Flux delivery
   ! ======================================================================

   pure subroutine cavity_flux_fill_impl(nx, ny, s_ice, active, melt, q_ocean, s_far, &
                                         heat_cavity, salt_cavity)
      !! Fill the two OWNED surface-flux components from the solved
      !! interface.  FULL OVERWRITE of the whole plane, never `+=` —
      !! the same rule the sea-ice coupler's fillers follow, and for the
      !! same reason (a `+=` ratchets across outer steps with no bound).
      !!
      !! See the module docstring for the sign derivation.  In one line:
      !! the ocean loses the heat the interface takes (`-q_ocean`), and
      !! the dilution by meltwater of salinity `s_ice` is emulated at
      !! fixed column mass by `-m_mass*(S_far - s_ice)`.
      !!
      !! No wet-mask factor here: `active` already carries it (the
      !! sampler composes it), and the assembler multiplies `Q_heat` /
      !! `Q_salt` by `ms%wet_mask` again.
      integer, intent(in) :: nx
         !! First dimension (ghosts included).
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: s_ice
         !! Ice salinity (g/kg).
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      real(wp), intent(in) :: melt(nx, ny)
         !! Melt mass flux (kg/m^2/s), > 0 melting.
      real(wp), intent(in) :: q_ocean(nx, ny)
         !! Ocean → interface heat flux (W/m^2).
      real(wp), intent(in) :: s_far(nx, ny)
         !! Far-field salinity the solve used (g/kg).
      real(wp), intent(out) :: heat_cavity(nx, ny)
         !! Heat component (W/m^2), positive down into the ocean.
      real(wp), intent(out) :: salt_cavity(nx, ny)
         !! Virtual salt component, positive salinifies.
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         if (active(i, j) > 0.5_wp) then
            heat_cavity(i, j) = -q_ocean(i, j)
            salt_cavity(i, j) = -melt(i, j)*(s_far(i, j) - s_ice)
         else
            heat_cavity(i, j) = 0.0_wp
            salt_cavity(i, j) = 0.0_wp
         end if
      end do
   end subroutine cavity_flux_fill_impl

   ! ======================================================================
   ! Status accounting
   ! ======================================================================

   pure subroutine cavity_status_counts_impl(nx, ny, active, status, n_nonfinite, &
                                             n_not_converged, n_no_root, n_other)
      !! Reduce the per-column status plane into four counts, ON DEVICE.
      !!
      !! `do concurrent ... reduce(+:)` and not `count(...)`: the array
      !! is device-resident under `mem:separate`, and a host intrinsic
      !! would silently read the stale host shadow — the same trap
      !! `continuity_t%n_limited_step` documents.
      !!
      !! Only ACTIVE columns are counted.  An inactive column always
      !! carries `CAVITY_MELT_OK`, but gating the count on the mask
      !! keeps that a property of the caller rather than of the kernel.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: active(nx, ny)
         !! Composed solve mask.
      integer, intent(in) :: status(nx, ny)
         !! `CAVITY_MELT_*` per column.
      integer, intent(out) :: n_nonfinite
         !! Count of `NONFINITE_INPUT` + `NONFINITE_STATE`.
      integer, intent(out) :: n_not_converged
         !! Count of `NOT_CONVERGED`.
      integer, intent(out) :: n_no_root
         !! Count of `NO_PHYSICAL_ROOT`.
      integer, intent(out) :: n_other
         !! Count of every other non-OK status.
      integer :: i, j, c_nf, c_nc, c_nr, c_ot

      c_nf = 0
      c_nc = 0
      c_nr = 0
      c_ot = 0
      do concurrent(j=1:ny, i=1:nx) reduce(+:c_nf, c_nc, c_nr, c_ot)
         if (active(i, j) > 0.5_wp) then
            select case (status(i, j))
            case (CAVITY_MELT_OK)
               continue
            case (CAVITY_MELT_NONFINITE_INPUT, CAVITY_MELT_NONFINITE_STATE)
               c_nf = c_nf + 1
            case (CAVITY_MELT_NOT_CONVERGED)
               c_nc = c_nc + 1
            case (CAVITY_MELT_NO_PHYSICAL_ROOT)
               c_nr = c_nr + 1
            case default
               c_ot = c_ot + 1
            end select
         end if
      end do
      n_nonfinite = c_nf
      n_not_converged = c_nc
      n_no_root = c_nr
      n_other = c_ot
   end subroutine cavity_status_counts_impl

   pure function cavity_melt_status_is_fatal(n_nonfinite) result(is_fatal)
      !! Is this step's status tally a FAIL-LOUD condition?
      !!
      !! A non-finite input or intermediate on a covered column means the
      !! state feeding the interface is already corrupt; the kernel's
      !! safe state (zero melt) would hide it behind a plausible run.
      !! `NOT_CONVERGED` and `NO_PHYSICAL_ROOT` are NOT fatal — they are
      !! counted and warned, and those columns take the documented
      !! zero-melt safe state.
      !!
      !! A `pure` predicate so the test suite can assert the DECISION
      !! without provoking `error stop`, the way the repo tests every
      !! other fail-loud rule.
      integer, intent(in) :: n_nonfinite
         !! `n_nonfinite_step` from `cavity_status_counts_impl`.
      logical :: is_fatal
      is_fatal = (n_nonfinite > 0)
   end function cavity_melt_status_is_fatal

   ! ======================================================================
   ! The driver
   ! ======================================================================

   subroutine ocean_cavity_flux_step(grid, cav, metrics, ms, eos, sf, active)
      !! One cavity basal-melt update: sample the far field, solve the
      !! three-equation interface on every covered column, deliver the
      !! two owned surface-flux components, and account for the solver
      !! status.
      !!
      !! CADENCE + PLACEMENT.  Called once per outer step at THERMO
      !! cadence from `engine_step_finalize`, immediately BEFORE
      !! `ocean_surface_flux_assemble` — the assembler must see the
      !! components this routine writes.  The tracers then integrate the
      !! assembled `Q_heat`/`Q_salt` on the NEXT outer step, which is the
      !! same one-step lag the sea-ice coupler documents.
      !!
      !! NON-`pure`, and it is the only non-pure procedure this module
      !! adds: it logs (`n_not_converged` warnings) and it FAILS LOUD on
      !! a non-finite column.  Its three kernels are all `pure`.
      !!
      !! `mem:separate`: every array read or written here is mapped by
      !! its owning slot's `enter_data`.  `ms%tracers(idx)%hTr` is
      !! dereferenced on the HOST before the kernels (the outer-shim
      !! rule for the array-of-derived-types registry).
      use pic_logger, only: global_logger
      use pic_strings, only: to_string
      use rdb_error_ring, only: fail
      use rdb_ocean_status, only: OCEAN_STATUS_ERR_SETUP
      type(hgrid_t), intent(in) :: grid
         !! Horizontal grid (for `nx_total`/`ny_total`).
      type(ocean_cavity_flux_t), intent(inout) :: cav
         !! The cavity-melt slot.
      type(ocean_metrics_t), intent(in) :: metrics
         !! Reads `cover_frac` only.
      type(multilayer_state_t), intent(in) :: ms
         !! Reads `h_layer`, the S/T tracer loads, the face velocities,
         !! `p_top` (THE interface pressure) and `wet_mask`.
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle — the liquidus.  Flat POD, by value.
      type(ocean_surface_flux_t), intent(inout) :: sf
         !! Writes `heat_cavity`/`salt_cavity` and latches
         !! `has_heat`/`has_salt` host-side.
      logical, intent(in), optional :: active
         !! Thermo-cadence gate.  Present-and-false ⇒ early return;
         !! absent ⇒ run (the `ocean_surface_flux_assemble` convention).
      integer :: nx, ny, nz, idx_t, idx_s

      if (.not. cav%is_init) return
      if (.not. cav%enable) return
      if (present(active)) then
         if (.not. active) return
      end if
      if (.not. allocated(ms%tracers)) return
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      if (idx_t <= 0 .or. idx_s <= 0) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call cavity_far_field_impl(nx, ny, nz, cav%far_field_depth, &
                                 metrics%cover_frac, ms%wet_mask, ms%h_layer, &
                                 ms%tracers(idx_t)%hTr, ms%tracers(idx_s)%hTr, &
                                 ms%u_face_x_layer, ms%v_face_y_layer, &
                                 cav%active, cav%t_far, cav%s_far, cav%u_far, cav%v_far)

      call cavity_melt_columns_2d(nx, ny, cav%active, cav%t_far, cav%s_far, ms%p_top, &
                                  cav%u_far, cav%v_far, cav%s_ice, cav%f_cor, &
                                  cav%cdrag_top, cav%u_tide, cav%ustar_min, &
                                  cav%par, cav%ice, eos, cav%const, &
                                  cav%ustar, cav%t_b, cav%s_b, cav%melt, cav%q_ocean, &
                                  cav%gamma_t, cav%gamma_s, cav%status)

      call cavity_flux_fill_impl(nx, ny, cav%s_ice, cav%active, cav%melt, cav%q_ocean, &
                                 cav%s_far, sf%heat_cavity, sf%salt_cavity)

      ! The filler's own obligations (PR-12 fill contract): latch the
      ! has_* flags HOST-side, never from a device reduction.  They stay
      ! latched for the rest of the run — the melt rate can legitimately
      ! pass through zero, and un-latching there would change the
      ! apply-path operand count mid-run.
      sf%has_heat = .true.
      sf%has_salt = .true.

      call cavity_status_counts_impl(nx, ny, cav%active, cav%status, &
                                     cav%n_nonfinite_step, cav%n_not_converged_step, &
                                     cav%n_no_root_step, cav%n_other_step)
      cav%n_not_converged_total = cav%n_not_converged_total &
                                  + int(cav%n_not_converged_step, int64)
      cav%n_no_root_total = cav%n_no_root_total + int(cav%n_no_root_step, int64)
      cav%n_other_total = cav%n_other_total + int(cav%n_other_step, int64)

      if (cavity_melt_status_is_fatal(cav%n_nonfinite_step)) then
         call fail("&ocean_cavity_melt_nml: "//to_string(cav%n_nonfinite_step)// &
                   " ice-covered column(s) fed the basal-melt solver a non-finite "// &
                   "far-field temperature, salinity, velocity or interface pressure "// &
                   "(CAVITY_MELT_NONFINITE_*).  The kernel's safe state is zero melt, "// &
                   "which would hide an already-corrupt column behind a plausible "// &
                   "run, so this is fatal.  NOT_CONVERGED / NO_PHYSICAL_ROOT are "// &
                   "counted and warned instead.", code=OCEAN_STATUS_ERR_SETUP)
      end if

      if (cav%n_not_converged_step > 0 .or. cav%n_no_root_step > 0 .or. &
          cav%n_other_step > 0) then
         call global_logger%warning("cavity melt: zero melt applied on "// &
                                    to_string(cav%n_not_converged_step)// &
                                    " non-converged, "//to_string(cav%n_no_root_step)// &
                                    " no-physical-root and "// &
                                    to_string(cav%n_other_step)// &
                                    " otherwise-refused column(s) this step")
      end if
   end subroutine ocean_cavity_flux_step

end module rdb_ocean_cavity_flux
