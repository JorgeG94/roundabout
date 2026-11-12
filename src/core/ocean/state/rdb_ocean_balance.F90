!! Balanced layer velocities from the current density field, using the model's own PGF.
module rdb_ocean_balance
   !! `ocean_balance_velocity` sets the layer velocities `u_face_x_layer`,
   !! `v_face_y_layer` in GEOSTROPHIC balance with the density field the state
   !! already carries.  The pressure gradient is not re-derived here: the
   !! model's own kernel (`ocean_pressure_force_compute`, whichever `&ocean_pgf_nml
   !! form` is configured) is run on the current `h_layer` / `rho_layer`, so the
   !! balance is exact for the discrete operator the dynamics will integrate
   !! (a balanced start does not radiate the discretisation's own imbalance).
   !!
   !! **Balance.**  With `a = (a_x, a_y)` the TOTAL horizontal pressure
   !! acceleration (the layer PGF plus whatever free-surface term the
   !! barotropic substep will add, see below), the steady state of
   !!
   !!    du/dt = +f v + a_x,        dv/dt = -f u + a_y
   !!
   !! is
   !!
   !!    u =  a_y / f ,             v = -a_x / f .
   !!
   !! On the C-grid `a_x` lives on u-faces and `a_y` on v-faces while `u` needs
   !! `a_y`, so each component takes the 4-point average of the other
   !! component's face accelerations (the two faces either side, on each of the
   !! two adjacent rows/columns):
   !!
   !!    u(i,j)   = (1/f_u(i,j)) * 1/4 [ a_y(i-1,j) + a_y(i,j) + a_y(i-1,j+1) + a_y(i,j+1) ]
   !!    v(i,j)   = -(1/f_v(i,j)) * 1/4 [ a_x(i,j-1) + a_x(i,j) + a_x(i+1,j-1) + a_x(i+1,j) ]
   !!
   !! `f_u = (f_corner(i,j) + f_corner(i,j+1))/2`, `f_v = (f_corner(i,j) +
   !! f_corner(i+1,j))/2` (corner `(i,j)` is the SW corner of cell `(i,j)`).
   !! These are the discrete operator `seed_balanced_vortex_ic` uses on its
   !! analytic potential: differencing a corner-averaged `phi` is identical to
   !! this 4-point average of the face differences.  At array edges where one
   !! of the four faces does not exist the available ones are averaged.
   !!
   !! **Free surface.**  The Montgomery / FV-lite / FV-Wright forms are measured
   !! from the free surface and carry no `-g grad(eta)`; the barotropic substep
   !! supplies it with `bt_work%g_bt`.  FV-MOM6 and gprime carry
   !! `pgf_free_surface_gravity(pgf)` of it themselves.  The net term added
   !! here is therefore `-(g_bt - g_pf) grad(eta)` with `eta = sum_k h - bt_H_ref`,
   !! so the balance sees exactly the acceleration the split step will apply.
   !!
   !! **Reference level.**  `ref_level = "keep_h"` trusts the thickness the
   !! caller supplied (their `eta`).  `"bottom_rest"` (default) first chooses the
   !! sea surface that leaves the BOTTOM pressure anomaly zero (no net load on
   !! the bed: no deep flow, a surface-intensified balanced eddy, as in the
   !! vortex context of McWilliams & Flierl, 1979):
   !!
   !!    eta(i,j) = -(1/rho_0) sum_k (rho_k - rho_ref,k) h_k
   !!
   !! with `rho_ref,k` the per-layer MEDIAN density over the wet interior
   !! (the far-field background of a compact anomaly; a plain mean would add a
   !! uniform offset equal to the anomaly's area-mean, which leaves the
   !! velocities unchanged but shifts the absolute sea level).  Each column is
   !! then stretched uniformly, `h_k <- h_k (b + eta)/sum h`, with every tracer's
   !! `c` divided by the same factor so the content is kept (a sigma-like
   !! choice, not an isopycnal one).
   !!
   !! **Residency (GPU `-gpu=mem:separate`).**  Must be called AFTER
   !! `ocean_state_enter_data` (`engine_enter_data`): the PGF and EOS kernels are
   !! `do concurrent` over device-mapped state and scratch.  The HOST arrays
   !! `h_layer` and every tracer `c` are the inputs and are pushed to the device
   !! first (`update device` on leaf arrays only); the host `u/v/hu/hv/h_layer/
   !! c/rho_layer` are the outputs and are pushed/pulled symmetrically, so on
   !! return host and device agree.  Never an aggregate derived-type update.
   !!
   !! **Scope.**  Geostrophic only (`BALANCE_GRADIENT` is reserved and fails
   !! loud).  Cartesian and spherical lon-lat grids; refuses
   !! tripolar / supergrid, porous barriers, wet/dry, ice-shelf cavities, and
   !! any face where `|f| < f_min` (the equatorial band has no geostrophic
   !! solution).  Every refusal returns an error status and leaves `u`, `v`
   !! (and `h`) untouched.
   !!
   !! **Decomposed domains (bitwise identical to one rank).**  Every stencil
   !! above is local; the only global quantities are (i) the per-layer MEDIAN
   !! `rho_ref,k`, computed on more than one rank by an exact bisection over
   !! the ordered bit patterns of the wet interior densities (global counts via
   !! an exact allreduce, so the result is the SAME element the one-rank
   !! quickselect returns), (ii) the `|f|` floor (global min) and (iii) the
   !! bottom-rest failure flag (global max, so every rank returns the same
   !! status).  The inputs are first made seam-valid (host exchange of `h`,
   !! tracers, and the prescribed `rho` when thermodynamics is off), and the
   !! produced `u`, `v`, `hu`, `hv` get a final host exchange, so a seam ghost
   !! is the neighbour's owned value rather than the one-sided edge average.
   !! On one rank the ghost-cell velocities come out of the same stencil on the
   !! ghost thickness/density; the first model step's halo exchange owns them.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, H_VANISHED
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_error_ring, only: fail
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IC_SEED, &
                               OCEAN_STATUS_ERR_NOT_INITIALISED
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_metrics, only: parse_grid_config, GRID_CONFIG_CARTESIAN, &
                                GRID_CONFIG_SPHERICAL
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pressure_force, only: ocean_pressure_force_compute
   use rdb_barotropic_coupling, only: derive_bt_from_layers, pgf_free_surface_gravity
   use rdb_ocean_halo, only: ocean_halo_is_decomposed, ocean_halo_centre, &
                             ocean_halo_face_x, ocean_halo_face_y
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_halo, only: halo_allreduce_min, halo_allreduce_max, halo_allreduce_sum_i8
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ocean_balance_velocity
   public :: BALANCE_NONE, BALANCE_GEOSTROPHIC, BALANCE_GRADIENT

   integer, parameter :: BALANCE_NONE = 0
      !! No balance requested: `ocean_balance_velocity` is a successful no-op.
   integer, parameter :: BALANCE_GEOSTROPHIC = 1
      !! Geostrophic balance with the model's own PGF.
   integer, parameter :: BALANCE_GRADIENT = 2
      !! Gradient-wind balance.  Reserved: not implemented yet (fails loud).

   real(wp), parameter :: F_MIN_DEFAULT = 1.0e-6_wp
      !! Default `|f|` floor (s^-1) below which geostrophy is refused.

contains

   subroutine ocean_balance_velocity(grid, cfg, state, mode, ierr, ref_level, f_min)
      !! Set `u_face_x_layer` / `v_face_y_layer` (and `hu`/`hv`, the barotropic
      !! state, and for `ref_level="bottom_rest"` the thicknesses) in balance
      !! with the density field; see the module header for the discretisation,
      !! sign convention and device contract.  Call after `enter_data`.
      !!
      !! `hu_face_x_layer` / `hv_face_y_layer` are `u * h_face` with `h_face` the
      !! mean of the two abutting cells (array-edge faces: the single cell), the
      !! convention of the balanced-vortex seed.  Faces closed by the land mask
      !! (`wet_u`/`wet_v` = 0) and faces abutting a vanished layer
      !! (`h <= H_VANISHED`) get zero velocity.
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
         !! Supplies the refusal gates (grid kind, wet/dry, thermodynamics).
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: mode
         !! `BALANCE_NONE` (no-op), `BALANCE_GEOSTROPHIC`; `BALANCE_GRADIENT`
         !! fails loud.
      integer, intent(out) :: ierr
         !! `OCEAN_STATUS_OK`, or `OCEAN_STATUS_ERR_IC_SEED` (refusal / bad
         !! argument; message in the error ring), or `..._NOT_INITIALISED`.
      character(len=*), intent(in), optional :: ref_level
         !! `"bottom_rest"` (default) or `"keep_h"`.
      real(wp), intent(in), optional :: f_min
         !! `|f|` floor, default `1e-6` s^-1.

      character(len=16) :: ref
      real(wp) :: fmin
      logical :: rest, rho_from_eos
      integer :: nx, ny, nz
      real(wp), allocatable :: scale_col(:, :), ax(:, :, :), ay(:, :, :)
      real(wp), allocatable :: u_new(:, :, :), v_new(:, :, :)

      ierr = OCEAN_STATUS_OK
      if (mode == BALANCE_NONE) return
      if (mode == BALANCE_GRADIENT) then
         call fail("ocean_balance_velocity: BALANCE_GRADIENT (gradient-wind) is not "// &
                   "implemented yet; use BALANCE_GEOSTROPHIC", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (mode /= BALANCE_GEOSTROPHIC) then
         call fail("ocean_balance_velocity: unknown balance mode "//to_string(mode), &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ref = "bottom_rest"
      if (present(ref_level)) ref = trim(adjustl(ref_level))
      if (ref /= "bottom_rest" .and. ref /= "keep_h") then
         call fail("ocean_balance_velocity: ref_level must be 'bottom_rest' or 'keep_h'; got '"// &
                   trim(ref)//"'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      rest = (ref == "bottom_rest")
      fmin = F_MIN_DEFAULT
      if (present(f_min)) fmin = f_min

      if (.not. state%is_init .or. .not. allocated(state%multilayer%h_layer) .or. &
          .not. allocated(state%pressure_force%dpdx_face%data) .or. &
          .not. allocated(state%coriolis_adv%f_corner)) then
         call fail("ocean_balance_velocity: ocean state is not configured "// &
                   "(call after engine_setup / engine_enter_data)", ierr, OCEAN_STATUS_ERR_NOT_INITIALISED)
         return
      end if
      if (.not. refusals_clear(grid, cfg, state, fmin, ierr)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = state%multilayer%nz_ml
      rho_from_eos = cfg%ocean%thermo%enable_thermodynamics .and. &
                     state%multilayer%idx_salinity > 0 .and. state%multilayer%idx_temperature > 0

      ! ---- decomposed: seam-valid host inputs (no-op on one rank) ----
      if (ocean_halo_is_decomposed()) then
         call ocean_halo_exchange_ml_state(state%multilayer, device_resident=.false.)
         if (.not. rho_from_eos) call ocean_halo_centre(state%multilayer%rho_layer, nz, .false.)
      end if

      ! ---- inputs host -> device (leaf arrays only) ----
      call push_thickness_tracers(state)
      call refresh_density(state, rho_from_eos)

      ! ---- step 1: reference level (SSH) and uniform column stretch ----
      if (rest) then
         allocate (scale_col(nx, ny))
         call bottom_rest_scale(grid, state, scale_col, ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call apply_column_scale(nx, ny, nz, scale_col, state%multilayer%h_layer)
         call push_thickness_tracers(state)
         call refresh_density(state, rho_from_eos)
         deallocate (scale_col)
      end if

      ! ---- step 3: the model's own PGF on the (stretched) columns ----
      call ocean_pressure_force_compute(grid, state%metrics, state%pressure_force, &
                                        state%multilayer, state%eos)
      call pull_pgf(state)

      ! ---- steps 4-5: total acceleration, f at faces, 4-point stencil (host) ----
      allocate (ax(nx + 1, ny, nz), ay(nx, ny + 1, nz), u_new(nx + 1, ny, nz), v_new(nx, ny + 1, nz))
      call total_acceleration_impl(nx, ny, nz, state%pressure_force%dpdx_face%data, &
                                   state%pressure_force%dpdy_face%data, &
                                   state%metrics%idxCu, state%metrics%idyCv, &
                                   state%multilayer%h_layer, balance_href(state), &
                                   state%dyn%bt_work%g_bt - pgf_free_surface_gravity(state%pressure_force), &
                                   ax, ay)
      call geostrophic_faces_impl(nx, ny, nz, ax, ay, state%coriolis_adv%f_corner, &
                                  state%metrics%wet_u, state%metrics%wet_v, &
                                  state%multilayer%h_layer, u_new, v_new)

      ! ---- step 7: commit u, v, hu, hv; push; barotropic state ----
      state%multilayer%u_face_x_layer = u_new
      state%multilayer%v_face_y_layer = v_new
      call mass_fluxes_impl(nx, ny, nz, u_new, v_new, state%multilayer%h_layer, &
                            state%multilayer%hu_face_x_layer, state%multilayer%hv_face_y_layer)
      if (ocean_halo_is_decomposed()) then
         call ocean_halo_face_x(state%multilayer%u_face_x_layer, nz, .false.)
         call ocean_halo_face_y(state%multilayer%v_face_y_layer, nz, .false.)
         call ocean_halo_face_x(state%multilayer%hu_face_x_layer, nz, .false.)
         call ocean_halo_face_y(state%multilayer%hv_face_y_layer, nz, .false.)
      end if
      if (allocated(state%barotropic%h)) call set_column_depth(state)
      call push_velocity(state)
      if (allocated(state%dyn%bt_work%bt_eta) .and. allocated(state%dyn%bt_work%bt_H_ref)) then
         call derive_bt_from_layers(grid, state%dyn%bt_work, state%multilayer, state%metrics)
      end if

      if (rest) then
         call logger%info("ocean_balance_velocity: geostrophic, ref_level=bottom_rest")
      else
         call logger%info("ocean_balance_velocity: geostrophic, ref_level=keep_h")
      end if
   end subroutine ocean_balance_velocity

   ! ------------------------------------------------------------------
   ! refusal gates
   ! ------------------------------------------------------------------

   function refusals_clear(grid, cfg, state, fmin, ierr) result(ok)
      !! `.true.` when no unsupported feature is active and `|f| >= fmin` on every
      !! open interior face; otherwise fails loud (status in `ierr`).  Pure
      !! host reads: mutates nothing.
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(in) :: state
      real(wp), intent(in) :: fmin
      integer, intent(inout) :: ierr
      logical :: ok

      integer :: gk, ng, nx, ny, i, j
      real(wp) :: fu, fv, fworst, fglob

      ok = .false.
      gk = parse_grid_config(trim(cfg%ocean%grid%grid_config))
      if (gk /= GRID_CONFIG_CARTESIAN .and. gk /= GRID_CONFIG_SPHERICAL) then
         call fail("ocean_balance_velocity: grid_config='"//trim(cfg%ocean%grid%grid_config)// &
                   "' is not supported (cartesian and spherical only; tripolar fold and "// &
                   "supergrid angle rotation are untested)", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (state%bc%north_fold) then
         call fail("ocean_balance_velocity: the tripolar north fold is not supported", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (state%metrics%use_porous) then
         call fail("ocean_balance_velocity: porous barriers are not supported", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (cfg%ocean%wetdry%enable) then
         call fail("ocean_balance_velocity: wet/dry is not supported", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (state%metrics%use_cavity) then
         call fail("ocean_balance_velocity: ice-shelf cavities are not supported", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ! Equatorial band: the smallest |f| on any open interior face.
      ng = grid%nghost
      nx = grid%nx_total
      ny = grid%ny_total
      fworst = huge(1.0_wp)
      do j = ng + 1, ny - ng
         do i = ng + 1, nx - ng + 1
            if (state%metrics%wet_u(i, j) < 0.5_wp) cycle
            fu = 0.5_wp*(state%coriolis_adv%f_corner(i, j) + state%coriolis_adv%f_corner(i, j + 1))
            fworst = min(fworst, abs(fu))
         end do
      end do
      do j = ng + 1, ny - ng + 1
         do i = ng + 1, nx - ng
            if (state%metrics%wet_v(i, j) < 0.5_wp) cycle
            fv = 0.5_wp*(state%coriolis_adv%f_corner(i, j) + state%coriolis_adv%f_corner(i + 1, j))
            fworst = min(fworst, abs(fv))
         end do
      end do
      ! Rank-uniform verdict: every rank sees the global worst face.
      if (ocean_halo_is_decomposed()) then
         call halo_allreduce_min(fworst, fglob)
         fworst = fglob
      end if
      if (fworst < fmin) then
         call fail("ocean_balance_velocity: |f| = "//to_string(fworst)//" s^-1 < f_min = "// &
                   to_string(fmin)//" s^-1 on an open face: no geostrophic solution at the "// &
                   "equator; move the domain off the equatorial band or lower f_min", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      ok = .true.
   end function refusals_clear

   ! ------------------------------------------------------------------
   ! device sync (leaf arrays only)
   ! ------------------------------------------------------------------

   subroutine push_thickness_tracers(state)
      !! Host `h_layer` and every tracer `c` -> device.
      type(ocean_state_t), intent(inout) :: state
      integer :: it
      associate (hl => state%multilayer%h_layer)
         !$acc update device(hl)
      end associate
      if (allocated(state%multilayer%tracers)) then
         do it = 1, size(state%multilayer%tracers)
            associate (htr => state%multilayer%tracers(it)%c)
               !$acc update device(htr)
            end associate
         end do
      end if
   end subroutine push_thickness_tracers

   subroutine refresh_density(state, from_eos)
      !! `rho_layer` from the EOS on the device, then host copy refreshed; or,
      !! when thermodynamics is off (prescribed density), push the host density.
      type(ocean_state_t), intent(inout) :: state
      logical, intent(in) :: from_eos
      if (from_eos) then
         call ocean_eos_compute(state%eos, state%multilayer, active=.true.)
         associate (rl => state%multilayer%rho_layer)
            !$acc update self(rl)
         end associate
      else
         associate (rl => state%multilayer%rho_layer)
            !$acc update device(rl)
         end associate
      end if
   end subroutine refresh_density

   subroutine pull_pgf(state)
      !! Face accelerations device -> host (the scratch buffers' leaf `data`).
      type(ocean_state_t), intent(inout) :: state
      associate (dpx => state%pressure_force%dpdx_face%data, &
                 dpy => state%pressure_force%dpdy_face%data)
         !$acc update self(dpx, dpy)
      end associate
   end subroutine pull_pgf

   subroutine push_velocity(state)
      !! Host `u, v, hu, hv` (and the legacy column depth) -> device.
      type(ocean_state_t), intent(inout) :: state
      associate (uf => state%multilayer%u_face_x_layer, vf => state%multilayer%v_face_y_layer, &
                 huf => state%multilayer%hu_face_x_layer, hvf => state%multilayer%hv_face_y_layer)
         !$acc update device(uf, vf, huf, hvf)
      end associate
      if (allocated(state%barotropic%h)) then
         associate (hb => state%barotropic%h)
            !$acc update device(hb) if_present
         end associate
      end if
   end subroutine push_velocity

   subroutine set_column_depth(state)
      !! Legacy `barotropic%h = b + eta = sum_k h_layer` (host; the seed's convention).
      type(ocean_state_t), intent(inout) :: state
      integer :: i, j, k
      do j = 1, size(state%barotropic%h, 2)
         do i = 1, size(state%barotropic%h, 1)
            state%barotropic%h(i, j) = 0.0_wp
            do k = 1, state%multilayer%nz_ml
               state%barotropic%h(i, j) = state%barotropic%h(i, j) + state%multilayer%h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine set_column_depth

   function balance_href(state) result(href)
      !! Reference depth `eta = sum h - href` is measured from: the barotropic
      !! workstate's `bt_H_ref` when set, else the bathymetry `b`.
      type(ocean_state_t), intent(in) :: state
      real(wp), allocatable :: href(:, :)
      if (allocated(state%dyn%bt_work%bt_H_ref)) then
         href = state%dyn%bt_work%bt_H_ref
      else
         href = state%barotropic%b
      end if
   end function balance_href

   ! ------------------------------------------------------------------
   ! step 1: reference level
   ! ------------------------------------------------------------------

   subroutine bottom_rest_scale(grid, state, scale_col, ierr)
      !! Per-column stretch factor `(b + eta)/sum h` with `eta` the zero
      !! bottom-pressure-anomaly sea level (module header).  Columns that are
      !! land (`wet_T < 0.5`, `b <= 0`) or empty keep scale 1.  Fails (status,
      !! nothing modified) if a wet column would be left with `b + eta <= 0`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(in) :: state
      real(wp), intent(out) :: scale_col(:, :)
      integer, intent(inout) :: ierr

      integer :: i, j, k, ng, nx, ny, nz, n
      integer(int64) :: n_global
      real(wp) :: rho0, eta, tot, bad, bad_g
      logical :: failed
      real(wp), allocatable :: rho_ref(:), buf(:)

      nx = grid%nx_total
      ny = grid%ny_total
      nz = state%multilayer%nz_ml
      ng = grid%nghost
      rho0 = state%pressure_force%rho0
      allocate (rho_ref(nz), buf((nx - 2*ng)*(ny - 2*ng)))
      do k = 1, nz
         n = 0
         do j = ng + 1, ny - ng
            do i = ng + 1, nx - ng
               if (state%metrics%wet_T(i, j) < 0.5_wp) cycle
               n = n + 1
               buf(n) = state%multilayer%rho_layer(i, j, k)
            end do
         end do
         if (ocean_halo_is_decomposed()) then
            call global_lower_median(buf, n, rho_ref(k), n_global)
         else
            n_global = n
            if (n > 0) call lower_median(buf, n, rho_ref(k))
         end if
         if (n_global == 0) then
            call fail("ocean_balance_velocity: no wet interior cells", ierr, OCEAN_STATUS_ERR_IC_SEED)
            return
         end if
      end do

      scale_col = 1.0_wp
      failed = .false.
      do j = 1, ny
         do i = 1, nx
            if (state%metrics%wet_T(i, j) < 0.5_wp) cycle
            if (state%barotropic%b(i, j) <= 0.0_wp) cycle
            eta = 0.0_wp
            tot = 0.0_wp
            do k = 1, nz
               eta = eta - (state%multilayer%rho_layer(i, j, k) - rho_ref(k)) &
                     *state%multilayer%h_layer(i, j, k)/rho0
               tot = tot + state%multilayer%h_layer(i, j, k)
            end do
            if (tot <= 0.0_wp) cycle
            if (state%barotropic%b(i, j) + eta <= 0.0_wp) then
               ! Ghost columns count too (the one-rank run stretches them with
               ! the same formula), so the verdict is the same on every rank.
               if (.not. failed) then
                  call fail("ocean_balance_velocity: bottom_rest sea level empties the column at ("// &
                            to_string(i)//","//to_string(j)//"): b + eta = "// &
                            to_string(state%barotropic%b(i, j) + eta)//" m", ierr, OCEAN_STATUS_ERR_IC_SEED)
               end if
               failed = .true.
               cycle
            end if
            scale_col(i, j) = (state%barotropic%b(i, j) + eta)/tot
         end do
      end do
      ! Collective: a rank that saw no bad column must still return the refusal.
      if (ocean_halo_is_decomposed()) then
         bad = merge(1.0_wp, 0.0_wp, failed)
         call halo_allreduce_max(bad, bad_g)
         if (bad_g > 0.5_wp .and. .not. failed) then
            call fail("ocean_balance_velocity: bottom_rest sea level empties a column on another rank", &
                      ierr, OCEAN_STATUS_ERR_IC_SEED)
         end if
      end if
   end subroutine bottom_rest_scale

   subroutine global_lower_median(a, n, m, n_global)
      !! Lower median over the union of every rank's `a(1:n)` (collective),
      !! identical to `lower_median` on the gathered array: the smallest value
      !! `x` with `count(a <= x) >= (N + 1)/2`, found by bisection on the
      !! monotone int64 image of the IEEE bit pattern.  Counts ride the exact
      !! real64 allreduce (`halo_allreduce_sum_i8`).  `n_global = N`; `m` is
      !! untouched when `N = 0`.  Non-finite entries are not supported.
      integer, intent(in) :: n
      real(wp), intent(in) :: a(:)
      real(wp), intent(out) :: m
      integer(int64), intent(out) :: n_global
      integer(int64) :: lo, hi, mid, cnt, kth, nl
      real(wp) :: vmin, vmax, gmin, gmax
      integer :: i

      nl = int(n, int64)
      call halo_allreduce_sum_i8(nl, n_global)
      if (n_global == 0) return
      vmin = huge(1.0_wp)
      vmax = -huge(1.0_wp)
      do i = 1, n
         vmin = min(vmin, a(i))
         vmax = max(vmax, a(i))
      end do
      call halo_allreduce_min(vmin, gmin)
      call halo_allreduce_max(vmax, gmax)
      kth = (n_global + 1_int64)/2_int64
      lo = ordered_key(gmin)
      hi = ordered_key(gmax)
      do while (lo < hi)
         mid = lo/2_int64 + hi/2_int64 + (mod(lo, 2_int64) + mod(hi, 2_int64))/2_int64
         mid = max(lo, min(mid, hi - 1_int64))
         nl = 0_int64
         do i = 1, n
            if (ordered_key(a(i)) <= mid) nl = nl + 1_int64
         end do
         call halo_allreduce_sum_i8(nl, cnt)
         if (cnt >= kth) then
            hi = mid
         else
            lo = mid + 1_int64
         end if
      end do
      m = key_value(lo)
   end subroutine global_lower_median

   pure function ordered_key(x) result(key)
      !! Monotone int64 image of a finite double (negative values flip).
      real(wp), intent(in) :: x
      integer(int64) :: key
      integer(int64) :: b
      b = transfer(x, b)
      if (b >= 0_int64) then
         key = b
      else
         key = -iand(b, huge(b))
      end if
   end function ordered_key

   pure function key_value(key) result(x)
      !! Inverse of `ordered_key`.
      integer(int64), intent(in) :: key
      real(wp) :: x
      integer(int64) :: b
      if (key >= 0_int64) then
         b = key
      else
         b = ior(-key, -huge(b) - 1_int64)
      end if
      x = transfer(b, x)
   end function key_value

   pure subroutine lower_median(a, n, m)
      !! Lower median `m` of `a(1:n)` by quickselect.  `a` is a scratch buffer:
      !! it is partially reordered.
      integer, intent(in) :: n
      real(wp), intent(inout) :: a(:)
      real(wp), intent(out) :: m
      integer :: lo, hi, kth, i, j
      real(wp) :: pivot, t

      kth = (n + 1)/2
      lo = 1
      hi = n
      do while (lo < hi)
         pivot = a((lo + hi)/2)
         i = lo
         j = hi
         do while (i <= j)
            do while (a(i) < pivot)
               i = i + 1
            end do
            do while (a(j) > pivot)
               j = j - 1
            end do
            if (i <= j) then
               t = a(i)
               a(i) = a(j)
               a(j) = t
               i = i + 1
               j = j - 1
            end if
         end do
         if (kth <= j) then
            hi = j
         else if (kth >= i) then
            lo = i
         else
            exit
         end if
      end do
      m = a(kth)
   end subroutine lower_median

   pure subroutine apply_column_scale(nx, ny, nz, scale_col, arr)
      !! `arr(i,j,k) *= scale_col(i,j)` (host).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: scale_col(nx, ny)
      real(wp), intent(inout) :: arr(nx, ny, nz)
      integer :: i, j, k
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               arr(i, j, k) = arr(i, j, k)*scale_col(i, j)
            end do
         end do
      end do
   end subroutine apply_column_scale

   ! ------------------------------------------------------------------
   ! steps 4-5 (host stencil kernels; explicit-shape)
   ! ------------------------------------------------------------------

   pure subroutine total_acceleration_impl(nx, ny, nz, dpdx, dpdy, idxCu, idyCv, &
                                           h_layer, href, g_net, ax, ay)
      !! `a = PGF - g_net * grad(eta)`, `eta = sum_k h - href`, differenced across
      !! each interior face with the face metric.  Array-edge faces carry the PGF
      !! unchanged (the kernel zeroes them).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dpdx(nx + 1, ny, nz), dpdy(nx, ny + 1, nz)
      real(wp), intent(in) :: idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: h_layer(nx, ny, nz), href(nx, ny)
      real(wp), intent(in) :: g_net
      real(wp), intent(out) :: ax(nx + 1, ny, nz), ay(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: s
      real(wp) :: eta(nx, ny)

      do j = 1, ny
         do i = 1, nx
            s = 0.0_wp
            do k = 1, nz
               s = s + h_layer(i, j, k)
            end do
            eta(i, j) = s - href(i, j)
         end do
      end do
      do k = 1, nz
         do j = 1, ny
            ax(1, j, k) = dpdx(1, j, k)
            ax(nx + 1, j, k) = dpdx(nx + 1, j, k)
            do i = 2, nx
               ax(i, j, k) = dpdx(i, j, k) - g_net*(eta(i, j) - eta(i - 1, j))*idxCu(i, j)
            end do
         end do
         do i = 1, nx
            ay(i, 1, k) = dpdy(i, 1, k)
            ay(i, ny + 1, k) = dpdy(i, ny + 1, k)
            do j = 2, ny
               ay(i, j, k) = dpdy(i, j, k) - g_net*(eta(i, j) - eta(i, j - 1))*idyCv(i, j)
            end do
         end do
      end do
   end subroutine total_acceleration_impl

   pure subroutine geostrophic_faces_impl(nx, ny, nz, ax, ay, f_corner, wet_u, wet_v, &
                                          h_layer, u, v)
      !! `u = a_y/f`, `v = -a_x/f` by the 4-point face averages of the module
      !! header; zero on closed (`wet_* < 0.5`) faces and faces next to a
      !! vanished layer.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: ax(nx + 1, ny, nz), ay(nx, ny + 1, nz)
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1)
      real(wp), intent(in) :: wet_u(nx + 1, ny), wet_v(nx, ny + 1)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(out) :: u(nx + 1, ny, nz), v(nx, ny + 1, nz)
      integer :: i, j, k, il, ir, jl, jr, n
      real(wp) :: f, s, hmin

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx + 1
               u(i, j, k) = 0.0_wp
               if (wet_u(i, j) < 0.5_wp) cycle
               il = max(1, i - 1)
               ir = min(nx, i)
               hmin = min(h_layer(il, j, k), h_layer(ir, j, k))
               ! skip, no geostrophic velocity next to a vanished layer
               if (hmin <= H_VANISHED) cycle
               s = 0.0_wp
               n = 0
               ! v-faces (i-1,j), (i,j), (i-1,j+1), (i,j+1): those that exist
               if (i - 1 >= 1) then
                  s = s + ay(i - 1, j, k) + ay(i - 1, j + 1, k)
                  n = n + 2
               end if
               if (i <= nx) then
                  s = s + ay(i, j, k) + ay(i, j + 1, k)
                  n = n + 2
               end if
               f = 0.5_wp*(f_corner(i, j) + f_corner(i, j + 1))
               if (abs(f) > 0.0_wp .and. n > 0) u(i, j, k) = (s/real(n, wp))/f
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               v(i, j, k) = 0.0_wp
               if (wet_v(i, j) < 0.5_wp) cycle
               jl = max(1, j - 1)
               jr = min(ny, j)
               hmin = min(h_layer(i, jl, k), h_layer(i, jr, k))
               ! skip, no geostrophic velocity next to a vanished layer
               if (hmin <= H_VANISHED) cycle
               s = 0.0_wp
               n = 0
               ! u-faces (i,j-1), (i,j), (i+1,j-1), (i+1,j): those that exist
               if (j - 1 >= 1) then
                  s = s + ax(i, j - 1, k) + ax(i + 1, j - 1, k)
                  n = n + 2
               end if
               if (j <= ny) then
                  s = s + ax(i, j, k) + ax(i + 1, j, k)
                  n = n + 2
               end if
               f = 0.5_wp*(f_corner(i, j) + f_corner(i + 1, j))
               if (abs(f) > 0.0_wp .and. n > 0) v(i, j, k) = -(s/real(n, wp))/f
            end do
         end do
      end do
   end subroutine geostrophic_faces_impl

   pure subroutine mass_fluxes_impl(nx, ny, nz, u, v, h_layer, hu, hv)
      !! `hu = u * h_face`, `hv = v * h_face`; `h_face` the mean of the two
      !! abutting cells, the single cell at an array edge.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: u(nx + 1, ny, nz), v(nx, ny + 1, nz), h_layer(nx, ny, nz)
      real(wp), intent(out) :: hu(nx + 1, ny, nz), hv(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: hf

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx + 1
               if (i == 1) then
                  hf = h_layer(1, j, k)
               else if (i == nx + 1) then
                  hf = h_layer(nx, j, k)
               else
                  hf = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
               end if
               hu(i, j, k) = u(i, j, k)*hf
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               if (j == 1) then
                  hf = h_layer(i, 1, k)
               else if (j == ny + 1) then
                  hf = h_layer(i, ny, k)
               else
                  hf = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
               end if
               hv(i, j, k) = v(i, j, k)*hf
            end do
         end do
      end do
   end subroutine mass_fluxes_impl

end module rdb_ocean_balance
