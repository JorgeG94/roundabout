!! Default fill routines for the ocean diagnostics registry.
module rdb_ocean_diag_fills
   !! Per-variable fill procedures the diag manager invokes on cadence-fire;
   !! the bridge between the diag manager and `ocean_state_t`.  Each fill
   !! `select type`-casts the `class(*)` state handle back to
   !! `ocean_state_t` and populates a pre-allocated buffer.  Face-staggered
   !! fields (u_face_x, v_face_y) are averaged to cell centres so every
   !! default variable is `(nx, ny, nz)`.
   !! Defaults: SSH (m), temperature (°C), salinity (PSU), u/v_centre (m/s),
   !! ke (m²/s²); conditionally ice_conc (1) / ice_thick (m) when
   !! `&ocean_ice_nml enable`.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, H_VANISHED, &
                            REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, H_VANISHED, &
                            REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4
#endif
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_diag, only: ocean_diag_t, DIAG_OP_INSTANT, DIAG_OP_MEAN, &
                             DIAG_OP_MAX, DIAG_OP_MIN, DIAG_OP_UNSET, &
                             diag_spec_t, diag_remap_proc, diag_fill_proc, &
                             DIAG_VGRID_LAYER, DIAG_VGRID_Z_FIXED, &
                             DIAG_VGRID_SIGMA, DIAG_VGRID_ZSTAR, &
                             DIAG_VGRID_DENSITY, DIAG_COORD_UNSET, &
                             DIAG_MISSING_VALUE
   use rdb_eos, only: eos_t, eos_density_point
   use rdb_ocean_vcoord, only: invert_density_targets
   use rdb_remap_column, only: remap_column
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ocean_pseudo_salt, only: ocean_pseudo_salt_deviation
   implicit none
   private
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: fill_ssh, fill_temperature, fill_salinity
   public :: fill_age
   public :: fill_pseudo_salt, fill_pseudo_salt_diff
   public :: fill_u_centre, fill_v_centre, fill_ke
   public :: fill_mld_epbl, fill_kd_epbl
   public :: fill_kd_kshear
   public :: fill_ice_conc, fill_ice_thick
   public :: remap_layer_to_z
   public :: remap_layer_to_sigma, remap_layer_to_zstar
   public :: remap_layer_to_density
   public :: coord_remap_proc
   public :: register_default_diags
   public :: is_canonical_diag_name
   public :: canonical_diag_gate_hint
   public :: canonical_diag_catalog_size, canonical_diag_catalog_name
   public :: set_diag_remap_method, parse_diag_remap_scheme
   public :: set_diag_mask_vanished, diag_mask_vanished_is_on

   integer, save :: diag_remap_method = REMAP_PPM
      !! Reconstruction for the conservative diagnostic vertical remap
      !! (`remap_column`).  Set from `&ocean_diag_nml diag_remap_scheme` via
      !! `set_diag_remap_method`; passed by value into the device `_impl`.
      !! PPM default; all schemes conservative (donor-cell overlap integral).

   logical, save :: diag_mask_vanished = .false.
      !! When `.true.`, the conservative remap fills target cells that
      !! overlap no water (below-bottom / pinched-out in a shallow column)
      !! with `DIAG_MISSING_VALUE` instead of 0.  Set from
      !! `&ocean_diag_nml mask_vanished_layers` at configure
      !! (`set_diag_mask_vanished`); read host-side and passed BY VALUE into
      !! the device `_impl`.  Default `.false.` => below-bottom cells read 0
      !! (bit-identical to the legacy remap).

   integer, parameter :: N_CANONICAL_DIAGS = 14
   character(len=16), parameter :: CANONICAL_DIAG_NAMES(N_CANONICAL_DIAGS) = &
                                   [character(len=16) :: &
                                    "SSH", "temperature", "salinity", "age", "u", "v", "KE", &
                                    "MLD_EPBL", "Kd_EPBL", "Kd_KSHEAR", &
                                    "ice_conc", "ice_thick", &
                                    "pseudo_salt", "pseudo_salt_diff"]
      !! The single source of truth for `is_canonical_diag_name` AND the
      !! P7 discoverability getters (`canonical_diag_catalog_size`/
      !! `_name`, wrapped by the C ABI as `rdb_ocean_canonical_*`) —
      !! one array, never a second hand-copied list that could drift.
      !! NOTE this is the catalog of NAMES `register_default_diags` may
      !! register, not what IS registered on a given live instance (some
      !! entries are gated — see `canonical_diag_gate_hint`); the live
      !! set is `ocean_diag_t%vars(1:nvars)%name` on an actual handle.

   real(wp), parameter :: VANISHED_TARGET_FLOOR = 1.0e-10_wp
      !! Target-cell thickness at/below which the cell is treated as
      !! overlapping no water (below-bottom / pinched-out).  Tiny absolute
      !! floor: catches the exactly-zero below-bottom cells (and z*/sigma
      !! cells in a dry column) without masking genuinely thin overlaps.

contains

   pure function parse_diag_remap_scheme(name) result(method)
      !! Map a `&ocean_diag_nml diag_remap_scheme` string to the `REMAP_*`
      !! enum; returns -1 for an unrecognised name (caller fails loud).
      character(len=*), intent(in) :: name
      integer :: method
      select case (trim(name))
      case ("pcm", "PCM")
         method = REMAP_PCM
      case ("plm", "PLM")
         method = REMAP_PLM
      case ("ppm", "PPM")
         method = REMAP_PPM
      case ("ppm_h4", "PPM_H4")
         method = REMAP_PPM_H4
      case default
         method = -1
      end select
   end function parse_diag_remap_scheme

   subroutine set_diag_remap_method(method)
      !! Seed the module-level diagnostic remap reconstruction (host only).
      integer, intent(in) :: method
      diag_remap_method = method
   end subroutine set_diag_remap_method

   subroutine set_diag_mask_vanished(flag)
      !! Enable / disable masking of vanished (no-water) remap target cells
      !! to `DIAG_MISSING_VALUE` (host only).  Default off.
      logical, intent(in) :: flag
      diag_mask_vanished = flag
   end subroutine set_diag_mask_vanished

   pure function diag_mask_vanished_is_on() result(on)
      !! Query the vanished-masking mode (used by the registration path to
      !! tag non-layer diagnostics with `has_missing` for the NetCDF writer).
      logical :: on
      on = diag_mask_vanished
   end function diag_mask_vanished_is_on

   ! Outer-shim + flat-impl: each fill recovers the concrete
   ! `ocean_state_t` via `select type` host-side, then forwards to a
   ! flat-impl `do concurrent` device kernel — keeps the `class(*)`
   ! polymorphic dispatch (blocks NVHPC device codegen) off the device.

   subroutine fill_ssh(state_handle, buf)
      !! SSH = total column thickness minus bathymetry depth, into the k=1
      !! plane of `buf`.  Multilayer or barotropic path per `use_multilayer`.
      !! Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (state%use_multilayer .and. allocated(state%multilayer%h_layer)) then
            call fill_ssh_ml_impl(state%multilayer%h_layer, state%barotropic%b, &
                                  state%multilayer%nz_ml, buf)
         else
            call fill_ssh_bt_impl(state%barotropic%h, state%barotropic%b, buf)
         end if
      end select
   end subroutine fill_ssh

   pure subroutine fill_ssh_ml_impl(h_layer, b, nz_ml, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h_layer(:, :, :)
      real(wp), intent(in)    :: b(:, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer, intent(in)    :: nz_ml
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny
      real(wp) :: col_sum
      nx = size(buf, 1)
      ny = size(buf, 2)
      do concurrent(j=1:ny, i=1:nx) local(col_sum, k)
         col_sum = 0.0_wp
         do k = 1, nz_ml
            col_sum = col_sum + h_layer(i, j, k)
         end do
         buf(i, j, 1) = col_sum - b(i, j)
      end do
   end subroutine fill_ssh_ml_impl

   pure subroutine fill_ssh_bt_impl(h, b, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h(:, :), b(:, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      nx = size(buf, 1)
      ny = size(buf, 2)
      do concurrent(j=1:ny, i=1:nx)
         buf(i, j, 1) = h(i, j) - b(i, j)
      end do
   end subroutine fill_ssh_bt_impl

   subroutine fill_temperature(state_handle, buf)
      !! Layer temperature = tracers(idx_temperature)%hTr / h_layer.
      !! Tracer-registry indirection dereferenced HOST-side before the
      !! flat-impl kernel (array-of-DT deep deref blocks NVHPC device
      !! codegen).  Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: it
      select type (state => state_handle)
      class is (ocean_state_t)
         it = state%multilayer%idx_temperature
         if (it <= 0) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_tracer_impl(state%multilayer%h_layer, &
                               state%multilayer%tracers(it)%hTr, buf)
      end select
   end subroutine fill_temperature

   subroutine fill_salinity(state_handle, buf)
      !! Layer salinity, same pattern as temperature.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: it
      select type (state => state_handle)
      class is (ocean_state_t)
         it = state%multilayer%idx_salinity
         if (it <= 0) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_tracer_impl(state%multilayer%h_layer, &
                               state%multilayer%tracers(it)%hTr, buf)
      end select
   end subroutine fill_salinity

   subroutine fill_age(state_handle, buf)
      !! Layer ideal age = tracers(idx_age)%hTr / h_layer (s).  Read-out only;
      !! registered when `&ocean_tracers_nml enable_ideal_age` (idx_age > 0).
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: it
      select type (state => state_handle)
      class is (ocean_state_t)
         it = state%multilayer%idx_age
         if (it <= 0) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_tracer_impl(state%multilayer%h_layer, &
                               state%multilayer%tracers(it)%hTr, buf)
      end select
   end subroutine fill_age

   subroutine fill_pseudo_salt(state_handle, buf)
      !! Layer pseudo-salt = tracers(idx_pseudo_salt)%hTr / h_layer (psu).
      !! Read-out only; registered when `&ocean_tracers_nml
      !! enable_pseudo_salt` (idx_pseudo_salt > 0).  Mirrors `fill_age`.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: it
      select type (state => state_handle)
      class is (ocean_state_t)
         it = state%multilayer%idx_pseudo_salt
         if (it <= 0) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_tracer_impl(state%multilayer%h_layer, &
                               state%multilayer%tracers(it)%hTr, buf)
      end select
   end subroutine fill_pseudo_salt

   subroutine fill_pseudo_salt_diff(state_handle, buf)
      !! Pseudo-salt deviation D = pseudo_salt - S (psu): a direct,
      !! measured proxy for how far the passive-tracer transport path
      !! has drifted from the active-tracer (salinity) path.  Gated on
      !! BOTH indices being registered.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: it_ps, it_s
      select type (state => state_handle)
      class is (ocean_state_t)
         it_ps = state%multilayer%idx_pseudo_salt
         it_s = state%multilayer%idx_salinity
         if (it_ps <= 0 .or. it_s <= 0) then
            call fill_zero_impl(buf)
            return
         end if
         call ocean_pseudo_salt_deviation(state%multilayer%h_layer, &
                                          state%multilayer%tracers(it_ps)%hTr, &
                                          state%multilayer%tracers(it_s)%hTr, buf, &
                                          size(buf, 1), size(buf, 2), size(buf, 3))
      end select
   end subroutine fill_pseudo_salt_diff

   pure subroutine fill_tracer_impl(h_layer, hTr, buf)
      !! Shared flat-impl for any tracer concentration field.
      !!
      !! P7 F5 fix: a vanishing layer used to write 0.0, and 0 degC / 0
      !! PSU are both LEGAL ocean values — a vanished (below-bottom or
      !! pinched-out) bed layer under ZSTAR_FULL therefore used to read
      !! back as plausible ice-point freshwater rather than as missing
      !! data (D3.2 of `06_python_surface_design.md`). Now: the divisor
      !! guard is `H_VANISHED` (the dynamic-vanish threshold, not the
      !! pure 1/0 armour `H_DIV_EPS` — see `rdb_constants`' D4 taxonomy),
      !! and the fill is IEEE NaN, matching the Python concentration
      !! accessor's own convention (D3.2) so an in-memory read
      !! (`model.diagnostic(...)`) and the NetCDF stream agree. This is
      !! independent of `&ocean_diag_nml mask_vanished_layers`: that knob
      !! only gates the conservative REMAP path's below-target-cell fill
      !! (`remap_column` via `register_one_canonical`'s `has_missing`);
      !! at the default LAYER vgrid there is no remap, so this fill is
      !! the only thing between "no water here" and a plausible-looking
      !! number, unconditionally.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: h_layer(:, :, :)
      real(wp), intent(in)    :: hTr(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz, nz_h
      real(wp) :: h, qnan
      nx = size(buf, 1)
      ny = size(buf, 2)
      nz = size(buf, 3)
      nz_h = size(h_layer, 3)
      ! Sentinel computed once, host-side, before the do-concurrent body
      ! (ieee_value is a host intrinsic — see rdb_ocean_ghost_poison's
      ! qnan for the same pattern): pass a plain IEEE bit pattern into
      ! the device kernel rather than calling the intrinsic per-cell.
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)
      do concurrent(k=1:min(nz, nz_h), j=1:ny, i=1:nx)
         h = h_layer(i, j, k)
         if (h > H_VANISHED) then
            buf(i, j, k) = hTr(i, j, k)/h
         else
            buf(i, j, k) = qnan
         end if
      end do
   end subroutine fill_tracer_impl

   subroutine fill_ice_conc(state_handle, buf)
      !! Total sea-ice concentration (0..1) at T-centres — the two-mode
      !! per-cell gather of `ice_cell_concentration_impl` (rdb_ice_state),
      !! inlined (see fill_ice_conc_thick_impl).  Registered by
      !! `register_default_diags` only when `&ocean_ice_nml enable`.
      !! Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%ice%is_init) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_ice_conc_thick_impl(state%ice%ncat, &
                                       size(state%ice%m_ice, 1), &
                                       size(state%ice%m_ice, 2), &
                                       state%metrics%wet_T, state%ice%part_size, &
                                       state%ice%m_ice, .false., buf)
      end select
   end subroutine fill_ice_conc

   subroutine fill_ice_thick(state_handle, buf)
      !! Grid-mean sea-ice thickness (m) at T-centres: mice/ICE_RHO_ICE
      !! (MOM6 effective-thickness convention) — the two-mode per-cell
      !! gather of `ice_cell_concentration_impl`, inlined (see
      !! fill_ice_conc_thick_impl).  Registered by `register_default_diags`
      !! only when `&ocean_ice_nml enable`.  Public only for the unit-test
      !! suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         if (.not. state%ice%is_init) then
            call fill_zero_impl(buf)
            return
         end if
         call fill_ice_conc_thick_impl(state%ice%ncat, &
                                       size(state%ice%m_ice, 1), &
                                       size(state%ice%m_ice, 2), &
                                       state%metrics%wet_T, state%ice%part_size, &
                                       state%ice%m_ice, .true., buf)
      end select
   end subroutine fill_ice_thick

   pure subroutine fill_ice_conc_thick_impl(ncat, nx, ny, wet_T, part_size, &
                                            m_ice, emit_thick, buf)
      !! Shared conc/thick device kernel.  Inlines the two-mode gather of
      !! `ice_cell_concentration_impl` (`rdb_ice_state` — convention of
      !! record; `test_ocean_ice_diags` pins the copies equal): ncat==1
      !! legacy lumped (per-CELL m_ice, ci = 0/1), ncat>1 SIS2 ITD
      !! (ci = min(1, Σ part_size), mice = Σ part_size·m_ice), land ⇒ 0.
      !! `emit_thick=.false.` ⇒ buf = ci; `.true.` ⇒ buf = mice/ICE_RHO_ICE
      !! (grid-mean thickness, m).  Scalar flag branch is constant-folded on
      !! the device — one kernel, no scratch companion.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      logical, intent(in) :: emit_thick
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, c, nxl, nyl
      real(wp) :: mice_val, ci_val, ci_sum

      nxl = min(nx, size(buf, 1))
      nyl = min(ny, size(buf, 2))

      if (ncat == 1) then
         do concurrent(j=1:nyl, i=1:nxl) local(mice_val, ci_val)
            if (wet_T(i, j) > 0.5_wp .and. m_ice(i, j, 1) > 0.0_wp) then
               mice_val = m_ice(i, j, 1)
               ci_val = 1.0_wp
            else
               mice_val = 0.0_wp
               ci_val = 0.0_wp
            end if
            buf(i, j, 1) = merge(mice_val/ICE_RHO_ICE, ci_val, emit_thick)
         end do
      else
         do concurrent(j=1:nyl, i=1:nxl) local(c, mice_val, ci_val, ci_sum)
            if (wet_T(i, j) > 0.5_wp) then
               mice_val = 0.0_wp
               ci_sum = 0.0_wp
               do c = 1, ncat
                  mice_val = mice_val + part_size(i, j, c)*m_ice(i, j, c)
                  ci_sum = ci_sum + part_size(i, j, c)
               end do
               ci_val = min(1.0_wp, ci_sum)
            else
               mice_val = 0.0_wp
               ci_val = 0.0_wp
            end if
            buf(i, j, 1) = merge(mice_val/ICE_RHO_ICE, ci_val, emit_thick)
         end do
      end if
   end subroutine fill_ice_conc_thick_impl

   pure subroutine fill_zero_impl(buf)
      !! Device-side zero of `buf` — used when a fill's input is missing
      !! (e.g. tracer not registered) so the fold reads a defined value.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: i, j, k, nx, ny, nz
      nx = size(buf, 1)
      ny = size(buf, 2)
      nz = size(buf, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = 0.0_wp
      end do
   end subroutine fill_zero_impl

   subroutine fill_u_centre(state_handle, buf)
      !! C-grid u-face → cell centre by simple 2-point average.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_u_centre_impl(state%multilayer%u_face_x_layer, buf)
      end select
   end subroutine fill_u_centre

   pure subroutine fill_u_centre_impl(u_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized u_face(:,:,:) has nx+1 first dim, incompatible with cell-sized buf.
      real(wp), intent(in)    :: u_face(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = min(size(buf, 1), size(u_face, 1) - 1)
      ny = min(size(buf, 2), size(u_face, 2))
      nz = min(size(buf, 3), size(u_face, 3))
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
      end do
   end subroutine fill_u_centre_impl

   subroutine fill_v_centre(state_handle, buf)
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_v_centre_impl(state%multilayer%v_face_y_layer, buf)
      end select
   end subroutine fill_v_centre

   pure subroutine fill_v_centre_impl(v_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized v_face(:,:,:) has ny+1 second dim, incompatible with cell-sized buf.
      real(wp), intent(in)    :: v_face(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = min(size(buf, 1), size(v_face, 1))
      ny = min(size(buf, 2), size(v_face, 2) - 1)
      nz = min(size(buf, 3), size(v_face, 3))
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
      end do
   end subroutine fill_v_centre_impl

   subroutine fill_ke(state_handle, buf)
      !! KE per cell on the C-grid: KE = 0.25·(u_W² + u_E² + v_S² + v_N²),
      !! the discrete C-grid KE-density (consistent with the Coriolis
      !! KE_ARAKAWA stencil).  Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_ke_impl(state%multilayer%u_face_x_layer, &
                           state%multilayer%v_face_y_layer, buf)
      end select
   end subroutine fill_ke

   pure subroutine fill_ke_impl(u_face, v_face, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! face-sized arrays have nx+1/ny+1 first/second dims, incompatible with buf.
      real(wp), intent(in)    :: u_face(:, :, :), v_face(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = min(size(buf, 1), size(u_face, 1) - 1, size(v_face, 1))
      ny = min(size(buf, 2), size(u_face, 2), size(v_face, 2) - 1)
      nz = min(size(buf, 3), size(u_face, 3), size(v_face, 3))
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = 0.25_wp*( &
                        u_face(i, j, k)*u_face(i, j, k) + &
                        u_face(i + 1, j, k)*u_face(i + 1, j, k) + &
                        v_face(i, j, k)*v_face(i, j, k) + &
                        v_face(i, j + 1, k)*v_face(i, j + 1, k))
      end do
   end subroutine fill_ke_impl

   subroutine remap_layer_to_z(state_handle, z_out, layer_buf, output_buf, is_extensive)
      !! Layer→fixed-z vertical remap, CONSERVATIVE (donor-cell overlap via
      !! `remap_column`, scheme `diag_remap_method`).  `is_extensive=.false.`
      !! ⇒ INTENSIVE (thickness-weighted average of overlapping source layers);
      !! `.true.` ⇒ EXTENSIVE (thickness-integrated field — column integral
      !! redistributed across targets, Σ preserved when the z-grid spans H).
      !! `z_out(:)` = target INTERFACE depths (m, positive-down, shallow→deep,
      !! implicit 0 surface); output cell m spans `[z_out(m-1), z_out(m)]`.
      !! Thicknesses clipped to column total H = Σ h_layer (exact
      !! conservation); below-seafloor cells read 0.  k=1 bed, k=nz surface.
      !! Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(in) :: z_out(:)
      real(wp), intent(in) :: layer_buf(:, :, :)
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in) :: is_extensive

      select type (state => state_handle)
      class is (ocean_state_t)
         call remap_layer_to_vcoord_impl(DIAG_VGRID_Z_FIXED, state%multilayer%h_layer, &
                                         z_out, layer_buf, output_buf, is_extensive, &
                                         diag_remap_method, diag_mask_vanished, &
                                         DIAG_MISSING_VALUE)
      end select
   end subroutine remap_layer_to_z

   subroutine remap_layer_to_sigma(state_handle, levels, layer_buf, output_buf, is_extensive)
      !! Layer→fixed-sigma vertical remap (terrain-following output grid).
      !! `levels(:)` are cumulative sigma fractions (0..1, shallow→deep); the
      !! m-th output cell spans `[levels(m-1), levels(m)] * col_h`.  Same
      !! conservative donor-cell overlap as `remap_layer_to_z`; see its
      !! docstring for the intensive/extensive contract.
      class(*), intent(in) :: state_handle
      real(wp), intent(in) :: levels(:)
      real(wp), intent(in) :: layer_buf(:, :, :)
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in) :: is_extensive

      select type (state => state_handle)
      class is (ocean_state_t)
         call remap_layer_to_vcoord_impl(DIAG_VGRID_SIGMA, state%multilayer%h_layer, &
                                         levels, layer_buf, output_buf, is_extensive, &
                                         diag_remap_method, diag_mask_vanished, &
                                         DIAG_MISSING_VALUE)
      end select
   end subroutine remap_layer_to_sigma

   subroutine remap_layer_to_zstar(state_handle, levels, layer_buf, output_buf, is_extensive)
      !! Layer→fixed-z* vertical remap (SSH-tracking stretched-depth output
      !! grid).  `levels(:)` are reference interface depths (m, positive-down,
      !! shallow→deep); the deepest is the reference total depth H_ref and the
      !! per-column grid is stretched by `col_h / H_ref`.  Identical to
      !! `remap_layer_to_sigma` under uniform levels — supply a non-uniform
      !! (fine-near-surface) reference for it to differ.  Same conservative
      !! donor-cell overlap; see `remap_layer_to_z` for the intensive/extensive
      !! contract.
      class(*), intent(in) :: state_handle
      real(wp), intent(in) :: levels(:)
      real(wp), intent(in) :: layer_buf(:, :, :)
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in) :: is_extensive

      select type (state => state_handle)
      class is (ocean_state_t)
         call remap_layer_to_vcoord_impl(DIAG_VGRID_ZSTAR, state%multilayer%h_layer, &
                                         levels, layer_buf, output_buf, is_extensive, &
                                         diag_remap_method, diag_mask_vanished, &
                                         DIAG_MISSING_VALUE)
      end select
   end subroutine remap_layer_to_zstar

   subroutine coord_remap_proc(coord, remap)
      !! Return the default remap procedure pointer for an output vgrid.
      !! Used by the registration path to attach the right conservative
      !! remap when a diagnostic selects a non-layer output coordinate.
      !! Null for LAYER / unknown (no remap needed).
      integer, intent(in) :: coord
      procedure(diag_remap_proc), pointer, intent(out) :: remap
      select case (coord)
      case (DIAG_VGRID_Z_FIXED)
         remap => remap_layer_to_z
      case (DIAG_VGRID_SIGMA)
         remap => remap_layer_to_sigma
      case (DIAG_VGRID_ZSTAR)
         remap => remap_layer_to_zstar
      case (DIAG_VGRID_DENSITY)
         remap => remap_layer_to_density
      case default
         remap => null()
      end select
   end subroutine coord_remap_proc

   pure subroutine remap_layer_to_vcoord_impl(coord_type, h_layer, levels, layer_buf, &
                                              output_buf, is_extensive, method, &
                                              mask_vanished, missing)
      !! Per-column conservative layer→output-coordinate remap as a
      !! `do concurrent` over (j, i).  Builds the source column TOP-DOWN
      !! (work index 1 = surface = state k=nz) and the target cells from the
      !! `levels` interface positions (implicit 0 surface), clips both to the
      !! column total, then runs the donor-cell overlap integral
      !! (`remap_column`, scheme `method`) on a common `n = max(nz, nz_out)`
      !! padded partition.  Intensive remaps the value directly; extensive
      !! divides in / multiplies out by thickness so the column integral
      !! redistributes (sum preserved).
      !!
      !! `coord_type` selects how `levels(m)` maps to a target interface
      !! depth (a per-column scale hoisted out of the inner loop):
      !!   * `Z_FIXED`: `levels` are absolute depths (m, positive-down) —
      !!     scale 1 (the bit-identical legacy fixed-z path).
      !!   * `SIGMA`:   `levels` are cumulative fractions (0..1) — depth =
      !!     `levels(m) * col_h` (terrain-following).
      !!   * `ZSTAR`:   `levels` are reference depths (deepest = H_ref) —
      !!     depth = `levels(m) * col_h / H_ref` (SSH-tracking: when
      !!     col_h == H_ref the grid is the reference grid).  Identical to
      !!     SIGMA under uniform levels; differs only with a non-uniform
      !!     (e.g. fine-near-surface) reference.
      !!
      !! Fixed-size `NZ_STACK_MAX` stack locals via `local(...)` — automatic
      !! arrays sized from a dummy crash NVHPC stdpar device codegen.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! size() used to derive loop bounds from the actual buffer dimensions.
      integer, intent(in)    :: coord_type
      real(wp), intent(in)    :: h_layer(:, :, :)
      real(wp), intent(in)    :: levels(:)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: layer_buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in)    :: is_extensive
      integer, intent(in)    :: method
      logical, intent(in)    :: mask_vanished
      real(wp), intent(in)    :: missing
      integer :: i, j, k, m, nx, ny, nz, nz_out, n
      real(wp) :: dz_old(NZ_STACK_MAX), dz_new(NZ_STACK_MAX)
      real(wp) :: q_old(NZ_STACK_MAX), q_new(NZ_STACK_MAX)
      real(wp) :: col_h, zf_prev, zf, dz, lvl_scale, href
      nx = size(layer_buf, 1)
      ny = size(layer_buf, 2)
      nz = size(layer_buf, 3)
      nz_out = size(levels)
      n = max(nz, nz_out)
      href = levels(nz_out)   ! ZSTAR reference total depth (deepest interface)
      do concurrent(j=1:ny, i=1:nx) &
         local(dz_old, dz_new, q_old, q_new, k, m, col_h, zf_prev, zf, dz, lvl_scale)
         ! --- source column TOP-DOWN: work index k = state index nz-k+1 ---
         col_h = 0.0_wp
         do k = 1, nz
            dz_old(k) = h_layer(i, j, nz - k + 1)
            q_old(k) = layer_buf(i, j, nz - k + 1)
            col_h = col_h + dz_old(k)
         end do
         if (is_extensive) then
            do k = 1, nz
               if (dz_old(k) > 1.0e-12_wp) then
                  q_old(k) = q_old(k)/dz_old(k)   ! integral -> concentration
               else
                  q_old(k) = 0.0_wp
               end if
            end do
         end if
         ! pad the source to n with zero-thickness layers at the deep end
         do k = nz + 1, n
            dz_old(k) = 0.0_wp
            q_old(k) = 0.0_wp
         end do
         ! --- per-column level->depth scale (loop-invariant, hoisted) ---
         select case (coord_type)
         case (DIAG_VGRID_SIGMA)
            lvl_scale = col_h
         case (DIAG_VGRID_ZSTAR)
            if (href > 1.0e-12_wp) then
               lvl_scale = col_h/href
            else
               lvl_scale = 1.0_wp
            end if
         case default   ! DIAG_VGRID_Z_FIXED: levels are absolute depths
            lvl_scale = 1.0_wp
         end select
         ! --- target cells from scaled interfaces (implicit 0 surface),
         !     clipped to the column total H so the totals match exactly ---
         zf_prev = 0.0_wp
         do m = 1, nz_out
            zf = min(max(levels(m)*lvl_scale, 0.0_wp), col_h)
            dz = zf - zf_prev
            if (dz < 0.0_wp) dz = 0.0_wp
            dz_new(m) = dz
            zf_prev = zf
         end do
         do m = nz_out + 1, n
            dz_new(m) = 0.0_wp
         end do
         call remap_column(method, n, dz_old, dz_new, q_old, q_new)
         do m = 1, nz_out
            if (mask_vanished .and. dz_new(m) <= VANISHED_TARGET_FLOOR) then
               ! Target cell overlaps no water (below-bottom / pinched-out):
               ! emit the missing sentinel instead of a misleading 0.
               output_buf(i, j, m) = missing
            else if (is_extensive) then
               output_buf(i, j, m) = q_new(m)*dz_new(m)
            else
               output_buf(i, j, m) = q_new(m)
            end if
         end do
      end do
   end subroutine remap_layer_to_vcoord_impl

   subroutine remap_layer_to_density(state_handle, z_out, layer_buf, output_buf, is_extensive)
      !! Conservative layer→DENSITY-space remap (`DIAG_VGRID_DENSITY`).
      !! `z_out(:)` = monotone-increasing target potential DENSITIES (kg/m³).
      !! Per column: layer potential density via device EOS at the diag
      !! reference pressure, invert the profile to target-interface depths
      !! (`invert_density_targets`), then remap.  Lightest target → surface.
      !! `is_extensive=.false.` ⇒ INTENSIVE (weighted average); `.true.` ⇒
      !! EXTENSIVE (column integral redistributed across bins, Σ preserved).
      !! Public only for the unit-test suite.
      class(*), intent(in) :: state_handle
      real(wp), intent(in) :: z_out(:)
      real(wp), intent(in) :: layer_buf(:, :, :)
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in) :: is_extensive
      integer :: it, is_
      real(wp) :: pref
      select type (state => state_handle)
      class is (ocean_state_t)
         it = state%multilayer%idx_temperature
         is_ = state%multilayer%idx_salinity
         pref = state%vcoord%rho_ref_pressure
         if (it <= 0 .or. is_ <= 0) then
            call fill_zero_impl(output_buf)
            return
         end if
         call remap_layer_to_density_impl(state%multilayer%h_layer, &
                                          state%multilayer%tracers(it)%hTr, &
                                          state%multilayer%tracers(is_)%hTr, &
                                          state%eos, pref, z_out, &
                                          layer_buf, output_buf, is_extensive, &
                                          diag_remap_method)
      end select
   end subroutine remap_layer_to_density

   pure subroutine remap_layer_to_density_impl(h_layer, hT, hS, eos, rho_ref_p, &
                                               rho_tgt, layer_buf, output_buf, &
                                               is_extensive, method)
      !! Per-column density-space remap, `do concurrent` over (j, i).
      !! Source column TOP-DOWN + layer potential density (EOS at
      !! `rho_ref_p`), invert profile to interface depths
      !! (`invert_density_targets`), then donor-cell remap.  Cells outside
      !! the column density range read 0.  Intensive/extensive as the z-remap.
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded);
      ! size() used to derive loop bounds from the actual buffer dimensions.
      real(wp), intent(in)    :: h_layer(:, :, :)
      real(wp), intent(in)    :: hT(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: hS(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      type(eos_t), intent(in) :: eos
      real(wp), intent(in)    :: rho_ref_p
      real(wp), intent(in)    :: rho_tgt(:)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(in)    :: layer_buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      real(wp), intent(inout) :: output_buf(:, :, :)
      logical, intent(in)    :: is_extensive
      integer, intent(in)    :: method
      integer :: i, j, k, m, nx, ny, nz, n_bin, n
      real(wp) :: dz_old(NZ_STACK_MAX), dz_new(NZ_STACK_MAX)
      real(wp) :: q_old(NZ_STACK_MAX), q_new(NZ_STACK_MAX)
      ! z_iface holds n_bin+2 entries (`invert_density_targets` writes
      ! z_new(1:n_int+2) with n_int = n_bin, and :798 reads z_iface(n_bin+2)),
      ! so it needs NZ_STACK_MAX+2 — at NZ_STACK_MAX+1 a density diagnostic
      ! with n_rho_out == NZ_STACK_MAX (which the output-level guard permits)
      ! wrote and read one element past the end.
      real(wp) :: rhoc(NZ_STACK_MAX), z_iface(NZ_STACK_MAX + 2)
      real(wp) :: rho_tgt_c(NZ_STACK_MAX)
      real(wp) :: hh, tt, ss
      nx = size(layer_buf, 1)
      ny = size(layer_buf, 2)
      nz = size(layer_buf, 3)
      n_bin = size(rho_tgt)      ! number of density bins (output cells)
      n = max(nz, n_bin)
      do concurrent(j=1:ny, i=1:nx) &
         local(dz_old, dz_new, q_old, q_new, rhoc, z_iface, rho_tgt_c, &
               k, m, hh, tt, ss)
         ! --- source column TOP-DOWN + layer potential density ---
         do k = 1, nz
            hh = h_layer(i, j, nz - k + 1)
            dz_old(k) = hh
            q_old(k) = layer_buf(i, j, nz - k + 1)
            if (hh > 1.0e-12_wp) then
               tt = hT(i, j, nz - k + 1)/hh
               ss = hS(i, j, nz - k + 1)/hh
            else
               tt = 0.0_wp
               ss = 0.0_wp
            end if
            rhoc(k) = eos_density_point(eos, tt, ss, rho_ref_p)
         end do
         if (is_extensive) then
            do k = 1, nz
               if (dz_old(k) > 1.0e-12_wp) then
                  q_old(k) = q_old(k)/dz_old(k)
               else
                  q_old(k) = 0.0_wp
               end if
            end do
         end if
         do k = nz + 1, n
            dz_old(k) = 0.0_wp
            q_old(k) = 0.0_wp
         end do
         ! --- invert instantaneous density to target-interface depths ---
         ! n_bin output cells; the LAST cell absorbs to the bed so
         ! Σ dz_new == H exactly (extensive conservation, densest mass kept).
         ! Copy assumed-shape rho_tgt(:) into a contiguous fixed-size local
         ! first: passing the assumed-shape actual to the explicit-shape
         ! dummy makes flang emit host copy helpers absent in the AMD device
         ! runtime → offload link fails.  Element reads off the descriptor OK.
         do m = 1, n_bin
            rho_tgt_c(m) = rho_tgt(m)
         end do
         call invert_density_targets(nz, dz_old, rhoc, n_bin, rho_tgt_c, z_iface)
         do m = 1, n_bin - 1
            dz_new(m) = z_iface(m + 1) - z_iface(m)
            if (dz_new(m) < 0.0_wp) dz_new(m) = 0.0_wp
         end do
         dz_new(n_bin) = z_iface(n_bin + 2) - z_iface(n_bin)
         if (dz_new(n_bin) < 0.0_wp) dz_new(n_bin) = 0.0_wp
         do m = n_bin + 1, n
            dz_new(m) = 0.0_wp
         end do
         call remap_column(method, n, dz_old, dz_new, q_old, q_new)
         do m = 1, n_bin
            if (is_extensive) then
               output_buf(i, j, m) = q_new(m)*dz_new(m)
            else
               output_buf(i, j, m) = q_new(m)
            end if
         end do
      end do
   end subroutine remap_layer_to_density_impl

   subroutine fill_mld_epbl(state_handle, buf)
      !! EPBL active-mixing-layer depth (m) into the k=1 plane.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_2d_impl(state%epbl%mld, buf)
      end select
   end subroutine fill_mld_epbl

   pure subroutine fill_2d_impl(field, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: field(:, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, nx, ny
      nx = size(buf, 1)
      ny = size(buf, 2)
      do concurrent(j=1:ny, i=1:nx)
         buf(i, j, 1) = field(i, j)
      end do
   end subroutine fill_2d_impl

   subroutine fill_kd_epbl(state_handle, buf)
      !! EPBL interface diffusivity (m^2/s).  Buffer is layer-shaped
      !! (nx, ny, nz); we emit the value at the BOTTOM interface of
      !! each layer (kd_int(:, :, k) convention), losing only the
      !! identically-zero surface interface.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_kd_epbl_impl(state%epbl%kd_int, buf)
      end select
   end subroutine fill_kd_epbl

   subroutine fill_kd_kshear(state_handle, buf)
      !! Kappa-shear interface diffusivity (m^2/s).  Buffer is
      !! layer-shaped (nx, ny, nz); we emit the value at the BOTTOM
      !! interface of each layer (kd_int(:, :, k) convention), losing
      !! only the identically-zero surface interface.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_kd_epbl_impl(state%kshear%kd_int, buf)
      end select
   end subroutine fill_kd_kshear

   pure subroutine fill_kd_epbl_impl(kd_int, buf)
      ! assumed-shape-ok: diag fill — fires once per output frame (cadence-bounded).
      real(wp), intent(in)    :: kd_int(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)  ! assumed-shape-ok: diag fill — cadence-bounded
      integer :: i, j, k, nx, ny, nz
      nx = size(buf, 1)
      ny = size(buf, 2)
      nz = size(buf, 3)
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = kd_int(i, j, k)
      end do
   end subroutine fill_kd_epbl_impl

   subroutine register_default_diags(state, dt_out, specs, default_coord)
      !! Register the canonical ocean diagnostic variable set.  Called once
      !! from setup (driver init or test) after `state%diag%init`.  When
      !! `specs` (parsed `&ocean_diag_nml diags`) is present, each canonical
      !! diagnostic consults it: a `:off` entry skips registration entirely,
      !! and `:cadence` / `:op` / `:coord` attributes override the defaults.
      !! `default_coord` (the global `&ocean_diag_nml vgrid`, default LAYER)
      !! sets the output vgrid for layered diagnostics absent a per-diag
      !! `:coord`.  Absent `specs` + LAYER default => the canonical defaults,
      !! bit-identical to the legacy behaviour.  The canonical set also
      !! gains `ice_conc` / `ice_thick` when `&ocean_ice_nml enable`.
      type(ocean_state_t), intent(inout), target :: state
      real(wp), intent(in), optional :: dt_out
      type(diag_spec_t), intent(in), optional :: specs(:)
      integer, intent(in), optional :: default_coord
      integer :: nz, dcoord
      real(wp) :: dtout

      dtout = 3600.0_wp
      if (present(dt_out)) dtout = dt_out
      dcoord = DIAG_VGRID_LAYER
      if (present(default_coord)) dcoord = default_coord
      nz = state%multilayer%nz_ml

      call register_one_canonical(state, specs, dcoord, dtout, "SSH", "m", fill_ssh, &
                                  1, DIAG_OP_INSTANT, "sea_surface_height_above_geoid", &
                                  "sea_surface_height")
      if (state%multilayer%idx_temperature > 0) then
         call register_one_canonical(state, specs, dcoord, dtout, "temperature", "degC", &
                                     fill_temperature, nz, DIAG_OP_MEAN, &
                                     "sea_water_potential_temperature", &
                                     "sea_water_potential_temperature")
      end if
      if (state%multilayer%idx_salinity > 0) then
         call register_one_canonical(state, specs, dcoord, dtout, "salinity", "psu", &
                                     fill_salinity, nz, DIAG_OP_MEAN, &
                                     "sea_water_salinity", "sea_water_salinity")
      end if
      if (state%multilayer%idx_age > 0) then
         call register_one_canonical(state, specs, dcoord, dtout, "age", "s", fill_age, &
                                     nz, DIAG_OP_MEAN, "ideal_age_of_sea_water", &
                                     "age_of_sea_water")
      end if
      if (state%multilayer%idx_pseudo_salt > 0) then
         call register_one_canonical(state, specs, dcoord, dtout, "pseudo_salt", "psu", &
                                     fill_pseudo_salt, nz, DIAG_OP_MEAN, &
                                     "pseudo_salt_passive_tracer")
         call register_one_canonical(state, specs, dcoord, dtout, "pseudo_salt_diff", "psu", &
                                     fill_pseudo_salt_diff, nz, DIAG_OP_MEAN, &
                                     "difference_between_pseudo_salt_and_salt")
      end if
      call register_one_canonical(state, specs, dcoord, dtout, "u", "m s-1", fill_u_centre, &
                                  nz, DIAG_OP_MEAN, "eastward_velocity_at_cell_centre")
      call register_one_canonical(state, specs, dcoord, dtout, "v", "m s-1", fill_v_centre, &
                                  nz, DIAG_OP_MEAN, "northward_velocity_at_cell_centre")
      call register_one_canonical(state, specs, dcoord, dtout, "KE", "m2 s-2", fill_ke, &
                                  nz, DIAG_OP_MEAN, "kinetic_energy_per_unit_mass")
      if (state%epbl%enable) then
         call register_one_canonical(state, specs, dcoord, dtout, "MLD_EPBL", "m", &
                                     fill_mld_epbl, 1, DIAG_OP_MEAN, &
                                     "epbl_active_mixing_layer_depth")
         call register_one_canonical(state, specs, dcoord, dtout, "Kd_EPBL", "m2 s-1", &
                                     fill_kd_epbl, nz, DIAG_OP_MEAN, &
                                     "epbl_diffusivity_at_layer_bottom_interface")
      end if
      if (state%kshear%enable) then
         call register_one_canonical(state, specs, dcoord, dtout, "Kd_KSHEAR", "m2 s-1", &
                                     fill_kd_kshear, nz, DIAG_OP_MEAN, &
                                     "kappa_shear_diffusivity_at_layer_bottom_interface")
      end if
      if (state%ice%enable) then
         call register_one_canonical(state, specs, dcoord, dtout, "ice_conc", "1", &
                                     fill_ice_conc, 1, DIAG_OP_MEAN, &
                                     "sea_ice_area_fraction", "sea_ice_area_fraction")
         ! ice_thick carries NO CF standard_name: the CF `sea_ice_thickness`
         ! conventionally means actual floe thickness, but this is the
         ! grid-mean (effective) thickness mice/ICE_RHO_ICE — so only the
         ! long_name identifies it (omit standard_name, as u/v/KE do).
         call register_one_canonical(state, specs, dcoord, dtout, "ice_thick", "m", &
                                     fill_ice_thick, 1, DIAG_OP_MEAN, &
                                     "grid_mean_sea_ice_thickness")
      end if
   end subroutine register_default_diags

   subroutine register_one_canonical(state, specs, default_coord, dtout, name, units, &
                                     fill, n3, def_op, long_name, standard_name)
      !! Register one canonical diagnostic, applying its `specs` entry:
      !! skip if `:off`; override cadence / time-op; pick the output vgrid
      !! (`:coord` override, else `default_coord`, else LAYER) and attach the
      !! matching conservative remap.  2D diagnostics (`n3 <= 1`) ignore any
      !! coord request (no vertical to remap).  All canonical diagnostics are
      !! INTENSIVE (thickness-weighted average on remap).
      type(ocean_state_t), intent(inout), target :: state
      type(diag_spec_t), intent(in), optional :: specs(:)
      integer, intent(in) :: default_coord, n3, def_op
      real(wp), intent(in) :: dtout
      character(len=*), intent(in) :: name, units, long_name
      character(len=*), intent(in), optional :: standard_name
      procedure(diag_fill_proc) :: fill
      integer :: nx, ny, op, coord
      real(wp) :: dto
      logical :: skip
      procedure(diag_remap_proc), pointer :: remap

      op = def_op
      dto = dtout
      coord = default_coord
      call resolve_canonical_spec(specs, name, op, dto, coord, skip)
      if (skip) return
      ! 2D fields have no vertical axis to remap onto.
      if (n3 <= 1) coord = DIAG_VGRID_LAYER
      call coord_remap_proc(coord, remap)
      nx = size(state%barotropic%h, 1)
      ny = size(state%barotropic%h, 2)
      if (associated(remap)) then
         call state%diag%register(name, units=units, fill=fill, n1=nx, n2=ny, n3=n3, &
                                  long_name=long_name, standard_name=standard_name, &
                                  time_op=op, dt_out=dto, output_vgrid=coord, &
                                  remap=remap, is_extensive=.false., &
                                  has_missing=(diag_mask_vanished .and. &
                                               coord /= DIAG_VGRID_DENSITY))
      else
         call state%diag%register(name, units=units, fill=fill, n1=nx, n2=ny, n3=n3, &
                                  long_name=long_name, standard_name=standard_name, &
                                  time_op=op, dt_out=dto)
      end if
   end subroutine register_one_canonical

   pure subroutine resolve_canonical_spec(specs, name, time_op, dt_out, coord, skip)
      !! Look up canonical diagnostic `name` in the parsed `specs` and apply
      !! its overrides: `skip=.true.` if the entry is `:off`; otherwise
      !! `time_op` / `dt_out` / `coord` are overwritten when the entry sets
      !! them.  No matching entry (or absent `specs`) leaves the caller's
      !! defaults untouched.
      type(diag_spec_t), intent(in), optional :: specs(:)
      character(len=*), intent(in) :: name
      integer, intent(inout) :: time_op
      real(wp), intent(inout) :: dt_out
      integer, intent(inout) :: coord
      logical, intent(out) :: skip
      integer :: i

      skip = .false.
      if (.not. present(specs)) return
      do i = 1, size(specs)
         if (trim(specs(i)%name) == trim(name)) then
            if (specs(i)%off) then
               skip = .true.
               return
            end if
            if (specs(i)%time_op /= DIAG_OP_UNSET) time_op = specs(i)%time_op
            if (specs(i)%dt_out > 0.0_wp) dt_out = specs(i)%dt_out
            if (specs(i)%coord /= DIAG_COORD_UNSET) coord = specs(i)%coord
            return
         end if
      end do
   end subroutine resolve_canonical_spec

   pure function is_canonical_diag_name(name) result(yes)
      !! `.true.` if `name` is one of the canonical default diagnostics
      !! registered by `register_default_diags` (including the
      !! conditionally-registered EPBL / kappa-shear / tracer / sea-ice
      !! diags).  Used by the selection orchestrator to route a spec entry
      !! to either the canonical-override path or the derived-catalog
      !! registration path.
      character(len=*), intent(in) :: name
      logical :: yes
      integer :: i

      yes = .false.
      do i = 1, N_CANONICAL_DIAGS
         if (trim(CANONICAL_DIAG_NAMES(i)) == trim(name)) then
            yes = .true.
            return
         end if
      end do
   end function is_canonical_diag_name

   pure function canonical_diag_catalog_size() result(n)
      !! Number of names in the canonical-diagnostic catalog (the static
      !! set `register_default_diags` MAY register — some entries are
      !! gated, see `canonical_diag_gate_hint`). Bounds for
      !! `canonical_diag_catalog_name`'s index argument. Public for the
      !! P7 discoverability C ABI (`rdb_ocean_canonical_catalog_size`)
      !! and the unit-test suite.
      integer :: n
      n = N_CANONICAL_DIAGS
   end function canonical_diag_catalog_size

   pure function canonical_diag_catalog_name(i) result(name)
      !! Name of canonical-catalog entry `i` (1-based). Public for the P7
      !! discoverability C ABI and the unit-test suite.
      integer, intent(in) :: i
      character(len=16) :: name
      name = CANONICAL_DIAG_NAMES(i)
   end function canonical_diag_catalog_name

   pure function canonical_diag_gate_hint(name) result(hint)
      !! The namelist gate whose closure suppresses canonical diagnostic
      !! `name`, for the "you asked for this and did not get it" warning.
      !! `""` for the four unconditional diagnostics (SSH/u/v/KE) — a hint
      !! of `""` means "this one should have registered; that is a bug,
      !! not a config".  MUST be kept in lock-step with
      !! `register_default_diags`'s gates; the count of non-empty hints
      !! must equal the number of conditionally-registered canonical
      !! diagnostics (`test_ocean_diag/diag_gate_hint_covers_every_gate` is
      !! the lock — any PR that adds a new gated canonical diagnostic must
      !! add its hint here or that test fails).  This duplicates the gate
      !! NAMES only, never the gate LOGIC — the single `if` in
      !! `register_default_diags` remains the one place the gate is
      !! evaluated.
      character(len=*), intent(in) :: name
      character(len=64) :: hint

      select case (trim(name))
      case ("temperature", "salinity")
         hint = "&ocean_thermo_nml enable_thermodynamics"
      case ("age")
         hint = "&ocean_tracers_nml enable_ideal_age"
      case ("MLD_EPBL", "Kd_EPBL")
         hint = "&ocean_epbl_nml enable"
      case ("Kd_KSHEAR")
         hint = "&ocean_kappa_shear_nml enable"
      case ("ice_conc", "ice_thick")
         hint = "&ocean_ice_nml enable"
      case default
         hint = ""
      end select
   end function canonical_diag_gate_hint

end module rdb_ocean_diag_fills
