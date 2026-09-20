!! Z-level T/S initial-condition overlay (capability A2).
module rdb_ocean_z_init
   !! Puts a GEOPOTENTIAL T(z)/S(z) profile onto the seeded layer
   !! centres (`&ocean_zinit_nml`) and writes the multilayer tracer slots
   !! as `hTr = value * h_layer`.  Two sources, `&ocean_zinit_nml source`:
   !!
   !!   * `"file"` (default) — reads pre-regridded T + S from a
   !!     model-grid NetCDF and interpolates linearly in depth.
   !!     File is PRE-REGRIDDED to `nx_phys x ny_phys` — no in-core
   !!     horizontal interpolation.  Vertical interp is
   !!     point-linear-in-depth, constant extrapolation beyond the source
   !!     range.  Dry columns (`wet_mask <= 0`) get the namelist
   !!     `land_fill_t/_s` constants.  Interior cells only — ghosts keep
   !!     whatever the analytical IC seeded, because the file carries no
   !!     data for them.
   !!   * `"linear"` — an ANALYTIC affine profile,
   !!     `T(z) = lin_t_ref + lin_dt_dz*z` (and the salinity twin) with
   !!     `z` the GEOPOTENTIAL height, positive UP, zero at the `z = 0`
   !!     datum, i.e. `z = -z_ctr`.  Needs no file, so it fills the FULL
   !!     array INCLUDING ghosts and land columns (the analytic formula
   !!     is defined everywhere — same reasoning as the formula
   !!     bathymetry setters, and it keeps a wall-adjacent EOS
   !!     evaluation off the `rho_0` fallback).  This is the profile an
   !!     idealised ice-shelf cavity needs: `&tracer_nml T_init_surface`
   !!     / `T_init_bottom` are linear in LAYER INDEX, which under a
   !!     terrain-following coordinate with a SLOPING lid tilts the
   !!     isopycnals with the coordinate and is therefore NOT a state of
   !!     rest.
   !!
   !! ### Depth is measured from `z = 0`, not from the column top
   !!
   !! `build_z_ctr` takes the geopotential depth of the column TOP
   !! (`z_top`, positive-down) and adds the thickness above each layer
   !! centre.  In the open ocean `z_top = 0` and this is the free-surface
   !! datum.  Under an ice shelf the water column starts `z_draft` metres
   !! down, so `z_top = metrics%z_draft(i,j)`; without it a geopotential
   !! profile lands `z_draft` metres too shallow on every shelf column
   !! and the isopycnals tilt with the ice base.  `z_draft` reaches here
   !! as an OPTIONAL ghosted `(nx_total, ny_total)` argument: absent (the
   !! no-cavity path) is `z_top = 0` and bit-identical to the original.
   !!
   !! Runs at seed time, HOST-side, BEFORE `ocean_state_enter_data`, so
   !! it uses plain host `do` loops (NOT `do concurrent`); the later
   !! enter_data copies the host-resident arrays up.  The module is only
   !! compiled with `RDB_ENABLE_NETCDF=ON` (the file reader needs it), so
   !! the `"linear"` source currently inherits that build requirement
   !! even though it opens nothing.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: ocean_zinit_config_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_var_3d, nc_get_var_1d
   use netcdf, only: nf90_noerr, nf90_inq_varid, nf90_inquire_variable, &
                     nf90_inquire_dimension
   use pic_logger, only: logger => global_logger
   use rdb_error_ring, only: fail
   use pic_strings, only: to_string
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IO
   implicit none
   private

   public :: seed_ts_from_zfile
   public :: seed_ts_linear_z
   public :: interp_column_linear_z
   public :: build_z_ctr
   public :: zinit_dims_ok

contains

   subroutine seed_ts_from_zfile(ms, grid, cfg, ierr, z_draft)
      !! Read T/S from the configured z-level NetCDF and overwrite the
      !! tracer slots with the depth-interpolated profiles.  Wet columns
      !! interpolate; dry columns get the namelist land-fill constants.
      !! Must be called AFTER bathymetry, the uniform `h_layer` seed, and
      !! `wet_mask` are in place.  Takes the multilayer slot directly (not
      !! the full ocean state) to avoid a circular dependency.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      type(ocean_zinit_config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a malformed/mismatched z-level IC file when
         !! present; absent behaves as today (`error stop`).
      real(wp), intent(in), optional :: z_draft(:, :)
         !! Ice-base depth (m, positive-down, `>= 0`), FULL ghosted
         !! `(nx_total, ny_total)` shape — `metrics%z_draft`.  The
         !! geopotential depth of every layer centre is measured from
         !! `z = 0`, so the column top sits at `z_draft` under a shelf.
         !! Absent ⇒ `z_top = 0` (open ocean), bit-identical to the
         !! pre-cavity behaviour.

      integer :: ncid, ng, nx, ny, nz_ml, nz_src, local_ierr
      integer :: idx_S, idx_T, i, j, k
      integer :: t_varid, s_varid, z_varid
      real(wp), allocatable :: t_src(:, :, :), s_src(:, :, :), z_src(:)
      real(wp) :: z_ctr(ms%nz_ml)
      real(wp) :: t_col(ms%nz_ml), s_col(ms%nz_ml)
      real(wp) :: z_top
      logical :: needs_transpose, have_draft

      have_draft = present(z_draft)
      if (have_draft) then
         if (.not. draft_shape_ok(z_draft, ms)) then
            call fail("ocean_zinit: z_draft must be the FULL ghosted "// &
                      "(nx_total, ny_total) array, matching h_layer", &
                      ierr, OCEAN_STATUS_ERR_IO)
            return
         end if
      end if
      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys
      nz_ml = ms%nz_ml
      idx_S = ms%idx_salinity
      idx_T = ms%idx_temperature

      call logger%info("Loading z-level T/S IC from: "//trim(cfg%file))

      call nc_open_read(trim(cfg%file), ncid, ierr=local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr)) return

      ! Locate the three variables: namelist override then name-fallbacks.
      ! ierr threaded down ONLY when THIS routine's own ierr is present:
      ! otherwise find_var/read_dims must keep reaching their own
      ! `error stop` (specific text) rather than the generic wrapper
      ! message below (P0.1 review F2).
      if (present(ierr)) then
         call find_var(ncid, cfg%t_var, &
                       [character(len=16) :: "temp", "T", "temperature"], "temperature", t_varid, &
                       ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call find_var(ncid, cfg%t_var, &
                       [character(len=16) :: "temp", "T", "temperature"], "temperature", t_varid)
      end if
      if (present(ierr)) then
         call find_var(ncid, cfg%s_var, &
                       [character(len=16) :: "salt", "S", "salinity"], "salinity", s_varid, &
                       ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call find_var(ncid, cfg%s_var, &
                       [character(len=16) :: "salt", "S", "salinity"], "salinity", s_varid)
      end if
      if (present(ierr)) then
         call find_var(ncid, cfg%z_var, &
                       [character(len=16) :: "z_src", "z", "depth", "lev"], "source axis", z_varid, &
                       ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call find_var(ncid, cfg%z_var, &
                       [character(len=16) :: "z_src", "z", "depth", "lev"], "source axis", z_varid)
      end if

      ! Validate dims + determine the storage-order permutation against
      ! the model grid; error-stops on mismatch.
      if (present(ierr)) then
         call read_dims(ncid, t_varid, grid, nz_src, needs_transpose, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call read_dims(ncid, t_varid, grid, nz_src, needs_transpose)
      end if

      ! Read the source axis + assert monotonic increase (positive-down).
      allocate (z_src(nz_src))
      call nc_get_var_1d(ncid, z_varid, z_src, ierr=local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      do k = 2, nz_src
         if (z_src(k) <= z_src(k - 1)) then
            call nc_close(ncid)
            call fail("ocean_zinit: z_src not monotonically increasing at level "// &
                      to_string(k)//" ("//to_string(z_src(k))//" <= "// &
                      to_string(z_src(k - 1))//")", ierr, OCEAN_STATUS_ERR_IO)
            return
         end if
      end do

      ! Read T/S into model-grid (x, y, z) interior arrays, undoing the
      ! C/Fortran dimension reversal when the file is C-ordered (z, y, x).
      allocate (t_src(nx, ny, nz_src), s_src(nx, ny, nz_src))
      call read_field_xyz(ncid, t_varid, t_src, nx, ny, nz_src, needs_transpose, &
                          ierr=local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      call read_field_xyz(ncid, s_varid, s_src, nx, ny, nz_src, needs_transpose, &
                          ierr=local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      call nc_close(ncid)

      ! Host-side per-column interpolation.  Plain do loops: this runs
      ! before enter_data, so the arrays are host-resident.
      do j = 1, ny
         do i = 1, nx
            ! Layer-centre GEOPOTENTIAL depths (positive-down from z = 0)
            ! from h_layer, offset by the depth of the column top.
            ! Bottom-up: k=1 bed (deepest), k=nz_ml surface (shallowest).
            z_top = 0.0_wp
            if (have_draft) z_top = z_draft(ng + i, ng + j)
            call build_z_ctr(ms%h_layer(ng + i, ng + j, :), nz_ml, z_top, z_ctr)

            if (ms%wet_mask(ng + i, ng + j) <= 0.0_wp) then
               ! Dry column — fill with the namelist land-fill constants.
               t_col = cfg%land_fill_t
               s_col = cfg%land_fill_s
            else
               ! Interp in DEPTH SPACE so the source's ascending
               ! positive-down order and the model's bottom-up order
               ! never need reconciling.
               call interp_column_linear_z(z_src, t_src(i, j, :), nz_src, z_ctr, nz_ml, t_col)
               call interp_column_linear_z(z_src, s_src(i, j, :), nz_src, z_ctr, nz_ml, s_col)
            end if

            ! Roundabout tracer convention: hTr = value * h_layer.
            if (idx_T > 0) then
               do k = 1, nz_ml
                  ms%tracers(idx_T)%hTr(ng + i, ng + j, k) = &
                     t_col(k)*ms%h_layer(ng + i, ng + j, k)
               end do
            end if
            if (idx_S > 0) then
               do k = 1, nz_ml
                  ms%tracers(idx_S)%hTr(ng + i, ng + j, k) = &
                     s_col(k)*ms%h_layer(ng + i, ng + j, k)
               end do
            end if
         end do
      end do

      deallocate (t_src, s_src, z_src)

      call logger%info("ocean_zinit: seeded T/S from "//to_string(nz_src)// &
                       " source z-levels onto "//to_string(nz_ml)//" model layers.")
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_ts_from_zfile

   pure subroutine build_z_ctr(h_layer, nz_ml, z_top, z_ctr)
      !! Layer-centre GEOPOTENTIAL depths (positive-down from the `z = 0`
      !! datum) from a column of layer thicknesses.  Bottom-up: k=1 bed,
      !! k=nz_ml surface.
      !! `z_ctr(k) = z_top + sum_{k'=k+1..nz_ml} h(k') + 0.5*h(k)`.
      !!
      !! `z_top` is the depth of the TOP of the water column: `0` in the
      !! open ocean (the free-surface datum), `z_draft(i,j)` under an ice
      !! shelf.  Passing `0` reproduces the pre-cavity arithmetic
      !! bit-for-bit — `above` starts at `z_top` and the very first
      !! addition is `0 + 0.5*h`, the same expression as before.
      integer, intent(in) :: nz_ml
      real(wp), intent(in) :: h_layer(nz_ml)
      real(wp), intent(in) :: z_top
      real(wp), intent(out) :: z_ctr(nz_ml)

      integer :: k
      real(wp) :: above

      ! Accumulate from the surface (k=nz_ml) downward.  `above` holds
      ! the depth of the top of layer k: the column-top depth plus the
      ! total thickness of all layers shallower than k.
      above = z_top
      do k = nz_ml, 1, -1
         z_ctr(k) = above + 0.5_wp*h_layer(k)
         above = above + h_layer(k)
      end do
   end subroutine build_z_ctr

   pure function draft_shape_ok(z_draft, ms) result(ok)
      !! True iff `z_draft` is the FULL ghosted `(nx_total, ny_total)`
      !! array the seeders index with `(ng+i, ng+j)` — i.e. it matches
      !! `h_layer`'s horizontal extent.  Guards the `(1, 1)` placeholder
      !! `ocean_metrics_t` allocates when the cavity is off from ever
      !! being read as a field.
      real(wp), intent(in) :: z_draft(:, :)
      type(multilayer_state_t), intent(in) :: ms
      logical :: ok

      ok = (size(z_draft, 1) == size(ms%h_layer, 1)) .and. &
           (size(z_draft, 2) == size(ms%h_layer, 2))
   end function draft_shape_ok

   pure elemental function linear_in_z(v_ref, dv_dz, z_depth) result(v)
      !! Affine profile `v(z) = v_ref + dv_dz*z` evaluated at
      !! geopotential height `z = -z_depth`, i.e.
      !! `v = v_ref - dv_dz*z_depth`.
      !!
      !! `z` is positive UP with the origin at the `z = 0` datum, which
      !! is the convention `&ocean_ic_nml eady_dT_dz` already uses.  A
      !! STABLE thermal column therefore has `dv_dz > 0` for temperature
      !! (warm on top) and `dv_dz < 0` for salinity (salty at the bed).
      real(wp), intent(in) :: v_ref, dv_dz, z_depth
      real(wp) :: v

      v = v_ref - dv_dz*z_depth
   end function linear_in_z

   subroutine seed_ts_linear_z(ms, cfg, ierr, z_draft)
      !! Seed T and S from the ANALYTIC affine geopotential profiles
      !! `T(z) = lin_t_ref + lin_dt_dz*z`, `S(z) = lin_s_ref +
      !! lin_ds_dz*z` (`z` positive UP, zero at the `z = 0` datum),
      !! sampled at each layer centre's TRUE geopotential depth and
      !! written as `hTr = value * h_layer`.
      !!
      !! Opens no file, so — unlike the NetCDF path — it fills the FULL
      !! array INCLUDING the ghost rows and the dry columns: the formula
      !! is defined everywhere, and a ghost column left on some other
      !! profile is the wall-adjacent spurious-density-jump gotcha
      !! (CLAUDE.md, "formula bathymetry setters must fill ghost rows").
      !! There is no `land_fill_t/_s` on this path for the same reason.
      !!
      !! NOT `pure`: it logs one line reporting what it seeded, matching
      !! `seed_ts_from_zfile`.  The arithmetic itself is the `pure
      !! elemental` `linear_in_z`.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_zinit_config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a malformed `z_draft` when present; absent
         !! behaves like the rest of the seed path (`error stop`).
      real(wp), intent(in), optional :: z_draft(:, :)
         !! Ice-base depth (m, positive-down), FULL ghosted
         !! `(nx_total, ny_total)` shape.  Absent ⇒ open ocean.

      integer :: nxt, nyt, nz_ml, idx_S, idx_T, i, j, k
      real(wp) :: z_ctr(ms%nz_ml)
      real(wp) :: z_top
      logical :: have_draft

      have_draft = present(z_draft)
      if (have_draft) then
         if (.not. draft_shape_ok(z_draft, ms)) then
            call fail("ocean_zinit: z_draft must be the FULL ghosted "// &
                      "(nx_total, ny_total) array, matching h_layer", &
                      ierr, OCEAN_STATUS_ERR_IO)
            return
         end if
      end if

      nxt = size(ms%h_layer, 1)
      nyt = size(ms%h_layer, 2)
      nz_ml = ms%nz_ml
      idx_S = ms%idx_salinity
      idx_T = ms%idx_temperature

      ! Host-side, before enter_data: plain do loops, NOT do concurrent.
      do j = 1, nyt
         do i = 1, nxt
            z_top = 0.0_wp
            if (have_draft) z_top = z_draft(i, j)
            call build_z_ctr(ms%h_layer(i, j, :), nz_ml, z_top, z_ctr)
            if (idx_T > 0) then
               do k = 1, nz_ml
                  ms%tracers(idx_T)%hTr(i, j, k) = &
                     linear_in_z(cfg%lin_t_ref, cfg%lin_dt_dz, z_ctr(k))* &
                     ms%h_layer(i, j, k)
               end do
            end if
            if (idx_S > 0) then
               do k = 1, nz_ml
                  ms%tracers(idx_S)%hTr(i, j, k) = &
                     linear_in_z(cfg%lin_s_ref, cfg%lin_ds_dz, z_ctr(k))* &
                     ms%h_layer(i, j, k)
               end do
            end if
         end do
      end do

      if (have_draft) then
         call logger%info("ocean_zinit: seeded analytic linear-in-z T/S onto "// &
                          to_string(nz_ml)//" model layers, measuring the "// &
                          "geopotential depth from the ice base (cavity draft).")
      else
         call logger%info("ocean_zinit: seeded analytic linear-in-z T/S onto "// &
                          to_string(nz_ml)//" model layers.")
      end if
      call logger%info("ocean_zinit:   T = "//to_string(cfg%lin_t_ref)//" + "// &
                       to_string(cfg%lin_dt_dz)//"*z,  S = "// &
                       to_string(cfg%lin_s_ref)//" + "// &
                       to_string(cfg%lin_ds_dz)//"*z   (z positive UP from z = 0).")
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_ts_linear_z

   pure subroutine interp_column_linear_z(z_src, v_src, nz_src, z_ctr, nz_ml, v_out)
      !! Linear-in-depth interpolation of source profile `v_src` on
      !! ascending depths `z_src` onto target depths `z_ctr`, with
      !! CONSTANT extrapolation beyond the source range.  Comparisons in
      !! depth space so index orientations never need reconciling.
      integer, intent(in) :: nz_src, nz_ml
      real(wp), intent(in) :: z_src(nz_src), v_src(nz_src)
      real(wp), intent(in) :: z_ctr(nz_ml)
      real(wp), intent(out) :: v_out(nz_ml)

      integer :: k, m
      real(wp) :: zt, w

      do k = 1, nz_ml
         zt = z_ctr(k)
         if (zt <= z_src(1)) then
            ! Shallower than the shallowest source level — constant.
            v_out(k) = v_src(1)
         else if (zt >= z_src(nz_src)) then
            ! Deeper than the deepest source level — constant.
            v_out(k) = v_src(nz_src)
         else
            ! Locate the bracketing source levels [m, m+1] and lerp.
            do m = 1, nz_src - 1
               if (zt >= z_src(m) .and. zt <= z_src(m + 1)) then
                  w = (zt - z_src(m))/(z_src(m + 1) - z_src(m))
                  v_out(k) = (1.0_wp - w)*v_src(m) + w*v_src(m + 1)
                  exit
               end if
            end do
         end if
      end do
   end subroutine interp_column_linear_z

   subroutine find_var(ncid, override, fallbacks, label, varid, ierr)
      !! Resolve a variable id: try the namelist `override` name first
      !! (when non-blank), then the documented `fallbacks` in order.
      !! Error-stops with a descriptive message when none is found.
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: override
      character(len=*), intent(in) :: fallbacks(:)
      character(len=*), intent(in) :: label
      integer, intent(out) :: varid
      integer, intent(out), optional :: ierr
         !! Non-zero when no matching variable is found, when present;
         !! absent behaves as today (`error stop`).

      integer :: status, n
      character(len=512) :: tried

      if (len_trim(override) > 0) then
         status = nf90_inq_varid(ncid, trim(override), varid)
         if (status == nf90_noerr) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if
         call nc_close(ncid)
         call fail("ocean_zinit: "//trim(label)//" variable override '"// &
                   trim(override)//"' not found in file", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      tried = ""
      do n = 1, size(fallbacks)
         status = nf90_inq_varid(ncid, trim(fallbacks(n)), varid)
         if (status == nf90_noerr) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if
         tried = trim(tried)//" "//trim(fallbacks(n))
      end do

      call nc_close(ncid)
      call fail("ocean_zinit: no "//trim(label)//" variable found (tried:"// &
                trim(tried)//")", ierr, OCEAN_STATUS_ERR_IO)
   end subroutine find_var

   subroutine read_dims(ncid, varid, grid, nz_src, needs_transpose, ierr)
      !! Validate the T/S variable's dims against the model grid and
      !! detect whether the file is C-ordered (z, y, x) =>
      !! `needs_transpose = .true.` (by the first Fortran dim name).
      !! Returns source z-level count `nz_src`; error-stops on mismatch.
      integer, intent(in) :: ncid, varid
      type(hgrid_t), intent(in) :: grid
      integer, intent(out) :: nz_src
      logical, intent(out) :: needs_transpose
      integer, intent(out), optional :: ierr
         !! Non-zero on a dimensionality/shape mismatch when present;
         !! absent behaves as today (`error stop`).

      integer :: var_ndims
      integer :: var_dimids(3)
      integer :: d1_len, d2_len, d3_len
      integer :: local_ierr
      character(len=64) :: d1_name
      logical :: ok

      call nc_check(nf90_inquire_variable(ncid, varid, ndims=var_ndims, dimids=var_dimids), &
                    "querying z-level T/S variable", local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      if (var_ndims /= 3) then
         call nc_close(ncid)
         call fail("ocean_zinit: T/S variable must be 3D (x,y,z); got "// &
                   to_string(var_ndims)//" dims", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      call nc_check(nf90_inquire_dimension(ncid, var_dimids(1), name=d1_name, len=d1_len), &
                    "querying T/S dim 1", local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(2), len=d2_len), &
                    "querying T/S dim 2", local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(3), len=d3_len), &
                    "querying T/S dim 3", local_ierr)
      if (.not. zinit_io_ok(local_ierr, ierr, ncid)) return

      ! C-ordered files store temp(x,y,z); Fortran's NetCDF lib reverses
      ! that to temp(z,y,x) so dim1 is "z".  Fortran-ordered files keep
      ! dim1 == "x".
      needs_transpose = (trim(d1_name) == "z" .or. trim(d1_name) == "depth" .or. &
                         trim(d1_name) == "lev" .or. trim(d1_name) == "z_src")

      if (needs_transpose) then
         ! Fortran storage (z, y, x): physical lengths are d3=x, d2=y, d1=z.
         ok = (d3_len == grid%nx_phys .and. d2_len == grid%ny_phys)
         nz_src = d1_len
      else
         ! Fortran storage (x, y, z): physical lengths are d1=x, d2=y, d3=z.
         ok = (d1_len == grid%nx_phys .and. d2_len == grid%ny_phys)
         nz_src = d3_len
      end if

      if (.not. ok) then
         call nc_close(ncid)
         call fail("ocean_zinit: T/S grid mismatch: file horizontal dims do not match "// &
                   "simulation "//to_string(grid%nx_phys)//" x "//to_string(grid%ny_phys), ierr, OCEAN_STATUS_ERR_IO)
         return
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine read_dims

   subroutine read_field_xyz(ncid, varid, dst, nx, ny, nz_src, needs_transpose, ierr)
      !! Read a 3D T/S variable into a model-grid `(nx, ny, nz_src)`
      !! interior array, permuting from the file's Fortran storage order.
      integer, intent(in) :: ncid, varid, nx, ny, nz_src
      logical, intent(in) :: needs_transpose
      real(wp), intent(out) :: dst(nx, ny, nz_src)
      integer, intent(out), optional :: ierr
         !! Non-zero on a read failure when present; absent behaves as
         !! today (`error stop`).

      integer :: i, j, k
      real(wp), allocatable :: buf(:, :, :)

      if (needs_transpose) then
         ! File is (z, y, x) in Fortran order.
         allocate (buf(nz_src, ny, nx))
         call nc_get_var_3d(ncid, varid, buf, ierr)
         if (present(ierr)) then
            if (ierr /= 0) then
               deallocate (buf)
               return
            end if
         end if
         do k = 1, nz_src
            do j = 1, ny
               do i = 1, nx
                  dst(i, j, k) = buf(k, j, i)
               end do
            end do
         end do
      else
         ! File is (x, y, z) in Fortran order — direct read.
         call nc_get_var_3d(ncid, varid, dst, ierr)
      end if

      if (allocated(buf)) deallocate (buf)
   end subroutine read_field_xyz

   function zinit_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` reader calls in `seed_ts_from_zfile`/`read_dims` into the
      !! caller's `ierr` contract: `.true.` on success; on failure,
      !! `.false.` with `ierr = OCEAN_STATUS_ERR_IO` when `ierr` is
      !! present (closing `ncid` first, when given, so a mid-read failure
      !! does not leak the file handle), or `error stop`s with the SAME
      !! generic text `nc_check` itself would have used had `ierr` never
      !! been threaded through — keeps the legacy (no `ierr`) behaviour
      !! byte-identical while unblocking the `ierr`-present return path
      !! (F1/F2 of the P0.1 review).
      integer, intent(in) :: local_ierr
      integer, intent(out), optional :: ierr
      integer, intent(in), optional :: ncid
      logical :: ok

      integer :: discard_ierr

      ok = (local_ierr == 0)
      if (ok) return

      if (present(ierr)) then
         if (present(ncid)) call nc_close(ncid, ierr=discard_ierr)
         ierr = OCEAN_STATUS_ERR_IO
         return
      end if

      error stop "NetCDF operation failed"
   end function zinit_io_ok

   function zinit_dims_ok(filename, t_var, nx_phys, ny_phys) result(ok)
      !! Validator: true iff the temperature variable's horizontal dims
      !! match `(nx_phys, ny_phys)`.  Returns `.false.` on mismatch or
      !! I/O error; does NOT `error stop` — safe from test code.
      use netcdf, only: nf90_open, nf90_nowrite, nf90_close
      character(len=*), intent(in) :: filename, t_var
      integer, intent(in) :: nx_phys, ny_phys
      logical :: ok

      integer :: ncid, varid, ierr, var_ndims
      integer :: var_dimids(3)
      integer :: d1_len, d2_len, d3_len
      character(len=64) :: d1_name
      logical :: needs_transpose

      ok = .false.
      ierr = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (ierr /= nf90_noerr) return

      ! Chain queries through `ierr`: first failure short-circuits the
      ! rest; the file is closed once at the end (no early return).
      ierr = nf90_inq_varid(ncid, trim(t_var), varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_inquire_variable(ncid, varid, ndims=var_ndims, dimids=var_dimids)
      end if
      if (ierr == nf90_noerr .and. var_ndims /= 3) ierr = -1
      if (ierr == nf90_noerr) then
         ierr = nf90_inquire_dimension(ncid, var_dimids(1), name=d1_name, len=d1_len)
      end if
      if (ierr == nf90_noerr) then
         ierr = nf90_inquire_dimension(ncid, var_dimids(2), len=d2_len)
      end if
      if (ierr == nf90_noerr) then
         ierr = nf90_inquire_dimension(ncid, var_dimids(3), len=d3_len)
      end if

      if (ierr == nf90_noerr) then
         needs_transpose = (trim(d1_name) == "z" .or. trim(d1_name) == "depth" .or. &
                            trim(d1_name) == "lev" .or. trim(d1_name) == "z_src")
         if (needs_transpose) then
            ok = (d3_len == nx_phys .and. d2_len == ny_phys)
         else
            ok = (d1_len == nx_phys .and. d2_len == ny_phys)
         end if
      end if

      ierr = nf90_close(ncid)
   end function zinit_dims_ok

end module rdb_ocean_z_init
