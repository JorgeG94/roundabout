!! Analytical tests for the linear-in-z SALINITY initial condition
!! (`&tracer_nml S_init_surface` / `S_init_bottom`).
module test_ocean_salinity_ic
   !! `S_init_surface`/`S_init_bottom` were registered, type-validated
   !! and then read NOWHERE — the temperature siblings
   !! (`T_init_surface`/`T_init_bottom`) had a seed path and salinity did
   !! not, so a namelist asking for a haline-stratified column silently
   !! got a uniform one.  (Two shipped double-diffusion cases,
   !! `validation_examples/ocean/seamount/seamount_ddiff*.nml`, ask for
   !! exactly that.)  This suite pins the seed that closes the gap:
   !!
   !!   * **Endpoints** — `S_init_bottom` lands at `k = 1` (the BED) and
   !!     `S_init_surface` at `k = nz_ml` (the SURFACE), matching the
   !!     repo's bottom-up vertical convention.  Note the stable haline
   !!     polarity is the INVERSE of temperature's: dense/salty water
   !!     belongs at the bed.
   !!   * **Linearity** — the profile is linear in layer index, which
   !!     under the sigma-style `h_layer = b/nz_ml` seed is linear in
   !!     layer-centre depth on every column.
   !!   * **Sloping bathymetry** — on a spoon basin (depth varies by
   !!     cell) the endpoint values and the linearity still hold per
   !!     column; only the physical gradient dS/dz scales with 1/depth.
   !!     This is the sigma-like behaviour the temperature path already
   !!     had, verified here so the two can't drift apart.
   !!   * **Ghosts + land** — every cell of the full local array,
   !!     ghosts included, carries `hTr = S(k)*h_layer`; the seed writes
   !!     the same expression everywhere, exactly as the temperature
   !!     seed does, so a wall-adjacent EOS evaluation sees a real
   !!     density rather than the `rho_0` vanished-layer fallback.
   !!   * **Pseudo-salt** — with `&ocean_tracers_nml enable_pseudo_salt`
   !!     the verification tracer is seeded from the FINAL salinity, so
   !!     its day-zero deviation is exactly zero on a stratified IC too.
   !!   * **Bit-identity** — both knobs at their `0.0` default keeps the
   !!     uniform `initial_salinity` column, byte-identical.
   !!
   !! Host-only: `ocean_state_seed_from_cfg` runs BEFORE
   !! `ocean_state_enter_data` by driver contract (it is deliberately a
   !! plain host loop — see `seed_tracer_stratified_impl`), so there is
   !! no device-resident array in the loop and no `mem:separate` mapping
   !! to arrange.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

   public :: collect_ocean_salinity_ic_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 6
   integer, parameter :: NZ = 10

   real(wp), parameter :: MAX_DEPTH = 1000.0_wp
   real(wp), parameter :: EDGE_DEPTH = 200.0_wp
   real(wp), parameter :: SLOPE_SCALE = 8000.0_wp
   real(wp), parameter :: DX = 2000.0_wp

   ! Stable haline column: FRESH at the surface, SALTY at the bed.
   real(wp), parameter :: S_SURF = 34.0_wp
   real(wp), parameter :: S_BED = 36.75_wp
   real(wp), parameter :: S_UNIFORM = 35.0_wp

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_ocean_salinity_ic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("salinity_ic_endpoints_and_linearity", test_endpoints_linear), &
                  new_unittest("salinity_ic_sloping_bathymetry", test_sloping_bathymetry), &
                  new_unittest("salinity_ic_defaults_bit_identical", test_defaults_uniform), &
                  new_unittest("salinity_ic_pseudo_salt_seeded_from_S", test_pseudo_salt_seed), &
                  new_unittest("salinity_ic_one_sided_fails_loud", test_one_sided_fails) &
                  ]
   end subroutine collect_ocean_salinity_ic_tests

   subroutine base_cfg(cfg, topo)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: topo

      cfg%sim_type = "ocean"
      cfg%nx = NX_PHYS
      cfg%ny = NY_PHYS
      cfg%nghost = NGHOST
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%ocean%topo%topo_config = topo
      cfg%ocean%topo%max_depth = MAX_DEPTH
      cfg%ocean%topo%edge_depth = EDGE_DEPTH
      cfg%ocean%topo%slope_scale = SLOPE_SCALE
      cfg%initial_salinity = S_UNIFORM
      cfg%initial_temperature = 10.0_wp
   end subroutine base_cfg

   subroutine build(cfg, grid, state)
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DX)
      call state%init_from_config(cfg, grid)
      call ocean_state_seed_from_cfg(state, grid, cfg)
   end subroutine build

   pure function column_s(state, iS, i, j, k) result(s)
      !! Seeded salinity CONCENTRATION at one cell, read back from the
      !! thickness-weighted prognostic.
      type(ocean_state_t), intent(in) :: state
      integer, intent(in) :: iS, i, j, k
      real(wp) :: s
      s = state%multilayer%tracers(iS)%hTr(i, j, k)/state%multilayer%h_layer(i, j, k)
   end function column_s

   pure function column_centre_span(state, i, j) result(dz)
      !! Vertical distance between the bed-layer and surface-layer
      !! CENTRES of one column, from the seeded thicknesses.
      type(ocean_state_t), intent(in) :: state
      integer, intent(in) :: i, j
      real(wp) :: dz
      integer :: k
      dz = 0.5_wp*(state%multilayer%h_layer(i, j, 1) &
                   + state%multilayer%h_layer(i, j, NZ))
      do k = 2, NZ - 1
         dz = dz + state%multilayer%h_layer(i, j, k)
      end do
   end function column_centre_span

   pure function s_expected(k) result(s)
      !! The analytical target: linear in layer index from the bed
      !! (k = 1) to the surface (k = NZ).
      integer, intent(in) :: k
      real(wp) :: s
      s = S_BED + (S_SURF - S_BED)*real(k - 1, wp)/real(NZ - 1, wp)
   end function s_expected

   subroutine test_endpoints_linear(error)
      !! Flat basin: S(k=1) = S_init_bottom, S(k=NZ) = S_init_surface,
      !! and every interior layer sits on the straight line between
      !! them — over the FULL local array, ghosts included.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: s_ijk, dz, dS_dz, z_centre
      integer :: i, j, k, iS

      call base_cfg(cfg, "flat")
      cfg%S_init_surface = S_SURF
      cfg%S_init_bottom = S_BED
      call build(cfg, grid, state)
      iS = state%multilayer%idx_salinity

      checks: block
         call check(error, iS > 0, "salinity must be registered")
         if (allocated(error)) exit checks

         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  s_ijk = state%multilayer%tracers(iS)%hTr(i, j, k) &
                          /state%multilayer%h_layer(i, j, k)
                  call check(error, abs(s_ijk - s_expected(k)) < TOL, &
                             "S(k) must be linear in layer index, ghosts included")
                  if (allocated(error)) exit checks
               end do
            end do
         end do

         ! Endpoint orientation, stated separately so a k-flip cannot
         ! hide inside the linearity loop.
         s_ijk = state%multilayer%tracers(iS)%hTr(NGHOST + 1, NGHOST + 1, 1) &
                 /state%multilayer%h_layer(NGHOST + 1, NGHOST + 1, 1)
         call check(error, abs(s_ijk - S_BED) < TOL, &
                    "k = 1 is the BED and must carry S_init_bottom")
         if (allocated(error)) exit checks
         s_ijk = state%multilayer%tracers(iS)%hTr(NGHOST + 1, NGHOST + 1, NZ) &
                 /state%multilayer%h_layer(NGHOST + 1, NGHOST + 1, NZ)
         call check(error, abs(s_ijk - S_SURF) < TOL, &
                    "k = NZ is the SURFACE and must carry S_init_surface")
         if (allocated(error)) exit checks

         ! Linear in LAYER-CENTRE DEPTH, stated in physical units: with
         ! the equal-thickness sigma seed the layer centres sit at
         ! z(k) = (k - 1/2)*dz above the bed, so S must be the affine
         ! function of z with slope (S_surf - S_bed)/((nz-1)*dz).
         dz = MAX_DEPTH/real(NZ, wp)
         dS_dz = (S_SURF - S_BED)/(real(NZ - 1, wp)*dz)
         do k = 1, NZ
            z_centre = (real(k, wp) - 0.5_wp)*dz
            s_ijk = state%multilayer%tracers(iS)%hTr(NGHOST + 1, NGHOST + 1, k) &
                    /state%multilayer%h_layer(NGHOST + 1, NGHOST + 1, k)
            call check(error, abs(s_ijk - (S_BED + dS_dz*(z_centre - 0.5_wp*dz))) &
                       < 1.0e-10_wp, &
                       "S must be affine in layer-centre depth with the expected slope")
            if (allocated(error)) exit checks
         end do
      end block checks

      call state%destroy()
   end subroutine test_endpoints_linear

   subroutine test_sloping_bathymetry(error)
      !! Spoon basin: the column depth varies cell to cell.  Endpoints
      !! and layer-index linearity must hold on EVERY column, and the
      !! physical dS/dz must scale as 1/depth — i.e. the profile is
      !! sigma-following, exactly like the temperature sibling.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: s_ijk, depth_a, depth_b, grad_a, grad_b, dz_a, dz_b
      integer :: i, j, k, iS, i_shal, i_deep

      call base_cfg(cfg, "spoon")
      cfg%S_init_surface = S_SURF
      cfg%S_init_bottom = S_BED
      call build(cfg, grid, state)
      iS = state%multilayer%idx_salinity

      checks: block
         ! Same profile on every wet column regardless of its depth.
         do k = 1, NZ
            do j = NGHOST + 1, NGHOST + NY_PHYS
               do i = NGHOST + 1, NGHOST + NX_PHYS
                  if (state%multilayer%h_layer(i, j, k) <= 0.0_wp) cycle
                  s_ijk = state%multilayer%tracers(iS)%hTr(i, j, k) &
                          /state%multilayer%h_layer(i, j, k)
                  call check(error, abs(s_ijk - s_expected(k)) < TOL, &
                             "sloping bed: S(k) must still be linear in layer index")
                  if (allocated(error)) exit checks
               end do
            end do
         end do

         ! Find two columns of genuinely different depth on the same row.
         i_shal = NGHOST + 1
         i_deep = NGHOST + NX_PHYS/2
         j = NGHOST + NY_PHYS/2
         depth_a = state%barotropic%b(i_shal, j)
         depth_b = state%barotropic%b(i_deep, j)
         call check(error, abs(depth_b - depth_a) > 1.0_wp, &
                    "spoon basin must give two columns of different depth")
         if (allocated(error)) exit checks

         ! Measured dS/dz between the bed and surface layer CENTRES,
         ! read back out of the SEEDED field and the seeded thicknesses
         ! (no analytical shortcut), one value per column.
         dz_a = column_centre_span(state, i_shal, j)
         dz_b = column_centre_span(state, i_deep, j)
         grad_a = (column_s(state, iS, i_shal, j, NZ) &
                   - column_s(state, iS, i_shal, j, 1))/dz_a
         grad_b = (column_s(state, iS, i_deep, j, NZ) &
                   - column_s(state, iS, i_deep, j, 1))/dz_b
         call check(error, grad_a /= grad_b, &
                    "two columns of different depth must have different dS/dz")
         if (allocated(error)) exit checks
         ! Sigma-following: the gradient scales exactly as 1/depth, so
         ! grad*depth is one and the same constant on both columns.
         call check(error, abs(grad_a*depth_a - grad_b*depth_b) < 1.0e-10_wp, &
                    "measured dS/dz must scale as 1/depth (sigma-following)")
         if (allocated(error)) exit checks
         call check(error, abs(grad_a*depth_a - (S_SURF - S_BED)*real(NZ, wp) &
                               /real(NZ - 1, wp)) < 1.0e-10_wp, &
                    "grad*depth must equal the analytical sigma value")
         if (allocated(error)) exit checks

         ! Endpoints on the shallow column too.
         s_ijk = state%multilayer%tracers(iS)%hTr(i_shal, j, 1) &
                 /state%multilayer%h_layer(i_shal, j, 1)
         call check(error, abs(s_ijk - S_BED) < TOL, &
                    "shallow column bed must still carry S_init_bottom")
         if (allocated(error)) exit checks
         s_ijk = state%multilayer%tracers(iS)%hTr(i_shal, j, NZ) &
                 /state%multilayer%h_layer(i_shal, j, NZ)
         call check(error, abs(s_ijk - S_SURF) < TOL, &
                    "shallow column surface must still carry S_init_surface")
      end block checks

      call state%destroy()
   end subroutine test_sloping_bathymetry

   subroutine test_defaults_uniform(error)
      !! Both knobs at their `0.0` default ⇒ the historical uniform
      !! `initial_salinity` column, byte-identical — asserted as the
      !! EXACT `hTr = S*h` product, not within a tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, k, iS

      call base_cfg(cfg, "flat")
      call build(cfg, grid, state)
      iS = state%multilayer%idx_salinity

      checks: block
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  call check(error, state%multilayer%tracers(iS)%hTr(i, j, k) &
                             == S_UNIFORM*state%multilayer%h_layer(i, j, k), &
                             "default S IC must stay the exact uniform hTr = S*h product")
                  if (allocated(error)) exit checks
               end do
            end do
         end do
      end block checks

      call state%destroy()
   end subroutine test_defaults_uniform

   subroutine test_pseudo_salt_seed(error)
      !! The pseudo-salt verification tracer is seeded from salinity
      !! AFTER every write to S, so a stratified S IC must give it the
      !! same stratified field — deviation exactly zero at t = 0.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, k, iS, iPS

      call base_cfg(cfg, "flat")
      cfg%S_init_surface = S_SURF
      cfg%S_init_bottom = S_BED
      cfg%ocean%tracers%enable_pseudo_salt = .true.
      call build(cfg, grid, state)
      iS = state%multilayer%idx_salinity
      iPS = state%multilayer%idx_pseudo_salt

      checks: block
         call check(error, iPS > 0, "pseudo-salt must be registered when enabled")
         if (allocated(error)) exit checks
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  call check(error, state%multilayer%tracers(iPS)%hTr(i, j, k) &
                             == state%multilayer%tracers(iS)%hTr(i, j, k), &
                             "pseudo-salt must be seeded from the FINAL stratified S")
                  if (allocated(error)) exit checks
               end do
            end do
         end do
      end block checks

      call state%destroy()
   end subroutine test_pseudo_salt_seed

   subroutine test_one_sided_fails(error)
      !! Setting only one end is the same silent-no-op class the dead
      !! knob was: the seed gates on BOTH being non-zero, so a lone
      !! `S_init_surface` would leave a uniform column.  Refuse it.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 "default config must validate cleanly (baseline)")
      if (allocated(error)) return

      cfg%S_init_surface = S_SURF
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "S_init_surface without S_init_bottom must fail loud")
      if (allocated(error)) return

      cfg%S_init_bottom = S_BED
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 "both ends set must validate")
   end subroutine test_one_sided_fails

end module test_ocean_salinity_ic
