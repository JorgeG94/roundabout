!! P1 diagnostic-remap tests — conservative z-remap + density-space remap.
!!
!! Exercises the upgraded diagnostic vertical remap (`rdb_ocean_diag_fills`):
!!   * conservative z-remap (intensive): uniform-in-z field is preserved;
!!     a linear-in-z field's column integral is preserved.
!!   * conservative z-remap (extensive): the column SUM is preserved to
!!     round-off.
!!   * density-space remap (intensive): a stratified column maps the source
!!     layer at each density into the matching bin; lightest → surface.
!!   * density-space remap (extensive): the column integral is conserved.
!!   * on-device (enter_data round-trip): the remap fires through the GPU
!!     pipeline and stays finite + conservative.
!!
!! Density predictability: the in-test EOS is forced linear with
!! `rho = rho0 - alpha_T·T` (beta_S = 0, T_ref = S_ref = 0), so a chosen
!! per-layer temperature gives a chosen potential density.  Warmer = lighter
!! ⇒ the surface (k=NZ) layer is the lightest, matching the lightest→surface
!! density-bin convention.
module test_ocean_diag_remap
   use rdb_constants, only: wp, REMAP_PCM, REMAP_PPM
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: DIAG_OP_INSTANT, DIAG_VGRID_Z_FIXED, DIAG_VGRID_DENSITY, &
                             DIAG_VGRID_SIGMA, DIAG_VGRID_ZSTAR, DIAG_MISSING_VALUE, &
                             DIAG_VGRID_LAYER, diag_remap_proc, diag_spec_t, parse_diag_spec
   use rdb_ocean_diag_fills, only: fill_temperature, &
                                   remap_layer_to_z, &
                                   remap_layer_to_sigma, remap_layer_to_zstar, &
                                   remap_layer_to_density, &
                                   set_diag_remap_method, set_diag_mask_vanished, &
                                   coord_remap_proc
   use rdb_driver, only: diag_vgrid_from_name
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_diag_remap_tests

   integer, parameter :: NX = 5, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1.0_wp
   real(wp), parameter :: TOL = 1.0e-9_wp

contains

   subroutine collect_ocean_diag_remap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("diag_remap_z_consistency", test_z_consistency), &
                  new_unittest("diag_remap_z_integral_preserved", test_z_integral), &
                  new_unittest("diag_remap_z_extensive_conserves", test_z_extensive), &
                  new_unittest("diag_remap_density_bins", test_density_bins), &
                  new_unittest("diag_remap_density_conserves", test_density_extensive), &
                  new_unittest("diag_remap_sigma_consistency", test_sigma_consistency), &
                  new_unittest("diag_remap_sigma_extensive_conserves", test_sigma_extensive), &
                  new_unittest("diag_remap_zstar_eq_sigma_uniform", test_zstar_eq_sigma), &
                  new_unittest("diag_remap_zstar_nonuniform_differs", test_zstar_nonuniform), &
                  new_unittest("diag_remap_mask_off_reads_zero", test_mask_off_zero), &
                  new_unittest("diag_remap_mask_on_reads_missing", test_mask_on_missing), &
                  new_unittest("diag_remap_ppm_higher_order", test_z_ppm_exact), &
                  new_unittest("diag_remap_on_device_finite", test_on_device), &
                  new_unittest("diag_dispatch_honors_is_extensive", test_dispatch_extensive), &
                  new_unittest("diag_density_vgrid_from_namelist", test_density_vgrid_from_namelist) &
                  ]
   end subroutine collect_ocean_diag_remap_tests

   subroutine test_z_ppm_exact(error)
      !! Selector / higher-order gate.  A field linear in depth (layer values
      !! = the midpoint of f(z)=z) remapped onto an OFFSET z-grid is
      !! reproduced essentially EXACTLY by PPM (a linear profile is in its
      !! reconstruction space) but STAIR-STEPPED by PCM (each layer is a
      !! constant = its mean).  Proves `diag_remap_scheme` actually switches
      !! the reconstruction — the diag remap used to be hard-wired to PCM.
      !! Source: H=30, 3 layers of 10 m; f=z midpoints (bottom-up) 25/15/5.
      !! Target cells [0,15],[15,30] ⇒ analytic averages 7.5 / 22.5.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 2)
      real(wp), parameter :: Z_IFACE(2) = [15.0_wp, 30.0_wp]
      real(wp), parameter :: A1 = 7.5_wp, A2 = 22.5_wp
      real(wp) :: err_ppm, err_pcm
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp
         layer_buf(:, :, 1) = 25.0_wp   ! bed   (z in [20,30], midpoint 25)
         layer_buf(:, :, 2) = 15.0_wp   !       (z in [10,20], midpoint 15)
         layer_buf(:, :, 3) = 5.0_wp    ! surf  (z in [0,10],  midpoint 5)

         call set_diag_remap_method(REMAP_PPM)
         call remap_layer_to_z(state, Z_IFACE, layer_buf, out_buf, .false.)
         err_ppm = max(abs(out_buf(2, 2, 1) - A1), abs(out_buf(2, 2, 2) - A2))

         call set_diag_remap_method(REMAP_PCM)
         call remap_layer_to_z(state, Z_IFACE, layer_buf, out_buf, .false.)
         err_pcm = max(abs(out_buf(2, 2, 1) - A1), abs(out_buf(2, 2, 2) - A2))

         call check(error, err_pcm > 0.5_wp, &
                    "PCM must stair-step a linear field on an offset grid (large error)")
         if (allocated(error)) exit checks
         call check(error, err_ppm < 1.0e-6_wp, &
                    "PPM must reproduce a linear-in-z field ~exactly (selector is live)")
         if (allocated(error)) exit checks
         call check(error, err_ppm < 0.01_wp*err_pcm, &
                    "PPM must be far more accurate than PCM here")
      end block checks
      call set_diag_remap_method(REMAP_PPM)   ! restore the production default
      call state%destroy()
   end subroutine test_z_ppm_exact

   subroutine setup_state(grid, state)
      !! NZ-layer cartesian ocean state with the linear EOS forced to the
      !! clean `rho = rho0 - alpha_T·T` form so layer densities are
      !! analytically predictable.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      ! Clean linear EOS: rho = rho0 - alpha_T*T (S term off).
      state%eos%rho0 = 1025.0_wp
      state%eos%alpha_T = 1.0_wp
      state%eos%beta_S = 0.0_wp
      state%eos%T_ref = 0.0_wp
      state%eos%S_ref = 0.0_wp
   end subroutine setup_state

   subroutine set_uniform_TS(state, temp_by_layer)
      !! Stamp h_layer = 10 and T (via hTr) so layer k has temperature
      !! `temp_by_layer(k)`; S held at S_ref (= 0) so density depends on T
      !! only.  rho(k) = 1025 - temp_by_layer(k).
      type(ocean_state_t), intent(inout) :: state
      real(wp), intent(in) :: temp_by_layer(:)
      integer :: it, is_, k
      it = state%multilayer%idx_temperature
      is_ = state%multilayer%idx_salinity
      state%multilayer%h_layer = 10.0_wp
      do k = 1, NZ
         state%multilayer%tracers(it)%hTr(:, :, k) = 10.0_wp*temp_by_layer(k)
         state%multilayer%tracers(is_)%hTr(:, :, k) = 0.0_wp
      end do
   end subroutine set_uniform_TS

   ! ------------------------------------------------------------------
   ! Conservative z-remap
   ! ------------------------------------------------------------------

   subroutine test_z_consistency(error)
      !! A uniform-in-z intensive field remapped to z returns the same
      !! value in every (non-empty) target cell.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 3)
      real(wp), parameter :: Z_IFACE(3) = [10.0_wp, 20.0_wp, 30.0_wp]
      real(wp), parameter :: VAL = 7.25_wp
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         layer_buf = VAL
         call remap_layer_to_z(state, Z_IFACE, layer_buf, out_buf, .false.)
         call check(error, abs(out_buf(2, 2, 1) - VAL) < TOL, "cell 1 should equal VAL")
         if (allocated(error)) exit checks
         call check(error, abs(out_buf(2, 2, 2) - VAL) < TOL, "cell 2 should equal VAL")
         if (allocated(error)) exit checks
         call check(error, abs(out_buf(2, 2, 3) - VAL) < TOL, "cell 3 should equal VAL")
      end block checks
      call state%destroy()
   end subroutine test_z_consistency

   subroutine test_z_integral(error)
      !! A linear-in-z intensive field: the thickness-weighted column
      !! integral of the remapped output equals the source integral when
      !! the target grid spans the whole column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 3)
      real(wp), parameter :: Z_IFACE(3) = [7.0_wp, 19.0_wp, 30.0_wp]
      real(wp) :: src_int, out_int, dz_new
      integer :: m
      real(wp) :: zf_prev, zf
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         ! intensive layer values (bottom-up): 1, 4, 9
         layer_buf(:, :, 1) = 1.0_wp
         layer_buf(:, :, 2) = 4.0_wp
         layer_buf(:, :, 3) = 9.0_wp
         call remap_layer_to_z(state, Z_IFACE, layer_buf, out_buf, .false.)
         ! source integral = Σ value·dz
         src_int = (1.0_wp + 4.0_wp + 9.0_wp)*10.0_wp
         ! output integral = Σ value·dz_new (cells clipped to H=30)
         out_int = 0.0_wp
         zf_prev = 0.0_wp
         do m = 1, 3
            zf = min(max(Z_IFACE(m), 0.0_wp), 30.0_wp)
            dz_new = zf - zf_prev
            out_int = out_int + out_buf(2, 2, m)*dz_new
            zf_prev = zf
         end do
         call check(error, abs(out_int - src_int) < 1.0e-8_wp, &
                    "remapped column integral should equal the source integral")
      end block checks
      call state%destroy()
   end subroutine test_z_integral

   subroutine test_z_extensive(error)
      !! An extensive field (thickness-integrated) remapped to z preserves
      !! the column SUM to round-off (the spanning target grid).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 4)
      real(wp), parameter :: Z_IFACE(4) = [5.0_wp, 13.0_wp, 22.0_wp, 30.0_wp]
      real(wp) :: src_sum, out_sum
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         ! extensive per-layer values (already integrated): 30, 50, 20
         layer_buf(:, :, 1) = 30.0_wp
         layer_buf(:, :, 2) = 50.0_wp
         layer_buf(:, :, 3) = 20.0_wp
         call remap_layer_to_z(state, Z_IFACE, layer_buf, out_buf, .true.)
         src_sum = 30.0_wp + 50.0_wp + 20.0_wp
         out_sum = sum(out_buf(2, 2, :))
         call check(error, abs(out_sum - src_sum) < 1.0e-8_wp, &
                    "extensive remap should preserve the column sum")
      end block checks
      call state%destroy()
   end subroutine test_z_extensive

   ! ------------------------------------------------------------------
   ! Conservative sigma / z* remap
   ! ------------------------------------------------------------------

   subroutine test_sigma_consistency(error)
      !! A uniform-in-z intensive field remapped to sigma returns the same
      !! value in every target cell (terrain-following grid spans the column).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 3)
      real(wp), parameter :: SIG(3) = [1.0_wp/3.0_wp, 2.0_wp/3.0_wp, 1.0_wp]
      real(wp), parameter :: VAL = 4.5_wp
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         layer_buf = VAL
         call remap_layer_to_sigma(state, SIG, layer_buf, out_buf, .false.)
         call check(error, abs(out_buf(2, 2, 1) - VAL) < TOL .and. &
                    abs(out_buf(2, 2, 2) - VAL) < TOL .and. &
                    abs(out_buf(2, 2, 3) - VAL) < TOL, &
                    "uniform field on sigma should be VAL in every cell")
      end block checks
      call state%destroy()
   end subroutine test_sigma_consistency

   subroutine test_sigma_extensive(error)
      !! An extensive field remapped to sigma preserves the column sum
      !! (sigma always spans the full column => no below-grid loss).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 4)
      real(wp), parameter :: SIG(4) = [0.2_wp, 0.5_wp, 0.8_wp, 1.0_wp]
      real(wp) :: out_sum
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         layer_buf(:, :, 1) = 30.0_wp
         layer_buf(:, :, 2) = 50.0_wp
         layer_buf(:, :, 3) = 20.0_wp
         call remap_layer_to_sigma(state, SIG, layer_buf, out_buf, .true.)
         out_sum = sum(out_buf(2, 2, :))
         call check(error, abs(out_sum - 100.0_wp) < 1.0e-8_wp, &
                    "extensive sigma remap should preserve the column sum")
      end block checks
      call state%destroy()
   end subroutine test_sigma_extensive

   subroutine test_zstar_eq_sigma(error)
      !! Under UNIFORM levels z* == sigma: zstar levels [10,20,30] (H_ref=30)
      !! on a 30 m column give the same interfaces as sigma [1/3,2/3,1].
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_s(NX, NY, 3), out_z(NX, NY, 3)
      real(wp), parameter :: SIG(3) = [1.0_wp/3.0_wp, 2.0_wp/3.0_wp, 1.0_wp]
      real(wp), parameter :: ZS(3) = [10.0_wp, 20.0_wp, 30.0_wp]
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30 = H_ref
         layer_buf(:, :, 1) = 1.0_wp
         layer_buf(:, :, 2) = 5.0_wp
         layer_buf(:, :, 3) = 9.0_wp
         call remap_layer_to_sigma(state, SIG, layer_buf, out_s, .false.)
         call remap_layer_to_zstar(state, ZS, layer_buf, out_z, .false.)
         call check(error, maxval(abs(out_s(2, 2, :) - out_z(2, 2, :))) < TOL, &
                    "uniform z* must equal sigma")
      end block checks
      call state%destroy()
   end subroutine test_zstar_eq_sigma

   subroutine test_zstar_nonuniform(error)
      !! A NON-uniform z* reference (fine near surface: [5,15,30]) differs
      !! from uniform sigma [1/3,2/3,1] for a depth-varying field.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_s(NX, NY, 3), out_z(NX, NY, 3)
      real(wp), parameter :: SIG(3) = [1.0_wp/3.0_wp, 2.0_wp/3.0_wp, 1.0_wp]
      real(wp), parameter :: ZS(3) = [5.0_wp, 15.0_wp, 30.0_wp]
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30 = H_ref
         layer_buf(:, :, 1) = 1.0_wp   ! bed
         layer_buf(:, :, 2) = 5.0_wp
         layer_buf(:, :, 3) = 9.0_wp   ! surface
         call remap_layer_to_sigma(state, SIG, layer_buf, out_s, .false.)
         call remap_layer_to_zstar(state, ZS, layer_buf, out_z, .false.)
         call check(error, maxval(abs(out_s(2, 2, :) - out_z(2, 2, :))) > TOL, &
                    "non-uniform z* must differ from sigma")
      end block checks
      call state%destroy()
   end subroutine test_zstar_nonuniform

   ! ------------------------------------------------------------------
   ! Vanished-target masking (B3)
   ! ------------------------------------------------------------------

   subroutine test_mask_off_zero(error)
      !! Masking OFF (default): a below-bottom z cell reads 0, not the
      !! sentinel — bit-identical to the legacy writer.  Column H=30, target
      !! depths [10,20,30,40] => cell 4 [30,40] is entirely below the bed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 4)
      real(wp), parameter :: ZI(4) = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         layer_buf = 7.0_wp
         call set_diag_mask_vanished(.false.)
         call remap_layer_to_z(state, ZI, layer_buf, out_buf, .false.)
         call check(error, abs(out_buf(2, 2, 1) - 7.0_wp) < TOL, "valid cell = 7")
         if (allocated(error)) exit checks
         call check(error, abs(out_buf(2, 2, 4)) < TOL, &
                    "below-bottom cell reads 0 when masking off")
      end block checks
      call state%destroy()
   end subroutine test_mask_off_zero

   subroutine test_mask_on_missing(error)
      !! Masking ON: the below-bottom cell carries DIAG_MISSING_VALUE while
      !! valid overlapping cells keep their data.  Resets the module flag so
      !! it does not leak into other tests in this binary.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 4)
      real(wp), parameter :: ZI(4) = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp   ! H = 30
         layer_buf = 7.0_wp
         call set_diag_mask_vanished(.true.)
         call remap_layer_to_z(state, ZI, layer_buf, out_buf, .false.)
         call check(error, abs(out_buf(2, 2, 1) - 7.0_wp) < TOL, &
                    "valid overlapping cell keeps its data (not masked)")
         if (allocated(error)) exit checks
         call check(error, abs(out_buf(2, 2, 4) - DIAG_MISSING_VALUE) < 1.0_wp, &
                    "below-bottom cell masked to missing sentinel")
      end block checks
      call set_diag_mask_vanished(.false.)   ! reset module state for other tests
      call state%destroy()
   end subroutine test_mask_on_missing

   ! ------------------------------------------------------------------
   ! Density-space remap
   ! ------------------------------------------------------------------

   subroutine test_density_bins(error)
      !! 3-layer column, densities (top-down) 1023 / 1024 / 1025 from the
      !! forced linear EOS (T = 2 / 1 / 0 at k=NZ / 2 / 1).  Field = layer
      !! temperature.  Target densities [1023.5, 1024.5] (n_bin = 2).
      !! The lightest source water (T=2, surface) must land in the first
      !! (lightest, surface) output cell; the densest (T=0, bed) in the
      !! last.  Checks lightest→surface ordering + that the cell values
      !! lie within the source range.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 2)
      real(wp), parameter :: RHO_TGT(2) = [1023.5_wp, 1024.5_wp]
      integer :: it, k, ng
      checks: block
         call setup_state(grid, state)
         ! T per layer (bottom-up): k=1 bed T=0, k=2 T=1, k=3 surface T=2.
         call set_uniform_TS(state, [0.0_wp, 1.0_wp, 2.0_wp])
         ! field = layer temperature (intensive), via fill_temperature path
         it = state%multilayer%idx_temperature
         ng = grid%nghost
         ! hTr/h_layer are ghost-inclusive (nx_total); layer_buf is interior
         ! (NX,NY) — slice to the interior so the shapes conform (the old whole
         ! -array assign was non-conformant: silently OOB without runtime checks).
         do k = 1, NZ
            layer_buf(:, :, k) = &
               state%multilayer%tracers(it)%hTr(ng + 1:ng + NX, ng + 1:ng + NY, k)/ &
               state%multilayer%h_layer(ng + 1:ng + NX, ng + 1:ng + NY, k)
         end do
         call remap_layer_to_density(state, RHO_TGT, layer_buf, out_buf, .false.)
         ! lightest bin (surface, cell 1) should be near the warmest (T=2)
         ! source; densest bin (cell 2) near the coldest (T=0).  Lightest →
         ! surface ⇒ out(.,.,1) > out(.,.,2).
         call check(error, out_buf(2, 2, 1) > out_buf(2, 2, 2) + 0.5_wp, &
                    "lightest density bin (surface) should hold the warmer water")
         if (allocated(error)) exit checks
         ! Strong value check (not just in-range): the lightest bin
         ! [0, ρ=1023.5 interface] coincides with the surface layer (T=2,
         ! ρ=1023), so the surface bin value must equal 2.0 to round-off.
         ! A mis-bracket that still landed in-range would fail this.
         call check(error, abs(out_buf(2, 2, 1) - 2.0_wp) < 1.0e-9_wp, &
                    "surface density bin should equal the surface-layer value (2.0)")
         if (allocated(error)) exit checks
         call check(error, out_buf(2, 2, 2) <= 1.0_wp + TOL .and. &
                    out_buf(2, 2, 2) >= 0.0_wp - TOL, &
                    "deep bin value should lie in the cold source range")
      end block checks
      call state%destroy()
   end subroutine test_density_bins

   subroutine test_density_extensive(error)
      !! Extensive field remapped to density bins preserves the column
      !! sum to round-off (the last bin absorbs down to the bed, so the
      !! density grid spans the whole column).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), out_buf(NX, NY, 2)
      real(wp), parameter :: RHO_TGT(2) = [1023.5_wp, 1024.5_wp]
      real(wp) :: src_sum, out_sum
      checks: block
         call setup_state(grid, state)
         call set_uniform_TS(state, [0.0_wp, 1.0_wp, 2.0_wp])
         ! extensive per-layer values (bottom-up): 12, 34, 56
         layer_buf(:, :, 1) = 12.0_wp
         layer_buf(:, :, 2) = 34.0_wp
         layer_buf(:, :, 3) = 56.0_wp
         call remap_layer_to_density(state, RHO_TGT, layer_buf, out_buf, .true.)
         src_sum = 12.0_wp + 34.0_wp + 56.0_wp
         out_sum = sum(out_buf(2, 2, :))
         call check(error, abs(out_sum - src_sum) < 1.0e-7_wp, &
                    "extensive density remap should preserve the column sum")
      end block checks
      call state%destroy()
   end subroutine test_density_extensive

   ! ------------------------------------------------------------------
   ! On-device (GPU) round-trip
   ! ------------------------------------------------------------------

   subroutine test_on_device(error)
      !! Register a Z_FIXED and a DENSITY temperature diagnostic, fire both
      !! through `step` under enter_data (device pipeline), confirm finite
      !! output + that the conservative z-remap preserved the column
      !! integral end-to-end.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: Z_IFACE(3) = [10.0_wp, 20.0_wp, 30.0_wp]
      real(wp), parameter :: RHO_TGT(2) = [1023.5_wp, 1024.5_wp]
      integer :: iz, ird, iv
      real(wp) :: out_int
      checks: block
         call setup_state(grid, state)
         call set_uniform_TS(state, [0.0_wp, 1.0_wp, 2.0_wp])

         call state%diag%set_output_z_levels(Z_IFACE)
         call state%diag%set_output_density_levels(RHO_TGT)
         call state%diag%register("T_z", units="degC", fill=fill_temperature, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z, &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call state%diag%register("T_rho", units="degC", fill=fill_temperature, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_DENSITY, &
                                  remap=remap_layer_to_density, &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)

         iz = 0
         ird = 0
         do iv = 1, state%diag%nvars
            if (state%diag%vars(iv)%name == "T_z") iz = iv
            if (state%diag%vars(iv)%name == "T_rho") ird = iv
         end do
         call check(error, iz > 0 .and. ird > 0, "both diags must register")
         if (allocated(error)) exit checks

         call check(error, size(state%diag%vars(iz)%output_buffer, 3) == 3, &
                    "z output_buffer should have nz_out=3")
         if (allocated(error)) exit checks
         call check(error, size(state%diag%vars(ird)%output_buffer, 3) == 2, &
                    "density output_buffer should have n_rho_out=2")
         if (allocated(error)) exit checks

         call check(error, all(ieee_is_finite_buf(state%diag%vars(iz)%output_buffer)), &
                    "z remap output must be finite")
         if (allocated(error)) exit checks
         call check(error, all(ieee_is_finite_buf(state%diag%vars(ird)%output_buffer)), &
                    "density remap output must be finite")
         if (allocated(error)) exit checks

         ! z grid spans the column (H=30) at 10/20/30 ⇒ intensive integral
         ! of T = (0+1+2)*10 = 30 must be preserved.
         out_int = sum(state%diag%vars(iz)%output_buffer(2, 2, :)*10.0_wp)
         call check(error, abs(out_int - 30.0_wp) < 1.0e-7_wp, &
                    "device z remap should preserve the T column integral (30)")
      end block checks
      call state%destroy()
   end subroutine test_on_device

   subroutine test_dispatch_extensive(error)
      !! The manager's `fill_and_remap` forwards `diag_var_t%is_extensive` to
      !! the remap proc (B1 wiring): an extensive-registered z diag takes the
      !! conservative-redistribution path, giving a DIFFERENT result from the
      !! intensive-averaged one for the same fill on a non-uniform profile —
      !! proving the flag drives dispatch, not just metadata.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      ! MISALIGNED target (2 cells over 3×10 m source layers): cell [0,15]
      ! straddles layers, so intensive (weighted average) and extensive
      ! (conservative redistribution) give genuinely different results — an
      ! aligned grid would make them coincide and the test vacuous.
      real(wp), parameter :: Z_IFACE(2) = [15.0_wp, 30.0_wp]
      integer :: i_int, i_ext, iv
      real(wp) :: dmax
      checks: block
         call setup_state(grid, state)
         call set_uniform_TS(state, [0.0_wp, 1.0_wp, 2.0_wp])
         call state%diag%set_output_z_levels(Z_IFACE)
         call state%diag%register("T_int", units="degC", fill=fill_temperature, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z, is_extensive=.false., &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call state%diag%register("T_ext", units="degC", fill=fill_temperature, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z, is_extensive=.true., &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)
         i_int = 0
         i_ext = 0
         do iv = 1, state%diag%nvars
            if (state%diag%vars(iv)%name == "T_int") i_int = iv
            if (state%diag%vars(iv)%name == "T_ext") i_ext = iv
         end do
         call check(error, i_int > 0 .and. i_ext > 0, "both diags must register")
         if (allocated(error)) exit checks
         call check(error, all(ieee_is_finite_buf(state%diag%vars(i_ext)%output_buffer)), &
                    "extensive remap output must be finite")
         if (allocated(error)) exit checks
         ! Non-uniform T (0,1,2): intensive average /= extensive redistribute,
         ! so the flag MUST change the result if dispatch honours it.
         dmax = maxval(abs(state%diag%vars(i_int)%output_buffer &
                           - state%diag%vars(i_ext)%output_buffer))
         call check(error, dmax > 1.0e-6_wp, &
                    "is_extensive must route to a different remap (intensive /= extensive)")
      end block checks
      call state%destroy()
   end subroutine test_dispatch_extensive

   subroutine test_density_vgrid_from_namelist(error)
      !! PR-9 §9.4 (namelist-reachability half — the numerics are covered
      !! by `diag_remap_density_bins`/`_conserves` above, which predate
      !! this PR).  Proves the "density" `&ocean_diag_nml vgrid` string
      !! now reaches the SAME remap procedure the state-level API uses:
      !!   1. `diag_vgrid_from_name("density")` (driver) -> DIAG_VGRID_DENSITY
      !!   2. `coord_remap_proc(DIAG_VGRID_DENSITY, ...)` -> remap_layer_to_density
      !!   3. `parse_diag_spec("temperature:density")` -> the same coord tag.
      !! Before this PR neither (1) nor (3) existed: `diag_vgrid_from_name`
      !! had no "density" case (silent fallback to LAYER) and the `:coord`
      !! attribute parser rejected `:density` outright.
      type(error_type), allocatable, intent(out) :: error
      procedure(diag_remap_proc), pointer :: remap
      type(diag_spec_t), allocatable :: s(:)
      integer :: vgrid
      checks: block
         vgrid = diag_vgrid_from_name("density")
         call check(error, vgrid == DIAG_VGRID_DENSITY, &
                    "diag_vgrid_from_name('density') must yield DIAG_VGRID_DENSITY")
         if (allocated(error)) exit checks

         call coord_remap_proc(vgrid, remap)
         call check(error, associated(remap, remap_layer_to_density), &
                    "coord_remap_proc(DIAG_VGRID_DENSITY) must yield remap_layer_to_density")
         if (allocated(error)) exit checks

         ! Sanity: an unrelated name still falls back to LAYER (no remap),
         ! so the new case did not widen the default.
         call check(error, diag_vgrid_from_name("not-a-real-vgrid") == DIAG_VGRID_LAYER, &
                    "unknown vgrid name must still fall back to LAYER")
         if (allocated(error)) exit checks

         s = parse_diag_spec("temperature:density")
         call check(error, s(1)%coord == DIAG_VGRID_DENSITY, &
                    "parse_diag_spec('temperature:density') must set coord=DIAG_VGRID_DENSITY")
      end block checks
   end subroutine test_density_vgrid_from_namelist

   pure elemental function ieee_is_finite_buf(x) result(ok)
      !! Local finite check (avoids importing ieee_arithmetic just for the
      !! buffer scan).  A value is finite iff it equals itself and its
      !! magnitude is below huge.
      real(wp), intent(in) :: x
      logical :: ok
      ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function ieee_is_finite_buf

end module test_ocean_diag_remap
