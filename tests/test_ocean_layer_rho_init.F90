!! Unit tests for `apply_layer_rho_init` — the direct per-layer
!! density-init helper (MOM6 `COORD_CONFIG="gprime"` analogue).
!!
!! `apply_layer_rho_init(ms, layer_rho_init, nz_ml)` writes
!! `ms%rho_layer(:,:,k) = layer_rho_init(k)` straight into the
!! multilayer C-grid state, bypassing the EOS path entirely.
!! Default sentinel is `-1.0` ⇒ no-op (leave `rho_layer` alone).
!! When set, the caller MUST pass exactly `nz_ml` non-sentinel
!! values in `k=1..nz_ml` order (k=1 bed, k=nz_ml surface — the
!! project-wide bottom-up convention).
!!
!! Cases:
!!
!!   * `sentinel_is_no_op` — `layer_rho_init = -1.0` everywhere
!!     leaves `rho_layer` at its alloc-time value (zero from
!!     `multilayer_state_init`).  Without this guarantee,
!!     setting the namelist default to "disabled" wouldn't actually
!!     keep the EOS path bit-identical.
!!
!!   * `writes_per_layer_density` — `layer_rho_init = [1036.0, 1035.0]`
!!     with `nz_ml = 2` writes 1036.0 to every interior + ghost cell
!!     of `rho_layer(:,:,1)` (bed) and 1035.0 to `rho_layer(:,:,2)`
!!     (surface).  Verifies (a) the values land in `rho_layer` exactly,
!!     (b) the k=1=bed / k=nz=surface convention is honoured (heavier
!!     value at bottom, lighter at top — gravitationally stable), and
!!     (c) every cell of `rho_layer` gets written (no ghost-row gap).
module test_ocean_layer_rho_init
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_state, only: apply_layer_rho_init, ocean_linear_layer_density
   implicit none
   private

   public :: collect_ocean_layer_rho_init_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_layer_rho_init_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sentinel_is_no_op", test_sentinel_is_no_op), &
                  new_unittest("writes_per_layer_density", test_writes_per_layer_density), &
                  new_unittest("linear_density_range", test_linear_density_range) &
                  ]
   end subroutine collect_ocean_layer_rho_init_tests

   subroutine make_ms(grid, ms, nz_ml)
      !! Helper: build a fresh `multilayer_state_t` ready for
      !! `apply_layer_rho_init`.  Caller must set `nz_ml` before
      !! `init` per the type's contract.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nz_ml
      ms%nz_ml = nz_ml
      call ms%init(grid)
   end subroutine make_ms

   subroutine test_sentinel_is_no_op(error)
      !! With every entry of `layer_rho_init` at the sentinel `-1.0`,
      !! `apply_layer_rho_init` returns without touching `rho_layer`.
      !! `multilayer_state_init` allocates `rho_layer` with
      !! `source = 0.0_wp`, so the post-call max abs value must still
      !! be zero (machine-epsilon strict — there's no floating-point
      !! arithmetic on this path).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NZ_ML = 2
      real(wp) :: layer_rho_init(NZ_ML)
      checks: block
         call grid%init(8, 8, NGHOST, 5000.0_wp, 5000.0_wp)
         call make_ms(grid, ms, NZ_ML)

         layer_rho_init = -1.0_wp
         call apply_layer_rho_init(ms, layer_rho_init, NZ_ML)

         call check(error, maxval(abs(ms%rho_layer)) < 1.0e-15_wp, &
                    "sentinel -1.0 should leave rho_layer at its zero IC")
      end block checks
      call ms%destroy()
   end subroutine test_sentinel_is_no_op

   subroutine test_writes_per_layer_density(error)
      !! With `layer_rho_init = [1036.0, 1035.0]` and `nz_ml = 2`,
      !! `apply_layer_rho_init` must:
      !!   (a) write 1036.0 to every cell of `rho_layer(:,:,1)`,
      !!   (b) write 1035.0 to every cell of `rho_layer(:,:,2)`,
      !!   (c) respect the bed-to-surface (k=1 bed, k=nz surface)
      !!       convention so the heavier value lands at the bed.
      !! Equality is exact (no FP arithmetic on the write path).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: NZ_ML = 2
      real(wp), parameter :: RHO_BED = 1036.0_wp
      real(wp), parameter :: RHO_SURF = 1035.0_wp
      real(wp) :: layer_rho_init(NZ_ML)
      checks: block
         call grid%init(8, 8, NGHOST, 5000.0_wp, 5000.0_wp)
         call make_ms(grid, ms, NZ_ML)

         layer_rho_init = [RHO_BED, RHO_SURF]   ! k=1 bed, k=2 surface
         call apply_layer_rho_init(ms, layer_rho_init, NZ_ML)

         call check(error, maxval(abs(ms%rho_layer(:, :, 1) - RHO_BED)) < 1.0e-15_wp, &
                    "k=1 (bed) layer should be RHO_BED everywhere, incl ghosts")
         if (allocated(error)) exit checks

         call check(error, maxval(abs(ms%rho_layer(:, :, 2) - RHO_SURF)) < 1.0e-15_wp, &
                    "k=2 (surface) layer should be RHO_SURF everywhere, incl ghosts")
         if (allocated(error)) exit checks

         ! Sanity: bed denser than surface (gravitationally stable).
         call check(error, ms%rho_layer(NGHOST + 1, NGHOST + 1, 1) > &
                    ms%rho_layer(NGHOST + 1, NGHOST + 1, 2), &
                    "k=1 (bed) must be denser than k=nz (surface) per project convention")
      end block checks
      call ms%destroy()
   end subroutine test_writes_per_layer_density

   subroutine test_linear_density_range(error)
      !! `ocean_linear_layer_density(rho_lightest, rho_range, nz_ml)` — the
      !! MOM6 `COORD_CONFIG="linear"` generator — must return `nz_ml`
      !! linearly-spaced densities, bed `k=1` heaviest → surface `k=nz`
      !! lightest, layer-centred:
      !!     ρ(k) = rho_lightest + rho_range·(nz_ml − k + 0.5)/nz_ml
      !! For rho_lightest = 1035, rho_range = 2.0, nz_ml = 10:
      !!   (a) surface ρ(10) = 1035 + 2·0.05  = 1035.1,
      !!   (b) bed     ρ(1)  = 1035 + 2·0.95  = 1036.9,
      !!   (c) uniform spacing ρ(k) − ρ(k+1) = rho_range/nz_ml = 0.2,
      !!   (d) monotone (bed denser than every shallower layer).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ_ML = 10
      real(wp), parameter :: RHO_LIGHTEST = 1035.0_wp
      real(wp), parameter :: RHO_RANGE = 2.0_wp
      real(wp), parameter :: DRHO = RHO_RANGE/real(NZ_ML, wp)   ! 0.2
      real(wp), parameter :: TOL = 1.0e-12_wp
      real(wp) :: rho(NZ_ML)
      integer :: k
      checks: block
         rho = ocean_linear_layer_density(RHO_LIGHTEST, RHO_RANGE, NZ_ML)

         ! (a) surface (k=nz) is the lightest, half a layer above the floor.
         call check(error, abs(rho(NZ_ML) - (RHO_LIGHTEST + 0.5_wp*DRHO)) < TOL, &
                    "surface layer density should be rho_lightest + 0.5*drho")
         if (allocated(error)) exit checks

         ! (b) bed (k=1) is the heaviest, half a layer below the ceiling.
         call check(error, abs(rho(1) - (RHO_LIGHTEST + RHO_RANGE - 0.5_wp*DRHO)) < TOL, &
                    "bed layer density should be rho_lightest + rho_range - 0.5*drho")
         if (allocated(error)) exit checks

         ! (c) uniform spacing, decreasing with k (bed→surface).
         do k = 1, NZ_ML - 1
            call check(error, abs((rho(k) - rho(k + 1)) - DRHO) < TOL, &
                       "adjacent layers must be uniformly spaced by rho_range/nz_ml")
            if (allocated(error)) exit checks
         end do

         ! (d) gravitationally stable: bed strictly denser than surface.
         call check(error, rho(1) > rho(NZ_ML) + (RHO_RANGE - DRHO) - TOL, &
                    "bed must be denser than surface by ~rho_range")
      end block checks
   end subroutine test_linear_density_range

end module test_ocean_layer_rho_init
