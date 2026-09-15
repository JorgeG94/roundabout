!! Unit tests for the ocean OBC type layer (PR 3 foundation).
module test_ocean_boundary
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, &
                                       ocean_bc_state_init, ocean_bc_state_destroy, &
                                       ocean_bc_type_from_string, &
                                       OBC_WALL, OBC_OPEN, OBC_TIDAL, OBC_NESTED, &
                                       OBC_INFLOW, OBC_DISCHARGE, OBC_CLAMPED, &
                                       OBC_SPONGE, OBC_CHAPMAN, OBC_INVALID
   use rdb_ocean_boundary_data, only: ocean_boundary_data_constant_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_boundary_tests

   integer, parameter :: NX = 8, NY = 6, NZ = 2, NGHOST = 2
   real(wp), parameter :: DX = 1.0_wp

contains

   subroutine collect_ocean_boundary_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bc_state_init_defaults_to_wall", test_bc_state_defaults), &
                  new_unittest("bc_state_destroy_clears_init", test_bc_state_destroy), &
                  new_unittest("bc_type_from_string_known_names", test_type_string_known), &
                  new_unittest("bc_type_from_string_unknown_falls_back", test_type_string_fallback), &
                  new_unittest("constant_data_source_writes_clamped_fields", &
                               test_constant_data_source) &
                  ]
   end subroutine collect_ocean_boundary_tests

   subroutine test_bc_state_defaults(error)
      !! Fresh `ocean_bc_state_t` has all four edges = OBC_WALL — the
      !! Phase 3 closed-wall convention.  Any silent regression here
      !! would mean kernel dispatch reaches a non-wall branch by
      !! default, which is exactly the kind of bug PR 3 must not
      !! introduce.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         call check(error, bc%is_init, "bc_state should be flagged is_init after init")
         if (allocated(error)) exit checks
         call check(error, bc%west%bc_type == OBC_WALL, "west default must be OBC_WALL")
         if (allocated(error)) exit checks
         call check(error, bc%east%bc_type == OBC_WALL, "east default must be OBC_WALL")
         if (allocated(error)) exit checks
         call check(error, bc%south%bc_type == OBC_WALL, "south default must be OBC_WALL")
         if (allocated(error)) exit checks
         call check(error, bc%north%bc_type == OBC_WALL, "north default must be OBC_WALL")
         if (allocated(error)) exit checks
         call check(error, bc%nx_total == grid%nx_total .and. bc%ny_total == grid%ny_total, &
                    "bc_state should cache the grid extents")
         if (allocated(error)) exit checks
         call check(error, bc%nz_ml == NZ .and. bc%n_tracers == 2, &
                    "bc_state should cache nz_ml and n_tracers")
      end block checks
      call ocean_bc_state_destroy(bc)
   end subroutine test_bc_state_defaults

   subroutine test_bc_state_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      call grid%init(NX, NY, NGHOST, DX, DX)
      call ocean_bc_state_init(bc, grid, nz_ml=NZ)
      call ocean_bc_state_destroy(bc)
      call check(error,.not. bc%is_init, &
                 "bc_state%is_init must clear after destroy")
   end subroutine test_bc_state_destroy

   subroutine test_type_string_known(error)
      !! Each known string parses to its corresponding integer tag.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, ocean_bc_type_from_string("wall") == OBC_WALL, "'wall' -> OBC_WALL")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("open") == OBC_OPEN, "'open' -> OBC_OPEN")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("tidal") == OBC_TIDAL, "'tidal' -> OBC_TIDAL")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("nested") == OBC_NESTED, "'nested' -> OBC_NESTED")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("inflow") == OBC_INFLOW, "'inflow' -> OBC_INFLOW")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("discharge") == OBC_DISCHARGE, "'discharge' -> OBC_DISCHARGE")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("clamped") == OBC_CLAMPED, "'clamped' -> OBC_CLAMPED")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("sponge") == OBC_SPONGE, "'sponge' -> OBC_SPONGE")
         if (allocated(error)) exit checks
         call check(error, ocean_bc_type_from_string("chapman") == OBC_CHAPMAN, "'chapman' -> OBC_CHAPMAN")
      end block checks
   end subroutine test_type_string_known

   subroutine test_type_string_fallback(error)
      !! PR-6 fail-loud: an unknown edge string returns OBC_INVALID (NOT
      !! the old OBC_WALL default).  A misspelled edge must not silently
      !! close the boundary to a wall — the user asked for a specific
      !! (often open) boundary, and a wall reflects every outgoing gravity
      !! wave.  `validate_config` rejects OBC_INVALID naming the edge.
      type(error_type), allocatable, intent(out) :: error
      call check(error, ocean_bc_type_from_string("not-a-real-bc") == OBC_INVALID, &
                 "unknown string must return OBC_INVALID (fail-loud, PR-6)")
   end subroutine test_type_string_fallback

   subroutine test_constant_data_source(error)
      !! The constant data source's `update(t, bc)` writes the
      !! configured scalars into `bc%<edge>%clamped_*` fields.  Proves
      !! the polymorphic dispatch is wired and the simplest backend
      !! works end-to-end.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      type(ocean_boundary_data_constant_t) :: src
      checks: block
         call grid%init(NX, NY, NGHOST, DX, DX)
         call ocean_bc_state_init(bc, grid, nz_ml=NZ)
         src%u_west = 0.1_wp
         src%v_north = -0.2_wp
         src%eta_east = 0.5_wp
         call src%update(t=100.0_wp, bc=bc)
         call check(error, abs(bc%west%clamped_u - 0.1_wp) < 1.0e-12_wp, &
                    "constant source should set bc%west%clamped_u")
         if (allocated(error)) exit checks
         call check(error, abs(bc%north%clamped_v - (-0.2_wp)) < 1.0e-12_wp, &
                    "constant source should set bc%north%clamped_v")
         if (allocated(error)) exit checks
         call check(error, abs(bc%east%clamped_eta - 0.5_wp) < 1.0e-12_wp, &
                    "constant source should set bc%east%clamped_eta")
         if (allocated(error)) exit checks
         call check(error, abs(bc%south%clamped_eta) < 1.0e-12_wp, &
                    "untouched constant source fields should default to zero")
      end block checks
      call src%destroy()
      call ocean_bc_state_destroy(bc)
   end subroutine test_constant_data_source

end module test_ocean_boundary
