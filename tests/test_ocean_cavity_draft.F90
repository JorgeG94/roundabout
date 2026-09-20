!! Analytical + contract tests for the static ice-shelf cavity GEOMETRY
!! and the barotropic DATUM that absorbs it (`&ocean_cavity_dyn_nml`).
module test_ocean_cavity_draft
   !! The whole capability is three equations:
   !!
   !! ```
   !! (D)  bt_H_ref  = b - z_draft
   !! (P)  p_ice_ref = (rho_ref*GRAVITY) * z_draft
   !! (I)  rho*g*z_draft + (bt_H_ref - b)*rho*g == 0   <=>   (D)
   !! ```
   !!
   !! and one decision (grounding: too little water under the ice ⇒ the
   !! column is LAND, through the SAME wet-mask seed the bathymetry uses).
   !! This suite pins each of them, plus the two things that make the
   !! geometry reach the kernels correctly at all — ghost filling and the
   !! metres→GRID-units conversion — plus the refusal matrix and the
   !! memory/restart bookkeeping.
   !!
   !! Host-only by construction: everything under test runs at SETUP,
   !! before `ocean_state_enter_data`, so there is no device-resident
   !! array in any loop here and no `mem:separate` mapping to arrange.
   !! (The three new arrays ARE mapped — by `metrics%enter_data`, with
   !! their `arr_bytes` terms asserted below — but nothing reads them on
   !! the device in this slice.)
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, LAND_DEPTH_THRESHOLD
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg, &
                              ocean_state_build_restart_registry
   use rdb_ocean_setup, only: configure_ocean_bt_split, configure_ocean_cavity, &
                              configure_ocean_pgf
   use rdb_ocean_cavity, only: set_draft_flat, set_draft_linear, &
                               parse_cavity_draft_config, parse_cavity_draft_source, &
                               CAVITY_DRAFT_NONE, CAVITY_DRAFT_FLAT, &
                               CAVITY_DRAFT_LINEAR, CAVITY_DRAFT_FILE, &
                               CAVITY_DRAFT_INVALID, CAVITY_SOURCE_DRAFT, &
                               CAVITY_SOURCE_THICKNESS, CAVITY_SOURCE_IN_SITU, &
                               CAVITY_SOURCE_INVALID, CAVITY_BOUND_INF, &
                               cavity_datum_residual
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

   public :: collect_ocean_cavity_draft_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 6
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: DY = 1000.0_wp
   real(wp), parameter :: BED = 1000.0_wp
   real(wp), parameter :: DRAFT = 300.0_wp

   ! Metres per degree of latitude at the shipped `rad_earth`, i.e. the
   ! factor `topo_length_to_grid_units` divides by on a spherical grid.
   real(wp), parameter :: RAD_EARTH = 6.378e6_wp
   real(wp), parameter :: PI_WP = 4.0_wp*atan(1.0_wp)
   real(wp), parameter :: M_PER_DEG = RAD_EARTH*PI_WP/180.0_wp

contains

   subroutine collect_ocean_cavity_draft_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_parse_enums", test_parse_enums), &
                  new_unittest("cavity_draft_flat_values_ghosts_front", test_flat_setter), &
                  new_unittest("cavity_draft_linear_profile_and_clip", test_linear_setter), &
                  new_unittest("cavity_draft_units_are_grid_units", test_units_are_grid_units), &
                  new_unittest("cavity_datum_and_init_balance", test_datum_and_init_balance), &
                  new_unittest("cavity_sloping_draft_bt_eta_zero", test_sloping_bt_eta_zero), &
                  new_unittest("cavity_grounding_marks_land", test_grounding_marks_land), &
                  new_unittest("cavity_no_ice_over_land", test_no_ice_over_land), &
                  new_unittest("cavity_grounded_fraction_fails_loud", test_grounded_fraction), &
                  new_unittest("cavity_p_ice_ref_is_same_product", test_p_ice_ref), &
                  new_unittest("cavity_off_bit_identical", test_off_bit_identical), &
                  new_unittest("cavity_bytes_accounting", test_bytes_accounting), &
                  new_unittest("cavity_restart_registry_excludes_draft", test_restart_registry), &
                  new_unittest("cavity_validate_config_matrix", test_validate_matrix) &
                  ]
   end subroutine collect_ocean_cavity_draft_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   function nml_case(cavity, topo, vcoord, bt, pgf, extra) result(nml)
      !! A minimal cavity-capable ocean namelist: flat 1000 m bed, sigma,
      !! FV-MOM6 PGF, split solver — i.e. inside the v1 envelope, so any
      !! refusal a test sees comes from the line that test changed.
      !!
      !! Each group that a test might need to REPLACE is an optional
      !! argument rather than something appended to a fixed prefix: a
      !! namelist carrying the same group twice is a trap (the second
      !! read wins silently, or does not, depending on the reader), and a
      !! test that accidentally relied on that would be testing the
      !! reader, not the cavity.
      character(len=*), intent(in) :: cavity
         !! The `&ocean_cavity_dyn_nml` body, WITHOUT the group name or
         !! the closing `/` — empty string for "no cavity group at all".
      character(len=*), intent(in), optional :: topo, vcoord, bt, pgf, extra
      character(len=:), allocatable :: nml
      character(len=:), allocatable :: topo_l, vcoord_l, bt_l, pgf_l, extra_l

      topo_l = "&ocean_topo_nml max_depth = 1000.0 /"
      if (present(topo)) topo_l = topo
      vcoord_l = "&vcoord_nml vcoord_type = 'sigma' /"
      if (present(vcoord)) vcoord_l = vcoord
      bt_l = "&ocean_bt_nml auto_n_inner = .false., n_inner = 8 /"
      if (present(bt)) bt_l = bt
      ! The production cavity spelling (P5.2): the isostatic load has to
      ! have a consumer, and `p_top_in_bc` is the only one.  Configure
      ! refuses a NON-UNIFORM draft without it, so the default carries it
      ! — the `expect_invalid` rows that probe the PGF envelope override
      ! this whole line anyway.
      pgf_l = "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"
      if (present(pgf)) pgf_l = pgf
      extra_l = ""
      if (present(extra)) extra_l = extra

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            topo_l//new_line("a")//pgf_l//new_line("a")//vcoord_l//new_line("a")// &
            bt_l//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
      if (len_trim(cavity) > 0) then
         nml = nml//"&ocean_cavity_dyn_nml "//cavity//" /"//new_line("a")
      end if
      if (len_trim(extra_l) > 0) nml = nml//extra_l//new_line("a")
   end function nml_case

   function cavity_on() result(body)
      !! The stock in-envelope cavity: a 300 m flat lid over the 1000 m bed.
      character(len=:), allocatable :: body
      body = "enable = .true., draft_config = 'flat', draft_depth = 300.0"
   end function cavity_on

   subroutine build_seeded(cfg, grid, state, nml, ierr)
      !! Production configure path down to the seeded state: parse ->
      !! `init_from_config` (which latches `metrics%use_cavity` before the
      !! allocations) -> `ocean_state_seed_from_cfg` (bathymetry, draft,
      !! grounding, layer split).
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state
      character(len=*), intent(in) :: nml
      integer, intent(out) :: ierr

      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      ! Take dx/dy FROM the namelist, never from the module constants:
      ! the spherical case ships degrees there, and a grid built in
      ! metres against a degrees namelist is exactly the units bug the
      ! `cavity_draft_units_are_grid_units` case exists to catch.
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, cfg%dx, cfg%dy)
      call state%init_from_config(cfg, grid)
      call ocean_state_seed_from_cfg(state, grid, cfg, ierr=ierr)
   end subroutine build_seeded

   subroutine build_configured(cfg, grid, state, nml, ierr)
      !! `build_seeded` plus the configure-time steps this slice touches:
      !! the PGF (which settles `rho_ref`), the `bt_H_ref` latch, and the
      !! cavity load + datum assertion.
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state
      character(len=*), intent(in) :: nml
      integer, intent(out) :: ierr

      call build_seeded(cfg, grid, state, nml, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call configure_ocean_pgf(cfg, state, 1, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call configure_ocean_bt_split(cfg, state, grid, 1)
      call configure_ocean_cavity(cfg, state, grid, 1, ierr=ierr)
   end subroutine build_configured

   pure function cell_x(i) result(x)
      !! Global physical x of cell-centre `i` (single rank, offset 0) —
      !! the same expression every formula setter uses.
      integer, intent(in) :: i
      real(wp) :: x
      x = (real(i - NGHOST, wp) - 0.5_wp)*DX
   end function cell_x

   ! ------------------------------------------------------------------
   ! Cases
   ! ------------------------------------------------------------------

   subroutine test_parse_enums(error)
      !! Both spellings parse, and an unknown one returns INVALID rather
      !! than a silent default — the caller then fails loud.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_cavity_draft_config("none") == CAVITY_DRAFT_NONE, "none")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_config("flat") == CAVITY_DRAFT_FLAT, "flat")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_config("linear") == CAVITY_DRAFT_LINEAR, "linear")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_config("file") == CAVITY_DRAFT_FILE, "file")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_config("wedge") == CAVITY_DRAFT_INVALID, &
                 "an unknown draft_config must be INVALID, never a silent default")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_source("draft") == CAVITY_SOURCE_DRAFT, "draft")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_source("thickness") == CAVITY_SOURCE_THICKNESS, &
                 "thickness")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_source("in_situ") == CAVITY_SOURCE_IN_SITU, &
                 "in_situ")
      if (allocated(error)) return
      call check(error, parse_cavity_draft_source("guess") == CAVITY_SOURCE_INVALID, &
                 "an unknown draft_source must be INVALID")
   end subroutine test_parse_enums

   subroutine test_flat_setter(error)
      !! `flat`: the draft inside the box, ZERO beyond the calving front,
      !! and GHOST ROWS FILLED BY THE FORMULA — not left at the alloc-time
      !! zero, which would put a phantom front one cell outside the wall
      !! (the formula-bathymetry ghost trap in its cavity form).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: z(NX_PHYS + 2*NGHOST, NY_PHYS + 2*NGHOST)
      integer :: i, j
      logical :: ghost_ok, front_ok, inside_ok

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      z = -1.0_wp   ! poison: an unwritten cell must fail, not read as 0
      ! Shelf unbounded to the west/south/north, calving front at 4000 m.
      call set_draft_flat(z, grid, DRAFT, -1.0e30_wp, 4000.0_wp, -1.0e30_wp, 1.0e30_wp)

      inside_ok = .true.
      front_ok = .true.
      do j = 1, NY_PHYS + 2*NGHOST
         do i = 1, NX_PHYS + 2*NGHOST
            if (cell_x(i) <= 4000.0_wp) then
               if (z(i, j) /= DRAFT) inside_ok = .false.
            else
               if (z(i, j) /= 0.0_wp) front_ok = .false.
            end if
         end do
      end do
      call check(error, inside_ok, "flat draft must equal draft_depth everywhere "// &
                 "inside the shelf box")
      if (allocated(error)) return
      call check(error, front_ok, "flat draft must be exactly 0 beyond the calving front")
      if (allocated(error)) return

      ! The west ghost columns (global x = -1500, -500) are INSIDE an
      ! unbounded-west box, so they carry the draft; the east ghosts
      ! (8500, 9500) are beyond the front and carry 0.  Either way they
      ! were WRITTEN — the poison value is gone.
      ghost_ok = .true.
      do j = 1, NY_PHYS + 2*NGHOST
         do i = 1, NGHOST
            if (z(i, j) /= DRAFT) ghost_ok = .false.
         end do
         do i = NX_PHYS + NGHOST + 1, NX_PHYS + 2*NGHOST
            if (z(i, j) /= 0.0_wp) ghost_ok = .false.
         end do
      end do
      call check(error, ghost_ok, "the setter must fill the FULL array including "// &
                 "ghost rows, by evaluating the formula at the ghost position")
      if (allocated(error)) return

      ! And the sentinel really is a sentinel: an unbounded front covers
      ! every cell, ghosts included.
      z = -1.0_wp
      call set_draft_flat(z, grid, DRAFT, -CAVITY_BOUND_INF, CAVITY_BOUND_INF, &
                          -CAVITY_BOUND_INF, CAVITY_BOUND_INF)
      call check(error, all(z == DRAFT), &
                 "a box at the CAVITY_BOUND_INF sentinel must cover the whole array")
   end subroutine test_flat_setter

   subroutine test_linear_setter(error)
      !! `linear`: `d0 + s*(x - x0)` inside the box, clipped at 0 from
      !! below (a formula that lifts the ice base above the sea surface is
      !! open water, not negative ice), 0 beyond the front.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: z(NX_PHYS + 2*NGHOST, NY_PHYS + 2*NGHOST)
      real(wp) :: expect, worst
      integer :: i, j

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      z = -1.0_wp
      ! 100 m at x = 0, deepening by 0.05 m/m eastward, front at 6000 m.
      call set_draft_linear(z, grid, 100.0_wp, 0.05_wp, 0.0_wp, 6000.0_wp, &
                            -1.0e30_wp, 1.0e30_wp)

      worst = 0.0_wp
      do j = 1, NY_PHYS + 2*NGHOST
         do i = 1, NX_PHYS + 2*NGHOST
            if (cell_x(i) >= 0.0_wp .and. cell_x(i) <= 6000.0_wp) then
               expect = 100.0_wp + 0.05_wp*(cell_x(i) - 0.0_wp)
            else
               expect = 0.0_wp
            end if
            worst = max(worst, abs(z(i, j) - expect))
         end do
      end do
      ! Bound: the setter evaluates the SAME expression the check does,
      ! but an FMA-contracting build may fuse `d0 + s*(x-x0)` in one of
      ! them and not the other, so this is a round-off bound on a value of
      ! order 400 m, not a bit-equality.
      call check(error, worst <= 8.0_wp*epsilon(1.0_wp)*400.0_wp, &
                 "linear draft must follow d0 + slope*(x - x0) inside the box "// &
                 "and be 0 outside it")
      if (allocated(error)) return

      ! Clip: a slope steep enough to drive the formula negative must give
      ! exactly 0 (open water), never a negative draft.
      z = -1.0_wp
      call set_draft_linear(z, grid, 100.0_wp, -1.0_wp, 0.0_wp, 1.0e30_wp, &
                            -1.0e30_wp, 1.0e30_wp)
      call check(error, all(z >= 0.0_wp), &
                 "the linear profile must clip at 0 — a negative draft is not ice")
      if (allocated(error)) return
      call check(error, z(NX_PHYS + NGHOST, 1) == 0.0_wp, &
                 "the clipped far end must be exactly 0")
   end subroutine test_linear_setter

   subroutine test_units_are_grid_units(error)
      !! THE UNITS TRAP.  The knobs are metres; the setters work in GRID
      !! coordinate units, which are DEGREES on a spherical grid.  A
      !! calving front at 200 km must therefore land at 200e3/M_PER_DEG
      !! degrees — about 1.8 deg — and NOT at "200000 degrees", which
      !! would cover the whole planet and hide the front entirely.
      !!
      !! This is the same failure mode `&ocean_topo_nml slope_scale`
      !! documents (a metres length against a degrees position collapsed
      !! the seamount to a flat basin); here it would silently delete the
      !! open ocean.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=:), allocatable :: nml
      integer :: ierr, i, ng, i_last_ice, i_first_open
      real(wp) :: front_deg

      ! 8 x 6 cells of 1 DEGREE each; the metres front at 200 km must
      ! land ~1.8 degrees in, i.e. inside cell 2.
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1.0, dy = 1.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, "// &
            "lat_south = 0.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 8 /"//new_line("a")// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', "// &
            "draft_depth = 300.0, draft_x1 = 200.0e3 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")

      call build_seeded(cfg, grid, state, nml, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "spherical cavity namelist must seed")
      if (allocated(error)) return

      front_deg = 200.0e3_wp/M_PER_DEG
      ng = NGHOST
      i_last_ice = 0
      i_first_open = 0
      do i = 1, NX_PHYS + 2*ng
         if (state%metrics%z_draft(i, ng + 1) > 0.0_wp) i_last_ice = i
         if (i_first_open == 0 .and. i > ng .and. &
             state%metrics%z_draft(i, ng + 1) == 0.0_wp) i_first_open = i
      end do
      ! Cell centres in DEGREES: (i - ng - 0.5)*1.0.
      call check(error, (real(i_last_ice - ng, wp) - 0.5_wp) <= front_deg, &
                 "the last ice column must sit west of the converted front")
      if (allocated(error)) return
      call check(error, (real(i_first_open - ng, wp) - 0.5_wp) > front_deg, &
                 "the first open column must sit east of the converted front")
      if (allocated(error)) return
      ! The witness that the conversion happened at all: unconverted, the
      ! front would be at 200000 "degrees" and there would be NO open water.
      call check(error, i_first_open <= NX_PHYS + ng, &
                 "an unconverted metres front would ice the whole domain — the "// &
                 "metres->grid-units conversion is missing")
   end subroutine test_units_are_grid_units

   subroutine test_datum_and_init_balance(error)
      !! (D) + (I) + the resting balance: `bt_H_ref = b - z_draft`,
      !! `sum(h_layer) = b - z_draft`, hence `bt_eta = 0` at t = 0.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr, i, j, nxt, nyt
      real(wp) :: resid, worst_sum, tol_datum, tol_sum

      call build_configured(cfg, grid, state, nml_case(cavity_on()), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "flat-lid cavity must configure")
      if (allocated(error)) return

      nxt = size(state%barotropic%b, 1)
      nyt = size(state%barotropic%b, 2)
      ! (I), with the common positive factor rho_ref*g divided out:
      ! max |bt_H_ref - (b - z_draft)|, in metres of reference depth.
      ! Bound: a few ulp of the largest depth in play — the two sides
      ! evaluate the same difference in two places and a contracting build
      ! may round them differently, so this is not asserted as bit-zero.
      resid = cavity_datum_residual(state%dyn%bt_work%bt_H_ref, state%barotropic%b, &
                                    state%metrics%z_draft, nxt, nyt)
      tol_datum = 8.0_wp*epsilon(1.0_wp)*BED
      call check(error, resid <= tol_datum, &
                 "the counted-once invariant bt_H_ref == b - z_draft must hold")
      if (allocated(error)) return
      call check(error, abs(state%dyn%bt_work%bt_H_ref(NGHOST + 1, NGHOST + 1) - &
                            (BED - DRAFT)) <= tol_datum, &
                 "bt_H_ref must be the WATER column (700 m), not the bed (1000 m)")
      if (allocated(error)) return

      ! sum_k h_layer == b - z_draft  =>  bt_eta == 0 at rest.  The sum of
      ! NZ equal pieces of (b - z_draft)/NZ carries at most NZ roundings.
      worst_sum = 0.0_wp
      do j = 1, nyt
         do i = 1, nxt
            worst_sum = max(worst_sum, abs(sum(state%multilayer%h_layer(i, j, :)) - &
                                           state%dyn%bt_work%bt_H_ref(i, j)))
         end do
      end do
      tol_sum = real(NZ + 1, wp)*epsilon(1.0_wp)*BED
      call check(error, worst_sum <= tol_sum, &
                 "bt_eta = sum(h_layer) - bt_H_ref must be zero at t = 0 under the ice")
      if (allocated(error)) return

      ! The barotropic water-column prognostic follows the same datum.
      call check(error, state%barotropic%h(NGHOST + 1, NGHOST + 1) == BED - DRAFT, &
                 "barotropic h must be seeded from b - z_draft")
   end subroutine test_datum_and_init_balance

   subroutine test_sloping_bt_eta_zero(error)
      !! The same resting balance under a SLOPING draft, where every
      !! column has a different water thickness — the case a constant
      !! `bt_H_ref` would pass by accident.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr, i, j, nxt, nyt
      real(wp) :: worst, spread, tol

      call build_configured(cfg, grid, state, nml_case( &
                            "enable = .true., draft_config = 'linear', "// &
                            "draft_depth = 100.0, draft_slope = 0.05, "// &
                            "draft_x0 = 0.0"), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "sloping-lid cavity must configure")
      if (allocated(error)) return

      nxt = size(state%barotropic%b, 1)
      nyt = size(state%barotropic%b, 2)
      worst = 0.0_wp
      do j = 1, nyt
         do i = 1, nxt
            worst = max(worst, abs(sum(state%multilayer%h_layer(i, j, :)) - &
                                   state%dyn%bt_work%bt_H_ref(i, j)))
         end do
      end do
      tol = real(NZ + 1, wp)*epsilon(1.0_wp)*BED
      call check(error, worst <= tol, &
                 "bt_eta must be zero at t = 0 under a SLOPING draft too")
      if (allocated(error)) return

      ! Non-vacuity: the draft really does vary across the domain, so the
      ! assertion above is not "constant == constant".
      spread = maxval(state%metrics%z_draft(NGHOST + 1:NGHOST + NX_PHYS, NGHOST + 1)) - &
               minval(state%metrics%z_draft(NGHOST + 1:NGHOST + NX_PHYS, NGHOST + 1))
      call check(error, spread > 100.0_wp, &
                 "the test draft must actually slope (else the balance is trivial)")
   end subroutine test_sloping_bt_eta_zero

   subroutine test_grounding_marks_land(error)
      !! GROUNDING: a column with less than `h_min_cavity` of water under
      !! the ice is LAND — `wet_mask == 0` — and it gets there through the
      !! SAME `seed_wet_mask_impl` the bathymetry uses, so the static
      !! metric-zeroing mask follows for free.  Never a thin film.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr, ng
      real(wp) :: water_west, water_east

      ! Seamount bathymetry: shallow in the middle, deep at the edges.
      ! A 700 m flat lid grounds the shallow centre and leaves the deep
      ! rim afloat.
      call build_configured(cfg, grid, state, nml_case( &
                            "enable = .true., draft_config = 'flat', "// &
                            "draft_depth = 700.0, h_min_cavity = 50.0, "// &
                            "grounded_max_frac = 0.9", &
                            topo="&ocean_topo_nml topo_config = 'seamount', "// &
                            "max_depth = 1000.0, edge_depth = 200.0, "// &
                            "slope_scale = 2000.0 /"), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "grounded cavity must configure")
      if (allocated(error)) return

      ng = NGHOST
      water_west = state%barotropic%b(ng + 1, ng + 1) - state%metrics%z_draft(ng + 1, ng + 1)
      water_east = state%barotropic%b(ng + NX_PHYS/2, ng + NY_PHYS/2) - &
                   state%metrics%z_draft(ng + NX_PHYS/2, ng + NY_PHYS/2)
      ! Sanity on the geometry the assertions below rest on.
      call check(error, water_west >= 50.0_wp, &
                 "the deep corner must stay afloat for this test to mean anything")
      if (allocated(error)) return
      call check(error, water_east < 50.0_wp, &
                 "the seamount centre must be grounded for this test to mean anything")
      if (allocated(error)) return

      call check(error, state%multilayer%wet_mask(ng + 1, ng + 1) == 1.0_wp, &
                 "an afloat column must stay wet")
      if (allocated(error)) return
      call check(error, state%multilayer%wet_mask(ng + NX_PHYS/2, ng + NY_PHYS/2) == 0.0_wp, &
                 "a GROUNDED column must be LAND (wet_mask == 0) — never a thin film "// &
                 "of water under grounded ice")
   end subroutine test_grounding_marks_land

   subroutine test_no_ice_over_land(error)
      !! NO ICE OVER LAND: where the bathymetry already says land, the
      !! draft is forced to 0, so a land column keeps the datum it always
      !! had (`bt_H_ref = b`) and the counted-once invariant stays exact
      !! on EVERY column — no `merge` anywhere downstream.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr, i, j, nxt, nyt
      logical :: clean
      integer :: n_land

      ! `island`: a central LAND square in a flat basin, under a lid that
      ! covers the whole domain.
      call build_configured(cfg, grid, state, nml_case( &
                            "enable = .true., draft_config = 'flat', "// &
                            "draft_depth = 300.0, grounded_max_frac = 0.9", &
                            topo="&ocean_topo_nml topo_config = 'island', "// &
                            "max_depth = 1000.0, slope_scale = 0.3 /"), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "island + cavity must configure")
      if (allocated(error)) return

      nxt = size(state%barotropic%b, 1)
      nyt = size(state%barotropic%b, 2)
      clean = .true.
      n_land = 0
      do j = 1, nyt
         do i = 1, nxt
            if (state%barotropic%b(i, j) < LAND_DEPTH_THRESHOLD) then
               n_land = n_land + 1
               if (state%metrics%z_draft(i, j) /= 0.0_wp) clean = .false.
            end if
         end do
      end do
      call check(error, n_land > 0, "the island bathymetry must actually contain land")
      if (allocated(error)) return
      call check(error, clean, "the draft must be zeroed on every pre-existing land "// &
                 "column (no ice over land)")
      if (allocated(error)) return
      call check(error, cavity_datum_residual(state%dyn%bt_work%bt_H_ref, &
                                              state%barotropic%b, &
                                              state%metrics%z_draft, nxt, nyt) &
                 <= 8.0_wp*epsilon(1.0_wp)*BED, &
                 "the datum invariant must hold on land columns too")
   end subroutine test_no_ice_over_land

   subroutine test_grounded_fraction(error)
      !! The sanity bound fails LOUD: a draft deeper than the whole basin
      !! grounds everything, and configure refuses rather than running a
      !! domain that is almost entirely land.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr

      call build_seeded(cfg, grid, state, nml_case( &
                        "enable = .true., draft_config = 'flat', "// &
                        "draft_depth = 995.0, h_min_cavity = 10.0, "// &
                        "grounded_max_frac = 0.5"), ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, &
                 "a draft that grounds the whole domain must FAIL configure, not "// &
                 "silently run a basin of land")
   end subroutine test_grounded_fraction

   subroutine test_p_ice_ref(error)
      !! (P): the load is `(rho_ref*GRAVITY)*z_draft` built as ONE
      !! product — the same one the FV_MOM6 surface BC forms, which is
      !! what lets `pa(nz+1) = rho_ref*g*(-z_draft) + p_ice_ref` cancel to
      !! bit-zero at rest once the load is wired.
      !!
      !! Asserted as BIT equality on purpose, and legitimately: both sides
      !! are the single product `(rho_ref*GRAVITY)*z_draft` with the same
      !! operands in the same order, and a lone multiply has no addend for
      !! an FMA to fuse.  (The moment a `+` appears in either side, this
      !! becomes a round-off bound instead.)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr, i, j, nxt, nyt
      real(wp) :: rho_g
      logical :: exact

      call build_configured(cfg, grid, state, nml_case( &
                            "enable = .true., draft_config = 'linear', "// &
                            "draft_depth = 100.0, draft_slope = 0.05, "// &
                            "draft_x0 = 0.0"), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "cavity must configure")
      if (allocated(error)) return

      nxt = size(state%metrics%z_draft, 1)
      nyt = size(state%metrics%z_draft, 2)
      rho_g = state%pressure_force%rho_ref*GRAVITY
      exact = .true.
      do j = 1, nyt
         do i = 1, nxt
            if (state%metrics%p_ice_ref(i, j) /= rho_g*state%metrics%z_draft(i, j)) then
               exact = .false.
            end if
         end do
      end do
      call check(error, exact, "p_ice_ref must be exactly (rho_ref*GRAVITY)*z_draft")
      if (allocated(error)) return
      call check(error, maxval(state%metrics%p_ice_ref) > 1.0e6_wp, &
                 "the test draft must produce a load worth checking (MPa scale)")
   end subroutine test_p_ice_ref

   subroutine test_off_bit_identical(error)
      !! Knob off ⇒ nothing moved: the three cavity arrays stay at their
      !! `(1,1)` placeholder, `bt_H_ref` is a BYTE copy of `b`, and the
      !! column thickness is the bed.  Structural, not an argument about
      !! `x + 0.0`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ierr

      call build_configured(cfg, grid, state, nml_case(""), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "default namelist must configure")
      if (allocated(error)) return
      call check(error,.not. state%metrics%use_cavity, "use_cavity must default off")
      if (allocated(error)) return
      call check(error, size(state%metrics%z_draft, 1) == 1 .and. &
                 size(state%metrics%z_draft, 2) == 1, &
                 "z_draft must stay at its (1,1) placeholder when the knob is off")
      if (allocated(error)) return
      call check(error, size(state%metrics%p_ice_ref, 1) == 1, &
                 "p_ice_ref must stay at its (1,1) placeholder when the knob is off")
      if (allocated(error)) return
      call check(error, all(state%dyn%bt_work%bt_H_ref == state%barotropic%b), &
                 "bt_H_ref must be a byte copy of b with the cavity off")
      if (allocated(error)) return
      call check(error, all(state%barotropic%h == state%barotropic%b), &
                 "barotropic h must be a byte copy of b with the cavity off")
   end subroutine test_off_bit_identical

   subroutine test_bytes_accounting(error)
      !! Every new array is counted in its owning type's `bytes()`: the
      !! cavity-on total must exceed the cavity-off total by exactly three
      !! full 2-D fields, less the three placeholders it replaced.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg_on, cfg_off
      type(hgrid_t) :: grid_on, grid_off
      type(ocean_state_t) :: on, off
      integer :: ierr
      integer :: nxt, nyt
      integer(kind=8) :: delta, expect

      call build_seeded(cfg_off, grid_off, off, nml_case(""), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "cavity-off state must seed")
      if (allocated(error)) return
      call build_seeded(cfg_on, grid_on, on, nml_case(cavity_on()), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "cavity-on state must seed")
      if (allocated(error)) return

      nxt = NX_PHYS + 2*NGHOST
      nyt = NY_PHYS + 2*NGHOST
      delta = on%metrics%bytes() - off%metrics%bytes()
      expect = 3_8*int(storage_size(1.0_wp)/8, 8)*(int(nxt, 8)*int(nyt, 8) - 1_8)
      call check(error, delta == expect, &
                 "metrics%bytes() must account for z_draft + cover_frac + p_ice_ref")
   end subroutine test_bytes_accounting

   subroutine test_restart_registry(error)
      !! The draft is PRESCRIBED STATIC geometry, rebuilt at configure
      !! from the namelist — like `b` and `bt_H_ref`, and excluded from
      !! the restart registry for the same reason: checkpointing it would
      !! create a file-vs-namelist ambiguity, and a resume whose draft
      !! disagreed with its datum would be silently wrong.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg_on, cfg_off
      type(hgrid_t) :: grid_on, grid_off
      type(ocean_state_t) :: on, off
      type(restart_registry_t) :: reg_on, reg_off
      integer :: ierr, e
      logical :: found

      call build_seeded(cfg_off, grid_off, off, nml_case(""), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "cavity-off state must seed")
      if (allocated(error)) return
      call build_seeded(cfg_on, grid_on, on, nml_case(cavity_on()), ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "cavity-on state must seed")
      if (allocated(error)) return

      call ocean_state_build_restart_registry(off, grid_off, reg_off)
      call ocean_state_build_restart_registry(on, grid_on, reg_on)
      call check(error, reg_on%n == reg_off%n, &
                 "enabling the cavity must not add a restart entry")
      if (allocated(error)) return
      found = .false.
      do e = 1, reg_on%n
         if (index(reg_on%entries(e)%tag, "draft") > 0 .or. &
             index(reg_on%entries(e)%tag, "cavity") > 0 .or. &
             index(reg_on%entries(e)%tag, "p_ice") > 0) found = .true.
      end do
      call check(error,.not. found, &
                 "no cavity field may be registered for restart (static, rebuilt "// &
                 "at configure)")
   end subroutine test_restart_registry

   subroutine test_validate_matrix(error)
      !! The refusal matrix, one row at a time.  Each row must FAIL
      !! validation; the base config must PASS, so a row that fails for an
      !! unrelated reason cannot hide.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: on

      on = cavity_on()

      call expect_valid(error, nml_case(on), "the v1 envelope itself must validate")
      if (allocated(error)) return

      ! --- geometry sources that are named but not implemented ---
      call expect_invalid(error, nml_case("enable = .true., draft_config = 'file'"), &
                          "draft_config='file'")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on//", draft_source = 'in_situ'"), &
                          "draft_source='in_situ'")
      if (allocated(error)) return
      call expect_invalid(error, nml_case("enable = .true., draft_config = 'linear', "// &
                                          "draft_depth = 100.0, draft_slope = 0.01"), &
                          "linear with no finite draft_x0 anchor")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on//", h_min_cavity = 0.0"), "h_min_cavity = 0")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on//", grounded_max_frac = 1.5"), &
                          "grounded_max_frac outside (0, 1]")
      if (allocated(error)) return

      ! --- the load must have a consumer once it has a gradient (P5.2) ---
      ! A draft that VARIES needs `&ocean_pgf_nml p_top_in_bc`: that knob
      ! is the only route by which rho_ref*g*z_draft reaches the FV_MOM6
      ! pa(nz+1) surface BC.  Refused, not auto-enabled.
      call expect_invalid(error, nml_case("enable = .true., draft_config = 'linear', "// &
                                          "draft_depth = 100.0, draft_slope = 0.01, "// &
                                          "draft_x0 = 0.0", &
                                          pgf="&ocean_pgf_nml form = 'fv_mom6' /"), &
                          "a sloping draft without p_top_in_bc")
      if (allocated(error)) return
      call expect_invalid(error, nml_case("enable = .true., draft_config = 'flat', "// &
                                          "draft_depth = 300.0, draft_x1 = 4000.0", &
                                          pgf="&ocean_pgf_nml form = 'fv_mom6' /"), &
                          "a flat draft with a CALVING FRONT (a step) and no p_top_in_bc")
      if (allocated(error)) return
      ! ... but a draft that is UNIFORM over the whole array is EXEMPT: a
      ! load with no gradient is bit-identically inert in the top BC, and
      ! that exemption is what keeps the unloaded flat-lid datum-
      ! equivalence gate expressible at all.
      call expect_valid(error, nml_case(on, pgf="&ocean_pgf_nml form = 'fv_mom6' /"), &
                        "a UNIFORM draft without p_top_in_bc (provably inert)")
      if (allocated(error)) return

      ! --- pressure-gradient envelope ---
      call expect_invalid(error, nml_case(on, pgf="&ocean_pgf_nml form = 'mont' /"), &
                          "PGF form /= fv_mom6")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, pgf="&ocean_pgf_nml form = 'fv_mom6', "// &
                                          "gfs_scale = 0.9 /"), "gfs_scale /= 1")
      if (allocated(error)) return

      ! --- vertical-coordinate envelope (every z-like family) ---
      call expect_invalid(error, nml_case(on, vcoord="&vcoord_nml vcoord_type = "// &
                                          "'zsigma' /"), "vcoord zsigma")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, vcoord="&vcoord_nml vcoord_type = "// &
                                          "'zstar_full' /"), "vcoord zstar_full")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, vcoord="&vcoord_nml vcoord_type = "// &
                                          "'z_fixed' /"), "vcoord z_fixed")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, vcoord="&vcoord_nml vcoord_type = "// &
                                          "'eulerian_z' /"), "vcoord eulerian_z")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, vcoord="&vcoord_nml vcoord_type = "// &
                                          "'sigma', thickness_config = 'uniform_z' /"), &
                          "thickness_config uniform_z")
      if (allocated(error)) return

      ! --- solver envelope ---
      call expect_invalid(error, nml_case(on, bt="&ocean_bt_nml auto_n_inner = "// &
                                          ".false., n_inner = 0 /"), "the unsplit driver")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, bt="&ocean_bt_nml auto_n_inner = "// &
                                          ".false., n_inner = 8, bt_halo = 8 /"), &
                          "bt_halo > 0")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, extra="&mpi_nml px = 2, py = 1 /"), &
                          "multi-rank")
      if (allocated(error)) return

      ! --- mutually exclusive capabilities ---
      call expect_invalid(error, nml_case(on, extra="&ocean_wetdry_nml enable = "// &
                                          ".true., dry_depth = 0.1, rewet_depth = 0.2 /"), &
                          "wet/dry")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, extra="&ocean_porous_nml enable = .true. /"), &
                          "porous barriers")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, extra="&ocean_ice_nml enable = .true. /"), &
                          "sea ice")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, extra="&ocean_tides_nml enable = .true., "// &
                                          "use_sal = .true. /"), "tidal SAL")
      if (allocated(error)) return
      call expect_invalid(error, nml_case(on, extra="&ocean_zinit_nml enable = .true., "// &
                                          "file = 'ts.nc' /"), "the z-level T/S IC")
      if (allocated(error)) return

      ! --- and the knob-off path stays acceptable everywhere ---
      call expect_valid(error, nml_case("", vcoord="&vcoord_nml vcoord_type = "// &
                                        "'zstar_full' /"), &
                        "with the cavity OFF every refused partner is fine again")
   end subroutine test_validate_matrix

   subroutine expect_valid(error, nml, what)
      !! `validate_config` reports through `ierr` when it is present
      !! (without it, it `error stop`s) — so every row here passes it.
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, what
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call check(error, .false., what//": namelist must parse")
         return
      end if
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, what)
   end subroutine expect_valid

   subroutine expect_invalid(error, nml, what)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, what
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return   ! a parse-level rejection is a refusal too
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "the cavity must be REFUSED with "//what)
   end subroutine expect_invalid

end module test_ocean_cavity_draft
