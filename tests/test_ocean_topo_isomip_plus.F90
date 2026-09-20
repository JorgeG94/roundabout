!! Analytic gates for the MISMIP+/ISOMIP+ bedrock formula bathymetry
!! (`&ocean_topo_nml topo_config = "isomip_plus"`).
!!
!! Source of truth: Asay-Davis et al. (2016), Geosci. Model Dev. 9,
!! 2471-2497.  The bed is their Eq. (1),
!!
!!     z_b(x,y) = max( Bx(x) + By(y), z_b,deep )
!!     Bx(x)    = B0 + B2*xt^2 + B4*xt^4 + B6*xt^6,   xt = x/x_bar   (2)
!!     By(y)    = d_c/(1 + exp(-2*(y - Ly/2 - w_c)/f_c))
!!              + d_c/(1 + exp( 2*(y - Ly/2 + w_c)/f_c))             (4)
!!
!! with Table 1: B0 = -150.0, B2 = -728.8, B4 = 343.91, B6 = -50.57 m,
!! x_bar = 300 km, d_c = 500 m, f_c = 4.0 km, w_c = 24.0 km,
!! z_b,deep = -720 m; and Table 3 for the ISOMIP+ ocean box,
!! x0 = 320 km, Lx = 480 km, Ly = 80 km.
!!
!! Every expected number below is derived from those printed constants
!! INDEPENDENTLY of the Fortran (the arithmetic is written out in the
!! comment that precedes each), so the test fails if the code's
!! coefficients, the polynomial nesting, the trough geometry, the deep
!! clip or the elevation->depth sign flip drifts.
module test_ocean_topo_isomip_plus
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: set_bathymetry_isomip_plus, &
                              isomip_plus_bx, isomip_plus_by, &
                              ISOMIP_ZB_DEEP, ISOMIP_DC, ISOMIP_WC
   implicit none
   private

   public :: collect_ocean_topo_isomip_plus_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: X0_ISOMIP = 320.0e3_wp
      !! Table 3 `x0` — southern (west, in model terms) edge of the
      !! ISOMIP+ ocean box on the MISMIP+ absolute x axis.
   real(wp), parameter :: DX_2KM = 2000.0_wp
      !! Table 4 COM resolution.
   real(wp), parameter :: MAX_DEPTH = 720.0_wp
      !! `-z_b,deep` — the Eq. (1) deep clip as a positive-down depth.
   real(wp), parameter :: TOL_M = 1.0e-6_wp
      !! Metres.  The expected values below are quoted to 8 decimals, so
      !! a micrometre tolerance is ~2 orders tighter than the quote and
      !! still far above double-precision round-off on ~1e3 m operands.

contains

   subroutine collect_ocean_topo_isomip_plus_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bx_matches_hand_evaluated_polynomial", test_bx_values), &
                  new_unittest("by_is_the_two_sided_logistic_trough", test_by_values), &
                  new_unittest("bed_matches_hand_evaluated_points", test_bed_points), &
                  new_unittest("deep_end_is_clipped_at_zb_deep", test_deep_clip), &
                  new_unittest("trough_is_deeper_than_side_walls", test_trough_deeper), &
                  new_unittest("fills_ghost_rows", test_ghost_fill), &
                  new_unittest("x_origin_shifts_the_profile", test_x_origin_shift), &
                  new_unittest("bed_above_sea_level_is_land", test_land_above_sea_level) &
                  ]
   end subroutine collect_ocean_topo_isomip_plus_tests

   function make_grid(nx, ny) result(g)
      !! Uniform 2 km Cartesian grid, the ISOMIP+ COM resolution.
      integer, intent(in) :: nx, ny
      type(hgrid_t) :: g
      call g%init(nx, ny, NGHOST, DX_2KM, DX_2KM)
   end function make_grid

   subroutine test_bx_values(error)
      !! `Bx` at four absolute x, each evaluated by hand from Eq. (2).
      !! With xt = x/300e3:
      !!
      !!  x = 320 km: xt = 16/15 = 1.0666666666666667
      !!    xt^2 = 1.1377777777777778
      !!    xt^4 = 1.2945393909465022
      !!    xt^6 = 1.4728803274907668
      !!    Bx = -150 - 728.8*1.1377777777777778
      !!             + 343.91*1.2945393909465022
      !!             -  50.57*1.4728803274907668
      !!       = -150 - 829.2522666666667 + 445.20609043...
      !!             -  74.48334216...
      !!       = -608.49218257 m
      !!
      !!  x = 450 km: xt = 1.5, xt^2 = 2.25, xt^4 = 5.0625, xt^6 = 11.390625
      !!    Bx = -150 - 1639.8 + 1741.0443750 - 576.02390625
      !!       = -624.77953125 m   (exact in binary64 to the digits shown)
      !!
      !!  x = 600 km: xt = 2, xt^2 = 4, xt^4 = 16, xt^6 = 64
      !!    Bx = -150 - 2915.2 + 5502.56 - 3236.48 = -799.12 m
      !!         (below z_b,deep, so Eq. (1) clips — see test_deep_clip)
      !!
      !!  x = 0:      Bx = B0 = -150 m
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, abs(isomip_plus_bx(0.0_wp) + 150.0_wp) < TOL_M, &
                    "Bx(0) must be B0 = -150 m")
         if (allocated(error)) exit checks
         call check(error, abs(isomip_plus_bx(320.0e3_wp) + 608.49218257_wp) < 1.0e-5_wp, &
                    "Bx(320 km) = -608.49218257 m")
         if (allocated(error)) exit checks
         call check(error, abs(isomip_plus_bx(450.0e3_wp) + 624.77953125_wp) < TOL_M, &
                    "Bx(450 km) = -624.77953125 m")
         if (allocated(error)) exit checks
         call check(error, abs(isomip_plus_bx(600.0e3_wp) + 799.12_wp) < 1.0e-9_wp, &
                    "Bx(600 km) = -799.12 m")
      end block checks
   end subroutine test_bx_values

   subroutine test_by_values(error)
      !! `By` from Eq. (4) on the prescribed Ly = 80 km box.
      !!
      !!  y = Ly/2 = 40 km (trough centre): yc = 0, so both terms are
      !!    d_c/(1 + exp(2*w_c/f_c)) = 500/(1 + exp(12)).
      !!    exp(12) = 162754.79141900392, so each term is
      !!    3.0720873011...e-3 and By = 6.144174602e-3 m — the trough
      !!    floor is `Bx` to within 6 mm, which is the point of Fig. 1b's
      !!    caveat that `By` is an offset, not a transect.
      !!
      !!  y = 0 (side wall): term 1 = 500/(1 + exp(32)) ~ 6.4e-12 ~ 0;
      !!    term 2 = 500/(1 + exp(-8)) = 500/1.00033546262790 = 499.83232493 m.
      !!    By(0) = 499.8323249348 m.
      !!
      !!  Symmetry: By(y) = By(Ly - y) exactly (the two terms swap).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: LY = 80.0e3_wp
      checks: block
         call check(error, abs(isomip_plus_by(40.0e3_wp, LY) - 6.144174602e-3_wp) < 1.0e-9_wp, &
                    "By(Ly/2) = 6.144174602e-3 m (trough floor ~ 0)")
         if (allocated(error)) exit checks
         call check(error, abs(isomip_plus_by(0.0_wp, LY) - 499.8323249348_wp) < 1.0e-8_wp, &
                    "By(0) = 499.8323249348 m (side wall, ~d_c)")
         if (allocated(error)) exit checks
         call check(error, abs(isomip_plus_by(0.0_wp, LY) - isomip_plus_by(LY, LY)) < 1.0e-12_wp, &
                    "By must be symmetric about Ly/2")
         if (allocated(error)) exit checks
         ! Far outside the trough the logistic saturates at d_c; well
         ! inside it it is ~0.  Both bracket the trough half-width w_c.
         call check(error, isomip_plus_by(0.5_wp*LY - ISOMIP_WC, LY) > 0.4_wp*ISOMIP_DC .and. &
                    isomip_plus_by(0.5_wp*LY - ISOMIP_WC, LY) < 0.6_wp*ISOMIP_DC, &
                    "By at |y-Ly/2| = w_c must be ~d_c/2 (logistic midpoint)")
      end block checks
   end subroutine test_by_values

   subroutine test_bed_points(error)
      !! The assembled setter, on the 240 x 40 ISOMIP+ COM grid with
      !! `x_origin = 320 km`, at three cell centres.  Cell centres are
      !! `x = x_origin + (i_phys - 0.5)*dx`, `y = (j_phys - 0.5)*dy`.
      !!
      !!  (i_phys, j_phys) = (1, 20): x = 321 km, y = 39 km
      !!    Bx(321 km) = -609.49919809 m  (xt = 1.07, xt^2 = 1.1449,
      !!      xt^4 = 1.31079601, xt^6 = 1.500580...;
      !!      -150 - 834.4055... + 450.8146... + ... = -609.49919809)
      !!    By(39 km) = 6.9283151e-3 m   (yc = -1 km: terms
      !!      500/(1+exp(12.5)) + 500/(1+exp(11.5)))
      !!    z_b = -609.49226978 m, above the clip, so
      !!    b = 609.49226978 m.
      !!
      !!  (i_phys, j_phys) = (1, 1): x = 321 km, y = 1 km
      !!    By(1 km) = 499.7236106815 m  (yc = -39 km)
      !!    z_b = -609.49919809 + 499.72361068 = -109.77558741
      !!    b = 109.77558741 m — the shallow side-wall shelf.
      !!
      !!  (i_phys, j_phys) = (66, 20): x = 451 km, y = 39 km
      !!    Bx(451 km) = -624.27112988 m
      !!    z_b = -624.26420157, b = 624.26420157 m.  451 km is one cell
      !!    past the prescribed steady-state grounding line at
      !!    x = 450 +/- 10 km (their Sect. 2.1).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng

      grid = make_grid(240, 40)
      ng = grid%nghost
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_isomip_plus(b, grid, MAX_DEPTH, X0_ISOMIP, 1.0_wp)

      checks: block
         call check(error, abs(b(ng + 1, ng + 20) - 609.49226978_wp) < 1.0e-5_wp, &
                    "b(x=321 km, y=39 km) = 609.49226978 m")
         if (allocated(error)) exit checks
         call check(error, abs(b(ng + 1, ng + 1) - 109.77558741_wp) < 1.0e-5_wp, &
                    "b(x=321 km, y=1 km) = 109.77558741 m (side-wall shelf)")
         if (allocated(error)) exit checks
         call check(error, abs(b(ng + 66, ng + 20) - 624.26420157_wp) < 1.0e-5_wp, &
                    "b(x=451 km, y=39 km) = 624.26420157 m")
      end block checks
   end subroutine test_bed_points

   subroutine test_deep_clip(error)
      !! Eq. (1)'s `max(..., z_b,deep)` must bite at the deep end.  At the
      !! trough centre, `Bx + By` drops below -720 m at x ~ 582.7 km; the
      !! ISOMIP+ box runs to 800 km, where Bx(800 km) = -6126.44 m, so the
      !! entire open-ocean end of the domain sits exactly at the clip.
      !! The clip value is the `max_depth` argument, NOT a hard-coded
      !! 720 — pass a different one and the plateau must follow it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng

      grid = make_grid(240, 40)
      ng = grid%nghost
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_isomip_plus(b, grid, MAX_DEPTH, X0_ISOMIP, 1.0_wp)

      checks: block
         ! i_phys = 240 -> x = 320 + 479 = 799 km, the last wet column.
         call check(error, abs(b(ng + 240, ng + 20) - MAX_DEPTH) < TOL_M, &
                    "deep end must sit exactly at the clip depth")
         if (allocated(error)) exit checks
         ! i_phys = 141 -> x = 320 + 281 = 601 km, already past 582.7 km.
         call check(error, abs(b(ng + 141, ng + 20) - MAX_DEPTH) < TOL_M, &
                    "x = 601 km at the trough centre is clipped")
         if (allocated(error)) exit checks
         ! Nothing anywhere may exceed the clip.
         call check(error, maxval(b) <= MAX_DEPTH + TOL_M, &
                    "no cell may be deeper than the clip")
         if (allocated(error)) exit checks
         ! i_phys = 66 -> x = 451 km: NOT clipped (624.26 m < 720 m).
         call check(error, b(ng + 66, ng + 20) < MAX_DEPTH - 1.0_wp, &
                    "the cavity end must NOT be clipped")
         if (allocated(error)) exit checks
         ! The clip really is the argument: halve it and the plateau moves.
         call set_bathymetry_isomip_plus(b, grid, 0.5_wp*MAX_DEPTH, X0_ISOMIP, 1.0_wp)
         call check(error, abs(maxval(b) - 0.5_wp*MAX_DEPTH) < TOL_M, &
                    "max_depth is the clip, not a hard-coded 720 m")
         if (allocated(error)) exit checks
         call check(error, abs(ISOMIP_ZB_DEEP + MAX_DEPTH) < TOL_M, &
                    "the exported z_b,deep constant must be -720 m")
      end block checks
   end subroutine test_deep_clip

   subroutine test_trough_deeper(error)
      !! The defining geometry: a central trough `d_c = 500 m` deeper
      !! than the side walls, half-width `w_c = 24 km` on an 80 km box.
      !! At a fixed x the depth must therefore be `~d_c` greater at
      !! y = Ly/2 than at y -> 0, and the profile must be symmetric in y.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng, j
      real(wp) :: drop

      grid = make_grid(240, 40)
      ng = grid%nghost
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_isomip_plus(b, grid, MAX_DEPTH, X0_ISOMIP, 1.0_wp)

      checks: block
         ! x = 321 km, unclipped there, so the full d_c offset survives:
         ! 609.49226978 - 109.77558741 = 499.71668237 m.
         drop = b(ng + 1, ng + 20) - b(ng + 1, ng + 1)
         call check(error, abs(drop - 499.71668237_wp) < 1.0e-5_wp, &
                    "trough is d_c (minus the logistic tails) deeper than the wall")
         if (allocated(error)) exit checks
         ! y symmetry across the 40-row box: j_phys and 41 - j_phys.
         do j = 1, 20
            if (abs(b(ng + 1, ng + j) - b(ng + 1, ng + 41 - j)) > 1.0e-9_wp) then
               call check(error, .false., "bed must be symmetric about y = Ly/2")
               exit checks
            end if
         end do
         ! Monotone from the wall into the trough (no spurious ridge).
         do j = 1, 19
            if (b(ng + 1, ng + j + 1) < b(ng + 1, ng + j) - 1.0e-9_wp) then
               call check(error, .false., "bed must deepen monotonically toward the trough")
               exit checks
            end if
         end do
      end block checks
   end subroutine test_trough_deeper

   subroutine test_ghost_fill(error)
      !! CLAUDE.md's standing rule: a formula setter must fill the FULL
      !! array, ghost rows included, or the EOS falls back to rho_0 there
      !! and a spurious density jump lands on every wall-adjacent face.
      !! With `x_origin = 320 km` the west ghost band sits at x = 317,
      !! 319 km, still well inside the region where the bed is below sea
      !! level, so every ghost cell must be strictly positive.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)

      grid = make_grid(60, 40)
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      call set_bathymetry_isomip_plus(b, grid, MAX_DEPTH, X0_ISOMIP, 1.0_wp)

      checks: block
         call check(error, minval(b) > 0.0_wp, &
                    "every cell including ghosts must be a positive depth")
         if (allocated(error)) exit checks
         call check(error, maxval(b) <= MAX_DEPTH + TOL_M, &
                    "every cell must respect the deep clip")
         if (allocated(error)) exit checks
         ! The west ghost column must be the formula at x = 319 km, not a
         ! copy of the first interior column at 321 km.  Bx(319 km) =
         ! -607.47101192, By(39 km) = 0.00692832 ->  607.46408360 m.
         call check(error, abs(b(grid%nghost, grid%nghost + 20) - 607.46408360_wp) < 1.0e-5_wp, &
                    "west ghost column evaluates the formula at x = 319 km")
      end block checks
   end subroutine test_ghost_fill

   subroutine test_x_origin_shift(error)
      !! `x_origin` translates the domain along the paper's absolute x
      !! axis and does nothing else.  Shifting the origin by exactly N
      !! cells must reproduce the same bed shifted by N columns, so a
      !! 40-cell (80 km) origin shift maps column i to column i + 40.
      !! `x_origin = 0` (the default) must reproduce the un-shifted
      !! MISMIP+ bed, which is what every other `topo_config` sees.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b_a(:, :), b_b(:, :)
      integer :: ng, i

      grid = make_grid(120, 40)
      ng = grid%nghost
      allocate (b_a(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (b_b(grid%nx_total, grid%ny_total), source=0.0_wp)

      call set_bathymetry_isomip_plus(b_a, grid, MAX_DEPTH, X0_ISOMIP, 1.0_wp)
      ! 40 cells x 2 km = 80 km further along the axis.
      call set_bathymetry_isomip_plus(b_b, grid, MAX_DEPTH, X0_ISOMIP + 80.0e3_wp, 1.0_wp)

      checks: block
         do i = 1, 80
            if (abs(b_a(ng + i + 40, ng + 20) - b_b(ng + i, ng + 20)) > 1.0e-9_wp) then
               call check(error, .false., "x_origin must translate the bed by whole cells")
               exit checks
            end if
         end do
         ! x_origin = 0 puts the west ghost band at negative x, where the
         ! polynomial is even in x and Bx(-x) = Bx(x): the bed there is
         ! B0 + By = -150 + 499.72 = +349.72 m of ELEVATION, i.e. land.
         call set_bathymetry_isomip_plus(b_a, grid, MAX_DEPTH, 0.0_wp, 1.0_wp)
         call check(error, b_a(ng + 1, ng + 1) == 0.0_wp, &
                    "x_origin = 0 puts the first side-wall column above sea level (land)")
      end block checks
   end subroutine test_x_origin_shift

   subroutine test_land_above_sea_level(error)
      !! A bed the formula puts ABOVE sea level must come back as
      !! `b = 0` (dry land, below `LAND_DEPTH_THRESHOLD = 2 m`, so
      !! `seed_wet_mask_impl` masks it), never as a negative depth — a
      !! negative `b` would seed a negative `h_layer` and NaN the column.
      !! At x = 0, y = 0: z_b = B0 + By(0) = -150 + 499.83 = +349.83 m.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng

      grid = make_grid(60, 40)
      ng = grid%nghost
      allocate (b(grid%nx_total, grid%ny_total), source=-1.0_wp)
      ! Origin at 0 km: the first ~70 columns of the side walls are land.
      call set_bathymetry_isomip_plus(b, grid, MAX_DEPTH, 0.0_wp, 1.0_wp)

      checks: block
         call check(error, minval(b) >= 0.0_wp, "no negative depth anywhere")
         if (allocated(error)) exit checks
         call check(error, b(ng + 1, ng + 1) == 0.0_wp, &
                    "bed above sea level reports b = 0 (land), not a negative depth")
         if (allocated(error)) exit checks
         ! ... while the trough at the same x is still deep water:
         ! Bx(1 km) = -150.00809774, By(39 km) = 0.00692832
         !   -> z_b = -150.00116942, b = 150.00116942 m.
         call check(error, abs(b(ng + 1, ng + 20) - 150.00116942_wp) < 1.0e-5_wp, &
                    "the trough at x = 1 km is 150.00116942 m deep")
      end block checks
   end subroutine test_land_above_sea_level

end module test_ocean_topo_isomip_plus
