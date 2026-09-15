!! Analytical tests for the Neverworld2 idealized-basin setters
!! (`set_bathymetry_neverworld2` + `ocean_surfstress_set_neverworld2`,
!! Marques et al. 2022 GMD; MOM6-inspired).
!!
!! The bathymetry/wind formulas are re-implemented independently here
!! (`ref_cosbell`/`ref_spike`/`ref_dfrac`/`ref_taux`) and compared cell-by-cell
!! to the production setters — a transcription error in any coefficient, sign,
!! or band breakpoint fails the checkvalue test.  Plus structural invariants:
!!   - ghost-row fill (sentinel survives nowhere) — the EOS-fallback footgun;
!!   - depth range [0, max_depth];
!!   - hard N/S walls + an OPEN re-entrant southern channel (the Drake analog);
!!   - aquaplanet (nl_continent_amp=0) has no interior land away from the walls;
!!   - wind τ_y≡0, ghost rows zero, polar τ_x=0, a peak westerly in the south.
module test_ocean_neverworld2
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: set_bathymetry_neverworld2
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   implicit none
   private

   public :: collect_ocean_neverworld2_tests

   integer, parameter :: NGHOST = 3
   real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
   real(wp), parameter :: MAX_DEPTH = 4000.0_wp

contains

   subroutine collect_ocean_neverworld2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nw2_bathy_fills_ghosts", test_bathy_ghost_fill), &
                  new_unittest("nw2_bathy_matches_formula", test_bathy_checkvalue), &
                  new_unittest("nw2_bathy_walls_and_channel", test_walls_channel), &
                  new_unittest("nw2_aquaplanet_no_interior_land", test_aquaplanet), &
                  new_unittest("nw2_min_depth_floor", test_min_depth_floor), &
                  new_unittest("nw2_wind_matches_formula", test_wind_checkvalue), &
                  new_unittest("nw2_wind_structure", test_wind_structure) &
                  ]
   end subroutine collect_ocean_neverworld2_tests

   ! ---- independent reference re-implementation of the documented formula ----

   pure function ref_cosbell(x, L) result(c)
      real(wp), intent(in) :: x, L
      real(wp) :: c
      c = 0.5_wp*(1.0_wp + cos(PI*min(abs(x/L), 1.0_wp)))
   end function ref_cosbell

   pure function ref_spike(x, L) result(s)
      real(wp), intent(in) :: x, L
      real(wp) :: s
      s = 1.0_wp - sin(PI*min(abs(x/L), 0.5_wp))
   end function ref_spike

   pure function ref_dfrac(x, y, ac, ar) result(d)
      real(wp), intent(in) :: x, y, ac, ar
      real(wp) :: d
      d = 1.0_wp &
          - 1.1_wp*ref_spike(y - 1.0_wp, 0.12_wp) &
          - 1.1_wp*ref_spike(y, 0.12_wp) &
          - ac*( &
          (1.2_wp*ref_spike(x, 0.2_wp) + 1.2_wp*ref_spike(x - 1.0_wp, 0.2_wp)) &
          *ref_spike(min(0.0_wp, y - 0.3_wp), 0.2_wp) &
          + 1.2_wp*ref_spike(x - 0.5_wp, 0.2_wp)*ref_spike(min(0.0_wp, y - 0.55_wp), 0.2_wp) &
          + 1.2_wp*(ref_spike(x, 0.12_wp) + ref_spike(x - 1.0_wp, 0.12_wp)) &
          *ref_spike(max(0.0_wp, y - 0.06_wp), 0.12_wp) &
          + 0.1_wp*(ref_cosbell(x, 0.1_wp) + ref_cosbell(x - 1.0_wp, 0.1_wp)) &
          + 0.5_wp*ref_cosbell(x - 0.16_wp, 0.05_wp)*(ref_cosbell(y - 0.18_wp, 0.13_wp)**0.4_wp) &
          + 0.4_wp*(ref_cosbell(x - 0.09_wp, 0.08_wp)**0.4_wp)*ref_cosbell(y - 0.26_wp, 0.05_wp) &
          + 0.4_wp*(ref_cosbell(x - 0.08_wp, 0.08_wp)**0.4_wp)*ref_cosbell(y - 0.1_wp, 0.05_wp)) &
          - ar*cos(14.0_wp*PI*x)*sin(14.0_wp*PI*y) &
          - ar*cos(20.0_wp*PI*x)*cos(20.0_wp*PI*y)
      if (d < 0.0_wp) d = 0.0_wp
   end function ref_dfrac

   pure function ref_taux(y, taux_mag) result(t)
      real(wp), intent(in) :: y, taux_mag
      real(wp) :: t
      real(wp), parameter :: OFF = 0.02_wp
      t = 0.0_wp
      if (y <= 0.29_wp) then
         t = taux_mag*((1.0_wp/0.29_wp)*y - (1.0_wp/(2.0_wp*PI))*sin(2.0_wp*PI*y/0.29_wp))
      else if (y <= 0.8_wp - OFF) then
         t = taux_mag*(0.35_wp + 0.65_wp*cos(PI*(y - 0.29_wp)/(0.51_wp - OFF)))
      else if (y <= 1.0_wp - OFF) then
         t = taux_mag*(1.5_wp*((y - 1.0_wp + OFF) - (0.1_wp/PI)*sin(10.0_wp*PI*(y - 0.8_wp + OFF))))
      end if
   end function ref_taux

   ! ---------------------------------- tests ----------------------------------

   subroutine test_bathy_ghost_fill(error)
      !! Every cell of `b` (interior + all four ghost bands) must be written.
      !! Neverworld2 produces legitimate zeros (land), so a `minval > 0` check
      !! is wrong; instead pre-fill with a negative sentinel and assert it
      !! survives nowhere — the real "ghost left at alloc-zero" guard.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: SENTINEL = -999.0_wp

      call grid%init(40, 70, NGHOST, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=SENTINEL)
      call set_bathymetry_neverworld2(b, grid, MAX_DEPTH, 1.0_wp, 0.05_wp, 0.0_wp)

      call check(error, all(b > SENTINEL + 1.0_wp), &
                 "set_bathymetry_neverworld2: a cell kept the sentinel (ghost not filled)")
      if (allocated(error)) return
      call check(error, minval(b) >= 0.0_wp, "nw2 bathy: negative depth (clamp failed)")
      if (allocated(error)) return
      ! MOM6 only clamps D = max(D,0) — there is NO upper clamp, so the
      ! roughness term (amplitude nl_roughness_amp on each of two cosines)
      ! can lift D_frac above 1 by up to 2*nl_roughness_amp.  Matching MOM6
      ! exactly (it is the ground truth), the setter does the same.
      call check(error, maxval(b) <= MAX_DEPTH*(1.0_wp + 2.0_wp*0.05_wp) + 1.0e-9_wp, &
                 "nw2 bathy: depth exceeds max_depth + roughness headroom")
      deallocate (b)
   end subroutine test_bathy_ghost_fill

   subroutine test_bathy_checkvalue(error)
      !! Cell-by-cell match to the independent reference formula over the full
      !! array — catches any coefficient/sign/term transcription error.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp), parameter :: AC = 1.0_wp, AR = 0.05_wp
      real(wp) :: x, y, expect, maxerr
      integer :: i, j, i_phys, j_phys, ng, nxp, nyp

      call grid%init(44, 66, NGHOST, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_neverworld2(b, grid, MAX_DEPTH, AC, AR, 0.0_wp)

      ng = grid%nghost
      nxp = grid%nx_phys
      nyp = grid%ny_phys
      maxerr = 0.0_wp
      do j = 1, grid%ny_total
         j_phys = j - ng
         y = (real(j_phys, wp) - 0.5_wp)/real(nyp, wp)
         do i = 1, grid%nx_total
            i_phys = i - ng
            x = (real(i_phys, wp) - 0.5_wp)/real(nxp, wp)
            expect = ref_dfrac(x, y, AC, AR)*MAX_DEPTH
            maxerr = max(maxerr, abs(b(i, j) - expect))
         end do
      end do
      call check(error, maxerr < 1.0e-9_wp, "nw2 bathy: deviates from reference formula")
      deallocate (b)
   end subroutine test_bathy_checkvalue

   subroutine test_walls_channel(error)
      !! The northern and southern boundary rows are walls (near-zero / land),
      !! the deep interior reaches ~max_depth, and the far-south interior has
      !! at least one OPEN (wet, deep) cell — the re-entrant Drake channel.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng, i, nxp, j_south
      real(wp) :: south_max

      call grid%init(60, 70, NGHOST, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_neverworld2(b, grid, MAX_DEPTH, 1.0_wp, 0.05_wp, 0.0_wp)
      ng = grid%nghost
      nxp = grid%nx_phys

      ! Deep interior exists somewhere (basin is not all land).
      call check(error, maxval(b) > 0.9_wp*MAX_DEPTH, "nw2: no deep interior")
      if (allocated(error)) return

      ! Southernmost physical row (Antarctica wall band): shallow everywhere.
      call check(error, maxval(b(:, ng + 1)) < 0.5_wp*MAX_DEPTH, &
                 "nw2: southern wall row not shallow")
      if (allocated(error)) return
      ! Northernmost physical row (great northern wall): shallow everywhere.
      call check(error, maxval(b(:, ng + grid%ny_phys)) < 0.5_wp*MAX_DEPTH, &
                 "nw2: northern wall row not shallow")
      if (allocated(error)) return

      ! Far-south interior band (a few rows north of the wall) must contain an
      ! open deep gap: the Drake channel. Scan row at ~y=0.07.
      j_south = ng + max(1, nint(0.07_wp*real(grid%ny_phys, wp)))
      south_max = 0.0_wp
      do i = ng + 1, ng + nxp
         south_max = max(south_max, b(i, j_south))
      end do
      call check(error, south_max > 0.5_wp*MAX_DEPTH, &
                 "nw2: southern channel fully blocked (no open Drake passage)")
      deallocate (b)
   end subroutine test_walls_channel

   subroutine test_aquaplanet(error)
      !! With nl_continent_amp = 0 the only land is the N/S wall bands; the
      !! mid-latitude interior is everywhere wet (no continents/ridges).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      integer :: ng, i, j, jlo, jhi
      logical :: any_dry

      call grid%init(60, 70, NGHOST, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_neverworld2(b, grid, MAX_DEPTH, 0.0_wp, 0.05_wp, 0.0_wp)
      ng = grid%nghost

      ! Interior band away from the walls: y in [0.2, 0.8].
      jlo = ng + nint(0.2_wp*real(grid%ny_phys, wp))
      jhi = ng + nint(0.8_wp*real(grid%ny_phys, wp))
      any_dry = .false.
      do j = jlo, jhi
         do i = ng + 1, ng + grid%nx_phys
            if (b(i, j) <= 1.0_wp) any_dry = .true.
         end do
      end do
      call check(error,.not. any_dry, "nw2 aquaplanet: interior land present (should be all wet)")
      deallocate (b)
   end subroutine test_aquaplanet

   subroutine test_min_depth_floor(error)
      !! The min_depth floor (MOM6 MINIMUM_DEPTH analogue) lifts every cell to
      !! at least min_depth — no land, no thin cells — while leaving deep cells
      !! (b > min_depth) at their unfloored formula value.  This keeps the
      !! C-grid dyn-core stable (true land + sub-metre layers blow up).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), b0(:, :)
      real(wp), parameter :: FLOOR = 500.0_wp
      integer :: i, j

      call grid%init(60, 70, NGHOST, 1.0_wp, 1.0_wp)
      allocate (b(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (b0(grid%nx_total, grid%ny_total), source=0.0_wp)
      call set_bathymetry_neverworld2(b0, grid, MAX_DEPTH, 1.0_wp, 0.05_wp, 0.0_wp)   ! unfloored
      call set_bathymetry_neverworld2(b, grid, MAX_DEPTH, 1.0_wp, 0.05_wp, FLOOR)     ! floored

      call check(error, minval(b) >= FLOOR - 1.0e-9_wp, "nw2 floor: a cell below min_depth")
      if (allocated(error)) return
      ! Where the unfloored depth already exceeds the floor, the value is
      ! unchanged; elsewhere it equals the floor.  b = max(b0, FLOOR).
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            if (abs(b(i, j) - max(b0(i, j), FLOOR)) > 1.0e-9_wp) then
               call check(error, .false., "nw2 floor: b /= max(unfloored, min_depth)")
               deallocate (b, b0)
               return
            end if
         end do
      end do
      call check(error, maxval(b) > 0.9_wp*MAX_DEPTH, "nw2 floor: deep interior lost")
      deallocate (b, b0)
   end subroutine test_min_depth_floor

   subroutine test_wind_checkvalue(error)
      !! τ_x matches the reference 3-band profile cell-by-cell; τ_y ≡ 0;
      !! ghost rows stay zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAUX = 0.2_wp
      real(wp) :: y, expect, maxerr
      integer :: i, j, ng, j_phys

      call grid%init(40, 70, NGHOST, 1.0_wp, 1.0_wp)
      call ss%init(grid, 1)
      call ss%set_wind_stress_neverworld2(grid, TAUX)
      ng = grid%nghost

      call check(error, all(abs(ss%tau_y) < 1.0e-15_wp), "nw2 wind: tau_y not zero")
      if (allocated(error)) return

      maxerr = 0.0_wp
      do j = ng + 1, ng + grid%ny_phys
         j_phys = j - ng
         y = (real(j_phys, wp) - 0.5_wp)/real(grid%ny_phys, wp)
         expect = ref_taux(y, TAUX)
         do i = 1, size(ss%tau_x, 1)
            maxerr = max(maxerr, abs(ss%tau_x(i, j) - expect))
         end do
      end do
      call check(error, maxerr < 1.0e-12_wp, "nw2 wind: tau_x deviates from reference profile")
      if (allocated(error)) return

      ! Ghost rows (south of ng+1, north of ng+ny_phys) stay zero.
      call check(error, all(abs(ss%tau_x(:, 1:ng)) < 1.0e-15_wp), &
                 "nw2 wind: south ghost rows not zero")
      if (allocated(error)) return
      call check(error, all(abs(ss%tau_x(:, ng + grid%ny_phys + 1:)) < 1.0e-15_wp), &
                 "nw2 wind: north ghost rows not zero")
      call ss%destroy()
   end subroutine test_wind_checkvalue

   subroutine test_wind_structure(error)
      !! Physical structure: a strong westerly (τ_x > 0) peak in the southern
      !! band, and τ_x → 0 at both polar edges.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAUX = 0.2_wp
      integer :: ng

      call grid%init(40, 70, NGHOST, 1.0_wp, 1.0_wp)
      call ss%init(grid, 1)
      call ss%set_wind_stress_neverworld2(grid, TAUX)
      ng = grid%nghost

      ! Peak westerly is a substantial fraction of TAUX (the band-2 plateau and
      ! band-1 ramp both peak near TAUX·1.0 around y≈0.25–0.29).
      call check(error, maxval(ss%tau_x) > 0.5_wp*TAUX, "nw2 wind: no strong westerly")
      if (allocated(error)) return
      ! Southernmost physical row (y≈0.007): essentially zero stress.
      call check(error, maxval(abs(ss%tau_x(:, ng + 1))) < 0.05_wp*TAUX, &
                 "nw2 wind: southern polar edge not near-zero")
      if (allocated(error)) return
      ! Northernmost physical row (y≈0.993 > 1-off): zero (polar cutoff).
      call check(error, maxval(abs(ss%tau_x(:, ng + grid%ny_phys))) < 1.0e-12_wp, &
                 "nw2 wind: northern polar cutoff not zero")
      call ss%destroy()
   end subroutine test_wind_structure

end module test_ocean_neverworld2
