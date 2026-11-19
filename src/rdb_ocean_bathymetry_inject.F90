!! Pre-create bathymetry ARRAY injection (Python runtime API plan, P2.5).
module rdb_ocean_bathymetry_inject
   !! Oceananigans-style geometry construction (`docs/ocean_python_api_plan.md`
   !! S5b, `06_python_surface_design.md` D6.2) hands bathymetry to the
   !! library as an interior-sized array instead of a `topo_config` formula
   !! or a NetCDF file. Two pieces live here, deliberately in a NetCDF-FREE
   !! module (unlike `rdb_bathymetry`, which requires
   !! `RDB_ENABLE_NETCDF=ON` purely because it also happens to contain the
   !! NetCDF reader): sign normalisation/validation, and ghost fill. Array
   !! injection is explicitly meant to need no filesystem/NetCDF at all — a
   !! caller building a grid in memory and handing it straight to `create()`
   !! must not be forced onto a NetCDF build just to get its ghosts filled.
   !!
   !! **Sign is the single most dangerous argument in this API (D6.2).**
   !! Roundabout's `%barotropic%b` is a POSITIVE-DOWN depth (a 4000 m-deep cell
   !! is `+4000`); GEBCO/ETOPO/Oceananigans `GridFittedBottom` ship a
   !! NEGATIVE-DOWN height (the same cell is `-4000`). Passing a GEBCO array
   !! straight through un-flipped makes every cell read as land (`b < 2.0 =
   !! LAND_DEPTH_THRESHOLD`), producing a clean, crash-free, entirely-wrong
   !! quiescent run — and a quiescent uniform-rho zero-velocity run is this
   !! project's OWN correctness check (memory `quiescent-ic`), so the sign
   !! error would look like a passing validation. The guard:
   !!   1. `convention` is REQUIRED — no default (`bathymetry_normalise_sign`
   !!      has no default branch; an unrecognised value is
   !!      `OCEAN_STATUS_ERR_SETUP`, not a silent fall-through).
   !!   2. after normalising to positive-down, the WET FRACTION is checked
   !!      (`b > LAND_DEPTH_THRESHOLD`) — zero wet cells is
   !!      `OCEAN_STATUS_ERR_BATHYMETRY_SIGN`, naming the array's median and
   !!      the convention that was asked for.
   !!   3. this is NOT a "no negatives" test: `b < 0` is legal under
   !!      wet/dry (`rdb_config.F90` `&ocean_wetdry_nml`), so the check is on
   !!      the wet fraction, never on the presence of a negative value.
   !!
   !! Ghost fill (`bathymetry_fill_ghosts_array`) mirrors
   !! `rdb_bathymetry::fill_bathymetry_ghosts_array` (constant extrapolation
   !! from the nearest interior cell) — CLAUDE.md's formula-bathymetry
   !! ghost-fill gotcha: an unfilled ghost row leaves `b=0` there, the EOS
   !! falls back to `rho=rho_0`, and the spurious density jump at the
   !! wall-adjacent face e-folds the basin in ~12 h. Deliberately
   !! duplicated (not `use`d from `rdb_bathymetry`) rather than pulled in
   !! through a NetCDF-gated module — see the module docstring above.
   use rdb_constants, only: wp, LAND_DEPTH_THRESHOLD
   use rdb_grid, only: hgrid_t
   use rdb_error_ring, only: fail
   use pic_strings, only: to_string
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP, &
                               OCEAN_STATUS_ERR_BATHYMETRY_SIGN
   implicit none
   private

   public :: bathymetry_normalise_sign
   public :: bathymetry_fill_ghosts_array
   public :: BATHY_CONVENTION_DEPTH_POSITIVE_DOWN
   public :: BATHY_CONVENTION_HEIGHT_POSITIVE_UP

   integer, parameter :: BATHY_CONVENTION_DEPTH_POSITIVE_DOWN = 1
      !! Caller's array already matches Roundabout's internal convention: a
      !! 4000 m-deep cell is `+4000`. No sign flip.
   integer, parameter :: BATHY_CONVENTION_HEIGHT_POSITIVE_UP = 2
      !! Caller's array is a GEBCO/ETOPO/Oceananigans-style height: a
      !! 4000 m-deep cell is `-4000`. Flipped (negated) to positive-down.

contains

   subroutine bathymetry_normalise_sign(b, convention, ierr)
      !! In-place sign-normalise `b` (any shape — interior or full array,
      !! the caller decides what it passes) to Roundabout's positive-down depth
      !! convention, then validate on the NORMALISED array: zero wet cells
      !! (`b > LAND_DEPTH_THRESHOLD`) is rejected as
      !! `OCEAN_STATUS_ERR_BATHYMETRY_SIGN`, naming the median depth and the
      !! convention that was requested, so the caller can see at a glance
      !! that the OTHER convention was probably meant. NOT a "no negatives"
      !! check (see module docstring) — a majority-negative but non-empty
      !! wet fraction is accepted (legal under wet/dry).
      real(wp), intent(inout) :: b(:, :)
      integer, intent(in) :: convention
         !! One of `BATHY_CONVENTION_DEPTH_POSITIVE_DOWN` /
         !! `_HEIGHT_POSITIVE_UP`. No default — any other value is
         !! `OCEAN_STATUS_ERR_SETUP`.
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP` on an unrecognised
         !! `convention`, `OCEAN_STATUS_ERR_BATHYMETRY_SIGN` on a
         !! zero-wet-cell array) when present; absent behaves as today
         !! (`error stop`).

      integer :: n_wet, n_total
      real(wp) :: wet_frac, median_b

      select case (convention)
      case (BATHY_CONVENTION_DEPTH_POSITIVE_DOWN)
         ! Already Roundabout's internal convention — no flip.
      case (BATHY_CONVENTION_HEIGHT_POSITIVE_UP)
         b = -b
      case default
         call fail("bathymetry_normalise_sign: unrecognised convention = "// &
                   to_string(convention)//" (must be BATHY_CONVENTION_DEPTH_POSITIVE_DOWN "// &
                   "= 1 or BATHY_CONVENTION_HEIGHT_POSITIVE_UP = 2 — this argument has no "// &
                   "default; see D6.2 in the Python API design notes)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end select

      n_total = size(b)
      n_wet = count(b > LAND_DEPTH_THRESHOLD)
      wet_frac = 0.0_wp
      if (n_total > 0) wet_frac = real(n_wet, wp)/real(n_total, wp)

      if (n_wet == 0) then
         median_b = bathymetry_median(b)
         call fail("bathymetry_normalise_sign: zero wet cells after sign "// &
                   "normalisation (convention = "//to_string(convention)// &
                   ", median depth = "//to_string(median_b)// &
                   " m) — every cell reads as land, which produces a clean, "// &
                   "crash-free, entirely wrong quiescent run. Check the convention "// &
                   "argument; the OTHER convention is probably intended.", &
                   ierr, OCEAN_STATUS_ERR_BATHYMETRY_SIGN)
         return
      end if

      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine bathymetry_normalise_sign

   pure function bathymetry_median(b) result(med)
      !! Approximate median via a full sort — `b` is a setup-time array
      !! (called once per `create()`, never per-step), so O(n log n) is
      !! fine; this exists purely to make the sign-error message name a
      !! representative depth rather than `minval`/`maxval` (which a single
      !! outlier cell would distort).
      real(wp), intent(in) :: b(:, :)
      real(wp) :: med
      real(wp), allocatable :: flat(:)
      integer :: n

      n = size(b)
      allocate (flat(n), source=reshape(b, [n]))
      call sort_real(flat)
      if (mod(n, 2) == 1) then
         med = flat((n + 1)/2)
      else
         med = 0.5_wp*(flat(n/2) + flat(n/2 + 1))
      end if
   end function bathymetry_median

   pure subroutine sort_real(a)
      !! Plain insertion sort. `a` is at most a few hundred cells in every
      !! realistic setup-time call (and correctness, not speed, is what
      !! matters for a diagnostic median) — no need for anything fancier.
      real(wp), intent(inout) :: a(:)
      integer :: i, j
      real(wp) :: key

      do i = 2, size(a)
         key = a(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= key) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = key
      end do
   end subroutine sort_real

   subroutine bathymetry_fill_ghosts_array(b, grid)
      !! Fill ghost-cell bathymetry by constant extrapolation from the
      !! nearest interior cell. Deliberate duplicate of
      !! `rdb_bathymetry::fill_bathymetry_ghosts_array` — see the module
      !! docstring for why this module cannot `use` that one.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid

      integer :: ng, i, j, nx, ny

      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys

      ! West and east ghost columns
      do j = 1, grid%ny_total
         do i = 1, ng
            b(i, j) = b(ng + 1, j)                  ! west
            b(ng + nx + i, j) = b(ng + nx, j)       ! east
         end do
      end do

      ! South and north ghost rows (corners already filled above)
      do j = 1, ng
         do i = 1, grid%nx_total
            b(i, j) = b(i, ng + 1)                  ! south
            b(i, ng + ny + j) = b(i, ng + ny)       ! north
         end do
      end do
   end subroutine bathymetry_fill_ghosts_array

end module rdb_ocean_bathymetry_inject
