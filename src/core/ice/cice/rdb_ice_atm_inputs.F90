!! Atmosphere inputs of the CICE-port ice thermodynamics and dynamics: the per-cell
!! forcing the vertical thermodynamics reads, the T-averaged ice velocity of the
!! high-frequency wind stress, and the T-point wind stress handed to the C-grid
!! dynamics at its faces.  Roundabout-original code (no upstream ancestry): the
!! formulas are the standard ones named below.
module rdb_ice_atm_inputs
   !! **What it builds** (`ice_atm_prepare_cell_inputs`, ONE `do concurrent (j, i)`
   !! launch over the owned cells, plain explicit-shape arrays).  From the nine-field
   !! prescribed-atmosphere state of `rdb_atm_state` (`tas` [K], `huss` [kg/kg], `uas, vas`
   !! [m/s, EARTH-relative], `psl` [Pa], `rlds` [W/m2], `prra`, `prsn` [kg/m2/s]; the
   !! caller passes the arrays, this module does not `use` the ocean) the inputs of
   !! `ice_step_therm1`:
   !!
   !!  * `potT = tas*(1e5/psl)^0.286` [K], the Poisson relation with `kappa = R/cp = 0.286`
   !!    and the 1000 hPa reference, taking `psl` as the pressure;
   !!  * `rhoa = psl/(287.058*(1 + 0.608*huss)*tas)` [kg/m3], the ideal-gas law with the
   !!    virtual-temperature factor (Large and Yeager 2009, the same law the ocean bulk
   !!    formulae use, but with the sea-ice constants set `R = 287.058`, `0.608`, which
   !!    differ from the ocean's `287.04`, `0.6078`: plan decision M-D8-b, "CICE's
   !!    potT/rhoa formulas for the ice");
   !!  * `uatm, vatm` = the earth-relative wind rotated onto the model axes with the
   !!    grid angle `angle_dx` exactly as the ocean bulk forcing does
   !!    (`u = cos(a) uas + sin(a) vas`, `v = -sin(a) uas + cos(a) vas`), `wind =
   !!    sqrt(uatm^2 + vatm^2)` (rotation-invariant), `Qa = huss`, `flw = rlds`, `fsnow =
   !!    prsn`, `frain = prra` (the forcing gives the rain/snow split; nothing is split
   !!    by air temperature).
   !!
   !! **Ocean and ice see slightly different air, on purpose.**  The ocean bulk forcing
   !! uses the air temperature `tas` directly with the ocean's own constants
   !! (`R = 287.04`, `0.6078`) and wind floor; the ice uses the potential temperature
   !! and the sea-ice constants set above (plan decision M-D8-b).  The pressure used
   !! here is the sea-level pressure `psl` of the forcing, where the reference uses
   !! the lowest-level pressure (equal only if the forcing provides the same field):
   !! an interim choice, logged in the PR report.
   !!
   !! **The grid angle** is static: `ice_atm_grid_rotation` computes `cos` and `sin`
   !! of it ONCE at setup into two arrays the caller keeps device-resident and passes
   !! to the builder (no transcendental per cell per step).
   !!
   !! Reference heights: `zlvl = z_wind` (the height of the wind and state data), `zlvs`
   !! a fixed 10 m (`ice_atm_heights`; interim: the data atmosphere never sends a scalar
   !! height).
   !!
   !! **Units, handled ONCE.**  The state is in kelvin; a source in degC passes
   !! `tair_offset = 273.15` (added to `tas` before everything else).  Nothing else in
   !! this module knows about degC.
   !!
   !! **Land and bad cells.**  A cell with `hm <= 0.5` is land: every output is the benign
   !! defined value (`ICE_ATM_BENIGN_*`: a positive temperature and density, zero wind,
   !! humidity, longwave and precipitation) and the inputs are never read (a NaN over
   !! land is harmless).  A wet cell whose inputs are non-finite or whose `tas` or
   !! `psl` is not positive (after the offset) gets the same benign outputs AND
   !! `ICE_ERR_NONFINITE` OR-ed into `ice_err(i, j)`: never a silent default.
   !! Only the owned range `is:ie, js:je` is written; the ring is untouched.
   !!
   !! **High-frequency velocity** (`ice_atm_velocity_T`).  `uvelT, vvelT` are the
   !! start-of-step ice velocity at the T points for the relative-wind stress of
   !! `rdb_ice_atmo` with `highfreq` (Roberts et al. 2014), in the form the C-grid
   !! reference uses: the face velocities are averaged to the four CORNERS of the cell
   !! with the dynamics' own corner average (`evp_wavg2` with the face areas, times the
   !! corner mask `uvm`: the expressions of `evp_average_impl`), and the four corners
   !! are averaged to T with the corner areas `uarea` as UNMASKED weights.  On a
   !! uniform grid that is a six-face stencil: the faces of the west and east columns
   !! with weights `(1/8)(u(j-1) + 2 u(j) + u(j+1))`.  Next to land the masked corners
   !! carry zero velocity but keep their weight, so the cell mean shrinks toward the
   !! coast exactly as upstream's does.  Closed faces must carry 0 (the dynamics'
   !! invariant; the face is never masked here).  ADDITION to the reference: a cell
   !! that is land, or whose start-of-step concentration is not above `amin`, gets
   !! exactly 0 (the thermodynamics reads the velocity only under ice).  A non-finite
   !! result flags the cell.  Needs one west and one south neighbour: `is, js >= 2`.
   !! Only called when `highfreq`; with it off `ice_atm_boundary_ice` ignores the ice
   !! velocity bitwise.
   !!
   !! **Wind stress to the faces -- the ONE place** (`ice_atm_stress_to_faces`).  The EVP
   !! preparation (`ice_evp_c_prep`) takes the wind stress ON THE FACES (`strairxE`,
   !! `strairyN`) and does not average anything itself, so this is the single T-to-face
   !! average: the flux map `S_E(i,j) = (S_T(i,j) A_T(i,j) + S_T(i+1,j) A_T(i+1,j)) /
   !! (2 A_E(i,j))` (and the y analogue to N), unmasked, area-weighted so that the
   !! force integrated over the two cells is the force integrated over the face
   !! cell; zero where the face area is zero (those faces drop their force; the identity
   !! `sum_E A_E S_E = sum_i A_i S_i` holds over the faces with positive area).  A land
   !! neighbour has zero stress and contributes zero: the face next to land gets half of
   !! the cell's stress scaled by `A_T/A_E`, with no renormalisation.  `strair_x_T, strair_y_T` are the
   !! merged T-point stress of `ice_step_therm1` (zero where there is no ice) and must carry a valid
   !! east and north neighbour ring (halo exchange before the call: the neighbour
   !! of an owned face is a ghost on a tile edge).
   !!
   !! Every kernel is a `do concurrent` over explicit-shape arrays; transcendentals
   !! are the reproducible `rdb_sm_*` bodies.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_evp_c, only: evp_wavg2
   implicit none
   private

   public :: ice_atm_prepare_cell_inputs, ice_atm_cell_inputs, ice_atm_heights
   public :: ice_atm_velocity_T, ice_atm_stress_to_faces, ice_atm_grid_rotation
   public :: ICE_ATM_PREF, ICE_ATM_KAPPA_AIR, ICE_ATM_RGAS, ICE_ATM_TVQ, ICE_ATM_ZLVS
   public :: ICE_ATM_BENIGN_POTT, ICE_ATM_BENIGN_RHOA

   real(wp), parameter :: ICE_ATM_PREF = 1.0e5_wp
      !! Reference pressure of the potential temperature [Pa].
   real(wp), parameter :: ICE_ATM_KAPPA_AIR = 0.286_wp
      !! `R/cp` of dry air in the potential-temperature exponent.
   real(wp), parameter :: ICE_ATM_RGAS = 287.058_wp
      !! Gas constant of dry air of the air density [J/kg/K].
   real(wp), parameter :: ICE_ATM_TVQ = 0.608_wp
      !! Virtual-temperature coefficient of the humidity.
   real(wp), parameter :: ICE_ATM_ZLVS = 10.0_wp
      !! Scalar reference height [m] (interim, fixed).
   real(wp), parameter :: ICE_ATM_BENIGN_POTT = 273.15_wp
      !! `potT` of a land or bad cell [K].
   real(wp), parameter :: ICE_ATM_BENIGN_RHOA = 1.0_wp
      !! `rhoa` of a land or bad cell [kg/m3].

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_pow.inc"
#include "rdb_safe_math_trig.inc"

   pure subroutine ice_atm_heights(z_wind, zlvl, zlvs, ok)
      !! The reference heights of the thermodynamics: `zlvl = z_wind`, `zlvs = 10`.
      !! `ok` is false for a non-finite or non-positive `z_wind` (the heights are then
      !! returned as zero).  Host setup.
      real(wp), intent(in) :: z_wind
      real(wp), intent(out) :: zlvl, zlvs
      logical, intent(out) :: ok

      ok = ieee_is_finite(z_wind)
      if (ok) ok = z_wind > 0.0_wp
      zlvl = 0.0_wp
      zlvs = 0.0_wp
      if (ok) then
         zlvl = z_wind
         zlvs = ICE_ATM_ZLVS
      end if
   end subroutine ice_atm_heights

   pure subroutine ice_atm_cell_inputs(tas, huss, uas, vas, psl, rlds, prra, prsn, ca, sa, tair_offset, &
                                       potT, rhoa, wind, uatm, vatm, Qa, flw, fsnow, frain, ok)
      !! One wet cell: the formulas of the module header on scalars.  `ok` is false on a
      !! non-finite input or a non-positive temperature/pressure, with every output the
      !! benign value.
      !$acc routine seq
      real(wp), intent(in) :: tas, huss, uas, vas, psl, rlds, prra, prsn, ca, sa, tair_offset
      real(wp), intent(out) :: potT, rhoa, wind, uatm, vatm, Qa, flw, fsnow, frain
      logical, intent(out) :: ok
      real(wp) :: tk

      potT = ICE_ATM_BENIGN_POTT
      rhoa = ICE_ATM_BENIGN_RHOA
      wind = 0.0_wp
      uatm = 0.0_wp
      vatm = 0.0_wp
      Qa = 0.0_wp
      flw = 0.0_wp
      fsnow = 0.0_wp
      frain = 0.0_wp
      ok = ieee_is_finite(tas) .and. ieee_is_finite(huss) .and. ieee_is_finite(uas) .and. ieee_is_finite(vas) &
           .and. ieee_is_finite(psl) .and. ieee_is_finite(rlds) .and. ieee_is_finite(prra) &
           .and. ieee_is_finite(prsn) .and. ieee_is_finite(ca) .and. ieee_is_finite(sa) &
           .and. ieee_is_finite(tair_offset)
      if (.not. ok) return
      tk = tas + tair_offset
      ok = tk > 0.0_wp .and. psl > 0.0_wp
      if (.not. ok) return
      uatm = ca*uas + sa*vas
      vatm = -sa*uas + ca*vas
      wind = sqrt(uatm*uatm + vatm*vatm)
      potT = tk*rdb_sm_pow(ICE_ATM_PREF/psl, ICE_ATM_KAPPA_AIR)
      rhoa = psl/(ICE_ATM_RGAS*(1.0_wp + ICE_ATM_TVQ*huss)*tk)
      Qa = huss
      flw = rlds
      fsnow = prsn
      frain = prra
      ok = ieee_is_finite(potT) .and. ieee_is_finite(rhoa) .and. rhoa > 0.0_wp
      if (.not. ok) then
         potT = ICE_ATM_BENIGN_POTT
         rhoa = ICE_ATM_BENIGN_RHOA
         wind = 0.0_wp
         uatm = 0.0_wp
         vatm = 0.0_wp
         Qa = 0.0_wp
         flw = 0.0_wp
         fsnow = 0.0_wp
         frain = 0.0_wp
      end if
   end subroutine ice_atm_cell_inputs

   subroutine ice_atm_grid_rotation(nx, ny, is, ie, js, je, angle, cosA, sinA)
      !! Setup: `cosA = cos(angle)`, `sinA = sin(angle)` on the owned cells (the reproducible
      !! `rdb_sm_*` bodies).  The caller keeps the two arrays device-resident.  A non-finite
      !! angle gives NaN, which the builder flags.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: angle(nx, ny)
      real(wp), intent(inout) :: cosA(nx, ny), sinA(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         cosA(i, j) = rdb_sm_cos(angle(i, j))
         sinA(i, j) = rdb_sm_sin(angle(i, j))
      end do
   end subroutine ice_atm_grid_rotation

   subroutine ice_atm_prepare_cell_inputs(nx, ny, is, ie, js, je, tair_offset, hm, cosA, sinA, &
                                          tas, huss, uas, vas, psl, rlds, prra, prsn, &
                                          potT, rhoa, wind, uatm, vatm, Qa, flw, fsnow, frain, ice_err)
      !! The per-cell atmosphere inputs of `ice_step_therm1` over the owned cells, one
      !! launch (module header).  `hm` is the wet-T mask (real 0/1), `angle` the grid
      !! angle [rad] of the rotation to the model axes (0 on a Cartesian grid), the
      !! eight forcing arrays are `atm_state_t`'s `tas ... prsn`.  Outputs on the
      !! owned cells; `ice_err` is OR-ed into.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: tair_offset
      real(wp), intent(in) :: hm(nx, ny), cosA(nx, ny), sinA(nx, ny)
      real(wp), intent(in) :: tas(nx, ny), huss(nx, ny), uas(nx, ny), vas(nx, ny), psl(nx, ny)
      real(wp), intent(in) :: rlds(nx, ny), prra(nx, ny), prsn(nx, ny)
      real(wp), intent(inout) :: potT(nx, ny), rhoa(nx, ny), wind(nx, ny), uatm(nx, ny), vatm(nx, ny)
      real(wp), intent(inout) :: Qa(nx, ny), flw(nx, ny), fsnow(nx, ny), frain(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j
      logical :: ok

      do concurrent(j=js:je, i=is:ie) local(ok)
         if (hm(i, j) > 0.5_wp) then
            call ice_atm_cell_inputs(tas(i, j), huss(i, j), uas(i, j), vas(i, j), psl(i, j), rlds(i, j), &
                                     prra(i, j), prsn(i, j), cosA(i, j), sinA(i, j), tair_offset, &
                                     potT(i, j), rhoa(i, j), wind(i, j), uatm(i, j), vatm(i, j), Qa(i, j), &
                                     flw(i, j), fsnow(i, j), frain(i, j), ok)
            if (.not. ok) ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
         else
            potT(i, j) = ICE_ATM_BENIGN_POTT
            rhoa(i, j) = ICE_ATM_BENIGN_RHOA
            wind(i, j) = 0.0_wp
            uatm(i, j) = 0.0_wp
            vatm(i, j) = 0.0_wp
            Qa(i, j) = 0.0_wp
            flw(i, j) = 0.0_wp
            fsnow(i, j) = 0.0_wp
            frain(i, j) = 0.0_wp
         end if
      end do
   end subroutine ice_atm_prepare_cell_inputs

   subroutine ice_atm_velocity_T(nx, ny, is, ie, js, je, amin, u, v, earea, narea, uarea, uvm, hm, aice_init, &
                                 uvelT, vvelT, ice_err)
      !! `uvelT, vvelT` at the owned T cells from the face velocities `u(nx+1, ny)`
      !! (east face of cell `(i, j)` at `u(i+1, j)`) and `v(nx, ny+1)` (north face at
      !! `v(i, j+1)`), one launch (module header): faces to the four corners with
      !! `evp_wavg2`, corners to T with the corner areas.  `earea, narea` are the face
      !! areas at the east and north faces of cell `(i, j)`, `uarea(i, j)` the area of
      !! its NE corner, `uvm` the corner mask (real 0/1), `hm` the wet-T mask,
      !! `aice_init` the start-of-step concentration, `amin` the ice threshold.
      !! Requires `is >= 2` and `js >= 2`.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: amin
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), uarea(nx, ny), uvm(nx, ny), hm(nx, ny), aice_init(nx, ny)
      real(wp), intent(inout) :: uvelT(nx, ny), vvelT(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j
      real(wp) :: wt, su, sv
      logical :: gate

      do concurrent(j=js:je, i=is:ie) local(wt, su, sv, gate)
         su = 0.0_wp
         sv = 0.0_wp
         gate = hm(i, j) > 0.5_wp
         if (gate) gate = aice_init(i, j) > amin
         if (gate) then
            wt = uarea(i, j) + uarea(i - 1, j) + uarea(i, j - 1) + uarea(i - 1, j - 1)
            if (wt /= 0.0_wp) then
               ! corners NE (i,j), NW (i-1,j), SE (i,j-1), SW (i-1,j-1): u from the east column of the
               ! corner's two faces, v from its two north faces (the expressions of `evp_average_impl`)
               su = (evp_wavg2(u(i + 1, j), u(i + 1, j + 1), earea(i, j), earea(i, j + 1))*uvm(i, j)*uarea(i, j) &
                     + evp_wavg2(u(i, j), u(i, j + 1), earea(i - 1, j), earea(i - 1, j + 1))*uvm(i - 1, j)*uarea(i - 1, j) &
                     + evp_wavg2(u(i + 1, j - 1), u(i + 1, j), earea(i, j - 1), earea(i, j))*uvm(i, j - 1)*uarea(i, j - 1) &
                     + evp_wavg2(u(i, j - 1), u(i, j), earea(i - 1, j - 1), earea(i - 1, j)) &
                     *uvm(i - 1, j - 1)*uarea(i - 1, j - 1))/wt
               sv = (evp_wavg2(v(i, j + 1), v(i + 1, j + 1), narea(i, j), narea(i + 1, j))*uvm(i, j)*uarea(i, j) &
                     + evp_wavg2(v(i - 1, j + 1), v(i, j + 1), narea(i - 1, j), narea(i, j))*uvm(i - 1, j)*uarea(i - 1, j) &
                     + evp_wavg2(v(i, j), v(i + 1, j), narea(i, j - 1), narea(i + 1, j - 1))*uvm(i, j - 1)*uarea(i, j - 1) &
                     + evp_wavg2(v(i - 1, j), v(i, j), narea(i - 1, j - 1), narea(i, j - 1)) &
                     *uvm(i - 1, j - 1)*uarea(i - 1, j - 1))/wt
            end if
            if (.not. (ieee_is_finite(su) .and. ieee_is_finite(sv))) then
               su = 0.0_wp
               sv = 0.0_wp
               ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
            end if
         end if
         uvelT(i, j) = su
         vvelT(i, j) = sv
      end do
   end subroutine ice_atm_velocity_T

   subroutine ice_atm_stress_to_faces(nx, ny, is, ie, js, je, strair_x_T, strair_y_T, tarea, earea, narea, &
                                      strairxE, strairyN)
      !! The wind stress on the ice at the east faces `strairxE(i, j)` (between T
      !! `(i, j)` and `(i+1, j)`) and north faces `strairyN(i, j)` of the owned cells:
      !! the area-weighted two-point flux map of the module header.  Needs the T stress
      !! at `ie + 1` and `je + 1` (halo).  One launch.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: strair_x_T(nx, ny), strair_y_T(nx, ny), tarea(nx, ny), earea(nx, ny), narea(nx, ny)
      real(wp), intent(inout) :: strairxE(nx, ny), strairyN(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         if (earea(i, j) > 0.0_wp) then
            strairxE(i, j) = 0.5_wp*(strair_x_T(i, j)*tarea(i, j) + strair_x_T(i + 1, j)*tarea(i + 1, j))/earea(i, j)
         else
            strairxE(i, j) = 0.0_wp
         end if
         if (narea(i, j) > 0.0_wp) then
            strairyN(i, j) = 0.5_wp*(strair_y_T(i, j)*tarea(i, j) + strair_y_T(i, j + 1)*tarea(i, j + 1))/narea(i, j)
         else
            strairyN(i, j) = 0.0_wp
         end if
      end do
   end subroutine ice_atm_stress_to_faces

end module rdb_ice_atm_inputs
