!! SIS2 -> Icepack restart conversion (D1.5): the array converter and the start-up read
!! path.  Closed-form expectations only (densities, the Bitz-Lipscomb enthalpy relation,
!! the SIS2 enthalpy relation, overlap weights); nothing is compared against the
!! converter's own output except the read-path test, which checks the file round trip
!! against a direct converter call on the same arrays.
!!
!! Tolerances are ~10x the MEASURED worst on gfortran (CPU): mass 1.7e-16 relative
!! (one division and one multiplication) -> 2e-15; closed-form heat sums 1.1e-15 -> 1e-14;
!! temperature round trip 5.2e-14 K (quadratic solve near the melting point) -> 5e-13;
!! layer-mean remap 2e-16 -> 2e-15.  The GPU build (FMA on in host code) must meet the same.  Everything runs on the host: the converter is set-up code.
module test_cice_restart_convert
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, &
                                   ocean_restart_metadata_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ice_q_from_T, &
                             ice_T_from_q, ice_Tsno_from_q, ice_salinity_profile, ice_qsno_init_from_T
   use rdb_ice_enthalpy, only: ice_enth_from_ts
   use rdb_ice_itd_params, only: ice_itd_bounds, ITD_BOUNDS_OK
   use rdb_ice_seed, only: ice_seed_uniform_impl, ICE_SEED_OK
   use rdb_ice_restart_convert, only: ice_restart_convert_sis2, ice_restart_remap_layers, &
                                      ice_convert_report_t, ICE_CONV_OK, ICE_CONV_NONFINITE, &
                                      ICE_CONV_NEGATIVE, ICE_CONV_BAD_EXTENT
   use rdb_ice_cice_state, only: ice_cice_state_t, ICE_STATE_VERSION
   use rdb_ice_cice_restart, only: ice_legacy_t, ice_cice_restart_resume, ICE_RESUME_V2, &
                                   ICE_RESUME_CONVERTED, ICE_RESUME_COLD, ICE_RESUME_FAILED
   implicit none
   private

   public :: collect_cice_restart_convert_tests

   integer, parameter :: NX = 6, NY = 4, NCAT = 5, NILYR = 7
   real(wp), parameter :: RHOI_OLD = 905.0_wp, RHOS = 330.0_wp, TF = -1.8_wp, SBULK = 4.0_wp
   real(wp), parameter :: MASS_TOL = 2.0e-15_wp
   real(wp), parameter :: HEAT_TOL = 1.0e-14_wp
   real(wp), parameter :: TEMP_TOL = 5.0e-13_wp
   real(wp), parameter :: REMAP_TOL = 2.0e-15_wp

   type :: sis2_t
      integer :: nc = 0, nk = 0
      real(wp), allocatable :: ps(:, :, :), mi(:, :, :), ms(:, :, :)
      real(wp), allocatable :: en(:, :, :, :), sa(:, :, :, :), es(:, :, :, :)
   end type sis2_t

   type :: newstate_t
      real(wp), allocatable :: a(:, :, :), v(:, :, :), vs(:, :, :), ts(:, :, :)
      real(wp), allocatable :: q(:, :, :, :), s(:, :, :, :), qs(:, :, :, :)
   end type newstate_t

contains

   subroutine collect_cice_restart_convert_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mass_exact_per_category", test_mass_exact), &
                  new_unittest("temperature_round_trip", test_t_round_trip), &
                  new_unittest("layer_reversal_gradient", test_reversal), &
                  new_unittest("resample_conserves_energy", test_resample), &
                  new_unittest("resample_through_converter", test_resample_converter), &
                  new_unittest("category_rebinning_merge", test_rebin), &
                  new_unittest("empty_land_tiny_cells", test_empty), &
                  new_unittest("lumped_single_category", test_lumped), &
                  new_unittest("area_above_one_rescaled", test_rescale), &
                  new_unittest("melting_point_cap_reported", test_cap), &
                  new_unittest("refuses_bad_input", test_refusals), &
                  new_unittest("lumped_thin_film_partial_cover", test_thin_film), &
                  new_unittest("dropped_mass_identity", test_dropped), &
                  new_unittest("tsfcn_under_snow", test_tsfc_snow), &
                  new_unittest("read_path_three_starts", test_read_path), &
                  new_unittest("read_path_refusals", test_read_refusals) &
                  ]
   end subroutine collect_cice_restart_convert_tests

   ! ------------------------------------------------------------------ helpers
   function params() result(p)
      type(ice_params_t) :: p
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, 1)
   end function params

   subroutine bounds(p, hb)
      type(ice_params_t), intent(in) :: p
      real(wp), intent(out) :: hb(0:NCAT)
      integer :: st
      call ice_itd_bounds(p, hb, st)
      if (st /= ITD_BOUNDS_OK) error stop "bounds"
   end subroutine bounds

   subroutine alloc_sis2(s, nc, nk)
      type(sis2_t), intent(out) :: s
      integer, intent(in) :: nc, nk
      s%nc = nc
      s%nk = nk
      allocate (s%ps(NX, NY, 0:nc), source=0.0_wp)
      s%ps(:, :, 0) = 1.0_wp
      allocate (s%mi(NX, NY, nc), s%ms(NX, NY, nc), source=0.0_wp)
      allocate (s%en(NX, NY, nc, nk), source=0.0_wp)
      allocate (s%sa(NX, NY, nc, nk), source=SBULK)
      allocate (s%es(NX, NY, nc, 1), source=0.0_wp)
   end subroutine alloc_sis2

   subroutine alloc_new(o)
      type(newstate_t), intent(out) :: o
      allocate (o%a(NX, NY, NCAT), o%v(NX, NY, NCAT), o%vs(NX, NY, NCAT), o%ts(NX, NY, NCAT), source=7.0_wp)
      allocate (o%q(NX, NY, NCAT, NILYR), o%s(NX, NY, NCAT, NILYR), source=7.0_wp)
      allocate (o%qs(NX, NY, NCAT, 1), source=7.0_wp)
   end subroutine alloc_new

   subroutine set_column(s, i, j, c, a, h, hs, Tb, Tt, Tsn)
      !! SIS2 category `c` of cell (i,j): area `a`, ice thickness `h` [m] at the SIS2
      !! density, snow depth `hs`, temperatures linear in depth from `Tb` (base, k=1) to `Tt`.
      type(sis2_t), intent(inout) :: s
      integer, intent(in) :: i, j, c
      real(wp), intent(in) :: a, h, hs, Tb, Tt, Tsn
      integer :: k
      s%ps(i, j, c) = a
      s%ps(i, j, 0) = max(1.0_wp - sum(s%ps(i, j, 1:s%nc)), 0.0_wp)
      s%mi(i, j, c) = RHOI_OLD*h
      s%ms(i, j, c) = RHOS*hs
      do k = 1, s%nk
         s%en(i, j, c, k) = ice_enth_from_ts(Tb + (Tt - Tb)*(real(k, wp) - 0.5_wp)/real(s%nk, wp), SBULK)
         s%sa(i, j, c, k) = SBULK
      end do
      s%es(i, j, c, 1) = ice_enth_from_ts(Tsn, 0.0_wp)
   end subroutine set_column

   subroutine convert(s, p, hb, o, rep, ierr)
      type(sis2_t), intent(in) :: s
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hb(0:NCAT)
      type(newstate_t), intent(inout) :: o
      type(ice_convert_report_t), intent(out) :: rep
      integer, intent(out) :: ierr
      call ice_restart_convert_sis2(NX, NY, s%nc, s%nk, NCAT, NILYR, 1, p, hb, TF, s%ps, s%mi, s%ms, &
                                    s%en, s%sa, s%es, o%a, o%v, o%vs, o%ts, o%q, o%s, o%qs, rep, ierr)
   end subroutine convert

   subroutine track(worst, x, bad)
      !! Accumulate the worst finite value; a non-finite one sets `bad` (max() would hide a NaN).
      real(wp), intent(inout) :: worst
      real(wp), intent(in) :: x
      logical, intent(inout) :: bad
      if (.not. ieee_is_finite(x)) then
         bad = .true.
      else
         worst = max(worst, x)
      end if
   end subroutine track

   function relerr(got, want) result(e)
      real(wp), intent(in) :: got, want
      real(wp) :: e
      e = abs(got - want)/max(abs(want), tiny(1.0_wp))
   end function relerr

   ! ------------------------------------------------------------------ tests
   subroutine test_mass_exact(error)
      !! Every old category lands in its own new bin (thickness at the bin middle): area is
      !! copied exactly, `rhoi*vicen = a*m_ice` and `rhos*vsnon = a*m_snow` per category
      !! and in total, with the thickness down by 905/917.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), h, worst_i, worst_s, a, u, hs, Tb, Tt, Tsn, Tl, sp(NILYR + 1), hin, hout, mi, ms, e
      logical :: bad
      integer(int64) :: seed
      integer :: i, j, c, ierr, k

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, NCAT, NILYR)
      call alloc_new(o)
      call ice_salinity_profile(NILYR, p%c%saltmax, p%c%min_salin, sp)
      hin = 0.0_wp
      hout = 0.0_wp
      seed = 424242_int64
      do j = 1, NY
         do i = 1, NX
            do c = 1, NCAT
               seed = lcg_next(seed)
               u = lcg_unit(seed)
               h = 0.5_wp*(hb(c - 1) + hb(min(c, NCAT - 1))) + 0.2_wp*(u - 0.5_wp)*(hb(min(c, NCAT - 1)) - hb(c - 1))
               if (c == NCAT) h = 5.0_wp + u
               a = 0.15_wp + 0.04_wp*u
               hs = 0.1_wp*u + 0.02_wp
               Tb = -2.0_wp - 20.0_wp*u
               Tt = -5.0_wp - 15.0_wp*u
               Tsn = -8.0_wp*u - 1.0_wp
               call set_column(s, i, j, c, a, h, hs, Tb, Tt, Tsn)
               ! closed-form heat of this category, both sides, from the relations alone
               mi = RHOI_OLD*h
               ms = RHOS*hs
               do k = 1, NILYR
                  hin = hin + a*mi*ice_enth_from_ts(Tb + (Tt - Tb)*(real(k, wp) - 0.5_wp)/real(NILYR, wp), SBULK) &
                        /real(NILYR, wp)
                  Tl = Tb + (Tt - Tb)*(real(NILYR + 1 - k, wp) - 0.5_wp)/real(NILYR, wp)
                  hout = hout + a*mi/p%c%rhoi*ice_q_from_T(Tl, sp(k), p%c, p%thermo%l_brine)/real(NILYR, wp)
               end do
               hin = hin + a*ms*ice_enth_from_ts(Tsn, 0.0_wp)
               hout = hout + a*ms/p%c%rhos*ice_qsno_init_from_T(Tsn, p%c)
            end do
         end do
      end do
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      worst_i = 0.0_wp
      worst_s = 0.0_wp
      bad = .false.
      do c = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               if (.not. (ieee_is_finite(o%v(i, j, c)) .and. ieee_is_finite(o%vs(i, j, c)))) then
                  call test_failed(error, "non-finite volume")
                  return
               end if
               if (o%a(i, j, c) /= s%ps(i, j, c)) then
                  call test_failed(error, "area copied exactly")
                  return
               end if
               e = relerr(p%c%rhoi*o%v(i, j, c), s%ps(i, j, c)*s%mi(i, j, c))
               call track(worst_i, e, bad)
               e = relerr(p%c%rhos*o%vs(i, j, c), s%ps(i, j, c)*s%ms(i, j, c))
               call track(worst_s, e, bad)
            end do
         end do
      end do
      call check(error,.not. bad, "mass errors finite")
      call check(error, worst_i <= MASS_TOL .and. worst_s <= MASS_TOL, "per-category mass exact to 2e-15")
      call check(error, relerr(rep%mass_ice_out, rep%mass_ice_in) <= 8.0_wp*MASS_TOL, "total ice mass exact (report)")
      call check(error, relerr(rep%mass_snow_out, rep%mass_snow_in) <= 8.0_wp*MASS_TOL, "total snow mass exact (report)")
      call check(error, relerr(o%v(1, 1, 1)/(o%a(1, 1, 1)), s%mi(1, 1, 1)/p%c%rhoi) < 1.0e-14_wp, "thickness = m/rhoi_new")
      call check(error, abs(p%c%rhoi*o%v(1, 1, 1)/o%a(1, 1, 1)/(s%mi(1, 1, 1)) - 1.0_wp) < 1.0e-14_wp, "mass per area")
      call check(error, rep%n_dropped == 0 .and. rep%n_rescaled == 0, "nothing dropped or rescaled")
      call check(error, rep%n_t_capped == 0, "no cap in this data set")
      call check(error, ieee_is_finite(rep%heat_rel_change), "heat change finite")
      if (allocated(error)) return
      call check(error, relerr(rep%heat_in, hin) < HEAT_TOL .and. relerr(rep%heat_out, hout) < HEAT_TOL, &
                 "heat before/after equal the closed-form category sums")
      call check(error, abs(rep%heat_rel_change - (hout - hin)/abs(hin)) < HEAT_TOL, &
                 "relative heat change equals the closed form")
   end subroutine test_mass_exact

   subroutine test_t_round_trip(error)
      !! Same layer count: the temperature recovered from the new `qicen`/`sicen` is the
      !! SIS2 one, layer by layer (top-down), snow too; `Tsfcn` is the top-layer value.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), Tb, Tt, Tsn, Tw, Tgot, worst, e, worst_s
      logical :: bad, nonfin
      integer :: i, j, k, ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, NCAT, NILYR)
      call alloc_new(o)
      do j = 1, NY
         do i = 1, NX
            call set_column(s, i, j, 2, 0.6_wp, 1.0_wp, 0.15_wp, -1.9_wp + 0.01_wp*i, -6.0_wp - 2.0_wp*j - i, &
                            -3.0_wp - 0.5_wp*i)
         end do
      end do
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      worst = 0.0_wp
      worst_s = 0.0_wp
      nonfin = .false.
      do j = 1, NY
         do i = 1, NX
            Tb = -1.9_wp + 0.01_wp*i
            Tt = -6.0_wp - 2.0_wp*j - i
            Tsn = -3.0_wp - 0.5_wp*i
            do k = 1, NILYR
               ! new layer k is the SIS2 layer NILYR+1-k: depth fraction from the base
               Tw = Tb + (Tt - Tb)*(real(NILYR + 1 - k, wp) - 0.5_wp)/real(NILYR, wp)
               call ice_T_from_q(o%q(i, j, 2, k), o%s(i, j, 2, k), p%c, p%thermo%l_brine, Tgot, bad)
               if (bad .or. .not. ieee_is_finite(Tgot)) then
                  call test_failed(error, "T from q failed")
                  return
               end if
               call track(worst, abs(Tgot - Tw), nonfin)
            end do
            call track(worst, abs(ice_Tsno_from_q(o%qs(i, j, 2, 1), p%c) - Tsn), nonfin)
            ! snow present: the surface temperature is the snow temperature
            e = abs(o%ts(i, j, 2) - Tsn)
            call track(worst_s, e, nonfin)
         end do
      end do
      call check(error,.not. nonfin, "round-trip errors finite")
      call check(error, worst < TEMP_TOL .and. worst_s < TEMP_TOL, "layer, snow and surface temperatures round-trip")
   end subroutine test_t_round_trip

   subroutine test_reversal(error)
      !! A strong vertical gradient (base -1.8, top -30 degC) must come out mirrored:
      !! `k = 1` is the COLDEST (top) layer, `k = nilyr` the warmest.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), T1, Tn
      logical :: bad
      integer :: ierr, k
      real(wp) :: Tprev, Tk

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, NCAT, NILYR)
      call alloc_new(o)
      call set_column(s, 3, 2, 3, 0.5_wp, 2.0_wp, 0.1_wp, -1.8_wp, -30.0_wp, -20.0_wp)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call ice_T_from_q(o%q(3, 2, 3, 1), o%s(3, 2, 3, 1), p%c, p%thermo%l_brine, T1, bad)
      call ice_T_from_q(o%q(3, 2, 3, NILYR), o%s(3, 2, 3, NILYR), p%c, p%thermo%l_brine, Tn, bad)
      call check(error, T1 < -25.0_wp .and. Tn > -4.0_wp, "k=1 (top) is the cold end, k=nilyr (base) the warm end")
      Tprev = -100.0_wp
      do k = 1, NILYR
         call ice_T_from_q(o%q(3, 2, 3, k), o%s(3, 2, 3, k), p%c, p%thermo%l_brine, Tk, bad)
         if (Tk <= Tprev) call test_failed(error, "temperature must increase monotonically downward")
         Tprev = Tk
      end do
      call check(error, abs(o%ts(3, 2, 3) + 20.0_wp) < 1.0e-9_wp, "Tsfcn is the snow temperature (snow present)")
   end subroutine test_reversal

   subroutine test_resample(error)
      !! `ice_restart_remap_layers`: hand values 2->4, a hand overlap row 4->7, the layer
      !! mean preserved (to rounding), a uniform column stays uniform, equal counts copy.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: x2(2), y4(4), x4(4), y7(7), u7(7), x7(7), z(7)
      integer :: nki, nko, i
      integer(int64) :: seed
      real(wp) :: worst
      logical :: nonfin

      x2 = [-3.0_wp, 5.0_wp]
      call ice_restart_remap_layers(2, 4, x2, y4)
      call check(error, all(y4 == [-3.0_wp, -3.0_wp, 5.0_wp, 5.0_wp]), "2 -> 4 piecewise constant")
      x4 = [1.0_wp, 2.0_wp, 4.0_wp, 8.0_wp]
      call ice_restart_remap_layers(4, 7, x4, y7)
      call check(error, abs(y7(1) - 1.0_wp) < 1.0e-15_wp .and. abs(y7(2) - (3.0_wp*1.0_wp + 2.0_wp)/4.0_wp) < 1.0e-15_wp, &
                 "4 -> 7 first overlap rows (x1, (3x1+x2)/4)")
      call ice_restart_remap_layers(4, 7, [3.5_wp, 3.5_wp, 3.5_wp, 3.5_wp], u7)
      call check(error, all(abs(u7 - 3.5_wp) < 1.0e-15_wp), "uniform stays uniform")
      x7 = [(0.5_wp*real(i, wp)**2, i=1, 7)]
      call ice_restart_remap_layers(7, 7, x7, z)
      call check(error, all(z == x7), "equal counts: exact copy")
      seed = 99_int64
      worst = 0.0_wp
      nonfin = .false.
      do nki = 1, 8
         do nko = 1, 8
            block
               real(wp) :: xi(nki), xo(nko)
               do i = 1, nki
                  seed = lcg_next(seed)
                  xi(i) = -3.0e8_wp*lcg_unit(seed)
               end do
               call ice_restart_remap_layers(nki, nko, xi, xo)
               call track(worst, abs(sum(xo)/real(nko, wp) - sum(xi)/real(nki, wp))/3.0e8_wp, nonfin)
            end block
         end do
      end do
      call check(error,.not. nonfin, "remap errors finite")
      call check(error, worst < REMAP_TOL, "layer-mean (energy) conserved for every 1..8 -> 1..8 pair")
   end subroutine test_resample

   subroutine test_resample_converter(error)
      !! nk 2 -> 7 and 4 -> 7 through the converter: mass exact, finite, ordered top-down,
      !! heat change of the same order as the equal-layer case (no energy invented).
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s2, s4, s7
      type(newstate_t) :: o2, o4, o7
      type(ice_convert_report_t) :: r2, r4, r7
      real(wp) :: hb(0:NCAT)
      integer :: ierr, i, j

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s2, NCAT, 2)
      call alloc_sis2(s4, NCAT, 4)
      call alloc_sis2(s7, NCAT, 7)
      call alloc_new(o2)
      call alloc_new(o4)
      call alloc_new(o7)
      do j = 1, NY
         do i = 1, NX
            call set_column(s2, i, j, 2, 0.7_wp, 1.1_wp, 0.1_wp, -2.0_wp, -12.0_wp, -9.0_wp)
            call set_column(s4, i, j, 2, 0.7_wp, 1.1_wp, 0.1_wp, -2.0_wp, -12.0_wp, -9.0_wp)
            call set_column(s7, i, j, 2, 0.7_wp, 1.1_wp, 0.1_wp, -2.0_wp, -12.0_wp, -9.0_wp)
         end do
      end do
      call convert(s2, p, hb, o2, r2, ierr)
      call check(error, ierr, ICE_CONV_OK, "2 -> 7")
      call convert(s4, p, hb, o4, r4, ierr)
      call check(error, ierr, ICE_CONV_OK, "4 -> 7")
      call convert(s7, p, hb, o7, r7, ierr)
      call check(error, ierr, ICE_CONV_OK, "7 -> 7")
      if (allocated(error)) return
      call check(error, relerr(r2%mass_ice_out, r2%mass_ice_in) < 1.0e-14_wp .and. &
                 relerr(r4%mass_ice_out, r4%mass_ice_in) < 1.0e-14_wp, "mass exact for resampled columns")
      call check(error, abs(r2%heat_rel_change - r7%heat_rel_change) < 0.1_wp .and. &
                 abs(r4%heat_rel_change - r7%heat_rel_change) < 0.03_wp, &
                 "resampling adds no energy beyond the layering error")
      call check(error, o2%q(1, 1, 2, 1) < o2%q(1, 1, 2, NILYR) .and. o4%q(1, 1, 2, 1) < o4%q(1, 1, 2, NILYR), &
                 "top layer colder (lower enthalpy) than the base layer after resampling")
   end subroutine test_resample_converter

   subroutine test_rebin(error)
      !! Re-binning: (a) two old categories whose thicknesses fall in ONE new bin merge with
      !! area/volume summed and the layer energy volume-weighted; (b) a thickness outside its
      !! old bin goes to the right new bin; (c) a slab above the top bound goes to the top category.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), q1, q2, vsum, want
      integer :: ierr, k

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 3, NILYR)
      call alloc_new(o)
      ! cell (1,1): categories 1 and 2 both thin (0.20 and 0.40 m SIS2-density thickness) -> new bin 1
      call set_column(s, 1, 1, 1, 0.2_wp, 0.20_wp, 0.0_wp, -2.0_wp, -4.0_wp, -3.0_wp)
      call set_column(s, 1, 1, 2, 0.3_wp, 0.40_wp, 0.05_wp, -2.0_wp, -10.0_wp, -3.0_wp)
      call set_column(s, 1, 1, 3, 0.4_wp, 1.00_wp, 0.0_wp, -2.0_wp, -8.0_wp, -3.0_wp)
      ! cell (2,1): old category 1 is really 0.9 m -> new bin 2 ; category 3 is 1500 m -> top bin
      call set_column(s, 2, 1, 1, 0.5_wp, 0.90_wp, 0.0_wp, -2.0_wp, -6.0_wp, -3.0_wp)
      call set_column(s, 2, 1, 3, 0.5_wp, 1500.0_wp, 0.0_wp, -2.0_wp, -6.0_wp, -3.0_wp)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      ! (a)
      call check(error, abs(o%a(1, 1, 1) - 0.5_wp) < 1.0e-15_wp, "merged area = 0.2 + 0.3")
      vsum = (0.2_wp*0.20_wp + 0.3_wp*0.40_wp)*RHOI_OLD/p%c%rhoi
      call check(error, relerr(o%v(1, 1, 1), vsum) < 1.0e-14_wp, "merged volume = sum of volumes")
      call check(error, abs(o%vs(1, 1, 1) - 0.3_wp*0.05_wp*RHOS/p%c%rhos) < 1.0e-15_wp, "snow volume carried")
      call check(error, o%a(1, 1, 3) == 0.0_wp .and. o%v(1, 1, 3) == 0.0_wp, "unused middle bin exactly empty")
      call check(error, abs(o%a(1, 1, 2) - 0.4_wp) < 1.0e-15_wp, "thick category stays in bin 2")
      ! volume-weighted layer energy (top-down layer k = SIS2 layer NILYR+1-k)
      do k = 1, NILYR
         call col_q(p, -2.0_wp, -4.0_wp, k, q1)
         call col_q(p, -2.0_wp, -10.0_wp, k, q2)
         want = (0.2_wp*0.20_wp*q1 + 0.3_wp*0.40_wp*q2)/(0.2_wp*0.20_wp + 0.3_wp*0.40_wp)
         if (relerr(o%q(1, 1, 1, k), want) > 1.0e-12_wp) then
            call test_failed(error, "merged layer enthalpy is the volume-weighted mean")
            exit
         end if
      end do
      ! (b), (c)
      call check(error, o%a(2, 1, 1) == 0.0_wp .and. abs(o%a(2, 1, 2) - 0.5_wp) < 1.0e-15_wp, "0.9 m goes to bin 2")
      call check(error, abs(o%a(2, 1, NCAT) - 0.5_wp) < 1.0e-15_wp .and. &
                 relerr(o%v(2, 1, NCAT), 0.5_wp*1500.0_wp*RHOI_OLD/p%c%rhoi) < 1.0e-14_wp, "above the top bound -> top category")
   end subroutine test_rebin

   subroutine col_q(p, Tb, Tt, knew, q)
      !! Expected new-layer enthalpy (top-down layer `knew`) of a SIS2 column linear in T.
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: Tb, Tt
      integer, intent(in) :: knew
      real(wp), intent(out) :: q
      real(wp) :: sp(NILYR + 1), T
      call ice_salinity_profile(NILYR, p%c%saltmax, p%c%min_salin, sp)
      T = Tb + (Tt - Tb)*(real(NILYR + 1 - knew, wp) - 0.5_wp)/real(NILYR, wp)
      q = ice_q_from_T(T, sp(knew), p%c, p%thermo%l_brine)
   end subroutine col_q

   subroutine test_empty(error)
      !! Land/open-water cells and sub-`puny` categories come out as exact zero with the SAME
      !! don't-care tracers as the uniform seed; a tiny occupied category is dropped and counted.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o, ref
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr, c, k, sst

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 3, NILYR)
      call alloc_new(o)
      call alloc_new(ref)
      call set_column(s, 2, 2, 2, 0.5_wp, 1.0_wp, 0.1_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      call set_column(s, 4, 3, 1, 1.0e-13_wp, 0.3_wp, 0.0_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call ice_seed_uniform_impl(NX, NY, NCAT, NILYR, 1, p, hb, 0.0_wp, 0.0_wp, 0.0_wp, min(TF, p%c%Timelt), TF, &
                                 ref%a, ref%v, ref%vs, ref%ts, ref%q, ref%s, ref%qs, sst)
      call check(error, sst, ICE_SEED_OK, "reference seed")
      do c = 1, NCAT
         if (c == 2) cycle
         call check(error, all(o%a(:, :, c) == 0.0_wp) .and. all(o%v(:, :, c) == 0.0_wp) .and. &
                    all(o%vs(:, :, c) == 0.0_wp), "unoccupied categories exactly zero")
         call check(error, all(o%ts(:, :, c) == ref%ts(:, :, c)), "empty Tsfcn as the seed")
         do k = 1, NILYR
            call check(error, all(o%q(:, :, c, k) == ref%q(:, :, c, k)) .and. all(o%s(:, :, c, k) == ref%s(:, :, c, k)), &
                       "empty layer tracers as the seed")
         end do
      end do
      call check(error, rep%n_dropped == 1, "tiny category dropped and counted")
      call check(error, o%a(2, 2, 2) == 0.5_wp, "the occupied cell is kept")
      call check(error, o%a(1, 1, 2) == 0.0_wp .and. o%vs(1, 1, 2) == 0.0_wp, "land cell empty")
   end subroutine test_empty

   subroutine test_lumped(error)
      !! Legacy `ncat = 1` lumped mode: concentration 1 where `m_ice > 0`, `vicen = m/rhoi`,
      !! the slab re-binned by its thickness; cells without ice stay empty.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 1, NILYR)
      call alloc_new(o)
      s%ps(:, :, 0) = 1.0_wp      ! never maintained in lumped mode
      call set_column(s, 3, 3, 1, 1.0_wp, 0.8_wp, 0.2_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      s%ps(3, 3, 0) = 1.0_wp
      s%ps(3, 3, 1) = 0.0_wp      ! garbage part_size must be ignored
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call check(error, o%a(3, 3, 2) == 1.0_wp .and. sum(o%a(3, 3, :)) == 1.0_wp, "concentration 1 in the bin of h = 0.79 m")
      call check(error, relerr(p%c%rhoi*o%v(3, 3, 2), RHOI_OLD*0.8_wp) <= MASS_TOL, "lumped mass exact")
      call check(error, sum(o%a(1, 1, :)) == 0.0_wp, "iceless cell empty")
   end subroutine test_lumped

   subroutine test_rescale(error)
      !! Total area 1.2 is rescaled to 1 keeping the volumes (mass exact), and counted.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 3, NILYR)
      call alloc_new(o)
      call set_column(s, 1, 1, 1, 0.6_wp, 0.3_wp, 0.0_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      call set_column(s, 1, 1, 2, 0.6_wp, 1.0_wp, 0.0_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call check(error, abs(sum(o%a(1, 1, :)) - 1.0_wp) < 4.0e-16_wp, "area rescaled to 1")
      call check(error, rep%n_rescaled == 1, "counted")
      call check(error, relerr(rep%mass_ice_out, rep%mass_ice_in) < 1.0e-14_wp, "mass unchanged by the rescale")
   end subroutine test_rescale

   subroutine test_cap(error)
      !! A layer so close to melting that the NEW brine profile (fresher top layers) has a
      !! lower melting point is capped at it, and the cap is counted.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), sp(NILYR + 1), T
      logical :: bad
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 2, NILYR)
      call alloc_new(o)
      call set_column(s, 1, 1, 2, 0.5_wp, 1.0_wp, 0.0_wp, -0.002_wp, -0.002_wp, -0.1_wp)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call ice_salinity_profile(NILYR, p%c%saltmax, p%c%min_salin, sp)
      call ice_T_from_q(o%q(1, 1, 2, 1), o%s(1, 1, 2, 1), p%c, p%thermo%l_brine, T, bad)
      call check(error, rep%n_t_capped > 0, "cap counted")
      call check(error, abs(T + p%c%depressT*sp(1)) < 1.0e-10_wp, "top layer sits at the new melting point")
      call check(error, ieee_is_finite(o%q(1, 1, 2, 1)) .and. o%q(1, 1, 2, 1) < 0.0_wp, "finite negative enthalpy")
   end subroutine test_cap

   subroutine test_refusals(error)
      !! NaN, negative and extent errors are refused with the outputs untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 3, NILYR)
      call alloc_new(o)
      call set_column(s, 2, 2, 2, 0.5_wp, 1.0_wp, 0.1_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      s%en(5, 4, 1, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_NONFINITE, "NaN enthalpy refused")
      call check(error, all(o%a == 7.0_wp) .and. all(o%q == 7.0_wp), "outputs untouched on refusal")
      s%en(5, 4, 1, 3) = -1.0e8_wp
      s%mi(1, 1, 1) = -1.0_wp
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_NEGATIVE, "negative mass refused")
      s%mi(1, 1, 1) = 0.0_wp
      s%ps(1, 1, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_NONFINITE, "NaN part_size refused")
      s%ps(1, 1, 1) = 0.0_wp
      call ice_restart_convert_sis2(NX, NY, s%nc, s%nk, NCAT, NILYR + 1, 1, p, hb, TF, s%ps, s%mi, s%ms, s%en, s%sa, &
                                    s%es, o%a, o%v, o%vs, o%ts, o%q(:, :, :, 1:NILYR), o%s(:, :, :, 1:NILYR), &
                                    o%qs, rep, ierr)
      call check(error, ierr, ICE_CONV_BAD_EXTENT, "layer count disagreeing with the parameters refused")
   end subroutine test_refusals

   subroutine test_thin_film(error)
      !! Lumped mode: a film thinner than `hi_min` is a PARTIAL cover of `hi_min` ice, mass
      !! exact; a slab thicker than `hi_min` is concentration 1.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT), hmin
      integer :: ierr

      p = params()
      hmin = p%thermo%hi_min
      call bounds(p, hb)
      call alloc_sis2(s, 1, NILYR)
      call alloc_new(o)
      call set_column(s, 1, 1, 1, 1.0_wp, 0.004_wp, 0.0_wp, -2.0_wp, -6.0_wp, -3.0_wp)   ! 4 mm film
      s%mi(2, 1, 1) = 1.0e-6_wp                                                            ! 1e-6 kg/m2
      s%en(2, 1, 1, :) = ice_enth_from_ts(-3.0_wp, SBULK)
      call set_column(s, 3, 1, 1, 1.0_wp, 0.8_wp, 0.0_wp, -2.0_wp, -6.0_wp, -3.0_wp)      ! thick slab (0.79 m: bin 2)
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call check(error, relerr(o%a(1, 1, 1), RHOI_OLD*0.004_wp/(p%c%rhoi*hmin)) < 1.0e-14_wp, "aicen = m/(rhoi*hi_min)")
      call check(error, relerr(o%v(1, 1, 1)/o%a(1, 1, 1), hmin) < 1.0e-13_wp, "thickness = hi_min")
      call check(error, relerr(p%c%rhoi*o%v(1, 1, 1), RHOI_OLD*0.004_wp) <= MASS_TOL, "film mass exact")
      call check(error, o%a(2, 1, 1) > 0.0_wp .and. o%a(2, 1, 1) < 1.0e-6_wp .and. &
                 relerr(p%c%rhoi*o%v(2, 1, 1), 1.0e-6_wp) <= MASS_TOL, "1e-6 kg/m2 film: tiny cover, mass exact")
      call check(error, o%a(3, 1, 2) == 1.0_wp, "thick slab keeps concentration 1")
      call check(error, relerr(rep%mass_ice_out + rep%mass_ice_dropped, rep%mass_ice_in) < 1.0e-14_wp, "mass identity")
   end subroutine test_thin_film

   subroutine test_dropped(error)
      !! in = out + dropped, for ice and snow, with snow on an ice-free category and a sub-`puny` category.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, 3, NILYR)
      call alloc_new(o)
      call set_column(s, 2, 2, 2, 0.5_wp, 1.0_wp, 0.1_wp, -2.0_wp, -9.0_wp, -4.0_wp)   ! kept
      call set_column(s, 3, 2, 1, 0.3_wp, 0.5_wp, 0.2_wp, -2.0_wp, -9.0_wp, -4.0_wp)
      s%mi(3, 2, 1) = 0.0_wp                                                             ! snow, no ice
      call set_column(s, 4, 3, 1, 1.0e-13_wp, 0.3_wp, 0.2_wp, -2.0_wp, -9.0_wp, -4.0_wp)  ! below puny
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call check(error, rep%mass_snow_dropped > 0.0_wp .and. rep%mass_ice_dropped >= 0.0_wp, "snow dropped")
      call check(error, relerr(rep%mass_snow_out + rep%mass_snow_dropped, rep%mass_snow_in) < 1.0e-14_wp, &
                 "snow: in = out + dropped")
      call check(error, relerr(rep%mass_ice_out + rep%mass_ice_dropped, rep%mass_ice_in) < 1.0e-14_wp, &
                 "ice: in = out + dropped")
      call check(error, abs(rep%mass_snow_dropped - (0.3_wp*RHOS*0.2_wp + 1.0e-13_wp*RHOS*0.2_wp)) &
                 < 1.0e-14_wp*rep%mass_snow_in, "dropped snow = the two lost categories")
   end subroutine test_dropped

   subroutine test_tsfc_snow(error)
      !! `Tsfcn` is the snow temperature where snow is present, the top-layer ice temperature otherwise.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(sis2_t) :: s
      type(newstate_t) :: o
      type(ice_convert_report_t) :: rep
      real(wp) :: hb(0:NCAT)
      integer :: ierr

      p = params()
      call bounds(p, hb)
      call alloc_sis2(s, NCAT, NILYR)
      call alloc_new(o)
      call set_column(s, 1, 1, 2, 0.5_wp, 1.0_wp, 0.1_wp, -2.0_wp, -12.0_wp, -6.0_wp)   ! snow at -6, top ice ~-11.6
      call set_column(s, 2, 1, 2, 0.5_wp, 1.0_wp, 0.0_wp, -2.0_wp, -12.0_wp, -6.0_wp)   ! no snow
      call convert(s, p, hb, o, rep, ierr)
      call check(error, ierr, ICE_CONV_OK, "converted")
      if (allocated(error)) return
      call check(error, abs(o%ts(1, 1, 2) + 6.0_wp) < 1.0e-9_wp, "under snow: snow temperature")
      call check(error, abs(o%ts(2, 1, 2) - (-12.0_wp + 10.0_wp*0.5_wp/real(NILYR, wp))) < 1.0e-9_wp, &
                 "no snow: top ice layer temperature")
   end subroutine test_tsfc_snow

   ! ------------------------------------------------------------------ read path
   subroutine write_reg(fn, reg, d, meta, ierr)
      character(len=*), intent(in) :: fn
      type(restart_registry_t), intent(inout) :: reg
      type(decomp_t), intent(in) :: d
      type(ocean_restart_metadata_t), intent(inout) :: meta
      integer, intent(out) :: ierr
      integer :: lun

      open (newunit=lun, file=fn, status="replace")
      close (lun, status="delete")
      call ocean_restart_write_local(fn, reg, d, meta, 1.0_wp, 1, 1, ierr)
   end subroutine write_reg

   subroutine read_start(fn, st, lg, d, meta, Tf_in, status, rep)
      !! The documented pattern: the slot in the run registry, the legacy in a SEPARATE registry.
      character(len=*), intent(in) :: fn
      type(ice_cice_state_t), target, intent(inout) :: st
      type(ice_legacy_t), target, intent(inout) :: lg
      type(decomp_t), intent(in) :: d
      type(ocean_restart_metadata_t), intent(inout) :: meta
      real(wp), intent(in) :: Tf_in
      integer, intent(out) :: status
      type(ice_convert_report_t), intent(out) :: rep
      type(restart_registry_t) :: rn, rl
      real(wp) :: t
      integer :: step, osc, ierr

      call st%register_restart(rn)
      call lg%register_restart(rl)
      call ocean_restart_read_local(fn, rn, d, meta, t, step, osc, ierr)
      if (ierr /= 0) error stop "read_start: run registry read"
      call ocean_restart_read_local(fn, rl, d, meta, t, step, osc, ierr)
      if (ierr /= 0) error stop "read_start: legacy registry read"
      call ice_cice_restart_resume(st, lg, rn, rl, Tf_in, status, rep)
   end subroutine read_start

   subroutine reg_sis2(rw, lw, skip)
      !! Write-side registration of the SIS2 tags (no probes), except the tag `skip`.
      type(restart_registry_t), intent(inout) :: rw
      type(ice_legacy_t), target, intent(in) :: lw
      character(len=*), intent(in) :: skip
      integer :: k

      if (skip /= "ice_part_size") call rw%register_3d("ice_part_size", lw%part_size, 0, lw%nx, lw%ny)
      if (skip /= "ice_m_ice") call rw%register_3d("ice_m_ice", lw%m_ice, 0, lw%nx, lw%ny)
      if (skip /= "ice_m_snow") call rw%register_3d("ice_m_snow", lw%m_snow, 0, lw%nx, lw%ny)
      do k = 1, lw%nk
         if (skip /= "ice_enth_ice_k"//char(48 + k)) then
            call rw%register_3d("ice_enth_ice_k"//char(48 + k), lw%enth_ice(:, :, :, k), 0, lw%nx, lw%ny)
         end if
         if (skip /= "ice_sal_ice_k"//char(48 + k)) then
            call rw%register_3d("ice_sal_ice_k"//char(48 + k), lw%sal_ice(:, :, :, k), 0, lw%nx, lw%ny)
         end if
      end do
      if (skip /= "ice_enth_snow_k1") call rw%register_3d("ice_enth_snow_k1", lw%enth_snow(:, :, :, 1), 0, lw%nx, lw%ny)
   end subroutine reg_sis2

   subroutine fill_legacy(lw, seedv)
      type(ice_legacy_t), intent(inout) :: lw
      integer(int64), intent(in) :: seedv
      integer(int64) :: seed
      integer :: i, j, c, k

      seed = seedv
      do j = 1, NY
         do i = 1, NX
            do c = 1, lw%ncat
               seed = lcg_next(seed)
               lw%part_size(i, j, c) = 0.25_wp*lcg_unit(seed)
               lw%part_size(i, j, 0) = 1.0_wp - sum(lw%part_size(i, j, 1:lw%ncat))
               lw%m_ice(i, j, c) = RHOI_OLD*(0.3_wp + 1.5_wp*c*lcg_unit(seed))
               lw%m_snow(i, j, c) = RHOS*0.1_wp*lcg_unit(seed)
               do k = 1, lw%nk
                  lw%enth_ice(i, j, c, k) = ice_enth_from_ts(-2.0_wp - 3.0_wp*k*lcg_unit(seed), SBULK)
               end do
               lw%enth_snow(i, j, c, 1) = ice_enth_from_ts(-5.0_wp, 0.0_wp)
            end do
         end do
      end do
   end subroutine fill_legacy

   subroutine test_read_path(error)
      !! The three starts, through a real restart file.  (1) A file with ONLY SIS2 tags is converted:
      !! the slot equals a direct converter call on the same arrays, and `state_version` is reset to 2.
      !! (2) A v2 file is read back bitwise and NOT converted even if SIS2 tags are present too.
      !! (3) A file without any ice tag cold-starts to open water.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 2, NXP = NX - 2*NG, NYP = NY - 2*NG
      type(ice_cice_state_t), target :: st
      type(ice_legacy_t), target :: lw, lr
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      type(restart_registry_t) :: rw
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      type(ice_convert_report_t) :: rep, rep2
      type(newstate_t) :: want
      real(wp), target :: dummy(NX, NY)
      real(wp) :: hb(0:NCAT)
      integer :: ierr, status, lun
      character(len=*), parameter :: FN = "test_cice_restart_convert.nc"

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = params()
      call st%init(g, p, ierr)
      call check(error, ierr, 0, "slot init")
      if (allocated(error)) return
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = NG
      dummy = 1.0_wp

      ! ---- (1) SIS2-only file ----
      call lw%init(NX, NY, 3, 4, NG)
      call fill_legacy(lw, 31337_int64)
      call reg_sis2(rw, lw, "")
      call rw%register_2d("dummy_marker", dummy, 0, NX, NY)
      call write_reg(FN, rw, d, meta, ierr)
      call check(error, ierr, 0, "write SIS2 file")
      if (allocated(error)) return
      call lr%init(NX, NY, 3, 4, NG)
      lr%part_size = -9.0_wp     ! junk: must be overwritten by the read
      st%state_version = 1.0_wp
      call read_start(FN, st, lr, d, meta, TF, status, rep)
      call check(error, status, ICE_RESUME_CONVERTED, "converted on resume")
      if (allocated(error)) return
      call check(error, st%state_version == ICE_STATE_VERSION, "version reset after conversion")
      call bounds(p, hb)
      call alloc_new(want)
      call ice_restart_convert_sis2(NX, NY, 3, 4, NCAT, NILYR, 1, p, hb, TF, lw%part_size, lw%m_ice, lw%m_snow, &
                                    lw%enth_ice, lw%sal_ice, lw%enth_snow, want%a, want%v, want%vs, want%ts, &
                                    want%q, want%s, want%qs, rep2, ierr)
      call check(error, ierr, ICE_CONV_OK, "reference conversion")
      if (allocated(error)) return
      call check(error, all(st%aicen == want%a) .and. all(st%vicen == want%v) .and. all(st%vsnon == want%vs), &
                 "slot a/v/vs equal the direct conversion")
      call check(error, all(st%qicen == want%q) .and. all(st%sicen == want%s) .and. all(st%qsnon == want%qs) &
                 .and. all(st%Tsfcn == want%ts), "slot tracers equal the direct conversion")
      call check(error, relerr(rep%mass_ice_out + rep%mass_ice_dropped, rep%mass_ice_in) < 1.0e-13_wp, &
                 "mass in the report is exact")
      call check(error, rep%heat_rel_change == rep2%heat_rel_change .and. ieee_is_finite(rep%heat_rel_change), &
                 "heat change reported")

      ! ---- (2) v2 file written from the converted slot; SIS2 tags present TOO: v2 must win ----
      call rw%clear()
      call st%register_restart(rw)
      call reg_sis2(rw, lw, "")
      call write_reg(FN, rw, d, meta, ierr)
      call check(error, ierr, 0, "write v2 file")
      block
         type(ice_cice_state_t), target :: st2
         type(ice_legacy_t), target :: l2
         call st2%init(g, p, ierr)
         call l2%init(NX, NY, 3, 4, NG)
         call read_start(FN, st2, l2, d, meta, TF, status, rep)
         call check(error, status, ICE_RESUME_V2, "v2 read, not converted")
         call check(error, all(st2%aicen == st%aicen) .and. all(st2%qicen == st%qicen) .and. &
                    all(st2%vsnon == st%vsnon) .and. all(st2%sicen == st%sicen), "v2 state bitwise")
         call st2%destroy()
         call l2%destroy()
      end block

      ! ---- (3) no ice tags at all ----
      call rw%clear()
      call rw%register_2d("dummy_marker", dummy, 0, NX, NY)
      call write_reg(FN, rw, d, meta, ierr)
      block
         type(ice_cice_state_t), target :: st3
         type(ice_legacy_t), target :: l3
         call st3%init(g, p, ierr)
         call l3%init(NX, NY, 3, 4, NG)
         st3%aicen = 0.9_wp
         st3%state_version = 1.0_wp
         call read_start(FN, st3, l3, d, meta, TF, status, rep)
         call check(error, status, ICE_RESUME_COLD, "cold start")
         call check(error, all(st3%aicen == 0.0_wp) .and. all(st3%vicen == 0.0_wp) .and. all(st3%aice0 == 1.0_wp), &
                    "cold start is open water")
         call check(error, st3%state_version == ICE_STATE_VERSION, "version reset after cold start")
         call st3%destroy()
         call l3%destroy()
      end block
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      call lw%destroy()
      call lr%destroy()
      call st%destroy()
   end subroutine test_read_path

   subroutine test_read_refusals(error)
      !! Partial or mismatched files are REFUSED with the state arrays untouched: SIS2 masses missing,
      !! one SIS2 layer missing, a file with MORE layers than the legacy was sized for, legacy extents
      !! (nghost) different from the slot's, a partial v2 set, a v2 structure mismatch, a newer version.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NG = 2, NXP = NX - 2*NG, NYP = NY - 2*NG
      type(ice_cice_state_t), target :: st, a
      type(ice_legacy_t), target :: lw, lw5, lr
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      type(restart_registry_t) :: rw
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      type(ice_convert_report_t) :: rep
      integer :: ierr, status, lun
      character(len=*), parameter :: FN = "test_cice_restart_refuse.nc"

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = params()
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = NG
      call lw%init(NX, NY, 3, 4, NG)
      call fill_legacy(lw, 5150_int64)
      call lw5%init(NX, NY, 3, 5, NG)
      call fill_legacy(lw5, 5151_int64)

      ! masses missing
      call rw%clear()
      call reg_sis2(rw, lw, "ice_m_ice")
      call write_reg(FN, rw, d, meta, ierr)
      call try_sis2(error, "masses missing")
      ! a layer missing
      call rw%clear()
      call reg_sis2(rw, lw, "ice_sal_ice_k2")
      call write_reg(FN, rw, d, meta, ierr)
      call try_sis2(error, "layer missing")
      ! five layers in the file, legacy sized for four
      call rw%clear()
      call reg_sis2(rw, lw5, "")
      call write_reg(FN, rw, d, meta, ierr)
      call try_sis2(error, "extra layer")
      ! legacy extents differ (nghost): good complete file, wrong nghost
      call rw%clear()
      call reg_sis2(rw, lw, "")
      call write_reg(FN, rw, d, meta, ierr)
      call st%init(g, p, ierr)
      call lr%init(NX, NY, 3, 4, NG + 1)
      call read_start(FN, st, lr, d, meta, TF, status, rep)
      call check(error, status, ICE_RESUME_FAILED, "legacy extents mismatch refused")
      call check(error, all(st%aicen == 0.0_wp), "state untouched on extents refusal")
      call st%destroy()
      call lr%destroy()

      ! partial v2: version + aicen only
      call a%init(g, p, ierr)
      call rw%clear()
      call rw%register_scalar("ice_state_version", a%state_version)
      call rw%register_3d("ice_aicen", a%aicen, 0, NX, NY)
      call write_reg(FN, rw, d, meta, ierr)
      call try_v2(error, "partial v2")
      ! structure mismatch
      a%nilyr_tag = 6.0_wp
      call rw%clear()
      call a%register_restart(rw)
      call write_reg(FN, rw, d, meta, ierr)
      call try_v2(error, "nilyr mismatch")
      a%nilyr_tag = real(NILYR, wp)
      ! newer version
      a%state_version = 3.0_wp
      call rw%clear()
      call a%register_restart(rw)
      call write_reg(FN, rw, d, meta, ierr)
      call try_v2(error, "version 3")
      call a%destroy()
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      call lw%destroy()
      call lw5%destroy()

   contains

      subroutine try_sis2(err, what)
         type(error_type), allocatable, intent(inout) :: err
         character(len=*), intent(in) :: what
         call st%init(g, p, ierr)
         call lr%init(NX, NY, 3, 4, NG)
         call read_start(FN, st, lr, d, meta, TF, status, rep)
         if (.not. allocated(err)) call check(err, status, ICE_RESUME_FAILED, what//": refused")
         if (.not. allocated(err)) call check(err, all(st%aicen == 0.0_wp) .and. all(st%vicen == 0.0_wp), &
                                              what//": state untouched")
         call st%destroy()
         call lr%destroy()
      end subroutine try_sis2

      subroutine try_v2(err, what)
         type(error_type), allocatable, intent(inout) :: err
         character(len=*), intent(in) :: what
         type(ice_cice_state_t), target :: b
         type(ice_legacy_t), target :: l
         call b%init(g, p, ierr)
         call l%init(NX, NY, 3, 4, NG)
         call read_start(FN, b, l, d, meta, TF, status, rep)
         if (.not. allocated(err)) call check(err, status, ICE_RESUME_FAILED, what//": refused")
         call b%destroy()
         call l%destroy()
      end subroutine try_v2

   end subroutine test_read_refusals

end module test_cice_restart_convert
