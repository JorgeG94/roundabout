!! Ice-only kernel benchmark (no ocean): the device-shaped column and cell kernels of the CICE port
!! on a synthetic, device-resident state (SoA, horizontal indices first).
!!
!! Usage: `bench_ice_kernels [--band] [nx ny [nwarm nrep [nvrows [tol]]]]` (defaults 512 512 10 50 96).
!! `--band` replaces the random ~30 %-ice-free mask by a spatially coherent one (ice in the first ~20 % of the rows).
!!
!! Every loop is ONE `do concurrent` launch (per (n, j, i) for column physics, per (j, i) for cell-coupled
!! steps), timed as the median of `nrep` launches after `nwarm` warm-ups with `!$acc wait` before each
!! timer read.  Before timing, each loop is VERIFIED once against a host evaluation of the same cell
!! kernel on the first `nvrows` rows (a fast wrong answer is worthless): bitwise on gfortran CPU builds, 1e-12
!! scale-aware on the offload and nvfortran builds (FMA contraction).  Any failure stops with a non-zero exit code.  The synthetic
!! state has ~30 % ice-free cells.  Loops that are not on this base yet are listed as STUB.
module bench_ice_kernels_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX, ICE_NTRCR_MAX
   use rdb_ice_params, only: ice_constants_t, ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, &
                             ICE_CONDUCT_BUBBLY, ice_q_from_T, ice_qsno_from_T, ice_ntrcr, ice_trcr_get, ice_trcr_set_mean
   use rdb_ice_bl99, only: ICE_NMAT_MAX, ice_conductivity_core, ice_matrix_calc_Tsfc_core, ice_tridiag_solver_core
   use rdb_ice_vprofile, only: ice_init_vertical_profile
   use rdb_ice_shortwave, only: ice_albedo_fhtan, ice_albedos_ccsm3, ice_absorbed_ccsm3
   use rdb_ice_itd_shift, only: itd_aggregate, itd_aggregate_cell
   use rdb_ice_thermo_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf, ICE_SFC_BULK
   implicit none
   private

   public :: agg_out_t, mean_get, bl99_cat, sw_cat, sfc_cat, agg_packed_cell, rnd_unit
   public :: bl99_impl, sw_impl, sfc_impl, agg_packed_impl, agg_cell_impl

   type :: agg_out_t
      !! Cell-mean outputs of one aggregate variant (component arrays are mapped one by one).
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), aice0(:, :)
      real(wp), allocatable :: Tsfc(:, :), qice(:, :, :), sice(:, :, :), qsno(:, :, :)
      integer(int32), allocatable :: ierr(:, :)
   end type agg_out_t

contains

   function rnd_unit(st) result(u)
      !! xorshift64 -> uniform in [0, 1); `st` is advanced.
      integer(int64), intent(inout) :: st
      real(wp) :: u

      st = ieor(st, ishft(st, 13))
      st = ieor(st, ishft(st, -7))
      st = ieor(st, ishft(st, 17))
      u = real(ishft(st, -11), wp)/9007199254740992.0_wp
   end function rnd_unit

   pure function mean_get(it, i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc) result(v)
      !! Cell-mean tracer `it` (packed order) read back from the cell-mean arrays.
      integer, intent(in) :: it, i, j, nx, ny, nilyr, nslyr
      real(wp), intent(in) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      real(wp) :: v

      if (it == 1) then
         v = Tc(i, j)
      else if (it <= 1 + nilyr) then
         v = qc(i, j, it - 1)
      else if (it <= 1 + nilyr + nslyr) then
         v = qsc(i, j, it - 1 - nilyr)
      else
         v = sc(i, j, it - 1 - nilyr - nslyr)
      end if
   end function mean_get

   ! ------------------------------------------------------------------ 1. BL99 per (n, j, i)
   pure subroutine bl99_cat(i, j, n, nx, ny, ncat, nilyr, nslyr, c, dt, aicen, vicen, vsnon, qicen, sicen, qsnon, &
                            Tsfcn, Tsf_o, Tin_o, ierr_o)
      !! One category column through the BL99 pieces, the surface iteration mimicked by 6 passes of
      !! profile -> conductivity -> matrix -> solve.  Gather once, scatter once; outputs per category.
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: Tsf_o(nx, ny, ncat), Tin_o(nx, ny, ncat, nilyr)
      integer(int32), intent(inout) :: ierr_o(nx, ny, ncat)
      integer :: k, it
      integer(int32) :: ierr
      real(wp) :: zq(ICE_NILYR_MAX), zqw(ICE_NILYR_MAX), zS(ICE_NILYR_MAX), zT(ICE_NILYR_MAX), Tm(ICE_NILYR_MAX)
      real(wp) :: zqs(ICE_NSLYR_MAX), zqsw(ICE_NSLYR_MAX), zTs(ICE_NSLYR_MAX)
      real(wp) :: kh(ICE_NMAT_MAX), etai(ICE_NILYR_MAX), Isw(ICE_NILYR_MAX), etas(ICE_NSLYR_MAX)
      real(wp) :: Ssw(ICE_NSLYR_MAX), sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX)
      real(wp) :: x(ICE_NMAT_MAX), hin, hilyr, hsn, hslyr, einit, ci, Tsf
      logical :: l_snow

      ierr = 0_int32
      do k = 1, nilyr
         zq(k) = qicen(i, j, n, k)
         zS(k) = sicen(i, j, n, k)
         Isw(k) = 0.5_wp
      end do
      do k = 1, nslyr
         zqs(k) = qsnon(i, j, n, k)
         Ssw(k) = 1.0_wp
      end do
      Tsf = Tsfcn(i, j, n)
      x = 0.0_wp
      do k = 1, nilyr
         zqw(k) = zq(k)
      end do
      do k = 1, nslyr
         zqsw(k) = zqs(k)
      end do
      call ice_init_vertical_profile(aicen(i, j, n), vicen(i, j, n), vsnon(i, j, n), nilyr, nslyr, c, .true., zS, &
                                     zqw, zqsw, zT, zTs, Tm, hin, hilyr, hsn, hslyr, einit, ierr)
      if (ierr == 0_int32) then
         l_snow = hslyr > c%hs_min
         do k = 1, nilyr
            ci = c%cp_ice - c%Lfresh*Tm(k)/(zT(k)*zT(k))
            etai(k) = dt/(c%rhoi*ci*hilyr)
         end do
         do k = 1, nslyr
            etas(k) = dt/(c%rhos*c%cp_ice*max(hslyr, c%hs_min))
         end do
         do it = 1, 6                                       ! the surface iteration: conductivity, matrix, solve
            call ice_conductivity_core(l_snow, nilyr, nslyr, hilyr, hslyr, zT, zS, c, ICE_CONDUCT_BUBBLY, kh)
            call ice_matrix_calc_Tsfc_core(l_snow, .true., nilyr, nslyr, Tsf, -1.8_wp, 20.0_wp - Tsf, -9.0_wp, zT, zTs, kh, &
                                           Ssw, Isw, etai, etas, sb, dg, sp, rh)
            call ice_tridiag_solver_core(nilyr + nslyr + 1, sb, dg, sp, rh, x, ierr)
            if (ierr /= 0_int32) exit
            Tsf = max(min(x(1), 0.0_wp), -60.0_wp)
         end do
      end if
      Tsf_o(i, j, n) = Tsf
      do k = 1, nilyr
         Tin_o(i, j, n, k) = x(1 + nslyr + k)
      end do
      ierr_o(i, j, n) = ierr
   end subroutine bl99_cat

   subroutine bl99_impl(nx, ny, ncat, nilyr, nslyr, c, dt, puny, aicen, vicen, vsnon, qicen, sicen, qsnon, &
                        Tsfcn, Tsf_o, Tin_o, ierr_o)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: dt, puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: Tsf_o(nx, ny, ncat), Tin_o(nx, ny, ncat, nilyr)
      integer(int32), intent(inout) :: ierr_o(nx, ny, ncat)
      integer :: i, j, n

      do concurrent(n=1:ncat, j=1:ny, i=1:nx)
         if (aicen(i, j, n) > puny) then
            call bl99_cat(i, j, n, nx, ny, ncat, nilyr, nslyr, c, dt, aicen, vicen, vsnon, qicen, sicen, qsnon, &
                          Tsfcn, Tsf_o, Tin_o, ierr_o)
         end if
      end do
   end subroutine bl99_impl

   ! ------------------------------------------------------------------ 2. shortwave per (n, j, i)
   pure subroutine sw_cat(i, j, n, nx, ny, ncat, nilyr, c, albicev, albicei, albsnowv, albsnowi, fhtan, aicen, vicen, &
                          vsnon, Tsfcn, swv, swn, fsfc_o, fint_o, fthru_o, isw_o)
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nilyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: albicev, albicei, albsnowv, albsnowi, fhtan
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(in) :: swv(nx, ny), swn(nx, ny)
      real(wp), intent(inout) :: fsfc_o(nx, ny, ncat), fint_o(nx, ny, ncat), fthru_o(nx, ny, ncat)
      real(wp), intent(inout) :: isw_o(nx, ny, ncat, nilyr)
      integer :: k
      real(wp) :: hi, hs, alv_i, ali_i, alv_s, ali_s, asnow, alv, ali, fs, fi, ft, iw(ICE_NILYR_MAX)

      hi = vicen(i, j, n)/aicen(i, j, n)
      hs = vsnon(i, j, n)/aicen(i, j, n)
      call ice_albedos_ccsm3(hi, hs, Tsfcn(i, j, n), albicev, albicei, albsnowv, albsnowi, fhtan, c%snowpatch, &
                             alv_i, ali_i, alv_s, ali_s, asnow, alv, ali)
      call ice_absorbed_ccsm3(nilyr, hi, asnow, swv(i, j), swn(i, j), alv_i, ali_i, alv_s, ali_s, fs, fi, ft, iw)
      fsfc_o(i, j, n) = fs
      fint_o(i, j, n) = fi
      fthru_o(i, j, n) = ft
      do k = 1, nilyr
         isw_o(i, j, n, k) = iw(k)
      end do
   end subroutine sw_cat

   subroutine sw_impl(nx, ny, ncat, nilyr, c, albicev, albicei, albsnowv, albsnowi, fhtan, puny, aicen, vicen, vsnon, &
                      Tsfcn, swv, swn, fsfc_o, fint_o, fthru_o, isw_o)
      integer, intent(in) :: nx, ny, ncat, nilyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: albicev, albicei, albsnowv, albsnowi, fhtan, puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(in) :: swv(nx, ny), swn(nx, ny)
      real(wp), intent(inout) :: fsfc_o(nx, ny, ncat), fint_o(nx, ny, ncat), fthru_o(nx, ny, ncat)
      real(wp), intent(inout) :: isw_o(nx, ny, ncat, nilyr)
      integer :: i, j, n

      do concurrent(n=1:ncat, j=1:ny, i=1:nx)
         if (aicen(i, j, n) > puny) then
            call sw_cat(i, j, n, nx, ny, ncat, nilyr, c, albicev, albicei, albsnowv, albsnowi, fhtan, aicen, vicen, &
                        vsnon, Tsfcn, swv, swn, fsfc_o, fint_o, fthru_o, isw_o)
         end if
      end do
   end subroutine sw_impl

   ! ------------------------------------------------------------------ 3. aggregate per (j, i)
   pure subroutine agg_packed_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                   Tf, puny, aice_o, vice_o, vsno_o, aice0_o, Tc, qc, sc, qsc, ierr_o)
      !! The packed-gather form: gathers the (ntrcr, ncat) block per cell, calls `itd_aggregate`, scatters.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: Tf, puny
      real(wp), intent(inout) :: aice_o(nx, ny), vice_o(nx, ny), vsno_o(nx, ny), aice0_o(nx, ny)
      real(wp), intent(inout) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      integer(int32), intent(inout) :: ierr_o(nx, ny)
      real(wp) :: a(ICE_NCAT_MAX), v(ICE_NCAT_MAX), s(ICE_NCAT_MAX), tr(ICE_NTRCR_MAX, ICE_NCAT_MAX), trc(ICE_NTRCR_MAX)
      real(wp) :: aice, vice, vsno, aice0
      integer(int32) :: ierr
      integer :: n, it

      do n = 1, ncat
         a(n) = aicen(i, j, n)
         v(n) = vicen(i, j, n)
         s(n) = vsnon(i, j, n)
         do it = 1, ice_ntrcr(nilyr, nslyr)
            tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
      ierr = 0_int32
      call itd_aggregate(ncat, ICE_NTRCR_MAX, nilyr, nslyr, a, v, s, tr, Tf, puny, aice, vice, vsno, aice0, trc, ierr)
      aice_o(i, j) = aice
      vice_o(i, j) = vice
      vsno_o(i, j) = vsno
      aice0_o(i, j) = aice0
      do it = 1, ice_ntrcr(nilyr, nslyr)
         call ice_trcr_set_mean(it, i, j, nx, ny, nilyr, nslyr, trc(it), Tc, qc, qsc, sc)
      end do
      ierr_o(i, j) = ierr
   end subroutine agg_packed_cell

   subroutine agg_packed_impl(nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, Tf, puny, &
                              aice_o, vice_o, vsno_o, aice0_o, Tc, qc, sc, qsc, ierr_o)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: Tf, puny
      real(wp), intent(inout) :: aice_o(nx, ny), vice_o(nx, ny), vsno_o(nx, ny), aice0_o(nx, ny)
      real(wp), intent(inout) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      integer(int32), intent(inout) :: ierr_o(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         call agg_packed_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, Tf, &
                              puny, aice_o, vice_o, vsno_o, aice0_o, Tc, qc, sc, qsc, ierr_o)
      end do
   end subroutine agg_packed_impl

   subroutine agg_cell_impl(nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, Tf, puny, &
                            aice_o, vice_o, vsno_o, aice0_o, Tc, qc, sc, qsc, ierr_o)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: Tf, puny
      real(wp), intent(inout) :: aice_o(nx, ny), vice_o(nx, ny), vsno_o(nx, ny), aice0_o(nx, ny)
      real(wp), intent(inout) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      integer(int32), intent(inout) :: ierr_o(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         call itd_aggregate_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, Tf, &
                                 puny, aice_o(i, j), vice_o(i, j), vsno_o(i, j), aice0_o(i, j), Tc, qc, sc, qsc, ierr_o)
      end do
   end subroutine agg_cell_impl

   ! ------------------------------------------------------------------ 6. surface heat flux per (n, j, i)
   pure subroutine sfc_cat(i, j, n, nx, ny, ncat, c, emis, aicen, Tsfcn, fsw, rhoa, flw, potT, Qa, shc, lhc, &
                           fsurf_o, dfsurf_o, flat_o, fsens_o, ierr_o)
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: emis
      real(wp), intent(in) :: aicen(nx, ny, ncat), Tsfcn(nx, ny, ncat), fsw(nx, ny, ncat)
      real(wp), intent(in) :: rhoa(nx, ny), flw(nx, ny), potT(nx, ny), Qa(nx, ny), shc(nx, ny), lhc(nx, ny)
      real(wp), intent(inout) :: fsurf_o(nx, ny, ncat), dfsurf_o(nx, ny, ncat), flat_o(nx, ny, ncat), fsens_o(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr_o(nx, ny, ncat)
      integer(int32) :: ierr
      real(wp) :: flwout, fsens, flat, fsurf, dfs, dflw, dfsens, dflat

      ierr = 0_int32
      call surface_heat_flux(Tsfcn(i, j, n), fsw(i, j, n), rhoa(i, j), flw(i, j), potT(i, j), Qa(i, j), shc(i, j), &
                             lhc(i, j), emis, c, ICE_SFC_BULK, 0.0_wp, 0.0_wp, flwout, fsens, flat, fsurf, ierr)
      call dsurface_heat_flux_dTsf(Tsfcn(i, j, n), rhoa(i, j), shc(i, j), lhc(i, j), emis, c, ICE_SFC_BULK, 0.0_wp, &
                                   dfs, dflw, dfsens, dflat, ierr)
      fsurf_o(i, j, n) = fsurf
      dfsurf_o(i, j, n) = dfs
      flat_o(i, j, n) = flat
      fsens_o(i, j, n) = fsens
      ierr_o(i, j, n) = ierr
   end subroutine sfc_cat

   subroutine sfc_impl(nx, ny, ncat, c, emis, puny, aicen, Tsfcn, fsw, rhoa, flw, potT, Qa, shc, lhc, &
                       fsurf_o, dfsurf_o, flat_o, fsens_o, ierr_o)
      integer, intent(in) :: nx, ny, ncat
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: emis, puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), Tsfcn(nx, ny, ncat), fsw(nx, ny, ncat)
      real(wp), intent(in) :: rhoa(nx, ny), flw(nx, ny), potT(nx, ny), Qa(nx, ny), shc(nx, ny), lhc(nx, ny)
      real(wp), intent(inout) :: fsurf_o(nx, ny, ncat), dfsurf_o(nx, ny, ncat), flat_o(nx, ny, ncat), fsens_o(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr_o(nx, ny, ncat)
      integer :: i, j, n

      do concurrent(n=1:ncat, j=1:ny, i=1:nx)
         if (aicen(i, j, n) > puny) then
            call sfc_cat(i, j, n, nx, ny, ncat, c, emis, aicen, Tsfcn, fsw, rhoa, flw, potT, Qa, shc, lhc, &
                         fsurf_o, dfsurf_o, flat_o, fsens_o, ierr_o)
         end if
      end do
   end subroutine sfc_impl

end module bench_ice_kernels_mod

program bench_ice_kernels
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit, error_unit
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ice_err_summary_t, ice_err_reduce, ice_err_fatal_set
   use rdb_ice_params, only: ice_constants_t, ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ice_q_from_T, &
                             ice_qsno_from_T, ice_ntrcr, ice_trcr_get
   use rdb_ice_shortwave, only: ice_albedo_fhtan
   use rdb_ice_lateral_melt, only: ice_lateral_melt_cell, ice_lateral_melt_step
   use rdb_ice_itd_shift, only: itd_aggregate_cell
   use bench_ice_kernels_mod
   implicit none

   integer, parameter :: NCAT = 5, NILYR = 7, NSLYR = 1
   integer :: nx, ny, nwarm, nrep, nvrows
   type(ice_params_t) :: p
   type(ice_constants_t) :: c
   real(wp) :: fhtan, dt, tf, emis
   real(wp), allocatable :: aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :)
   real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
   real(wp), allocatable :: swv(:, :), swn(:, :), rhoa(:, :), flw(:, :), potT(:, :), Qa(:, :), shc(:, :), lhc(:, :)
   real(wp), allocatable :: rside(:, :), fresh(:, :), fsalt(:, :), fhocn(:, :), meltl(:, :)
   integer(int32), allocatable :: ierr_cell(:, :), ierr_cat(:, :, :)
   ! outputs (device) and host references
   real(wp), allocatable :: Tsf_o(:, :, :), Tin_o(:, :, :, :), fsfc(:, :, :), fint(:, :, :), fthru(:, :, :), isw(:, :, :, :)
   real(wp), allocatable :: fsurf(:, :, :), dfsurf(:, :, :), flat(:, :, :), fsens(:, :, :)
   real(wp), allocatable :: Tsf_r(:, :, :), Tin_r(:, :, :, :), fsfc_r(:, :, :), fint_r(:, :, :), fthru_r(:, :, :)
   real(wp), allocatable :: isw_r(:, :, :, :), fsurf_r(:, :, :), dfsurf_r(:, :, :), flat_r(:, :, :), fsens_r(:, :, :)
   integer(int32), allocatable :: ierr_r(:, :, :)
   type(agg_out_t) :: op, oc, orf
   real(wp), allocatable :: a0h(:, :, :), v0h(:, :, :), vs0h(:, :, :)
   real(wp), allocatable :: fresh_r(:, :), fsalt_r(:, :), fhocn_r(:, :), meltl_r(:, :)
   integer(int32), allocatable :: ie2_r(:, :)
   type(ice_err_summary_t) :: summary
   real(wp) :: worst_rel = 0.0_wp
   real(wp) :: tol, ms(8), per_cell(8), t_ext(8), nitems(8)
   character(len=40) :: names(8)
   logical :: ok, verified(8)
   integer :: i, j, k, n, nv, loop, order(7), il
   integer(int64) :: st
   real(wp) :: u, Tsf, Si, Tm
   character(len=32) :: arg
   integer :: nargs, npos
   character(len=32) :: pos(8)
   logical :: band = .false.
   integer(int64) :: nact

   nx = 512; ny = 512; nwarm = 10; nrep = 50; nvrows = 96
   nargs = command_argument_count()
   npos = 0
   do k = 1, nargs                                       ! "--band" anywhere; the rest are positional
      call get_command_argument(k, arg)
      if (trim(arg) == "--band") then
         band = .true.
      else
         npos = npos + 1
         pos(min(npos, 8)) = arg
      end if
   end do
   if (npos >= 2) then
      read (pos(1), *) nx
      read (pos(2), *) ny
   end if
   if (npos >= 4) then
      read (pos(3), *) nwarm
      read (pos(4), *) nrep
   end if
   if (npos >= 5) read (pos(5), *) nvrows
   nvrows = min(nvrows, ny)
#if defined(RDB_GPU_OFFLOAD) || defined(__NVCOMPILER)
   tol = 1.0e-12_wp   ! FMA contraction (device, and host code compiled by nvfortran -fast): measured 5e-16
#else
   tol = 0.0_wp
#endif
   if (npos >= 6) read (pos(6), *) tol                  ! optional override of the verification tolerance

   p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, NSLYR)
   if (p%bad) error stop "bench_ice_kernels: params refused"
   c = p%c
   fhtan = ice_albedo_fhtan(p%sw%ahmax)
   dt = 3600.0_wp
   tf = -1.8_wp
   emis = 0.95_wp

   allocate (aicen(nx, ny, NCAT), vicen(nx, ny, NCAT), vsnon(nx, ny, NCAT), Tsfcn(nx, ny, NCAT))
   allocate (qicen(nx, ny, NCAT, NILYR), sicen(nx, ny, NCAT, NILYR), qsnon(nx, ny, NCAT, NSLYR))
   allocate (swv(nx, ny), swn(nx, ny), rhoa(nx, ny), flw(nx, ny), potT(nx, ny), Qa(nx, ny), shc(nx, ny), lhc(nx, ny))
   allocate (rside(nx, ny), fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny))
   allocate (ierr_cell(nx, ny), ierr_cat(nx, ny, NCAT))
   allocate (Tsf_o(nx, ny, NCAT), Tin_o(nx, ny, NCAT, NILYR), fsfc(nx, ny, NCAT), fint(nx, ny, NCAT), fthru(nx, ny, NCAT))
   allocate (isw(nx, ny, NCAT, NILYR), fsurf(nx, ny, NCAT), dfsurf(nx, ny, NCAT), flat(nx, ny, NCAT), fsens(nx, ny, NCAT))
   allocate (Tsf_r(nx, ny, NCAT), Tin_r(nx, ny, NCAT, NILYR), fsfc_r(nx, ny, NCAT), fint_r(nx, ny, NCAT), fthru_r(nx, ny, NCAT))
   allocate (isw_r(nx, ny, NCAT, NILYR), fsurf_r(nx, ny, NCAT), dfsurf_r(nx, ny, NCAT), flat_r(nx, ny, NCAT))
   allocate (fsens_r(nx, ny, NCAT), ierr_r(nx, ny, NCAT))
   call alloc_agg(op); call alloc_agg(oc); call alloc_agg(orf)
   allocate (a0h(nx, ny, NCAT), v0h(nx, ny, NCAT), vs0h(nx, ny, NCAT))
   allocate (fresh_r(nx, ny), fsalt_r(nx, ny), fhocn_r(nx, ny), meltl_r(nx, ny), ie2_r(nx, ny))

   ! ---- synthetic state: ~30 % ice-free cells, plausible profiles (enthalpy from a temperature profile)
   st = 88172645463325252_int64
   do j = 1, ny
      do i = 1, nx
         u = rnd_unit(st)
         do n = 1, NCAT
            if (merge(real(j, wp) > 0.2_wp*real(ny, wp), u < 0.3_wp, band)) then   ! ice-free cell
               aicen(i, j, n) = 0.0_wp; vicen(i, j, n) = 0.0_wp; vsnon(i, j, n) = 0.0_wp
               Tsfcn(i, j, n) = 0.0_wp
            else
               aicen(i, j, n) = 0.02_wp + 0.16_wp*rnd_unit(st)
               vicen(i, j, n) = aicen(i, j, n)*(0.1_wp + 0.9_wp*real(n, wp) + 0.8_wp*rnd_unit(st))
               vsnon(i, j, n) = merge(aicen(i, j, n)*0.3_wp*rnd_unit(st), 0.0_wp, rnd_unit(st) > 0.3_wp)
               Tsfcn(i, j, n) = -2.0_wp - 28.0_wp*rnd_unit(st)
            end if
            Tsf = Tsfcn(i, j, n)
            do k = 1, NILYR
               Si = 0.5_wp + 5.0_wp*rnd_unit(st)
               Tm = -c%depressT*Si
               sicen(i, j, n, k) = Si
               qicen(i, j, n, k) = ice_q_from_T(min(Tsf + (tf - Tsf)*real(k, wp)/real(NILYR + 1, wp), Tm - 0.01_wp), Si, c, .true.)
            end do
            qsnon(i, j, n, 1) = ice_qsno_from_T(min(Tsf, -0.5_wp), c)
         end do
         swv(i, j) = 40.0_wp + 120.0_wp*rnd_unit(st)
         swn(i, j) = 30.0_wp + 80.0_wp*rnd_unit(st)
         rhoa(i, j) = 1.25_wp + 0.05_wp*rnd_unit(st)
         flw(i, j) = 180.0_wp + 80.0_wp*rnd_unit(st)
         potT(i, j) = 245.0_wp + 25.0_wp*rnd_unit(st)
         Qa(i, j) = 1.0e-3_wp + 1.0e-3_wp*rnd_unit(st)
         shc(i, j) = 5.0_wp + 10.0_wp*rnd_unit(st)
         lhc(i, j) = 1.0e4_wp + 2.0e4_wp*rnd_unit(st)
         rside(i, j) = 1.0e-4_wp*rnd_unit(st)
      end do
   end do
   fresh = 0.0_wp; fsalt = 0.0_wp; fhocn = 0.0_wp; meltl = 0.0_wp
   ierr_cell = 0_int32; ierr_cat = 0_int32
   Tsf_o = 0.0_wp; Tin_o = 0.0_wp; fsfc = 0.0_wp; fint = 0.0_wp; fthru = 0.0_wp; isw = 0.0_wp
   fsurf = 0.0_wp; dfsurf = 0.0_wp; flat = 0.0_wp; fsens = 0.0_wp
   Tsf_r = 0.0_wp; Tin_r = 0.0_wp; fsfc_r = 0.0_wp; fint_r = 0.0_wp; fthru_r = 0.0_wp; isw_r = 0.0_wp
   fsurf_r = 0.0_wp; dfsurf_r = 0.0_wp; flat_r = 0.0_wp; fsens_r = 0.0_wp; ierr_r = 0_int32
   a0h = aicen; v0h = vicen; vs0h = vsnon
   fresh_r = 0.0_wp; fsalt_r = 0.0_wp; fhocn_r = 0.0_wp; meltl_r = 0.0_wp

   !$acc enter data copyin(c, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, swv, swn, rhoa, flw, potT, Qa, shc, lhc) &
   !$acc copyin(rside, fresh, fsalt, fhocn, meltl, ierr_cell, ierr_cat, Tsf_o, Tin_o, fsfc, fint, fthru, isw) &
   !$acc copyin(fsurf, dfsurf, flat, fsens)
   call map_agg(op); call map_agg(oc)

   nact = count(aicen > c%puny)
   write (output_unit, '(a,a,a,i0,a,f5.1,a)') "ice mask: ", merge("coherent band (first 20 % of rows)", "random (~30 % ice-free cells)     ", band), &
      "; active category columns ", nact, " = ", 100.0_wp*real(nact, wp)/(real(nx, wp)*ny*NCAT), " % of nx*ny*ncat"
   write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0)') "bench_ice_kernels: nx=", nx, " ny=", ny, " ncat=", NCAT, " nilyr=", NILYR, &
      " nslyr=", NSLYR
#ifdef RDB_GPU_OFFLOAD
   write (output_unit, '(a)') "build: offload (device)"
#else
   write (output_unit, '(a)') "build: host (do concurrent on the CPU)"
#endif
   write (output_unit, '(a,i0,a,i0,a,i0,a)') "timing: ", nwarm, " warm-up launches, median of ", nrep, &
      " launches; verification on rows 1:", nvrows, " against the host cell kernel"

   names = [character(len=40) :: "bl99 (6 passes) per (n,j,i)", "shortwave per (n,j,i)", "aggregate packed-gather per (j,i)", &
            "aggregate in-place per (j,i)", "lateral_melt_step per (j,i)", "ice_err_reduce (clean)", &
            "surface heat flux+deriv per (n,j,i)", "STUB atmosphere/strength (not on base)"]
   nitems = [real(nx, wp)*ny*NCAT, real(nx, wp)*ny*NCAT, real(nx, wp)*ny, real(nx, wp)*ny, real(nx, wp)*ny, &
             real(nx, wp)*ny, real(nx, wp)*ny*NCAT, 0.0_wp]
   ms = -1.0_wp
   verified = .false.
   ok = .true.
   nv = nvrows

   ! execution order: lateral melt (mutates the state) last
   order = [1, 2, 3, 4, 7, 6, 5]
   do il = 1, 7
      loop = order(il)
      call run_loop(loop)
   end do

   write (output_unit, '(a)') ""
   write (output_unit, '(a)') "loop                                   ms/launch   Mcell-cat/s  ns/(i,j,n)  ms/25km step (1440x1152)  verified"
   do loop = 1, 8
      if (ms(loop) < 0.0_wp) then
         write (output_unit, '(a40,a)') names(loop), "   -- not benchmarked --"
         cycle
      end if
      per_cell(loop) = ms(loop)*1.0e6_wp/(real(nx, wp)*ny*NCAT)
      t_ext(loop) = ms(loop)*(1440.0_wp*1152.0_wp)/(real(nx, wp)*ny)
      write (output_unit, '(a40,f12.4,f14.2,f12.2,f18.2,l12)') names(loop), ms(loop), nitems(loop)/(ms(loop)*1.0e-3_wp)/1.0e6_wp, &
         per_cell(loop), t_ext(loop), verified(loop)
      if (any(loop == [1, 2, 7])) write (output_unit, '(a40,f12.2,a)') "", ms(loop)*1.0e6_wp/real(max(nact, 1_int64), wp), &
         " ns per ACTIVE cell-category"

   end do
   write (output_unit, '(a,es9.2,a,es9.2)') "verification: tolerance ", tol, ", largest relative deviation seen ", worst_rel
   write (output_unit, '(a)') "(Mcell-cat/s = launch items per second: nx*ny*ncat for per-(n,j,i) loops, nx*ny for per-(j,i) loops;"
   write (output_unit, '(a)') " ns/(i,j,n) divides every loop by nx*ny*ncat so the loops are comparable.)"

   !$acc exit data delete(c, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, swv, swn, rhoa, flw, potT, Qa, shc, lhc) &
   !$acc delete(rside, fresh, fsalt, fhocn, meltl, ierr_cell, ierr_cat, Tsf_o, Tin_o, fsfc, fint, fthru, isw) &
   !$acc delete(fsurf, dfsurf, flat, fsens)
   call unmap_agg(op); call unmap_agg(oc)
   if (.not. ok) then
      write (error_unit, '(a)') "bench_ice_kernels: VERIFICATION FAILED"
      error stop 1
   end if

contains

   subroutine alloc_agg(o)
      type(agg_out_t), intent(inout) :: o

      allocate (o%aice(nx, ny), o%vice(nx, ny), o%vsno(nx, ny), o%aice0(nx, ny))
      allocate (o%Tsfc(nx, ny), o%qice(nx, ny, NILYR), o%sice(nx, ny, NILYR), o%qsno(nx, ny, NSLYR))
      allocate (o%ierr(nx, ny))
      o%aice = 0.0_wp; o%vice = 0.0_wp; o%vsno = 0.0_wp; o%aice0 = 0.0_wp
      o%Tsfc = 0.0_wp; o%qice = 0.0_wp; o%sice = 0.0_wp; o%qsno = 0.0_wp
      o%ierr = 0_int32
   end subroutine alloc_agg

   subroutine map_agg(o)
      type(agg_out_t), intent(inout) :: o

      !$acc enter data copyin(o%aice, o%vice, o%vsno, o%aice0, o%Tsfc, o%qice, o%sice, o%qsno, o%ierr)
   end subroutine map_agg

   subroutine unmap_agg(o)
      type(agg_out_t), intent(inout) :: o

      !$acc exit data delete(o%aice, o%vice, o%vsno, o%aice0, o%Tsfc, o%qice, o%sice, o%qsno, o%ierr)
   end subroutine unmap_agg

   subroutine launch(which)
      !! One launch of loop `which` on the mapped arrays, followed by a device wait.
      integer, intent(in) :: which

      select case (which)
      case (1)
         call bl99_impl(nx, ny, NCAT, NILYR, NSLYR, c, dt, c%puny, aicen, vicen, vsnon, qicen, sicen, qsnon, Tsfcn, Tsf_o, &
                        Tin_o, ierr_cat)
      case (2)
         call sw_impl(nx, ny, NCAT, NILYR, c, p%sw%albicev, p%sw%albicei, p%sw%albsnowv, p%sw%albsnowi, fhtan, c%puny, &
                      aicen, vicen, vsnon, Tsfcn, swv, swn, fsfc, fint, fthru, isw)
      case (3)
         call agg_packed_impl(nx, ny, NCAT, NILYR, NSLYR, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, tf, c%puny, &
                              op%aice, op%vice, op%vsno, op%aice0, op%Tsfc, op%qice, op%sice, op%qsno, op%ierr)
      case (4)
         call agg_cell_impl(nx, ny, NCAT, NILYR, NSLYR, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, tf, c%puny, &
                            oc%aice, oc%vice, oc%vsno, oc%aice0, oc%Tsfc, oc%qice, oc%sice, oc%qsno, oc%ierr)
      case (5)
         call ice_lateral_melt_step(nx, ny, NCAT, NILYR, NSLYR, 1, nx, 1, ny, dt, p, rside, aicen, vicen, vsnon, qicen, &
                                    qsnon, fresh, fsalt, fhocn, meltl, ierr_cell)
      case (6)
         call ice_err_reduce(ierr_cell, nx, ny, 1, nx, 1, ny, 0, 0, nx, ice_err_fatal_set(.true.), summary)
      case (7)
         call sfc_impl(nx, ny, NCAT, c, emis, c%puny, aicen, Tsfcn, fsfc, rhoa, flw, potT, Qa, shc, lhc, fsurf, dfsurf, &
                       flat, fsens, ierr_cat)
      end select
      !$acc wait
   end subroutine launch

   subroutine time_loop(which)
      integer, intent(in) :: which
      real(wp) :: t(nrep), tmp
      integer(int64) :: c0, c1, rate
      integer :: r, a, b

      call system_clock(count_rate=rate)
      do r = 1, nwarm
         call launch(which)
      end do
      do r = 1, nrep
         call system_clock(c0)
         call launch(which)
         call system_clock(c1)
         t(r) = real(c1 - c0, wp)/real(rate, wp)
      end do
      do a = 2, nrep                                       ! insertion sort
         tmp = t(a)
         b = a - 1
         do while (b >= 1)
            if (t(b) <= tmp) exit
            t(b + 1) = t(b)
            b = b - 1
         end do
         t(b + 1) = tmp
      end do
      ms(which) = 1.0e3_wp*t((nrep + 1)/2)
   end subroutine time_loop

   logical function close_to(d, h)
      real(wp), intent(in) :: d, h

      if (.not. (ieee_is_finite(d) .and. ieee_is_finite(h))) then
         close_to = .false.
         return
      end if
      worst_rel = max(worst_rel, abs(d - h)/max(abs(d), abs(h), 1.0_wp))
      if (tol == 0.0_wp) then
         close_to = (d == h)
      else
         close_to = abs(d - h) <= tol*max(abs(d), abs(h), 1.0_wp)
      end if
   end function close_to

   subroutine fail(msg)
      character(len=*), intent(in) :: msg

      write (error_unit, '(a)') "VERIFY FAIL: "//msg
      ok = .false.
   end subroutine fail

   subroutine run_loop(which)
      integer, intent(in) :: which
      integer :: ii, jj, nn, kk, nbad

      ! ---- one verification launch, outputs pulled to the host, compared with the host cell kernel on rows 1:nv
      call launch(which)
      nbad = 0
      select case (which)
      case (1)
         !$acc update self(Tsf_o, Tin_o, ierr_cat)
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
                  if (aicen(ii, jj, nn) > c%puny) &
                     call bl99_cat(ii, jj, nn, nx, ny, NCAT, NILYR, NSLYR, c, dt, aicen, vicen, vsnon, qicen, sicen, qsnon, &
                                   Tsfcn, Tsf_r, Tin_r, ierr_r)
               end do
            end do
         end do
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
                  if (.not. close_to(Tsf_o(ii, jj, nn), Tsf_r(ii, jj, nn))) nbad = nbad + 1
                  do kk = 1, NILYR
                     if (.not. close_to(Tin_o(ii, jj, nn, kk), Tin_r(ii, jj, nn, kk))) nbad = nbad + 1
                  end do
                  if (ierr_cat(ii, jj, nn) /= ierr_r(ii, jj, nn)) nbad = nbad + 1
               end do
            end do
         end do
         if (any(ierr_cat /= 0_int32)) call fail("bl99 set an error bit")
         ! the device result must also not be the initial (unwritten) state
         if (count(Tin_o(:, 1:nv, :, 1) /= 0.0_wp) == 0) call fail("bl99 wrote nothing")
      case (2)
         !$acc update self(fsfc, fint, fthru, isw)
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
                  if (aicen(ii, jj, nn) > c%puny) &
                     call sw_cat(ii, jj, nn, nx, ny, NCAT, NILYR, c, p%sw%albicev, p%sw%albicei, p%sw%albsnowv, &
                                 p%sw%albsnowi, fhtan, aicen, vicen, vsnon, Tsfcn, swv, swn, fsfc_r, fint_r, fthru_r, isw_r)
               end do
            end do
         end do
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
                  if (.not. (close_to(fsfc(ii, jj, nn), fsfc_r(ii, jj, nn)) .and. close_to(fint(ii, jj, nn), fint_r(ii, jj, nn)) &
                             .and. close_to(fthru(ii, jj, nn), fthru_r(ii, jj, nn)))) nbad = nbad + 1
                  do kk = 1, NILYR
                     if (.not. close_to(isw(ii, jj, nn, kk), isw_r(ii, jj, nn, kk))) nbad = nbad + 1
                  end do
               end do
            end do
         end do
      case (3, 4)
         if (which == 3) then
            call launch(4)                                  ! both variants have run before the comparison
            call verify_agg(nbad)
         else
            nbad = merge(0, 1, verified(3))                 ! compared together with the packed variant
         end if
      case (5)
         ! host reference on the INITIAL state (the verification launch above applied exactly one melt)
         fresh_r = 0.0_wp; fsalt_r = 0.0_wp; fhocn_r = 0.0_wp; meltl_r = 0.0_wp; ie2_r = 0_int32
         !$acc update self(aicen, vicen, vsnon, fresh, fsalt, fhocn, meltl, ierr_cell)
         do jj = 1, nv
            do ii = 1, nx
               call ice_lateral_melt_cell(ii, jj, nx, ny, NCAT, NILYR, NSLYR, rside(ii, jj), dt, c%rhoi, c%rhos, &
                                          p%thermo%ice_ref_salinity, c%puny, a0h, v0h, vs0h, qicen, qsnon, &
                                          fresh_r, fsalt_r, fhocn_r, meltl_r, ie2_r)
            end do
         end do
         do jj = 1, nv
            do ii = 1, nx
               if (.not. (close_to(fresh(ii, jj), fresh_r(ii, jj)) .and. close_to(fsalt(ii, jj), fsalt_r(ii, jj)) .and. &
                          close_to(fhocn(ii, jj), fhocn_r(ii, jj)) .and. close_to(meltl(ii, jj), meltl_r(ii, jj)))) nbad = nbad + 1
               do nn = 1, NCAT
                  if (.not. (close_to(aicen(ii, jj, nn), a0h(ii, jj, nn)) .and. close_to(vicen(ii, jj, nn), v0h(ii, jj, nn)) &
                             .and. close_to(vsnon(ii, jj, nn), vs0h(ii, jj, nn)))) nbad = nbad + 1
               end do
            end do
         end do
         if (count(fresh(:, 1:nv) /= 0.0_wp) == 0) call fail("lateral melt wrote nothing")
         if (any(ierr_cell /= 0_int32)) call fail("lateral melt set an error bit")
      case (6)
         if (summary%mask /= 0_int32 .or. summary%first_idx /= 0_int64) nbad = nbad + 1
      case (7)
         !$acc update self(fsurf, dfsurf, flat, fsens, ierr_cat)
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
                  if (aicen(ii, jj, nn) > c%puny) &
                     call sfc_cat(ii, jj, nn, nx, ny, NCAT, c, emis, aicen, Tsfcn, fsfc, rhoa, flw, potT, Qa, shc, lhc, &
                                  fsurf_r, dfsurf_r, flat_r, fsens_r, ierr_r)
               end do
            end do
         end do
         do nn = 1, NCAT
            do jj = 1, nv
               do ii = 1, nx
              if (.not. (close_to(fsurf(ii, jj, nn), fsurf_r(ii, jj, nn)) .and. close_to(dfsurf(ii, jj, nn), dfsurf_r(ii, jj, nn)) &
                             .and. close_to(flat(ii, jj, nn), flat_r(ii, jj, nn)) .and. close_to(fsens(ii, jj, nn), &
                                                                                              fsens_r(ii, jj, nn)))) nbad = nbad + 1
               end do
            end do
         end do
         if (any(ierr_cat /= 0_int32)) call fail("surface flux set an error bit")
      end select
      verified(which) = (nbad == 0)
      if (nbad /= 0) then
         write (error_unit, '(a,a,a,i0,a)') "VERIFY FAIL: loop '", trim(names(which)), "': ", nbad, " mismatching values"
         ok = .false.
      end if
      ! ---- timing (after verification, outside it)
      call time_loop(which)
   end subroutine run_loop

   subroutine verify_agg(nbad)
      !! Aggregate loops: the device result of the variant `which` against the host in-place kernel on rows 1:nv,
      !! and (once both exist) the packed variant against the in-place variant on EVERY cell.
      integer, intent(inout) :: nbad
      integer :: ii, jj, ntr, it
      real(wp) :: r(4)

      !$acc update self(op%aice, op%vice, op%vsno, op%aice0, op%Tsfc, op%qice, op%sice, op%qsno, op%ierr)
      !$acc update self(oc%aice, oc%vice, oc%vsno, oc%aice0, oc%Tsfc, oc%qice, oc%sice, oc%qsno, oc%ierr)
      orf%ierr = 0_int32
      ntr = ice_ntrcr(NILYR, NSLYR)
      do jj = 1, nv
         do ii = 1, nx
            call itd_aggregate_cell(ii, jj, nx, ny, NCAT, NILYR, NSLYR, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, tf, &
                                    c%puny, r(1), r(2), r(3), r(4), orf%Tsfc, orf%qice, orf%sice, orf%qsno, orf%ierr)
            orf%aice(ii, jj) = r(1); orf%vice(ii, jj) = r(2); orf%vsno(ii, jj) = r(3); orf%aice0(ii, jj) = r(4)
         end do
      end do
      do jj = 1, nv
         do ii = 1, nx
            if (.not. (close_to(oc%aice(ii, jj), orf%aice(ii, jj)) .and. close_to(oc%vice(ii, jj), orf%vice(ii, jj)) .and. &
                       close_to(oc%vsno(ii, jj), orf%vsno(ii, jj)) .and. close_to(oc%aice0(ii, jj), orf%aice0(ii, jj)))) &
               nbad = nbad + 1
            do it = 1, ntr
               if (.not. close_to(mean_get(it, ii, jj, nx, ny, NILYR, NSLYR, oc%Tsfc, oc%qice, oc%sice, oc%qsno), &
                                  mean_get(it, ii, jj, nx, ny, NILYR, NSLYR, orf%Tsfc, orf%qice, orf%sice, &
                                           orf%qsno))) nbad = nbad + 1
            end do
         end do
      end do
      ! packed vs in-place, every cell and tracer (device vs device)
      do jj = 1, ny
         do ii = 1, nx
            if (.not. (close_to(oc%aice(ii, jj), op%aice(ii, jj)) .and. close_to(oc%vice(ii, jj), op%vice(ii, jj)) .and. &
                       close_to(oc%vsno(ii, jj), op%vsno(ii, jj)) .and. close_to(oc%aice0(ii, jj), op%aice0(ii, jj)))) &
               nbad = nbad + 1
            do it = 1, ntr
               if (.not. close_to(mean_get(it, ii, jj, nx, ny, NILYR, NSLYR, oc%Tsfc, oc%qice, oc%sice, oc%qsno), &
                                  ice_trcr_get(it, ii, jj, 1, nx, ny, 1, NILYR, NSLYR, op%Tsfc, op%qice, op%sice, &
                                               op%qsno))) nbad = nbad + 1
            end do
         end do
      end do
      if (any(oc%ierr /= 0_int32) .or. any(op%ierr /= 0_int32)) call fail("aggregate set an error bit")
   end subroutine verify_agg

end program bench_ice_kernels
