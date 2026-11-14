!! `rdb_ice_rad_state` and `rdb_ice_radiation` (D7.3): the radiation slot, the band split, the category-parallel
!! step and the cell aggregate of the CCSM3 shortwave (Briegleb et al. 2004).
!!
!! * `slot_*`              -- lifecycle, `bytes()` against its closed form, device canary on EVERY array (the
!!                            host values are copied in, a kernel adds a per-element canary, every element comes back);
!! * `band_split_*`        -- visible + near-IR add back to the input, the visible share is `sw_vis_frac`;
!! * `shim_fails_loud`     -- uninitialised / unconfigured slot and shape mismatches return `ierr`;
!! * `closure_*`           -- per category, `sw = reflected + fswsfc + fswint + fswthru` and `sum_k Iswabs = fswint`,
!!                            for visible-only, near-IR-only and mixed forcing;
!! * `aggregate_*`         -- the grid-box albedos equal the hand area-weighted sum; ice-free cell = exact 0;
!! * `defaults_*`          -- empty categories get the upstream defaults; zero shortwave gives exactly zero absorption;
!! * `nan_*`               -- a non-finite input flags its own category and cell only;
!! * `device_vs_host_*`    -- EVERY output of every (i, j, n[, k]) of a random interior range against an independent
!!                            host reference (the scalar kernels called directly), ring untouched;
!! * `oracle_*`            -- 10 multi-category cells (nilyr 7) against upstream `icepack_step_radiation`;
!! * `driver_throughput`   -- 512 x 512 x 5 x 7 cell-categories per second.
!!
!! Tolerances are measured: the worst relative difference is printed by each comparison.  On the GPU build FMA
!! contraction applies to host code too, so those comparisons are scale-aware ~1e-12 instead of CPU-tight.
module test_cice_radiation
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_sw_params_t
   use rdb_ice_rad_state, only: ice_rad_state_t, ICE_RAD_INIT_OK, ICE_RAD_INIT_BAD_EXTENT
   use rdb_ice_radiation, only: ice_sw_band_split, ice_step_radiation, ICE_RAD_OK, ICE_RAD_ERR_NOT_INIT, &
                                ICE_RAD_ERR_NOT_CONFIGURED, ICE_RAD_ERR_SHAPE
   use rdb_ice_shortwave, only: ice_albedo_fhtan, ice_albedos_ccsm3, ice_absorbed_ccsm3, ICE_ALBOCN, ICE_SW_PUNY
   implicit none
   private

   public :: collect_cice_radiation_tests

   integer, parameter :: NC = 5, NL = 7
   real(wp), parameter :: SNOWPATCH = 0.02_wp
   real(wp), parameter :: SENT = -777.0_wp
   integer(int32), parameter :: ESENT = 1073741824_int32
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: RTOL = 1.0e-12_wp
      !! device/host scale-aware tolerance (FMA on the GPU build)
   real(wp), parameter :: CTOL = 1.0e-13_wp
      !! closure / aggregate tolerance relative to the incoming flux
#else
   real(wp), parameter :: RTOL = 1.0e-15_wp
      !! host: measured worst difference is exactly 0 (same statements, same order); 1e-15 leaves one ulp
   real(wp), parameter :: CTOL = 2.0e-15_wp
      !! measured worst closure residual 2.1e-16 of the incoming flux (about 10x)
#endif

   type :: fix_t
      !! Inputs of one radiation run plus the slot.
      integer :: nxp = 0, nyp = 0, ng = 0, nx = 0, ny = 0
      real(wp), allocatable :: aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), tsfcn(:, :, :), sw(:, :)
      integer(int32), allocatable :: ice_err(:, :)
      type(ice_rad_state_t) :: rad
      type(ice_sw_params_t) :: swp
   end type fix_t

contains

   subroutine collect_cice_radiation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("slot_lifecycle_bytes", test_slot), &
                  new_unittest("slot_device_canary_every_array", test_canary), &
                  new_unittest("shim_fails_loud", test_fail_loud), &
                  new_unittest("band_split_sums_to_input", test_band_split), &
                  new_unittest("closure_per_category_each_band", test_closure), &
                  new_unittest("aggregate_is_hand_weighted_mean", test_aggregate), &
                  new_unittest("defaults_and_zero_shortwave", test_defaults), &
                  new_unittest("nan_flags_own_cell_only", test_nan), &
                  new_unittest("device_vs_host_every_output", test_device_vs_host), &
                  new_unittest("oracle_icepack_step_radiation", test_oracle), &
                  new_unittest("driver_throughput", test_throughput) &
                  ]
   end subroutine collect_cice_radiation_tests

   ! ------------------------------------------------------------------ fixture

   subroutine fix_make(f, nxp, nyp, ng, seed, vis_frac)
      !! Random category state and shortwave over the whole grid (ring included, so the ring has
      !! inputs too), slot initialised and its host arrays filled with the sentinels.
      type(fix_t), intent(inout) :: f
      integer, intent(in) :: nxp, nyp, ng
      integer(int64), intent(in) :: seed
      real(wp), intent(in) :: vis_frac
      integer(int64) :: s
      integer :: i, j, n, ierr
      real(wp) :: u, h

      f%nxp = nxp
      f%nyp = nyp
      f%ng = ng
      f%nx = nxp + 2*ng
      f%ny = nyp + 2*ng
      allocate (f%aicen(f%nx, f%ny, NC), f%vicen(f%nx, f%ny, NC), f%vsnon(f%nx, f%ny, NC), f%tsfcn(f%nx, f%ny, NC))
      allocate (f%sw(f%nx, f%ny), f%ice_err(f%nx, f%ny))
      f%ice_err = 0_int32
      s = lcg_next(seed)
      do j = 1, f%ny
         do i = 1, f%nx
            s = lcg_next(s)
            u = lcg_unit(s)
            f%sw(i, j) = merge(0.0_wp, 400.0_wp*lcg_unit(lcg_next(s)), u < 0.08_wp)
            do n = 1, NC
               s = lcg_next(s)
               u = lcg_unit(s)
               s = lcg_next(s)
               if (u < 0.2_wp) then
                  f%aicen(i, j, n) = 0.0_wp
                  f%vicen(i, j, n) = 0.0_wp
                  f%vsnon(i, j, n) = 0.0_wp
               else if (u < 0.25_wp) then
                  f%aicen(i, j, n) = 5.0e-12_wp   ! below puny: counts as empty
                  f%vicen(i, j, n) = 1.0e-12_wp
                  f%vsnon(i, j, n) = 0.0_wp
               else
                  f%aicen(i, j, n) = 0.02_wp + 0.2_wp*lcg_unit(s)
                  s = lcg_next(s)
                  h = 0.05_wp + 5.0_wp*lcg_unit(s)**2
                  f%vicen(i, j, n) = f%aicen(i, j, n)*h
                  s = lcg_next(s)
                  f%vsnon(i, j, n) = merge(0.0_wp, f%aicen(i, j, n)*0.6_wp*lcg_unit(s), lcg_unit(lcg_next(s)) < 0.33_wp)
               end if
               s = lcg_next(s)
               f%tsfcn(i, j, n) = -25.0_wp*lcg_unit(s)*merge(1.0_wp, 0.1_wp, lcg_unit(lcg_next(s)) < 0.7_wp)
            end do
         end do
      end do
      f%swp = ice_sw_params_t(albicev=0.86_wp, albicei=0.44_wp, albsnowv=0.98_wp, albsnowi=0.70_wp, ahmax=0.1_wp, &
                              sw_redist=.true., sw_frac=0.9_wp, sw_dtemp=0.02_wp, sw_vis_frac=vis_frac)
      call f%rad%init(f%nx, f%ny, NC, NL, ierr)
      call f%rad%configure(f%swp%ahmax, SNOWPATCH)
      f%rad%alvn = SENT; f%rad%alin = SENT
      f%rad%fswsfcn = SENT; f%rad%fswintn = SENT; f%rad%fswthrun = SENT
      f%rad%Iswabsn = SENT
      f%rad%alv_ai = SENT; f%rad%ali_ai = SENT
      f%rad%errn = ESENT
   end subroutine fix_make

   subroutine fix_run(f)
      !! Map inputs and slot, run the driver over the interior, bring every output back, unmap.
      type(fix_t), intent(inout) :: f
      integer :: is, ie, js, je, ierr
      is = f%ng + 1
      ie = f%ng + f%nxp
      js = f%ng + 1
      je = f%ng + f%nyp
      !$acc enter data copyin(f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err)
      call f%rad%enter_data()
      call ice_step_radiation(f%rad, f%swp, is, ie, js, je, f%nx, f%ny, NC, NL, &
                              f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err, ierr)
      if (ierr /= ICE_RAD_OK) error stop "fix_run: ice_step_radiation refused"
      !$acc update self(f%rad%alvn, f%rad%alin)
      !$acc update self(f%rad%fswsfcn, f%rad%fswintn, f%rad%fswthrun, f%rad%Iswabsn)
      !$acc update self(f%rad%alv_ai, f%rad%ali_ai, f%rad%errn, f%ice_err)
      call f%rad%exit_data()
      !$acc exit data delete(f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err)
   end subroutine fix_run

   subroutine fix_free(f)
      type(fix_t), intent(inout) :: f
      call f%rad%destroy()
      deallocate (f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err)
   end subroutine fix_free

   subroutine ref_cell(f, i, j, n, alvn, alin, fsfc, fint, fthru, isw)
      !! Independent host reference of one category: the scalar kernels called directly.
      type(fix_t), intent(in) :: f
      integer, intent(in) :: i, j, n
      real(wp), intent(out) :: alvn, alin, fsfc, fint, fthru, isw(NL)
      real(wp) :: a, hi, hs, vi, ii, vs, is_, as, swv, swn
      alvn = ICE_ALBOCN
      alin = ICE_ALBOCN
      fsfc = 0.0_wp
      fint = 0.0_wp
      fthru = 0.0_wp
      isw = 0.0_wp
      a = f%aicen(i, j, n)
      if (a <= ICE_SW_PUNY) return
      hi = f%vicen(i, j, n)/a
      hs = f%vsnon(i, j, n)/a
      call ice_albedos_ccsm3(hi, hs, f%tsfcn(i, j, n), f%swp%albicev, f%swp%albicei, f%swp%albsnowv, f%swp%albsnowi, &
                             ice_albedo_fhtan(f%swp%ahmax), SNOWPATCH, vi, ii, vs, is_, as, alvn, alin)
      call ice_sw_band_split(f%sw(i, j), f%swp%sw_vis_frac, swv, swn)
      call ice_absorbed_ccsm3(NL, hi, as, swv, swn, vi, ii, vs, is_, fsfc, fint, fthru, isw)
   end subroutine ref_cell

   subroutine accum(w, x)
      !! `w = max(w, x)` that cannot hide a NaN: a non-finite `x` poisons `w` with `huge`.
      real(wp), intent(inout) :: w
      real(wp), intent(in) :: x
      if (ieee_is_finite(x)) then
         w = max(w, x)
      else
         w = huge(1.0_wp)
      end if
   end subroutine accum

   logical function close_to(x, ref, tol)
      !! NaN-proof relative comparison with a floor of 1 in the scale.
      real(wp), intent(in) :: x, ref, tol
      close_to = .false.
      if (.not. (ieee_is_finite(x) .and. ieee_is_finite(ref))) return
      close_to = abs(x - ref) <= tol*max(abs(ref), 1.0_wp)
   end function close_to

   ! ------------------------------------------------------------------ slot

   subroutine test_slot(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_rad_state_t) :: r
      integer, parameter :: NX = 9, NY = 6
      integer :: ierr
      integer(int64) :: expect

      call r%init(0, NY, NC, NL, ierr)
      call check(error, ierr, ICE_RAD_INIT_BAD_EXTENT, "nx = 0 refused")
      if (allocated(error)) return
      call r%init(NX, NY, NC, 9, ierr)
      call check(error, ierr, ICE_RAD_INIT_BAD_EXTENT, "nilyr beyond the bound refused")
      if (allocated(error)) return
      call r%init(NX, NY, 9, NL, ierr)
      call check(error, ierr, ICE_RAD_INIT_BAD_EXTENT, "ncat beyond the bound refused")
      if (allocated(error)) return
      call check(error, r%is_init, .false., "refused init allocates nothing")
      if (allocated(error)) return
      call check(error, r%bytes(), 0_int64, "bytes of an uninitialised slot")
      if (allocated(error)) return
      call r%init(NX, NY, NC, NL, ierr)
      call check(error, ierr, ICE_RAD_INIT_OK, "init ok")
      if (allocated(error)) return
      ! 2 albedos + 3 fluxes + nilyr layer arrays per category, 2 aggregates, all doubles; errn int32
      expect = 8_int64*(int(NX*NY*NC*(5 + NL), int64) + int(2*NX*NY, int64)) + 4_int64*int(NX*NY*NC, int64)
      write (output_unit, '(a,i0,a)') "  slot size at 9 x 6 x 5 x 7: ", expect, " bytes"
      call check(error, r%bytes(), expect, "bytes against the closed form")
      if (allocated(error)) return
      if (any(r%alvn /= 0.0_wp) .or. any(r%Iswabsn /= 0.0_wp) .or. any(r%ali_ai /= 0.0_wp) .or. any(r%errn /= 0)) then
         call test_failed(error, "init must zero every array")
         return
      end if
      call r%init(NX + 1, NY, NC, NL, ierr)   ! re-init destroys first
      call check(error, r%bytes() > expect, .true., "re-init reallocates")
      if (allocated(error)) return
      call r%destroy()
      call check(error, r%bytes(), 0_int64, "destroy releases everything")
   end subroutine test_slot

   subroutine test_canary(error)
      !! Host values go up with `copyin`; a kernel adds a per-element canary to EVERY array; every element
      !! (ghost ring included) comes back and is compared.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 9, NY = 6
      type(ice_rad_state_t) :: r
      integer :: ierr, i, j, n, k, bad
      real(wp) :: h0

      call r%init(NX, NY, NC, NL, ierr)
      h0 = 0.5_wp
      r%alvn = h0; r%alin = h0
      r%fswsfcn = h0; r%fswintn = h0; r%fswthrun = h0
      r%Iswabsn = h0
      r%alv_ai = h0; r%ali_ai = h0
      r%errn = 7_int32
      call r%enter_data()
      call canary3(NX, NY, NC, r%alvn, 1)
      call canary3(NX, NY, NC, r%alin, 2)
      call canary3(NX, NY, NC, r%fswsfcn, 5)
      call canary3(NX, NY, NC, r%fswintn, 6)
      call canary3(NX, NY, NC, r%fswthrun, 7)
      call canary4(NX, NY, NC, NL, r%Iswabsn, 8)
      call canary2(NX, NY, r%alv_ai, 10)
      call canary2(NX, NY, r%ali_ai, 12)
      call canary3i(NX, NY, NC, r%errn)
      !$acc update self(r%alvn, r%alin, r%fswsfcn, r%fswintn, r%fswthrun)
      !$acc update self(r%Iswabsn, r%alv_ai, r%ali_ai, r%errn)
      call r%exit_data()
      bad = 0
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               if (r%alvn(i, j, n) /= h0 + c3(1, i, j, n)) bad = bad + 1
               if (r%alin(i, j, n) /= h0 + c3(2, i, j, n)) bad = bad + 1
               if (r%fswsfcn(i, j, n) /= h0 + c3(5, i, j, n)) bad = bad + 1
               if (r%fswintn(i, j, n) /= h0 + c3(6, i, j, n)) bad = bad + 1
               if (r%fswthrun(i, j, n) /= h0 + c3(7, i, j, n)) bad = bad + 1
               if (r%errn(i, j, n) /= 7_int32 + int(i + 10*j + 100*n, int32)) bad = bad + 1
               do k = 1, NL
                  if (r%Iswabsn(i, j, n, k) /= h0 + c4(8, i, j, n, k)) bad = bad + 1
               end do
            end do
         end do
      end do
      do j = 1, NY
         do i = 1, NX
            if (r%alv_ai(i, j) /= h0 + c2(10, i, j)) bad = bad + 1
            if (r%ali_ai(i, j) /= h0 + c2(12, i, j)) bad = bad + 1
         end do
      end do
      call check(error, bad, 0, "every element of every array carries its canary (copyin + device write + update self)")
      call r%destroy()
   end subroutine test_canary

   pure function c2(id, i, j) result(v)
      integer, intent(in) :: id, i, j
      real(wp) :: v
      v = real(id*100000 + i + 10*j, wp)
   end function c2

   pure function c3(id, i, j, n) result(v)
      integer, intent(in) :: id, i, j, n
      real(wp) :: v
      v = real(id*100000 + i + 10*j + 100*n, wp)
   end function c3

   pure function c4(id, i, j, n, k) result(v)
      integer, intent(in) :: id, i, j, n, k
      real(wp) :: v
      v = real(id*100000 + i + 10*j + 100*n + 1000*k, wp)
   end function c4

   subroutine canary2(nx, ny, a, id)
      integer, intent(in) :: nx, ny, id
      real(wp), intent(inout) :: a(nx, ny)
      integer :: i, j
      !$acc data present(a)
      do concurrent(j=1:ny, i=1:nx)
         a(i, j) = a(i, j) + real(id*100000 + i + 10*j, wp)
      end do
      !$acc end data
   end subroutine canary2

   subroutine canary3(nx, ny, nc_, a, id)
      integer, intent(in) :: nx, ny, nc_, id
      real(wp), intent(inout) :: a(nx, ny, nc_)
      integer :: i, j, n
      !$acc data present(a)
      do concurrent(n=1:nc_, j=1:ny, i=1:nx)
         a(i, j, n) = a(i, j, n) + real(id*100000 + i + 10*j + 100*n, wp)
      end do
      !$acc end data
   end subroutine canary3

   subroutine canary4(nx, ny, nc_, nk, a, id)
      integer, intent(in) :: nx, ny, nc_, nk, id
      real(wp), intent(inout) :: a(nx, ny, nc_, nk)
      integer :: i, j, n, k
      !$acc data present(a)
      do concurrent(k=1:nk, n=1:nc_, j=1:ny, i=1:nx)
         a(i, j, n, k) = a(i, j, n, k) + real(id*100000 + i + 10*j + 100*n + 1000*k, wp)
      end do
      !$acc end data
   end subroutine canary4

   subroutine canary3i(nx, ny, nc_, a)
      integer, intent(in) :: nx, ny, nc_
      integer(int32), intent(inout) :: a(nx, ny, nc_)
      integer :: i, j, n
      !$acc data present(a)
      do concurrent(n=1:nc_, j=1:ny, i=1:nx)
         a(i, j, n) = a(i, j, n) + int(i + 10*j + 100*n, int32)
      end do
      !$acc end data
   end subroutine canary3i

   ! ------------------------------------------------------------------ band split

   subroutine test_band_split(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: sw, vf, a, b, worst, e
      integer :: m
      integer(int64) :: s
      s = 99_int64
      worst = 0.0_wp
      do m = 1, 4000
         s = lcg_next(s)
         sw = 1000.0_wp*lcg_unit(s)
         s = lcg_next(s)
         vf = lcg_unit(s)
         call ice_sw_band_split(sw, vf, a, b)
         e = abs((a + b) - sw)/max(sw, tiny(1.0_wp))
         if (.not. ieee_is_finite(e)) then
            call test_failed(error, "non-finite band split")
            return
         end if
         worst = max(worst, e)
         if (.not. (a >= 0.0_wp .and. b >= 0.0_wp)) then
            call test_failed(error, "negative band")
            return
         end if
         if (a /= sw*vf) then
            call test_failed(error, "visible total is not sw_vis_frac*sw")
            return
         end if
      end do
      write (output_unit, '(a,es10.2)') "  band split: worst |swv + swn - sw|/sw = ", worst
      call check(error, worst <= epsilon(1.0_wp), "the two bands add back to the input within one eps")
      if (allocated(error)) return
      call ice_sw_band_split(100.0_wp, 1.0_wp, a, b)
      call check(error, b == 0.0_wp .and. a == 100.0_wp, .true., "all-visible forcing has no near-IR")
      if (allocated(error)) return
      call ice_sw_band_split(100.0_wp, 0.0_wp, a, b)
      call check(error, a == 0.0_wp .and. b == 100.0_wp, .true., "all-near-IR forcing has no visible")
   end subroutine test_band_split

   subroutine test_fail_loud(error)
      !! The host shim refuses (ierr /= 0, nothing runs, outputs untouched) an uninitialised slot, an
      !! unconfigured slot and extents that differ from the slot's.
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      type(ice_rad_state_t) :: dead
      integer :: ierr
      integer(int32) :: e0(4, 3)
      real(wp) :: a(4, 3, NC), sw(4, 3)

      a = 0.1_wp
      sw = 100.0_wp
      e0 = 0_int32
      call ice_step_radiation(dead, f%swp, 1, 4, 1, 3, 4, 3, NC, NL, a, a, a, a, sw, e0, ierr)
      call check(error, ierr, ICE_RAD_ERR_NOT_INIT, "uninitialised slot refused")
      if (allocated(error)) return
      call fix_make(f, 5, 4, 1, 3_int64, 0.52_wp)
      f%rad%configured = .false.
      call ice_step_radiation(f%rad, f%swp, 2, 6, 2, 5, f%nx, f%ny, NC, NL, f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, &
                              f%ice_err, ierr)
      call check(error, ierr, ICE_RAD_ERR_NOT_CONFIGURED, "unconfigured slot refused")
      if (allocated(error)) return
      call f%rad%configure(f%swp%ahmax, SNOWPATCH)
      block
         real(wp), allocatable :: b(:, :, :), w(:, :)
         integer(int32), allocatable :: ee(:, :)
         allocate (b(f%nx + 1, f%ny, NC), w(f%nx + 1, f%ny), ee(f%nx + 1, f%ny))
         b = 0.1_wp; w = 1.0_wp; ee = 0_int32
         call ice_step_radiation(f%rad, f%swp, 2, 6, 2, 5, f%nx + 1, f%ny, NC, NL, b, b, b, b, w, ee, ierr)
         call check(error, ierr, ICE_RAD_ERR_SHAPE, "nx different from the slot refused")
         if (allocated(error)) return
      end block
      call ice_step_radiation(f%rad, f%swp, 2, 6, 2, 5, f%nx, f%ny, NC, NL - 1, f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, &
                              f%ice_err, ierr)
      call check(error, ierr, ICE_RAD_ERR_SHAPE, "nilyr different from the slot refused")
      if (allocated(error)) return
      call ice_step_radiation(f%rad, f%swp, 2, 9, 2, 5, f%nx, f%ny, NC, NL, f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, &
                              f%ice_err, ierr)
      call check(error, ierr, ICE_RAD_ERR_SHAPE, "cell range outside the grid refused")
      if (allocated(error)) return
      call check(error, all(f%rad%alvn == SENT) .and. all(f%rad%errn == ESENT), .true., "a refused call leaves the slot untouched")
      call fix_free(f)
   end subroutine test_fail_loud

   ! ------------------------------------------------------------------ closure

   subroutine test_closure(error)
      !! `sw = reflected + fswsfc + fswint + fswthru` per category, with `reflected` the band-weighted albedo
      !! reflection, and `sum_k Iswabs = fswint`; plus the band statements: near-IR only is absorbed at the surface
      !! (`fswint = fswthru = 0`, `fswsfc = sw (1 - alin)`); visible only splits into the three.
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      real(wp), parameter :: VF(3) = [0.52_wp, 1.0_wp, 0.0_wp]
      real(wp) :: swv, swn, refl, resid, wl, wb, w_close, w_layer, w_band, sc
      integer :: iv, i, j, n, nice

      w_close = 0.0_wp
      w_layer = 0.0_wp
      w_band = 0.0_wp
      nice = 0
      do iv = 1, 3
         call fix_make(f, 6, 5, 2, 31_int64 + iv, VF(iv))
         call fix_run(f)
         do n = 1, NC
            do j = f%ng + 1, f%ng + f%nyp
               do i = f%ng + 1, f%ng + f%nxp
                  if (f%aicen(i, j, n) <= ICE_SW_PUNY) cycle
                  nice = nice + 1
                  call ice_sw_band_split(f%sw(i, j), VF(iv), swv, swn)
                  refl = swv*f%rad%alvn(i, j, n) + swn*f%rad%alin(i, j, n)
                  sc = max(f%sw(i, j), 1.0_wp)
                  resid = f%sw(i, j) - refl - f%rad%fswsfcn(i, j, n) - f%rad%fswintn(i, j, n) - f%rad%fswthrun(i, j, n)
                  if (.not. ieee_is_finite(resid)) then
                     call test_failed(error, "non-finite closure residual")
                     return
                  end if
                  w_close = max(w_close, abs(resid)/sc)  ! resid checked finite above
                  wl = abs(sum(f%rad%Iswabsn(i, j, n, :)) - f%rad%fswintn(i, j, n))/sc
                  if (.not. ieee_is_finite(wl)) then
                     call test_failed(error, "non-finite layer residual")
                     return
                  end if
                  w_layer = max(w_layer, wl)
                  if (VF(iv) == 0.0_wp) then
                     if (f%rad%fswintn(i, j, n) /= 0.0_wp .or. f%rad%fswthrun(i, j, n) /= 0.0_wp) then
                        call test_failed(error, "near-IR-only forcing must not penetrate")
                        return
                     end if
                     wb = abs(f%rad%fswsfcn(i, j, n) - f%sw(i, j)*(1.0_wp - f%rad%alin(i, j, n)))/sc
                     if (.not. ieee_is_finite(wb)) then
                        call test_failed(error, "non-finite band residual")
                        return
                     end if
                     w_band = max(w_band, wb)
                  end if
                  if (VF(iv) == 1.0_wp) then
                     wb = abs(f%rad%fswsfcn(i, j, n) + f%rad%fswintn(i, j, n) + f%rad%fswthrun(i, j, n) &
                              - f%sw(i, j)*(1.0_wp - f%rad%alvn(i, j, n)))/sc
                     if (.not. ieee_is_finite(wb)) then
                        call test_failed(error, "non-finite band residual")
                        return
                     end if
                     w_band = max(w_band, wb)
                  end if
                  if (f%rad%fswsfcn(i, j, n) < 0.0_wp .or. f%rad%fswintn(i, j, n) < 0.0_wp &
                      .or. f%rad%fswthrun(i, j, n) < 0.0_wp) then
                     call test_failed(error, "negative absorbed flux")
                     return
                  end if
               end do
            end do
         end do
         call fix_free(f)
      end do
   write (output_unit, '(a,i0,a,3es10.2)') "  closure over ", nice, " ice categories: worst residual/sw (total, layers, band) = ", &
         w_close, w_layer, w_band
      call check(error, nice > 200, .true., "enough ice categories exercised")
      if (allocated(error)) return
      call check(error, w_close <= CTOL, .true., "sw = reflected + fswsfc + fswint + fswthru")
      if (allocated(error)) return
      call check(error, w_layer <= CTOL, .true., "sum_k Iswabs = fswint")
      if (allocated(error)) return
      call check(error, w_band <= CTOL, .true., "single-band closure")
   end subroutine test_closure

   ! ------------------------------------------------------------------ aggregate

   subroutine test_aggregate(error)
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      real(wp) :: v1, v3, w1, w3, worst
      integer :: i, j, n, nempty, ncell

      call fix_make(f, 8, 6, 2, 77_int64, 0.52_wp)
      ! force three whole cells ice-free
      f%aicen(3, 3, :) = 0.0_wp
      f%aicen(4, 4, :) = 0.0_wp
      f%aicen(5, 5, :) = 3.0e-12_wp
      call fix_run(f)
      worst = 0.0_wp
      nempty = 0
      ncell = 0
      do j = f%ng + 1, f%ng + f%nyp
         do i = f%ng + 1, f%ng + f%nxp
            ncell = ncell + 1
            v1 = 0.0_wp; v3 = 0.0_wp
            do n = 1, NC
               if (f%aicen(i, j, n) > ICE_SW_PUNY) then
                  v1 = v1 + f%rad%alvn(i, j, n)*f%aicen(i, j, n)
                  v3 = v3 + f%rad%alin(i, j, n)*f%aicen(i, j, n)
               end if
            end do
            if (all(f%aicen(i, j, :) <= ICE_SW_PUNY)) then
               nempty = nempty + 1
               if (f%rad%alv_ai(i, j) /= 0.0_wp .or. f%rad%ali_ai(i, j) /= 0.0_wp) then
                  call test_failed(error, "ice-free cell must aggregate to exactly 0")
                  return
               end if
            end if
            w1 = abs(f%rad%alv_ai(i, j) - v1)
            w3 = abs(f%rad%ali_ai(i, j) - v3)
            if (.not. (ieee_is_finite(w1) .and. ieee_is_finite(w3))) then
               call test_failed(error, "non-finite aggregate")
               return
            end if
            worst = max(worst, w1, w3)
            if (f%rad%alv_ai(i, j) > 1.0_wp .or. f%rad%alv_ai(i, j) < 0.0_wp) then
               call test_failed(error, "aggregate albedo outside [0,1]")
               return
            end if
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,es10.2)') "  aggregate over ", ncell, " cells (", nempty, " ice-free): worst |driver - hand| = ", worst
      call check(error, nempty >= 3, .true., "ice-free cells exercised")
      if (allocated(error)) return
      call check(error, worst <= RTOL, .true., "aggregate equals the hand area-weighted sum")
      call fix_free(f)
   end subroutine test_aggregate

   ! ------------------------------------------------------------------ defaults and zero SW

   subroutine test_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      integer :: i, j, n, k, nemp, bad
      real(wp) :: rv, ri, rf, rn, rt, ris(NL)

      call fix_make(f, 7, 5, 1, 5_int64, 0.52_wp)
      call fix_run(f)
      nemp = 0
      bad = 0
      do n = 1, NC
         do j = f%ng + 1, f%ng + f%nyp
            do i = f%ng + 1, f%ng + f%nxp
               if (f%aicen(i, j, n) > ICE_SW_PUNY) cycle
               nemp = nemp + 1
               if (f%rad%alvn(i, j, n) /= ICE_ALBOCN .or. f%rad%alin(i, j, n) /= ICE_ALBOCN) bad = bad + 1
               if (f%rad%fswsfcn(i, j, n) /= 0.0_wp .or. f%rad%fswintn(i, j, n) /= 0.0_wp &
                   .or. f%rad%fswthrun(i, j, n) /= 0.0_wp) bad = bad + 1
               do k = 1, NL
                  if (f%rad%Iswabsn(i, j, n, k) /= 0.0_wp) bad = bad + 1
               end do
               if (f%rad%errn(i, j, n) /= 0_int32) bad = bad + 1
            end do
         end do
      end do
      call check(error, nemp > 20, .true., "empty categories exercised")
      if (allocated(error)) return
      call check(error, bad, 0, "empty categories: albedo ICE_ALBOCN, zero fluxes, no error")
      if (allocated(error)) return
      call fix_free(f)

      ! zero incoming shortwave: every absorbed flux exactly zero, albedos still those of the state
      call fix_make(f, 7, 5, 1, 6_int64, 0.52_wp)
      f%sw = 0.0_wp
      call fix_run(f)
      bad = 0
      do n = 1, NC
         do j = f%ng + 1, f%ng + f%nyp
            do i = f%ng + 1, f%ng + f%nxp
               if (f%rad%fswsfcn(i, j, n) /= 0.0_wp .or. f%rad%fswintn(i, j, n) /= 0.0_wp &
                   .or. f%rad%fswthrun(i, j, n) /= 0.0_wp) bad = bad + 1
               do k = 1, NL
                  if (f%rad%Iswabsn(i, j, n, k) /= 0.0_wp) bad = bad + 1
               end do
               call ref_cell(f, i, j, n, rv, ri, rf, rn, rt, ris)
               if (.not. close_to(f%rad%alvn(i, j, n), rv, RTOL) .or. .not. close_to(f%rad%alin(i, j, n), ri, RTOL)) &
                  bad = bad + 1
            end do
         end do
      end do
      call check(error, bad, 0, "zero shortwave: all absorbed fluxes exactly zero, albedos unchanged")
      call fix_free(f)
   end subroutine test_defaults

   ! ------------------------------------------------------------------ NaN

   subroutine test_nan(error)
      !! Four poisoned inputs, each in its own cell: aicen (category 2), vicen (3), Tsfcn (1), sw (all categories
      !! of the cell).  Only those cells are flagged, only those categories get the word, their outputs are the
      !! finite defaults, and every other cell is bit-identical to the clean run.
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      real(wp), allocatable :: base_f(:, :, :), base_a(:, :)
      integer(int32), allocatable :: base_e(:, :, :)
      real(wp) :: nan
      integer :: i, j, n, bad, nflag
      logical :: p_cell

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call fix_make(f, 7, 6, 2, 8_int64, 0.52_wp)
      call fix_run(f)
      allocate (base_f(f%nx, f%ny, NC), base_a(f%nx, f%ny), base_e(f%nx, f%ny, NC))
      base_f = f%rad%fswsfcn
      base_a = f%rad%alv_ai
      base_e = f%rad%errn
      call check(error, all(f%ice_err(f%ng + 1:f%ng + f%nxp, f%ng + 1:f%ng + f%nyp) == 0_int32), .true., "clean run flags nothing")
      if (allocated(error)) return
      f%aicen(4, 4, 2) = nan
      f%vicen(5, 3, 3) = nan
      f%tsfcn(6, 5, 1) = nan
      f%sw(7, 6) = nan
      f%rad%fswsfcn = SENT
      f%rad%alv_ai = SENT
      f%rad%errn = ESENT
      call fix_run(f)
      bad = 0
      nflag = 0
      do j = f%ng + 1, f%ng + f%nyp
         do i = f%ng + 1, f%ng + f%nxp
            p_cell = (i == 4 .and. j == 4) .or. (i == 5 .and. j == 3) .or. (i == 6 .and. j == 5) .or. (i == 7 .and. j == 6)
            if (p_cell) then
               if (iand(f%ice_err(i, j), ICE_ERR_NONFINITE) == 0_int32) bad = bad + 1
               nflag = nflag + 1
            else
               if (f%ice_err(i, j) /= 0_int32) bad = bad + 1
               if (f%rad%alv_ai(i, j) /= base_a(i, j)) bad = bad + 1
            end if
            do n = 1, NC
               if (p_cell) then
                  if (.not. (ieee_is_finite(f%rad%fswsfcn(i, j, n)) .and. ieee_is_finite(f%rad%alvn(i, j, n)))) bad = bad + 1
                  if ((i == 4 .and. j == 4 .and. n == 2) .or. (i == 5 .and. j == 3 .and. n == 3) &
                      .or. (i == 6 .and. j == 5 .and. n == 1) .or. (i == 7 .and. j == 6)) then
                     if (f%rad%errn(i, j, n) /= ICE_ERR_NONFINITE) bad = bad + 1
                     if (f%rad%fswsfcn(i, j, n) /= 0.0_wp .or. f%rad%alvn(i, j, n) /= ICE_ALBOCN) bad = bad + 1
                  else
                     if (f%rad%errn(i, j, n) /= 0_int32) bad = bad + 1
                  end if
               else
                  if (f%rad%errn(i, j, n) /= base_e(i, j, n)) bad = bad + 1
                  if (f%rad%fswsfcn(i, j, n) /= base_f(i, j, n)) bad = bad + 1
               end if
            end do
         end do
      end do
      call check(error, nflag, 4, "four poisoned cells")
      if (allocated(error)) return
      call check(error, bad, 0, "NaN flags its own category/cell only; the rest is bit-identical to the clean run")
      call fix_free(f)
   end subroutine test_nan

   ! ------------------------------------------------------------------ device vs host

   subroutine test_device_vs_host(error)
      !! Interior 9 x 7 of an 13 x 11 grid, ncat 5, nilyr 7: EVERY output of every (i, j, n[, k]) against the
      !! independent host reference; the ring keeps its sentinels; the interior is fully overwritten.
      type(error_type), allocatable, intent(out) :: error
      type(fix_t) :: f
      real(wp) :: av, ai, fsfc, fint, fthru, isw(NL), v1, v3, worst, x
      integer :: i, j, n, k, bad, nint_
      logical :: inside

      call fix_make(f, 9, 7, 2, 1234_int64, 0.52_wp)
      call fix_run(f)
      bad = 0
      nint_ = 0
      worst = 0.0_wp
      do n = 1, NC
         do j = 1, f%ny
            do i = 1, f%nx
               inside = i > f%ng .and. i <= f%ng + f%nxp .and. j > f%ng .and. j <= f%ng + f%nyp
               if (.not. inside) then
                  if (f%rad%alvn(i, j, n) /= SENT .or. f%rad%alin(i, j, n) /= SENT .or. f%rad%fswsfcn(i, j, n) /= SENT &
                      .or. f%rad%fswintn(i, j, n) /= SENT .or. f%rad%fswthrun(i, j, n) /= SENT &
                      .or. f%rad%errn(i, j, n) /= ESENT) bad = bad + 1
                  if (any(f%rad%Iswabsn(i, j, n, :) /= SENT)) bad = bad + 1
                  if (n == 1) then
                     if (f%rad%alv_ai(i, j) /= SENT .or. f%rad%ali_ai(i, j) /= SENT .or. f%ice_err(i, j) /= 0_int32) &
                        bad = bad + 1
                  end if
                  cycle
               end if
               nint_ = nint_ + 1
               call ref_cell(f, i, j, n, av, ai, fsfc, fint, fthru, isw)
               if (.not. (close_to(f%rad%alvn(i, j, n), av, RTOL) .and. close_to(f%rad%alin(i, j, n), ai, RTOL) &
                          .and. close_to(f%rad%fswsfcn(i, j, n), fsfc, RTOL) .and. close_to(f%rad%fswintn(i, j, n), fint, RTOL) &
                          .and. close_to(f%rad%fswthrun(i, j, n), fthru, RTOL))) bad = bad + 1
               call accum(worst, abs(f%rad%alvn(i, j, n) - av))
               call accum(worst, abs(f%rad%alin(i, j, n) - ai))
               call accum(worst, abs(f%rad%fswsfcn(i, j, n) - fsfc)/max(abs(fsfc), 1.0_wp))
               call accum(worst, abs(f%rad%fswintn(i, j, n) - fint)/max(abs(fint), 1.0_wp))
               call accum(worst, abs(f%rad%fswthrun(i, j, n) - fthru)/max(abs(fthru), 1.0_wp))
               do k = 1, NL
                  if (.not. close_to(f%rad%Iswabsn(i, j, n, k), isw(k), RTOL)) bad = bad + 1
                  x = abs(f%rad%Iswabsn(i, j, n, k) - isw(k))/max(abs(isw(k)), 1.0_wp)
                  call accum(worst, x)
               end do
               if (f%rad%errn(i, j, n) /= 0_int32) bad = bad + 1
               if (n == NC) then
                  v1 = 0.0_wp; v3 = 0.0_wp
                  do k = 1, NC
                     if (f%aicen(i, j, k) > ICE_SW_PUNY) then
                        call ref_cell(f, i, j, k, av, ai, fsfc, fint, fthru, isw)
                        v1 = v1 + av*f%aicen(i, j, k)
                        v3 = v3 + ai*f%aicen(i, j, k)
                     end if
                  end do
                  if (.not. (close_to(f%rad%alv_ai(i, j), v1, RTOL) .and. close_to(f%rad%ali_ai(i, j), v3, RTOL) &
                             .and. f%ice_err(i, j) == 0_int32)) bad = bad + 1
               end if
            end do
         end do
      end do
      write (output_unit, '(a,i0,a,es10.2)') "  device vs host over ", nint_, " cell-categories: worst scaled difference ", worst
      call check(error, bad, 0, "every output of every interior cell-category matches the host reference; ring untouched")
      call fix_free(f)
   end subroutine test_device_vs_host

   ! ------------------------------------------------------------------ oracle

   subroutine load_gold(gi, gc)
      real(wp), intent(out) :: gi(88, 10), gc(88, 10)
      real(wp) :: GOLD_ICEPACK(88, 10), GOLD_CESM(88, 10)
#include "test_cice_radiation_golden.inc"
      gi = GOLD_ICEPACK
      gc = GOLD_CESM
   end subroutine load_gold

   subroutine test_oracle(error)
      !! Ten multi-category cells (ncat 5, nilyr 7, OM3 and code-default albedo sets, empty and below-puny
      !! categories, melting, thin and thick ice, snow, night) through the FULL driver against upstream
      !! `icepack_step_radiation` + the CICE area-weighted aggregate, for the upstream and the CESM-coupled
      !! `snowpatch` (0.02 / 0.005).  NOT bitwise on gfortran, for two reasons: the port forms the band
      !! totals as `sw*0.52` and `sw - sw*0.52` where upstream sums per-band products (`sw*0.28 + sw*0.24`, ...),
      !! so the absorbed fluxes differ by an ulp; and the reproducible `atan`/`exp` library bodies (<= 2.5 ulp)
      !! differ from libm.  Measured worst scaled difference 2.6e-15 (CPU), 1.8e-15 (GPU); tolerance 3e-14.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: gi(88, 10), gc(88, 10), g(88), got(88), worst, x
      type(fix_t) :: f
      integer :: set, c, m, nbad, n, k
      real(wp) :: sp

      call load_gold(gi, gc)
      nbad = 0
      worst = 0.0_wp
      do set = 1, 2
         do c = 1, 10
            g = merge(gi(:, c), gc(:, c), set == 1)
            sp = merge(0.02_wp, 0.005_wp, set == 1)
            call fix_make(f, 1, 1, 0, 1_int64, 0.52_wp)
            f%swp%albicev = g(1); f%swp%albicei = g(2); f%swp%albsnowv = g(3); f%swp%albsnowi = g(4); f%swp%ahmax = g(5)
            f%sw(1, 1) = g(6)
            f%aicen(1, 1, :) = g(7:11)
            f%vicen(1, 1, :) = g(12:16)
            f%vsnon(1, 1, :) = g(17:21)
            f%tsfcn(1, 1, :) = g(22:26)
            call f%rad%configure(f%swp%ahmax, sp)
            call fix_run(f)
            do n = 1, NC
               got(27 + n - 1) = f%rad%alvn(1, 1, n)
               got(32 + n - 1) = f%rad%alin(1, 1, n)
               got(37 + n - 1) = f%rad%fswsfcn(1, 1, n)
               got(42 + n - 1) = f%rad%fswintn(1, 1, n)
               got(47 + n - 1) = f%rad%fswthrun(1, 1, n)
               do k = 1, NL
                  got(52 + (n - 1)*NL + k - 1) = f%rad%Iswabsn(1, 1, n, k)
               end do
            end do
            got(87) = f%rad%alv_ai(1, 1)
            got(88) = f%rad%ali_ai(1, 1)
            do m = 27, 88
               if (.not. close_to(got(m), g(m), 3.0e-14_wp)) then
                  nbad = nbad + 1
               end if
               x = abs(got(m) - g(m))/max(abs(g(m)), 1.0_wp)
               call accum(worst, x)
            end do
            call fix_free(f)
         end do
      end do
  write (output_unit, '(a,es10.2,a,i0)') "  driver vs upstream step_radiation: worst scaled difference ", worst, ", failures ", nbad
      call check(error, nbad, 0, "driver disagrees with icepack_step_radiation")
   end subroutine test_oracle

   ! ------------------------------------------------------------------ throughput

   subroutine test_throughput(error)
      !! 512 x 512 x 5 categories x 7 layers; sampled host comparison; prints cell-categories/s.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NREP = 10
      type(fix_t) :: f
      integer(int64) :: t0, t1, rate
      real(wp) :: secs, av, ai, fsfc, fint, fthru, isw(NL)
      integer :: k, i, j, n, r, bad, ierr

      call fix_make(f, 512, 512, 0, 4242_int64, 0.52_wp)
      !$acc enter data copyin(f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err)
      call f%rad%enter_data()
      call ice_step_radiation(f%rad, f%swp, 1, f%nx, 1, f%ny, f%nx, f%ny, NC, NL, &
                              f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err, ierr)
      if (ierr /= ICE_RAD_OK) error stop "throughput: refused"
      call system_clock(t0, rate)
      do r = 1, NREP
         call ice_step_radiation(f%rad, f%swp, 1, f%nx, 1, f%ny, f%nx, f%ny, NC, NL, &
                                 f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err, ierr)
      end do
      !$acc update self(f%rad%fswsfcn, f%rad%fswintn, f%rad%fswthrun, f%rad%alvn, f%rad%alin, f%rad%Iswabsn)
      call system_clock(t1)
      call f%rad%exit_data()
      !$acc exit data delete(f%aicen, f%vicen, f%vsnon, f%tsfcn, f%sw, f%ice_err)
      secs = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,f10.4,a,es11.3,a)') "radiation driver (step + aggregate): ", secs/real(NREP, wp)*1.0e3_wp, &
         " ms/step, ", real(f%nx*f%ny*NC, wp)*real(NREP, wp)/secs, " cell-categories/s"
      bad = 0
      do k = 1, f%nx*f%ny*NC, 997
         n = (k - 1)/(f%nx*f%ny) + 1
         j = (k - 1 - (n - 1)*f%nx*f%ny)/f%nx + 1
         i = k - (n - 1)*f%nx*f%ny - (j - 1)*f%nx
         call ref_cell(f, i, j, n, av, ai, fsfc, fint, fthru, isw)
         if (.not. (close_to(f%rad%fswsfcn(i, j, n), fsfc, RTOL) .and. close_to(f%rad%fswthrun(i, j, n), fthru, RTOL) &
                    .and. close_to(f%rad%alvn(i, j, n), av, RTOL) .and. close_to(f%rad%Iswabsn(i, j, n, NL), isw(NL), RTOL))) &
            bad = bad + 1
      end do
      call check(error, bad, 0, "sampled cell-categories match the host reference")
      call fix_free(f)
   end subroutine test_throughput

end module test_cice_radiation
