!! Streaming tracer remap (`ocean_remap_tracer_field_stream`) vs the
!! stack-based reference (`ocean_remap_tracer_field`), both in CONCENTRATION
!! form (White & Adcroft 2008): `c_old` -> PPM on `c` -> overlap integral ->
!! bounded `c_new` -> `hTr = c_new*h_new`, with no pool merge and no defect
!! fold.
!!
!! The stream is the DEFAULT kernel (`&vcoord_nml remap_kernel = "stream"`).
!! On a field of adversarial columns (uniform, 1e3-stretched, bed / top /
!! isolated / consecutive filler runs, source filler under a live target and
!! vice versa, sawtooth / step / spike / constant profiles, dead columns,
!! exact-zero target layers, a column total that differs by a few ulp, coarse
!! and fine targets) at `nz in {3, 4, 7, 50, 128}` it is checked that
!!   * stream and reference agree (hTr, budget, c_new, bounding amount) to
!!     `STREAM_REF_ULPS` ulp of the column scale, with and without the
!!     optional `budget` -- BITWISE with FMA contraction off for the ALE remap
!!     files, a few ulp without (see the constant);
!!   * on every column, both kernels (a) conserve the column content to
!!     round-off once the RECORDED bounding amount is added back, (b) keep
!!     `c_new` inside the range of the source column, and (c) give an exactly
!!     empty target (`h_new = 0`) the mean of the source cell containing its
!!     bottom interface.
!! Fillers enter with invariant I1' (`hTr = h*c_donor`), which is what the tail
!! `enforce_vanished_content` guarantees in production.
!!
!! Columns 27-33 are built from dyadic thicknesses so that source and target
!! interfaces coincide EXACTLY (`z_old(ko) == z_new(k)`): identical grids,
!! targets spanning many source cells, many targets inside one source cell,
!! zero-thickness source / target cells interleaved with live ones, and
!! random dyadic grids with zeros.
!!
!! `nz = 200` exceeds `NZ_STACK_MAX`, where the reference cannot run: the
!! stream alone is checked against (a)-(c).  `nz = 1, 2` and every
!! non-default option must dispatch to the reference (`remap_use_stream`).
!!
!! Device residency (`mem:separate`): every array a kernel touches, scratch
!! included, is mapped here and the results are pulled back with
!! `update self` on the arrays (never on a derived type).
module test_ocean_remap_stream
   use rdb_constants, only: wp, REMAP_PPM, REMAP_PLM, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_remap, only: ocean_remap_tracer_field, ocean_remap_tracer_field_stream, &
                              remap_use_stream
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_remap_stream_tests

   integer, parameter :: NC = 33
      !! Columns (= nx) in the adversarial field.
   real(wp), parameter :: H_TOT = 100.0_wp
   real(wp), parameter :: STREAM_REF_ULPS = 16.0_wp
      !! Allowed stream-vs-reference disagreement, in units of
      !! `epsilon * (|hTr| scale)`.  The two kernels evaluate textually
      !! identical expressions in different files and loop shapes; a compiler
      !! that contracts `a*b + c` into an FMA decides per expression context,
      !! so they differ by a few ulp on a default (`-march=native`) build and
      !! agree BITWISE when FMA contraction is off for the ALE remap files
      !! (`-ffp-contract=off` / `-Mnofma`, measured: 0 differing cells at
      !! nz = 3, 4, 7, 50, 128 with the flag).  Set to `0` once the build
      !! pins that.

contains

   subroutine collect_ocean_remap_stream_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("stream_matches_reference_nz3", test_nz3), &
                  new_unittest("stream_matches_reference_nz4", test_nz4), &
                  new_unittest("stream_matches_reference_nz7", test_nz7), &
                  new_unittest("stream_matches_reference_nz50", test_nz50), &
                  new_unittest("stream_matches_reference_nz128", test_nz128), &
                  new_unittest("stream_nz200_conserves_bounded_fallback", test_nz200), &
                  new_unittest("dispatch_falls_back_to_reference", test_dispatch) &
                  ]
   end subroutine collect_ocean_remap_stream_tests

   subroutine test_nz3(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_kernels(error, 3)
   end subroutine test_nz3

   subroutine test_nz4(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_kernels(error, 4)
   end subroutine test_nz4

   subroutine test_nz7(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_kernels(error, 7)
   end subroutine test_nz7

   subroutine test_nz50(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_kernels(error, 50)
   end subroutine test_nz50

   subroutine test_nz128(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_kernels(error, 128)
   end subroutine test_nz128

   ! ---------------------------------------------------------------- helpers

   real(wp) function lcg(state)
      !! Deterministic uniform (0,1) generator (no dependence on the
      !! compiler's `random_number`).
      integer, intent(inout) :: state
      state = int(mod(1103515245_8*int(state, 8) + 12345_8, 2147483648_8))
      lcg = (real(state, wp) + 0.5_wp)/2147483648.0_wp
   end function lcg

   subroutine layer_thickness(nz, kind, filler, total, h)
      !! `filler(k)` layers get a sub-threshold thickness; the live layers
      !! share `total` (uniform, geometric bed-thin, geometric bed-thick or
      !! alternating coarse/fine).
      integer, intent(in) :: nz, kind
      logical, intent(in) :: filler(nz)
      real(wp), intent(in) :: total
      real(wp), intent(out) :: h(nz)
      real(wp) :: w(nz), ws, r
      integer :: k, nlive
      nlive = count(.not. filler)
      r = 1.0_wp
      if (nz > 1) r = 1.0e3_wp**(1.0_wp/real(max(nz - 1, 1), wp))
      do k = 1, nz
         select case (kind)
         case (1)
            w(k) = r**(k - 1)          ! thin bed, thick top
         case (2)
            w(k) = r**(nz - k)         ! thick bed, thin top
         case (3)
            if (mod(k, 2) == 0) then
               w(k) = 0.2_wp
            else
               w(k) = 1.8_wp
            end if
         case default
            w(k) = 1.0_wp
         end select
      end do
      ws = 0.0_wp
      do k = 1, nz
         if (.not. filler(k)) ws = ws + w(k)
      end do
      do k = 1, nz
         if (filler(k)) then
            h(k) = 5.0e-5_wp
         else
            h(k) = total*w(k)/ws
         end if
      end do
      if (nlive == 0) h = 5.0e-5_wp
   end subroutine layer_thickness

   subroutine build_column(ic, nz, h_old, h_new, q)
      !! Column `ic` of the adversarial field.
      integer, intent(in) :: ic, nz
      real(wp), intent(out) :: h_old(nz), h_new(nz), q(nz)
      logical :: fo(nz), fn(nz), direct
      integer :: k, kind_o, kind_n, rs, run
      real(wp) :: c(nz), tot_n

      fo = .false.
      fn = .false.
      direct = .false.
      kind_o = 0
      kind_n = 0
      tot_n = H_TOT
      do k = 1, nz
         c(k) = 10.0_wp + 5.0_wp*sin(0.7_wp*real(k, wp))
      end do
      select case (ic)
      case (1)   ! uniform / uniform
         continue
      case (2)   ! stretched old (thin bed) / uniform new
         kind_o = 1
      case (3)   ! uniform old / stretched new
         kind_n = 1
      case (4)   ! bed filler run of 1
         fo(1) = .true.
         fn(1) = .true.
      case (5)   ! bed filler run of 3
         run = min(3, nz - 2)
         fo(1:run) = .true.
         fn(1:run) = .true.
      case (6)   ! bed filler run of many
         run = max(1, nz/2)
         fo(1:run) = .true.
         fn(1:run) = .true.
      case (7)   ! top filler run (ice shelf)
         run = min(2, nz - 2)
         fo(nz - run + 1:nz) = .true.
         fn(nz - run + 1:nz) = .true.
      case (8)   ! isolated mid filler
         fo(nz/2 + 1) = .true.
         fn(nz/2 + 1) = .true.
      case (9)   ! consecutive mid fillers
         if (nz >= 5) then
            fo(2:min(4, nz - 1)) = .true.
            fn(2:min(4, nz - 1)) = .true.
         else
            fo(2) = .true.
            fn(2) = .true.
         end if
      case (10)  ! source filler, target all live
         fo(1) = .true.
      case (11)  ! source live, target bed filler run
         run = min(3, nz - 2)
         fn(1:run) = .true.
      case (12)  ! source mid filler, target top filler
         fo(nz/2) = .true.
         fn(nz) = .true.
      case (13)  ! sawtooth, uniform grids
         do k = 1, nz
            if (mod(k, 2) == 0) then
               c(k) = 20.0_wp
            else
               c(k) = 1.0_wp
            end if
         end do
         kind_n = 3
      case (14)  ! step
         do k = 1, nz
            if (k > nz/2) then
               c(k) = 30.0_wp
            else
               c(k) = 5.0_wp
            end if
         end do
         kind_o = 2
      case (15)  ! single spike
         c = 2.0_wp
         c(nz/2 + 1) = 50.0_wp
         kind_n = 1
      case (16)  ! constant
         c = 7.25_wp
         kind_o = 1
         kind_n = 2
      case (17)  ! dead source column
         fo = .true.
      case (18)  ! target with exact zeros (even layers)
         do k = 2, nz, 2
            fn(k) = .true.
         end do
      case (19)  ! column totals differ by a few ulp, fillers present
         fo(1) = .true.
         fn(1) = .true.
         tot_n = H_TOT*(1.0_wp + 3.0_wp*epsilon(1.0_wp))
      case (20)  ! few ulp, no fillers (fast path, no fold)
         tot_n = H_TOT*(1.0_wp - 3.0_wp*epsilon(1.0_wp))
         kind_o = 1
      case (21)  ! coarse / fine alternation target
         kind_n = 3
      case (22)  ! thick-bed source, uniform target (in-place hazard)
         kind_o = 2
      case (23)  ! dead target column
         fn = .true.
      case (24)  ! both dead
         fo = .true.
         fn = .true.
      case (27:33)  ! exact-tie / zero-thickness / many-to-one grids, set directly
         direct = .true.
      case default  ! 25, 26: random thickness, random fillers, random content
         rs = 977*ic
         do k = 1, nz
            fo(k) = lcg(rs) < 0.2_wp
            fn(k) = lcg(rs) < 0.2_wp
         end do
         fo(nz/2 + 1) = .false.
         fn(nz/2 + 1) = .false.
      end select

      if (direct) then
         call direct_grids(ic, nz, h_old, h_new)
         call content_with_i1p(nz, h_old, c, q)
         return
      end if
      call layer_thickness(nz, kind_o, fo, H_TOT, h_old)
      call layer_thickness(nz, kind_n, fn, tot_n, h_new)
      if (ic == 18) then
         do k = 2, nz, 2
            h_new(k) = 0.0_wp
         end do
         call rescale_live(nz, h_new, tot_n)
      end if
      if (ic >= 25) then
         rs = 31*ic
         do k = 1, nz
            if (.not. fo(k)) h_old(k) = h_old(k)*(0.2_wp + 1.6_wp*lcg(rs))
            if (.not. fn(k)) h_new(k) = h_new(k)*(0.2_wp + 1.6_wp*lcg(rs))
            c(k) = 4.0_wp*(lcg(rs) - 0.5_wp)
         end do
         call rescale_live(nz, h_old, H_TOT)
         call rescale_live(nz, h_new, H_TOT)
      end if
      call content_with_i1p(nz, h_old, c, q)
   end subroutine build_column

   subroutine donor_conc(nz, h, c_live, c_don)
      !! Concentration every layer reads: its own on a live layer, the nearest
      !! live layer ABOVE's on a filler (the topmost live layer's for the run
      !! above the water), `0` on a column with no live layer.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(nz), c_live(nz)
      real(wp), intent(out) :: c_don(nz)
      integer :: k, kt
      real(wp) :: cd
      kt = 0
      do k = nz, 1, -1
         if (h(k) > H_VANISHED) then
            kt = k
            exit
         end if
      end do
      c_don = 0.0_wp
      if (kt == 0) return
      cd = c_live(kt)
      do k = nz, 1, -1
         if (h(k) > H_VANISHED) cd = c_live(k)
         c_don(k) = cd
      end do
   end subroutine donor_conc

   subroutine content_with_i1p(nz, h_old, c, q)
      !! Content `h*c` on a live layer, `h*c_donor` on a filler (invariant I1').
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old(nz), c(nz)
      real(wp), intent(out) :: q(nz)
      real(wp) :: cd(nz)
      call donor_conc(nz, h_old, c, cd)
      q = h_old*cd
   end subroutine content_with_i1p

   subroutine direct_grids(ic, nz, h_old, h_new)
      !! Dyadic thicknesses chosen so that source and target interfaces
      !! coincide EXACTLY (`z_old(ko) == z_new(k)`), with zero-thickness
      !! layers interleaved among live ones, one target spanning many source
      !! cells, and many targets inside one source cell.
      integer, intent(in) :: ic, nz
      real(wp), intent(out) :: h_old(nz), h_new(nz)
      integer :: k, rs, na, nb
      real(wp) :: so, sn
      na = max(1, nz/3)
      nb = max(1, nz/3)
      h_old = 0.0_wp
      h_new = 0.0_wp
      select case (ic)
      case (27)   ! identical stretched grids: every interface is a tie
         do k = 1, nz
            h_old(k) = real(k, wp)*0.25_wp
         end do
         h_new = h_old
      case (28)   ! unit source cells, three thick targets spanning many of them
         h_old = 1.0_wp
         h_new(1) = real(na, wp)
         h_new(2) = real(nb, wp)
         h_new(nz) = real(nz - na - nb, wp)
      case (29)   ! three thick source cells, many unit targets inside each
         h_new = 1.0_wp
         h_old(1) = real(na, wp)
         h_old(2) = real(nb, wp)
         h_old(nz) = real(nz - na - nb, wp)
      case (30)   ! zero-thickness source cells interleaved with live ones
         do k = 2, nz, 2
            h_old(k) = 2.0_wp
         end do
         h_new(1:2*(nz/2)) = 1.0_wp
      case (31)   ! zero-thickness target cells interleaved with live ones
         h_old(1:2*(nz/2)) = 1.0_wp
         do k = 2, nz, 2
            h_new(k) = 2.0_wp
         end do
      case default   ! 32, 33: random dyadic thicknesses with random zeros
         rs = 4099*ic + nz
         so = 0.0_wp
         sn = 0.0_wp
         do k = 1, nz
            if (lcg(rs) > 0.25_wp) h_old(k) = real(1 + int(4.0_wp*lcg(rs)), wp)*0.5_wp
            if (lcg(rs) > 0.25_wp) h_new(k) = real(1 + int(4.0_wp*lcg(rs)), wp)*0.5_wp
            so = so + h_old(k)
            sn = sn + h_new(k)
         end do
         h_old(nz) = h_old(nz) + 1.0_wp
         so = so + 1.0_wp
         h_new(nz) = h_new(nz) + 1.0_wp
         sn = sn + 1.0_wp
         ! equalise the column totals on the top layer (dyadic, exact)
         if (so > sn) then
            h_new(nz) = h_new(nz) + (so - sn)
         else
            h_old(nz) = h_old(nz) + (sn - so)
         end if
      end select
   end subroutine direct_grids

   subroutine rescale_live(nz, h, total)
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h(nz)
      real(wp), intent(in) :: total
      real(wp) :: s
      integer :: k
      s = 0.0_wp
      do k = 1, nz
         if (h(k) > H_VANISHED) s = s + h(k)
      end do
      if (s > 0.0_wp) then
         do k = 1, nz
            if (h(k) > H_VANISHED) h(k) = h(k)*(total/s)
         end do
      end if
   end subroutine rescale_live

   subroutine fill_field(nz, h_old, h_new, hTr, budget)
      integer, intent(in) :: nz
      real(wp), intent(out) :: h_old(NC, 1, nz), h_new(NC, 1, nz)
      real(wp), intent(out) :: hTr(NC, 1, nz), budget(NC, 1, nz)
      real(wp) :: ho(nz), hn(nz), q(nz)
      integer :: ic, k
      do ic = 1, NC
         call build_column(ic, nz, ho, hn, q)
         do k = 1, nz
            h_old(ic, 1, k) = ho(k)
            h_new(ic, 1, k) = hn(k)
            hTr(ic, 1, k) = q(k)
            budget(ic, 1, k) = 1.0e-3_wp*real(ic, wp) + 1.0e-6_wp*real(k, wp)
         end do
      end do
   end subroutine fill_field

   ! ------------------------------------------------------------------ tests

   logical function column_has_filler(nz, h_old, h_new)
      !! Any layer of either grid at or below the vanished marker.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old(nz), h_new(nz)
      column_has_filler = any(h_old <= H_VANISHED) .or. any(h_new <= H_VANISHED)
   end function column_has_filler

   subroutine column_properties(nz, h_old, h_new, hTr0, hTr1, c_new, clamp, n_cons, n_bound, n_fb, n_clamp)
      !! The three claims on ONE remapped column, in concentration form.
      !!  * conservation: `sum(hTr1) - sum(hTr0) - sum(clamp)` is round-off,
      !!    skipped where the column totals differ by more than 1e-9*H_TOT (a
      !!    dead target column keeps no water, by contract);
      !!  * bounding: on a matched column the recorded amount is itself
      !!    round-off (on a mismatched one it also carries the dilution of the
      !!    uncovered sliver, `~c*dh`, by design), and every target `c_new`
      !!    lies inside the range of the source column;
      !!  * fallback: an exactly empty target holds the mean of the source
      !!    cell containing its bottom interface (the topmost cell if the
      !!    interface lies above the column).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old(nz), h_new(nz), hTr0(nz), hTr1(nz), c_new(nz), clamp(nz)
      integer, intent(inout) :: n_cons, n_bound, n_fb, n_clamp
      real(wp) :: c_live(nz), c_old(nz), z_old(nz), zn_lo, s_abs, tol, cmin, cmax
      integer :: k, m

      c_live = 0.0_wp
      do k = 1, nz
         if (h_old(k) > H_VANISHED) c_live(k) = hTr0(k)/h_old(k)
      end do
      call donor_conc(nz, h_old, c_live, c_old)
      s_abs = max(sum(abs(hTr0)), 1.0_wp)
      tol = 1.0e-13_wp*s_abs
      if (abs(sum(h_new) - sum(h_old)) <= 1.0e-9_wp*H_TOT) then
         if (abs(sum(hTr1) - sum(hTr0) - sum(clamp)) > tol) n_cons = n_cons + 1
         if (abs(sum(clamp)) > tol) n_clamp = n_clamp + 1
      end if
      cmin = minval(c_old)
      cmax = maxval(c_old)
      z_old(1) = h_old(1)
      do k = 2, nz
         z_old(k) = z_old(k - 1) + h_old(k)
      end do
      zn_lo = 0.0_wp
      do k = 1, nz
         if (c_new(k) < cmin .or. c_new(k) > cmax) n_bound = n_bound + 1
         if (h_new(k) == 0.0_wp) then
            m = nz
            do while (m > 1)
               if (z_old(m - 1) <= zn_lo) exit
               m = m - 1
            end do
            if (z_old(nz) <= zn_lo) m = nz
            if (c_new(k) /= c_old(m)) n_fb = n_fb + 1
         end if
         zn_lo = zn_lo + h_new(k)
      end do
   end subroutine column_properties

   subroutine compare_kernels(error, nz)
      !! Reference vs stream on the full adversarial field, with AND without
      !! the optional `budget`, to `STREAM_REF_ULPS`; conservation, bounding
      !! and the empty-target fallback on every column of BOTH.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nz
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable :: tr_ref(:, :, :), tr_str(:, :, :), tr0(:, :, :)
      real(wp), allocatable :: bud_ref(:, :, :), bud_str(:, :, :), bud0(:, :, :)
      real(wp), allocatable :: cn_ref(:, :, :), cn_str(:, :, :), cl_ref(:, :, :), cl_str(:, :, :)
      real(wp), allocatable :: c_old(:, :, :)
      real(wp), allocatable :: z_old(:, :, :), z_new(:, :, :), q_lft(:, :, :), q_rgt(:, :, :), q_six(:, :, :)
      integer :: ic, k, n_free_cols
      real(wp) :: worst, worst2, worst_free, scale, u
      integer :: nc_ref, nb_ref, nf_ref, nk_ref, nc_str, nb_str, nf_str, nk_str
      logical :: free
      character(len=160) :: msg

      allocate (h_old(NC, 1, nz), h_new(NC, 1, nz), tr0(NC, 1, nz), bud0(NC, 1, nz))
      allocate (tr_ref(NC, 1, nz), tr_str(NC, 1, nz), bud_ref(NC, 1, nz), bud_str(NC, 1, nz))
      allocate (cn_ref(NC, 1, nz), cn_str(NC, 1, nz), cl_ref(NC, 1, nz), cl_str(NC, 1, nz))
      allocate (c_old(NC, 1, nz), z_old(NC, 1, nz), z_new(NC, 1, nz))
      allocate (q_lft(NC, 1, nz), q_rgt(NC, 1, nz), q_six(NC, 1, nz))
      c_old = 0.0_wp
      z_old = 0.0_wp
      z_new = 0.0_wp
      q_lft = 0.0_wp
      q_rgt = 0.0_wp
      q_six = 0.0_wp
      cn_ref = 0.0_wp
      cn_str = 0.0_wp
      cl_ref = 0.0_wp
      cl_str = 0.0_wp
      call fill_field(nz, h_old, h_new, tr0, bud0)
      tr_ref = tr0
      tr_str = tr0
      bud_ref = bud0
      bud_str = bud0

      !$acc enter data copyin(h_old, h_new, tr_ref, tr_str, bud_ref, bud_str, cn_ref, cn_str, cl_ref, cl_str)
      !$acc enter data copyin(c_old, z_old, z_new, q_lft, q_rgt, q_six)
      call ocean_remap_tracer_field(NC, 1, nz, h_old, h_new, tr_ref, REMAP_PPM, &
                                    .false., .false., budget=bud_ref, &
                                    c_new_out=cn_ref, clamp_out=cl_ref)
      call ocean_remap_tracer_field_stream(NC, 1, nz, h_old, h_new, tr_str, &
                                           c_old, z_old, z_new, q_lft, q_rgt, q_six, &
                                           budget=bud_str, c_new_out=cn_str, clamp_out=cl_str)
      !$acc update self(tr_ref, tr_str, bud_ref, bud_str, cn_ref, cn_str, cl_ref, cl_str)

      worst = 0.0_wp
      worst_free = 0.0_wp
      n_free_cols = 0
      nc_ref = 0
      nb_ref = 0
      nf_ref = 0
      nk_ref = 0
      nc_str = 0
      nb_str = 0
      nf_str = 0
      nk_str = 0
      do ic = 1, NC
         free = .not. column_has_filler(nz, h_old(ic, 1, :), h_new(ic, 1, :))
         if (free) n_free_cols = n_free_cols + 1
         do k = 1, nz
            scale = max(abs(tr0(ic, 1, k)), abs(tr_ref(ic, 1, k)), maxval(abs(tr0(ic, 1, :))), 1.0e-300_wp)
            u = max(abs(tr_ref(ic, 1, k) - tr_str(ic, 1, k)), abs(bud_ref(ic, 1, k) - bud_str(ic, 1, k)), &
                    abs(cl_ref(ic, 1, k) - cl_str(ic, 1, k)))/(epsilon(1.0_wp)*scale)
            ! The target concentration is compared as CONTENT (times h_new): a
            ! few-ulp-thin target amplifies the integral's round-off by
            ! h_src/h_new in c (which is why the bounding exists), not in hTr.
            u = max(u, abs(cn_ref(ic, 1, k) - cn_str(ic, 1, k))*h_new(ic, 1, k)/(epsilon(1.0_wp)*scale))
            worst = max(worst, u)
            if (free) worst_free = max(worst_free, u)
         end do
         call column_properties(nz, h_old(ic, 1, :), h_new(ic, 1, :), tr0(ic, 1, :), tr_ref(ic, 1, :), &
                                cn_ref(ic, 1, :), cl_ref(ic, 1, :), nc_ref, nb_ref, nf_ref, nk_ref)
         call column_properties(nz, h_old(ic, 1, :), h_new(ic, 1, :), tr0(ic, 1, :), tr_str(ic, 1, :), &
                                cn_str(ic, 1, :), cl_str(ic, 1, :), nc_str, nb_str, nf_str, nk_str)
      end do
      ! Without a budget and without the diagnostic outputs (the third
      ! tracer-type call shape).
      tr_ref = tr0
      tr_str = tr0
      !$acc update device(tr_ref, tr_str)
      call ocean_remap_tracer_field(NC, 1, nz, h_old, h_new, tr_ref, REMAP_PPM, &
                                    .false., .false.)
      call ocean_remap_tracer_field_stream(NC, 1, nz, h_old, h_new, tr_str, &
                                           c_old, z_old, z_new, q_lft, q_rgt, q_six)
      !$acc update self(tr_ref, tr_str)
      worst2 = 0.0_wp
      do ic = 1, NC
         do k = 1, nz
            scale = max(abs(tr0(ic, 1, k)), abs(tr_ref(ic, 1, k)), maxval(abs(tr0(ic, 1, :))), 1.0e-300_wp)
            worst2 = max(worst2, abs(tr_ref(ic, 1, k) - tr_str(ic, 1, k))/(epsilon(1.0_wp)*scale))
         end do
      end do
      ! Unmap BEFORE any check can return, so a failing assertion cannot
      ! leave the next test a partially present table.
      !$acc exit data delete(h_old, h_new, tr_ref, tr_str, bud_ref, bud_str, cn_ref, cn_str, cl_ref, cl_str)
      !$acc exit data delete(c_old, z_old, z_new, q_lft, q_rgt, q_six)

      call check(error, n_free_cols > 4, .true., "adversarial field must hold filler-free columns")
      if (allocated(error)) return
      write (msg, '(a,i0,a,es9.2,a,es9.2,a)') "nz=", nz, ": stream vs reference differ by ", worst, &
         " ulp (filler-free columns: ", worst_free, ")"
      call check(error, worst <= STREAM_REF_ULPS, trim(msg))
      if (allocated(error)) return
      write (msg, '(a,i0,a)') "nz=", nz, ": reference column not conserved"
      call check(error, nc_ref, 0, trim(msg))
      if (allocated(error)) return
      call check(error, nc_str, 0, "stream column not conserved")
      if (allocated(error)) return
      call check(error, nb_ref + nb_str, 0, "c_new outside the source range")
      if (allocated(error)) return
      call check(error, nf_ref + nf_str, 0, "empty target does not hold the source-cell mean")
      if (allocated(error)) return
      call check(error, nk_ref + nk_str, 0, "recorded bounding amount is not round-off")
      if (allocated(error)) return

      write (msg, '(a,i0,a,es9.2,a)') "nz=", nz, ": hTr (no budget) differs by ", worst2, " ulp"
      call check(error, worst2 <= STREAM_REF_ULPS, trim(msg))
   end subroutine compare_kernels

   subroutine test_nz200(error)
      !! `nz = 200 > NZ_STACK_MAX`: only the stream can run.  Conservation
      !! (with the recorded bounding amount), bounding and the empty-target
      !! fallback hold on every column.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 200
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :), hTr(:, :, :), tr0(:, :, :)
      real(wp), allocatable :: bud(:, :, :), cn(:, :, :), cl(:, :, :), c_old(:, :, :)
      real(wp), allocatable :: z_old(:, :, :), z_new(:, :, :), q_lft(:, :, :), q_rgt(:, :, :), q_six(:, :, :)
      real(wp) :: s_bud
      integer :: ic, k, n_cons, n_bound, n_fb, n_clamp, n_bud

      allocate (h_old(NC, 1, NZ), h_new(NC, 1, NZ), hTr(NC, 1, NZ), tr0(NC, 1, NZ))
      allocate (bud(NC, 1, NZ), cn(NC, 1, NZ), cl(NC, 1, NZ), c_old(NC, 1, NZ))
      allocate (z_old(NC, 1, NZ), z_new(NC, 1, NZ), q_lft(NC, 1, NZ), q_rgt(NC, 1, NZ), q_six(NC, 1, NZ))
      z_old = 0.0_wp
      z_new = 0.0_wp
      q_lft = 0.0_wp
      q_rgt = 0.0_wp
      q_six = 0.0_wp
      c_old = 0.0_wp
      cn = 0.0_wp
      cl = 0.0_wp
      call fill_field(NZ, h_old, h_new, tr0, bud)
      hTr = tr0
      bud = 0.0_wp
      !$acc enter data copyin(h_old, h_new, hTr, bud, cn, cl, c_old, z_old, z_new, q_lft, q_rgt, q_six)
      call ocean_remap_tracer_field_stream(NC, 1, NZ, h_old, h_new, hTr, &
                                           c_old, z_old, z_new, q_lft, q_rgt, q_six, &
                                           budget=bud, c_new_out=cn, clamp_out=cl)
      !$acc update self(hTr, bud, cn, cl)
      !$acc exit data delete(h_old, h_new, hTr, bud, cn, cl, c_old, z_old, z_new, q_lft, q_rgt, q_six)

      n_cons = 0
      n_bound = 0
      n_fb = 0
      n_clamp = 0
      n_bud = 0
      do ic = 1, NC
         call column_properties(NZ, h_old(ic, 1, :), h_new(ic, 1, :), tr0(ic, 1, :), hTr(ic, 1, :), &
                                cn(ic, 1, :), cl(ic, 1, :), n_cons, n_bound, n_fb, n_clamp)
         ! The budget telescopes to the content change.
         s_bud = 0.0_wp
         do k = 1, NZ
            s_bud = s_bud + bud(ic, 1, k)
         end do
         if (abs(s_bud - (sum(hTr(ic, 1, :)) - sum(tr0(ic, 1, :)))) > &
             1.0e-13_wp*max(sum(abs(tr0(ic, 1, :))), 1.0_wp)) n_bud = n_bud + 1
      end do
      call check(error, n_cons, 0, "nz=200: column content not conserved (with the bounding amount)")
      if (allocated(error)) return
      call check(error, n_bud, 0, "nz=200: budget does not telescope to the content change")
      if (allocated(error)) return
      call check(error, n_bound, 0, "nz=200: c_new outside the source range")
      if (allocated(error)) return
      call check(error, n_fb, 0, "nz=200: empty target does not hold the source-cell mean")
      if (allocated(error)) return
      call check(error, n_clamp, 0, "nz=200: recorded bounding amount is not round-off")
   end subroutine test_nz200

   subroutine test_dispatch(error)
      !! Host dispatch: the stream only for the default PPM column with
      !! `nz >= 3`; everything else keeps the reference.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_vcoord_t) :: vc
      integer :: nz

      call grid%init(4, 4, 1, 1.0_wp, 1.0_wp)
      do nz = 1, 4
         call vc%init(grid, nz_ml=nz)
         vc%remap_stream = .true.
         if (nz <= 2) then
            call check(error, remap_use_stream(vc, nz, REMAP_PPM, .false., .false.), .false., &
                       "nz <= 2 must use the reference")
         else
            call check(error, remap_use_stream(vc, nz, REMAP_PPM, .false., .false.), .true., &
                       "nz >= 3 default PPM must use the stream")
         end if
         if (allocated(error)) return
         if (nz >= 3) then
            call check(error, remap_use_stream(vc, nz, REMAP_PLM, .false., .false.), .false., &
                       "non-PPM method must use the reference")
            if (allocated(error)) return
            call check(error, remap_use_stream(vc, nz, REMAP_PPM, .true., .false.), .false., &
                       "bnd_extrap must use the reference")
            if (allocated(error)) return
            call check(error, remap_use_stream(vc, nz, REMAP_PPM, .false., .true.), .false., &
                       "nonunif must use the reference")
            if (allocated(error)) return
            vc%remap_stream = .false.
            call check(error, remap_use_stream(vc, nz, REMAP_PPM, .false., .false.), .false., &
                       "remap_kernel = reference must use the reference")
            if (allocated(error)) return
         end if
         call vc%destroy()
      end do
   end subroutine test_dispatch

end module test_ocean_remap_stream
