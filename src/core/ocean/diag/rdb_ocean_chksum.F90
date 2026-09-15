!! MOM6-style per-phase field checksums ("DEBUG=True" analogue).
module rdb_ocean_chksum
   !! Windowed, greppable per-field statistics printed at the split-RK2
   !! phase seams, so (a) two runs can be diffed log-to-log to find the
   !! FIRST diverging operator, and (b) a non-finite value is attributed
   !! to the phase that MINTED it, not the phase that noticed it (the
   !! nan-catch in the truncation / BT fold names the messenger only).
   !!
   !! Row format (fixed-width, grep "CHKSUM"):
   !!   CHKSUM <step> s<stage> <label> <field> <sum> <min> <max> <nonfin> <bits>
   !!
   !! `bits` is a decomposition-INVARIANT POPCNT reduction (int64 sum of
   !! the population count of each element's IEEE-754 bit pattern) —
   !! identical across rank counts / loop orders, so a serial log diffs
   !! against a decomposed run and the first mismatching `bits` names the
   !! diverging operator.  sum/min/max stay per-rank.
   !!
   !! All reductions run device-side (`!$acc parallel loop reduction`,
   !! inert on host builds) — no D->H field copies, so a death-window
   !! sample cadence is affordable on GPU.  NaN semantics: min/max
   !! reductions are NaN-blind (comparisons with NaN are FALSE), so a
   !! contaminated field can print plausible min/max — `nonfin` is the
   !! authoritative corruption signal; `sum` goes NaN with the field.
   !!
   !! Gated by `&ocean_debug_nml chksum` with an optional
   !! `[chksum_start_step, chksum_end_step]` outer-step window (mirrors
   !! debug_ke_attr).  Default off => bit-identical (untaken branches).
   !! Single-rank semantics: sums are per-rank (no halo reduction);
   !! multi-rank runs print one row set per rank.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_halo, only: halo_allreduce_sum_i8
   implicit none
   private

   public :: chksum_probe_t, chksum_stats_t
   public :: chksum_active, chksum_stats_3d, chksum_stats_2d
   public :: chksum_state, chksum_bt, chksum_hotface, chksum_argmax
   public :: rdb_debug_chksum, chksum_loc_extents
   public :: LOC_H, LOC_U, LOC_V, LOC_Q

   integer, parameter :: LOC_H = 1
      !! Tracer / thickness cell centre — extents `(nx_total, ny_total)`.
   integer, parameter :: LOC_U = 2
      !! Arakawa-C u-face (west/east) — extents `(nx_total+1, ny_total)`.
   integer, parameter :: LOC_V = 3
      !! Arakawa-C v-face (south/north) — extents `(nx_total, ny_total+1)`.
   integer, parameter :: LOC_Q = 4
      !! Vorticity / corner point — extents `(nx_total+1, ny_total+1)`.

   interface rdb_debug_chksum
      !! Grid-location-aware checksum: derives the array extents from
      !! `grid` + `loc` (no hand-written bounds at the call site), runs
      !! the device-side reductions, and emits one CHKSUM row.  Resolves
      !! by rank — 3D field trio / 2D BT-work variant.
      module procedure rdb_debug_chksum_3d
      module procedure rdb_debug_chksum_2d
   end interface rdb_debug_chksum

   type :: chksum_probe_t
      !! Config + one-shot header state.  Lives on `ocean_dyn_t`;
      !! host-only (never device-mapped).
      logical :: enable = .false.
         !! Master gate — `&ocean_debug_nml chksum`.
      integer :: start_step = 0
         !! First outer step to sample (0 = from the start).
      integer :: end_step = 0
         !! Last outer step to sample (0 = no upper bound).
      logical :: header_done = .false.
         !! Column-header row emitted.
      logical :: interior = .false.
         !! Restrict every reduction to PHYSICAL cells (exclude the whole
         !! ghost ring).  Default `.false.` = legacy whole-array
         !! behaviour, so existing logs are unchanged.
         !!
         !! Set this to make `bits` a valid 1-rank-vs-N-rank gate.  With
         !! ghosts included the reduced SET of values differs between
         !! decompositions even for correct code (ghosts hold BC junk on
         !! one layout and neighbour data on another), so a `bits` diff
         !! is meaningless.  Restricted to physical cells the set is
         !! identical by construction, and since the POPCNT reduction is
         !! exact integer addition — associative and commutative — a
         !! correct run gives the SAME `bits` on any decomposition.
         !! A stage-by-stage diff then localises a decomposition bug to
         !! the first seam whose `bits` disagree.
   end type chksum_probe_t

   type :: chksum_stats_t
      !! One field's reduction bundle.
      real(wp) :: total = 0.0_wp
         !! Plain (non-compensated) sum — bitwise-stable per build, so
         !! log-diffable between two runs of the SAME binary; NaN when
         !! the field is contaminated.
      real(wp) :: minv = 0.0_wp
      real(wp) :: maxv = 0.0_wp
         !! NaN-blind extrema (see module docstring).
      integer :: nonfin = 0
         !! Count of non-finite entries — the corruption signal.
      integer(int64) :: bits = 0_int64
         !! Decomposition-invariant bitcount: the wrapped int64 sum of
         !! POPCNT over the IEEE-754 bit pattern of every reduced element.
         !! Unlike the FP `total` (only stable within one binary),
         !! integer addition is EXACTLY associative + commutative, so this
         !! is identical across rank counts and loop orders — the
         !! authoritative log-diff signal for a decomposition / halo bug.
   end type chksum_stats_t

contains

   pure function chksum_active(probe, step) result(active)
      !! Gate: enabled AND inside the step window.
      type(chksum_probe_t), intent(in) :: probe
      integer, intent(in) :: step
      logical :: active
      active = probe%enable
      if (active .and. probe%start_step > 0) active = step >= probe%start_step
      if (active .and. probe%end_step > 0) active = step <= probe%end_step
   end function chksum_active

   subroutine chksum_stats_3d(arr, nx, ny, nz, stats, i0, i1, j0, j1)
      !! Device-side sum/min/max/nonfinite over an explicit-shape 3D
      !! field.  No `present` clause: present_or_copyin reads the device
      !! copy in production and copies-in host data in unmapped unit
      !! tests (same convention as the BT fold's loud count).
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: arr(nx, ny, nz)
      type(chksum_stats_t), intent(out) :: stats
      integer, intent(in), optional :: i0, i1, j0, j1
         !! Optional index window (default: the whole array).  Used by the
         !! probe's `interior` mode to exclude the ghost ring.
      real(wp) :: s, mn, mx
      integer(int64) :: bsum
      integer :: nf, i, j, k, ia, ib, ja, jb
      ia = 1
      ib = nx
      ja = 1
      jb = ny
      if (present(i0)) ia = i0
      if (present(i1)) ib = i1
      if (present(j0)) ja = j0
      if (present(j1)) jb = j1

      s = 0.0_wp
      mn = huge(1.0_wp)
      mx = -huge(1.0_wp)
      nf = 0
      bsum = 0_int64
      ! GPU-offload risk: transfer()/popcnt() inside a device region has
      ! NO precedent in this tree.  The intended design accumulates `bsum`
      ! here (device-side, in the same reduction).  If nvfortran cannot
      ! offload transfer/popcnt, this ONE accumulation moves to a
      ! window-cadence host fallback (`!$acc update self(arr)` then a plain
      ! host loop) — a localized one-block change; the FP reductions stay.
      !$acc parallel loop collapse(3) reduction(+:s, nf, bsum) &
      !$acc   reduction(min:mn) reduction(max:mx)
      do k = 1, nz
         do j = ja, jb
            do i = ia, ib
               s = s + arr(i, j, k)
               mn = min(mn, arr(i, j, k))
               mx = max(mx, arr(i, j, k))
               if (.not. ieee_is_finite(arr(i, j, k))) nf = nf + 1
               bsum = bsum + int(popcnt(transfer(arr(i, j, k), 0_int64)), int64)
            end do
         end do
      end do
      stats%total = s
      stats%minv = mn
      stats%maxv = mx
      stats%nonfin = nf
      stats%bits = bsum
   end subroutine chksum_stats_3d

   subroutine chksum_stats_2d(arr, nx, ny, stats, i0, i1, j0, j1)
      !! 2D twin of chksum_stats_3d (BT work fields).
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: arr(nx, ny)
      type(chksum_stats_t), intent(out) :: stats
      integer, intent(in), optional :: i0, i1, j0, j1
         !! Optional index window (default: the whole array) — see the 3D twin.
      real(wp) :: s, mn, mx
      integer(int64) :: bsum
      integer :: nf, i, j, ia, ib, ja, jb
      ia = 1
      ib = nx
      ja = 1
      jb = ny
      if (present(i0)) ia = i0
      if (present(i1)) ib = i1
      if (present(j0)) ja = j0
      if (present(j1)) jb = j1

      s = 0.0_wp
      mn = huge(1.0_wp)
      mx = -huge(1.0_wp)
      nf = 0
      bsum = 0_int64
      ! GPU-offload risk on transfer/popcnt in a device region — see the
      ! matching note in chksum_stats_3d for the host-fallback recipe.
      !$acc parallel loop collapse(2) reduction(+:s, nf, bsum) &
      !$acc   reduction(min:mn) reduction(max:mx)
      do j = ja, jb
         do i = ia, ib
            s = s + arr(i, j)
            mn = min(mn, arr(i, j))
            mx = max(mx, arr(i, j))
            if (.not. ieee_is_finite(arr(i, j))) nf = nf + 1
            bsum = bsum + int(popcnt(transfer(arr(i, j), 0_int64)), int64)
         end do
      end do
      stats%total = s
      stats%minv = mn
      stats%maxv = mx
      stats%nonfin = nf
      stats%bits = bsum
   end subroutine chksum_stats_2d

   subroutine chksum_row(probe, label, stage, step, field, stats)
      !! Emit one CHKSUM row (write(*,...) like the KE_ATTR probe —
      !! probe output bypasses the logger by design: greppable, no
      !! prefix, survives logger-level filtering).
      type(chksum_probe_t), intent(inout) :: probe
      character(len=*), intent(in) :: label, field
      integer, intent(in) :: stage, step
      type(chksum_stats_t), intent(in) :: stats
      integer(int64) :: bits_global

      if (.not. probe%header_done) then
         write (*, "(a)") "# CHKSUM columns: step stage phase field "// &
            "sum min max nonfin bits"
         write (*, "(a)") "# nonfin > 0 = the phase just run minted non-finite values"
         write (*, "(a)") "# bits = decomposition-invariant POPCNT reduction "// &
            "(sum/min/max are per-rank)"
         probe%header_done = .true.
      end if
      ! `bits` is the ONLY column reduced across ranks — its whole purpose
      ! is being identical on 1 rank vs N.  sum/min/max stay per-rank
      ! (documented).  Single-rank ⇒ the facade wrapper is an identity copy.
      call halo_allreduce_sum_i8(stats%bits, bits_global)
      write (*, "(a,1x,i6,1x,a,i1,1x,a16,1x,a12,3(1x,es20.12),1x,i9,1x,i20)") &
         "CHKSUM", step, "s", stage, adjustl(label), adjustl(field), &
         stats%total, stats%minv, stats%maxv, stats%nonfin, bits_global
   end subroutine chksum_row

   pure subroutine chksum_loc_extents(grid, loc, nx, ny)
      !! Map a grid-location tag (`LOC_{H,U,V,Q}`) to the field's 2D
      !! extents on the Arakawa-C grid — u/v faces carry the extra
      !! wall-normal row/column, the corner both.  This is the single
      !! place the C-grid extent convention lives; call sites pass `loc`.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: loc
      integer, intent(out) :: nx, ny

      select case (loc)
      case (LOC_U)
         nx = grid%nx_total + 1
         ny = grid%ny_total
      case (LOC_V)
         nx = grid%nx_total
         ny = grid%ny_total + 1
      case (LOC_Q)
         nx = grid%nx_total + 1
         ny = grid%ny_total + 1
      case default
         ! LOC_H (and any unknown tag) => cell-centre extents.
         nx = grid%nx_total
         ny = grid%ny_total
      end select
   end subroutine chksum_loc_extents

   pure subroutine chksum_loc_interior(grid, loc, i0, i1, j0, j1)
      !! Physical-cell index range for a grid-location tag: the ghost
      !! ring excluded.  Face locations carry one MORE physical entry
      !! than cell centres in their staggered direction (an x-face array
      !! spans `nx_phys + 1` faces), which is why `LOC_U`/`LOC_Q` extend
      !! `i1` by one and `LOC_V`/`LOC_Q` extend `j1`.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: loc
      integer, intent(out) :: i0, i1, j0, j1

      i0 = grid%nghost + 1
      j0 = grid%nghost + 1
      i1 = grid%nghost + grid%nx_phys
      j1 = grid%nghost + grid%ny_phys
      select case (loc)
      case (LOC_U)
         i1 = i1 + 1
      case (LOC_V)
         j1 = j1 + 1
      case (LOC_Q)
         i1 = i1 + 1
         j1 = j1 + 1
      case default
         ! LOC_H (and any unknown tag): plain cell-centre range, already set.
      end select
   end subroutine chksum_loc_interior

   subroutine rdb_debug_chksum_3d(grid, arr, label, field, stage, step, probe, loc)
      !! Location-aware 3D checksum: derives `(nx, ny)` from `grid` + `loc`
      !! and `nz` from the array, so a new instrumentation site is one
      !! line with no hand-written bounds.  Gated by the probe window.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: arr(:, :, :)
      ! assumed-shape-ok: no device loop here; forwards a contiguous
      ! actual to the explicit-shape chksum_stats_3d kernel.
      character(len=*), intent(in) :: label, field
      integer, intent(in) :: stage, step, loc
      type(chksum_probe_t), intent(inout) :: probe
      integer :: nx, ny, nz, ia, ib, ja, jb
      type(chksum_stats_t) :: st

      if (.not. chksum_active(probe, step)) return
      call chksum_loc_extents(grid, loc, nx, ny)
      nz = size(arr, 3)
      if (probe%interior) then
         call chksum_loc_interior(grid, loc, ia, ib, ja, jb)
         call chksum_stats_3d(arr, nx, ny, nz, st, ia, ib, ja, jb)
      else
         call chksum_stats_3d(arr, nx, ny, nz, st)
      end if
      call chksum_row(probe, label, stage, step, field, st)
   end subroutine rdb_debug_chksum_3d

   subroutine rdb_debug_chksum_2d(grid, arr, label, field, stage, step, probe, loc)
      !! 2D twin of rdb_debug_chksum_3d (BT-work fields).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: arr(:, :)
      ! assumed-shape-ok: no device loop here; forwards a contiguous
      ! actual to the explicit-shape chksum_stats_2d kernel.
      character(len=*), intent(in) :: label, field
      integer, intent(in) :: stage, step, loc
      type(chksum_probe_t), intent(inout) :: probe
      integer :: nx, ny, ia, ib, ja, jb
      type(chksum_stats_t) :: st

      if (.not. chksum_active(probe, step)) return
      call chksum_loc_extents(grid, loc, nx, ny)
      if (probe%interior) then
         call chksum_loc_interior(grid, loc, ia, ib, ja, jb)
         call chksum_stats_2d(arr, nx, ny, st, ia, ib, ja, jb)
      else
         call chksum_stats_2d(arr, nx, ny, st)
      end if
      call chksum_row(probe, label, stage, step, field, st)
   end subroutine rdb_debug_chksum_2d

   subroutine chksum_state(grid, ms, probe, label, stage, step)
      !! Sample the prognostic trio (h, u, v) at a phase seam.  Call
      !! AFTER the phase named by `label`; drains device queues first so
      !! the async apply chain has landed.  Bounds come from the
      !! location-aware API (LOC_H / LOC_U / LOC_V), not hand-written.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(in) :: ms
      type(chksum_probe_t), intent(inout) :: probe
      character(len=*), intent(in) :: label
      integer, intent(in) :: stage, step

      if (.not. chksum_active(probe, step)) return
      !$acc wait
      call rdb_debug_chksum(grid, ms%h_layer, label, "h_layer", stage, step, probe, LOC_H)
      call rdb_debug_chksum(grid, ms%u_face_x_layer, label, "u_face", stage, step, probe, LOC_U)
      call rdb_debug_chksum(grid, ms%v_face_y_layer, label, "v_face", stage, step, probe, LOC_V)
   end subroutine chksum_state

   subroutine chksum_bt(grid, bt_work, probe, label, stage, step)
      !! Sample the BT substep's 2D in/out fields at the fold seam:
      !! the fold delta is `bt_ubt_end - ubt_at_n - dt*F_bt_u`, so these
      !! four (+ eta) name which INPUT went non-finite when the fold's
      !! loud count fires — the "what fed it" the nan-catch cannot see.
      type(hgrid_t), intent(in) :: grid
      type(barotropic_workstate_t), intent(in) :: bt_work
      type(chksum_probe_t), intent(inout) :: probe
      character(len=*), intent(in) :: label
      integer, intent(in) :: stage, step

      if (.not. chksum_active(probe, step)) return
      !$acc wait
      call rdb_debug_chksum(grid, bt_work%bt_eta, label, "bt_eta", stage, step, probe, LOC_H)
      call rdb_debug_chksum(grid, bt_work%bt_ubt_end, label, "bt_ubt_end", stage, step, &
                            probe, LOC_U)
      call rdb_debug_chksum(grid, bt_work%bt_vbt_end, label, "bt_vbt_end", stage, step, &
                            probe, LOC_V)
      call rdb_debug_chksum(grid, bt_work%ubt_at_n, label, "ubt_at_n", stage, step, &
                            probe, LOC_U)
      call rdb_debug_chksum(grid, bt_work%F_bt_u, label, "F_bt_u", stage, step, probe, LOC_U)
      call rdb_debug_chksum(grid, bt_work%F_bt_v, label, "F_bt_v", stage, step, probe, LOC_V)
   end subroutine chksum_bt

   subroutine chksum_hotface(grid, ms, visc_rem_u, visc_rem_v, probe, label, stage, step)
      !! Hot-face anatomy at a phase seam: interior argmax |u| and |v|
      !! with the local thickness pair (donor/receiver cells), the
      !! per-face viscous remnant, and the column context (thickness of
      !! the layer below/above at the max face).  The forensic question
      !! this answers: WHICH face takes the explicit dt·F kick, is it an
      !! outcrop edge (massive|vanished thickness pair), and is visc_rem
      !! actually small there (i.e. would MOM6's attenuation have caught
      !! it)?  Row format (grep "HOTFACE"):
      !!   HOTFACE <step> s<stage> <label> u|v (i,j,k) val hL hR rem h_dn h_up
      !! Debug-window only (chksum_active gate); D→H of the two face
      !! fields per sample — affordable at window cadence, not at
      !! production cadence.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: visc_rem_u(:, :, :), visc_rem_v(:, :, :)
      ! assumed-shape-ok: host-side debug probe, no device loops here.
      type(chksum_probe_t), intent(inout) :: probe
      character(len=*), intent(in) :: label
      integer, intent(in) :: stage, step

      integer :: ig, i0, i1, j0, j1, im, jm, km
      real(wp) :: amax

      if (.not. chksum_active(probe, step)) return
      !$acc wait
      !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(visc_rem_u, visc_rem_v)

      ig = grid%nghost
      i0 = ig + 1
      i1 = grid%nx_total - ig
      j0 = ig + 1
      j1 = grid%ny_total - ig

      ! ---- u argmax over interior u-faces ----
      call chksum_argmax(ms%u_face_x_layer, grid%nx_total + 1, grid%ny_total, &
                         ms%nz_ml, i0, i1 + 1, j0, j1, im, jm, km, amax)
      call hotface_row(label, stage, step, "u", im, jm, km, &
                       ms%u_face_x_layer(im, jm, km), &
                       ms%h_layer(im - 1, jm, km), ms%h_layer(im, jm, km), &
                       visc_rem_u(im, jm, km), &
                       ms%h_layer(im - 1, jm, max(km - 1, 1)), &
                       ms%h_layer(im - 1, jm, min(km + 1, ms%nz_ml)))

      ! ---- v argmax over interior v-faces ----
      call chksum_argmax(ms%v_face_y_layer, grid%nx_total, grid%ny_total + 1, &
                         ms%nz_ml, i0, i1, j0, j1 + 1, im, jm, km, amax)
      call hotface_row(label, stage, step, "v", im, jm, km, &
                       ms%v_face_y_layer(im, jm, km), &
                       ms%h_layer(im, jm - 1, km), ms%h_layer(im, jm, km), &
                       visc_rem_v(im, jm, km), &
                       ms%h_layer(im, jm - 1, max(km - 1, 1)), &
                       ms%h_layer(im, jm - 1, min(km + 1, ms%nz_ml)))
   end subroutine chksum_hotface

   pure subroutine chksum_argmax(arr, n1, n2, n3, i0, i1, j0, j1, im, jm, km, amax)
      !! Interior argmax of |arr|: first-encountered strict maximum over
      !! `i ∈ [i0, i1]`, `j ∈ [j0, j1]`, all k (the ghost ring is excluded
      !! by the caller's bounds).  Host-side (debug-window cadence only).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: arr(n1, n2, n3)
      integer, intent(in) :: i0, i1, j0, j1
      integer, intent(out) :: im, jm, km
      real(wp), intent(out) :: amax
      integer :: i, j, k

      amax = -1.0_wp
      im = 0
      jm = 0
      km = 0
      do k = 1, n3
         do j = j0, j1
            do i = i0, i1
               if (abs(arr(i, j, k)) > amax) then
                  amax = abs(arr(i, j, k))
                  im = i
                  jm = j
                  km = k
               end if
            end do
         end do
      end do
   end subroutine chksum_argmax

   subroutine hotface_row(label, stage, step, comp, i, j, k, val, &
                          h_left, h_right, rem, h_dn, h_up)
      !! One HOTFACE row (see chksum_hotface docstring for columns).
      character(len=*), intent(in) :: label, comp
      integer, intent(in) :: stage, step, i, j, k
      real(wp), intent(in) :: val, h_left, h_right, rem, h_dn, h_up

      write (*, "(a,1x,i6,1x,a,i1,1x,a16,1x,a1,' (',i4,',',i4,',',i3,')',6(1x,es13.5))") &
         "HOTFACE", step, "s", stage, adjustl(label), comp, i, j, k, &
         val, h_left, h_right, rem, h_dn, h_up
   end subroutine hotface_row

end module rdb_ocean_chksum
