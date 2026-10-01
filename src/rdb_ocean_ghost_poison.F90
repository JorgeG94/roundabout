!! Ghost-band sentinel-NaN fill for exchange-coverage regression testing.
!!
!! Provides `ocean_poison_ghost_bands`, called at the TOP of each outer step
!! when `dyn%poison_ghosts = .true.` (driven by `&ocean_mpi_nml poison_ghosts`).
!! The knob is default `.false.` — the branch is never taken in production,
!! so the cost is one logical test per outer step.
!!
!! Band-selection rule (avoids false positives on wall runs):
!!   An edge is poisoned iff it is exchange-owned:
!!     MPI seam:   `.not. bc%has_*`  — the neighbour rank fills it each step.
!!     Periodic:   `bc%periodic_x/y` — the exchange/wrap fills it each step
!!                 on BOTH single-rank (local wrap) and multi-rank (message).
!!   Physical wall/OBC edges (has_* = .true., non-periodic) are NEVER poisoned:
!!   those ghosts are BC-owned (zeros from init), read-and-discarded by stencils,
!!   and never refilled by exchanges — poisoning them = false-positive NaN.
!!
!! Cell-centred band limits (nx_total = nxl + 2*ng):
!!   west:  i = 1..ng
!!   east:  i = ng+nxl+1..nx_total   (= ng+nx_phys+1..nx_total)
!!   south: j = 1..ng
!!   north: j = ng+nyl+1..ny_total
!!
!! X-face band limits (nx_face = nx_total+1 = nxl + 2*ng + 1):
!!   west ghost faces:  i = 1..ng
!!   east ghost faces:  i = ng+nxl+2..nx_face   (> ng+nx_phys+1)
!!   NEVER poison i = ng+1 (west seam face) or i = ng+nxl+1 (east seam face):
!!   the two seam-face copies carry ownership/wrap-direction asymmetries that
!!   make it unsafe to overwrite them — see exchange contract §1.1 property (b).
!!   Y-bands of an x-face array are centre-type (j = 1..ng and j = ng+nyl+1..ny_total).
!!
!! Y-face band limits (ny_face = ny_total+1):
!!   NEVER poison j = ng+1 (south seam face) or j = ng+nyl+1 (north seam face).
!!   X-bands of a y-face array are centre-type.
module rdb_ocean_ghost_poison
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   implicit none
   private

   public :: ocean_poison_ghost_bands

contains

   subroutine ocean_poison_ghost_bands(grid, ms, bt_work, ss, &
                                       poison_w, poison_e, poison_s, poison_n)
      !! Sentinel-NaN the exchange-covered ghost bands of all multilayer
      !! prognostic fields and the BT workstate fields.  Called at the TOP
      !! of each outer step (before any exchange or kernel) so that any
      !! kernel consuming an unexchanged ghost produces a loud NaN at the
      !! offending step.
      !!
      !! `poison_w/e/s/n` encode the band-selection rule: pass
      !! `(.not. bc%has_*) .or. bc%periodic_*`.
      !! When all four are `.false.` (e.g. a single-rank all-wall run with no
      !! periodic axes) this routine is a no-op — every ghost is BC-owned.
      !! Do not allocate anything: device-resident state arrays are written
      !! via plain `do concurrent` band-fill kernels, identical to the
      !! periodic-wrap pattern.
      type(hgrid_t), intent(in) :: grid
         !! Horizontal grid (carries nx_total, ny_total, nx_phys, ny_phys, nghost).
      type(multilayer_state_t), intent(inout) :: ms
         !! Multilayer C-grid state whose exchange-covered ghost bands are poisoned.
      type(barotropic_workstate_t), intent(inout) :: bt_work
         !! BT workstate whose exchange-covered ghost bands are poisoned.
      type(ocean_surface_stress_t), intent(inout) :: ss
         !! Surface-stress slot.  Its `tau_x`/`tau_y` bands became
         !! exchange-owned with the forcing-halo fix, so they are poisoned
         !! here too — see the block at the end of the body.
      logical, intent(in) :: poison_w
         !! Poison west ghost band.  `.true.` iff the west edge is exchange-owned
         !! (MPI seam or periodic-x).
      logical, intent(in) :: poison_e
         !! Poison east ghost band.  `.true.` iff the east edge is exchange-owned.
      logical, intent(in) :: poison_s
         !! Poison south ghost band.  `.true.` iff the south edge is exchange-owned.
      logical, intent(in) :: poison_n
         !! Poison north ghost band.  `.true.` iff the north edge is exchange-owned.

      integer :: nxl, nyl, ng, nzt, it
      real(wp) :: qnan

      nxl = grid%nx_phys
      nyl = grid%ny_phys
      ng = grid%nghost
      nzt = ms%nz_ml

      ! Compute the sentinel value on the host — ieee_value is a host intrinsic.
      ! Pass the scalar into the do-concurrent kernels as a local value (device-safe:
      ! plain IEEE representation, no intrinsic call inside the DC body).
      qnan = ieee_value(0.0_wp, ieee_quiet_nan)

      ! ------------------------------------------------------------------
      ! ML h_layer: centre (nx_total, ny_total, nz_ml)
      ! ------------------------------------------------------------------
      call poison_centre_3d(ms%h_layer, nxl, nyl, ng, nzt, qnan, &
                            poison_w, poison_e, poison_s, poison_n)

      ! ------------------------------------------------------------------
      ! ML u_face_x_layer: face-x (nx_total+1, ny_total, nz_ml)
      ! ------------------------------------------------------------------
      call poison_face_x_3d(ms%u_face_x_layer, nxl, nyl, ng, nzt, qnan, &
                            poison_w, poison_e, poison_s, poison_n)

      ! ------------------------------------------------------------------
      ! ML v_face_y_layer: face-y (nx_total, ny_total+1, nz_ml)
      ! ------------------------------------------------------------------
      call poison_face_y_3d(ms%v_face_y_layer, nxl, nyl, ng, nzt, qnan, &
                            poison_w, poison_e, poison_s, poison_n)

      ! ------------------------------------------------------------------
      ! Tracers: per-tracer hTr — centre (nx_total, ny_total, nz_ml).
      ! Outer-shim rule: array-of-DT cannot be dereferenced on-device;
      ! the loop is on the HOST, each element passed by explicit-shape helper.
      ! ------------------------------------------------------------------
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (allocated(ms%tracers(it)%hTr)) then
               call poison_centre_3d(ms%tracers(it)%hTr, nxl, nyl, ng, nzt, qnan, &
                                     poison_w, poison_e, poison_s, poison_n)
            end if
         end do
      end if

      ! ------------------------------------------------------------------
      ! BT workstate: bt_eta (centre), bt_ubt (face-x), bt_vbt (face-y).
      ! These are re-derived/exchanged early in each stage.  Poisoning
      ! them proves the derive/coupling+exchange chain covers the bands.
      ! Guard with allocated() — the BT workstate is only allocated when
      ! ocean_dyn_init has been called with nz_ml present (split path).
      ! ------------------------------------------------------------------
      if (allocated(bt_work%bt_eta)) then
         call poison_centre_2d(bt_work%bt_eta, nxl, nyl, ng, qnan, &
                               poison_w, poison_e, poison_s, poison_n)
      end if
      if (allocated(bt_work%bt_ubt)) then
         call poison_face_x_2d(bt_work%bt_ubt, nxl, nyl, ng, qnan, &
                               poison_w, poison_e, poison_s, poison_n)
      end if
      if (allocated(bt_work%bt_vbt)) then
         call poison_face_y_2d(bt_work%bt_vbt, nxl, nyl, ng, qnan, &
                               poison_w, poison_e, poison_s, poison_n)
      end if

      ! ------------------------------------------------------------------
      ! Surface stress (face-x / face-y).  Exchange-owned as of the
      ! forcing-halo fix: `ocean_seam_refresh_surface_stress` fills these
      ! bands at configure and after every file blend.  Poisoning them is
      ! what proves that — several kernels read one cell beyond their own
      ! (`stress_mag` face-to-centre average feeding KPP/EPBL, and the MLE
      ! corner average), so an unexchanged band shows up as NaN rather
      ! than as a quiet decomposition-dependent answer.
      ! ------------------------------------------------------------------
      if (allocated(ss%tau_x)) then
         call poison_face_x_2d(ss%tau_x, nxl, nyl, ng, qnan, &
                               poison_w, poison_e, poison_s, poison_n)
      end if
      if (allocated(ss%tau_y)) then
         call poison_face_y_2d(ss%tau_y, nxl, nyl, ng, qnan, &
                               poison_w, poison_e, poison_s, poison_n)
      end if

   end subroutine ocean_poison_ghost_bands

   ! ===========================================================================
   ! Private flat-impl helpers — explicit-shape so the dc-assumed-shape
   ! pre-commit hook passes; device-resident arrays are written via
   ! plain do concurrent (device-safe on GPU builds).
   ! ===========================================================================

   subroutine poison_centre_2d(fld, nxl, nyl, ng, qnan, pw, pe, ps, pn)
      !! Poison centre-type 2D ghost bands.
      integer, intent(in) :: nxl, nyl, ng
      real(wp), intent(inout) :: fld(nxl + 2*ng, nyl + 2*ng)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxt, nyt, i, j

      nxt = nxl + 2*ng
      nyt = nyl + 2*ng

      ! West band: i = 1..ng
      if (pw) then
         do concurrent(j=1:nyt, i=1:ng)
            fld(i, j) = qnan
         end do
      end if
      ! East band: i = ng+nxl+1..nxt
      if (pe) then
         do concurrent(j=1:nyt, i=ng + nxl + 1:nxt)
            fld(i, j) = qnan
         end do
      end if
      ! South band: j = 1..ng
      if (ps) then
         do concurrent(j=1:ng, i=1:nxt)
            fld(i, j) = qnan
         end do
      end if
      ! North band: j = ng+nyl+1..nyt
      if (pn) then
         do concurrent(j=ng + nyl + 1:nyt, i=1:nxt)
            fld(i, j) = qnan
         end do
      end if
   end subroutine poison_centre_2d

   subroutine poison_centre_3d(fld, nxl, nyl, ng, nz, qnan, pw, pe, ps, pn)
      !! Poison centre-type 3D ghost bands.
      integer, intent(in) :: nxl, nyl, ng, nz
      real(wp), intent(inout) :: fld(nxl + 2*ng, nyl + 2*ng, nz)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxt, nyt, i, j, k

      nxt = nxl + 2*ng
      nyt = nyl + 2*ng

      ! West band: i = 1..ng
      if (pw) then
         do concurrent(k=1:nz, j=1:nyt, i=1:ng)
            fld(i, j, k) = qnan
         end do
      end if
      ! East band: i = ng+nxl+1..nxt
      if (pe) then
         do concurrent(k=1:nz, j=1:nyt, i=ng + nxl + 1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
      ! South band: j = 1..ng
      if (ps) then
         do concurrent(k=1:nz, j=1:ng, i=1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
      ! North band: j = ng+nyl+1..nyt
      if (pn) then
         do concurrent(k=1:nz, j=ng + nyl + 1:nyt, i=1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
   end subroutine poison_centre_3d

   subroutine poison_face_x_2d(fld, nxl, nyl, ng, qnan, pw, pe, ps, pn)
      !! Poison x-face 2D ghost bands.
      !! fld has shape (nx_total+1, ny_total) = (nxl+2*ng+1, nyl+2*ng).
      !! NEVER poison i = ng+1 (west seam face) or i = ng+nxl+1 (east seam face):
      !! these carry ownership/wrap-direction asymmetries that make it unsafe
      !! to overwrite them.  West ghost faces: i = 1..ng.
      !! East ghost faces: i = ng+nxl+2..nx_face.  Y-bands centre-type.
      integer, intent(in) :: nxl, nyl, ng
      real(wp), intent(inout) :: fld(nxl + 2*ng + 1, nyl + 2*ng)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxf, nyt, i, j

      nxf = nxl + 2*ng + 1
      nyt = nyl + 2*ng

      ! West ghost faces: i = 1..ng  (NOT i = ng+1 — that is the seam face)
      if (pw) then
         do concurrent(j=1:nyt, i=1:ng)
            fld(i, j) = qnan
         end do
      end if
      ! East ghost faces: i = ng+nxl+2..nxf  (NOT i = ng+nxl+1 — seam face)
      if (pe) then
         do concurrent(j=1:nyt, i=ng + nxl + 2:nxf)
            fld(i, j) = qnan
         end do
      end if
      ! South band (centre-type): j = 1..ng
      if (ps) then
         do concurrent(j=1:ng, i=1:nxf)
            fld(i, j) = qnan
         end do
      end if
      ! North band (centre-type): j = ng+nyl+1..nyt
      if (pn) then
         do concurrent(j=ng + nyl + 1:nyt, i=1:nxf)
            fld(i, j) = qnan
         end do
      end if
   end subroutine poison_face_x_2d

   subroutine poison_face_x_3d(fld, nxl, nyl, ng, nz, qnan, pw, pe, ps, pn)
      !! Poison x-face 3D ghost bands.  Seam-face exclusion: see poison_face_x_2d.
      integer, intent(in) :: nxl, nyl, ng, nz
      real(wp), intent(inout) :: fld(nxl + 2*ng + 1, nyl + 2*ng, nz)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxf, nyt, i, j, k

      nxf = nxl + 2*ng + 1
      nyt = nyl + 2*ng

      if (pw) then
         do concurrent(k=1:nz, j=1:nyt, i=1:ng)
            fld(i, j, k) = qnan
         end do
      end if
      if (pe) then
         do concurrent(k=1:nz, j=1:nyt, i=ng + nxl + 2:nxf)
            fld(i, j, k) = qnan
         end do
      end if
      if (ps) then
         do concurrent(k=1:nz, j=1:ng, i=1:nxf)
            fld(i, j, k) = qnan
         end do
      end if
      if (pn) then
         do concurrent(k=1:nz, j=ng + nyl + 1:nyt, i=1:nxf)
            fld(i, j, k) = qnan
         end do
      end if
   end subroutine poison_face_x_3d

   subroutine poison_face_y_2d(fld, nxl, nyl, ng, qnan, pw, pe, ps, pn)
      !! Poison y-face 2D ghost bands.
      !! fld has shape (nx_total, ny_total+1) = (nxl+2*ng, nyl+2*ng+1).
      !! NEVER poison j = ng+1 (south seam face) or j = ng+nyl+1 (north seam face).
      !! South ghost faces: j = 1..ng.  North ghost faces: j = ng+nyl+2..ny_face.
      !! X-bands are centre-type: i = 1..ng and i = ng+nxl+1..nx_total.
      integer, intent(in) :: nxl, nyl, ng
      real(wp), intent(inout) :: fld(nxl + 2*ng, nyl + 2*ng + 1)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxt, nyf, i, j

      nxt = nxl + 2*ng
      nyf = nyl + 2*ng + 1

      ! West band (centre-type): i = 1..ng
      if (pw) then
         do concurrent(j=1:nyf, i=1:ng)
            fld(i, j) = qnan
         end do
      end if
      ! East band (centre-type): i = ng+nxl+1..nxt
      if (pe) then
         do concurrent(j=1:nyf, i=ng + nxl + 1:nxt)
            fld(i, j) = qnan
         end do
      end if
      ! South ghost faces: j = 1..ng  (NOT j = ng+1 — seam face)
      if (ps) then
         do concurrent(j=1:ng, i=1:nxt)
            fld(i, j) = qnan
         end do
      end if
      ! North ghost faces: j = ng+nyl+2..nyf  (NOT j = ng+nyl+1 — seam face)
      if (pn) then
         do concurrent(j=ng + nyl + 2:nyf, i=1:nxt)
            fld(i, j) = qnan
         end do
      end if
   end subroutine poison_face_y_2d

   subroutine poison_face_y_3d(fld, nxl, nyl, ng, nz, qnan, pw, pe, ps, pn)
      !! Poison y-face 3D ghost bands.  Seam-face exclusion: see poison_face_y_2d.
      integer, intent(in) :: nxl, nyl, ng, nz
      real(wp), intent(inout) :: fld(nxl + 2*ng, nyl + 2*ng + 1, nz)
      real(wp), intent(in) :: qnan
      logical, intent(in) :: pw, pe, ps, pn

      integer :: nxt, nyf, i, j, k

      nxt = nxl + 2*ng
      nyf = nyl + 2*ng + 1

      if (pw) then
         do concurrent(k=1:nz, j=1:nyf, i=1:ng)
            fld(i, j, k) = qnan
         end do
      end if
      if (pe) then
         do concurrent(k=1:nz, j=1:nyf, i=ng + nxl + 1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
      if (ps) then
         do concurrent(k=1:nz, j=1:ng, i=1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
      if (pn) then
         do concurrent(k=1:nz, j=ng + nyl + 2:nyf, i=1:nxt)
            fld(i, j, k) = qnan
         end do
      end if
   end subroutine poison_face_y_3d

end module rdb_ocean_ghost_poison
