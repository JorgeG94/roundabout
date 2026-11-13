!! Massless-layer merge column helper (ocean dyn-core, D4).
module rdb_massless
   !! Pure column helpers that merge vanished (sub-`H_VANISHED`) layers into a
   !! "massive" grid (`nzc <= nz`) so the column solver never divides by a
   !! vanished thickness: thickness-weighted means onto the merged grid, plus
   !! the inverse interpolation of an interface quantity back to the original
   !! interfaces. Reference: Jackson, Hallberg & Legg (2008), JPO 38, 1033.
   !!
   !! ORDERING (load-bearing): LOCAL SURFACE-DOWN indices — local `k=1` =
   !! surface, `k=nz` = bed; interfaces `K=1` (surface) ... `K=nz+1` (bed).
   !! The consumer does the global(bed-up) <-> local(surface-down) flip; the
   !! helper never sees global indices.
   !!
   !! MERGE DIRECTION (do not invent a symmetric rule): walking surface-down, a
   !! merged layer opens only when the current cluster already has mass AND the
   !! incoming layer is massive (> h_min). So a vanished layer folds into the
   !! preceding (surfaceward) massive layer; a leading surface run folds into
   !! the first massive layer beneath; a trailing bed run folds into the last
   !! massive layer above.
   !!
   !! GPU NOTE: `!$acc routine seq`, callable cross-module from the kappa-shear
   !! `do concurrent` column kernel. NVHPC will NOT inline across the module
   !! boundary, but these run once per column per compute (thermo cadence), so
   !! the divergence cost is acceptable. All arrays are fixed-size
   !! `NZ_STACK_MAX`(+1) explicit-shape locals (no descriptor walk, no alloc).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, H_DIV_EPS
#endif
   implicit none
   private
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: massless_build_maps
   public :: massless_merge_fields
   public :: massless_interp_back

   integer, parameter :: NZL = NZ_STACK_MAX
      !! Layer-array dimension (must match the consumer's NZL).
   integer, parameter :: NZLI = NZ_STACK_MAX + 1
      !! Interface-array dimension (= NZL + 1).

contains

   pure subroutine massless_build_maps(h, nz, h_min, nzc, hc, kc, kf)
      !! Build the merged massive-layer grid and the maps `kc`, `kf`.
      !! `kc(k)` (k=1..nz) = merged-layer index original layer `k` folds into;
      !! `kc(nz+1) = nzc+1` is the bed sentinel. `kf(K)` = fractional position
      !! of original interface `K` inside merged layer `kc(K)` (0 = coincides
      !! with a merged interface -> pure lookup; (0,1) = interior -> linear
      !! blend in `massless_interp_back`).
      !! Identity case (no layer < h_min): `nzc==nz`, `kc(k)==k`, `kf(k)==0`,
      !! `hc==h` bit-for-bit.
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(NZL)
         !! Local surface-down layer thicknesses (m), >= 0.
      real(wp), intent(in) :: h_min
         !! Merge threshold (H_VANISHED).
      integer, intent(out) :: nzc
         !! Number of merged massive layers (<= nz).
      real(wp), intent(out) :: hc(NZL)
         !! Merged thicknesses (m), valid 1..nzc.
      integer, intent(out) :: kc(NZLI)
         !! Layer/interface -> merged index, valid 1..nz+1 (sentinel at nz+1).
      real(wp), intent(out) :: kf(NZLI)
         !! Fractional interface weight, valid 1..nz+1.

      integer :: k
      real(wp) :: dz_in

      ! Build loop (surface-down, 1-based).
      nzc = 1
      hc(1) = 0.0_wp
      do k = 1, nz
         ! Open a new merged layer iff the current cluster has mass AND
         ! this layer is massive.
         if (hc(nzc) > 0.0_wp .and. h(k) > h_min) then
            nzc = nzc + 1
            hc(nzc) = 0.0_wp
         end if
         kc(k) = nzc
         hc(nzc) = hc(nzc) + h(k)
      end do
      kc(nz + 1) = nzc + 1   ! bed-interface sentinel

      ! --- kf: interface interpolation weights ---
      kf(1) = 0.0_wp
      dz_in = h(1)
      do k = 2, nz
         if (kc(k) > kc(k - 1)) then
            ! Interface k sits ON a merged interface.
            kf(k) = 0.0_wp
            dz_in = h(k)
         else
            ! Interface k is interior to merged layer kc(k).
            kf(k) = dz_in/max(hc(kc(k)), H_DIV_EPS)
            dz_in = dz_in + h(k)
         end if
      end do
      kf(nz + 1) = 0.0_wp
   end subroutine massless_build_maps

   pure subroutine massless_merge_fields(h, kc, nz, nzc, &
                                         u, v, t, s, uc, vc, tc, sc)
      !! Thickness-weighted merged means for u, v, T, S — the solver receives
      !! MEANS, not integrals, and must not re-divide. Accumulation runs in the
      !! same surface-down `k` order as `massless_build_maps` so the
      !! column-integral conservation holds to round-off.
      !$acc routine seq
      integer, intent(in) :: nz, nzc
      real(wp), intent(in) :: h(NZL)
      integer, intent(in) :: kc(NZLI)
      real(wp), intent(in) :: u(NZL), v(NZL), t(NZL), s(NZL)
      real(wp), intent(out) :: uc(NZL), vc(NZL), tc(NZL), sc(NZL)

      integer :: k, kk
      real(wp) :: hc_acc(NZL)
      real(wp) :: denom

      do kk = 1, nzc
         hc_acc(kk) = 0.0_wp
         uc(kk) = 0.0_wp
         vc(kk) = 0.0_wp
         tc(kk) = 0.0_wp
         sc(kk) = 0.0_wp
      end do

      do k = 1, nz
         kk = kc(k)
         hc_acc(kk) = hc_acc(kk) + h(k)
         uc(kk) = uc(kk) + u(k)*h(k)
         vc(kk) = vc(kk) + v(k)*h(k)
         tc(kk) = tc(kk) + t(k)*h(k)
         sc(kk) = sc(kk) + s(k)*h(k)
      end do

      ! Finalise means.  hc_acc(kk) > 0 for every kk (a cluster opens only
      ! on a massive layer or holds an absorbed leading run + its first
      ! massive layer), so H_DIV_EPS is armour, not a clamp.
      do kk = 1, nzc
         denom = max(hc_acc(kk), H_DIV_EPS)
         uc(kk) = uc(kk)/denom
         vc(kk) = vc(kk)/denom
         tc(kk) = tc(kk)/denom
         sc(kk) = sc(kk)/denom
      end do
   end subroutine massless_merge_fields

   pure subroutine massless_interp_back(qc, kc, kf, nz, q)
      !! Inverse map: interpolate an interface quantity `qc(1:nzc+1)` on the
      !! merged grid back to the original `nz+1` interfaces. `kf==0` takes a
      !! pure lookup; interior interfaces take the linear blend between merged
      !! interfaces `kc(K)` and `kc(K)+1` (in range since `kf>0` only when
      !! `kc(K) < nzc+1`). Identity case: `kf==0` everywhere -> `q==qc`.
      !$acc routine seq
      integer, intent(in) :: nz
      real(wp), intent(in) :: qc(NZLI)
         !! Interface quantity on the merged grid (valid 1..nzc+1).
      integer, intent(in) :: kc(NZLI)
      real(wp), intent(in) :: kf(NZLI)
      real(wp), intent(out) :: q(NZLI)
         !! Interface quantity on the original grid (filled 1..nz+1).

      integer :: kk

      do kk = 1, nz + 1
         if (kf(kk) == 0.0_wp) then
            q(kk) = qc(kc(kk))
         else
            q(kk) = (1.0_wp - kf(kk))*qc(kc(kk)) + kf(kk)*qc(kc(kk) + 1)
         end if
      end do
   end subroutine massless_interp_back

end module rdb_massless
