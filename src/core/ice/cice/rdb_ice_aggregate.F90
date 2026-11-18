!! Cell aggregates of the CICE-port ice state (roundabout-original).
module rdb_ice_aggregate
   !! In-place aggregate of the per-category state into the per-cell totals
   !! `aice`, `vice`, `vsno` and open water `aice0` (the state update of the
   !! ice step; the tracer aggregate is imported with the ITD PR).
   !!
   !! GPU shape (MASTER_PLAN C.1, brief "cell-coupled kernels"): ONE
   !! `do concurrent` over cells, the category loop INSIDE the cell, reading
   !! the SoA arrays `aicen(nx,ny,ncat)` in place through `(i,j,n)` indexing.
   !! No column is gathered into a packed block, nothing is allocated, nothing
   !! is accumulated across different `n` into a shared cell (each work item
   !! owns exactly one cell), so the summation order is the fixed order
   !! `n = 1..ncat` and the result is reproducible.
   !!
   !! The range `is:ie, js:je` is the caller's: the owned cells (a reduction or
   !! diagnostic consumer) or the whole ghosted array (a ghost cell then gets
   !! the aggregate of its own ghost state, which is exactly right once the
   !! halo of the state has been refreshed).
   !!
   !! Non-finite input.  A cell whose `aice + vice + vsno` is not finite gets
   !! non-finite outputs (no `max`/`min` clamp may launder it: `aice0` takes
   !! `aice` itself) and `ICE_ERR_NONFINITE` ORed into `ice_err(i,j)`.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   implicit none
   private

   public :: ice_aggregate_state_impl

contains

   subroutine ice_aggregate_state_impl(nx, ny, ncat, is, ie, js, je, aicen, vicen, vsnon, &
                                       aice, vice, vsno, aice0, ice_err)
      !! `aice = sum_n aicen`, `vice = sum_n vicen`, `vsno = sum_n vsnon`,
      !! `aice0 = max(1 - aice, 0)` over the cells `is:ie, js:je`.  All arrays
      !! must be device-present on an offload build.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)

      integer :: i, j, n
      real(wp) :: a, v, s

      !$acc data present(aicen, vicen, vsnon, aice, vice, vsno, aice0, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, a, v, s)
         a = 0.0_wp
         v = 0.0_wp
         s = 0.0_wp
         do n = 1, ncat
            a = a + aicen(i, j, n)
            v = v + vicen(i, j, n)
            s = s + vsnon(i, j, n)
         end do
         aice(i, j) = a
         vice(i, j) = v
         vsno(i, j) = s
         if (ieee_is_finite(a + v + s)) then
            aice0(i, j) = max(1.0_wp - a, 0.0_wp)
         else
            aice0(i, j) = a + v + s   ! non-finite stays non-finite
            ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
         end if
      end do
      !$acc end data
   end subroutine ice_aggregate_state_impl

end module rdb_ice_aggregate
