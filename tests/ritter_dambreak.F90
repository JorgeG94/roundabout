!! Exact Riemann solver for the 1D shallow water dam break problem
module ritter_dambreak
   !! Provides the exact solution for the dam break (Riemann) problem
   !! following Toro, "Riemann Solvers and Numerical Methods for Fluid
   !! Dynamics", 3rd edition, Chapter 5.
   !!
   !! Valid for wet-bed dam break with zero initial velocities (u_L = u_R = 0).
   !! The left wave is a rarefaction, the right wave is a shock.
   use rdb_constants, only: wp, GRAVITY
   implicit none
   private

   public :: exact_riemann_dambreak
   public :: exact_riemann_sample

contains

   subroutine exact_riemann_dambreak(h_l, h_r, h_star, u_star)
      !! Solve the exact Riemann problem for h_star and u_star
      !!
      !! Uses Newton iteration on f(h) = f_L(h) + f_R(h) = 0
      !! where f_K is the wave function for the K-th wave family.
      real(wp), intent(in) :: h_l
         !! Left water depth (m)
      real(wp), intent(in) :: h_r
         !! Right water depth (m)
      real(wp), intent(out) :: h_star
         !! Star-region water depth (m)
      real(wp), intent(out) :: u_star
         !! Star-region velocity (m/s)

      real(wp) :: a_l, a_r
      real(wp) :: h_guess, f_l, f_r, df_l, df_r, delta
      real(wp) :: q_k
      integer :: iter
      integer, parameter :: MAX_ITER = 50
      real(wp), parameter :: TOL = 1.0e-12_wp

      a_l = sqrt(GRAVITY*h_l)
      a_r = sqrt(GRAVITY*h_r)

      ! Two-rarefaction approximation as initial guess
      h_guess = ((a_l + a_r)/2.0_wp)**2/GRAVITY

      do iter = 1, MAX_ITER
         ! Left wave function and derivative
         if (h_guess <= h_l) then
            ! Rarefaction
            f_l = 2.0_wp*(sqrt(GRAVITY*h_guess) - a_l)
            df_l = sqrt(GRAVITY/h_guess)
         else
            ! Shock
            q_k = sqrt(0.5_wp*GRAVITY*(h_guess + h_l)/(h_guess*h_l))
            f_l = (h_guess - h_l)*q_k
            df_l = q_k - 0.25_wp*GRAVITY*(h_guess - h_l) &
                   /(q_k*h_guess*h_guess)
         end if

         ! Right wave function and derivative
         if (h_guess <= h_r) then
            ! Rarefaction
            f_r = 2.0_wp*(sqrt(GRAVITY*h_guess) - a_r)
            df_r = sqrt(GRAVITY/h_guess)
         else
            ! Shock
            q_k = sqrt(0.5_wp*GRAVITY*(h_guess + h_r)/(h_guess*h_r))
            f_r = (h_guess - h_r)*q_k
            df_r = q_k - 0.25_wp*GRAVITY*(h_guess - h_r) &
                   /(q_k*h_guess*h_guess)
         end if

         ! Newton update: f_L + f_R + (u_R - u_L) = 0, with u_L = u_R = 0
         delta = (f_l + f_r)/(df_l + df_r)
         h_guess = h_guess - delta
         h_guess = max(h_guess, 1.0e-15_wp)

         if (abs(delta) < TOL*h_guess) exit
      end do

      h_star = h_guess
      ! u* = 0.5*(u_L + u_R) + 0.5*(f_R - f_L), with u_L = u_R = 0
      u_star = 0.5_wp*(f_r - f_l)

   end subroutine exact_riemann_dambreak

   subroutine exact_riemann_sample(h_l, h_r, h_star, u_star, &
                                   x_dam, t, x, h_out, u_out)
      !! Sample the exact solution at position x and time t
      !!
      !! Determines which wave region (x,t) falls in and returns the
      !! exact depth and velocity.
      real(wp), intent(in) :: h_l
         !! Left water depth (m)
      real(wp), intent(in) :: h_r
         !! Right water depth (m)
      real(wp), intent(in) :: h_star
         !! Star-region depth from exact_riemann_dambreak
      real(wp), intent(in) :: u_star
         !! Star-region velocity from exact_riemann_dambreak
      real(wp), intent(in) :: x_dam
         !! Dam position (m)
      real(wp), intent(in) :: t
         !! Time (s)
      real(wp), intent(in) :: x
         !! Spatial coordinate (m)
      real(wp), intent(out) :: h_out
         !! Exact water depth at (x, t)
      real(wp), intent(out) :: u_out
         !! Exact velocity at (x, t)

      real(wp) :: a_l, a_r, a_star, xi
      real(wp) :: s_head_l, s_tail_l, s_shock_r

      if (t < 1.0e-15_wp) then
         if (x < x_dam) then
            h_out = h_l
         else
            h_out = h_r
         end if
         u_out = 0.0_wp
         return
      end if

      a_l = sqrt(GRAVITY*h_l)
      a_r = sqrt(GRAVITY*h_r)
      a_star = sqrt(GRAVITY*h_star)

      ! Similarity variable
      xi = (x - x_dam)/t

      ! Left wave (rarefaction for h_star < h_l)
      if (h_star <= h_l) then
         s_head_l = -a_l
         s_tail_l = u_star - a_star

         if (xi <= s_head_l) then
            ! Undisturbed left state
            h_out = h_l
            u_out = 0.0_wp
            return
         else if (xi <= s_tail_l) then
            ! Inside rarefaction fan
            h_out = ((2.0_wp*a_l - xi)/3.0_wp)**2/GRAVITY
            u_out = (2.0_wp/3.0_wp)*(a_l + xi)
            return
         end if
      else
         ! Left shock (shouldn't happen for h_L > h_R dam break, but handle it)
         ! Toro Eq. 10.21: S_L = u_L - a_L * sqrt(h*(h*+h_L)/(2*h_L^2))
         s_head_l = -a_l*sqrt(h_star*(h_star + h_l)/(2.0_wp*h_l*h_l))
         if (xi <= s_head_l) then
            h_out = h_l
            u_out = 0.0_wp
            return
         end if
      end if

      ! Right wave (shock for h_star > h_r)
      if (h_star > h_r) then
         ! Rankine-Hugoniot shock speed (Toro Eq. 10.22)
         s_shock_r = a_r*sqrt(h_star*(h_star + h_r)/(2.0_wp*h_r*h_r))

         if (xi >= s_shock_r) then
            ! Undisturbed right state
            h_out = h_r
            u_out = 0.0_wp
            return
         end if
      else
         ! Right rarefaction (shouldn't happen for standard dam break)
         s_tail_l = u_star + a_star
         if (xi >= s_tail_l) then
            h_out = h_r
            u_out = 0.0_wp
            return
         end if
      end if

      ! Star region
      h_out = h_star
      u_out = u_star

   end subroutine exact_riemann_sample

end module ritter_dambreak
