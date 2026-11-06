!! Generic 5-point Boole density integrals for the FV pressure-gradient force.
module rdb_ocean_pgf_reconstruct
   !! The 5-point Boole-quadrature density integrals, through the generic
   !! `eos_t` handle, that feed the finite-volume pressure-gradient force
   !! (FV_MOM6 variant of `rdb_ocean_pressure_force`): the in-layer
   !! vertical rule over a PCM / PLM / PPM sub-layer T/S profile
   !! (`boole_dpa_intz_layer`) and the cross-face rules (`boole_dpa_face`,
   !! `boole_dpa_face_pcm`).  They are the REFERENCE the per-EOS fast paths
   !! are tested against, and the rule the kernel runs for an EOS without a
   !! twin (the linear EOS).
   !!
   !! The sub-layer PLM / PPM edge reconstruction (`plm_edges_layer`,
   !! `ppm_edges_layer`, with the linear-exact boundary pair
   !! `boundary_edges_linear`) and the Wright / Roquet twins of these rules
   !! live in `rdb_ocean_pressure_force`, next to the kernel that calls
   !! them, so they inline into it.
   !!
   !! Algorithm references (cite the PAPER, not other codebases):
   !!   * Adcroft, Hallberg & Harrison (2008), Ocean Modelling 24, 1-2 —
   !!     analytic finite-volume pressure gradient; the layer-integrated
   !!     pressure anomaly `dpa` and its first moment `intz_dpa`.
   !!   * White, Adcroft & Hallberg (2009), J. Comput. Phys. 228 —
   !!     high-order in-layer T/S reconstruction for the density integral.
   !!
   !! Why reconstruct: the PCM (layer-mean) density integral is exact only
   !! for in-layer-uniform density. On thick, sloped layers the moment
   !! error does NOT cancel in the horizontal PGF difference (spurious
   !! acceleration); a monotone PLM/PPM sub-layer profile removes it to
   !! high order.
   !!
   !! Conventions (load-bearing): bottom-up k (k=1 bed, k=nz surface);
   !! within a layer the `_t` (top) edge is SHALLOWER (toward k+1), the
   !! `_b` (bottom) edge is DEEPER (toward k).  z is surface-relative and
   !! negative below the surface; the Boussinesq hydrostatic pressure
   !! estimate at a sub-point is `p = -g*rho0*z = g*rho0*depth`.
   use rdb_constants, only: wp, GRAVITY
   use rdb_eos, only: eos_t, eos_density_point
   implicit none
   private

   public :: boole_dpa_intz_layer
   public :: boole_dpa_face
   public :: boole_dpa_face_pcm

   ! Reconstruction-scheme tags (mirror MOM6 Recon_Scheme; only consulted
   ! when reconstruct_for_pressure is on).
   integer, parameter, public :: PGF_RECON_PLM = 1
      !! Piecewise-linear sub-layer T/S (two-stage h-weighted slope).
   integer, parameter, public :: PGF_RECON_PPM = 2
      !! Piecewise-parabolic sub-layer T/S (implicit-h4 edges + CW limiter
      !! + parabolic s6 integrand).

   integer, parameter :: N_BOOLE = 5
      !! Number of evenly-spaced sub-points in the closed Newton-Cotes
      !! (Boole, 5th-order) quadrature of the in-layer density anomaly.

contains

   pure subroutine boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                        e_top, dz, &
                                        t_t, t_b, t_mean, &
                                        s_t, s_b, s_mean, &
                                        parabolic, dpa, intz_dpa)
      !$omp declare target
      !! 5-point Boole-quadrature density-anomaly integral over one layer
      !! (Adcroft, Hallberg & Harrison 2008; White, Adcroft & Hallberg
      !! 2009).  Returns the layer-integrated pressure-anomaly increment
      !! `dpa = g * int rho' dz` and its first moment (from the layer
      !! TOP edge inward) `intz_dpa = 0.5 * g * dz^2 * bracket`.
      !!
      !! Sub-layer T/S vary between the top (shallower) edge `_t` and the
      !! bottom (deeper) edge `_b`.  When `parabolic` is .false. (PLM) the
      !! profile is linear in the fractional depth; when .true. (PPM) the
      !! curvature `q6 = 3*(2*q_mean - (q_t + q_b))` is added so the
      !! profile is the in-layer parabola through `_t`, `_b` and the layer
      !! mean.
      !!
      !! The density anomaly at each of the 5 evenly-spaced sub-points
      !! (top -> bottom) is rho'(z) = EOS(T,S,p) - rho_ref with the
      !! Boussinesq hydrostatic pressure estimate p = -g*rho0*z.  rho_ref
      !! is subtracted from the absolute EOS density per sub-point.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in)  :: rho0
         !! Boussinesq reference density used in the pressure estimate.
      real(wp), intent(in)  :: rho_ref
         !! Anomaly reference subtracted from the EOS density.
      real(wp), intent(in)  :: e_top
         !! Surface-relative height of the SHALLOWER interface (<= 0).
      real(wp), intent(in)  :: dz
         !! Layer thickness (m), dz >= 0.  Marches down from e_top.
      real(wp), intent(in)  :: t_t, t_b, t_mean
         !! Temperature: top edge, bottom edge, layer mean.
      real(wp), intent(in)  :: s_t, s_b, s_mean
         !! Salinity: top edge, bottom edge, layer mean.
      logical, intent(in)  :: parabolic
         !! .true. -> add the PPM curvature (s6/t6) term.
      real(wp), intent(out) :: dpa
         !! g * int rho' dz over the layer (Pa).
      real(wp), intent(out) :: intz_dpa
         !! 0.5 * g * dz^2 * bracket (Pa*m), first moment from the top.

      real(wp) :: gxrho, wt_t, wt_b, t6, s6
      real(wp) :: t5, s5, z5, p5, rho_anom
      real(wp) :: r5(N_BOOLE)
      integer  :: n

      ! The sub-points and weights are written out in place.  The per-EOS
      ! twins of this rule (`boole_dpa_intz_layer_wright`,
      ! `roquet_recon_dpa_intz`) live in `rdb_ocean_pressure_force`, next to
      ! the kernel that calls them, so they inline into it.
      gxrho = GRAVITY*rho0

      ! PPM curvature (zero for PLM).
      t6 = 0.0_wp
      s6 = 0.0_wp
      if (parabolic) then
         t6 = 3.0_wp*(2.0_wp*t_mean - (t_t + t_b))
         s6 = 3.0_wp*(2.0_wp*s_mean - (s_t + s_b))
      end if

      do n = 1, N_BOOLE
         wt_t = 0.25_wp*real(N_BOOLE - n, wp)   ! 1, .75, .5, .25, 0
         wt_b = 1.0_wp - wt_t
         ! Linear blend + parabolic correction.  At wt_t in [0,1] the
         ! parabola through (_t at wt_t=1, _b at wt_t=0, mean) is
         !   q(wt_t) = wt_t*q_t + wt_b*q_b + q6*wt_t*wt_b.
         t5 = wt_t*t_t + wt_b*t_b + t6*wt_t*wt_b
         s5 = wt_t*s_t + wt_b*s_b + s6*wt_t*wt_b
         z5 = e_top - 0.25_wp*real(n - 1, wp)*dz   ! marches DOWN from top
         p5 = -gxrho*z5
         r5(n) = eos_density_point(eos, t5, s5, p5) - rho_ref
      end do

      rho_anom = (1.0_wp/90.0_wp)*(7.0_wp*(r5(1) + r5(5)) &
                                   + 32.0_wp*(r5(2) + r5(4)) + 12.0_wp*r5(3))
      dpa = GRAVITY*dz*rho_anom
      intz_dpa = 0.5_wp*GRAVITY*dz*dz*(rho_anom &
                                       - (1.0_wp/90.0_wp)*(16.0_wp*(r5(4) - r5(2)) + 7.0_wp*(r5(5) - r5(1))))
   end subroutine boole_dpa_intz_layer

   pure subroutine boole_dpa_face(eos, rho0, rho_ref, &
                                  e_top_l, e_top_r, dz_l, dz_r, &
                                  t_t_l, t_b_l, t_m_l, t_t_r, t_b_r, t_m_r, &
                                  s_t_l, s_b_l, s_m_l, s_t_r, s_b_r, s_m_r, &
                                  dpa_l, dpa_r, parabolic, dpa_face)
      !$omp declare target
      !! HORIZONTAL (cross-face) Boole quadrature of the layer pressure
      !! increment `dpa = g * int rho' dz` — the face integral the FV
      !! pressure-gradient contour needs (Adcroft, Hallberg & Harrison
      !! 2008 §3; Yung, Hallberg, Adcroft & Morrison 2026 §2.4).
      !!
      !! WHY A QUADRATURE AND NOT A MEAN.  The FV assembly evaluates the
      !! top/bottom edges of the control volume as `Delta_e * pbar`, where
      !! `pbar` is the mean pressure ALONG that edge; `pbar` is marched
      !! down from the surface by adding this routine's result layer by
      !! layer.  `Delta_e * pbar` is the exact `int p dz` along the edge
      !! only when `pbar` is the true along-face mean.  Replacing it with
      !! the two-column average `0.5*(dpa_L + dpa_R)` is a TRAPEZOID: it is
      !! exact only if `p` is linear in x along the edge.  With a tilted
      !! interface, `z` is linear in x but `p` is QUADRATIC in z under a
      !! linear stratification, so the trapezoid leaves a curvature
      !! residual `g*(-drho/dz)*Delta_e^2/12` at every interface — the
      !! terrain-following "pressure gradient error of the second kind"
      !! (Haney 1991; Mellor, Ezer & Oey 1994), reported for the sloping
      !! ice-shelf surface as the LINEAR pressure reconstruction by Yung
      !! et al. (2026) §3.1 and cured there by the same device.
      !!
      !! THE FIX.  Sample five evenly-spaced sub-columns across the face.
      !! At fraction `w` from the left column, interpolate the interface
      !! height, the thickness, and the T/S edge + mean triples linearly,
      !! then run the same in-layer vertical Boole quadrature.  Both
      !! interpolations are in the SAME parameter, so a sub-column's
      !! profile is the true profile at the interpolated depth: for T/S
      !! linear in z the sub-column reconstruction is exact and `dpa(w)`
      !! is a quadratic in `w`, which the 5-point closed Newton-Cotes rule
      !! integrates exactly.  The whole PGF then vanishes to round-off on
      !! a resting linear-EOS/linear-stratification column under ANY
      !! layer geometry — the algorithm's defining property.
      !!
      !! THE END POINTS ARE THE COLUMNS' OWN INTEGRALS.  At `w = 0` and
      !! `w = 1` the sub-column IS the left / right column, so its `dpa` is
      !! the one Pass 1 already integrated; the caller passes it in
      !! (`dpa_l` / `dpa_r`, recovered from the `pa` stack as MOM6
      !! `int_density_dz_generic_plm` does) and only the three interior
      !! sub-columns are integrated here.
      !!
      !! COST: 3 sub-columns x 5 sub-points = 15 EOS evaluations per face
      !! per layer (25 before the end points were reused), against 0 for the
      !! trapezoid.  `reconstruct_for_pressure` is opt-in and already the
      !! expensive branch.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in)  :: rho0
         !! Boussinesq reference density used in the pressure estimate.
      real(wp), intent(in)  :: rho_ref
         !! Anomaly reference subtracted from the EOS density.
      real(wp), intent(in)  :: e_top_l, e_top_r
         !! Shallower-interface heights of the layer in the LEFT and
         !! RIGHT columns (m, surface-relative, <= 0).
      real(wp), intent(in)  :: dz_l, dz_r
         !! Layer thicknesses in the left / right columns (m, >= 0).
      real(wp), intent(in)  :: t_t_l, t_b_l, t_m_l
         !! Left column temperature: top edge, bottom edge, layer mean.
      real(wp), intent(in)  :: t_t_r, t_b_r, t_m_r
         !! Right column temperature triple.
      real(wp), intent(in)  :: s_t_l, s_b_l, s_m_l
         !! Left column salinity triple.
      real(wp), intent(in)  :: s_t_r, s_b_r, s_m_r
         !! Right column salinity triple.
      real(wp), intent(in)  :: dpa_l, dpa_r
         !! The left / right columns' own `g * int rho' dz` over the layer
         !! (Pa) -- the `w = 0` / `w = 1` end points of the rule.
      logical, intent(in)  :: parabolic
         !! .true. -> the sub-column profiles carry the PPM curvature.
      real(wp), intent(out) :: dpa_face
         !! Along-face mean of `g * int rho' dz` over the layer (Pa).

      real(wp) :: wr, wl, dpa_m, intz_m, acc
      integer  :: m
      real(wp), parameter :: BOOLE_W(N_BOOLE) = &
                             [7.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp]

      acc = BOOLE_W(1)*dpa_l + BOOLE_W(N_BOOLE)*dpa_r
      do m = 2, N_BOOLE - 1
         wr = 0.25_wp*real(m - 1, wp)   ! 0 at the left column .. 1 at the right
         wl = 1.0_wp - wr
         call boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                   wl*e_top_l + wr*e_top_r, &
                                   wl*dz_l + wr*dz_r, &
                                   wl*t_t_l + wr*t_t_r, &
                                   wl*t_b_l + wr*t_b_r, &
                                   wl*t_m_l + wr*t_m_r, &
                                   wl*s_t_l + wr*s_t_r, &
                                   wl*s_b_l + wr*s_b_r, &
                                   wl*s_m_l + wr*s_m_r, &
                                   parabolic, dpa_m, intz_m)
         acc = acc + BOOLE_W(m)*dpa_m
      end do
      dpa_face = acc/90.0_wp
   end subroutine boole_dpa_face

   pure subroutine boole_dpa_face_pcm(eos, rho0, rho_ref, &
                                      e_top_l, e_top_r, dz_l, dz_r, &
                                      t_l, t_r, s_l, s_r, dpa_l, dpa_r, &
                                      hwt_ll, hwt_lr, hwt_rr, hwt_rl, dpa_face)
      !$omp declare target
      !! The cross-face Boole quadrature of `boole_dpa_face` for a
      !! CONSTANT-BY-LAYER (PCM) T/S column, with MOM6's near-bottom
      !! mass-weighting of the interpolated T/S (MOM6
      !! `int_density_dz_generic_pcm`, `intx_dpa`).
      !!
      !! The two end points are the columns' own vertical integrals
      !! `dpa_l` / `dpa_r` (the caller's `boole_dpa_intz_layer` results).
      !! The three interior sub-columns interpolate the interface height
      !! and thickness LINEARLY in the cross-face fraction, and T/S with
      !! the mass-weighted fractions
      !! `wtT_L = wl*hwt_ll + wr*hwt_rl`, `wtT_R = wl*hwt_lr + wr*hwt_rr`;
      !! `hwt_ll = hwt_rr = 1`, `hwt_lr = hwt_rl = 0` is plain linear
      !! interpolation (no mass weighting).  Each sub-column is integrated
      !! in the vertical at its own IN-SITU pressure `p = -g*rho0*z`.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in)  :: rho0
         !! Boussinesq reference density used in the pressure estimate.
      real(wp), intent(in)  :: rho_ref
         !! Anomaly reference subtracted from the EOS density.
      real(wp), intent(in)  :: e_top_l, e_top_r
         !! Height of the SHALLOWER interface of the layer in the left /
         !! right column (m, geopotential, negative below the datum).
      real(wp), intent(in)  :: dz_l, dz_r
         !! Layer thicknesses in the left / right columns (m, >= 0).
      real(wp), intent(in)  :: t_l, t_r, s_l, s_r
         !! Layer-mean temperature / salinity in the left / right column.
      real(wp), intent(in)  :: dpa_l, dpa_r
         !! The columns' own `g * int rho' dz` over the layer (Pa).
      real(wp), intent(in)  :: hwt_ll, hwt_lr, hwt_rr, hwt_rl
         !! MOM6 `hWt_LL/LR/RR/RL` mass-weighting fractions.
      real(wp), intent(out) :: dpa_face
         !! Along-face mean of `g * int rho' dz` over the layer (Pa).

      real(wp) :: wr, wl, wtt_l, wtt_r, tm, sm, dpa_m, intz_m, acc
      integer  :: m
      real(wp), parameter :: BOOLE_W(N_BOOLE) = &
                             [7.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp]

      acc = BOOLE_W(1)*dpa_l + BOOLE_W(N_BOOLE)*dpa_r
      do m = 2, N_BOOLE - 1
         wr = 0.25_wp*real(m - 1, wp)   ! 0 at the left column .. 1 at the right
         wl = 1.0_wp - wr
         wtt_l = wl*hwt_ll + wr*hwt_rl
         wtt_r = wl*hwt_lr + wr*hwt_rr
         tm = wtt_l*t_l + wtt_r*t_r
         sm = wtt_l*s_l + wtt_r*s_r
         call boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                   wl*e_top_l + wr*e_top_r, &
                                   wl*dz_l + wr*dz_r, &
                                   tm, tm, tm, sm, sm, sm, &
                                   .false., dpa_m, intz_m)
         acc = acc + BOOLE_W(m)*dpa_m
      end do
      dpa_face = acc/90.0_wp
   end subroutine boole_dpa_face_pcm

end module rdb_ocean_pgf_reconstruct
