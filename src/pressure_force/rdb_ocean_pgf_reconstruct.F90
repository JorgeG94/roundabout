!! Sub-layer T/S edge reconstruction for the FV pressure-gradient force.
module rdb_ocean_pgf_reconstruct
   !! Per-column PLM / PPM edge-value reconstruction of a layer-mean
   !! tracer profile, plus the 5-point Boole-quadrature density integral
   !! that feeds the finite-volume pressure-gradient force (FV_MOM6
   !! variant of `rdb_ocean_pressure_force`).
   !!
   !! Algorithm references (cite the PAPER, not other codebases):
   !!   * Adcroft, Hallberg & Harrison (2008), Ocean Modelling 24, 1-2 —
   !!     analytic finite-volume pressure gradient; the layer-integrated
   !!     pressure anomaly `dpa` and its first moment `intz_dpa`.
   !!   * White, Adcroft & Hallberg (2009), J. Comput. Phys. 228 —
   !!     high-order in-layer T/S reconstruction for the density integral.
   !!   * Colella & Woodward (1984), J. Comput. Phys. 54 — the PPM
   !!     parabola limiter reused for the monotone edge values.
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
   !!
   !! Boundary layers (k=1, k=nz) take a LINEAR-EXACT one-sided edge pair
   !! from the single interior neighbour (`boundary_edges_linear`), not a
   !! PCM flatten.  That matters: under a terrain-following coordinate the
   !! layers adjacent to the tilted boundary are exactly where the
   !! sigma pressure-gradient truncation error lives, and flattening them
   !! to PCM left the FULL error there while the interior was corrected.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, GRAVITY
#endif
   use rdb_eos, only: eos_t, eos_density_point
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: plm_edges_column
   public :: ppm_edges_column
   public :: boole_dpa_intz_layer
   public :: boole_dpa_face

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

   pure subroutine boundary_edges_linear(h_self, h_nbr, q_self, dq_up, q_t, q_b)
      !$acc routine seq
      !! Linear-exact one-sided edge pair for a BOUNDARY layer (k=1 or
      !! k=nz), where a centred slope has no second neighbour.
      !!
      !! `dq_up` is the layer-mean increment toward the SURFACE across the
      !! two cell centres (`q(2)-q(1)` at the bed, `q(nz)-q(nz-1)` at the
      !! surface).  The centres are `(h_self + h_nbr)/2` apart, so the
      !! per-metre slope is `dq_up/((h_self+h_nbr)/2)` and the half-jump
      !! across this layer is
      !!
      !!     d = dq_up * h_self / (h_self + h_nbr)
      !!
      !! giving `q_t = q + d` (shallower edge) and `q_b = q - d`.  For a
      !! profile that is linear in z this reproduces the true edge values
      !! EXACTLY, for any thickness pair — which is the property the FV
      !! pressure-gradient quadrature needs (Adcroft, Hallberg & Harrison
      !! 2008; White, Adcroft & Hallberg 2009 §2): a PCM flatten here
      !! leaves the full terrain-following truncation error in the layers
      !! next to the tilted boundary.
      !!
      !! Limiter: `|d| <= |dq_up|`, i.e. the edge never leaves the
      !! interval the two cell means span on the other side.  Since
      !! `h_self/(h_self+h_nbr) < 1` it never bites on a real thickness
      !! pair — it is armour against a degenerate `h_nbr <= 0`, and it
      !! keeps the extrapolation from manufacturing a density inversion.
      real(wp), intent(in)  :: h_self
         !! Thickness of the boundary layer itself (m).
      real(wp), intent(in)  :: h_nbr
         !! Thickness of its single interior neighbour (m).
      real(wp), intent(in)  :: q_self
         !! Layer mean of the boundary layer.
      real(wp), intent(in)  :: dq_up
         !! Layer-mean increment toward the surface, neighbour -> self at
         !! the surface layer, self -> neighbour at the bed layer.
      real(wp), intent(out) :: q_t
         !! Top (shallower) edge value.
      real(wp), intent(out) :: q_b
         !! Bottom (deeper) edge value.

      real(wp), parameter :: H_TINY = 1.0e-30_wp
      real(wp) :: d

      d = dq_up*h_self/max(h_self + h_nbr, H_TINY)
      d = sign(min(abs(d), abs(dq_up)), d)
      q_t = q_self + d
      q_b = q_self - d
   end subroutine boundary_edges_linear

   pure subroutine plm_edges_column(nz, h, q, q_t, q_b)
      !$acc routine seq
      !! Per-column PLM top/bottom edge values of a layer-mean field `q`,
      !! via a two-stage h-weighted van-Leer slope (White, Adcroft &
      !! Hallberg 2009 §2). Returns the SHALLOWER edge in `q_t` (toward
      !! k+1) and the DEEPER edge in `q_b` (toward k-1), bottom-up.
      !! Boundary layers (k=1, k=nz) -> `boundary_edges_linear`, the
      !! linear-exact one-sided pair.
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz)
         !! Layer thicknesses (m), k=1 bed .. k=nz surface.
      real(wp), intent(in)  :: q(nz)
         !! Layer-mean scalar (T or S).
      real(wp), intent(out) :: q_t(nz)
         !! Top (shallower) edge value per layer.
      real(wp), intent(out) :: q_b(nz)
         !! Bottom (deeper) edge value per layer.

      real(wp) :: slp(NZ_STACK_MAX)
      real(wp) :: h_l, h_c, h_r, sig_c, sig_l, sig_r, slp_max
      real(wp) :: e_t, e_b, q_lo, q_hi
      integer  :: k

      ! Single-layer column: PCM is the only option (no neighbour).
      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if

      ! ---- Stage 1: h-weighted limited central slope per interior layer ----
      ! sig_c is the change ACROSS the layer measured deeper->shallower:
      ! positive sig_c means q increases toward the surface (k+1).  In
      ! bottom-up indexing the shallower neighbour is k+1, the deeper is
      ! k-1.  h-weighted central slope (van-Leer / White-Adcroft-Hallberg):
      !   sig_c = (q(k+1)-q(k-1)) * h(k) / (h(k-1)+2 h(k)+h(k+1))  * 2
      ! then limited to 2*min(|q(k)-q_deeper|,|q_shallower-q(k)|), zeroed
      ! at extrema.
      slp(1) = 0.0_wp
      slp(nz) = 0.0_wp
      do k = 2, nz - 1
         h_l = h(k - 1)   ! deeper neighbour thickness
         h_c = h(k)
         h_r = h(k + 1)   ! shallower neighbour thickness
         sig_l = q(k) - q(k - 1)     ! deeper one-sided (toward k-1)
         sig_r = q(k + 1) - q(k)     ! shallower one-sided (toward k+1)
         if (sig_l*sig_r <= 0.0_wp) then
            slp(k) = 0.0_wp          ! local extremum -> flatten
         else
            sig_c = 2.0_wp*(q(k + 1) - q(k - 1))*h_c/(h_l + 2.0_wp*h_c + h_r)
            slp_max = 2.0_wp*min(abs(sig_l), abs(sig_r))
            slp(k) = sign(min(abs(sig_c), slp_max), sig_c)
         end if
      end do

      ! ---- Stage 2: monotonized edges bounded against neighbour edges ----
      ! Build the raw edges from the limited slope, then clamp each edge
      ! so it lies between the cell mean and the adjacent cell mean
      ! (White, Adcroft & Hallberg 2009 §2 monotonization — prevents the
      ! reconstructed edge from over/undershooting the neighbour mean,
      ! which would manufacture a density inversion under the EOS).
      call boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call boundary_edges_linear(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                                 q_t(nz), q_b(nz))
      do k = 2, nz - 1
         e_t = q(k) + 0.5_wp*slp(k)   ! shallower edge (toward k+1)
         e_b = q(k) - 0.5_wp*slp(k)   ! deeper edge (toward k-1)
         ! Bound the shallower edge between q(k) and q(k+1).
         q_lo = min(q(k), q(k + 1))
         q_hi = max(q(k), q(k + 1))
         q_t(k) = max(q_lo, min(q_hi, e_t))
         ! Bound the deeper edge between q(k) and q(k-1).
         q_lo = min(q(k), q(k - 1))
         q_hi = max(q(k), q(k - 1))
         q_b(k) = max(q_lo, min(q_hi, e_b))
      end do
   end subroutine plm_edges_column

   pure subroutine ppm_edges_column(nz, h, q, q_t, q_b)
      !$acc routine seq
      !! Per-column PPM top/bottom edge values via the implicit-h4
      !! (thickness-weighted, exactly 4th-order on non-uniform layers)
      !! interface estimate + Colella & Woodward (1984) parabola limiter.
      !! Returns the SHALLOWER edge in `q_t`, the DEEPER edge in `q_b`
      !! (bottom-up). The h4 interior edges reuse the White & Adcroft
      !! (2008) non-uniform stencil (matching `remap_column_ppm_h4` Step 1
      !! to round-off).
      !!
      !! NOTE the PGF integrand on top of these edges is PARABOLIC (see
      !! `boole_dpa_intz_layer`): q6 = 3*(2*q_mean - (q_t + q_b)) is the
      !! in-layer curvature — why PPM differs from PLM at the density
      !! integral even at identical edge values.
      !!
      !! Boundary layers (k=1, k=nz) -> `boundary_edges_linear`.  The pair
      !! is symmetric about the layer mean, so `q6 = 3*(2*q - (q_t+q_b))`
      !! is identically zero there: the boundary layer carries a straight
      !! line, which is the exact profile whenever `q(z)` is linear.
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz)
      real(wp), intent(in)  :: q(nz)
      real(wp), intent(out) :: q_t(nz)
      real(wp), intent(out) :: q_b(nz)

      real(wp) :: edge(NZ_STACK_MAX)
         !! edge(k) = value at the SHALLOWER interface of layer k (between
         !! layer k and k+1) = the DEEPER interface of layer k+1.  Defined
         !! for k = 1 .. nz-1.
      real(wp) :: q_lo, q_hi, ql, qr, qm, dq, dq_l, dq_r, q6
      real(wp) :: h0, h1, h2, h3, hf, h_sum
      real(wp) :: h01, h12, h23, h012, h123, h0123
      real(wp) :: f1, f2, f3, et1, et2, et3
      real(wp), parameter :: H_NEGLECT = 1.0e-30_wp
      real(wp), parameter :: H_MIN_FRAC = 1.0e-5_wp
      integer  :: k

      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if
      if (nz == 2) then
         call boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
         call boundary_edges_linear(h(2), h(1), q(2), q(2) - q(1), q_t(2), q_b(2))
         return
      end if

      ! ---- Interior interface estimates (implicit-h4 explicit form) ----
      ! edge(k) sits between layer k (deeper) and k+1 (shallower); the
      ! four-cell stencil is h(k-1..k+2).
      do k = 2, nz - 2
         h0 = h(k - 1)
         h1 = h(k)
         h2 = h(k + 1)
         h3 = h(k + 2)
         h_sum = h0 + h1 + h2 + h3
         if (h0 + h1 <= 0.0_wp .or. h1 + h2 <= 0.0_wp .or. h2 + h3 <= 0.0_wp) then
            hf = H_MIN_FRAC*max(H_NEGLECT, h_sum)
            h0 = max(h0, hf)
            h1 = max(h1, hf)
            h2 = max(h2, hf)
            h3 = max(h3, hf)
         end if
         h01 = h0 + h1
         h12 = h1 + h2
         h23 = h2 + h3
         h012 = h0 + h1 + h2
         h123 = h1 + h2 + h3
         h0123 = h0 + h1 + h2 + h3
         f1 = h01*h23/h12
         f2 = h2*q(k) + h1*q(k + 1)
         f3 = 1.0_wp/h012 + 1.0_wp/h123
         et1 = f1*f2*f3
         et2 = (h2*h23/(h012*h01))*((h0 + 2.0_wp*h1)*q(k) - h1*q(k - 1))
         et3 = (h1*h01/(h123*h23))*((2.0_wp*h2 + h3)*q(k + 1) - h2*q(k + 2))
         edge(k) = (et1 + et2 + et3)/h0123
      end do
      ! Near-boundary interior interfaces: thickness-weighted (h2)
      ! interface estimate — the linear-exact value at the shared face of
      ! two piecewise-linear cells, q_f = (q_k h_{k+1} + q_{k+1} h_k) /
      ! (h_k + h_{k+1}).  These touch interior layer 2 (interface 1|2) and
      ! nz-1 (interface nz-1|nz); the bounding layers 1 and nz take the
      ! one-sided linear pair, so only these two matter.  (A plain mean
      ! biases the thick interior layer's edge on non-uniform thicknesses.)
      edge(1) = (q(1)*h(2) + q(2)*h(1))/(h(1) + h(2))
      edge(nz - 1) = (q(nz - 1)*h(nz) + q(nz)*h(nz - 1))/(h(nz - 1) + h(nz))

      ! ---- Per-layer edges + Colella-Woodward parabola limiter ----
      ! Boundary layers: linear-exact one-sided pair (q6 == 0 there).
      call boundary_edges_linear(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call boundary_edges_linear(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                                 q_t(nz), q_b(nz))
      do k = 2, nz - 1
         qm = q(k)
         ql = edge(k - 1)   ! deeper interface  -> bottom edge
         qr = edge(k)       ! shallower interface -> top edge
         ! Clip both edges into the local monotone bounds [min,max] of the
         ! three adjacent means.
         q_lo = min(q(k - 1), q(k), q(k + 1))
         q_hi = max(q(k - 1), q(k), q(k + 1))
         ql = max(q_lo, min(q_hi, ql))
         qr = max(q_lo, min(q_hi, qr))
         dq = qr - ql
         dq_l = qm - ql
         dq_r = qr - qm
         if (dq_l*dq_r <= 0.0_wp) then
            ! Local extremum -> flatten to PCM.
            ql = qm
            qr = qm
         else
            q6 = 6.0_wp*qm - 3.0_wp*(ql + qr)
            if (abs(q6) > abs(dq)) then
               if (q6*dq > 0.0_wp) then
                  ql = 3.0_wp*qm - 2.0_wp*qr
               else
                  qr = 3.0_wp*qm - 2.0_wp*ql
               end if
            end if
         end if
         q_b(k) = ql     ! deeper edge
         q_t(k) = qr     ! shallower edge
      end do
   end subroutine ppm_edges_column

   pure subroutine boole_dpa_intz_layer(eos, rho0, rho_ref, &
                                        e_top, dz, &
                                        t_t, t_b, t_mean, &
                                        s_t, s_b, s_mean, &
                                        parabolic, dpa, intz_dpa)
      !$acc routine seq
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
                                  parabolic, dpa_face)
      !$acc routine seq
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
      !! COST: 5 sub-columns x 5 sub-points = 25 EOS evaluations per face
      !! per layer, against 0 for the trapezoid.  `reconstruct_for_pressure`
      !! is opt-in and already the expensive branch.
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
      logical, intent(in)  :: parabolic
         !! .true. -> the sub-column profiles carry the PPM curvature.
      real(wp), intent(out) :: dpa_face
         !! Along-face mean of `g * int rho' dz` over the layer (Pa).

      real(wp) :: wr, wl, dpa_m, intz_m, acc
      integer  :: m
      real(wp), parameter :: BOOLE_W(N_BOOLE) = &
                             [7.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp]

      acc = 0.0_wp
      do m = 1, N_BOOLE
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

end module rdb_ocean_pgf_reconstruct
