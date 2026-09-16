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
   !! Boundary layers (k=1, k=nz) fall back to PCM edges
   !! (`Q_t = Q_b = Q`) — no boundary extrapolation, matching the
   !! conservative default.
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

   pure subroutine plm_edges_column(nz, h, q, q_t, q_b)
      !$acc routine seq
      !! Per-column PLM top/bottom edge values of a layer-mean field `q`,
      !! via a two-stage h-weighted van-Leer slope (White, Adcroft &
      !! Hallberg 2009 §2). Returns the SHALLOWER edge in `q_t` (toward
      !! k+1) and the DEEPER edge in `q_b` (toward k-1), bottom-up.
      !! Boundary layers (k=1, k=nz) -> PCM edges (q_t=q_b=q).
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

      ! Single / two-layer columns: all PCM (every layer is a boundary).
      if (nz <= 2) then
         do k = 1, nz
            q_t(k) = q(k)
            q_b(k) = q(k)
         end do
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
      q_t(1) = q(1)
      q_b(1) = q(1)
      q_t(nz) = q(nz)
      q_b(nz) = q(nz)
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
      !! Boundary layers (k=1, k=nz) -> PCM edges.
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

      if (nz <= 2) then
         do k = 1, nz
            q_t(k) = q(k)
            q_b(k) = q(k)
         end do
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
      ! nz-1 (interface nz-1|nz); the bounding layers 1 and nz fall back to
      ! PCM, so only these two matter.  (A plain mean biases the thick
      ! interior layer's edge on non-uniform thicknesses.)
      edge(1) = (q(1)*h(2) + q(2)*h(1))/(h(1) + h(2))
      edge(nz - 1) = (q(nz - 1)*h(nz) + q(nz)*h(nz - 1))/(h(nz - 1) + h(nz))

      ! ---- Per-layer edges + Colella-Woodward parabola limiter ----
      ! Boundary layers: PCM.
      q_t(1) = q(1)
      q_b(1) = q(1)
      q_t(nz) = q(nz)
      q_b(nz) = q(nz)
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

end module rdb_ocean_pgf_reconstruct
