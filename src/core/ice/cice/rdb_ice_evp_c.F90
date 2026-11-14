module rdb_ice_evp_c
   !! Sea-ice C-grid elastic-viscous-plastic dynamics: the subcycle kernels,
   !! the strain-rate diagnostics (`divu`, `shear`, `vort`, `rdg_conv`,
   !! `rdg_shear`) and the driver `ice_evp_c_step`.
   !!
   !! Metric terms of the C-grid strain rates: Hunke and Dukowicz (2002).
   !! Rheology: Hibler (1979) viscous-plastic law with the elliptical yield
   !! curve and plastic potential, Konig Beatty and Holland (2010) tensile
   !! strength (`Ktens`), regularised by `deltaminEVP` (Hunke 2001) and the
   !! Hibler or Kreyscher et al. (2000) capping of the viscosities; elastic
   !! relaxation of the stresses following Hunke and Dukowicz (1997), in the
   !! classic form of Hunke (2001) or the revised form of Bouillon et al.
   !! (2013) (`revised_evp`).  Discretisation: Arakawa C grid, with the
   !! normal stresses (`sigma_1 = sigma_11 + sigma_22`, `sigma_2 = sigma_11 -
   !! sigma_22`) at T and the shear stress `sigma_12` at the corners (Kimmritz
   !! et al. 2016; Lemieux et al. 2024).  The momentum equation is solved for
   !! `u` at E and `v` at N, the Coriolis term using the other component
   !! interpolated to the point, the ice-ocean drag semi-implicitly in the
   !! new velocity (Hunke and Dukowicz 1997, sec. 3); no seabed stress.
   !!
   !! Roundabout port of the CICE6 `grid_ice = 'C'` dynamics (BSD-3, see
   !! `NOTICE`).  The index convention of `rdb_ice_evp_c_prep` applies
   !! (`E(i,j) = u(i+1,j)`, `N(i,j) = v(i,j+1)`, `U(i,j) = stress12U(i+1,j+1)`).
   !!
   !! **One subcycle** (4 launches, all `do concurrent`, explicit-shape
   !! dummies, scalars by value, no host data motion):
   !!
   !! | launch | kernel | writes |
   !! |--------|--------|--------|
   !! | K1 | corner strain rate (`shearU`, `deltaU` for `avg_strength`) | `shearU`, `deltaU` |
   !! | K2 | T strain rates, viscosities, replacement pressure, T stresses | `zetax2T`, `etax2T`, `stresspT`, `stressmT` |
   !! | K3 | corner viscosity (T average or strength average), corner stress | `stress12U` |
   !! | K4 | stress divergence + momentum at E and N | `un`, `vn` |
   !!
   !! The velocity is PING-PONGED between two face pairs: every kernel reads the
   !! old pair `(uo, vo)` and K4 writes the new pair `(un, vn)` at the ice faces
   !! only.  The tangential and corner velocity components (formerly a fifth
   !! launch, `evp_average_impl`, and four arrays) are interpolated INLINE, with
   !! the same expressions, by K1 (corner) and K4 (face), so the old pair is
   !! the only velocity read.  The ice masks are packed in one byte per point
   !! (`imask`, built once per step), the ice-free corners' `shearU`/`deltaU`
   !! are zeroed once per step in the preparation, and the kernels skip every
   !! point off its mask.  A step is the preparation, 2 face copies, `4 ndte`
   !! launches (+ the wraps), the averages once and the finish kernel.
   !!
   !! **Halo (multi-rank PR).**  Between K4 and the next K1 of EVERY subcycle,
   !! and once after `ice_evp_c_prep`, the face pair holding the latest velocity
   !! must be refreshed with `ice_halo_faces` (`rdb_ice_cgrid_metrics`: exchange
   !! on split axes, periodic wrap, tripolar fold) -- nothing else is exchanged:
   !! the T stresses and viscosities and the corner shear and stress are
   !! recomputed in the ghost rings (width 3, see `rdb_ice_evp_c_prep`) from the
   !! haloed faces.  Under the fold every stress component is invariant (a 180
   !! degree rotation flips both axes), so the ghost value computed from the
   !! folded velocities is the rotated image.  `ice_evp_c_step` shows the
   !! sequence (wrap, `ice_evp_c_copy_faces`, then `ice_evp_c_substep` with the
   !! two pairs swapped each subcycle, a wrap of the new pair after each, and a
   !! copy back when `ndte` is odd); a multi-rank caller replaces the wrap by
   !! the halo exchange.  `ice_evp_c_average` (the interpolation, now off the
   !! subcycle path) runs once after the loop because `ice_evp_c_finish` reads
   !! its `uvelN`/`vvelE`.  Single rank periodic: `wrap_x`, `wrap_y`.
   !!
   !! **Errors.**  A face speed that is non-finite or above `vel_cap` raises
   !! `ICE_ERR_DYN_VEL` in the error word of the owning cell, after the loop;
   !! nothing is clamped.  Non-finite input propagates (every max/min on the
   !! viscosity path is `ieee_is_finite`-guarded so a NaN is never turned into
   !! a plausible number).
   use, intrinsic :: iso_fortran_env, only: int8, int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_DYN_VEL
   use rdb_ice_params, only: ice_constants_t, ice_ridge_params_t
   use rdb_ice_evp_c_prep, only: ice_evp_c_const_t, ice_evp_c_work_t, ice_evp_c_prep, EVP_BIT_T, EVP_BIT_U, &
                                 EVP_BIT_E, EVP_BIT_N
   implicit none
   private

   public :: ice_evp_c_step, ice_evp_c_substep, ice_evp_c_average, ice_evp_c_finish
   public :: ice_evp_c_wrap_faces, ice_evp_c_div_stress, ice_evp_c_copy_faces
   public :: evp_strain_u_impl, evp_stress_t_impl, evp_stress_u_impl, evp_momentum_impl
   public :: evp_average_impl, evp_copy_faces_impl, evp_finish_impl, evp_div_stress_impl, evp_wrap_u_impl, evp_wrap_v_impl
   public :: evp_div_x, evp_div_y, evp_visc_replpress

   real(wp), parameter :: EPS13 = 1.0e-13_wp

contains

   ! ------------------------------------------------------------------
   ! scalar helpers (shared by the kernels and the diagnostic routines)
   ! ------------------------------------------------------------------

   pure subroutine evp_visc_replpress(strength, dmin, delta, cap_max, ktens, epp2i, zetax2, etax2, rep_prs)
      !! Viscosities and replacement pressure of the elliptical yield curve
      !! (Hibler 1979; Konig Beatty and Holland 2010; Kreyscher et al. 2000):
      !! `tmp = P / max(Delta, Dmin)` (Hibler cap) or `P / (Delta + Dmin)`;
      !! `zeta2 = (1 + Ktens) tmp`, `P_rep = (1 - Ktens) tmp Delta`,
      !! `eta2 = zeta2 / e_plasticpot^2`.  A non-finite `Delta` is not
      !! capped (the cap would turn NaN into `Dmin`); it propagates.
      !$acc routine seq
      real(wp), value :: strength, dmin, delta, ktens, epp2i
      logical, value :: cap_max
      real(wp), intent(out) :: zetax2, etax2, rep_prs
      real(wp) :: tmp, dd

      if (cap_max) then
         if (ieee_is_finite(delta)) then
            dd = max(delta, dmin)
         else
            dd = delta
         end if
         tmp = strength/dd
      else
         tmp = strength/(delta + dmin)
      end if
      zetax2 = (1.0_wp + ktens)*tmp
      rep_prs = (1.0_wp - ktens)*tmp*delta
      etax2 = epp2i*zetax2
   end subroutine evp_visc_replpress

   pure function evp_div_x(dxE, dyE, dyT1, dyT0, dxU1, dxU0, sp1, sp0, sm1, sm0, s121, s120, arear, rheo) &
      result(f)
      !! x stress divergence at E(i,j) [N/m2] (Lemieux et al. 2024):
      !! `T(i+1)` is index 1, `T(i)` index 0, `U(i,j)` index 1, `U(i,j-1)` 0.
      !$acc routine seq
      real(wp), value :: dxE, dyE, dyT1, dyT0, dxU1, dxU0, sp1, sp0, sm1, sm0, s121, s120, arear, rheo
      real(wp) :: f
      f = rheo*arear*(0.5_wp*dyE*(sp1 - sp0) &
                      + (0.5_wp/dyE)*(dyT1*dyT1*sm1 - dyT0*dyT0*sm0) &
                      + (1.0_wp/dxE)*(dxU1*dxU1*s121 - dxU0*dxU0*s120))
   end function evp_div_x

   pure function evp_div_y(dxN, dyN, dxT1, dxT0, dyU1, dyU0, sp1, sp0, sm1, sm0, s121, s120, arear, rheo) &
      result(f)
      !! y stress divergence at N(i,j) [N/m2]: `T(j+1)` index 1, `T(j)` 0,
      !! `U(i,j)` 1, `U(i-1,j)` 0.
      !$acc routine seq
      real(wp), value :: dxN, dyN, dxT1, dxT0, dyU1, dyU0, sp1, sp0, sm1, sm0, s121, s120, arear, rheo
      real(wp) :: f
      f = rheo*arear*(0.5_wp*dxN*(sp1 - sp0) &
                      - (0.5_wp/dxN)*(dxT1*dxT1*sm1 - dxT0*dxT0*sm0) &
                      + (1.0_wp/dyN)*(dyU1*dyU1*s121 - dyU0*dyU0*s120))
   end function evp_div_y

   pure function evp_wavg4(a1, a2, a3, a4, w1, w2, w3, w4) result(s)
      !! Area-weighted mean of four face values (0 when the weights sum to 0):
      !! the tangential velocity interpolation of the C grid, written once so the
      !! kernels that need it inline compute it with the same operands, same order.
      !$acc routine seq
      real(wp), value :: a1, a2, a3, a4, w1, w2, w3, w4
      real(wp) :: s
      real(wp) :: wt
      wt = w1 + w2 + w3 + w4
      s = 0.0_wp
      if (wt /= 0.0_wp) s = (a1*w1 + a2*w2 + a3*w3 + a4*w4)/wt
   end function evp_wavg4

   pure function evp_wavg2(a1, a2, w1, w2) result(s)
      !! Two-point version of `evp_wavg4` (corner velocities).
      !$acc routine seq
      real(wp), value :: a1, a2, w1, w2
      real(wp) :: s
      real(wp) :: wt
      wt = w1 + w2
      s = 0.0_wp
      if (wt /= 0.0_wp) s = (a1*w1 + a2*w2)/wt
   end function evp_wavg2

   ! ------------------------------------------------------------------
   ! K1: corner strain rate
   ! ------------------------------------------------------------------

   subroutine evp_strain_u_impl(nx, ny, is, ie, js, je, e_factor, avg_strength, u, v, earea, narea, &
                                dxE, dyN, dxU, dyU, ratiodxN, ratiodxNr, ratiodyE, ratiodyEr, epm, npm, uvm, &
                                imask, shearU, deltaU)
      !! K1.  Shear strain rate x corner area, `shearU`, at
      !! `U(is-2:ie+1, js-2:je+1)` where the corner carries ice (the ice-free
      !! corners were zeroed once in the preparation and are skipped);
      !! with `avg_strength` also `deltaU = sqrt(div^2 + e_factor (tension^2
      !! + shear^2))`.  The tangential velocities at the corner (and, for
      !! `deltaU`, at the four surrounding faces) are interpolated here from
      !! `(u, v)` with the same expressions as `evp_average_impl`, so no
      !! averaged array is read.  The C-grid lateral boundary condition is the
      !! `ratiodxN` / `ratiodyE` extrapolation of the tangential face
      !! velocities (a closed face contributes the mirrored (sign-reversed)
      !! value of the open one, so the wall-parallel velocity vanishes at a
      !! straight wall: NO slip, Hunke and Dukowicz (2002) / Lemieux et al.
      !! (2024)).  It acts at the corners where at least two of the four T
      !! cells are wet.
      integer, value :: nx, ny, is, ie, js, je
      real(wp), value :: e_factor
      logical, value :: avg_strength
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny)
      real(wp), intent(in) :: dxE(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: ratiodxN(nx, ny), ratiodxNr(nx, ny), ratiodyE(nx, ny), ratiodyEr(nx, ny)
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny), uvm(nx, ny)
      integer(int8), intent(in) :: imask(nx, ny)
      real(wp), intent(inout) :: shearU(nx, ny), deltaU(nx, ny)
      integer :: i, j
      real(wp) :: uEijp1, uEij, vNip1j, vNij, uNip1j, uNij, vEijp1, vEij, shr, dvg, ten, uvelU, vvelU
      real(wp) :: uvelNij, uvelNip1j, vvelEij, vvelEjp1

      do concurrent(j=js - 2:je + 1, i=is - 2:ie + 1) local(uEijp1, uEij, vNip1j, vNij, uNip1j, uNij, &
                                                            vEijp1, vEij, shr, dvg, ten, uvelU, vvelU, &
                                                            uvelNij, uvelNip1j, vvelEij, vvelEjp1)
         if (iand(imask(i, j), EVP_BIT_U) /= 0_int8) then
            uvelU = evp_wavg2(u(i + 1, j), u(i + 1, j + 1), earea(i, j), earea(i, j + 1))*uvm(i, j)
            vvelU = evp_wavg2(v(i, j + 1), v(i + 1, j + 1), narea(i, j), narea(i + 1, j))*uvm(i, j)
            uEijp1 = u(i + 1, j + 1)*epm(i, j + 1) &
                     + (epm(i, j) - epm(i, j + 1))*epm(i, j)*ratiodyE(i, j)*u(i + 1, j)
            uEij = u(i + 1, j)*epm(i, j) &
                   + (epm(i, j + 1) - epm(i, j))*epm(i, j + 1)*ratiodyEr(i, j)*u(i + 1, j + 1)
            vNip1j = v(i + 1, j + 1)*npm(i + 1, j) &
                     + (npm(i, j) - npm(i + 1, j))*npm(i, j)*ratiodxN(i, j)*v(i, j + 1)
            vNij = v(i, j + 1)*npm(i, j) &
                   + (npm(i + 1, j) - npm(i, j))*npm(i + 1, j)*ratiodxNr(i, j)*v(i + 1, j + 1)
            shr = dxU(i, j)*(uEijp1 - uEij) - uvelU*(dxE(i, j + 1) - dxE(i, j)) &
                  + dyU(i, j)*(vNip1j - vNij) - vvelU*(dyN(i + 1, j) - dyN(i, j))
            shearU(i, j) = shr
            if (avg_strength) then
               uvelNij = evp_wavg4(u(i, j), u(i + 1, j), u(i, j + 1), u(i + 1, j + 1), &
                                   earea(i - 1, j), earea(i, j), earea(i - 1, j + 1), earea(i, j + 1))*npm(i, j)
               uvelNip1j = evp_wavg4(u(i + 1, j), u(i + 2, j), u(i + 1, j + 1), u(i + 2, j + 1), &
                                     earea(i, j), earea(i + 1, j), earea(i, j + 1), earea(i + 1, j + 1))*npm(i + 1, j)
               vvelEij = evp_wavg4(v(i, j), v(i + 1, j), v(i, j + 1), v(i + 1, j + 1), &
                                   narea(i, j - 1), narea(i + 1, j - 1), narea(i, j), narea(i + 1, j))*epm(i, j)
               vvelEjp1 = evp_wavg4(v(i, j + 1), v(i + 1, j + 1), v(i, j + 2), v(i + 1, j + 2), &
                                    narea(i, j), narea(i + 1, j), narea(i, j + 1), narea(i + 1, j + 1))*epm(i, j + 1)
               uNip1j = uvelNip1j*npm(i + 1, j) &
                        + (npm(i, j) - npm(i + 1, j))*npm(i, j)*ratiodxN(i, j)*uvelNij
               uNij = uvelNij*npm(i, j) &
                      + (npm(i + 1, j) - npm(i, j))*npm(i + 1, j)*ratiodxNr(i, j)*uvelNip1j
               vEijp1 = vvelEjp1*epm(i, j + 1) &
                        + (epm(i, j) - epm(i, j + 1))*epm(i, j)*ratiodyE(i, j)*vvelEij
               vEij = vvelEij*epm(i, j) &
                      + (epm(i, j + 1) - epm(i, j))*epm(i, j + 1)*ratiodyEr(i, j)*vvelEjp1
               dvg = dyU(i, j)*(uNip1j - uNij) + uvelU*(dyN(i + 1, j) - dyN(i, j)) &
                     + dxU(i, j)*(vEijp1 - vEij) + vvelU*(dxE(i, j + 1) - dxE(i, j))
               ten = dyU(i, j)*(uNip1j - uNij) - uvelU*(dyN(i + 1, j) - dyN(i, j)) &
                     - dxU(i, j)*(vEijp1 - vEij) + vvelU*(dxE(i, j + 1) - dxE(i, j))
               deltaU(i, j) = sqrt(dvg*dvg + e_factor*(ten*ten + shr*shr))
            end if
         end if
      end do
   end subroutine evp_strain_u_impl

   ! ------------------------------------------------------------------
   ! K2: T viscosities and stresses
   ! ------------------------------------------------------------------

   subroutine evp_stress_t_impl(nx, ny, is, ie, js, je, arlx1i, revp, denom1, epp2i, e_factor, ktens, deltamin, &
                                cap_max, u, v, dxN, dyE, dxT, dyT, uarea, tarea, strength, shearU, imask, &
                                zetax2T, etax2T, stresspT, stressmT)
      !! K2.  Divergence and tension at T from the faces, the shear squared
      !! as the area-weighted mean of the four corner values, Delta, the
      !! viscosities and replacement pressure, and the elastic update of
      !! `stresspT`/`stressmT` (Hunke and Dukowicz 1997; Bouillon et al. 2013)
      !! over `T(is-1:ie+1, js-1:je+1)` where `iceTmask`.
      integer, value :: nx, ny, is, ie, js, je
      real(wp), value :: arlx1i, revp, denom1, epp2i, e_factor, ktens, deltamin
      logical, value :: cap_max
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: dxN(nx, ny), dyE(nx, ny), dxT(nx, ny), dyT(nx, ny)
      real(wp), intent(in) :: uarea(nx, ny), tarea(nx, ny), strength(nx, ny), shearU(nx, ny)
      integer(int8), intent(in) :: imask(nx, ny)
      real(wp), intent(inout) :: zetax2T(nx, ny), etax2T(nx, ny), stresspT(nx, ny), stressmT(nx, ny)
      integer :: i, j
      real(wp) :: uE1, uE0, vN1, vN0, divT, tenT, w1, w2, w3, w4, ssq, delta, zeta2, eta2, rep, keep

      keep = 1.0_wp - arlx1i*revp
      do concurrent(j=js - 1:je + 1, i=is - 1:ie + 1) local(uE1, uE0, vN1, vN0, divT, tenT, w1, w2, w3, w4, &
                                                            ssq, delta, zeta2, eta2, rep)
         if (iand(imask(i, j), EVP_BIT_T) /= 0_int8) then
            uE1 = u(i + 1, j)
            uE0 = u(i, j)
            vN1 = v(i, j + 1)
            vN0 = v(i, j)
            divT = dyE(i, j)*uE1 - dyE(i - 1, j)*uE0 + dxN(i, j)*vN1 - dxN(i, j - 1)*vN0
            tenT = dyT(i, j)*dyT(i, j)*(uE1/dyE(i, j) - uE0/dyE(i - 1, j)) &
                   - dxT(i, j)*dxT(i, j)*(vN1/dxN(i, j) - vN0/dxN(i, j - 1))
            w1 = uarea(i, j)
            w2 = uarea(i, j - 1)
            w3 = uarea(i - 1, j - 1)
            w4 = uarea(i - 1, j)
            ssq = (shearU(i, j)*shearU(i, j)*w1 + shearU(i, j - 1)*shearU(i, j - 1)*w2 &
                   + shearU(i - 1, j - 1)*shearU(i - 1, j - 1)*w3 + shearU(i - 1, j)*shearU(i - 1, j)*w4) &
                  /(w1 + w2 + w3 + w4)
            delta = sqrt(divT*divT + e_factor*(tenT*tenT + ssq))
            call evp_visc_replpress(strength(i, j), deltamin*tarea(i, j), delta, cap_max, ktens, epp2i, &
                                    zeta2, eta2, rep)
            zetax2T(i, j) = zeta2
            etax2T(i, j) = eta2
            stresspT(i, j) = (stresspT(i, j)*keep + arlx1i*(zeta2*divT - rep))*denom1
            stressmT(i, j) = (stressmT(i, j)*keep + arlx1i*eta2*tenT)*denom1
         end if
      end do
   end subroutine evp_stress_t_impl

   ! ------------------------------------------------------------------
   ! K3: corner stress
   ! ------------------------------------------------------------------

   subroutine evp_stress_u_impl(nx, ny, is, ie, js, je, arlx1i, revp, denom1, epp2i, ktens, deltamin, cap_max, &
                                avg_strength, hm, tarea, uarea, etax2T, strengthU, deltaU, shearU, imask, &
                                stress12U)
      !! K3.  Corner shear stress over `U(is-1:ie, js-1:je)` where
      !! `iceUmask`.  `avg_zeta` (CICE default): the shear viscosity at the
      !! corner is the masked, area-weighted mean of the four T viscosities;
      !! `avg_strength`: it is rebuilt from the averaged strength and the
      !! corner Delta (`deltaminEVP * uarea` regularisation).
      integer, value :: nx, ny, is, ie, js, je
      real(wp), value :: arlx1i, revp, denom1, epp2i, ktens, deltamin
      logical, value :: cap_max, avg_strength
      real(wp), intent(in) :: hm(nx, ny), tarea(nx, ny), uarea(nx, ny), etax2T(nx, ny)
      real(wp), intent(in) :: strengthU(nx, ny), deltaU(nx, ny), shearU(nx, ny)
      integer(int8), intent(in) :: imask(nx, ny)
      real(wp), intent(inout) :: stress12U(nx + 1, ny + 1)
      integer :: i, j
      real(wp) :: w1, w2, w3, w4, wt, etaU, zeta2, rep, keep

      keep = 1.0_wp - arlx1i*revp
      do concurrent(j=js - 1:je, i=is - 1:ie) local(w1, w2, w3, w4, wt, etaU, zeta2, rep)
         if (iand(imask(i, j), EVP_BIT_U) /= 0_int8) then
            if (avg_strength) then
               call evp_visc_replpress(strengthU(i, j), deltamin*uarea(i, j), deltaU(i, j), cap_max, ktens, &
                                       epp2i, zeta2, etaU, rep)
            else
               w1 = hm(i, j)*tarea(i, j)
               w2 = hm(i + 1, j)*tarea(i + 1, j)
               w3 = hm(i, j + 1)*tarea(i, j + 1)
               w4 = hm(i + 1, j + 1)*tarea(i + 1, j + 1)
               wt = w1 + w2 + w3 + w4
               etaU = 0.0_wp
               if (wt /= 0.0_wp) then
                  etaU = (w1*etax2T(i, j) + w2*etax2T(i + 1, j) + w3*etax2T(i, j + 1) &
                          + w4*etax2T(i + 1, j + 1))/wt
               end if
            end if
            stress12U(i + 1, j + 1) = (stress12U(i + 1, j + 1)*keep + arlx1i*0.5_wp*etaU*shearU(i, j))*denom1
         end if
      end do
   end subroutine evp_stress_u_impl

   ! ------------------------------------------------------------------
   ! K4: stress divergence + momentum
   ! ------------------------------------------------------------------

   subroutine evp_momentum_impl(nx, ny, is, ie, js, je, revp, brlx, rhow, cw, cosw, sinw, free_drift, rheo_amin, &
                                uo, vo, un, vn, earea, narea, epm, npm, dxE, dyE, dxN, dyN, dxU, dyU, dxT, dyT, &
                                earear, narear, stresspT, stressmT, stress12U, imask, aiE, aiN, emassdti, &
                                nmassdti, fmE, fmN, forcexE, forceyN, uinit, vinit, uocnE, vocnE, uocnN, vocnN)
      !! K4.  For the owned faces `E(is:ie, js:je)` and `N(is:ie, js:je)`
      !! where they carry ice: the stress divergence (a local scalar, with
      !! the rheology dropped where `free_drift` and the concentration is
      !! below `rheo_amin`) and the momentum update
      !!   `u = (ccb v_t + cc1)/cca`, `v = (-ccb u_t + cc2)/cca`,
      !!   `cca = (brlx + revp) m/dt + vrel cosw`, `ccb = f m + sgn(f) vrel sinw`,
      !!   `cc1 = strint + force + vrel water + (m/dt)(brlx u + revp u_init)`
      !! with `vrel = a rho_w C_w |u_ocn - u|` evaluated at the old velocity
      !! and `(u_t, v_t)` the tangential components interpolated from the OLD
      !! velocity here (same expressions as `evp_average_impl`).  Ping-pong: the
      !! old pair `(uo, vo)` is only read, the new pair `(un, vn)` only written at
      !! the ice faces; at every other face `(un, vn)` must already hold the old
      !! value (the caller copies the pair once per step) and the ghost ring is the
      !! caller's (wrap or halo exchange of `(un, vn)`, then swap).
      integer, value :: nx, ny, is, ie, js, je
      real(wp), value :: revp, brlx, rhow, cw, cosw, sinw, rheo_amin
      logical, value :: free_drift
      real(wp), intent(in) :: uo(nx + 1, ny), vo(nx, ny + 1)
      real(wp), intent(inout) :: un(nx + 1, ny), vn(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), epm(nx, ny), npm(nx, ny)
      real(wp), intent(in) :: dxE(nx, ny), dyE(nx, ny), dxN(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), earear(nx, ny), narear(nx, ny)
      real(wp), intent(in) :: stresspT(nx, ny), stressmT(nx, ny), stress12U(nx + 1, ny + 1)
      integer(int8), intent(in) :: imask(nx, ny)
      real(wp), intent(in) :: aiE(nx, ny), aiN(nx, ny), emassdti(nx, ny), nmassdti(nx, ny)
      real(wp), intent(in) :: fmE(nx, ny), fmN(nx, ny), forcexE(nx, ny), forceyN(nx, ny)
      real(wp), intent(in) :: uinit(nx, ny), vinit(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)
      integer :: i, j
      real(wp) :: rheo, strint, uold, vold, vrel, sg, cca, ccb, cc, wat
      integer(int8) :: mk

      do concurrent(j=js:je, i=is:ie) local(rheo, strint, uold, vold, vrel, sg, cca, ccb, cc, wat, mk)
         mk = imask(i, j)
         ! ---- u at E(i,j) ----
         if (iand(mk, EVP_BIT_E) /= 0_int8) then
            rheo = 1.0_wp
            if (free_drift .and. .not. (aiE(i, j) > rheo_amin)) rheo = 0.0_wp
            strint = evp_div_x(dxE(i, j), dyE(i, j), dyT(i + 1, j), dyT(i, j), dxU(i, j), dxU(i, j - 1), &
                               stresspT(i + 1, j), stresspT(i, j), stressmT(i + 1, j), stressmT(i, j), &
                               stress12U(i + 1, j + 1), stress12U(i + 1, j), earear(i, j), rheo)
            uold = uo(i + 1, j)
            vold = evp_wavg4(vo(i, j), vo(i + 1, j), vo(i, j + 1), vo(i + 1, j + 1), &
                             narea(i, j - 1), narea(i + 1, j - 1), narea(i, j), narea(i + 1, j))*epm(i, j)
            vrel = aiE(i, j)*rhow*cw*sqrt((uocnE(i, j) - uold)**2 + (vocnE(i, j) - vold)**2)
            sg = sign(1.0_wp, fmE(i, j))
            wat = uocnE(i, j)*cosw - vocnE(i, j)*sinw*sg
            cca = (brlx + revp)*emassdti(i, j) + vrel*cosw
            ccb = fmE(i, j) + sg*vrel*sinw
            cc = strint + forcexE(i, j) + vrel*wat + emassdti(i, j)*(brlx*uold + revp*uinit(i, j))
            un(i + 1, j) = (ccb*vold + cc)/cca
         end if
         ! ---- v at N(i,j) ----
         if (iand(mk, EVP_BIT_N) /= 0_int8) then
            rheo = 1.0_wp
            if (free_drift .and. .not. (aiN(i, j) > rheo_amin)) rheo = 0.0_wp
            strint = evp_div_y(dxN(i, j), dyN(i, j), dxT(i, j + 1), dxT(i, j), dyU(i, j), dyU(i - 1, j), &
                               stresspT(i, j + 1), stresspT(i, j), stressmT(i, j + 1), stressmT(i, j), &
                               stress12U(i + 1, j + 1), stress12U(i, j + 1), narear(i, j), rheo)
            uold = evp_wavg4(uo(i, j), uo(i + 1, j), uo(i, j + 1), uo(i + 1, j + 1), &
                             earea(i - 1, j), earea(i, j), earea(i - 1, j + 1), earea(i, j + 1))*npm(i, j)
            vold = vo(i, j + 1)
            vrel = aiN(i, j)*rhow*cw*sqrt((uocnN(i, j) - uold)**2 + (vocnN(i, j) - vold)**2)
            sg = sign(1.0_wp, fmN(i, j))
            wat = vocnN(i, j)*cosw + uocnN(i, j)*sinw*sg
            cca = (brlx + revp)*nmassdti(i, j) + vrel*cosw
            ccb = fmN(i, j) + sg*vrel*sinw
            cc = strint + forceyN(i, j) + vrel*wat + nmassdti(i, j)*(brlx*vold + revp*vinit(i, j))
            vn(i, j + 1) = (-ccb*uold + cc)/cca
         end if
      end do
   end subroutine evp_momentum_impl

   subroutine evp_div_stress_impl(nx, ny, is, ie, js, je, free_drift, rheo_amin, dxE, dyE, dxN, dyN, dxU, dyU, &
                                  dxT, dyT, earear, narear, stresspT, stressmT, stress12U, iceEmask, iceNmask, &
                                  aiE, aiN, strintxE, strintyN)
      !! The stress divergence of K4 written to arrays (diagnostics and
      !! tests; the solver never calls it).  Same helpers as K4.
      integer, intent(in) :: nx, ny, is, ie, js, je
      logical, intent(in) :: free_drift
      real(wp), intent(in) :: rheo_amin
      real(wp), intent(in) :: dxE(nx, ny), dyE(nx, ny), dxN(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), earear(nx, ny), narear(nx, ny)
      real(wp), intent(in) :: stresspT(nx, ny), stressmT(nx, ny), stress12U(nx + 1, ny + 1)
      logical, intent(in) :: iceEmask(nx, ny), iceNmask(nx, ny)
      real(wp), intent(in) :: aiE(nx, ny), aiN(nx, ny)
      real(wp), intent(inout) :: strintxE(nx, ny), strintyN(nx, ny)
      integer :: i, j
      real(wp) :: rheo

      do concurrent(j=js:je, i=is:ie) local(rheo)
         strintxE(i, j) = 0.0_wp
         strintyN(i, j) = 0.0_wp
         if (iceEmask(i, j)) then
            rheo = 1.0_wp
            if (free_drift .and. .not. (aiE(i, j) > rheo_amin)) rheo = 0.0_wp
            strintxE(i, j) = evp_div_x(dxE(i, j), dyE(i, j), dyT(i + 1, j), dyT(i, j), dxU(i, j), dxU(i, j - 1), &
                                       stresspT(i + 1, j), stresspT(i, j), stressmT(i + 1, j), stressmT(i, j), &
                                       stress12U(i + 1, j + 1), stress12U(i + 1, j), earear(i, j), rheo)
         end if
         if (iceNmask(i, j)) then
            rheo = 1.0_wp
            if (free_drift .and. .not. (aiN(i, j) > rheo_amin)) rheo = 0.0_wp
            strintyN(i, j) = evp_div_y(dxN(i, j), dyN(i, j), dxT(i, j + 1), dxT(i, j), dyU(i, j), dyU(i - 1, j), &
                                       stresspT(i, j + 1), stresspT(i, j), stressmT(i, j + 1), stressmT(i, j), &
                                       stress12U(i + 1, j + 1), stress12U(i, j + 1), narear(i, j), rheo)
         end if
      end do
   end subroutine evp_div_stress_impl

   ! ------------------------------------------------------------------
   ! K5: averages (once per step, for the finish kernel)
   ! ------------------------------------------------------------------

   subroutine evp_average_impl(nx, ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm, uvelN, vvelE, uvelU, vvelU)
      !! K5.  The tangential and corner velocity components, masked and
      !! area-weighted (`earea` for `u`, `narea` for `v`), over
      !! `(is-2:ie+2, js-2:je+2)`:
      !!   `uvelN(i,j)` = mean of `E(i-1:i, j:j+1)`, times `npm`;
      !!   `vvelE(i,j)` = mean of `N(i:i+1, j-1:j)`, times `epm`;
      !!   `uvelU(i,j)` = mean of `E(i, j:j+1)`, `vvelU(i,j)` = mean of
      !!   `N(i:i+1, j)`, times `uvm`.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), epm(nx, ny), npm(nx, ny), uvm(nx, ny)
      real(wp), intent(inout) :: uvelN(nx, ny), vvelE(nx, ny), uvelU(nx, ny), vvelU(nx, ny)
      integer :: i, j
      real(wp) :: w1, w2, w3, w4, wt, s

      do concurrent(j=js - 2:je + 2, i=is - 2:ie + 2) local(w1, w2, w3, w4, wt, s)
         ! u at N: E(i-1,j), E(i,j), E(i-1,j+1), E(i,j+1)
         w1 = earea(i - 1, j)
         w2 = earea(i, j)
         w3 = earea(i - 1, j + 1)
         w4 = earea(i, j + 1)
         wt = w1 + w2 + w3 + w4
         s = 0.0_wp
         if (wt /= 0.0_wp) s = (u(i, j)*w1 + u(i + 1, j)*w2 + u(i, j + 1)*w3 + u(i + 1, j + 1)*w4)/wt
         uvelN(i, j) = s*npm(i, j)
         ! v at E: N(i,j-1), N(i+1,j-1), N(i,j), N(i+1,j)
         w1 = narea(i, j - 1)
         w2 = narea(i + 1, j - 1)
         w3 = narea(i, j)
         w4 = narea(i + 1, j)
         wt = w1 + w2 + w3 + w4
         s = 0.0_wp
         if (wt /= 0.0_wp) s = (v(i, j)*w1 + v(i + 1, j)*w2 + v(i, j + 1)*w3 + v(i + 1, j + 1)*w4)/wt
         vvelE(i, j) = s*epm(i, j)
         ! u at U: E(i,j), E(i,j+1)
         w1 = earea(i, j)
         w2 = earea(i, j + 1)
         wt = w1 + w2
         s = 0.0_wp
         if (wt /= 0.0_wp) s = (u(i + 1, j)*w1 + u(i + 1, j + 1)*w2)/wt
         uvelU(i, j) = s*uvm(i, j)
         ! v at U: N(i,j), N(i+1,j)
         w1 = narea(i, j)
         w2 = narea(i + 1, j)
         wt = w1 + w2
         s = 0.0_wp
         if (wt /= 0.0_wp) s = (v(i, j + 1)*w1 + v(i + 1, j + 1)*w2)/wt
         vvelU(i, j) = s*uvm(i, j)
      end do
   end subroutine evp_average_impl

   ! ------------------------------------------------------------------
   ! finish: deformation rates, ice-ocean stress, velocity guard
   ! ------------------------------------------------------------------

   subroutine evp_finish_impl(nx, ny, is, ie, js, je, e_factor, epp2i, rhow, cw, cosw, sinw, vel_cap, u, v, &
                              uvelN, vvelE, shearU, dxN, dyE, dxT, dyT, uarea, tarear, imask, &
                              aiE, aiN, fmE, fmN, uocnE, vocnE, uocnN, vocnN, divu, shear, vort, &
                              rdg_conv, rdg_shear, strocnxE, strocnyN, ice_err)
      !! Over the owned cells: the strain-rate diagnostics at T
      !! (`deformationsC_T`: `divu`, `shear`, `vort`, and the ridging
      !! convergence and shear `rdg_conv = -min(divu, 0)`, `rdg_shear =
      !! (Delta_pp^2/Delta - |divu|)/2`, Lemieux et al. 2016 / Flato and
      !! Hibler 1995), the ice-ocean stress on the ice at E (x) and N (y),
      !! and the velocity guard (`ICE_ERR_DYN_VEL`).  `shearU` is the one of
      !! the last subcycle, as in the reference scheme.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: e_factor, epp2i, rhow, cw, cosw, sinw, vel_cap
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1), uvelN(nx, ny), vvelE(nx, ny), shearU(nx, ny)
      real(wp), intent(in) :: dxN(nx, ny), dyE(nx, ny), dxT(nx, ny), dyT(nx, ny), uarea(nx, ny), tarear(nx, ny)
      integer(int8), intent(in) :: imask(nx, ny)
      real(wp), intent(in) :: aiE(nx, ny), aiN(nx, ny), fmE(nx, ny), fmN(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)
      real(wp), intent(inout) :: divu(nx, ny), shear(nx, ny), vort(nx, ny), rdg_conv(nx, ny), rdg_shear(nx, ny)
      real(wp), intent(inout) :: strocnxE(nx, ny), strocnyN(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j
      real(wp) :: uE1, uE0, vN1, vN0, divT, tenT, shrT, w1, w2, w3, w4, ssq, dlt, dpp2, tmp, dv, vrel, sg
      logical :: bad, nonfin

      do concurrent(j=js:je, i=is:ie) local(uE1, uE0, vN1, vN0, divT, tenT, shrT, w1, w2, w3, w4, ssq, dlt, &
                                            dpp2, tmp, dv, vrel, sg, bad, nonfin)
         ! ---- strain-rate diagnostics at T ----
         nonfin = .false.
         divu(i, j) = 0.0_wp
         shear(i, j) = 0.0_wp
         vort(i, j) = 0.0_wp
         rdg_conv(i, j) = 0.0_wp
         rdg_shear(i, j) = 0.0_wp
         if (iand(imask(i, j), EVP_BIT_T) /= 0_int8) then
            uE1 = u(i + 1, j)
            uE0 = u(i, j)
            vN1 = v(i, j + 1)
            vN0 = v(i, j)
            divT = dyE(i, j)*uE1 - dyE(i - 1, j)*uE0 + dxN(i, j)*vN1 - dxN(i, j - 1)*vN0
            tenT = dyT(i, j)*dyT(i, j)*(uE1/dyE(i, j) - uE0/dyE(i - 1, j)) &
                   - dxT(i, j)*dxT(i, j)*(vN1/dxN(i, j) - vN0/dxN(i, j - 1))
            shrT = dxT(i, j)*dxT(i, j)*(uvelN(i, j)/dxN(i, j) - uvelN(i, j - 1)/dxN(i, j - 1)) &
                   + dyT(i, j)*dyT(i, j)*(vvelE(i, j)/dyE(i, j) - vvelE(i - 1, j)/dyE(i - 1, j))
            w1 = uarea(i, j)
            w2 = uarea(i, j - 1)
            w3 = uarea(i - 1, j - 1)
            w4 = uarea(i - 1, j)
            ssq = (shearU(i, j)*shearU(i, j)*w1 + shearU(i, j - 1)*shearU(i, j - 1)*w2 &
                   + shearU(i - 1, j - 1)*shearU(i - 1, j - 1)*w3 + shearU(i - 1, j)*shearU(i - 1, j)*w4) &
                  /(w1 + w2 + w3 + w4)
            dlt = sqrt(divT*divT + e_factor*(tenT*tenT + ssq))
            dpp2 = divT*divT + epp2i*(tenT*tenT + ssq)
            dv = divT*tarear(i, j)
            divu(i, j) = dv
            if (ieee_is_finite(dlt)) then
               tmp = (dpp2/max(dlt, EPS13))*tarear(i, j)
            else
               tmp = dlt
            end if
            if (ieee_is_finite(dv)) then
               rdg_conv(i, j) = -min(dv, 0.0_wp)
            else
               rdg_conv(i, j) = dv            ! min() would launder NaN to 0 under relaxed FP
               nonfin = .true.
            end if
            rdg_shear(i, j) = 0.5_wp*(tmp - abs(dv))
            shear(i, j) = tarear(i, j)*sqrt(tenT*tenT + shrT*shrT)
            vort(i, j) = tarear(i, j)*((dyE(i, j)*vvelE(i, j) - dyE(i - 1, j)*vvelE(i - 1, j)) &
                                       - (dxN(i, j)*uvelN(i, j) - dxN(i, j - 1)*uvelN(i, j - 1)))
         end if
         ! ---- ice-ocean stress on the ice (force per unit cell area) ----
         strocnxE(i, j) = 0.0_wp
         strocnyN(i, j) = 0.0_wp
         bad = nonfin
         if (iand(imask(i, j), EVP_BIT_E) /= 0_int8) then
            vrel = aiE(i, j)*rhow*cw*sqrt((uocnE(i, j) - u(i + 1, j))**2 + (vocnE(i, j) - vvelE(i, j))**2)
            sg = sign(1.0_wp, fmE(i, j))
            strocnxE(i, j) = vrel*((uocnE(i, j) - u(i + 1, j))*cosw - (vocnE(i, j) - vvelE(i, j))*sinw*sg)
            bad = .not. ieee_is_finite(u(i + 1, j))
            if (.not. bad) bad = abs(u(i + 1, j)) > vel_cap
         end if
         if (iand(imask(i, j), EVP_BIT_N) /= 0_int8) then
            vrel = aiN(i, j)*rhow*cw*sqrt((uocnN(i, j) - uvelN(i, j))**2 + (vocnN(i, j) - v(i, j + 1))**2)
            sg = sign(1.0_wp, fmN(i, j))
            strocnyN(i, j) = vrel*((vocnN(i, j) - v(i, j + 1))*cosw + (uocnN(i, j) - uvelN(i, j))*sinw*sg)
            if (.not. ieee_is_finite(v(i, j + 1))) bad = .true.
            if (abs(v(i, j + 1)) > vel_cap) bad = .true.
         end if
         if (bad) ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_DYN_VEL)
      end do
   end subroutine evp_finish_impl

   ! ------------------------------------------------------------------
   ! local periodic wrap (single rank)
   ! ------------------------------------------------------------------

   subroutine evp_wrap_u_impl(nx, ny, is, ie, js, je, ng, wrap_x, wrap_y, u)
      !! Periodic wrap of the x-face array: x seam (`u(is) = u(ie+1)`) then
      !! the ghost columns, then y.  Single rank only; the multi-rank and
      !! fold cases use `ice_halo_faces`.
      integer, intent(in) :: nx, ny, is, ie, js, je, ng
      logical, intent(in) :: wrap_x, wrap_y
      real(wp), intent(inout) :: u(nx + 1, ny)
      integer :: i, j, kk

      if (wrap_x) then
         do concurrent(j=1:ny, kk=0:ng)
            u(is - kk, j) = u(ie + 1 - kk, j)
         end do
         do concurrent(j=1:ny, kk=1:ng)
            u(ie + 1 + kk, j) = u(is + kk, j)
         end do
      end if
      if (wrap_y) then
         do concurrent(kk=1:ng, i=1:nx + 1)
            u(i, js - kk) = u(i, je + 1 - kk)
            u(i, je + kk) = u(i, js - 1 + kk)
         end do
      end if
   end subroutine evp_wrap_u_impl

   subroutine evp_wrap_v_impl(nx, ny, is, ie, js, je, ng, wrap_x, wrap_y, v)
      !! Periodic wrap of the y-face array (see `evp_wrap_u_impl`).
      integer, intent(in) :: nx, ny, is, ie, js, je, ng
      logical, intent(in) :: wrap_x, wrap_y
      real(wp), intent(inout) :: v(nx, ny + 1)
      integer :: i, kk

      if (wrap_x) then
         do concurrent(kk=1:ng, i=1:ny + 1)
            v(is - kk, i) = v(ie + 1 - kk, i)
            v(ie + kk, i) = v(is - 1 + kk, i)
         end do
      end if
      if (wrap_y) then
         do concurrent(i=1:nx, kk=0:ng)
            v(i, js - kk) = v(i, je + 1 - kk)
         end do
         do concurrent(i=1:nx, kk=1:ng)
            v(i, je + 1 + kk) = v(i, js + kk)
         end do
      end if
   end subroutine evp_wrap_v_impl

   ! ------------------------------------------------------------------
   ! host shims
   ! ------------------------------------------------------------------

   subroutine ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, wrap_x, wrap_y, u, v)
      !! Single-rank periodic wrap of the face pair (`wrap_x`, `wrap_y`
      !! false: no-op).  The ghost width is `nghost = is - 1`.
      integer, intent(in) :: nx, ny, is, ie, js, je
      logical, intent(in) :: wrap_x, wrap_y
      real(wp), intent(inout) :: u(nx + 1, ny), v(nx, ny + 1)
      if (.not. (wrap_x .or. wrap_y)) return
      call evp_wrap_u_impl(nx, ny, is, ie, js, je, is - 1, wrap_x, wrap_y, u)
      call evp_wrap_v_impl(nx, ny, is, ie, js, je, is - 1, wrap_x, wrap_y, v)
   end subroutine ice_evp_c_wrap_faces

   subroutine evp_copy_faces_impl(nx, ny, us, vs, ud, vd)
      !! `(ud, vd) = (us, vs)` on the whole face pair (one launch).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: us(nx + 1, ny), vs(nx, ny + 1)
      real(wp), intent(inout) :: ud(nx + 1, ny), vd(nx, ny + 1)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx + 1)
         ud(i, j) = us(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         vd(i, j) = vs(i, j)
      end do
   end subroutine evp_copy_faces_impl

   subroutine ice_evp_c_copy_faces(nx, ny, us, vs, ud, vd)
      !! Copy the face pair `(us, vs)` to `(ud, vd)` on the device (the
      !! per-step set-up of the ping-pong buffers, and the copy back).
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: us(nx + 1, ny), vs(nx, ny + 1)
      real(wp), intent(inout) :: ud(nx + 1, ny), vd(nx, ny + 1)
      call evp_copy_faces_impl(nx, ny, us, vs, ud, vd)
   end subroutine ice_evp_c_copy_faces

   subroutine ice_evp_c_average(w, nx, ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm)
      !! K5 on the work arrays: fills `uvelN`, `vvelE`, `uvelU`, `vvelU`
      !! from the face pair.  The subcycles no longer need them (K1 and K4
      !! interpolate inline); `ice_evp_c_step` calls it once after the last
      !! subcycle for `ice_evp_c_finish`, and diagnostics may call it any time
      !! after the face halo of `(u, v)`.
      type(ice_evp_c_work_t), intent(inout) :: w
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), epm(nx, ny), npm(nx, ny), uvm(nx, ny)
      call evp_average_impl(nx, ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm, &
                            w%uvelN, w%vvelE, w%uvelU, w%vvelU)
   end subroutine ice_evp_c_average

   subroutine ice_evp_c_substep(w, k, nx, ny, is, ie, js, je, uo, vo, un, vn, &
                                dxE, dyE, dxN, dyN, dxU, dyU, dxT, dyT, uarea, tarea, earea, narea, earear, narear, &
                                ratiodxN, ratiodxNr, ratiodyE, ratiodyEr, epm, npm, uvm, hm, &
                                uocnE, vocnE, uocnN, vocnN)
      !! K1-K4 of one subcycle (4 launches, no halo): the piece a multi-rank
      !! caller runs between two face exchanges.  Ping-pong: reads `(uo, vo)`,
      !! writes the owned ice faces of `(un, vn)`; at every other face `(un, vn)`
      !! must already equal `(uo, vo)` (copy the pair once per step with
      !! `ice_evp_c_copy_faces`, after the first halo).  The caller then
      !! refreshes the ghost ring of `(un, vn)` (`ice_halo_faces`, or
      !! `ice_evp_c_wrap_faces` on one periodic rank) and swaps the roles of
      !! the two pairs for the next subcycle.
      type(ice_evp_c_work_t), intent(inout) :: w
      type(ice_evp_c_const_t), intent(in) :: k
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: uo(nx + 1, ny), vo(nx, ny + 1)
      real(wp), intent(inout) :: un(nx + 1, ny), vn(nx, ny + 1)
      real(wp), intent(in) :: dxE(nx, ny), dyE(nx, ny), dxN(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), uarea(nx, ny), tarea(nx, ny), earea(nx, ny), narea(nx, ny)
      real(wp), intent(in) :: earear(nx, ny), narear(nx, ny)
      real(wp), intent(in) :: ratiodxN(nx, ny), ratiodxNr(nx, ny), ratiodyE(nx, ny), ratiodyEr(nx, ny)
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny), uvm(nx, ny), hm(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)

      call evp_strain_u_impl(nx, ny, is, ie, js, je, k%e_factor, k%avg_strength, uo, vo, earea, narea, &
                             dxE, dyN, dxU, dyU, ratiodxN, ratiodxNr, ratiodyE, ratiodyEr, &
                             epm, npm, uvm, w%imask, w%shearU, w%deltaU)
      call evp_stress_t_impl(nx, ny, is, ie, js, je, k%arlx1i, k%revp, k%denom1, k%epp2i, k%e_factor, k%Ktens, &
                             k%deltamin, k%cap_max, uo, vo, dxN, dyE, dxT, dyT, uarea, tarea, w%strength, &
                             w%shearU, w%imask, w%zetax2T, w%etax2T, w%stresspT, w%stressmT)
      call evp_stress_u_impl(nx, ny, is, ie, js, je, k%arlx1i, k%revp, k%denom1, k%epp2i, k%Ktens, k%deltamin, &
                             k%cap_max, k%avg_strength, hm, tarea, uarea, w%etax2T, w%strengthU, w%deltaU, &
                             w%shearU, w%imask, w%stress12U)
      call evp_momentum_impl(nx, ny, is, ie, js, je, k%revp, k%brlx, k%rhow, k%cw, k%cosw, k%sinw, &
                             k%free_drift, k%rheo_amin, uo, vo, un, vn, earea, narea, epm, npm, &
                             dxE, dyE, dxN, dyN, dxU, dyU, dxT, dyT, earear, narear, w%stresspT, w%stressmT, &
                             w%stress12U, w%imask, w%aiE, w%aiN, w%emassdti, w%nmassdti, w%fmE, w%fmN, &
                             w%forcexE, w%forceyN, w%uinit, w%vinit, uocnE, vocnE, uocnN, vocnN)
   end subroutine ice_evp_c_substep

   subroutine ice_evp_c_finish(w, k, nx, ny, is, ie, js, je, u, v, dxN, dyE, dxT, dyT, uarea, tarear, &
                               uocnE, vocnE, uocnN, vocnN, ice_err)
      !! The after-loop kernel on the work arrays: `divu`, `shear`, `vort`,
      !! `rdg_conv`, `rdg_shear`, `strocnxE`, `strocnyN` and the velocity
      !! guard.
      type(ice_evp_c_work_t), intent(inout) :: w
      type(ice_evp_c_const_t), intent(in) :: k
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: u(nx + 1, ny), v(nx, ny + 1)
      real(wp), intent(in) :: dxN(nx, ny), dyE(nx, ny), dxT(nx, ny), dyT(nx, ny), uarea(nx, ny), tarear(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)

      call evp_finish_impl(nx, ny, is, ie, js, je, k%e_factor, k%epp2i, k%rhow, k%cw, k%cosw, k%sinw, k%vel_cap, &
                           u, v, w%uvelN, w%vvelE, w%shearU, dxN, dyE, dxT, dyT, uarea, tarear, w%imask, &
                           w%aiE, w%aiN, w%fmE, w%fmN, uocnE, vocnE, uocnN, vocnN, &
                           w%divu, w%shear, w%vort, w%rdg_conv, w%rdg_shear, w%strocnxE, w%strocnyN, ice_err)
   end subroutine ice_evp_c_finish

   subroutine ice_evp_c_div_stress(w, k, nx, ny, is, ie, js, je, dxE, dyE, dxN, dyN, dxU, dyU, dxT, dyT, &
                                   earear, narear)
      !! Diagnostic: the stress divergence at E and N into
      !! `w%strintxE`/`w%strintyN` for the current stresses and masks.
      type(ice_evp_c_work_t), intent(inout) :: w
      type(ice_evp_c_const_t), intent(in) :: k
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dxE(nx, ny), dyE(nx, ny), dxN(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), earear(nx, ny), narear(nx, ny)
      call evp_div_stress_impl(nx, ny, is, ie, js, je, k%free_drift, k%rheo_amin, dxE, dyE, dxN, dyN, dxU, dyU, &
                               dxT, dyT, earear, narear, w%stresspT, w%stressmT, w%stress12U, w%iceEmask, &
                               w%iceNmask, w%aiE, w%aiN, w%strintxE, w%strintyN)
   end subroutine ice_evp_c_div_stress

   subroutine ice_evp_c_step(w, k, nx, ny, ncat, is, ie, js, je, rp, c, wrap_x, wrap_y, &
                             dxE, dyE, dxN, dyN, dxU, dyU, dxT, dyT, uarea, tarea, tarear, earea, narea, &
                             earear, narear, ratiodxN, ratiodxNr, ratiodyE, ratiodyEr, epm, npm, uvm, hm, &
                             aice, vice, vsno, aice_init, aice0, aicen, vicen, &
                             fcorE, fcorN, strairxE, strairyN, ss_tltxE, ss_tltyN, &
                             uocnE, vocnE, uocnN, vocnN, u, v, ice_err)
      !! One dynamics step: preparation, `ndte` EVP subcycles, deformation
      !! rates and ice-ocean stress.  Single rank: the periodic wrap of the
      !! face pair is local (`wrap_x`, `wrap_y`); the fold and the multi-rank
      !! exchange are the caller's (`ice_evp_c_prep`, `ice_evp_c_substep`,
      !! `ice_evp_c_average`, `ice_evp_c_finish` with `ice_halo_faces` in
      !! between; see the module header).  The velocity, the five stresses
      !! and the ice masks in `w` are the state carried between steps.
      !! Arguments: see `ice_evp_c_prep` for the T-cell inputs and the forcing
      !! at the faces; the metric arrays are the ice names of
      !! `ice_cgrid_metrics_t` plus `dxT, dyT, tarea, tarear, hm`.
      !! `uocnE/vocnE` and `uocnN/vocnN` are the two ocean-velocity
      !! components at E and at N.
      type(ice_evp_c_work_t), intent(inout) :: w
      type(ice_evp_c_const_t), intent(in) :: k
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      type(ice_ridge_params_t), intent(in) :: rp
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: wrap_x, wrap_y
      real(wp), intent(in) :: dxE(nx, ny), dyE(nx, ny), dxN(nx, ny), dyN(nx, ny), dxU(nx, ny), dyU(nx, ny)
      real(wp), intent(in) :: dxT(nx, ny), dyT(nx, ny), uarea(nx, ny), tarea(nx, ny), tarear(nx, ny)
      real(wp), intent(in) :: earea(nx, ny), narea(nx, ny), earear(nx, ny), narear(nx, ny)
      real(wp), intent(in) :: ratiodxN(nx, ny), ratiodxNr(nx, ny), ratiodyE(nx, ny), ratiodyEr(nx, ny)
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny), uvm(nx, ny), hm(nx, ny)
      real(wp), intent(in) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice_init(nx, ny), aice0(nx, ny)
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat)
      real(wp), intent(in) :: fcorE(nx, ny), fcorN(nx, ny)
      real(wp), intent(in) :: strairxE(nx, ny), strairyN(nx, ny), ss_tltxE(nx, ny), ss_tltyN(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)
      real(wp), intent(inout) :: u(nx + 1, ny), v(nx, ny + 1)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: ksub

      call ice_evp_c_prep(w, k, nx, ny, ncat, is, ie, js, je, rp, c, hm, aice, vice, vsno, aice_init, aice0, &
                          aicen, vicen, tarea, epm, npm, fcorE, fcorN, strairxE, strairyN, ss_tltxE, &
                          ss_tltyN, uocnE, vocnN, u, v, ice_err)
      call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, wrap_x, wrap_y, u, v)
      call ice_evp_c_copy_faces(nx, ny, u, v, w%u2, w%v2)
      do ksub = 1, k%ndte
         if (mod(ksub, 2) == 1) then          ! old pair (u, v) -> new pair (u2, v2)
            call ice_evp_c_substep(w, k, nx, ny, is, ie, js, je, u, v, w%u2, w%v2, dxE, dyE, dxN, dyN, dxU, dyU, &
                                   dxT, dyT, uarea, tarea, earea, narea, earear, narear, ratiodxN, ratiodxNr, &
                                   ratiodyE, ratiodyEr, epm, npm, uvm, hm, uocnE, vocnE, uocnN, vocnN)
            call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, wrap_x, wrap_y, w%u2, w%v2)
         else                                 ! old pair (u2, v2) -> new pair (u, v)
            call ice_evp_c_substep(w, k, nx, ny, is, ie, js, je, w%u2, w%v2, u, v, dxE, dyE, dxN, dyN, dxU, dyU, &
                                   dxT, dyT, uarea, tarea, earea, narea, earear, narear, ratiodxN, ratiodxNr, &
                                   ratiodyE, ratiodyEr, epm, npm, uvm, hm, uocnE, vocnE, uocnN, vocnN)
            call ice_evp_c_wrap_faces(nx, ny, is, ie, js, je, wrap_x, wrap_y, u, v)
         end if
      end do
      if (mod(k%ndte, 2) == 1) call ice_evp_c_copy_faces(nx, ny, w%u2, w%v2, u, v)
      call ice_evp_c_average(w, nx, ny, is, ie, js, je, u, v, earea, narea, epm, npm, uvm)
      call ice_evp_c_finish(w, k, nx, ny, is, ie, js, je, u, v, dxN, dyE, dxT, dyT, uarea, tarear, &
                            uocnE, vocnE, uocnN, vocnN, ice_err)
   end subroutine ice_evp_c_step

end module rdb_ice_evp_c
