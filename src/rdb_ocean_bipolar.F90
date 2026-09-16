!! Analytic TRIPOLAR (Murray 1996) bipolar-cap coordinate map.
module rdb_ocean_bipolar
   !! Pure corner-coordinate generator for the bipolar Arctic cap of a
   !! TRIPOLAR ocean grid (Murray 1996, J. Comput. Phys. 126, 251-273).
   !!
   !! Ordinary lon-lat for `lat <= phi_join`; a conformal BIPOLAR cap for
   !! `lat > phi_join`. Two grid poles sit ON the join latitude at antipodal
   !! longitudes `lon_pole` and `lon_pole+180`, so the join circle is a
   !! coordinate line and the grid is C0-continuous across it.
   !!
   !! Map: stereographic projection from the south pole
   !! (Z = tan((90-lat)/2)*exp(i*lon)), de-rotate by `lon_pole`, Mobius
   !! w = (Zp-A)/(Zp+A) with A = tan((90-phi_join)/2). Conformal coords
   !! (mu, xi) = (log|w|, arg w) are orthogonal by construction. The logical
   !! parameterisation makes the join a coordinate line and continuous with
   !! the lon-lat grid: i sets mu = log|tan((lam-lon_pole)/2)| (lon = lam on
   !! the join ring); j (s in [0,1]) rotates xi = sgn*(pi/2 + s*pi/2) from
   !! the join (s=0) to the self-conjugate fold line (s=1). Inverse:
   !! w = exp(mu+i*xi), Zp = A*(1+w)/(1-w), Z = Zp*exp(i*lon_pole).
   use rdb_constants, only: wp, PI, DEG2RAD, RAD2DEG
   implicit none
   private

   public :: bipolar_corner_latlon
   public :: bipolar_pole_lat

contains

   pure subroutine bipolar_corner_latlon(lam_deg, s, phi_join, lon_pole, lat, lon)
      !! Map a logical cap location to geographic (lat, lon) in degrees.
      !!   `lam_deg` : pseudo-longitude (geographic lon at the join ring),
      !!               any real (wrapped internally); the i-direction.
      !!   `s`       : cap-row fraction, 0 = join ring (`lat = phi_join`),
      !!               1 = fold line; the j-direction.
      !!   `phi_join`: join latitude (deg); cap covers lat > phi_join.
      !!   `lon_pole`: longitude (deg) of the first cap pole; partner at
      !!               `lon_pole + 180`.
      real(wp), intent(in) :: lam_deg, s, phi_join, lon_pole
      real(wp), intent(out) :: lat, lon

      real(wp) :: a_focus, half, t, mu, xi, sgn
      real(wp) :: wr, wi, denr, deni, dnorm, zpr, zpi
      real(wp) :: lpr, zr, zi, rmag

      a_focus = tan(DEG2RAD*(90.0_wp - phi_join)*0.5_wp)

      ! mu from the join-continuity relation: lon_join = lon_pole + 2 atan(exp(mu)).
      half = DEG2RAD*(lam_deg - lon_pole)*0.5_wp
      t = tan(half)
      if (t > 0.0_wp) then
         sgn = 1.0_wp
         mu = log(t)
      else if (t < 0.0_wp) then
         sgn = -1.0_wp
         mu = log(-t)
      else
         ! lam == lon_pole (mod 360): a cap pole.  mu -> -inf; clamp finite
         ! so the corner sits arbitrarily close to the pole without an Inf.
         sgn = 1.0_wp
         mu = -50.0_wp
      end if

      ! xi rotates from the join (sgn*pi/2) to the fold (sgn*pi) as s: 0 -> 1.
      xi = sgn*(PI*0.5_wp + s*PI*0.5_wp)

      ! w = exp(mu) * exp(i*xi)
      wr = exp(mu)*cos(xi)
      wi = exp(mu)*sin(xi)

      ! Zp = A * (1 + w)/(1 - w)
      denr = 1.0_wp - wr
      deni = -wi
      dnorm = denr*denr + deni*deni
      if (dnorm <= 0.0_wp) then
         ! w == 1 (mu -> +inf): the other pole.  Place at phi_join on the
         ! partner meridian.
         lat = phi_join
         lon = lon_pole + 180.0_wp
         if (lon >= 180.0_wp) lon = lon - 360.0_wp
         return
      end if
      ! (1+w)/(1-w)
      zpr = ((1.0_wp + wr)*denr + wi*deni)/dnorm
      zpi = (wi*denr - (1.0_wp + wr)*deni)/dnorm
      zpr = a_focus*zpr
      zpi = a_focus*zpi

      ! Z = Zp * exp(i*lon_pole)
      lpr = DEG2RAD*lon_pole
      zr = zpr*cos(lpr) - zpi*sin(lpr)
      zi = zpr*sin(lpr) + zpi*cos(lpr)

      ! stereographic inverse
      rmag = sqrt(zr*zr + zi*zi)
      lat = 90.0_wp - RAD2DEG*2.0_wp*atan(rmag)
      if (zr == 0.0_wp .and. zi == 0.0_wp) then
         lon = lon_pole + 180.0_wp
      else
         lon = RAD2DEG*atan2(zi, zr)
      end if
      ! wrap lon to [-180, 180)
      lon = modulo(lon + 180.0_wp, 360.0_wp) - 180.0_wp
   end subroutine bipolar_corner_latlon

   pure function bipolar_pole_lat(phi_join) result(lat)
      !! Latitude (deg) of the two cap poles — they sit ON the join
      !! latitude, so this is simply `phi_join`.
      real(wp), intent(in) :: phi_join
      real(wp) :: lat
      lat = phi_join
   end function bipolar_pole_lat

end module rdb_ocean_bipolar
