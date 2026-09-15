!! Analytical + device tests for the prognostic mesoscale eddy kinetic
!! energy slot (capability [5], `rdb_ocean_meke`).  All cases RUN THE DEVICE
!! KERNELS via `meke_step` (the production entry).  A synthetic GM slot
!! supplies `gm%gm_src` (the PE-release source) directly so the tests do not
!! depend on the full slope/EOS chain.
!!
!!  1. meke_equilibrium       — constant gm_src + linear damping, no transport
!!                              ⇒ E relaxes to the analytic src/damp.
!!  2. meke_drag_stable       — huge sdt*damp_rate ⇒ E stays positive + decays
!!                              (backward-Euler implicit drag).
!!  3. meke_diffusion_conserves — closed domain, Sum E*area conserved under
!!                              the harmonic-mass diffusion stage.
!!  4. meke_kh_closure        — given E, geometry ⇒ kh = khcoeff*sqrt(2*gt2*E)*L
!!                              for a single (grid) length scale (hand-computed).
!!  5. meke_feeds_khth        — khth_fac>0 ⇒ varmix%khth picks up the geom-mean;
!!                              khth_fac=0 ⇒ varmix%khth unchanged (seam identity).
module test_ocean_meke
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_gm, only: ocean_gm_t
   use rdb_ocean_varmix, only: ocean_varmix_t
   use rdb_ocean_meke, only: ocean_meke_t, meke_step
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_meke_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp
   real(wp), parameter :: DZ = 50.0_wp

contains

   subroutine collect_ocean_meke_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("meke_equilibrium", test_equilibrium), &
                  new_unittest("meke_drag_stable", test_drag_stable), &
                  new_unittest("meke_diffusion_conserves", test_diffusion_conserves), &
                  new_unittest("meke_kh_closure", test_kh_closure), &
                  new_unittest("meke_feeds_khth", test_feeds_khth), &
                  new_unittest("meke_biharmonic_conserves", test_biharmonic_conserves), &
                  new_unittest("meke_rhines_length", test_rhines_length), &
                  new_unittest("meke_advection_conserves", test_advection_conserves), &
                  new_unittest("meke_advection_inert", test_advection_inert), &
                  new_unittest("meke_bbl_drag", test_bbl_drag), &
                  new_unittest("meke_frictional_source", test_frictional_source), &
                  new_unittest("meke_drag_equilibrium_cdrag", test_drag_equilibrium_cdrag) &
                  ]
   end subroutine collect_ocean_meke_tests

   ! ------------------------------------------------------------------
   ! Frictional-source seam: a NEGATIVE hvisc KE-dissipation rate (ke_diss)
   ! with frcoeff>=0 injects eddy energy (`src -= frcoeff·i_mass·ke_diss`,
   ! and ke_diss<0 ⇒ positive source).  frcoeff<0 ignores ke_diss entirely;
   ! ke_diss=0 with frcoeff on is a no-op (bit-identical seam).
   ! ------------------------------------------------------------------
   subroutine run_frict(frcoeff, mom_val, e_out)
      real(wp), intent(in) :: frcoeff, mom_val
      real(wp), intent(out) :: e_out
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp, DZL = 50.0_wp
      real(wp), allocatable :: mom(:, :)
      call make_grid(grid, NX, NY, DX)
      call setup_ms(ms, grid, NZ, DZL)
      call make_cartesian_metrics(metrics, grid)
      call setup_gm(gm, grid, NZ, 0.0_wp)
      call meke%init(grid, nz_ml=NZ)
      meke%enable = .true.
      meke%gmcoeff = -1.0_wp   ! GM source off
      meke%frcoeff = frcoeff
      meke%bgsrc = 0.0_wp
      meke%damping = 0.0_wp
      meke%cdrag = 0.0_wp      ! no drag ⇒ E changes only via the frictional source
      meke%khcoeff = -1.0_wp
      meke%kh = -1.0_wp
      meke%k4 = -1.0_wp
      meke%meke = 1.0_wp
      allocate (mom(grid%nx_total, grid%ny_total), source=mom_val)
      call map_in(ms, gm, meke)
      !$acc enter data copyin(mom)
      call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT, ke_diss_ext=mom)
      !$acc exit data delete(mom)
      call map_out(ms, gm, meke)
      e_out = meke%meke(NGHOST + 1, NGHOST + 1)
      deallocate (mom)
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_frict

   subroutine test_frictional_source(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: e_off, e_on, e_on0
      checks: block
         call run_frict(-1.0_wp, -100.0_wp, e_off)   ! frcoeff off ⇒ mom ignored
         call run_frict(0.5_wp, -100.0_wp, e_on)     ! frcoeff on, mom<0 ⇒ source
         call run_frict(0.5_wp, 0.0_wp, e_on0)       ! frcoeff on but mom=0
         call check(error, e_on > e_off + 1.0e-6_wp, &
                    "negative ke_diss + frcoeff>=0 must inject eddy energy (E_on > E_off)")
         if (allocated(error)) exit checks
         call check(error, abs(e_off - 1.0_wp) < 1.0e-12_wp, &
                    "frcoeff<0 ⇒ no source, no drag ⇒ E unchanged from IC")
         if (allocated(error)) exit checks
         call check(error, abs(e_on0 - e_off) < 1.0e-12_wp, &
                    "ke_diss=0 ⇒ frictional source inert (bit-identical)")
      end block checks
   end subroutine test_frictional_source

   ! ------------------------------------------------------------------
   ! BBL-drag seam: the resolved bed-layer eddy velocity (|u_bed|²) adds to
   ! the MEKE bottom-drag rate when `use_bbl_drag` is on, so a non-zero bed
   ! velocity dissipates MORE eddy energy; a zero bed velocity reproduces
   ! the prior (off) drag exactly (bit-identical seam).
   ! ------------------------------------------------------------------
   subroutine run_bbl(use_bbl, ubed, e_out, bf2_out)
      logical, intent(in) :: use_bbl
      real(wp), intent(in) :: ubed
      real(wp), intent(out) :: e_out
      real(wp), intent(out), optional :: bf2_out
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp, DZL = 10.0_wp
      call make_grid(grid, NX, NY, DX)
      call setup_ms(ms, grid, NZ, DZL)
      ms%u_face_x_layer = ubed        ! uniform bed-layer (k=1) velocity
      ms%v_face_y_layer = 0.0_wp
      call make_cartesian_metrics(metrics, grid)
      call setup_gm(gm, grid, NZ, 0.0_wp)
      call meke%init(grid, nz_ml=NZ)
      meke%enable = .true.
      meke%gmcoeff = -1.0_wp   ! GM source off
      meke%bgsrc = 0.0_wp
      meke%damping = 0.0_wp    ! isolate the bottom-drag term
      meke%khcoeff = -1.0_wp
      meke%kh = -1.0_wp
      meke%k4 = -1.0_wp
      meke%cdrag = 0.1_wp      ! large ⇒ measurable drag
      meke%cd_scale = 0.0_wp
      meke%use_bbl_drag = use_bbl
      meke%meke = 1.0_wp
      call map_in(ms, gm, meke)
      call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
      call map_out(ms, gm, meke)
      e_out = meke%meke(NGHOST + 1, NGHOST + 1)
      if (present(bf2_out)) bf2_out = meke%bottom_fac2(NGHOST + 1, NGHOST + 1)
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_bbl

   subroutine test_bbl_drag(error)
      !! Test 1 (PR-5): keep the three seam-ordering assertions, then pin the
      !! MAGNITUDE of the post-drag E to the closed-form backward-Euler
      !! update -- the only assertion that can see a missing/spurious rho0
      !! factor in `drag_rate` (a positive multiplicative constant preserves
      !! positivity, ordering, and the u_bbl2=0 bit-identity, but NOT the
      !! magnitude).  Mirrors run_bbl's local constants exactly:
      !!   NZ=1, DZL=10, DT=1800, cdrag=0.1, cd_scale=0, damping=0,
      !!   uscale=0 (default), dtscale=1 (default), RHO0=1025, E0=1.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: e_off, e_on, e_on0, bf2
      real(wp) :: i_mass, u_bbl2, drag_rate, damp_rate, e_expect
      real(wp), parameter :: DZL = 10.0_wp, DT = 1800.0_wp, CDRAG = 0.1_wp
      real(wp), parameter :: UBED = 2.0_wp, E0 = 1.0_wp
      real(wp), parameter :: USCALE = 0.0_wp, DTSCALE = 1.0_wp, DAMPING = 0.0_wp
      integer, parameter :: NZ = 1
      checks: block
         call run_bbl(.false., UBED, e_off)          ! BBL off (bed vel ignored)
         call run_bbl(.true., UBED, e_on, bf2)        ! BBL on, |u_bed|=2
         call run_bbl(.true., 0.0_wp, e_on0)          ! BBL on but bed vel = 0
         call check(error, e_on > 0.0_wp .and. e_off > 0.0_wp, "drag keeps MEKE positive")
         if (allocated(error)) exit checks
         call check(error, e_on < e_off - 1.0e-6_wp, &
                    "BBL eddy velocity must add bottom drag (E_on < E_off)")
         if (allocated(error)) exit checks
         call check(error, abs(e_on0 - e_off) < 1.0e-12_wp, &
                    "zero bed velocity ⇒ bit-identical to BBL-off")
         if (allocated(error)) exit checks

         ! bottomFac2 must be exactly 1 (cd_scale=0, no wavespeed ⇒
         ! ldeform=0) for the closed form below to be valid; a future
         ! default change (e.g. PR-1's wavespeed wire) must fail HERE,
         ! loudly, rather than corrupt the magnitude assertion silently.
         call check(error, abs(bf2 - 1.0_wp) < 1.0e-12_wp, &
                    "bottom_fac2 must be exactly 1 for the closed-form shortcut to hold")
         if (allocated(error)) exit checks

         ! Closed-form backward-Euler post-drag E (PR-5: drag_rate carries
         ! the rho0 factor -- without it this assertion fails by ~43x).
         i_mass = 1.0_wp/(RHO0*real(NZ, wp)*DZL)
         u_bbl2 = UBED*UBED
         drag_rate = (RHO0*i_mass)*sqrt(CDRAG*CDRAG*(2.0_wp*1.0_wp*E0 + u_bbl2 + USCALE*USCALE))
         damp_rate = DAMPING + drag_rate*1.0_wp
         e_expect = E0/(1.0_wp + DT*DTSCALE*damp_rate)
         call check(error, abs(e_on - e_expect)/abs(e_expect) < 1.0e-12_wp, &
                    "BBL-on E must match the closed-form drag_rate=rho0*i_mass*sqrt(...) update")
      end block checks
   end subroutine test_bbl_drag

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine setup_ms(ms, grid, nz, dz)
      !! Uniform-thickness, single-layer-ish multilayer state (no rho_layer ⇒
      !! mass = rho0*Sum h via the gm%rho0 fallback).
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = dz
      ms%rho_layer = RHO0   ! uniform density ⇒ mass = RHO0*Sum h
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_ms

   subroutine setup_gm(gm, grid, nz, gm_src_val)
      !! A synthetic GM slot: init + enable + a uniform gm_src.
      type(ocean_gm_t), intent(inout) :: gm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: gm_src_val
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%rho0 = RHO0
      gm%gm_src = gm_src_val
   end subroutine setup_gm

   subroutine map_in(ms, gm, meke)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_gm_t), intent(inout) :: gm
      type(ocean_meke_t), intent(inout) :: meke
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(gm)
      call gm%enter_data()
      !$acc enter data copyin(meke)
      call meke%enter_data()
   end subroutine map_in

   subroutine map_out(ms, gm, meke)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_gm_t), intent(inout) :: gm
      type(ocean_meke_t), intent(inout) :: meke
      !$acc update self(meke%meke, meke%kh_diff, meke%le, meke%bottom_fac2)
      call meke%exit_data()
      !$acc exit data delete(meke)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! Test 1: E -> src/damp under constant source + linear damping
   ! ------------------------------------------------------------------
   subroutine test_equilibrium(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp
      real(wp), parameter :: DAMP = 1.0e-6_wp, DT = 1800.0_wp
      real(wp), parameter :: GMCOEFF = 0.2_wp, GMSRC = 1.0e-3_wp
      real(wp) :: i_mass, src_per_mass, e_eq, e_mid, relerr
      integer :: it
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, GMSRC)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%damping = DAMP
         meke%gmcoeff = GMCOEFF
         meke%khcoeff = -1.0_wp  ! closure off (focus on the budget)
         meke%kh = -1.0_wp       ! no diffusion ⇒ damp_step=1 (single drag)
         meke%k4 = -1.0_wp
         meke%cd_scale = 0.0_wp  ! bottomFac2 floor = min_gamma2 (tiny)
         meke%uscale = 0.0_wp
         meke%cdrag = 0.0_wp     ! drag_rate=0 ⇒ damp_rate = damping only

         call map_in(ms, gm, meke)
         do it = 1, 200000
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         call map_out(ms, gm, meke)

         ! Analytic equilibrium: E_eq = src/damp, src = gmcoeff*I_mass*gm_src.
         i_mass = 1.0_wp/(RHO0*real(NZ, wp)*DZ)
         src_per_mass = GMCOEFF*i_mass*GMSRC
         e_eq = src_per_mass/DAMP
         e_mid = meke%meke(NGHOST + 1, NGHOST + 1)
         relerr = abs(e_mid - e_eq)/abs(e_eq)
         call check(error, relerr < 1.0e-4_wp, &
                    "MEKE must relax to src/damp under constant source + linear damping")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_equilibrium

   ! ------------------------------------------------------------------
   ! Test (PR-5): analytic E^{3/2} equilibrium with cdrag /= 0.
   !
   ! With cdrag=0 (test_equilibrium above) the sink is LINEAR in E and the
   ! fixed point is the textbook src/damp.  With cdrag/=0 the bottom-drag
   ! sink is `drag_rate*bf2*E ~ E^{3/2}` (drag_rate itself ~ sqrt(E)), so
   ! the fixed point of `dE/dt = bgsrc - drag_rate*bf2*E` is NOT src/damp:
   !
   !   bf2=1 (cd_scale=0, no wavespeed/varmix ⇒ ldeform=0)
   !   rho0*i_mass = 1/depth_tot exactly (uniform rho_layer=RHO0=gm%rho0)
   !   S = bgsrc = drag_rate*E = (1/H)*cdrag*sqrt(2*E)*E = (cdrag/H)*sqrt(2)*E^1.5
   !   ⇒ E_eq = ( S*H / (cdrag*sqrt(2)) )^(2/3)
   !
   ! This is PR-5's bug signature to the power 2/3: a missing/spurious rho0
   ! factor moves E_eq by rho0^(2/3) (~102x for RHO0=1025), not linearly.
   ! ------------------------------------------------------------------
   subroutine test_drag_equilibrium_cdrag(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 4.0_wp
      real(wp), parameter :: SRC = 1.0e-6_wp, CDRAG = 5.0e-3_wp
      integer, parameter :: N_ITER = 100000
      real(wp) :: h_tot, e_eq, e_final, e_prev100, relerr, drift100
      integer :: it
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)   ! gm_src unused (gmcoeff<0)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp    ! GM source off
         meke%frcoeff = -1.0_wp    ! frictional source off
         meke%bgsrc = SRC          ! constant background source (S)
         meke%damping = 0.0_wp     ! linear damping off ⇒ drag is the ONLY sink
         meke%uscale = 0.0_wp
         meke%cdrag = CDRAG
         meke%cd_scale = 0.0_wp    ! ⇒ bottomFac2 = 1 (no wavespeed ⇒ ldeform=0)
         meke%use_bbl_drag = .false.  ! u_bbl2 stays 0
         meke%khcoeff = -1.0_wp    ! closure off
         meke%kh = -1.0_wp         ! no diffusion ⇒ damp_step=1 (single drag half)
         meke%k4 = -1.0_wp
         meke%advection_factor = 0.0_wp

         call map_in(ms, gm, meke)
         do it = 1, N_ITER - 100
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         !$acc update self(meke%meke)
         e_prev100 = meke%meke(NGHOST + 1, NGHOST + 1)
         do it = N_ITER - 99, N_ITER
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         call map_out(ms, gm, meke)
         e_final = meke%meke(NGHOST + 1, NGHOST + 1)

         ! Convergence guard: the last-100-step change must be tiny, so a
         ! non-converged run fails as "did not reach equilibrium" rather
         ! than accidentally passing near the initial condition.
         drift100 = abs(e_final - e_prev100)/max(abs(e_final), 1.0e-30_wp)
         call check(error, drift100 < 1.0e-6_wp, &
                    "MEKE must have converged to a fixed point (last-100-step drift too large)")
         if (allocated(error)) exit checks

         h_tot = real(NZ, wp)*DZ
         e_eq = (SRC*h_tot/(CDRAG*sqrt(2.0_wp)))**(2.0_wp/3.0_wp)
         relerr = abs(e_final - e_eq)/abs(e_eq)
         call check(error, relerr < 1.0e-4_wp, &
                    "MEKE must relax to the E^1.5 equilibrium E_eq=(S*H/(cdrag*sqrt2))^(2/3) &
                    &under a constant source + cdrag/=0 bottom drag")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_drag_equilibrium_cdrag

   ! ------------------------------------------------------------------
   ! Test 2: implicit drag stays positive + decays for huge stiffness
   ! ------------------------------------------------------------------
   subroutine test_drag_stable(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      real(wp) :: e0, e1
      integer :: i, j
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)   ! no source
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp   ! GM source off
         meke%bgsrc = 0.0_wp
         meke%damping = 1.0e3_wp  ! HUGE linear damping ⇒ stiff sink
         meke%khcoeff = -1.0_wp
         meke%kh = -1.0_wp
         meke%k4 = -1.0_wp
         meke%cdrag = 0.0_wp

         ! Seed a positive eddy energy everywhere.
         e0 = 5.0_wp
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               meke%meke(i, j) = e0
            end do
         end do

         call map_in(ms, gm, meke)
         call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         call map_out(ms, gm, meke)

         e1 = meke%meke(NGHOST + 1, NGHOST + 1)
         call check(error, e1 > 0.0_wp, "implicit drag must keep MEKE positive")
         if (allocated(error)) exit checks
         call check(error, e1 < e0, "implicit drag must decay MEKE under a stiff sink")
         if (allocated(error)) exit checks
         call check(error, e1 == e1 .and. e1 < huge(1.0_wp), "MEKE must be finite (no NaN/Inf)")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_drag_stable

   ! ------------------------------------------------------------------
   ! Test 3: harmonic-mass diffusion conserves Sum E*area (closed domain)
   ! ------------------------------------------------------------------
   subroutine test_diffusion_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 8, NY = 8, NZ = 1
      real(wp), parameter :: DX = 5000.0_wp, DT = 1800.0_wp
      real(wp) :: area, tot0, tot, drift
      integer :: i, j, it, ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp   ! no source
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp    ! no sink
         meke%cdrag = 0.0_wp
         meke%khcoeff = -1.0_wp   ! no closure feedback into kh_diff
         meke%khmeke_fac = 0.0_wp ! Kh_u = const meke%kh only
         meke%kh = 300.0_wp       ! constant lateral diffusion ⇒ Strang split
         meke%k4 = -1.0_wp

         ! A blob of eddy energy in the physical interior.
         ic = NGHOST + NX/2
         jc = NGHOST + NY/2
         meke%meke = 0.0_wp
         meke%meke(ic, jc) = 1.0_wp
         meke%meke(ic + 1, jc) = 1.0_wp

         ! Uniform area ⇒ harmonic mass is uniform; Sum E conserved.
         area = DX*DX

         call map_in(ms, gm, meke)
         tot0 = sum_phys(meke, grid)
         do it = 1, 30
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         tot = sum_phys(meke, grid)
         call map_out(ms, gm, meke)

         drift = abs(tot - tot0)/abs(tot0)
         call check(error, drift < 1.0e-10_wp, &
                    "harmonic-mass diffusion must conserve Sum E over a closed domain")
         if (allocated(error)) exit checks
         ! Non-trivial: the blob must have spread (peak dropped).
         call check(error, meke%meke(ic, jc) < 1.0_wp, &
                    "diffusion must spread the MEKE blob (peak drops)")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_diffusion_conserves

   ! Biharmonic (k4>0) lateral diffusion conserves Sum E over a closed domain
   ! AND exercises the damp_step=0.5 two-half-drag Strang branch (k4>=0 ⇒ the
   ! drag is split; with damping=cdrag=0 the halves are no-ops, isolating the
   ! biharmonic flux).  Both paths are untested by test_diffusion_conserves
   ! (which uses the Laplacian kh, not k4).
   subroutine test_biharmonic_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 8, NY = 8, NZ = 1
      real(wp), parameter :: DX = 5000.0_wp, DT = 1800.0_wp
      real(wp) :: tot0, tot, drift, peak0
      integer :: it, ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp    ! no source
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp     ! no sink (isolate biharmonic; drag halves no-op)
         meke%cdrag = 0.0_wp
         meke%khcoeff = -1.0_wp
         meke%khmeke_fac = 0.0_wp
         meke%kh = -1.0_wp         ! Laplacian OFF
         meke%k4 = 1.0e8_wp        ! biharmonic ON ⇒ damp_step=0.5 (two half-drags)

         ic = NGHOST + NX/2
         jc = NGHOST + NY/2
         meke%meke = 0.0_wp
         meke%meke(ic, jc) = 1.0_wp
         meke%meke(ic + 1, jc) = 1.0_wp
         peak0 = 1.0_wp

         call map_in(ms, gm, meke)
         tot0 = sum_phys(meke, grid)
         do it = 1, 30
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         tot = sum_phys(meke, grid)
         call map_out(ms, gm, meke)

         drift = abs(tot - tot0)/abs(tot0)
         call check(error, drift < 1.0e-10_wp, &
                    "biharmonic diffusion must conserve Sum E over a closed domain")
         if (allocated(error)) exit checks
         call check(error, meke%meke(ic, jc) < peak0, &
                    "biharmonic must spread the MEKE blob (peak drops)")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_biharmonic_conserves

   ! ------------------------------------------------------------------
   ! Rhines length: set a beta-plane f_centre, run with ONLY alpha_rhines
   ! active ⇒ Lmix == Lrhines = sqrt(Ueddy/beta) (hand-computed).
   ! ------------------------------------------------------------------
   subroutine test_rhines_length(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 7, NY = 7, NZ = 1
      real(wp), parameter :: DX = 4000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: F0 = 1.0e-4_wp, BETA = 2.0e-11_wp, E_SET = 0.05_wp
      real(wp) :: f_centre(NX + 2*NGHOST, NY + 2*NGHOST)
      real(wp) :: ueddy, beta_grad, l_expect, l_got, relerr, y
      integer :: i, j, ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp     ! no source
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp      ! E untouched by drag
         meke%cdrag = 0.0_wp        ! Lfrict=0 ⇒ barotrFac2 = 1
         meke%kh = -1.0_wp          ! no diffusion ⇒ damp_step=1, E untouched
         meke%k4 = -1.0_wp
         meke%khcoeff = -1.0_wp     ! focus on le (Lmix), not kh
         ! ONLY the Rhines scale active.
         meke%alpha_grid = 0.0_wp
         meke%alpha_deform = 0.0_wp
         meke%alpha_rhines = 1.0_wp
         meke%alpha_eady = 0.0_wp
         meke%alpha_frict = 0.0_wp
         meke%min_gamma2 = 1.0e-4_wp

         ! beta-plane |f| at cell centres: |F0 + BETA*(y - 0)|, y = (j-ng-0.5)*dy.
         do j = 1, grid%ny_total
            y = (real(j - grid%nghost, wp) - 0.5_wp)*grid%dy
            do i = 1, grid%nx_total
               f_centre(i, j) = abs(F0 + BETA*y)
            end do
         end do
         call meke%set_f_centre(grid, f_centre)
         meke%meke = E_SET

         call map_in(ms, gm, meke)
         call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         call map_out(ms, gm, meke)

         ! barotrFac2 = 1 (cdrag=0) ⇒ Ueddy = sqrt(2*E).  The beta-plane has
         ! d|f|/dy = BETA (uniform, F0+BETA*y > 0), d|f|/dx = 0
         ! ⇒ |grad f| = BETA.  Lrhines = sqrt(Ueddy/beta); alpha_rhines=1.
         ueddy = sqrt(2.0_wp*E_SET)
         beta_grad = BETA
         l_expect = sqrt(ueddy/beta_grad)
         ic = NGHOST + NX/2 + 1
         jc = NGHOST + NY/2 + 1
         l_got = meke%le(ic, jc)
         call check(error, l_got > 0.0_wp, "Lrhines must be finite/positive with alpha_rhines>0 + beta-plane f")
         if (allocated(error)) exit checks
         relerr = abs(l_got - l_expect)/abs(l_expect)
         call check(error, relerr < 1.0e-6_wp, &
                    "Lmix must equal Lrhines = sqrt(Ueddy/beta) when only alpha_rhines is active")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_rhines_length

   ! ------------------------------------------------------------------
   ! Advection conserves Sum E (uniform mass flux, closed domain) AND
   ! transports the blob; factor=0 leaves E unchanged (inertness).
   ! ------------------------------------------------------------------
   subroutine test_advection_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 10, NY = 6, NZ = 1
      real(wp), parameter :: DX = 5000.0_wp, DT = 600.0_wp
      real(wp) :: tot0, tot, drift, e_left0, e_left1
      integer :: it, ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp
         meke%cdrag = 0.0_wp
         meke%khcoeff = -1.0_wp
         meke%kh = -1.0_wp          ! no diffusion ⇒ damp_step=1
         meke%k4 = -1.0_wp
         meke%advection_factor = 0.5_wp  ! advection ON

         ! Uniform eastward volume transport on interior u-faces (m^3/s).
         ! mass_flux_x_layer(i,j,k) = u*h_face*dy ; set a small CFL-safe value.
         ms%mass_flux_x_layer = 50.0_wp   ! all faces; edge faces zeroed by kernel gate
         ms%mass_flux_y_layer = 0.0_wp

         ! A blob in the interior.
         ic = NGHOST + 3
         jc = NGHOST + NY/2
         meke%meke = 0.0_wp
         meke%meke(ic, jc) = 1.0_wp
         meke%meke(ic + 1, jc) = 1.0_wp

         call map_in(ms, gm, meke)
         ! mass_flux_*_layer is `create` (not copyin) in ms enter_data ⇒ push
         ! the host values onto the device explicitly.
         !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
         tot0 = sum_phys(meke, grid)
         !$acc update self(meke%meke)
         e_left0 = meke%meke(ic, jc)
         do it = 1, 20
            call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         end do
         tot = sum_phys(meke, grid)
         !$acc update self(meke%meke)
         e_left1 = meke%meke(ic, jc)
         call map_out(ms, gm, meke)

         drift = abs(tot - tot0)/abs(tot0)
         call check(error, drift < 1.0e-10_wp, &
                    "upwind advection must conserve Sum E over a closed domain")
         if (allocated(error)) exit checks
         ! Eastward flow drains the upstream cell (peak at ic drops).
         call check(error, e_left1 < e_left0, &
                    "eastward advection must transport E downstream (upstream cell drains)")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_advection_conserves

   subroutine test_advection_inert(error)
      !! advection_factor = 0 ⇒ the advection stage is a no-op (E unchanged by
      !! transport even with a non-zero mass flux present).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 8, NY = 6, NZ = 1
      real(wp), parameter :: DX = 5000.0_wp, DT = 600.0_wp
      real(wp) :: e_before, e_after
      integer :: ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp
         meke%cdrag = 0.0_wp
         meke%khcoeff = -1.0_wp
         meke%kh = -1.0_wp
         meke%k4 = -1.0_wp
         meke%advection_factor = 0.0_wp  ! advection OFF (default)

         ms%mass_flux_x_layer = 50.0_wp  ! non-zero flux present but must be ignored
         ms%mass_flux_y_layer = 0.0_wp

         ic = NGHOST + 3
         jc = NGHOST + NY/2
         meke%meke = 0.0_wp
         meke%meke(ic, jc) = 1.0_wp

         call map_in(ms, gm, meke)
         e_before = 1.0_wp
         call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         !$acc update self(meke%meke)
         e_after = meke%meke(ic, jc)
         call map_out(ms, gm, meke)

         call check(error, abs(e_after - e_before) < 1.0e-14_wp, &
                    "advection_factor=0 ⇒ advection stage is inert (E unchanged)")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_advection_inert

   function sum_phys(meke, grid) result(s)
      !! Sum of E over the FULL array (uniform area+mass ⇒ proportional to
      !! Sum E*area*mass).  Flux-form diffusion vanishes only at the ARRAY
      !! edges (where the face metrics zero), so the conserved quantity is the
      !! full-array sum, not the physical-interior sum (the interior exchanges
      !! with the ghost halo across its inner faces).  Pulls device first.
      type(ocean_meke_t), intent(in) :: meke
      type(hgrid_t), intent(in) :: grid
      real(wp) :: s
      integer :: i, j
      !$acc update self(meke%meke)
      s = 0.0_wp
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            s = s + meke%meke(i, j)
         end do
      end do
   end function sum_phys

   ! ------------------------------------------------------------------
   ! Test 4: kh closure = khcoeff*sqrt(2*gamma_t2*E)*Lmix (single scale)
   ! ------------------------------------------------------------------
   subroutine test_kh_closure(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 4000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: KHCOEFF = 0.3_wp, E_SET = 0.05_wp
      real(wp) :: lgrid, gt2, ueddy, kh_expect, kh_got, relerr
      integer :: i, j, ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp    ! no source bump
         meke%bgsrc = 0.0_wp
         meke%damping = 0.0_wp     ! no drag change to E
         meke%cdrag = 0.0_wp       ! ⇒ Lfrict=0 ⇒ gamma factors = floor/1
         meke%kh = -1.0_wp         ! no diffusion ⇒ damp_step=1, E untouched
         meke%k4 = -1.0_wp
         meke%khcoeff = KHCOEFF
         ! ONLY the grid length scale is active: alpha_grid>0, rest 0.
         meke%alpha_grid = 1.0_wp
         meke%alpha_deform = 0.0_wp
         meke%alpha_rhines = 0.0_wp
         meke%alpha_eady = 0.0_wp
         meke%alpha_frict = 0.0_wp
         ! cdrag=0 ⇒ Lfrict=0 ⇒ bottomFac2 floor, barotrFac2 = 1 (Ct*Lfrict=0).
         meke%min_gamma2 = 1.0e-4_wp

         meke%meke = E_SET

         call map_in(ms, gm, meke)
         call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
         call map_out(ms, gm, meke)

         ! With cdrag=0 ⇒ Lfrict=0 ⇒ Ct*Lfrict=0 ⇒ barotrFac2 = 1 (no floor
         ! needed since 1 > min_gamma2).  Single length scale = Lgrid =
         ! sqrt(area) = DX (square cells), alpha_grid=1 ⇒ Lmix = Lgrid.
         lgrid = sqrt(DX*DX)
         gt2 = 1.0_wp
         ueddy = sqrt(2.0_wp*gt2*E_SET)
         kh_expect = KHCOEFF*ueddy*lgrid
         ic = NGHOST + 1
         jc = NGHOST + 1
         kh_got = meke%kh_diff(ic, jc)
         relerr = abs(kh_got - kh_expect)/abs(kh_expect)
         call check(error, relerr < 1.0e-9_wp, &
                    "kh = khcoeff*sqrt(2*gamma_t2*E)*Lgrid (single-scale closure)")
         if (allocated(error)) exit checks
         call check(error, abs(meke%le(ic, jc) - lgrid) < 1.0e-6_wp, &
                    "Lmix must equal Lgrid when only alpha_grid is active")
      end block checks
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_kh_closure

   ! ------------------------------------------------------------------
   ! Test 5: feedback seam adds geom-mean kh into varmix%khth (fac>0);
   !         fac=0 ⇒ varmix%khth unchanged (bit-identity of the seam).
   ! ------------------------------------------------------------------
   subroutine test_feeds_khth(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      type(ocean_varmix_t) :: vm
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 4000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: KHTH_FAC = 0.5_wp, BASE = 100.0_wp
      real(wp) :: khu_before, khu_after, expect
      integer :: ic, jc
      checks: block
         call make_grid(grid, NX, NY, DX)
         call setup_ms(ms, grid, NZ, DZ)
         call make_cartesian_metrics(metrics, grid)
         call setup_gm(gm, grid, NZ, 0.0_wp)
         call vm%init(grid, nz_ml=NZ)
         vm%enable = .true.
         vm%khth_u = BASE
         vm%khth_v = BASE
         vm%khtr_u = BASE
         vm%khtr_v = BASE
         ! VarMix SN faces are zero (not set) ⇒ Eady scale inert; fine.
         vm%sn_u = 0.0_wp
         vm%sn_v = 0.0_wp

         call meke%init(grid, nz_ml=NZ)
         meke%enable = .true.
         meke%gmcoeff = -1.0_wp
         meke%damping = 0.0_wp
         meke%cdrag = 0.0_wp
         meke%kh = -1.0_wp
         meke%k4 = -1.0_wp
         meke%khcoeff = 0.3_wp
         meke%alpha_grid = 1.0_wp
         meke%khth_fac = 0.0_wp   ! FIRST: seam inert.
         meke%meke = 0.05_wp

         ! ---- Pass A: khth_fac = 0 ⇒ varmix%khth unchanged. ----
         ! NOTE: the enter_data() TBP MUST be on its own line — a `; call`
         ! after an !$acc directive is swallowed into the directive comment
         ! and never executes (the component arrays then never attach).
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(gm)
         call gm%enter_data()
         !$acc enter data copyin(vm)
         call vm%enter_data()
         !$acc enter data copyin(meke)
         call meke%enter_data()
         call meke_step(grid, metrics, meke, gm, varmix=vm, ms=ms, dt=DT)
         !$acc update self(vm%khth_u, vm%khth_v, meke%kh_diff)
         ic = NGHOST + 2
         jc = NGHOST + 1
         khu_before = vm%khth_u(ic, jc)
         call check(error, abs(khu_before - BASE) < 1.0e-9_wp, &
                    "khth_fac=0 ⇒ MEKE feedback seam leaves varmix%khth unchanged")
         if (allocated(error)) then
            call teardown_A(ms, gm, vm, meke); exit checks
         end if

         ! ---- Pass B: khth_fac > 0 ⇒ varmix%khth picks up geom-mean kh. ----
         meke%khth_fac = KHTH_FAC
         !$acc update device(meke%khth_fac)
         ! reset varmix base to BASE on device.
         vm%khth_u = BASE; vm%khth_v = BASE
         !$acc update device(vm%khth_u, vm%khth_v)
         call meke_step(grid, metrics, meke, gm, varmix=vm, ms=ms, dt=DT)
         !$acc update self(vm%khth_u, meke%kh_diff)
         khu_after = vm%khth_u(ic, jc)
         ! geom mean of kh_diff(ic-1,jc) and kh_diff(ic,jc) (uniform field).
         expect = BASE + KHTH_FAC*sqrt(meke%kh_diff(ic - 1, jc)*meke%kh_diff(ic, jc))
         call teardown_A(ms, gm, vm, meke)

         call check(error, khu_after > khu_before + 1.0e-9_wp, &
                    "khth_fac>0 ⇒ MEKE feedback must raise varmix%khth")
         if (allocated(error)) exit checks
         call check(error, abs(khu_after - expect) < 1.0e-6_wp, &
                    "MEKE seam must add khth_fac*sqrt(kh_i*kh_{i+1}) into varmix%khth")
      end block checks
      call meke%destroy()
      call vm%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_feeds_khth

   subroutine teardown_A(ms, gm, vm, meke)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_gm_t), intent(inout) :: gm
      type(ocean_varmix_t), intent(inout) :: vm
      type(ocean_meke_t), intent(inout) :: meke
      call meke%exit_data()
      !$acc exit data delete(meke)
      call vm%exit_data()
      !$acc exit data delete(vm)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine teardown_A

end module test_ocean_meke
