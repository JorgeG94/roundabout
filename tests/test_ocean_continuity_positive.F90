!! Unit tests for the positive-definite continuity limiter (P2/P3).
!!
!! Exercises the production split entry `continuity_tracer_step_split` on the
!! GPU (`-gpu=mem:separate`) build — the whole point is device data motion, so
!! every state object AND `ct` is mapped, and every host read is preceded by an
!! `!$acc update self` of the COMPONENT arrays (never the aggregate derived
!! type — that would clobber the host descriptors).  Directives are inert
!! no-ops on the host/multicore build, so they are unconditional.
!!
!! Three properties (plan P4):
!!   (1) overdraw_trap — a prescribed `uhbt` demands far more transport than an
!!       interior column holds; through the renorm path the limiter keeps every
!!       h_layer >= h_lim, conserves the DOMAIN mass to round-off (closed
!!       walls), engages (`n_limited_step > 0`), and a uniform tracer stays
!!       uniform (CWC under limiting).
!!   (2) theta_one_bitident — the SAME benign setup, knob OFF vs knob ON with
!!       gentle fluxes that limit nothing: h_layer + hTr must be BIT-identical
!!       (`==`) and `n_limited_step == 0` on the ON run (guards against an
!!       always-on perturbation of the untriggered path).
!!   (3) below_floor_freeze — one interior column initialised AT h_lim with an
!!       outgoing flux: its outflux is zeroed (theta = 0), its h stays >= h_lim,
!!       and its inflow face (donor healthy, theta = 1) is untouched — inflow is
!!       never blocked.
module test_ocean_continuity_positive
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, continuity_tracer_step_split, &
                             continuity_zonal_flux, pd_limit_zonal_impl
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_continuity_positive_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: DT = 100.0_wp

contains

   subroutine collect_ocean_continuity_positive_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overdraw_trap", test_overdraw_trap), &
                  new_unittest("theta_one_bitident", test_theta_one_bitident), &
                  new_unittest("below_floor_freeze", test_below_floor_freeze), &
                  new_unittest("ucor_rematch", test_ucor_rematch) &
                  ]
   end subroutine collect_ocean_continuity_positive_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, DX, DX)
   end subroutine make_grid

   ! ------------------------------------------------------------------
   ! (1) Overdraw trap through the uhbt renorm path
   ! ------------------------------------------------------------------
   subroutine test_overdraw_trap(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), allocatable :: uhbt(:, :)
      real(wp), parameter :: H0 = 1.0_wp
      real(wp), parameter :: HLIM = 0.9_wp
      real(wp), parameter :: T0 = 10.0_wp
      real(wp), parameter :: TOL = 1.0e-11_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: nx, ny, it, jt, kt
      real(wp) :: mass0, mass1, hmin_out, tr_dev, tr

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ! Positive-definite knob ON with a high floor so even a CFL-bracketed
         ! renorm velocity over-drains (avail = H0 - HLIM = 0.1).
         ct%positive_definite = .true.
         ct%h_lim = HLIM

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = T0*H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0

         ! uhbt: a huge zonal transport demand at every east face — the renorm
         ! saturates each layer at the CFL cap; the limiter then caps realised
         ! outflow at avail.  Wall faces are skipped by the renorm.
         allocate (uhbt(nx + 1, ny), source=1.0e9_wp)

         ! Domain mass BEFORE (interior cells only; areaT = DX*DX uniform).
         mass0 = sum(ms%h_layer(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*DX*DX

         !$acc enter data copyin(ms, ct)
         call ms%enter_data()
         call ct%enter_data()
         !$acc enter data copyin(uhbt)

         call continuity_tracer_step_split(grid, metrics, ct, ms, DT, uhbt=uhbt)

         !$acc update self(ms%h_layer, ms%tracers(ms%idx_temperature)%hTr)
         !$acc exit data delete(uhbt)
         call ct%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ct)

         ! Domain mass AFTER.
         mass1 = sum(ms%h_layer(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*DX*DX

         ! (a) positivity: every interior layer stays >= h_lim.
         hmin_out = minval(ms%h_layer(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))
         call check(error, hmin_out >= HLIM - TOL, &
                    "overdraw: every h_layer must stay >= h_lim after limiting")
         if (allocated(error)) exit checks

         ! (b) zero mass created: domain total conserved to round-off (walls closed).
         call check(error, abs(mass1 - mass0) <= TOL*mass0, &
                    "overdraw: domain mass must be conserved to round-off")
         if (allocated(error)) exit checks

         ! (c) the limiter actually engaged.
         call check(error, ct%n_limited_step > 0, &
                    "overdraw: n_limited_step must be > 0 (limiter engaged)")
         if (allocated(error)) exit checks

         ! (d) CWC: a uniform tracer stays uniform through the limited step.
         tr_dev = 0.0_wp
         do kt = 1, NZ
            do jt = NGHOST + 1, ny - NGHOST
               do it = NGHOST + 1, nx - NGHOST
                  tr = ms%tracers(ms%idx_temperature)%hTr(it, jt, kt)/ms%h_layer(it, jt, kt)
                  tr_dev = max(tr_dev, abs(tr - T0))
               end do
            end do
         end do
         call check(error, tr_dev <= TOL*T0, &
                    "overdraw: uniform tracer must stay uniform (CWC under limiting)")
      end block checks

      if (allocated(uhbt)) deallocate (uhbt)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_overdraw_trap

   ! ------------------------------------------------------------------
   ! (2) theta == 1 non-perturbation: knob OFF == knob ON, bit-identical
   ! ------------------------------------------------------------------
   subroutine test_theta_one_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      real(wp), parameter :: H0 = 1.0_wp
      real(wp), parameter :: T0 = 10.0_wp
      real(wp), parameter :: UGENTLE = 0.1_wp   !! 0.01*DX/DT — demand << avail
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: nx, ny
      logical :: h_ident, tr_ident

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total

         ms_a%nz_ml = NZ
         ms_b%nz_ml = NZ
         call ms_a%init(grid)
         call ms_b%init(grid)
         call ct_a%init(grid, nz_ml=NZ)
         call ct_b%init(grid, nz_ml=NZ)

         ! A = knob OFF (reference), B = knob ON, h_lim = 0.  Identical ICs.
         ct_a%positive_definite = .false.
         ct_b%positive_definite = .true.
         ct_b%h_lim = 0.0_wp

         ms_a%h_layer = H0
         ms_a%u_face_x_layer = UGENTLE
         ms_a%v_face_y_layer = UGENTLE
         ms_a%tracers(ms_a%idx_salinity)%hTr = T0*H0
         ms_a%tracers(ms_a%idx_temperature)%hTr = T0*H0
         ms_b%h_layer = H0
         ms_b%u_face_x_layer = UGENTLE
         ms_b%v_face_y_layer = UGENTLE
         ms_b%tracers(ms_b%idx_salinity)%hTr = T0*H0
         ms_b%tracers(ms_b%idx_temperature)%hTr = T0*H0

         !$acc enter data copyin(ms_a, ms_b, ct_a, ct_b)
         call ms_a%enter_data()
         call ms_b%enter_data()
         call ct_a%enter_data()
         call ct_b%enter_data()

         call continuity_tracer_step_split(grid, metrics, ct_a, ms_a, DT)
         call continuity_tracer_step_split(grid, metrics, ct_b, ms_b, DT)

         !$acc update self(ms_a%h_layer, ms_a%tracers(ms_a%idx_temperature)%hTr)
         !$acc update self(ms_b%h_layer, ms_b%tracers(ms_b%idx_temperature)%hTr)
         call ct_a%exit_data()
         call ct_b%exit_data()
         call ms_a%exit_data()
         call ms_b%exit_data()
         !$acc exit data delete(ms_a, ms_b, ct_a, ct_b)

         ! Bit-identity (==, not tolerance): the knob-ON untriggered path must
         ! not perturb a single bit relative to knob-OFF.
         h_ident = all(ms_a%h_layer == ms_b%h_layer)
         tr_ident = all(ms_a%tracers(ms_a%idx_temperature)%hTr == &
                        ms_b%tracers(ms_b%idx_temperature)%hTr)

         call check(error, ct_b%n_limited_step == 0, &
                    "theta==1: n_limited_step must be 0 on the benign ON run")
         if (allocated(error)) exit checks
         call check(error, h_ident, &
                    "theta==1: h_layer must be BIT-identical OFF vs ON")
         if (allocated(error)) exit checks
         call check(error, tr_ident, &
                    "theta==1: hTr must be BIT-identical OFF vs ON")
      end block checks

      call ct_a%destroy()
      call ct_b%destroy()
      call ms_a%destroy()
      call ms_b%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_theta_one_bitident

   ! ------------------------------------------------------------------
   ! (3) below-floor freeze: at-floor donor zeroes its outflux, inflow untouched
   ! ------------------------------------------------------------------
   subroutine test_below_floor_freeze(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), parameter :: H0 = 1.0_wp
      real(wp), parameter :: HLIM = 0.5_wp
      real(wp), parameter :: T0 = 10.0_wp
      real(wp), parameter :: UEAST = 4.0_wp   !! 0.4*DX/DT — demand 0.4 < avail 0.5 (healthy)
      real(wp), parameter :: TOL = 1.0e-11_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: nx, ny, i0, j0
      real(wp) :: flux_out, flux_in, h_frozen

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ct%positive_definite = .true.
         ct%h_lim = HLIM

         ! Uniform h, uniform eastward flow; one interior column pinned AT the
         ! floor (avail = 0) so its outflux must be frozen.
         i0 = NGHOST + 4
         j0 = NGHOST + 3
         ms%h_layer = H0
         ms%h_layer(i0, j0, :) = HLIM
         ms%u_face_x_layer = UEAST
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = T0*ms%h_layer
         ms%tracers(ms%idx_temperature)%hTr = T0*ms%h_layer

         !$acc enter data copyin(ms, ct)
         call ms%enter_data()
         call ct%enter_data()

         call continuity_tracer_step_split(grid, metrics, ct, ms, DT)

         !$acc update self(ms%h_layer, ms%mass_flux_x_layer)
         call ct%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ct)

         ! Frozen column's OUTflux (east face i0+1, donor = frozen cell, theta=0).
         flux_out = maxval(abs(ms%mass_flux_x_layer(i0 + 1, j0, :)))
         ! Frozen column's INflux (west face i0, donor = healthy west cell, theta=1).
         flux_in = minval(ms%mass_flux_x_layer(i0, j0, :))
         h_frozen = minval(ms%h_layer(i0, j0, :))

         call check(error, ct%n_limited_step > 0, &
                    "freeze: n_limited_step must be > 0 (at-floor donor limited)")
         if (allocated(error)) exit checks
         call check(error, flux_out <= TOL, &
                    "freeze: at-floor column's outflux must be zeroed (theta=0)")
         if (allocated(error)) exit checks
         call check(error, flux_in > 0.0_wp, &
                    "freeze: inflow face (healthy donor) must NOT be blocked")
         if (allocated(error)) exit checks
         call check(error, h_frozen >= HLIM - TOL, &
                    "freeze: at-floor column's h must stay >= h_lim")
      end block checks

      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_below_floor_freeze

   ! ------------------------------------------------------------------
   ! (4) v1.1 u_cor re-matching: the corrector velocity is scaled by the
   !     SAME per-face θ as the mass flux (flux/velocity stay consistent).
   ! ------------------------------------------------------------------
   subroutine test_ucor_rematch(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), allocatable :: uhbt(:, :), u_cor(:, :, :)
      real(wp), allocatable :: flux_pre(:, :, :), flux_post(:, :, :)
      real(wp), allocatable :: ucor_pre(:, :, :), ucor_post(:, :, :)
      real(wp), parameter :: H0 = 1.0_wp
      real(wp), parameter :: HLIM = 0.9_wp
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: nx, ny, i, j, k
      real(wp) :: rf, ru, max_ratio_dev, min_ratio_flux

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ct%positive_definite = .true.
         ct%h_lim = HLIM

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%tracers(ms%idx_salinity)%hTr = 10.0_wp*H0
         ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*H0

         allocate (uhbt(nx + 1, ny), source=1.0e9_wp)
         allocate (u_cor(nx + 1, ny, NZ), source=0.0_wp)
         allocate (flux_pre(nx + 1, ny, NZ), flux_post(nx + 1, ny, NZ))
         allocate (ucor_pre(nx + 1, ny, NZ), ucor_post(nx + 1, ny, NZ))

         !$acc enter data copyin(ms, ct)
         call ms%enter_data()
         call ct%enter_data()
         !$acc enter data copyin(uhbt)
         !$acc enter data create(u_cor)

         ! Renorm fill: captures the transport-matched u_cor + renormalised
         ! mass_flux_x_layer BEFORE the limiter (the mom6-scheme capture point).
         call continuity_zonal_flux(grid, metrics, ct, ms, DT, uhbt=uhbt, u_cor=u_cor)
         !$acc update self(ms%mass_flux_x_layer, u_cor)
         flux_pre = ms%mass_flux_x_layer
         ucor_pre = u_cor

         ! The P2 limiter, WITH u_cor forwarded — scales both by the same θ.
         call pd_limit_zonal_impl(nx, ny, NZ, DT, ct%h_lim, metrics%iareaT, &
                                  ms%h_layer, ms%mass_flux_x_layer, &
                                  ct%pd_theta%data, ct%n_limited_step, u_cor=u_cor)
         !$acc update self(ms%mass_flux_x_layer, u_cor)
         flux_post = ms%mass_flux_x_layer
         ucor_post = u_cor

         !$acc exit data delete(u_cor, uhbt)
         call ct%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, ct)

         ! On every touched face (flux_pre /= 0), the velocity ratio must equal
         ! the flux ratio (both = the same θ_face); walls (flux_pre = 0) skipped.
         max_ratio_dev = 0.0_wp
         min_ratio_flux = 1.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 2, nx
                  if (abs(flux_pre(i, j, k)) > 1.0e-30_wp .and. &
                      abs(ucor_pre(i, j, k)) > 1.0e-30_wp) then
                     rf = flux_post(i, j, k)/flux_pre(i, j, k)
                     ru = ucor_post(i, j, k)/ucor_pre(i, j, k)
                     max_ratio_dev = max(max_ratio_dev, abs(ru - rf))
                     min_ratio_flux = min(min_ratio_flux, rf)
                  end if
               end do
            end do
         end do

         call check(error, ct%n_limited_step > 0, &
                    "ucor: limiter must engage (n_limited_step > 0)")
         if (allocated(error)) exit checks
         call check(error, min_ratio_flux < 1.0_wp, &
                    "ucor: at least one face must actually be scaled (θ < 1)")
         if (allocated(error)) exit checks
         call check(error, max_ratio_dev <= TOL, &
                    "ucor: u_cor must scale by the SAME θ as the mass flux on every face")
      end block checks

      if (allocated(uhbt)) deallocate (uhbt)
      if (allocated(u_cor)) deallocate (u_cor)
      if (allocated(flux_pre)) deallocate (flux_pre, flux_post, ucor_pre, ucor_post)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_ucor_rematch

end module test_ocean_continuity_positive
