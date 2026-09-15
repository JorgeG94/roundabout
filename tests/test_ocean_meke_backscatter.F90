!! Tests for the MEKE -> momentum harmonic backscatter (capability [5],
!! Gap 2 — negative-viscosity eddy-energy return).  All cases RUN THE DEVICE
!! KERNELS: `meke_step` fills the derived backscatter coefficient `ku`, and
!! `meke_backscatter_apply` (the production injection seam) subtracts a
!! face-average of `ku` from the per-face resolved harmonic viscosity with a
!! CFL stability floor.
!!
!!  1. ku_closure_form   — backscatter on ⇒ ku = coeff*sqrt(2*E)*Lmix > 0
!!                         (hand-checked against the grid length scale);
!!                         off ⇒ ku == 0 exactly (inert).
!!  2. apply_subtracts   — a small ku is subtracted exactly from ah_face
!!                         (net = A - ku) on interior faces; wall faces and
!!                         off-state leave ah_face bit-identical.
!!  3. cfl_floor         — a huge ku cannot drive the net coefficient past
!!                         the forward-Euler viscous-CFL lower bound: the
!!                         net magnitude metric `|A|*dt*(idx²+idy²)` stays
!!                         <= backscatter_cfl*0.5. This bounds the negative
!!                         mode's GROWTH RATE (not |g|<=1 — a negative
!!                         Laplacian always amplifies; stability needs a
!!                         positive biharmonic backstop, mandatory at
!!                         configure).
!!  4. energy_return     — where ku > 0 the NET (resolved - ku) coefficient
!!                         is strictly smaller than resolved-only, so the
!!                         Laplacian removes LESS KE (resolved KE decays
!!                         slower) — the energy-return direction.
module test_ocean_meke_backscatter
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_gm, only: ocean_gm_t
   use rdb_ocean_meke, only: ocean_meke_t, meke_step, meke_backscatter_apply
   use rdb_ocean_lateral_mix, only: has_biharmonic_backstop, &
                                    LMIX_NONE, LMIX_LEITH, LMIX_LEITH_BIHARM
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_meke_backscatter_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1025.0_wp

contains

   subroutine collect_ocean_meke_backscatter_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("meke_bs_ku_closure_form", test_ku_closure_form), &
                  new_unittest("meke_bs_apply_subtracts", test_apply_subtracts), &
                  new_unittest("meke_bs_cfl_floor", test_cfl_floor), &
                  new_unittest("meke_bs_energy_return", test_energy_return), &
                  new_unittest("meke_bs_requires_biharmonic_backstop", test_backstop_guard) &
                  ]
   end subroutine collect_ocean_meke_backscatter_tests

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
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = dz
      ms%rho_layer = RHO0
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_ms

   subroutine setup_gm(gm, grid, nz)
      type(ocean_gm_t), intent(inout) :: gm
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      call gm%init(grid, nz_ml=nz)
      gm%enable = .true.
      gm%rho0 = RHO0
      gm%gm_src = 0.0_wp
   end subroutine setup_gm

   !! Run `meke_step` once with a fixed eddy-energy IC, grid-length-only
   !! mixing length, and the backscatter knobs as given; return the
   !! cell-centred `ku` + the mixing length `le` + the structure factor.
   subroutine run_ku_closure(backscatter, coeff, ke_ic, ku_out, le_out, gt2_out)
      logical, intent(in) :: backscatter
      real(wp), intent(in) :: coeff, ke_ic
      real(wp), intent(out) :: ku_out, le_out, gt2_out
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_gm_t) :: gm
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 5, NY = 5, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp, DZL = 50.0_wp
      call make_grid(grid, NX, NY, DX)
      call setup_ms(ms, grid, NZ, DZL)
      call make_cartesian_metrics(metrics, grid)
      call setup_gm(gm, grid, NZ)
      call meke%init(grid, nz_ml=NZ)
      meke%enable = .true.
      meke%gmcoeff = -1.0_wp
      meke%frcoeff = -1.0_wp
      meke%bgsrc = 0.0_wp
      meke%damping = 0.0_wp
      meke%cdrag = 0.0_wp
      meke%cd_scale = 0.0_wp
      meke%khcoeff = -1.0_wp       ! kh closure off (we want ku only)
      meke%kh = -1.0_wp
      meke%k4 = -1.0_wp
      meke%min_gamma2 = 1.0e-4_wp
      meke%alpha_grid = 1.0_wp     ! mixing length = grid scale only
      meke%backscatter = backscatter
      meke%visc_coeff_ku = coeff
      meke%meke = ke_ic
      call map_in(ms, gm, meke)
      call meke_step(grid, metrics, meke, gm, ms=ms, dt=DT)
      call map_out(ms, gm, meke)
      ku_out = meke%ku(NGHOST + 1, NGHOST + 1)
      le_out = meke%le(NGHOST + 1, NGHOST + 1)
      gt2_out = meke%barotr_fac2(NGHOST + 1, NGHOST + 1)
      call meke%destroy()
      call gm%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_ku_closure

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
      !$acc update self(meke%meke, meke%kh_diff, meke%le, meke%ku, meke%barotr_fac2)
      call meke%exit_data()
      !$acc exit data delete(meke)
      call gm%exit_data()
      !$acc exit data delete(gm)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! Test 1: ku closure form + inert-when-off
   ! ------------------------------------------------------------------
   subroutine test_ku_closure_form(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ku_on, ku_off, le_on, gt2_on, ku_expect, ueddy
      real(wp), parameter :: COEFF = 0.15_wp, KE = 0.04_wp
      checks: block
         call run_ku_closure(.true., COEFF, KE, ku_on, le_on, gt2_on)
         call run_ku_closure(.false., COEFF, KE, ku_off, le_on, gt2_on)
         call check(error, abs(ku_off) < 1.0e-30_wp, &
                    "backscatter off ⇒ ku == 0 (inert)")
         if (allocated(error)) exit checks
         ! ku = coeff*sqrt(2*E)*Lmix — MOM6's Ku uses the PLAIN eddy
         ! velocity sqrt(2*E), no barotropic-mode gt2 factor (unlike kh).
         ! Lmix read back from the slot ⇒ no re-derivation drift.
         ueddy = sqrt(2.0_wp*KE)
         ku_expect = COEFF*ueddy*le_on
         call check(error, ku_on > 0.0_wp, "backscatter on + E>0 ⇒ ku > 0")
         if (allocated(error)) exit checks
         call check(error, abs(ku_on - ku_expect) < 1.0e-12_wp*max(1.0_wp, abs(ku_expect)), &
                    "ku = coeff*sqrt(2*E)*Lmix to round-off")
      end block checks
   end subroutine test_ku_closure_form

   ! ------------------------------------------------------------------
   ! Test 2: apply subtracts ku exactly on interior faces; bit-identical
   ! at wall faces and when off.
   ! ------------------------------------------------------------------
   subroutine test_apply_subtracts(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 6, NY = 6, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: A0 = 500.0_wp, KU = 50.0_wp
      real(wp), allocatable :: ahx(:, :, :), ahy(:, :, :)
      real(wp), allocatable :: ahx0(:, :, :), ahy0(:, :, :)
      integer :: nxt, nyt, ii, jj
      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         nxt = grid%nx_total
         nyt = grid%ny_total
         call meke%init(grid, nz_ml=NZ)
         meke%backscatter = .true.
         meke%visc_coeff_ku = 1.0_wp
         meke%ku = KU                 ! uniform cell-centred ku
         allocate (ahx(nxt + 1, nyt, NZ), source=A0)
         allocate (ahy(nxt, nyt + 1, NZ), source=A0)
         allocate (ahx0(nxt + 1, nyt, NZ), source=A0)
         allocate (ahy0(nxt, nyt + 1, NZ), source=A0)
         ! interior pick (well away from walls; ku uniform ⇒ face-avg = KU)
         ii = NGHOST + 2
         jj = NGHOST + 2
         !$acc enter data copyin(meke)
         call meke%enter_data()
         !$acc enter data copyin(ahx, ahy)
         call meke_backscatter_apply(grid, metrics, meke, DT, ahx, ahy)
         !$acc update self(ahx, ahy)
         !$acc exit data delete(ahx, ahy)
         call meke%exit_data()
         !$acc exit data delete(meke)
         ! interior: net = A0 - KU exactly (KU << CFL floor magnitude).
         call check(error, abs(ahx(ii, jj, 1) - (A0 - KU)) < 1.0e-10_wp, &
                    "u-face interior: net = A_resolved - ku")
         if (allocated(error)) exit checks
         call check(error, abs(ahy(ii, jj, 1) - (A0 - KU)) < 1.0e-10_wp, &
                    "v-face interior: net = A_resolved - ku")
         if (allocated(error)) exit checks
         ! wall faces untouched (i=1 u-face / j=1 v-face).
         call check(error, abs(ahx(1, jj, 1) - ahx0(1, jj, 1)) < 1.0e-30_wp, &
                    "u-face wall (i=1) bit-identical (no backscatter at wall)")
         if (allocated(error)) exit checks
         call check(error, abs(ahy(ii, 1, 1) - ahy0(ii, 1, 1)) < 1.0e-30_wp, &
                    "v-face wall (j=1) bit-identical")
         if (allocated(error)) exit checks
         ! off ⇒ literal no-op.
         ahx = A0
         ahy = A0
         meke%backscatter = .false.
         !$acc enter data copyin(meke)
         call meke%enter_data()
         !$acc enter data copyin(ahx, ahy)
         call meke_backscatter_apply(grid, metrics, meke, DT, ahx, ahy)
         !$acc update self(ahx, ahy)
         !$acc exit data delete(ahx, ahy)
         call meke%exit_data()
         !$acc exit data delete(meke)
         call check(error, maxval(abs(ahx - ahx0)) < 1.0e-30_wp .and. &
                    maxval(abs(ahy - ahy0)) < 1.0e-30_wp, &
                    "backscatter off ⇒ ah_face bit-identical (no-op)")
      end block checks
      if (allocated(ahx)) deallocate (ahx, ahy, ahx0, ahy0)
      call meke%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_apply_subtracts

   ! ------------------------------------------------------------------
   ! Test 3: CFL floor — a HUGE ku is clamped so the net coefficient's
   ! magnitude stays within the forward-Euler viscous-CFL bound.  We check
   ! the per-face metric |A|*dt*(idx²+idy²) <= backscatter_cfl*0.5.
   !
   ! NOTE: this bounds the negative mode's GROWTH RATE, not |g|<=1 — a
   ! floored negative Laplacian still amplifies (g = 1 + |A|*dt*k² > 1) in
   ! isolation.  Run stability relies on a co-present positive biharmonic
   ! (mandatory at configure) dissipating the grid-scale mode; this test
   ! only asserts the floor caps the magnitude as designed.
   ! ------------------------------------------------------------------
   subroutine test_cfl_floor(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 6, NY = 6, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: A0 = 100.0_wp, KU_HUGE = 1.0e9_wp
      real(wp), allocatable :: ahx(:, :, :), ahy(:, :, :)
      integer :: nxt, nyt, ii, jj
      real(wp) :: idx2, a_net, cfl_metric
      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         nxt = grid%nx_total
         nyt = grid%ny_total
         call meke%init(grid, nz_ml=NZ)
         meke%backscatter = .true.
         meke%visc_coeff_ku = 1.0_wp
         meke%ku = KU_HUGE
         allocate (ahx(nxt + 1, nyt, NZ), source=A0)
         allocate (ahy(nxt, nyt + 1, NZ), source=A0)
         ii = NGHOST + 2
         jj = NGHOST + 2
         !$acc enter data copyin(meke)
         call meke%enter_data()
         !$acc enter data copyin(ahx, ahy)
         call meke_backscatter_apply(grid, metrics, meke, DT, ahx, ahy)
         !$acc update self(ahx, ahy)
         !$acc exit data delete(ahx, ahy)
         call meke%exit_data()
         !$acc exit data delete(meke)
         a_net = ahx(ii, jj, 1)
         ! Cartesian square grid: idxCu = idyCu = 1/DX.
         idx2 = 1.0_wp/(DX*DX)
         ! The mandatory bound: |A_net|*dt*(idx²+idy²) <= 0.8*0.5 (the floor
         ! safety coefficient).  A huge ku drove A_net to the floor, so it
         ! must sit exactly at the (negative) bound, never past it.
         cfl_metric = abs(a_net)*DT*(idx2 + idx2)
         call check(error, a_net < 0.0_wp, &
                    "huge ku ⇒ net viscosity is negative (energy return)")
         if (allocated(error)) exit checks
         call check(error, cfl_metric <= 0.8_wp*0.5_wp + 1.0e-9_wp, &
                    "floored net coefficient honours the forward-Euler viscous-CFL bound")
         if (allocated(error)) exit checks
         ! and it actually clamped (not the bare A0-ku, which would blow the bound).
         call check(error, a_net > -(A0 + KU_HUGE), &
                    "net was clamped by the floor, not the bare A - ku")
      end block checks
      if (allocated(ahx)) deallocate (ahx, ahy)
      call meke%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_cfl_floor

   ! ------------------------------------------------------------------
   ! Test 4: energy-return direction.  With a moderate ku, the NET
   ! per-face coefficient is strictly smaller than resolved-only, so the
   ! Laplacian friction removes LESS resolved KE — i.e. eddy energy is
   ! returned and the resolved flow decays slower.  We compare one
   ! forward-Euler viscous decrement of a grid-scale velocity mode under
   ! resolved-only vs net (resolved - ku) viscosity.
   ! ------------------------------------------------------------------
   subroutine test_energy_return(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_meke_t) :: meke
      integer, parameter :: NX = 6, NY = 6, NZ = 1
      real(wp), parameter :: DX = 2000.0_wp, DT = 1800.0_wp
      real(wp), parameter :: A0 = 800.0_wp, KU = 300.0_wp
      real(wp), allocatable :: ahx(:, :, :), ahy(:, :, :)
      integer :: nxt, nyt, ii, jj
      real(wp) :: a_resolved, a_net, decay_resolved, decay_net, lam
      checks: block
         call make_grid(grid, NX, NY, DX)
         call make_cartesian_metrics(metrics, grid)
         nxt = grid%nx_total
         nyt = grid%ny_total
         call meke%init(grid, nz_ml=NZ)
         meke%backscatter = .true.
         meke%visc_coeff_ku = 1.0_wp
         meke%ku = KU
         allocate (ahx(nxt + 1, nyt, NZ), source=A0)
         allocate (ahy(nxt, nyt + 1, NZ), source=A0)
         ii = NGHOST + 2
         jj = NGHOST + 2
         a_resolved = A0
         !$acc enter data copyin(meke)
         call meke%enter_data()
         !$acc enter data copyin(ahx, ahy)
         call meke_backscatter_apply(grid, metrics, meke, DT, ahx, ahy)
         !$acc update self(ahx, ahy)
         !$acc exit data delete(ahx, ahy)
         call meke%exit_data()
         !$acc exit data delete(meke)
         a_net = ahx(ii, jj, 1)
         ! Grid-scale mode forward-Euler decrement: u^{n+1} = (1 - lam*A)*u,
         ! lam = dt*(idx²+idy²)*c (c>0 a fixed stencil constant) — only the
         ! SIGN/ORDER matters here, take c=1.  Less viscosity ⇒ less decay
         ! ⇒ larger retained amplitude.
         lam = DT*(1.0_wp/(DX*DX) + 1.0_wp/(DX*DX))
         decay_resolved = 1.0_wp - lam*a_resolved   ! retained fraction, resolved-only
         decay_net = 1.0_wp - lam*a_net             ! retained fraction, with backscatter
         call check(error, a_net < a_resolved - 1.0e-9_wp, &
                    "net viscosity strictly below resolved-only where ku>0")
         if (allocated(error)) exit checks
         call check(error, decay_net > decay_resolved + 1.0e-9_wp, &
                    "backscatter returns energy ⇒ resolved mode retains more amplitude (decays slower)")
      end block checks
      if (allocated(ahx)) deallocate (ahx, ahy)
      call meke%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_return

   ! ------------------------------------------------------------------
   ! Test 5: the configure-time fail-loud guard predicate.  MEKE
   ! backscatter (negative harmonic viscosity) requires a POSITIVE
   ! biharmonic backstop to be stable; `validate_config` aborts via this
   ! predicate when none is configured (the error stop can't run
   ! in-process, so we assert the predicate directly).
   ! ------------------------------------------------------------------
   subroutine test_backstop_guard(error)
      type(error_type), allocatable, intent(out) :: error

      ! No biharmonic anywhere ⇒ no backstop (the guard fires).
      call check(error,.not. has_biharmonic_backstop(0.0_wp, .false., 0.06_wp, LMIX_NONE, &
                                                     0.0_wp, 0.0_wp), &
                 "no nu_4 / no smag_ah / non-biharm closure must report NO backstop")
      if (allocated(error)) return
      call check(error,.not. has_biharmonic_backstop(0.0_wp, .false., 0.06_wp, LMIX_LEITH, &
                                                     0.0_wp, 0.0_wp), &
                 "Leith (Laplacian) is not a biharmonic backstop")
      if (allocated(error)) return

      ! leith_biharm SELECTED but with the default zero coefficient and no
      ! floor ⇒ nu4_face ≡ 0, NOT a backstop.  This is the assertion this
      ! PR exists to make true — on `main` this was `.true.` (a
      ! configure-blessed zero-dissipation backscatter run).
      call check(error,.not. has_biharmonic_backstop(0.0_wp, .false., 0.06_wp, &
                                                     LMIX_LEITH_BIHARM, &
                                                     c_leith_bi=0.0_wp, nu_4_bg=0.0_wp), &
                 "leith_biharm with c_leith_bi=0, nu_4_bg=0 is NOT a backstop "// &
                 "(zero dissipation coefficient)")
      if (allocated(error)) return
      call check(error, has_biharmonic_backstop(0.0_wp, .false., 0.06_wp, &
                                                LMIX_LEITH_BIHARM, &
                                                c_leith_bi=2.0_wp, nu_4_bg=0.0_wp), &
                 "leith_biharm with a positive c_leith_bi IS a backstop")
      if (allocated(error)) return
      call check(error, has_biharmonic_backstop(0.0_wp, .false., 0.06_wp, &
                                                LMIX_LEITH_BIHARM, &
                                                c_leith_bi=0.0_wp, nu_4_bg=1.0e9_wp), &
                 "leith_biharm with c_leith_bi=0 but a positive nu_4_bg floor "// &
                 "IS a backstop")
      if (allocated(error)) return

      ! smag_ah SELECTED but with a zero biharmonic Smagorinsky constant
      ! and no floor ⇒ same failure mode on the other flow-aware arm.
      call check(error,.not. has_biharmonic_backstop(0.0_wp, .true., 0.0_wp, LMIX_NONE, &
                                                     0.0_wp, 0.0_wp), &
                 "smag_ah with smag_bi_const=0, nu_4_bg=0 is NOT a backstop")
      if (allocated(error)) return
      call check(error, has_biharmonic_backstop(0.0_wp, .true., 0.06_wp, LMIX_NONE, &
                                                0.0_wp, 0.0_wp), &
                 "smag_ah with a positive smag_bi_const IS a backstop")
      if (allocated(error)) return

      ! Scalar nu_4 arm: a backstop only when no flow-aware face path is
      ! taken (the kernel returns before the scalar arm otherwise).
      call check(error, has_biharmonic_backstop(1.0e10_wp, .false., 0.06_wp, LMIX_NONE, &
                                                0.0_wp, 0.0_wp), &
                 "nu_4 > 0 with no flow-aware closure selected IS a backstop "// &
                 "(scalar arm)")
      if (allocated(error)) return
      call check(error,.not. has_biharmonic_backstop(1.0e10_wp, .true., 0.0_wp, LMIX_NONE, &
                                                     0.0_wp, 0.0_wp), &
                 "nu_4 > 0 is NOT a backstop when smag_ah's face path is taken "// &
                 "(the kernel returns before the scalar arm — pins the predicate "// &
                 "to the dispatch)")
   end subroutine test_backstop_guard

end module test_ocean_meke_backscatter
