!! Analytical tests for the law-of-the-wall bottom-boundary-layer diapycnal
!! mixing (`rdb_ocean_bbl_mixing`).
!!
!! Every expected value is derived independently of the kernel under test:
!!   * (a) the classical wall law `kappa u* z (1 - z/D)` (f = 0);
!!   * (b) the rotation-limited value `kappa u*^2 / |f|` at mid-depth;
!!   * (c) the energy ledger `Sum Kd dz max(N^2, Omega^2) = effic P_bbl`;
!!   * (d) the Joule decay `exp(-Sum h / lambda)` of the power offered;
!!   * (e) the unstratified fill to the cap;
!!   * (f) the face -> cell friction velocity and KE loss of a uniform flow;
!!   * (g) fail-loud configuration;
!!   * (h) the additive, kt-only merge;
!!   * (i) the column kernel's interface N^2 against the linear-EOS value.
!! Bottom-up indexing: k=1 bed, k=nz surface; `kd(K)` at the bottom
!! interface of layer K.
module test_ocean_bbl_mixing
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY, OMEGA_EARTH
   use rdb_config, only: config_t, read_config_from_string
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_bbl_mixing, only: ocean_bbl_mixing_t, bbl_column_sweep, &
                                   bbl_mixing_set_tke_impl, bbl_mixing_compute_impl, &
                                   bbl_mixing_merge_into_kt
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_bbl_mixing
   use rdb_error_ring, only: error_ring_clear, error_ring_get
   implicit none
   private

   public :: collect_ocean_bbl_mixing_tests

   real(wp), parameter :: VK = 0.41_wp
   real(wp), parameter :: HUGE_KD = 1.0e30_wp
   real(wp), parameter :: OM2 = OMEGA_EARTH**2

contains

   subroutine collect_ocean_bbl_mixing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wall_law_f0", test_wall_law_f0), &
                  new_unittest("rotation_limit", test_rotation_limit), &
                  new_unittest("energy_budget", test_energy_budget), &
                  new_unittest("joule_decay", test_joule_decay), &
                  new_unittest("unstratified_fill", test_unstratified_fill), &
                  new_unittest("face_to_cell", test_face_to_cell), &
                  new_unittest("fail_loud", test_fail_loud), &
                  new_unittest("merge_kt_only_additive", test_merge), &
                  new_unittest("column_kernel_n2", test_column_kernel_n2) &
                  ]
   end subroutine collect_ocean_bbl_mixing_tests

   ! ------------------------------------------------------------------
   ! (a) f = 0, abundant power: the classical wall law with surface closure
   ! ------------------------------------------------------------------
   subroutine test_wall_law_f0(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 12
      real(wp) :: h(NZ), n2(NZ), kd(NZ + 1), z, dtot, expect, ustar, worst
      integer :: k

      ustar = 0.01_wp
      do k = 1, NZ
         h(k) = 5.0_wp + 3.0_wp*real(k, wp)   ! non-uniform
      end do
      n2 = 0.0_wp                              ! floor Omega^2 only
      dtot = sum(h)
      call bbl_column_sweep(NZ, 1, h, n2, ustar, 0.0_wp, HUGE_KD, 200.0_wp, .true., &
                            OM2, VK, HUGE_KD, kd)
      worst = 0.0_wp
      z = 0.0_wp
      do k = 1, NZ - 1
         z = z + h(k)
         expect = VK*ustar*z*(1.0_wp - z/dtot)
         worst = max(worst, abs(kd(k + 1) - expect)/expect)
      end do
      call check(error, worst < 1.0e-12_wp, "f=0 wall law kappa u* z (1-z/D): rel err too big")
      if (allocated(error)) return
      call check(error, kd(1) == 0.0_wp .and. kd(NZ + 1) == 0.0_wp, "bed/surface interfaces must be 0")
   end subroutine test_wall_law_f0

   ! ------------------------------------------------------------------
   ! (b) strong rotation: Kd -> kappa u*^2 / |f| away from the ends
   ! ------------------------------------------------------------------
   subroutine test_rotation_limit(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 20
      real(wp) :: h(NZ), n2(NZ), kd(NZ + 1), kd0(NZ + 1), ustar, f, z, dtot, expect
      integer :: k, kmid

      ustar = 0.01_wp
      f = 1.0e-1_wp
      h = 5.0_wp
      n2 = 0.0_wp
      dtot = sum(h)
      ! tke = 1e300: the rotation-shortened decay (lambda = u*/|f| = 0.1 m)
      ! still leaves the mid-column ledger far above its Omega^2 cost.
      call bbl_column_sweep(NZ, 1, h, n2, ustar, f, 1.0e300_wp, 1.0e30_wp, .true., OM2, VK, HUGE_KD, kd)
      call bbl_column_sweep(NZ, 1, h, n2, ustar, 0.0_wp, 1.0e300_wp, 1.0e30_wp, .true., OM2, VK, HUGE_KD, kd0)
      kmid = NZ/2
      z = sum(h(1:kmid))
      ! Exact formula, evaluated independently.
      expect = VK*ustar**2*z*(dtot - z)/(ustar*dtot + f*z*(dtot - z))
      call check(error, abs(kd(kmid + 1) - expect) < 1.0e-12_wp*expect, &
                 "rotation: exact wall-law-with-rotation value")
      if (allocated(error)) return
      ! Limit: within 1 % of kappa u*^2/|f| at mid-depth (u* D/(f z (D-z)) = 0.4 %).
      call check(error, abs(kd(kmid + 1)*f/(VK*ustar**2) - 1.0_wp) < 1.0e-2_wp, &
                 "rotation: mid-depth Kd must approach kappa u*^2/|f|")
      if (allocated(error)) return
      do k = 2, NZ
         if (kd(k) >= kd0(k)) then
            call check(error, .false., "rotation must reduce Kd below the f=0 wall law at every interface")
            return
         end if
      end do
   end subroutine test_rotation_limit

   ! ------------------------------------------------------------------
   ! (c) energy budget: Sum Kd dz max(N^2, Omega^2) <= effic P, equality
   !     when the ledger is exhausted (no decay, f = 0)
   ! ------------------------------------------------------------------
   subroutine test_energy_budget(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 16
      real(wp) :: h(NZ), n2(NZ), kd(NZ + 1), tke, spent, dzk
      integer :: k

      do k = 1, NZ
         h(k) = 8.0_wp + real(mod(k, 3), wp)
         n2(k) = 1.0e-5_wp*(1.0_wp + real(k, wp))
      end do
      n2(NZ) = 0.0_wp
      ! Scarce power: ledger exhausted at some interface.
      tke = 2.0e-9_wp
      call bbl_column_sweep(NZ, 1, h, n2, 0.01_wp, 0.0_wp, tke, 1.0e30_wp, .true., OM2, VK, 0.1_wp, kd)
      spent = 0.0_wp
      do k = 1, NZ - 1
         dzk = 0.5_wp*(h(k) + h(k + 1))
         spent = spent + kd(k + 1)*dzk*max(n2(k), OM2)
      end do
      call check(error, abs(spent - tke) < 1.0e-12_wp*tke, &
                 "exhausted ledger: Sum Kd dz N^2 must equal effic*P_bbl")
      if (allocated(error)) return
      ! Abundant power: strictly under budget.
      tke = 1.0_wp
      call bbl_column_sweep(NZ, 1, h, n2, 0.01_wp, 0.0_wp, tke, 1.0e30_wp, .true., OM2, VK, 0.1_wp, kd)
      spent = 0.0_wp
      do k = 1, NZ - 1
         dzk = 0.5_wp*(h(k) + h(k + 1))
         spent = spent + kd(k + 1)*dzk*max(n2(k), OM2)
      end do
      call check(error, spent < tke .and. spent > 0.0_wp, "abundant power: spent must stay below the offer")
   end subroutine test_energy_budget

   ! ------------------------------------------------------------------
   ! (d) Joule decay of the power with height above the bed
   ! ------------------------------------------------------------------
   subroutine test_joule_decay(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: h(NZ), n2(NZ), kd(NZ + 1), tke, lam, expect, dzk, f, ustar
      integer :: k

      do k = 1, NZ
         h(k) = 10.0_wp + 5.0_wp*real(k, wp)
      end do
      ! Unstratified below (free: costs no power), a huge N^2 at the top
      ! interior interface so it uses everything offered there.
      n2 = 0.0_wp
      n2(NZ - 1) = 1.0e3_wp
      tke = 3.0e-6_wp
      ustar = 0.01_wp

      ! (i) f = 0: lambda = max_decay.
      lam = 120.0_wp
      call bbl_column_sweep(NZ, 1, h, n2, ustar, 0.0_wp, tke, lam, .false., OM2, VK, 0.1_wp, kd)
      dzk = 0.5_wp*(h(NZ - 1) + h(NZ))
      expect = tke*exp(-sum(h(1:NZ - 1))/lam)/(dzk*n2(NZ - 1))
      call check(error, abs(kd(NZ) - expect) < 1.0e-12_wp*expect, &
                 "top interface Kd must carry tke*exp(-Sum h/lambda)")
      if (allocated(error)) return
      call check(error, all(kd(2:NZ - 1) == 0.1_wp), "free (unstratified) interfaces below fill to the cap")
      if (allocated(error)) return

      ! (ii) rotation shortens lambda: 1/lambda = f/u* > 1/max_decay.
      f = 1.0e-4_wp                       ! f/u* = 1e-2 m^-1 -> lambda = 100 m < 1e30
      call bbl_column_sweep(NZ, 1, h, n2, ustar, f, tke, 1.0e30_wp, .false., OM2, VK, 0.1_wp, kd)
      expect = tke*exp(-sum(h(1:NZ - 1))*f/ustar)/(dzk*n2(NZ - 1))
      call check(error, abs(kd(NZ) - expect) < 1.0e-12_wp*expect, &
                 "rotation-shortened decay length u*/|f|")
      if (allocated(error)) return

      ! (iii) max_decay caps lambda from above: u*/|f| = 1e4 m > 120 m.
      f = 1.0e-6_wp
      call bbl_column_sweep(NZ, 1, h, n2, ustar, f, tke, lam, .false., OM2, VK, 0.1_wp, kd)
      expect = tke*exp(-sum(h(1:NZ - 1))/lam)/(dzk*n2(NZ - 1))
      call check(error, abs(kd(NZ) - expect) < 1.0e-12_wp*expect, "lambda = min(max_decay, u*/|f|)")
   end subroutine test_joule_decay

   ! ------------------------------------------------------------------
   ! (e) unstratified branch
   ! ------------------------------------------------------------------
   subroutine test_unstratified_fill(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 9
      real(wp) :: h(NZ), n2(NZ), kd(NZ + 1)

      h = 20.0_wp
      n2 = 0.0_wp
      call bbl_column_sweep(NZ, 1, h, n2, 0.01_wp, 0.0_wp, 1.0e-9_wp, 200.0_wp, .false., OM2, VK, 0.05_wp, kd)
      call check(error, all(kd(2:NZ) == 0.05_wp) .and. kd(1) == 0.0_wp .and. kd(NZ + 1) == 0.0_wp, &
                 "N^2=0 and no Omega floor: Kd = K_max at every interior interface while power remains")
      if (allocated(error)) return
      call bbl_column_sweep(NZ, 1, h, n2, 0.01_wp, 0.0_wp, 0.0_wp, 200.0_wp, .false., OM2, VK, 0.05_wp, kd)
      call check(error, all(kd == 0.0_wp), "no power left: unstratified Kd must be 0")
      if (allocated(error)) return
      ! A bed above k=1 (z_fixed filler): nothing below kbed.
      call bbl_column_sweep(NZ, 4, h, n2, 0.01_wp, 0.0_wp, 1.0e-9_wp, 200.0_wp, .false., OM2, VK, 0.05_wp, kd)
      call check(error, all(kd(1:4) == 0.0_wp) .and. all(kd(5:NZ) == 0.05_wp), &
                 "sweep starts at the live bed layer kbed")
   end subroutine test_unstratified_fill

   ! ------------------------------------------------------------------
   ! (f) face -> cell: uniform flow
   ! ------------------------------------------------------------------
   subroutine test_face_to_cell(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 4, NZ = 4
      real(wp), parameter :: CD = 3.0e-3_wp
      real(wp), parameter :: THICK = 3.0_wp
      real(wp) :: u(NX + 1, NY, NZ), v(NX, NY + 1, NZ), h(NX, NY, NZ), wet(NX, NY)
      real(wp) :: kvu(NX + 1, NY), kvv(NX, NY + 1), thu(NX + 1, NY), thv(NX, NY + 1)
      real(wp) :: fc(NX + 1, NY + 1), at(NX, NY), acu(NX + 1, NY), acv(NX, NY + 1)
      real(wp) :: ust(NX, NY), ke(NX, NY), fa(NX, NY)
      real(wp) :: sq, ubar, ustar_f, u1, u2

      sq = sqrt(CD)
      u1 = 0.3_wp
      u2 = -0.1_wp
      h = 2.0_wp                          ! layer 1 fully inside the BBL; layer 2 partly
      wet = 1.0_wp
      u = u1
      u(:, :, 2:NZ) = u2
      v = 0.0_wp
      fc = -1.0e-4_wp
      at = 100.0_wp
      acu = 100.0_wp
      acv = 100.0_wp
      ! BBL-mean along-face u over THICK = 3 m: (2*u1 + 1*u2)/3.
      ubar = (2.0_wp*u1 + 1.0_wp*u2)/3.0_wp
      ! Pick the face friction velocity as sqrt(c_d)*U (any value works).
      ustar_f = sq*0.25_wp
      thu = THICK
      thv = THICK
      kvu = sq*ustar_f*THICK
      kvv = sq*ustar_f*THICK

      !$acc enter data copyin(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv) create(ust, ke, fa)
      call bbl_mixing_set_tke_impl(NX, NY, NZ, CD, u, v, h, wet, kvu, kvv, thu, thv, fc, &
                                   at, acu, acv, ust, ke, fa)
      !$acc update self(ust, ke, fa)
      !$acc exit data delete(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv, ust, ke, fa)

      call check(error, abs(ust(2, 2) - sqrt(2.0_wp)*ustar_f) < 1.0e-12_wp*ustar_f, &
                 "u*_h = sqrt(1/2 * 4 u*_f^2) = sqrt(2) u*_f for four equal faces")
      if (allocated(error)) return
      ! v = 0 on the v-faces: only the two u-faces carry KE; 4-face sum, no 1/2.
      call check(error, abs(ke(2, 2) - 2.0_wp*sq*ustar_f*ubar**2) < 1.0e-12_wp*2.0_wp*sq*ustar_f*ubar**2, &
                 "P_bbl = 2 sqrt(c_d) u* ubar^2 (two u-faces, v = 0)")
      if (allocated(error)) return
      call check(error, abs(fa(2, 2) - 1.0e-4_wp) < 1.0e-16_wp, "|f|_h = mean of |f| at the corners")
      if (allocated(error)) return

      ! Dry cell is zero.
      wet(3, 3) = 0.0_wp
      !$acc enter data copyin(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv) create(ust, ke, fa)
      call bbl_mixing_set_tke_impl(NX, NY, NZ, CD, u, v, h, wet, kvu, kvv, thu, thv, fc, &
                                   at, acu, acv, ust, ke, fa)
      !$acc update self(ust, ke, fa)
      !$acc exit data delete(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv, ust, ke, fa)
      call check(error, ust(3, 3) == 0.0_wp .and. ke(3, 3) == 0.0_wp, "dry cell: u*_h = P = 0")
      if (allocated(error)) return

      ! A face with zero BBL thickness / viscosity (land face) contributes nothing.
      wet(3, 3) = 1.0_wp
      kvu(2, 2) = 0.0_wp
      !$acc enter data copyin(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv) create(ust, ke, fa)
      call bbl_mixing_set_tke_impl(NX, NY, NZ, CD, u, v, h, wet, kvu, kvv, thu, thv, fc, &
                                   at, acu, acv, ust, ke, fa)
      !$acc update self(ust, ke, fa)
      !$acc exit data delete(u, v, h, wet, kvu, kvv, thu, thv, fc, at, acu, acv, ust, ke, fa)
      ! Cell (2,2): west u-face (2,2) dead -> one u-face + two v-faces.
      call check(error, abs(ke(2, 2) - sq*ustar_f*ubar**2) < 1.0e-12_wp*sq*ustar_f*ubar**2, &
                 "a dead face drops out of the KE sum")
      if (allocated(error)) return
      call check(error, abs(ust(2, 2) - sqrt(1.5_wp)*ustar_f) < 1.0e-12_wp*ustar_f, &
                 "dead face: u*_h = sqrt(1/2 * 3 u*_f^2)")
   end subroutine test_face_to_cell

   ! ------------------------------------------------------------------
   ! (g) fail-loud configuration
   ! ------------------------------------------------------------------
   subroutine test_fail_loud(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ierr
      character(len=:), allocatable :: msg

      ! No per-face BBL (default slot).
      call configure_case(.false., 3.0e-3_wp, ierr, msg)
      call check(error, ierr /= 0 .and. index(msg, "per-face") > 0, &
                 "enable without the per-face BBL must fail loud; got: "//msg)
      if (allocated(error)) return
      ! Per-face BBL but c_d = 0.
      call configure_case(.true., 0.0_wp, ierr, msg)
      call check(error, ierr /= 0 .and. index(msg, "c_d") > 0, &
                 "enable with c_d = 0 must fail loud; got: "//msg)
      if (allocated(error)) return
      ! Valid configuration passes and allocates.
      call configure_case(.true., 3.0e-3_wp, ierr, msg)
      call check(error, ierr == 0, "valid BBL-mixing configuration must pass")
   end subroutine test_fail_loud

   subroutine configure_case(per_face, cd, ierr, msg)
      logical, intent(in) :: per_face
      real(wp), intent(in) :: cd
      integer, intent(out) :: ierr
      character(len=:), allocatable, intent(out) :: msg
      type(config_t) :: cfg
      type(ocean_state_t) :: st
      character(len=:), allocatable :: nml

      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 4, ny = 4, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_vmix_nml use_closure = .true. /"//new_line("a")// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//new_line("a")// &
            "&ocean_bbl_mixing_nml enable = .true., effic = 0.01 /"//new_line("a")
      call read_config_from_string(nml, cfg)
      allocate (st%multilayer%h_layer(4, 4, 3), source=1.0_wp)
      st%vdiff%bbl_glue = per_face
      st%vdiff%bbl_per_face = per_face
      st%vdiff%bbl_cd = cd
      call error_ring_clear()
      call configure_ocean_bbl_mixing(cfg, st, 1, ierr=ierr)
      msg = trim(error_ring_get(0))
      if (ierr == 0) then
         if (.not. (st%vdiff%bbl_mix%enable .and. allocated(st%vdiff%bbl_mix%kd_int) .and. &
                    abs(st%vdiff%bbl_mix%effic - 0.01_wp) < 1.0e-15_wp .and. &
                    abs(st%vdiff%bbl_mix%c_d - cd) < 1.0e-15_wp)) then
            ierr = -1
            msg = "configure did not latch the knobs / allocate"
         end if
      end if
   end subroutine configure_case

   ! ------------------------------------------------------------------
   ! (h) merge: kt only, additive, interior interfaces only
   ! ------------------------------------------------------------------
   subroutine test_merge(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 2, NZ = 5
      type(ocean_bbl_mixing_t) :: bm
      real(wp) :: kt(NX, NY, NZ + 1), kv(NX, NY, NZ + 1)
      integer :: k

      call bm%init(NX, NY, NZ)
      do k = 1, NZ + 1
         bm%kd_int(:, :, k) = 1.0e-3_wp*real(k, wp)
      end do
      kt = 2.0_wp
      kv = 5.0_wp
      !$acc enter data copyin(bm%kd_int, kt)
      call bbl_mixing_merge_into_kt(bm, NX, NY, NZ + 1, kt)
      !$acc update self(kt)
      !$acc exit data delete(bm%kd_int, kt)
      do k = 2, NZ
         if (abs(kt(1, 1, k) - (2.0_wp + 1.0e-3_wp*real(k, wp))) > 1.0e-15_wp) then
            call check(error, .false., "merge must add kd_int to kt at interior interfaces")
            return
         end if
      end do
      call check(error, kt(1, 1, 1) == 2.0_wp .and. kt(1, 1, NZ + 1) == 2.0_wp, &
                 "bed/surface kt end-caps untouched")
      if (allocated(error)) return
      ! The merge has no viscosity argument: kv is structurally untouched.
      call check(error, all(kv == 5.0_wp), "kv untouched (tracer diffusivity only)")
      call bm%destroy()
   end subroutine test_merge

   ! ------------------------------------------------------------------
   ! (i) column kernel: linear-EOS N^2 -> energy ledger equality
   ! ------------------------------------------------------------------
   subroutine test_column_kernel_n2(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 2, NY = 2, NZ = 10
      type(eos_t) :: eos
      real(wp) :: h(NX, NY, NZ), ht(NX, NY, NZ), hs(NX, NY, NZ), wet(NX, NY)
      integer :: kb(NX, NY)
      real(wp) :: ust(NX, NY), ke(NX, NY), fa(NX, NY), kd(NX, NY, NZ + 1)
      real(wp) :: tk(NZ), dz, dbuoy, n2k, spent, effic, tke, dzk
      real(wp) :: rho0_l
      integer :: k

      eos%variant = EOS_VARIANT_LINEAR
      eos%alpha_T = 0.2_wp
      eos%beta_S = 1035.0_wp*7.6e-4_wp
      eos%rho0 = 1035.0_wp
      dz = 20.0_wp
      do k = 1, NZ
         tk(k) = 10.0_wp + 0.05_wp*real(k*k, wp)    ! stable, k-dependent N^2
      end do
      do k = 1, NZ
         h(:, :, k) = dz
         ht(:, :, k) = tk(k)*dz
         hs(:, :, k) = 35.0_wp*dz
      end do
      wet = 1.0_wp
      kb = 1
      ust = 0.01_wp
      ke = 1.0e-4_wp
      fa = 0.0_wp
      effic = 0.5_wp
      tke = effic*ke(1, 1)
      rho0_l = eos%rho0

      !$acc enter data copyin(h, ht, hs, wet, kb, ust, ke, fa) create(kd)
      call bbl_mixing_compute_impl(NX, NY, NZ, h, ht, hs, wet, kb, eos, rho0_l, ust, ke, fa, &
                                   effic, 1.0e30_wp, .true., OM2, VK, 0.1_wp, kd)
      !$acc update self(kd)
      !$acc exit data delete(h, ht, hs, wet, kb, ust, ke, fa, kd)

      dbuoy = GRAVITY*eos%alpha_T/eos%rho0
      spent = 0.0_wp
      do k = 1, NZ - 1
         n2k = dbuoy*(tk(k + 1) - tk(k))/dz
         dzk = dz
         spent = spent + kd(1, 1, k + 1)*dzk*max(n2k, OM2)
      end do
      call check(error, abs(spent - tke) < 1.0e-9_wp*tke, &
                 "column kernel: Sum Kd dz max(N^2,Omega^2) must equal effic*P (linear-EOS N^2)")
      if (allocated(error)) return
      call check(error, kd(2, 2, 1) == 0.0_wp .and. kd(2, 2, NZ + 1) == 0.0_wp, "end-caps zero")
   end subroutine test_column_kernel_n2

end module test_ocean_bbl_mixing
