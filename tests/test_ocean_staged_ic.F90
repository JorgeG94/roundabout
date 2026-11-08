!! Staged initial state (`rdb_ocean_ic_inject`): the apply point in
!! `engine_setup` sits between the IC seed and the sponge reference snapshot.
!!
!! Cases:
!!   (A1) staged T reaches the sponge reference (ghost columns included, on a
!!        periodic axis they are the wrapped image) and staged u reaches `u_ref`.
!!   (A2) staged S seeds pseudo-salt bit for bit; staging pseudo-salt wins.
!!   (A4) staged h then T: the concentration read back is the staged one within
!!        a few ulp, whichever order the caller staged them in.
!!   (F)  h <= 0 on a wet cell, an unknown tracer name and a wrong shape all fail
!!        loud with their own status code.
!!   (B)  nothing staged: `ocean_ic_apply_staged` is a no-op (state unchanged).
!!
!! Host-only by construction: every check reads the state straight after
!! `engine_setup`, before `engine_enter_data`, so there are no device copies to
!! mirror and the test is identical on the GPU build (`mem:separate` never
!! enters).
module test_ocean_staged_ic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_teardown
   use rdb_ocean_ic_inject, only: staged_tracer_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IC_SEED, &
                               OCEAN_STATUS_ERR_NOT_FOUND, OCEAN_STATUS_ERR_BAD_SHAPE
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_staged_ic_tests

   integer, parameter :: NX = 10, NY = 8, NZ = 3, NG = 3
   real(wp), parameter :: DEPTH = 300.0_wp

   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_staged_ic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("staged_T_and_u_reach_the_sponge_reference", test_sponge_reference), &
                  new_unittest("staged_S_seeds_pseudo_salt_and_staging_it_wins", test_pseudo_salt), &
                  new_unittest("staged_h_then_T_gives_the_staged_concentration", test_h_then_t), &
                  new_unittest("staged_ic_fails_loud", test_fail_loud) &
                  ]
   end subroutine collect_ocean_staged_ic_tests

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function make_nml(pseudo) result(nml)
      logical, intent(in) :: pseudo
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny, sdepth

      write (snx, '(i0)') NX
      write (sny, '(i0)') NY
      write (sdepth, '(f0.1)') DEPTH
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
            "dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = "//trim(sdepth)//" /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'sponge', sponge_width = 2, sponge_strength = 1.0e-4, "// &
            "sponge_relax_tracers = .true. /"//NL// &
            "&ocean_sponge_nml enable = .true., damp_source = 'band', "// &
            "target_source = 'ic', relax_uv = .true., relax_tracers = .true. /"//NL
      if (pseudo) nml = nml//"&ocean_tracers_nml enable_pseudo_salt = .true. /"//NL
      nml = nml//"&output_nml output_to_file = .false. /"//NL
   end function make_nml

   pure function conc_field(base, di, dj, dk) result(f)
      real(wp), intent(in) :: base, di, dj, dk
      real(wp) :: f(NX, NY, NZ)
      integer :: i, j, k
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = base + di*real(i, wp) + dj*real(j, wp) + dk*real(k, wp)
            end do
         end do
      end do
   end function conc_field

   pure function thickness_field() result(h)
      real(wp) :: h(NX, NY, NZ)
      integer :: k
      ! Weights 1:2:3 bottom-up; the column sums to DEPTH.
      do k = 1, NZ
         h(:, :, k) = DEPTH*real(k, wp)/6.0_wp
      end do
   end function thickness_field

   subroutine setup(engine, cfg, nml, ierr)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      character(len=*), intent(in) :: nml
      integer, intent(out) :: ierr
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
   end subroutine setup

   subroutine stage_tracer(engine, name, conc)
      type(ocean_engine_t), intent(inout) :: engine
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: conc(NX, NY, NZ)
      type(staged_tracer_t), allocatable :: grown(:)
      integer :: n, i
      n = 0
      if (allocated(engine%staged_tracers)) n = size(engine%staged_tracers)
      allocate (grown(n + 1))
      do i = 1, n
         grown(i)%name = engine%staged_tracers(i)%name
         call move_alloc(engine%staged_tracers(i)%conc, grown(i)%conc)
      end do
      grown(n + 1)%name = name
      allocate (grown(n + 1)%conc(NX, NY, NZ), source=conc)
      call move_alloc(grown, engine%staged_tracers)
   end subroutine stage_tracer

   subroutine test_sponge_reference(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      real(wp) :: tc(NX, NY, NZ), u(NX + 1, NY, NZ), want, got, worst
      integer :: ierr, i, j, k, iw, jw, it

      call read_and_stage_sponge(engine, cfg, tc, u, ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "engine_setup with a staged IC failed")
      if (allocated(error)) return

      it = engine%state%multilayer%idx_temperature
      worst = 0.0_wp
      ! Whole array incl. ghosts: x is periodic (wrapped image), y zero-gradient.
      do k = 1, NZ
         do j = 1, NY + 2*NG
            do i = 1, NX + 2*NG
               iw = modulo(i - NG - 1, NX) + 1
               jw = min(max(j - NG, 1), NY)
               want = tc(iw, jw, k)
               got = engine%state%sponge%ref_tracer(i, j, k, it)
               worst = max(worst, abs(got - want)/abs(want))
            end do
         end do
      end do
      call check(error, worst < 1.0e-14_wp, "sponge ref_tracer(T) is not the staged field "// &
                 "(interior + wrapped/zero-gradient ghosts)")
      if (allocated(error)) return

      worst = 0.0_wp
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX + 1
               worst = max(worst, abs(engine%state%sponge%u_ref(NG + i, NG + j, k) - u(i, j, k)))
            end do
         end do
      end do
      call check(error, worst == 0.0_wp, "sponge u_ref is not the staged u")
      call engine_teardown(engine)
   end subroutine test_sponge_reference

   subroutine read_and_stage_sponge(engine, cfg, tc, u, ierr)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(out) :: tc(NX, NY, NZ), u(NX + 1, NY, NZ)
      integer, intent(out) :: ierr
      integer :: i, j, k

      tc = conc_field(8.0_wp, 0.1_wp, 0.03_wp, 1.7_wp)
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX + 1
               u(i, j, k) = 0.01_wp*real(i, wp) - 0.002_wp*real(j, wp) + 0.1_wp*real(k, wp)
            end do
         end do
      end do
      allocate (engine%staged_h(NX, NY, NZ), source=thickness_field())
      allocate (engine%staged_u(NX + 1, NY, NZ), source=u)
      call stage_tracer(engine, "temperature", tc)
      call setup(engine, cfg, make_nml(.false.), ierr)
   end subroutine read_and_stage_sponge

   subroutine test_pseudo_salt(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: e1, e2
      type(config_t) :: c1, c2
      real(wp) :: sc(NX, NY, NZ), pc(NX, NY, NZ), c
      integer :: ierr, ip, is, i, j, k

      sc = conc_field(33.0_wp, 0.05_wp, 0.01_wp, 0.4_wp)
      allocate (e1%staged_h(NX, NY, NZ), source=thickness_field())
      call stage_tracer(e1, "salinity", sc)
      call setup(e1, c1, make_nml(.true.), ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "setup (S staged) failed")
      if (allocated(error)) return
      ip = e1%state%multilayer%idx_pseudo_salt
      is = e1%state%multilayer%idx_salinity
      call check(error, ip > 0, "pseudo-salt was not registered")
      if (allocated(error)) return
      call check(error, all(e1%state%multilayer%tracers(ip)%c == &
                            e1%state%multilayer%tracers(is)%c), &
                 "pseudo-salt is not bit-for-bit the staged S")
      if (allocated(error)) return
      call engine_teardown(e1)

      pc = conc_field(20.0_wp, 0.0_wp, 0.0_wp, 1.0_wp)
      allocate (e2%staged_h(NX, NY, NZ), source=thickness_field())
      call stage_tracer(e2, "salinity", sc)
      call stage_tracer(e2, "pseudo_salt", pc)
      call setup(e2, c2, make_nml(.true.), ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "setup (S and pseudo-salt staged) failed")
      if (allocated(error)) return
      ip = e2%state%multilayer%idx_pseudo_salt
      is = e2%state%multilayer%idx_salinity
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               c = e2%state%multilayer%tracers(ip)%c(NG + i, NG + j, k)
               call check(error, abs(c - pc(i, j, k)) <= 8.0_wp*epsilon(1.0_wp)*abs(pc(i, j, k)), &
                          "staged pseudo-salt was overwritten by the S re-seed")
               if (allocated(error)) return
            end do
         end do
      end do
      call check(error, any(e2%state%multilayer%tracers(ip)%c /= &
                            e2%state%multilayer%tracers(is)%c), &
                 "pseudo-salt still equals S although it was staged separately")
      call engine_teardown(e2)
   end subroutine test_pseudo_salt

   subroutine test_h_then_t(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      real(wp) :: tc(NX, NY, NZ), h(NX, NY, NZ), c
      integer :: ierr, it, i, j, k

      tc = conc_field(8.0_wp, 0.1_wp, 0.03_wp, 1.7_wp)
      h = thickness_field()
      h(3, 4, :) = h(3, 4, :)*[1.2_wp, 0.9_wp, 0.95_wp]   ! non-uniform split, same order
      h(3, 4, :) = h(3, 4, :)*DEPTH/sum(h(3, 4, :))
      ! T is staged BEFORE h in the engine's lists on purpose: order is irrelevant.
      call stage_tracer(engine, "temperature", tc)
      allocate (engine%staged_h(NX, NY, NZ), source=h)
      call setup(engine, cfg, make_nml(.false.), ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "setup failed")
      if (allocated(error)) return
      it = engine%state%multilayer%idx_temperature
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               c = engine%state%multilayer%tracers(it)%c(NG + i, NG + j, k)
               call check(error, abs(c - tc(i, j, k)) <= 4.0_wp*epsilon(1.0_wp)*abs(tc(i, j, k)), &
                          "T concentration is not the staged one within a few ulp")
               if (allocated(error)) return
               call check(error, engine%state%multilayer%h_layer(NG + i, NG + j, k) == h(i, j, k), &
                          "h_layer is not the staged h")
               if (allocated(error)) return
            end do
         end do
      end do
      ! S was not staged: its concentration must survive the h change.
      it = engine%state%multilayer%idx_salinity
      c = engine%state%multilayer%tracers(it)%c(NG + 3, NG + 4, 2)
      call check(error, abs(c - 35.0_wp) <= 8.0_wp*epsilon(1.0_wp)*35.0_wp, &
                 "unstaged S lost its concentration across a staged h")
      call engine_teardown(engine)
   end subroutine test_h_then_t

   subroutine test_fail_loud(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: e
      type(config_t) :: cfg
      real(wp) :: h(NX, NY, NZ), c(NX, NY, NZ)
      integer :: ierr

      ! h <= 0 on a wet cell
      h = thickness_field()
      h(2, 2, 2) = 0.0_wp
      allocate (e%staged_h(NX, NY, NZ), source=h)
      call setup(e, cfg, make_nml(.false.), ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "h <= 0 must be ERR_IC_SEED")
      if (allocated(error)) return
      call engine_teardown(e)

      ! column sum wildly off (h in the wrong units)
      h = thickness_field()*0.01_wp
      allocate (e%staged_h(NX, NY, NZ), source=h)
      call setup(e, cfg, make_nml(.false.), ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "inconsistent column sum must be ERR_IC_SEED")
      if (allocated(error)) return
      call engine_teardown(e)

      ! unknown tracer
      c = 1.0_wp
      call stage_tracer(e, "no_such_tracer", c)
      call setup(e, cfg, make_nml(.false.), ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_NOT_FOUND, "unknown tracer must be ERR_NOT_FOUND")
      if (allocated(error)) return
      call engine_teardown(e)

      ! wrong shape (u missing a column)
      allocate (e%staged_u(NX, NY, NZ), source=0.0_wp)
      call setup(e, cfg, make_nml(.false.), ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_BAD_SHAPE, "wrong u shape must be ERR_BAD_SHAPE")
      call engine_teardown(e)
   end subroutine test_fail_loud

end module test_ocean_staged_ic
