!! Velocity-form BIHARMONIC viscosity under `&vcoord_nml
!! zfixed_closed_faces`: a closed face-layer (`metrics%open_u/open_v == 0`)
!! is a FREE-SLIP wall for both biharmonic paths — the scalar `nu_4`
!! (`hvisc_compute_biharmonic_impl`) and the flow-aware `nu4_face_*`
!! (`hvisc_compute_biharmonic_face_impl`, `smag_ah` / `leith_biharm`).
!!
!! The stencil rule (every face-to-face difference of BOTH chained
!! Laplacians weighted by `open(a)*open(b)`, tendency x own flag) is
!! derived and checked in the `biharm_zfixed` prototype; this suite is its
!! Fortran gate.  Two layers: layer 1 (bed) carries the staircase, layer 2
!! is all open and must behave exactly as the ungated kernel does.
!!
!! Covers, each for BOTH paths:
!!   * `free_slip_along_wall_{scalar,face}` — a uniform along-wall u beside
!!     a horizontal closed step, and a uniform along-wall v beside a
!!     vertical one, receive ZERO biharmonic tendency.  This is the case
!!     that FAILS on the pre-fix kernel (it read the zero stored at the
!!     closed face as a Dirichlet value: max|du| = O(nu4*U/dx^4)).
!!   * `zero_on_closed_{scalar,face}` — exactly zero tendency on every
!!     closed face-layer for a generic field on a staircase + pillar.
!!   * `energy_{scalar,face}` — with a constant nu4 on a region bounded by
!!     closed faces, `dE/dt = Sum A u du + Sum A v dv <= 0` for several
!!     fields, and E is non-increasing over a short forward-Euler march.
!!   * `open_layer_unchanged` — the all-open layer 2 is bit-identical to
!!     the ungated kernel's result (the gate is a multiply by exactly 1).
!!   * `refusals` — `stress_tensor` is still REFUSED under closed faces;
!!     `nu_4 > 0` now passes the hvisc gate (the refusal reached is the
!!     later `correction_bc_pgf` one).
module test_ocean_hvisc_biharm_zfixed
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, LMIX_NONE
   use rdb_ocean_horizontal_viscosity, only: &
      ocean_horizontal_viscosity_t, &
      ocean_horizontal_viscosity_compute_tendencies
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_closed_faces
   use rdb_ocean_status, only: OCEAN_STATUS_ERR_SETUP
   use rdb_ocean_vcoord, only: VCOORD_Z_FIXED
   use rdb_error_ring, only: error_ring_clear, error_ring_get
   implicit none
   private

   public :: collect_ocean_hvisc_biharm_zfixed_tests

   integer, parameter :: NXP = 16, NYP = 14, NZ = 2, NG = 2
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: NU4 = 1.0e9_wp
      !! m^4/s.  CFL ceiling at dx = 1 km, dt = 60 s, bound_coef = 0.8:
      !! 0.8*2/(60*(2e-6)^2) ~ 6.7e9 -- NU4 sits well below it, so the
      !! per-face clamp is inactive and nu4 is spatially constant.
   real(wp), parameter :: DT = 60.0_wp
   integer, parameter :: GEOM_HWALL = 1, GEOM_VWALL = 2, GEOM_STAIR = 3, GEOM_RING = 4

contains

   subroutine collect_ocean_hvisc_biharm_zfixed_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("free_slip_along_wall_scalar", test_free_slip_scalar), &
                  new_unittest("free_slip_along_wall_face", test_free_slip_face), &
                  new_unittest("zero_on_closed_scalar", test_zero_closed_scalar), &
                  new_unittest("zero_on_closed_face", test_zero_closed_face), &
                  new_unittest("energy_scalar", test_energy_scalar), &
                  new_unittest("energy_face", test_energy_face), &
                  new_unittest("open_layer_unchanged", test_open_layer_unchanged), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_ocean_hvisc_biharm_zfixed_tests

   ! ------------------------------------------------------------------
   ! Harness
   ! ------------------------------------------------------------------

   subroutine filler_of(geom, nx, ny, fil)
      !! T-cell filler pattern of the BED layer (layer 1) for `geom`
      !! (total indices, 1-based).
      integer, intent(in) :: geom, nx, ny
      logical, intent(out) :: fil(nx, ny)
      integer :: i
      fil = .false.
      select case (geom)
      case (GEOM_HWALL)
         ! horizontal step: filler for j >= ny-4, a lower tread for i <= 6
         fil(:, ny - 4:) = .true.
         fil(1:6, ny - 6:) = .true.
      case (GEOM_VWALL)
         fil(nx - 4:, :) = .true.
         fil(nx - 6:, ny - 5:) = .true.
      case (GEOM_STAIR, GEOM_RING)
         do i = 1, nx
            fil(i, ny - 3 - (4*(i - 1))/nx:) = .true.
         end do
         fil(6:7, 5:6) = .true.     ! an isolated pillar (seamount top)
         if (geom == GEOM_RING) then
            ! the open region is bounded by closed faces only, so the
            ! array-edge mirror rows (mask-independent) are never reached
            fil(1:2, :) = .true.
            fil(nx - 1:, :) = .true.
            fil(:, 1:2) = .true.
         end if
      end select
   end subroutine filler_of

   subroutine set_masks(metrics, geom, nx, ny)
      !! Layer 1: a face is closed iff the layer is a filler on EITHER
      !! side (array-edge faces stay open, as `ocean_vcoord_closed_face_masks`
      !! builds them).  Layer 2: all open.  Host edit + device push (the
      !! slot was sized via `nz_closed` BEFORE the map).
      type(ocean_metrics_t), intent(inout) :: metrics
      integer, intent(in) :: geom, nx, ny
      logical :: fil(nx, ny)
      integer :: i, j
      call filler_of(geom, nx, ny, fil)
      metrics%open_u = 1.0_wp
      metrics%open_v = 1.0_wp
      do j = 1, ny
         do i = 2, nx
            if (fil(i - 1, j) .or. fil(i, j)) metrics%open_u(i, j, 1) = 0.0_wp
         end do
      end do
      do j = 2, ny
         do i = 1, nx
            if (fil(i, j - 1) .or. fil(i, j)) metrics%open_v(i, j, 1) = 0.0_wp
         end do
      end do
      !$acc update device(metrics%open_u, metrics%open_v)
      metrics%use_closed_faces = .true.
   end subroutine set_masks

   subroutine run_biharm(geom, use_face, closed, u_in, v_in, du, dv, open_u, open_v)
      !! One biharmonic-only tendency (`nu_h = 0`) on a fresh state.
      !! `closed = .false.` runs the UNGATED kernel on the same masked
      !! velocities (what the pre-fix code did with closed faces on).
      integer, intent(in) :: geom
      logical, intent(in) :: use_face, closed
      real(wp), intent(in) :: u_in(:, :, :), v_in(:, :, :)
      real(wp), allocatable, intent(out) :: du(:, :, :), dv(:, :, :)
      real(wp), allocatable, intent(out) :: open_u(:, :, :), open_v(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_lateral_mix_t) :: lmix
      integer :: nx, ny

      call grid%init(NXP, NYP, NG, DX, DX)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 100.0_wp
      ms%u_face_x_layer = u_in
      ms%v_face_y_layer = v_in
      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = 0.0_wp
      hv%nu_4 = merge(0.0_wp, NU4, use_face)
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_NONE
      lmix%smag_ah_active = use_face
      lmix%nu4_face_x = NU4
      lmix%nu4_face_y = NU4

      call make_cartesian_metrics(metrics, grid, nz_closed=NZ)
      call set_masks(metrics, geom, nx, ny)
      allocate (open_u, source=metrics%open_u)
      allocate (open_v, source=metrics%open_v)
      if (.not. closed) metrics%use_closed_faces = .false.
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()

      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, &
                                                         lateral_mix=lmix, dt=DT)
      !$acc update self(hv%du_visc%data, hv%dv_visc%data)
      allocate (du, source=hv%du_visc%data)
      allocate (dv, source=hv%dv_visc%data)

      call hv%exit_data()
      !$acc exit data delete(hv)
      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call hv%destroy()
      call lmix%destroy()
      call ms%destroy()
   end subroutine run_biharm

   subroutine masks_only(geom, open_u, open_v)
      !! The layer masks for `geom`, without running anything.
      integer, intent(in) :: geom
      real(wp), allocatable, intent(out) :: open_u(:, :, :), open_v(:, :, :)
      logical, allocatable :: fil(:, :)
      integer :: nx, ny, i, j
      nx = NXP + 2*NG
      ny = NYP + 2*NG
      allocate (fil(nx, ny))
      call filler_of(geom, nx, ny, fil)
      allocate (open_u(nx + 1, ny, NZ), open_v(nx, ny + 1, NZ), source=1.0_wp)
      do j = 1, ny
         do i = 2, nx
            if (fil(i - 1, j) .or. fil(i, j)) open_u(i, j, 1) = 0.0_wp
         end do
      end do
      do j = 2, ny
         do i = 1, nx
            if (fil(i, j - 1) .or. fil(i, j)) open_v(i, j, 1) = 0.0_wp
         end do
      end do
   end subroutine masks_only

   subroutine generic_field(seed, open_u, open_v, u, v)
      !! A smooth-plus-grid-scale field, zeroed on closed faces (what
      !! `mask_layer_velocities` leaves).
      integer, intent(in) :: seed
      real(wp), intent(in) :: open_u(:, :, :), open_v(:, :, :)
      real(wp), allocatable, intent(out) :: u(:, :, :), v(:, :, :)
      integer :: i, j, k
      real(wp) :: s
      s = real(seed, wp)
      allocate (u, mold=open_u)
      allocate (v, mold=open_v)
      do k = 1, size(u, 3)
         do j = 1, size(u, 2)
            do i = 1, size(u, 1)
               u(i, j, k) = (sin(0.37_wp*s*real(i, wp) + 0.11_wp*real(j*k, wp)) + &
                             0.3_wp*cos(1.7_wp*real(i + 2*j, wp) + s))*open_u(i, j, k)
            end do
         end do
         do j = 1, size(v, 2)
            do i = 1, size(v, 1)
               v(i, j, k) = (cos(0.23_wp*real(i, wp) - 0.41_wp*s*real(j, wp)) - &
                             0.4_wp*sin(2.3_wp*real(3*i + j, wp) - s))*open_v(i, j, k)
            end do
         end do
      end do
   end subroutine generic_field

   ! ------------------------------------------------------------------
   ! Free slip
   ! ------------------------------------------------------------------

   subroutine test_free_slip_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      call check_free_slip(error, .false.)
   end subroutine test_free_slip_scalar

   subroutine test_free_slip_face(error)
      type(error_type), allocatable, intent(out) :: error
      call check_free_slip(error, .true.)
   end subroutine test_free_slip_face

   subroutine check_free_slip(error, use_face)
      !! Uniform along-wall flow beside a closed step must receive ZERO
      !! biharmonic tendency (a free-slip wall does not drag it).  The
      !! ungated kernel on the same field is required to drag it — that
      !! is the pre-fix defect this suite exists to pin.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_face
      real(wp), parameter :: U0 = 0.3_wp
      real(wp), allocatable :: ou(:, :, :), ov(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: du(:, :, :), dv(:, :, :), du0(:, :, :), dv0(:, :, :)
      real(wp), allocatable :: o2u(:, :, :), o2v(:, :, :)
      real(wp) :: gated, ungated
      character(len=8) :: tag

      tag = merge("face  ", "scalar", use_face)
      ! (1) u along a horizontal step
      call masks_only(GEOM_HWALL, ou, ov)
      u = U0*ou
      allocate (v, mold=ov)
      v = 0.0_wp
      call run_biharm(GEOM_HWALL, use_face, .true., u, v, du, dv, o2u, o2v)
      call run_biharm(GEOM_HWALL, use_face, .false., u, v, du0, dv0, o2u, o2v)
      gated = maxval(abs(du))
      ungated = maxval(abs(du0*ou))
      call check(error, gated <= 1.0e-12_wp*NU4*U0/DX**4, &
                 trim(tag)//": uniform u along a closed step must be untouched "// &
                 "(free slip); max|du| = "//r2s(gated))
      if (allocated(error)) return
      call check(error, ungated > 1.0e-3_wp*NU4*U0/DX**4, &
                 trim(tag)//": the UNGATED kernel must drag the along-wall flow "// &
                 "(else this test cannot see the defect); max|du| = "//r2s(ungated))
      if (allocated(error)) return
      deallocate (u, v, du, dv, du0, dv0)

      ! (2) v along a vertical step
      call masks_only(GEOM_VWALL, ou, ov)
      allocate (u, mold=ou)
      u = 0.0_wp
      v = -U0*ov
      call run_biharm(GEOM_VWALL, use_face, .true., u, v, du, dv, o2u, o2v)
      call run_biharm(GEOM_VWALL, use_face, .false., u, v, du0, dv0, o2u, o2v)
      gated = maxval(abs(dv))
      ungated = maxval(abs(dv0*ov))
      call check(error, gated <= 1.0e-12_wp*NU4*U0/DX**4, &
                 trim(tag)//": uniform v along a closed step must be untouched "// &
                 "(free slip); max|dv| = "//r2s(gated))
      if (allocated(error)) return
      call check(error, ungated > 1.0e-3_wp*NU4*U0/DX**4, &
                 trim(tag)//": the UNGATED kernel must drag the along-wall v; "// &
                 "max|dv| = "//r2s(ungated))
   end subroutine check_free_slip

   ! ------------------------------------------------------------------
   ! Zero on closed faces
   ! ------------------------------------------------------------------

   subroutine test_zero_closed_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      call check_zero_closed(error, .false.)
   end subroutine test_zero_closed_scalar

   subroutine test_zero_closed_face(error)
      type(error_type), allocatable, intent(out) :: error
      call check_zero_closed(error, .true.)
   end subroutine test_zero_closed_face

   subroutine check_zero_closed(error, use_face)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_face
      real(wp), allocatable :: ou(:, :, :), ov(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: du(:, :, :), dv(:, :, :), o2u(:, :, :), o2v(:, :, :)
      real(wp) :: on_closed, on_open
      call masks_only(GEOM_STAIR, ou, ov)
      call generic_field(1, ou, ov, u, v)
      call run_biharm(GEOM_STAIR, use_face, .true., u, v, du, dv, o2u, o2v)
      on_closed = max(maxval(abs(du*(1.0_wp - ou))), maxval(abs(dv*(1.0_wp - ov))))
      on_open = max(maxval(abs(du*ou)), maxval(abs(dv*ov)))
      call check(error, on_closed == 0.0_wp, &
                 merge("face  ", "scalar", use_face)// &
                 ": closed face-layers must get EXACTLY zero biharmonic tendency; "// &
                 "max = "//r2s(on_closed))
      if (allocated(error)) return
      call check(error, on_open > 0.0_wp, "the generic field must be damped somewhere")
   end subroutine check_zero_closed

   ! ------------------------------------------------------------------
   ! Energy
   ! ------------------------------------------------------------------

   subroutine test_energy_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      call check_energy(error, .false.)
   end subroutine test_energy_scalar

   subroutine test_energy_face(error)
      type(error_type), allocatable, intent(out) :: error
      call check_energy(error, .true.)
   end subroutine test_energy_face

   subroutine check_energy(error, use_face)
      !! Constant nu4, open region bounded by closed faces only: the gated
      !! operator is `-nu4 A^-1 L A^-1 L` with `L` symmetric, so
      !! `dE/dt = Sum u du + Sum v dv` (uniform face area) is <= 0 for
      !! every field, and a forward-Euler march at dt well inside the CFL
      !! bound never raises E.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_face
      integer, parameter :: NSTEP = 8
      real(wp), allocatable :: ou(:, :, :), ov(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: du(:, :, :), dv(:, :, :), o2u(:, :, :), o2v(:, :, :)
      real(wp) :: dedt, e_old, e_new
      integer :: seed, n
      character(len=6) :: tag

      tag = merge("face  ", "scalar", use_face)
      call masks_only(GEOM_RING, ou, ov)
      do seed = 1, 3
         call generic_field(seed, ou, ov, u, v)
         call run_biharm(GEOM_RING, use_face, .true., u, v, du, dv, o2u, o2v)
         dedt = sum(u*du) + sum(v*dv)
         call check(error, dedt <= 0.0_wp, &
                    trim(tag)//": the gated biharmonic must not create energy; "// &
                    "dE/dt = "//r2s(dedt))
         if (allocated(error)) return
         call check(error, dedt < 0.0_wp, trim(tag)//": the field must be damped")
         if (allocated(error)) return
         deallocate (u, v, du, dv)
      end do

      call generic_field(4, ou, ov, u, v)
      e_old = 0.5_wp*(sum(u*u) + sum(v*v))
      do n = 1, NSTEP
         call run_biharm(GEOM_RING, use_face, .true., u, v, du, dv, o2u, o2v)
         u = (u + DT*du)*ou
         v = (v + DT*dv)*ov
         e_new = 0.5_wp*(sum(u*u) + sum(v*v))
         call check(error, e_new <= e_old, &
                    trim(tag)//": E rose under the gated biharmonic at step "// &
                    i2s(n)//": "//r2s(e_old)//" -> "//r2s(e_new))
         if (allocated(error)) return
         e_old = e_new
         deallocate (du, dv)
      end do
   end subroutine check_energy

   ! ------------------------------------------------------------------
   ! The all-open layer is untouched by the gate
   ! ------------------------------------------------------------------

   subroutine test_open_layer_unchanged(error)
      !! Layer 2 is all open: the gated result there must be bit-identical
      !! to the ungated kernel's (each gate multiplies by exactly 1.0).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ou(:, :, :), ov(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: du(:, :, :), dv(:, :, :), du0(:, :, :), dv0(:, :, :)
      real(wp), allocatable :: o2u(:, :, :), o2v(:, :, :)
      logical :: face
      integer :: p
      call masks_only(GEOM_STAIR, ou, ov)
      call generic_field(2, ou, ov, u, v)
      do p = 1, 2
         face = p == 2
         call run_biharm(GEOM_STAIR, face, .true., u, v, du, dv, o2u, o2v)
         call run_biharm(GEOM_STAIR, face, .false., u, v, du0, dv0, o2u, o2v)
         call check(error, all(du(:, :, 2) == du0(:, :, 2)) .and. &
                    all(dv(:, :, 2) == dv0(:, :, 2)), &
                    "the all-open layer must be bit-identical to the ungated kernel")
         if (allocated(error)) return
         deallocate (du, dv, du0, dv0)
      end do
   end subroutine test_open_layer_unchanged

   ! ------------------------------------------------------------------
   ! Configure gate
   ! ------------------------------------------------------------------

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: msg
      integer :: ierr

      call configure_case("&ocean_hvisc_nml stress_tensor = .true. /", ierr, msg)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP .and. &
                 index(msg, "stress_tensor") > 0, &
                 "stress_tensor must still be REFUSED under zfixed_closed_faces; got: "//msg)
      if (allocated(error)) return

      ! nu_4 > 0 together with a LATER refusal: the hvisc gate must let
      ! nu_4 through, so the refusal reached is correction_bc_pgf's.
      call configure_case("&ocean_hvisc_nml nu_4 = 1.0e9 /"//new_line("a")// &
                          "&ocean_bt_nml correction_bc_pgf = .true. /", ierr, msg)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP .and. &
                 index(msg, "correction_bc_pgf") > 0 .and. index(msg, "nu_4") == 0, &
                 "nu_4 must now pass the closed-face hvisc gate; got: "//msg)
   end subroutine test_refusals

   subroutine configure_case(extra, ierr, msg)
      character(len=*), intent(in) :: extra
      integer, intent(out) :: ierr
      character(len=:), allocatable, intent(out) :: msg
      type(config_t) :: cfg
      type(ocean_state_t) :: st
      type(hgrid_t) :: grid
      character(len=:), allocatable :: nml
      st%multilayer%is_init = .true.
      st%vcoord%coord_type = VCOORD_Z_FIXED
      st%vcoord%z_fixed_h_ref = 400.0_wp
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            '&vcoord_nml vcoord_type = "z_fixed", zfixed_closed_faces = .true. /'// &
            new_line("a")//extra//new_line("a")
      call read_config_from_string(nml, cfg)
      call error_ring_clear()
      call configure_ocean_closed_faces(cfg, st, grid, 1, ierr=ierr)
      msg = trim(error_ring_get(0))
   end subroutine configure_case

   pure function r2s(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(es12.4)') x
      s = trim(adjustl(buf))
   end function r2s

   pure function i2s(n) result(s)
      integer, intent(in) :: n
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write (buf, '(i0)') n
      s = trim(buf)
   end function i2s

end module test_ocean_hvisc_biharm_zfixed
