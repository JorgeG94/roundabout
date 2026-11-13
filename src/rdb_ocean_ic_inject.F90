!! Pre-create initial-STATE injection (h / u / v / tracer concentrations).
module rdb_ocean_ic_inject
   !! A caller building an initial condition in memory (Python's
   !! `Model.stage_h` / `stage_u` / `stage_v` / `stage_tracer`, or a Fortran
   !! caller filling `ocean_engine_t%staged_*` directly) hands the library
   !! interior-sized arrays in the pending-create window.  `engine_setup`
   !! calls `ocean_ic_apply_staged` ONCE, right after the built-in IC seed and
   !! BEFORE the sponge reference snapshot, so every setup-time snapshot and
   !! derivation (sponge `ref_tracer`/`u_ref`/`v_ref`, land-face masking,
   !! tidal-mixing fold, wet/dry mask, the init-time seam exchange and
   !! `enter_data`) sees the user's fields instead of the built-in ones.
   !!
   !! Everything here is HOST-side plain `do` loops: the call precedes
   !! `ocean_state_enter_data`, so no device copies exist yet and none are
   !! made (`mem:separate` safe by construction).
   !!
   !! ## Order (load-bearing)
   !! 1. `h` first — validated (finite, `h > 0` on wet cells; column sum
   !!    within `IC_H_COLUMN_RTOL` of the water column).  Tracers NOT staged
   !!    keep their CONCENTRATION (stored directly; untouched).
   !! 2. Staged tracers: `c = conc` (the value `rdb_ocean_set_tracer`
   !!    forms), against the NEW `h`, so staging order is irrelevant.
   !! 3. `u` / `v` on the interior faces.
   !! 4. Ghosts: periodic wrap on a periodic axis, zero-gradient at every
   !!    other edge (never `h = 0` ghosts: the EOS would fall back to rho0 and
   !!    manufacture a density jump at wall-adjacent faces).
   !! 5. Pseudo-salt re-seeded from S unless pseudo-salt itself was staged.
   !! 6. EOS, gated by `thermo_active` exactly as the step gates it, so
   !!    `rho_layer` is current for the sponge snapshot while a prescribed
   !!    `layer_rho_init` is left alone when thermodynamics are off.
   !!
   !! Nothing staged => the first statement returns: bit-identical to the
   !! ordinary create.
   !!
   !! Single-rank only: staged arrays are the whole (global) interior.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_error_ring, only: fail
   use pic_strings, only: to_string
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pseudo_salt, only: ocean_pseudo_salt_seed
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IC_SEED, &
                               OCEAN_STATUS_ERR_BAD_SHAPE, OCEAN_STATUS_ERR_NOT_FOUND
   implicit none
   private

   public :: staged_tracer_t
   public :: ocean_ic_apply_staged
   public :: ocean_ic_fill_state_ghosts
   public :: IC_H_COLUMN_RTOL

   real(wp), parameter :: IC_H_COLUMN_RTOL = 0.25_wp
      !! Largest accepted `|sum_k h - D| / D` on a wet column, `D` the water
      !! column (`b`, or `b - z_draft` under an ice shelf).  A real free
      !! surface moves a column by O(1 m) against kilometres of water, so 25 %
      !! admits every physical sea-surface height and every thin-shelf
      !! column, while catching what staging mistakes actually look like: `h`
      !! in the wrong units (cm, fractions of the depth), the wrong
      !! `nz`-ordering of a layer-thickness profile, or an array that never
      !! summed to the bathymetry that was staged with it.

   type :: staged_tracer_t
      !! One staged tracer: registry name plus the CONCENTRATION field.
      character(len=32) :: name = ""
      real(wp), allocatable :: conc(:, :, :)
         !! `(nx_phys, ny_phys, nz)`, k = 1 the bed.
   end type staged_tracer_t

contains

   subroutine ocean_ic_apply_staged(state, grid, staged_h, staged_u, staged_v, &
                                    staged_tracers, per_x, per_y, thermo_active, ierr)
      !! Apply the staged initial state to the just-seeded `state`.  See the
      !! module docstring for the ordering and the failure modes.  Every
      !! failure sets `ierr` (`OCEAN_STATUS_ERR_BAD_SHAPE`, `_NOT_FOUND`, or
      !! `_IC_SEED`) and leaves `state` partially written: the caller aborts
      !! the create.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in), optional :: staged_h(:, :, :)
         !! `(nx_phys, ny_phys, nz)` layer thickness [m], k = 1 the bed.
      real(wp), intent(in), optional :: staged_u(:, :, :)
         !! `(nx_phys+1, ny_phys, nz)` x-face velocity [m/s].
      real(wp), intent(in), optional :: staged_v(:, :, :)
         !! `(nx_phys, ny_phys+1, nz)` y-face velocity [m/s].
      type(staged_tracer_t), intent(in), optional :: staged_tracers(:)
      logical, intent(in) :: per_x
         !! West/east edges are periodic.
      logical, intent(in) :: per_y
         !! South/north edges are periodic.
      logical, intent(in) :: thermo_active
         !! `&ocean_thermo_nml enable_thermodynamics`.
      integer, intent(out), optional :: ierr

      integer :: nx, ny, nz, ng, nxt, nyt, it, is, ntr, ierr_local
      logical :: any_h, any_u, any_v, any_tr, pseudo_staged
      integer, allocatable :: tr_idx(:)
      logical, allocatable :: is_staged(:)

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      any_h = present(staged_h)
      any_u = present(staged_u)
      any_v = present(staged_v)
      any_tr = .false.
      if (present(staged_tracers)) any_tr = size(staged_tracers) > 0
      if (.not. (any_h .or. any_u .or. any_v .or. any_tr)) return

      nx = grid%nx_phys
      ny = grid%ny_phys
      nz = state%multilayer%nz_ml
      ng = grid%nghost
      nxt = grid%nx_total
      nyt = grid%ny_total

      ! ---- shapes ----
      if (any_h) then
         if (.not. shape_ok(staged_h, nx, ny, nz)) then
            call shape_fail("stage_h", shape(staged_h), [nx, ny, nz], ierr)
            return
         end if
      end if
      if (any_u) then
         if (.not. shape_ok(staged_u, nx + 1, ny, nz)) then
            call shape_fail("stage_u", shape(staged_u), [nx + 1, ny, nz], ierr)
            return
         end if
      end if
      if (any_v) then
         if (.not. shape_ok(staged_v, nx, ny + 1, nz)) then
            call shape_fail("stage_v", shape(staged_v), [nx, ny + 1, nz], ierr)
            return
         end if
      end if

      ! ---- tracer names against the registry ----
      ntr = 0
      if (allocated(state%multilayer%tracers)) ntr = size(state%multilayer%tracers)
      allocate (is_staged(max(ntr, 1)), source=.false.)
      pseudo_staged = .false.
      if (any_tr) then
         allocate (tr_idx(size(staged_tracers)))
         do is = 1, size(staged_tracers)
            tr_idx(is) = 0
            do it = 1, ntr
               if (trim(state%multilayer%tracers(it)%name) == trim(staged_tracers(is)%name)) then
                  tr_idx(is) = it
                  exit
               end if
            end do
            if (tr_idx(is) == 0) then
               call fail("stage_tracer: no tracer named '"//trim(staged_tracers(is)%name)// &
                         "' in the registry (S, T, and any tracer enabled by the namelist "// &
                         "or registered in the pending window are stageable)", &
                         ierr_local, OCEAN_STATUS_ERR_NOT_FOUND)
               if (present(ierr)) ierr = ierr_local
               return
            end if
            if (.not. allocated(staged_tracers(is)%conc)) then
               call fail("stage_tracer '"//trim(staged_tracers(is)%name)//"': no data", &
                         ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
               if (present(ierr)) ierr = ierr_local
               return
            end if
            if (.not. shape_ok(staged_tracers(is)%conc, nx, ny, nz)) then
               call shape_fail("stage_tracer '"//trim(staged_tracers(is)%name)//"'", &
                               shape(staged_tracers(is)%conc), [nx, ny, nz], ierr)
               return
            end if
            if (any(.not. ieee_is_finite(staged_tracers(is)%conc))) then
               call fail("stage_tracer '"//trim(staged_tracers(is)%name)// &
                         "': non-finite concentration", ierr_local, OCEAN_STATUS_ERR_IC_SEED)
               if (present(ierr)) ierr = ierr_local
               return
            end if
            is_staged(tr_idx(is)) = .true.
            if (tr_idx(is) == state%multilayer%idx_pseudo_salt) pseudo_staged = .true.
         end do
      end if
      if (any_u) then
         if (any(.not. ieee_is_finite(staged_u))) then
            call fail("stage_u: non-finite velocity", ierr_local, OCEAN_STATUS_ERR_IC_SEED)
            if (present(ierr)) ierr = ierr_local
            return
         end if
      end if
      if (any_v) then
         if (any(.not. ieee_is_finite(staged_v))) then
            call fail("stage_v: non-finite velocity", ierr_local, OCEAN_STATUS_ERR_IC_SEED)
            if (present(ierr)) ierr = ierr_local
            return
         end if
      end if

      ! ---- 1. h (validated), unstaged tracers keep their concentration ----
      if (any_h) then
         call apply_h(state, staged_h, is_staged, nx, ny, nz, ng, ierr_local)
         if (ierr_local /= OCEAN_STATUS_OK) then
            if (present(ierr)) ierr = ierr_local
            return
         end if
      end if

      ! ---- 2. tracers: c = conc ----
      if (any_tr) then
         do is = 1, size(staged_tracers)
            call set_conc(state%multilayer%tracers(tr_idx(is))%c, &
                          state%multilayer%wet_mask, staged_tracers(is)%conc, nx, ny, nz, ng, &
                          nxt, nyt)
         end do
      end if

      ! ---- 3. velocities on interior faces ----
      if (any_u) then
         state%multilayer%u_face_x_layer(ng + 1:ng + nx + 1, ng + 1:ng + ny, 1:nz) = staged_u
      end if
      if (any_v) then
         state%multilayer%v_face_y_layer(ng + 1:ng + nx, ng + 1:ng + ny + 1, 1:nz) = staged_v
      end if

      ! ---- 4. ghosts ----
      call ocean_ic_fill_state_ghosts(state, grid, per_x, per_y, fill_h=any_h, &
                                      fill_u=any_u, fill_v=any_v, tracer_mask=is_staged, &
                                      all_tracers=any_h)

      ! ---- 5. pseudo-salt follows S unless the caller staged it ----
      if (.not. pseudo_staged) call ocean_pseudo_salt_seed(state%multilayer)

      ! ---- 6. rho_layer from the final T/S/h ----
      call ocean_eos_compute(state%eos, state%multilayer, active=thermo_active)
   end subroutine ocean_ic_apply_staged

   subroutine ocean_ic_fill_state_ghosts(state, grid, per_x, per_y, fill_h, fill_u, fill_v, &
                                         tracer_mask, all_tracers)
      !! Fill the ghost bands of the host multilayer state from its interior:
      !! periodic image on a periodic axis, zero-gradient elsewhere (never
      !! `h = 0` ghosts, see the module docstring).  HOST plain loops, so it
      !! is correct before `enter_data` (the staged path) AND on the live
      !! handle (`ocean_rebaseline_apply` pushes the result to the device
      !! afterwards).  With every optional absent it fills `h`, `u`, `v` and
      !! every tracer.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: per_x
      logical, intent(in) :: per_y
      logical, intent(in), optional :: fill_h
      logical, intent(in), optional :: fill_u
      logical, intent(in), optional :: fill_v
      logical, intent(in), optional :: tracer_mask(:)
         !! Registry-indexed: tracers to fill (absent: all).
      logical, intent(in), optional :: all_tracers
         !! `.true.`: fill every tracer regardless of `tracer_mask`.

      integer :: nx, ny, nz, ng, nxt, nyt, it
      logical :: do_h, do_u, do_v, do_all

      do_h = .true.
      do_u = .true.
      do_v = .true.
      if (present(fill_h)) do_h = fill_h
      if (present(fill_u)) do_u = fill_u
      if (present(fill_v)) do_v = fill_v
      do_all = .not. present(tracer_mask)
      if (present(all_tracers)) do_all = do_all .or. all_tracers

      nx = grid%nx_phys
      ny = grid%ny_phys
      nz = state%multilayer%nz_ml
      ng = grid%nghost
      nxt = grid%nx_total
      nyt = grid%ny_total

      if (do_h) call fill_ghosts_centre(state%multilayer%h_layer, nxt, nyt, nz, nx, ny, ng, &
                                        per_x, per_y)
      if (allocated(state%multilayer%tracers)) then
         do it = 1, size(state%multilayer%tracers)
            if (.not. allocated(state%multilayer%tracers(it)%c)) cycle
            if (.not. do_all) then
               if (it > size(tracer_mask)) cycle
               if (.not. tracer_mask(it)) cycle
            end if
            call fill_ghosts_centre(state%multilayer%tracers(it)%c, nxt, nyt, nz, &
                                    nx, ny, ng, per_x, per_y)
         end do
      end if
      if (do_u) call fill_ghosts_face_x(state%multilayer%u_face_x_layer, nxt + 1, nyt, nz, &
                                        nx, ny, ng, per_x, per_y)
      if (do_v) call fill_ghosts_face_y(state%multilayer%v_face_y_layer, nxt, nyt + 1, nz, &
                                        nx, ny, ng, per_x, per_y)
   end subroutine ocean_ic_fill_state_ghosts

   pure function shape_ok(a, n1, n2, n3) result(ok)
      real(wp), intent(in) :: a(:, :, :)
      integer, intent(in) :: n1, n2, n3
      logical :: ok
      ok = size(a, 1) == n1 .and. size(a, 2) == n2 .and. size(a, 3) == n3
   end function shape_ok

   subroutine shape_fail(what, got, want, ierr)
      character(len=*), intent(in) :: what
      integer, intent(in) :: got(3), want(3)
      integer, intent(out), optional :: ierr
      integer :: ierr_local
      call fail(what//": shape ("//to_string(got(1))//", "//to_string(got(2))//", "// &
                to_string(got(3))//") does not match the grid, expected ("// &
                to_string(want(1))//", "//to_string(want(2))//", "//to_string(want(3))// &
                ") with k = 1 the bed", ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
      if (present(ierr)) ierr = ierr_local
   end subroutine shape_fail

   subroutine apply_h(state, staged_h, is_staged, nx, ny, nz, ng, ierr)
      !! Validate `staged_h` on wet cells and write it; unstaged tracers are
      !! rescaled so their concentration survives the thickness change.
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: nx, ny, nz, ng
      real(wp), intent(in) :: staged_h(nx, ny, nz)
      logical, intent(in) :: is_staged(:)
      integer, intent(out) :: ierr

      integer :: i, j, k, it, ierr_local
      real(wp) :: col, water, hold
      logical :: wet

      ierr = OCEAN_STATUS_OK
      do j = 1, ny
         do i = 1, nx
            wet = state%multilayer%wet_mask(ng + i, ng + j) > 0.5_wp
            if (.not. wet) cycle
            col = 0.0_wp
            do k = 1, nz
               if (.not. (staged_h(i, j, k) > 0.0_wp)) then
                  call fail("stage_h: h <= 0 on a wet cell at (i,j,k) = ("//to_string(i)//", "// &
                            to_string(j)//", "//to_string(k)//"), h = "// &
                            to_string(staged_h(i, j, k))//" (wet = bathymetry above the land "// &
                            "threshold; shape is (nx, ny, nz), k = 1 the bed)", &
                            ierr_local, OCEAN_STATUS_ERR_IC_SEED)
                  ierr = ierr_local
                  return
               end if
               col = col + staged_h(i, j, k)
            end do
            water = state%barotropic%b(ng + i, ng + j)
            if (state%metrics%use_cavity) water = water - state%metrics%z_draft(ng + i, ng + j)
            if (abs(col - water) > IC_H_COLUMN_RTOL*water) then
               call fail("stage_h: column at (i,j) = ("//to_string(i)//", "//to_string(j)// &
                         ") sums to "//to_string(col)//" m but the water column is "// &
                         to_string(water)//" m (allowed |sum(h) - D| <= "// &
                         to_string(IC_H_COLUMN_RTOL)//" D); check units and that k = 1 is "// &
                         "the bed", ierr_local, OCEAN_STATUS_ERR_IC_SEED)
               ierr = ierr_local
               return
            end if
         end do
      end do

      ! Unstaged tracers keep their concentration: `c` is stored directly, so
      ! restating `h` leaves it untouched.
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (state%multilayer%wet_mask(ng + i, ng + j) <= 0.5_wp) cycle
               state%multilayer%h_layer(ng + i, ng + j, k) = staged_h(i, j, k)
            end do
         end do
      end do
   end subroutine apply_h

   pure subroutine set_conc(c, wet_mask, conc, nx, ny, nz, ng, nxt, nyt)
      !! `c = conc` on wet interior cells (land cells keep the seed).
      integer, intent(in) :: nx, ny, nz, ng, nxt, nyt
      real(wp), intent(inout) :: c(nxt, nyt, nz)
      real(wp), intent(in) :: wet_mask(nxt, nyt)
      real(wp), intent(in) :: conc(nx, ny, nz)
      integer :: i, j, k
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (wet_mask(ng + i, ng + j) <= 0.5_wp) cycle
               c(ng + i, ng + j, k) = conc(i, j, k)
            end do
         end do
      end do
   end subroutine set_conc

   pure subroutine fill_ghosts_centre(f, nxt, nyt, nz, nx, ny, ng, per_x, per_y)
      !! Ghosts of a cell-centred field: periodic image on a periodic axis,
      !! else zero-gradient (nearest interior cell).  x first, then y over
      !! the full x range, so corners are consistent.
      integer, intent(in) :: nxt, nyt, nz, nx, ny, ng
      real(wp), intent(inout) :: f(nxt, nyt, nz)
      logical, intent(in) :: per_x, per_y
      integer :: i, j, k
      do k = 1, nz
         do j = ng + 1, ng + ny
            do i = 1, ng
               if (per_x) then
                  f(i, j, k) = f(i + nx, j, k)
                  f(ng + nx + i, j, k) = f(ng + i, j, k)
               else
                  f(i, j, k) = f(ng + 1, j, k)
                  f(ng + nx + i, j, k) = f(ng + nx, j, k)
               end if
            end do
         end do
         do i = 1, nxt
            do j = 1, ng
               if (per_y) then
                  f(i, j, k) = f(i, j + ny, k)
                  f(i, ng + ny + j, k) = f(i, ng + j, k)
               else
                  f(i, j, k) = f(i, ng + 1, k)
                  f(i, ng + ny + j, k) = f(i, ng + ny, k)
               end if
            end do
         end do
      end do
   end subroutine fill_ghosts_centre

   pure subroutine fill_ghosts_face_x(f, nxf, nyt, nz, nx, ny, ng, per_x, per_y)
      !! Ghosts of an x-face field (interior faces ng+1 .. ng+nx+1).
      integer, intent(in) :: nxf, nyt, nz, nx, ny, ng
      real(wp), intent(inout) :: f(nxf, nyt, nz)
      logical, intent(in) :: per_x, per_y
      integer :: i, j, k
      do k = 1, nz
         do j = ng + 1, ng + ny
            do i = 1, ng
               if (per_x) then
                  f(i, j, k) = f(i + nx, j, k)
                  f(ng + nx + 1 + i, j, k) = f(ng + 1 + i, j, k)
               else
                  f(i, j, k) = f(ng + 1, j, k)
                  f(ng + nx + 1 + i, j, k) = f(ng + nx + 1, j, k)
               end if
            end do
         end do
         do i = 1, nxf
            do j = 1, ng
               if (per_y) then
                  f(i, j, k) = f(i, j + ny, k)
                  f(i, ng + ny + j, k) = f(i, ng + j, k)
               else
                  f(i, j, k) = f(i, ng + 1, k)
                  f(i, ng + ny + j, k) = f(i, ng + ny, k)
               end if
            end do
         end do
      end do
   end subroutine fill_ghosts_face_x

   pure subroutine fill_ghosts_face_y(f, nxt, nyf, nz, nx, ny, ng, per_x, per_y)
      !! Ghosts of a y-face field (interior faces ng+1 .. ng+ny+1).
      integer, intent(in) :: nxt, nyf, nz, nx, ny, ng
      real(wp), intent(inout) :: f(nxt, nyf, nz)
      logical, intent(in) :: per_x, per_y
      integer :: i, j, k
      do k = 1, nz
         do j = ng + 1, ng + ny + 1
            do i = 1, ng
               if (per_x) then
                  f(i, j, k) = f(i + nx, j, k)
                  f(ng + nx + i, j, k) = f(ng + i, j, k)
               else
                  f(i, j, k) = f(ng + 1, j, k)
                  f(ng + nx + i, j, k) = f(ng + nx, j, k)
               end if
            end do
         end do
         do i = 1, nxt
            do j = 1, ng
               if (per_y) then
                  f(i, j, k) = f(i, j + ny, k)
                  f(i, ng + ny + 1 + j, k) = f(i, ng + 1 + j, k)
               else
                  f(i, j, k) = f(i, ng + 1, k)
                  f(i, ng + ny + 1 + j, k) = f(i, ng + ny + 1, k)
               end if
            end do
         end do
      end do
   end subroutine fill_ghosts_face_y

end module rdb_ocean_ic_inject
