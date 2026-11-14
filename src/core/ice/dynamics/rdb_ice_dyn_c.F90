module rdb_ice_dyn_c
   !! Engine driver of the CICE6 C-grid EVP dynamics under the SIS2-derived ice
   !! state (`&ocean_ice_dyn_nml dyn_scheme = "cice_c"`, CICE6 port D9.2): the
   !! coupling kernels that feed `rdb_ice_evp_c` from the ocean state and the
   !! ice aggregates, the subcycle loop with the face halo, and the hand-off of
   !! the ice-ocean stress.  The solver is Lemieux et al. (2024) (see
   !! `rdb_ice_evp_c`); the coupling follows the C-grid design of the port plan
   !! (sea-surface tilt from the ocean surface height at the faces, ocean
   !! velocity at the faces, `strocnT` from the face stress per ice area).
   !!
   !! Roundabout-original (no upstream code).
   !!
   !! **Index map** (see `rdb_ice_cgrid_metrics`): `E(i,j) = u(i+1,j)` is the
   !! east face of `T(i,j)`, `N(i,j) = v(i,j+1)` its north face.  The ice
   !! velocity is `ocean_sea_ice_t%u_ice/v_ice`, the arrays the category
   !! transport already reads, so the transport follows the new velocities with
   !! no further wiring.
   !!
   !! **Ice state seen by the solver.**  The scheme-neutral aggregates of the
   !! ice slot (`aice`, `vice`, `vsno`, `aice0`, filled by `ice_aggregates_fill`
   !! after every thermodynamics window; ghosts valid) with the SIS2 densities,
   !! so the mass the dynamics moves is the mass of the state.  The strength is
   !! the Hibler (1979) strength of the aggregates (`kstrength = 0`); the start
   !! of the step concentration is the current one (the dynamics runs on the
   !! state it finds).
   !!
   !! **One step** (`ice_dyn_c_step`), all kernels device-resident:
   !!
   !!   1. `ssh = sum_k h_layer - b` at T (one launch);
   !!   2. the inputs at E and N in ONE launch over the owned cells: Coriolis
   !!      (mean of the two end-point corner values of the face), sea-surface
   !!      slope `(ssh(i+1) - ssh(i))/dxCu` (and `dyCv` in y), wind stress on the
   !!      ice (the snapshot `tau_a` weighted by the face ice concentration, the
   !!      concentration weighting the drag uses, so the momentum budget of the
   !!      pair ice + ocean closes at fractional cover), the ocean velocity
   !!      component along each face from the surface layer of the face and the
   !!      tangential one as the area-weighted mean of the four surrounding faces
   !!      (the same weights the solver uses for the ice velocity).  The
   !!      surface layer of a face is the first LIVE layer at or below
   !!      `k_top_u/v` of both neighbouring columns;
   !!   3. `ice_evp_c_prep`, the face halo, `ndte` subcycles of
   !!      `ice_evp_c_substep` with `ice_halo_faces` after each, the average and
   !!      `ice_evp_c_finish` (the strain-rate diagnostics, `strocnxE/yN`, the
   !!      velocity guard);
   !!   4. the face stress per ice area, the mirrors of the two face masks and
   !!      `strocnT` (per ice area at T: `0.5 (w_E(i-1) A_E(i-1) + w_E(i) A_E(i))
   !!      / A_T`, `w = strocn/ai` where `ai > puny`; the division by the face
   !!      area comes BEFORE the average to T);
   !!   5. one reduction of the error words; a fatal bit (`ICE_ERR_DYN_VEL` for
   !!      a non-finite or above-cap velocity, the strength guard, ...) aborts
   !!      through `ice_err_enforce` -- nothing is ever clamped.
   !!
   !! The hand-off of `-strocn` to the ocean surface stress is
   !! `ice_ocean_stress_flux_c` (`rdb_ice_ocean_coupler`).
   !!
   !! Single rank, no tripolar fold (refused at configure; D9.3 and the
   !! multi-rank PR lift them).  The face halo of every subcycle is
   !! `ice_halo_faces` (`rdb_ice_cgrid_metrics`), which on a periodic axis keeps
   !! the east/north seam face the solver computes and copies it onto the
   !! duplicate (the ocean's own wrap goes the other way and would clobber it).
   !! For the same reason the blended stress at the west/south seam face is
   !! copied from the east/north one before the ocean's seam refresh; see
   !! `ice_ocean_stress_flux_c`.
   !!
   !! Open ocean boundaries are land for the ice (the dynamics envelope allows
   !! wall or periodic edges only).  `strocnT_x/y` is restart state with no
   !! consumer until the basal-heat wiring.  When the ice mass loads the
   !! surface the tilt must use the pressure-loaded surface, not `sum h - b`.
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d
   use rdb_ice_column, only: ICE_RHO_ICE, ICE_RHO_SNOW
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_errors, only: ICE_ERR_FATAL_MASK, ice_err_summary_t, ice_err_clear, ice_err_reduce, &
                             ice_err_enforce
   use rdb_ice_evp_c_prep, only: ice_evp_c_const_t, ice_evp_c_constants, ice_evp_c_prep, ice_evp_c_ghost_ok
   use rdb_ice_evp_c, only: ice_evp_c_substep, ice_evp_c_average, ice_evp_c_finish, ice_evp_c_copy_faces
   use rdb_ice_cgrid_metrics, only: ice_cgrid_metrics_t, ice_halo_faces
   implicit none
   private

   public :: ice_dyn_c_metrics_build, ice_dyn_c_step
   public :: dync_ssh_impl, dync_inputs_impl, dync_post_impl, dync_strocnt_impl

contains

   subroutine ice_dyn_c_metrics_build(grid, metrics, bc, ice, cm, ierr)
      !! Build the C-grid metric set `cm` (`ocean_state_t%ice_cgrid`) from the
      !! FINAL ocean metrics (land mask applied) when the C-grid slot is live,
      !! and map it; the ice slot is mapped by the state's own `enter_data`,
      !! which must come after this call.  Host, configure time.  `ierr /= 0`:
      !! the window is too narrow for the solver (`nghost < 3`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_bc_state_t), intent(in) :: bc
      type(ocean_sea_ice_t), intent(in) :: ice
      type(ice_cgrid_metrics_t), intent(inout) :: cm
      integer, intent(out) :: ierr
      integer :: ng

      ierr = 0
      if (.not. ice%dync%is_init) return
      ng = grid%nghost
      if (.not. ice_evp_c_ghost_ok(grid%nx_total, grid%ny_total, ng + 1, ng + grid%nx_phys, &
                                   ng + 1, ng + grid%ny_phys)) then
         ierr = 1
         return
      end if
      call cm%init(grid, metrics, bc, ierr)
   end subroutine ice_dyn_c_metrics_build

   subroutine ice_dyn_c_step(grid, metrics, f_corner, b, cm, ice, ms, bc, dt, err_summary)
      !! One outer dynamics step (see the module note).  Aborts on a fatal
      !! error bit.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), intent(in), contiguous :: f_corner(:, :)
         !! Coriolis parameter at the SW corners, `(nx+1, ny+1)`.
      real(wp), intent(in), contiguous :: b(:, :)
         !! Bathymetry depth at T, positive down [m] (`barotropic%b`).
      type(ice_cgrid_metrics_t), intent(in) :: cm
         !! The C-grid metric set (`ocean_state_t%ice_cgrid`), built and mapped.
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(in) :: dt
         !! Outer step [s].
      type(ice_err_summary_t), intent(out), optional :: err_summary
         !! When present the reduced error words are returned and NOTHING aborts
         !! (the test seam of the velocity guard); absent: a fatal bit stops the
         !! run through `ice_err_enforce`.

      type(ice_evp_c_const_t) :: k
      type(ice_err_summary_t) :: summ
      integer :: nx, ny, nz, ng, is, ie, js, je, ksub
      logical :: odd

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      ng = grid%nghost
      is = ng + 1
      ie = ng + grid%nx_phys
      js = ng + 1
      je = ng + grid%ny_phys

      k = ice_evp_c_constants(ice%dync%p, dt)
      k%rhoi = ICE_RHO_ICE
      k%rhos = ICE_RHO_SNOW
      odd = mod(k%ndte, 2) == 1

      call ice_err_clear(ice%dync%ice_err, nx, ny)
      call dync_ssh_impl(nx, ny, nz, ms%h_layer, b, metrics%wet_T, ice%dync%ssh)
      call dync_inputs_impl(nx, ny, nz, is, ie, js, je, ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
                            ms%k_top_u, ms%k_top_v, f_corner, ice%dync%ssh, metrics%idxCu, metrics%idyCv, &
                            cm%hm, cm%tarea, cm%earea, cm%narea, &
                            cm%epm, cm%npm, ice%aice, ice%tau_a_x, ice%tau_a_y, &
                            ice%dync%fcorE, ice%dync%fcorN, ice%dync%strairxE, ice%dync%strairyN, &
                            ice%dync%ss_tltxE, ice%dync%ss_tltyN, ice%dync%uocnE, ice%dync%vocnE, &
                            ice%dync%uocnN, ice%dync%vocnN)
      call ice_evp_c_prep(ice%dync%w, k, nx, ny, 1, is, ie, js, je, ice%dync%p%ridge, ice%dync%p%c, &
                          cm%hm, ice%aice, ice%vice, ice%vsno, ice%aice, ice%aice0, ice%aice, &
                          ice%vice, cm%tarea, cm%epm, cm%npm, ice%dync%fcorE, &
                          ice%dync%fcorN, ice%dync%strairxE, ice%dync%strairyN, ice%dync%ss_tltxE, &
                          ice%dync%ss_tltyN, ice%dync%uocnE, ice%dync%vocnN, ice%u_ice, ice%v_ice, &
                          ice%dync%ice_err)
      call ice_halo_faces(grid, bc, ice%u_ice, ice%v_ice)
      call ice_evp_c_copy_faces(nx, ny, ice%u_ice, ice%v_ice, ice%dync%w%u2, ice%dync%w%v2)
      do ksub = 1, k%ndte
         if (mod(ksub, 2) == 1) then
            call ice_evp_c_substep(ice%dync%w, k, nx, ny, is, ie, js, je, ice%u_ice, ice%v_ice, &
                                   ice%dync%w%u2, ice%dync%w%v2, &
                                   cm%dxE, cm%dyE, cm%dxN, cm%dyN, &
                                   cm%dxU, cm%dyU, cm%dxT, cm%dyT, &
                                   cm%uarea, cm%tarea, cm%earea, cm%narea, &
                                   cm%earear, cm%narear, cm%ratiodxN, &
                                   cm%ratiodxNr, cm%ratiodyE, cm%ratiodyEr, &
                                   cm%epm, cm%npm, cm%uvm, cm%hm, &
                                   ice%dync%uocnE, ice%dync%vocnE, ice%dync%uocnN, ice%dync%vocnN)
            call ice_halo_faces(grid, bc, ice%dync%w%u2, ice%dync%w%v2)
         else
            call ice_evp_c_substep(ice%dync%w, k, nx, ny, is, ie, js, je, ice%dync%w%u2, ice%dync%w%v2, &
                                   ice%u_ice, ice%v_ice, &
                                   cm%dxE, cm%dyE, cm%dxN, cm%dyN, &
                                   cm%dxU, cm%dyU, cm%dxT, cm%dyT, &
                                   cm%uarea, cm%tarea, cm%earea, cm%narea, &
                                   cm%earear, cm%narear, cm%ratiodxN, &
                                   cm%ratiodxNr, cm%ratiodyE, cm%ratiodyEr, &
                                   cm%epm, cm%npm, cm%uvm, cm%hm, &
                                   ice%dync%uocnE, ice%dync%vocnE, ice%dync%uocnN, ice%dync%vocnN)
            call ice_halo_faces(grid, bc, ice%u_ice, ice%v_ice)
         end if
      end do
      if (odd) call ice_evp_c_copy_faces(nx, ny, ice%dync%w%u2, ice%dync%w%v2, ice%u_ice, ice%v_ice)
      call ice_evp_c_average(ice%dync%w, nx, ny, is, ie, js, je, ice%u_ice, ice%v_ice, cm%earea, &
                             cm%narea, cm%epm, cm%npm, cm%uvm)
      call ice_evp_c_finish(ice%dync%w, k, nx, ny, is, ie, js, je, ice%u_ice, ice%v_ice, cm%dxN, &
                            cm%dyE, cm%dxT, cm%dyT, cm%uarea, &
                            cm%tarear, ice%dync%uocnE, ice%dync%vocnE, ice%dync%uocnN, &
                            ice%dync%vocnN, ice%dync%ice_err)

      call dync_post_impl(nx, ny, is, ie, js, je, ice%dync%p%c%puny, ice%dync%w%strocnxE, &
                          ice%dync%w%strocnyN, ice%dync%w%aiE, ice%dync%w%aiN, ice%dync%w%iceEmask, &
                          ice%dync%w%iceNmask, ice%dync%wE, ice%dync%wN, ice%dync%emask_r, ice%dync%nmask_r)
      call ocean_periodic_wrap_centre_2d(ice%dync%wE, nx, ny, grid%nx_phys, grid%ny_phys, ng, &
                                         bc%periodic_x, bc%periodic_y)
      call ocean_periodic_wrap_centre_2d(ice%dync%wN, nx, ny, grid%nx_phys, grid%ny_phys, ng, &
                                         bc%periodic_x, bc%periodic_y)
      call dync_strocnt_impl(nx, ny, is, ie, js, je, ice%dync%wE, ice%dync%wN, cm%earea, &
                             cm%narea, cm%tarear, ice%dync%strocnT_x, ice%dync%strocnT_y)

      call ice_err_reduce(ice%dync%ice_err, nx, ny, is, ie, js, je, 0, 0, grid%nx_phys, &
                          ICE_ERR_FATAL_MASK, summ)
      if (present(err_summary)) then
         err_summary = summ
      else
         call ice_err_enforce(summ, grid%nx_phys, ICE_ERR_FATAL_MASK)
      end if
   end subroutine ice_dyn_c_step

   ! ------------------------------------------------------------------
   ! kernels (explicit-shape, flat; every array device-present)
   ! ------------------------------------------------------------------

   pure function wavg4(a1, a2, a3, a4, w1, w2, w3, w4) result(s)
      !! Area-weighted mean of four face values (0 when the weights sum to 0):
      !! the operands and their order of the solver's tangential interpolation.
      !$acc routine seq
      real(wp), value :: a1, a2, a3, a4, w1, w2, w3, w4
      real(wp) :: s
      real(wp) :: wt
      wt = w1 + w2 + w3 + w4
      s = 0.0_wp
      if (wt /= 0.0_wp) s = (a1*w1 + a2*w2 + a3*w3 + a4*w4)/wt
   end function wavg4

   pure function live_u(nx, ny, nz, h_layer, i, j, k0) result(kl)
      !! Surface layer of the x-face `(i,j)` (between the columns `i-1` and
      !! `i`): the first layer counting down from `k0` that is live (strictly
      !! above `H_VANISHED`) in BOTH columns; `k0` itself when none is.
      !$acc routine seq
      integer, value :: nx, ny, nz, i, j, k0
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      integer :: kl
      integer :: kk
      kl = k0
      do kk = k0, 1, -1
         if (h_layer(i - 1, j, kk) > H_VANISHED .and. h_layer(i, j, kk) > H_VANISHED) then
            kl = kk
            exit
         end if
      end do
   end function live_u

   pure function live_v(nx, ny, nz, h_layer, i, j, k0) result(kl)
      !! Surface layer of the y-face `(i,j)` (between the rows `j-1` and `j`),
      !! as `live_u`.
      !$acc routine seq
      integer, value :: nx, ny, nz, i, j, k0
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      integer :: kl
      integer :: kk
      kl = k0
      do kk = k0, 1, -1
         if (h_layer(i, j - 1, kk) > H_VANISHED .and. h_layer(i, j, kk) > H_VANISHED) then
            kl = kk
            exit
         end if
      end do
   end function live_v

   subroutine dync_ssh_impl(nx, ny, nz, h_layer, b, wet_T, ssh)
      !! Sea-surface height at T, `ssh = sum_k h_layer - b` on wet cells (the
      !! identity of the SSH diagnostic), 0 on land.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz), b(nx, ny), wet_T(nx, ny)
      real(wp), intent(inout) :: ssh(nx, ny)
      integer :: i, j, k
      real(wp) :: col

      do concurrent(j=1:ny, i=1:nx) local(k, col)
         col = 0.0_wp
         if (wet_T(i, j) > 0.5_wp) then
            do k = 1, nz
               col = col + h_layer(i, j, k)
            end do
            col = col - b(i, j)
         end if
         ssh(i, j) = col
      end do
   end subroutine dync_ssh_impl

   subroutine dync_inputs_impl(nx, ny, nz, is, ie, js, je, h_layer, u_x, v_y, k_top_u, k_top_v, f_corner, &
                               ssh, idxCu, idyCv, hm, tarea, earea, narea, epm, npm, aice, tau_a_x, tau_a_y, &
                               fcorE, fcorN, strairxE, strairyN, ss_tltxE, ss_tltyN, uocnE, vocnE, uocnN, vocnN)
      !! The forcing of the solver at the owned faces `E(is:ie, js:je)` and
      !! `N(is:ie, js:je)`, one launch.  `E(i,j)` is the x-face `i+1`, `N(i,j)`
      !! the y-face `j+1`.
      !!   * `fcorE`, `fcorN`: mean of the Coriolis parameter at the two ends of
      !!     the face;
      !!   * `ss_tltxE = (ssh(i+1,j) - ssh(i,j)) idxCu(i+1,j)` and
      !!     `ss_tltyN = (ssh(i,j+1) - ssh(i,j)) idyCv(i,j+1)`, zero on a closed
      !!     face;
      !!   * `strairxE`, `strairyN`: the wind stress snapshot times the face ice
      !!     concentration `a = (hm A_T aice)` two-point mean (the expression of
      !!     the preparation, so the ice and the ocean weigh the same `a`);
      !!   * the ocean velocity: the component along the face from the face's
      !!     surface layer, the tangential one as the mean of the four
      !!     surrounding faces weighted by their face area (the weights of the
      !!     solver's interpolation), masked by the face mask.
      integer, intent(in) :: nx, ny, nz, is, ie, js, je
      real(wp), intent(in) :: h_layer(nx, ny, nz), u_x(nx + 1, ny, nz), v_y(nx, ny + 1, nz)
      integer, intent(in) :: k_top_u(nx + 1, ny), k_top_v(nx, ny + 1)
      real(wp), intent(in) :: f_corner(nx + 1, ny + 1), ssh(nx, ny), idxCu(nx + 1, ny), idyCv(nx, ny + 1)
      real(wp), intent(in) :: hm(nx, ny), tarea(nx, ny), earea(nx, ny), narea(nx, ny), epm(nx, ny), npm(nx, ny)
      real(wp), intent(in) :: aice(nx, ny), tau_a_x(nx + 1, ny), tau_a_y(nx, ny + 1)
      real(wp), intent(inout) :: fcorE(nx, ny), fcorN(nx, ny), strairxE(nx, ny), strairyN(nx, ny)
      real(wp), intent(inout) :: ss_tltxE(nx, ny), ss_tltyN(nx, ny)
      real(wp), intent(inout) :: uocnE(nx, ny), vocnE(nx, ny), uocnN(nx, ny), vocnN(nx, ny)
      integer :: i, j
      integer :: kA, kB, kC, kD
      real(wp) :: w1, w2, wt, a

      do concurrent(j=js:je, i=is:ie) local(kA, kB, kC, kD, w1, w2, wt, a)
         ! ---- Coriolis and sea-surface slope ----
         fcorE(i, j) = 0.5_wp*(f_corner(i + 1, j) + f_corner(i + 1, j + 1))
         fcorN(i, j) = 0.5_wp*(f_corner(i, j + 1) + f_corner(i + 1, j + 1))
         ss_tltxE(i, j) = (ssh(i + 1, j) - ssh(i, j))*idxCu(i + 1, j)*epm(i, j)
         ss_tltyN(i, j) = (ssh(i, j + 1) - ssh(i, j))*idyCv(i, j + 1)*npm(i, j)
         ! ---- wind stress on the ice: the snapshot times the face concentration ----
         w1 = hm(i, j)*tarea(i, j)
         w2 = hm(i + 1, j)*tarea(i + 1, j)
         wt = w1 + w2
         a = 0.0_wp
         if (wt /= 0.0_wp) a = (w1*aice(i, j) + w2*aice(i + 1, j))/wt
         strairxE(i, j) = a*tau_a_x(i + 1, j)
         w2 = hm(i, j + 1)*tarea(i, j + 1)
         wt = w1 + w2
         a = 0.0_wp
         if (wt /= 0.0_wp) a = (w1*aice(i, j) + w2*aice(i, j + 1))/wt
         strairyN(i, j) = a*tau_a_y(i, j + 1)
         ! ---- ocean velocity at E: u along, v from the four surrounding y-faces ----
         kA = live_u(nx, ny, nz, h_layer, i + 1, j, k_top_u(i + 1, j))
         uocnE(i, j) = u_x(i + 1, j, kA)*epm(i, j)
         kA = live_v(nx, ny, nz, h_layer, i, j, k_top_v(i, j))
         kB = live_v(nx, ny, nz, h_layer, i + 1, j, k_top_v(i + 1, j))
         kC = live_v(nx, ny, nz, h_layer, i, j + 1, k_top_v(i, j + 1))
         kD = live_v(nx, ny, nz, h_layer, i + 1, j + 1, k_top_v(i + 1, j + 1))
         vocnE(i, j) = wavg4(v_y(i, j, kA), v_y(i + 1, j, kB), v_y(i, j + 1, kC), v_y(i + 1, j + 1, kD), &
                             narea(i, j - 1), narea(i + 1, j - 1), narea(i, j), narea(i + 1, j))*epm(i, j)
         ! ---- ocean velocity at N: v along, u from the four surrounding x-faces ----
         kA = live_v(nx, ny, nz, h_layer, i, j + 1, k_top_v(i, j + 1))
         vocnN(i, j) = v_y(i, j + 1, kA)*npm(i, j)
         kA = live_u(nx, ny, nz, h_layer, i, j, k_top_u(i, j))
         kB = live_u(nx, ny, nz, h_layer, i + 1, j, k_top_u(i + 1, j))
         kC = live_u(nx, ny, nz, h_layer, i, j + 1, k_top_u(i, j + 1))
         kD = live_u(nx, ny, nz, h_layer, i + 1, j + 1, k_top_u(i + 1, j + 1))
         uocnN(i, j) = wavg4(u_x(i, j, kA), u_x(i + 1, j, kB), u_x(i, j + 1, kC), u_x(i + 1, j + 1, kD), &
                             earea(i - 1, j), earea(i, j), earea(i - 1, j + 1), earea(i, j + 1))*npm(i, j)
      end do
   end subroutine dync_inputs_impl

   subroutine dync_post_impl(nx, ny, is, ie, js, je, puny, strocnxE, strocnyN, aiE, aiN, iceEmask, iceNmask, &
                             wE, wN, emask_r, nmask_r)
      !! Over the owned cells: the face stress per ICE area `w = strocn/ai`
      !! (0 where `ai <= puny`) and the 0/1 mirrors of the two face ice masks.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: strocnxE(nx, ny), strocnyN(nx, ny), aiE(nx, ny), aiN(nx, ny)
      logical, intent(in) :: iceEmask(nx, ny), iceNmask(nx, ny)
      real(wp), intent(inout) :: wE(nx, ny), wN(nx, ny), emask_r(nx, ny), nmask_r(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         wE(i, j) = 0.0_wp
         wN(i, j) = 0.0_wp
         if (aiE(i, j) > puny) wE(i, j) = strocnxE(i, j)/aiE(i, j)
         if (aiN(i, j) > puny) wN(i, j) = strocnyN(i, j)/aiN(i, j)
         emask_r(i, j) = merge(1.0_wp, 0.0_wp, iceEmask(i, j))
         nmask_r(i, j) = merge(1.0_wp, 0.0_wp, iceNmask(i, j))
      end do
   end subroutine dync_post_impl

   subroutine dync_strocnt_impl(nx, ny, is, ie, js, je, wE, wN, earea, narea, tarear, strocnT_x, strocnT_y)
      !! `strocnT` per ice area at the owned T cells: the flux-conserving
      !! two-face mean of the face stress per ice area, `0.5 (w(i-1) A(i-1) +
      !! w(i) A(i)) / A_T` (x from E, y from N).  The faces `E(is-1)`,
      !! `N(js-1)` must hold their wrapped values (or 0 at a wall).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: wE(nx, ny), wN(nx, ny), earea(nx, ny), narea(nx, ny), tarear(nx, ny)
      real(wp), intent(inout) :: strocnT_x(nx, ny), strocnT_y(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         strocnT_x(i, j) = 0.5_wp*(wE(i - 1, j)*earea(i - 1, j) + wE(i, j)*earea(i, j))*tarear(i, j)
         strocnT_y(i, j) = 0.5_wp*(wN(i, j - 1)*narea(i, j - 1) + wN(i, j)*narea(i, j))*tarear(i, j)
      end do
   end subroutine dync_strocnt_impl

end module rdb_ice_dyn_c
