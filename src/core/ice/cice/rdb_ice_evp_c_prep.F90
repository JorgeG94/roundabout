module rdb_ice_evp_c_prep
   !! Sea-ice C-grid EVP dynamics: the per-step preparation (masses, ice-free
   !! masks, forcing terms at the E/N points, ice strength), the all-scalar
   !! constants derived from `ice_dyn_params_t`, and the persistent work/state
   !! type of the solver in `rdb_ice_evp_c`.  Hunke and Dukowicz (1997), Hunke
   !! (2001), Bouillon et al. (2013) for the revised EVP, Lemieux et al. (2024)
   !! for the C-grid discretisation.
   !!
   !! Roundabout port of the CICE6 `grid_ice = 'C'` preparation
   !! (BSD-3, see `NOTICE`).
   !!
   !! **Index convention.**  Everything is addressed by the ice names of
   !! `rdb_ice_cgrid_metrics`: an `(nx, ny)` array holds, at `(i, j)`, the value
   !! at the T cell `T(i,j)`, its EAST face `E(i,j)`, its NORTH face `N(i,j)` or
   !! its NE corner `U(i,j)`.  Only three arrays are shifted, because they ARE
   !! the ocean's staggered arrays: `u(nx+1, ny)` (`E(i,j) = u(i+1,j)`),
   !! `v(nx, ny+1)` (`N(i,j) = v(i,j+1)`) and the persistent corner stress
   !! `stress12U(nx+1, ny+1)` (`U(i,j) = stress12U(i+1,j+1)`).
   !!
   !! **Ranges.**  The owned T cells are `is:ie x js:je`.  The solver needs the
   !! ghost ring of the T-cell state to be valid out to 3 cells and the face
   !! pair `(u, v)` to be valid out to 3 faces beyond the owned ones (the
   !! corner stress `stress12U(U(is-1,.))` needs the T viscosity of `T(is-1,.)`,
   !! whose Delta needs the shear of the corner `U(is-2,.)`): `nghost >= 3`
   !! (`is >= 4`), which is also the multi-rank minimum of the ocean.
   !! Every derived quantity of the ghost rings is recomputed locally from the
   !! haloed faces; nothing but `(u, v)` is exchanged, once per subcycle
   !! (see `rdb_ice_evp_c`).
   !!
   !! **Masks and means.**  Concentration, mass, Coriolis-scaled mass and the
   !! drag factor at E, N and U are the masked, area-weighted means of the
   !! neighbouring T cells (weights `hm * tarea`); the ice-free mask is
   !! `mask .and. a > dyn_area_min .and. m > dyn_mass_min`; the T mask is the
   !! 3x3 dilation of the cells that pass it (so a stress is computed one cell
   !! beyond the ice edge).  A face that becomes ice-covered starts at the
   !! ocean velocity.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_params_t, ice_constants_t, ice_ridge_params_t, &
                             ICE_VISC_AVG_STRENGTH
   use rdb_ice_ridge, only: ice_compute_strength
   implicit none
   private

   public :: ice_evp_c_const_t, ice_evp_c_work_t
   public :: ice_evp_c_constants, ice_evp_c_prep, ice_evp_c_ghost_ok

   type :: ice_evp_c_const_t
      !! Everything the kernels read as scalars (derived once per step on the
      !! host).  No default initialisers.
      integer :: ndte                 !! subcycles
      real(wp) :: dt                  !! dynamics time step [s]
      real(wp) :: revp                !! 0 classic EVP, 1 revised EVP
      real(wp) :: arlx1i              !! 1/alpha of the stress relaxation
      real(wp) :: denom1              !! 1/(1+arlx1i) classic, 1 revised
      real(wp) :: brlx                !! beta of the momentum relaxation
      real(wp) :: epp2i               !! 1/e_plasticpot^2
      real(wp) :: e_factor            !! e_yieldcurve^2 / e_plasticpot^4
      real(wp) :: Ktens               !! tensile strength fraction
      real(wp) :: deltamin            !! deltaminEVP
      logical :: cap_max              !! Hibler cap (max) rather than the sum
      real(wp) :: amin                !! dyn_area_min
      real(wp) :: mmin                !! dyn_mass_min
      real(wp) :: rheo_amin           !! rheo_area_min
      logical :: free_drift           !! free_drift_approx
      logical :: avg_strength         !! visc_method = avg_strength
      real(wp) :: rhoi, rhos, rhow    !! densities [kg/m3]
      real(wp) :: gravit              !! gravity [m/s2]
      real(wp) :: cw                  !! ice-ocean drag coefficient
      real(wp) :: cosw, sinw          !! turning angle
      real(wp) :: vel_cap             !! face-speed ceiling [m/s]
   end type ice_evp_c_const_t

   type :: ice_evp_c_work_t
      !! Persistent state and work arrays of the C-grid EVP; allocated once
      !! (`init`), mapped once (`enter_data`).  `(nx, ny)` unless noted.
      integer :: nx = 0
      integer :: ny = 0
      logical :: mapped = .false.
      ! ---- persistent state ----
      ! Restart (coupling PR): stresspT, stressmT, stress12U and iceEmask/iceNmask are CARRIED state
      ! (the masks decide the start velocity of newly covered faces, so a bit-exact restart needs
      ! them).  zetax2T/etax2T are recomputed every step (zeroed in prep): not restart state.  The
      ! ghost-ring stresses evolve locally and are never exchanged: after reading a restart they
      ! must be restored from the checkpoint or exchanged once.
      real(wp), allocatable :: stresspT(:, :)   !! sigma11 + sigma22 at T [N/m]
      real(wp), allocatable :: stressmT(:, :)   !! sigma11 - sigma22 at T [N/m]
      real(wp), allocatable :: zetax2T(:, :)    !! 2 * bulk viscosity at T [N/m s]
      real(wp), allocatable :: etax2T(:, :)     !! 2 * shear viscosity at T [N/m s]
      real(wp), allocatable :: stress12U(:, :)  !! sigma12 at the corner, `(nx+1, ny+1)`
      logical, allocatable :: iceEmask(:, :)    !! last step's ice mask at E (new-ice start)
      logical, allocatable :: iceNmask(:, :)    !! last step's ice mask at N
      ! ---- per-step preparation ----
      real(wp), allocatable :: tmass(:, :)      !! ice + snow mass at T [kg/m2]
      logical, allocatable :: tmphm(:, :)       !! cell passes the thresholds
      logical, allocatable :: iceTmask(:, :)    !! 3x3 dilation of `tmphm`
      logical, allocatable :: iceUmask(:, :)    !! ice mask at the corner
      real(wp), allocatable :: strength(:, :)   !! ice strength at T [N/m]
      real(wp), allocatable :: strengthU(:, :)  !! strength at the corner (avg_strength only)
      real(wp), allocatable :: aiE(:, :)        !! concentration at E
      real(wp), allocatable :: aiN(:, :)        !! concentration at N
      real(wp), allocatable :: emassdti(:, :)   !! mass at E over dt [kg/m2/s]
      real(wp), allocatable :: nmassdti(:, :)   !! mass at N over dt [kg/m2/s]
      real(wp), allocatable :: fmE(:, :)        !! Coriolis * mass at E [kg/s]
      real(wp), allocatable :: fmN(:, :)        !! Coriolis * mass at N [kg/s]
      real(wp), allocatable :: forcexE(:, :)    !! wind stress + sea-surface tilt, x at E [N/m2]
      real(wp), allocatable :: forceyN(:, :)    !! wind stress + sea-surface tilt, y at N [N/m2]
      real(wp), allocatable :: uinit(:, :)      !! u at E at the start of the step [m/s]
      real(wp), allocatable :: vinit(:, :)      !! v at N at the start of the step [m/s]
      ! ---- per-subcycle ----
      real(wp), allocatable :: uvelN(:, :)      !! u interpolated to N [m/s]
      real(wp), allocatable :: vvelE(:, :)      !! v interpolated to E [m/s]
      real(wp), allocatable :: uvelU(:, :)      !! u interpolated to the corner [m/s]
      real(wp), allocatable :: vvelU(:, :)      !! v interpolated to the corner [m/s]
      real(wp), allocatable :: shearU(:, :)     !! shear strain rate x area at the corner [m2/s]
      real(wp), allocatable :: deltaU(:, :)     !! Delta x area at the corner (avg_strength only)
      real(wp), allocatable :: strintxE(:, :)   !! x stress divergence at E (written by the diagnostic routine only)
      real(wp), allocatable :: strintyN(:, :)   !! y stress divergence at N (idem)
      ! ---- results of the step ----
      real(wp), allocatable :: divu(:, :)       !! divergence [1/s]
      real(wp), allocatable :: shear(:, :)      !! shear strain rate [1/s]
      real(wp), allocatable :: vort(:, :)       !! vorticity [1/s]
      real(wp), allocatable :: rdg_conv(:, :)   !! ridging convergence [1/s]
      real(wp), allocatable :: rdg_shear(:, :)  !! ridging shear [1/s]
      real(wp), allocatable :: strocnxE(:, :)   !! x ice-ocean stress at E, force ON the ice [N/m2]
      real(wp), allocatable :: strocnyN(:, :)   !! y ice-ocean stress at N, force ON the ice [N/m2]
   contains
      procedure :: init => work_init
      procedure :: destroy => work_destroy
      procedure :: enter_data => work_enter_data
      procedure :: exit_data => work_exit_data
      procedure :: bytes => work_bytes
   end type ice_evp_c_work_t

contains

   pure function ice_evp_c_ghost_ok(nx, ny, is, ie, js, je) result(ok)
      !! One-time validity check of the owned window (call at configure or
      !! driver setup; the kernel path has no check): the solver reaches
      !! three cells into the ghost rings (`is >= 4`, `ie <= nx - 3`, likewise
      !! in y).  The width is set by `visc_method = "avg_zeta"` (the corner
      !! shear viscosity averages the T viscosity of the ghost ring, whose
      !! Delta needs the shear of the corner two cells out); `avg_strength`
      !! needs only two in principle, but the kernels run one range for both,
      !! so 3 is required for either.
      integer, intent(in) :: nx, ny, is, ie, js, je
      logical :: ok
      ok = is >= 4 .and. js >= 4 .and. ie <= nx - 3 .and. je <= ny - 3 .and. ie >= is .and. je >= js
   end function ice_evp_c_ghost_ok

   pure function ice_evp_c_constants(p, dt) result(k)
      !! The all-scalar constants of one dynamics step.  Classic EVP
      !! (Hunke 2001): `arlx = 2 E ndte`, `brlx = ndte`, `denom1 = 1/(1 +
      !! 1/arlx)`; revised EVP (Bouillon et al. 2013): `arlx`, `brlx` from the
      !! knobs, `denom1 = 1`, stress term `(1 - 1/arlx)`.  Does not validate
      !! (`ice_dyn_params_check` does).
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      type(ice_evp_c_const_t) :: k
      real(wp) :: arlx

      k%ndte = p%dyn%ndte
      k%dt = dt
      k%epp2i = 1.0_wp/(p%dyn%e_plasticpot*p%dyn%e_plasticpot)
      k%e_factor = p%dyn%e_yieldcurve*p%dyn%e_yieldcurve &
                   /(p%dyn%e_plasticpot**4)
      if (p%dyn%revised_evp) then
         k%revp = 1.0_wp
         k%denom1 = 1.0_wp
         arlx = p%dyn%arlx
         k%arlx1i = 1.0_wp/arlx
         k%brlx = p%dyn%brlx
      else
         k%revp = 0.0_wp
         arlx = 2.0_wp*p%dyn%elasticDamp*real(p%dyn%ndte, wp)
         k%arlx1i = 1.0_wp/arlx
         k%brlx = real(p%dyn%ndte, wp)
         k%denom1 = 1.0_wp/(1.0_wp + k%arlx1i)
      end if
      k%Ktens = p%dyn%Ktens
      k%deltamin = p%dyn%deltaminEVP
      k%cap_max = p%dyn%capping_max
      k%amin = p%dyn%dyn_area_min
      k%mmin = p%dyn%dyn_mass_min
      k%rheo_amin = p%dyn%rheo_area_min
      k%free_drift = p%dyn%free_drift_approx
      k%avg_strength = p%dyn%visc_method == ICE_VISC_AVG_STRENGTH
      k%rhoi = p%c%rhoi
      k%rhos = p%c%rhos
      k%rhow = p%c%rhow
      k%gravit = p%c%gravit
      k%cw = p%dyn%cdn_ocn
      k%cosw = p%dyn%cosw
      k%sinw = p%dyn%sinw
      k%vel_cap = p%dyn%vel_cap
   end function ice_evp_c_constants

   ! ------------------------------------------------------------------
   ! work type
   ! ------------------------------------------------------------------

   subroutine work_init(this, nx, ny)
      !! Allocate and zero every array.  Re-initialising is refused.
      class(ice_evp_c_work_t), intent(inout) :: this
      integer, intent(in) :: nx, ny

      if (allocated(this%stresspT)) error stop "ice_evp_c_work_t%init: already initialised"
      this%nx = nx
      this%ny = ny
      allocate (this%stresspT(nx, ny), this%stressmT(nx, ny), this%zetax2T(nx, ny), this%etax2T(nx, ny))
      allocate (this%stress12U(nx + 1, ny + 1))
      allocate (this%iceEmask(nx, ny), this%iceNmask(nx, ny))
      allocate (this%tmass(nx, ny), this%tmphm(nx, ny), this%iceTmask(nx, ny), this%iceUmask(nx, ny))
      allocate (this%strength(nx, ny), this%strengthU(nx, ny))
      allocate (this%aiE(nx, ny), this%aiN(nx, ny), this%emassdti(nx, ny), this%nmassdti(nx, ny))
      allocate (this%fmE(nx, ny), this%fmN(nx, ny), this%forcexE(nx, ny), this%forceyN(nx, ny))
      allocate (this%uinit(nx, ny), this%vinit(nx, ny))
      allocate (this%uvelN(nx, ny), this%vvelE(nx, ny), this%uvelU(nx, ny), this%vvelU(nx, ny))
      allocate (this%shearU(nx, ny), this%deltaU(nx, ny), this%strintxE(nx, ny), this%strintyN(nx, ny))
      allocate (this%divu(nx, ny), this%shear(nx, ny), this%vort(nx, ny))
      allocate (this%rdg_conv(nx, ny), this%rdg_shear(nx, ny))
      allocate (this%strocnxE(nx, ny), this%strocnyN(nx, ny))
      this%stresspT = 0.0_wp
      this%stressmT = 0.0_wp
      this%zetax2T = 0.0_wp
      this%etax2T = 0.0_wp
      this%stress12U = 0.0_wp
      this%iceEmask = .false.
      this%iceNmask = .false.
      this%tmass = 0.0_wp
      this%tmphm = .false.
      this%iceTmask = .false.
      this%iceUmask = .false.
      this%strength = 0.0_wp
      this%strengthU = 0.0_wp
      this%aiE = 0.0_wp
      this%aiN = 0.0_wp
      this%emassdti = 0.0_wp
      this%nmassdti = 0.0_wp
      this%fmE = 0.0_wp
      this%fmN = 0.0_wp
      this%forcexE = 0.0_wp
      this%forceyN = 0.0_wp
      this%uinit = 0.0_wp
      this%vinit = 0.0_wp
      this%uvelN = 0.0_wp
      this%vvelE = 0.0_wp
      this%uvelU = 0.0_wp
      this%vvelU = 0.0_wp
      this%shearU = 0.0_wp
      this%deltaU = 0.0_wp
      this%strintxE = 0.0_wp
      this%strintyN = 0.0_wp
      this%divu = 0.0_wp
      this%shear = 0.0_wp
      this%vort = 0.0_wp
      this%rdg_conv = 0.0_wp
      this%rdg_shear = 0.0_wp
      this%strocnxE = 0.0_wp
      this%strocnyN = 0.0_wp
   end subroutine work_init

   subroutine work_destroy(this)
      class(ice_evp_c_work_t), intent(inout) :: this
      if (this%mapped) error stop "ice_evp_c_work_t%destroy: exit_data first"
      if (.not. allocated(this%stresspT)) return
      deallocate (this%stresspT, this%stressmT, this%zetax2T, this%etax2T, this%stress12U)
      deallocate (this%iceEmask, this%iceNmask, this%tmass, this%tmphm, this%iceTmask, this%iceUmask)
      deallocate (this%strength, this%strengthU, this%aiE, this%aiN, this%emassdti, this%nmassdti)
      deallocate (this%fmE, this%fmN, this%forcexE, this%forceyN, this%uinit, this%vinit)
      deallocate (this%uvelN, this%vvelE, this%uvelU, this%vvelU, this%shearU, this%deltaU)
      deallocate (this%strintxE, this%strintyN, this%divu, this%shear, this%vort)
      deallocate (this%rdg_conv, this%rdg_shear, this%strocnxE, this%strocnyN)
      this%nx = 0
      this%ny = 0
   end subroutine work_destroy

   pure function arr_bytes(a) result(nb)
      real(wp), allocatable, intent(in) :: a(:, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function arr_bytes

   pure function larr_bytes(a) result(nb)
      logical, allocatable, intent(in) :: a(:, :)
      integer(int64) :: nb
      nb = 0_int64
      if (allocated(a)) nb = int(size(a), int64)*int(storage_size(a)/8, int64)
   end function larr_bytes

   pure function work_bytes(this) result(nb)
      !! Allocated footprint in bytes.
      class(ice_evp_c_work_t), intent(in) :: this
      integer(int64) :: nb
      nb = arr_bytes(this%stresspT) &
           + arr_bytes(this%stressmT) &
           + arr_bytes(this%zetax2T) &
           + arr_bytes(this%etax2T) &
           + arr_bytes(this%stress12U) &
           + arr_bytes(this%tmass) &
           + arr_bytes(this%strength) &
           + arr_bytes(this%strengthU) &
           + arr_bytes(this%aiE) &
           + arr_bytes(this%aiN) &
           + arr_bytes(this%emassdti) &
           + arr_bytes(this%nmassdti) &
           + arr_bytes(this%fmE) &
           + arr_bytes(this%fmN) &
           + arr_bytes(this%forcexE) &
           + arr_bytes(this%forceyN) &
           + arr_bytes(this%uinit) &
           + arr_bytes(this%vinit) &
           + arr_bytes(this%uvelN) &
           + arr_bytes(this%vvelE) &
           + arr_bytes(this%uvelU) &
           + arr_bytes(this%vvelU) &
           + arr_bytes(this%shearU) &
           + arr_bytes(this%deltaU) &
           + arr_bytes(this%strintxE) &
           + arr_bytes(this%strintyN) &
           + arr_bytes(this%divu) &
           + arr_bytes(this%shear) &
           + arr_bytes(this%vort) &
           + arr_bytes(this%rdg_conv) &
           + arr_bytes(this%rdg_shear) &
           + arr_bytes(this%strocnxE) &
           + arr_bytes(this%strocnyN) &
           + larr_bytes(this%iceEmask) &
           + larr_bytes(this%iceNmask) &
           + larr_bytes(this%tmphm) &
           + larr_bytes(this%iceTmask) &
           + larr_bytes(this%iceUmask)
   end function work_bytes

   subroutine work_enter_data(this)
      !! Map every array to the device (`copyin`: the zero state).
      class(ice_evp_c_work_t), intent(inout) :: this
      !$acc enter data copyin(this%stresspT, this%stressmT, this%zetax2T, this%etax2T, this%stress12U)
      !$acc enter data copyin(this%iceEmask, this%iceNmask, this%tmass, this%tmphm, this%iceTmask, this%iceUmask)
      !$acc enter data copyin(this%strength, this%strengthU, this%aiE, this%aiN, this%emassdti, this%nmassdti)
      !$acc enter data copyin(this%fmE, this%fmN, this%forcexE, this%forceyN, this%uinit, this%vinit)
      !$acc enter data copyin(this%uvelN, this%vvelE, this%uvelU, this%vvelU, this%shearU, this%deltaU)
      !$acc enter data copyin(this%strintxE, this%strintyN, this%divu, this%shear, this%vort)
      !$acc enter data copyin(this%rdg_conv, this%rdg_shear, this%strocnxE, this%strocnyN)
      this%mapped = .true.
   end subroutine work_enter_data

   subroutine work_exit_data(this)
      class(ice_evp_c_work_t), intent(inout) :: this
      !$acc exit data delete(this%stresspT, this%stressmT, this%zetax2T, this%etax2T, this%stress12U)
      !$acc exit data delete(this%iceEmask, this%iceNmask, this%tmass, this%tmphm, this%iceTmask, this%iceUmask)
      !$acc exit data delete(this%strength, this%strengthU, this%aiE, this%aiN, this%emassdti, this%nmassdti)
      !$acc exit data delete(this%fmE, this%fmN, this%forcexE, this%forceyN, this%uinit, this%vinit)
      !$acc exit data delete(this%uvelN, this%vvelE, this%uvelU, this%vvelU, this%shearU, this%deltaU)
      !$acc exit data delete(this%strintxE, this%strintyN, this%divu, this%shear, this%vort)
      !$acc exit data delete(this%rdg_conv, this%rdg_shear, this%strocnxE, this%strocnyN)
      this%mapped = .false.
   end subroutine work_exit_data

   ! ------------------------------------------------------------------
   ! preparation
   ! ------------------------------------------------------------------

   subroutine ice_evp_c_prep(w, k, nx, ny, ncat, is, ie, js, je, rp, c, &
                             hm, aice, vice, vsno, aice_init, aice0, aicen, vicen, &
                             tarea, epm, npm, &
                             fcorE, fcorN, strairxE, strairyN, ss_tltxE, ss_tltyN, &
                             uocnE, vocnN, u, v, ice_err)
      !! Once-per-step preparation: masses, ice-free masks, strength, the
      !! forcing at E and N, the new-ice start velocity.  Call BEFORE the face
      !! halo of `(u, v)`; the solver then needs `ice_evp_c_average` once.
      !!
      !! Inputs at T (`hm` = the wet-T mask, real 0/1, ghost band included,
      !! `aice_init` = concentration at the START of the thermodynamic step,
      !! the one the CICE dynamics averages; `aice` the current one, which
      !! decides the dilated mask): see the argument comments.  At E and N:
      !! `epm/npm` the open masks (the corner ice mask is built from `hm`, see `prep_corner_impl`), `fcorE/fcorN` the Coriolis parameter,
      !! `strairxE` the x wind stress on the ice at E, `strairyN` the y one at
      !! N, `ss_tltxE/ss_tltyN` the sea-surface slope along the face,
      !! `uocnE` (x) / `vocnN` (y) the ocean surface velocity component ALONG
      !! the face (the new-ice start; the tangential components enter only the
      !! drag and are read by `ice_evp_c_substep`).  `u`, `v` are updated in place (zero
      !! where the face carries no ice).  Strength is `ice_compute_strength`
      !! (Rothrock 1975 / Hibler 1979, D6.1) over the dilated mask.
      !! The sea-surface tilt force is `-gravit * mass * slope` (`ssh_stress =
      !! "coupled"`).
      type(ice_evp_c_work_t), intent(inout) :: w
      type(ice_evp_c_const_t), intent(in) :: k
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      type(ice_ridge_params_t), intent(in) :: rp
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: hm(nx, ny)
         !! Wet-T mask (0/1).
      real(wp), intent(in) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice_init(nx, ny), aice0(nx, ny)
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat)
      real(wp), intent(in) :: tarea(nx, ny)
         !! T-cell area [m2].
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny)
      real(wp), intent(in) :: fcorE(nx, ny), fcorN(nx, ny)
      real(wp), intent(in) :: strairxE(nx, ny), strairyN(nx, ny), ss_tltxE(nx, ny), ss_tltyN(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnN(nx, ny)
      real(wp), intent(inout) :: u(nx + 1, ny), v(nx, ny + 1)
      integer(int32), intent(inout) :: ice_err(nx, ny)

      call prep_tmass_impl(nx, ny, is, ie, js, je, hm, vice, vsno, aice, k%rhoi, k%rhos, k%amin, k%mmin, &
                           w%tmass, w%tmphm)
      call prep_icetmask_impl(nx, ny, is, ie, js, je, hm, w%tmphm, w%iceTmask, &
                              w%stresspT, w%stressmT, w%zetax2T, w%etax2T)
      call prep_face_impl(nx, ny, is, ie, js, je, k%dt, k%amin, k%mmin, k%gravit, hm, w%tmass, aice_init, &
                          tarea, epm, npm, fcorE, fcorN, strairxE, strairyN, ss_tltxE, ss_tltyN, &
                          uocnE, vocnN, u, v, w%iceEmask, w%iceNmask, w%aiE, w%aiN, w%emassdti, &
                          w%nmassdti, w%fmE, w%fmN, w%forcexE, w%forceyN, w%uinit, w%vinit)
      call prep_corner_impl(nx, ny, is, ie, js, je, k%amin, k%mmin, hm, w%tmass, aice_init, tarea, &
                            w%iceUmask, w%stress12U)
      call ice_compute_strength(nx, ny, ncat, is - 1, ie + 1, js - 1, je + 1, rp, c, w%iceTmask, &
                                aice, vice, aice0, aicen, vicen, w%strength, ice_err)
      if (k%avg_strength) then
         call prep_strength_u_impl(nx, ny, is, ie, js, je, hm, tarea, w%strength, w%strengthU)
      end if
   end subroutine ice_evp_c_prep

   subroutine prep_tmass_impl(nx, ny, is, ie, js, je, hm, vice, vsno, aice, rhoi, rhos, amin, mmin, &
                              tmass, tmphm)
      !! `tmass = rhoi*vice + rhos*vsno` (0 on land) and the threshold mask,
      !! over the owned cells + 2 rings.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: hm(nx, ny), vice(nx, ny), vsno(nx, ny), aice(nx, ny)
      real(wp), intent(in) :: rhoi, rhos, amin, mmin
      real(wp), intent(inout) :: tmass(nx, ny)
      logical, intent(inout) :: tmphm(nx, ny)
      integer :: i, j
      real(wp) :: m

      do concurrent(j=js - 2:je + 2, i=is - 2:ie + 2) local(m)
         if (hm(i, j) > 0.5_wp) then
            m = rhoi*vice(i, j) + rhos*vsno(i, j)
         else
            m = 0.0_wp
         end if
         tmass(i, j) = m
         tmphm(i, j) = hm(i, j) > 0.5_wp .and. aice(i, j) > amin .and. m > mmin
      end do
   end subroutine prep_tmass_impl

   subroutine prep_icetmask_impl(nx, ny, is, ie, js, je, hm, tmphm, iceTmask, stresspT, stressmT, zetax2T, etax2T)
      !! 3x3 dilation of the threshold mask over the T cells the solver
      !! touches (owned + 1 ring); the stresses and viscosities are reset
      !! (the viscosities every step, the stresses off the ice mask).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: hm(nx, ny)
      logical, intent(in) :: tmphm(nx, ny)
      logical, intent(inout) :: iceTmask(nx, ny)
      real(wp), intent(inout) :: stresspT(nx, ny), stressmT(nx, ny), zetax2T(nx, ny), etax2T(nx, ny)
      integer :: i, j
      logical :: any9

      do concurrent(j=js - 1:je + 1, i=is - 1:ie + 1) local(any9)
         any9 = tmphm(i - 1, j + 1) .or. tmphm(i, j + 1) .or. tmphm(i + 1, j + 1) .or. &
                tmphm(i - 1, j) .or. tmphm(i, j) .or. tmphm(i + 1, j) .or. &
                tmphm(i - 1, j - 1) .or. tmphm(i, j - 1) .or. tmphm(i + 1, j - 1)
         iceTmask(i, j) = any9 .and. hm(i, j) > 0.5_wp
         zetax2T(i, j) = 0.0_wp
         etax2T(i, j) = 0.0_wp
         if (.not. iceTmask(i, j)) then
            stresspT(i, j) = 0.0_wp
            stressmT(i, j) = 0.0_wp
         end if
      end do
   end subroutine prep_icetmask_impl

   subroutine prep_face_impl(nx, ny, is, ie, js, je, dt, amin, mmin, gravit, hm, tmass, aice_init, &
                             tarea, epm, npm, fcorE, fcorN, strairxE, strairyN, ss_tltxE, ss_tltyN, &
                             uocnE, vocnN, u, v, iceEmask, iceNmask, aiE, aiN, emassdti, nmassdti, &
                             fmE, fmN, forcexE, forceyN, uinit, vinit)
      !! Masses, concentration and forcing at E and N (two-point masked,
      !! area-weighted means of the T cells), the ice masks, the new-ice start
      !! and the start-of-step velocity, for the owned faces `E(is:ie, js:je)`,
      !! `N(is:ie, js:je)`.  `u(i+1,j)` / `v(i,j+1)` are written by the thread
      !! that owns the face only.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: dt, amin, mmin, gravit
      real(wp), intent(in) :: hm(nx, ny), tmass(nx, ny), aice_init(nx, ny), tarea(nx, ny)
      real(wp), intent(in) :: epm(nx, ny), npm(nx, ny), fcorE(nx, ny), fcorN(nx, ny)
      real(wp), intent(in) :: strairxE(nx, ny), strairyN(nx, ny), ss_tltxE(nx, ny), ss_tltyN(nx, ny)
      real(wp), intent(in) :: uocnE(nx, ny), vocnN(nx, ny)
      real(wp), intent(inout) :: u(nx + 1, ny), v(nx, ny + 1)
      logical, intent(inout) :: iceEmask(nx, ny), iceNmask(nx, ny)
      real(wp), intent(inout) :: aiE(nx, ny), aiN(nx, ny), emassdti(nx, ny), nmassdti(nx, ny)
      real(wp), intent(inout) :: fmE(nx, ny), fmN(nx, ny), forcexE(nx, ny), forceyN(nx, ny)
      real(wp), intent(inout) :: uinit(nx, ny), vinit(nx, ny)
      integer :: i, j
      real(wp) :: w1, w2, wt, m, a
      logical :: ice, old

      do concurrent(j=js:je, i=is:ie) local(w1, w2, wt, m, a, ice, old)
         ! ---- E(i,j): between T(i,j) and T(i+1,j) ----
         w1 = hm(i, j)*tarea(i, j)
         w2 = hm(i + 1, j)*tarea(i + 1, j)
         wt = w1 + w2
         m = 0.0_wp
         a = 0.0_wp
         if (wt /= 0.0_wp) then
            m = (w1*tmass(i, j) + w2*tmass(i + 1, j))/wt
            a = (w1*aice_init(i, j) + w2*aice_init(i + 1, j))/wt
         end if
         ice = epm(i, j) > 0.5_wp .and. a > amin .and. m > mmin
         old = iceEmask(i, j)
         iceEmask(i, j) = ice
         if (ice) then
            if (.not. old) u(i + 1, j) = uocnE(i, j)
            aiE(i, j) = a
            emassdti(i, j) = m/dt
            fmE(i, j) = fcorE(i, j)*m
            forcexE(i, j) = strairxE(i, j) - gravit*m*ss_tltxE(i, j)
            uinit(i, j) = u(i + 1, j)
         else
            u(i + 1, j) = 0.0_wp
            aiE(i, j) = 0.0_wp
            emassdti(i, j) = 0.0_wp
            fmE(i, j) = 0.0_wp
            forcexE(i, j) = 0.0_wp
            uinit(i, j) = 0.0_wp
         end if
         ! ---- N(i,j): between T(i,j) and T(i,j+1) ----
         w1 = hm(i, j)*tarea(i, j)
         w2 = hm(i, j + 1)*tarea(i, j + 1)
         wt = w1 + w2
         m = 0.0_wp
         a = 0.0_wp
         if (wt /= 0.0_wp) then
            m = (w1*tmass(i, j) + w2*tmass(i, j + 1))/wt
            a = (w1*aice_init(i, j) + w2*aice_init(i, j + 1))/wt
         end if
         ice = npm(i, j) > 0.5_wp .and. a > amin .and. m > mmin
         old = iceNmask(i, j)
         iceNmask(i, j) = ice
         if (ice) then
            if (.not. old) v(i, j + 1) = vocnN(i, j)
            aiN(i, j) = a
            nmassdti(i, j) = m/dt
            fmN(i, j) = fcorN(i, j)*m
            forceyN(i, j) = strairyN(i, j) - gravit*m*ss_tltyN(i, j)
            vinit(i, j) = v(i, j + 1)
         else
            v(i, j + 1) = 0.0_wp
            aiN(i, j) = 0.0_wp
            nmassdti(i, j) = 0.0_wp
            fmN(i, j) = 0.0_wp
            forceyN(i, j) = 0.0_wp
            vinit(i, j) = 0.0_wp
         end if
      end do
   end subroutine prep_face_impl

   subroutine prep_corner_impl(nx, ny, is, ie, js, je, amin, mmin, hm, tmass, aice_init, tarea, &
                               iceUmask, stress12U)
      !! Ice mask at the corners `U(is-2:ie+1, js-2:je+1)` (four-point masked
      !! area-weighted means of T) and the corner stress reset off the mask.
      !! The corner is a stress point when AT LEAST TWO of its four T cells are
      !! wet (`sum(hm) > 1.5`, the reference scheme's C-grid corner mask), not
      !! all four: at a coastal corner the free-slip-looking extrapolation of
      !! `evp_strain_u_impl` then gives a no-slip wall (the mirrored tangential
      !! velocity makes the wall-parallel flow vanish at the wall) and a
      !! non-zero wall shear stress.  The corner VELOCITIES keep the all-wet
      !! `uvm` mask (`evp_average_impl`).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: amin, mmin
      real(wp), intent(in) :: hm(nx, ny), tmass(nx, ny), aice_init(nx, ny), tarea(nx, ny)
      logical, intent(inout) :: iceUmask(nx, ny)
      real(wp), intent(inout) :: stress12U(nx + 1, ny + 1)
      integer :: i, j
      real(wp) :: w1, w2, w3, w4, wt, m, a
      logical :: ice

      do concurrent(j=js - 2:je + 1, i=is - 2:ie + 1) local(w1, w2, w3, w4, wt, m, a, ice)
         w1 = hm(i, j)*tarea(i, j)
         w2 = hm(i + 1, j)*tarea(i + 1, j)
         w3 = hm(i, j + 1)*tarea(i, j + 1)
         w4 = hm(i + 1, j + 1)*tarea(i + 1, j + 1)
         wt = w1 + w2 + w3 + w4
         m = 0.0_wp
         a = 0.0_wp
         if (wt /= 0.0_wp) then
            m = (w1*tmass(i, j) + w2*tmass(i + 1, j) + w3*tmass(i, j + 1) + w4*tmass(i + 1, j + 1))/wt
            a = (w1*aice_init(i, j) + w2*aice_init(i + 1, j) + w3*aice_init(i, j + 1) &
                 + w4*aice_init(i + 1, j + 1))/wt
         end if
         ice = (hm(i, j) + hm(i + 1, j) + hm(i, j + 1) + hm(i + 1, j + 1)) > 1.5_wp .and. a > amin .and. m > mmin
         iceUmask(i, j) = ice
         if (.not. ice) stress12U(i + 1, j + 1) = 0.0_wp
      end do
   end subroutine prep_corner_impl

   subroutine prep_strength_u_impl(nx, ny, is, ie, js, je, hm, tarea, strength, strengthU)
      !! `strengthU = T2US(strength)` over `U(is-1:ie, js-1:je)` (masked,
      !! area-weighted four-point mean); static for the step.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: hm(nx, ny), tarea(nx, ny), strength(nx, ny)
      real(wp), intent(inout) :: strengthU(nx, ny)
      integer :: i, j
      real(wp) :: w1, w2, w3, w4, wt

      do concurrent(j=js - 1:je, i=is - 1:ie) local(w1, w2, w3, w4, wt)
         w1 = hm(i, j)*tarea(i, j)
         w2 = hm(i + 1, j)*tarea(i + 1, j)
         w3 = hm(i, j + 1)*tarea(i, j + 1)
         w4 = hm(i + 1, j + 1)*tarea(i + 1, j + 1)
         wt = w1 + w2 + w3 + w4
         strengthU(i, j) = 0.0_wp
         if (wt /= 0.0_wp) then
            strengthU(i, j) = (w1*strength(i, j) + w2*strength(i + 1, j) + w3*strength(i, j + 1) &
                               + w4*strength(i + 1, j + 1))/wt
         end if
      end do
   end subroutine prep_strength_u_impl

end module rdb_ice_evp_c_prep
