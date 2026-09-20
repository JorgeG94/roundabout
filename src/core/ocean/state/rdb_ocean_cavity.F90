module rdb_ocean_cavity
   !! Static ice-shelf cavity GEOMETRY: the prescribed ice draft
   !! `z_draft(i,j)` and the arithmetic that turns it into a water column,
   !! a grounding decision and an isostatic load.
   !!
   !! ### The datum (the whole design in three lines)
   !!
   !! With `b` the bed depth (m, positive down) and `z_draft >= 0` the
   !! ice-base depth (m, positive down):
   !!
   !! ```
   !! (D)  datum :  bt_H_ref = b - z_draft            (was: bt_H_ref = b)
   !!               ... and exactly 0 on a GROUNDED column, which has no
   !!               water column at all -- see `cavity_datum_impl`
   !! (P)  load  :  p_ice_ref = rho_ref*GRAVITY*z_draft
   !! (I)  invariant :  rho_ref*g*z_draft + (bt_H_ref - b)*rho_ref*g == 0
   !!               on every WET column (a grounded one carries no
   !!               barotropic momentum equation, so no load to count)
   !! ```
   !!
   !! (I) says the load is counted exactly ONCE: whatever the datum
   !! absorbs must not also be handed to the `eta_forcing` seam.  The
   !! two halves of "once" are:
   !!
   !!   * BAROTROPIC — the datum (D), and nothing else.  The split solver
   !!     replaces the depth mean of the layer PGF with the barotropic
   !!     solution, so the column-integrated pressure force is discarded
   !!     and `-G*grad(eta - eta_forcing)` is the only barotropic term
   !!     there is.  `p_ice_ref` therefore never joins `sf%p_surf`, out of
   !!     which `eta_ib` — and hence `eta_forcing` — is built.
   !!   * PRESSURE — `p_ice_ref` (P), assembled into
   !!     `multilayer_state_t%p_top = p_ice_ref + sf%p_surf` and read by
   !!     the FV_MOM6 `pa(nz+1)` top BC and the in-situ EOS.  Only the
   !!     load ANOMALY (total minus what the datum carries, i.e. exactly
   !!     `sf%p_surf` under the Boussinesq-isostatic convention) reaches
   !!     the seam, so an inverse-barometer run with no cavity is
   !!     bit-identical and a cavity with no `p_surf` sends the seam
   !!     nothing at all.
   !!
   !! The
   !! consequence of (D) is that the free-surface anomaly
   !! `bt_eta = sum(h_layer) - bt_H_ref` is ZERO under the shelf at rest,
   !! so every consumer of the water-column thickness
   !! `D = bt_H_ref + bt_eta` — the barotropic continuity face thickness,
   !! the Chapman phase speed, the ALE `remap_h_ref`, the wet/dry depth —
   !! is already correct with no cavity branch of its own.  That is
   !! Losch (2008) §2.1's own convention ("the 'sea-surface height' eta is
   !! the deviation from the 'reference' ice-shelf draft h"), and it is
   !! why the draft is NOT carried in `eta`.
   !!
   !! The geopotential interface stack is a DIFFERENT datum and stays
   !! absolute: the FV PGF builds `e_face(1) = -b` from the true bed and
   !! stacks upward, so the column top lands at `-z_draft + eta` on its
   !! own.  Two datums, kept apart on purpose.
   !!
   !! ### Grounding
   !!
   !! A column whose water thickness `b - z_draft` falls below
   !! `h_min_cavity` is LAND — it goes through the same
   !! `seed_wet_mask_impl` the bathymetry uses, so the static
   !! metric-zeroing land mask and the finite land-state seeding follow
   !! for free.  There is deliberately NO thin film of water under
   !! grounded ice.
   !!
   !! ### Units
   !!
   !! Every knob in `&ocean_cavity_dyn_nml` is in METRES; the setters
   !! below work in GRID coordinate units (metres on a Cartesian grid,
   !! DEGREES on spherical/curvilinear).  `cavity_fill_draft` converts at
   !! the dispatch with `topo_length_to_grid_units`, exactly as
   !! `&ocean_topo_nml slope_scale` is converted for the formula
   !! bathymetry — without it a metres position against a degrees grid
   !! puts the shelf outside the domain.
   !!
   !! Host-only by construction: every routine here runs at SETUP, before
   !! `ocean_state_enter_data`, on plain `do` loops (a `do concurrent`
   !! here would make `-stdpar=gpu` round-trip unmapped arrays through
   !! the device once per loop — the documented ic-seed trap).
   use rdb_constants, only: wp, LAND_DEPTH_THRESHOLD
   use rdb_grid, only: hgrid_t
   implicit none
   private

   public :: CAVITY_DRAFT_NONE, CAVITY_DRAFT_FLAT, CAVITY_DRAFT_LINEAR
   public :: CAVITY_DRAFT_FILE, CAVITY_DRAFT_INVALID
   public :: CAVITY_SOURCE_DRAFT, CAVITY_SOURCE_THICKNESS
   public :: CAVITY_SOURCE_IN_SITU, CAVITY_SOURCE_INVALID
   public :: parse_cavity_draft_config, parse_cavity_draft_source
   public :: parse_cavity_draft_sign
   public :: CAVITY_SIGN_DEPTH, CAVITY_SIGN_ELEVATION, CAVITY_SIGN_INVALID
   public :: cavity_draft_apply_sign
   public :: set_draft_flat, set_draft_linear
   public :: cavity_water_column_impl
   public :: cavity_apply_land_exclusion
   public :: cavity_count_grounded
   public :: cavity_fill_cover_frac
   public :: cavity_fill_p_ice_ref
   public :: cavity_draft_is_finite_nonneg
   public :: cavity_datum_impl
   public :: cavity_datum_residual
   public :: CAVITY_BOUND_INF

   real(wp), parameter :: CAVITY_BOUND_INF = 1.0e29_wp
      !! "No limit on this side" sentinel for the shelf-box corners.  A
      !! bound whose magnitude is at or above this is IGNORED, which is
      !! what the `&ocean_cavity_dyn_nml draft_x0/x1/y0/y1` defaults of
      !! +/-1e30 mean.  It is a sentinel rather than a sign convention
      !! because a NEGATIVE position is legitimate: the ghost band sits at
      !! negative global x/y, and a shelf that reaches the west wall must
      !! cover those ghost columns too — a box that stops at x = 0 would
      !! put a phantom calving front one cell outside the wall, which is
      !! the formula-bathymetry ghost trap in its cavity form.

   ! ---- draft_config enum ----
   integer, parameter :: CAVITY_DRAFT_NONE = 0
      !! `z_draft = 0` everywhere (the identity, even when enabled).
   integer, parameter :: CAVITY_DRAFT_FLAT = 1
      !! Uniform draft inside the shelf box.
   integer, parameter :: CAVITY_DRAFT_LINEAR = 2
      !! Linear-in-x draft inside the shelf box.
   integer, parameter :: CAVITY_DRAFT_FILE = 3
      !! Static 2-D NetCDF draft — deferred (fails loud at configure).
   integer, parameter :: CAVITY_DRAFT_INVALID = -1
      !! Unrecognised spelling.

   ! ---- draft_source enum ----
   integer, parameter :: CAVITY_SOURCE_DRAFT = 0
      !! The formula gives the ice-base DEPTH directly.
   integer, parameter :: CAVITY_SOURCE_THICKNESS = 1
      !! The formula gives an ice THICKNESS; `z_draft = rho_ice*h/rho_0`.
   integer, parameter :: CAVITY_SOURCE_IN_SITU = 2
      !! True isostasy against the in-situ column — deferred.
   integer, parameter :: CAVITY_SOURCE_INVALID = -1

   ! ---- draft_sign enum (`draft_config="file"` only) ----
   integer, parameter :: CAVITY_SIGN_DEPTH = 0
      !! The file variable IS the ice-base DEPTH, positive down and `>= 0`
      !! — Roundabout's own `z_draft` convention, so the values pass
      !! through unchanged.
   integer, parameter :: CAVITY_SIGN_ELEVATION = 1
      !! The file variable is the ice-base ELEVATION `z_d`, positive UP
      !! and therefore `<= 0` under a floating shelf — the convention the
      !! ISOMIP+ geometry file uses ("iceDraft ... the elevation of the
      !! ice-ocean interface (z_d)", Asay-Davis et al. 2016 Sect. 3.3).
      !! Values are NEGATED on load.
   integer, parameter :: CAVITY_SIGN_INVALID = -1
      !! Unrecognised spelling.  There is deliberately NO default that
      !! guesses from the data: a draft file whose sign is inferred from
      !! `minval < 0` would silently flip an all-zero (open-ocean)
      !! or partially-calved field, and the two conventions differ by the
      !! entire ice load.

contains

   pure function parse_cavity_draft_config(name) result(code)
      !! `&ocean_cavity_dyn_nml draft_config` -> `CAVITY_DRAFT_*`.
      !! Returns `CAVITY_DRAFT_INVALID` on an unrecognised spelling — the
      !! caller fails loud; there is no silent default.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("none", "")
         code = CAVITY_DRAFT_NONE
      case ("flat")
         code = CAVITY_DRAFT_FLAT
      case ("linear", "linear_x")
         code = CAVITY_DRAFT_LINEAR
      case ("file")
         code = CAVITY_DRAFT_FILE
      case default
         code = CAVITY_DRAFT_INVALID
      end select
   end function parse_cavity_draft_config

   pure function parse_cavity_draft_source(name) result(code)
      !! `&ocean_cavity_dyn_nml draft_source` -> `CAVITY_SOURCE_*`.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("draft", "")
         code = CAVITY_SOURCE_DRAFT
      case ("thickness")
         code = CAVITY_SOURCE_THICKNESS
      case ("in_situ")
         code = CAVITY_SOURCE_IN_SITU
      case default
         code = CAVITY_SOURCE_INVALID
      end select
   end function parse_cavity_draft_source

   pure function parse_cavity_draft_sign(name) result(code)
      !! `&ocean_cavity_dyn_nml draft_sign` -> `CAVITY_SIGN_*`.
      !! Returns `CAVITY_SIGN_INVALID` on an unrecognised spelling — the
      !! caller fails loud.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("depth", "positive_down", "")
         code = CAVITY_SIGN_DEPTH
      case ("elevation", "positive_up")
         code = CAVITY_SIGN_ELEVATION
      case default
         code = CAVITY_SIGN_INVALID
      end select
   end function parse_cavity_draft_sign

   pure subroutine cavity_draft_apply_sign(z_draft, nx, ny, sign_code)
      !! Normalise a freshly loaded draft field onto Roundabout's
      !! convention (DEPTH, positive down, `>= 0`).
      !!
      !!   * `CAVITY_SIGN_DEPTH` — identity (bit-for-bit).
      !!   * `CAVITY_SIGN_ELEVATION` — `z_draft = -z_d`.
      !!
      !! Open water in the ISOMIP+ file is `z_d = 0`, which negates to
      !! `-0.0`.  Left as the literal negation would make a downstream
      !! `z_draft /= 0` test (`cavity_apply_land_exclusion`) fire on a
      !! cell with no ice, so the zero is re-normalised explicitly.
      !! No clipping otherwise: a POSITIVE elevation (ice base above sea
      !! level, i.e. grounded ice or a file in the wrong convention)
      !! negates to a negative depth and is caught fail-loud by
      !! `cavity_draft_is_finite_nonneg` — which is the point of having
      !! an explicit knob instead of a guess.
      integer, intent(in) :: nx, ny, sign_code
      real(wp), intent(inout) :: z_draft(nx, ny)
      integer :: i, j
      if (sign_code /= CAVITY_SIGN_ELEVATION) return
      do j = 1, ny
         do i = 1, nx
            if (z_draft(i, j) == 0.0_wp) then
               z_draft(i, j) = 0.0_wp
            else
               z_draft(i, j) = -z_draft(i, j)
            end if
         end do
      end do
   end subroutine cavity_draft_apply_sign

   pure subroutine set_draft_flat(z_draft, grid, draft, x0, x1, y0, y1)
      !! Uniform draft `draft` inside the shelf box `[x0,x1] x [y0,y1]`,
      !! zero outside it (open ocean, including everything beyond the
      !! calving front at `x1`).  Any bound at or beyond
      !! `CAVITY_BOUND_INF` is ignored — the shelf is then open on that
      !! side, which is what the namelist defaults ask for.
      !!
      !! Positions are in GRID coordinate units and are compared against
      !! the CELL-CENTRE position, built the same way every formula
      !! bathymetry setter builds it:
      !! `x = (i - nghost + i_offset_global - 0.5)*dx`.  Because the
      !! offsets and the global extents come off `grid`, each rank fills
      !! its own window of ONE global shelf; on a single rank the offsets
      !! are 0 and this is the undecomposed formula.
      !!
      !! FILLS THE FULL ARRAY INCLUDING GHOSTS, by evaluating the formula
      !! at the ghost index — the same rule the bathymetry setters follow.
      !! A ghost row left at its alloc-time zero would put a phantom
      !! calving front one cell outside every wall.
      real(wp), intent(inout) :: z_draft(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: draft
         !! Draft depth (m, positive down) inside the box.
      real(wp), intent(in) :: x0, x1, y0, y1
         !! Shelf box in GRID coordinate units; a bound at or beyond
         !! `CAVITY_BOUND_INF` is ignored (open on that side).
      integer :: i, j, ng, nxt, nyt, ioff, joff
      real(wp) :: x_phys, y_phys

      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxt = size(z_draft, 1)
      nyt = size(z_draft, 2)

      do j = 1, nyt
         y_phys = (real(j - ng + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, nxt
            x_phys = (real(i - ng + ioff, wp) - 0.5_wp)*grid%dx
            if (in_shelf_box(x_phys, y_phys, x0, x1, y0, y1)) then
               z_draft(i, j) = draft
            else
               z_draft(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine set_draft_flat

   pure subroutine set_draft_linear(z_draft, grid, draft0, slope, x0, x1, y0, y1)
      !! Linear-in-x draft `z_draft = draft0 + slope*(x - x0)` inside the
      !! shelf box, clipped at 0 from below (a formula that would lift the
      !! ice base above the sea surface is open water, not negative ice),
      !! and zero outside the box.
      !!
      !! `slope` is in draft-metres per GRID unit; the dispatch converts
      !! the dimensionless namelist slope.  Ghost rule as `set_draft_flat`.
      real(wp), intent(inout) :: z_draft(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: draft0
         !! Draft at `x = x0` (m, positive down).
      real(wp), intent(in) :: slope
         !! d(draft)/dx in draft-metres per grid unit.
      real(wp), intent(in) :: x0, x1, y0, y1
      integer :: i, j, ng, nxt, nyt, ioff, joff
      real(wp) :: x_phys, y_phys, d_local

      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxt = size(z_draft, 1)
      nyt = size(z_draft, 2)

      do j = 1, nyt
         y_phys = (real(j - ng + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, nxt
            x_phys = (real(i - ng + ioff, wp) - 0.5_wp)*grid%dx
            if (in_shelf_box(x_phys, y_phys, x0, x1, y0, y1)) then
               d_local = draft0 + slope*(x_phys - x0)
               if (d_local < 0.0_wp) d_local = 0.0_wp
               z_draft(i, j) = d_local
            else
               z_draft(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine set_draft_linear

   pure function in_shelf_box(x, y, x0, x1, y0, y1) result(inside)
      !! `.true.` inside the shelf box.  Each of the four bounds is
      !! applied only when it is FINITE in the `CAVITY_BOUND_INF` sense,
      !! so the namelist defaults (+/-1e30) describe an unbounded shelf
      !! and a user who only wants a calving front in x need not describe
      !! y at all.
      real(wp), intent(in) :: x, y, x0, x1, y0, y1
      logical :: inside
      inside = .true.
      if (abs(x0) < CAVITY_BOUND_INF) inside = inside .and. (x >= x0)
      if (abs(x1) < CAVITY_BOUND_INF) inside = inside .and. (x <= x1)
      if (abs(y0) < CAVITY_BOUND_INF) inside = inside .and. (y >= y0)
      if (abs(y1) < CAVITY_BOUND_INF) inside = inside .and. (y <= y1)
   end function in_shelf_box

   pure subroutine cavity_water_column_impl(water, b, z_draft, nx, ny)
      !! `water = b - z_draft` — the reference water-column thickness the
      !! datum, the layer split and the wet-mask seed all work on.  Kept
      !! as one named routine so the three call sites cannot drift.
      !!
      !! Explicit-shape by the house rule; host-only (setup).
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: water(nx, ny)
      real(wp), intent(in) :: b(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      integer :: i, j
      do j = 1, ny
         do i = 1, nx
            water(i, j) = b(i, j) - z_draft(i, j)
         end do
      end do
   end subroutine cavity_water_column_impl

   pure subroutine cavity_apply_land_exclusion(z_draft, b, nx, ny, n_over_land)
      !! NO ICE OVER LAND (design rule R2): force `z_draft = 0` on every
      !! column the bathymetry already calls land (`b <
      !! LAND_DEPTH_THRESHOLD`), and report how many were touched so the
      !! caller can log it.
      !!
      !! Why zero rather than refuse: a formula shelf box drawn over a
      !! continent is a normal thing to write, and the physical meaning is
      !! unambiguous (there is no ocean there).  Zeroing it keeps land
      !! columns byte-identical to a cavity-free run AND keeps the datum
      !! invariant (I) exact everywhere — `bt_H_ref = b - z_draft` and
      !! `p_ice_ref = rho_ref*g*z_draft` are then consistent on every
      !! column, land included, with no `merge` anywhere downstream.
      !!
      !! GROUNDED columns (wet bed, but too little water under the ice)
      !! are NOT touched here: they keep their draft and are removed by
      !! the wet mask instead, which is what makes the metric-zeroing land
      !! mask do the rest.
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: z_draft(nx, ny)
      real(wp), intent(in) :: b(nx, ny)
      integer, intent(out) :: n_over_land
         !! Count of columns whose draft was zeroed (ghosts included).
      integer :: i, j
      n_over_land = 0
      do j = 1, ny
         do i = 1, nx
            if (b(i, j) < LAND_DEPTH_THRESHOLD .and. z_draft(i, j) /= 0.0_wp) then
               z_draft(i, j) = 0.0_wp
               n_over_land = n_over_land + 1
            end if
         end do
      end do
   end subroutine cavity_apply_land_exclusion

   pure subroutine cavity_count_grounded(b, z_draft, h_min, ng, nx_phys, ny_phys, &
                                         nx, ny, n_grounded, n_interior)
      !! Count the INTERIOR columns the cavity grounds: wet bed
      !! (`b >= LAND_DEPTH_THRESHOLD`) but water thickness
      !! `b - z_draft < h_min`.  Interior-only, because the ghost band
      !! carries extrapolated bathymetry and would bias the fraction the
      !! `grounded_max_frac` sanity bound is taken against.
      integer, intent(in) :: ng, nx_phys, ny_phys, nx, ny
      real(wp), intent(in) :: b(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      real(wp), intent(in) :: h_min
      integer, intent(out) :: n_grounded
      integer, intent(out) :: n_interior
      integer :: i, j
      n_grounded = 0
      n_interior = nx_phys*ny_phys
      do j = ng + 1, ng + ny_phys
         do i = ng + 1, ng + nx_phys
            if (b(i, j) >= LAND_DEPTH_THRESHOLD .and. &
                b(i, j) - z_draft(i, j) < h_min) then
               n_grounded = n_grounded + 1
            end if
         end do
      end do
   end subroutine cavity_count_grounded

   pure subroutine cavity_fill_cover_frac(cover_frac, z_draft, nx, ny)
      !! v1 ice-cover fraction: BINARY, `1` wherever there is any draft.
      !! An area-blended calving front (a partially covered cell) is a
      !! melt-physics decision, not a geometry one, so it is deliberately
      !! left for the slice that needs it — the field exists now only so
      !! that slice does not have to re-open the metrics lifecycle.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: cover_frac(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      integer :: i, j
      do j = 1, ny
         do i = 1, nx
            if (z_draft(i, j) > 0.0_wp) then
               cover_frac(i, j) = 1.0_wp
            else
               cover_frac(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine cavity_fill_cover_frac

   pure subroutine cavity_fill_p_ice_ref(p_ice_ref, z_draft, rho_g, nx, ny)
      !! `p_ice_ref = (rho_ref*GRAVITY) * z_draft` (Pa).
      !!
      !! `rho_g` is passed as ONE pre-multiplied scalar on purpose: the
      !! FV_MOM6 surface BC forms `rho_ref*GRAVITY*eta_geo` with the same
      !! product, so building the load this way makes
      !! `pa(nz+1) = rho_ref*g*(-z_draft) + p_ice_ref` cancel to bit-zero
      !! at rest instead of merely to a small number.  Splitting it into
      !! two multiplies would give up that exactness for nothing.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: p_ice_ref(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      real(wp), intent(in) :: rho_g
         !! `rho_ref*GRAVITY` (kg m^-2 s^-2), pre-multiplied by the caller.
      integer :: i, j
      do j = 1, ny
         do i = 1, nx
            p_ice_ref(i, j) = rho_g*z_draft(i, j)
         end do
      end do
   end subroutine cavity_fill_p_ice_ref

   pure function cavity_draft_is_finite_nonneg(z_draft, nx, ny) result(ok)
      !! Configure-time guard: every draft entry is finite and `>= 0`.
      !!
      !! The `>= 0` test is written as `.not. (z >= 0)` so a NaN FAILS it:
      !! every comparison with NaN is false, so the naive `z < 0` test
      !! would wave a NaN straight through — the same trap the
      !! NaN-laundering clamp gotcha describes, in its cheap form.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: z_draft(nx, ny)
      logical :: ok
      integer :: i, j
      ok = .true.
      do j = 1, ny
         do i = 1, nx
            if (.not. (z_draft(i, j) >= 0.0_wp)) ok = .false.
            if (.not. (z_draft(i, j) < huge(1.0_wp))) ok = .false.
         end do
      end do
   end function cavity_draft_is_finite_nonneg

   pure subroutine cavity_datum_impl(H_ref, b, z_draft, h_min, nx, ny)
      !! The BAROTROPIC DATUM: `bt_H_ref = b - z_draft` on a column that
      !! has water under the ice, and exactly `0` on one that is
      !! GROUNDED (`b - z_draft < h_min`, i.e. land by the very rule
      !! `seed_wet_mask_impl` applies).
      !!
      !! ### Why the grounded branch is not `b - z_draft`
      !!
      !! `bt_H_ref` is the reference WATER-COLUMN thickness.  A grounded
      !! column has no water column, and `b - z_draft` there is not a
      !! small thickness — it is NEGATIVE, by hundreds of metres for a
      !! real draft over a real bed.  Carrying that number had three
      !! consequences, none of them wanted:
      !!
      !!   * `bt_eta = sum h_layer - bt_H_ref` came out at `+|b -
      !!     z_draft|` on every grounded column — a phantom few-hundred-
      !!     metre free surface, masked out of the dynamics but visible
      !!     in `eta` min/max and in the `ssh` diagnostic;
      !!   * the ALE target on that column is built as
      !!     `(remap_h_ref + bt_eta)*dsig` with `remap_h_ref = total_h -
      !!     bt_eta`, so recovering the land column's `nz*H_VANISHED`
      !!     total meant CANCELLING two numbers of order the draft.  The
      !!     land thickness then jittered at `eps*z_draft` — round-off of
      !!     the wrong quantity — instead of sitting bit-stably at
      !!     `H_VANISHED`;
      !!   * it made a grounded column distinguishable from an ordinary
      !!     land column (`bt_H_ref = b`), for no dynamical reason: every
      !!     face metric on a land cell is zeroed, so nothing downstream
      !!     reads either value.
      !!
      !! Zero is the value that says "no water column" and the value that
      !! makes a grounded column arithmetically indistinguishable from
      !! the ordinary land it IS.
      !!
      !! ### Why this does not un-count the ice load
      !!
      !! Invariant (I) — the load reaching the barotropic mode exactly
      !! once, through the datum — is a statement about the BAROTROPIC
      !! MOMENTUM EQUATION, and that equation exists only on wet columns:
      !! on a land column every face metric is zero, `-G*grad(eta -
      !! eta_forcing)` is multiplied by nothing, and there is no load to
      !! count once or twice.  So (I) is asserted, by
      !! `cavity_datum_residual`, over the WET columns — which is where
      !! it is a physical statement rather than a bookkeeping one.
      !!
      !! Explicit-shape by the house rule; host-only (setup).
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: H_ref(nx, ny)
      real(wp), intent(in) :: b(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      real(wp), intent(in) :: h_min
         !! `&ocean_cavity_dyn_nml h_min_cavity` (m, validated `> 0`) —
         !! the SAME grounding cutoff the wet-mask seed applies, passed
         !! rather than re-spelled so the two decisions cannot drift.
      real(wp) :: water
      integer :: i, j
      do j = 1, ny
         do i = 1, nx
            water = b(i, j) - z_draft(i, j)
            if (water < h_min) then
               H_ref(i, j) = 0.0_wp
            else
               H_ref(i, j) = water
            end if
         end do
      end do
   end subroutine cavity_datum_impl

   pure function cavity_datum_residual(bt_H_ref, b, z_draft, h_min, nx, ny) result(resid)
      !! Max violation of the counted-once invariant (I) in METRES of
      !! reference depth, over the WET columns:
      !! `max |bt_H_ref - (b - z_draft)|` where `b - z_draft >= h_min`.
      !!
      !! (I) itself is `rho_ref*g*z_draft + (bt_H_ref - b)*rho_ref*g == 0`;
      !! dividing out the common positive factor `rho_ref*g` leaves
      !! exactly this length, which is the form worth asserting — it is
      !! scale-free and it does not fabricate a product that the code
      !! never forms.
      !!
      !! GROUNDED columns are excluded on purpose, and `cavity_datum_impl`
      !! carries the argument: they are land, they carry no barotropic
      !! momentum equation, and their datum is deliberately `0` rather
      !! than a negative water column.  Including them would assert a
      !! load-counting statement where there is no load being counted.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: bt_H_ref(nx, ny)
      real(wp), intent(in) :: b(nx, ny)
      real(wp), intent(in) :: z_draft(nx, ny)
      real(wp), intent(in) :: h_min
      real(wp) :: resid
      real(wp) :: water
      integer :: i, j
      resid = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            water = b(i, j) - z_draft(i, j)
            if (water < h_min) cycle
            resid = max(resid, abs(bt_H_ref(i, j) - water))
         end do
      end do
   end function cavity_datum_residual

end module rdb_ocean_cavity
