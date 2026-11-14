!! Ocean -> ice import of the CICE-port sea ice (roundabout-original).
module rdb_ice_import_ocean
   !! `ice_import_ocean` builds, once per ICE step, the per-cell ocean inputs the
   !! ice thermodynamics reads (MASTER_PLAN C.4, "Imports every ice step"):
   !!
   !! | field | unit | meaning |
   !! |---|---|---|
   !! | `sst` | degC | surface temperature of the FIRST LIVE layer (ice-boundary convention) |
   !! | `sss` | ppt | surface salinity of that layer (ice-boundary convention) |
   !! | `sss_ocn` | psu, ocean convention | the same salinity as the OCEAN's tracer carries it: the reference of the virtual-salt adapter |
   !! | `Tf` | degC | the ICE-side freezing point `-depressT*sss` (`ice_sea_freezing_temperature`, so `tfrz_option` is honoured) |
   !! | `frzmlt` | W m-2 | the melt-potential provider's signed field, copied (see below) |
   !!
   !! **First live layer.**  Under the bottom-up convention the surface is `k_top`
   !! (`ms%k_top`, `nz` off a rigid top) and the scan moves DOWN from it
   !! (`rdb_vl_first_live`), so a top layer that collapsed at run time does not
   !! hand the ice the water it no longer holds.  Under `cons_abs` the ocean
   !! carries Conservative Temperature and Reference Salinity; `sst`/`sss` are
   !! converted to potential temperature and Practical Salinity exactly as the
   !! basal-flux seam does (`pt_from_ct`, `sp_from_sr`), `sss_ocn` is NOT.
   !!
   !! **Two freezing points, kept apart.**  `Tf` is the ice's linear liquidus.
   !! The ocean's frazil clamp and the melt potential use the EOS liquidus; the
   !! import does not reconcile them (the ice side clamps).
   !!
   !! **`frzmlt` cadence and dt contract.**  `frzmlt` is a RATE only over the
   !! window the provider divided by: `frzmlt = bank/dt` with `dt` the interval
   !! the frazil bank accumulated over.  The `scheme = "cice"` ice steps EVERY
   !! outer step and `add_new_ice` zeroes the bank where `frzmlt > 0`, so the
   !! provider (`ice_melt_potential`) must be called EVERY ice step, immediately
   !! before this import, with `dt` = the ice step; a provider left at thermo
   !! cadence would hand `add_new_ice` the same bank again on every ice step of
   !! the window (the ice formed twice).  The import only copies: it does not
   !! know the cadence, so the call order is the wiring PR's obligation
   !! (provider -> import -> ice step).  The provider's ghosts are not filled and
   !! neither are the import's: every field here is valid on OWNED cells only.
   !!
   !! **Melt-potential cap over a thermo window (interim).**  The ice steps every outer
   !! step but the ocean is debited only at the thermo hand-off, so with N ice steps per
   !! window the same upper-ocean heat would be offered N times (up to N x the available
   !! heat could be drawn).  When the averager's OPEN window is passed (`avg`, `dt`) the
   !! melt side is capped: `MP_eff = max(0, MP + t_acc*min(acc_fhocn, 0))` with
   !! `MP = -frzmlt*dt` the offered potential [J m-2] (`dt` = the provider's step = the
   !! ice step), `t_acc` the seconds already accumulated and `acc_fhocn` the window-mean
   !! heat flux so far, so heat the ice already took this window is no longer on offer;
   !! `frzmlt_eff = -MP_eff/dt`.  The freezing side (`frzmlt > 0`, the bank) is never
   !! touched.  At the first step of a window `t_acc = 0`, the cap is skipped and the
   !! result is BITWISE the copied `frzmlt`.  The plan is silent on this; the reference
   !! coupled model has the same exposure.  A non-finite `acc_fhocn` poisons the cell.
   !!
   !! **Engine contract (scheme select).**  The engine must `select case` on the scheme:
   !! the SIS2 branch keeps its thermo-cadence provider call and `ice_frazil_uptake`, the
   !! CICE branch calls provider -> import every ice step.  The two branches never run
   !! together (two spenders of the one frazil bank, two providers of `frzmlt`).  The
   !! frazil half is consistent as built: the bank fills every outer step, the provider
   !! reads it every ice step with the ice `dt`, and the new-ice driver zeroes it.
   !!
   !! **Land, dry, vanished.**  A land cell, or a wet cell with no live layer,
   !! gets `sst = sss = sss_ocn = Tf = frzmlt = 0` (no drive: `sst - Tf = 0`).
   !!
   !! **Non-finite input.**  A non-finite `h` at the surface start layer, a
   !! non-finite `T`/`S` of the chosen layer, a non-finite `frzmlt` on a wet cell,
   !! or an unknown `tfrz_option` sets `ICE_ERR_NONFINITE` in the cell's ice error
   !! word and writes NaN to every output of that cell (never a plausible
   !! value).
   !!
   !! **Not imported here.**  The ocean surface velocity and the sea-surface tilt
   !! the dynamics needs live on FACES and are consumed by the EVP prep (`uocnE/N`
   !! in `rdb_ice_evp_c`), which the dynamics coupling PR owns; the thermodynamics
   !! needs only the five fields above.
   !!
   !! **Lifecycle.**  `ice_ocean_in_t%init` allocates (zero); `enter_data` maps the
   !! COMPONENT arrays with `copyin`; nothing here is restart state (recomputed
   !! at the start of every ice step).  The slot exists only under the new scheme;
   !! under `scheme = "sis2"` nothing is allocated, mapped or counted.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp, H_VANISHED, NZ_STACK_MAX
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, TS_CONS_ABS
   use rdb_teos10_convert, only: pt_from_ct, sp_from_sr
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ICE_TFRZ_LINEAR_SALT, ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_avg, only: ice_flux_avg_t
   implicit none
   private

   public :: ice_ocean_in_t, ice_import_ocean, ice_import_ocean_impl, ice_melt_cap
   public :: ICE_IMPORT_OK, ICE_IMPORT_ERR_SHAPE

   integer, parameter :: ICE_IMPORT_OK = 0
      !! `ice_import_ocean` ran.
   integer, parameter :: ICE_IMPORT_ERR_SHAPE = 1
      !! Slots disagree on shape, T/S unregistered, or not initialised; nothing ran.

   type :: ice_ocean_in_t
      !! Per-cell ocean inputs of the ice (see the module docstring).
      logical :: is_init = .false.
         !! True between `init` and `destroy`.
      logical :: mapped = .false.
         !! True between `enter_data` and `exit_data`.
      integer :: nx = 0
         !! Total extent in x including ghosts.
      integer :: ny = 0
         !! Total extent in y including ghosts.
      integer :: nghost = 0
         !! Ghost width.
      real(wp), allocatable :: sst(:, :)
         !! [degC] surface temperature, ice-boundary convention.
      real(wp), allocatable :: sss(:, :)
         !! [ppt] surface salinity, ice-boundary convention.
      real(wp), allocatable :: sss_ocn(:, :)
         !! [ocean salinity units] surface salinity as the ocean tracer carries it.
      real(wp), allocatable :: Tf(:, :)
         !! [degC] ice-side freezing point of `sss`.
      real(wp), allocatable :: frzmlt(:, :)
         !! [W m-2] signed freezing (+) / melt (-) potential, copied from the provider.
   contains
      procedure, non_overridable :: init => in_init
      procedure, non_overridable :: destroy => in_destroy
      procedure, non_overridable :: enter_data => in_enter_data
      procedure, non_overridable :: exit_data => in_exit_data
      procedure, non_overridable :: bytes => in_bytes
   end type ice_ocean_in_t

contains

   subroutine in_init(this, grid)
      !! Allocate the five arrays at the ghosted extent, zeroed.
      class(ice_ocean_in_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      if (this%is_init) call this%destroy()
      this%nx = grid%nx_total
      this%ny = grid%ny_total
      this%nghost = grid%nghost
      allocate (this%sst(this%nx, this%ny), source=0.0_wp)
      allocate (this%sss(this%nx, this%ny), source=0.0_wp)
      allocate (this%sss_ocn(this%nx, this%ny), source=0.0_wp)
      allocate (this%Tf(this%nx, this%ny), source=0.0_wp)
      allocate (this%frzmlt(this%nx, this%ny), source=0.0_wp)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine in_init

   subroutine in_destroy(this)
      !! Reverse of `init` (unmaps first).
      class(ice_ocean_in_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%sst)) deallocate (this%sst)
      if (allocated(this%sss)) deallocate (this%sss)
      if (allocated(this%sss_ocn)) deallocate (this%sss_ocn)
      if (allocated(this%Tf)) deallocate (this%Tf)
      if (allocated(this%frzmlt)) deallocate (this%frzmlt)
      this%nx = 0
      this%ny = 0
      this%nghost = 0
   end subroutine in_destroy

   subroutine in_enter_data(this)
      !! `copyin` the component arrays (host zeros).
      class(ice_ocean_in_t), intent(inout) :: this
      select type (this)
      type is (ice_ocean_in_t)
         call in_enter_data_impl(this)
      end select
   end subroutine in_enter_data

   subroutine in_enter_data_impl(this)
      type(ice_ocean_in_t), intent(inout) :: this
      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%sst, this%sss, this%sss_ocn, this%Tf, this%frzmlt)
      this%mapped = .true.
   end subroutine in_enter_data_impl

   subroutine in_exit_data(this)
      !! `delete` the component arrays (per-step scratch).
      class(ice_ocean_in_t), intent(inout) :: this
      select type (this)
      type is (ice_ocean_in_t)
         call in_exit_data_impl(this)
      end select
   end subroutine in_exit_data

   subroutine in_exit_data_impl(this)
      type(ice_ocean_in_t), intent(inout) :: this
      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%sst, this%sss, this%sss_ocn, this%Tf, this%frzmlt)
      this%mapped = .false.
   end subroutine in_exit_data_impl

   pure function in_bytes(this) result(nbytes)
      !! Counted allocatable footprint (one `arr_bytes` term per array).
      class(ice_ocean_in_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = 0_int64 &
               + arr_bytes(this%sst) &
               + arr_bytes(this%sss) &
               + arr_bytes(this%sss_ocn) &
               + arr_bytes(this%Tf) &
               + arr_bytes(this%frzmlt) &
               + 0_int64
   end function in_bytes

   subroutine ice_import_ocean(grid, eos, ms, frzmlt_src, cice, oin, ierr, avg, dt)
      !! Outer shim: dereference the tracer registry on the host and forward bare
      !! explicit-shape arrays to `ice_import_ocean_impl`.  `frzmlt_src(nx, ny)` is
      !! the provider's signed field (`ocean_sea_ice_t%frzmlt`, see the dt
      !! contract above); `cice%p` supplies `depressT`, `Tocnfrz`, `tfrz_option`;
      !! `cice%ice_err` receives the non-finite bit.  `ierr = ICE_IMPORT_ERR_SHAPE`
      !! (nothing ran) when a slot is uninitialised, the shapes disagree or T/S
      !! are unregistered.  All slots must be mapped on an offload build.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: frzmlt_src(:, :)
      type(ice_cice_state_t), intent(inout) :: cice
      type(ice_ocean_in_t), intent(inout) :: oin
      integer, intent(out) :: ierr
      type(ice_flux_avg_t), intent(in), optional :: avg
         !! The OPEN averaging window: enables the melt cap (module docstring).
      real(wp), intent(in), optional :: dt
         !! The ice step [s] the provider divided by; required with `avg`.

      integer :: idx_T, idx_S
      logical :: use_cap
      real(wp) :: t_acc, dt_c

      ierr = ICE_IMPORT_OK
      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (.not. (oin%is_init .and. cice%is_init) .or. idx_T <= 0 .or. idx_S <= 0 .or. &
          oin%nx /= grid%nx_total .or. oin%ny /= grid%ny_total .or. &
          cice%nx /= grid%nx_total .or. cice%ny /= grid%ny_total .or. &
          size(frzmlt_src, 1) /= grid%nx_total .or. size(frzmlt_src, 2) /= grid%ny_total) then
         ierr = ICE_IMPORT_ERR_SHAPE
         return
      end if
      use_cap = present(avg) .and. present(dt)
      t_acc = 0.0_wp
      dt_c = 1.0_wp
      if (use_cap) then
         if (.not. avg%is_init .or. avg%nx /= grid%nx_total .or. avg%ny /= grid%ny_total .or. &
             .not. (dt > 0.0_wp)) then
            ierr = ICE_IMPORT_ERR_SHAPE
            return
         end if
         t_acc = avg%t_acc
         dt_c = dt
      end if
      if (use_cap) then
         call ice_import_ocean_impl(ms%tracers(idx_T)%c, ms%tracers(idx_S)%c, ms%h_layer, &
                                    ms%wet_mask, ms%k_top, frzmlt_src, eos%ts_convention == TS_CONS_ABS, &
                                    cice%p%c%depressT, cice%p%c%Tocnfrz, cice%p%thermo%tfrz_option, &
                                    ieee_value(1.0_wp, ieee_quiet_nan), grid%nghost, ms%nz_ml, &
                                    grid%nx_total, grid%ny_total, &
                                    oin%sst, oin%sss, oin%sss_ocn, oin%Tf, oin%frzmlt, cice%ice_err, &
                                    .true., avg%acc_fhocn, t_acc, dt_c)
         return
      end if
      call ice_import_ocean_impl(ms%tracers(idx_T)%c, ms%tracers(idx_S)%c, ms%h_layer, &
                                 ms%wet_mask, ms%k_top, frzmlt_src, eos%ts_convention == TS_CONS_ABS, &
                                 cice%p%c%depressT, cice%p%c%Tocnfrz, cice%p%thermo%tfrz_option, &
                                 ieee_value(1.0_wp, ieee_quiet_nan), grid%nghost, ms%nz_ml, &
                                 grid%nx_total, grid%ny_total, &
                                 oin%sst, oin%sss, oin%sss_ocn, oin%Tf, oin%frzmlt, cice%ice_err, &
                                 .false., frzmlt_src, 0.0_wp, 1.0_wp)
   end subroutine ice_import_ocean

   pure subroutine ice_import_ocean_impl(c_T, c_S, h_layer, wet_mask, k_top, frzmlt_src, cons, &
                                         depressT, Tocnfrz, tfrz_option, nan, nghost, nz, nx, ny, &
                                         sst, sss, sss_ocn, Tf, frzmlt, ice_err, &
                                         use_cap, acc_fhocn, t_acc, dt)
      !! Device kernel over the OWNED cells, one thread per column, ONE launch;
      !! ghosts untouched.  Every output of a cell is written (the five arrays
      !! need no prior zeroing).
      integer, intent(in) :: nghost, nz, nx, ny
      real(wp), intent(in) :: c_T(nx, ny, nz)
      real(wp), intent(in) :: c_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      integer, intent(in) :: k_top(nx, ny)
      real(wp), intent(in) :: frzmlt_src(nx, ny)
      logical, intent(in) :: cons
      real(wp), intent(in) :: depressT, Tocnfrz, nan
      integer, intent(in) :: tfrz_option
      real(wp), intent(inout) :: sst(nx, ny), sss(nx, ny), sss_ocn(nx, ny), Tf(nx, ny), frzmlt(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      logical, intent(in) :: use_cap
         !! Apply `ice_melt_cap` with the open window `(acc_fhocn, t_acc)` and step `dt`.
      real(wp), intent(in) :: acc_fhocn(nx, ny)
      real(wp), intent(in) :: t_acc, dt

      integer :: i, j, kt, k, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h, t_k, s_k, t_ice, s_ice, tf_k
      logical :: bad, live

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(kt, k, h, t_k, s_k, t_ice, s_ice, tf_k, bad, live)
         t_ice = 0.0_wp
         s_ice = 0.0_wp
         s_k = 0.0_wp
         tf_k = 0.0_wp
         bad = .false.
         live = .false.
         sst(i, j) = 0.0_wp
         sss(i, j) = 0.0_wp
         sss_ocn(i, j) = 0.0_wp
         Tf(i, j) = 0.0_wp
         frzmlt(i, j) = 0.0_wp
         if (wet_mask(i, j) > 0.5_wp) then
            kt = k_top(i, j)
            if (.not. ieee_is_finite(h_layer(i, j, kt))) bad = .true.
            k = rdb_vl_first_live(nx, ny, nz, h_layer, i, j, kt)
            h = h_layer(i, j, k)
            if (.not. ieee_is_finite(h)) bad = .true.
            if (.not. ieee_is_finite(frzmlt_src(i, j))) bad = .true.
            live = rdb_vl_is_live(h)
            if (live .and. .not. bad) then
               t_k = c_T(i, j, k)
               s_k = c_S(i, j, k)
               if (.not. (ieee_is_finite(t_k) .and. ieee_is_finite(s_k))) then
                  bad = .true.
               else
                  if (cons) then
                     t_ice = pt_from_ct(s_k, t_k)
                     s_ice = sp_from_sr(s_k)
                  else
                     t_ice = t_k
                     s_ice = s_k
                  end if
                  ! Same three options as `ice_sea_freezing_temperature`, written inline, so
                  ! the kernel calls nothing but the vanished-layer helpers.
                  select case (tfrz_option)
                  case (ICE_TFRZ_LINEAR_SALT)
                     tf_k = -depressT*s_ice
                  case (ICE_TFRZ_CONSTANT)
                     tf_k = Tocnfrz
                  case (ICE_TFRZ_MINUS1P8)
                     tf_k = -1.8_wp
                  case default
                     tf_k = 0.0_wp
                     bad = .true.
                  end select
               end if
            end if
            if (bad) then
               sst(i, j) = nan
               sss(i, j) = nan
               sss_ocn(i, j) = nan
               Tf(i, j) = nan
               frzmlt(i, j) = nan
               ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
            else if (live) then
               sst(i, j) = t_ice
               sss(i, j) = s_ice
               sss_ocn(i, j) = s_k
               Tf(i, j) = tf_k
               frzmlt(i, j) = frzmlt_src(i, j)
               if (use_cap) frzmlt(i, j) = ice_melt_cap(frzmlt_src(i, j), acc_fhocn(i, j), t_acc, dt, nan)
            end if
         end if
      end do
   end subroutine ice_import_ocean_impl

   pure function ice_melt_cap(frzmlt, acc_fhocn, t_acc, dt, nan) result(eff)
      !! `frzmlt_eff = -max(0, MP + t_acc*min(acc_fhocn, 0))/dt` with `MP = -frzmlt*dt`
      !! on the melt side (`frzmlt < 0`); the freezing side and a window that has taken
      !! no heat (`t_acc*min(acc_fhocn, 0) == 0`, in particular the first step) return
      !! `frzmlt` itself, bit for bit.  A non-finite `acc_fhocn` returns `nan`.
      !$acc routine seq
      real(wp), intent(in) :: frzmlt, acc_fhocn, t_acc, dt, nan
      real(wp) :: eff
      real(wp) :: used, mp

      eff = frzmlt
      if (.not. ieee_is_finite(acc_fhocn)) then
         eff = nan
         return
      end if
      used = t_acc*min(acc_fhocn, 0.0_wp)
      if (frzmlt < 0.0_wp .and. used /= 0.0_wp) then
         mp = -frzmlt*dt
         eff = -max(0.0_wp, mp + used)/dt
      end if
   end function ice_melt_cap

#include "rdb_vanished_layer.inc"

end module rdb_ice_import_ocean
