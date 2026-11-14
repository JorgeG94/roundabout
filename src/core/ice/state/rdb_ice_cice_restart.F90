!! Start-up read path of the CICE-port ice state: v2 read, SIS2 convert, or cold start.
module rdb_ice_cice_restart
   !! The three ways a `scheme = "cice"` run can start, decided from which
   !! restart variables the file really had (`restart_registry_t%entry_found`),
   !! never from "which values look non-zero":
   !!
   !!  1. EVERY v2 tag found (`ice_state_version == ICE_STATE_VERSION`, the
   !!     structure scalars `ice_ncat/nilyr/nslyr` equal to the slot's, and all
   !!     arrays) -> already read into the slot (`ICE_RESUME_V2`).
   !!  2. no v2 tag at all and EVERY SIS2 tag found (`ice_part_size`, `ice_m_ice`,
   !!     `ice_m_snow`, every `ice_enth_ice_k<k>`/`ice_sal_ice_k<k>` for the
   !!     `legacy` layer count, `ice_enth_snow_k1`, and NO layer beyond that count)
   !!     -> the SIS2 state read into `ice_legacy_t` is converted ONCE by
   !!     `ice_restart_convert_sis2` (`ICE_RESUME_CONVERTED`), mass (exact, with
   !!     the dropped mass) and heat (relative change) logged.  One-way: every
   !!     later write carries version 2 only (`st%state_version` is reset).
   !!     Decision M-D1-c: this module is removed one release after the flip.
   !!  3. no ice tag of either kind -> cold start: the empty-ice seed
   !!     (`ICE_RESUME_COLD`), logged.
   !!  4. anything else -- a PARTIAL set of either kind, a version other than
   !!     `ICE_STATE_VERSION`, a structure mismatch, legacy extents (`nx`, `ny`,
   !!     `nghost`) different from the slot's -- is REFUSED (`ICE_RESUME_FAILED`,
   !!     reason logged); a half-read file is never silently accepted.
   !!
   !! HOST code on host arrays before `enter_data` (it runs once, before the
   !! device map; the copyin then carries the result up).  Usage by the engine:
   !! `legacy%init(...)`; register `st%register_restart(reg)` in the RUN registry
   !! and `legacy%register_restart(reg_legacy)` in a SEPARATE, short-lived,
   !! read-only registry; read the file into both; `ice_cice_restart_resume`;
   !! then `legacy%destroy()` and drop `reg_legacy` (the legacy tags must never
   !! enter the run registry: their pointers dangle after `destroy` and every
   !! later checkpoint would carry stale SIS2 tags); `enter_data`; and
   !! `ice_update_state(.true.)` (the aggregates are not restart state).
   !!
   !! Wired into the engine's restart read (`resume_cice_ice`): the SIS2 counts are
   !! `&ocean_ice_nml ncat` / `nk_ice` (a SIS2 file must have been written with the
   !! same ones), and `Tf = -depressT*S_ref` of the EOS reference salinity.
   use rdb_constants, only: wp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_restart_convert, only: ice_restart_convert_sis2, ice_convert_report_t, ICE_CONV_OK
   use rdb_ice_seed, only: ICE_SEED_OK
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_state_seed_uniform, ICE_STATE_VERSION
   implicit none
   private

   public :: ice_legacy_t, ice_cice_restart_resume
   public :: ICE_RESUME_V2, ICE_RESUME_CONVERTED, ICE_RESUME_COLD, ICE_RESUME_FAILED

   integer, parameter :: ICE_RESUME_V2 = 0
      !! The v2 tags were read; the slot is already the restart state.
   integer, parameter :: ICE_RESUME_CONVERTED = 1
      !! A SIS2 checkpoint was converted.
   integer, parameter :: ICE_RESUME_COLD = 2
      !! No ice state in the file: empty cold start.
   integer, parameter :: ICE_RESUME_FAILED = 3
      !! Refused (partial or mismatched file, bad version, extent mismatch,
      !! non-finite or negative input).  The caller must abort.

   type :: ice_legacy_t
      !! Short-lived host-only copy of the SIS2 restart variables (same tags and
      !! array conventions as `rdb_ice_state`; `scheme = "cice"` allocates no SIS2
      !! slot, so the file is read into this).  Defaults mirror the SIS2 init
      !! (open water, zero ice, bulk salinity 4 ppt) so an absent variable is inert.
      !! Register it in a registry of its OWN, never the run registry.
      integer :: nx = 0, ny = 0, ncat = 0, nk = 0, nghost = 0
      real(wp), allocatable :: part_size(:, :, :)    !! (nx,ny,0:ncat)
      real(wp), allocatable :: m_ice(:, :, :)        !! (nx,ny,ncat) kg/m2
      real(wp), allocatable :: m_snow(:, :, :)       !! (nx,ny,ncat) kg/m2
      real(wp), allocatable :: enth_ice(:, :, :, :)  !! (nx,ny,ncat,nk) J/kg, bottom-up
      real(wp), allocatable :: sal_ice(:, :, :, :)   !! (nx,ny,ncat,nk) ppt, bottom-up
      real(wp), allocatable :: enth_snow(:, :, :, :)  !! (nx,ny,ncat,1) J/kg
      real(wp), allocatable :: probe(:, :, :)        !! (nx,ny,ncat) sink for tags BEYOND `nk`
   contains
      procedure :: init => legacy_init
      procedure :: destroy => legacy_destroy
      procedure :: register_restart => legacy_register_restart
   end type ice_legacy_t

contains

   subroutine legacy_init(this, nx, ny, ncat, nk, nghost)
      !! Allocate the legacy arrays for the SIS2 structure `ncat` x `nk` on an
      !! `nx` x `ny` tile (ghosts included, `nghost` wide).
      class(ice_legacy_t), intent(inout) :: this
      integer, intent(in) :: nx, ny, ncat, nk, nghost

      call this%destroy()
      this%nx = nx
      this%ny = ny
      this%ncat = ncat
      this%nk = nk
      this%nghost = nghost
      allocate (this%part_size(nx, ny, 0:ncat), source=0.0_wp)
      this%part_size(:, :, 0) = 1.0_wp
      allocate (this%m_ice(nx, ny, ncat), source=0.0_wp)
      allocate (this%m_snow(nx, ny, ncat), source=0.0_wp)
      allocate (this%enth_ice(nx, ny, ncat, nk), source=0.0_wp)
      allocate (this%sal_ice(nx, ny, ncat, nk), source=4.0_wp)
      allocate (this%enth_snow(nx, ny, ncat, 1), source=0.0_wp)
      allocate (this%probe(nx, ny, ncat), source=0.0_wp)
   end subroutine legacy_init

   subroutine legacy_destroy(this)
      class(ice_legacy_t), intent(inout) :: this

      if (allocated(this%part_size)) deallocate (this%part_size)
      if (allocated(this%m_ice)) deallocate (this%m_ice)
      if (allocated(this%m_snow)) deallocate (this%m_snow)
      if (allocated(this%enth_ice)) deallocate (this%enth_ice)
      if (allocated(this%sal_ice)) deallocate (this%sal_ice)
      if (allocated(this%enth_snow)) deallocate (this%enth_snow)
      if (allocated(this%probe)) deallocate (this%probe)
      this%nx = 0
      this%ny = 0
      this%ncat = 0
      this%nk = 0
      this%nghost = 0
   end subroutine legacy_destroy

   subroutine legacy_register_restart(this, reg)
      !! Register the SIS2 tags (`ice_part_size`, `ice_m_ice`, `ice_m_snow`,
      !! `ice_enth_ice_k<n>`, `ice_sal_ice_k<n>`, `ice_enth_snow_k1`), all
      !! `optional`, `device_mapped = .false.`, in the caller's LEGACY registry,
      !! plus the PROBE tags `ice_enth_ice_k<nk+1>` / `ice_sal_ice_k<nk+1>` whose
      !! only job is to reveal a file with more layers than `nk`.
      class(ice_legacy_t), target, intent(in) :: this
      type(restart_registry_t), intent(inout) :: reg
      integer :: k

      if (.not. allocated(this%part_size)) return
      call reg%register_3d("ice_part_size", this%part_size, 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
      call reg%register_3d("ice_m_ice", this%m_ice, 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
      call reg%register_3d("ice_m_snow", this%m_snow, 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
      do k = 1, this%nk
         call reg%register_3d("ice_enth_ice_k"//to_string(k), this%enth_ice(:, :, :, k), 0, &
                              this%nx, this%ny, optional=.true., device_mapped=.false.)
         call reg%register_3d("ice_sal_ice_k"//to_string(k), this%sal_ice(:, :, :, k), 0, &
                              this%nx, this%ny, optional=.true., device_mapped=.false.)
      end do
      call reg%register_3d("ice_enth_snow_k1", this%enth_snow(:, :, :, 1), 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
      call reg%register_3d("ice_enth_ice_k"//to_string(this%nk + 1), this%probe, 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
      call reg%register_3d("ice_sal_ice_k"//to_string(this%nk + 1), this%probe, 0, this%nx, this%ny, &
                           optional=.true., device_mapped=.false.)
   end subroutine legacy_register_restart

   subroutine count_v2(st, reg, nfound, ntotal)
      !! How many of the slot's v2 tags (version, structure scalars, arrays) the
      !! file had.
      type(ice_cice_state_t), intent(in) :: st
      type(restart_registry_t), intent(in) :: reg
      integer, intent(out) :: nfound, ntotal
      integer :: k

      nfound = 0
      ntotal = 0
      call tally("ice_state_version")
      call tally("ice_ncat")
      call tally("ice_nilyr")
      call tally("ice_nslyr")
      call tally("ice_aicen")
      call tally("ice_vicen")
      call tally("ice_vsnon")
      call tally("ice_Tsfcn")
      do k = 1, st%nilyr
         call tally("ice_qicen_k"//to_string(k))
         call tally("ice_sicen_k"//to_string(k))
      end do
      do k = 1, st%nslyr
         call tally("ice_qsnon_k"//to_string(k))
      end do
   contains
      subroutine tally(tag)
         character(len=*), intent(in) :: tag
         ntotal = ntotal + 1
         if (reg%entry_found(tag)) nfound = nfound + 1
      end subroutine tally
   end subroutine count_v2

   subroutine count_sis2(legacy, reg_legacy, nfound, ntotal, beyond)
      !! How many SIS2 tags the file had, and whether it has a layer beyond `nk`.
      type(ice_legacy_t), intent(in) :: legacy
      type(restart_registry_t), intent(in) :: reg_legacy
      integer, intent(out) :: nfound, ntotal
      logical, intent(out) :: beyond
      integer :: k

      nfound = 0
      ntotal = 0
      call tally("ice_part_size")
      call tally("ice_m_ice")
      call tally("ice_m_snow")
      call tally("ice_enth_snow_k1")
      do k = 1, legacy%nk
         call tally("ice_enth_ice_k"//to_string(k))
         call tally("ice_sal_ice_k"//to_string(k))
      end do
      beyond = reg_legacy%entry_found("ice_enth_ice_k"//to_string(legacy%nk + 1)) .or. &
               reg_legacy%entry_found("ice_sal_ice_k"//to_string(legacy%nk + 1))
   contains
      subroutine tally(tag)
         character(len=*), intent(in) :: tag
         ntotal = ntotal + 1
         if (reg_legacy%entry_found(tag)) nfound = nfound + 1
      end subroutine tally
   end subroutine count_sis2

   subroutine ice_cice_restart_resume(st, legacy, reg, reg_legacy, Tf, status, rep)
      !! After `reg` (the run registry holding `st%register_restart`) and
      !! `reg_legacy` (a separate registry holding `legacy%register_restart`) have
      !! been READ from the same file, decide which start applies (module
      !! docstring), do it, and log it.  `Tf` [degC] is the ocean freezing point
      !! for empty categories / cold start.  `status` is `ICE_RESUME_*`; `rep` is
      !! filled on a conversion.  Host arrays, before `enter_data`.  On a refusal
      !! raised BEFORE the conversion (structure, extents, partial SIS2 set) the
      !! state arrays are untouched.
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_legacy_t), intent(in) :: legacy
      type(restart_registry_t), intent(in) :: reg, reg_legacy
      real(wp), intent(in) :: Tf
      integer, intent(out) :: status
      type(ice_convert_report_t), intent(out) :: rep

      integer :: ierr, nf2, nt2, nf1, nt1
      logical :: beyond, have_legacy

      rep = ice_convert_report_t()
      status = ICE_RESUME_FAILED
      if (.not. st%is_init) return
      call count_v2(st, reg, nf2, nt2)
      have_legacy = allocated(legacy%part_size)
      nf1 = 0
      nt1 = 0
      beyond = .false.
      if (have_legacy) call count_sis2(legacy, reg_legacy, nf1, nt1, beyond)

      if (nf2 > 0) then
         if (nf2 /= nt2) then
            call logger%error("ice restart: PARTIAL CICE-port tag set ("//to_string(nf2)//" of "// &
                              to_string(nt2)//" found); refusing")
            return
         end if
         if (st%state_version /= ICE_STATE_VERSION) then
            call logger%error("ice restart: ice_state_version "//to_string(st%state_version)// &
                              " is not the supported "//to_string(ICE_STATE_VERSION)//"; refusing")
            return
         end if
         if (st%ncat_tag /= real(st%ncat, wp) .or. st%nilyr_tag /= real(st%nilyr, wp) .or. &
             st%nslyr_tag /= real(st%nslyr, wp)) then
            call logger%error("ice restart: file structure (ncat, nilyr, nslyr) = "// &
                              to_string(st%ncat_tag)//", "//to_string(st%nilyr_tag)//", "// &
                              to_string(st%nslyr_tag)//" differs from the run's "//to_string(st%ncat)//", "// &
                              to_string(st%nilyr)//", "//to_string(st%nslyr)//"; refusing")
            return
         end if
         status = ICE_RESUME_V2
         call logger%info("ice restart: CICE-port state (version "//to_string(int(st%state_version))//") read")
         return
      end if

      if (nf1 > 0) then
         if (nf1 /= nt1) then
            call logger%error("ice restart: PARTIAL SIS2 tag set ("//to_string(nf1)//" of "// &
                              to_string(nt1)//" found; masses, a layer or part_size missing); refusing")
            return
         end if
         if (beyond) then
            call logger%error("ice restart: the SIS2 file has more than "//to_string(legacy%nk)// &
                              " ice layers; refusing")
            return
         end if
         if (legacy%nx /= st%nx .or. legacy%ny /= st%ny .or. legacy%nghost /= st%nghost) then
            call logger%error("ice restart: legacy extents (nx, ny, nghost) = "//to_string(legacy%nx)//", "// &
                              to_string(legacy%ny)//", "//to_string(legacy%nghost)//" differ from the slot's "// &
                              to_string(st%nx)//", "//to_string(st%ny)//", "//to_string(st%nghost)//"; refusing")
            return
         end if
         call ice_restart_convert_sis2(st%nx, st%ny, legacy%ncat, legacy%nk, st%ncat, st%nilyr, st%nslyr, &
                                       st%p, st%hin_max, Tf, legacy%part_size, legacy%m_ice, legacy%m_snow, &
                                       legacy%enth_ice, legacy%sal_ice, legacy%enth_snow, &
                                       st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon, &
                                       rep, ierr)
         if (ierr /= ICE_CONV_OK) then
            call logger%error("ice restart: SIS2 -> CICE conversion refused, code "//to_string(ierr))
            return
         end if
         st%state_version = ICE_STATE_VERSION
         status = ICE_RESUME_CONVERTED
         call logger%info("ice restart: converted a SIS2 checkpoint (one-way). ice mass in/out/dropped "// &
                          "[kg/m2 summed] "//to_string(rep%mass_ice_in)//" / "//to_string(rep%mass_ice_out)// &
                          " / "//to_string(rep%mass_ice_dropped)//", snow mass in/out/dropped "// &
                          to_string(rep%mass_snow_in)//" / "//to_string(rep%mass_snow_out)//" / "// &
                          to_string(rep%mass_snow_dropped)//"; heat relative change "// &
                          to_string(rep%heat_rel_change)//"; cells rescaled "//to_string(rep%n_rescaled)// &
                          ", categories dropped "//to_string(rep%n_dropped)//", layers T-capped "// &
                          to_string(rep%n_t_capped))
         return
      end if
      call ice_state_seed_uniform(st, 0.0_wp, 0.0_wp, 0.0_wp, min(Tf, st%p%c%Timelt), Tf, ierr)
      if (ierr /= ICE_SEED_OK) then
         call logger%error("ice restart: cold-start seed refused, code "//to_string(ierr))
         return
      end if
      st%state_version = ICE_STATE_VERSION
      status = ICE_RESUME_COLD
      call logger%info("ice restart: no ice state in the file; cold start with open water")
   end subroutine ice_cice_restart_resume

end module rdb_ice_cice_restart
