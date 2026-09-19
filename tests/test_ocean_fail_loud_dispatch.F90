!! PR-6 fail-loud dispatch pack — pure-predicate tests.
!!
!! `validate_config` `error stop`s and cannot be unit-tested; the `pure`
!! predicates it consumes ARE the testable surface (that is why they are
!! pure).  Each test pins one predicate: valid tags accepted, the
!! typo/reserved/inert cases rejected, and — critically for the parser
!! flips — every VALID tag still round-trips (the regression net for the
!! flipped `case default`s).
module test_ocean_fail_loud_dispatch
   use rdb_constants, only: wp
   use rdb_coriolis_adv, only: parse_pv_variant, pv_variant_is_implemented, &
                               PV_VARIANT_SADOURNY, PV_VARIANT_SADOURNY_HK, &
                               PV_VARIANT_SADOURNY_ENERGY, PV_VARIANT_AL81, &
                               PV_VARIANT_INVALID
   use rdb_ocean_bottom_drag, only: parse_bdrag_variant, bdrag_variant_is_implemented, &
                                    BDRAG_LINEAR, BDRAG_QUADRATIC, BDRAG_INVALID
   use rdb_ocean_vmix, only: vmix_interior_closure_is_implemented, &
                             VMIX_INTERIOR_PP81, VMIX_INTERIOR_LARGE94, VMIX_INTERIOR_CVMIX
   use rdb_ocean_boundary_types, only: ocean_bc_type_from_string, &
                                       OBC_WALL, OBC_OPEN, OBC_TIDAL, OBC_NESTED, &
                                       OBC_INFLOW, OBC_DISCHARGE, OBC_CLAMPED, OBC_SPONGE, &
                                       OBC_CHAPMAN, OBC_PERIODIC, OBC_TRIPOLAR_FOLD, OBC_INVALID
   use rdb_ocean_pressure_force, only: gprime_nz_is_supported, &
                                       OPGF_VARIANT_GPRIME, OPGF_VARIANT_FV_LITE
   use rdb_ocean_lateral_mix, only: leith_biharm_is_inert, &
                                    LMIX_LEITH, LMIX_LEITH_BIHARM
   use rdb_ocean_tidal_mixing, only: tidal_mixing_is_inert
   use rdb_vcoord, only: vcoord_h_min_role, vcoord_h_min_is_coherent, &
                         VCOORD_HMIN_INERT, VCOORD_HMIN_KEEPALIVE, VCOORD_HMIN_UNUSED
   use rdb_constants, only: H_VANISHED, VCOORD_ZSTAR_FULL, VCOORD_Z_FIXED, &
                            VCOORD_RHO, VCOORD_HYCOM, VCOORD_SIGMA, VCOORD_ZSTAR
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_fail_loud_dispatch_tests

contains

   subroutine collect_ocean_fail_loud_dispatch_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("pv_variant_predicate", test_pv_variant_predicate), &
                  new_unittest("pv_parse_invalid", test_pv_parse_invalid), &
                  new_unittest("bdrag_variant_predicate", test_bdrag_variant_predicate), &
                  new_unittest("vmix_interior_predicate", test_vmix_interior_predicate), &
                  new_unittest("obc_parse_invalid", test_obc_parse_invalid), &
                  new_unittest("gprime_nz_predicate", test_gprime_nz_predicate), &
                  new_unittest("leith_biharm_inert_predicate", test_leith_biharm_inert_predicate), &
                  new_unittest("tidal_mixing_inert_predicate", test_tidal_mixing_inert_predicate), &
                  new_unittest("vcoord_h_min_role_split", test_vcoord_h_min_role_split), &
                  new_unittest("vcoord_h_min_coherence", test_vcoord_h_min_coherence) &
                  ]
   end subroutine collect_ocean_fail_loud_dispatch_tests

   subroutine test_pv_variant_predicate(error)
      !! AL81 promises simultaneous energy+enstrophy conservation; the
      !! Sadourny kernel delivers enstrophy only.  `.false.` for AL81 is
      !! the assertion that the model refuses to substitute one
      !! conservation law for another.
      type(error_type), allocatable, intent(out) :: error
      call check(error, pv_variant_is_implemented(PV_VARIANT_SADOURNY), "SADOURNY implemented")
      if (allocated(error)) return
      call check(error, pv_variant_is_implemented(PV_VARIANT_SADOURNY_HK), "SADOURNY_HK implemented")
      if (allocated(error)) return
      call check(error, pv_variant_is_implemented(PV_VARIANT_SADOURNY_ENERGY), &
                 "SADOURNY_ENERGY implemented")
      if (allocated(error)) return
      call check(error,.not. pv_variant_is_implemented(PV_VARIANT_AL81), &
                 "AL81 is reserved-but-unwired — must be rejected")
      if (allocated(error)) return
      call check(error,.not. pv_variant_is_implemented(PV_VARIANT_INVALID), &
                 "INVALID must be rejected")
   end subroutine test_pv_variant_predicate

   subroutine test_pv_parse_invalid(error)
      !! The parser must DISTINGUISH a typo (INVALID) from the
      !! reserved-but-unwired al81 (AL81): both abort, different messages,
      !! and a future AL81 PR flips one predicate line not the parser.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_pv_variant("nonsense") == PV_VARIANT_INVALID, &
                 "typo -> PV_VARIANT_INVALID")
      if (allocated(error)) return
      call check(error, parse_pv_variant("al81") == PV_VARIANT_AL81, "al81 -> PV_VARIANT_AL81")
      if (allocated(error)) return
      call check(error, parse_pv_variant("sadourny") == PV_VARIANT_SADOURNY, "sadourny round-trips")
      if (allocated(error)) return
      call check(error, parse_pv_variant("sadourny_hk") == PV_VARIANT_SADOURNY_HK, "hk round-trips")
      if (allocated(error)) return
      call check(error, parse_pv_variant("sadourny_energy") == PV_VARIANT_SADOURNY_ENERGY, &
                 "energy round-trips")
   end subroutine test_pv_parse_invalid

   subroutine test_bdrag_variant_predicate(error)
      !! Quadratic (τ=ρ·C_d·|u|·u) and linear (τ=ρ·r·u) have different
      !! coefficient dimensions and energy-decay laws; a typo must not
      !! choose between them.
      type(error_type), allocatable, intent(out) :: error
      call check(error, bdrag_variant_is_implemented(BDRAG_LINEAR), "LINEAR implemented")
      if (allocated(error)) return
      call check(error, bdrag_variant_is_implemented(BDRAG_QUADRATIC), "QUADRATIC implemented")
      if (allocated(error)) return
      call check(error,.not. bdrag_variant_is_implemented(BDRAG_INVALID), "INVALID rejected")
      if (allocated(error)) return
      call check(error, parse_bdrag_variant("qudratic") == BDRAG_INVALID, "typo -> BDRAG_INVALID")
      if (allocated(error)) return
      call check(error, parse_bdrag_variant("linear") == BDRAG_LINEAR, "linear round-trips")
      if (allocated(error)) return
      call check(error, parse_bdrag_variant("quadratic") == BDRAG_QUADRATIC, "quadratic round-trips")
   end subroutine test_bdrag_variant_predicate

   subroutine test_vmix_interior_predicate(error)
      !! Selecting a kernel-less interior closure leaves kv/kt stale — no
      !! interior mixing at all.  Only PP81 has a kernel.
      type(error_type), allocatable, intent(out) :: error
      call check(error, vmix_interior_closure_is_implemented(VMIX_INTERIOR_PP81), &
                 "PP81 implemented")
      if (allocated(error)) return
      call check(error,.not. vmix_interior_closure_is_implemented(VMIX_INTERIOR_LARGE94), &
                 "LARGE94 reserved-but-unwired — rejected")
      if (allocated(error)) return
      call check(error,.not. vmix_interior_closure_is_implemented(VMIX_INTERIOR_CVMIX), &
                 "CVMIX reserved-but-unwired — rejected")
   end subroutine test_vmix_interior_predicate

   subroutine test_obc_parse_invalid(error)
      !! A typo can no longer close a boundary and reflect every outgoing
      !! gravity wave.  The round-trip half is the regression net for the
      !! flipped `case default` that 11 valid tags fall past.
      type(error_type), allocatable, intent(out) :: error
      call check(error, ocean_bc_type_from_string("opne") == OBC_INVALID, "typo -> OBC_INVALID")
      if (allocated(error)) return
      ! every valid tag still round-trips
      call check(error, ocean_bc_type_from_string("wall") == OBC_WALL, "wall")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("open") == OBC_OPEN, "open")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("tidal") == OBC_TIDAL, "tidal")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("nested") == OBC_NESTED, "nested")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("inflow") == OBC_INFLOW, "inflow")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("discharge") == OBC_DISCHARGE, "discharge")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("clamped") == OBC_CLAMPED, "clamped")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("sponge") == OBC_SPONGE, "sponge")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("chapman") == OBC_CHAPMAN, "chapman")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("periodic") == OBC_PERIODIC, "periodic")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("tripolar_fold") == OBC_TRIPOLAR_FOLD, &
                 "tripolar_fold")
      if (allocated(error)) return
      ! case-insensitivity preserved by the flip
      call check(error, ocean_bc_type_from_string("OPEN") == OBC_OPEN, "OPEN case-insensitive")
      if (allocated(error)) return
      call check(error, ocean_bc_type_from_string("Open") == OBC_OPEN, "Open case-insensitive")
   end subroutine test_obc_parse_invalid

   subroutine test_gprime_nz_predicate(error)
      !! gprime hard-writes only top+bottom, so nz/=2 is unsupported; the
      !! FV_LITE case asserts the guard has not over-fired onto the
      !! general-nz variants.
      type(error_type), allocatable, intent(out) :: error
      call check(error, gprime_nz_is_supported(OPGF_VARIANT_GPRIME, 2), "gprime nz=2 supported")
      if (allocated(error)) return
      call check(error,.not. gprime_nz_is_supported(OPGF_VARIANT_GPRIME, 3), "gprime nz=3 rejected")
      if (allocated(error)) return
      call check(error,.not. gprime_nz_is_supported(OPGF_VARIANT_GPRIME, 50), &
                 "gprime nz=50 rejected (the live in-tree bug)")
      if (allocated(error)) return
      call check(error,.not. gprime_nz_is_supported(OPGF_VARIANT_GPRIME, 1), "gprime nz=1 rejected")
      if (allocated(error)) return
      call check(error, gprime_nz_is_supported(OPGF_VARIANT_FV_LITE, 50), &
                 "fv_lite nz=50 supported — guard must not over-fire")
   end subroutine test_gprime_nz_predicate

   subroutine test_leith_biharm_inert_predicate(error)
      !! ν₄ is linear in c_leith_bi ⇒ c=0 is provably a no-op.  The third
      !! case asserts the guard does not fire on a closure that does not
      !! read the knob.
      type(error_type), allocatable, intent(out) :: error
      call check(error, leith_biharm_is_inert(LMIX_LEITH_BIHARM, 0.0_wp), &
                 "leith_biharm + c=0 is inert")
      if (allocated(error)) return
      call check(error,.not. leith_biharm_is_inert(LMIX_LEITH_BIHARM, 0.06_wp), &
                 "leith_biharm + c>0 is not inert")
      if (allocated(error)) return
      call check(error,.not. leith_biharm_is_inert(LMIX_LEITH, 0.0_wp), &
                 "plain leith does not read c_leith_bi — guard must not fire")
   end subroutine test_leith_biharm_inert_predicate

   subroutine test_tidal_mixing_inert_predicate(error)
      !! Kd ∝ E ⇒ E=0 ⇒ Kd ≡ 0.  The e_compute=.true. case is the one an
      !! over-eager guard breaks; a disabled closure is not "inert".
      type(error_type), allocatable, intent(out) :: error
      call check(error, tidal_mixing_is_inert(.true., 0.0_wp, .false.), &
                 "enabled, E=0, no e_compute -> inert")
      if (allocated(error)) return
      call check(error,.not. tidal_mixing_is_inert(.true., 0.0_wp, .true.), &
                 "e_compute supplies E -> not inert")
      if (allocated(error)) return
      call check(error,.not. tidal_mixing_is_inert(.true., 1.0e-3_wp, .false.), &
                 "prescribed E>0 -> not inert")
      if (allocated(error)) return
      call check(error,.not. tidal_mixing_is_inert(.false., 0.0_wp, .false.), &
                 "disabled is not 'inert'")
   end subroutine test_tidal_mixing_inert_predicate

   subroutine test_vcoord_h_min_role_split(error)
      !! `zstar_h_min` spells TWO opposite contracts and the coordinate
      !! family — not the value — picks which.  Pinning the split is the
      !! point: the geometric families floor BELOW-BED filler that must stay
      !! vanished, the density families floor REAL collapsed layers that must
      !! stay alive (and therefore get `max(zstar_h_min, 2*H_VANISHED)`
      !! instead).  A future family added to the wrong bucket silently
      !! reclassifies "thin".
      type(error_type), allocatable, intent(out) :: error
      call check(error, vcoord_h_min_role(VCOORD_ZSTAR_FULL) == VCOORD_HMIN_INERT, &
                 "ZSTAR_FULL floors below-bed filler -> INERT role")
      if (allocated(error)) return
      call check(error, vcoord_h_min_role(VCOORD_Z_FIXED) == VCOORD_HMIN_INERT, &
                 "Z_FIXED floors above-column filler -> INERT role")
      if (allocated(error)) return
      call check(error, vcoord_h_min_role(VCOORD_RHO) == VCOORD_HMIN_KEEPALIVE, &
                 "RHO inflates tracer-carrying layers -> KEEPALIVE role")
      if (allocated(error)) return
      call check(error, vcoord_h_min_role(VCOORD_HYCOM) == VCOORD_HMIN_KEEPALIVE, &
                 "HYCOM runs the RHO inversion -> KEEPALIVE role")
      if (allocated(error)) return
      call check(error, vcoord_h_min_role(VCOORD_SIGMA) == VCOORD_HMIN_UNUSED, &
                 "sigma never reads zstar_h_min")
      if (allocated(error)) return
      call check(error, vcoord_h_min_role(VCOORD_ZSTAR) == VCOORD_HMIN_UNUSED, &
                 "zstar-lite never reads zstar_h_min")
   end subroutine test_vcoord_h_min_role_split

   subroutine test_vcoord_h_min_coherence(error)
      !! The guard `validate_config` consumes.  The predicate states the RULE;
      !! the call site picks the severity (non-positive aborts, above-marker
      !! warns today because the Python worked example sits in that band), so
      !! promoting that to fail-loud never re-derives the rule — it edits one
      !! line.  Two rejections here, and — the part an over-eager guard
      !! breaks — every shipped value still accepted:
      !! the 1.0e-4 type default AND the 1.5e-4 five shipped namelists set,
      !! which sits exactly ON `H_VANISHED` and is legal because every
      !! downstream vanish gate is a strict `>`.
      type(error_type), allocatable, intent(out) :: error
      call check(error, vcoord_h_min_is_coherent(VCOORD_ZSTAR_FULL, 1.0e-4_wp), &
                 "type default 1.0e-4 accepted under ZSTAR_FULL")
      if (allocated(error)) return
      call check(error, vcoord_h_min_is_coherent(VCOORD_ZSTAR_FULL, H_VANISHED), &
                 "shipped 1.5e-4 == H_VANISHED accepted (gates are strict '>')")
      if (allocated(error)) return
      call check(error, vcoord_h_min_is_coherent(VCOORD_Z_FIXED, H_VANISHED), &
                 "shipped 1.5e-4 accepted under Z_FIXED too")
      if (allocated(error)) return
      call check(error,.not. vcoord_h_min_is_coherent(VCOORD_ZSTAR_FULL, 1.0e-3_wp), &
                 "above H_VANISHED under ZSTAR_FULL -> filler goes live, rejected")
      if (allocated(error)) return
      call check(error,.not. vcoord_h_min_is_coherent(VCOORD_Z_FIXED, 2.0_wp*H_VANISHED), &
                 "2*H_VANISHED is the RHO keep-alive floor, NOT legal for Z_FIXED")
      if (allocated(error)) return
      ! The density families read the same knob under the opposite contract:
      ! a large value there is a meaningful pre-compaction strip threshold,
      ! and the regrid lifts its own floor to max(h_min, 2*H_VANISHED).
      call check(error, vcoord_h_min_is_coherent(VCOORD_RHO, 1.0e-3_wp), &
                 "RHO keep-alive role permits a floor above H_VANISHED")
      if (allocated(error)) return
      call check(error, vcoord_h_min_is_coherent(VCOORD_HYCOM, 1.0e-2_wp), &
                 "HYCOM keep-alive role permits a floor above H_VANISHED")
      if (allocated(error)) return
      ! Non-positive defeats the knob's one purpose (never an exactly-zero
      ! target_h) on EVERY family, including the ones that ignore it.
      call check(error,.not. vcoord_h_min_is_coherent(VCOORD_ZSTAR_FULL, 0.0_wp), &
                 "zero floor rejected — target_h would be exactly 0")
      if (allocated(error)) return
      call check(error,.not. vcoord_h_min_is_coherent(VCOORD_RHO, -1.0e-4_wp), &
                 "negative floor rejected on the density families too")
      if (allocated(error)) return
      call check(error,.not. vcoord_h_min_is_coherent(VCOORD_SIGMA, 0.0_wp), &
                 "non-positive rejected even where the knob is unused")
   end subroutine test_vcoord_h_min_coherence

end module test_ocean_fail_loud_dispatch
