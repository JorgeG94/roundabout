!! Sea-ice incremental-remapping reconstruction (D10.4): masks, limited
!! gradient and `construct_fields` of `rdb_ice_remap_recon`.
module test_cice_remap_recon
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Tests written for the roundabout port of the CICE incremental-remapping reconstruction
   ! (CICE Consortium, https://github.com/CICE-Consortium); the independent host oracle
   ! follows the published scheme (Dukowicz and Baumgardner 2000; Lipscomb and Hunke 2004),
   ! not the port.  Modified by the roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   !! **What is pinned.**  (1) limiter bounds: for a random positive state
   !! the reconstructed value at the four corners of every cell (mass and every
   !! tracer, type 1 and type 2) lies within the min/max of the 3x3
   !! neighbourhood, the worst violation measured; the limiter is active in a
   !! large fraction of the cells; (2) mean preservation, by a 2x2 Gauss
   !! quadrature that shares nothing with the module: `<m> = aim`,
   !! `<m T> = aim tm` (type 1), `<m T_parent T> = aim tm_parent tm` (type 2),
   !! on a random masked state; (3) a field linear in the index coordinates is
   !! reconstructed exactly (limiter inactive); (4) a constant tracer has a
   !! zero gradient whatever the mass; (5) mask substitution against two hand
   !! computations (a land neighbour for the mass, an ice-free neighbour for a
   !! tracer, each compared with what NOT substituting would give); (6)
   !! stretched metrics: the reconstruction is metric-free, so a field linear
   !! in physical distance on a geometrically stretched grid carries the
   !! closed-form error `(r-1)^2/(4r)`, and a linearly varying width is
   !! exact; (7) non-finite input raises the counter; (8) the library against
   !! an independent host implementation (per-neighbour loops, stored
   !! centroids) on EVERY output of EVERY cell and category of random masked
   !! states, the ring outside the range untouched; (9) the numbers of a
   !! stdlib-only Python transcript; (10) configuration refusals; (11) an empty
   !! domain and an all-land mask; (12) throughput.
   use, intrinsic :: iso_fortran_env, only: int8, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_transport_pack, only: transport_workspace_t, transport_tracer_table_t, transport_tracer_table_build, &
                                     transport_tracer_table_from_depend, TT_OK
   use rdb_ice_remap_recon, only: remap_recon_t, remap_reconstruct, RECON_CFG_OK, RECON_CFG_TABLE, RECON_CFG_WORKSPACE, &
                                  RECON_CFG_RANGE, RECON_CFG_MASK, RECON_XXAV
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_remap_recon_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: TOL = 1.0e-12_wp
      !! Device vs host / identity tolerance relative to the field scale: FMA contraction is on.
#else
   real(wp), parameter :: TOL = 1.0e-14_wp
      !! CPU builds contract nothing: the expressions agree to a few ulp.
#endif
   real(wp), parameter :: GAUSS = 0.5_wp/1.7320508075688772935_wp
      !! Abscissa of the 2-point Gauss rule on the unit cell, `+-1/(2 sqrt 3)`.

   real(wp), parameter :: PY_AIM(9, 4) = reshape([ &
                                               2.99999999999999989e-01_wp, 3.49999999999999978e-01_wp, 4.19999999999999984e-01_wp, &
                                               2.80000000000000027e-01_wp, 4.00000000000000022e-01_wp, 5.50000000000000044e-01_wp, &
                                               3.30000000000000016e-01_wp, 4.69999999999999973e-01_wp, 5.00000000000000000e-01_wp, &
                                               2.00000000000000011e-01_wp, 2.99999999999999989e-01_wp, 0.00000000000000000e+00_wp, &
                                               7.70000000000000018e-01_wp, 4.50000000000000011e-01_wp, 5.00000000000000000e-01_wp, &
                                               1.00000000000000006e-01_wp, 5.99999999999999978e-01_wp, 5.50000000000000044e-01_wp, &
                                               5.00000000000000000e-01_wp, 5.00000000000000000e-01_wp, 5.00000000000000000e-01_wp, &
                                               4.50000000000000011e-01_wp, 5.00000000000000000e-01_wp, 9.00000000000000022e-01_wp, &
                                               5.00000000000000000e-01_wp, 5.00000000000000000e-01_wp, 5.00000000000000000e-01_wp, &
                                               5.00000000000000028e-02_wp, 2.50000000000000000e-01_wp, 5.99999999999999978e-01_wp, &
                                               5.99999999999999978e-02_wp, 2.99999999999999989e-01_wp, 8.00000000000000044e-01_wp, &
                                                7.00000000000000067e-02_wp, 3.49999999999999978e-01_wp, 9.49999999999999956e-01_wp &
                                                 ], [9, 4])
   real(wp), parameter :: PY_HM(9, 4) = reshape([ &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               0.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               0.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                                1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp &
                                                ], [9, 4])
   real(wp), parameter :: PY_HI(9, 4) = reshape([ &
                                               1.10000000000000009e+00_wp, 1.30000000000000004e+00_wp, 1.19999999999999996e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.50000000000000000e+00_wp, 1.89999999999999991e+00_wp, &
                                               1.39999999999999991e+00_wp, 1.60000000000000009e+00_wp, 1.69999999999999996e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.19999999999999996e+00_wp, 0.00000000000000000e+00_wp, &
                                               0.00000000000000000e+00_wp, 1.39999999999999991e+00_wp, 1.60000000000000009e+00_wp, &
                                               6.99999999999999956e-01_wp, 1.80000000000000004e+00_wp, 1.10000000000000009e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               9.49999999999999956e-01_wp, 1.00000000000000000e+00_wp, 1.39999999999999991e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, 1.00000000000000000e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.19999999999999996e+00_wp, 1.39999999999999991e+00_wp, &
                                               1.00000000000000000e+00_wp, 1.30000000000000004e+00_wp, 1.60000000000000009e+00_wp, &
                                                0.00000000000000000e+00_wp, 1.39999999999999991e+00_wp, 1.80000000000000004e+00_wp &
                                                ], [9, 4])
   real(wp), parameter :: PY_HS(9, 4) = reshape([ &
                                               1.00000000000000006e-01_wp, 1.19999999999999996e-01_wp, 2.00000000000000011e-01_wp, &
                                               5.00000000000000028e-02_wp, 1.49999999999999994e-01_wp, 2.20000000000000001e-01_wp, &
                                               8.00000000000000017e-02_wp, 8.99999999999999967e-02_wp, 2.99999999999999989e-01_wp, &
                                               1.00000000000000006e-01_wp, 2.00000000000000011e-01_wp, 0.00000000000000000e+00_wp, &
                                               0.00000000000000000e+00_wp, 1.00000000000000006e-01_wp, 2.99999999999999989e-01_wp, &
                                               5.00000000000000028e-02_wp, 1.49999999999999994e-01_wp, 1.00000000000000006e-01_wp, &
                                               1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, &
                                               1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, &
                                               1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, 1.00000000000000006e-01_wp, &
                                               0.00000000000000000e+00_wp, 1.00000000000000006e-01_wp, 2.00000000000000011e-01_wp, &
                                               0.00000000000000000e+00_wp, 2.00000000000000011e-01_wp, 2.99999999999999989e-01_wp, &
                                                1.00000000000000006e-01_wp, 2.99999999999999989e-01_wp, 4.00000000000000022e-01_wp &
                                                ], [9, 4])
   real(wp), parameter :: PY_T3(9, 4) = reshape([ &
                                            -5.00000000000000000e+00_wp, -6.00000000000000000e+00_wp, -7.50000000000000000e+00_wp, &
                                            -4.00000000000000000e+00_wp, -8.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, &
                                            -3.00000000000000000e+00_wp, -7.00000000000000000e+00_wp, -6.50000000000000000e+00_wp, &
                                             -3.00000000000000000e+00_wp, -4.00000000000000000e+00_wp, 0.00000000000000000e+00_wp, &
                                             0.00000000000000000e+00_wp, -5.00000000000000000e+00_wp, -6.00000000000000000e+00_wp, &
                                            -2.00000000000000000e+00_wp, -7.00000000000000000e+00_wp, -8.00000000000000000e+00_wp, &
                                            -9.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, &
                                            -8.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, -3.00000000000000000e+00_wp, &
                                            -9.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, -9.00000000000000000e+00_wp, &
                                            -1.00000000000000000e+00_wp, -2.00000000000000000e+00_wp, -3.00000000000000000e+00_wp, &
                                            -1.00000000000000000e+00_wp, -2.50000000000000000e+00_wp, -4.00000000000000000e+00_wp, &
                                             -1.00000000000000000e+00_wp, -3.00000000000000000e+00_wp, -5.00000000000000000e+00_wp &
                                                ], [9, 4])
   real(wp), parameter :: PY_T4(9, 4) = reshape([ &
                                            -2.10000000000000000e+08_wp, -2.30000000000000000e+08_wp, -2.00000000000000000e+08_wp, &
                                            -2.50000000000000000e+08_wp, -2.20000000000000000e+08_wp, -2.60000000000000000e+08_wp, &
                                            -2.40000000000000000e+08_wp, -2.00000000000000000e+08_wp, -2.20000000000000000e+08_wp, &
                                             -2.00000000000000000e+08_wp, -2.20000000000000000e+08_wp, 0.00000000000000000e+00_wp, &
                                             0.00000000000000000e+00_wp, -2.40000000000000000e+08_wp, -2.10000000000000000e+08_wp, &
                                            -2.30000000000000000e+08_wp, -2.50000000000000000e+08_wp, -2.20000000000000000e+08_wp, &
                                            -2.00000000000000000e+08_wp, -2.00000000000000000e+08_wp, -2.00000000000000000e+08_wp, &
                                            -2.00000000000000000e+08_wp, -2.00000000000000000e+08_wp, -1.00000000000000000e+08_wp, &
                                            -2.00000000000000000e+08_wp, -2.00000000000000000e+08_wp, -2.00000000000000000e+08_wp, &
                                            -2.00000000000000000e+08_wp, -2.10000000000000000e+08_wp, -2.20000000000000000e+08_wp, &
                                            -2.00000000000000000e+08_wp, -2.20000000000000000e+08_wp, -2.40000000000000000e+08_wp, &
                                              0.00000000000000000e+00_wp, -2.30000000000000000e+08_wp, -2.60000000000000000e+08_wp &
                                                ], [9, 4])
   real(wp), parameter :: PY_OUT(17, 4) = reshape([ &
                                               4.00000000000000022e-01_wp, 1.35000000000000009e-01_wp, 5.99999999999999978e-02_wp, &
                                               1.48546874999999989e+00_wp, 4.49999999999999956e-01_wp, 1.50000000000000022e-01_wp, &
                                               5.28525390624999952e-02_wp, 2.07122395833333299e-02_wp, 1.47796874999999994e-01_wp, &
                                             8.49999999999999922e-02_wp, -1.49999999999999994e-02_wp, -7.94621295279912143e+00_wp, &
                                            -1.75631174533479717e+00_wp, -3.51262349066959412e-01_wp, -2.20046420898437500e+08_wp, &
                                              -5.00000000000000000e+06_wp, 1.50000000000000000e+07_wp, 4.50000000000000011e-01_wp, &
                                               2.49999999999999944e-02_wp, 1.49999999999999994e-01_wp, 1.39104991394148025e+00_wp, &
                                               5.94836488812392594e-01_wp, 2.23063683304647264e-01_wp, 4.00069666420785272e-02_wp, &
                                               4.08777969018932902e-02_wp, 1.00000000000000006e-01_wp, 1.49999999999999994e-01_wp, &
                                            -2.50000000000000083e-02_wp, -4.94444444444444375e+00_wp, -3.00000000000000000e+00_wp, &
                                             -1.50000000000000000e+00_wp, -2.39991284107338786e+08_wp, 1.00087158926612176e+07_wp, &
                                              -1.00087158926612176e+07_wp, 5.00000000000000000e-01_wp, 9.99999999999999778e-02_wp, &
                                               0.00000000000000000e+00_wp, 9.98387096774193572e-01_wp, 9.67741935483871912e-02_wp, &
                                               0.00000000000000000e+00_wp, 2.47043010752688204e-02_wp, 0.00000000000000000e+00_wp, &
                                               1.00000000000000006e-01_wp, 0.00000000000000000e+00_wp, 0.00000000000000000e+00_wp, &
                                              -9.00000000000000000e+00_wp, 0.00000000000000000e+00_wp, 0.00000000000000000e+00_wp, &
                                              -2.00000000000000000e+08_wp, 0.00000000000000000e+00_wp, 0.00000000000000000e+00_wp, &
                                               2.99999999999999989e-01_wp, 3.69999999999999996e-01_wp, 4.99999999999999889e-02_wp, &
                                               1.26777777777777789e+00_wp, 3.00000000000000044e-01_wp, 9.99999999999999778e-02_wp, &
                                               1.19461063627730277e-01_wp, 1.99548907882241174e-02_wp, 1.83194444444444465e-01_wp, &
                                              1.49999999999999994e-01_wp, 9.99999999999999917e-02_wp, -2.33888888888888857e+00_wp, &
                                            -1.50000000000000000e+00_wp, -5.00000000000000000e-01_wp, -2.17411229819563150e+08_wp, &
                                                  -2.00000000000000000e+07_wp, -1.00000000000000000e+07_wp &
                                                  ], [17, 4])

   type :: rig_t
      !! A packed workspace, the reconstruction workspace and the land mask, all device-mapped.
      integer :: nx = 0, ny = 0, ncat = 0, ntrace = 0
      type(transport_workspace_t) :: ws
      type(remap_recon_t) :: rw
      real(wp), allocatable :: hm(:, :)
   end type rig_t

   type :: orc_t
      !! Output of the independent host oracle.
      integer(int8), allocatable :: mmask(:, :, :), tmask(:, :, :, :)
      real(wp), allocatable :: mc(:, :, :), mx(:, :, :), my(:, :, :)
      real(wp), allocatable :: tc(:, :, :, :), tx(:, :, :, :), ty(:, :, :, :)
   end type orc_t

contains

   subroutine collect_cice_remap_recon_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("limiter_bounds_corners", test_bounds), &
                  new_unittest("mean_preservation_quadrature", test_mean), &
                  new_unittest("linear_field_exact", test_linear), &
                  new_unittest("constant_tracer_zero_gradient", test_constant), &
                  new_unittest("mask_substitution_hand_computed", test_masks_hand), &
                  new_unittest("stretched_metrics", test_stretched), &
                  new_unittest("type2_mask_substitution_hand_computed", test_type2_mask_hand), &
                  new_unittest("land_mask_must_be_0_or_1", test_hm_contract), &
                  new_unittest("nonfinite_input_counted", test_nan), &
                  new_unittest("oracle_every_output_every_cell", test_oracle), &
                  new_unittest("python_transcript", test_python), &
                  new_unittest("configuration_refusals", test_cfg), &
                  new_unittest("empty_and_all_land", test_empty), &
                  new_unittest("throughput", test_throughput) &
                  ]
   end subroutine collect_cice_remap_recon_tests

   ! ------------------------------------------------------------------
   ! rig
   ! ------------------------------------------------------------------

   subroutine rig_init(r, nx, ny, ncat, tab)
      type(rig_t), intent(inout) :: r
      integer, intent(in) :: nx, ny, ncat
      type(transport_tracer_table_t), intent(in) :: tab
      integer :: ierr

      r%nx = nx; r%ny = ny; r%ncat = ncat; r%ntrace = tab%ntrace
      call r%ws%init(nx, ny, ncat, tab, ierr)
      call r%rw%init(nx, ny, ncat, tab%ntrace, ierr)
      allocate (r%hm(nx, ny))
      r%hm = 1.0_wp
      call r%ws%enter_data()
      call r%rw%enter_data()
      !$acc enter data copyin(r%hm)
   end subroutine rig_init

   subroutine rig_destroy(r)
      type(rig_t), intent(inout) :: r
      !$acc exit data delete(r%hm)
      call r%rw%destroy()
      call r%ws%destroy()
      if (allocated(r%hm)) deallocate (r%hm)
   end subroutine rig_destroy

   subroutine rig_push(r)
      !! Host to device: the packed inputs, the mask and (sentinel-fillable) outputs.
      type(rig_t), intent(inout) :: r
      !$acc update device(r%ws%aim, r%ws%trm, r%hm)
      !$acc update device(r%rw%mmask, r%rw%tmask, r%rw%mc, r%rw%mx, r%rw%my, r%rw%tc, r%rw%tx, r%rw%ty)
   end subroutine rig_push

   subroutine rig_pull(r)
      type(rig_t), intent(inout) :: r
      !$acc update self(r%rw%mmask, r%rw%tmask, r%rw%mc, r%rw%mx, r%rw%my, r%rw%tc, r%rw%tx, r%rw%ty)
   end subroutine rig_pull

   subroutine rig_run(r, is, ie, js, je, nbad, cfg)
      type(rig_t), intent(inout) :: r
      integer, intent(in) :: is, ie, js, je
      integer, intent(out) :: nbad, cfg
      call remap_reconstruct(r%ws, r%rw, is, ie, js, je, PUNY, r%hm, nbad, cfg)
   end subroutine rig_run

   subroutine table_nil(nil, tab)
      integer, intent(in) :: nil
      type(transport_tracer_table_t), intent(out) :: tab
      call transport_tracer_table_build(nil, 1, tab)
   end subroutine table_nil

   subroutine sentinel_outputs(r)
      type(rig_t), intent(inout) :: r
      r%rw%mmask = 5_int8
      r%rw%tmask = 5_int8
      r%rw%mc = 7.0_wp; r%rw%mx = 7.0_wp; r%rw%my = 7.0_wp
      r%rw%tc = 7.0_wp; r%rw%tx = 7.0_wp; r%rw%ty = 7.0_wp
   end subroutine sentinel_outputs

   subroutine fill_random(r, seed, nil)
      !! Random masked packed state: ~25 % ice-free cells, ~5 % sub-puny areas, ~8 % land, parent
      !! tracers (`hi`, `hs`) zero in ~10 % of the ice cells; packed like the pack (tracers zero
      !! where `aim <= puny`).  The land mask is independent of the ice (a stronger test).
      type(rig_t), intent(inout) :: r
      integer(int64), intent(in) :: seed
      integer, intent(in) :: nil
      integer(int64) :: st
      integer :: i, j, n, nt
      real(wp) :: u

      st = lcg_next(seed)
      do n = 0, r%ncat
         do j = 1, r%ny
            do i = 1, r%nx
               st = lcg_next(st); u = lcg_unit(st)
               if (u < 0.25_wp) then
                  r%ws%aim(i, j, n) = 0.0_wp
               else if (u < 0.30_wp) then
                  st = lcg_next(st)
                  r%ws%aim(i, j, n) = 0.5_wp*PUNY*lcg_unit(st)
               else
                  st = lcg_next(st)
                  r%ws%aim(i, j, n) = 0.02_wp + 0.98_wp*lcg_unit(st)
               end if
            end do
         end do
      end do
      do n = 1, r%ncat
         do nt = 1, r%ntrace
            do j = 1, r%ny
               do i = 1, r%nx
                  st = lcg_next(st); u = lcg_unit(st)
                  if (r%ws%aim(i, j, n) <= PUNY) then
                     r%ws%trm(i, j, nt, n) = 0.0_wp
                  else if (nt <= 2 .and. u < 0.10_wp) then
                     r%ws%trm(i, j, nt, n) = 0.0_wp
                  else if (nt <= 2) then
                     r%ws%trm(i, j, nt, n) = 0.05_wp + 3.0_wp*u
                  else
                     r%ws%trm(i, j, nt, n) = -(1.0_wp + 9.0_wp*u)*real(nt, wp)
                  end if
               end do
            end do
         end do
      end do
      do j = 1, r%ny
         do i = 1, r%nx
            st = lcg_next(st)
            r%hm(i, j) = 1.0_wp
            if (lcg_unit(st) < 0.08_wp) r%hm(i, j) = 0.0_wp
         end do
      end do
      if (nil < 0) r%hm = 1.0_wp
   end subroutine fill_random

   subroutine fill_smooth_positive(r, seed)
      !! Fully iced random positive state (every `tmask = 1`): rough enough that the limiter is active.
      type(rig_t), intent(inout) :: r
      integer(int64), intent(in) :: seed
      integer(int64) :: st
      integer :: i, j, n, nt

      st = lcg_next(seed)
      do n = 0, r%ncat
         do j = 1, r%ny
            do i = 1, r%nx
               st = lcg_next(st)
               r%ws%aim(i, j, n) = 0.1_wp + 0.9_wp*lcg_unit(st)
            end do
         end do
      end do
      do n = 1, r%ncat
         do nt = 1, r%ntrace
            do j = 1, r%ny
               do i = 1, r%nx
                  st = lcg_next(st)
                  r%ws%trm(i, j, nt, n) = real(nt, wp)*(0.5_wp + lcg_unit(st))
               end do
            end do
         end do
      end do
      r%hm = 1.0_wp
   end subroutine fill_smooth_positive

   ! ------------------------------------------------------------------
   ! independent oracle
   ! ------------------------------------------------------------------

   subroutine orc_grad(win, ok, cnx, cny, gx, gy)
      !! Limited gradient written per neighbour: substitute by loops, extremes by loops, the four
      !! corners by a loop over the sign pairs.
      real(wp), intent(in) :: win(-1:1, -1:1), cnx, cny
      logical, intent(in) :: ok(-1:1, -1:1)
      real(wp), intent(out) :: gx, gy
      real(wp) :: v(-1:1, -1:1), lo, hi, gxt, gyt, qmn, qmx, dev, a1, a2, alpha
      integer :: di, dj, sx, sy

      gx = 0.0_wp; gy = 0.0_wp
      if (.not. ok(0, 0)) return
      do dj = -1, 1
         do di = -1, 1
            if (ok(di, dj)) then
               v(di, dj) = win(di, dj)
            else
               v(di, dj) = win(0, 0)
            end if
         end do
      end do
      gxt = (v(1, 0) - v(-1, 0))*0.5_wp
      gyt = (v(0, 1) - v(0, -1))*0.5_wp
      lo = v(0, 0); hi = v(0, 0)
      do dj = -1, 1
         do di = -1, 1
            lo = min(lo, v(di, dj)); hi = max(hi, v(di, dj))
         end do
      end do
      lo = lo - win(0, 0); hi = hi - win(0, 0)
      qmn = huge(1.0_wp); qmx = -huge(1.0_wp)
      do sy = -1, 1, 2
         do sx = -1, 1, 2
            dev = (0.5_wp*real(sx, wp) - cnx)*gxt + (0.5_wp*real(sy, wp) - cny)*gyt
            qmn = min(qmn, dev); qmx = max(qmx, dev)
         end do
      end do
      a1 = 1.0_wp; a2 = 1.0_wp
      if (abs(qmn) > abs(lo)) a1 = max(0.0_wp, lo/qmn)
      if (abs(qmx) > abs(hi)) a2 = max(0.0_wp, hi/qmx)
      alpha = min(a1, a2)
      gx = alpha*gxt
      gy = alpha*gyt
   end subroutine orc_grad

   subroutine oracle(r, is, ie, js, je, o)
      type(rig_t), intent(in) :: r
      integer, intent(in) :: is, ie, js, je
      type(orc_t), intent(out) :: o
      integer :: i, j, n, nt, di, dj, p
      real(wp) :: win(-1:1, -1:1), gx, gy, mxav, myav, mtx(r%ntrace), mty(r%ntrace), mmv
      logical :: ok(-1:1, -1:1)
      logical :: hasdep(r%ntrace)

      allocate (o%mmask(r%nx, r%ny, 0:r%ncat), o%tmask(r%nx, r%ny, r%ntrace, r%ncat))
      allocate (o%mc(r%nx, r%ny, 0:r%ncat), o%mx(r%nx, r%ny, 0:r%ncat), o%my(r%nx, r%ny, 0:r%ncat))
      allocate (o%tc(r%nx, r%ny, r%ntrace, r%ncat), o%tx(r%nx, r%ny, r%ntrace, r%ncat), o%ty(r%nx, r%ny, r%ntrace, r%ncat))
      o%mmask = 0_int8; o%tmask = 0_int8
      o%mc = 0.0_wp; o%mx = 0.0_wp; o%my = 0.0_wp; o%tc = 0.0_wp; o%tx = 0.0_wp; o%ty = 0.0_wp
      hasdep = r%ws%tab%has_dependents(1:r%ntrace)
      do n = 0, r%ncat
         do j = js - 1, je + 1
            do i = is - 1, ie + 1
               if (r%ws%aim(i, j, n) > PUNY) o%mmask(i, j, n) = 1_int8
               if (n >= 1) then
                  do nt = 1, r%ntrace
                     if (hasdep(nt) .and. o%mmask(i, j, n) == 1_int8 .and. abs(r%ws%trm(i, j, nt, n)) > PUNY) &
                        o%tmask(i, j, nt, n) = 1_int8
                  end do
               end if
            end do
         end do
      end do
      do n = 0, r%ncat
         do j = js, je
            do i = is, ie
               do dj = -1, 1
                  do di = -1, 1
                     win(di, dj) = r%ws%aim(i + di, j + dj, n)
                     ok(di, dj) = r%hm(i + di, j + dj) > PUNY
                  end do
               end do
               call orc_grad(win, ok, 0.0_wp, 0.0_wp, gx, gy)
               o%mx(i, j, n) = gx; o%my(i, j, n) = gy
               mmv = r%ws%aim(i, j, n)
               if (mmv > PUNY) o%mc(i, j, n) = mmv
               if (n == 0 .or. .not. (mmv > PUNY)) cycle
               mxav = gx*RECON_XXAV/mmv
               myav = gy*RECON_XXAV/mmv
               mtx = 0.0_wp; mty = 0.0_wp
               do nt = 1, r%ntrace
                  do dj = -1, 1
                     do di = -1, 1
                        win(di, dj) = r%ws%trm(i + di, j + dj, nt, n)
                     end do
                  end do
                  if (r%ws%tab%tracer_type(nt) == 1) then
                     do dj = -1, 1
                        do di = -1, 1
                           ok(di, dj) = o%mmask(i + di, j + dj, n) == 1_int8
                        end do
                     end do
                     call orc_grad(win, ok, mxav, myav, gx, gy)
                     o%tx(i, j, nt, n) = gx; o%ty(i, j, nt, n) = gy
                     o%tc(i, j, nt, n) = win(0, 0) - gx*mxav - gy*myav
                     if (hasdep(nt) .and. o%tmask(i, j, nt, n) == 1_int8) then
                        mtx(nt) = (o%mc(i, j, n)*gx + o%mx(i, j, n)*o%tc(i, j, nt, n))*RECON_XXAV &
                                  *(1.0_wp/(mmv*win(0, 0)))
                        mty(nt) = (o%mc(i, j, n)*gy + o%my(i, j, n)*o%tc(i, j, nt, n))*RECON_XXAV &
                                  *(1.0_wp/(mmv*win(0, 0)))
                     end if
                  else
                     p = r%ws%tab%depend(nt)
                     do dj = -1, 1
                        do di = -1, 1
                           ok(di, dj) = o%tmask(i + di, j + dj, p, n) == 1_int8
                        end do
                     end do
                     call orc_grad(win, ok, mtx(p), mty(p), gx, gy)
                     o%tx(i, j, nt, n) = gx; o%ty(i, j, nt, n) = gy
                     o%tc(i, j, nt, n) = win(0, 0) - gx*mtx(p) - gy*mty(p)
                  end if
               end do
            end do
         end do
      end do
   end subroutine oracle

   function maxdiff3(got, want, scale) result(d)
      !! Largest `|got - want|/scale`; `huge` for a non-finite entry.
      real(wp), intent(in) :: got(:, :, :), want(:, :, :), scale
      real(wp) :: d
      integer :: i, j, k
      d = 0.0_wp
      do k = 1, size(got, 3)
         do j = 1, size(got, 2)
            do i = 1, size(got, 1)
               if (.not. (ieee_is_finite(got(i, j, k)) .and. ieee_is_finite(want(i, j, k)))) then
                  d = huge(1.0_wp)
               else if (abs(got(i, j, k) - want(i, j, k)) > d*scale) then
                  d = abs(got(i, j, k) - want(i, j, k))/scale
               end if
            end do
         end do
      end do
   end function maxdiff3

   function maxdiff4(got, want, scale) result(d)
      real(wp), intent(in) :: got(:, :, :, :), want(:, :, :, :), scale
      real(wp) :: d
      integer :: l
      d = 0.0_wp
      do l = 1, size(got, 4)
         call acc_max(d, maxdiff3(got(:, :, :, l), want(:, :, :, l), scale))
      end do
   end function maxdiff4

   ! ------------------------------------------------------------------
   ! tests
   ! ------------------------------------------------------------------

   subroutine acc_max(w, x)
      !! `w = max(w, x)`, but a non-finite `x` poisons `w` to `huge` (a NaN must never hide in a `max`).
      real(wp), intent(inout) :: w
      real(wp), intent(in) :: x
      if (.not. ieee_is_finite(x)) then
         w = huge(1.0_wp)
      else
         w = max(w, x)
      end if
   end subroutine acc_max

   subroutine corner_bounds(r, worst_m, worst_t, ncell, nlim)
      !! Worst violation of the 3x3 bounds AFTER substitution, host side: the mass over cells with
      !! `hm > puny` and area (neighbours with `hm <= puny` replaced by the cell's value); a type-1
      !! tracer over cells with area (neighbours without area replaced); a type-2 tracer over cells
      !! whose PARENT mask holds (neighbours whose parent mask is zero replaced).
      type(rig_t), intent(in) :: r
      real(wp), intent(out) :: worst_m, worst_t
      integer, intent(out) :: ncell, nlim
      integer :: i, j, n, nt, p, di, dj, sx, sy
      real(wp) :: win(-1:1, -1:1), lo, hi, sc, v, gxt
      logical :: ok(-1:1, -1:1)

      worst_m = 0.0_wp; worst_t = 0.0_wp; ncell = 0; nlim = 0
      do n = 0, r%ncat
         sc = max(maxval(abs(r%ws%aim(:, :, n))), 1.0e-300_wp)
         do j = 2, r%ny - 1
            do i = 2, r%nx - 1
               if (.not. (r%hm(i, j) > PUNY .and. r%ws%aim(i, j, n) > PUNY)) cycle
               do dj = -1, 1
                  do di = -1, 1
                     ok(di, dj) = r%hm(i + di, j + dj) > PUNY
                     win(di, dj) = merge(r%ws%aim(i + di, j + dj, n), r%ws%aim(i, j, n), ok(di, dj))
                  end do
               end do
               lo = minval(win); hi = maxval(win)
               do sy = -1, 1, 2
                  do sx = -1, 1, 2
                     v = r%rw%mc(i, j, n) + 0.5_wp*real(sx, wp)*r%rw%mx(i, j, n) + 0.5_wp*real(sy, wp)*r%rw%my(i, j, n)
                     call acc_max(worst_m, max(v - hi, lo - v)/sc)
                  end do
               end do
               gxt = 0.5_wp*(win(1, 0) - win(-1, 0))
               ncell = ncell + 1
               if (abs(r%rw%mx(i, j, n) - gxt) > 1.0e-12_wp*sc) nlim = nlim + 1
            end do
         end do
      end do
      do n = 1, r%ncat
         do nt = 1, r%ntrace
            sc = max(maxval(abs(r%ws%trm(:, :, nt, n))), 1.0e-300_wp)
            p = r%ws%tab%depend(nt)
            do j = 2, r%ny - 1
               do i = 2, r%nx - 1
                  if (.not. r%ws%aim(i, j, n) > PUNY) cycle
                  if (r%ws%tab%tracer_type(nt) == 2) then
                     if (.not. abs(r%ws%trm(i, j, p, n)) > PUNY) cycle
                  end if
                  do dj = -1, 1
                     do di = -1, 1
                        ok(di, dj) = r%ws%aim(i + di, j + dj, n) > PUNY
                        if (r%ws%tab%tracer_type(nt) == 2) ok(di, dj) = ok(di, dj) .and. abs(r%ws%trm(i + di, j + dj, p, n)) > PUNY
                        win(di, dj) = merge(r%ws%trm(i + di, j + dj, nt, n), r%ws%trm(i, j, nt, n), ok(di, dj))
                     end do
                  end do
                  lo = minval(win); hi = maxval(win)
                  do sy = -1, 1, 2
                     do sx = -1, 1, 2
                        v = r%rw%tc(i, j, nt, n) + 0.5_wp*real(sx, wp)*r%rw%tx(i, j, nt, n) &
                            + 0.5_wp*real(sy, wp)*r%rw%ty(i, j, nt, n)
                        call acc_max(worst_t, max(v - hi, lo - v)/sc)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine corner_bounds

   subroutine test_bounds(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, ncell, nlim, k
      real(wp) :: worst_t, worst_m

      call table_nil(4, tab)
      call rig_init(r, 40, 36, 2, tab)
      do k = 1, 2
         if (k == 1) then
            call fill_smooth_positive(r, 4242_int64)
         else
            call fill_random(r, 515_int64, 4)
         end if
         call rig_push(r)
         call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
         call rig_pull(r)
         call check(error, cfg, RECON_CFG_OK)
         if (allocated(error)) return
         call check(error, nbad, 0)
         if (allocated(error)) return
         call corner_bounds(r, worst_m, worst_t, ncell, nlim)
         call check(error, ncell > 100 .and. real(nlim, wp) > 0.1_wp*real(ncell, wp), "the limiter is active in many mass cells")
         if (allocated(error)) return
         write (error_unit, "(a,a,a,es10.3,a,es10.3,a,f6.1,a)") "  [recon] corner bound violation (", &
            merge("smooth state", "masked state", k == 1), "): mass ", worst_m, ", tracers ", worst_t, &
            " (relative); limiter active in ", 100.0_wp*real(nlim, wp)/real(max(ncell, 1), wp), " % of mass cells"
         call check(error, worst_m <= TOL .and. worst_t <= TOL, "corner values within the substituted 3x3 bounds")
         if (allocated(error)) return
      end do
      call rig_destroy(r)
   end subroutine test_bounds

   subroutine test_type2_mask_hand(error)
      !! A type-2 tracer (`qice`-like child of `hi`) in a cell whose WEST neighbour has ice (mass mask 1)
      !! but zero parent thickness (parent mask 0) and a finite, very different child value: the
      !! child gradient substitutes that neighbour by the cell's own value, scalar hand computation.
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg
      real(wp) :: mm, hi0, mtx, mty, tcp, txp, typ, gx, gy, pmn, pmx, qmn, qmx, a1, a2, alpha, etx, ety, etc, c0
      real(wp) :: massmask_tx

      call table_nil(2, tab)
      call rig_init(r, 3, 3, 1, tab)
      r%ws%aim = 0.0_wp
      r%ws%trm = 0.0_wp
      r%ws%aim(:, :, 1) = 0.4_wp
      r%ws%aim(3, 2, 1) = 0.6_wp
      r%ws%trm(:, :, 1, 1) = 1.0_wp
      r%ws%trm(1, 2, 1, 1) = 0.0_wp      ! west: ice, but no parent thickness
      r%ws%trm(3, 2, 1, 1) = 1.3_wp
      c0 = -2.0e8_wp
      r%ws%trm(:, :, 4, 1) = c0
      r%ws%trm(1, 2, 4, 1) = -9.0e8_wp   ! west child value, finite and very different
      r%ws%trm(3, 2, 4, 1) = -2.4e8_wp   ! east
      r%ws%trm(2, 3, 4, 1) = -2.2e8_wp   ! north
      r%ws%trm(2, 1, 4, 1) = -1.8e8_wp   ! south
      r%hm = 1.0_wp
      call rig_push(r)
      call rig_run(r, 2, 2, 2, 2, nbad, cfg)
      call rig_pull(r)
      call check(error, cfg, RECON_CFG_OK)
      if (allocated(error)) return
      call check(error, int(r%rw%tmask(1, 2, 1, 1)) == 0 .and. int(r%rw%tmask(2, 2, 1, 1)) == 1 .and. &
                 int(r%rw%mmask(1, 2, 1)) == 1, "west: mass mask 1, parent mask 0")
      if (allocated(error)) return
      ! centroid of area x parent in the cell, from the parent's own outputs
      mm = 0.4_wp; hi0 = 1.0_wp
      tcp = r%rw%tc(2, 2, 1, 1); txp = r%rw%tx(2, 2, 1, 1); typ = r%rw%ty(2, 2, 1, 1)
      mtx = (r%rw%mc(2, 2, 1)*txp + r%rw%mx(2, 2, 1)*tcp)*RECON_XXAV/(mm*hi0)
      mty = (r%rw%mc(2, 2, 1)*typ + r%rw%my(2, 2, 1)*tcp)*RECON_XXAV/(mm*hi0)
      ! the west value is replaced by the centre: gx = (e - c)/2, gy = (n - s)/2
      gx = 0.5_wp*(-2.4e8_wp - c0)
      gy = 0.5_wp*(-2.2e8_wp - (-1.8e8_wp))
      pmn = (-2.4e8_wp) - c0
      pmx = (-1.8e8_wp) - c0
      qmn = min((0.5_wp - mtx)*gx + (0.5_wp - mty)*gy, (0.5_wp - mtx)*gx - (0.5_wp + mty)*gy, &
                -(0.5_wp + mtx)*gx - (0.5_wp + mty)*gy, (0.5_wp - mty)*gy - (0.5_wp + mtx)*gx)
      qmx = max((0.5_wp - mtx)*gx + (0.5_wp - mty)*gy, (0.5_wp - mtx)*gx - (0.5_wp + mty)*gy, &
                -(0.5_wp + mtx)*gx - (0.5_wp + mty)*gy, (0.5_wp - mty)*gy - (0.5_wp + mtx)*gx)
      a1 = 1.0_wp; a2 = 1.0_wp
      if (abs(qmn) > abs(pmn)) a1 = max(0.0_wp, pmn/qmn)
      if (abs(qmx) > abs(pmx)) a2 = max(0.0_wp, pmx/qmx)
      alpha = min(a1, a2)
      etx = alpha*gx; ety = alpha*gy
      etc = c0 - etx*mtx - ety*mty
      call check(error, abs(r%rw%tx(2, 2, 4, 1) - etx) <= 1.0e-6_wp .and. abs(r%rw%ty(2, 2, 4, 1) - ety) <= 1.0e-6_wp &
                 .and. abs(r%rw%tc(2, 2, 4, 1) - etc) <= 1.0e-6_wp, "child gradient and centre value (hand computation)")
      if (allocated(error)) return
      ! what the MASS mask would give: the west value -9e8 used, gx = (e - w)/2 = +3.3e8
      massmask_tx = 0.5_wp*(-2.4e8_wp - (-9.0e8_wp))
      call check(error, abs(r%rw%tx(2, 2, 4, 1) - massmask_tx) > 1.0e8_wp, "differs from masking with the mass mask")
      call rig_destroy(r)
   end subroutine test_type2_mask_hand

   subroutine test_hm_contract(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg

      call table_nil(2, tab)
      call rig_init(r, 8, 8, 1, tab)
      call fill_smooth_positive(r, 5_int64)
      r%hm(1, 8) = 0.5_wp      ! on the read ring corner of the range 2:7
      call rig_push(r)
      r%rw%mc = 7.0_wp
      !$acc update device(r%rw%mc)
      call rig_run(r, 2, 7, 2, 7, nbad, cfg)
      !$acc update self(r%rw%mc)
      call check(error, cfg, RECON_CFG_MASK, "a land mask value of 0.5 is refused")
      if (allocated(error)) return
      call check(error, all(r%rw%mc == 7.0_wp), "nothing reconstructed on a refusal")
      if (allocated(error)) return
      r%hm(1, 8) = 1.0_wp
      call rig_push(r)
      call rig_run(r, 2, 7, 2, 7, nbad, cfg)
      call check(error, cfg, RECON_CFG_OK, "0/1 accepted")
      call rig_destroy(r)
   end subroutine test_hm_contract

   subroutine test_mean(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, i, j, n, nt, p, sx, sy, nchk1, nchk2
      real(wp) :: m, t, tp, s1, s2, a1, a2, xx, yy, mm, w1, w2, res_m, res1, res2

      call table_nil(4, tab)
      call rig_init(r, 34, 30, 3, tab)
      call fill_random(r, 99_int64, 4)
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      call check(error, cfg, RECON_CFG_OK)
      if (allocated(error)) return
      res_m = 0.0_wp; res1 = 0.0_wp; res2 = 0.0_wp; nchk1 = 0; nchk2 = 0
      do n = 1, r%ncat
         do j = 2, r%ny - 1
            do i = 2, r%nx - 1
               mm = r%ws%aim(i, j, n)
               if (.not. mm > PUNY) cycle
               s1 = 0.0_wp
               do sy = -1, 1, 2
                  do sx = -1, 1, 2
                     xx = real(sx, wp)*GAUSS; yy = real(sy, wp)*GAUSS
                     s1 = s1 + 0.25_wp*(r%rw%mc(i, j, n) + r%rw%mx(i, j, n)*xx + r%rw%my(i, j, n)*yy)
                  end do
               end do
               call acc_max(res_m, abs(s1 - mm)/mm)
               do nt = 1, r%ntrace
                  t = r%ws%trm(i, j, nt, n)
                  s1 = 0.0_wp; a1 = 0.0_wp; s2 = 0.0_wp; a2 = 0.0_wp
                  p = r%ws%tab%depend(nt)
                  do sy = -1, 1, 2
                     do sx = -1, 1, 2
                        xx = real(sx, wp)*GAUSS; yy = real(sy, wp)*GAUSS
                        m = r%rw%mc(i, j, n) + r%rw%mx(i, j, n)*xx + r%rw%my(i, j, n)*yy
                        w1 = r%rw%tc(i, j, nt, n) + r%rw%tx(i, j, nt, n)*xx + r%rw%ty(i, j, nt, n)*yy
                        if (r%ws%tab%tracer_type(nt) == 1) then
                           s1 = s1 + 0.25_wp*m*w1
                           a1 = a1 + 0.25_wp*abs(m*w1)
                        else if (r%ws%aim(i, j, n) > PUNY .and. abs(r%ws%trm(i, j, p, n)) > PUNY) then
                           w2 = r%rw%tc(i, j, p, n) + r%rw%tx(i, j, p, n)*xx + r%rw%ty(i, j, p, n)*yy
                           s2 = s2 + 0.25_wp*m*w2*w1
                           a2 = a2 + 0.25_wp*abs(m*w2*w1)
                        end if
                     end do
                  end do
                  if (r%ws%tab%tracer_type(nt) == 1) then
                     call acc_max(res1, abs(s1 - mm*t)/max(a1, tiny(1.0_wp)))
                     nchk1 = nchk1 + 1
                  else if (abs(r%ws%trm(i, j, p, n)) > PUNY) then
                     tp = r%ws%trm(i, j, p, n)
                     call acc_max(res2, abs(s2 - mm*tp*t)/max(a2, tiny(1.0_wp)))
                     nchk2 = nchk2 + 1
                  end if
               end do
            end do
         end do
      end do
      write (error_unit, "(a,es10.3,a,es10.3,a,es10.3,a,i0,a,i0,a)") "  [recon] mean residuals: <m>=aim ", res_m, &
         ", <mT>=aim tm ", res1, ", <m Tp T>=aim tp t ", res2, " (", nchk1, " type-1 and ", nchk2, " type-2 samples)"
      call check(error, nchk1 > 1000 .and. nchk2 > 1000, "enough samples")
      if (allocated(error)) return
      call check(error, res_m <= TOL, "mass mean")
      if (allocated(error)) return
      call check(error, res1 <= TOL, "type-1 mean")
      if (allocated(error)) return
      call check(error, res2 <= TOL, "type-2 mean")
      call rig_destroy(r)
   end subroutine test_mean

   subroutine test_linear(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, i, j, nt
      real(wp) :: ax, ay, bx(12), by(12), worst, mxav, myav, wc
      real(wp) :: wmx, wtx

      call table_nil(4, tab)
      call rig_init(r, 20, 18, 1, tab)
      ax = 0.011_wp; ay = -0.017_wp
      do nt = 1, r%ntrace
         bx(nt) = 0.013_wp*real(nt, wp)*(-1.0_wp)**nt
         by(nt) = 0.007_wp*real(nt + 1, wp)*(-1.0_wp)**(nt/2)
      end do
      do j = 1, r%ny
         do i = 1, r%nx
            r%ws%aim(i, j, 0) = 0.4_wp
            r%ws%aim(i, j, 1) = 0.5_wp + ax*real(i, wp) + ay*real(j, wp)
            do nt = 1, r%ntrace
               r%ws%trm(i, j, nt, 1) = real(nt, wp) + 2.0_wp + bx(nt)*real(i, wp) + by(nt)*real(j, wp)
            end do
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      call check(error, cfg, RECON_CFG_OK)
      if (allocated(error)) return
      wmx = 0.0_wp; wtx = 0.0_wp; worst = 0.0_wp
      do j = 2, r%ny - 1
         do i = 2, r%nx - 1
            call acc_max(wmx, abs(r%rw%mx(i, j, 1) - ax))
            call acc_max(wmx, abs(r%rw%my(i, j, 1) - ay))
            mxav = r%rw%mx(i, j, 1)*RECON_XXAV/r%ws%aim(i, j, 1)
            myav = r%rw%my(i, j, 1)*RECON_XXAV/r%ws%aim(i, j, 1)
            do nt = 1, r%ntrace
               call acc_max(wtx, abs(r%rw%tx(i, j, nt, 1) - bx(nt)))
               call acc_max(wtx, abs(r%rw%ty(i, j, nt, 1) - by(nt)))
               if (r%ws%tab%tracer_type(nt) == 1) then
                  wc = r%ws%trm(i, j, nt, 1) - bx(nt)*mxav - by(nt)*myav
                  call acc_max(worst, abs(r%rw%tc(i, j, nt, 1) - wc))
               end if
            end do
         end do
      end do
      write (error_unit, "(a,es10.3,a,es10.3,a,es10.3)") "  [recon] linear field: mass gradient error ", wmx, &
         ", tracer gradient error ", wtx, ", type-1 centre error ", worst
      call check(error, wmx <= 10.0_wp*TOL*abs(ax) .and. wtx <= 10.0_wp*TOL, "gradients exact")
      if (allocated(error)) return
      call check(error, worst <= 10.0_wp*TOL, "centre values")
      call rig_destroy(r)
   end subroutine test_linear

   subroutine test_constant(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, i, j, n, nt, nbadcell

      call table_nil(4, tab)
      call rig_init(r, 30, 26, 2, tab)
      call fill_smooth_positive(r, 7_int64)
      do n = 1, r%ncat
         do nt = 1, r%ntrace
            r%ws%trm(:, :, nt, n) = 0.37_wp*real(nt, wp) + 1.0_wp
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      nbadcell = 0
      do n = 1, r%ncat
         do nt = 1, r%ntrace
            do j = 2, r%ny - 1
               do i = 2, r%nx - 1
                  if (r%rw%tx(i, j, nt, n) /= 0.0_wp .or. r%rw%ty(i, j, nt, n) /= 0.0_wp) nbadcell = nbadcell + 1
                  if (r%rw%tc(i, j, nt, n) /= r%ws%trm(i, j, nt, n)) nbadcell = nbadcell + 1
               end do
            end do
         end do
      end do
      call check(error, nbadcell, 0, "constant tracers: zero gradient, centre = mean, exactly")
      call rig_destroy(r)
   end subroutine test_constant

   subroutine test_masks_hand(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg
      real(wp) :: alpha, gx, gy, cnx, cny, qmn, expect_tx, expect_ty, expect_tc, nosub_mx

      ! A: the WEST neighbour is land (hm = 0) holding a wild area: the mass gradient substitutes the
      ! cell's own value 0.4, so gx = (0.6 - 0.4)/2 = 0.1, gy = (0.5 - 0.3)/2 = 0.1; pmn = -0.1,
      ! pmx = 0.7 - 0.4 = 0.3; corners +-0.05 +-0.05 so qmn = -0.1 and |qmn| = |pmn|: no limiting.
      call table_nil(4, tab)
      call rig_init(r, 3, 3, 1, tab)
      r%ws%aim = 0.0_wp
      r%ws%trm = 0.0_wp
      r%ws%aim(:, :, 1) = 0.4_wp
      r%ws%aim(3, 2, 1) = 0.6_wp   ! east
      r%ws%aim(2, 3, 1) = 0.5_wp   ! north
      r%ws%aim(2, 1, 1) = 0.3_wp   ! south
      r%ws%aim(3, 3, 1) = 0.7_wp   ! north-east
      r%ws%aim(1, 2, 1) = 9.9_wp   ! west: land with junk
      r%hm = 1.0_wp
      r%hm(1, 2) = 0.0_wp
      call rig_push(r)
      call rig_run(r, 2, 2, 2, 2, nbad, cfg)
      call rig_pull(r)
      call check(error, cfg, RECON_CFG_OK)
      if (allocated(error)) return
      nosub_mx = 0.5_wp*(0.6_wp - 9.9_wp)
      call check(error, abs(r%rw%mx(2, 2, 1) - 0.1_wp) < 1.0e-15_wp .and. abs(r%rw%my(2, 2, 1) - 0.1_wp) < 1.0e-15_wp, &
                 "land neighbour replaced by the cell's own value (hand computation)")
      if (allocated(error)) return
      call check(error, abs(nosub_mx - r%rw%mx(2, 2, 1)) > 1.0_wp, "...which is very different from using the junk value")
      call rig_destroy(r)
      if (allocated(error)) return

      ! B: the WEST neighbour is an ice-free OCEAN cell (area 0, tracer 0): for the mass it counts as
      ! zero: gx = (0.6 - 0)/2 = 0.3, gy = 0.1, no limiting (corners +-0.15 +-0.05 inside [-0.4, 0.3]
      ! only the negative one... qmn = -0.2 < |pmn| = 0.4, qmx = 0.2 < pmx = 0.3).  The ice thickness
      ! hi (a type-1 tracer with dependents) is masked by the AREA mask: the west cell carries no ice, so
      ! its hi is replaced by the cell's own 1.0.  East 1.4, north 1.2, south 0.8: gx = gy = 0.2,
      ! pmn = -0.2, pmx = +0.4.  The centre of mass is (mx, my)/(12 mm) and the limiter
      ! alpha = pmn/qmn with qmn = -(0.5 + cnx) gx - (0.5 + cny) gy.
      call rig_init(r, 3, 3, 1, tab)
      r%ws%aim = 0.0_wp
      r%ws%trm = 0.0_wp
      r%ws%aim(:, :, 1) = 0.4_wp
      r%ws%aim(3, 2, 1) = 0.6_wp
      r%ws%aim(2, 3, 1) = 0.5_wp
      r%ws%aim(2, 1, 1) = 0.3_wp
      r%ws%aim(3, 3, 1) = 0.7_wp
      r%ws%aim(1, 2, 1) = 0.0_wp
      r%ws%trm(:, :, 1, 1) = 1.0_wp
      r%ws%trm(1, 2, 1, 1) = 0.0_wp
      r%ws%trm(3, 2, 1, 1) = 1.4_wp
      r%ws%trm(2, 3, 1, 1) = 1.2_wp
      r%ws%trm(2, 1, 1, 1) = 0.8_wp
      r%hm = 1.0_wp
      call rig_push(r)
      call rig_run(r, 2, 2, 2, 2, nbad, cfg)
      call rig_pull(r)
      call check(error, abs(r%rw%mx(2, 2, 1) - 0.3_wp) < 1.0e-15_wp .and. abs(r%rw%my(2, 2, 1) - 0.1_wp) < 1.0e-15_wp, &
                 "ice-free ocean neighbour counts as zero for the mass")
      if (allocated(error)) return
      cnx = 0.3_wp/(12.0_wp*0.4_wp)
      cny = 0.1_wp/(12.0_wp*0.4_wp)
      gx = 0.2_wp; gy = 0.2_wp
      qmn = -(0.5_wp + cnx)*gx - (0.5_wp + cny)*gy
      alpha = 0.2_wp/abs(qmn)
      expect_tx = alpha*gx
      expect_ty = alpha*gy
      expect_tc = 1.0_wp - expect_tx*cnx - expect_ty*cny
      call check(error, alpha < 1.0_wp .and. alpha > 0.9_wp, "hand case limits (alpha in (0.9, 1))")
      if (allocated(error)) return
      call check(error, abs(r%rw%tx(2, 2, 1, 1) - expect_tx) < 1.0e-14_wp .and. &
                 abs(r%rw%ty(2, 2, 1, 1) - expect_ty) < 1.0e-14_wp .and. abs(r%rw%tc(2, 2, 1, 1) - expect_tc) < 1.0e-14_wp, &
                 "tracer: ice-free neighbour replaced by the cell's own value, limited (hand computation)")
      if (allocated(error)) return
      ! Without the substitution the west value 0 would give gx = 0.7: a different gradient
      call check(error, abs(r%rw%tx(2, 2, 1, 1) - 0.5_wp*(1.4_wp - 0.0_wp)) > 0.1_wp, "...unlike using the empty cell's zero")
      call rig_destroy(r)
   end subroutine test_masks_hand

   subroutine test_stretched(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, i, j
      real(wp) :: rr, dxi, e_geo, worst_lin, e_pred, worst_geo, ycen, dyj, ey
      real(wp), allocatable :: xcen(:), wid(:)

      call table_nil(2, tab)
      call rig_init(r, 24, 16, 1, tab)
      allocate (xcen(r%nx), wid(r%nx))
      ! (a) linearly varying widths dx_i = 1 + 0.03 i: the centre positions are quadratic in i, the
      ! centred difference of the scaled field equals the cell width EXACTLY.
      ! centres of cells whose edges are X(i) = i + 0.015 i (i+1): c(i) = (X(i-1) + X(i))/2
      do i = 1, r%nx
         xcen(i) = 0.5_wp*((real(i - 1, wp) + 0.015_wp*real((i - 1)*i, wp)) + (real(i, wp) + 0.015_wp*real(i*(i + 1), wp)))
         wid(i) = (real(i, wp) + 0.015_wp*real(i*(i + 1), wp)) - (real(i - 1, wp) + 0.015_wp*real((i - 1)*i, wp))
      end do
      r%ws%aim = 0.5_wp
      do j = 1, r%ny
         do i = 1, r%nx
            r%ws%trm(i, j, 1, 1) = 3.0_wp + xcen(i)
            r%ws%trm(i, j, 2, 1) = 1.0_wp
            r%ws%trm(i, j, 3:, 1) = 2.0_wp
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      worst_lin = 0.0_wp
      do j = 2, r%ny - 1
         do i = 2, r%nx - 1
            call acc_max(worst_lin, abs(r%rw%tx(i, j, 1, 1)/wid(i) - 1.0_wp))
         end do
      end do
      write (error_unit, "(a,es10.3)") "  [recon] linearly varying widths: worst |d phi/dx - 1| = ", worst_lin
      call check(error, worst_lin <= 10.0_wp*TOL, "linearly varying widths: physical gradient exact")
      if (allocated(error)) return

      ! (b) geometric stretching r = 1.05: edges X(i) = r^i, centres r^i (r+1)/2, widths r^i (r-1); the
      ! scaled centred difference over the width is (r+1)^2/(4 r) -> error (r-1)^2/(4 r)
      rr = 1.05_wp
      do i = 1, r%nx
         xcen(i) = rr**(i - 1)*(rr + 1.0_wp)*0.5_wp
         wid(i) = rr**(i - 1)*(rr - 1.0_wp)
      end do
      do j = 1, r%ny
         do i = 1, r%nx
            r%ws%trm(i, j, 1, 1) = 3.0_wp + xcen(i)
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      e_pred = (rr - 1.0_wp)**2/(4.0_wp*rr)
      worst_geo = 0.0_wp
      do j = 2, r%ny - 1
         do i = 2, r%nx - 1
            dxi = wid(i)
            e_geo = r%rw%tx(i, j, 1, 1)/dxi - 1.0_wp
            call acc_max(worst_geo, abs(e_geo - e_pred))
         end do
      end do
      write (error_unit, "(a,es10.3,a,es10.3)") "  [recon] geometric stretching 1.05: predicted error (r-1)^2/(4r) = ", e_pred, &
         ", worst departure from it ", worst_geo
      call check(error, worst_geo <= 1.0e-12_wp, "stretched grid: the closed-form second-order error")
      if (allocated(error)) return

      ! (c) latitude rows of a sector (uniform dlat, dx ~ cos(lat) does not enter): a field linear in
      ! the row index is exact in y
      do j = 1, r%ny
         ycen = 10.0_wp + 0.25_wp*real(j, wp)
         dyj = 0.25_wp
         do i = 1, r%nx
            r%ws%trm(i, j, 1, 1) = ycen
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call rig_pull(r)
      ey = 0.0_wp
      do j = 2, r%ny - 1
         do i = 2, r%nx - 1
            call acc_max(ey, abs(r%rw%ty(i, j, 1, 1)/dyj - 1.0_wp))
         end do
      end do
      call check(error, ey <= 10.0_wp*TOL, "uniform latitude spacing: exact")
      call rig_destroy(r)
   end subroutine test_stretched

   subroutine test_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, k
      real(wp) :: bad

      call table_nil(2, tab)
      call rig_init(r, 12, 12, 2, tab)
      do k = 1, 5
         call fill_smooth_positive(r, 31_int64)
         bad = ieee_value(1.0_wp, ieee_quiet_nan)
         if (k == 5) bad = ieee_value(1.0_wp, ieee_positive_inf)
         select case (k)
         case (1)
            r%ws%aim(11, 6, 1) = bad          ! ring of the mask range
         case (2)
            r%ws%trm(6, 7, 2, 2) = bad
         case (3)
            r%hm(5, 5) = bad
         case (4)
            r%ws%aim(6, 6, 0) = bad
         case (5)
            r%ws%trm(3, 3, 1, 1) = bad
         end select
         call rig_push(r)
         call rig_run(r, 2, 11, 2, 11, nbad, cfg)
         call check(error, cfg, RECON_CFG_OK)
         if (allocated(error)) return
         call check(error, nbad > 0, "non-finite input is counted")
         if (allocated(error)) return
      end do
      call fill_smooth_positive(r, 31_int64)
      call rig_push(r)
      call rig_run(r, 2, 11, 2, 11, nbad, cfg)
      call check(error, nbad, 0, "clean input counts nothing")
      call rig_destroy(r)
   end subroutine test_nan

   subroutine test_oracle(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      type(orc_t) :: o
      integer :: nbad, cfg, is, ie, js, je, ic, nmis
      integer(int64) :: seed
      real(wp) :: dmc, dmx, dtc, dtx, sc_m, sc_t, sc_g, wm, wt
      logical :: ring_ok

      call table_nil(4, tab)
      call rig_init(r, 29, 25, 3, tab)
      wm = 0.0_wp; wt = 0.0_wp
      do ic = 1, 4
         seed = 1000_int64 + 17_int64*int(ic, int64)
         call fill_random(r, seed, 4)
         if (ic >= 3) call fill_random(r, seed + 5_int64, 4)
         is = 2; ie = r%nx - 1; js = 2; je = r%ny - 1
         if (ic == 2 .or. ic == 4) then
            is = 4; ie = r%nx - 4; js = 3; je = r%ny - 5
         end if
         call sentinel_outputs(r)
         call rig_push(r)
         call rig_run(r, is, ie, js, je, nbad, cfg)
         call rig_pull(r)
         call check(error, cfg, RECON_CFG_OK)
         if (allocated(error)) return
         call check(error, nbad, 0, "no non-finite input")
         if (allocated(error)) return
         call oracle(r, is, ie, js, je, o)
         ! decisions: the masks, exactly, on their range
         nmis = count(r%rw%mmask(is - 1:ie + 1, js - 1:je + 1, :) /= o%mmask(is - 1:ie + 1, js - 1:je + 1, :)) &
                + count(r%rw%tmask(is - 1:ie + 1, js - 1:je + 1, :, :) /= o%tmask(is - 1:ie + 1, js - 1:je + 1, :, :))
         call check(error, nmis, 0, "masks equal the oracle exactly")
         if (allocated(error)) return
         sc_m = max(maxval(abs(o%mc)), 1.0_wp)
         sc_t = max(maxval(abs(r%ws%trm)), 1.0_wp)
         sc_g = sc_t
         dmc = maxdiff3(r%rw%mc(is:ie, js:je, :), o%mc(is:ie, js:je, :), sc_m)
         dmx = max(maxdiff3(r%rw%mx(is:ie, js:je, :), o%mx(is:ie, js:je, :), sc_m), &
                   maxdiff3(r%rw%my(is:ie, js:je, :), o%my(is:ie, js:je, :), sc_m))
         dtc = maxdiff4(r%rw%tc(is:ie, js:je, :, :), o%tc(is:ie, js:je, :, :), sc_t)
         dtx = max(maxdiff4(r%rw%tx(is:ie, js:je, :, :), o%tx(is:ie, js:je, :, :), sc_g), &
                   maxdiff4(r%rw%ty(is:ie, js:je, :, :), o%ty(is:ie, js:je, :, :), sc_g))
         call acc_max(wm, dmc)
         call acc_max(wm, dmx)
         call acc_max(wt, dtc)
         call acc_max(wt, dtx)
         call check(error, max(dmc, dmx) <= TOL, "mass outputs equal the oracle")
         if (allocated(error)) return
         call check(error, max(dtc, dtx) <= TOL, "tracer outputs equal the oracle")
         if (allocated(error)) return
         ! the ring outside the range keeps its sentinel
         ring_ok = count(r%rw%mc /= 7.0_wp) == count(r%rw%mc(is:ie, js:je, :) /= 7.0_wp) &
                   .and. count(r%rw%mx /= 7.0_wp) == count(r%rw%mx(is:ie, js:je, :) /= 7.0_wp) &
                   .and. count(r%rw%my /= 7.0_wp) == count(r%rw%my(is:ie, js:je, :) /= 7.0_wp) &
                   .and. count(r%rw%tc /= 7.0_wp) == count(r%rw%tc(is:ie, js:je, :, :) /= 7.0_wp) &
                   .and. count(r%rw%tx /= 7.0_wp) == count(r%rw%tx(is:ie, js:je, :, :) /= 7.0_wp) &
                   .and. count(r%rw%ty /= 7.0_wp) == count(r%rw%ty(is:ie, js:je, :, :) /= 7.0_wp) &
                   .and. count(r%rw%mmask /= 5_int8) == count(r%rw%mmask(is - 1:ie + 1, js - 1:je + 1, :) /= 5_int8) &
                   .and. count(r%rw%tmask /= 5_int8) == count(r%rw%tmask(is - 1:ie + 1, js - 1:je + 1, :, :) /= 5_int8)
         call check(error, ring_ok, "outputs outside the range are untouched")
         if (allocated(error)) return
         call check(error, count(r%rw%mc(is:ie, js:je, :) == 7.0_wp) + count(r%rw%tx(is:ie, js:je, :, :) == 7.0_wp), 0, &
                    "every output on the range is written")
         if (allocated(error)) return
      end do
      write (error_unit, "(a,es10.3,a,es10.3,a)") "  [recon] oracle, 4 random masked states: worst relative difference mass ", &
         wm, ", tracers ", wt, " (tolerance "//merge("1e-12", "1e-14", TOL > 1.0e-13_wp)//")"
      call rig_destroy(r)
   end subroutine test_oracle

   subroutine test_python(error)
      !! The numbers of /python_prototypes/cice_remap/recon_transcript.py (a stdlib-only
      !! transcript written corner by corner): four cells, five fields each.
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, ic, k, i, j, c
      real(wp) :: got(17), worst, sc
      integer :: dep(2)

      dep = [0, 1]
      call transport_tracer_table_from_depend(2, dep, 0, 0, tab)
      call check(error, tab%ierr, TT_OK)
      if (allocated(error)) return
      worst = 0.0_wp
      do ic = 1, 4
         call rig_init(r, 3, 3, 1, tab)
         r%ws%aim = 0.0_wp
         r%ws%trm = 0.0_wp
         do j = 1, 3
            do i = 1, 3
               k = (j - 1)*3 + i
               r%ws%aim(i, j, 1) = PY_AIM(k, ic)
               r%hm(i, j) = PY_HM(k, ic)
               r%ws%trm(i, j, 1, 1) = PY_HI(k, ic)
               r%ws%trm(i, j, 2, 1) = PY_HS(k, ic)
               r%ws%trm(i, j, 3, 1) = PY_T3(k, ic)
               r%ws%trm(i, j, 4, 1) = PY_T4(k, ic)
            end do
         end do
         call rig_push(r)
         call rig_run(r, 2, 2, 2, 2, nbad, cfg)
         call rig_pull(r)
         call check(error, cfg, RECON_CFG_OK)
         if (allocated(error)) return
         got = [r%rw%mc(2, 2, 1), r%rw%mx(2, 2, 1), r%rw%my(2, 2, 1), &
                r%rw%tc(2, 2, 1, 1), r%rw%tx(2, 2, 1, 1), r%rw%ty(2, 2, 1, 1), 0.0_wp, 0.0_wp, &
                r%rw%tc(2, 2, 2, 1), r%rw%tx(2, 2, 2, 1), r%rw%ty(2, 2, 2, 1), &
                r%rw%tc(2, 2, 3, 1), r%rw%tx(2, 2, 3, 1), r%rw%ty(2, 2, 3, 1), &
                r%rw%tc(2, 2, 4, 1), r%rw%tx(2, 2, 4, 1), r%rw%ty(2, 2, 4, 1)]
         ! the parent centroid (entries 7, 8) is not an output: re-form it from the outputs
         if (r%rw%tmask(2, 2, 1, 1) /= 0_int8) then
            got(7) = (r%rw%mc(2, 2, 1)*r%rw%tx(2, 2, 1, 1) + r%rw%mx(2, 2, 1)*r%rw%tc(2, 2, 1, 1))*RECON_XXAV &
                     /(r%ws%aim(2, 2, 1)*r%ws%trm(2, 2, 1, 1))
            got(8) = (r%rw%mc(2, 2, 1)*r%rw%ty(2, 2, 1, 1) + r%rw%my(2, 2, 1)*r%rw%tc(2, 2, 1, 1))*RECON_XXAV &
                     /(r%ws%aim(2, 2, 1)*r%ws%trm(2, 2, 1, 1))
         end if
         do c = 1, 17
            sc = max(abs(PY_OUT(c, ic)), 1.0e-3_wp)
            if (c >= 15) sc = 3.0e8_wp
            if (c >= 12 .and. c <= 14) sc = 10.0_wp
            if (.not. ieee_is_finite(got(c))) then
               worst = huge(1.0_wp)
            else
               call acc_max(worst, abs(got(c) - PY_OUT(c, ic))/sc)
            end if
         end do
         call rig_destroy(r)
      end do
      write (error_unit, "(a,es10.3)") "  [recon] Python transcript, 4 cells x 17 numbers: worst scaled difference ", worst
      call check(error, worst <= TOL, "Python transcript")
   end subroutine test_python

   subroutine test_cfg(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(remap_recon_t) :: rw2, rw0
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, ierr, keep

      call table_nil(2, tab)
      call rig_init(r, 10, 10, 2, tab)
      call fill_smooth_positive(r, 3_int64)
      call rig_push(r)
      call rig_run(r, 1, 9, 2, 9, nbad, cfg)
      call check(error, cfg, RECON_CFG_RANGE, "is = 1 leaves no input ring")
      if (allocated(error)) return
      call rig_run(r, 2, 10, 2, 9, nbad, cfg)
      call check(error, cfg, RECON_CFG_RANGE, "ie = nx leaves no input ring")
      if (allocated(error)) return
      call rig_run(r, 5, 4, 2, 9, nbad, cfg)
      call check(error, cfg, RECON_CFG_RANGE, "empty range")
      if (allocated(error)) return
      call rw2%init(10, 11, 2, tab%ntrace, ierr)
      call remap_reconstruct(r%ws, rw2, 2, 9, 2, 9, PUNY, r%hm, nbad, cfg)
      call check(error, cfg, RECON_CFG_WORKSPACE, "extent mismatch")
      if (allocated(error)) return
      call remap_reconstruct(r%ws, rw0, 2, 9, 2, 9, PUNY, r%hm, nbad, cfg)
      call check(error, cfg, RECON_CFG_WORKSPACE, "uninitialised workspace")
      if (allocated(error)) return
      call rw0%init(2, 2, 1, 1, ierr)
      call check(error, ierr /= 0, "too small an extent is refused")
      if (allocated(error)) return
      keep = r%ws%tab%tracer_type(3)
      r%ws%tab%tracer_type(3) = 3
      call rig_run(r, 2, 9, 2, 9, nbad, cfg)
      call check(error, cfg, RECON_CFG_TABLE, "type 3 is refused")
      if (allocated(error)) return
      r%ws%tab%tracer_type(3) = keep
      keep = r%ws%tab%depend(5)
      r%ws%tab%depend(5) = 5
      call rig_run(r, 2, 9, 2, 9, nbad, cfg)
      call check(error, cfg, RECON_CFG_TABLE, "a parent that is not earlier is refused")
      r%ws%tab%depend(5) = keep
      call check(error, r%rw%bytes() > 0_int64, "bytes")
      call rw2%destroy()
      call rig_destroy(r)
   end subroutine test_cfg

   subroutine test_empty(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nbad, cfg, nnz

      call table_nil(3, tab)
      call rig_init(r, 14, 12, 2, tab)
      r%ws%aim = 0.0_wp
      r%ws%trm = 0.0_wp
      call sentinel_outputs(r)
      call rig_push(r)
      call rig_run(r, 2, 13, 2, 11, nbad, cfg)
      call rig_pull(r)
      nnz = count(r%rw%mc(2:13, 2:11, :) /= 0.0_wp) + count(r%rw%mx(2:13, 2:11, :) /= 0.0_wp) &
            + count(r%rw%tc(2:13, 2:11, :, :) /= 0.0_wp) + count(r%rw%tx(2:13, 2:11, :, :) /= 0.0_wp) &
            + count(r%rw%ty(2:13, 2:11, :, :) /= 0.0_wp) + count(r%rw%mmask(1:14, 1:12, :) /= 0_int8) &
            + count(r%rw%tmask(1:14, 1:12, :, :) /= 0_int8)
      call check(error, nnz, 0, "an empty domain reconstructs to zeros and no mask")
      if (allocated(error)) return
      ! all land: ice present but no cell is ocean: the mass gradient is zero, the centre values stay
      call fill_smooth_positive(r, 77_int64)
      r%hm = 0.0_wp
      call rig_push(r)
      call rig_run(r, 2, 13, 2, 11, nbad, cfg)
      call rig_pull(r)
      call check(error, count(r%rw%mx(2:13, 2:11, :) /= 0.0_wp) + count(r%rw%my(2:13, 2:11, :) /= 0.0_wp), 0, &
                 "all land: zero mass gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(r%rw%mc(2:13, 2:11, 1:2) - r%ws%aim(2:13, 2:11, 1:2))) == 0.0_wp, "mc = aim")
      call rig_destroy(r)
   end subroutine test_empty

   subroutine test_throughput(error)
      type(error_type), allocatable, intent(out) :: error
      type(rig_t) :: r
      type(transport_tracer_table_t) :: tab
      integer :: nxp, nyp, rep, nbad, cfg, i, j, n, nt
      integer, parameter :: NREP = 5
      integer(int64) :: t0, t1, rate
      real(wp) :: tc, cells

#ifdef RDB_GPU_OFFLOAD
      nxp = 1024; nyp = 1024
#else
      nxp = 512; nyp = 512
#endif
      call table_nil(7, tab)
      call rig_init(r, nxp + 4, nyp + 4, 5, tab)
      do n = 0, 5
         do j = 1, r%ny
            do i = 1, r%nx
               r%ws%aim(i, j, n) = 0.04_wp + 0.01_wp*real(mod(i + 2*j + n, 7), wp)
            end do
         end do
      end do
      do n = 1, 5
         do nt = 1, r%ntrace
            do j = 1, r%ny
               do i = 1, r%nx
                  r%ws%trm(i, j, nt, n) = real(nt, wp)*(1.0_wp + 0.1_wp*real(mod(i + j + n + nt, 5), wp))
               end do
            end do
         end do
      end do
      call rig_push(r)
      call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      call system_clock(t0, rate)
      do rep = 1, NREP
         call rig_run(r, 2, r%nx - 1, 2, r%ny - 1, nbad, cfg)
      end do
      call system_clock(t1)
      tc = real(t1 - t0, wp)/real(rate, wp)/real(NREP, wp)
      cells = real(r%nx, wp)*real(r%ny, wp)*5.0_wp
      write (error_unit, "(a,i0,a,i0,a,i0,a,es10.3,a,f8.2,a,f7.2,a)") "  [recon] ", nxp, "x", nyp, "x5 cats, ntrace ", r%ntrace, &
         ": ", cells/tc, " cat-cells/s (", tc*1000.0_wp, " ms for 2 launches), buffers ", &
         real(r%rw%bytes(), wp)/1.0e9_wp, " GB"
      call check(error, nbad == 0 .and. cfg == RECON_CFG_OK .and. ieee_is_finite(cells/tc), "clean")
      call rig_destroy(r)
   end subroutine test_throughput

end module test_cice_remap_recon
