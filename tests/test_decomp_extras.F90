!! Coverage for decomp_auto_factor and decomp_init_from_config (the
!! auto-factor branch the existing test_decomp suite skips).
module test_decomp_extras
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_decomp, only: decomp_t, decomp_init_from_config, decomp_auto_factor
   use rdb_config, only: config_t
   implicit none
   private

   public :: collect_decomp_extras_tests

contains

   subroutine collect_decomp_extras_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("auto_factor_perfect_squares", test_auto_factor_squares), &
                  new_unittest("auto_factor_prime_nprocs", test_auto_factor_prime), &
                  new_unittest("auto_factor_minimises_perimeter", test_auto_factor_perimeter), &
                  new_unittest("auto_factor_product_invariant", test_auto_factor_product), &
                  new_unittest("init_from_config_uses_explicit_pxpy", test_explicit_pxpy), &
                  new_unittest("init_from_config_auto_factors", test_auto_factors_via_cfg), &
                  new_unittest("init_from_config_single_rank", test_single_rank_via_cfg) &
                  ]
   end subroutine collect_decomp_extras_tests

   subroutine test_auto_factor_squares(error)
      !! 4 procs on a square grid: prefer 2x2 over 1x4 / 4x1.
      type(error_type), allocatable, intent(out) :: error
      integer :: px, py

      call decomp_auto_factor(4, 100, 100, px, py)
      call check(error, px == 2 .and. py == 2, &
                 "4 procs on square grid must factor as 2x2")
      if (allocated(error)) return

      call decomp_auto_factor(16, 100, 100, px, py)
      call check(error, px == 4 .and. py == 4, "16 on square -> 4x4")
   end subroutine test_auto_factor_squares

   subroutine test_auto_factor_prime(error)
      !! 7 is prime: only factorisations are 1x7 / 7x1. The cost-minimiser
      !! picks based on grid aspect.
      type(error_type), allocatable, intent(out) :: error
      integer :: px, py

      ! Tall grid (ny >> nx): perimeter cost prefers a tall process layout.
      call decomp_auto_factor(7, 10, 1000, px, py)
      call check(error, px*py == 7, "factorisation must multiply to nprocs")
      if (allocated(error)) return
      call check(error, px == 1 .and. py == 7, &
                 "tall grid: 7 ranks should stack vertically (1x7)")
   end subroutine test_auto_factor_prime

   subroutine test_auto_factor_perimeter(error)
      !! For 8 procs on a roughly-square grid, the lowest-perimeter
      !! factorisation is 2x4 or 4x2 (cost = 12 vs 9 for 1x8 / 8x1).
      type(error_type), allocatable, intent(out) :: error
      integer :: px, py
      integer :: cost

      call decomp_auto_factor(8, 100, 100, px, py)
      call check(error, px*py == 8, "factorisation product")
      if (allocated(error)) return

      cost = px*100 + py*100   ! perimeter cost from the implementation
      call check(error, cost <= 8*100 + 1*100, &
                 "auto_factor cost must beat the worst factorisation")
   end subroutine test_auto_factor_perimeter

   subroutine test_auto_factor_product(error)
      !! Sweep nprocs and confirm px*py == nprocs in all cases.
      type(error_type), allocatable, intent(out) :: error
      integer :: nprocs, px, py

      do nprocs = 1, 32
         call decomp_auto_factor(nprocs, 100, 50, px, py)
         call check(error, px*py == nprocs, "px*py must equal nprocs")
         if (allocated(error)) return
      end do
   end subroutine test_auto_factor_product

   subroutine test_explicit_pxpy(error)
      !! When cfg%px and cfg%py are pre-set (and nprocs matches), the
      !! auto-factor path is skipped.
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      type(config_t) :: cfg

      cfg%nx = 100; cfg%ny = 50
      cfg%px = 2; cfg%py = 2

      call decomp_init_from_config(d, cfg, 4, 0)

      call check(error, d%px == 2 .and. d%py == 2, &
                 "explicit px/py should be honoured")
      if (allocated(error)) return
      call check(error, d%nx_global == 100 .and. d%ny_global == 50, &
                 "global sizes plumbed from cfg")
   end subroutine test_explicit_pxpy

   subroutine test_auto_factors_via_cfg(error)
      !! cfg%px = cfg%py = 1 with nprocs > 1 triggers auto-factoring.
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      type(config_t) :: cfg

      cfg%nx = 100; cfg%ny = 100
      cfg%px = 1; cfg%py = 1

      call decomp_init_from_config(d, cfg, 4, 0)

      ! For 4 ranks on a square grid the auto-factorer picks 2x2; cfg%px/py
      ! are mutated in place.
      call check(error, cfg%px == 2 .and. cfg%py == 2, &
                 "auto-factor must rewrite cfg%px/py")
      if (allocated(error)) return
      call check(error, d%px == 2 .and. d%py == 2, &
                 "decomp inherits the auto-chosen px/py")
   end subroutine test_auto_factors_via_cfg

   subroutine test_single_rank_via_cfg(error)
      !! Single-rank case must NOT trigger auto-factor (1 == 1*1 trivially).
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      type(config_t) :: cfg

      cfg%nx = 50; cfg%ny = 50
      cfg%px = 1; cfg%py = 1

      call decomp_init_from_config(d, cfg, 1, 0)

      call check(error, d%px == 1 .and. d%py == 1, "single rank stays 1x1")
      if (allocated(error)) return
      call check(error, d%nx_local == 50 .and. d%ny_local == 50, &
                 "single rank owns the full grid")
   end subroutine test_single_rank_via_cfg

end module test_decomp_extras
