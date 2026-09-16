!! Central Roundabout namelist schema (strict validation + parameter dumps).
module rdb_config_schema
   !! Re-export shim for Roundabout's strict namelist schema. Re-exports
   !! `build_rdb_schema` (defined in `rdb_config`, kept here to avoid a
   !! circular dependency) and owns the `nml_dirname` helper.
   use rdb_config, only: build_rdb_schema
   implicit none
   private

   public :: build_rdb_schema
   public :: nml_dirname

contains

   pure function nml_dirname(path) result(dir)
      !! Directory part of `path` with a trailing '/', or '' if none.
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: dir
      integer :: i
      dir = ""
      do i = len_trim(path), 1, -1
         if (path(i:i) == "/") then
            dir = path(1:i)
            return
         end if
      end do
   end function nml_dirname

end module rdb_config_schema
