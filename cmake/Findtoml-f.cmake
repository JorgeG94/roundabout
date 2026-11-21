#
# toml-f: the TOML parser behind the configuration front end
# (src/core/config/rdb_toml_read.F90).
#
# Pinned to JorgeG94/toml-f's chore/get_to_work_nvhpc (e0bc8c6), commits on
# upstream f9ff547 that an nvfortran build needs: -Mbackslash/-Mallocatable=03
# in toml-f's CMake, convert_ucs fixed for nvfortran (a four-byte \U escape came
# back padded with blanks), and a cmake_minimum_required policy range. Move to
# an upstream tag once those are merged.
#
# Built STATIC and without its tests: toml-f's library type follows
# BUILD_SHARED_LIBS, and its tests follow BUILD_TESTING, which would register
# toml-f's own checks in our ctest. Both are set only for this fetch and put
# back afterwards. It inherits CMAKE_POSITION_INDEPENDENT_CODE, which the top
# level CMakeLists.txt turns ON before the dependencies are found whenever the
# shared core (RDB_BUILD_SHARED / RDB_TESTS_LINK_SHARED) is requested.
#
set(_lib "toml-f")
set(_url "https://github.com/JorgeG94/toml-f.git")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
set(_rev "e0bc8c69985adcae9b61af45b8dd18ceea54ef2c")

set(_rdb_saved_shared "${BUILD_SHARED_LIBS}")
set(_rdb_saved_testing "${BUILD_TESTING}")
set(BUILD_SHARED_LIBS OFF)
set(BUILD_TESTING OFF)

# Two nvfortran (26.5) defects make the unmodified library crash inside
# `toml_load` on VALID input; gfortran is correct on both.
#
# UNINITIALISED DESCRIPTORS, at every optimisation level. nvfortran hands
# `move_alloc` / intrinsic assignment a descriptor it never initialised when the
# source is a BLOCK-local allocatable, an allocatable component of a local, or a
# derived-type function result, and its runtime acts on whatever the stack slot
# held. `cmake/patches/toml-f-nvfortran-descriptors.patch` (its header has the
# detail) rewrites the three statements on the load path. It changes nothing on
# gfortran and is applied for every compiler, so there is one toml-f source in
# the tree. The marker below is a string the patch adds.
#
# A DROPPED ARGUMENT, at -O2 and above. In toml-f's `add_array_to_table`
# (src/tomlf/type.f90) the call `table%get(key, tmp)` that follows
# `table%push_back(val, istat)` is emitted WITHOUT loading the hidden length of
# `key` into its argument register, so `get` compares the key at whatever length
# `push_back` left there. When `push_back` grew the table's list -- the 17th,
# 26th, 39th, 59th, 89th ... key of one table -- the lookup fails, `add_array`
# hands the parser a null pointer and `parse_inline_array` segfaults. It is
# register allocation, not source: there is nothing to patch and nothing says it
# stops at that one procedure, so nvfortran builds the whole library at -O1.
# toml-f only parses the input once per run; -O1 costs nothing measurable.
# `-fast` has to go with `-O3` (left on the command line it forces the level
# back up), which is why the per-config flags are replaced for the fetch rather
# than `-O1` appended to the target.
#
# Gates, both in `tests/test_toml_config.F90`: `parse_on_dirty_stack` for the
# descriptors, `key_position_sweep` for the dropped argument.
set(_rdb_tomlf_patch
    "${CMAKE_CURRENT_LIST_DIR}/patches/toml-f-nvfortran-descriptors.patch")
set(_rdb_tomlf_patch_marker "subroutine new_parser(parser, root, config)")
set(_rdb_tomlf_cfgs RELEASE RELWITHDEBINFO MINSIZEREL)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  foreach(_cfg IN LISTS _rdb_tomlf_cfgs)
    set(_rdb_saved_flags_${_cfg} "${CMAKE_Fortran_FLAGS_${_cfg}}")
    set(CMAKE_Fortran_FLAGS_${_cfg} "-O1")
  endforeach()
  message(STATUS "toml-f: built at -O1 (nvfortran miscompiles it at -O2+)")
endif()

my_fetch_package("${_lib}" "${_url}" "${_rev}" "${_rdb_tomlf_patch}"
                 "${_rdb_tomlf_patch_marker}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  foreach(_cfg IN LISTS _rdb_tomlf_cfgs)
    set(CMAKE_Fortran_FLAGS_${_cfg} "${_rdb_saved_flags_${_cfg}}")
    unset(_rdb_saved_flags_${_cfg})
  endforeach()
endif()
set(BUILD_SHARED_LIBS "${_rdb_saved_shared}")
set(BUILD_TESTING "${_rdb_saved_testing}")

unset(_rdb_tomlf_patch)
unset(_rdb_tomlf_patch_marker)
unset(_rdb_tomlf_cfgs)
unset(_cfg)
unset(_rdb_saved_shared)
unset(_rdb_saved_testing)
unset(_lib)
unset(_url)
unset(_rev)
