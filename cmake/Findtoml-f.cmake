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
my_fetch_package("${_lib}" "${_url}" "${_rev}")
set(BUILD_SHARED_LIBS "${_rdb_saved_shared}")
set(BUILD_TESTING "${_rdb_saved_testing}")

unset(_rdb_saved_shared)
unset(_rdb_saved_testing)
unset(_lib)
unset(_url)
unset(_rev)
