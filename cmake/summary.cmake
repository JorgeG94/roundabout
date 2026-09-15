# Configure-time build summary. A function: it only READS variables, so its own
# scope costs nothing. The "GPU offload:" line is load-bearing — it reports what
# the guards in compiler_flags.cmake resolved RDB_ENABLE_GPU to.
function(rdb_print_build_summary)
  message(STATUS "")
  message(STATUS "=== Roundabout Build Configuration ===")
  message(
    STATUS
      "Compiler:       ${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER})"
  )
  message(STATUS "Build type:     ${CMAKE_BUILD_TYPE}")
  message(STATUS "Fortran flags:  ${CMAKE_Fortran_FLAGS}")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug flags:    ${CMAKE_Fortran_FLAGS_DEBUG}")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Release flags:  ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
  message(STATUS "GPU offload:    ${RDB_ENABLE_GPU}")
  message(STATUS "Double prec:    ${RDB_ENABLE_DOUBLE}")
  message(STATUS "NetCDF I/O:     ${RDB_ENABLE_NETCDF}")
  message(STATUS "Tests:          ${RDB_ENABLE_TESTING}")
  message(STATUS "Benchmarks:     ${RDB_BUILD_BENCHMARKS}")
  message(STATUS "=================================")
  message(STATUS "")
endfunction()
