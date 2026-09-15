# gcov/lcov coverage report target (GNU/gfortran only). A function: the only
# things that escape its scope are the `coverage` target and the find_program()
# cache entries, both global.
function(rdb_add_coverage_target)
  find_program(LCOV_EXECUTABLE lcov)
  find_program(GENHTML_EXECUTABLE genhtml)
  if(LCOV_EXECUTABLE AND GENHTML_EXECUTABLE)
    # `--ignore-errors mismatch,unused` were added in lcov 2.x; on lcov 1.x
    # those keys are unknown arguments and the whole capture step fails. Detect
    # the major version once and gate the flags accordingly.
    execute_process(
      COMMAND ${LCOV_EXECUTABLE} --version
      OUTPUT_VARIABLE _lcov_version_raw
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "([0-9]+)\\.([0-9]+)" _lcov_version_match
                 "${_lcov_version_raw}")
    set(_lcov_major "${CMAKE_MATCH_1}")
    if(_lcov_major GREATER_EQUAL 2)
      set(_lcov_capture_ignore --ignore-errors mismatch,gcov,source,unused)
      set(_lcov_remove_ignore --ignore-errors unused)
      message(STATUS "Coverage: lcov ${CMAKE_MATCH_0} (2.x flags)")
    else()
      set(_lcov_capture_ignore --ignore-errors gcov,source)
      set(_lcov_remove_ignore "")
      message(STATUS "Coverage: lcov ${CMAKE_MATCH_0} (1.x flags)")
    endif()

    add_custom_target(
      coverage
      COMMAND ${CMAKE_COMMAND} -E remove -f coverage.info coverage_filtered.info
      COMMAND ${LCOV_EXECUTABLE} --directory . --zerocounters
      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
      COMMAND ${LCOV_EXECUTABLE} --directory . --capture --output-file
              coverage.info ${_lcov_capture_ignore}
      COMMAND
        ${LCOV_EXECUTABLE} --remove coverage.info '/usr/*' '/opt/*' '/Library/*'
        '*/build/*' '*/tests/*' --output-file coverage_filtered.info
        ${_lcov_remove_ignore}
      COMMAND ${GENHTML_EXECUTABLE} coverage_filtered.info --output-directory
              coverage_report --ignore-errors source
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Running tests and generating gcov/lcov coverage report")
    message(STATUS "Coverage target enabled: build target `coverage`")
  else()
    message(
      WARNING
        "RDB_ENABLE_COVERAGE=ON but lcov/genhtml not found on PATH "
        "(lcov=${LCOV_EXECUTABLE}, genhtml=${GENHTML_EXECUTABLE}). "
        "Install with `brew install lcov` (macOS) or `apt install lcov` (Linux)."
    )
  endif()
endfunction()
