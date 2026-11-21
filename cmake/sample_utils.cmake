macro("my_fetch_package" package url)
  # Optional third argument for git reference (branch, tag, or commit)
  set(_git_ref "${ARGV2}")
  if(NOT _git_ref)
    set(_git_ref "HEAD")
  endif()
  # Optional fourth and fifth arguments: a patch file (`patch -p1`) applied to
  # the fetched source right after it is cloned -- for a fix that has not landed
  # upstream -- and a string the patch adds under src/, which tells an
  # already-patched checkout from one that isn't on a reconfigure.
  set(_patch_file "${ARGV3}")
  set(_patch_marker "${ARGV4}")

  string(TOLOWER "${package}" _pkg_lc)
  string(TOUPPER "${package}" _pkg_uc)

  message(STATUS "Retrieving ${package} from ${url} (ref: ${_git_ref})")
  include(FetchContent)

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(FETCHCONTENT_QUIET FALSE)
  endif()

  if(_patch_file)
    if(NOT _patch_marker)
      message(
        FATAL_ERROR
          "my_fetch_package(${package} ... ${_patch_file}): needs a 5th "
          "argument, a string the patch adds under src/, to tell an "
          "already-patched checkout from one that isn't")
    endif()
    if(NOT EXISTS "${_patch_file}")
      message(
        FATAL_ERROR "my_fetch_package(${package}): no patch ${_patch_file}")
    endif()
    find_program(_patch_exe patch REQUIRED)
    find_program(_sh_exe sh REQUIRED)
    # Idempotent: a reconfigure can re-run the patch step on a source tree the
    # first configure already patched, and `patch` refuses a hunk it finds
    # applied. Grepping for a string the patch adds is the check: present,
    # nothing to do; absent, apply for real and let a failure stop the
    # configure.
    FetchContent_Declare(
      "${_pkg_lc}"
      GIT_REPOSITORY "${url}"
      GIT_TAG "${_git_ref}"
      PATCH_COMMAND
        "${_sh_exe}" -c
        "grep -qrF '${_patch_marker}' src/ || '${_patch_exe}' -p1 -i '${_patch_file}'"
    )
  else()
    FetchContent_Declare(
      "${_pkg_lc}"
      GIT_REPOSITORY "${url}"
      GIT_TAG "${_git_ref}")
  endif()
  FetchContent_MakeAvailable("${_pkg_lc}")

  # Only create the imported interface if the target doesn't already exist This
  # allows the fetched package to create its own namespaced target
  if(NOT TARGET "${package}::${package}")
    add_library("${package}::${package}" INTERFACE IMPORTED)
    target_link_libraries("${package}::${package}" INTERFACE "${package}")
  endif()

  if(NOT EXISTS "${${_pkg_lc}_BINARY_DIR}/include")
    file(MAKE_DIRECTORY "${${_pkg_lc}_BINARY_DIR}/include")
  endif()

  unset(_pkg_lc)
  unset(_pkg_uc)
  unset(_git_ref)
  unset(_patch_file)
  unset(_patch_marker)

  # sanity check
  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "Could not find dependency ${package}")
  endif()
endmacro()
