set(_lib "pic-mpi")
set(_url "https://github.com/JorgeG94/pic-mpi/")

# Pass PIC_ENABLE_MPI to the fetched package: OFF builds pic-mpi's serial
# backend, which is what makes a no-MPI rdb build link at all.
if(DEFINED PIC_ENABLE_MPI)
  set(PIC_ENABLE_MPI
      ${PIC_ENABLE_MPI}
      CACHE BOOL "Link pic-mpi against MPI" FORCE)
endif()

# Pass PIC_USE_LEGACY_MPI option to the fetched package if set
if(DEFINED PIC_USE_LEGACY_MPI)
  set(PIC_USE_LEGACY_MPI
      ${PIC_USE_LEGACY_MPI}
      CACHE BOOL "Use legacy MPI module" FORCE)
endif()

# Pass PIC_USE_VAPAA option to the fetched package if set (also set transitively
# from RDB_USE_VAPAA in the top-level CMakeLists).
if(DEFINED PIC_USE_VAPAA)
  set(PIC_USE_VAPAA
      ${PIC_USE_VAPAA}
      CACHE BOOL "Use vapaa to link to a C MPI library" FORCE)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Use the tagged release
set(_rev "main")
my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_url)
unset(_rev)
