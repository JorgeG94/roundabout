set(_lib "test-drive")
set(_pkg "test-drive")
set(_url "https://github.com/fortran-lang/test-drive")
# NOT v0.6.1, the newest tag: its quad-precision complex checks
# (`check_complex_xdp`) need `_FortranACAbsF128`, which the LLVM Flang runtime
# does not provide -- every test binary then fails to link with "undefined
# reference to _FortranACAbsF128". `main` is AHEAD of v0.6.1 and does not have
# the problem.
#
# So this pins the exact commit verified here at 174/174 on gfortran, ifx,
# nvfortran and flang. A SHA is a stronger pin than a tag anyway (tags can be
# moved); revisit when test-drive cuts a release newer than v0.6.1. (metalquicha
# pins v0.6.1 -- it does not build with Flang, so it never hits this.)
set(_rev "791484a79a24cbc00be6f74e53c6e030cb75efe2")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
