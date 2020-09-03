# automatically enable doctest to generate ctest targets
include(../cmake/doctest.cmake)

add_library(doctest-main STATIC doctest_main.cpp)
target_link_libraries(doctest-main PUBLIC project_options mpfr-cxx)

add_executable(test-mpfr-fn mpfr_fn.cpp)
target_link_libraries(test-mpfr-fn PUBLIC doctest-main)

include_directories(../include)

doctest_discover_tests(test-mpfr-fn)
