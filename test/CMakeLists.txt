# automatically enable doctest to generate ctest targets
include(../cmake/doctest.cmake)

add_library(backward OBJECT backward.cpp)
target_link_libraries(backward PUBLIC CONAN_PKG::backward-cpp)
add_library(doctest_main STATIC doctest_main.cpp)
target_link_libraries(doctest_main PUBLIC project_options mpfr-cxx)

set(testlibs backward doctest_main)

add_executable(test_mpfr_layer mpfr_layer.cpp)
target_link_libraries(test_mpfr_layer PUBLIC ${testlibs})

add_executable(test_math math.cpp)
target_link_libraries(test_math PUBLIC ${testlibs})

include_directories(../include)

doctest_discover_tests(test_mpfr_layer)
doctest_discover_tests(test_math)
