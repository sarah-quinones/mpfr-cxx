# automatically enable catch2 to generate ctest targets
include(../cmake/doctest.cmake)

find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_library(nanobench-main STATIC nanobench_main.cpp)
target_link_libraries(nanobench-main PRIVATE project_options)

add_executable(bench-pow2 pow2.cpp)
target_link_libraries(
  bench-pow2 PRIVATE ${Boost_LIBRARIES} nanobench-main mpfr-cxx
)

include_directories(../include)
