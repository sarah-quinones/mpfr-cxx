find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_library(nanobench-main STATIC nanobench_main.cpp)
target_link_libraries(
  nanobench-main PUBLIC project_options mpfr-cxx ${Boost_LIBRARIES}
)

add_executable(bench-pow2 pow2.cpp)
target_link_libraries(bench-pow2 PRIVATE nanobench-main)

add_executable(bench-std-vec std_vector.cpp)
target_link_libraries(bench-std-vec PRIVATE nanobench-main)

add_executable(bench-operations operations.cpp)
target_link_libraries(bench-operations PRIVATE nanobench-main)

include_directories(../include)
