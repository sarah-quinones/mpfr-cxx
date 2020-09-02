# automatically enable catch2 to generate ctest targets
include(../cmake/doctest.cmake)

add_library(doctest-main STATIC doctest_main.cpp)
target_link_libraries(doctest-main PRIVATE project_options)

include_directories(../include)
