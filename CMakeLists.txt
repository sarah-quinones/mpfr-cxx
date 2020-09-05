cmake_minimum_required(VERSION 3.16)
project(mpfr-cxx CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(ENABLE_DOCS_ONLY)
  add_subdirectory("docs")
  return()
endif()

include(cmake/standard_project_settings.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link this 'library' to use the warnings specified in compiler_warnings.cmake
add_library(project_warnings INTERFACE)

include(cmake/compiler_warnings.cmake)
set_project_warnings(project_warnings)

include(cmake/sanitizers.cmake)
enable_sanitizers(project_options)

include(cmake/static_analyzers.cmake)

# Conan dependencies
set(CONAN_REQUIRES
    # Boost Software License 1.0
    boost/1.74.0
    # MIT License
    fmt/7.0.3
    backward-cpp/1.5
    # LGPL v3.0
    mpfr/4.0.2
)
target_compile_definitions(
  project_options INTERFACE SPDLOG_FMT_EXTERNAL gsl_CONFIG_DEFAULTS_VERSION=1
)
set(CONAN_OPTIONS)

include(cmake/conan.cmake)
run_conan()

include_directories(include)

add_library(mpfr-cxx INTERFACE)
target_link_libraries(
  mpfr-cxx
  INTERFACE project_options
            project_warnings
            CONAN_PKG::boost
            CONAN_PKG::fmt
            CONAN_PKG::backward-cpp
            CONAN_PKG::mpfr
)

if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

if(ENABLE_BENCHMARK)
  add_subdirectory(benchmark)
endif()

if(ENABLE_DOCS)
  add_subdirectory ("docs")
endif()
