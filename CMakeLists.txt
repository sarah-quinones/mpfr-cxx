cmake_minimum_required(VERSION 3.16)
project(mpfr-cxx CXX)

include(cmake/standard_project_settings.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link this 'library' to use the warnings specified in compiler_warnings.cmake
add_library(project_warnings INTERFACE)

# standard compiler warnings
include(cmake/compiler_warnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/sanitizers.cmake)
enable_sanitizers(project_options)

# enable doxygen
include(cmake/doxygen.cmake)
enable_doxygen()

# allow for static analysis options
include(cmake/static_analyzers.cmake)

option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
  target_precompile_headers(
    project_options
    INTERFACE
    <vector>
    <string>
    <map>
    <utility>
  )
endif()

# Conan dependencies
set(CONAN_REQUIRES
    # Boost Software License 1.0
    boost/1.73.0
    # MIT License
    fmt/6.2.1
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

if(ENABLE_TESTING)
  add_subdirectory(benchmark)
endif()
