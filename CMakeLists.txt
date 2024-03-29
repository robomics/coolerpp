# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.20)
cmake_policy(VERSION 3.20...3.26)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE
    OFF
    CACHE BOOL "Enable 'developer mode'")

include(FetchContent)
# cmake-format: off
FetchContent_Declare(
        _coolerpp_project_options
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.29.0.tar.xz
        URL_HASH SHA512=5cef5c4aab91328d41d0c846db326075eebab0bdcd53919d0bde288eb3ff4544ec51116234af25cb93475bb1ed806c51cef757868a56ceafdcf08a553afb6bf5
)
# cmake-format: on

FetchContent_MakeAvailable(_coolerpp_project_options)
include(${_coolerpp_project_options_SOURCE_DIR}/Index.cmake)

include(cmake/Versioning.cmake)
project(
  coolerpp
  LANGUAGES CXX
  VERSION ${COOLERPP_PROJECT_VERSION_MAJOR}.${COOLERPP_PROJECT_VERSION_MINOR}.${COOLERPP_PROJECT_VERSION_PATCH}
  HOMEPAGE_URL https://github.com/robomics/coolerpp
  DESCRIPTION "High-performance C++ library to read and write .cool files.")

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(CMAKE_CONFIGURATION_TYPES
        Debug Release RelWithDebInfo
        CACHE STRING "Enabled build types" FORCE)
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES
        ${CMAKE_BUILD_TYPE}
        CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

include(${_coolerpp_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake)

# dynamic_project_options sets recommended defaults and provides user and developer modes and full GUI support for
# choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode set(<feature_name>_USER_DEFAULT
# <value>) - set default for user mode

# Initialize project_options variable related to this project This overwrites `project_options` and sets
# `project_warnings` uncomment the options to enable them:

set(ENABLE_CACHE_DEFAULT ON)
set(ENABLE_COMPILE_COMMANDS_SYMLINK_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_USER OFF)
set(ENABLE_DOXYGEN_DEVELOPER ON)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_USER OFF)
set(ENABLE_SANITIZER_ADDRESS_DEVELOPER ON)
set(ENABLE_SANITIZER_LEAK_USER OFF)
set(ENABLE_SANITIZER_LEAK_DEVELOPER ON)
set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR_USER OFF)
set(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR_DEVELOPER ON)

dynamic_project_options(
  PREFIX
  coolerpp
  # CLANG_WARNINGS # Override the defaults for the CLANG warnings GCC_WARNINGS   # Override the defaults for the GCC
  CPPCHECK_OPTIONS
  --enable=performance,portability,style,warning
  --inline-suppr
  # We cannot act on a bug/missing feature of cppcheck
  --suppress=internalAstError
  # if a file does not have an internalAstError, we get an unmatchedSuppression error
  --suppress=unmatchedSuppression
  --suppress=passedByValue
  --inconclusive
  MSVC_WARNINGS
  ${MSVC_WARNINGS}
  CLANG_WARNINGS
  ${CLANG_WARNINGS}
  GCC_WARNINGS
  ${GCC_WARNINGS}
  CUDA_WARNINGS
  ${CUDA_WARNINGS})

target_compile_features(coolerpp_project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

# Tweak fmt
target_compile_definitions(coolerpp_project_options INTERFACE FMT_ENFORCE_COMPILE_STRING)

option(COOLERPP_ENABLE_TESTING "Build unit tests" ON)
option(COOLERPP_BUILD_EXAMPLES "Build examples" ON)
option(COOLERPP_BUILD_BENCHMARKS "Build benchmarks" OFF)

add_subdirectory(src)

option(COOLERPP_ENABLE_TESTING "Build unit tests" ON)
option(COOLERPP_BUILD_EXAMPLES "Build examples" ON)

if(COOLERPP_ENABLE_TESTING)
  enable_testing()
  message("-- Building unit tests.")
  target_compile_definitions(coolerpp_project_options INTERFACE COOLERPP_ENABLE_TESTING)
  add_subdirectory(test)
endif()

if(COOLERPP_BUILD_EXAMPLES)
  message("-- Building examples.")
  add_subdirectory(examples)
endif()

if(COOLERPP_BUILD_BENCHMARKS)
  message("Building benchmarks.")
  add_subdirectory(benchmarks)
endif()
