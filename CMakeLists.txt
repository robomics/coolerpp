# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE
    OFF
    CACHE BOOL "Enable 'developer mode'")

include(FetchContent)
# cmake-format: off
FetchContent_Declare(
        _project_options
        URL https://github.com/cpp-best-practices/project_options/archive/refs/tags/v0.24.0.tar.gz
        URL_HASH SHA512=f3afb9506724a2235a9f23243a78e73ce82017aaa84ba892c5d83e0a4f26c0771963c5e73c1190d8bfbf79a4d35539738e096e6182e19c7a78f6b44598ee9155
)
# cmake-format: on

FetchContent_MakeAvailable(_project_options)
include(${_project_options_SOURCE_DIR}/Index.cmake)

project(
  coolerpp
  LANGUAGES CXX
  VERSION 0.0.1
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

include(${_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake)
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
set(ENABLE_CONAN_DEFAULT ON)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_USER OFF)
set(ENABLE_DOXYGEN_DEVELOPER ON)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

dynamic_project_options(
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

target_compile_features(project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

# Tweak fmt
target_compile_definitions(project_options INTERFACE FMT_ENFORCE_COMPILE_STRING)

add_subdirectory(src)

if(ENABLE_TESTING)
  enable_testing()
  message("Building unit tests.")
  target_compile_definitions(project_options INTERFACE ENABLE_TESTING)
  add_subdirectory(test)
endif()

option(BUILD_EXAMPLES "Build examples" ON)
if(BUILD_EXAMPLES)
  message("Building examples.")
  add_subdirectory(examples)
endif()
