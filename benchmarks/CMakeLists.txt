# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 CONFIG REQUIRED)
find_package(Filesystem REQUIRED)

include(CTest)
include(Catch)

add_executable(coolerpp_benchmark_main)

target_sources(
  coolerpp_benchmark_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/micro/iteration_bench.cpp)

target_include_directories(coolerpp_benchmark_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/micro/include/)

target_link_libraries(
  coolerpp_benchmark_main
  PRIVATE coolerpp_project_warnings coolerpp_project_options
  PUBLIC Coolerpp::Coolerpp)

target_link_system_libraries(
  coolerpp_benchmark_main
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  coolerpp_benchmark_main
  TEST_SPEC
  "[bench]"
  TEST_SUFFIX
  " - BENCH"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Benchmarking/
  EXTRA_ARGS
  --success)
