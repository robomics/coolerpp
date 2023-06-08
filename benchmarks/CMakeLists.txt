# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 CONFIG REQUIRED)
find_package(Filesystem REQUIRED)

include(CTest)
include(Catch)

add_executable(coolerpp_benchmark)

target_sources(
  coolerpp_benchmark
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/micro/bin_table_bench.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/micro/chromosome_bench.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/micro/iteration_bench.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/micro/pixel_bench.cpp)

target_include_directories(coolerpp_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/micro/include/)

target_link_libraries(
  coolerpp_benchmark
  PRIVATE coolerpp_project_warnings coolerpp_project_options
  PUBLIC Coolerpp::Coolerpp)

target_link_system_libraries(
  coolerpp_benchmark
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem)
