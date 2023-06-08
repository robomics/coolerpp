# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(coolerpp_tools)

target_sources(
  coolerpp_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dump.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/load.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/merge.cpp)

target_link_libraries(
  coolerpp_tools
  PRIVATE coolerpp_project_warnings coolerpp_project_options
  PUBLIC CLI11::CLI11
         Coolerpp::Coolerpp
         fmt::fmt
         std::filesystem)

target_include_directories(coolerpp_tools PRIVATE include/)
