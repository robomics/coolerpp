# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat CONFIG REQUIRED)
find_package(Filesystem REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)
find_package(tsl-hopscotch-map CONFIG REQUIRED)
find_package(tsl-ordered-map CONFIG REQUIRED)

# Cooler
add_library(coolerpp)
add_library(Coolerpp::Coolerpp ALIAS coolerpp)
target_sources(
  coolerpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/attribute_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bin_table.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bin_table_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/chromosome.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/index.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/uri.cpp)

target_include_directories(coolerpp PUBLIC include/)

target_link_libraries(coolerpp PRIVATE project_warnings project_options)

target_link_system_libraries(
  coolerpp
  PUBLIC
  FastFloat::fast_float
  fmt::fmt
  HighFive::HighFive
  tsl::hopscotch_map
  tsl::ordered_map
  std::filesystem)

install(TARGETS coolerpp INCLUDES)
