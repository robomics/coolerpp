# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat CONFIG REQUIRED)
find_package(Filesystem REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)
find_package(tsl-hopscotch-map CONFIG REQUIRED)
find_package(tsl-ordered-map CONFIG REQUIRED)

add_library(coolerpp)
add_library(Coolerpp::Coolerpp ALIAS coolerpp)
target_sources(
  coolerpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/attribute_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bin_table_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/chromosome_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_accessors_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_read_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_standard_attr_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_validation_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_write_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_accessors_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_iterator_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_read_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dataset_write_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/index.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/pixel.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/pixel_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/uri.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp
          ${CMAKE_CURRENT_SOURCE_DIR}/variant_buff_impl.hpp)

target_include_directories(coolerpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)

# cmake-format: off
target_link_libraries(
  coolerpp
  PRIVATE coolerpp_project_warnings coolerpp_project_options
  PUBLIC fmt::fmt # For some reason linking fmt as sys library raises a
                  # "Found relative path while evaluating include directories"
                  # error
)
# cmake-format: on

target_link_system_libraries(
  coolerpp
  PUBLIC
  FastFloat::fast_float
  HighFive::HighFive
  tsl::hopscotch_map
  tsl::ordered_map
  std::filesystem)

install(TARGETS coolerpp INCLUDES)
