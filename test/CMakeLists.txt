# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 CONFIG REQUIRED)
find_package(Filesystem REQUIRED)

include(CTest)
include(Catch)

add_executable(test_main)

target_sources(
  test_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/attribute_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/bin_table_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/chromosome_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/cooler_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/dataset_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/index_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/pixel_selector_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/variant_buff_test.cpp)

target_include_directories(test_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/include/)

target_link_libraries(
  test_main
  PRIVATE project_warnings project_options
  PUBLIC Coolerpp::Coolerpp)

target_link_system_libraries(
  test_main
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  test_main
  TEST_SPEC
  "[short]"
  TEST_SUFFIX
  " - SHORT"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  EXTRA_ARGS
  --success)

catch_discover_tests(
  test_main
  TEST_SPEC
  "[medium]"
  TEST_SUFFIX
  " - MEDIUM"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  EXTRA_ARGS
  --success)

catch_discover_tests(
  test_main
  TEST_SPEC
  "[long]"
  TEST_SUFFIX
  " - LONG"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  REPORTER
  compact # This dramatically speeds up the test suit
  EXTRA_ARGS
  --success)
