# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(FMT CONFIG REQUIRED)

set(example_targets "coolerpp_dump;coolerpp_load")

foreach(target_name ${example_targets})
  add_executable(${target_name} ${target_name}_example.cpp)
  target_link_libraries(
    ${target_name}
    PRIVATE coolerpp_project_warnings coolerpp_project_options
    PUBLIC Coolerpp::Coolerpp fmt::fmt std::filesystem)
endforeach()
