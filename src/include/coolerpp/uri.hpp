// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <string_view>
#include <utility>

namespace coolerpp {

struct CoolerURI {
  std::string file_path;
  std::string group_path;

  CoolerURI() = default;
  CoolerURI(std::string_view p2, std::string_view p1);
  CoolerURI(std::string p1, std::string p2);
  CoolerURI(std::pair<std::string_view, std::string_view> paths);
  CoolerURI(std::pair<std::string, std::string> paths);
};

[[nodiscard]] CoolerURI parse_cooler_uri(std::string_view uri);
}  // namespace coolerpp
